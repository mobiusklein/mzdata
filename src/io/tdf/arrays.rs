use std::{
    iter::FromIterator,
    ops::{Range, RangeBounds},
};

use mzpeaks::{feature::Feature, IonMobility, MZPeakSetType, MZ};
use timsrust::core::{Converter, Frame, Im, ScanIndex};

use crate::{
    mzpeaks::{CentroidPeak, PeakSet},
    params::Unit,
    prelude::*,
    spectrum::{
        bindata::{ArrayRetrievalError, BinaryArrayMap3D},
        ArrayType, BinaryArrayMap, BinaryDataArrayType, DataArray,
    },
};

use mzsignal::feature_mapping::{FeatureGraphBuilder, IMMSMapExtracter};

pub struct FrameToArraysMapper<'a> {
    frame: &'a Frame,
    handle: &'a timsrust::TimsTofPath
}

impl<'a> FrameToArraysMapper<'a> {
    pub fn new(frame: &'a timsrust::core::Frame, handle: &'a timsrust::TimsTofPath) -> Self {
        Self { frame, handle }
    }

    pub fn process_3d_slice(&self, iv: impl RangeBounds<usize>) -> BinaryArrayMap3D {
        let ions = self.frame.ions();
        let scan_offsets = ions.scan_offsets();

        // `scan_offsets` is a cumulative prefix-sum of peak OFFSETS into
        // tof_indices/intensities (length `n_scan_rows + 1`, last entry = total
        // peak count); scan `k`'s peaks are tof_indices[scan_offsets[k] ..
        // scan_offsets[k + 1]]. #50 corrected the read offset so peak counts are
        // right, but the mobility was still taken from the loop position
        // (`convert(i + first_scan)`), tagging each scan's peaks with the *next*
        // scan's mobility. Index scans directly so each is read once and labelled
        // with its own mobility. (The end-bound mapping is also restored to the
        // conventional half-open meaning; the selected peak set is unchanged.)

        let n_scan_rows = scan_offsets.len() - 1;

        let first_scan = match iv.start_bound() {
            std::ops::Bound::Included(i) => *i,
            std::ops::Bound::Excluded(i) => *i + 1,
            std::ops::Bound::Unbounded => 0,
        }
        .min(n_scan_rows);

        let final_scan = match iv.end_bound() {
            std::ops::Bound::Included(i) => *i + 1,
            std::ops::Bound::Excluded(i) => *i,
            std::ops::Bound::Unbounded => n_scan_rows,
        }
        .min(n_scan_rows);

        let span = final_scan.saturating_sub(first_scan);
        let mut im_dimension = Vec::with_capacity(span);
        let mut arrays = Vec::with_capacity(span);
        let tof_indices = ions.tof_indices();
        let intensities = ions.intensities();
        let mz_converter = self.handle.mz_converter().unwrap();
        let im_converter = self.handle.im_converter().unwrap();

        for scan in first_scan..final_scan {
            let begin = scan_offsets[scan];
            let end = scan_offsets[scan + 1].min(tof_indices.len());
            let width = end.saturating_sub(begin);

            let mut mz_array_bytes: Vec<u8> =
                Vec::with_capacity(width * BinaryDataArrayType::Float64.size_of());
            let mut intensity_array_bytes: Vec<u8> =
                Vec::with_capacity(width * BinaryDataArrayType::Float32.size_of());

            tof_indices[begin..end].iter().for_each(|tof_idx| {
                let v: f64 = mz_converter.convert(*tof_idx).into();
                mz_array_bytes
                    .extend_from_slice(&v.to_le_bytes())
            });
            intensities[begin..end].iter().for_each(|intensity| {
                intensity_array_bytes.extend_from_slice(&(u32::from(*intensity) as f32).to_le_bytes())
            });
            let scan_idx = ScanIndex::from(scan as u16);
            let drift: Im = im_converter.convert(scan_idx);
            im_dimension.push(drift.into());

            let mut mz_array = DataArray::wrap(
                &ArrayType::MZArray,
                BinaryDataArrayType::Float64,
                mz_array_bytes,
            );
            mz_array.unit = Unit::MZ;
            let mut intensity_array = DataArray::wrap(
                &ArrayType::IntensityArray,
                BinaryDataArrayType::Float32,
                intensity_array_bytes,
            );
            intensity_array.unit = Unit::DetectorCounts;

            let mut arrays_at = BinaryArrayMap::new();
            arrays_at.add(mz_array);
            arrays_at.add(intensity_array);
            arrays_at.sort_by_array(&ArrayType::MZArray).unwrap();
            arrays.push(arrays_at);
        }

        // We read out the IM dimension in descending order, here we reverse it to be in ascending
        // IM order
        im_dimension.reverse();
        arrays.reverse();

        BinaryArrayMap3D::from_ion_mobility_dimension_and_arrays(
            im_dimension,
            ArrayType::MeanInverseReducedIonMobilityArray,
            Unit::VoltSecondPerSquareCentimeter,
            arrays,
        )
    }
}

pub fn consolidate_peaks<CP: CentroidLike + From<CentroidPeak>>(
    arrays: &BinaryArrayMap3D,
    scan_range: &Range<u32>,
    handle: &timsrust::TimsTofPath,
    error_tolerance: Tolerance,
) -> Result<MZPeakSetType<CP>, ArrayRetrievalError> {
    let im_converter = handle.im_converter().unwrap();
    let peaks: Result<Vec<_>, ArrayRetrievalError> = scan_range
        .clone()
        .rev()
        .map(|i| -> Result<(f64, PeakSet), ArrayRetrievalError> {
            let idx = ScanIndex::from(i as u16);
            let im: Im = im_converter.convert(idx);
            let im = im.into();
            if let Some(arrays_point) = arrays.get_ion_mobility(im) {
                let mzs = arrays_point.mzs()?;
                let intens = arrays_point.intensities()?;
                let peaks: PeakSet = mzs
                    .iter()
                    .copied()
                    .zip(intens.iter().copied())
                    .map(|(mz, i)| CentroidPeak::new(mz, i, 0))
                    .collect();
                Ok((im, peaks))
            } else {
                Ok((im, PeakSet::empty()))
            }
        })
        .collect();

    let peaks = peaks?;
    if peaks.is_empty() {
        return Ok(MZPeakSetType::empty());
    }

    if peaks.len() == 1 {
        return Ok(peaks
            .into_iter()
            .next()
            .unwrap()
            .1
            .into_iter()
            .map(|p| p.into())
            .collect());
    }

    let mut extracter = IMMSMapExtracter::from_iter(peaks);
    let features = extracter.extract_features(error_tolerance, 2, 0.01);
    let merger = mzsignal::feature_mapping::FeatureMerger::<
        MZ,
        IonMobility,
        Feature<MZ, IonMobility>,
    >::default();
    let features = merger
        .bridge_feature_gaps(features, error_tolerance, f64::INFINITY)
        .features;

    let peaks: MZPeakSetType<CP> = features
        .iter()
        .map(|f| CentroidPeak::new(f.mz(), f.intensity(), 0).into())
        .collect();

    Ok(peaks)
}
