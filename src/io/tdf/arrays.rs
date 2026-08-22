use std::{
    iter::FromIterator,
    ops::{Range, RangeBounds},
};

use mzpeaks::{feature::Feature, IonMobility, MZPeakSetType, MZ};
use timsrust::{converters::ConvertableDomain, Metadata};

use crate::{
    io::tdf::{CalibrationParameters, sql::SQLFrame}, mzpeaks::{CentroidPeak, PeakSet}, params::Unit, prelude::*, spectrum::{
        ArrayType, BinaryArrayMap, BinaryDataArrayType, DataArray, bindata::{ArrayRetrievalError, BinaryArrayMap3D},
    },
};

use mzsignal::feature_mapping::{FeatureGraphBuilder, IMMSMapExtracter};

pub struct FrameToArraysMapper<'a> {
    frame: &'a timsrust::Frame,
    calibration_models: &'a CalibrationParameters,
    frame_meta: &'a SQLFrame
}

impl<'a> FrameToArraysMapper<'a> {
    pub fn new(frame: &'a timsrust::Frame, calibration_models: &'a CalibrationParameters, frame_meta: &'a SQLFrame) -> Self {
        Self { frame, calibration_models, frame_meta }
    }

    fn find_calibration_models(&self) -> (super::calibration::MzCalibrationModel, super::calibration::TimsCalibrationModel) {
        let mz_model = match self.calibration_models.find_mz_model_for_frame(self.frame_meta) {
            Ok(model) => model,
            Err(e) => {
                log::warn!("Falling back to basic m/z model: {e}");
                self.calibration_models.basic_mz_model.into()
            }
        };

        let im_model = match self.calibration_models.find_tims_model_for_frame(self.frame_meta) {
            Ok(model) => model,
            Err(e) => {
                log::warn!("Falling back to basic ion mobility model: {e}");
                self.calibration_models.basic_im_model.into()
            },
        };
        (mz_model, im_model)
    }

    pub fn process_3d_slice(&self, iv: impl RangeBounds<usize>) -> BinaryArrayMap3D {
        // `scan_offsets` is a cumulative prefix-sum of peak OFFSETS into
        // tof_indices/intensities (length `n_scan_rows + 1`, last entry = total
        // peak count); scan `k`'s peaks are tof_indices[scan_offsets[k] ..
        // scan_offsets[k + 1]].
        let n_scan_rows = self.frame.scan_offsets.len() - 1;

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

        let (mz_converter,
             im_converter) = self.find_calibration_models();

        for scan in first_scan..final_scan {
            let begin = self.frame.scan_offsets[scan];
            let end = self.frame.scan_offsets[scan + 1].min(self.frame.tof_indices.len());
            let width = end.saturating_sub(begin);

            let mut mz_array_bytes: Vec<u8> =
                Vec::with_capacity(width * BinaryDataArrayType::Float64.size_of());
            let mut intensity_array_bytes: Vec<u8> =
                Vec::with_capacity(width * BinaryDataArrayType::Float32.size_of());

            self.frame.tof_indices[begin..end].iter().for_each(|tof_idx| {
                mz_array_bytes
                    .extend_from_slice(&mz_converter.convert(*tof_idx).to_le_bytes())
            });
            self.frame.intensities[begin..end].iter().for_each(|intensity| {
                intensity_array_bytes.extend_from_slice(&((*intensity as u64) as f32).to_le_bytes())
            });
            let drift = im_converter.convert(scan as u32);
            im_dimension.push(drift);

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
    metadata: &Metadata,
    error_tolerance: Tolerance,
) -> Result<MZPeakSetType<CP>, ArrayRetrievalError> {
    let peaks: Result<Vec<_>, ArrayRetrievalError> = scan_range
        .clone()
        .rev()
        .map(|i| -> Result<(f64, PeakSet), ArrayRetrievalError> {
            let im = metadata.im_converter.convert(i);
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
