use std::{iter::FromIterator, ops::{Range, RangeBounds}};

use mzpeaks::MZPeakSetType;
use timsrust::{converters::ConvertableDomain, Metadata};

use crate::{
    mzpeaks::{CentroidPeak, PeakSet},
    params::Unit,
    prelude::*,
    spectrum::{
        bindata::{ArrayRetrievalError, BinaryArrayMap3D},
        ArrayType, BinaryArrayMap, BinaryDataArrayType, DataArray,
    },
};

use mzsignal::feature_mapping::IMMSMapExtracter;


#[derive(Debug, Clone, Copy, PartialEq)]
struct NonNaNF64(f64);

impl From<f64> for NonNaNF64 {
    fn from(value: f64) -> Self {
        Self::wrap(value)
            .unwrap_or_else(|| panic!("Expected an order-able f64 value, but found {}", value))
    }
}

impl NonNaNF64 {
    fn wrap(value: f64) -> Option<Self> {
        if value.is_nan() {
            None
        } else {
            Some(Self(value))
        }
    }
}

impl PartialOrd for NonNaNF64 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for NonNaNF64 {}

impl Ord for NonNaNF64 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}

#[allow(unused)]
pub(crate) struct CoSorter2D<'a, Y> {
    dim1: &'a mut Vec<f64>,
    dim2: &'a mut Vec<Y>,
}

#[allow(unused)]
impl<'a, Y> CoSorter2D<'a, Y> {
    pub fn new(dim1: &'a mut Vec<f64>, dim2: &'a mut Vec<Y>) -> Self {
        Self { dim1, dim2 }
    }

    pub fn sort(self) {
        let n = self.dim1.len();
        let mut mask: Vec<_> = (0..n).into_iter().collect();
        mask.sort_by_key(|i| NonNaNF64(self.dim1[*i]));

        const TOMBSTONE: usize = usize::MAX;

        for idx in 0..n {
            if mask[idx] != TOMBSTONE {
                let mut current_idx = idx;
                loop {
                    let next_idx = mask[current_idx];
                    mask[current_idx] = TOMBSTONE;
                    if mask[next_idx] == TOMBSTONE {
                        break;
                    }
                    self.dim1.swap(current_idx, next_idx);
                    self.dim2.swap(current_idx, next_idx);
                    current_idx = next_idx;
                }
            }
        }
    }
}

pub struct FrameToArraysMapper<'a> {
    frame: &'a timsrust::Frame,
    metadata: &'a Metadata,
}

impl<'a> FrameToArraysMapper<'a> {
    pub fn new(frame: &'a timsrust::Frame, metadata: &'a Metadata) -> Self {
        Self { frame, metadata }
    }

    pub fn process_3d_slice(&self, iv: impl RangeBounds<usize>) -> BinaryArrayMap3D {
        let n_scans = self.frame.scan_offsets.len();

        let first_scan = match iv.start_bound() {
            std::ops::Bound::Included(i) => *i,
            std::ops::Bound::Excluded(i) => *i,
            std::ops::Bound::Unbounded => 0,
        };

        let final_scan = match iv.end_bound() {
            std::ops::Bound::Included(i) => *i,
            std::ops::Bound::Excluded(i) => *i + 1,
            std::ops::Bound::Unbounded => n_scans,
        }
        .min(n_scans);

        let mut im_dimension = Vec::with_capacity(final_scan - first_scan + 1);
        let mut arrays = Vec::with_capacity(final_scan - first_scan + 1);

        let mut scan_begin = first_scan;
        for (i, scan_end) in self.frame.scan_offsets[first_scan..final_scan]
            .iter()
            .copied()
            .enumerate()
        {
            let width = scan_end.saturating_sub(scan_begin);

            let mut mz_array_bytes: Vec<u8> =
                Vec::with_capacity(width * BinaryDataArrayType::Float64.size_of());
            let mut intensity_array_bytes: Vec<u8> =
                Vec::with_capacity(width * BinaryDataArrayType::Float32.size_of());

            self.frame.tof_indices[scan_begin..scan_begin + width]
                .iter()
                .for_each(|tof_idx| {
                    mz_array_bytes.extend_from_slice(
                        &self.metadata.mz_converter.convert(*tof_idx).to_le_bytes(),
                    )
                });
            (scan_begin..(scan_begin + width))
                .into_iter()
                .for_each(|idx| {
                    intensity_array_bytes.extend_from_slice(
                        &((self.frame.intensities[idx] as u64) as f32).to_le_bytes(),
                    );
                });
            let drift = self.metadata.im_converter.convert((i + first_scan) as u32) as f64;
            im_dimension.push(drift);

            let mz_array = DataArray::wrap(
                &ArrayType::MZArray,
                BinaryDataArrayType::Float64,
                mz_array_bytes,
            );
            let intensity_array = DataArray::wrap(
                &ArrayType::IntensityArray,
                BinaryDataArrayType::Float32,
                intensity_array_bytes,
            );

            let mut arrays_at = BinaryArrayMap::new();
            arrays_at.add(mz_array);
            arrays_at.add(intensity_array);
            arrays.push(arrays_at);
            scan_begin = scan_end;
        }

        // We read out the IM dimension in descending order, here we reverse it to be in ascending
        // IM order
        im_dimension.reverse();
        arrays.reverse();

        let arrays = BinaryArrayMap3D::from_ion_mobility_dimension_and_arrays(
            im_dimension,
            ArrayType::IonMobilityArray,
            Unit::VoltSecondPerSquareCentimeter,
            arrays,
        );
        arrays
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
        .into_iter().rev()
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
        return Ok(peaks.into_iter().next().unwrap().1.into_iter().map(|p| p.into()).collect());
    }

    let mut extracter = IMMSMapExtracter::from_iter(peaks);
    let features = extracter.extract_features(error_tolerance, 2, f64::INFINITY);

    let peaks: MZPeakSetType<CP> = features
        .iter()
        .map(|f| CentroidPeak::new(f.mz(), f.intensity(), 0).into())
        .collect();

    Ok(peaks)
}
