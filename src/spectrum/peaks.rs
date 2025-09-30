use std::borrow::Cow;
use std::iter::FusedIterator;

use mzpeaks::prelude::*;
use mzpeaks::{
    IndexType, MZPeakSetType, MassPeakSetType,
    peak::{CentroidPeak, DeconvolutedPeak, MZPoint},
    peak_set::PeakSetVec,
};

use super::BinaryArrayMap;
use super::bindata::ArrayRetrievalError;
use crate::prelude::BuildFromArrayMap;
use crate::spectrum::bindata::ArraysAvailable;
use crate::utils::mass_charge_ratio;

trait SummaryOps {
    /// Compute the base peak of a spectrum
    fn base_peak(&self) -> CentroidPeak;

    /// Find the minimum and maximum m/z values of a spectrum
    fn mz_range(&self) -> (f64, f64);

    /// Compute the total ion current for a spectrum
    fn tic(&self) -> f32;

    fn search(&self, query: f64, error_tolerance: Tolerance) -> Option<usize>;

    /// Find the number of points in a profile spectrum, or the number of peaks
    /// for a centroid spectrum
    fn len(&self) -> usize;

    /// Get the `i`th point in the peak data.
    ///
    /// **NOTE**: There is no guarantee that values are stored in m/z order.
    fn get(&self, i: usize) -> Option<MZPoint>;

    /// Compute summary information about the peak data's m/z range, base peak, total ion current,
    /// and number of data points
    fn fetch_summaries(&self) -> SpectrumSummary;

    #[allow(unused)]
    /// Check if the collection is empty
    fn is_empty(&self) -> bool;
}

impl SummaryOps for BinaryArrayMap {
    fn base_peak(&self) -> CentroidPeak {
        let mut peak = CentroidPeak::new(0.0, 0.0, 0);
        if let Ok(intensities) = self.intensities() {
            let result = intensities
                .iter()
                .enumerate()
                .max_by(|ia, ib| ia.1.total_cmp(ib.1));
            if let Ok(mzs) = self.mzs() {
                if let Some((i, inten)) = result {
                    peak = CentroidPeak::new(mzs[i], *inten, i as IndexType)
                }
            }
        }
        peak
    }

    fn mz_range(&self) -> (f64, f64) {
        if let Ok(mzs) = self.mzs() {
            if mzs.is_empty() {
                (0.0, 0.0)
            } else {
                (*mzs.first().unwrap(), *mzs.last().unwrap())
            }
        } else {
            (0.0, 0.0)
        }
    }

    fn tic(&self) -> f32 {
        if let Ok(intensities) = self.intensities() {
            intensities.iter().sum()
        } else {
            0.0
        }
    }

    fn search(&self, query: f64, error_tolerance: Tolerance) -> Option<usize> {
        self.search(query, error_tolerance)
    }

    fn len(&self) -> usize {
        self.mzs().map(|arr| arr.len()).unwrap_or_default()
    }

    fn get(&self, i: usize) -> Option<MZPoint> {
        match (self.mzs(), self.intensities()) {
            (Ok(mzs), Ok(intensities)) => match (mzs.get(i), intensities.get(i)) {
                (Some(mz), Some(intensity)) => Some(MZPoint::new(*mz, *intensity)),
                _ => None,
            },
            _ => None,
        }
    }

    fn fetch_summaries(&self) -> SpectrumSummary {
        let (mzs, intensities) = match (self.mzs(), self.intensities()) {
            (Ok(mzs), Ok(intensities)) => (mzs, intensities),
            (_, _) => (Cow::Owned(Vec::new()), Cow::Owned(Vec::new())),
        };

        let (tic, (bpmz, bpint, bpidx)) = mzs.iter().zip(intensities.iter()).enumerate().fold(
            (0.0, (0.0, 0.0f32, 0)),
            |(mut tic, (mut bpmz, mut bpint, mut bpidx)), (idx, (mz, int))| {
                tic += int;
                if *int > bpint {
                    bpint = *int;
                    bpmz = *mz;
                    bpidx = idx;
                }

                (tic, (bpmz, bpint, bpidx))
            },
        );
        let mz_range = if mzs.is_empty() {
            (0.0, 0.0)
        } else {
            (mzs.first().copied().unwrap(), mzs.last().copied().unwrap())
        };
        SpectrumSummary::new(
            tic,
            CentroidPeak::new(bpmz, bpint, bpidx as u32),
            mz_range,
            mzs.len(),
        )
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<C: CentroidLike> SummaryOps for MZPeakSetType<C> {
    fn base_peak(&self) -> CentroidPeak {
        let result = self
            .iter()
            .enumerate()
            .max_by(|ia, ib| ia.1.intensity().partial_cmp(&ib.1.intensity()).unwrap());
        if let Some((i, peak)) = result {
            CentroidPeak::new(peak.coordinate(), peak.intensity(), i as IndexType)
        } else {
            CentroidPeak::new(0.0, 0.0, 0)
        }
    }

    fn mz_range(&self) -> (f64, f64) {
        if PeakCollection::is_empty(self) {
            (0.0, 0.0)
        } else {
            (
                self.first().as_ref().unwrap().mz(),
                self.last().as_ref().unwrap().mz(),
            )
        }
    }

    fn tic(&self) -> f32 {
        self.iter().map(|p| p.intensity()).sum()
    }

    fn search(&self, query: f64, error_tolerance: Tolerance) -> Option<usize> {
        PeakCollection::search(self, query, error_tolerance)
    }

    fn len(&self) -> usize {
        PeakCollection::len(self)
    }

    fn get(&self, i: usize) -> Option<MZPoint> {
        self.peaks
            .get(i)
            .map(|p| MZPoint::new(p.mz(), p.intensity()))
    }

    fn fetch_summaries(&self) -> SpectrumSummary {
        let state = self.iter().fold(
            (
                0.0f32,
                CentroidPeak::default(),
                (f64::INFINITY, f64::NEG_INFINITY),
            ),
            |(mut tic, mut bp, (mut mz_min, mut mz_max)), p| {
                tic += p.intensity();
                if p.intensity() > bp.intensity {
                    bp = p.as_centroid();
                }
                let mz = p.mz();
                mz_min = mz_min.min(mz);
                mz_max = mz_max.max(mz);
                (tic, bp, (mz_min, mz_max))
            },
        );

        SpectrumSummary::new(state.0, state.1, state.2, SummaryOps::len(self))
    }

    fn is_empty(&self) -> bool {
        PeakCollection::is_empty(self)
    }
}

impl<D: DeconvolutedCentroidLike> SummaryOps for MassPeakSetType<D> {
    fn base_peak(&self) -> CentroidPeak {
        let result = self
            .iter()
            .enumerate()
            .max_by(|ia, ib| ia.1.intensity().partial_cmp(&ib.1.intensity()).unwrap());
        if let Some((i, peak)) = result {
            CentroidPeak::new(
                mass_charge_ratio(peak.coordinate(), peak.charge()),
                peak.intensity(),
                i as IndexType,
            )
        } else {
            CentroidPeak::new(0.0, 0.0, 0)
        }
    }

    fn mz_range(&self) -> (f64, f64) {
        if PeakCollection::is_empty(self) {
            (0.0, 0.0)
        } else {
            self.iter()
                .map(|p| {
                    let m = p.neutral_mass();
                    let z = p.charge();
                    mass_charge_ratio(m, z)
                })
                .fold((f64::MAX, f64::MIN), |state, mz| {
                    (state.0.min(mz), state.1.max(mz))
                })
        }
    }

    fn tic(&self) -> f32 {
        self.iter().map(|p| p.intensity()).sum()
    }

    fn search(&self, query: f64, error_tolerance: Tolerance) -> Option<usize> {
        PeakCollection::search(self, query, error_tolerance)
    }

    fn len(&self) -> usize {
        PeakCollection::len(self)
    }

    fn get(&self, i: usize) -> Option<MZPoint> {
        self.peaks.get(i).map(|p| {
            MZPoint::new(
                mass_charge_ratio(p.neutral_mass(), p.charge()),
                p.intensity(),
            )
        })
    }

    fn fetch_summaries(&self) -> SpectrumSummary {
        let state = self.iter().fold(
            (
                0.0f32,
                CentroidPeak::default(),
                (f64::INFINITY, f64::NEG_INFINITY),
            ),
            |(mut tic, mut bp, (mut mz_min, mut mz_max)), p| {
                tic += p.intensity();
                let mz = mass_charge_ratio(p.neutral_mass(), p.charge());
                if p.intensity() > bp.intensity {
                    bp.intensity = p.intensity();
                    bp.index = p.get_index();
                    bp.mz = mz;
                }
                mz_min = mz_min.min(mz);
                mz_max = mz_max.max(mz);
                (tic, bp, (mz_min, mz_max))
            },
        );

        SpectrumSummary::new(state.0, state.1, state.2, PeakCollection::len(self))
    }

    fn is_empty(&self) -> bool {
        PeakCollection::is_empty(self)
    }
}

/// Represent an owned representation of one the kinds of peak data that a [`SpectrumLike`](crate::spectrum::SpectrumLike) instance
/// might otherwise carry.
///
/// # See also
/// [`RefPeakDataLevel`] for borrowed version.
#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum PeakDataLevel<
    C: CentroidLike = CentroidPeak,
    D: DeconvolutedCentroidLike = DeconvolutedPeak,
> {
    Missing,
    RawData(BinaryArrayMap),
    Centroid(MZPeakSetType<C>),
    Deconvoluted(MassPeakSetType<D>),
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> PeakDataLevel<C, D> {
    /// Iterate over [`MZPoint`](mzpeaks::peak::MZPoint) data encoded in the peak data.
    ///
    /// **NOTE**: Values are produced in the order they are stored, so the data are not guaranteed
    /// to be ordered by m/z.
    pub fn iter(&self) -> PeakDataIterDispatch<'_, C, D> {
        match self {
            PeakDataLevel::Missing => PeakDataIterDispatch::PeakData(PeakDataIter::new(self)),
            PeakDataLevel::RawData(a) => {
                PeakDataIterDispatch::RawData(RawIter::from_binary_array_map(a).unwrap())
            }
            PeakDataLevel::Centroid(_) => PeakDataIterDispatch::PeakData(PeakDataIter::new(self)),
            PeakDataLevel::Deconvoluted(_) => {
                PeakDataIterDispatch::PeakData(PeakDataIter::new(self))
            }
        }
    }

    /// Attempt to reconstruct one of the peak levels based upon the available data arrays
    pub(crate) fn try_build_peaks(arrays: &BinaryArrayMap) -> Result<Self, ArrayRetrievalError>
    where
        C: BuildFromArrayMap,
        D: BuildFromArrayMap,
    {
        {
            if let ArraysAvailable::Ok = D::has_arrays_for(arrays) {
                let peaks = D::try_from_arrays(arrays)?.into();
                Ok(Self::Deconvoluted(peaks))
            } else if let ArraysAvailable::Ok = C::has_arrays_for(arrays) {
                let peaks = C::try_from_arrays(arrays)?.into();
                Ok(Self::Centroid(peaks))
            } else {
                Ok(Self::Missing)
            }
        }
    }
}

impl<C: CentroidLike + BuildFromArrayMap, D: DeconvolutedCentroidLike + BuildFromArrayMap>
    TryFrom<&BinaryArrayMap> for PeakDataLevel<C, D>
{
    type Error = ArrayRetrievalError;

    fn try_from(value: &BinaryArrayMap) -> Result<Self, Self::Error> {
        Self::try_build_peaks(value)
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> PeakDataLevel<C, D> {
    /// Compute the base peak of a spectrum
    pub fn base_peak(&self) -> CentroidPeak {
        match self {
            Self::Missing => CentroidPeak::new(0.0, 0.0, 0),
            Self::RawData(arrays) => arrays.base_peak(),
            Self::Centroid(peaks) => SummaryOps::base_peak(peaks),
            Self::Deconvoluted(peaks) => SummaryOps::base_peak(peaks),
        }
    }

    /// Find the minimum and maximum m/z values of a spectrum
    pub fn mz_range(&self) -> (f64, f64) {
        match self {
            Self::Missing => (0.0, 0.0),
            Self::RawData(arrays) => arrays.mz_range(),
            Self::Centroid(peaks) => peaks.mz_range(),
            Self::Deconvoluted(peaks) => peaks.mz_range(),
        }
    }

    /// Compute the total ion current for a spectrum
    pub fn tic(&self) -> f32 {
        match self {
            Self::Missing => 0.0,
            Self::RawData(arrays) => arrays.tic(),
            Self::Centroid(peaks) => peaks.tic(),
            Self::Deconvoluted(peaks) => peaks.tic(),
        }
    }

    /// Iterate over [`MZPoint`](mzpeaks::peak::MZPoint) data encoded in the peak data.
    ///
    /// **NOTE**: Values are produced in the order they are stored, so the data are not guaranteed
    /// to be ordered by m/z.
    pub fn search(&self, query: f64, error_tolerance: Tolerance) -> Option<usize> {
        match self {
            Self::Missing => None,
            Self::RawData(arrays) => arrays.search(query, error_tolerance),
            Self::Centroid(peaks) => SummaryOps::search(peaks, query, error_tolerance),
            Self::Deconvoluted(peaks) => SummaryOps::search(peaks, query, error_tolerance),
        }
    }

    /// Find the number of points in a profile spectrum, or the number of peaks
    /// for a centroid spectrum
    pub fn len(&self) -> usize {
        match self {
            Self::Missing => 0,
            Self::RawData(arrays) => SummaryOps::len(arrays),
            Self::Centroid(peaks) => SummaryOps::len(peaks),
            Self::Deconvoluted(peaks) => SummaryOps::len(peaks),
        }
    }

    /// Get the `i`th point in the peak data.
    ///
    /// **NOTE**: There is no guarantee that values are stored in m/z order.
    pub fn get(&self, i: usize) -> Option<MZPoint> {
        match self {
            PeakDataLevel::Missing => None,
            PeakDataLevel::RawData(arr) => SummaryOps::get(arr, i),
            PeakDataLevel::Centroid(peaks) => peaks.get(i),
            PeakDataLevel::Deconvoluted(peaks) => peaks.get(i),
        }
    }

    /// Compute summary information about the peak data's m/z range, base peak, total ion current,
    /// and number of data points
    pub fn fetch_summaries(&self) -> SpectrumSummary {
        match self {
            Self::Missing => SpectrumSummary::new(0.0, CentroidPeak::default(), (0.0, 0.0), 0),
            Self::RawData(arrays) => arrays.fetch_summaries(),
            Self::Centroid(peaks) => peaks.fetch_summaries(),
            Self::Deconvoluted(peaks) => peaks.fetch_summaries(),
        }
    }

    /// Check if the collection is empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// An [`MZPoint`] iterator over paired, potentially borrowed, arrays of m/z and intensity data that may not have been
/// extracted into discrete peaks.
///
/// This type should not need to be used directly in most cases. Instead see [`PeakDataIterDispatch`]
pub struct RawIter<'a> {
    mz_array: Cow<'a, [f64]>,
    intensity_array: Cow<'a, [f32]>,
    n: usize,
    i: usize,
}

impl<'a> RawIter<'a> {
    fn new(mz_array: Cow<'a, [f64]>, intensity_array: Cow<'a, [f32]>) -> Self {
        let n = mz_array.len();
        let i = 0;
        Self {
            mz_array,
            intensity_array,
            n,
            i,
        }
    }

    fn from_binary_array_map(arrays: &'a BinaryArrayMap) -> Result<Self, ArrayRetrievalError> {
        let mz_array = arrays.mzs()?;
        let intensity_array = arrays.intensities()?;

        Ok(Self::new(mz_array, intensity_array))
    }
}

impl Iterator for RawIter<'_> {
    type Item = MZPoint;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i >= self.n {
            None
        } else {
            let i = self.i;
            let mz = self.mz_array[i];
            let intens = self.intensity_array[i];
            self.i += 1;
            Some(MZPoint::new(mz, intens))
        }
    }
}

/// An [`MZPoint`] iterator over owned peak collections abstracted over by [`PeakDataLevel`]
///
/// This type should not need to be used directly in most cases. Instead see [`PeakDataIterDispatch`]
pub struct PeakDataIter<'a, C: CentroidLike, D: DeconvolutedCentroidLike> {
    peaks: &'a PeakDataLevel<C, D>,
    i: usize,
    n: usize,
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> Iterator for PeakDataIter<'_, C, D> {
    type Item = MZPoint;

    fn next(&mut self) -> Option<Self::Item> {
        let i = self.i;
        if i >= self.n {
            None
        } else {
            self.i += 1;
            self.get(i)
        }
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> ExactSizeIterator for PeakDataIter<'_, C, D> {
    fn len(&self) -> usize {
        self.n
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> FusedIterator for PeakDataIter<'_, C, D> {}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> DoubleEndedIterator for PeakDataIter<'_, C, D> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let n = self.n;
        if n < 1 {
            None
        } else {
            let pt = self.get(n - 1);
            self.n = self.n.saturating_sub(1);
            pt
        }
    }
}

impl<'a, C: CentroidLike, D: DeconvolutedCentroidLike> PeakDataIter<'a, C, D> {
    pub fn new(peaks: &'a PeakDataLevel<C, D>) -> Self {
        let n = peaks.len();
        Self { peaks, i: 0, n }
    }

    pub fn get(&self, i: usize) -> Option<MZPoint> {
        self.peaks.get(i)
    }
}

/// An [`MZPoint`] iterator over peak collections abstracted by source type and ownership status.
///
/// This iterator strategy has a slight amount of overhead on each operation as it dispatches to the
/// appropriate underlying type. It is necessary to allow [`PeakDataLevel`] and [`RefPeakDataLevel`]
/// to have a single type to return for an iterator.
pub enum PeakDataIterDispatch<'a, C: CentroidLike, D: DeconvolutedCentroidLike> {
    PeakData(PeakDataIter<'a, C, D>),
    RawData(RawIter<'a>),
    RefPeakData(RefPeakDataIter<'a, C, D>),
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> Iterator for PeakDataIterDispatch<'_, C, D> {
    type Item = MZPoint;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            PeakDataIterDispatch::PeakData(a) => a.next(),
            PeakDataIterDispatch::RawData(a) => a.next(),
            PeakDataIterDispatch::RefPeakData(a) => a.next(),
        }
    }
}

#[derive(Debug)]
/// An variant for dispatching to different strategies of computing
/// common statistics of different levels of peak data.
pub enum RefPeakDataLevel<'a, C: CentroidLike, D: DeconvolutedCentroidLike> {
    Missing,
    RawData(&'a BinaryArrayMap),
    Centroid(&'a MZPeakSetType<C>),
    Deconvoluted(&'a MassPeakSetType<D>),
}

/// A set of common summary metrics describing a mass spectrum
#[derive(Debug, Default, Clone)]
pub struct SpectrumSummary {
    /// The total ion current for a spectrum
    pub tic: f32,
    /// The base peak, the most intense peak in the spectrum
    pub base_peak: CentroidPeak,
    /// The minimum and maximum m/z observed
    pub mz_range: (f64, f64),
    /// The number of peaks or data points in the spectrum
    pub count: usize,
}

impl SpectrumSummary {
    pub fn new(tic: f32, base_peak: CentroidPeak, mz_range: (f64, f64), count: usize) -> Self {
        Self {
            tic,
            base_peak,
            mz_range,
            count,
        }
    }

    pub fn len(&self) -> usize {
        self.count
    }

    pub fn is_empty(&self) -> bool {
        self.count == 0
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> RefPeakDataLevel<'_, C, D> {
    /// Compute the base peak of a spectrum
    pub fn base_peak(&self) -> CentroidPeak {
        match self {
            RefPeakDataLevel::Missing => CentroidPeak::new(0.0, 0.0, 0),
            RefPeakDataLevel::RawData(arrays) => arrays.base_peak(),
            RefPeakDataLevel::Centroid(peaks) => SummaryOps::base_peak(*peaks),
            RefPeakDataLevel::Deconvoluted(peaks) => SummaryOps::base_peak(*peaks),
        }
    }

    /// Find the minimum and maximum m/z values of a spectrum
    pub fn mz_range(&self) -> (f64, f64) {
        match self {
            Self::Missing => (0.0, 0.0),
            Self::RawData(arrays) => arrays.mz_range(),
            Self::Centroid(peaks) => peaks.mz_range(),
            Self::Deconvoluted(peaks) => peaks.mz_range(),
        }
    }

    /// Compute the total ion current for a spectrum
    pub fn tic(&self) -> f32 {
        match self {
            Self::Missing => 0.0,
            Self::RawData(arrays) => arrays.tic(),
            Self::Centroid(peaks) => peaks.tic(),
            Self::Deconvoluted(peaks) => peaks.tic(),
        }
    }

    pub fn search(&self, query: f64, error_tolerance: Tolerance) -> Option<usize> {
        match self {
            Self::Missing => None,
            Self::RawData(arrays) => arrays.search(query, error_tolerance),
            Self::Centroid(peaks) => SummaryOps::search(*peaks, query, error_tolerance),
            Self::Deconvoluted(peaks) => SummaryOps::search(*peaks, query, error_tolerance),
        }
    }

    /// Find the number of points in a profile spectrum, or the number of peaks
    /// for a centroid spectrum
    pub fn len(&self) -> usize {
        match self {
            RefPeakDataLevel::Missing => 0,
            RefPeakDataLevel::RawData(arrays) => SummaryOps::len(*arrays),
            RefPeakDataLevel::Centroid(peaks) => SummaryOps::len(*peaks),
            RefPeakDataLevel::Deconvoluted(peaks) => SummaryOps::len(*peaks),
        }
    }

    /// Check if the collection is empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Iterate over [`MZPoint`](mzpeaks::peak::MZPoint) data encoded in the peak data.
    ///
    /// **NOTE**: Values are produced in the order they are stored, so the data are not guaranteed
    /// to be ordered by m/z.
    pub fn iter(&self) -> PeakDataIterDispatch<'_, C, D> {
        match self {
            RefPeakDataLevel::Missing => {
                PeakDataIterDispatch::RefPeakData(RefPeakDataIter::new(self))
            }
            RefPeakDataLevel::RawData(a) => {
                PeakDataIterDispatch::RawData(RawIter::from_binary_array_map(a).unwrap())
            }
            RefPeakDataLevel::Centroid(_) => {
                PeakDataIterDispatch::RefPeakData(RefPeakDataIter::new(self))
            }
            RefPeakDataLevel::Deconvoluted(_) => {
                PeakDataIterDispatch::RefPeakData(RefPeakDataIter::new(self))
            }
        }
    }

    /// Get the `i`th point in the peak data.
    ///
    /// **NOTE**: There is no guarantee that values are stored in m/z order.
    pub fn get(&self, i: usize) -> Option<MZPoint> {
        match self {
            RefPeakDataLevel::Missing => None,
            RefPeakDataLevel::RawData(arr) => SummaryOps::get(*arr, i),
            RefPeakDataLevel::Centroid(peaks) => peaks.get(i),
            RefPeakDataLevel::Deconvoluted(peaks) => peaks.get(i),
        }
    }

    /// Compute summary information about the peak data's m/z range, base peak, total ion current,
    /// and number of data points
    pub fn fetch_summaries(&self) -> SpectrumSummary {
        match self {
            RefPeakDataLevel::Missing => {
                SpectrumSummary::new(0.0, CentroidPeak::default(), (0.0, 0.0), 0)
            }
            RefPeakDataLevel::RawData(arrays) => arrays.fetch_summaries(),
            RefPeakDataLevel::Centroid(peaks) => peaks.fetch_summaries(),
            RefPeakDataLevel::Deconvoluted(peaks) => peaks.fetch_summaries(),
        }
    }
}

impl<C: CentroidLike + Clone, D: DeconvolutedCentroidLike + Clone> RefPeakDataLevel<'_, C, D> {
    pub fn cloned(&self) -> PeakDataLevel<C, D> {
        match self {
            RefPeakDataLevel::Missing => PeakDataLevel::Missing,
            RefPeakDataLevel::RawData(a) => PeakDataLevel::RawData(BinaryArrayMap::clone(a)),
            RefPeakDataLevel::Centroid(a) => PeakDataLevel::Centroid(PeakSetVec::clone(a)),
            RefPeakDataLevel::Deconvoluted(a) => PeakDataLevel::Deconvoluted(PeakSetVec::clone(a)),
        }
    }
}

/// An [`MZPoint`] iterator over borrowed peak collections abstracted over by [`RefPeakDataLevel`]
///
/// This type should not need to be used directly in most cases. Instead see [`PeakDataIterDispatch`]
pub struct RefPeakDataIter<'a, C: CentroidLike, D: DeconvolutedCentroidLike> {
    peaks: &'a RefPeakDataLevel<'a, C, D>,
    i: usize,
    n: usize,
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> Iterator for RefPeakDataIter<'_, C, D> {
    type Item = MZPoint;

    fn next(&mut self) -> Option<Self::Item> {
        let i = self.i;
        if i >= self.n {
            None
        } else {
            self.i += 1;
            self.get(i)
        }
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> ExactSizeIterator for RefPeakDataIter<'_, C, D> {
    fn len(&self) -> usize {
        self.n
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> FusedIterator for RefPeakDataIter<'_, C, D> {}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> DoubleEndedIterator
    for RefPeakDataIter<'_, C, D>
{
    fn next_back(&mut self) -> Option<Self::Item> {
        let n = self.n;
        if n < 1 {
            None
        } else {
            let pt = self.get(n - 1);
            self.n = self.n.saturating_sub(1);
            pt
        }
    }
}

impl<'a, C: CentroidLike, D: DeconvolutedCentroidLike> RefPeakDataIter<'a, C, D> {
    pub fn new(peaks: &'a RefPeakDataLevel<'a, C, D>) -> Self {
        let n = peaks.len();
        Self { peaks, i: 0, n }
    }

    pub fn get(&self, i: usize) -> Option<MZPoint> {
        self.peaks.get(i)
    }
}
