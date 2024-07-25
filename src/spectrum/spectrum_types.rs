use std::convert::TryFrom;
use std::{borrow::Cow, ops::Index};

use mzpeaks::Mass;
use thiserror::Error;

use mzpeaks::{
    peak_set::PeakSetVec, prelude::*, CentroidLike, CentroidPeak, DeconvolutedCentroidLike,
    DeconvolutedPeak, IndexType, MZPeakSetType, MassPeakSetType, PeakCollection, PeakSet,
    Tolerance, MZ,
};

#[cfg(feature = "mzsignal")]
use mzsignal::{
    denoise::{denoise, DenoisingError},
    peak_picker::{PeakFitType, PeakPicker, PeakPickerError},
    reprofile::{self, PeakShape, PeakShapeModel},
    FittedPeak,
};

use crate::params::{ParamDescribed, ParamList, Unit, Value};
#[allow(unused)]
use crate::spectrum::bindata::{ArrayType, BinaryArrayMap, BinaryDataArrayType};
use crate::spectrum::scan_properties::{
    Acquisition, IonMobilityMeasure, Precursor, ScanPolarity, SignalContinuity, SpectrumDescription,
};
use crate::utils::mass_charge_ratio;

use super::bindata::{ArrayRetrievalError, ArraysAvailable, BuildArrayMapFrom, BuildFromArrayMap};
#[allow(unused)]
use super::DataArray;

/// A blanket trait that ties together all the assumed behaviors of an m/z coordinate centroid peak
pub trait CentroidPeakAdapting: CentroidLike + Default + From<CentroidPeak> {}
impl<C: CentroidLike + Default + From<CentroidPeak>> CentroidPeakAdapting for C {}

/// A blanket trait that ties together all the assumed behaviors of an neutral mass coordinate centroid peak
pub trait DeconvolutedPeakAdapting:
    DeconvolutedCentroidLike + Default + From<DeconvolutedPeak>
{
}
impl<D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak>> DeconvolutedPeakAdapting
    for D
{
}

/// Represent an owned representation of one the kinds of peak data that a [`SpectrumLike`] instance
/// might otherwise carry.
///
/// # See also
/// [`RefPeakDataLevel`] for borrowed version.
#[derive(Debug)]
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
    /// Compute the base peak of a spectrum
    pub fn base_peak(&self) -> CentroidPeak {
        match self {
            Self::Missing => CentroidPeak::new(0.0, 0.0, 0),
            Self::RawData(arrays) => {
                let mut peak = CentroidPeak::new(0.0, 0.0, 0);
                if let Ok(intensities) = arrays.intensities() {
                    let result = intensities
                        .iter()
                        .enumerate()
                        .max_by(|ia, ib| ia.1.total_cmp(ib.1));
                    if let Ok(mzs) = arrays.mzs() {
                        if let Some((i, inten)) = result {
                            peak = CentroidPeak::new(mzs[i], *inten, i as IndexType)
                        }
                    }
                }
                peak
            }
            Self::Centroid(peaks) => {
                let result = peaks
                    .iter()
                    .enumerate()
                    .max_by(|ia, ib| ia.1.intensity().partial_cmp(&ib.1.intensity()).unwrap());
                if let Some((i, peak)) = result {
                    CentroidPeak::new(peak.coordinate(), peak.intensity(), i as IndexType)
                } else {
                    CentroidPeak::new(0.0, 0.0, 0)
                }
            }
            Self::Deconvoluted(peaks) => {
                let result = peaks
                    .iter()
                    .enumerate()
                    .max_by(|ia, ib| ia.1.intensity().partial_cmp(&ib.1.intensity()).unwrap());
                if let Some((i, peak)) = result {
                    CentroidPeak::new(
                        crate::utils::mass_charge_ratio(peak.coordinate(), peak.charge()),
                        peak.intensity(),
                        i as IndexType,
                    )
                } else {
                    CentroidPeak::new(0.0, 0.0, 0)
                }
            }
        }
    }

    /// Find the minimum and maximum m/z values of a spectrum
    pub fn mz_range(&self) -> (f64, f64) {
        match self {
            Self::Missing => (0.0, 0.0),
            Self::RawData(arrays) => {
                if let Ok(mzs) = arrays.mzs() {
                    if mzs.len() == 0 {
                        (0.0, 0.0)
                    } else {
                        (*mzs.first().unwrap(), *mzs.last().unwrap())
                    }
                } else {
                    (0.0, 0.0)
                }
            }
            Self::Centroid(peaks) => {
                if peaks.is_empty() {
                    (0.0, 0.0)
                } else {
                    (
                        peaks.first().as_ref().unwrap().mz(),
                        peaks.last().as_ref().unwrap().mz(),
                    )
                }
            }
            Self::Deconvoluted(peaks) => {
                if peaks.len() == 0 {
                    (0.0, 0.0)
                } else {
                    peaks
                        .iter()
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
        }
    }

    /// Compute the total ion current for a spectrum
    pub fn tic(&self) -> f32 {
        match self {
            Self::Missing => 0.0,
            Self::RawData(arrays) => {
                if let Ok(intensities) = arrays.intensities() {
                    intensities.iter().sum()
                } else {
                    0.0
                }
            }
            Self::Centroid(peaks) => peaks.iter().map(|p| p.intensity()).sum(),
            Self::Deconvoluted(peaks) => peaks.iter().map(|p| p.intensity()).sum(),
        }
    }

    pub fn search(&self, query: f64, error_tolerance: Tolerance) -> Option<usize> {
        match self {
            Self::Missing => None,
            Self::RawData(arrays) => arrays.search(query, error_tolerance),
            Self::Centroid(peaks) => peaks.search(query, error_tolerance),
            Self::Deconvoluted(peaks) => peaks.search(query, error_tolerance),
        }
    }

    /// Find the number of points in a profile spectrum, or the number of peaks
    /// for a centroid spectrum
    pub fn len(&self) -> usize {
        match self {
            Self::Missing => 0,
            Self::RawData(arrays) => arrays.mzs().map(|arr| arr.len()).unwrap_or_default(),
            Self::Centroid(peaks) => peaks.len(),
            Self::Deconvoluted(peaks) => peaks.len(),
        }
    }

    /// Check if the collection is empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
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

#[derive(Debug, Default, Clone)]
pub struct SpectrumSummary {
    pub tic: f32,
    pub bp: CentroidPeak,
    pub mz_range: (f64, f64),
    pub count: usize,
}

impl SpectrumSummary {
    pub fn new(tic: f32, bp: CentroidPeak, mz_range: (f64, f64), count: usize) -> Self {
        Self {
            tic,
            bp,
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

impl<'a, C: CentroidLike, D: DeconvolutedCentroidLike> RefPeakDataLevel<'a, C, D> {
    /// Compute the base peak of a spectrum
    pub fn base_peak(&self) -> CentroidPeak {
        match self {
            RefPeakDataLevel::Missing => CentroidPeak::new(0.0, 0.0, 0),
            RefPeakDataLevel::RawData(arrays) => {
                let mut peak = CentroidPeak::new(0.0, 0.0, 0);
                if let Ok(intensities) = arrays.intensities() {
                    let result = intensities
                        .iter()
                        .enumerate()
                        .max_by(|ia, ib| ia.1.total_cmp(ib.1));
                    if let Ok(mzs) = arrays.mzs() {
                        if let Some((i, inten)) = result {
                            peak = CentroidPeak::new(mzs[i], *inten, i as IndexType)
                        }
                    }
                }
                peak
            }
            RefPeakDataLevel::Centroid(peaks) => {
                let result = peaks
                    .iter()
                    .enumerate()
                    .max_by(|ia, ib| ia.1.intensity().partial_cmp(&ib.1.intensity()).unwrap());
                if let Some((i, peak)) = result {
                    CentroidPeak::new(peak.coordinate(), peak.intensity(), i as IndexType)
                } else {
                    CentroidPeak::new(0.0, 0.0, 0)
                }
            }
            RefPeakDataLevel::Deconvoluted(peaks) => {
                let result = peaks
                    .iter()
                    .enumerate()
                    .max_by(|ia, ib| ia.1.intensity().partial_cmp(&ib.1.intensity()).unwrap());
                if let Some((i, peak)) = result {
                    CentroidPeak::new(
                        crate::utils::mass_charge_ratio(peak.coordinate(), peak.charge()),
                        peak.intensity(),
                        i as IndexType,
                    )
                } else {
                    CentroidPeak::new(0.0, 0.0, 0)
                }
            }
        }
    }

    /// Find the minimum and maximum m/z values of a spectrum
    pub fn mz_range(&self) -> (f64, f64) {
        match self {
            Self::Missing => (0.0, 0.0),
            Self::RawData(arrays) => {
                if let Ok(mzs) = arrays.mzs() {
                    if mzs.len() == 0 {
                        (0.0, 0.0)
                    } else {
                        eprintln!("Size: {} | {:?} - {:?}", mzs.len(), mzs.first(), mzs.last());
                        (*mzs.first().unwrap(), *mzs.last().unwrap())
                    }
                } else {
                    (0.0, 0.0)
                }
            }
            Self::Centroid(peaks) => {
                if peaks.is_empty() {
                    (0.0, 0.0)
                } else {
                    (
                        peaks.first().as_ref().unwrap().mz(),
                        peaks.last().as_ref().unwrap().mz(),
                    )
                }
            }
            Self::Deconvoluted(peaks) => {
                if peaks.len() == 0 {
                    (0.0, 0.0)
                } else {
                    peaks
                        .iter()
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
        }
    }

    /// Compute the total ion current for a spectrum
    pub fn tic(&self) -> f32 {
        match self {
            Self::Missing => 0.0,
            Self::RawData(arrays) => {
                if let Ok(intensities) = arrays.intensities() {
                    intensities.iter().sum()
                } else {
                    0.0
                }
            }
            Self::Centroid(peaks) => peaks.iter().map(|p| p.intensity()).sum(),
            Self::Deconvoluted(peaks) => peaks.iter().map(|p| p.intensity()).sum(),
        }
    }

    pub fn search(&self, query: f64, error_tolerance: Tolerance) -> Option<usize> {
        match self {
            Self::Missing => None,
            Self::RawData(arrays) => arrays.search(query, error_tolerance),
            Self::Centroid(peaks) => peaks.search(query, error_tolerance),
            Self::Deconvoluted(peaks) => peaks.search(query, error_tolerance),
        }
    }

    /// Find the number of points in a profile spectrum, or the number of peaks
    /// for a centroid spectrum
    pub fn len(&self) -> usize {
        match self {
            RefPeakDataLevel::Missing => 0,
            RefPeakDataLevel::RawData(arrays) => {
                arrays.mzs().map(|arr| arr.len()).unwrap_or_default()
            }
            RefPeakDataLevel::Centroid(peaks) => peaks.len(),
            RefPeakDataLevel::Deconvoluted(peaks) => peaks.len(),
        }
    }

    /// Check if the collection is empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn fetch_summaries(&self) -> SpectrumSummary {
        match self {
            RefPeakDataLevel::Missing => {
                SpectrumSummary::new(0.0, CentroidPeak::default(), (0.0, 0.0), 0)
            }
            RefPeakDataLevel::RawData(arrays) => {
                let mzs = arrays.mzs().unwrap();
                let intensities = arrays.intensities().unwrap();
                let (tic, (bpmz, bpint, bpidx)) =
                    mzs.iter().zip(intensities.iter()).enumerate().fold(
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
            RefPeakDataLevel::Centroid(peaks) => {
                let state = peaks.iter().fold(
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

                SpectrumSummary::new(state.0, state.1, state.2, peaks.len())
            }
            RefPeakDataLevel::Deconvoluted(peaks) => {
                let state = peaks.iter().fold(
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

                SpectrumSummary::new(state.0, state.1, state.2, peaks.len())
            }
        }
    }
}

impl<'a, C: CentroidLike + Clone, D: DeconvolutedCentroidLike + Clone> RefPeakDataLevel<'a, C, D> {
    pub fn cloned(&self) -> PeakDataLevel<C, D> {
        match self {
            RefPeakDataLevel::Missing => PeakDataLevel::Missing,
            RefPeakDataLevel::RawData(a) => PeakDataLevel::RawData(BinaryArrayMap::clone(a)),
            RefPeakDataLevel::Centroid(a) => PeakDataLevel::Centroid(PeakSetVec::clone(a)),
            RefPeakDataLevel::Deconvoluted(a) => PeakDataLevel::Deconvoluted(PeakSetVec::clone(a)),
        }
    }
}

/// A trait for providing a uniform delegated access to spectrum metadata
pub trait SpectrumLike<
    C: CentroidLike = CentroidPeak,
    D: DeconvolutedCentroidLike = DeconvolutedPeak,
>
{
    /// The method to access the spectrum description itself, which supplies
    /// the data for most other methods on this trait.
    fn description(&self) -> &SpectrumDescription;

    /// The method to access the spectrum descript itself, mutably.
    fn description_mut(&mut self) -> &mut SpectrumDescription;

    /// Access the acquisition information for this spectrum.
    #[inline]
    fn acquisition(&self) -> &Acquisition {
        &self.description().acquisition
    }

    /// Access the precursor information, if it exists.
    #[inline]
    fn precursor(&self) -> Option<&Precursor> {
        let desc = self.description();
        if let Some(precursor) = &desc.precursor {
            Some(precursor)
        } else {
            None
        }
    }

    /// Iterate over all precursors of the spectrum
    fn precursor_iter(&self) -> impl Iterator<Item = &Precursor> {
        let desc = self.description();
        desc.precursor.iter()
    }

    /// Mutably access the precursor information, if it exists
    fn precursor_mut(&mut self) -> Option<&mut Precursor> {
        let desc = self.description_mut();
        if let Some(precursor) = desc.precursor.as_mut() {
            Some(precursor)
        } else {
            None
        }
    }

    /// Iterate over all precursors of the spectrum mutably
    fn precursor_iter_mut(&mut self) -> impl Iterator<Item = &mut Precursor> {
        let desc = self.description_mut();
        desc.precursor.iter_mut()
    }

    /// A shortcut method to retrieve the scan start time of a spectrum
    #[inline]
    fn start_time(&self) -> f64 {
        let acq = self.acquisition();
        match acq.scans.first() {
            Some(evt) => evt.start_time,
            None => 0.0,
        }
    }

    /// Access the MS exponentiation level
    #[inline]
    fn ms_level(&self) -> u8 {
        self.description().ms_level
    }

    /// Access the native ID string for the spectrum
    #[inline]
    fn id(&self) -> &str {
        &self.description().id
    }

    /// Access the index of the spectrum in the source file
    #[inline]
    fn index(&self) -> usize {
        self.description().index
    }

    /// Access a description of how raw the signal is, whether a
    /// profile spectrum is available or only centroids are present.
    #[inline]
    fn signal_continuity(&self) -> SignalContinuity {
        self.description().signal_continuity
    }

    #[inline]
    fn polarity(&self) -> ScanPolarity {
        self.description().polarity
    }

    #[inline]
    fn params(&self) -> &ParamList {
        &self.description().params
    }

    /// Access the point measure of ion mobility associated with the scan if present. This is distinct from
    /// having a frame-level scan across the ion mobility dimension.
    fn ion_mobility(&self) -> Option<f64> {
        self.acquisition()
            .iter()
            .flat_map(|s| s.ion_mobility())
            .next()
    }

    /// Check if this spectrum has a point measure of ion mobility. This is distinct from
    /// having a frame-level scan across the ion mobility dimension.
    fn has_ion_mobility(&self) -> bool {
        self.ion_mobility().is_some()
    }

    /// Retrieve the most processed representation of the mass spectrum's
    /// signal
    fn peaks(&'_ self) -> RefPeakDataLevel<'_, C, D>;

    fn into_peaks_and_description(self) -> (PeakDataLevel<C, D>, SpectrumDescription);

    fn raw_arrays(&'_ self) -> Option<&'_ BinaryArrayMap>;

    /// Compute and update the the total ion current, base peak, and m/z range for
    /// the spectrum based upon its current peak data.
    ///
    /// Uses [`RefPeakDataLevel::fetch_summaries`]
    fn update_summaries(&mut self) {
        let SpectrumSummary {
            tic,
            bp,
            mz_range,
            count: _,
        } = self.peaks().fetch_summaries();

        let desc = self.description_mut();
        let params = desc.params_mut();
        let tic_curie = curie!(MS:1000285);

        if let Some(p) = params.iter_mut().find(|p| **p == tic_curie) {
            p.value = Value::Float(tic as f64);
        } else {
            let mut p = tic_curie.as_param();
            p.name = "total ion current".to_string();
            p.value = Value::Float(tic as f64);
            p.unit = Unit::DetectorCounts;
            params.push(p)
        }

        let bpmz_curie = curie!(MS:1000504);
        if let Some(p) = params.iter_mut().find(|p| **p == bpmz_curie) {
            p.value = Value::Float(bp.mz());
        } else {
            let mut p = bpmz_curie.as_param();
            p.name = "base peak m/z".to_string();
            p.value = Value::Float(bp.mz());
            p.unit = Unit::MZ;
            params.push(p);
        }

        let lowest_mz_curie = curie!(MS:1000528);
        if let Some(p) = params.iter_mut().find(|p| **p == lowest_mz_curie) {
            p.value = Value::Float(mz_range.0);
        } else {
            let mut p = lowest_mz_curie.as_param();
            p.name = "lowest observed m/z".to_string();
            p.value = Value::Float(mz_range.0);
            p.unit = Unit::MZ;
            params.push(p);
        }

        let highest_mz_curie = curie!(MS:1000527);
        if let Some(p) = params.iter_mut().find(|p| **p == highest_mz_curie) {
            p.value = Value::Float(mz_range.0);
        } else {
            let mut p = highest_mz_curie.as_param();
            p.name = "highest observed m/z".to_string();
            p.value = Value::Float(mz_range.1);
            p.unit = Unit::MZ;
            params.push(p);
        }

        let bpint_curie = curie!(MS:1000505);
        if let Some(p) = params.iter_mut().find(|p| **p == bpint_curie) {
            p.value = Value::Float(bp.intensity() as f64);
        } else {
            let mut p = bpint_curie.as_param();
            p.name = "base peak intensity".to_string();
            p.value = Value::Float(bp.intensity() as f64);
            p.unit = Unit::DetectorCounts;
            params.push(p);
        }
    }
}

#[derive(Default, Debug, Clone)]
/// Represents a spectrum that hasn't been processed yet, with only
/// data arrays, potentially no discrete peaks. A raw spectrum may still
/// be centroided, but the peaks still need to be decoded.
pub struct RawSpectrum {
    /// The spectrum metadata describing acquisition conditions and details.
    pub description: SpectrumDescription,
    /// The data arrays describing the m/z, intensity, and potentially other
    /// measured properties
    pub arrays: BinaryArrayMap,
}


impl ParamDescribed for RawSpectrum {
    fn params(&self) -> &[crate::params::Param] {
        <SpectrumDescription as ParamDescribed>::params(&self.description)
    }

    fn params_mut(&mut self) -> &mut crate::params::ParamList {
        <SpectrumDescription as ParamDescribed>::params_mut(&mut self.description)
    }
}

/// Errors that may arise when converting between different [`SpectrumLike`] types
#[derive(Debug, Clone, PartialEq, Error)]
pub enum SpectrumConversionError {
    #[error("m/z array does not match size of intensity array")]
    MZIntensityArraySizeMismatch,
    #[error("Operation expected charge-deconvolved data but did not find it")]
    NotDeconvoluted,
    #[error("Operation expected centroided data but did not find it")]
    NotCentroided,
    #[error("No peak data of any kind was found")]
    NoPeakData,
    #[error("An error occurred while accessing raw data arrays: {0}")]
    ArrayRetrievalError(
        #[from]
        #[source]
        ArrayRetrievalError,
    ),
}

/// Errors that may arise when performing signal processing or other data transformation
#[derive(Debug, Clone, thiserror::Error)]
pub enum SpectrumProcessingError {
    #[cfg(feature = "mzsignal")]
    #[error("An error occurred while denoising: {0:?}")]
    DenoisingError(
        #[from]
        #[source]
        DenoisingError,
    ),

    #[cfg(feature = "mzsignal")]
    #[error("An error occurred while peak picking: {0:?}")]
    PeakPickerError(
        #[from]
        #[source]
        PeakPickerError,
    ),

    #[error("An error occurred while trying to convert spectrum types: {0}")]
    SpectrumConversionError(
        #[from]
        #[source]
        SpectrumConversionError,
    ),

    #[error("An error occurred while accessing raw data arrays: {0}")]
    ArrayRetrievalError(
        #[from]
        #[source]
        ArrayRetrievalError,
    ),
}

impl<'transient, 'lifespan: 'transient> RawSpectrum {
    pub fn new(description: SpectrumDescription, arrays: BinaryArrayMap) -> Self {
        Self {
            description,
            arrays,
        }
    }

    /// Convert a spectrum into a [`CentroidSpectrumType`].
    ///
    /// # Errors
    /// This operation returns [`SpectrumConversionError::NotCentroided`] if the
    /// [`SpectrumLike::signal_continuity`] != [`SignalContinuity::Centroid`].
    ///
    /// # See also
    /// To pick peaks from profile data, see [`RawSpectrum::pick_peaks_into`]
    pub fn into_centroid<C: CentroidLike + Default>(
        self,
    ) -> Result<CentroidSpectrumType<C>, SpectrumConversionError>
    where
        C: BuildFromArrayMap,
    {
        if !matches!(
            self.description.signal_continuity,
            SignalContinuity::Centroid
        ) {
            return Err(SpectrumConversionError::NotCentroided);
        }

        let peaks = match C::try_from_arrays(&self.arrays) {
            Ok(peaks) => peaks.into(),
            Err(e) => return Err(e.into()),
        };
        let mut centroid = CentroidSpectrumType::<C> {
            description: self.description,
            peaks,
        };
        centroid.description.signal_continuity = SignalContinuity::Centroid;

        Ok(centroid)
    }

    /// Access the m/z array.
    ///
    /// # Panics
    /// This function will panic if [`ArrayType::MZArray`] is not
    /// present in [`RawSpectrum::arrays`]
    pub fn mzs(&'lifespan self) -> Cow<'transient, [f64]> {
        self.arrays.mzs().unwrap()
    }

    /// Access the intensity array.
    ///
    /// # Panics
    /// This function will panic if [`ArrayType::IntensityArray`] is not
    /// present in [`RawSpectrum::arrays`]
    pub fn intensities(&'lifespan self) -> Cow<'transient, [f32]> {
        self.arrays.intensities().unwrap()
    }

    pub fn mzs_mut(&mut self) -> Result<&mut [f64], ArrayRetrievalError> {
        self.arrays.mzs_mut()
    }

    pub fn intensities_mut(&mut self) -> Result<&mut [f32], ArrayRetrievalError> {
        self.arrays.intensities_mut()
    }

    /// Explicitly decode any [`DataArray`] that is encoded or compressed still
    /// so that they are ready for use.
    ///
    /// # See also
    /// [`BinaryArrayMap::decode_all_arrays`]
    pub fn decode_all_arrays(&mut self) -> Result<(), ArrayRetrievalError> {
        self.arrays.decode_all_arrays()
    }

    /// Convert a spectrum into a [`MultiLayerSpectrum`].
    pub fn into_spectrum<C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default>(
        self,
    ) -> Result<MultiLayerSpectrum<C, D>, SpectrumConversionError>
    where
        C: BuildFromArrayMap,
    {
        Ok(MultiLayerSpectrum::<C, D> {
            arrays: Some(self.arrays),
            description: self.description,
            ..Default::default()
        })
    }

    /// Apply a local denoising algorithm to the signal that is suitable for profile mode
    /// data.
    #[cfg(feature = "mzsignal")]
    pub fn denoise(&mut self, scale: f32) -> Result<(), SpectrumProcessingError> {
        let mut intensities_copy = self.arrays.intensities()?.into_owned();
        let mz_array = self.arrays.mzs()?;
        match denoise(&mz_array, &mut intensities_copy, scale) {
            Ok(_) => {
                let view = self.arrays.get_mut(&ArrayType::IntensityArray).unwrap();
                view.store_as(BinaryDataArrayType::Float32)
                    .expect("Failed to reformat intensity array");
                view.update_buffer(&intensities_copy)
                    .expect("Failed to update intensity array buffer");
                Ok(())
            }
            Err(err) => Err(SpectrumProcessingError::DenoisingError(err)),
        }
    }

    /// pick peaks with `peak_picker` and convert this spectrum into a [`MultiLayerSpectrum`] with a centroid peak list
    /// as well as raw data arrays.
    ///
    /// # See also
    /// [`MultiLayerSpectrum::pick_peaks_with`]
    #[cfg(feature = "mzsignal")]
    pub fn pick_peaks_with_into<
        C: CentroidPeakAdapting + BuildFromArrayMap + From<FittedPeak>,
        D: DeconvolutedPeakAdapting + BuildFromArrayMap,
    >(
        self,
        peak_picker: &PeakPicker,
    ) -> Result<MultiLayerSpectrum<C, D>, SpectrumProcessingError> {
        let mut result = self.into_spectrum()?;
        match result.pick_peaks_with(peak_picker) {
            Ok(_) => Ok(result),
            Err(err) => Err(err),
        }
    }

    /// Pick peaks with a minimum signal-to-noise threshold and convert this spectrum into a [`MultiLayerSpectrum`] with
    /// a centroid peak list as well as raw data arrays.
    ///
    /// # See also
    /// [`RawSpectrum::pick_peaks_with_into`]
    #[cfg(feature = "mzsignal")]
    pub fn pick_peaks_into<
        C: CentroidPeakAdapting + BuildFromArrayMap + From<FittedPeak>,
        D: DeconvolutedPeakAdapting + BuildFromArrayMap,
    >(
        self,
        signal_to_noise_threshold: f32,
    ) -> Result<MultiLayerSpectrum<C, D>, SpectrumProcessingError> {
        let peak_picker = PeakPicker {
            fit_type: PeakFitType::Quadratic,
            signal_to_noise_threshold,
            ..Default::default()
        };
        self.pick_peaks_with_into(&peak_picker)
    }

    /// Pick peaks in the specified m/z ranges with a minimum signal-to-noise threshold and convert
    /// this spectrum into a [`MultiLayerSpectrum`] with a centroid peak list as well as raw data arrays.
    ///
    /// # See also
    /// [`MultiLayerSpectrum::pick_peaks_in_intervals_into`]
    #[cfg(feature = "mzsignal")]
    pub fn pick_peaks_in_intervals_into<
        C: CentroidPeakAdapting + BuildFromArrayMap + From<FittedPeak>,
        D: DeconvolutedPeakAdapting + BuildFromArrayMap,
    >(
        self,
        signal_to_noise_threshold: f32,
        intervals: &[(f64, f64)],
    ) -> Result<MultiLayerSpectrum<C, D>, SpectrumProcessingError> {
        let mut result = self.into_spectrum()?;
        match result.pick_peaks_in_intervals(signal_to_noise_threshold, intervals) {
            Ok(_) => Ok(result),
            Err(err) => Err(err),
        }
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> SpectrumLike<C, D> for RawSpectrum {
    #[inline]
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }

    fn description_mut(&mut self) -> &mut SpectrumDescription {
        &mut self.description
    }

    fn peaks(&'_ self) -> RefPeakDataLevel<'_, C, D> {
        RefPeakDataLevel::RawData(&self.arrays)
    }

    fn raw_arrays(&'_ self) -> Option<&'_ BinaryArrayMap> {
        Some(&self.arrays)
    }

    fn into_peaks_and_description(
        self,
    ) -> (
        PeakDataLevel<C, D>,
        SpectrumDescription,
    ) {
        (PeakDataLevel::RawData(self.arrays), self.description)
    }
}

#[derive(Default, Debug, Clone)]
/// Represents a spectrum that has been centroided into discrete m/z points, a
/// process also called "peak picking".
///
/// This type of spectrum represents data in exactly one format.
pub struct CentroidSpectrumType<C: CentroidLike + Default> {
    /// The spectrum metadata describing acquisition conditions and details.
    pub description: SpectrumDescription,
    /// The picked centroid peaks, sorted by m/z in a fast searchable structure.
    pub peaks: MZPeakSetType<C>,
}

#[cfg(feature = "mzsignal")]
impl<C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap> CentroidSpectrumType<C> {
    /// Convert the centroid peaks in theoretical profile signal, producing a [`MultiLayerSpectrum`]
    /// which contains both the `peaks` as well as a raw data array in profile mode.
    ///
    /// # Arguments
    /// - `dx`: The m/z spacing to reprofile over. The smaller the value, the higher the resolution, but also
    ///    the more memory required. For high resolution instruments, a value between `0.005` and `0.001` is suitable.
    /// - `fwhm`: The full width at half max to assume when reprofiling each peak. This affects the peak shape width.
    pub fn reprofile_with_shape_into(
        self,
        dx: f64,
        fwhm: f32,
    ) -> Result<MultiLayerSpectrum<C, DeconvolutedPeak>, SpectrumProcessingError> {
        let mut spectrum = self.into_spectrum()?;
        spectrum.reprofile_with_shape(dx, fwhm)?;
        Ok(spectrum)
    }
}

impl<C: CentroidLike + Default> ParamDescribed for CentroidSpectrumType<C> {
    fn params(&self) -> &[crate::params::Param] {
        <SpectrumDescription as ParamDescribed>::params(&self.description)
    }

    fn params_mut(&mut self) -> &mut crate::params::ParamList {
        <SpectrumDescription as ParamDescribed>::params_mut(&mut self.description)
    }
}

impl<C: CentroidLike + Default> SpectrumLike<C> for CentroidSpectrumType<C> {
    #[inline]
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }

    fn description_mut(&mut self) -> &mut SpectrumDescription {
        &mut self.description
    }

    fn peaks(&'_ self) -> RefPeakDataLevel<'_, C, DeconvolutedPeak> {
        RefPeakDataLevel::Centroid(&self.peaks)
    }

    fn raw_arrays(&'_ self) -> Option<&'_ BinaryArrayMap> {
        None
    }

    fn into_peaks_and_description(
        self,
    ) -> (PeakDataLevel<C, DeconvolutedPeak>, SpectrumDescription) {
        (PeakDataLevel::Centroid(self.peaks), self.description)
    }
}

impl<C: CentroidLike + Default> CentroidSpectrumType<C> {
    pub fn new(description: SpectrumDescription, peaks: MZPeakSetType<C>) -> Self {
        Self { description, peaks }
    }

    /// Convert a spectrum into a [`MultiLayerSpectrum`]
    pub fn into_spectrum<D>(self) -> Result<MultiLayerSpectrum<C, D>, SpectrumConversionError>
    where
        D: DeconvolutedCentroidLike + Default,
    {
        let val = MultiLayerSpectrum::<C, D> {
            peaks: Some(self.peaks),
            description: self.description,
            ..Default::default()
        };
        Ok(val)
    }
}

pub type CentroidSpectrum = CentroidSpectrumType<CentroidPeak>;

impl<C: CentroidPeakAdapting> Index<usize> for CentroidSpectrumType<C> {
    type Output = <MZPeakSetType<C> as Index<usize>>::Output;

    fn index(&self, index: usize) -> &Self::Output {
        self.peaks.index(index)
    }
}

impl<C: CentroidPeakAdapting> PeakCollection<C, MZ> for CentroidSpectrumType<C>
where
    MZPeakSetType<C>: PeakCollection<C, MZ>,
    <MZPeakSetType<C> as Index<usize>>::Output: CoordinateLike<MZ>,
{
    fn len(&self) -> usize {
        self.peaks.len()
    }

    fn get_item(&self, i: usize) -> &C {
        self.peaks.get_item(i)
    }

    fn get_slice(&self, i: std::ops::Range<usize>) -> &[C] {
        self.peaks.get_slice(i)
    }

    fn search_by(&self, query: f64) -> Result<usize, usize> {
        self.peaks.search_by(query)
    }

    fn iter(&self) -> impl Iterator<Item = &C>
    where
        C: 'static,
    {
        self.peaks.iter()
    }
}

/// Represents a spectrum that has been centroided, deisotoped, and charge state deconvolved.
///
/// This type of spectrum represents data in exactly one format.
#[derive(Default, Debug, Clone)]
pub struct DeconvolutedSpectrumType<D: DeconvolutedCentroidLike + Default> {
    /// The spectrum metadata describing acquisition conditions and details.
    pub description: SpectrumDescription,
    /// The deisotoped and charge state deconvolved peaks, sorted by neutral mass
    /// in a fast searchable structure.
    pub deconvoluted_peaks: MassPeakSetType<D>,
}

impl<D: DeconvolutedCentroidLike + Default> DeconvolutedSpectrumType<D> {
    pub fn new(description: SpectrumDescription, deconvoluted_peaks: MassPeakSetType<D>) -> Self {
        Self {
            description,
            deconvoluted_peaks,
        }
    }

    /// Convert a spectrum into a [`MultiLayerSpectrum`]
    pub fn into_spectrum<C>(self) -> Result<MultiLayerSpectrum<C, D>, SpectrumConversionError>
    where
        C: CentroidLike + Default,
    {
        let val = MultiLayerSpectrum::<C, D> {
            deconvoluted_peaks: Some(self.deconvoluted_peaks),
            description: self.description,
            ..Default::default()
        };
        Ok(val)
    }
}

impl<D: DeconvolutedCentroidLike + Default> ParamDescribed for DeconvolutedSpectrumType<D> {
    fn params(&self) -> &[crate::params::Param] {
        <SpectrumDescription as ParamDescribed>::params(&self.description)
    }

    fn params_mut(&mut self) -> &mut crate::params::ParamList {
        <SpectrumDescription as ParamDescribed>::params_mut(&mut self.description)
    }
}

impl<D: DeconvolutedCentroidLike + Default> SpectrumLike<CentroidPeak, D>
    for DeconvolutedSpectrumType<D>
{
    #[inline]
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }

    fn description_mut(&mut self) -> &mut SpectrumDescription {
        &mut self.description
    }

    fn peaks(&'_ self) -> RefPeakDataLevel<'_, CentroidPeak, D> {
        RefPeakDataLevel::Deconvoluted(&self.deconvoluted_peaks)
    }

    fn raw_arrays(&'_ self) -> Option<&'_ BinaryArrayMap> {
        None
    }

    fn into_peaks_and_description(self) -> (PeakDataLevel<CentroidPeak, D>, SpectrumDescription) {
        (
            PeakDataLevel::Deconvoluted(self.deconvoluted_peaks),
            self.description,
        )
    }
}

pub type DeconvolutedSpectrum = DeconvolutedSpectrumType<DeconvolutedPeak>;

impl<D: DeconvolutedPeakAdapting> Index<usize> for DeconvolutedSpectrumType<D> {
    type Output = <MassPeakSetType<D> as Index<usize>>::Output;

    fn index(&self, index: usize) -> &Self::Output {
        self.deconvoluted_peaks.index(index)
    }
}

impl<D: DeconvolutedPeakAdapting> PeakCollection<D, Mass> for DeconvolutedSpectrumType<D>
where
    MassPeakSetType<D>: PeakCollection<D, Mass>,
    <MassPeakSetType<D> as Index<usize>>::Output: CoordinateLike<Mass>,
{
    fn len(&self) -> usize {
        self.deconvoluted_peaks.len()
    }

    fn get_item(&self, i: usize) -> &D {
        self.deconvoluted_peaks.get_item(i)
    }

    fn get_slice(&self, i: std::ops::Range<usize>) -> &[D] {
        self.deconvoluted_peaks.get_slice(i)
    }

    fn search_by(&self, query: f64) -> Result<usize, usize> {
        self.deconvoluted_peaks.search_by(query)
    }

    fn iter(&self) -> impl Iterator<Item = &D>
    where
        D: 'static,
    {
        self.deconvoluted_peaks.iter()
    }
}


#[derive(Default, Debug, Clone)]
/// Represent a spectrum with multiple layers of representation of the
/// peak data.
///
/// This type is useful because it permits us to hold spectra in incrementally
pub struct MultiLayerSpectrum<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
> {
    /// The spectrum metadata describing acquisition conditions and details.
    pub description: SpectrumDescription,

    /// The (potentially absent) data arrays describing the m/z, intensity,
    /// and potentially other measured properties
    pub arrays: Option<BinaryArrayMap>,

    // The (potentially absent) centroid peaks
    pub peaks: Option<MZPeakSetType<C>>,

    // The (potentially absent) deconvoluted peaks
    pub deconvoluted_peaks: Option<MassPeakSetType<D>>,
}

impl<C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default> ParamDescribed
    for MultiLayerSpectrum<C, D>
{
    fn params(&self) -> &[crate::params::Param] {
        <SpectrumDescription as ParamDescribed>::params(&self.description)
    }

    fn params_mut(&mut self) -> &mut crate::params::ParamList {
        <SpectrumDescription as ParamDescribed>::params_mut(&mut self.description)
    }
}

impl<C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default> SpectrumLike<C, D>
    for MultiLayerSpectrum<C, D>
{
    #[inline]
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }

    fn description_mut(&mut self) -> &mut SpectrumDescription {
        &mut self.description
    }

    fn peaks(&'_ self) -> RefPeakDataLevel<'_, C, D> {
        if let Some(peaks) = &self.deconvoluted_peaks {
            RefPeakDataLevel::Deconvoluted(peaks)
        } else if let Some(peaks) = &self.peaks {
            RefPeakDataLevel::Centroid(peaks)
        } else if let Some(arrays) = &self.arrays {
            RefPeakDataLevel::RawData(arrays)
        } else {
            RefPeakDataLevel::Missing
        }
    }

    fn raw_arrays(&'_ self) -> Option<&'_ BinaryArrayMap> {
        self.arrays.as_ref()
    }

    fn into_peaks_and_description(self) -> (PeakDataLevel<C, D>, SpectrumDescription) {
        if let Some(peaks) = self.deconvoluted_peaks {
            (PeakDataLevel::Deconvoluted(peaks), self.description)
        } else if let Some(peaks) = self.peaks {
            (PeakDataLevel::Centroid(peaks), self.description)
        } else if let Some(arrays) = self.arrays {
            (PeakDataLevel::RawData(arrays), self.description)
        } else {
            (PeakDataLevel::Missing, self.description)
        }
    }
}

impl<C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default> MultiLayerSpectrum<C, D>
where
    C: BuildFromArrayMap + BuildArrayMapFrom,
    D: BuildFromArrayMap + BuildArrayMapFrom,
{
    pub fn new(
        description: SpectrumDescription,
        arrays: Option<BinaryArrayMap>,
        peaks: Option<MZPeakSetType<C>>,
        deconvoluted_peaks: Option<MassPeakSetType<D>>,
    ) -> Self {
        Self {
            description,
            arrays,
            peaks,
            deconvoluted_peaks,
        }
    }

    /// Convert a spectrum into a [`CentroidSpectrumType`]
    pub fn into_centroid(self) -> Result<CentroidSpectrumType<C>, SpectrumConversionError> {
        if let Some(peaks) = self.peaks {
            let mut result = CentroidSpectrumType::<C> {
                peaks,
                description: self.description,
            };
            result.description.signal_continuity = SignalContinuity::Centroid;
            return Ok(result);
        } else if self.signal_continuity() == SignalContinuity::Centroid {
            if let Some(arrays) = &self.arrays {
                let peaks = C::try_from_arrays(arrays)?.into();
                let mut centroid = CentroidSpectrumType::<C> {
                    description: self.description,
                    peaks,
                };
                centroid.description.signal_continuity = SignalContinuity::Centroid;
                return Ok(centroid);
            } else {
                let mut result = CentroidSpectrumType::<C> {
                    peaks: MZPeakSetType::<C>::empty(),
                    description: self.description,
                };
                result.description.signal_continuity = SignalContinuity::Centroid;
                return Ok(result);
            }
        }
        Err(SpectrumConversionError::NotCentroided)
    }

    /// Convert a spectrum into a [`RawSpectrum`]
    pub fn into_raw(self) -> Result<RawSpectrum, SpectrumConversionError> {
        if let Some(arrays) = self.arrays {
            Ok(RawSpectrum {
                arrays,
                description: self.description,
            })
        } else if let Some(peaks) = self.deconvoluted_peaks {
            let arrays = D::as_arrays(&peaks[0..]);
            let mut result = RawSpectrum {
                arrays,
                description: self.description,
            };
            result.description.signal_continuity = SignalContinuity::Centroid;
            Ok(result)
        } else if let Some(peaks) = self.peaks {
            let arrays = C::as_arrays(&peaks[0..]);

            let mut result = RawSpectrum {
                arrays,
                description: self.description,
            };
            result.description.signal_continuity = SignalContinuity::Centroid;
            Ok(result)
        } else {
            let arrays = BinaryArrayMap::from(&PeakSet::empty());
            Ok(RawSpectrum {
                arrays,
                description: self.description,
            })
        }
    }

    pub fn from_arrays_and_description(
        arrays: BinaryArrayMap,
        description: SpectrumDescription,
    ) -> Self {
        Self {
            description,
            arrays: Some(arrays),
            ..Default::default()
        }
    }

    pub fn from_description(description: SpectrumDescription) -> Self {
        Self {
            description,
            ..Default::default()
        }
    }

    pub fn from_peaks_data_levels_and_description(
        peaks: PeakDataLevel<C, D>,
        description: SpectrumDescription,
    ) -> Self {
        match peaks {
            PeakDataLevel::Missing => Self::new(description, None, None, None),
            PeakDataLevel::RawData(arrays) => Self::new(description, Some(arrays), None, None),
            PeakDataLevel::Centroid(peaks) => Self::new(description, None, Some(peaks), None),
            PeakDataLevel::Deconvoluted(peaks) => Self::new(description, None, None, Some(peaks)),
        }
    }

    pub fn from_spectrum_like<S: SpectrumLike<C, D>>(spectrum: S) -> Self {
        let (peaks, description) = spectrum.into_peaks_and_description();
        Self::from_peaks_data_levels_and_description(peaks, description)
    }

    #[cfg(feature = "mzsignal")]
    pub fn denoise(&mut self, scale: f32) -> Result<(), SpectrumProcessingError> {
        match &mut self.arrays {
            Some(arrays) => {
                let mut intensities_copy = arrays.intensities()?.into_owned();
                let mz_array = arrays.mzs()?;
                match denoise(&mz_array, &mut intensities_copy, scale) {
                    Ok(_) => {
                        let view = arrays.get_mut(&ArrayType::IntensityArray).unwrap();
                        view.store_as(BinaryDataArrayType::Float32)?;
                        view.update_buffer(&intensities_copy)?;
                        Ok(())
                    }
                    Err(err) => Err(SpectrumProcessingError::DenoisingError(err)),
                }
            }
            None => Err(SpectrumProcessingError::SpectrumConversionError(
                SpectrumConversionError::NoPeakData,
            )),
        }
    }

    #[cfg(feature = "mzsignal")]
    pub fn reprofile_with_shape(
        &mut self,
        dx: f64,
        fwhm: f32,
    ) -> Result<(), SpectrumProcessingError> {
        if let Some(peaks) = self.peaks.as_ref() {
            let reprofiler = reprofile::PeakSetReprofiler::new(
                peaks.first().map(|p| p.mz() - 1.0).unwrap_or_default(),
                peaks.last().map(|p| p.mz() + 1.0).unwrap_or_default(),
                dx,
            );

            let models: Vec<_> = peaks
                .iter()
                .map(|p| {
                    PeakShapeModel::from_centroid(p.mz(), p.intensity(), fwhm, PeakShape::Gaussian)
                })
                .collect();

            let pair = reprofiler.reprofile_from_models(&models);

            let arrays = match self.arrays.as_mut() {
                Some(arrays) => arrays,
                None => {
                    self.arrays = Some(BinaryArrayMap::new());
                    self.arrays.as_mut().unwrap()
                }
            };
            arrays.add(DataArray::wrap(
                &ArrayType::MZArray,
                BinaryDataArrayType::Float64,
                pair.mz_array
                    .iter()
                    .map(|i| i.to_le_bytes())
                    .flatten()
                    .collect(),
            ));

            arrays.add(DataArray::wrap(
                &ArrayType::IntensityArray,
                BinaryDataArrayType::Float32,
                pair.intensity_array
                    .iter()
                    .map(|i| i.to_le_bytes())
                    .flatten()
                    .collect(),
            ));
            Ok(())
        } else {
            Err(SpectrumConversionError::NoPeakData.into())
        }
    }
}

impl<C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default> MultiLayerSpectrum<C, D>
where
    C: BuildFromArrayMap,
    D: BuildFromArrayMap,
{
    /// Attempt to reconstruct one of the peak layers based upon the available data arrays
    /// if this spectrum is in centroid mode.
    ///
    /// # Errors
    /// If constructing the peak list encounters an error, or if data array access/decoding fails,
    /// a [`SpectrumConversionError`] is returned.
    pub fn try_build_peaks(&mut self) -> Result<RefPeakDataLevel<'_, C, D>, SpectrumConversionError> {
        if matches!(self.signal_continuity(), SignalContinuity::Centroid) {
            if let Some(arrays) = self.arrays.as_ref() {
                if let ArraysAvailable::Ok = D::has_arrays_for(arrays) {
                    let peaks = D::try_from_arrays(arrays)?.into();
                    self.deconvoluted_peaks = Some(peaks);
                    return Ok(self.peaks())
                }

                if let ArraysAvailable::Ok = C::has_arrays_for(arrays) {
                    let peaks = C::try_from_arrays(arrays)?.into();
                    self.peaks = Some(peaks);
                    return Ok(self.peaks())
                }
                return Ok(RefPeakDataLevel::Missing)
            } else {
                return Ok(self.peaks())
            }
        } else {
            return Ok(self.peaks())
        }
    }

    pub fn try_build_centroids(&mut self) -> Result<&MZPeakSetType<C>, SpectrumConversionError> {
        if self.peaks.is_some() {
            Ok(self.peaks.as_ref().unwrap())
        } else if let Some(arrays) = self.arrays.as_ref() {
            match self.signal_continuity() {
                SignalContinuity::Centroid => {
                    let peaks = C::try_from_arrays(arrays)?.into();
                    self.peaks = Some(peaks);
                    Ok(self.peaks.as_ref().unwrap())
                }
                _ => Err(SpectrumConversionError::NotCentroided),
            }
        } else {
            Err(SpectrumConversionError::NoPeakData)
        }
    }

    pub fn try_build_deconvoluted_centroids(
        &mut self,
    ) -> Result<&MassPeakSetType<D>, SpectrumConversionError> {
        if let Some(ref peaks) = self.deconvoluted_peaks {
            Ok(peaks)
        } else if let Some(ref arrays) = self.arrays {
            match self.signal_continuity() {
                SignalContinuity::Centroid => {
                    let peaks = D::try_from_arrays(arrays)?.into();
                    self.deconvoluted_peaks = Some(peaks);
                    Ok(self.deconvoluted_peaks.as_ref().unwrap())
                }
                _ => Err(SpectrumConversionError::NotCentroided),
            }
        } else {
            Err(SpectrumConversionError::NoPeakData)
        }
    }
}

#[cfg(feature = "mzsignal")]
impl<C: CentroidLike + Default + From<FittedPeak>, D: DeconvolutedCentroidLike + Default>
    MultiLayerSpectrum<C, D>
{
    /// Using a pre-configured [`mzsignal::PeakPicker`](mzsignal::peak_picker::PeakPicker) to pick peaks
    /// from `arrays`'s signal, populating the `peaks` field.
    ///
    /// If [`SpectrumLike::signal_continuity`] returns [`SignalContinuity::Centroid`], then no filtering is
    /// performed and all points are directly converted into picked peaks.
    pub fn pick_peaks_with(
        &mut self,
        peak_picker: &PeakPicker,
    ) -> Result<(), SpectrumProcessingError> {
        if let Some(arrays) = &self.arrays {
            let mz_array = arrays.mzs()?;
            let intensity_array = arrays.intensities()?;

            if matches!(self.signal_continuity(), SignalContinuity::Centroid) {
                let mut peaks: MZPeakSetType<C> = mz_array
                    .iter()
                    .zip(intensity_array.iter())
                    .map(|(mz, inten)| FittedPeak::new(*mz, *inten, 0, 0.0, 0.0).into())
                    .collect();
                peaks.sort();
                self.peaks = Some(peaks);
                Ok(())
            } else {
                let mut acc = Vec::new();
                match peak_picker.discover_peaks(&mz_array, &intensity_array, &mut acc) {
                    Ok(_) => {
                        let peaks: MZPeakSetType<C> = acc.into_iter().map(|p| C::from(p)).collect();
                        self.peaks = Some(peaks);
                        Ok(())
                    }
                    Err(err) => Err(SpectrumProcessingError::PeakPickerError(err)),
                }
            }
        } else {
            Err(SpectrumProcessingError::SpectrumConversionError(
                SpectrumConversionError::NoPeakData,
            ))
        }
    }

    /// Pick peaks with a minimum signal-to-noise threshold, populating [`MultiLayerSpectrum::peaks`] with
    /// an [`mzpeaks::MZPeakSetType`] over `C`.
    ///
    /// # See also
    /// [`MultiLayerSpectrum::pick_peaks_with`]
    pub fn pick_peaks(
        &mut self,
        signal_to_noise_threshold: f32,
    ) -> Result<(), SpectrumProcessingError> {
        let peak_picker = PeakPicker {
            fit_type: PeakFitType::Quadratic,
            signal_to_noise_threshold,
            ..Default::default()
        };
        self.pick_peaks_with(&peak_picker)
    }

    /// Pick peaks in the specified m/z intervals with a minimum signal-to-noise threshold,
    /// populating [`MultiLayerSpectrum::peaks`] with an [`mzpeaks::MZPeakSetType`] over `C`
    /// from those intervals.
    ///
    /// If [`SpectrumLike::signal_continuity`] returns [`SignalContinuity::Centroid`], then
    /// array data is converted to peaks as in [`MultiLayerSpectrum::pick_peaks`] and then
    /// only those peaks in the requested intervals are kept.
    pub fn pick_peaks_in_intervals(
        &mut self,
        signal_to_noise_threshold: f32,
        intervals: &[(f64, f64)],
    ) -> Result<(), SpectrumProcessingError> {
        let peak_picker = PeakPicker {
            fit_type: PeakFitType::Quadratic,
            signal_to_noise_threshold,
            ..Default::default()
        };

        if matches!(self.signal_continuity(), SignalContinuity::Centroid) {
            self.pick_peaks_with(&peak_picker)?;
            let peaks = self
                .peaks
                .take()
                .unwrap()
                .into_iter()
                .filter(|p| {
                    let mz = p.mz();
                    intervals.iter().any(|(lo, hi)| *lo <= mz && mz < *hi)
                })
                .collect();
            self.peaks = Some(peaks);
            Ok(())
        } else if let Some(arrays) = &self.arrays {
            let mz_array = arrays.mzs()?;
            let intensity_array = arrays.intensities()?;
            let mut acc = Vec::new();
            for (mz_start, mz_end) in intervals {
                match peak_picker.discover_peaks_in_interval(
                    &mz_array,
                    &intensity_array,
                    &mut acc,
                    *mz_start,
                    *mz_end,
                ) {
                    Ok(_) => {}
                    Err(err) => return Err(SpectrumProcessingError::PeakPickerError(err)),
                }
            }
            let peaks: MZPeakSetType<C> = acc.into_iter().map(|p| C::from(p)).collect();
            self.peaks = Some(peaks);
            return Ok(());
        } else {
            return Err(SpectrumProcessingError::SpectrumConversionError(
                SpectrumConversionError::NoPeakData,
            ));
        }
    }
}

impl<C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default>
    TryFrom<MultiLayerSpectrum<C, D>> for CentroidSpectrumType<C>
where
    C: BuildFromArrayMap + BuildArrayMapFrom,
    D: BuildFromArrayMap + BuildArrayMapFrom,
{
    type Error = SpectrumConversionError;

    fn try_from(
        spectrum: MultiLayerSpectrum<C, D>,
    ) -> Result<CentroidSpectrumType<C>, Self::Error> {
        spectrum.into_centroid()
    }
}

/// A ready-to-use parameterized type for representing a spectrum
/// in multiple overlapping layers.
pub type Spectrum = MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak>;

impl<C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default> From<CentroidSpectrumType<C>>
    for MultiLayerSpectrum<C, D>
{
    fn from(spectrum: CentroidSpectrumType<C>) -> MultiLayerSpectrum<C, D> {
        spectrum.into_spectrum().unwrap()
    }
}

impl<C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default> From<RawSpectrum>
    for MultiLayerSpectrum<C, D>
where
    C: BuildFromArrayMap,
    D: BuildFromArrayMap,
{
    fn from(spectrum: RawSpectrum) -> MultiLayerSpectrum<C, D> {
        spectrum.into_spectrum().unwrap()
    }
}

impl<C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default>
    From<MultiLayerSpectrum<C, D>> for RawSpectrum
where
    C: BuildFromArrayMap + BuildArrayMapFrom,
    D: BuildFromArrayMap + BuildArrayMapFrom,
{
    fn from(spectrum: MultiLayerSpectrum<C, D>) -> RawSpectrum {
        spectrum.into_raw().unwrap()
    }
}

impl<C: CentroidLike + Default> TryFrom<RawSpectrum> for CentroidSpectrumType<C>
where
    C: BuildFromArrayMap,
{
    type Error = SpectrumConversionError;

    fn try_from(spectrum: RawSpectrum) -> Result<CentroidSpectrumType<C>, Self::Error> {
        spectrum.into_centroid()
    }
}

impl<C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default>
    TryFrom<MultiLayerSpectrum<C, D>> for DeconvolutedSpectrumType<D>
where
    C: BuildFromArrayMap + BuildArrayMapFrom,
    D: BuildFromArrayMap + BuildArrayMapFrom,
{
    type Error = SpectrumConversionError;

    fn try_from(value: MultiLayerSpectrum<C, D>) -> Result<Self, Self::Error> {
        match value.peaks() {
            RefPeakDataLevel::Deconvoluted(_) => {
                let peaks = value.deconvoluted_peaks.unwrap();
                Ok(DeconvolutedSpectrumType {
                    description: value.description,
                    deconvoluted_peaks: peaks,
                })
            }
            RefPeakDataLevel::RawData(arrays) => {
                let peaks = match D::try_from_arrays(arrays) {
                    Ok(peaks) => peaks,
                    Err(e) => return Err(SpectrumConversionError::ArrayRetrievalError(e)),
                };
                Ok(DeconvolutedSpectrumType {
                    description: value.description,
                    deconvoluted_peaks: PeakSetVec::new(peaks),
                })
            }
            _ => Err(SpectrumConversionError::NotDeconvoluted),
        }
    }
}

#[cfg(test)]
mod test {

    #[cfg(feature = "mzsignal")]
    #[test_log::test]
    fn test_profile_read() {
        use super::*;
        use crate::io::mzml::MzMLReader;
        use crate::prelude::*;
        let mut reader = MzMLReader::open_path("./test/data/three_test_scans.mzML")
            .expect("Failed to open test file");
        reader.reset();
        let mut scan = reader.next().expect("Failed to read spectrum");
        assert_eq!(scan.signal_continuity(), SignalContinuity::Profile);
        assert_eq!(scan.ms_level(), 1);
        assert_eq!(scan.polarity(), ScanPolarity::Positive);
        assert!(scan.precursor().is_none());

        if let Err(err) = scan.pick_peaks(1.0) {
            panic!("Should not have an error! {}", err);
        }

        if let Some(peaks) = &scan.peaks {
            assert_eq!(peaks.len(), 2107);
            let n = peaks.len();
            for i in 1..n {
                let diff = peaks[i].get_index() - peaks[i - 1].get_index();
                assert_eq!(diff, 1);
            }

            peaks
                .has_peak(562.741, Tolerance::PPM(3f64))
                .expect("Expected to find peak");
            peaks
                .has_peak(563.240, Tolerance::PPM(3f64))
                .expect("Expected to find peak");
            let p = peaks
                .has_peak(563.739, Tolerance::PPM(1f64))
                .expect("Expected to find peak");
            assert!((p.mz() - 563.739).abs() < 1e-3)
        }
    }
}
