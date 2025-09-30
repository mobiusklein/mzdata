use std::convert::TryFrom;
use std::{borrow::Cow, ops::Index};

use mzpeaks::Mass;
use thiserror::Error;

use mzpeaks::{
    CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak, MZ, MZPeakSetType,
    MassPeakSetType, PeakCollection, PeakSet, peak_set::PeakSetVec, prelude::*,
};

#[cfg(feature = "mzsignal")]
use mzsignal::{
    FittedPeak,
    denoise::{DenoisingError, denoise},
    peak_picker::{PeakFitType, PeakPicker, PeakPickerError},
    reprofile::{self, PeakShape, PeakShapeModel},
};

use crate::params::{ParamDescribed, Unit, Value};
#[allow(unused)]
use crate::spectrum::bindata::{ArrayType, BinaryArrayMap, BinaryDataArrayType};
use crate::spectrum::peaks::{PeakDataLevel, RefPeakDataLevel, SpectrumSummary};
use crate::spectrum::scan_properties::{
    Acquisition, IonMobilityMeasure, Precursor, ScanPolarity, SignalContinuity, SpectrumDescription,
};

#[allow(unused)]
use super::DataArray;
use super::HasIonMobility;
use super::bindata::{ArrayRetrievalError, BuildArrayMapFrom, BuildFromArrayMap};

/// A blanket trait that ties together all the assumed behaviors of an m/z coordinate centroid peak
pub trait CentroidPeakAdapting: CentroidLike + From<CentroidPeak> {}
impl<C: CentroidLike + From<CentroidPeak>> CentroidPeakAdapting for C {}

/// A blanket trait that ties together all the assumed behaviors of an neutral mass coordinate centroid peak
pub trait DeconvolutedPeakAdapting: DeconvolutedCentroidLike + From<DeconvolutedPeak> {}
impl<D: DeconvolutedCentroidLike + From<DeconvolutedPeak>> DeconvolutedPeakAdapting for D {}

/// A trait for providing a uniform delegated access to spectrum metadata
pub trait SpectrumLike<
    C: CentroidLike = CentroidPeak,
    D: DeconvolutedCentroidLike = DeconvolutedPeak,
>: ParamDescribed
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

    /// Access the (first) precursor information, if it exists.
    #[inline]
    fn precursor(&self) -> Option<&Precursor> {
        let desc = self.description();
        desc.precursor.first()
    }

    /// Iterate over all precursors of the spectrum
    fn precursor_iter(&self) -> impl Iterator<Item = &Precursor> {
        let desc = self.description();
        desc.precursor.iter()
    }

    /// Mutably access the (first) precursor information, if it exists
    fn precursor_mut(&mut self) -> Option<&mut Precursor> {
        let desc = self.description_mut();
        desc.precursor.first_mut()
    }

    /// Iterate over all precursors of the spectrum mutably
    fn precursor_iter_mut(&mut self) -> impl Iterator<Item = &mut Precursor> {
        let desc = self.description_mut();
        desc.precursor.iter_mut()
    }

    /// Add a precursor to the list of precursors for this spectrum.
    ///
    /// Precursors beyond the first one correspond to lower exponentiated spectra, e.g. for an MS3 spectrum
    /// the first precursor is the MS2 product ion that was selected, and the second precursor corresponds
    /// to the original MS1 ion that was chosen for MS2.
    ///
    /// Higher order precursors may be accessed with [`SpectrumLike::precursor_iter`].
    fn add_precursor(&mut self, precursor: Precursor) {
        self.description_mut().precursor.push(precursor);
    }

    /// Remove the precursor entry at `index` in the precursor list.
    ///
    /// Care should be taken if you are attempting to re-order precursors.
    /// It may be simpler to use [`SpectrumLike::precursor_iter_mut`] to
    /// re-arrange entries.
    fn remove_precursor(&mut self, index: usize) -> Precursor {
        self.description_mut().precursor.remove(index)
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

    /// Access a description of the spectrum polarity
    #[inline]
    fn polarity(&self) -> ScanPolarity {
        self.description().polarity
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

    /// Retrieve the most processed representation of the mass spectrum's signal
    fn peaks(&'_ self) -> RefPeakDataLevel<'_, C, D>;

    /// Consume the spectrum, decomposing it into the [`SpectrumDescription`] and an owning
    /// [`PeakDataLevel`].
    ///
    /// # See also
    /// [`SpectrumLike::peaks`]
    /// [`SpectrumLike::description`]
    fn into_peaks_and_description(self) -> (PeakDataLevel<C, D>, SpectrumDescription);

    /// Obtain a reference to the [`BinaryArrayMap`] if one is available for the peak
    /// information. This may not be the most refined version of the peak signal if
    /// it has been processed further.
    fn raw_arrays(&'_ self) -> Option<&'_ BinaryArrayMap>;

    /// Check if this spectrum has an ion mobility dimension/array. This is distinct from
    /// having a scan-level point measure of ion mobility.
    fn has_ion_mobility_dimension(&self) -> bool {
        self.raw_arrays()
            .map(|a| a.has_ion_mobility())
            .unwrap_or_default()
    }

    /// A single point of entry for checking if ion mobility is present.
    fn has_ion_mobility_class(&self) -> HasIonMobility {
        if self.has_ion_mobility_dimension() {
            HasIonMobility::Dimension
        } else if self.has_ion_mobility() {
            HasIonMobility::Point
        } else {
            HasIonMobility::None
        }
    }

    /// Find the type of spectrum described.
    ///
    /// A spectrum in `mzdata` is *usually* a mass spectrum of some sort,
    /// but that's not guaranteed to be the case. `mzdata` can handle non-MS spectra,
    /// but little of the signal processing machinery it provides currently supports
    /// those other kinds of data.
    fn spectrum_type(&self) -> Option<crate::meta::SpectrumType> {
        self.description().spectrum_type()
    }

    /// Set the kind of spectrum represented.
    fn set_spectrum_type(&mut self, spectrum_type: crate::meta::SpectrumType) {
        self.description_mut().set_spectrum_type(spectrum_type);
    }

    /// Compute and update the the total ion current, base peak, and m/z range for
    /// the spectrum based upon its current peak data.
    ///
    /// Uses [`RefPeakDataLevel::fetch_summaries`]
    fn update_summaries(&mut self) {
        let SpectrumSummary {
            tic,
            base_peak: bp,
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
    pub fn into_centroid<C>(self) -> Result<CentroidSpectrumType<C>, SpectrumConversionError>
    where
        C: BuildFromArrayMap + CentroidLike,
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
    pub fn into_spectrum<C, D>(self) -> Result<MultiLayerSpectrum<C, D>, SpectrumConversionError>
    where
        C: BuildFromArrayMap,
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
    {
        Ok(MultiLayerSpectrum::<C, D> {
            arrays: Some(self.arrays),
            description: self.description,
            ..Default::default()
        })
    }

    /// Apply a local denoising algorithm to the signal that is suitable for profile mode
    /// data using [`mzsignal::denoise`].
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
    /// [`MultiLayerSpectrum::pick_peaks_in_intervals`]
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

impl RawSpectrum {
    /// The method to access the spectrum description itself, which supplies
    /// the data for most other methods on this trait.
    pub fn description(&self) -> &SpectrumDescription {
        <RawSpectrum as SpectrumLike<CentroidPeak, DeconvolutedPeak>>::description(self)
    }

    /// The method to access the spectrum descript itself, mutably.
    pub fn description_mut(&mut self) -> &mut SpectrumDescription {
        <RawSpectrum as SpectrumLike<CentroidPeak, DeconvolutedPeak>>::description_mut(self)
    }

    /// Access the acquisition information for this spectrum.
    #[inline]
    pub fn acquisition(&self) -> &Acquisition {
        &self.description().acquisition
    }

    /// Access the precursor information, if it exists.
    #[inline]
    pub fn precursor(&self) -> Option<&Precursor> {
        let desc = self.description();
        desc.precursor.first()
    }

    /// Iterate over all precursors of the spectrum
    pub fn precursor_iter(&self) -> impl Iterator<Item = &Precursor> {
        let desc = self.description();
        desc.precursor.iter()
    }

    /// Mutably access the (first) precursor information, if it exists
    pub fn precursor_mut(&mut self) -> Option<&mut Precursor> {
        let desc = self.description_mut();
        desc.precursor.first_mut()
    }

    /// Iterate over all precursors of the spectrum mutably
    pub fn precursor_iter_mut(&mut self) -> impl Iterator<Item = &mut Precursor> {
        let desc = self.description_mut();
        desc.precursor.iter_mut()
    }

    /// A shortcut method to retrieve the scan start time of a spectrum
    #[inline]
    pub fn start_time(&self) -> f64 {
        let acq = self.acquisition();
        match acq.scans.first() {
            Some(evt) => evt.start_time,
            None => 0.0,
        }
    }

    /// Access the MS exponentiation level
    #[inline]
    pub fn ms_level(&self) -> u8 {
        self.description().ms_level
    }

    /// Access the native ID string for the spectrum
    #[inline]
    pub fn id(&self) -> &str {
        &self.description().id
    }

    /// Access the index of the spectrum in the source file
    #[inline]
    pub fn index(&self) -> usize {
        self.description().index
    }

    /// Access a description of how raw the signal is, whether a
    /// profile spectrum is available or only centroids are present.
    #[inline]
    pub fn signal_continuity(&self) -> SignalContinuity {
        self.description().signal_continuity
    }

    /// Access a description of the spectrum polarity
    #[inline]
    pub fn polarity(&self) -> ScanPolarity {
        self.description().polarity
    }

    /// Access the point measure of ion mobility associated with the scan if present. This is distinct from
    /// having a frame-level scan across the ion mobility dimension.
    pub fn ion_mobility(&self) -> Option<f64> {
        <RawSpectrum as SpectrumLike<CentroidPeak, DeconvolutedPeak>>::ion_mobility(self)
    }

    /// Check if this spectrum has a point measure of ion mobility. This is distinct from
    /// having a frame-level scan across the ion mobility dimension.
    pub fn has_ion_mobility(&self) -> bool {
        <RawSpectrum as SpectrumLike<CentroidPeak, DeconvolutedPeak>>::ion_mobility(self).is_some()
    }

    /// Retrieve the most processed representation of the mass spectrum's signal
    pub fn peaks(&'_ self) -> RefPeakDataLevel<'_, CentroidPeak, DeconvolutedPeak> {
        <RawSpectrum as SpectrumLike<CentroidPeak, DeconvolutedPeak>>::peaks(self)
    }

    pub fn into_peaks_and_description(
        self,
    ) -> (
        PeakDataLevel<CentroidPeak, DeconvolutedPeak>,
        SpectrumDescription,
    ) {
        <RawSpectrum as SpectrumLike<CentroidPeak, DeconvolutedPeak>>::into_peaks_and_description(
            self,
        )
    }

    /// Obtain a reference to the [`BinaryArrayMap`] if one is available for the peak
    /// information. This may not be the most refined version of the peak signal if
    /// it has been processed further.
    pub fn raw_arrays(&'_ self) -> Option<&'_ BinaryArrayMap> {
        <RawSpectrum as SpectrumLike<CentroidPeak, DeconvolutedPeak>>::raw_arrays(self)
    }

    /// Check if this spectrum has an ion mobility dimension/array. This is distinct from
    /// having a scan-level point measure of ion mobility.
    pub fn has_ion_mobility_dimension(&self) -> bool {
        self.raw_arrays()
            .map(|a| a.has_ion_mobility())
            .unwrap_or_default()
    }

    /// Compute and update the the total ion current, base peak, and m/z range for
    /// the spectrum based upon its current peak data.
    ///
    /// Uses [`RefPeakDataLevel::fetch_summaries`]
    pub fn update_summaries(&mut self) {
        <RawSpectrum as SpectrumLike<CentroidPeak, DeconvolutedPeak>>::update_summaries(self)
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

    fn into_peaks_and_description(self) -> (PeakDataLevel<C, D>, SpectrumDescription) {
        (PeakDataLevel::RawData(self.arrays), self.description)
    }
}

#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
/// Represents a spectrum that has been centroided into discrete m/z points, a
/// process also called "peak picking".
///
/// This type of spectrum represents data in exactly one format.
pub struct CentroidSpectrumType<C: CentroidLike> {
    /// The spectrum metadata describing acquisition conditions and details.
    pub description: SpectrumDescription,
    /// The picked centroid peaks, sorted by m/z in a fast searchable structure.
    pub peaks: MZPeakSetType<C>,
}

impl<C: CentroidLike> Default for CentroidSpectrumType<C> {
    fn default() -> Self {
        Self {
            description: Default::default(),
            peaks: PeakSetVec::empty(),
        }
    }
}

/// When the `mzsignal` feature is enabled, [`CentroidSpectrumType`] can be converted back into a profile-like
/// state using a peak shape model.
#[cfg(feature = "mzsignal")]
impl<C: CentroidLike + BuildArrayMapFrom + BuildFromArrayMap> CentroidSpectrumType<C> {
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

impl<C: CentroidLike> ParamDescribed for CentroidSpectrumType<C> {
    fn params(&self) -> &[crate::params::Param] {
        <SpectrumDescription as ParamDescribed>::params(&self.description)
    }

    fn params_mut(&mut self) -> &mut crate::params::ParamList {
        <SpectrumDescription as ParamDescribed>::params_mut(&mut self.description)
    }
}

/// [`CentroidSpectrumType`] implements [`SpectrumLike`] for `C`.
impl<C: CentroidLike> SpectrumLike<C> for CentroidSpectrumType<C> {
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

impl<C: CentroidLike> CentroidSpectrumType<C> {
    pub fn new(description: SpectrumDescription, peaks: MZPeakSetType<C>) -> Self {
        Self { description, peaks }
    }

    /// Convert a spectrum into a [`MultiLayerSpectrum`] on `C`, but preserve all information
    /// associated with this spectrum as-is.
    pub fn into_spectrum<D>(self) -> Result<MultiLayerSpectrum<C, D>, SpectrumConversionError>
    where
        D: DeconvolutedCentroidLike,
    {
        let val = MultiLayerSpectrum::<C, D> {
            peaks: Some(self.peaks),
            description: self.description,
            ..Default::default()
        };
        Ok(val)
    }

    /// Convert a [`CentroidSpectrumType`] into a [`MultiLayerSpectrum`] over any peak type, but
    /// recode the peak list as [`DataArray`], potentially losing information.
    pub fn reinterpret<C1, D1>(self) -> Result<MultiLayerSpectrum<C1, D1>, SpectrumConversionError>
    where
        C1: CentroidLike + BuildArrayMapFrom + BuildFromArrayMap,
        D1: DeconvolutedCentroidLike + BuildArrayMapFrom + BuildFromArrayMap,
        C: BuildArrayMapFrom,
    {
        let arrays = C::as_arrays(self.peaks.as_slice());
        Ok(MultiLayerSpectrum::new(
            self.description,
            Some(arrays),
            None,
            None,
        ))
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

    fn iter<'a>(&'a self) -> impl Iterator<Item = &'a C>
    where
        C: 'a,
    {
        self.peaks.iter()
    }
}

/// Represents a spectrum that has been centroided, deisotoped, and charge state deconvolved.
///
/// This type of spectrum represents data in exactly one format.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DeconvolutedSpectrumType<D: DeconvolutedCentroidLike> {
    /// The spectrum metadata describing acquisition conditions and details.
    pub description: SpectrumDescription,
    /// The deisotoped and charge state deconvolved peaks, sorted by neutral mass
    /// in a fast searchable structure.
    pub deconvoluted_peaks: MassPeakSetType<D>,
}

impl<D: DeconvolutedCentroidLike> Default for DeconvolutedSpectrumType<D> {
    fn default() -> Self {
        Self {
            description: Default::default(),
            deconvoluted_peaks: PeakSetVec::empty(),
        }
    }
}

impl<D: DeconvolutedCentroidLike> DeconvolutedSpectrumType<D> {
    pub fn new(description: SpectrumDescription, deconvoluted_peaks: MassPeakSetType<D>) -> Self {
        Self {
            description,
            deconvoluted_peaks,
        }
    }

    /// Convert a spectrum into a [`MultiLayerSpectrum`] on `D`, but preserve all information
    /// associated with this spectrum as-is.
    pub fn into_spectrum<C>(self) -> Result<MultiLayerSpectrum<C, D>, SpectrumConversionError>
    where
        C: CentroidLike,
    {
        let val = MultiLayerSpectrum::<C, D> {
            deconvoluted_peaks: Some(self.deconvoluted_peaks),
            description: self.description,
            ..Default::default()
        };
        Ok(val)
    }

    /// Convert a [`DeconvolutedSpectrumType`] into a [`MultiLayerSpectrum`] over any peak type, but
    /// recode the peak list as [`DataArray`], potentially losing information.
    pub fn reinterpret<C1, D1>(self) -> MultiLayerSpectrum<C1, D1>
    where
        C1: CentroidLike + BuildArrayMapFrom + BuildFromArrayMap,
        D1: DeconvolutedCentroidLike + BuildArrayMapFrom + BuildFromArrayMap,
        D: BuildArrayMapFrom,
    {
        let arrays = D::as_arrays(self.deconvoluted_peaks.as_slice());
        MultiLayerSpectrum::new(self.description, Some(arrays), None, None)
    }
}

impl<D: DeconvolutedCentroidLike> ParamDescribed for DeconvolutedSpectrumType<D> {
    fn params(&self) -> &[crate::params::Param] {
        <SpectrumDescription as ParamDescribed>::params(&self.description)
    }

    fn params_mut(&mut self) -> &mut crate::params::ParamList {
        <SpectrumDescription as ParamDescribed>::params_mut(&mut self.description)
    }
}

/// [`DeconvolutedSpectrumType`] implements [`SpectrumLike`] for `D`.
impl<D: DeconvolutedCentroidLike> SpectrumLike<CentroidPeak, D> for DeconvolutedSpectrumType<D> {
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

    fn iter<'a>(&'a self) -> impl Iterator<Item = &'a D>
    where
        D: 'a,
    {
        self.deconvoluted_peaks.iter()
    }
}

#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
/// Represent a spectrum with multiple layers of representation of the
/// peak data.
///
/// This type is useful because it permits us to hold spectra in incrementally
/// more processed representations without loss of information.
pub struct MultiLayerSpectrum<
    C: CentroidLike = CentroidPeak,
    D: DeconvolutedCentroidLike = DeconvolutedPeak,
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

impl<C: CentroidLike, D: DeconvolutedCentroidLike> Default for MultiLayerSpectrum<C, D> {
    fn default() -> Self {
        Self {
            description: Default::default(),
            arrays: None,
            peaks: None,
            deconvoluted_peaks: None,
        }
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> ParamDescribed for MultiLayerSpectrum<C, D> {
    fn params(&self) -> &[crate::params::Param] {
        <SpectrumDescription as ParamDescribed>::params(&self.description)
    }

    fn params_mut(&mut self) -> &mut crate::params::ParamList {
        <SpectrumDescription as ParamDescribed>::params_mut(&mut self.description)
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> SpectrumLike<C, D> for MultiLayerSpectrum<C, D> {
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

impl<C: CentroidLike, D: DeconvolutedCentroidLike> MultiLayerSpectrum<C, D>
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

    /// Construct a [`MultiLayerSpectrum`] with `description` and the provided raw data arrays
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

    /// Construct a [`MultiLayerSpectrum`] with `description` and *no* peak data of any kind
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

    /// Converts an arbitrary [`SpectrumLike`] type into a [`MultiLayerSpectrum`] with the same
    /// peak types.
    ///
    /// This is a thin wrapper around [`SpectrumLike::into_peaks_and_description`] and [`MultiLayerSpectrum::from_peaks_data_levels_and_description`]
    pub fn from_spectrum_like<S: SpectrumLike<C, D>>(spectrum: S) -> Self {
        let (peaks, description) = spectrum.into_peaks_and_description();
        Self::from_peaks_data_levels_and_description(peaks, description)
    }

    /// Convert a [`MultiLayerSpectrum`] with one set of peak types to another with any peak type, but
    /// recode the peak list as [`DataArray`] as an intermediary, potentially losing information.
    ///
    /// ## Note
    /// Peak data is not actually reconstructed from the intermediate data arrays. If that is desired,
    /// call [`MultiLayerSpectrum::try_build_peaks`]
    pub fn reinterpret<C1, D1>(self) -> MultiLayerSpectrum<C1, D1>
    where
        C1: CentroidLike + BuildArrayMapFrom + BuildFromArrayMap,
        D1: DeconvolutedCentroidLike + BuildArrayMapFrom + BuildFromArrayMap,
        C: BuildArrayMapFrom,
        D: BuildArrayMapFrom,
    {
        let arrays = if let Some(peaks) = self.deconvoluted_peaks {
            D::as_arrays(peaks.as_slice())
        } else if let Some(peaks) = self.peaks {
            C::as_arrays(peaks.as_slice())
        } else if let Some(arrays) = self.arrays {
            arrays
        } else {
            BinaryArrayMap::new()
        };

        MultiLayerSpectrum::new(self.description, Some(arrays), None, None)
    }

    #[cfg(feature = "mzsignal")]
    /// Apply a local denoising algorithm to the signal that is suitable for profile mode
    /// data using [`mzsignal::denoise`].
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
    /// When `mzsignal` is available, use [`mzsignal::reprofile::PeakSetReprofiler`] to create
    /// new profile m/z and intensity [`DataArray`]s from [`MultiLayerSpectrum::peaks`], updating
    /// [`MultiLayerSpectrum::arrays`].
    ///
    /// **NOTE**: This does not modify [`SpectrumLike::signal_continuity`], and if it is desired
    /// you should update this manually to [`SignalContinuity::Profile`].
    ///
    /// # Arguments
    /// - `dx`: The m/z grid spacing that the profile signal will be projected onto. The
    ///   best value to use here depends upon the spectrum's resolution, but when resources
    ///   are not an issue between `0.005` and `0.001` are is a choice. This may produce a lot
    ///   of wasted space however.
    /// - `fwhm`: The uniform peak full-width-at-half-max to assume for all peaks. It controls
    ///   the broadness of the theoretical gaussian profile created around each centroid peak.
    ///   If non-uniform shapes are needed, directly use [`PeakSetReprofiler`](mzsignal::reprofile::PeakSetReprofiler).
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
                pair.mz_array.iter().flat_map(|i| i.to_le_bytes()).collect(),
            ));

            arrays.add(DataArray::wrap(
                &ArrayType::IntensityArray,
                BinaryDataArrayType::Float32,
                pair.intensity_array
                    .iter()
                    .flat_map(|i| i.to_le_bytes())
                    .collect(),
            ));
            Ok(())
        } else {
            Err(SpectrumConversionError::NoPeakData.into())
        }
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> MultiLayerSpectrum<C, D>
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
    pub fn try_build_peaks(
        &mut self,
    ) -> Result<RefPeakDataLevel<'_, C, D>, SpectrumConversionError> {
        if self.peaks.is_some() || self.deconvoluted_peaks.is_some() {
            return Ok(self.peaks());
        }
        if matches!(self.signal_continuity(), SignalContinuity::Centroid) {
            if let Some(arrays) = self.arrays.as_ref() {
                let peak_data: PeakDataLevel<C, D> = PeakDataLevel::try_from(arrays)?;
                match peak_data {
                    PeakDataLevel::Missing => Ok(RefPeakDataLevel::Missing),
                    PeakDataLevel::RawData(_) => panic!("not possible"),
                    PeakDataLevel::Centroid(peak_set_vec) => {
                        self.peaks = Some(peak_set_vec);
                        Ok(self.peaks())
                    }
                    PeakDataLevel::Deconvoluted(peak_set_vec) => {
                        self.deconvoluted_peaks = Some(peak_set_vec);
                        Ok(self.peaks())
                    }
                }
            } else {
                Ok(self.peaks())
            }
        } else {
            Ok(self.peaks())
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

/// When [`mzsignal`] is available, [`MultiLayerSpectrum`] supports in-place signal processing operations.
///
/// The peak picking steps need to convert an [`mzsignal::FittedPeak`] into `C`. This is trivial for [`CentroidPeak`]
/// and [`mzsignal::FittedPeak`] itself.
#[cfg(feature = "mzsignal")]
impl<C: CentroidLike + From<FittedPeak>, D: DeconvolutedCentroidLike> MultiLayerSpectrum<C, D> {
    /// Using a pre-configured [`mzsignal::PeakPicker`](mzsignal::peak_picker::PeakPicker) to pick peaks
    /// from `arrays`'s signal, populating the `peaks` field.
    ///
    /// If [`SpectrumLike::signal_continuity`] returns [`SignalContinuity::Centroid`], then no filtering is
    /// performed and all points are directly converted into picked peaks.
    ///
    /// **NOTE**: This does not modify [`SpectrumLike::signal_continuity`], and if it is desired
    /// you should update this manually to [`SignalContinuity::Centroid`].
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
    ///
    /// **NOTE**: This does not modify [`SpectrumLike::signal_continuity`], and if it is desired
    /// you should update this manually to [`SignalContinuity::Centroid`].
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

impl<C: CentroidLike, D: DeconvolutedCentroidLike> TryFrom<MultiLayerSpectrum<C, D>>
    for CentroidSpectrumType<C>
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

impl<C: CentroidLike, D: DeconvolutedCentroidLike> From<CentroidSpectrumType<C>>
    for MultiLayerSpectrum<C, D>
{
    fn from(spectrum: CentroidSpectrumType<C>) -> MultiLayerSpectrum<C, D> {
        spectrum.into_spectrum().unwrap()
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> From<RawSpectrum> for MultiLayerSpectrum<C, D>
where
    C: BuildFromArrayMap,
    D: BuildFromArrayMap,
{
    fn from(spectrum: RawSpectrum) -> MultiLayerSpectrum<C, D> {
        spectrum.into_spectrum().unwrap()
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> From<MultiLayerSpectrum<C, D>> for RawSpectrum
where
    C: BuildFromArrayMap + BuildArrayMapFrom,
    D: BuildFromArrayMap + BuildArrayMapFrom,
{
    fn from(spectrum: MultiLayerSpectrum<C, D>) -> RawSpectrum {
        spectrum.into_raw().unwrap()
    }
}

impl<C: CentroidLike> TryFrom<RawSpectrum> for CentroidSpectrumType<C>
where
    C: BuildFromArrayMap,
{
    type Error = SpectrumConversionError;

    fn try_from(spectrum: RawSpectrum) -> Result<CentroidSpectrumType<C>, Self::Error> {
        spectrum.into_centroid()
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> TryFrom<MultiLayerSpectrum<C, D>>
    for DeconvolutedSpectrumType<D>
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

#[cfg(all(test, feature = "mzml"))]
mod test {
    use std::io;

    use super::*;
    use crate::io::DetailLevel;
    use crate::io::mzml::MzMLReader;
    use crate::prelude::*;

    #[test_log::test]
    fn test_peakdata_lazy() -> io::Result<()> {
        let mut reader = MzMLReader::open_path("./test/data/small.mzML")?;
        reader.detail_level = DetailLevel::Lazy;
        let spec = reader.get_spectrum_by_index(0).unwrap();

        let peaks = spec.peaks();
        let n = peaks.len();
        assert_eq!(n, 19913);

        let iter = peaks.iter();
        let data: Vec<_> = iter.collect();

        let n2 = data.len();
        assert_eq!(n2, 19913);

        let p1 = peaks.get(5000);
        let p2 = data.get(5000).cloned();
        assert_eq!(p1, p2);
        Ok(())
    }

    macro_rules! behaviors {
        ($spec:ident) => {
            assert_eq!($spec.id(), "controllerType=0 controllerNumber=1 scan=10014");
            assert_eq!($spec.start_time(), 22.12829);
            assert_eq!($spec.ms_level(), 1);
            assert_eq!($spec.polarity(), ScanPolarity::Positive);
            assert!($spec.precursor().is_none());
            assert_eq!($spec.precursor_iter().count(), 0);
            assert!(!$spec.has_ion_mobility());
            assert!(!$spec.has_ion_mobility_dimension());
            assert!(!$spec.params().is_empty());
        };
    }

    #[allow(unused)]
    fn test_spectrum_behavior<T: SpectrumLike>(spec: &T) {
        assert_eq!(
            spec.spectrum_type(),
            Some(crate::meta::SpectrumType::MS1Spectrum)
        );
        behaviors!(spec);
    }

    #[test_log::test]
    fn test_peakdata() -> io::Result<()> {
        let mut reader = MzMLReader::open_path("./test/data/small.mzML")?;
        let spec = reader.get_spectrum_by_index(0).unwrap();

        let peaks = spec.peaks();
        let n = peaks.len();
        assert_eq!(n, 19913);

        let iter = peaks.iter();
        let data: Vec<_> = iter.collect();

        let n2 = data.len();
        assert_eq!(n2, 19913);

        let p1 = peaks.get(5000);
        let p2 = data.get(5000).cloned();
        assert_eq!(p1, p2);
        Ok(())
    }

    #[cfg(feature = "mzsignal")]
    #[test_log::test]
    fn test_profile_read() {
        let mut reader = MzMLReader::open_path("./test/data/three_test_scans.mzML")
            .expect("Failed to open test file");
        reader.reset();
        let mut scan = reader.next().expect("Failed to read spectrum");
        assert_eq!(scan.signal_continuity(), SignalContinuity::Profile);
        test_spectrum_behavior(&scan);
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

        scan.pick_peaks_with(&PeakPicker::default()).unwrap();
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

        let mut cent_scan = scan.clone().into_centroid().unwrap();
        assert_eq!(cent_scan.signal_continuity(), SignalContinuity::Centroid);
        assert_eq!(SpectrumLike::index(&cent_scan), 0);
        behaviors!(cent_scan);
        cent_scan.update_summaries();

        let tmp: Spectrum = cent_scan.into_spectrum().unwrap();
        let raw_scan = tmp.into_raw().unwrap();
        behaviors!(raw_scan);
        assert_eq!(raw_scan.index(), 0);
        assert_eq!(raw_scan.signal_continuity(), SignalContinuity::Centroid);
    }

    #[cfg(feature = "mzsignal")]
    #[test]
    fn test_reprofile() {
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

        let mut duplicate = scan.clone();
        duplicate.reprofile_with_shape(0.001, 0.01).unwrap();
        duplicate.peaks = None;
        duplicate.pick_peaks(1.0).unwrap();
        let peak = duplicate.peaks.as_ref().unwrap().base_peak().unwrap();
        eprintln!("{}", peak);
    }
}
