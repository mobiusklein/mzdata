//! Represent the collection of attributes and data that compose a single mass spectrum.
//!
//! Because a mass spectrum may be obtained from sources with varying levels of detail,
//! several alternative structures are provided with a common set of trait-based methods
//! to unify access:
//!
//! 1. [`RawSpectrum`] for representing a spectrum that has not been decoded into distinct
//!        peaks yet, but whose data may be continuous or discrete.
//! 2. [`CentroidSpectrum`] for representing spectra from sources which are guaranteed to
//!        be pre-centroided, like those from MGF files or other simple text representations.
//! 3. [`Spectrum`] for representing a multi-layer representation of a spectrum where both
//!        raw data and a distinct peak list are available.
//!
//! These structures all implement the [`SpectrumBehavior`] trait
use std::borrow;
use std::convert::TryFrom;

use mzpeaks::prelude::*;
use mzpeaks::{
    CentroidLike, DeconvolutedCentroidLike, DeconvolutedPeakSet, MZPeakSetType, MassPeakSetType,
    PeakSet,
};
use mzpeaks::{CentroidPeak, DeconvolutedPeak, MassErrorType};

use mzsignal::denoise::{denoise, DenoisingError};
use mzsignal::peak_picker::{PeakFitType, PeakPicker, PeakPickerError};

use crate::params::ParamList;
use crate::spectrum::scan_properties::{
    Acquisition, Precursor, ScanPolarity, SignalContinuity, SpectrumDescription,
};
use crate::spectrum::signal::{ArrayType, BinaryArrayMap, BinaryDataArrayType};

pub trait CentroidPeakAdapting: CentroidLike + Default + From<CentroidPeak> {}
impl<C: CentroidLike + Default + From<CentroidPeak>> CentroidPeakAdapting for C {}
pub trait DeconvolutedPeakAdapting:
    DeconvolutedCentroidLike + Default + From<DeconvolutedPeak>
{
}
impl<D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak>> DeconvolutedPeakAdapting
    for D
{
}

#[derive(Debug)]
/// An variant for dispatching to different strategies of computing
/// common statistics of different levels of peak data.
pub enum PeakDataLevel<'lifespan, C: CentroidLike, D: DeconvolutedCentroidLike> {
    Missing,
    RawData(&'lifespan BinaryArrayMap),
    Centroid(&'lifespan MZPeakSetType<C>),
    Deconvoluted(&'lifespan MassPeakSetType<D>),
}

impl<'lifespan, C: CentroidLike, D: DeconvolutedCentroidLike> PeakDataLevel<'lifespan, C, D> {
    /// Compute the base peak of a spectrum
    pub fn base_peak(&self) -> (usize, f64, f32) {
        match self {
            PeakDataLevel::Missing => (0, 0.0, 0.0),
            PeakDataLevel::RawData(arrays) => {
                let intensities = arrays.intensities();
                let result = intensities
                    .iter()
                    .enumerate()
                    .max_by(|ia, ib| ia.1.partial_cmp(ib.1).unwrap());
                if let Some((i, inten)) = result {
                    (i, arrays.mzs()[i], *inten)
                } else {
                    (0, 0.0, 0.0)
                }
            }
            PeakDataLevel::Centroid(peaks) => {
                let result = peaks
                    .iter()
                    .enumerate()
                    .max_by(|ia, ib| ia.1.intensity().partial_cmp(&ib.1.intensity()).unwrap());
                if let Some((i, peak)) = result {
                    (i, peak.coordinate(), peak.intensity())
                } else {
                    (0, 0.0, 0.0)
                }
            }
            PeakDataLevel::Deconvoluted(peaks) => {
                let result = peaks
                    .iter()
                    .enumerate()
                    .max_by(|ia, ib| ia.1.intensity().partial_cmp(&ib.1.intensity()).unwrap());
                if let Some((i, peak)) = result {
                    (i, peak.coordinate(), peak.intensity())
                } else {
                    (0, 0.0, 0.0)
                }
            }
        }
    }

    pub fn tic(&self) -> f32 {
        match self {
            PeakDataLevel::Missing => 0.0,
            PeakDataLevel::RawData(arrays) => {
                let intensities = arrays.intensities();
                intensities.iter().sum()
            }
            PeakDataLevel::Centroid(peaks) => peaks.iter().map(|p| p.intensity()).sum(),
            PeakDataLevel::Deconvoluted(peaks) => peaks.iter().map(|p| p.intensity()).sum(),
        }
    }

    pub fn search(
        &self,
        query: f64,
        error_tolerance: f64,
        error_type: MassErrorType,
    ) -> Option<usize> {
        match self {
            PeakDataLevel::Missing => None,
            PeakDataLevel::RawData(arrays) => arrays.search(query, error_tolerance, error_type),
            PeakDataLevel::Centroid(peaks) => peaks.search(query, error_tolerance, error_type),
            PeakDataLevel::Deconvoluted(peaks) => peaks.search(query, error_tolerance, error_type),
        }
    }
}

/// A trait for providing a uniform delegated access to spectrum metadata
pub trait SpectrumBehavior<
    C: CentroidLike = CentroidPeak,
    D: DeconvolutedCentroidLike = DeconvolutedPeak,
>
{
    /// The method to access the spectrum description itself, which supplies
    /// the data for most other methods on this trait.
    fn description(&self) -> &SpectrumDescription;

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

    /// A shortcut method to retrieve the scan start time
    /// of a spectrum.
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

    /// Retrieve the most processed representation of the mass spectrum's
    /// signal
    fn peaks(&'_ self) -> PeakDataLevel<'_, C, D>;
}

#[derive(Default, Debug, Clone)]
/// Represents a spectrum that hasn't been processed yet, with only
/// data arrays, potentially no discrete peaks.
pub struct RawSpectrum {
    /// The spectrum metadata describing acquisition conditions and details.
    pub description: SpectrumDescription,
    /// The data arrays describing the m/z, intensity, and potentially other
    /// measured properties
    pub arrays: BinaryArrayMap,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpectrumConversionError {
    MZIntensityArraySizeMismatch,
    NotDeconvoluted,
    NotCentroided,
    NoPeakData,
}

impl std::fmt::Display for SpectrumConversionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::error::Error for SpectrumConversionError {}


#[derive(Debug, Clone)]
pub enum SpectrumProcessingError {
    DenoisingError(DenoisingError),
    PeakPickerError(PeakPickerError),
    SpectrumConversionError(SpectrumConversionError),
}

impl std::fmt::Display for SpectrumProcessingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::error::Error for SpectrumProcessingError {}

impl<'transient, 'lifespan: 'transient> RawSpectrum {
    /// Convert a spectrum into a [`CentroidSpectrum`]
    pub fn into_centroid(self) -> Result<CentroidSpectrum, SpectrumConversionError> {
        if !matches!(
            self.description.signal_continuity,
            SignalContinuity::Centroid
        ) {
            return Err(SpectrumConversionError::NotCentroided);
        }

        let mz_array = self.mzs();
        let intensity_array = self.intensities();

        if mz_array.len() != intensity_array.len() {
            return Err(SpectrumConversionError::MZIntensityArraySizeMismatch);
        }

        let mut peaks = PeakSet::empty();
        for (mz, intensity) in mz_array.iter().zip(intensity_array.iter()) {
            peaks.push(CentroidPeak {
                mz: *mz,
                intensity: *intensity,
                ..CentroidPeak::default()
            });
        }
        let mut centroid = CentroidSpectrum {
            description: self.description,
            peaks,
        };
        centroid.description.signal_continuity = SignalContinuity::Centroid;

        Ok(centroid)
    }

    pub fn mzs(&'lifespan self) -> borrow::Cow<'transient, [f64]> {
        self.arrays.mzs()
    }

    pub fn intensities(&'lifespan self) -> borrow::Cow<'transient, [f32]> {
        self.arrays.intensities()
    }

    /// Convert a spectrum into a [`Spectrum`]
    pub fn into_spectrum<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>(
        self,
    ) -> Result<MultiLayerSpectrum<C, D>, SpectrumConversionError> {
        Ok(MultiLayerSpectrum::<C, D> {
            arrays: Some(self.arrays),
            description: self.description,
            ..Default::default()
        })
    }

    pub fn denoise(&mut self, scale: f32) -> Result<(), SpectrumProcessingError> {
        let mut intensities_copy = self.arrays.intensities().into_owned();
        let mz_array = self.arrays.mzs();
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

    pub fn pick_peaks_with_into(
        self,
        peak_picker: &PeakPicker,
    ) -> Result<MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak>, SpectrumProcessingError> {
        let mut result = self.into_spectrum().unwrap();
        match result.pick_peaks_with(peak_picker) {
            Ok(_) => Ok(result),
            Err(err) => Err(err),
        }
    }

    pub fn pick_peaks_into(
        self,
        signal_to_noise_threshold: f32,
        fit_type: PeakFitType,
    ) -> Result<MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak>, SpectrumProcessingError> {
        let peak_picker = PeakPicker {
            fit_type,
            signal_to_noise_threshold,
            ..Default::default()
        };
        self.pick_peaks_with_into(&peak_picker)
    }
}

impl<'lifespan> SpectrumBehavior for RawSpectrum {
    #[inline]
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }

    fn peaks(&'_ self) -> PeakDataLevel<'_, CentroidPeak, DeconvolutedPeak> {
        PeakDataLevel::RawData(&self.arrays)
    }
}

#[derive(Default, Debug, Clone)]
/// Represents a spectrum that has been centroided
pub struct CentroidSpectrumType<C: CentroidLike + Default> {
    /// The spectrum metadata describing acquisition conditions and details.
    pub description: SpectrumDescription,
    /// The picked centroid peaks
    pub peaks: MZPeakSetType<C>,
}

impl<'lifespan, C: CentroidLike + Default> SpectrumBehavior<C> for CentroidSpectrumType<C> {
    #[inline]
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }

    fn peaks(&'_ self) -> PeakDataLevel<'_, C, DeconvolutedPeak> {
        PeakDataLevel::Centroid(&self.peaks)
    }
}

impl<C: CentroidLike + Default> CentroidSpectrumType<C> {
    /// Convert a spectrum into a [`Spectrum`]
    pub fn into_spectrum<D: DeconvolutedPeakAdapting>(
        self,
    ) -> Result<MultiLayerSpectrum<C, D>, SpectrumConversionError> {
        let val = MultiLayerSpectrum::<C, D> {
            peaks: Some(self.peaks),
            description: self.description,
            ..Default::default()
        };
        Ok(val)
    }
}

pub type CentroidSpectrum = CentroidSpectrumType<CentroidPeak>;

#[derive(Default, Debug, Clone)]
pub struct DeconvolutedSpectrumType<D: DeconvolutedCentroidLike + Default> {
    /// The spectrum metadata describing acquisition conditions and details.
    pub description: SpectrumDescription,
    /// The deisotoped and charge state deconvolved peaks
    pub deconvoluted_peaks: MassPeakSetType<D>,
}

impl<'lifespan, D: DeconvolutedCentroidLike + Default> SpectrumBehavior<CentroidPeak, D>
    for DeconvolutedSpectrumType<D>
{
    #[inline]
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }

    fn peaks(&'_ self) -> PeakDataLevel<'_, CentroidPeak, D> {
        PeakDataLevel::Deconvoluted(&self.deconvoluted_peaks)
    }
}

pub type DeconvolutedSpectrum = DeconvolutedSpectrumType<DeconvolutedPeak>;

#[derive(Default, Debug, Clone)]
/// Represent a spectrum with multiple layers of representation of the
/// peak data.
pub struct MultiLayerSpectrum<C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default> {
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

impl<'lifespan, C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default>
    SpectrumBehavior<C, D> for MultiLayerSpectrum<C, D>
{
    #[inline]
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }

    fn peaks(&'_ self) -> PeakDataLevel<'_, C, D> {
        if let Some(peaks) = &self.peaks {
            PeakDataLevel::Centroid(peaks)
        } else if let Some(arrays) = &self.arrays {
            PeakDataLevel::RawData(arrays)
        } else {
            PeakDataLevel::Missing
        }
    }
}

impl<'lifespan, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MultiLayerSpectrum<C, D> {
    /// Convert a spectrum into a [`CentroidSpectrumType<C>`]
    pub fn into_centroid(self) -> Result<CentroidSpectrumType<C>, SpectrumConversionError> {
        if let Some(peaks) = self.peaks {
            let mut result = CentroidSpectrumType::<C> {
                peaks,
                description: self.description,
            };
            result.description.signal_continuity = SignalContinuity::Centroid;
            return Ok(result);
        } else {
            if self.signal_continuity() == SignalContinuity::Centroid {
                if let Some(arrays) = &self.arrays {
                    let peaks = MZPeakSetType::<C>::from(arrays);
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
            let arrays = BinaryArrayMap::from(DeconvolutedPeakSet::wrap(
                peaks.into_iter().map(|p| p.as_centroid()).collect(),
            ));

            let mut result = RawSpectrum {
                arrays,
                description: self.description,
            };
            result.description.signal_continuity = SignalContinuity::Centroid;
            Ok(result)
        } else if let Some(peaks) = self.peaks {
            let arrays = BinaryArrayMap::from(PeakSet::wrap(
                peaks.into_iter().map(|p| p.as_centroid()).collect(),
            ));

            let mut result = RawSpectrum {
                arrays,
                description: self.description,
            };
            result.description.signal_continuity = SignalContinuity::Centroid;
            Ok(result)
        } else {
            let arrays = BinaryArrayMap::from(PeakSet::empty());
            Ok(RawSpectrum {
                arrays,
                description: self.description,
            })
        }
    }

    pub fn pick_peaks_with(
        &mut self,
        peak_picker: &PeakPicker,
    ) -> Result<(), SpectrumProcessingError> {
        if let Some(arrays) = &self.arrays {
            let mz_array = arrays.mzs();
            let intensity_array = arrays.intensities();
            let mut acc = Vec::new();
            match peak_picker.discover_peaks(&mz_array, &intensity_array, &mut acc) {
                Ok(_) => {
                    let peaks: MZPeakSetType<C> =
                        acc.into_iter().map(|p| C::from(p.into())).collect();
                    self.peaks = Some(peaks);
                    Ok(())
                }
                Err(err) => Err(SpectrumProcessingError::PeakPickerError(err)),
            }
        } else {
            Err(SpectrumProcessingError::SpectrumConversionError(
                SpectrumConversionError::NoPeakData,
            ))
        }
    }

    pub fn pick_peaks(
        &mut self,
        signal_to_noise_threshold: f32,
        fit_type: PeakFitType,
    ) -> Result<(), SpectrumProcessingError> {
        let peak_picker = PeakPicker {
            fit_type,
            signal_to_noise_threshold,
            ..Default::default()
        };
        self.pick_peaks_with(&peak_picker)
    }

    pub fn denoise(&mut self, scale: f32) -> Result<(), SpectrumProcessingError> {
        match &mut self.arrays {
            Some(arrays) => {
                let mut intensities_copy = arrays.intensities().into_owned();
                let mz_array = arrays.mzs();
                match denoise(&mz_array, &mut intensities_copy, scale) {
                    Ok(_) => {
                        let view = arrays.get_mut(&ArrayType::IntensityArray).unwrap();
                        view.store_as(BinaryDataArrayType::Float32)
                            .expect("Failed to reformat intensity array");
                        view.update_buffer(&intensities_copy)
                            .expect("Failed to update intensity array buffer");
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
}

impl<'lifespan, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    TryFrom<MultiLayerSpectrum<C, D>> for CentroidSpectrumType<C>
{
    type Error = SpectrumConversionError;

    fn try_from(
        spectrum: MultiLayerSpectrum<C, D>,
    ) -> Result<CentroidSpectrumType<C>, Self::Error> {
        spectrum.into_centroid()
    }
}

pub type Spectrum = MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak>;

impl<C: CentroidPeakAdapting> From<CentroidSpectrumType<C>>
    for MultiLayerSpectrum<C, DeconvolutedPeak>
{
    fn from(spectrum: CentroidSpectrumType<C>) -> MultiLayerSpectrum<C, DeconvolutedPeak> {
        spectrum.into_spectrum().unwrap()
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> From<RawSpectrum>
    for MultiLayerSpectrum<C, D>
{
    fn from(spectrum: RawSpectrum) -> MultiLayerSpectrum<C, D> {
        spectrum.into_spectrum().unwrap()
    }
}

impl From<Spectrum> for RawSpectrum {
    fn from(spectrum: Spectrum) -> RawSpectrum {
        spectrum.into_raw().unwrap()
    }
}

impl TryFrom<RawSpectrum> for CentroidSpectrum {
    type Error = SpectrumConversionError;

    fn try_from(spectrum: RawSpectrum) -> Result<CentroidSpectrum, Self::Error> {
        spectrum.into_centroid()
    }
}

impl TryFrom<Spectrum> for DeconvolutedSpectrum {
    type Error = SpectrumConversionError;

    fn try_from(spectrum: Spectrum) -> Result<Self, Self::Error> {
        if spectrum.signal_continuity() == SignalContinuity::Profile {
            Err(SpectrumConversionError::NotCentroided)
        } else if let Some(arrays) = &spectrum.arrays {
            if arrays.has_array(&ArrayType::ChargeArray) {
                let peaks: DeconvolutedPeakSet = DeconvolutedPeakSet::from(arrays);
                return Ok(DeconvolutedSpectrum {
                    description: spectrum.description,
                    deconvoluted_peaks: peaks,
                });
            }
            Err(Self::Error::NotDeconvoluted)
        } else {
            Err(Self::Error::NoPeakData)
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::io::mzml::MzMLReader;
    use crate::io::prelude::*;

    use mzsignal::arrayops::ArrayPair;
    use mzsignal::plot::*;
    use mzsignal::plot::{SVGBuilder, RED};

    #[test]
    fn test_profile_read() {
        let mut reader = MzMLReader::open_path("./test/data/three_test_scans.mzML")
            .expect("Failed to open test file");
        let mut scan = reader.next().unwrap();
        assert_eq!(scan.signal_continuity(), SignalContinuity::Profile);
        assert_eq!(scan.ms_level(), 1);
        assert_eq!(scan.polarity(), ScanPolarity::Positive);
        assert!(matches!(scan.precursor(), None));

        match scan.pick_peaks(1.0, PeakFitType::Quadratic) {
            Err(err) => {
                panic!("Should not have an error! {:?}", err);
            }
            Ok(_) => {}
        }

        if let Some(peaks) = &scan.peaks {
            assert_eq!(peaks.len(), 2097);
            let n = peaks.len();
            for i in 1..n {
                let diff = peaks[i].get_index() - peaks[i - 1].get_index();
                assert_eq!(diff, 1);
            }

            let mut builder = SVGBuilder::default();
            let arrays = &scan.arrays.unwrap();
            let ser = ArrayPair::new(arrays.mzs(), arrays.intensities());
            let mut ser2 = SpectrumSeries::from(peaks.iter());
            ser2.color(RED.mix(1.0));
            builder
                .path("./test/data/0.svg")
                .size(1028, 512)
                .add_series(ser2)
                .add_series(&ser)
                .xlim(562f64, 565f64)
                .draw()
                .expect("Failed to draw");

            peaks
                .has_peak(562.741, 3f64, MassErrorType::PPM)
                .expect("Expected to find peak");
            peaks
                .has_peak(563.240, 3f64, MassErrorType::PPM)
                .expect("Expected to find peak");
            let p = peaks
                .has_peak(563.739, 1f64, MassErrorType::PPM)
                .expect("Expected to find peak");
            assert!((p.mz() - 563.739).abs() < 1e-3)
        }
    }
}
