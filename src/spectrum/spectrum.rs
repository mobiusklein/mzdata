/*!
Represent the collection of attributes and data that compose a single mass spectrum.

Because a mass spectrum may be obtained from sources with varying levels of detail,
several alternative structures are provided with a common set of trait-based methods
to unify access:

1. [`RawSpectrum`] for representing a spectrum that has not been decoded into distinct
       peaks yet, but whose data may be continuous or discrete.
2. [`CentroidSpectrum`] for representing spectra from sources which are guaranteed to
       be pre-centroided, like those from MGF files or other simple text representations.
3. [`Spectrum`] for representing a multi-layer representation of a spectrum where both
       raw data and a distinct peak list are available.

These structures all implement the [`SpectrumBehavior`] trait

The [`SpectrumBehavior`] trait is included in the crate prelude, and gives the caller
read-only access to components that describe a spectrum's metadata, and
```rust
use std::fs::File;
use mzdata::prelude::*;
use mzpeaks::{Tolerance, prelude::*};
use mzdata::io::MzMLReader;
use mzdata::spectrum::{SignalContinuity};

let reader = MzMLReader::new(File::open("./test/data/small.mzML").unwrap());
for spectrum in reader {
    println!("Scan {} => BP {}", spectrum.id(), spectrum.peaks().base_peak().mz);
    if spectrum.signal_continuity() < SignalContinuity::Profile {
        let peak_picked = spectrum.into_centroid().unwrap();
        println!("Matches for 579.155: {:?}",
                 peak_picked.peaks.all_peaks_for(
                     579.155, Tolerance::Da(0.02)));
    }
}
```

*/
use std::borrow::Cow;
use std::convert::TryFrom;

use thiserror::Error;

use mzpeaks::{prelude::*, IndexType};
use mzpeaks::{
    CentroidLike, DeconvolutedCentroidLike, DeconvolutedPeakSet, MZPeakSetType, MassPeakSetType,
    PeakSet,
};
use mzpeaks::{CentroidPeak, DeconvolutedPeak, Tolerance};

#[cfg(feature = "mzsignal")]
use mzsignal::denoise::{denoise, DenoisingError};
#[cfg(feature = "mzsignal")]
use mzsignal::peak_picker::{PeakFitType, PeakPicker, PeakPickerError};
#[cfg(feature = "mzsignal")]
use mzsignal::FittedPeak;

use crate::params::ParamList;
use crate::spectrum::scan_properties::{
    Acquisition, Precursor, ScanPolarity, SignalContinuity, SpectrumDescription,
};
use crate::spectrum::bindata::{ArrayType, BinaryArrayMap, BinaryDataArrayType};
use crate::utils::mass_charge_ratio;

use super::bindata::{ArrayRetrievalError, BuildArrayMapFrom, BuildFromArrayMap};

pub trait CentroidPeakAdapting: CentroidLike + Default {}
impl<C: CentroidLike + Default> CentroidPeakAdapting for C {}
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
    pub fn base_peak(&self) -> CentroidPeak {
        match self {
            PeakDataLevel::Missing => CentroidPeak::new(0.0, 0.0, 0),
            PeakDataLevel::RawData(arrays) => {
                let intensities = arrays.intensities().unwrap();
                let result = intensities
                    .iter()
                    .enumerate()
                    .max_by(|ia, ib| ia.1.partial_cmp(ib.1).unwrap());
                if let Some((i, inten)) = result {
                    CentroidPeak::new(arrays.mzs().unwrap()[i], *inten, i as IndexType)
                } else {
                    CentroidPeak::new(0.0, 0.0, 0)
                }
            }
            PeakDataLevel::Centroid(peaks) => {
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
            PeakDataLevel::Deconvoluted(peaks) => {
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

    pub fn mz_range(&self) -> (f64, f64) {
        match self {
            PeakDataLevel::Missing => (0.0, 0.0),
            PeakDataLevel::RawData(arrays) => {
                let mzs = arrays.mzs().unwrap();
                if mzs.len() == 0 {
                    (0.0, 0.0)
                } else {
                    (*mzs.first().unwrap(), *mzs.last().unwrap())
                }
            },
            PeakDataLevel::Centroid(peaks) => {
                if peaks.len() == 0 {
                    (0.0, 0.0)
                } else {
                    (peaks[0].mz(), peaks[peaks.len() - 1].mz())
                }
            },
            PeakDataLevel::Deconvoluted(peaks) => {
                if peaks.len() == 0 {
                    (0.0, 0.0)
                } else {
                    peaks.iter().map(|p| {
                        let m = p.neutral_mass();
                        let z = p.charge();
                        mass_charge_ratio(m, z)
                    }).fold((f64::MAX, f64::MIN), |state, mz| {
                        (state.0.min(mz), state.1.max(mz))
                    })
                }
            },
        }
    }

    pub fn tic(&self) -> f32 {
        match self {
            PeakDataLevel::Missing => 0.0,
            PeakDataLevel::RawData(arrays) => {
                let intensities = arrays.intensities();
                intensities.unwrap().iter().sum()
            }
            PeakDataLevel::Centroid(peaks) => peaks.iter().map(|p| p.intensity()).sum(),
            PeakDataLevel::Deconvoluted(peaks) => peaks.iter().map(|p| p.intensity()).sum(),
        }
    }

    pub fn search(&self, query: f64, error_tolerance: Tolerance) -> Option<usize> {
        match self {
            PeakDataLevel::Missing => None,
            PeakDataLevel::RawData(arrays) => arrays.search(query, error_tolerance),
            PeakDataLevel::Centroid(peaks) => peaks.search(query, error_tolerance),
            PeakDataLevel::Deconvoluted(peaks) => peaks.search(query, error_tolerance),
        }
    }

    pub fn len(&self) -> usize {
        match self {
            PeakDataLevel::Missing => 0,
            PeakDataLevel::RawData(arrays) => arrays
                .mzs()
                .expect("Failed to locate m/z array while computing size")
                .len(),
            PeakDataLevel::Centroid(peaks) => peaks.len(),
            PeakDataLevel::Deconvoluted(peaks) => peaks.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
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

    fn precursor_mut(&mut self) -> Option<&mut Precursor> {
        let desc = self.description_mut();
        if let Some(precursor) = desc.precursor.as_mut() {
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

/// Errors that may arise when converting between different SpectrumBehavior-like types
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

#[derive(Debug, Clone, Error)]
pub enum SpectrumProcessingError {
    #[cfg(feature = "mzsignal")]
    #[error("An error occurred while denoising: {0:?}")]
    DenoisingError(DenoisingError),

    #[cfg(feature = "mzsignal")]
    #[error("An error occurred while peak picking: {0:?}")]
    PeakPickerError(PeakPickerError),

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
    /// Convert a spectrum into a [`CentroidSpectrumType`]
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

    pub fn mzs(&'lifespan self) -> Cow<'transient, [f64]> {
        self.arrays.mzs().unwrap()
    }

    pub fn intensities(&'lifespan self) -> Cow<'transient, [f32]> {
        self.arrays.intensities().unwrap()
    }

    pub fn mzs_mut(&mut self) -> Result<&mut [f64], ArrayRetrievalError> {
        self.arrays.mzs_mut()
    }

    pub fn intensities_mut(&mut self) -> Result<&mut [f32], ArrayRetrievalError> {
        self.arrays.intensities_mut()
    }

    pub fn decode_all_arrays(&mut self) -> Result<(), ArrayRetrievalError> {
        self.arrays.decode_all_arrays()
    }

    /// Convert a spectrum into a [`Spectrum`]
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

    #[cfg(feature = "mzsignal")]
    pub fn pick_peaks_with_into(
        self,
        peak_picker: &PeakPicker,
    ) -> Result<MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak>, SpectrumProcessingError> {
        let mut result = self.into_spectrum()?;
        match result.pick_peaks_with(peak_picker) {
            Ok(_) => Ok(result),
            Err(err) => Err(err),
        }
    }

    #[cfg(feature = "mzsignal")]
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

impl SpectrumLike for RawSpectrum {
    #[inline]
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }

    fn description_mut(&mut self) -> &mut SpectrumDescription {
        &mut self.description
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

impl<C: CentroidLike + Default> SpectrumLike<C> for CentroidSpectrumType<C> {
    #[inline]
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }

    fn description_mut(&mut self) -> &mut SpectrumDescription {
        &mut self.description
    }

    fn peaks(&'_ self) -> PeakDataLevel<'_, C, DeconvolutedPeak> {
        PeakDataLevel::Centroid(&self.peaks)
    }
}

impl<C: CentroidLike + Default> CentroidSpectrumType<C> {
    /// Convert a spectrum into a [`Spectrum`]
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

#[derive(Default, Debug, Clone)]
pub struct DeconvolutedSpectrumType<D: DeconvolutedCentroidLike + Default> {
    /// The spectrum metadata describing acquisition conditions and details.
    pub description: SpectrumDescription,
    /// The deisotoped and charge state deconvolved peaks
    pub deconvoluted_peaks: MassPeakSetType<D>,
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

    fn peaks(&'_ self) -> PeakDataLevel<'_, CentroidPeak, D> {
        PeakDataLevel::Deconvoluted(&self.deconvoluted_peaks)
    }
}

pub type DeconvolutedSpectrum = DeconvolutedSpectrumType<DeconvolutedPeak>;

#[derive(Default, Debug, Clone)]
/// Represent a spectrum with multiple layers of representation of the
/// peak data.
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

    fn peaks(&'_ self) -> PeakDataLevel<'_, C, D> {
        if let Some(peaks) = &self.deconvoluted_peaks {
            PeakDataLevel::Deconvoluted(peaks)
        } else if let Some(peaks) = &self.peaks {
            PeakDataLevel::Centroid(peaks)
        } else if let Some(arrays) = &self.arrays {
            PeakDataLevel::RawData(arrays)
        } else {
            PeakDataLevel::Missing
        }
    }
}

impl<C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default> MultiLayerSpectrum<C, D>
where
    C: BuildFromArrayMap + BuildArrayMapFrom,
    D: BuildFromArrayMap,
{
    pub fn try_build_centroids(&mut self) -> Result<&MZPeakSetType<C>, SpectrumConversionError> {
        if self.peaks.is_some() {
            Ok(self.peaks.as_ref().unwrap())
        } else if let Some(arrays) = self.arrays.as_ref() {
            match self.signal_continuity() {
                SignalContinuity::Centroid => {
                    let peaks = C::try_from_arrays(arrays)?.into();
                    self.peaks = Some(peaks);
                    Ok(self.peaks.as_ref().unwrap())
                },
                _ => {
                    Err(SpectrumConversionError::NotCentroided)
                }
            }
        } else {
            Err(SpectrumConversionError::NoPeakData)
        }
    }

    pub fn try_build_deconvoluted_centroids(&mut self) -> Result<&MassPeakSetType<D>, SpectrumConversionError> {
        if let Some(ref peaks) = self.deconvoluted_peaks {
            Ok(peaks)
        } else if let Some(ref arrays) = self.arrays {
            match self.signal_continuity() {
                SignalContinuity::Centroid => {
                    let peaks = D::try_from_arrays(arrays)?.into();
                    self.deconvoluted_peaks = Some(peaks);
                    Ok(self.deconvoluted_peaks.as_ref().unwrap())
                },
                _ => {
                    Err(SpectrumConversionError::NotCentroided)
                }
            }
        } else {
            Err(SpectrumConversionError::NoPeakData)
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
            let arrays = BinaryArrayMap::from(&DeconvolutedPeakSet::wrap(
                peaks.into_iter().map(|p| p.as_centroid()).collect(),
            ));

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
}

#[cfg(feature = "mzsignal")]
impl<C: CentroidLike + Default + From<FittedPeak>, D: DeconvolutedCentroidLike + Default>
    MultiLayerSpectrum<C, D>
{
    pub fn pick_peaks_with(
        &mut self,
        peak_picker: &PeakPicker,
    ) -> Result<(), SpectrumProcessingError> {
        if let Some(arrays) = &self.arrays {
            let mz_array = arrays.mzs()?;
            let intensity_array = arrays.intensities()?;
            let mut acc = Vec::new();
            match peak_picker.discover_peaks(&mz_array, &intensity_array, &mut acc) {
                Ok(_) => {
                    let peaks: MZPeakSetType<C> = acc.into_iter().map(|p| C::from(p)).collect();
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
    use crate::prelude::*;

    #[test_log::test]
    fn test_profile_read() {
        let mut reader = MzMLReader::open_path("./test/data/three_test_scans.mzML")
            .expect("Failed to open test file");
        reader.reset();
        let mut scan = reader.next().expect("Failed to read spectrum");
        assert_eq!(scan.signal_continuity(), SignalContinuity::Profile);
        assert_eq!(scan.ms_level(), 1);
        assert_eq!(scan.polarity(), ScanPolarity::Positive);
        assert!(scan.precursor().is_none());

        if let Err(err) = scan.pick_peaks(1.0, PeakFitType::Quadratic) {
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
