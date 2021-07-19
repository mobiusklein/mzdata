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

use crate::peaks::CentroidPeak;
use crate::peaks::{PeakCollection, PeakSet};
use crate::spectrum::scan_properties::{
    Acquisition, Precursor, SignalContinuity, SpectrumDescription,
};
use crate::spectrum::signal::BinaryArrayMap;


#[derive(Debug)]
/// An variant for dispatching to different strategies of computing
/// common statistics of different levels of peak data.
pub enum PeakDataLevel<'lifespan> {
    Missing,
    RawData(&'lifespan BinaryArrayMap),
    Centroid(&'lifespan PeakSet),
}

impl<'lifespan> PeakDataLevel<'lifespan> {
    /// Compute the base peak of a spectrum
    pub fn base_peak(&self) -> (usize, f64, f32) {
        match self {
            PeakDataLevel::Missing => {
                (0, 0.0, 0.0)
            },
            PeakDataLevel::RawData(arrays) => {
                let intensities = arrays.intensities();
                let result = intensities.iter().enumerate().max_by(|ia, ib| ia.1.partial_cmp(&ib.1).unwrap());
                if let Some((i, inten)) = result {
                    (i, arrays.mzs()[i], *inten)
                } else {
                    (0, 0.0, 0.0)
                }
            },
            PeakDataLevel::Centroid(peaks) => {
                let result = peaks.iter().enumerate().max_by(|ia, ib| ia.1.intensity.partial_cmp(&ib.1.intensity).unwrap());
                if let Some((i, peak)) = result {
                    (i, peak.mz, peak.intensity)
                } else {
                    (0, 0.0, 0.0)
                }
            }
        }
    }
}

/// A trait for providing a uniform delegated access to spectrum metadata
pub trait SpectrumBehavior {
    /// The method to access the spectrum description itself, which supplies
    /// the data for most other methods on this trait.
    fn description(&self) -> &SpectrumDescription;

    /// Access the acquisition information for this spectrum.
    fn acquisition(&self) -> &Acquisition {
        &self.description().acquisition
    }

    /// Access the precursor information, if it exists.
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
    fn start_time(&self) -> f64 {
        let acq = self.acquisition();
        match acq.scans.first() {
            Some(evt) => evt.start_time,
            None => 0.0,
        }
    }

    /// Access the MS exponentiation level
    fn ms_level(&self) -> u8 {
        self.description().ms_level
    }

    /// Access the native ID string for the spectrum
    fn id(&self) -> &str {
        &self.description().id
    }

    /// Access the index of the spectrum in the source file
    fn index(&self) -> usize {
        self.description().index
    }

    /// Retrieve the most processed representation of the mass spectrum's
    /// signal
    fn peaks(&'_ self) -> PeakDataLevel<'_>;
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
    NotCentroided,
    NoPeakData,
}

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
            return Err(SpectrumConversionError::NotCentroided);
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
    pub fn into_spectrum(self) -> Result<Spectrum, SpectrumConversionError> {
        Ok(Spectrum {
            arrays: Some(self.arrays),
            description: self.description,
            ..Default::default()
        })
    }
}

impl<'lifespan> SpectrumBehavior for RawSpectrum {
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }

    fn peaks(&'_ self) -> PeakDataLevel<'_> {
        PeakDataLevel::RawData(&self.arrays)
    }
}

#[derive(Default, Debug, Clone)]
/// Represents a spectrum that has been centroided
pub struct CentroidSpectrum {
    /// The spectrum metadata describing acquisition conditions and details.
    pub description: SpectrumDescription,
    /// The picked centroid peaks
    pub peaks: PeakSet,
}

impl<'lifespan> SpectrumBehavior for CentroidSpectrum {
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }

    fn peaks(&'_ self) -> PeakDataLevel<'_> {
        PeakDataLevel::Centroid(&self.peaks)
    }
}

impl CentroidSpectrum {
    /// Convert a spectrum into a [`Spectrum`]
    pub fn into_spectrum(self) -> Result<Spectrum, SpectrumConversionError> {
        Ok(Spectrum {
            peaks: Some(self.peaks),
            description: self.description,
            ..Default::default()
        })
    }
}

#[derive(Default, Debug, Clone)]
/// Represent a spectrum with multiple layers of representation of the
/// peak data.
///
/// **Not directly used at this time**
pub struct Spectrum {
    /// The spectrum metadata describing acquisition conditions and details.
    pub description: SpectrumDescription,
    /// The (potentially absent) data arrays describing the m/z, intensity,
    /// and potentially other measured properties
    pub arrays: Option<BinaryArrayMap>,
    // The (potentially absent) centroid peaks
    pub peaks: Option<PeakSet>,
}

impl<'lifespan> SpectrumBehavior for Spectrum {
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }

    fn peaks(&'_ self) -> PeakDataLevel<'_> {
        if let Some(peaks) = &self.peaks {
            PeakDataLevel::Centroid(peaks)
        } else if let Some(arrays) = &self.arrays {
            PeakDataLevel::RawData(arrays)
        } else {
            PeakDataLevel::Missing
        }
    }

}

impl Spectrum {
    /// Convert a spectrum into a [`CentroidSpectrum`]
    pub fn into_centroid(self) -> Result<CentroidSpectrum, SpectrumConversionError> {
        if let Some(peaks) = self.peaks {
            let mut result = CentroidSpectrum {
                peaks,
                description: self.description,
            };
            result.description.signal_continuity = SignalContinuity::Centroid;
            return Ok(result);
        } else if let Some(arrays) = &self.arrays {
            let peaks = PeakSet::from(arrays);
            let mut centroid = CentroidSpectrum {
                description: self.description,
                peaks,
            };
            centroid.description.signal_continuity = SignalContinuity::Centroid;
            return Ok(centroid);
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
        } else if let Some(peaks) = self.peaks {
            let arrays = BinaryArrayMap::from(peaks);

            let mut result = RawSpectrum {
                arrays,
                description: self.description,
            };
            result.description.signal_continuity = SignalContinuity::Centroid;
            Ok(result)
        } else {
            Err(SpectrumConversionError::NoPeakData)
        }
    }
}
