use std::borrow;

use crate::peaks::CentroidPeak;
use crate::peaks::{PeakCollection, PeakSet};
use crate::spectrum::scan_properties::{
    Acquisition, Precursor, SignalContinuity, SpectrumDescription,
};
use crate::spectrum::signal::BinaryArrayMap;

pub trait SpectrumBehavior {
    fn description(&self) -> &SpectrumDescription;

    fn acquisition(&self) -> &Acquisition {
        &self.description().acquisition
    }

    fn precursor(&self) -> Option<&Precursor> {
        let desc = self.description();
        if let Some(precursor) = &desc.precursor {
            Some(precursor)
        } else {
            None
        }
    }

    fn start_time(&self) -> f64 {
        let acq = self.acquisition();
        match acq.scans.first() {
            Some(evt) => evt.start_time,
            None => 0.0,
        }
    }

    fn ms_level(&self) -> u8 {
        self.description().ms_level
    }

    fn id(&self) -> &str {
        &self.description().id
    }

    fn index(&self) -> usize {
        self.description().index
    }
}

#[derive(Default, Debug, Clone)]
pub struct RawSpectrum {
    pub description: SpectrumDescription,
    pub arrays: BinaryArrayMap,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpectrumConversionError {
    MZIntensityArraySizeMismatch,
    NotCentroided,
    NoPeakData,
}

impl<'transient, 'lifespan: 'transient> RawSpectrum {
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

    pub fn into_spectrum(self) -> Result<Spectrum, SpectrumConversionError> {
        Ok(Spectrum {
            arrays: Some(self.arrays),
            description: self.description,
            ..Default::default()
        })
    }
}

impl SpectrumBehavior for RawSpectrum {
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }
}

#[derive(Default, Debug, Clone)]
pub struct CentroidSpectrum {
    pub description: SpectrumDescription,
    pub peaks: PeakSet,
}

impl SpectrumBehavior for CentroidSpectrum {
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }
}

impl CentroidSpectrum {
    pub fn into_spectrum(self) -> Result<Spectrum, SpectrumConversionError> {
        Ok(Spectrum {
            peaks: Some(self.peaks),
            description: self.description,
            ..Default::default()
        })
    }
}

#[derive(Default, Debug, Clone)]
pub struct Spectrum {
    pub description: SpectrumDescription,
    pub arrays: Option<BinaryArrayMap>,
    pub peaks: Option<PeakSet>,
}

impl SpectrumBehavior for Spectrum {
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }
}

impl Spectrum {
    pub fn into_centroid(self) -> Result<CentroidSpectrum, SpectrumConversionError> {
        if let Some(peaks) = self.peaks {
            let mut result = CentroidSpectrum {
                peaks,
                description: self.description,
            };
            result.description.signal_continuity = SignalContinuity::Centroid;
            return Ok(result);
        } else if let Some(arrays) = &self.arrays {
            let mz_array = arrays.mzs();
            let intensity_array = arrays.intensities();
            let mut peaks = PeakSet::empty();

            if mz_array.len() != intensity_array.len() {
                return Err(SpectrumConversionError::NotCentroided);
            }

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
            return Ok(centroid);
        }
        Err(SpectrumConversionError::NotCentroided)
    }

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
