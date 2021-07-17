use std::borrow;

use crate::peaks::CentroidPeak;
use crate::peaks::{PeakCollection, PeakSet};
use crate::spectrum::scan_properties::{
    Acquisition, Precursor, ScanSiganlContinuity, SpectrumDescription,
};
use crate::spectrum::signal::{ArrayType, BinaryArrayMap};

pub trait SpectrumBehavior {
    fn description(&self) -> &SpectrumDescription;

    fn acquisition(&self) -> &Acquisition {
        &self.description().acquisition
    }

    fn precursor(&self) -> &Option<Precursor> {
        let desc = self.description();
        &desc.precursor
    }

    fn start_time(&self) -> f64 {
        let acq = self.acquisition();
        match acq.scans.first() {
            Some(evt) => evt.start_time,
            None => 0.0,
        }
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
}

impl RawSpectrum {
    pub fn into_centroid(self) -> Result<CentroidSpectrum, SpectrumConversionError> {
        if !matches!(self.description.is_profile, ScanSiganlContinuity::Centroid) {
            return Err(SpectrumConversionError::NotCentroided);
        }

        let mut centroid = CentroidSpectrum {
            description: self.description,
            peaks: PeakSet::empty(),
        };

        let mz_array = self
            .arrays
            .get(&ArrayType::MZArray)
            .expect("Did not find m/z array")
            .to_f64()
            .expect("Failed to decode m/z array");
        let intensity_array = self
            .arrays
            .get(&ArrayType::IntensityArray)
            .expect("Did not find intensity array")
            .to_f32()
            .expect("Failed to decode intensity array");

        if mz_array.len() != intensity_array.len() {
            return Err(SpectrumConversionError::NotCentroided);
        }

        for (mz, intensity) in mz_array.iter().zip(intensity_array.iter()) {
            centroid.peaks.push(CentroidPeak {
                mz: *mz,
                intensity: *intensity,
                ..CentroidPeak::default()
            });
        }
        centroid.description.is_profile = ScanSiganlContinuity::Centroid;

        Ok(centroid)
    }

    pub fn mzs(&self) -> borrow::Cow<[f64]> {
        let mz_array = self
            .arrays
            .get(&ArrayType::MZArray)
            .expect("Did not find m/z array")
            .to_f64()
            .expect("Failed to decode m/z array");
        mz_array
    }

    pub fn intensities(&self) -> borrow::Cow<[f32]> {
        let intensity_array = self
            .arrays
            .get(&ArrayType::IntensityArray)
            .expect("Did not find intensity array")
            .to_f32()
            .expect("Failed to decode intensity array");
        intensity_array
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
