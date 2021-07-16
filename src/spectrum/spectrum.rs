use std::borrow;

use crate::spectrum::scan_properties::{SpectrumDescription, Acquisition, Precursor, ScanSiganlContinuity};
use crate::spectrum::signal::{BinaryArrayMap, ArrayType};
use crate::peaks::CentroidPeak;
use crate::peaks::{PeakSet, PeakCollection};


pub trait SpectrumBehavior {
    fn get_description(&self) -> &SpectrumDescription;

    fn get_acquisition(&self) -> &Acquisition {
        &self.get_description().acquisition
    }

    fn get_precursor(&self) -> &Option<Precursor> {
        let desc = self.get_description();
        &desc.precursor
    }

    fn get_start_time(&self) -> f64 {
        let acq = self.get_acquisition();
        match acq.scans.first() {
            Some(evt) => evt.start_time,
            None => 0.0
        }
    }

    fn get_id(&self) -> &str {
        &self.get_description().id
    }

    fn get_index(&self) -> usize {
        self.get_description().index
    }
}


#[derive(Default, Debug, Clone)]
pub struct RawSpectrum {
    pub description: SpectrumDescription,
    pub arrays: BinaryArrayMap
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

        let mz_array = self.arrays.get(&ArrayType::MZArray
            ).expect("Did not find m/z array").to_f64().expect(
                "Failed to decode m/z array");
        let intensity_array = self.arrays.get(
            &ArrayType::IntensityArray).expect(
                "Did not find intensity array").to_f32().expect(
                    "Failed to decode intensity array");

        if mz_array.len() != intensity_array.len() {
            return Err(SpectrumConversionError::NotCentroided);
        }

        for (mz, intensity) in mz_array.iter().zip(intensity_array.iter()) {
            centroid.peaks.push(
                CentroidPeak {
                    mz: *mz,
                    intensity: *intensity,
                    ..CentroidPeak::default() });
        }
        centroid.description.is_profile = ScanSiganlContinuity::Centroid;

        Ok(centroid)
    }

    pub fn get_mzs(&self) -> borrow::Cow<[f64]> {
        let mz_array = self.arrays.get(&ArrayType::MZArray
                    ).expect("Did not find m/z array").to_f64().expect(
                        "Failed to decode m/z array");
        mz_array
    }

    pub fn get_intensities(&self) -> borrow::Cow<[f32]> {
        let intensity_array = self.arrays.get(
            &ArrayType::IntensityArray).expect(
                "Did not find intensity array").to_f32().expect(
                    "Failed to decode intensity array");
        intensity_array
    }
}

impl SpectrumBehavior for RawSpectrum {
    fn get_description(&self) -> &SpectrumDescription {
        &self.description
    }
}


#[derive(Default, Debug, Clone)]
pub struct CentroidSpectrum {
    pub description: SpectrumDescription,
    pub peaks: PeakSet
}


impl SpectrumBehavior for CentroidSpectrum {
    fn get_description(&self) -> &SpectrumDescription {
        &self.description
    }
}
