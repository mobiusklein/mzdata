
use crate::spectrum::scan_properties::{SpectrumDescription, Acquisition, Precursor};
use crate::spectrum::signal::{BinaryArrayMap};
use crate::peak_set::PeakSet;


pub trait SpectrumBehavior {
    fn get_description(&self) -> &SpectrumDescription;

    fn get_acquisition(&self) -> &Acquisition {
        &self.get_description().acquisition
    }

    fn get_precursor(&self) -> &Option<Precursor> {
        let desc = self.get_description();
        &desc.precursor
    }

    fn get_scan_time(&self) -> f64 {
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
