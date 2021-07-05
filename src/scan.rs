use std::collections::HashMap;

use crate::peak_set::{PeakSet};

#[derive(Debug, Clone, Default)]
pub struct Precursor {
    pub mz: f64,
    pub intensity: f32,
    pub charge: i32,
    pub precursor_id: String,
}


#[derive(Debug, Clone, Default)]
pub struct ScanDescription {
    pub id: String,
    pub index: usize,
    pub ms_level: i8,
    pub time: f32,
    pub precursor_information: Option<Precursor>,
    pub annotations: HashMap<String, String>,
}


pub struct RawDataArrays {
    pub mz: Vec<f64>,
    pub intensity: Vec<f64>
}


pub struct RawScan {
    pub description: ScanDescription,
    pub arrays: RawDataArrays
}


pub struct CentroidScan {
    pub description: ScanDescription,
    pub peaks: PeakSet
}

