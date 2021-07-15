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
    pub polarity: i8,

    pub precursor_information: Option<Precursor>,
    pub annotations: HashMap<String, String>,
}


type I32ArrayMap = HashMap<String, Vec<i32>>;
type F64ArrayMap = HashMap<String, Vec<f64>>;


#[derive(Debug, Default)]
pub struct BinaryArrayMap {
    pub double_arrays: F64ArrayMap,
    pub int_arrays: I32ArrayMap
}


pub struct RawDataArrays {
    pub arrays: BinaryArrayMap
}


pub struct RawScan {
    pub description: ScanDescription,
    pub arrays: RawDataArrays
}


pub struct CentroidScan {
    pub description: ScanDescription,
    pub peaks: PeakSet
}

