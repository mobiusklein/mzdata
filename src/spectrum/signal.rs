use std::collections::HashMap;

use super::params::{ParamList};

type Bytes = Vec<u8>;


#[derive(Debug, Clone, Copy, PartialEq, Hash)]
pub enum BinaryDataArrayType {
    Unknown,
    Float64,
    Float32,
    Int64,
    Int32,
    ASCII,
}

impl Default for BinaryDataArrayType {
    fn default() -> BinaryDataArrayType {
        BinaryDataArrayType::Unknown
    }
}


#[derive(Debug, Clone, Copy, PartialEq, Hash)]
pub enum BinaryCompressionType {
    NoCompression,
    Zlib,
    NumpressLinear,
    NumpressSLOF,
    NumpressPIC,
}

impl Default for BinaryCompressionType {
    fn default() -> BinaryCompressionType {
        BinaryCompressionType::NoCompression
    }
}


#[derive(Debug, Clone, PartialEq, Hash, Eq)]
pub enum ArrayType {
    Unknown,
    MZArray,
    IntensityArray,
    ChargeArray,
    TimeArray,
    WavelengthArray,
    IonMobilityArray,
    MeanIonMobilityArray,
    RawIonMobilityArray,
    DeconvolutedIonMobilityArray,
    NonStandardDataArray {name: String},
}

impl Default for ArrayType {
    fn default() -> ArrayType {
        ArrayType::Unknown
    }
}


#[derive(Debug, Default, Clone)]
pub struct DataArray {
    pub data: Bytes,
    pub dtype: BinaryDataArrayType,
    pub compression: BinaryCompressionType,
    pub name: ArrayType,
    pub params: ParamList
}


impl DataArray {
    pub fn new() -> DataArray {
        DataArray {
            ..Default::default()
        }
    }

    pub fn clear(&mut self) {
        self.data.clear();
        self.params.clear();
    }
}


#[derive(Debug, Default, Clone)]
pub struct BinaryArrayMap {
    pub byte_buffer_map: HashMap<ArrayType, DataArray>
}

impl BinaryArrayMap {
    pub fn new() -> BinaryArrayMap {
        BinaryArrayMap {
            ..Default::default()
        }
    }

    pub fn add(&mut self, array: DataArray) {
        self.byte_buffer_map.insert(array.name.clone(), array);
    }

    pub fn clear(&mut self) {
        self.byte_buffer_map.clear();
    }
}
