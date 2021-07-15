use std::fmt;
use std::fmt::{Formatter, Debug, Display};
use std::collections::HashMap;
use std::io::prelude::*;

use flate2::write::ZlibDecoder;
use flate2::Compression;
use base64;

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


#[derive(Default, Clone)]
pub struct DataArray {
    pub data: Bytes,
    pub dtype: BinaryDataArrayType,
    pub compression: BinaryCompressionType,
    pub name: ArrayType,
    pub params: ParamList
}

impl Debug for DataArray {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_struct("DataArray")
            .field("name", &self.name)
            .field("data size", &self.data.len())
            .field("dtype", &self.dtype)
            .field("compression", &self.compression)
            .field("params", &self.params)
            .finish()
    }
}

impl DataArray {
    pub fn new() -> DataArray {
        DataArray {
            ..Default::default()
        }
    }

    fn _decompress(&self, bytestring: &Bytes) -> Bytes {
        let result = Bytes::new();
        let mut decompressor = ZlibDecoder::new(result);
        decompressor.write_all(&bytestring[..]);
        let result = decompressor.finish().expect("Error decompression");
        result
    }

    pub fn decode(&self) -> Bytes {
        let bytestring = base64::decode(&self.data).expect("Failed to decode base64 array");
        println!("Compression: {:?}, {} Bytes Input", self.compression, bytestring.len());
        match self.compression {
            BinaryCompressionType::NoCompression => {
                return bytestring;
            },
            BinaryCompressionType::Zlib => {
                return self._decompress(&bytestring);
            },
            _ => {
                panic!("Could not decompress array {:?}", self.compression);
            }
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
