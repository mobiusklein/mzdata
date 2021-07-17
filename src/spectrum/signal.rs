use std::borrow::Cow;
use std::collections::HashMap;
use std::fmt;
use std::fmt::{Debug, Formatter};
use std::io::prelude::*;
use std::mem;
use std::slice;

use num_traits::cast::AsPrimitive;
use num_traits::Num;

use base64;
use flate2::write::ZlibDecoder;

use super::params::ParamList;
use crate::peaks::{PeakCollection, PeakSet};

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

impl BinaryDataArrayType {
    pub const fn size_of(&self) -> usize {
        match self {
            BinaryDataArrayType::Unknown | BinaryDataArrayType::ASCII => 1,
            BinaryDataArrayType::Float32 | BinaryDataArrayType::Int32 => 4,
            BinaryDataArrayType::Float64 | BinaryDataArrayType::Int64 => 8,
        }
    }
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
    Decoded,
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
    NonStandardDataArray { name: String },
}

impl Default for ArrayType {
    fn default() -> ArrayType {
        ArrayType::Unknown
    }
}

impl ArrayType {
    pub const fn preferred_dtype(&self) -> BinaryDataArrayType {
        match self {
            ArrayType::MZArray => BinaryDataArrayType::Float64,
            ArrayType::IntensityArray => BinaryDataArrayType::Float32,
            ArrayType::ChargeArray => BinaryDataArrayType::Int32,
            _ => BinaryDataArrayType::Float32,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum ArrayRetrievalError {
    NotFound,
    DecompressionError,
    DecodeError,
    DataTypeSizeMismatch,
}

#[derive(Default, Clone)]
pub struct DataArray {
    pub data: Bytes,
    pub dtype: BinaryDataArrayType,
    pub compression: BinaryCompressionType,
    pub name: ArrayType,
    pub params: ParamList,
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

/// A type to represent a base64-encoded, possibly compressed data
/// array of a fixed size, usually numeric, type. It can be decoded,
/// and it can
impl<'transient, 'lifespan: 'transient> DataArray {
    pub fn new() -> DataArray {
        DataArray {
            ..Default::default()
        }
    }

    pub fn from_name(name: &ArrayType) -> DataArray {
        DataArray {
            dtype: name.preferred_dtype(),
            name: name.clone(),
            compression: BinaryCompressionType::Decoded,
            ..Default::default()
        }
    }

    pub fn from_name_and_type(name: &ArrayType, dtype: BinaryDataArrayType) -> DataArray {
        DataArray {
            dtype,
            name: name.clone(),
            compression: BinaryCompressionType::Decoded,
            ..Default::default()
        }
    }

    pub fn from_name_type_size(
        name: &ArrayType,
        dtype: BinaryDataArrayType,
        size: usize,
    ) -> DataArray {
        DataArray {
            dtype,
            name: name.clone(),
            data: Bytes::with_capacity(size),
            compression: BinaryCompressionType::Decoded,
            ..Default::default()
        }
    }

    fn decompres_zlib(&self, bytestring: &[u8]) -> Bytes {
        let result = Bytes::new();
        let mut decompressor = ZlibDecoder::new(result);
        decompressor
            .write_all(bytestring)
            .expect("Error decompression");
        decompressor.finish().expect("Error decompression")
    }

    pub fn decode_and_store(&mut self) -> Result<BinaryCompressionType, ArrayRetrievalError> {
        match self.decode() {
            Ok(data) => {
                match data {
                    // The only time this is a borrow is when the data are already
                    // decoded.
                    Cow::Borrowed(_view) => Ok(self.compression),
                    Cow::Owned(buffer) => {
                        self.data = buffer;
                        self.compression = BinaryCompressionType::Decoded;
                        Ok(self.compression)
                    }
                }
            }
            Err(err) => Err(err),
        }
    }

    pub fn decode(&'lifespan self) -> Result<Cow<'lifespan, Bytes>, ArrayRetrievalError> {
        match self.compression {
            BinaryCompressionType::Decoded => Ok(Cow::Borrowed(&self.data)),
            BinaryCompressionType::NoCompression => {
                let bytestring = base64::decode(&self.data).expect("Failed to decode base64 array");
                Ok(Cow::Owned(bytestring))
            }
            BinaryCompressionType::Zlib => {
                let bytestring = base64::decode(&self.data).expect("Failed to decode base64 array");
                Ok(Cow::Owned(self.decompres_zlib(&bytestring)))
            }
            _ => Err(ArrayRetrievalError::DecompressionError),
        }
    }

    pub fn coerce_from<T: Clone + Sized>(
        &'lifespan self,
        buffer: Cow<'transient, Bytes>,
    ) -> Result<Cow<'transient, [T]>, ArrayRetrievalError> {
        let n = buffer.len();
        let z = mem::size_of::<T>();
        if n % z != 0 {
            return Err(ArrayRetrievalError::DataTypeSizeMismatch);
        }
        let m = n / z;
        unsafe {
            Ok(Cow::Borrowed(slice::from_raw_parts(
                buffer.as_ptr() as *const T,
                m,
            )))
        }
    }

    pub fn coerce<T: Clone + Sized>(
        &'lifespan self,
    ) -> Result<Cow<'transient, [T]>, ArrayRetrievalError> {
        match self.decode() {
            Ok(data) => self.coerce_from(data),
            Err(err) => Err(err),
        }
    }

    pub fn clear(&mut self) {
        self.data.clear();
        self.params.clear();
    }

    fn transmute<S: Num + Clone + AsPrimitive<D>, D: Num + Clone + Copy + 'static>(
        &'lifespan self,
    ) -> Result<Cow<'transient, [D]>, ArrayRetrievalError> {
        match self.coerce::<S>() {
            Ok(view) => {
                let owned = view.into_owned();
                let res = owned.iter().map(|a| a.as_()).collect();
                Ok(Cow::Owned(res))
            }
            Err(err) => Err(err),
        }
    }

    pub fn to_f32(&'lifespan self) -> Result<Cow<'transient, [f32]>, ArrayRetrievalError> {
        type D = f32;
        match self.dtype {
            BinaryDataArrayType::Float32 => self.coerce::<D>(),
            BinaryDataArrayType::Float64 => {
                type S = f64;
                self.transmute::<S, D>()
            }
            BinaryDataArrayType::Int32 => {
                type S = i32;
                self.transmute::<S, D>()
            }
            BinaryDataArrayType::Int64 => {
                type S = i64;
                self.transmute::<S, D>()
            }
            _ => Err(ArrayRetrievalError::DataTypeSizeMismatch),
        }
    }

    pub fn to_f64(&'lifespan self) -> Result<Cow<'transient, [f64]>, ArrayRetrievalError> {
        type D = f64;
        match self.dtype {
            BinaryDataArrayType::Float32 => {
                type S = f32;
                self.transmute::<S, D>()
            }
            BinaryDataArrayType::Float64 => self.coerce(),
            BinaryDataArrayType::Int32 => {
                type S = i32;
                self.transmute::<S, D>()
            }
            BinaryDataArrayType::Int64 => {
                type S = i64;
                self.transmute::<S, D>()
            }
            _ => Err(ArrayRetrievalError::DataTypeSizeMismatch),
        }
    }

    pub fn to_i32(&'lifespan self) -> Result<Cow<'transient, [i32]>, ArrayRetrievalError> {
        type D = i32;
        match self.dtype {
            BinaryDataArrayType::Float32 => {
                type S = f32;
                self.transmute::<S, D>()
            }
            BinaryDataArrayType::Float64 => {
                type S = f64;
                self.transmute::<S, D>()
            }
            BinaryDataArrayType::Int32 => self.coerce::<D>(),
            BinaryDataArrayType::Int64 => {
                type S = i64;
                self.transmute::<S, D>()
            }
            _ => Err(ArrayRetrievalError::DataTypeSizeMismatch),
        }
    }
}

#[derive(Debug, Default, Clone)]
pub struct BinaryArrayMap {
    pub byte_buffer_map: HashMap<ArrayType, DataArray>,
}

impl<'transient, 'lifespan: 'transient> BinaryArrayMap {
    pub fn new() -> BinaryArrayMap {
        BinaryArrayMap {
            ..Default::default()
        }
    }

    pub fn add(&mut self, array: DataArray) {
        self.byte_buffer_map.insert(array.name.clone(), array);
    }

    pub fn get(&self, array_type: &ArrayType) -> Option<&DataArray> {
        self.byte_buffer_map.get(array_type)
    }

    pub fn clear(&mut self) {
        self.byte_buffer_map.clear();
    }

    pub fn mzs(&'lifespan self) -> Cow<'transient, [f64]> {
        let mz_array = self
            .get(&ArrayType::MZArray)
            .expect("Did not find m/z array")
            .to_f64()
            .expect("Failed to decode m/z array");
        mz_array
    }

    pub fn intensities(&'lifespan self) -> Cow<'transient, [f32]> {
        let intensities = self
            .get(&ArrayType::IntensityArray)
            .expect("Did not find intensity array")
            .to_f32()
            .expect("Failed to decode intensity array");
        intensities
    }

    pub fn from(peaks: PeakSet) -> BinaryArrayMap {
        let mut arrays = BinaryArrayMap::new();

        let mut mz_array = DataArray::from_name_type_size(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            peaks.len() * BinaryDataArrayType::Float64.size_of(),
        );

        let mut intensity_array = DataArray::from_name_type_size(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            peaks.len() * BinaryDataArrayType::Float32.size_of(),
        );

        mz_array.compression = BinaryCompressionType::Decoded;
        intensity_array.compression = BinaryCompressionType::Decoded;

        for p in peaks.iter() {
            let mz: f64 = p.mz;
            let inten: f32 = p.intensity;

            let raw_bytes: [u8; mem::size_of::<f64>()] = unsafe { mem::transmute(mz) };
            mz_array.data.extend(raw_bytes);

            let raw_bytes: [u8; mem::size_of::<f32>()] = unsafe { mem::transmute(inten) };
            intensity_array.data.extend(raw_bytes);
        }

        arrays.add(mz_array);
        arrays.add(intensity_array);
        arrays
    }
}
