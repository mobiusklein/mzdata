use std::borrow::Cow;
use std::collections::HashMap;
use std::fmt;
use std::fmt::{Debug, Formatter};
use std::io::prelude::*;
use std::mem;
use std::ops::Mul;
use bytemuck::{self, Pod};
use std::slice;

use num_traits::cast::AsPrimitive;
use num_traits::{Num, Float};

use log::warn;

use base64;
use base64::{engine::general_purpose::STANDARD as Base64Std, Engine as _};
use flate2::write::{ZlibDecoder, ZlibEncoder};
use flate2::Compression;
use mzpeaks::prelude::*;
use mzpeaks::{CentroidPeak, DeconvolutedPeak, DeconvolutedPeakSet, MZPeakSetType, PeakSet};

#[cfg(feature = "numpress")]
use numpress;

use crate::params::{ParamList, Unit};
use crate::utils::neutral_mass;

type Bytes = Vec<u8>;

pub fn to_bytes<T: Pod>(data: &[T]) -> Bytes {
    bytemuck::cast_slice(data).to_vec()
}

pub fn as_bytes<T: Pod>(data: &[T]) -> &[u8] {
    bytemuck::cast_slice(data)
}


pub fn vec_as_bytes<T: Pod>(data: Vec<T>) -> Bytes {
    bytemuck::cast_vec(data)
}


pub fn linear_prediction_decoding<F: Float + Mul<F>>(values: &mut [F]) -> &mut [F] {
    let two = F::from(2.0).unwrap();

    for i in 0..values.len() {
        if i < 2 {
            continue;
        }
        let v = values[i] + two * values[i - 1] - values[i - 2] - values[1];
        values[i] = v;
    }
    values
}


pub fn linear_prediction_encoding<F: Float + Mul<F>>(values: &mut [F]) -> &mut [F] {
    let n = values.len();
    if n < 3 {
        return values
    }
    let offset = values[1];
    let mut prev2 = values[0];
    let mut prev1 = values[1];
    let two = F::from(2.0).unwrap();
    for i in 0..values.len() {
        values[i] = offset + values[i] - two * prev1 + prev2;
        let tmp = prev1;
        prev1 = values[i] + two * prev1 - prev2 - offset;
        prev2 = tmp;
    }
    values
}


pub fn delta_decoding<F: Float + Mul<F>>(values: &mut [F]) -> &mut [F] {
    for i in 0..values.len() {
        if i < 2 {
            continue;
        }
        values[i] = values[i] + values[i - 1] - values[0];
    }
    values
}


pub fn delta_encoding<F: Float + Mul<F>>(values: &mut [F]) -> &mut [F] {
    let n = values.len();
    if n < 2 {
        return values
    }
    let mut prev = values[0];
    let offset = values[0];
    for i in 1..n {
        let tmp = values[i];
        values[i] = offset + values[i] - prev;
        prev = tmp;
    }
    values
}


/// The canonical primitive data types found in MS data file formats
/// supported by the PSI-MS controlled vocabulary
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


/// The range of compression and encoding states that a raw byte buffer
/// might be in during different stages of decoding. Other than `[BinaryCompressionType::Decoded]`,
/// these states may or may not include intermediate base64 encoding.
#[derive(Debug, Clone, Copy, PartialEq, Hash)]
pub enum BinaryCompressionType {
    NoCompression,
    Zlib,
    NumpressLinear,
    NumpressSLOF,
    NumpressPIC,
    NumpressLinearZlib,
    NumpressSLOFZlib,
    NumpressPICZlib,
    LinearPrediction,
    DeltaPrediction,
    Decoded,
}

impl BinaryCompressionType {
    /// Generate a user-understandable message about why a compression conversion operation failed
    pub fn unsupported_msg(&self, context: Option<&str>) -> String {
        match context {
            Some(ctx) => format!("Cannot decode array compressed with {:?} ({})", self, ctx),
            None => format!("Cannot decode array compressed with {:?}", self)
        }
    }
}

impl Default for BinaryCompressionType {
    fn default() -> BinaryCompressionType {
        BinaryCompressionType::NoCompression
    }
}


/// The kinds of data arrays found in mass spectrometry data files governed
/// by the PSI-MS controlled vocabulary.
#[derive(Debug, Clone, PartialEq, Hash, Eq, PartialOrd, Ord)]
pub enum ArrayType {
    Unknown,
    MZArray,
    IntensityArray,
    ChargeArray,
    SignalToNoiseArray,
    TimeArray,
    WavelengthArray,
    IonMobilityArray,
    MeanIonMobilityArray,
    RawIonMobilityArray,
    DeconvolutedIonMobilityArray,
    NonStandardDataArray { name: Box<String> },
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


/// A high level set of failure modes that an operation to retrieve a typed memory buffer
/// from a `[BinaryArrayMap]` might encounter. May also be used to represented conversion
/// during reading or writing.
#[derive(Debug, Clone)]
pub enum ArrayRetrievalError {
    NotFound,
    DecompressionError(String),
    DecodeError,
    DataTypeSizeMismatch,
}

impl std::fmt::Display for ArrayRetrievalError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::error::Error for ArrayRetrievalError {}

#[cfg(feature = "numpress")]
impl From<numpress::Error> for ArrayRetrievalError {
    fn from(value: numpress::Error) -> Self {
        ArrayRetrievalError::DecompressionError(value.to_string())
    }
}


/// Represents a data array
#[derive(Default, Clone)]
pub struct DataArray {
    pub data: Bytes,
    pub dtype: BinaryDataArrayType,
    pub compression: BinaryCompressionType,
    pub name: ArrayType,
    pub params: Option<Box<ParamList>>,
    pub unit: Unit,
}

impl Debug for DataArray {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_struct("DataArray")
            .field("name", &self.name)
            .field("data size", &self.data.len())
            .field("dtype", &self.dtype)
            .field("compression", &self.compression)
            .field("params", &self.params)
            .field("unit", &self.unit)
            .finish()
    }
}

pub trait ByteArrayView<'transient, 'lifespan: 'transient> {
    fn view(&'lifespan self) -> Result<Cow<'lifespan, [u8]>, ArrayRetrievalError>;

    fn coerce_from<T: Pod>(
        buffer: Cow<'transient, [u8]>,
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

    fn coerce<T: Pod>(
        &'lifespan self,
    ) -> Result<Cow<'transient, [T]>, ArrayRetrievalError> {
        match self.view() {
            Ok(data) => Self::coerce_from(data),
            Err(err) => Err(err),
        }
    }

    fn convert<S: Num + Clone + AsPrimitive<D> + Pod, D: Num + Clone + Copy + 'static>(
        &'lifespan self,
    ) -> Result<Cow<'transient, [D]>, ArrayRetrievalError> {
        match self.coerce::<S>() {
            Ok(view) => {
                match view {
                    Cow::Borrowed(view) => {
                        Ok(Cow::Owned(view.into_iter().map(|a| a.as_()).collect()))
                    }
                    Cow::Owned(owned) => {
                        let res = owned.iter().map(|a| a.as_()).collect();
                        Ok(Cow::Owned(res))
                    }
                }
            }
            Err(err) => Err(err),
        }
    }

    fn dtype(&self) -> BinaryDataArrayType;

    fn to_f32(&'lifespan self) -> Result<Cow<'transient, [f32]>, ArrayRetrievalError> {
        type D = f32;
        match self.dtype() {
            BinaryDataArrayType::Float32 | BinaryDataArrayType::ASCII => self.coerce::<D>(),
            BinaryDataArrayType::Float64 => {
                type S = f64;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Int32 => {
                type S = i32;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Int64 => {
                type S = i64;
                self.convert::<S, D>()
            }
            _ => Err(ArrayRetrievalError::DataTypeSizeMismatch),
        }
    }

    fn to_f64(&'lifespan self) -> Result<Cow<'transient, [f64]>, ArrayRetrievalError> {
        type D = f64;
        match self.dtype() {
            BinaryDataArrayType::Float32 => {
                type S = f32;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Float64 | BinaryDataArrayType::ASCII => self.coerce(),
            BinaryDataArrayType::Int32 => {
                type S = i32;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Int64 => {
                type S = i64;
                self.convert::<S, D>()
            }
            _ => Err(ArrayRetrievalError::DataTypeSizeMismatch),
        }
    }

    fn to_i32(&'lifespan self) -> Result<Cow<'transient, [i32]>, ArrayRetrievalError> {
        type D = i32;
        match self.dtype() {
            BinaryDataArrayType::Float32 => {
                type S = f32;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Float64 => {
                type S = f64;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Int32 | BinaryDataArrayType::ASCII => self.coerce::<D>(),
            BinaryDataArrayType::Int64 => {
                type S = i64;
                self.convert::<S, D>()
            }
            _ => Err(ArrayRetrievalError::DataTypeSizeMismatch),
        }
    }

    fn to_i64(&'lifespan self) -> Result<Cow<'transient, [i64]>, ArrayRetrievalError> {
        type D = i64;
        match self.dtype() {
            BinaryDataArrayType::Float32 => {
                type S = f32;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Float64 => {
                type S = f64;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Int64 | BinaryDataArrayType::ASCII => self.coerce::<D>(),
            BinaryDataArrayType::Int32 => {
                type S = i32;
                self.convert::<S, D>()
            }
            _ => Err(ArrayRetrievalError::DataTypeSizeMismatch),
        }
    }
}

pub trait ByteArrayViewMut<'transient, 'lifespan: 'transient>:
    ByteArrayView<'transient, 'lifespan>
{
    fn view_mut(&'transient mut self) -> Result<&'transient mut Bytes, ArrayRetrievalError>;

    fn coerce_from_mut<T: Clone + Sized>(
        buffer: &mut [u8],
    ) -> Result<&'transient mut [T], ArrayRetrievalError> {
        let n = buffer.len();
        let z = mem::size_of::<T>();
        if n % z != 0 {
            return Err(ArrayRetrievalError::DataTypeSizeMismatch);
        }
        let m = n / z;
        unsafe { Ok(slice::from_raw_parts_mut(buffer.as_ptr() as *mut T, m)) }
    }

    fn coerce_mut<T: Clone + Sized>(
        &'lifespan mut self,
    ) -> Result<&'transient mut [T], ArrayRetrievalError> {
        let view = match self.view_mut() {
            Ok(data) => data,
            Err(err) => return Err(err),
        };
        Self::coerce_from_mut(view)
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

    pub fn slice(&self, start: usize, end: usize) -> Result<DataArray, ArrayRetrievalError> {
        if end < start || (end - start) % self.dtype.size_of() != 0 {
            Err(ArrayRetrievalError::DataTypeSizeMismatch)
        } else {
            let data = self.decode()?;
            let slice = data[start..end].to_vec();
            let subset = Self::wrap(&self.name, self.dtype, slice);
            Ok(subset)
        }
    }

    pub fn slice_buffer(&self, start: usize, end: usize) -> Result<Cow<'_, [u8]>, ArrayRetrievalError> {
        if end < start || (end - start) % self.dtype.size_of() != 0 {
            Err(ArrayRetrievalError::DataTypeSizeMismatch)
        } else {
            let data = self.decode()?;
            match data {
                Cow::Borrowed(view) => {
                    Ok(Cow::Borrowed(&view[start..end]))
                },
                Cow::Owned(view) => {
                    let data = &view[start..end];
                    Ok(Cow::Owned(data.to_owned()))
                }
            }
        }
    }

    pub fn wrap(name: &ArrayType, dtype: BinaryDataArrayType, data: Bytes) -> DataArray {
        DataArray {
            dtype,
            name: name.clone(),
            data,
            compression: BinaryCompressionType::Decoded,
            ..Default::default()
        }
    }

    pub fn update_buffer<T: Pod>(&mut self, data_buffer: &[T]) -> Result<usize, ArrayRetrievalError> {
        if self.dtype.size_of() != mem::size_of::<T>() {
            Err(ArrayRetrievalError::DataTypeSizeMismatch)
        } else {
            self.data = to_bytes(data_buffer);
            Ok(self.data.len())
        }
    }

    pub fn encode_bytestring(&self, compression: BinaryCompressionType) -> Bytes {
        let bytestring = match self.compression {
            BinaryCompressionType::Decoded => Cow::Borrowed(self.data.as_slice()),
            _ => self.decode().expect("Failed to decode binary data"),
        };
        match compression {
            BinaryCompressionType::Zlib => {
                let compressed = Self::compress_zlib(&bytestring);
                Base64Std.encode(compressed).into()
            }
            BinaryCompressionType::NoCompression => Base64Std.encode(bytestring.as_ref()).into(),
            BinaryCompressionType::Decoded => panic!("Should never happen"),
            _ => {
                panic!("Compresion type {:?} is unsupported", compression)
            }
        }
    }

    pub fn encode(&self, compression: BinaryCompressionType) -> DataArray {
        let mut dup = self.clone();
        let bytestring = match self.compression {
            BinaryCompressionType::Decoded => dup.data,
            _ => {
                dup.decode_and_store()
                    .expect("Failed to decode binary data array");
                dup.data
            }
        };
        dup.compression = BinaryCompressionType::NoCompression;

        dup.data = match compression {
            BinaryCompressionType::Zlib => {
                let compressed = Self::compress_zlib(&bytestring);
                Base64Std.encode(compressed).into()
            }
            BinaryCompressionType::NoCompression => Base64Std.encode(bytestring).into(),
            BinaryCompressionType::Decoded => panic!("Should never happen"),
            _ => {
                panic!("Compresion type {:?} is unsupported", compression)
            }
        };
        dup.compression = compression;
        dup
    }

    pub fn compress_zlib(bytestring: &[u8]) -> Bytes {
        let result = Bytes::new();
        let mut compressor = ZlibEncoder::new(result, Compression::best());
        compressor.write_all(bytestring).expect("Error compressing");
        compressor.finish().expect("Error compressing")
    }

    pub fn decompres_zlib(bytestring: &[u8]) -> Bytes {
        let result = Bytes::new();
        let mut decompressor = ZlibDecoder::new(result);
        decompressor
            .write_all(bytestring)
            .expect("Error decompression");
        decompressor.finish().expect("Error decompression")
    }

    #[cfg(feature = "numpress")]
    pub fn compress_numpress_linear(data: &[f64]) -> Result<Bytes, ArrayRetrievalError> {
        let scaling = numpress::optimal_scaling(data);
        match numpress::numpress_compress(data, scaling) {
            Ok(data) => Ok(data),
            Err(e) => {
                Err(ArrayRetrievalError::DecompressionError(e.to_string()))
            },
        }
    }

    #[cfg(feature = "numpress")]
    pub fn decompres_numpress_linear(data: &[u8]) -> Result<Vec<f64>, ArrayRetrievalError> {
        match numpress::numpress_decompress(data) {
            Ok(data) => Ok(data),
            Err(e) => {
                Err(ArrayRetrievalError::DecompressionError(e.to_string()))
            },
        }
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

    pub fn decode(&'lifespan self) -> Result<Cow<'lifespan, [u8]>, ArrayRetrievalError> {
        match self.compression {
            BinaryCompressionType::Decoded => Ok(Cow::Borrowed(&self.data.as_slice())),
            BinaryCompressionType::NoCompression => {
                let bytestring = Base64Std
                    .decode(&self.data)
                    .expect("Failed to decode base64 array");
                Ok(Cow::Owned(bytestring))
            }
            BinaryCompressionType::Zlib => {
                let bytestring = Base64Std
                    .decode(&self.data)
                    .expect("Failed to decode base64 array");
                Ok(Cow::Owned(Self::decompres_zlib(&bytestring)))
            },
            #[cfg(feature = "numpress")]
            BinaryCompressionType::NumpressLinear => {
                match self.dtype {
                    BinaryDataArrayType::Float64 => {
                        let mut bytestring = Base64Std
                            .decode(&self.data)
                            .expect("Failed to decode base64 array");
                        let decoded = Self::decompres_numpress_linear(&mut bytestring)?;
                        let view = vec_as_bytes(decoded);
                        Ok(Cow::Owned(view))
                    }
                    _ => {
                        Err(ArrayRetrievalError::DecompressionError(
                            self.compression.unsupported_msg(Some(
                                format!("Not compatible with {:?}", self.dtype).as_str(),
                            )),
                        ))
                    }
                }
            }
            mode => Err(ArrayRetrievalError::DecompressionError(format!("Cannot decode array encoded with {:?}", mode))),
        }
    }

    pub(crate) fn decoded_slice(&'lifespan self, start: usize, end: usize) -> Result<Cow<'lifespan, [u8]>, ArrayRetrievalError> {
        if start > end || (end - start) % self.dtype.size_of() != 0 {
            return Err(ArrayRetrievalError::DataTypeSizeMismatch)
        }
        match self.compression {
            BinaryCompressionType::Decoded => {
                Ok(Cow::Borrowed(&self.data.as_slice()[start..end]))
            },
            BinaryCompressionType::NoCompression => {
                let bytestring = Base64Std
                    .decode(&self.data)
                    .expect("Failed to decode base64 array");
                Ok(Cow::Owned(bytestring[start..end].to_vec()))
            }
            BinaryCompressionType::Zlib => {
                let bytestring = Base64Std
                    .decode(&self.data)
                    .expect("Failed to decode base64 array");
                Ok(Cow::Owned(Self::decompres_zlib(&bytestring)[start..end].to_vec()))
            }
            mode => Err(ArrayRetrievalError::DecompressionError(format!("Cannot decode array compressed with {:?}", mode))),
        }
    }

    pub fn decode_mut(&'transient mut self) -> Result<&'transient mut Bytes, ArrayRetrievalError> {
        match self.compression {
            BinaryCompressionType::Decoded => Ok(&mut self.data),
            BinaryCompressionType::NoCompression => {
                let bytestring = Base64Std
                    .decode(&self.data)
                    .expect("Failed to decode base64 array");
                self.data = bytestring;
                self.compression = BinaryCompressionType::Decoded;
                Ok(&mut self.data)
            }
            BinaryCompressionType::Zlib => {
                let bytestring = Base64Std
                    .decode(&self.data)
                    .expect("Failed to decode base64 array");
                self.data = bytestring;
                self.compression = BinaryCompressionType::Decoded;
                Ok(&mut self.data)
            }
            mode => Err(ArrayRetrievalError::DecompressionError(format!("Cannot decode array compressed with {:?}", mode))),
        }
    }

    pub fn clear(&mut self) {
        self.data.clear();
        self.params = None;
    }

    pub fn store_as(&mut self, dtype: BinaryDataArrayType) -> Result<usize, ArrayRetrievalError> {
        if self.dtype == dtype {
            return Ok(self.data.len());
        }
        let result = match dtype {
            BinaryDataArrayType::Float32 => {
                let view = match self.to_f32() {
                    Ok(view) => view,
                    Err(err) => return Err(err),
                };
                let recast = to_bytes(&view);
                self.update_buffer(&recast)
            }
            BinaryDataArrayType::Float64 => {
                let view = match self.to_f64() {
                    Ok(view) => view,
                    Err(err) => return Err(err),
                };
                let recast = to_bytes(&view);
                self.update_buffer(&recast)
            }
            BinaryDataArrayType::Int32 => {
                let view = match self.to_i32() {
                    Ok(view) => view,
                    Err(err) => return Err(err),
                };
                let recast = to_bytes(&view);
                self.update_buffer(&recast)
            }
            BinaryDataArrayType::Int64 => {
                let view = match self.to_i64() {
                    Ok(view) => view,
                    Err(err) => return Err(err),
                };
                let recast = to_bytes(&view);
                self.update_buffer(&recast)
            }
            _ => Ok(0),
        };
        self.dtype = dtype;
        result
    }
}

impl<'transient, 'lifespan: 'transient> ByteArrayView<'transient, 'lifespan> for DataArray {
    fn view(&'lifespan self) -> Result<Cow<'lifespan, [u8]>, ArrayRetrievalError> {
        self.decode()
    }

    fn dtype(&self) -> BinaryDataArrayType {
        self.dtype
    }
}

impl<'transient, 'lifespan: 'transient> ByteArrayViewMut<'transient, 'lifespan> for DataArray {
    fn view_mut(&'transient mut self) -> Result<&'transient mut Bytes, ArrayRetrievalError> {
        self.decode_mut()
    }
}



impl_param_described_deferred!(DataArray);


#[derive(Clone, Debug)]
pub struct DataArraySlice<'a> {
    source: &'a DataArray,
    pub start: usize,
    pub end: usize,
}

impl<'a> DataArraySlice<'a> {
    pub fn new(source: &'a DataArray, mut start: usize, mut end: usize) -> Self {
        if start > end {
            mem::swap(&mut start, &mut end);
        }
        Self { source, start, end }
    }


    pub fn decode(&'a self) -> Result<Cow<'a, [u8]>, ArrayRetrievalError> {
        self.source.decoded_slice(self.start, self.end)
    }

}

impl<'transient, 'lifespan: 'transient> ByteArrayView<'transient, 'lifespan>
    for DataArraySlice<'lifespan>
{
    fn view(&'lifespan self) -> Result<Cow<'lifespan, [u8]>, ArrayRetrievalError> {
        self.decode()
    }

    fn dtype(&self) -> BinaryDataArrayType {
        self.source.dtype()
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

    pub fn len(&self) -> usize {
        self.byte_buffer_map.len()
    }

    pub fn is_empty(&self) -> bool {
        self.byte_buffer_map.is_empty()
    }

    pub fn iter(&self) -> std::collections::hash_map::Iter<ArrayType, DataArray> {
        self.byte_buffer_map.iter()
    }

    pub fn iter_mut(&mut self) -> std::collections::hash_map::IterMut<ArrayType, DataArray> {
        self.byte_buffer_map.iter_mut()
    }

    pub fn decode(&mut self) -> Result<(), ArrayRetrievalError> {
        for (_key, value) in self.iter_mut() {
            match value.compression {
                BinaryCompressionType::Decoded => {}
                _ => {
                    value.decode_and_store()?;
                }
            }
        }
        Ok(())
    }

    pub fn add(&mut self, array: DataArray) {
        self.byte_buffer_map.insert(array.name.clone(), array);
    }

    pub fn get(&'transient self, array_type: &ArrayType) -> Option<&'transient DataArray> {
        self.byte_buffer_map.get(array_type)
    }

    pub fn get_mut(&mut self, array_type: &ArrayType) -> Option<&mut DataArray> {
        self.byte_buffer_map.get_mut(array_type)
    }

    pub fn has_array(&self, array_type: &ArrayType) -> bool {
        self.byte_buffer_map.contains_key(array_type)
    }

    pub fn clear(&mut self) {
        self.byte_buffer_map.clear();
    }

    pub fn search(&self, query: f64, error_tolerance: Tolerance) -> Option<usize> {
        let mzs = self.mzs();
        let (lower, _upper) = error_tolerance.bounds(query);
        match mzs[..].binary_search_by(|m| m.partial_cmp(&lower).unwrap()) {
            Ok(i) => {
                let mut best_error = error_tolerance.call(query, mzs[i]).abs();
                let mut best_index = i;
                let mut index = i + 1;
                while index < mzs.len() {
                    let error = error_tolerance.call(query, mzs[index]).abs();
                    if error < best_error {
                        best_index = index;
                        best_error = error;
                    }
                    index += 1;
                }
                if best_error < error_tolerance.tol() {
                    return Some(best_index);
                }
                None
            }
            Err(_err) => None,
        }
    }

    pub fn mzs(&'transient self) -> Cow<'transient, [f64]> {
        let mz_array = self
            .get(&ArrayType::MZArray)
            .expect("Did not find m/z array")
            .to_f64()
            .expect("Failed to decode m/z array");
        mz_array
    }

    pub fn intensities(&'transient self) -> Cow<'transient, [f32]> {
        let intensities = self
            .get(&ArrayType::IntensityArray)
            .expect("Did not find intensity array")
            .to_f32()
            .expect("Failed to decode intensity array");
        intensities
    }

    pub fn charges(&'lifespan self) -> Option<Cow<'transient, [i32]>> {
        match self.get(&ArrayType::ChargeArray) {
            Some(data_array) => match data_array.to_i32() {
                Ok(array) => Some(array),
                Err(err) => {
                    warn!("Failed to decode charge state array: {:?}", err);
                    None
                }
            },
            None => None,
        }
    }
}

impl From<&PeakSet> for BinaryArrayMap {
    fn from(peaks: &PeakSet) -> BinaryArrayMap {
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
            let mz: f64 = p.coordinate();
            let inten: f32 = p.intensity();

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

impl<C: CentroidLike + From<CentroidPeak>> From<BinaryArrayMap> for MZPeakSetType<C> {
    fn from(arrays: BinaryArrayMap) -> MZPeakSetType<C> {
        (&arrays).into()
    }
}

impl<C: CentroidLike + From<CentroidPeak>> From<&BinaryArrayMap> for MZPeakSetType<C> {
    fn from(arrays: &BinaryArrayMap) -> MZPeakSetType<C> {
        let mz_array = arrays.mzs();
        let intensity_array = arrays.intensities();
        let mut peaks = Vec::with_capacity(mz_array.len());

        for (i, (mz, intensity)) in mz_array.iter().zip(intensity_array.iter()).enumerate() {
            peaks.push(
                CentroidPeak {
                    mz: *mz,
                    intensity: *intensity,
                    index: i as u32,
                }
                .into(),
            )
        }

        MZPeakSetType::<C>::new(peaks)
    }
}

impl From<&BinaryArrayMap> for DeconvolutedPeakSet {
    fn from(arrays: &BinaryArrayMap) -> DeconvolutedPeakSet {
        let mz_array = arrays.mzs();
        let intensity_array = arrays.intensities();
        let charge_array = arrays
            .charges()
            .expect("Charge state array is required for deconvoluted peaks");
        let mut peaks = Vec::with_capacity(mz_array.len());
        for (i, ((mz, intensity), charge)) in mz_array
            .iter()
            .zip(intensity_array.iter())
            .zip(charge_array.iter())
            .enumerate()
        {
            peaks.push(DeconvolutedPeak {
                neutral_mass: neutral_mass(*mz, *charge),
                intensity: *intensity,
                charge: *charge,
                index: i as u32,
            })
        }

        DeconvolutedPeakSet::new(peaks)
    }
}

impl From<&DeconvolutedPeakSet> for BinaryArrayMap {
    fn from(peaks: &DeconvolutedPeakSet) -> BinaryArrayMap {
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

        let mut charge_array = DataArray::from_name_type_size(
            &ArrayType::ChargeArray,
            BinaryDataArrayType::Int32,
            peaks.len() * BinaryDataArrayType::Int32.size_of(),
        );

        mz_array.compression = BinaryCompressionType::Decoded;
        intensity_array.compression = BinaryCompressionType::Decoded;
        charge_array.compression = BinaryCompressionType::Decoded;

        for p in peaks.iter() {
            let mz: f64 = p.mz();
            let inten: f32 = p.intensity();
            let charge = p.charge();

            let raw_bytes: [u8; mem::size_of::<f64>()] = unsafe { mem::transmute(mz) };
            mz_array.data.extend(raw_bytes);

            let raw_bytes: [u8; mem::size_of::<f32>()] = unsafe { mem::transmute(inten) };
            intensity_array.data.extend(raw_bytes);

            let raw_bytes: [u8; mem::size_of::<i32>()] = unsafe { mem::transmute(charge) };
            charge_array.data.extend(raw_bytes);
        }

        arrays.add(mz_array);
        arrays.add(intensity_array);
        arrays.add(charge_array);
        arrays
    }
}
