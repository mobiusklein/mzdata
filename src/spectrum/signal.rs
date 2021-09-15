use std::borrow::Cow;
use std::collections::HashMap;
use std::fmt;
use std::fmt::{Debug, Formatter};
use std::io::prelude::*;
use std::mem;
use std::slice;

use num_traits::cast::AsPrimitive;
use num_traits::Num;

use log::warn;

use base64;
use flate2::write::{ZlibDecoder, ZlibEncoder};
use flate2::Compression;
use mzpeaks::prelude::*;
use mzpeaks::{
    CentroidPeak, DeconvolutedPeak, DeconvolutedPeakSet, MZPeakSetType, MassErrorType, PeakSet,
};

use crate::params::{ParamList, Unit};
use crate::utils::neutral_mass;

type Bytes = Vec<u8>;

pub fn to_bytes<T>(data: &[T]) -> Bytes {
    let n = data.len();
    let z = mem::size_of::<T>();
    let m = n * z;
    unsafe {
        let byte_buffer = slice::from_raw_parts(data.as_ptr() as *const u8, m);
        let mut result = Bytes::new();
        result.copy_from_slice(byte_buffer);
        result
    }
}

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

#[derive(Debug, Clone, PartialEq, Hash, Eq, PartialOrd, Ord)]
pub enum ArrayType {
    Unknown,
    MZArray,
    IntensityArray,
    ChargeArray,
    TimeArray,
    WavelengthArray,
    IonMobilityArray(Unit),
    MeanIonMobilityArray(Unit),
    RawIonMobilityArray(Unit),
    DeconvolutedIonMobilityArray(Unit),
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

impl std::fmt::Display for ArrayRetrievalError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::error::Error for ArrayRetrievalError {}

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

    pub fn wrap(name: &ArrayType, dtype: BinaryDataArrayType, data: Bytes) -> DataArray {
        DataArray {
            dtype,
            name: name.clone(),
            data,
            compression: BinaryCompressionType::Decoded,
            ..Default::default()
        }
    }

    pub fn update_buffer<T>(&mut self, data_buffer: &[T]) -> Result<usize, ArrayRetrievalError> {
        if self.dtype.size_of() != mem::size_of::<T>() {
            Err(ArrayRetrievalError::DataTypeSizeMismatch)
        } else {
            self.data = to_bytes(data_buffer);
            Ok(self.data.len())
        }
    }

    pub fn encode_bytestring(&self, compression: BinaryCompressionType) -> Bytes {
        let bytestring = match self.compression {
            BinaryCompressionType::Decoded => Cow::Borrowed(&self.data),
            _ => self.decode().expect("Failed to decode binary data"),
        };
        match compression {
            BinaryCompressionType::Zlib => {
                let compressed = Self::compress_zlib(&bytestring);
                base64::encode(compressed).into()
            }
            BinaryCompressionType::NoCompression => base64::encode(bytestring.as_ref()).into(),
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
                base64::encode(compressed).into()
            }
            BinaryCompressionType::NoCompression => base64::encode(bytestring).into(),
            BinaryCompressionType::Decoded => panic!("Should never happen"),
            _ => {
                panic!("Compresion type {:?} is unsupported", compression)
            }
        };
        dup.compression = compression;
        dup
    }

    fn compress_zlib(bytestring: &[u8]) -> Bytes {
        let result = Bytes::new();
        let mut compressor = ZlibEncoder::new(result, Compression::best());
        compressor.write_all(bytestring).expect("Error compressing");
        compressor.finish().expect("Error compressing")
    }

    fn decompres_zlib(bytestring: &[u8]) -> Bytes {
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
                Ok(Cow::Owned(Self::decompres_zlib(&bytestring)))
            }
            _ => Err(ArrayRetrievalError::DecompressionError),
        }
    }

    pub fn decode_mut(&'transient mut self) -> Result<&'transient mut Bytes, ArrayRetrievalError> {
        match self.compression {
            BinaryCompressionType::Decoded => Ok(&mut self.data),
            BinaryCompressionType::NoCompression => {
                let bytestring = base64::decode(&self.data).expect("Failed to decode base64 array");
                self.data = bytestring;
                self.compression = BinaryCompressionType::Decoded;
                Ok(&mut self.data)
            }
            BinaryCompressionType::Zlib => {
                let bytestring = base64::decode(&self.data).expect("Failed to decode base64 array");
                self.data = bytestring;
                self.compression = BinaryCompressionType::Decoded;
                Ok(&mut self.data)
            }
            _ => Err(ArrayRetrievalError::DecompressionError),
        }
    }

    pub fn coerce_from<T: Clone + Sized>(
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

    pub fn coerce_from_mut<T: Clone + Sized>(
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

    pub fn coerce_mut<T: Clone + Sized>(
        &'lifespan mut self,
    ) -> Result<&'transient mut [T], ArrayRetrievalError> {
        let view = match self.decode_mut() {
            Ok(data) => data,
            Err(err) => return Err(err),
        };
        Self::coerce_from_mut(view)
    }

    pub fn coerce<T: Clone + Sized>(
        &'lifespan self,
    ) -> Result<Cow<'transient, [T]>, ArrayRetrievalError> {
        match self.decode() {
            Ok(data) => Self::coerce_from(data),
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

    pub fn to_i64(&'lifespan self) -> Result<Cow<'transient, [i64]>, ArrayRetrievalError> {
        type D = i64;
        match self.dtype {
            BinaryDataArrayType::Float32 => {
                type S = f32;
                self.transmute::<S, D>()
            }
            BinaryDataArrayType::Float64 => {
                type S = f64;
                self.transmute::<S, D>()
            }
            BinaryDataArrayType::Int64 => self.coerce::<D>(),
            BinaryDataArrayType::Int32 => {
                type S = i32;
                self.transmute::<S, D>()
            }
            _ => Err(ArrayRetrievalError::DataTypeSizeMismatch),
        }
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

    pub fn search(
        &self,
        query: f64,
        error_tolerance: f64,
        error_type: MassErrorType,
    ) -> Option<usize> {
        let mzs = self.mzs();
        let lower = error_type.lower_bound(query, error_tolerance);
        match mzs[..].binary_search_by(|m| m.partial_cmp(&lower).unwrap()) {
            Ok(i) => {
                let mut best_error = error_type.call(query, mzs[i]).abs();
                let mut best_index = i;
                let mut index = i + 1;
                while index < mzs.len() {
                    let error = error_type.call(query, mzs[index]).abs();
                    if error < best_error {
                        best_index = index;
                        best_error = error;
                    }
                    index += 1;
                }
                if best_error < error_tolerance {
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

impl From<PeakSet> for BinaryArrayMap {
    fn from(peaks: PeakSet) -> BinaryArrayMap {
        (&peaks).into()
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

impl From<DeconvolutedPeakSet> for BinaryArrayMap {
    fn from(peaks: DeconvolutedPeakSet) -> BinaryArrayMap {
        (&peaks).into()
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
