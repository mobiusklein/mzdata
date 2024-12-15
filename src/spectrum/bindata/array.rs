use std::borrow::Cow;
use std::fmt::{self, Formatter};
use std::io::prelude::*;
use std::mem;

use base64_simd;
use bytemuck::Pod;
use flate2::write::{ZlibDecoder, ZlibEncoder};
use flate2::Compression;

use crate::params::{ParamList, Unit};

use super::encodings::{
    to_bytes, ArrayRetrievalError, ArrayType, BinaryCompressionType,
    BinaryDataArrayType, Bytes,
};
use super::traits::{ByteArrayView, ByteArrayViewMut};
#[allow(unused)]
use super::vec_as_bytes;

/// Represents a data array that holds a byte buffer that may be compressed, base64 encoded,
/// or raw little endian bytes, and provides views of those bytes as a small range of supported
/// types.
///
/// This type is modeled after the `<binaryDataArray>` element in mzML.
///
/// # Note
/// This type tries to walk a fine line between convenience and performance and as such is easy
/// to misuse. All operations that view the byte buffer as arbitrary data need that data to be
/// decoded and decompressed in order to borrow it. If the byte buffer is not already stored
/// decoded, the operation will copy and decode the buffer in its entirety before performing any
/// other operation, so repeated method calls may incur excessive overhead. If this is happening,
/// please use [`DataArray::decode_and_store`] to store the decoded representation explicitly.
///
/// Normally, [`SpectrumSource`](crate::io::SpectrumSource)-implementing file readers will eagerly decode all arrays
/// as soon as they are ready. If they are operating in lazy mode, the buffers will need to be decoded
/// explicitly, again using [`DataArray::decode_and_store`] or operations should make as much use of the
/// copied arrays as possible instead.
#[derive(Default, Clone)]
pub struct DataArray {
    pub data: Bytes,
    pub dtype: BinaryDataArrayType,
    pub compression: BinaryCompressionType,
    pub name: ArrayType,
    pub params: Option<Box<ParamList>>,
    pub unit: Unit,
    item_count: Option<usize>
}

impl core::fmt::Debug for DataArray {
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

const EMPTY_BUFFER: [u8; 0] = [];

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

    pub fn slice_buffer(
        &self,
        start: usize,
        end: usize,
    ) -> Result<Cow<'_, [u8]>, ArrayRetrievalError> {
        if end < start || (end - start) % self.dtype.size_of() != 0 {
            Err(ArrayRetrievalError::DataTypeSizeMismatch)
        } else {
            let data = self.decode()?;
            match data {
                Cow::Borrowed(view) => Ok(Cow::Borrowed(&view[start..end])),
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

    fn set_buffer_of_type(&mut self, data_buffer: Vec<u8>) -> Result<usize, ArrayRetrievalError> {
        self.item_count = Some(data_buffer.len() / self.dtype().size_of());
        self.data = data_buffer;
        Ok(self.data.len())
    }

    pub fn update_buffer<T: Pod>(
        &mut self,
        data_buffer: &[T],
    ) -> Result<usize, ArrayRetrievalError> {
        if self.dtype.size_of() != mem::size_of::<T>() {
            Err(ArrayRetrievalError::DataTypeSizeMismatch)
        } else {
            self.item_count = Some(data_buffer.len());
            self.data = to_bytes(data_buffer);
            Ok(self.data.len())
        }
    }

    pub fn push<T: Pod>(&mut self, value: T) -> Result<(), ArrayRetrievalError> {
        if !matches!(self.compression, BinaryCompressionType::Decoded) {
            self.decode_and_store()?;
        };
        if self.dtype.size_of() != mem::size_of::<T>() {
            Err(ArrayRetrievalError::DataTypeSizeMismatch)
        } else {
            let data = bytemuck::bytes_of(&value);
            self.data.extend(data.iter());
            self.item_count = self.item_count.map(|i| i + 1);
            Ok(())
        }
    }

    pub fn extend<T: Pod>(&mut self, values: &[T]) -> Result<(), ArrayRetrievalError> {
        if !matches!(self.compression, BinaryCompressionType::Decoded) {
            self.decode_and_store()?;
        };
        if self.dtype.size_of() != mem::size_of::<T>() {
            Err(ArrayRetrievalError::DataTypeSizeMismatch)
        } else {
            self.item_count = self.item_count.map(|i| i + values.len());
            let data = bytemuck::cast_slice(values);
            self.data.extend(data.iter());
            Ok(())
        }
    }

    pub fn encode_bytestring(&self, compression: BinaryCompressionType) -> Bytes {
        if self.compression == compression {
            log::debug!("Fast-path encoding {}:{}", self.name, self.dtype);
            return self.data.clone()
        }
        let bytestring = match self.compression {
            BinaryCompressionType::Decoded => Cow::Borrowed(self.data.as_slice()),
            _ => self.decode().expect("Failed to decode binary data"),
        };
        match compression {
            BinaryCompressionType::Zlib => {
                let compressed = Self::compress_zlib(&bytestring);
                base64_simd::STANDARD.encode_type::<Bytes>(&compressed)
            }
            BinaryCompressionType::NoCompression => base64_simd::STANDARD.encode_type::<Bytes>(bytestring.as_ref()),
            BinaryCompressionType::Decoded => panic!("Should never happen"),
            _ => {
                panic!("Compresion type {:?} is unsupported", compression)
            }
        }
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
            .unwrap_or_else(|e| panic!("Decompression error: {}", e));
        decompressor.finish().unwrap_or_else(|e| panic!("Decompression error: {}", e))
    }

    #[cfg(feature = "numpress")]
    pub fn compress_numpress_linear(data: &[f64]) -> Result<Bytes, ArrayRetrievalError> {
        let scaling = numpress::optimal_scaling(data);
        match numpress::numpress_compress(data, scaling) {
            Ok(data) => Ok(data),
            Err(e) => Err(ArrayRetrievalError::DecompressionError(e.to_string())),
        }
    }

    #[cfg(feature = "numpress")]
    pub fn decompres_numpress_linear(data: &[u8]) -> Result<Vec<f64>, ArrayRetrievalError> {
        match numpress::numpress_decompress(data) {
            Ok(data) => Ok(data),
            Err(e) => Err(ArrayRetrievalError::DecompressionError(e.to_string())),
        }
    }

    /// Decode the compressed data, if needed, and store that buffer in `self.data`. After
    /// decoding `self.compression` will always be [`BinaryCompressionType::Decoded`].
    ///
    /// The return value is the content of `self.compression` after decoding.
    ///
    /// This may fail if the decoding fails for any reason.
    pub fn decode_and_store(&mut self) -> Result<BinaryCompressionType, ArrayRetrievalError> {
        match self.decode() {
            Ok(data) => {
                match data {
                    // The only time this is a borrow is when the data are already
                    // decoded.
                    Cow::Borrowed(_view) => Ok(self.compression),
                    Cow::Owned(buffer) => {
                        self.item_count = Some(buffer.len() / self.dtype.size_of());
                        self.data = buffer;
                        self.compression = BinaryCompressionType::Decoded;
                        Ok(self.compression)
                    }
                }
            }
            Err(err) => Err(err),
        }
    }

    /// Decompress and base64-decode encoded bytes, and return the data.
    ///
    /// If the data were already decoded, the existing bytes are returned. Otherwise one or
    /// more buffers may be allocated to hold the decompressed and decoded bytes.
    pub fn decode(&'lifespan self) -> Result<Cow<'lifespan, [u8]>, ArrayRetrievalError> {
        if self.data.is_empty() {
            return Ok(Cow::Borrowed(&EMPTY_BUFFER))
        }
        match self.compression {
            BinaryCompressionType::Decoded => Ok(Cow::Borrowed(self.data.as_slice())),
            BinaryCompressionType::NoCompression => {
                let bytestring = base64_simd::STANDARD.decode_type::<Bytes>(&self.data)
                    .unwrap_or_else(|e| panic!("Failed to decode base64 array: {}", e));
                Ok(Cow::Owned(bytestring))
            }
            BinaryCompressionType::Zlib => {
                let bytestring = base64_simd::STANDARD.decode_type::<Bytes>(&self.data)
                    .unwrap_or_else(|e| panic!("Failed to decode base64 array: {}", e));
                Ok(Cow::Owned(Self::decompres_zlib(&bytestring)))
            }
            #[cfg(feature = "numpress")]
            BinaryCompressionType::NumpressLinear => match self.dtype {
                BinaryDataArrayType::Float64 => {
                    let mut bytestring = base64_simd::STANDARD.decode_type::<Bytes>(&self.data)
                        .unwrap_or_else(|e| panic!("Failed to decode base64 array: {}", e));
                    let decoded = Self::decompres_numpress_linear(&mut bytestring)?;
                    let view = vec_as_bytes(decoded);
                    Ok(Cow::Owned(view))
                }
                _ => Err(ArrayRetrievalError::DecompressionError(
                    self.compression.unsupported_msg(Some(
                        format!("Not compatible with {:?}", self.dtype).as_str(),
                    )),
                )),
            },
            mode => Err(ArrayRetrievalError::DecompressionError(format!(
                "Cannot decode array encoded with {:?}",
                mode
            ))),
        }
    }

    pub(crate) fn decoded_slice(
        &'lifespan self,
        start: usize,
        end: usize,
    ) -> Result<Cow<'lifespan, [u8]>, ArrayRetrievalError> {
        if start > end || (end - start) % self.dtype.size_of() != 0 {
            return Err(ArrayRetrievalError::DataTypeSizeMismatch);
        }
        match self.compression {
            BinaryCompressionType::Decoded => Ok(Cow::Borrowed(&self.data.as_slice()[start..end])),
            BinaryCompressionType::NoCompression => {
                let bytestring = base64_simd::STANDARD.decode_type::<Bytes>(&self.data)
                    .unwrap_or_else(|e| panic!("Failed to decode base64 array: {}", e));
                Ok(Cow::Owned(bytestring[start..end].to_vec()))
            }
            BinaryCompressionType::Zlib => {
                let bytestring = base64_simd::STANDARD.decode_type::<Bytes>(&self.data)
                    .unwrap_or_else(|e| panic!("Failed to decode base64 array: {}", e));
                Ok(Cow::Owned(
                    Self::decompres_zlib(&bytestring)[start..end].to_vec(),
                ))
            }
            mode => Err(ArrayRetrievalError::DecompressionError(format!(
                "Cannot decode array slice compressed with {:?}",
                mode
            ))),
        }
    }

    pub fn decode_mut(&'transient mut self) -> Result<&'transient mut Bytes, ArrayRetrievalError> {
        if self.data.is_empty() {
            return Ok(&mut self.data)
        }
        match self.compression {
            BinaryCompressionType::Decoded => Ok(&mut self.data),
            BinaryCompressionType::NoCompression => {
                let bytestring = base64_simd::STANDARD.decode_type::<Bytes>(&self.data)
                    .unwrap_or_else(|e| panic!("Failed to decode base64 array: {}", e));
                self.data = bytestring;
                self.compression = BinaryCompressionType::Decoded;
                Ok(&mut self.data)
            }
            BinaryCompressionType::Zlib => {
                let bytestring = base64_simd::STANDARD.decode_type::<Bytes>(&self.data)
                    .unwrap_or_else(|e| panic!("Failed to decode base64 array: {}", e));
                self.data = bytestring;
                self.compression = BinaryCompressionType::Decoded;
                Ok(&mut self.data)
            },
            #[cfg(feature = "numpress")]
            BinaryCompressionType::NumpressLinear => match self.dtype {
                BinaryDataArrayType::Float64 => {
                    let mut bytestring = base64_simd::STANDARD.decode_type::<Bytes>(&self.data)
                        .unwrap_or_else(|e| panic!("Failed to decode base64 array: {}", e));
                    let decoded = Self::decompres_numpress_linear(&mut bytestring)?;
                    let view = vec_as_bytes(decoded);
                    self.data = view;
                    Ok(&mut self.data)
                }
                _ => Err(ArrayRetrievalError::DecompressionError(
                    self.compression.unsupported_msg(Some(
                        format!("Not compatible with {:?}", self.dtype).as_str(),
                    )),
                )),
            },
            mode => Err(ArrayRetrievalError::DecompressionError(format!(
                "Cannot decode array compressed with {:?}",
                mode
            ))),
        }
    }

    pub fn clear(&mut self) {
        self.data.clear();
        self.params = None;
        self.item_count = None;
    }

    /// The reverse of [`DataArray::decode_and_store`], this method compresses `self.data` to the desired
    /// compression method and stores that buffer as `self.data`.
    pub fn store_compressed(&mut self, compression: BinaryCompressionType) -> Result<(), ArrayRetrievalError> {
        if self.compression == compression {
            Ok(())
        } else {
            let bytes = self.encode_bytestring(compression);
            self.data = bytes;
            self.compression = compression;
            Ok(())
        }
    }

    /// Recode the stored data as the requested binary data type.
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
                self.dtype = dtype;
                self.set_buffer_of_type(recast)
            }
            BinaryDataArrayType::Float64 => {
                let view = match self.to_f64() {
                    Ok(view) => view,
                    Err(err) => return Err(err),
                };
                let recast = to_bytes(&view);
                self.dtype = dtype;
                self.set_buffer_of_type(recast)
            }
            BinaryDataArrayType::Int32 => {
                let view = match self.to_i32() {
                    Ok(view) => view,
                    Err(err) => return Err(err),
                };
                let recast = to_bytes(&view);
                self.dtype = dtype;
                self.set_buffer_of_type(recast)
            }
            BinaryDataArrayType::Int64 => {
                let view = match self.to_i64() {
                    Ok(view) => view,
                    Err(err) => return Err(err),
                };
                let recast = to_bytes(&view);
                self.dtype = dtype;
                self.set_buffer_of_type(recast)
            }
            _ => Ok(0),
        };
        result
    }

    /// Test if the the array describes an ion mobility quantity.
    ///
    /// # See also
    /// [`ArrayType::is_ion_mobility`]
    pub const fn is_ion_mobility(&self) -> bool {
        self.name.is_ion_mobility()
    }
}

impl<'transient, 'lifespan: 'transient> ByteArrayView<'transient, 'lifespan> for DataArray {
    fn view(&'lifespan self) -> Result<Cow<'lifespan, [u8]>, ArrayRetrievalError> {
        self.decode()
    }

    fn dtype(&self) -> BinaryDataArrayType {
        self.dtype
    }

    fn data_len(&'lifespan self) -> Result<usize, ArrayRetrievalError> {
        if let Some(z) = self.item_count {
            Ok(z)
        } else {
            let view = self.view()?;
            let n = view.len();
            Ok(n / self.dtype().size_of())
        }
    }

    fn unit(&self) -> Unit {
        self.unit
    }

    fn name(&self) -> &ArrayType {
        &self.name
    }
}

impl<'transient, 'lifespan: 'transient> ByteArrayViewMut<'transient, 'lifespan> for DataArray {
    fn view_mut(&'transient mut self) -> Result<&'transient mut Bytes, ArrayRetrievalError> {
        self.decode_mut()
    }

    fn unit_mut(&mut self) -> &mut Unit {
        &mut self.unit
    }
}

impl_param_described_deferred!(DataArray);


/// Represent a slice of a [`DataArray`] that manages offsets and decoding automatically.
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

    pub const fn is_ion_mobility(&self) -> bool {
        self.source.is_ion_mobility()
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

    fn unit(&self) -> Unit {
        self.source.unit
    }

    fn name(&self) -> &ArrayType {
        self.source.name()
    }
}


#[cfg(test)]
mod test {
    use super::*;
    use std::io;
    use std::fs;

    use super::DataArray;

    fn make_array_from_file() -> io::Result<DataArray> {
        let mut fh = fs::File::open("./test/data/mz_f64_zlib_bas64.txt")?;
        let mut buf = String::new();
        fh.read_to_string(&mut buf)?;
        let bytes: Vec<u8> = buf.into();
        let mut da = DataArray::wrap(&ArrayType::MZArray, BinaryDataArrayType::Float64, bytes);
        da.compression = BinaryCompressionType::Zlib;
        *da.unit_mut() = Unit::MZ;
        assert_eq!(da.unit(), Unit::MZ);
        assert!(!da.is_ion_mobility());
        assert_eq!(da.name(), &ArrayType::MZArray);
        Ok(da)
    }

    #[test]
    fn test_decode() -> io::Result<()> {
        let mut da = make_array_from_file()?;
        da.decode_and_store()?;
        let view = da.to_f64()?;
        assert_eq!(view.len(), 19800);
        Ok(())
    }

    #[test]
    fn test_decode_store() -> io::Result<()> {
        let mut da = make_array_from_file()?;
        da.decode_and_store()?;
        let back = da.clone();
        da.store_as(BinaryDataArrayType::Float32)?;
        let view = da.to_f64()?;
        assert_eq!(view.len(), 19800);
        for (a, b) in back.iter_f64()?.zip(view.iter().copied()) {
            let err= (a as f64 - b).abs();
            assert!((a as f64 - b).abs() < 1e-3, "{} - {} = {}", a, b, err);
        }
        for (a, b) in back.iter_f64()?.zip(da.iter_f32()?.map(|x| x as f64)) {
            let err= (a as f64 - b).abs();
            assert!((a as f64 - b).abs() < 1e-3, "{} - {} = {}", a, b, err);
        }
        da.store_as(BinaryDataArrayType::Float64)?;
        let view = da.to_f64()?;
        assert_eq!(view.len(), 19800);
        for (a, b) in back.iter_f64()?.zip(view.iter().copied()) {
            let err= (a as f64 - b).abs();
            assert!((a as f64 - b).abs() < 1e-3, "{} - {} = {}", a, b, err);
        }
        Ok(())
    }

    #[test]
    fn test_decode_roundtrip() -> io::Result<()> {
        let mut da = make_array_from_file()?;
        let compressed_size = da.data_len()?;
        da.decode_and_store()?;
        let view = da.to_f64()?;
        assert_eq!(view.len(), 19800);
        drop(view);
        da.store_compressed(BinaryCompressionType::Zlib)?;
        assert_eq!(da.compression, BinaryCompressionType::Zlib);
        assert_eq!(da.data_len()?, compressed_size);
        let view = da.to_f64()?;
        assert_eq!(view.len(), 19800);
        Ok(())
    }

    #[test]
    fn test_decode_empty() {
        let mut da = DataArray::wrap(&ArrayType::MZArray, BinaryDataArrayType::Float64, Vec::new());
        da.compression = BinaryCompressionType::Zlib;

        assert_eq!(da.data.len(), 0);
        assert_eq!(da.data_len().unwrap(), 0);
        assert_eq!(da.decode().unwrap().len(), 0);
        assert_eq!(da.to_f64().unwrap().len(), 0);
    }
}