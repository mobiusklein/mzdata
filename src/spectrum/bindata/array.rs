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

    pub fn update_buffer<T: Pod>(
        &mut self,
        data_buffer: &[T],
    ) -> Result<usize, ArrayRetrievalError> {
        if self.dtype.size_of() != mem::size_of::<T>() {
            Err(ArrayRetrievalError::DataTypeSizeMismatch)
        } else {
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
            let data = bytemuck::cast_slice(values);
            self.data.extend(data.iter());
            Ok(())
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
            .expect("Error decompression");
        decompressor.finish().expect("Error decompression")
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
        if self.data.is_empty() {
            return Ok(Cow::Borrowed(&EMPTY_BUFFER))
        }
        match self.compression {
            BinaryCompressionType::Decoded => Ok(Cow::Borrowed(self.data.as_slice())),
            BinaryCompressionType::NoCompression => {
                let bytestring = base64_simd::STANDARD.decode_type::<Bytes>(&self.data)
                    .expect("Failed to decode base64 array");
                Ok(Cow::Owned(bytestring))
            }
            BinaryCompressionType::Zlib => {
                let bytestring = base64_simd::STANDARD.decode_type::<Bytes>(&self.data)
                    .expect("Failed to decode base64 array");
                Ok(Cow::Owned(Self::decompres_zlib(&bytestring)))
            }
            #[cfg(feature = "numpress")]
            BinaryCompressionType::NumpressLinear => match self.dtype {
                BinaryDataArrayType::Float64 => {
                    let mut bytestring = base64_simd::STANDARD.decode_type::<Bytes>(&self.data)
                        .expect("Failed to decode base64 array");
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
                    .expect("Failed to decode base64 array");
                Ok(Cow::Owned(bytestring[start..end].to_vec()))
            }
            BinaryCompressionType::Zlib => {
                let bytestring = base64_simd::STANDARD.decode_type::<Bytes>(&self.data)
                    .expect("Failed to decode base64 array");
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
                    .expect("Failed to decode base64 array");
                self.data = bytestring;
                self.compression = BinaryCompressionType::Decoded;
                Ok(&mut self.data)
            }
            BinaryCompressionType::Zlib => {
                let bytestring = base64_simd::STANDARD.decode_type::<Bytes>(&self.data)
                    .expect("Failed to decode base64 array");
                self.data = bytestring;
                self.compression = BinaryCompressionType::Decoded;
                Ok(&mut self.data)
            },
            #[cfg(feature = "numpress")]
            BinaryCompressionType::NumpressLinear => match self.dtype {
                BinaryDataArrayType::Float64 => {
                    let mut bytestring = base64_simd::STANDARD.decode_type::<Bytes>(&self.data)
                        .expect("Failed to decode base64 array");
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