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
    to_bytes, ArrayRetrievalError, ArrayType, BinaryCompressionType, BinaryDataArrayType,
    Bytes,
};
use super::traits::{ByteArrayView, ByteArrayViewMut};
#[allow(unused)]
use super::vec_as_bytes;

#[allow(unused)]
use super::encodings::{
    reverse_transpose_f32, reverse_transpose_f64, reverse_transpose_i32, reverse_transpose_i64,
    transpose_f32, transpose_f64, transpose_i32, transpose_i64,
};

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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DataArray {
    pub data: Bytes,
    pub dtype: BinaryDataArrayType,
    pub compression: BinaryCompressionType,
    pub name: ArrayType,
    pub params: Option<Box<ParamList>>,
    pub unit: Unit,
    item_count: Option<usize>,
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
            self.data.extend_from_slice(data);
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
            self.data.extend_from_slice(data);
            Ok(())
        }
    }

    pub fn encode_bytestring(&self, compression: BinaryCompressionType) -> Bytes {
        if self.compression == compression {
            log::debug!("Fast-path encoding {}:{}", self.name, self.dtype);
            return self.data.clone();
        }
        let bytestring = match self.compression {
            BinaryCompressionType::Decoded => Cow::Borrowed(self.data.as_slice()),
            _ => self.decode().expect("Failed to decode binary data"),
        };
        match compression {
            BinaryCompressionType::Decoded => panic!("Should never happen"),
            BinaryCompressionType::Zlib => {
                let compressed = Self::compress_zlib(&bytestring);
                base64_simd::STANDARD.encode_type::<Bytes>(&compressed)
            }
            BinaryCompressionType::NoCompression => {
                base64_simd::STANDARD.encode_type::<Bytes>(bytestring.as_ref())
            }
            #[cfg(feature = "numpress")]
            BinaryCompressionType::NumpressLinear => {
                if self.dtype != BinaryDataArrayType::Float64 {
                    panic!("Cannot Numpress non-float64 data!");
                }
                let compressed =
                    Self::compress_numpress_linear(bytemuck::cast_slice(&bytestring)).unwrap();
                base64_simd::STANDARD.encode_type::<Bytes>(&compressed)
            }
            #[cfg(feature = "numpress")]
            BinaryCompressionType::NumpressSLOF => {
                let compressed = match self.dtype {
                    BinaryDataArrayType::Float32 => {
                        Self::compress_numpress_slof(bytemuck::cast_slice::<u8, f32>(&bytestring)).unwrap()
                    },
                    BinaryDataArrayType::Float64 => {
                        Self::compress_numpress_slof(bytemuck::cast_slice::<u8, f64>(&bytestring)).unwrap()
                    },
                    _ => {
                        panic!("Cannot Numpress non-float data!");
                    }
                };
                base64_simd::STANDARD.encode_type::<Bytes>(&compressed)
            }
            #[cfg(feature = "numpress")]
            BinaryCompressionType::NumpressLinearZlib => {
                if self.dtype != BinaryDataArrayType::Float64 {
                    panic!("Cannot Numpress non-float64 data!");
                }
                let compressed = Self::compress_numpress_linear(bytemuck::cast_slice(&bytestring))
                    .inspect_err(|e| {
                        log::error!("Failed to compress buffer with numpress: {e}");
                    })
                    .unwrap();
                let compressed = Self::compress_zlib(&compressed);
                base64_simd::STANDARD.encode_type::<Bytes>(&compressed)
            }
            #[cfg(all(feature = "numpress", feature = "zstd"))]
            BinaryCompressionType::NumpressLinearZstd => {
                if self.dtype != BinaryDataArrayType::Float64 {
                    panic!("Cannot Numpress non-float64 data!");
                }
                let compressed = Self::compress_numpress_linear(bytemuck::cast_slice(&bytestring))
                    .inspect_err(|e| {
                        log::error!("Failed to compress buffer with numpress: {e}");
                    })
                    .unwrap();
                let compressed = Self::compress_zstd(&compressed, BinaryDataArrayType::Unknown, false);
                base64_simd::STANDARD.encode_type::<Bytes>(&compressed)
            }
            #[cfg(feature = "numpress")]
            BinaryCompressionType::NumpressSLOFZlib => {
                let compressed = match self.dtype {
                    BinaryDataArrayType::Float32 => {
                        Self::compress_numpress_slof(bytemuck::cast_slice::<u8, f32>(&bytestring)).unwrap()
                    },
                    BinaryDataArrayType::Float64 => {
                        Self::compress_numpress_slof(bytemuck::cast_slice::<u8, f64>(&bytestring)).unwrap()
                    },
                    _ => {
                        panic!("Cannot Numpress non-float data!");
                    }
                };
                let bytestring = Self::compress_zlib(&compressed);
                base64_simd::STANDARD.encode_type::<Bytes>(&bytestring)
            }
            #[cfg(all(feature = "numpress", feature = "zstd"))]
            BinaryCompressionType::NumpressSLOFZstd => {
                let compressed = match self.dtype {
                    BinaryDataArrayType::Float32 => {
                        Self::compress_numpress_slof(bytemuck::cast_slice::<u8, f32>(&bytestring)).unwrap()
                    },
                    BinaryDataArrayType::Float64 => {
                        Self::compress_numpress_slof(bytemuck::cast_slice::<u8, f64>(&bytestring)).unwrap()
                    },
                    _ => {
                        panic!("Cannot Numpress non-float data!");
                    }
                };
                let bytestring = Self::compress_zstd(&compressed, BinaryDataArrayType::Unknown, false);
                base64_simd::STANDARD.encode_type::<Bytes>(&bytestring)
            }
            #[cfg(feature = "zstd")]
            BinaryCompressionType::Zstd => {
                let compressed = Self::compress_zstd(&bytestring, self.dtype, false);
                base64_simd::STANDARD.encode_type::<Bytes>(&compressed)
            }
            #[cfg(feature = "zstd")]
            BinaryCompressionType::ShuffleZstd => {
                let compressed = Self::compress_zstd(&bytestring, self.dtype, true);
                base64_simd::STANDARD.encode_type::<Bytes>(&compressed)
            }
            #[cfg(feature = "zstd")]
            BinaryCompressionType::DeltaShuffleZstd => {
                let compressed = Self::compress_delta_zstd(&bytestring, self.dtype, true);
                base64_simd::STANDARD.encode_type::<Bytes>(&compressed)
            }
            #[cfg(feature = "zstd")]
            BinaryCompressionType::ZstdDict => {
                let compressed = Self::compress_dict_zstd(&bytestring, self.dtype);
                base64_simd::STANDARD.encode_type::<Bytes>(&compressed)
            }
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

    pub fn decompress_zlib(bytestring: &[u8]) -> Bytes {
        let result = Bytes::new();
        let mut decompressor = ZlibDecoder::new(result);
        decompressor
            .write_all(bytestring)
            .unwrap_or_else(|e| panic!("Decompression error: {}", e));
        decompressor
            .finish()
            .unwrap_or_else(|e| panic!("Decompression error: {}", e))
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
    pub fn compress_numpress_slof<T: numpress::AsFloat64>(data: &[T]) -> Result<Bytes, ArrayRetrievalError> {
        let scaling = numpress::optimal_slof_fixed_point(data);
        let mut buf = Bytes::new();
        match numpress::encode_slof(data, &mut buf, scaling) {
            Ok(_) => Ok(buf),
            Err(e) => Err(ArrayRetrievalError::DecompressionError(e.to_string())),
        }
    }

    #[cfg(feature = "numpress")]
    pub fn decompress_numpress_linear(data: &[u8]) -> Result<Vec<f64>, ArrayRetrievalError> {
        match numpress::numpress_decompress(data) {
            Ok(data) => Ok(data),
            Err(e) => Err(ArrayRetrievalError::DecompressionError(e.to_string())),
        }
    }

    #[cfg(feature = "numpress")]
    pub fn decompress_numpress_slof(data: &[u8], dtype: BinaryDataArrayType) -> Result<Cow<'static, [u8]>, ArrayRetrievalError> {
        let mut buf = Vec::new();

        let decoded = match numpress::decode_slof(data, &mut buf) {
            Ok(_) => buf,
            Err(e) => return Err(ArrayRetrievalError::DecompressionError(e.to_string())),
        };

        match dtype {
            BinaryDataArrayType::Float64 => {
                let view = vec_as_bytes(decoded);
                Ok(Cow::Owned(view))
            },
            BinaryDataArrayType::Float32 => {
                let n = decoded.len() * 4;
                let mut view: Vec<u8> = Vec::with_capacity(n);
                view = decoded.into_iter().map(|v| v as f32).fold(view, |mut view, val| {
                    view.extend_from_slice(bytemuck::bytes_of(&val));
                    view
                });
                Ok(Cow::Owned(view))
            },
            _ => {
                Err(ArrayRetrievalError::DecompressionError(
                    BinaryCompressionType::NumpressSLOF.unsupported_msg(Some(
                        format!("Not compatible with {:?}", dtype).as_str(),
                    )),
                ))
            }
        }
    }

    #[cfg(feature = "zstd")]
    pub(crate) fn compress_zstd(
        bytestring: &[u8],
        dtype: BinaryDataArrayType,
        shuffle: bool,
    ) -> Bytes {
        let level: i32 = std::env::var("MZDATA_ZSTD_LEVEL")
            .map(|v| v.parse())
            .unwrap_or(Ok(zstd::DEFAULT_COMPRESSION_LEVEL))
            .unwrap_or(zstd::DEFAULT_COMPRESSION_LEVEL);
        if !shuffle {
            return zstd::bulk::compress(bytestring, level).unwrap();
        }
        match dtype {
            BinaryDataArrayType::Unknown | BinaryDataArrayType::ASCII => {
                zstd::bulk::compress(bytestring, level).unwrap()
            }
            BinaryDataArrayType::Float64 => {
                zstd::bulk::compress(&transpose_f64(bytemuck::cast_slice(bytestring)), level)
                    .unwrap()
            }
            BinaryDataArrayType::Float32 => {
                zstd::bulk::compress(&transpose_f32(bytemuck::cast_slice(bytestring)), level)
                    .unwrap()
            }
            BinaryDataArrayType::Int64 => {
                zstd::bulk::compress(&transpose_i64(bytemuck::cast_slice(bytestring)), level)
                    .unwrap()
            }
            BinaryDataArrayType::Int32 => {
                zstd::bulk::compress(&transpose_i32(bytemuck::cast_slice(bytestring)), level)
                    .unwrap()
            }
        }
    }

    #[cfg(feature = "zstd")]
    pub(crate) fn compress_delta_zstd(
        bytestring: &[u8],
        dtype: BinaryDataArrayType,
        shuffle: bool,
    ) -> Bytes {
        use bytemuck::cast_slice;

        use super::delta_encoding;

        match dtype {
            BinaryDataArrayType::Unknown | BinaryDataArrayType::ASCII => {
                Self::compress_zstd(bytestring, dtype, shuffle)
            }
            BinaryDataArrayType::Float64 => {
                let mut buf = cast_slice::<_, f64>(bytestring).to_vec();
                delta_encoding(&mut buf);
                Self::compress_zstd(cast_slice(&buf), dtype, shuffle)
            }
            BinaryDataArrayType::Float32 => {
                let mut buf = cast_slice::<_, f32>(bytestring).to_vec();
                delta_encoding(&mut buf);
                Self::compress_zstd(cast_slice(&buf), dtype, shuffle)
            }
            BinaryDataArrayType::Int64 => {
                let mut buf = cast_slice::<_, i64>(bytestring).to_vec();
                delta_encoding(&mut buf);
                Self::compress_zstd(cast_slice(&buf), dtype, shuffle)
            }
            BinaryDataArrayType::Int32 => {
                let mut buf = cast_slice::<_, i32>(bytestring).to_vec();
                delta_encoding(&mut buf);
                Self::compress_zstd(cast_slice(&buf), dtype, shuffle)
            }
        }
    }

    #[cfg(feature = "zstd")]
    pub fn compress_dict_zstd(bytestring: &[u8], dtype: BinaryDataArrayType) -> Bytes {
        use super::encodings::dictionary_encoding;
        log::trace!("Dictionary encoding {} bytes as {dtype}", bytestring.len());
        if bytestring.is_empty() {
            return Self::compress_zstd(&bytestring, dtype, false);
        }
        match dtype {
            BinaryDataArrayType::Float64 => {
                let compressed =
                    dictionary_encoding(bytemuck::cast_slice::<u8, f64>(&bytestring))
                        .unwrap();
                let compressed = Self::compress_zstd(&compressed, dtype, false);
                compressed
            }
            BinaryDataArrayType::Float32 => {
                let compressed =
                    dictionary_encoding(bytemuck::cast_slice::<u8, f32>(&bytestring))
                        .unwrap();
                let compressed = Self::compress_zstd(&compressed, dtype, false);
                compressed
            }
            BinaryDataArrayType::Int64 => {
                let compressed =
                    dictionary_encoding(bytemuck::cast_slice::<u8, i64>(&bytestring))
                        .unwrap();
                let compressed = Self::compress_zstd(&compressed, dtype, false);
                compressed
            }
            BinaryDataArrayType::Int32 => {
                let compressed =
                    dictionary_encoding(bytemuck::cast_slice::<u8, i32>(&bytestring))
                        .unwrap();
                let compressed = Self::compress_zstd(&compressed, dtype, false);
                compressed
            }
            _ => {
                let compressed =
                    dictionary_encoding(bytemuck::cast_slice::<u8, u8>(&bytestring))
                        .unwrap();
                let compressed = Self::compress_zstd(&compressed, dtype, false);
                compressed
            }
        }
    }

    #[cfg(feature = "zstd")]
    pub(crate) fn decompress_zstd(data: &[u8], dtype: BinaryDataArrayType, shuffle: bool) -> Bytes {
        let mut decoder = zstd::Decoder::new(std::io::Cursor::new(data)).unwrap();
        let mut buf = Vec::new();
        decoder.read_to_end(&mut buf).unwrap();
        if !shuffle {
            return buf;
        }
        match dtype {
            BinaryDataArrayType::Unknown | BinaryDataArrayType::ASCII => buf,
            BinaryDataArrayType::Float64 => reverse_transpose_f64(&buf),
            BinaryDataArrayType::Float32 => reverse_transpose_f32(&buf),
            BinaryDataArrayType::Int64 => reverse_transpose_i64(&buf),
            BinaryDataArrayType::Int32 => reverse_transpose_i32(&buf),
        }
    }

    #[cfg(feature = "zstd")]
    pub(crate) fn decompress_delta_zstd(
        data: &[u8],
        dtype: BinaryDataArrayType,
        shuffle: bool,
    ) -> Bytes {
        use super::delta_decoding;

        let mut delta = Self::decompress_zstd(data, dtype, shuffle);
        match dtype {
            BinaryDataArrayType::Unknown | BinaryDataArrayType::ASCII => delta,
            BinaryDataArrayType::Float64 => {
                let buf = bytemuck::cast_slice_mut::<_, f64>(&mut delta);
                delta_decoding(buf);
                delta
            }
            BinaryDataArrayType::Float32 => {
                let buf = bytemuck::cast_slice_mut::<_, f32>(&mut delta);
                delta_decoding(buf);
                delta
            }
            BinaryDataArrayType::Int64 => {
                let buf = bytemuck::cast_slice_mut::<_, i64>(&mut delta);
                delta_decoding(buf);
                delta
            }
            BinaryDataArrayType::Int32 => {
                let buf = bytemuck::cast_slice_mut::<_, i32>(&mut delta);
                delta_decoding(buf);
                delta
            }
        }
    }

    #[cfg(feature = "zstd")]
    pub(crate) fn decompress_dict_zstd(bytestring: &[u8], dtype: BinaryDataArrayType) -> Bytes {
        use super::encodings::dictionary_decoding;

        let data = Self::decompress_zstd(bytestring, dtype, false);
        match dtype {
            BinaryDataArrayType::ASCII | BinaryDataArrayType::Unknown => dictionary_decoding(&data).unwrap(),
            BinaryDataArrayType::Float64 => {
                to_bytes(&dictionary_decoding::<f64>(&data).unwrap())
            },
            BinaryDataArrayType::Float32 => {
                to_bytes(&dictionary_decoding::<f32>(&data).unwrap())
            },
            BinaryDataArrayType::Int64 => {
                to_bytes(&dictionary_decoding::<i64>(&data).unwrap())
            },
            BinaryDataArrayType::Int32 => {
                to_bytes(&dictionary_decoding::<i32>(&data).unwrap())
            },
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
            return Ok(Cow::Borrowed(&EMPTY_BUFFER));
        }

        macro_rules! base64_decode {
            () => {
                base64_simd::STANDARD
                    .decode_type::<Bytes>(&self.data)
                    .unwrap_or_else(|e| panic!("Failed to decode base64 array: {}", e))
            };
        }

        match self.compression {
            BinaryCompressionType::Decoded => Ok(Cow::Borrowed(self.data.as_slice())),
            BinaryCompressionType::NoCompression => {
                let bytestring = base64_decode!();
                Ok(Cow::Owned(bytestring))
            }
            BinaryCompressionType::Zlib => {
                let bytestring = base64_decode!();
                Ok(Cow::Owned(Self::decompress_zlib(&bytestring)))
            }
            #[cfg(feature = "zstd")]
            BinaryCompressionType::Zstd => {
                let bytestring = base64_decode!();
                Ok(Cow::Owned(Self::decompress_zstd(
                    &bytestring,
                    self.dtype,
                    false,
                )))
            }
            #[cfg(feature = "zstd")]
            BinaryCompressionType::ShuffleZstd => {
                let bytestring = base64_decode!();
                Ok(Cow::Owned(Self::decompress_zstd(
                    &bytestring,
                    self.dtype,
                    true,
                )))
            }
            #[cfg(feature = "zstd")]
            BinaryCompressionType::DeltaShuffleZstd => {
                let bytestring = base64_decode!();
                Ok(Cow::Owned(Self::decompress_delta_zstd(
                    &bytestring,
                    self.dtype,
                    true,
                )))
            }

            #[cfg(feature = "zstd")]
            BinaryCompressionType::ZstdDict => {
                let bytestring = base64_decode!();
                Ok(Cow::Owned(Self::decompress_dict_zstd(
                    &bytestring,
                    self.dtype,
                )))
            }
            #[cfg(feature = "numpress")]
            BinaryCompressionType::NumpressLinear => match self.dtype {
                BinaryDataArrayType::Float64 => {
                    let bytestring = base64_decode!();
                    let decoded = Self::decompress_numpress_linear(&bytestring)?;
                    let view = vec_as_bytes(decoded);
                    Ok(Cow::Owned(view))
                }
                _ => Err(ArrayRetrievalError::DecompressionError(
                    self.compression.unsupported_msg(Some(
                        format!("Not compatible with {:?}", self.dtype).as_str(),
                    )),
                )),
            },

            #[cfg(feature = "numpress")]
            BinaryCompressionType::NumpressSLOF => {
                let bytestring = base64_decode!();
                Self::decompress_numpress_slof(&bytestring, self.dtype)
            }

            #[cfg(feature = "numpress")]
            BinaryCompressionType::NumpressLinearZlib => match self.dtype {
                BinaryDataArrayType::Float64 => {
                    let bytestring = base64_decode!();
                    let bytestring = Self::decompress_zlib(bytemuck::cast_slice(&bytestring));
                    let decoded = Self::decompress_numpress_linear(&bytestring)?;
                    let view = vec_as_bytes(decoded);
                    Ok(Cow::Owned(view))
                }
                _ => Err(ArrayRetrievalError::DecompressionError(
                    self.compression.unsupported_msg(Some(
                        format!("Not compatible with {:?}", self.dtype).as_str(),
                    )),
                )),
            },

            #[cfg(feature = "numpress")]
            BinaryCompressionType::NumpressSLOFZlib => {
                let bytestring = base64_decode!();
                let bytestring = Self::decompress_zlib(bytemuck::cast_slice(&bytestring));
                Self::decompress_numpress_slof(&bytestring, self.dtype)
            }

            #[cfg(all(feature = "numpress", feature = "zstd"))]
            BinaryCompressionType::NumpressLinearZstd => match self.dtype {
                BinaryDataArrayType::Float64 => {
                    let bytestring = base64_decode!();
                    let bytestring = Self::decompress_zstd(
                        bytemuck::cast_slice(&bytestring),
                        BinaryDataArrayType::Unknown,
                        false,
                    );
                    let decoded = Self::decompress_numpress_linear(&bytestring)?;
                    let view = vec_as_bytes(decoded);
                    Ok(Cow::Owned(view))
                }
                _ => Err(ArrayRetrievalError::DecompressionError(
                    self.compression.unsupported_msg(Some(
                        format!("Not compatible with {:?}", self.dtype).as_str(),
                    )),
                )),
            },

            #[cfg(all(feature = "numpress", feature = "zstd"))]
            BinaryCompressionType::NumpressSLOFZstd => {
                let bytestring = base64_decode!();
                let bytestring = Self::decompress_zstd(bytemuck::cast_slice(&bytestring), BinaryDataArrayType::Unknown, false);
                Self::decompress_numpress_slof(&bytestring, self.dtype)
            }

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
            _ => {
                Ok(Cow::Owned(self.decode()?[start..end].to_vec()))
            }
        }
    }

    pub fn decode_mut(&'transient mut self) -> Result<&'transient mut Bytes, ArrayRetrievalError> {
        if self.data.is_empty() {
            return Ok(&mut self.data);
        }

        if matches!(self.compression, BinaryCompressionType::Decoded) {
            return Ok(&mut self.data)
        }

        match self.decode()? {
            Cow::Borrowed(_) => {
                return Ok(&mut self.data)
            },
            Cow::Owned(owned) => {
                self.data = owned;
                self.compression = BinaryCompressionType::Decoded;
                return Ok(&mut self.data)
            },
        }
    }

    pub fn clear(&mut self) {
        self.data.clear();
        self.params = None;
        self.item_count = None;
    }

    /// The reverse of [`DataArray::decode_and_store`], this method compresses `self.data` to the desired
    /// compression method and stores that buffer as `self.data`.
    pub fn store_compressed(
        &mut self,
        compression: BinaryCompressionType,
    ) -> Result<(), ArrayRetrievalError> {
        if self.compression == compression {
            Ok(())
        } else {
            self.item_count = self.data_len().ok();
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
        match dtype {
            BinaryDataArrayType::Float32 => {
                let view = self.to_f32()?;
                let recast = to_bytes(&view);
                self.dtype = dtype;
                self.set_buffer_of_type(recast)
            }
            BinaryDataArrayType::Float64 => {
                let view = self.to_f64()?;
                let recast = to_bytes(&view);
                self.dtype = dtype;
                self.set_buffer_of_type(recast)
            }
            BinaryDataArrayType::Int32 => {
                let view = self.to_i32()?;
                let recast = to_bytes(&view);
                self.dtype = dtype;
                self.set_buffer_of_type(recast)
            }
            BinaryDataArrayType::Int64 => {
                let view = self.to_i64()?;
                let recast = to_bytes(&view);
                self.dtype = dtype;
                self.set_buffer_of_type(recast)
            }
            _ => Ok(0),
        }
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
    use std::fs;
    use std::io;

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

    fn make_array_from_file_im_zstd() -> io::Result<DataArray> {
        let mut fh = fs::File::open("test/data/im_f64_zstd_base64.txt")?;
        let mut buf = String::new();
        fh.read_to_string(&mut buf)?;
        let bytes: Vec<u8> = buf.into();
        let mut da = DataArray::wrap(&ArrayType::MeanInverseReducedIonMobilityArray, BinaryDataArrayType::Float64, bytes);
        da.compression = BinaryCompressionType::Zstd;
        *da.unit_mut() = Unit::VoltSecondPerSquareCentimeter;
        assert_eq!(da.unit(), Unit::VoltSecondPerSquareCentimeter);
        assert!(da.is_ion_mobility());
        assert_eq!(da.name(), &ArrayType::MeanInverseReducedIonMobilityArray);
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
            let err = (a - b).abs();
            assert!((a - b).abs() < 1e-3, "{} - {} = {}", a, b, err);
        }
        for (a, b) in back.iter_f64()?.zip(da.iter_f32()?.map(|x| x as f64)) {
            let err = (a - b).abs();
            assert!((a - b).abs() < 1e-3, "{} - {} = {}", a, b, err);
        }
        da.store_as(BinaryDataArrayType::Float64)?;
        let view = da.to_f64()?;
        assert_eq!(view.len(), 19800);
        for (a, b) in back.iter_f64()?.zip(view.iter().copied()) {
            let err = (a - b).abs();
            assert!((a - b).abs() < 1e-3, "{} - {} = {}", a, b, err);
        }
        Ok(())
    }

    #[cfg(feature = "zstd")]
    #[test]
    fn test_decode_delta_zstd() {
        let points: Vec<f64> = (0..200_000usize).map(|i| i as f64 * 0.01 + 100.0).collect();
        let mut da = DataArray::from_name(&ArrayType::MZArray);
        da.extend(&points).unwrap();

        let decoded_len = da.data.len();
        da.store_compressed(BinaryCompressionType::ShuffleZstd)
            .unwrap();
        let zstd_len = da.data.len();

        da.store_compressed(BinaryCompressionType::Zlib).unwrap();
        let zlib_len = da.data.len();

        da.decode_and_store().unwrap();

        da.store_compressed(BinaryCompressionType::DeltaShuffleZstd)
            .unwrap();
        let delta_zstd_len = da.data.len();
        eprintln!("decoded: {decoded_len};\nzlib: {zlib_len};\nzstd: {zstd_len};\ndelta-zstd: {delta_zstd_len}");
        da.decode_and_store().unwrap();
        let view = da.to_f64().unwrap();
        let err: f64 = points
            .iter()
            .zip(view.iter())
            .map(|(a, b)| {
                assert!((a - b).abs() < 1e-3, "{a} - {b} = {}", a - b);
                (a - b).abs()
            })
            .sum();
        let mean_err = err / (points.len() as f64);
        eprintln!("mean abs error: {mean_err:0.8}")
    }

    #[cfg(feature = "zstd")]
    #[test]
    fn test_decode_zstd() -> io::Result<()> {
        let mut da = make_array_from_file()?;
        let zlib_len = da.data.len();
        da.decode_and_store()?;

        let decoded_len = da.data.len();

        da.store_compressed(BinaryCompressionType::ShuffleZstd)?;

        let zstd_len = da.data.len();

        eprintln!("zlib: {zlib_len};\ndecoded: {decoded_len};\nzstd: {zstd_len}");
        da.decode_and_store()?;

        let mut da_ref = make_array_from_file()?;
        da_ref.decode_and_store()?;
        assert_eq!(da.data, da_ref.data);
        Ok(())
    }

    #[cfg(feature = "numpress")]
    #[test]
    fn test_numpress_linear() -> io::Result<()> {
        let mut da = make_array_from_file()?;
        let zlib_len = da.data.len();
        da.decode_and_store()?;

        let decoded_len = da.data.len();

        da.store_compressed(BinaryCompressionType::NumpressLinear)?;
        let numpress_len = da.data.len();

        eprintln!("zlib: {zlib_len};\ndecoded: {decoded_len};\nnumpress: {numpress_len}");
        da.decode_and_store()?;

        let mut da_ref = make_array_from_file()?;
        da_ref.decode_and_store()?;
        for (a, b) in da.iter_f64()?.zip(da_ref.iter_f64()?) {
            assert!((a - b).abs() < 1e-3, "{a} - {b} = {} which is too large a deviation", (a - b).abs())
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
        let mut da = DataArray::wrap(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            Vec::new(),
        );
        da.compression = BinaryCompressionType::Zlib;

        assert_eq!(da.data.len(), 0);
        assert_eq!(da.data_len().unwrap(), 0);
        assert_eq!(da.decode().unwrap().len(), 0);
        assert_eq!(da.to_f64().unwrap().len(), 0);
    }


    #[cfg(feature = "zstd")]
    #[test]
    fn test_dict_from_base64() -> io::Result<()> {
        let mut da = make_array_from_file_im_zstd()?;

        da.decode_and_store()?;
        assert_eq!(da.data_len()?, 221);

        da.store_compressed(BinaryCompressionType::ZstdDict)?;
        assert_eq!(da.data_len()?, 221);

        Ok(())
    }
}
