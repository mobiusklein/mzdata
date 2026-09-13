//! The `DataArray` implementation for data buffer manipulation used for, amongst other things,
//! mapping numeric data to and from optionally compressed little endian base-64 byte arrays with
//! as few data copies as possible, algebraically.
//!
//! Also covers the `ArrayType` enum that covers the standardized array flavors telling us which
//! dimension an array describes, `BinaryArrayMap` collections of (`ArrayType`, `DataArray`), and
//! traits for converting to and from `BinaryArrayMap` based upon the set of available `ArrayType`.
//! This includes a generalization for `BinaryArrayMap3D` when ion mobility is available.
//!
//! ## Compression
//!
//! `mzdata` supports several kinds of data compression methods with tradeoffs
//!
//! ### Lossless methods
//! - `zlib` compression [`BinaryCompressionType::Zlib`]: A well known, widely available general-purpose compressor codec with many implementations,
//!     with varying speed and size tradeoffs taken from [`flate2`].
//! - `zstd` compression [`BinaryCompressionType::Zstd`]: A newer general-purpose codec that is not widely used in mzML but very popular in other applications.
//!     Much faster and often better compression compared to `zlib`.
//!     Requires feature `zstd`
//! - `byte-shuffled zstd` compression [`BinaryCompressionType::ShuffleZstd`]: This first applies a byte shuffle transform to the data to compress
//!     the `n`th byte of each value contiguously, and then apply `zstd`. This makes sorted data compress more effectively.
//!     Requires feature `zstd`
//! - `dictionary encoded zstd` [`BinaryCompressionType::ZstdDict`]: Another elaboration of `byte-shuffled zstd` where the values are first dictionary-encoded
//!     so that repeated values are denoted by indices rather than the full size value itself. The dictionary is sorted and then both the dictionary and the
//!     indices are byte shuffle encoded and then `zstd` compressed. This is even more effective than previous codecs for data with repeated, semi-sorted values
//!     like m/z and ion mobility.
//!     Requires feature `zstd`
//!
//! ### Lossy methods
//! - `numpress-linear` compression [`BinaryCompressionType::NumpressLinear`]: A combination simple linear model fit and byte packing codec that stores the model
//!     coefficients as float64 and the residuals as unsigned 16 or 32-bit integers. With an optimal fixed point, this will produce a median loss of accuracy below
//!     of about 2e-8 absolute units, but may vary widely from dataset to dataset. This codec **MUST** be used with sorted data only, and is only suitable for
//!     float64 data. It may be combined with `zlib` [`BinaryCompressionType::NumpressLinearZlib`] or `zstd` [`BinaryCompressionType::NumpressLinearZstd`] for even
//!     greater space savings. See details at <https://doi.org/10.1074/mcp.O114.037879>. /// Requires feature `numpress`
//! - `numpress-slof` compression [`BinaryCompressionType::NumpressSLOF`]: A lossy numerical approximation that log transforms and scales values to make them very
//!     close to 16-bit integers. This codec is lossier the larger the value is, and is only appropriate for count data like intensities. It can also combine with
//!     `zlib` [`BinaryCompressionType::NumpressSLOFZlib`] and `zstd` [`BinaryCompressionType::NumpressSLOFZstd`]. See details at <https://doi.org/10.1074/mcp.O114.037879>
//!     Requires feature `numpress`
//! - `numpress-pic` compression [`BinaryCompressionType::NumpressPIC`]: This simply rounds a floating point value and then packs it
//!     into a truncated integer. May be suitable for intensity values, but not consistently better than [`BinaryCompressionType::NumpressSLOF`].
//!     See details at <https://doi.org/10.1074/mcp.O114.037879>
//!     Requires feature `numpress`
//! - Linear prediction compression [`BinaryCompressionType::LinearPrediction`]: The same concept as [`BinaryCompressionType::NumpressLinear`] without the truncation
//!     All of the same concerns apply though it loses less precision. This does not perform any form of compression. It was uses as part of mzMLb which applies
//!     container-level compression. Only suitable for float64 data.
//! - Delta encoding [`BinaryCompressionType::DeltaPrediction`]: This computes first order difference of a sorted array of values and stores the difference
//!     and the starting point. It is sometimes more accurate than [`BinaryCompressionType::LinearPrediction`]. Like [`BinaryCompressionType::LinearPrediction`],
//!     it does not do any actual compression. It was meant for use with mzMLb, which compresses the entire container. Only suitable for float64 data.
//! - Delta-shuffled `zstd` [`BinaryCompressionType::DeltaShuffleZstd`]: Combines [`BinaryCompressionType::DeltaPrediction`] with [`BinaryCompressionType::ShuffleZstd`].
//!     Requires feature `zstd`
//!
//!
//! ## Supported Data types
//!
//! `mzdata` supports five basic binary data types via the [`BinaryDataArrayType`] enum. Data that doesn't follow any of these types
//! may still be transported with the [`BinaryDataArrayType::Unknown`] variant, but it will not be directly usable.
//!
//! These values may be decoded on the fly or stored for faster repeated access. Coercion with [`DataArray::coerce`] relies on
//! [`bytemuck`] to do the heavy lifting. Conversion through [`DataArray::convert`] maps between different types but makes copies.
//!
//!
//!

mod encodings;
mod array;
mod conversion;
mod map;
mod traits;

pub mod utils;


pub use encodings::{
    ArrayRetrievalError, ArrayType, BinaryCompressionType, BinaryDataArrayType, Bytes, as_bytes,
    delta_decoding, delta_encoding, linear_prediction_decoding, linear_prediction_encoding,
    to_bytes, vec_as_bytes,
};

pub use array::{DataArray, DataArraySlice};
pub use map::{BinaryArrayMap, BinaryArrayMap3D};
pub use traits::{ByteArrayView, ByteArrayViewMut};

pub use conversion::{
    ArraysAvailable, BuildArrayMap3DFrom, BuildArrayMapFrom, BuildFromArrayMap, BuildFromArrayMap3D,
};