//! The `DataArray` implementation for data buffer manipulation used for, amongst other things,
//! mapping numeric data to and from optionally compressed little endian base-64 byte arrays with
//! as few data copies as possible, algebraically.
//!
//! Also covers the `ArrayType` enum that covers the standardized array flavors telling us which
//! dimension an array describes, `BinaryArrayMap` collections of (`ArrayType`, `DataArray`), and
//! traits for converting to and from `BinaryArrayMap` based upon the set of available `ArrayType`.
//! This includes a generalization for `BinaryArrayMap3D` when ion mobility is available.

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