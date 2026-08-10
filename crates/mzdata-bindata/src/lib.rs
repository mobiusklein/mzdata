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