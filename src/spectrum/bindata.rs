mod array;
mod conversion;
mod encodings;
mod map;
mod traits;

pub use array::{DataArray, DataArraySlice};
pub use encodings::{
    as_bytes, delta_decoding, delta_encoding, linear_prediction_decoding,
    linear_prediction_encoding, to_bytes, vec_as_bytes, ArrayRetrievalError, ArrayType,
    BinaryCompressionType, BinaryDataArrayType, Bytes,
};
pub use conversion::{BuildArrayMapFrom, BuildFromArrayMap, ArraysAvailable};
pub use map::BinaryArrayMap;
pub use traits::{ByteArrayView, ByteArrayViewMut};
