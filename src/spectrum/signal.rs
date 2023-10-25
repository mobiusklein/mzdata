mod bindata;
mod conversion;
mod encodings;
mod map;
mod traits;

pub use bindata::{DataArray, DataArraySlice};
pub use encodings::{
    as_bytes, delta_decoding, delta_encoding, linear_prediction_decoding,
    linear_prediction_encoding, to_bytes, vec_as_bytes, ArrayRetrievalError, ArrayType,
    BinaryCompressionType, BinaryDataArrayType, Bytes,
};
pub use map::BinaryArrayMap;
pub use traits::{ByteArrayView, ByteArrayViewMut};
