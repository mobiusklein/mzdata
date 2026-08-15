//! The `DataArray` implementation for data buffer manipulation used for, amongst other things,
//! mapping numeric data to and from optionally compressed little endian base-64 byte arrays with
//! as few data copies as possible, algebraically.
//!
//! Also covers the `ArrayType` enum that covers the standardized array flavors telling us which
//! dimension an array describes, `BinaryArrayMap` collections of (`ArrayType`, `DataArray`), and
//! traits for converting to and from `BinaryArrayMap` based upon the set of available `ArrayType`.
//! This includes a generalization for `BinaryArrayMap3D` when ion mobility is available.

pub use mzdata_bindata::*;