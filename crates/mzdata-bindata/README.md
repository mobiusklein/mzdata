# mzdata-bindata

[![Latest Version](https://img.shields.io/crates/v/mzdata?style=for-the-badge&color=mediumpurple&logo=rust)](https://crates.io/crates/mzdata)
[![docs.rs](https://img.shields.io/docsrs/mzdata?style=for-the-badge&logo=docs.rs&color=mediumseagreen)](https://docs.rs/mzdata/latest/mzdata/)

This is a component of [`mzdata`](https://crates.io/crates/mzdata).

The `DataArray` implementation for data buffer manipulation used for, amongst other things,
mapping numeric data to and from optionally compressed little endian base-64 byte arrays with
as few data copies as possible, algebraically.

Also covers the `ArrayType` enum that covers the standardized array flavors telling us which
dimension an array describes, `BinaryArrayMap` collections of (`ArrayType`, `DataArray`), and
traits for converting to and from `BinaryArrayMap` based upon the set of available `ArrayType`.
This includes a generalization for `BinaryArrayMap3D` when ion mobility is available.
