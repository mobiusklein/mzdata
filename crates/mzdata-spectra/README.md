# mzdata-spectra

[![Latest Version](https://img.shields.io/crates/v/mzdata?style=for-the-badge&color=mediumpurple&logo=rust)](https://crates.io/crates/mzdata-spectra)
[![docs.rs](https://img.shields.io/docsrs/mzdata?style=for-the-badge&logo=docs.rs&color=mediumseagreen)](https://docs.rs/mzdata/latest/mzdata/)


This re-exports a subset of the [`mzdata`](https://github.com/mobiusklein/mzdata), specifically
by not enabling the default features that include the mzML and MGF reading and writing components.

This cuts the minimum number of dependencies considerably, but retains the metadata mappings and
spectrum data model, as well as all the traits and "helper" types under `mzdata::io`. For more detail,
please see `mzdata`'s source code and documentation.