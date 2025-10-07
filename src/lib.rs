//! `mzdata` provides basic access to raw and processed mass spectrometry data formats in
//! Rust.
//!
//! For a guide, see the [tutorial] section.
//!
//! The library currently supports reading:
//!   1. MGF files using [`MGFReader`] in [`mzdata::io::mgf`](crate::io::mgf)
//!   2. mzML & indexedmzML files using [`MzMLReader`] in [`mzdata::io::mzml`](crate::io::mzml)
//!   3. mzMLb files using [`MzMLbReader`] in [`mzdata::io::mzmlb`](crate::io::mzmlb), if the `mzmlb` feature is enabled
//!   4. Thermo RAW files using [`ThermoRawReader`](crate::io::thermo::ThermoRawReader) in [`mzdata::io::thermo`](crate::io::thermo), if the `thermo` feature is enabled
//!   5. Bruker TDF files using [`TDFSpectrumReader`](crate::io::tdf::TDFSpectrumReader) in [`mzdata::io::tdf`](crate::io::tdf), if the `bruker_tdf` feature is enabled
//!
//! and writing:
//!   1. MGF files using [`MGFWriter`] in [`mzdata::io::mgf`](crate::io::mgf)
//!   2. mzML & indexedmzML files using [`MzMLWriter`] in [`mzdata::io::mzml`](crate::io::mzml)
//!   3. mzMLb files using [`MzMLbWriter`] in [`mzdata::io::mzmlb`](crate::io::mzmlb), if the `mzmlb` feature is enabled
//!
//! This menagerie of different formats and gzip compression or not can be inferred from a path or [`io::Read`](std::io::Read) using [`io::infer_format`] and [`io::infer_from_stream`].
//! Conventional dispatch is possible through [`MZReader`]. The [`mz_read`] macro provides a convenient means of working with
//! a value with zero added overhead, but with a limited scope. The [`mz_write`] macro is the equivalent for opening a writer.
//! There are additional tools for dealing with file format dispatch in [`MassSpectrometryReadWriteProcess`](crate::io::MassSpectrometryReadWriteProcess).
//!
//! It also includes a set of representation layers for spectra in [`mzdata::spectrum`](crate::spectrum)
//!
//! ## Example
//! ```rust
//! use std::fs;
//! use mzdata::prelude::*;
//! use mzpeaks::Tolerance;
//! use mzdata::MZReader;
//! use mzdata::spectrum::SignalContinuity;
//!
//! let reader = MZReader::open_path("./test/data/small.mzML").unwrap();
//! for spectrum in reader {
//!     println!("Scan {} => BP {}", spectrum.id(), spectrum.peaks().base_peak().mz);
//!
//!     if spectrum.signal_continuity() == SignalContinuity::Centroid {
//!         let peak_picked = spectrum.into_centroid().unwrap();
//!         println!("Matches for 579.155: {:?}",
//!                  peak_picked.peaks.all_peaks_for(
//!                     579.155, Tolerance::Da(0.02)
//!                 )
//!         );
//!     }
//! }
//! ```
//!
//! It uses [`mzpeaks`] to represent peaks and peak lists, and re-exports the basic types. While the high-level
//! types are templated on simple peak types, more complex, application-specific peak types can be substituted.
//! See [`mzdata::spectrum::bindata`](crate::spectrum::bindata) for more information about how to directly convert
//! data arrays to peak lists.
//!
//!
//! ## Traits
//! The library makes heavy use of traits to abstract over the implementation details of different file formats.
//! These traits are included in [`mzdata::prelude`](crate::prelude). It also imports [`mzpeaks::prelude`].
//!
//!
//! ## Features
//! [`mzdata`](crate) provides many optional features, some of which are self-contained, while others layer funcitonality.
//!
//! ### `mzsignal`-related
//!
//! TLDR: Unless you are already using `ndarray-linalg` in your dependency graph, you should enable `mzsignal` + `nalgebra`.
//!
//! The [`mzsignal`](mzsignal) crate provides signal processing, peak picking, and feature finding funcitonality. Part of this behavior requires a linear algebra implementation. [`mzsignal`](mzsignal) is flexible. It can use either `nalgebra`, a pure Rust library that is self-contained but optimized for small matrices, or `ndarray-linalg` which requires an external LAPACK library be available either at build time or run time, all of which are outside the basic Rust ecosystem. Enabling the `mzsignal` feature **requires** one of the following features:
//!
//! - `nalgebra` - No external dependencies.
//! - `openblas` - Requires OpenBlas (see https://crates.io/crates/ndarray-linalg)
//! - `intel-mkl` - Requires Intel's Math Kernel Library (see https://crates.io/crates/ndarray-linalg)
//! - `netlib` - Requires the NETLIB (see https://crates.io/crates/ndarray-linalg)
//!
//! ### File Formats
//!
//! [`mzdata`](crate) supports reading several file formats, some of which add large dependencies and can be opted into or out of.
//!
//! |   Feature    | File Format                                                     | Dependency                                                                                                    |
//! | :----------: | :-------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------ |
//! |   `mzmlb`    | [mzMLb](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00192) | HDF5 C shared library at runtime or statically linked with `hdf5-rs`, possibly a C compiler                   |
//! |   `thermo`   | Thermo-Fisher RAW Format                                        | .NET runtime at build time _and_ runtime, possibly a C compiler                                               |
//! | `bruker_tdf` | Bruker TDF Format                                               | SQLite3 C library at runtime or statically linked with `rusqlite`, requires `mzsignal` for flattening spectra |
//!
//! Additionally, mzML and MGF are supported by default, but they can be disabled by skipping default features and not enabling the `mzml` and `mgf` features.
//!
//! To complicate matters the `hdf5_static` feature combined with `mzmlb` handles statically linking the HDF5 C library and `zlib` together to avoid symbol collision with other compression libraries used by `mzdata`.
//!
//! ### Compression
//!
//! `mzdata` uses [`flate2`](https://crates.io/crates/flate2) to compress and decompress `zlib`-type compressed streams, but there are three different backends available with different tradeoffs in speed and build convenience:
//!
//! - `zlib` - The historical implementation. Faster than `miniz_oxide` and consistently produces the largest compression. Requires a nearly ubiquitous C library at build time.
//! - `zlib-ng-compat` - The fasted, often nearly largest if not largest compression and decompression. Requires a C library or a C compiler at build time.
//! - `zlib-ng` - C library dependency, I encountered build errors but your mileage may vary. Requires a C library or a C compiler at build time.
//! - `miniz_oxide` - Pure Rust backend, the slowest in practice.
//!
//! `mzdata` was also a test-bed for some experimental compression techniques.
//!
//! - `zstd` - Enables layered Zstandard and byte shuffling + dictionary encoding methods.
//!
//! ### Async I/O
//!
//! `mzdata` uses synchronous I/O by default, but includes code for some async options:
//!
//! - `async_partial` - Implements trait-level asynchronous versions of the spectrum reading traits and implementations for mzML, MGF, and Thermo RAW files using [`tokio`](https://crates.io/crates/tokio), but doesn't enable the `tokio/fs` module which carries additional requirements which is not compatible with all platforms.
//! - `async` - Enables `async_partial` and `tokio/fs`.
//!
//! ### PROXI
//!
//! `mzdata` includes [PROXI](https://www.psidev.info/proxi) clients for fetching spectra from supporting servers on the internet using [USIs](https://www.psidev.info/usi).
//!
//! - `proxi` - Provides a synchronous client in `mzdata::io::proxi` and adds `mzdata::io::usi::USI::download_spectrum_blocking`
//! - `async-proxi` - Provides an asynchronous client in `mzdata::io::proxi` and adds `mzdata::io::usi::USI::download_spectrum_async`
//!
//! ### Other
//!
//! - `serde` - Enables `serde` serialization and deserialization for most library types that aren't directly connected to an I/O device.
//! - `parallelism` - Enables `rayon` parallel iterators on a small number of internal operations to speed up some operations relating to decompression signal processing. This is unlikely to be notice-able in most cases. More benefit is had by simply processing multiple spectra in parallel using `rayon`'s bridging adapters.
//!
pub mod io;
pub mod meta;
#[macro_use]
pub mod params;
pub mod prelude;
pub mod spectrum;
pub mod utils;

pub use crate::io::{MZReader, MZReaderBuilder};
#[cfg(feature = "mgf")]
pub use crate::io::mgf::{MGFReader, MGFWriter};
#[cfg(feature = "mzml")]
pub use crate::io::mzml::{MzMLReader, MzMLWriter};

#[cfg(feature = "mzmlb")]
pub use crate::io::mzmlb::{
    MzMLbReader, MzMLbWriter, MzMLbWriterBuilder,
};

#[cfg(feature = "thermo")]
pub use crate::io::thermo::ThermoRawReader;

pub use crate::params::{Param, ParamList};

pub use crate::spectrum::{CentroidSpectrum, RawSpectrum, Spectrum};

#[cfg(doc)]
pub mod tutorial;

pub use mzpeaks;

#[cfg(feature = "mzsignal")]
pub use mzsignal;