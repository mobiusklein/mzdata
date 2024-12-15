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
//! # Example
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
pub mod io;
pub mod meta;
#[macro_use]
pub mod params;
pub mod prelude;
pub mod spectrum;
pub mod utils;

pub use crate::io::{MZReader, MZReaderBuilder};
pub use crate::io::mgf::{MGFReader, MGFWriter};
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