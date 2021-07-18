//! mzdata provides basic access to raw and processed mass spectrometry data formats in
//! Rust.
//!
//! The library currently supports reading:
//!   1. MGF files using [`MGFReader`]
//!   2. mzML files using [`MzMLReader`]
//!
//! It also provides a sorted data structure for representing peak lists, [`PeakSet`]
//! and a trait implementing the majority of the logic, [`PeakCollection`].
//!
//!
pub mod io;
pub mod mass_error;
pub mod peaks;
pub mod spectrum;
pub mod meta;
pub mod params;

pub use crate::peaks::coordinate::{CoordinateDimension, CoordinateLike, Mass, MZ};

pub use crate::mass_error::MassErrorType;
pub use crate::peaks::peak::{CentroidPeak, DeconvolutedPeak};
pub use crate::peaks::{PeakCollection, PeakSet};

pub use crate::io::mgf::MGFReader;
pub use crate::io::mzml::MzMLReader;

pub use crate::spectrum::{CentroidSpectrum, RawSpectrum, SpectrumBehavior};
