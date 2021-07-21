//! The data structures and components that represent a mass spectrum and how
//! we access their data.
//!
//! [`crate::spectrum::spectrum`] contains high level structures for individual spectra,
//! [`crate::spectrum::scan_properties`] contains structures for describing different facets
//! of mass spectrum acquisition and other such metadata, and [`crate::spectrum::signal`]
//! includes structures for dealing with raw signal data arrays that may or may not
//! be byte-encoded but not strongly typed, though it does not include signal processing
//! as that is outside the scope of this crate.

pub mod scan_properties;
pub mod signal;
pub mod spectrum;

pub use crate::spectrum::scan_properties::*;
pub use crate::spectrum::signal::{ArrayType, BinaryArrayMap, BinaryDataArrayType, DataArray};
pub use crate::spectrum::spectrum::{CentroidSpectrum, RawSpectrum, SpectrumBehavior, Spectrum};
