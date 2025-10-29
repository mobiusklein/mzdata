//! Implements a parser for the imzML (Imaging Mass Spectrometry Markup Language) format
//! for representing mass spectrometry imaging data.
//!
//! **Requires the `imzml` feature**
//!
//! imzML is based on mzML but stores binary data in an external `.ibd` file instead of 
//! embedding it as base64 encoded data in the XML. This allows for more efficient storage
//! of large imaging datasets.
//!
//! The format consists of:
//! - `.imzml` file: XML metadata based on mzML schema with imaging-specific CV terms
//! - `.ibd` file: Binary data file containing mass spectra data
//!
//! Data can be stored in two modes:
//! - **Continuous**: All spectra share the same m/z values (more efficient for similar spectra)
//! - **Processed**: Each spectrum has its own m/z and intensity arrays
//!
//! See: <https://www.ms-imaging.org/imzml/>

#![cfg(feature = "imzml")]

pub mod reader;

pub use reader::{ImzMLReaderType, ImzMLReader, is_imzml};

// Re-export UUID for convenience
pub use uuid::Uuid;

#[cfg(test)]
mod tests;
