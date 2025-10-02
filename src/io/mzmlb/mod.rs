//! Implements a parser for the mzMLb file format
//! for representing raw and processed mass spectra, providing a
//! [`RandomAccessSpectrumIterator`](crate::io::traits::RandomAccessSpectrumIterator)
//! interface for reading, and [`SpectrumWriter`](crate::io::traits::SpectrumWriter)
//! interface for writing.
//!
//! **Requires the `mzmlb` feature**
//!
//! The mzMLb format embeds a variant of the mzML format within an HDF5 file, storing
//! the spectrum metadata in XML and data arrays in separate datasets in the same file.
//!
//! This crate supports both reading and writing (indexed) mzML documents with spectra
//! of varying degrees of complexity (raw profiles, centroids, processed centroids), though
//! extensive customization of the coercion process relies on the [`BuildFromArrayMap`](crate::spectrum::bindata::BuildFromArrayMap) and
//! [`BuildArrayMapFrom`](crate::spectrum::bindata::BuildArrayMapFrom) traits
//! for reading and writing conversion to [`BinaryArrayMap`](crate::spectrum::bindata::BinaryArrayMap).

mod reader;
mod common;
mod writer;

pub use reader::{MzMLbReader, MzMLbError, MzMLbReaderType, MzMLbSpectrumBuilder};
pub use writer::{MzMLbWriterType, MzMLbWriterError, MzMLbWriterBuilder, MzMLbWriter};