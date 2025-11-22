//! Reader implementation for Bruker's TDF data files, [`TDFFrameReaderType`] for ion mobility frames
//! and [`TDFSpectrumReaderType`] for summed or sliced spectra.
//!
//! **Requires the `bruker_tdf` feature**
//!
//! Depends upon the [`timsrust`] library, a cross-platform, pure Rust implementation of the Bruker-specifc
//! file reading behaviors and [`rusqlite`] for reading the SQLite3 .tdf files.
mod constants;
mod arrays;
mod sql;
mod reader;

pub use reader::{TDFFrameReader, TDFFrameReaderType, TDFSpectrumReader, TDFSpectrumReaderType, is_tdf};
pub use sql::{ChromatographyData, SQLTrace};