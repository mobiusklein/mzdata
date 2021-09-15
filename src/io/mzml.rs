//! Implements a parser for the PSI-MS mzML and indexedmzML XML file formats
//! for representing raw and processed mass spectra.

pub mod reader;
pub mod writer;

pub use crate::io::mzml::reader::{
    FileMetadataBuilder, MzMLParserError, MzMLParserState, MzMLReader, MzMLReaderType,
    ParserResult, SpectrumBuilding,
};

pub use crate::io::mzml::writer::{MzMLWriter, MzMLWriterState, MzMLWriterType, XMLResult};
