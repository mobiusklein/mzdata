/*!
Implements a parser for the PSI-MS mzML and indexedmzML XML file formats
for representing raw and processed mass spectra, providing a
[`RandomAccessSpectrumIterator`](crate::io::traits::RandomAccessSpectrumIterator)
interface for reading, and [`SpectrumWriter`](crate::io::traits::SpectrumWriter)
interface for writing.

The mzML format is standardized by the Proteomics Standards Initiative (PSI), with
a formal schema defined at <https://www.psidev.info/mzML>.

This crate supports both reading and writing (indexed) mzML documents with spectra
of varying degrees of complexity (raw profiles, centroids, processed centroids), though
extensive customization of the coercion process relies on the [`BuildFromArrayMap`](crate::spectrum::bindata::BuildFromArrayMap) and
[`BuildArrayMapFrom`](crate::spectrum::bindata::BuildArrayMapFrom) traits
for reading and writing conversion to [`BinaryArrayMap`](crate::spectrum::bindata::BinaryArrayMap).

*/

mod reader;
mod reading_shared;
mod writer;

#[cfg(feature = "async")]
mod async_reader;

pub use reading_shared::{
    CVParamParse, MzMLParserError, MzMLParserState, MzMLSAX, XMLParseBase,
    FileMetadataBuilder, EntryType
};

#[allow(unused)]
pub(crate) use reading_shared::{IncrementingIdMap, ParserResult};

pub use crate::io::mzml::reader::{
    MzMLReader, MzMLReaderType, MzMLSpectrumBuilder,
    SpectrumBuilding,
};

pub(crate) use crate::io::mzml::reader::is_mzml;

pub use crate::io::mzml::writer::{MzMLWriter, MzMLWriterState, MzMLWriterType, MzMLWriterError};

#[cfg(feature = "async")]
pub use crate::io::mzml::async_reader::{
    MzMLReader as AsyncMzMLReader, MzMLReaderType as AsyncMzMLReaderType,
};
