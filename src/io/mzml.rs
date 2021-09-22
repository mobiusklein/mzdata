/*!
Implements a parser for the PSI-MS mzML and indexedmzML XML file formats
for representing raw and processed mass spectra, providing a
[`RandomAccessSpectrumIterator`](crate::io::traits::RandomAccessSpectrumIterator)
interface for reading, and a simple [`ScanWriter`](crate::io::traits::ScanWriter)
interface for writing.

The mzML format is standardized by the Proteomics Standards Initiative (PSI), with
a formal schema defined at <https://www.psidev.info/mzML>.

This crate supports both reading and writing (indexed) mzML documents with spectra
of varying degrees of complexity (raw profiles, centroids, processed centroids), though
extensive customization of the coercion process relies on the [`std::convert::From`] trait
templated on [`&BinaryArrayMap`](crate::spectrum::signal::BinaryArrayMap) for reading
and a complementary reverse conversion for writing.

*/

mod reader;
mod writer;

pub use crate::io::mzml::reader::{
    CVParamParse, FileMetadataBuilder, MzMLParserError, MzMLParserState, MzMLReader,
    MzMLReaderType, ParserResult, SpectrumBuilding, XMLParseBase,
};

pub use crate::io::mzml::writer::{MzMLWriter, MzMLWriterState, MzMLWriterType, WriterResult};
