//! Reading and writing mass spectrometry data file formats and abstractions over them.
//!
//! There are many data file formats for recording mass spectrometry data.
//!

mod infer_format;
pub mod mgf;
pub mod mzml;
#[cfg(feature = "mzmlb")]
pub mod mzmlb;
mod offset_index;
pub(crate) mod traits;
mod utils;
mod shorthand;

pub(crate) mod compression;

pub use crate::io::infer_format::{
    infer_format, infer_from_path, infer_from_stream, MassSpectrometryFormat,
    MassSpectrometryReadWriteProcess, Sink, Source, MZReader, MZReaderType
};
pub use crate::io::mgf::{MGFError, MGFReader, MGFWriter};
#[cfg(feature = "async")]
pub use crate::io::mzml::AsyncMzMLReader;
pub use crate::io::mzml::{MzMLParserError, MzMLReader, MzMLWriter};
#[cfg(feature = "mzmlb")]
pub use crate::io::mzmlb::{MzMLbError, MzMLbReader};
pub use crate::io::offset_index::OffsetIndex;
pub use crate::io::traits::{
    MZFileReader, RandomAccessSpectrumIterator, SpectrumSource, SpectrumWriter, SpectrumAccessError,
    SpectrumGrouping, SpectrumIterator, StreamingSpectrumIterator, SpectrumReceiver, RandomAccessSpectrumSource,
    SpectrumSourceWithMetadata, MemorySpectrumSource, RandomAccessSpectrumGroupingIterator,
};
pub use crate::io::utils::{DetailLevel, PreBufferedStream, checksum_file};
pub use compression::RestartableGzDecoder;

#[cfg(feature = "thermorawfilereader")]
pub mod thermo;
#[cfg(feature = "thermorawfilereader")]
pub use thermo::ThermoRawReader;