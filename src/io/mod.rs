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
#[cfg(feature = "proxi")]
pub mod proxi;
mod shorthand;
pub(crate) mod traits;
mod utils;

pub(crate) mod compression;

pub use crate::io::infer_format::{
    infer_format, infer_from_path, infer_from_stream, IMMZReaderType, MZReader, MZReaderBuilder,
    MZReaderType, MassSpectrometryFormat, MassSpectrometryReadWriteProcess, Sink, Source,
};

#[cfg(feature = "mgf")]
pub use crate::io::mgf::{MGFError, MGFReader, MGFWriter};

#[cfg(all(feature = "async", feature = "mzml"))]
pub use crate::io::mzml::AsyncMzMLReader;

#[cfg(feature = "mzml")]
pub use crate::io::mzml::{MzMLParserError, MzMLReader, MzMLWriter};

#[cfg(feature = "mzmlb")]
pub use crate::io::mzmlb::{MzMLbError, MzMLbReader};

pub use crate::io::offset_index::OffsetIndex;
pub use crate::io::traits::{
    BorrowedGeneric3DIonMobilityFrameSource, ChromatogramIterator, ChromatogramSource,
    Generic3DIonMobilityFrameSource, IntoIonMobilityFrameSource, IntoIonMobilityFrameSourceError,
    IonMobilityFrameAccessError, IonMobilityFrameGrouping, IonMobilityFrameIterator,
    IonMobilityFrameSource, MZFileReader, MemorySpectrumSource,
    RandomAccessIonMobilityFrameGroupingIterator, RandomAccessIonMobilityFrameIterator,
    RandomAccessSpectrumGroupingIterator, RandomAccessSpectrumIterator, RandomAccessSpectrumSource,
    SpectrumAccessError, SpectrumGrouping, SpectrumIterator, SpectrumReceiver, SpectrumSource,
    SpectrumSourceWithMetadata, SpectrumWriter, StreamingSpectrumIterator,
};

#[cfg(feature = "async_partial")]
pub use crate::io::traits::AsyncSpectrumSource;

#[cfg(feature = "async")]
pub use crate::io::traits::AsyncMZFileReader;

pub use crate::io::utils::{DetailLevel, PreBufferedStream};

#[cfg(feature = "checksum")]
pub use crate::io::utils::checksum_file;

pub use compression::RestartableGzDecoder;

#[cfg(any(feature = "thermo", feature = "doc-only"))]
pub mod thermo;

#[cfg(any(feature = "thermo", feature = "doc-only"))]
pub use thermo::ThermoRawReader;

#[cfg(all(
    feature = "async_partial",
    any(feature = "thermo", feature = "doc-only")
))]
pub use thermo::AsyncThermoRawReader;

#[cfg(feature = "async_partial")]
pub use infer_format::{AsyncMZReader, AsyncMZReaderBuilder, AsyncMZReaderType};

#[cfg(feature = "bruker_tdf")]
pub mod tdf;

pub mod usi;
