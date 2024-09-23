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
mod shorthand;
pub(crate) mod traits;
mod utils;

pub(crate) mod compression;

pub use crate::io::infer_format::{
    infer_format, infer_from_path, infer_from_stream, MZReader, MZReaderType,
    MassSpectrometryFormat, MassSpectrometryReadWriteProcess, Sink, Source,
};
pub use crate::io::mgf::{MGFError, MGFReader, MGFWriter};
#[cfg(feature = "async")]
pub use crate::io::mzml::AsyncMzMLReader;
pub use crate::io::mzml::{MzMLParserError, MzMLReader, MzMLWriter};
#[cfg(feature = "mzmlb")]
pub use crate::io::mzmlb::{MzMLbError, MzMLbReader};
pub use crate::io::offset_index::OffsetIndex;
pub use crate::io::traits::{
    BorrowedGeneric3DIonMobilityFrameSource, ChromatogramIterator, ChromatogramSource,
    Generic3DIonMobilityFrameSource, IonMobilityFrameAccessError, IonMobilityFrameGrouping,
    IonMobilityFrameIterator, IonMobilityFrameSource, MZFileReader, MemorySpectrumSource,
    RandomAccessIonMobilityFrameIterator, RandomAccessSpectrumGroupingIterator,
    RandomAccessSpectrumIterator, RandomAccessSpectrumSource, SpectrumAccessError,
    SpectrumGrouping, SpectrumIterator, SpectrumReceiver, SpectrumSource,
    SpectrumSourceWithMetadata, SpectrumWriter, StreamingSpectrumIterator,
};
pub use crate::io::utils::{checksum_file, DetailLevel, PreBufferedStream};
pub use compression::RestartableGzDecoder;

#[cfg(feature = "thermo")]
pub mod thermo;
#[cfg(feature = "thermo")]
pub use thermo::ThermoRawReader;

pub mod usi;
pub mod proxi;
