mod infer_format;
pub mod mgf;
pub mod mzml;
#[cfg(feature = "mzmlb")]
pub mod mzmlb;
mod offset_index;
pub mod traits;
mod utils;

pub(crate) mod compression;

pub use crate::io::infer_format::{
    infer_format, infer_from_path, infer_from_stream, open_file, MassSpectrometryFormat,
};
pub use crate::io::mgf::{MGFError, MGFReader, MGFWriter};
#[cfg(feature = "async")]
pub use crate::io::mzml::AsyncMzMLReaderType;
pub use crate::io::mzml::{MzMLParserError, MzMLReader, MzMLWriter};
#[cfg(feature = "mzmlb")]
pub use crate::io::mzmlb::{MzMLbError, MzMLbReader};
pub use crate::io::offset_index::OffsetIndex;
pub use crate::io::traits::{
    RandomAccessSpectrumIterator, ScanAccessError, ScanSource, SpectrumIterator,
};
pub use crate::io::utils::{DetailLevel, PreBufferedStream};
