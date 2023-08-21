pub mod mgf;
pub mod mzml;
pub mod offset_index;
pub mod prelude;
pub mod traits;
pub mod utils;
pub mod infer_format;

pub(crate) mod compression;

pub use crate::io::mgf::MGFReader;
pub use crate::io::mzml::MzMLReader;
pub use crate::io::offset_index::OffsetIndex;
pub use crate::io::traits::{RandomAccessSpectrumIterator, ScanAccessError, SpectrumIterator, ScanSource};
pub use crate::io::infer_format::open_file;
