//! A set of foundational traits used throughout the library.
pub use crate::io::traits::{
    MZFileReader, RandomAccessSpectrumGroupingIterator, RandomAccessSpectrumIterator, ScanSource,
    ScanWriter, SeekRead, SpectrumAccessError, SpectrumGrouping, SpectrumIterator,
};

pub use crate::meta::MSDataFileMetadata;
pub use crate::params::{ParamDescribed, ParamLike};
pub use crate::spectrum::bindata::{
    BuildArrayMapFrom, BuildFromArrayMap, ByteArrayView, ByteArrayViewMut,
};
pub use crate::spectrum::{IonProperties, PrecursorSelection, SpectrumLike};

#[cfg(feature = "mzsignal")]
pub use crate::spectrum::group::SpectrumGroupAveraging;

#[doc(hidden)]
pub use std::convert::TryInto;
#[doc(hidden)]
pub use std::io::prelude::*;
