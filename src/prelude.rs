//! A set of foundational traits used throughout the library.
pub use crate::io::traits::{
    MZFileReader, RandomAccessSpectrumIterator, SpectrumAccessError, ScanSource, ScanWriter, SeekRead,
    SpectrumGrouping, SpectrumIterator, RandomAccessSpectrumGroupingIterator,
};

pub use crate::meta::MSDataFileMetadata;
pub use crate::params::{ParamDescribed, ParamLike};
pub use crate::spectrum::{IonProperties, PrecursorSelection, SpectrumLike};
pub use crate::spectrum::bindata::{ByteArrayView, ByteArrayViewMut, BuildArrayMapFrom, BuildFromArrayMap};

#[cfg(feature = "mzsignal")]
pub use crate::spectrum::group::SpectrumGroupAveraging;

#[doc(hidden)]
pub use std::convert::TryInto;
#[doc(hidden)]
pub use std::io::prelude::*;
