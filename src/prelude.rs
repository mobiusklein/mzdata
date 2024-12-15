//! A set of foundational traits used throughout the library.
pub use crate::io::traits::{
    MZFileReader, RandomAccessSpectrumGroupingIterator, RandomAccessSpectrumIterator,
    RandomAccessSpectrumSource as _, SpectrumSourceWithMetadata as _, SpectrumSource,
    SpectrumWriter, SeekRead, SpectrumAccessError, IonMobilityFrameSource,
    RandomAccessIonMobilityFrameIterator, ChromatogramSource,
    IonMobilityFrameWriter, RandomAccessIonMobilityFrameGroupingIterator,
};

pub use crate::meta::MSDataFileMetadata;
pub use crate::params::{ParamDescribed, ParamLike, ParamValue, ParamDescribedRead};
pub use crate::spectrum::bindata::{
    BuildArrayMapFrom, BuildFromArrayMap, ByteArrayView, ByteArrayViewMut,
};
pub use crate::spectrum::{
    IonProperties, PrecursorSelection, SpectrumLike, IonMobilityMeasure, IonMobilityFrameLike,
    SpectrumGrouping, IonMobilityFrameGrouping,
};

#[cfg(feature = "mzsignal")]
pub use crate::spectrum::group::SpectrumGroupAveraging;

#[doc(hidden)]
pub use std::convert::TryInto;
#[doc(hidden)]
pub use std::io::prelude::*;
#[doc(hidden)]
pub use mzpeaks::prelude::*;