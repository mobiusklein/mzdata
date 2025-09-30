//! A set of foundational traits used throughout the library.
pub use crate::io::traits::{
    ChromatogramSource, IntoIonMobilityFrameSource, IonMobilityFrameSource, IonMobilityFrameWriter,
    MZFileReader, RandomAccessIonMobilityFrameGroupingIterator,
    RandomAccessIonMobilityFrameIterator, RandomAccessSpectrumGroupingIterator,
    RandomAccessSpectrumIterator, RandomAccessSpectrumSource as _, SeekRead, SpectrumAccessError,
    SpectrumSource, SpectrumSourceWithMetadata as _, SpectrumWriter,
};

#[cfg(feature = "async_partial")]
pub use crate::io::traits::AsyncSpectrumSource;

pub use crate::meta::MSDataFileMetadata;
pub use crate::params::{ParamDescribed, ParamDescribedRead, ParamLike, ParamValue};
pub use crate::spectrum::bindata::{
    BuildArrayMap3DFrom, BuildArrayMapFrom, BuildFromArrayMap, BuildFromArrayMap3D, ByteArrayView,
    ByteArrayViewMut,
};
pub use crate::spectrum::{
    ChromatogramLike, IonMobilityFrameGrouping, IonMobilityFrameLike, IonMobilityMeasure,
    IonProperties, PrecursorSelection, SpectrumGrouping, SpectrumLike,
};

#[cfg(feature = "mzsignal")]
pub use crate::spectrum::group::SpectrumGroupAveraging;

#[doc(hidden)]
pub use mzpeaks::prelude::*;
#[doc(hidden)]
pub use std::convert::TryInto;
#[doc(hidden)]
pub use std::io::prelude::*;
