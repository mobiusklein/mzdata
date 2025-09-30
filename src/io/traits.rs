mod chromatogram;
mod frame;
mod spectrum;
mod util;

pub use spectrum::{
    MZFileReader, MemorySpectrumSource, RandomAccessSpectrumGroupingIterator,
    RandomAccessSpectrumIterator, RandomAccessSpectrumSource, SpectrumAccessError,
    SpectrumIterator, SpectrumReceiver, SpectrumSource, SpectrumSourceWithMetadata, SpectrumWriter,
    StreamingSpectrumIterator,
};
pub use util::SeekRead;

pub use frame::{
    BorrowedGeneric3DIonMobilityFrameSource, Generic3DIonMobilityFrameSource,
    IntoIonMobilityFrameSource, IntoIonMobilityFrameSourceError, IonMobilityFrameAccessError,
    IonMobilityFrameIterator, IonMobilityFrameSource, IonMobilityFrameWriter,
    RandomAccessIonMobilityFrameGroupingIterator, RandomAccessIonMobilityFrameIterator,
};

pub use chromatogram::{ChromatogramIterator, ChromatogramSource};

pub use crate::spectrum::group::{IonMobilityFrameGrouping, SpectrumGrouping};

#[cfg(feature = "async_partial")]
pub use spectrum::{AsyncRandomAccessSpectrumIterator, AsyncSpectrumSource, SpectrumStream};

#[cfg(feature = "async")]
pub use spectrum::AsyncMZFileReader;

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_object_safe() {
        // If `SpectrumSource` were not object safe, this code
        // couldn't compile.
        let _f = |_x: &dyn SpectrumSource| {};
    }
}
