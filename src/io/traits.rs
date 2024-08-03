mod chromatogram;
mod frame;
mod spectrum;
mod util;

pub use spectrum::{
    MZFileReader, MemorySpectrumSource, RandomAccessSpectrumGroupingIterator,
    RandomAccessSpectrumIterator, RandomAccessSpectrumSource, SpectrumAccessError,
    SpectrumGrouping, SpectrumIterator, SpectrumReceiver, SpectrumSource,
    SpectrumSourceWithMetadata, SpectrumWriter, StreamingSpectrumIterator,
};
pub use util::SeekRead;

pub use frame::{
    BorrowedGeneric3DIonMobilityFrameSource, Generic3DIonMobilityFrameSource,
    IonMobilityFrameAccessError, IonMobilityFrameGrouping, IonMobilityFrameIterator,
    IonMobilityFrameSource, IonMobilityFrameWriter, RandomAccessIonMobilityFrameIterator,
};

pub use chromatogram::{ChromatogramIterator, ChromatogramSource};

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
