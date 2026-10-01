mod chromatogram;
mod frame;
mod spectrum;
mod util;

pub use spectrum::{
    MZFileReader, MemorySpectrumSource, RandomAccessSpectrumGroupingIterator,
    RandomAccessSpectrumIterator, RandomAccessSpectrumSource, SpectrumAccessError,
    SpectrumIterator, SpectrumReceiver, SpectrumSource,
    SpectrumSourceWithMetadata, SpectrumWriter, StreamingSpectrumIterator,
};
pub use util::SeekRead;

pub use frame::{
    BorrowedGeneric3DIonMobilityFrameSource, Generic3DIonMobilityFrameSource,
    IonMobilityFrameAccessError, IonMobilityFrameIterator,
    IonMobilityFrameSource, IonMobilityFrameWriter, RandomAccessIonMobilityFrameIterator,
    RandomAccessIonMobilityFrameGroupingIterator,
    IntoIonMobilityFrameSourceError,
    IntoIonMobilityFrameSource
};

pub use chromatogram::{ChromatogramIterator, ChromatogramSource};

pub use crate::spectrum::group::{SpectrumGrouping, IonMobilityFrameGrouping};

#[cfg(feature = "async_partial")]
pub use spectrum::{AsyncSpectrumSource, AsyncRandomAccessSpectrumIterator, SpectrumStream};

#[cfg(feature = "async_partial")]
pub use chromatogram::AsyncChromatogramSource;

#[cfg(feature = "async_partial")]
pub use frame::{
    AsyncGeneric3DIonMobilityFrameSource, AsyncIntoIonMobilityFrameSource,
    AsyncIonMobilityFrameSource, AsyncRandomAccessIonMobilityFrameIterator, IonMobilityFrameStream,
};

#[cfg(feature = "mzsignal")]
pub use spectrum::PeakPicking;

#[cfg(feature = "async")]
pub use spectrum::AsyncMZFileReader;

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_empty_spectrum_time_seek() {
        let mut source: MemorySpectrumSource = MemorySpectrumSource::default();
        let mut iter = source.iter();
        for time in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(matches!(
                iter.start_from_time(time),
                Err(SpectrumAccessError::SpectrumNotFound)
            ));
            assert!(iter.next().is_none());
            assert!(iter.next_back().is_none());
        }
        let mut scan = crate::spectrum::MultiLayerSpectrum::default();
        scan.description.id = "first".into();
        let mut source: MemorySpectrumSource = std::collections::VecDeque::from([scan]).into();
        let mut iter = source.iter();
        assert_eq!(iter.next().unwrap().description.id, "first");
        assert!(matches!(
            iter.start_from_time(f64::NAN),
            Err(SpectrumAccessError::IOError(None))
        ));
        assert!(iter.next().is_none());
        iter.start_from_time(0.0).unwrap();
        assert_eq!(iter.next().unwrap().description.id, "first");
    }

    #[test]
    fn test_empty_frame_time_seek() {
        let source: MemorySpectrumSource = MemorySpectrumSource::default();
        let mut source: Generic3DIonMobilityFrameSource<_, _, _> =
            Generic3DIonMobilityFrameSource::new(source);
        let mut iter = source.iter();
        for time in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(matches!(
                iter.start_from_time(time),
                Err(IonMobilityFrameAccessError::FrameNotFound)
            ));
            assert!(iter.next().is_none());
            assert!(iter.next_back().is_none());
        }
        let scan = crate::spectrum::MultiLayerSpectrum::default();
        let source: MemorySpectrumSource = std::collections::VecDeque::from([scan]).into();
        let mut source: Generic3DIonMobilityFrameSource<_, _, _> =
            Generic3DIonMobilityFrameSource::new(source);
        let mut iter = source.iter();
        assert!(matches!(
            iter.start_from_time(0.0),
            Err(IonMobilityFrameAccessError::IOError(None))
        ));
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_object_safe() {
        // If `SpectrumSource` were not object safe, this code
        // couldn't compile.
        let _f = |_x: &dyn SpectrumSource| {};
    }
}
