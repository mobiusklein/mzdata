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
#[cfg(any(feature = "mgf", feature = "mzml"))]
pub(crate) use util::DetailLevelGuard;
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

    #[cfg(feature = "mzml")]
    #[test]
    fn test_failed_detail_checks() -> std::io::Result<()> {
        use crate::io::{DetailLevel, OffsetIndex};
        for level in [
            DetailLevel::Full,
            DetailLevel::Lazy,
            DetailLevel::MetadataOnly,
        ] {
            let mut reader = crate::MzMLReader::open_path("test/data/small.mzML")?;
            reader.set_detail_level(level);
            let mut index = OffsetIndex::new("spectrum".into());
            index.insert("invalid offset", u64::MAX);
            index.init = true;
            reader.set_index(index);
            assert!(reader.get_spectrum_by_time(1.0).is_none());
            assert_eq!(*reader.detail_level(), level);
            assert!(reader.has_ion_mobility().is_none());
            assert_eq!(*reader.detail_level(), level);
            reader.set_index(OffsetIndex::new("spectrum".into()));
            assert!(reader.get_spectrum_by_time(1.0).is_none());
            assert_eq!(
                reader.has_ion_mobility(),
                Some(crate::spectrum::HasIonMobility::None)
            );
            assert_eq!(*reader.detail_level(), level);
        }
        Ok(())
    }

    #[test]
    fn test_object_safe() {
        // If `SpectrumSource` were not object safe, this code
        // couldn't compile.
        let _f = |_x: &dyn SpectrumSource| {};
    }
}
