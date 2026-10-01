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
    fn test_object_safe() {
        // If `SpectrumSource` were not object safe, this code
        // couldn't compile.
        let _f = |_x: &dyn SpectrumSource| {};
    }

    fn spectra(n: usize) -> MemorySpectrumSource {
        use crate::spectrum::{
            bindata::{ArrayType, BinaryArrayMap, BinaryDataArrayType, DataArray},
            MultiLayerSpectrum,
        };
        MemorySpectrumSource::new(
            (0..n)
                .map(|index| {
                    let mut spectrum = MultiLayerSpectrum::default();
                    spectrum.description.index = index;
                    spectrum.description.id = index.to_string();
                    let mut arrays = BinaryArrayMap::new();
                    arrays.add(DataArray::wrap(
                        &ArrayType::MeanInverseReducedIonMobilityArray,
                        BinaryDataArrayType::Float64,
                        Vec::new(),
                    ));
                    spectrum.arrays = Some(arrays);
                    spectrum
                })
                .collect(),
        )
    }

    fn check_remaining(
        mut it: impl ExactSizeIterator<Item = usize> + DoubleEndedIterator,
        n: usize,
    ) {
        let mut expected = 0..n;
        while !expected.is_empty() {
            assert_eq!(it.len(), expected.len());
            assert_eq!(it.next(), expected.next());
            assert_eq!(it.len(), expected.len());
            assert_eq!(it.next_back(), expected.next_back());
        }
        assert_eq!(it.len(), 0);
        assert_eq!(it.next(), None);
        assert_eq!(it.next_back(), None);
    }

    #[test]
    fn test_spectrum_iterator_remaining() {
        use crate::spectrum::SpectrumLike;
        for n in 0..6 {
            let mut source = spectra(n);
            let mut it = source.iter();
            check_remaining(it.by_ref().map(|s| s.index()), n);
            assert_eq!(SpectrumSource::len(&it), n);
            it.reset();
            assert_eq!(ExactSizeIterator::len(&it), n);
            assert!(it.nth(usize::MAX).is_none());
            assert!(it.next_back().is_none());
            if n > 0 {
                it.start_from_index(n - 1).unwrap();
                assert_eq!(ExactSizeIterator::len(&it), 1);
                assert_eq!(it.next().unwrap().index(), n - 1);
                it.start_from_id("0").unwrap();
                assert_eq!(it.next().unwrap().index(), 0);
                assert!(it.nth(usize::MAX).is_none());
            }
            assert!(it.start_from_index(n).is_err());
        }
    }

    #[test]
    fn test_frame_iterator_remaining() {
        use crate::spectrum::IonMobilityFrameLike;
        use mzpeaks::{CentroidPeak, DeconvolutedPeak};
        for n in 0..6 {
            let mut source =
                Generic3DIonMobilityFrameSource::<CentroidPeak, DeconvolutedPeak, _>::new(spectra(
                    n,
                ));
            let mut it = source.iter();
            check_remaining(it.by_ref().map(|s| s.index()), n);
            assert_eq!(IonMobilityFrameSource::len(&it), n);
            it.reset();
            assert_eq!(ExactSizeIterator::len(&it), n);
            assert!(it.nth(usize::MAX).is_none());
            assert!(it.next_back().is_none());
            if n > 0 {
                it.start_from_index(n - 1).unwrap();
                assert_eq!(ExactSizeIterator::len(&it), 1);
                assert_eq!(it.next().unwrap().index(), n - 1);
                it.start_from_id("0").unwrap();
                assert_eq!(it.next().unwrap().index(), 0);
                assert!(it.nth(usize::MAX).is_none());
            }
            assert!(it.start_from_index(n).is_err());
        }
    }

    #[test]
    fn test_chromatogram_iterator_remaining() {
        use crate::spectrum::Chromatogram;
        struct Chromatograms {
            count: usize,
            fail: bool,
            calls: usize,
        }
        impl ChromatogramSource for Chromatograms {
            fn get_chromatogram_by_id(&mut self, _: &str) -> Option<Chromatogram> {
                None
            }
            fn get_chromatogram_by_index(&mut self, index: usize) -> Option<Chromatogram> {
                self.calls += 1;
                let fail = std::mem::take(&mut self.fail);
                (index < self.count && !fail).then(Chromatogram::default)
            }
            fn count_chromatograms(&self) -> usize {
                self.count
            }
        }
        for (count, fail) in [(0, false), (1, false), (3, false), (3, true)] {
            let mut source = Chromatograms {
                count,
                fail,
                calls: 0,
            };
            let mut it = source.iter_chromatograms().fuse();
            for remaining in (1..=count).rev() {
                assert_eq!(it.len(), remaining);
                if fail {
                    break;
                }
                assert!(it.next().is_some());
            }
            assert!(it.next().is_none());
            assert_eq!(it.len(), 0);
            assert!(it.next().is_none());
            assert_eq!(source.calls, if fail { 1 } else { count });
        }
    }
}
