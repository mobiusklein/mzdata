use std::iter::{FusedIterator, ExactSizeIterator};

use crate::spectrum::Chromatogram;


/// A trait that for retrieving [`Chromatogram`]s from a source.
pub trait ChromatogramSource {
    /// Get a [`Chromatogram`] by its identifier, if it exists.
    fn get_chromatogram_by_id(&mut self, id: &str) -> Option<Chromatogram>;

    /// Get a [`Chromatogram`] by its index, if it exists.
    fn get_chromatogram_by_index(&mut self, index: usize) -> Option<Chromatogram>;

    /// Iterate over [`Chromatogram`]s with a [`ChromatogramIterator`]
    fn iter_chromatograms(&mut self) -> ChromatogramIterator<'_, Self>
    where
        Self: Sized,
    {
        ChromatogramIterator::new(self)
    }

    /// Get the number of pre-defined chromatograms available in this source
    fn count_chromatograms(&self) -> usize;
}

/// A facade for a [`ChromatogramSource`] that is an [`Iterator`] over [`Chromatogram`] instances
/// using [`ChromatogramSource::get_chromatogram_by_index`]
#[derive(Debug)]
pub struct ChromatogramIterator<'a, R: ChromatogramSource> {
    source: &'a mut R,
    index: usize,
    done: bool,
}

impl<'a, R: ChromatogramSource> ChromatogramIterator<'a, R> {
    pub fn new(source: &'a mut R) -> Self {
        Self {
            source,
            index: 0,
            done: false,
        }
    }
}

impl<R: ChromatogramSource> Iterator for ChromatogramIterator<'_, R> {
    type Item = Chromatogram;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done || self.index >= self.source.count_chromatograms() {
            self.done = true;
            return None;
        }
        if let Some(chrom) = self.source.get_chromatogram_by_index(self.index) {
            self.index += 1;
            Some(chrom)
        } else {
            self.done = true;
            None
        }
    }
}

impl<R: ChromatogramSource> FusedIterator for ChromatogramIterator<'_, R> {}

impl<R: ChromatogramSource> ExactSizeIterator for ChromatogramIterator<'_, R> {
    fn len(&self) -> usize {
        if self.done {
            0
        } else {
            self.source.count_chromatograms().saturating_sub(self.index)
        }
    }
}

#[cfg(feature = "async_partial")]
mod async_impl {
    use super::*;

    use futures::{
        stream,
        Stream,
    };

    pub trait AsyncChromatogramSource: Send {
        /// Get a [`Chromatogram`] by its identifier, if it exists.
        fn get_chromatogram_by_id(&mut self, id: &str) -> impl std::future::Future<Output = Option<Chromatogram>>;

        /// Get a [`Chromatogram`] by its index, if it exists.
        fn get_chromatogram_by_index(&mut self, index: usize) -> impl std::future::Future<Output = Option<Chromatogram>>;

        /// Wrap this source in a [`Stream`] over its chromatograms
        ///
        /// The returned stream is [`Unpin`], so it can be driven directly with
        /// [`StreamExt::next`](futures::StreamExt::next) without pinning it first.
        fn as_stream(&mut self) -> impl Stream<Item = Chromatogram> + Unpin + '_ {
            let it = 0..self.count_chromatograms();
            Box::pin(stream::unfold((self, it), |(reader, mut rng)| async {
                let i = rng.next()?;

                let spec = reader.get_chromatogram_by_index(i);
                spec.await.map(|val| (val, (reader, rng)))
            }))
        }

        /// Get the number of pre-defined chromatograms available in this source
        fn count_chromatograms(&self) -> usize;
    }
}

#[cfg(feature = "async_partial")]
pub use async_impl::AsyncChromatogramSource;
