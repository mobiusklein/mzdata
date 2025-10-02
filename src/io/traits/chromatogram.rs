use std::iter::FusedIterator;

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
}

/// A facade for a [`ChromatogramSource`] that is an [`Iterator`] over [`Chromatogram`] instances
/// using [`ChromatogramSource::get_chromatogram_by_index`]
#[derive(Debug)]
pub struct ChromatogramIterator<'a, R: ChromatogramSource> {
    source: &'a mut R,
    index: usize,
}

impl<'a, R: ChromatogramSource> ChromatogramIterator<'a, R> {
    pub fn new(source: &'a mut R) -> Self {
        Self { source, index: 0 }
    }
}

impl<R: ChromatogramSource> Iterator for ChromatogramIterator<'_, R> {
    type Item = Chromatogram;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(chrom) = self.source.get_chromatogram_by_index(self.index) {
            self.index += 1;
            Some(chrom)
        } else {
            None
        }
    }
}

impl<R: ChromatogramSource> FusedIterator for ChromatogramIterator<'_, R> {}