use std::convert::TryFrom;
use std::io;
use std::marker::PhantomData;

use log::warn;
use thiserror::Error;

use mzpeaks::{
    feature::{ChargedFeature, Feature, FeatureLike},
    CentroidLike, DeconvolutedCentroidLike, IonMobility, KnownCharge, Mass, MZ,
};

use crate::spectrum::group::IonMobilityFrameGroupingIterator;
use crate::spectrum::spectrum_types::MultiLayerSpectrum;
use crate::spectrum::{IonMobilityFrameLike, MultiLayerIonMobilityFrame};
use crate::{
    io::{DetailLevel, OffsetIndex},
    prelude::{MSDataFileMetadata, SpectrumLike},
    spectrum::{HasIonMobility, IonMobilityFrameGroup},
};

use super::{
    IonMobilityFrameGrouping, RandomAccessSpectrumIterator, SpectrumAccessError, SpectrumSource,
};

/// An analog of [`SpectrumSource`] for [`IonMobilityFrameLike`] producing types
pub trait IonMobilityFrameSource<
    C: FeatureLike<MZ, IonMobility> = Feature<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge = ChargedFeature<Mass, IonMobility>,
    S: IonMobilityFrameLike<C, D> = MultiLayerIonMobilityFrame<C, D>,
>: Iterator<Item = S>
{
    /// Rewind the current position of the source to the beginning
    fn reset(&mut self);

    /// Get the [`DetailLevel`] the reader currently uses
    fn detail_level(&self) -> &DetailLevel;

    /// Set the [`DetailLevel`] for the reader, changing
    /// the amount of work done immediately on loading a
    /// spectrum.
    ///
    /// # Note
    /// Not all readers support all detail levels, and the
    /// behavior when requesting one of those levels will
    /// depend upon the underlying reader.
    fn set_detail_level(&mut self, detail_level: DetailLevel);

    /// Retrieve a frame by it's native ID
    fn get_frame_by_id(&mut self, id: &str) -> Option<S>;

    /// Retrieve a frame by it's integer index
    fn get_frame_by_index(&mut self, index: usize) -> Option<S>;

    /// Retrieve a frame by its scan start time
    /// Considerably more complex than seeking by ID or index, this involves
    /// a binary search over the frame index and assumes that frames are stored
    /// in chronological order.
    fn get_frame_by_time(&mut self, time: f64) -> Option<S> {
        let n = self.len();
        let mut lo: usize = 0;
        let mut hi: usize = n;

        let mut best_error: f64 = f64::INFINITY;
        let mut best_match: Option<S> = None;

        if lo == hi {
            return None;
        }
        while hi != lo {
            let mid = (hi + lo) / 2;
            let scan = self.get_frame_by_index(mid)?;
            let scan_time = scan.start_time();
            let err = (scan_time - time).abs();

            if err < best_error {
                best_error = err;
                best_match = Some(scan);
            }
            if hi.saturating_sub(1) == lo {
                return best_match;
            } else if scan_time > time {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        best_match
    }

    /// Retrieve the number of frames in source file
    fn len(&self) -> usize {
        self.get_index().len()
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Access the frame offset index to enumerate all frames by ID
    fn get_index(&self) -> &OffsetIndex;

    /// Set the frame offset index. This method shouldn't be needed if not writing
    /// a new adapter
    fn set_index(&mut self, index: OffsetIndex);

    /// Helper method to support seeking to an ID
    fn _offset_of_id(&self, id: &str) -> Option<u64> {
        self.get_index().get(id)
    }

    /// Helper method to support seeking to an index
    fn _offset_of_index(&self, index: usize) -> Option<u64> {
        self.get_index()
            .get_index(index)
            .map(|(_id, offset)| offset)
    }

    /// Helper method to support seeking to a specific time.
    /// Considerably more complex than seeking by ID or index.
    fn _offset_of_time(&mut self, time: f64) -> Option<u64> {
        match self.get_frame_by_time(time) {
            Some(scan) => self._offset_of_index(scan.index()),
            None => None,
        }
    }

    /// Get the nth [`IonMobilityFrameGroup`] from this source
    fn get_group_by_index(&mut self, index: usize) -> Option<IonMobilityFrameGroup<C, D, S>>
    where
        Self: Sized,
    {
        self.groups().nth(index)
    }

    /// Open a new iterator over this stream
    fn iter(&mut self) -> IonMobilityFrameIterator<'_, C, D, S, Self>
    where
        Self: Sized,
    {
        IonMobilityFrameIterator::new(self)
    }

    /// Create a new [`IonMobilityFrameIterator`] over `self` and use that state to drive a [`IonMobilityFrameGroupingIterator`]
    fn groups(
        &mut self,
    ) -> IonMobilityFrameGroupingIterator<IonMobilityFrameIterator<'_, C, D, S, Self>, C, D, S>
    where
        Self: Sized,
    {
        IonMobilityFrameGroupingIterator::new(self.iter())
    }

    /// Consume `self` to create a [`IonMobilityFrameGroupingIterator`]. This is ideal for non-rewindable streams
    /// like [`io::stdin`] which don't implement [`io::Seek`]
    fn into_groups(self) -> IonMobilityFrameGroupingIterator<Self, C, D, S>
    where
        Self: Sized,
    {
        IonMobilityFrameGroupingIterator::new(self)
    }
}

/// A generic iterator over a [`IonMobilityFrameSource`] implementer that assumes the
/// source has already been indexed. Otherwise, the source's own iterator
/// behavior should be used.
pub struct IonMobilityFrameIterator<
    'lifespan,
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
    S: IonMobilityFrameLike<C, D>,
    R: IonMobilityFrameSource<C, D, S>,
> {
    source: &'lifespan mut R,
    frame_type: PhantomData<S>,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
    index: usize,
    back_index: usize,
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        R: IonMobilityFrameSource<C, D, S>,
    > IonMobilityFrameIterator<'_, C, D, S, R>
{
    pub fn new(source: &mut R) -> IonMobilityFrameIterator<'_, C, D, S, R> {
        IonMobilityFrameIterator::<C, D, S, R> {
            source,
            index: 0,
            back_index: 0,
            frame_type: PhantomData,
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
        }
    }
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        R: IonMobilityFrameSource<C, D, S>,
    > Iterator for IonMobilityFrameIterator<'_, C, D, S, R>
{
    type Item = S;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index + self.back_index >= self.len() {
            return None;
        }
        let result = self.source.get_frame_by_index(self.index);
        self.index += 1;
        result
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.index += n;
        self.next()
    }
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        R: IonMobilityFrameSource<C, D, S>,
    > ExactSizeIterator for IonMobilityFrameIterator<'_, C, D, S, R>
{
    fn len(&self) -> usize {
        self.source.len()
    }
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        R: IonMobilityFrameSource<C, D, S>,
    > DoubleEndedIterator for IonMobilityFrameIterator<'_, C, D, S, R>
{
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.index + self.back_index >= self.len() {
            return None;
        };
        let i = self.len() - (self.back_index + 1);
        let result = self.source.get_frame_by_index(i);
        self.back_index += 1;
        result
    }
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        R: IonMobilityFrameSource<C, D, S>,
    > IonMobilityFrameSource<C, D, S> for IonMobilityFrameIterator<'_, C, D, S, R>
{
    fn reset(&mut self) {
        self.index = 0;
        self.back_index = 0;
    }

    fn get_frame_by_id(&mut self, id: &str) -> Option<S> {
        self.source.get_frame_by_id(id)
    }

    fn get_frame_by_index(&mut self, index: usize) -> Option<S> {
        self.source.get_frame_by_index(index)
    }

    fn get_frame_by_time(&mut self, time: f64) -> Option<S> {
        self.source.get_frame_by_time(time)
    }

    fn get_index(&self) -> &OffsetIndex {
        self.source.get_index()
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.source.set_index(index);
    }

    fn detail_level(&self) -> &DetailLevel {
        self.source.detail_level()
    }

    fn set_detail_level(&mut self, detail_level: DetailLevel) {
        self.source.set_detail_level(detail_level);
    }
}

/// If the underlying iterator implements [`MSDataFileMetadata`] then [`IonMobilityFrameIterator`] will
/// forward that implementation, assuming it is available.
impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        R: IonMobilityFrameSource<C, D, S>,
    > MSDataFileMetadata for IonMobilityFrameIterator<'_, C, D, S, R>
where
    R: MSDataFileMetadata,
{
    crate::delegate_impl_metadata_trait!(source);
}

/// Adapt a [`SpectrumSource`] that contains spectra with a non-scalar ion mobility
/// dimension to a [`IonMobilityFrameSource`].
#[derive(Debug)]
pub struct Generic3DIonMobilityFrameSource<
    CP: CentroidLike,
    DP: DeconvolutedCentroidLike,
    R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
    C: FeatureLike<MZ, IonMobility> = Feature<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge = ChargedFeature<Mass, IonMobility>,
> {
    source: R,
    _cp: PhantomData<CP>,
    _dp: PhantomData<DP>,
    _c: PhantomData<C>,
    _d: PhantomData<D>,
}

impl<
        CP: CentroidLike,
        DP: DeconvolutedCentroidLike,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > MSDataFileMetadata for Generic3DIonMobilityFrameSource<CP, DP, R, C, D>
where
    R: MSDataFileMetadata,
{
    crate::delegate_impl_metadata_trait!(source);
}

impl<
        CP: CentroidLike,
        DP: DeconvolutedCentroidLike,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > IonMobilityFrameSource<C, D, MultiLayerIonMobilityFrame<C, D>>
    for Generic3DIonMobilityFrameSource<CP, DP, R, C, D>
{
    fn detail_level(&self) -> &DetailLevel {
        self.source.detail_level()
    }

    fn set_detail_level(&mut self, detail_level: DetailLevel) {
        self.source.set_detail_level(detail_level);
    }

    fn reset(&mut self) {
        self.source.reset()
    }

    fn get_frame_by_id(&mut self, id: &str) -> Option<MultiLayerIonMobilityFrame<C, D>> {
        self.source.get_spectrum_by_id(id).and_then(|s| {
            MultiLayerIonMobilityFrame::try_from(s).map_or_else(
                |err| {
                    warn!("Failed to convert {id} to MultiLayerIonMobilityFrame: {err}");
                    None
                },
                Some,
            )
        })
    }

    fn get_frame_by_index(&mut self, index: usize) -> Option<MultiLayerIonMobilityFrame<C, D>> {
        self.source.get_spectrum_by_index(index).and_then(|s| {
            MultiLayerIonMobilityFrame::try_from(s).map_or_else(
                |err| {
                    warn!("Failed to convert {index} to MultiLayerIonMobilityFrame: {err}");
                    None
                },
                Some,
            )
        })
    }

    fn get_frame_by_time(&mut self, time: f64) -> Option<MultiLayerIonMobilityFrame<C, D>> {
        self.source.get_spectrum_by_time(time).and_then(|s| {
            MultiLayerIonMobilityFrame::try_from(s).map_or_else(
                |err| {
                    warn!("Failed to convert {time} to MultiLayerIonMobilityFrame: {err}");
                    None
                },
                Some,
            )
        })
    }

    fn get_index(&self) -> &OffsetIndex {
        self.source.get_index()
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.source.set_index(index)
    }
}

impl<
        CP: CentroidLike,
        DP: DeconvolutedCentroidLike,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > Iterator for Generic3DIonMobilityFrameSource<CP, DP, R, C, D>
{
    type Item = MultiLayerIonMobilityFrame<C, D>;

    fn next(&mut self) -> Option<Self::Item> {
        self.source
            .by_ref()
            .filter_map(|s| MultiLayerIonMobilityFrame::try_from(s).ok())
            .next()
    }
}

impl<
        CP: CentroidLike,
        DP: DeconvolutedCentroidLike,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > Generic3DIonMobilityFrameSource<CP, DP, R, C, D>
{
    pub fn new(source: R) -> Self {
        Self {
            source,
            _cp: PhantomData,
            _dp: PhantomData,
            _c: PhantomData,
            _d: PhantomData,
        }
    }

    pub fn get_inner(&self) -> &R {
        &self.source
    }

    pub fn get_mut(&mut self) -> &mut R {
        &mut self.source
    }

    pub fn into_inner(self) -> R {
        self.source
    }
}

/// Errors that may occur when reading a spectrum from a [`RandomAccessIonMobilityFrameIterator`]
#[derive(Debug, Error)]
pub enum IonMobilityFrameAccessError {
    /// An undetermined error failing to locate the requested frame
    #[error("The requested frame was not found")]
    FrameNotFound,

    /// An error resolving a frame by it's native ID
    #[error("The requested frame native ID {0} was not found")]
    FrameIdNotFound(String),

    /// An error resolving a frame by it's index
    #[error("The requested frame index {0} was not found")]
    FrameIndexNotFound(usize),

    /// An I/O error prevented reading the frame, even if it could be found.
    #[error("I/O error occurred while reading: {0:?}")]
    IOError(#[source] Option<io::Error>),
}

impl From<IonMobilityFrameAccessError> for io::Error {
    fn from(value: IonMobilityFrameAccessError) -> Self {
        let s = value.to_string();
        match value {
            IonMobilityFrameAccessError::FrameNotFound => {
                io::Error::new(io::ErrorKind::NotFound, s)
            }
            IonMobilityFrameAccessError::FrameIdNotFound(_) => {
                io::Error::new(io::ErrorKind::NotFound, s)
            }
            IonMobilityFrameAccessError::FrameIndexNotFound(_) => {
                io::Error::new(io::ErrorKind::NotFound, s)
            }
            IonMobilityFrameAccessError::IOError(e) => match e {
                Some(e) => e,
                None => io::Error::other(s),
            },
        }
    }
}

impl From<IonMobilityFrameAccessError> for SpectrumAccessError {
    fn from(value: IonMobilityFrameAccessError) -> Self {
        match value {
            IonMobilityFrameAccessError::FrameNotFound => SpectrumAccessError::SpectrumNotFound,
            IonMobilityFrameAccessError::FrameIdNotFound(id) => {
                SpectrumAccessError::SpectrumIdNotFound(id)
            }
            IonMobilityFrameAccessError::FrameIndexNotFound(i) => {
                SpectrumAccessError::SpectrumIndexNotFound(i)
            }
            IonMobilityFrameAccessError::IOError(error) => SpectrumAccessError::IOError(error),
        }
    }
}

impl From<SpectrumAccessError> for IonMobilityFrameAccessError {
    fn from(value: SpectrumAccessError) -> Self {
        match value {
            SpectrumAccessError::SpectrumNotFound => IonMobilityFrameAccessError::FrameNotFound,
            SpectrumAccessError::SpectrumIdNotFound(id) => {
                IonMobilityFrameAccessError::FrameIdNotFound(id)
            }
            SpectrumAccessError::SpectrumIndexNotFound(i) => {
                IonMobilityFrameAccessError::FrameIndexNotFound(i)
            }
            SpectrumAccessError::IOError(error) => IonMobilityFrameAccessError::IOError(error),
        }
    }
}

/// An extension of [`IonMobilityFrameSource`] that supports relocatable iteration relative to a
/// specific spectrum coordinate or identifier.
pub trait RandomAccessIonMobilityFrameIterator<
    C: FeatureLike<MZ, IonMobility> = Feature<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge = ChargedFeature<Mass, IonMobility>,
    S: IonMobilityFrameLike<C, D> = MultiLayerIonMobilityFrame<C, D>,
>: IonMobilityFrameSource<C, D, S>
{
    /// Start iterating from the frame whose native ID matches `id`
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, IonMobilityFrameAccessError>;

    /// Start iterating from the frame whose index is `index`
    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, IonMobilityFrameAccessError>;

    /// Start iterating from the frame starting closest to `time`
    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, IonMobilityFrameAccessError>;
}

impl<
        CP: CentroidLike,
        DP: DeconvolutedCentroidLike,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > RandomAccessIonMobilityFrameIterator<C, D, MultiLayerIonMobilityFrame<C, D>>
    for Generic3DIonMobilityFrameSource<CP, DP, R, C, D>
where
    R: RandomAccessSpectrumIterator<CP, DP, MultiLayerSpectrum<CP, DP>>,
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, IonMobilityFrameAccessError> {
        match self.source.start_from_id(id) {
            Ok(_) => Ok(self),
            Err(e) => Err(IonMobilityFrameAccessError::from(e)),
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, IonMobilityFrameAccessError> {
        match self.source.start_from_index(index) {
            Ok(_) => Ok(self),
            Err(e) => Err(IonMobilityFrameAccessError::from(e)),
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, IonMobilityFrameAccessError> {
        match self.source.start_from_time(time) {
            Ok(_) => Ok(self),
            Err(e) => Err(IonMobilityFrameAccessError::from(e)),
        }
    }
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        R: IonMobilityFrameSource<C, D, S>,
    > RandomAccessIonMobilityFrameIterator<C, D, S> for IonMobilityFrameIterator<'_, C, D, S, R>
{
    /// Start iterating from the spectrum whose native ID matches `id`
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, IonMobilityFrameAccessError> {
        if let Some(scan) = self.get_frame_by_id(id) {
            self.index = scan.index();
            self.back_index = 0;
            Ok(self)
        } else if self.get_index().contains_key(id) {
            Err(IonMobilityFrameAccessError::IOError(None))
        } else {
            Err(IonMobilityFrameAccessError::FrameIdNotFound(id.to_string()))
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, IonMobilityFrameAccessError> {
        if index < self.len() {
            self.index = index;
            self.back_index = 0;
            Ok(self)
        } else {
            Err(IonMobilityFrameAccessError::FrameIndexNotFound(index))
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, IonMobilityFrameAccessError> {
        if let Some(scan) = self.get_frame_by_time(time) {
            self.index = scan.index();
            self.back_index = 0;
            Ok(self)
        } else if self
            .get_frame_by_index(self.len() - 1)
            .expect("Failed to fetch spectrum for boundary testing")
            .start_time()
            < time
        {
            Err(IonMobilityFrameAccessError::FrameNotFound)
        } else {
            Err(IonMobilityFrameAccessError::IOError(None))
        }
    }
}

/// Common interface for ion mobility frame writing.
pub trait IonMobilityFrameWriter<
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
>
{
    /// Write out a single frame
    fn write_frame<S: IonMobilityFrameLike<C, D> + 'static>(
        &mut self,
        spectrum: &S,
    ) -> io::Result<usize>;

    /// Write out a single owned frame.
    ///
    /// This may produce fewer copies for some implementations, but the default implementation
    /// just delegates to [`IonMobilityFrameWriter::write_frame`]
    fn write_frame_owned<S: IonMobilityFrameLike<C, D> + 'static>(
        &mut self,
        spectrum: S,
    ) -> io::Result<usize> {
        self.write_frame(&spectrum)
    }

    /// As [`std::io::Write::flush`]
    fn flush_frame(&mut self) -> io::Result<()>;

    /// Consume an [`Iterator`] over [`IonMobilityFrameLike`] references
    fn write_all_frames<'b, S: IonMobilityFrameLike<C, D> + 'static, T: Iterator<Item = &'b S>>(
        &mut self,
        iterator: T,
    ) -> io::Result<usize> {
        let mut n = 0;
        for spectrum in iterator {
            n += self.write_frame(spectrum)?;
        }
        Ok(n)
    }

    /// Consume an [`Iterator`] over [`IonMobilityFrameLike`]
    fn write_all_frames_owned<S: IonMobilityFrameLike<C, D> + 'static, T: Iterator<Item = S>>(
        &mut self,
        iterator: T,
    ) -> io::Result<usize> {
        let mut n = 0;
        for spectrum in iterator {
            n += self.write_frame_owned(spectrum)?;
        }
        Ok(n)
    }

    /// Write a [`IonMobilityFrameGrouping`] out in order
    fn write_frame_group<
        S: IonMobilityFrameLike<C, D> + 'static,
        G: IonMobilityFrameGrouping<C, D, S> + 'static,
    >(
        &mut self,
        group: &G,
    ) -> io::Result<usize> {
        let mut n = 0;
        if let Some(precursor) = group.precursor() {
            n += self.write_frame(precursor)?;
        }
        for product in group.products() {
            n += self.write_frame(product)?;
        }
        Ok(n)
    }

    /// Write an owned [`IonMobilityFrameGrouping`] out in order
    ///
    /// This may produce fewer copies for some implementations.
    fn write_frame_group_owned<
        S: IonMobilityFrameLike<C, D> + 'static,
        G: IonMobilityFrameGrouping<C, D, S> + 'static,
    >(
        &mut self,
        group: G,
    ) -> io::Result<usize> {
        let (precursor, products) = group.into_parts();
        let mut n = 0;
        if let Some(precursor) = precursor {
            n += self.write_frame_owned(precursor)?;
        }
        for product in products {
            n += self.write_frame_owned(product)?;
        }
        Ok(n)
    }

    /// Consume an [`Iterator`] over [`IonMobilityFrameGrouping`] references
    fn write_all_frame_groups<
        'b,
        S: IonMobilityFrameLike<C, D> + 'static,
        G: IonMobilityFrameGrouping<C, D, S> + 'static,
        T: Iterator<Item = &'b G>,
    >(
        &mut self,
        iterator: T,
    ) -> io::Result<usize> {
        let mut n = 0;
        for group in iterator {
            n += self.write_frame_group(group)?;
        }
        Ok(n)
    }

    /// Consume an [`Iterator`] over [`IonMobilityFrameGrouping`]
    fn write_all_frame_groups_owned<
        S: IonMobilityFrameLike<C, D> + 'static,
        G: IonMobilityFrameGrouping<C, D, S> + 'static,
        T: Iterator<Item = G>,
    >(
        &mut self,
        iterator: T,
    ) -> io::Result<usize> {
        let mut n = 0;
        for group in iterator {
            n += self.write_frame_group_owned(group)?;
        }
        Ok(n)
    }

    /// Completes the data file format, preventing new data from being able incorporate additional
    /// data. Does not formally close the underlying writing stream.
    fn close_frames(&mut self) -> io::Result<()>;
}

/// Adapt a [`SpectrumSource`] that contains spectra with a non-scalar ion mobility
/// dimension to a [`IonMobilityFrameSource`].
#[derive(Debug)]
pub struct BorrowedGeneric3DIonMobilityFrameSource<
    'a,
    CP: CentroidLike,
    DP: DeconvolutedCentroidLike,
    R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
    C: FeatureLike<MZ, IonMobility> = Feature<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge = ChargedFeature<Mass, IonMobility>,
> {
    source: &'a mut R,
    _cp: PhantomData<CP>,
    _dp: PhantomData<DP>,
    _c: PhantomData<C>,
    _d: PhantomData<D>,
}

impl<
        CP: CentroidLike,
        DP: DeconvolutedCentroidLike,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > MSDataFileMetadata for BorrowedGeneric3DIonMobilityFrameSource<'_, CP, DP, R, C, D>
where
    R: MSDataFileMetadata,
{
    crate::delegate_impl_metadata_trait!(source);
}

impl<
        CP: CentroidLike,
        DP: DeconvolutedCentroidLike,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > IonMobilityFrameSource<C, D, MultiLayerIonMobilityFrame<C, D>>
    for BorrowedGeneric3DIonMobilityFrameSource<'_, CP, DP, R, C, D>
{
    fn detail_level(&self) -> &DetailLevel {
        self.source.detail_level()
    }

    fn set_detail_level(&mut self, detail_level: DetailLevel) {
        self.source.set_detail_level(detail_level);
    }

    fn reset(&mut self) {
        self.source.reset()
    }

    fn get_frame_by_id(&mut self, id: &str) -> Option<MultiLayerIonMobilityFrame<C, D>> {
        self.source.get_spectrum_by_id(id).and_then(|s| {
            MultiLayerIonMobilityFrame::try_from(s).map_or_else(
                |err| {
                    warn!("Failed to convert {id} to MultiLayerIonMobilityFrame: {err}");
                    None
                },
                Some,
            )
        })
    }

    fn get_frame_by_index(&mut self, index: usize) -> Option<MultiLayerIonMobilityFrame<C, D>> {
        self.source.get_spectrum_by_index(index).and_then(|s| {
            MultiLayerIonMobilityFrame::try_from(s).map_or_else(
                |err| {
                    warn!("Failed to convert {index} to MultiLayerIonMobilityFrame: {err}");
                    None
                },
                Some,
            )
        })
    }

    fn get_frame_by_time(&mut self, time: f64) -> Option<MultiLayerIonMobilityFrame<C, D>> {
        self.source.get_spectrum_by_time(time).and_then(|s| {
            MultiLayerIonMobilityFrame::try_from(s).map_or_else(
                |err| {
                    warn!("Failed to convert {time} to MultiLayerIonMobilityFrame: {err}");
                    None
                },
                Some,
            )
        })
    }

    fn get_index(&self) -> &OffsetIndex {
        self.source.get_index()
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.source.set_index(index)
    }
}

impl<
        CP: CentroidLike,
        DP: DeconvolutedCentroidLike,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > Iterator for BorrowedGeneric3DIonMobilityFrameSource<'_, CP, DP, R, C, D>
{
    type Item = MultiLayerIonMobilityFrame<C, D>;

    fn next(&mut self) -> Option<Self::Item> {
        self.source
            .by_ref()
            .filter_map(|s| MultiLayerIonMobilityFrame::try_from(s).ok())
            .next()
    }
}

impl<
        'a,
        CP: CentroidLike,
        DP: DeconvolutedCentroidLike,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > BorrowedGeneric3DIonMobilityFrameSource<'a, CP, DP, R, C, D>
{
    #[allow(unused)]
    pub fn new(source: &'a mut R) -> Self {
        Self {
            source,
            _cp: PhantomData,
            _dp: PhantomData,
            _c: PhantomData,
            _d: PhantomData,
        }
    }
}

impl<
        CP: CentroidLike,
        DP: DeconvolutedCentroidLike,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > RandomAccessIonMobilityFrameIterator<C, D, MultiLayerIonMobilityFrame<C, D>>
    for BorrowedGeneric3DIonMobilityFrameSource<'_, CP, DP, R, C, D>
where
    R: RandomAccessSpectrumIterator<CP, DP, MultiLayerSpectrum<CP, DP>>,
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, IonMobilityFrameAccessError> {
        match self.source.start_from_id(id) {
            Ok(_) => Ok(self),
            Err(e) => Err(match e {
                SpectrumAccessError::SpectrumIdNotFound(id) => {
                    IonMobilityFrameAccessError::FrameIdNotFound(id)
                }
                SpectrumAccessError::SpectrumNotFound => IonMobilityFrameAccessError::FrameNotFound,
                SpectrumAccessError::IOError(e) => IonMobilityFrameAccessError::IOError(e),
                _ => todo!(),
            }),
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, IonMobilityFrameAccessError> {
        match self.source.start_from_index(index) {
            Ok(_) => Ok(self),
            Err(e) => Err(match e {
                SpectrumAccessError::SpectrumIndexNotFound(i) => {
                    IonMobilityFrameAccessError::FrameIndexNotFound(i)
                }
                SpectrumAccessError::SpectrumNotFound => IonMobilityFrameAccessError::FrameNotFound,
                SpectrumAccessError::IOError(e) => IonMobilityFrameAccessError::IOError(e),
                _ => todo!(),
            }),
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, IonMobilityFrameAccessError> {
        match self.source.start_from_time(time) {
            Ok(_) => Ok(self),
            Err(e) => Err(match e {
                SpectrumAccessError::SpectrumNotFound => IonMobilityFrameAccessError::FrameNotFound,
                SpectrumAccessError::IOError(e) => IonMobilityFrameAccessError::IOError(e),
                _ => todo!(),
            }),
        }
    }
}

/// Analogous to to [`RandomAccessIonMobilityFrameIterator`], but for [`IonMobilityFrameGrouping`] implementations.
pub trait RandomAccessIonMobilityFrameGroupingIterator<
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
    S: IonMobilityFrameLike<C, D> = MultiLayerIonMobilityFrame<C, D>,
    G: IonMobilityFrameGrouping<C, D, S> = IonMobilityFrameGroup<
        C,
        D,
        MultiLayerIonMobilityFrame<C, D>,
    >,
>: Iterator<Item = G>
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, IonMobilityFrameAccessError>;
    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, IonMobilityFrameAccessError>;
    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, IonMobilityFrameAccessError>;
    fn reset_state(&mut self);
}

#[derive(Debug, Error)]
pub enum IntoIonMobilityFrameSourceError {
    #[error("No ion mobility frames were found")]
    NoIonMobilityFramesFound,

    #[error("Cannot convert to requested frame type")]
    ConversionNotPossible,
}

/// Convert a [`SpectrumSource`] to an [`IonMobilityFrameSource`] if it detects ion mobility frames (3D spectra)
pub trait IntoIonMobilityFrameSource<C: CentroidLike, D: DeconvolutedCentroidLike>:
    SpectrumSource<C, D, MultiLayerSpectrum<C, D>> + Sized
{
    /// The [`IonMobilityFrameSource`]-implementing type for this [`SpectrumSource`].
    ///
    /// When another type isn't available, [`Generic3DIonMobilityFrameSource`]
    type IonMobilityFrameSource<CF: FeatureLike<MZ, IonMobility>, DF: FeatureLike<Mass, IonMobility> + KnownCharge>: IonMobilityFrameSource<CF, DF, MultiLayerIonMobilityFrame<CF, DF>>;

    /// Attempt to convert the [`SpectrumSource`] into an [`IonMobilityFrameSource`], returning [`IntoIonMobilityFrameSourceError`]
    /// if it is not possible
    fn try_into_frame_source<
        CF: FeatureLike<MZ, IonMobility>,
        DF: FeatureLike<Mass, IonMobility> + KnownCharge,
    >(
        self,
    ) -> Result<Self::IonMobilityFrameSource<CF, DF>, IntoIonMobilityFrameSourceError>;

    /// Call [`IntoIonMobilityFrameSource::try_into_frame_source`], panicking if an error is returned.
    fn into_frame_source<
        CF: FeatureLike<MZ, IonMobility>,
        DF: FeatureLike<Mass, IonMobility> + KnownCharge,
    >(
        self,
    ) -> Self::IonMobilityFrameSource<CF, DF> {
        self.try_into_frame_source().unwrap()
    }

    /// Reads a sparse 1% of the entries from the [`SpectrumSource`], testing
    /// for the presence of ion mobility data.
    fn has_ion_mobility(&mut self) -> Option<HasIonMobility> {
        let details = *self.detail_level();
        self.set_detail_level(DetailLevel::Lazy);
        let n = self.len();
        let mut handle = self.iter();
        let mut status = HasIonMobility::None;
        let step_size = if n > 100 { n / 100 } else { n };
        for i in (0..n).step_by(step_size) {
            let spec = handle.get_spectrum_by_index(i)?;
            let cls = spec.has_ion_mobility_class();
            status = status.max(cls);
            if status > HasIonMobility::None {
                break;
            }
        }
        self.set_detail_level(details);
        Some(status)
    }
}
