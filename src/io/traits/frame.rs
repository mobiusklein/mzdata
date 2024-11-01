use std::convert::TryFrom;
use std::io;
use std::marker::PhantomData;

use thiserror::Error;
use log::warn;

use mzpeaks::{
    IonMobility,
    KnownCharge, Mass, MZ,
    feature::{ChargedFeature, Feature, FeatureLike}
};

use crate::{io::OffsetIndex, prelude::MSDataFileMetadata};
use crate::spectrum::group::IonMobilityFrameGroupingIterator;
use crate::spectrum::spectrum_types::MultiLayerSpectrum;
use crate::spectrum::{
    CentroidPeakAdapting, DeconvolutedPeakAdapting, IonMobilityFrameLike,
    MultiLayerIonMobilityFrame,
};

use super::{RandomAccessSpectrumIterator, SpectrumAccessError, SpectrumSource};



/// An analog of [`SpectrumSource`] for [`IonMobilityFrameLike`] producing types
pub trait IonMobilityFrameSource<
    C: FeatureLike<MZ, IonMobility> = Feature<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge = ChargedFeature<Mass, IonMobility>,
    S: IonMobilityFrameLike<C, D> = MultiLayerIonMobilityFrame<C, D>,
>: Iterator<Item = S>
{
    fn reset(&mut self);

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

    fn iter(&mut self) -> IonMobilityFrameIterator<C, D, S, Self>
    where
        Self: Sized,
    {
        IonMobilityFrameIterator::new(self)
    }

    /// Create a new [`IonMobilityFrameIterator`] over `self` and use that state to drive a [`IonMobilityFrameGroupingIterator`]
    fn groups(&mut self) -> IonMobilityFrameGroupingIterator<IonMobilityFrameIterator<'_, C, D, S, Self>, C, D, S>
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
        'lifespan,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        R: IonMobilityFrameSource<C, D, S>,
    > IonMobilityFrameIterator<'lifespan, C, D, S, R>
{
    pub fn new(source: &mut R) -> IonMobilityFrameIterator<C, D, S, R> {
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
        'lifespan,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        R: IonMobilityFrameSource<C, D, S>,
    > Iterator for IonMobilityFrameIterator<'lifespan, C, D, S, R>
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
}

impl<
        'lifespan,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        R: IonMobilityFrameSource<C, D, S>,
    > ExactSizeIterator for IonMobilityFrameIterator<'lifespan, C, D, S, R>
{
    fn len(&self) -> usize {
        self.source.len()
    }
}

impl<
        'lifespan,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        R: IonMobilityFrameSource<C, D, S>,
    > DoubleEndedIterator for IonMobilityFrameIterator<'lifespan, C, D, S, R>
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
        'lifespan,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        R: IonMobilityFrameSource<C, D, S>,
    > IonMobilityFrameSource<C, D, S> for IonMobilityFrameIterator<'lifespan, C, D, S, R>
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
}

/// If the underlying iterator implements [`MSDataFileMetadata`] then [`IonMobilityFrameIterator`] will
/// forward that implementation, assuming it is available.
impl<
        'lifespan,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        R: IonMobilityFrameSource<C, D, S>,
    > MSDataFileMetadata for IonMobilityFrameIterator<'lifespan, C, D, S, R>
where
    R: MSDataFileMetadata,
{
    crate::delegate_impl_metadata_trait!(source);
}

/// Adapt a [`SpectrumSource`] that contains spectra with a non-scalar ion mobility
/// dimension to a [`IonMobilityFrameSource`].
#[derive(Debug)]
pub struct Generic3DIonMobilityFrameSource<
    CP: CentroidPeakAdapting,
    DP: DeconvolutedPeakAdapting,
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
        CP: CentroidPeakAdapting,
        DP: DeconvolutedPeakAdapting,
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
        CP: CentroidPeakAdapting,
        DP: DeconvolutedPeakAdapting,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > IonMobilityFrameSource<C, D, MultiLayerIonMobilityFrame<C, D>>
    for Generic3DIonMobilityFrameSource<CP, DP, R, C, D>
{
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
                |val| Some(val),
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
                |val| Some(val),
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
                |val| Some(val),
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
        CP: CentroidPeakAdapting,
        DP: DeconvolutedPeakAdapting,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > Iterator for Generic3DIonMobilityFrameSource<CP, DP, R, C, D>
{
    type Item = MultiLayerIonMobilityFrame<C, D>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(s) = self.source.next() {
            Some(
                MultiLayerIonMobilityFrame::try_from(s)
                    .expect("Failed to convert spectrum into frame"),
            )
        } else {
            None
        }
    }
}

impl<
        CP: CentroidPeakAdapting,
        DP: DeconvolutedPeakAdapting,
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

/// Errors that may occur when reading a spectrum from a [`RandomAccessSpectrumIterator`]
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
                None => io::Error::new(io::ErrorKind::Other, s),
            },
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
        CP: CentroidPeakAdapting,
        DP: DeconvolutedPeakAdapting,
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

impl<
        'lifespan,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        R: IonMobilityFrameSource<C, D, S>,
    > RandomAccessIonMobilityFrameIterator<C, D, S>
    for IonMobilityFrameIterator<'lifespan, C, D, S, R>
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

/// An abstraction over [`IonMobilityFrameGroup`](crate::spectrum::IonMobilityFrameGroup)'s interface.
pub trait IonMobilityFrameGrouping<
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
    S: IonMobilityFrameLike<C, D> = MultiLayerIonMobilityFrame<C, D>,
>: Default
{
    /// Get the precursor spectrum, which may be absent
    fn precursor(&self) -> Option<&S>;
    /// Get a mutable reference to the precursor spectrum, which may be absent
    fn precursor_mut(&mut self) -> Option<&mut S>;
    /// Explicitly set the precursor spectrum directly.
    fn set_precursor(&mut self, prec: S);

    /// Get a reference to the collection of product frames
    fn products(&self) -> &[S];

    /// Get a mutable reference to the collection of product frames
    fn products_mut(&mut self) -> &mut Vec<S>;

    /// The total number of frames in the group
    fn total_frames(&self) -> usize {
        self.precursor().is_some() as usize + self.products().len()
    }

    /// The frame that occurred first chronologically
    fn earliest_frame(&self) -> Option<&S> {
        self.precursor().or_else(|| {
            self.products().iter().min_by(|a, b| {
                a.acquisition()
                    .start_time()
                    .total_cmp(&b.acquisition().start_time())
            })
        })
    }

    /// The frame that occurred last chronologically
    fn latest_frame(&self) -> Option<&S> {
        self.precursor().or_else(|| {
            self.products().iter().max_by(|a, b| {
                a.acquisition()
                    .start_time()
                    .total_cmp(&b.acquisition().start_time())
            })
        })
    }

    /// The lowest MS level in the group
    fn lowest_ms_level(&self) -> Option<u8> {
        let prec_level = self.precursor().map(|p| p.ms_level()).unwrap_or(u8::MAX);
        let val = self
            .products()
            .iter()
            .fold(prec_level, |state, s| state.min(s.ms_level()));
        if val > 0 {
            Some(val)
        } else {
            None
        }
    }

    /// The highest MS level in the group
    fn highest_ms_level(&self) -> Option<u8> {
        let prec_level = self
            .precursor()
            .map(|p| p.ms_level())
            .unwrap_or_else(|| u8::MIN);
        let val = self
            .products()
            .iter()
            .fold(prec_level, |state, s| state.max(s.ms_level()));
        if val > 0 {
            Some(val)
        } else {
            None
        }
    }

    /// Decompose the group into its components, discarding any additional metrics
    fn into_parts(self) -> (Option<S>, Vec<S>);
}


/// Common interface for ion mobility frame writing
pub trait IonMobilityFrameWriter<
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
>
{
    /// Write out a single frame
    fn write_frame<S: IonMobilityFrameLike<C, D> + 'static>(&mut self, spectrum: &S) -> io::Result<usize>;

    /// Write out a single owned frame.
    ///
    /// This may produce fewer copies for some implementations.
    fn write_frame_owned<S: IonMobilityFrameLike<C, D> + 'static>(&mut self, spectrum: S) -> io::Result<usize> {
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
    fn write_all_frames_owned<'b, S: IonMobilityFrameLike<C, D> + 'static, T: Iterator<Item = S>>(
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
    fn write_frame_group<S: IonMobilityFrameLike<C, D> + 'static, G: IonMobilityFrameGrouping<C, D, S> + 'static>(
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
        'b,
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
    CP: CentroidPeakAdapting,
    DP: DeconvolutedPeakAdapting,
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
        'a,
        CP: CentroidPeakAdapting,
        DP: DeconvolutedPeakAdapting,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > MSDataFileMetadata for BorrowedGeneric3DIonMobilityFrameSource<'a, CP, DP, R, C, D>
where
    R: MSDataFileMetadata,
{
    crate::delegate_impl_metadata_trait!(source);
}

impl<
        'a,
        CP: CentroidPeakAdapting,
        DP: DeconvolutedPeakAdapting,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > IonMobilityFrameSource<C, D, MultiLayerIonMobilityFrame<C, D>>
    for BorrowedGeneric3DIonMobilityFrameSource<'a, CP, DP, R, C, D>
{
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
                |val| Some(val),
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
                |val| Some(val),
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
                |val| Some(val),
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
        'a,
        CP: CentroidPeakAdapting,
        DP: DeconvolutedPeakAdapting,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > Iterator for BorrowedGeneric3DIonMobilityFrameSource<'a, CP, DP, R, C, D>
{
    type Item = MultiLayerIonMobilityFrame<C, D>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(s) = self.source.next() {
            Some(
                MultiLayerIonMobilityFrame::try_from(s)
                    .expect("Failed to convert spectrum into frame"),
            )
        } else {
            None
        }
    }
}

impl<
        'a,
        CP: CentroidPeakAdapting,
        DP: DeconvolutedPeakAdapting,
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
        'a,
        CP: CentroidPeakAdapting,
        DP: DeconvolutedPeakAdapting,
        R: SpectrumSource<CP, DP, MultiLayerSpectrum<CP, DP>>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
    > RandomAccessIonMobilityFrameIterator<C, D, MultiLayerIonMobilityFrame<C, D>>
    for BorrowedGeneric3DIonMobilityFrameSource<'a, CP, DP, R, C, D>
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


