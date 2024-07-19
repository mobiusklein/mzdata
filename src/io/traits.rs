use std::collections::{HashMap, VecDeque};
use std::convert::TryFrom;
use std::fs;
use std::io;
use std::iter::FusedIterator;
use std::marker::PhantomData;
use std::ops::Index;
use std::path::{self, PathBuf};
use std::sync::mpsc::Receiver;

use thiserror::Error;
use log::warn;

use mzpeaks::{
    CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak, IonMobility,
    KnownCharge, Mass, MZ,
    feature::{ChargedFeature, Feature, FeatureLike}
};

use crate::meta::{
    DataProcessing, FileDescription, InstrumentConfiguration, MassSpectrometryRun, Software,
};
use crate::prelude::MSDataFileMetadata;
use crate::spectrum::group::{IonMobilityFrameGroupingIterator, SpectrumGroup, SpectrumGroupingIterator};
use crate::spectrum::spectrum_types::{MultiLayerSpectrum, SpectrumLike};
use crate::spectrum::{
    CentroidPeakAdapting, Chromatogram, DeconvolutedPeakAdapting, IonMobilityFrameLike,
    MultiLayerIonMobilityFrame,
};

use super::utils::FileSource;
use super::OffsetIndex;

pub trait SeekRead: io::Read + io::Seek {}
impl<T: io::Read + io::Seek> SeekRead for T {}

/// A base trait defining the behaviors of a source of spectra.
///
/// A [`SpectrumSource`]
pub trait SpectrumSource<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumLike<C, D> = MultiLayerSpectrum<C, D>,
>: Iterator<Item = S>
{
    fn reset(&mut self);

    /// Retrieve a spectrum by it's native ID
    fn get_spectrum_by_id(&mut self, id: &str) -> Option<S>;

    /// Retrieve a spectrum by it's integer index
    fn get_spectrum_by_index(&mut self, index: usize) -> Option<S>;

    /// Retrieve a spectrum by its scan start time
    /// Considerably more complex than seeking by ID or index, this involves
    /// a binary search over the spectrum index and assumes that spectra are stored
    /// in chronological order.
    fn get_spectrum_by_time(&mut self, time: f64) -> Option<S> {
        let n = self.len();
        if n == 0 {
            if !self.get_index().init {
                warn!("Attempting to use `get_spectrum_by_time` when the spectrum index has not been initialized.");
                return None;
            }
        }
        let mut lo: usize = 0;
        let mut hi: usize = n;

        let mut best_error: f64 = f64::INFINITY;
        let mut best_match: Option<S> = None;

        if lo == hi {
            return None;
        }
        while hi != lo {
            let mid = (hi + lo) / 2;
            let scan = self.get_spectrum_by_index(mid)?;
            let scan_time = scan.start_time();
            let err = (scan_time - time).abs();

            if err < best_error {
                best_error = err;
                best_match = Some(scan);
            }
            if hi.saturating_sub(1) == lo {
                return best_match
            }
            else if scan_time > time {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        best_match
    }

    /// Retrieve the number of spectra in source file, usually by getting
    /// the length of the index. If the index isn't initialized, this will
    /// be 0.
    fn len(&self) -> usize {
        self.get_index().len()
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Access the spectrum offset index to enumerate all spectra by ID
    fn get_index(&self) -> &OffsetIndex;

    /// Set the spectrum offset index. This method shouldn't be needed if not writing
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
        match self.get_spectrum_by_time(time) {
            Some(scan) => self._offset_of_index(scan.index()),
            None => None,
        }
    }

    /// Open a new iterator over this stream.
    fn iter(&mut self) -> SpectrumIterator<C, D, S, Self>
    where
        Self: Sized,
    {
        SpectrumIterator::new(self)
    }

    /// Create a new `SpectrumIterator` over `self` and use that state to drive a `SpectrumGroupIterator`
    fn groups(&mut self) -> SpectrumGroupingIterator<SpectrumIterator<'_, C, D, S, Self>, C, D, S>
    where
        Self: Sized,
    {
        SpectrumGroupingIterator::new(self.iter())
    }

    /// Consume `self` to create a `SpectrumGroupIterator`. This is ideal for non-rewindable streams
    /// like [`io::stdin`] which don't implement [`io::Seek`]
    fn into_groups(self) -> SpectrumGroupingIterator<Self, C, D, S>
    where
        Self: Sized,
    {
        SpectrumGroupingIterator::new(self)
    }
}

/// A generic iterator over a [`SpectrumSource`] implementer that assumes the
/// source has already been indexed. Otherwise, the source's own iterator
/// behavior should be used.
pub struct SpectrumIterator<
    'lifespan,
    C: CentroidLike + Default,
    D: DeconvolutedCentroidLike + Default,
    S: SpectrumLike<C, D>,
    R: SpectrumSource<C, D, S>,
> {
    source: &'lifespan mut R,
    spectrum_type: PhantomData<S>,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
    index: usize,
    back_index: usize,
}

impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        R: SpectrumSource<C, D, S>,
        S: SpectrumLike<C, D>,
    > SpectrumIterator<'lifespan, C, D, S, R>
{
    pub fn new(source: &mut R) -> SpectrumIterator<C, D, S, R> {
        SpectrumIterator::<C, D, S, R> {
            source,
            index: 0,
            back_index: 0,
            spectrum_type: PhantomData,
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
        }
    }
}

impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        R: SpectrumSource<C, D, S>,
        S: SpectrumLike<C, D>,
    > Iterator for SpectrumIterator<'lifespan, C, D, S, R>
{
    type Item = S;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index + self.back_index >= self.len() {
            return None;
        }
        let result = self.source.get_spectrum_by_index(self.index);
        self.index += 1;
        result
    }
}

impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        R: SpectrumSource<C, D, S>,
        S: SpectrumLike<C, D>,
    > ExactSizeIterator for SpectrumIterator<'lifespan, C, D, S, R>
{
    fn len(&self) -> usize {
        self.source.len()
    }
}

impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        R: SpectrumSource<C, D, S>,
        S: SpectrumLike<C, D>,
    > DoubleEndedIterator for SpectrumIterator<'lifespan, C, D, S, R>
{
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.index + self.back_index >= self.len() {
            return None;
        };
        let i = self.len() - (self.back_index + 1);
        let result = self.source.get_spectrum_by_index(i);
        self.back_index += 1;
        result
    }
}

impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        R: SpectrumSource<C, D, S>,
    > SpectrumSource<C, D, S> for SpectrumIterator<'lifespan, C, D, S, R>
{
    fn reset(&mut self) {
        self.index = 0;
        self.back_index = 0;
    }

    fn get_spectrum_by_id(&mut self, id: &str) -> Option<S> {
        self.source.get_spectrum_by_id(id)
    }

    fn get_spectrum_by_index(&mut self, index: usize) -> Option<S> {
        self.source.get_spectrum_by_index(index)
    }

    fn get_spectrum_by_time(&mut self, time: f64) -> Option<S> {
        self.source.get_spectrum_by_time(time)
    }

    fn get_index(&self) -> &OffsetIndex {
        self.source.get_index()
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.source.set_index(index);
    }
}

/// If the underlying iterator implements [`MSDataFileMetadata`] then [`SpectrumIterator`] will
/// forward that implementation, assuming it is available.
impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        R: SpectrumSource<C, D, S>,
        S: SpectrumLike<C, D>,
    > MSDataFileMetadata for SpectrumIterator<'lifespan, C, D, S, R>
where
    R: MSDataFileMetadata,
{
    crate::delegate_impl_metadata_trait!(source);
}

/// A trait defining some helper methods to make efficient use of indices
/// automatic when opening a file from a path-like object.
pub trait MZFileReader<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumLike<C, D> = MultiLayerSpectrum<C, D>,
>: SpectrumSource<C, D, S> + Sized
{
    /// An on-trait method of constructing an index. Assumed
    /// to be a trivial wrapper.
    fn construct_index_from_stream(&mut self) -> u64;

    /// Re-construct an offset index from this readable object, assuming
    /// it is a JSON stream over the serialized index.
    fn read_index(&mut self, reader: Box<dyn io::Read>) -> Result<&mut Self, serde_json::Error> {
        match OffsetIndex::from_reader(reader) {
            Ok(index) => {
                self.set_index(index);
                Ok(self)
            }
            Err(err) => Err(err),
        }
    }

    fn write_index(&self, writer: Box<dyn io::Write>) -> Result<&Self, serde_json::Error> {
        match self.get_index().to_writer(writer) {
            Ok(_) => Ok(self),
            Err(err) => Err(err),
        }
    }

    /// The preferred method of opening a file from a path-like object.
    /// This method will open the file at the provided path, test whether
    /// there is an accompanied index file next to it on the file system,
    /// and if not, build one and save it or otherwise read in the index.
    ///
    /// The index building process is usually neglible on "regular" IO file
    /// systems.
    fn open_path<P>(path: P) -> io::Result<Self>
    where
        P: Into<path::PathBuf> + Clone,
    {
        let source: FileSource<fs::File> = FileSource::from(path.clone());
        let index_file_name = source.index_file_name();

        match fs::File::open(path.into()) {
            Ok(file) => {
                let mut reader = Self::open_file(file)?;
                if let Some(index_path) = &index_file_name {
                    if index_path.exists() {
                        let index_stream = fs::File::open(index_path)?;
                        match reader.read_index(Box::new(io::BufReader::new(index_stream))) {
                            Ok(_) => {}
                            Err(_err) => {
                                reader.construct_index_from_stream();
                            }
                        }
                    } else {
                        reader.construct_index_from_stream();
                    }
                }
                Ok(reader)
            }
            Err(err) => Err(err),
        }
    }

    /// Given a regular file, construct a new instance without indexing.
    fn open_file(source: fs::File) -> io::Result<Self>;
}

fn _save_index<
    C: CentroidLike + Default,
    D: DeconvolutedCentroidLike + Default,
    S: SpectrumLike<C, D>,
>(
    index_path: &PathBuf,
    reader: &impl MZFileReader<C, D, S>,
) -> io::Result<()> {
    let index_stream = fs::File::create(index_path)?;
    match reader.write_index(Box::new(io::BufWriter::new(index_stream))) {
        Ok(_) => {}
        Err(err) => {
            warn!(
                "Failed to write index to {} because {:?}",
                index_path.display(),
                err
            );
        }
    }
    Ok(())
}

/// Errors that may occur when reading a spectrum from a [`RandomAccessSpectrumIterator`]
#[derive(Debug, Error)]
pub enum SpectrumAccessError {
    /// An undetermined error failing to locate the requested spectrum
    #[error("The requested spectrum was not found")]
    SpectrumNotFound,
    /// An error resolving a spectrum by it's native ID
    #[error("The requested spectrum native ID {0} was not found")]
    SpectrumIdNotFound(String),
    /// An error resolving a spectrum by it's index
    #[error("The requested spectrum index {0} was not found")]
    SpectrumIndexNotFound(usize),
    /// An I/O error prevented reading the spectrum, even if it could be found.
    #[error("I/O error occurred while reading: {0:?}")]
    IOError(#[source] Option<io::Error>),
}

impl From<SpectrumAccessError> for io::Error {
    fn from(value: SpectrumAccessError) -> Self {
        let s = value.to_string();
        match value {
            SpectrumAccessError::SpectrumNotFound => io::Error::new(io::ErrorKind::NotFound, s),
            SpectrumAccessError::SpectrumIdNotFound(_) => {
                io::Error::new(io::ErrorKind::NotFound, s)
            }
            SpectrumAccessError::SpectrumIndexNotFound(_) => {
                io::Error::new(io::ErrorKind::NotFound, s)
            }
            SpectrumAccessError::IOError(e) => match e {
                Some(e) => e,
                None => io::Error::new(io::ErrorKind::Other, s),
            },
        }
    }
}

/// An extension of [`SpectrumSource`] that supports relocatable iteration relative to a
/// specific spectrum coordinate or identifier.
pub trait RandomAccessSpectrumIterator<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumLike<C, D> = MultiLayerSpectrum<C, D>,
>: SpectrumSource<C, D, S>
{
    /// Start iterating from the spectrum whose native ID matches `id`
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, SpectrumAccessError>;

    /// Start iterating from the spectrum whose index is `index`
    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, SpectrumAccessError>;

    /// Start iterating from the spectrum starting closest to `time`
    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, SpectrumAccessError>;
}

impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        R: SpectrumSource<C, D, S>,
    > RandomAccessSpectrumIterator<C, D, S> for SpectrumIterator<'lifespan, C, D, S, R>
{
    /// Start iterating from the spectrum whose native ID matches `id`
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, SpectrumAccessError> {
        if let Some(scan) = self.get_spectrum_by_id(id) {
            self.index = scan.index();
            self.back_index = 0;
            Ok(self)
        } else if self.get_index().contains_key(id) {
            Err(SpectrumAccessError::IOError(None))
        } else {
            Err(SpectrumAccessError::SpectrumIdNotFound(id.to_string()))
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, SpectrumAccessError> {
        if index < self.len() {
            self.index = index;
            self.back_index = 0;
            Ok(self)
        } else {
            Err(SpectrumAccessError::SpectrumIndexNotFound(index))
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, SpectrumAccessError> {
        if let Some(scan) = self.get_spectrum_by_time(time) {
            self.index = scan.index();
            self.back_index = 0;
            Ok(self)
        } else if self
            .get_spectrum_by_index(self.len() - 1)
            .expect("Failed to fetch spectrum for boundary testing")
            .start_time()
            < time
        {
            Err(SpectrumAccessError::SpectrumNotFound)
        } else {
            Err(SpectrumAccessError::IOError(None))
        }
    }
}

/// A union trait for [`SpectrumSource`] and [`RandomAccessSpectrumIterator`]
pub trait RandomAccessSpectrumSource<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumLike<C, D> = MultiLayerSpectrum<C, D>,
>: SpectrumSource<C, D, S> + RandomAccessSpectrumIterator<C, D, S>
{
}

impl<
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        T: SpectrumSource<C, D, S> + RandomAccessSpectrumIterator<C, D, S>,
    > RandomAccessSpectrumSource<C, D, S> for T
{
}

pub trait SpectrumSourceWithMetadata<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumLike<C, D> = MultiLayerSpectrum<C, D>,
>: SpectrumSource<C, D, S> + MSDataFileMetadata
{
}

impl<
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        T: SpectrumSource<C, D, S> + MSDataFileMetadata,
    > SpectrumSourceWithMetadata<C, D, S> for T
{
}

/// An alternative implementation of [`SpectrumSource`] for non-rewindable underlying iterators.
///
/// When the source doesn't support [`io::Seek`](std::io::Seek), most reader types don't
/// even implement [`SpectrumSource`], although they still implement
/// [`Iterator`]. The [`StreamingSpectrumIterator`]
/// wrapper does implement parts of [`SpectrumSource`] using less efficient
/// mechanism, but in situations where it cannot satisfy the request, it will `panic` instead. It also,
/// naturally doesn't support reading spectra that have already been seen as the stream cannot be reversed.
pub struct StreamingSpectrumIterator<
    C: CentroidLike + Default,
    D: DeconvolutedCentroidLike + Default,
    S: SpectrumLike<C, D>,
    I: Iterator<Item = S>,
> {
    source: I,
    buffer: VecDeque<S>,
    _index: OffsetIndex,
    _c: PhantomData<C>,
    _d: PhantomData<D>,
}

impl<
        C: CentroidLike + Default + Send,
        D: DeconvolutedCentroidLike + Default + Send,
        S: SpectrumLike<C, D> + Send,
    > From<SpectrumReceiver<C, D, S>>
    for StreamingSpectrumIterator<C, D, S, SpectrumReceiver<C, D, S>>
{
    fn from(value: SpectrumReceiver<C, D, S>) -> Self {
        Self::new(value)
    }
}

impl<
        C: CentroidLike + Default + Send,
        D: DeconvolutedCentroidLike + Default + Send,
        S: SpectrumLike<C, D> + Send,
    > From<Receiver<S>> for StreamingSpectrumIterator<C, D, S, SpectrumReceiver<C, D, S>>
{
    fn from(value: Receiver<S>) -> Self {
        Self::new(value.into())
    }
}

impl<
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        I: Iterator<Item = S>,
    > SpectrumSource<C, D, S> for StreamingSpectrumIterator<C, D, S, I>
{
    fn reset(&mut self) {
        panic!("Cannot reset StreamingSpectrumIterator")
    }

    fn iter(&mut self) -> SpectrumIterator<C, D, S, Self>
    where
        Self: Sized,
    {
        panic!(
            "Cannot create a wrapping iterator for StreamingSpectrumIterator, just use it directly"
        )
    }

    fn get_spectrum_by_id(&mut self, id: &str) -> Option<S> {
        self.find(|s| s.id() == id)
    }

    fn get_spectrum_by_index(&mut self, index: usize) -> Option<S> {
        self.find(|s| s.index() == index)
    }

    fn get_spectrum_by_time(&mut self, time: f64) -> Option<S> {
        let mut placeholder: Option<S> = None;
        let mut delta = f64::INFINITY;
        while let Some(s) = self.next() {
            let new_delta = (s.start_time() - time).abs();
            if s.start_time() < time {
                placeholder = Some(s);
                delta = new_delta;
            } else if s.start_time() >= time {
                if new_delta < delta {
                    return Some(s);
                } else {
                    self.push_front(s);
                    return placeholder;
                }
            }
        }
        None
    }

    fn get_index(&self) -> &OffsetIndex {
        &self._index
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self._index = index
    }
}

impl<
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        I: Iterator<Item = S>,
    > Iterator for StreamingSpectrumIterator<C, D, S, I>
{
    type Item = S;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.buffer.is_empty() {
            self.buffer.pop_front()
        } else {
            self.source.next()
        }
    }
}

impl<
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        I: Iterator<Item = S>,
    > StreamingSpectrumIterator<C, D, S, I>
{
    pub fn new(source: I) -> Self {
        Self {
            source,
            buffer: VecDeque::new(),
            _index: OffsetIndex::new("spectrum".to_string()),
            _c: PhantomData,
            _d: PhantomData,
        }
    }

    fn push_front(&mut self, spectrum: S) {
        self.buffer.push_front(spectrum);
    }
}

/// [`StreamingSpectrumIterator`] implements [`RandomAccessSpectrumIterator`] in a limited fashion
/// by reading through successive spectra until the target spectrum is found. This will exhaust the
/// underlying iterator if the requested coordinate is not found.
impl<
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        I: Iterator<Item = S>,
    > RandomAccessSpectrumIterator<C, D, S> for StreamingSpectrumIterator<C, D, S, I>
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, SpectrumAccessError> {
        match self.get_spectrum_by_id(id) {
            Some(s) => {
                self.push_front(s);
                Ok(self)
            }
            None => Err(SpectrumAccessError::SpectrumIdNotFound(id.to_string())),
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, SpectrumAccessError> {
        match self.get_spectrum_by_index(index) {
            Some(s) => {
                self.push_front(s);
                Ok(self)
            }
            None => Err(SpectrumAccessError::SpectrumIndexNotFound(index)),
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, SpectrumAccessError> {
        match self.get_spectrum_by_time(time) {
            Some(s) => {
                self.push_front(s);
                Ok(self)
            }
            None => Err(SpectrumAccessError::SpectrumNotFound),
        }
    }
}

/// If the underlying iterator implements [`MSDataFileMetadata`] then [`StreamingSpectrumIterator`] will
/// forward that implementation, assuming it is available.
impl<
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        I: Iterator<Item = S>,
    > MSDataFileMetadata for StreamingSpectrumIterator<C, D, S, I>
where
    I: MSDataFileMetadata,
{
    crate::delegate_impl_metadata_trait!(source);
}

/// An in-memory communication, non-rewindable channel carrying spectra
/// with associated metadata.
///
/// This type is meant to be wrapped in a [`StreamingSpectrumIterator`] for
/// compatibility with other interfaces.
pub struct SpectrumReceiver<
    C: CentroidLike + Default + Send,
    D: DeconvolutedCentroidLike + Default + Send,
    S: SpectrumLike<C, D> + Send,
> {
    receiver: Receiver<S>,

    pub(crate) file_description: FileDescription,
    /// A mapping of different instrument configurations (source, analyzer, detector) components
    /// by ID string.
    pub(crate) instrument_configurations: HashMap<u32, InstrumentConfiguration>,
    /// The different software components that were involved in the processing and creation of this
    /// file.
    pub(crate) softwares: Vec<Software>,
    /// The data processing and signal transformation operations performed on the raw data in previous
    /// source files to produce this file's contents.
    pub(crate) data_processings: Vec<DataProcessing>,
    // SpectrumList attributes
    pub(crate) run: MassSpectrometryRun,
    num_spectra: Option<u64>,
    _c: PhantomData<C>,
    _d: PhantomData<D>,
}

impl<
        C: CentroidLike + Default + Send,
        D: DeconvolutedCentroidLike + Default + Send,
        S: SpectrumLike<C, D> + Send,
    > Iterator for SpectrumReceiver<C, D, S>
{
    type Item = S;

    fn next(&mut self) -> Option<Self::Item> {
        match self.receiver.recv() {
            Ok(s) => Some(s),
            Err(e) => {
                log::warn!("Failed to receive spectrum: {}", e);
                None
            }
        }
    }
}

impl<
        C: CentroidLike + Default + Send,
        D: DeconvolutedCentroidLike + Default + Send,
        S: SpectrumLike<C, D> + Send,
    > From<Receiver<S>> for SpectrumReceiver<C, D, S>
{
    fn from(value: Receiver<S>) -> Self {
        Self {
            receiver: value,
            file_description: Default::default(),
            instrument_configurations: Default::default(),
            softwares: Default::default(),
            data_processings: Default::default(),
            run: Default::default(),
            num_spectra: Default::default(),
            _c: PhantomData,
            _d: PhantomData,
        }
    }
}

impl<
        C: CentroidLike + Default + Send,
        D: DeconvolutedCentroidLike + Default + Send,
        S: SpectrumLike<C, D> + Send,
    > MSDataFileMetadata for SpectrumReceiver<C, D, S>
{
    crate::impl_metadata_trait!();

    fn spectrum_count_hint(&self) -> Option<u64> {
        self.num_spectra
    }

    fn run_description(&self) -> Option<&MassSpectrometryRun> {
        Some(&self.run)
    }

    fn run_description_mut(&mut self) -> Option<&mut MassSpectrometryRun> {
        Some(&mut self.run)
    }
}

impl<
        C: CentroidLike + Default + Send,
        D: DeconvolutedCentroidLike + Default + Send,
        S: SpectrumLike<C, D> + Send,
    > SpectrumReceiver<C, D, S>
{
    #[allow(unused)]
    pub fn new(
        receiver: Receiver<S>,
        file_description: FileDescription,
        instrument_configurations: HashMap<u32, InstrumentConfiguration>,
        softwares: Vec<Software>,
        data_processings: Vec<DataProcessing>,
        run: MassSpectrometryRun,
        num_spectra: Option<u64>,
    ) -> Self {
        Self {
            receiver,
            file_description,
            instrument_configurations,
            softwares,
            data_processings,
            run,
            num_spectra,
            _c: PhantomData,
            _d: PhantomData,
        }
    }
}

/// An abstraction over [`SpectrumGroup`](crate::spectrum::SpectrumGroup)'s interface.
pub trait SpectrumGrouping<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumLike<C, D> = MultiLayerSpectrum<C, D>,
>: Default
{
    /// Get the precursor spectrum, which may be absent
    fn precursor(&self) -> Option<&S>;
    /// Get a mutable reference to the precursor spectrum, which may be absent
    fn precursor_mut(&mut self) -> Option<&mut S>;
    /// Explicitly set the precursor spectrum directly.
    fn set_precursor(&mut self, prec: S);

    /// Get a reference to the collection of product spectra
    fn products(&self) -> &[S];

    /// Get a mutable reference to the collection of product spectra
    fn products_mut(&mut self) -> &mut Vec<S>;

    /// The total number of spectra in the group
    fn total_spectra(&self) -> usize {
        self.precursor().is_some() as usize + self.products().len()
    }

    /// The spectrum that occurred first chronologically
    fn earliest_spectrum(&self) -> Option<&S> {
        self.precursor().or_else(|| {
            self.products().iter().min_by(|a, b| {
                a.acquisition()
                    .start_time()
                    .total_cmp(&b.acquisition().start_time())
            })
        })
    }

    /// The spectrum that occurred last chronologically
    fn latest_spectrum(&self) -> Option<&S> {
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

/// Analogous to to [`RandomAccessSpectrumIterator`], but for [`SpectrumGrouping`] implementations.
pub trait RandomAccessSpectrumGroupingIterator<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumLike<C, D> = MultiLayerSpectrum<C, D>,
    G: SpectrumGrouping<C, D, S> = SpectrumGroup<C, D, S>,
>: Iterator<Item = G>
{
    fn start_from_id(&mut self, id: &str) -> Result<&Self, SpectrumAccessError>;
    fn start_from_index(&mut self, index: usize) -> Result<&Self, SpectrumAccessError>;
    fn start_from_time(&mut self, time: f64) -> Result<&Self, SpectrumAccessError>;
    fn reset_state(&mut self);
}

/// A collection of spectra held in memory but providing an interface
/// identical to a data file. This structure owns its data, so in order
/// to yield ownership for [`SpectrumSource`], they are cloned
#[derive(Debug, Default)]
pub struct MemorySpectrumSource<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumLike<C, D> = MultiLayerSpectrum<C, D>,
> {
    spectra: VecDeque<S>,
    position: usize,
    offsets: OffsetIndex,
    _c: PhantomData<C>,
    _d: PhantomData<D>,
}

impl<
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D> + Clone,
    > MemorySpectrumSource<C, D, S>
{
    pub fn new(spectra: VecDeque<S>) -> Self {
        let mut offsets = OffsetIndex::new("spectrum".to_string());
        spectra.iter().enumerate().for_each(|(i, s)| {
            offsets.insert(s.id().to_string(), i as u64);
        });

        Self {
            spectra,
            position: 0,
            offsets,
            _c: PhantomData,
            _d: PhantomData,
        }
    }
}

impl<
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D> + Clone,
    > Iterator for MemorySpectrumSource<C, D, S>
{
    type Item = S;

    fn next(&mut self) -> Option<Self::Item> {
        if self.position < self.spectra.len() {
            let idx = self.position;
            self.position += 1;
            let value = self.spectra.index(idx);
            Some(value.clone())
        } else {
            None
        }
    }
}

impl<
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D> + Clone,
    > SpectrumSource<C, D, S> for MemorySpectrumSource<C, D, S>
{
    fn reset(&mut self) {
        self.position = 0;
    }

    fn get_spectrum_by_id(&mut self, id: &str) -> Option<S> {
        self.offsets.get(id).map(|i| {
            let value = &self.spectra[i as usize];
            value.clone()
        })
    }

    fn get_spectrum_by_index(&mut self, index: usize) -> Option<S> {
        if index < self.len() {
            Some(self.spectra.index(index).clone())
        } else {
            None
        }
    }

    fn get_index(&self) -> &OffsetIndex {
        &self.offsets
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.offsets = index
    }
}

impl<
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D> + Clone,
    > RandomAccessSpectrumIterator<C, D, S> for MemorySpectrumSource<C, D, S>
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, SpectrumAccessError> {
        match self.offsets.get(id) {
            Some(offset) => {
                self.position = offset as usize;
                Ok(self)
            }
            None => Err(SpectrumAccessError::SpectrumNotFound),
        }
    }

    fn start_from_index(&mut self, id: usize) -> Result<&mut Self, SpectrumAccessError> {
        match self.offsets.get_index(id) {
            Some((_, offset)) => {
                self.position = offset as usize;
                Ok(self)
            }
            None => Err(SpectrumAccessError::SpectrumNotFound),
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, SpectrumAccessError> {
        if let Some(scan) = self.get_spectrum_by_time(time) {
            self.position = scan.index();
            Ok(self)
        } else {
            Err(SpectrumAccessError::SpectrumNotFound)
        }
    }
}

impl<
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D> + Clone,
    > From<VecDeque<S>> for MemorySpectrumSource<C, D, S>
{
    fn from(value: VecDeque<S>) -> Self {
        Self::new(value)
    }
}

/// Common interface for spectrum writing
pub trait SpectrumWriter<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
>
{
    /// Write out a single spectrum
    fn write<S: SpectrumLike<C, D> + 'static>(&mut self, spectrum: &S) -> io::Result<usize>;

    /// Write out a single owned spectrum.
    ///
    /// This may produce fewer copies for some implementations.
    fn write_owned<S: SpectrumLike<C, D> + 'static>(&mut self, spectrum: S) -> io::Result<usize> {
        self.write(&spectrum)
    }

    /// As [`std::io::Write::flush`]
    fn flush(&mut self) -> io::Result<()>;

    /// Consume an [`Iterator`] over [`MultiLayerSpectrum`] references
    fn write_all<'b, S: SpectrumLike<C, D> + 'static, T: Iterator<Item = &'b S>>(
        &mut self,
        iterator: T,
    ) -> io::Result<usize> {
        let mut n = 0;
        for spectrum in iterator {
            n += self.write(spectrum)?;
        }
        Ok(n)
    }

    /// Consume an [`Iterator`] over [`MultiLayerSpectrum`]
    fn write_all_owned<'b, S: SpectrumLike<C, D> + 'static, T: Iterator<Item = S>>(
        &mut self,
        iterator: T,
    ) -> io::Result<usize> {
        let mut n = 0;
        for spectrum in iterator {
            n += self.write_owned(spectrum)?;
        }
        Ok(n)
    }

    /// Write a [`SpectrumGroup`] out in order
    fn write_group<S: SpectrumLike<C, D> + 'static, G: SpectrumGrouping<C, D, S> + 'static>(
        &mut self,
        group: &G,
    ) -> io::Result<usize> {
        let mut n = 0;
        if let Some(precursor) = group.precursor() {
            n += self.write(precursor)?;
        }
        for product in group.products() {
            n += self.write(product)?;
        }
        Ok(n)
    }

    /// Write an owned [`SpectrumGroup`] out in order
    ///
    /// This may produce fewer copies for some implementations.
    fn write_group_owned<
        S: SpectrumLike<C, D> + 'static,
        G: SpectrumGrouping<C, D, S> + 'static,
    >(
        &mut self,
        group: G,
    ) -> io::Result<usize> {
        let (precursor, products) = group.into_parts();
        let mut n = 0;
        if let Some(precursor) = precursor {
            n += self.write_owned(precursor)?;
        }
        for product in products {
            n += self.write_owned(product)?;
        }
        Ok(n)
    }

    /// Consume an [`Iterator`] over [`SpectrumGroup`] references
    fn write_all_groups<
        'b,
        S: SpectrumLike<C, D> + 'static,
        G: SpectrumGrouping<C, D, S> + 'static,
        T: Iterator<Item = &'b G>,
    >(
        &mut self,
        iterator: T,
    ) -> io::Result<usize> {
        let mut n = 0;
        for group in iterator {
            n += self.write_group(group)?;
        }
        Ok(n)
    }

    /// Consume an [`Iterator`] over [`SpectrumGroup`]
    fn write_all_groups_owned<
        'b,
        S: SpectrumLike<C, D> + 'static,
        G: SpectrumGrouping<C, D, S> + 'static,
        T: Iterator<Item = G>,
    >(
        &mut self,
        iterator: T,
    ) -> io::Result<usize> {
        let mut n = 0;
        for group in iterator {
            n += self.write_group_owned(group)?;
        }
        Ok(n)
    }

    /// Completes the data file format, preventing new data from being able incorporate additional
    /// data. Does not formally close the underlying writing stream.
    fn close(&mut self) -> io::Result<()>;
}

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

    /// Create a new `SpectrumIterator` over `self` and use that state to drive a `SpectrumGroupIterator`
    fn groups(&mut self) -> IonMobilityFrameGroupingIterator<IonMobilityFrameIterator<'_, C, D, S, Self>, C, D, S>
    where
        Self: Sized,
    {
        IonMobilityFrameGroupingIterator::new(self.iter())
    }

    /// Consume `self` to create a `SpectrumGroupIterator`. This is ideal for non-rewindable streams
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

/// If the underlying iterator implements [`MSDataFileMetadata`] then [`SpectrumIterator`] will
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

/// A trait that for retrieving [`Chromatogram`]s from a source.
pub trait ChromatogramSource {
    /// Get a [`Chromatogram`] by its identifier, if it exists.
    fn get_chromatogram_by_id(&mut self, id: &str) -> Option<Chromatogram>;

    /// Get a [`Chromatogram`] by its index, if it exists.
    fn get_chromatogram_by_index(&mut self, index: usize) -> Option<Chromatogram>;

    /// Iterate over [`Chromatogram`]s with a [`ChromatogramIter`]
    fn iter_chromatograms(&mut self) -> ChromatogramIterator<'_, Self>
    where
        Self: Sized,
    {
        ChromatogramIterator::new(self)
    }
}

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

impl<'a, R: ChromatogramSource> Iterator for ChromatogramIterator<'a, R> {
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

impl<'a, R: ChromatogramSource> FusedIterator for ChromatogramIterator<'a, R> {}

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
