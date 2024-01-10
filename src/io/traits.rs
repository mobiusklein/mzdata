use log::warn;
use std::collections::VecDeque;
use std::fs;
use std::io;
use std::marker::PhantomData;
use std::ops::Index;
use std::path::{self, PathBuf};
use thiserror::Error;

use mzpeaks::{CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak};

use crate::prelude::MSDataFileMetadata;
use crate::spectrum::group::SpectrumGroupingIterator;
use crate::spectrum::spectrum::{MultiLayerSpectrum, SpectrumLike};
use crate::spectrum::SpectrumGroup;

use super::utils::FileSource;
use super::OffsetIndex;

pub trait SeekRead: io::Read + io::Seek {}
impl<T: io::Read + io::Seek> SeekRead for T {}

/// A base trait defining the behaviors of a source of spectra.
///
/// A [`ScanSource`]
pub trait ScanSource<
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
            } else if (scan_time - time).abs() < 1e-3 {
                return Some(scan);
            } else if scan_time > time {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        best_match
    }

    /// Retrieve the number of spectra in source file
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

/// A generic iterator over a [`ScanSource`] implementer that assumes the
/// source has already been indexed. Otherwise, the source's own iterator
/// behavior should be used.
pub struct SpectrumIterator<
    'lifespan,
    C: CentroidLike + Default,
    D: DeconvolutedCentroidLike + Default,
    S: SpectrumLike<C, D>,
    R: ScanSource<C, D, S>,
> {
    source: &'lifespan mut R,
    phantom: PhantomData<S>,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
    index: usize,
    back_index: usize,
}

impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        R: ScanSource<C, D, S>,
        S: SpectrumLike<C, D>,
    > SpectrumIterator<'lifespan, C, D, S, R>
{
    pub fn new(source: &mut R) -> SpectrumIterator<C, D, S, R> {
        SpectrumIterator::<C, D, S, R> {
            source,
            index: 0,
            back_index: 0,
            phantom: PhantomData,
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
        }
    }
}

impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        R: ScanSource<C, D, S>,
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
        R: ScanSource<C, D, S>,
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
        R: ScanSource<C, D, S>,
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
        R: ScanSource<C, D, S>,
    > ScanSource<C, D, S> for SpectrumIterator<'lifespan, C, D, S, R>
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

/// A trait defining some helper methods to make efficient use of indices
/// automatic when opening a file from a path-like object.
pub trait MZFileReader<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumLike<C, D> = MultiLayerSpectrum<C, D>,
>: ScanSource<C, D, S> + Sized
{
    /// An on-trait method of constructing an index. Assumed
    /// to be a trivial wrapper.
    fn construct_index_from_stream(&mut self) -> u64;

    /// Re-construct an offset index from this readable object, assuming
    /// it is a JSON stream over the serialized index.
    fn read_index(&mut self, reader: Box<dyn io::Read>) -> Result<&Self, serde_json::Error> {
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
                let mut reader = Self::open_file(file);
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
    fn open_file(source: fs::File) -> Self;
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
            SpectrumAccessError::SpectrumIdNotFound(_) => io::Error::new(io::ErrorKind::NotFound, s),
            SpectrumAccessError::SpectrumIndexNotFound(_) => io::Error::new(io::ErrorKind::NotFound, s),
            SpectrumAccessError::IOError(e) => {
                match e {
                    Some(e) => {
                        e
                    },
                    None => {
                        io::Error::new(io::ErrorKind::Other, s)
                    },
                }
            },
        }
    }
}

/// An extension of [`ScanSource`] that supports relocatable iteration relative to a
/// specific spectrum coordinate or identifier.
pub trait RandomAccessSpectrumIterator<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumLike<C, D> = MultiLayerSpectrum<C, D>,
>: ScanSource<C, D, S>
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
        R: ScanSource<C, D, S>,
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

/// An alternative implementation of [`ScanSource`] for non-rewindable underlying streams
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
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        I: Iterator<Item = S>,
    > ScanSource<C, D, S> for StreamingSpectrumIterator<C, D, S, I>
{
    fn reset(&mut self) {
        panic!("Cannot reset StreamingSpectrumIterator")
    }

    fn get_spectrum_by_id(&mut self, id: &str) -> Option<S> {
        while let Some(s) = self.next() {
            if s.id() == id {
                return Some(s);
            }
        }
        None
    }

    fn get_spectrum_by_index(&mut self, index: usize) -> Option<S> {
        while let Some(s) = self.next() {
            if s.index() == index {
                return Some(s);
            }
        }
        None
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

impl<
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        I: Iterator<Item = S>,
    > MSDataFileMetadata for StreamingSpectrumIterator<C, D, S, I> where I: MSDataFileMetadata {

    fn data_processings(&self) -> &Vec<crate::meta::DataProcessing> {
        self.source.data_processings()
    }

    fn instrument_configurations(&self) -> &std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
        self.source.instrument_configurations()
    }

    fn file_description(&self) -> &crate::meta::FileDescription {
        self.source.file_description()
    }

    fn softwares(&self) -> &Vec<crate::meta::Software> {
        self.source.softwares()
    }

    fn data_processings_mut(&mut self) -> &mut Vec<crate::meta::DataProcessing> {
        self.source.data_processings_mut()
    }

    fn instrument_configurations_mut(&mut self) -> &mut std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
        self.source.instrument_configurations_mut()
    }

    fn file_description_mut(&mut self) -> &mut crate::meta::FileDescription {
        self.source.file_description_mut()
    }

    fn softwares_mut(&mut self) -> &mut Vec<crate::meta::Software> {
        self.source.softwares_mut()
    }

    fn spectrum_count_hint(&self) -> Option<u64> {
        self.source.spectrum_count_hint()
    }

    fn run_description(&self) -> Option<&crate::meta::MassSpectrometryRun> {
        self.source.run_description()
    }

    fn run_description_mut(&mut self) -> Option<&mut crate::meta::MassSpectrometryRun> {
        self.source.run_description_mut()
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
        let prec_level = self
            .precursor()
            .map(|p| p.ms_level())
            .unwrap_or(u8::MAX);
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
/// to yield ownership for [`ScanSource`], they are cloned
#[derive(Debug, Default)]
pub struct MemoryScanSource<
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
    > MemoryScanSource<C, D, S>
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
    > Iterator for MemoryScanSource<C, D, S>
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
    > ScanSource<C, D, S> for MemoryScanSource<C, D, S>
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
    > RandomAccessSpectrumIterator<C, D, S> for MemoryScanSource<C, D, S>
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
    > From<VecDeque<S>> for MemoryScanSource<C, D, S>
{
    fn from(value: VecDeque<S>) -> Self {
        Self::new(value)
    }
}

/// Common interface for spectrum writing
pub trait ScanWriter<
    'a,
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
>
{
    /// Write out a single spectrum
    fn write<S: SpectrumLike<C, D> + 'static>(&mut self, spectrum: &S) -> io::Result<usize>;

    /// As [`std::io::Write::flush`]
    fn flush(&mut self) -> io::Result<()>;

    /// Consume an [`Iterator`] over [`Spectrum`](crate::spectrum::MultiLayerSpectrum) references
    fn write_all<S: SpectrumLike<C, D> + 'static, T: Iterator<Item = &'a S>>(
        &mut self,
        iterator: T,
    ) -> io::Result<usize> {
        let mut n = 0;
        for spectrum in iterator {
            n += self.write(spectrum)?;
        }
        Ok(n)
    }

    /// Write a [`SpectrumGroup`](crate::spectrum::SpectrumGroup) out in order
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

    /// Consume an [`Iterator`] over [`SpectrumGroup`](crate::spectrum::SpectrumGroup) references
    fn write_all_groups<
        S: SpectrumLike<C, D> + 'static,
        G: SpectrumGrouping<C, D, S> + 'static,
        T: Iterator<Item = &'a G>,
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

    fn close(&mut self) -> io::Result<()>;
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_object_safe() {
        // If `ScanSource` were not object safe, this code
        // couldn't compile.
        let _f = |_x: &dyn ScanSource| {};
    }
}
