use std::fs;
use std::io;
use std::marker::PhantomData;
use std::path;

use log::warn;

use mzpeaks::{CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak};

use crate::spectrum::spectrum::{
    MultiLayerSpectrum, SpectrumBehavior,
};

use super::utils::FileSource;
use super::OffsetIndex;

pub trait SeekRead: io::Read + io::Seek {}
impl<T: io::Read + io::Seek> SeekRead for T {}

trait IndexedScanSource<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumBehavior<C, D> = MultiLayerSpectrum<C, D>,
> {
    /// Retrieve the number of spectra in source file
    fn len(&self) -> usize {
        self.get_index().len()
    }

    fn get_index(&self) -> &OffsetIndex;

    fn set_index(&mut self, index: OffsetIndex);

    /// Helper method to support seeking to an ID
    fn _offset_of_id(&self, id: &str) -> Option<u64> {
        self.get_index().get(id).map(|offset| offset)
    }

    /// Helper method to support seeking to an index
    fn _offset_of_index(&self, index: usize) -> Option<u64> {
        self.get_index()
            .get_index(index)
            .map(|(_id, offset)| offset)
    }

    /// Helper method to support seeking to a specific time.
    /// Considerably more complex than seeking by ID or index.
    fn _offset_of_time(&mut self, time: f64) -> Option<u64>; /*{
                                                                 match self.get_spectrum_by_time(time) {
                                                                     Some(scan) => self._offset_of_index(scan.index()),
                                                                     None => None,
                                                                 }
                                                             }*/

    /// Re-construct an offset index from this readable object, assuming
    /// it is a JSON stream over the serialized index.
    fn read_index<R: io::Read>(&mut self, reader: R) -> Result<&Self, serde_json::Error> {
        match OffsetIndex::from_reader(reader) {
            Ok(index) => {
                self.set_index(index);
                Ok(self)
            }
            Err(err) => Err(err),
        }
    }

    fn write_index<W: io::Write>(&self, writer: W) -> Result<&Self, serde_json::Error> {
        match self.get_index().to_writer(writer) {
            Ok(_) => Ok(self),
            Err(err) => Err(err),
        }
    }
}

pub trait RandomAccessScanSource<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumBehavior<C, D> = MultiLayerSpectrum<C, D>,
> {
    fn len(&self) -> usize;

    /// Retrieve a spectrum by it's native ID
    fn get_spectrum_by_id(&mut self, id: &str) -> Option<S>;

    /// Retrieve a spectrum by it's integer index
    fn get_spectrum_by_index(&mut self, index: usize) -> Option<S>;

    /// Retrieve a spectrum by its scan start time
    /// Considerably more complex than seeking by ID or index.
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
}

/// A base trait defining the behaviors of a source of spectra.
pub trait ScanSource<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumBehavior<C, D> = MultiLayerSpectrum<C, D>,
>: Iterator<Item = S> + Sized
{
    fn reset(&mut self) -> &Self;

    /// Retrieve a spectrum by it's native ID
    fn get_spectrum_by_id(&mut self, id: &str) -> Option<S>;

    /// Retrieve a spectrum by it's integer index
    fn get_spectrum_by_index(&mut self, index: usize) -> Option<S>;

    /// Retrieve a spectrum by its scan start time
    /// Considerably more complex than seeking by ID or index.
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

    fn get_index(&self) -> &OffsetIndex;

    fn set_index(&mut self, index: OffsetIndex);

    /// Helper method to support seeking to an ID
    fn _offset_of_id(&self, id: &str) -> Option<u64> {
        self.get_index().get(id).map(|offset| offset)
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

    /// Re-construct an offset index from this readable object, assuming
    /// it is a JSON stream over the serialized index.
    fn read_index<R: io::Read>(&mut self, reader: R) -> Result<&Self, serde_json::Error> {
        match OffsetIndex::from_reader(reader) {
            Ok(index) => {
                self.set_index(index);
                Ok(self)
            }
            Err(err) => Err(err),
        }
    }

    fn write_index<W: io::Write>(&self, writer: W) -> Result<&Self, serde_json::Error> {
        match self.get_index().to_writer(writer) {
            Ok(_) => Ok(self),
            Err(err) => Err(err),
        }
    }

    /// Open a new iterator over this stream.
    fn iter(&mut self) -> ScanIterator<C, D, S, Self> {
        ScanIterator::new(self)
    }
}



/// A generic iterator over a [`ScanSource`] implementer that assumes the
/// source has already been indexed. Otherwise, the source's own iterator
/// behavior should be used.
pub struct ScanIterator<
    'lifespan,
    C: CentroidLike + Default,
    D: DeconvolutedCentroidLike + Default,
    S: SpectrumBehavior<C, D>,
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
        S: SpectrumBehavior<C, D>,
    > ScanIterator<'lifespan, C, D, S, R>
{
    pub fn new(source: &mut R) -> ScanIterator<C, D, S, R> {
        ScanIterator::<C, D, S, R> {
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
        S: SpectrumBehavior<C, D>,
    > Iterator for ScanIterator<'lifespan, C, D, S, R>
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

impl<'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        R: ScanSource<C, D, S>,
        S: SpectrumBehavior<C, D>,> ExactSizeIterator
    for ScanIterator<'lifespan, C, D, S, R>
{
    fn len(&self) -> usize {
        self.source.len()
    }
}

impl<'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        R: ScanSource<C, D, S>,
        S: SpectrumBehavior<C, D>> DoubleEndedIterator
    for ScanIterator<'lifespan, C, D, S, R>
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
        S: SpectrumBehavior<C, D>,
        R: ScanSource<C, D, S>,
    > ScanSource<C, D, S> for ScanIterator<'lifespan, C, D, S, R>
{
    fn reset(&mut self) -> &Self {
        self.index = 0;
        self.back_index = 0;
        self
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
    C: CentroidLike + Default,
    D: DeconvolutedCentroidLike + Default,
    S: SpectrumBehavior<C, D>,
>: ScanSource<C, D, S> + Sized
{
    /// An on-trait method of constructing an index. Assumed
    /// to be a trivial wrapper.
    fn construct_index_from_stream(&mut self) -> u64;

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

        match fs::File::open(path.clone().into()) {
            Ok(file) => {
                let mut reader = Self::open_file(file);
                if let Some(index_path) = &index_file_name {
                    if index_path.exists() {
                        let index_stream = fs::File::open(index_path).expect(&format!(
                            "Failed to open index file {}",
                            index_path.display()
                        ));
                        match reader.read_index(io::BufReader::new(&index_stream)) {
                            Ok(_) => {}
                            Err(_err) => {
                                reader.construct_index_from_stream();
                                match reader.write_index(io::BufWriter::new(index_stream)) {
                                    Ok(_) => {}
                                    Err(err) => {
                                        warn!(
                                            "Failed to write index to {} because {:?}",
                                            index_path.display(),
                                            err
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        reader.construct_index_from_stream();
                        let index_stream = fs::File::create(index_path).expect(&format!(
                            "Failed to create index file {}",
                            index_path.display()
                        ));
                        match reader.write_index(io::BufWriter::new(index_stream)) {
                            Ok(_) => {}
                            Err(err) => {
                                warn!(
                                    "Failed to write index to {} because {:?}",
                                    index_path.display(),
                                    err
                                );
                            }
                        }
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

#[derive(Debug)]
pub enum ScanAccessError {
    ScanNotFound,
    IOError(Option<io::Error>),
}

pub trait RandomAccessScanIterator<C: CentroidLike + Default=CentroidPeak, D: DeconvolutedCentroidLike + Default=DeconvolutedPeak, S: SpectrumBehavior<C,D>=MultiLayerSpectrum<C,D>>: ScanSource<C, D, S> {
    fn start_from_id(&mut self, id: &str) -> Result<&Self, ScanAccessError>;
    fn start_from_index(&mut self, id: usize) -> Result<&Self, ScanAccessError>;
    fn start_from_time(&mut self, time: f64) -> Result<&Self, ScanAccessError>;
}

impl<'lifespan, C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default, S: SpectrumBehavior<C,D>, R: ScanSource<C, D, S>> RandomAccessScanIterator<C, D, S>
    for ScanIterator<'lifespan, C, D, S, R>
{
    fn start_from_id(&mut self, id: &str) -> Result<&Self, ScanAccessError> {
        if let Some(scan) = self.get_spectrum_by_id(id) {
            self.index = scan.index();
            self.back_index = 0;
            Ok(self)
        } else {
            if self.get_index().contains_key(id) {
                Err(ScanAccessError::IOError(None))
            } else {
                Err(ScanAccessError::ScanNotFound)
            }
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&Self, ScanAccessError> {
        if index < self.len() {
            self.index = index;
            self.back_index = 0;
            Ok(self)
        } else {
            Err(ScanAccessError::ScanNotFound)
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&Self, ScanAccessError> {
        if let Some(scan) = self.get_spectrum_by_time(time) {
            self.index = scan.index();
            self.back_index = 0;
            Ok(self)
        } else {
            if self
                .get_spectrum_by_index(self.len() - 1)
                .expect("Failed to fetch spectrum for boundary testing")
                .start_time()
                < time
            {
                Err(ScanAccessError::ScanNotFound)
            } else {
                Err(ScanAccessError::IOError(None))
            }
        }
    }
}
