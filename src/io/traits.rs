use std::fs;
use std::io;
use std::marker::PhantomData;
use std::path;

use log::warn;

use crate::spectrum::SpectrumBehavior;

use super::utils::FileSource;
use super::OffsetIndex;

pub trait SeekRead: io::Read + io::Seek {}
impl<T: io::Read + io::Seek> SeekRead for T {}

pub trait ScanSource<S: SpectrumBehavior>: Iterator<Item = S> + Sized {
    fn reset(&mut self) -> &Self;

    /// Retrieve the number of spectra in source file
    fn len(&self) -> usize {
        self.get_index().len()
    }

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

    /// Retrieve a spectrum by it's native ID
    fn get_spectrum_by_id(&mut self, id: &str) -> Option<S>;

    /// Retrieve a spectrum by it's integer index
    fn get_spectrum_by_index(&mut self, index: usize) -> Option<S>;

    fn get_index(&self) -> &OffsetIndex;

    fn set_index(&mut self, index: OffsetIndex);

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

    /// Open a new iterator over this stream.
    fn iter(&mut self) -> ScanIterator<Self, S> {
        ScanIterator::new(self)
    }
}

pub trait MZFileReader<S: SpectrumBehavior>: ScanSource<S> + Sized {
    fn construct_index_from_stream(&mut self) -> u64;

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
                        let index_stream = fs::File::open(index_path).unwrap();
                        match reader.read_index(&index_stream) {
                            Ok(_) => {}
                            Err(_err) => {
                                reader.construct_index_from_stream();
                                match reader.write_index(index_stream) {
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
                        let index_stream = fs::File::create(index_path).unwrap();
                        match reader.write_index(index_stream) {
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

    fn open_file(source: fs::File) -> Self;
}

pub struct ScanIterator<'lifespan, R: ScanSource<S>, S: SpectrumBehavior> {
    source: &'lifespan mut R,
    phantom: PhantomData<S>,
    index: usize,
    back_index: usize,
}

impl<'lifespan, R: ScanSource<S>, S: SpectrumBehavior> ScanIterator<'lifespan, R, S> {
    pub fn new(source: &mut R) -> ScanIterator<R, S> {
        ScanIterator {
            source,
            index: 0,
            back_index: 0,
            phantom: PhantomData,
        }
    }
}

impl<'lifespan, R: ScanSource<S>, S: SpectrumBehavior> Iterator for ScanIterator<'lifespan, R, S> {
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

impl<'lifespan, R: ScanSource<S>, S: SpectrumBehavior> ExactSizeIterator
    for ScanIterator<'lifespan, R, S>
{
    fn len(&self) -> usize {
        self.source.len()
    }
}

impl<'lifespan, R: ScanSource<S>, S: SpectrumBehavior> DoubleEndedIterator
    for ScanIterator<'lifespan, R, S>
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

impl<'lifespan, R: ScanSource<S>, S: SpectrumBehavior> ScanSource<S>
    for ScanIterator<'lifespan, R, S>
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

#[derive(Debug)]
pub enum ScanAccessError {
    ScanNotFound,
    IOError(Option<io::Error>),
}

pub trait RandomAccessScanIterator<S: SpectrumBehavior>: ScanSource<S> {
    fn start_from_id(&mut self, id: &str) -> Result<&Self, ScanAccessError>;
    fn start_from_index(&mut self, id: usize) -> Result<&Self, ScanAccessError>;
    fn start_from_time(&mut self, time: f64) -> Result<&Self, ScanAccessError>;
}

impl<'lifespan, R: ScanSource<S>, S: SpectrumBehavior> RandomAccessScanIterator<S>
    for ScanIterator<'lifespan, R, S>
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
