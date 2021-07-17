use std::io;

use indexmap::IndexMap;

use crate::spectrum::SpectrumBehavior;

pub trait ScanIterator<S: SpectrumBehavior>: Iterator<Item = S> {}

pub trait ScanSource<S: SpectrumBehavior>: ScanIterator<S> {
    fn reset(&mut self) -> &Self;

    /// Retrieve the number of spectra in source file
    fn len(&self) -> usize {
        self.get_index().len()
    }

    fn _offset_of_id(&self, id: &str) -> Option<u64> {
        self.get_index().get(id).map(|offset| *offset)
    }

    fn _offset_of_index(&self, index: usize) -> Option<u64> {
        self.get_index()
            .get_index(index)
            .map(|(_id, offset)| *offset)
    }

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

    fn get_index(&self) -> &IndexMap<String, u64>;

    /// Retrieve a spectrum by its scan start time
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

#[derive(Debug)]
pub enum ScanAccessError {
    ScanNotFound,
    IOError(io::Error),
}

pub trait RandomAccessScanIterator<S: SpectrumBehavior>: ScanSource<S> {
    fn start_from_id(&mut self, id: &str) -> Result<&Self, ScanAccessError>;
    fn start_from_index(&mut self, id: usize) -> Result<&Self, ScanAccessError>;
    fn start_from_time(&mut self, time: f64) -> Result<&Self, ScanAccessError>;
}
