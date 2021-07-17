use crate::spectrum::{SpectrumBehavior};

pub trait ScanIterator<S: SpectrumBehavior> : Iterator<Item=S> {}

pub trait ScanSource<S: SpectrumBehavior> : ScanIterator<S> {
    fn reset(&mut self) -> &Self;
    fn get_spectrum_by_id(&mut self, id: &str) -> Option<S>;
    fn get_spectrum_by_index(&mut self, index: usize) -> Option<S>;
    fn len(&self) -> usize;

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
            } else if scan_time == time {
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
