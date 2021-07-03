use std::fmt;
use std::ops;
use std::option::Option;

use crate::peak::{Peak};
use crate::mass_error::{MassErrorType};
use crate::coordinate::{CoordinateLike, MZ};


pub trait PeakCollection<T: CoordinateLike<C>, C> : ops::Index<usize>
    where <Self as ops::Index<usize>>::Output : CoordinateLike<C>
    {
    fn sort(&mut self);

    fn len(&self) -> usize;
    fn get_item(&self, i: usize) -> &T;
    fn get_slice(&self, i: ops::Range<usize>) -> &[T];

    fn _search_by(&self, query: f64) -> Result<usize, usize>;

    fn _closest_peak(&self, query: f64, error_tolerance: f64, i: usize, error_type: MassErrorType) -> Option<usize> {
        let mut j = i;
        let mut best = j;
        let mut best_err = error_type.call(
            self.get_item(j).get_coordinate(), query).abs();
        let n = self.len();
        while j < n {
            let err = error_type.call(
                self.get_item(j).get_coordinate(), query).abs();
            if err < best_err && err < error_tolerance{
                best_err = err;
                best = j;
            }
            j += 1;
        }
        if best_err > error_tolerance {
            return None
        }
        return Some(best);
    }

    fn search(&self, query: f64, error_tolerance: f64, error_type: MassErrorType) -> Option<usize> {
        let lower_bound = error_type.lower_bound(query, error_tolerance);
        let i = match self._search_by(lower_bound) {
            Ok(j) => self._closest_peak(query, error_tolerance, j, error_type),
            Err(j) => self._closest_peak(query, error_tolerance, j, error_type),
        };
        return i;
    }

    fn has_peak(&self, query: f64, error_tolerance: f64, error_type: MassErrorType) -> Option<&T> {
        return match self.search(query, error_tolerance, error_type) {
            Some(j) => Some(self.get_item(j)),
            None => None
        }
    }

    fn between(&self, low: f64, high: f64, error_tolerance: f64, error_type: MassErrorType) -> &[T] {
        let lower_bound = error_type.lower_bound(low, error_tolerance);
        let upper_bound = error_type.upper_bound(high, error_tolerance);

        let n = self.len();
        if n == 0 {
            return &self.get_slice(0..0);
        }

        let mut lower_index = match self._search_by(lower_bound) {
            Ok(j) => j,
            Err(j) => j,
        };

        let mut upper_index = match self._search_by(upper_bound) {
            Ok(j) => j,
            Err(j) => j,
        };

        if lower_index < n {
            if self[lower_index].get_coordinate() < lower_bound {
                lower_index += 1;
            }
        }

        if upper_index < n {
            if self[upper_index].get_coordinate() > upper_bound {
                upper_index -= 1;
            }
        }

        let subset = &self.get_slice(lower_index..upper_index + 1);
        return subset
    }

    fn all_peaks_for(&self, query: f64, error_tolerance: f64, error_type: MassErrorType) -> &[T] {
        let lower_bound = error_type.lower_bound(query, error_tolerance);
        let upper_bound = error_type.upper_bound(query, error_tolerance);

        let n = self.len();
        if n == 0 {
            return &self.get_slice(0..0);
        }

        let mut lower_index = match self._search_by(lower_bound) {
            Ok(j) => j,
            Err(j) => j,
        };
        if lower_index < n {
            if self[lower_index].get_coordinate() < lower_bound {
                lower_index += 1;
            }
        }

        let mut upper_index = lower_index;

        for i in lower_index + 1..n {
            if self[i].get_coordinate() >= upper_bound {
                break
            } else {
                upper_index = i;
            }
        }

        return &self.get_slice(lower_index..upper_index + 1);
    }
}


#[derive(Clone, Debug)]
pub struct PeakSet {
    pub peaks: Vec<Peak>,
}

impl PeakSet {

    pub fn new(mut peaks: Vec<Peak>) -> Self {
        Self::_sort(&mut peaks);
        let inst: PeakSet = PeakSet { peaks: peaks, };
        return inst;
    }

    pub fn from<V: Iterator>(peaks: V, sort: bool) -> PeakSet where V: Iterator<Item=Peak> {
        let peaks: Vec<Peak> = peaks.collect();
        if sort {
            return Self::new(peaks)
        } else {
            return Self::wrap(peaks)
        }
    }

    pub fn wrap(peaks: Vec<Peak>) -> PeakSet{
        let inst: PeakSet = PeakSet { peaks: peaks };
        return inst;
    }

    fn _sort(peaks: &mut Vec<Peak>) {
        peaks.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (i, p) in peaks.iter_mut().enumerate() {
            p.index = i as u32;
        }
    }
}

impl PeakCollection<Peak, MZ> for PeakSet {
    fn sort(&mut self) {
        Self::_sort(&mut self.peaks);
    }

    fn len(&self) -> usize {
        return self.peaks.len();
    }

    fn get_item(&self, i: usize) -> &Peak {
        return &self[i]
    }

    fn get_slice(&self, i: ops::Range<usize>) -> &[Peak] {
        return &self.peaks[i]
    }

    fn _search_by(&self, query: f64) -> Result<usize, usize> {
        self.peaks.binary_search_by(|peak| peak.get_coordinate().partial_cmp(&query).unwrap())
    }
}


impl fmt::Display for PeakSet {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "PeakSet([").expect("Write failed");
        let n = self.len() - 1;
        for (i, peak) in self.into_iter().enumerate() {
            write!(f, "{}", peak).expect("Write failed");
            if i != n {
                write!(f, ",").expect("Write failed");
            }
        }
        write!(f, "])").expect("Write failed");
        return Ok(())
    }
}

impl ops::Index<usize> for PeakSet {
    type Output = Peak;

    fn index(&self, i: usize) -> &Self::Output {
        return &(self.peaks[i]);
    }
}

impl<'a> IntoIterator for &'a PeakSet {
    type Item = &'a Peak;
    type IntoIter = PeakSetIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        return PeakSetIter::new(self);
    }
}

// Iterators

pub struct PeakSetIter<'a> {
    iter: std::slice::Iter<'a, Peak>
}

impl<'a> PeakSetIter<'a> {
    fn new(peaks: &'a PeakSet) -> PeakSetIter<'a> {
        return PeakSetIter {
            iter: peaks.peaks.iter()
        }
    }
}

impl<'a> Iterator for PeakSetIter<'a> {
    type Item = &'a Peak;

    fn next(&mut self) -> Option<Self::Item> {
        return self.iter.next()
    }
}
