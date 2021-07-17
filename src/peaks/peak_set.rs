use super::coordinate::{CoordinateLike, IndexedCoordinate, MZ};
use super::peak::CentroidPeak;
use crate::mass_error::MassErrorType;
use std::fmt;
use std::marker;
use std::ops;

/// A trait for an ordered container of mass spectral peaks. The trait
/// interoperates with [`CoordinateLike`] to make searching efficient.
pub trait PeakCollection<T: CoordinateLike<C>, C>: ops::Index<usize>
where
    <Self as ops::Index<usize>>::Output: CoordinateLike<C>,
{
    fn push(&mut self, peak: T);
    fn sort(&mut self);

    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    fn get_item(&self, i: usize) -> &T;
    fn get_slice(&self, i: ops::Range<usize>) -> &[T];

    fn _search_by(&self, query: f64) -> Result<usize, usize>;

    fn _closest_peak(
        &self,
        query: f64,
        error_tolerance: f64,
        i: usize,
        error_type: MassErrorType,
    ) -> Option<usize> {
        let mut j = i;
        let mut best = j;
        let mut best_err = error_type
            .call(self.get_item(j).get_coordinate(), query)
            .abs();
        let n = self.len();
        while j < n {
            let err = error_type
                .call(self.get_item(j).get_coordinate(), query)
                .abs();
            if err < best_err && err < error_tolerance {
                best_err = err;
                best = j;
            }
            j += 1;
        }
        if best_err > error_tolerance {
            return None;
        }
        Some(best)
    }

    fn search(&self, query: f64, error_tolerance: f64, error_type: MassErrorType) -> Option<usize> {
        let lower_bound = error_type.lower_bound(query, error_tolerance);
        match self._search_by(lower_bound) {
            Ok(j) => self._closest_peak(query, error_tolerance, j, error_type),
            Err(j) => self._closest_peak(query, error_tolerance, j, error_type),
        }
    }

    fn has_peak(&self, query: f64, error_tolerance: f64, error_type: MassErrorType) -> Option<&T> {
        return match self.search(query, error_tolerance, error_type) {
            Some(j) => Some(self.get_item(j)),
            None => None,
        };
    }

    fn between(
        &self,
        low: f64,
        high: f64,
        error_tolerance: f64,
        error_type: MassErrorType,
    ) -> &[T] {
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
        subset
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
                break;
            } else {
                upper_index = i;
            }
        }

        return &self.get_slice(lower_index..upper_index + 1);
    }
}

/// Represent a sorted list of processed mass spectral peaks. It is a
/// concrete implementation of [`PeakCollection`] based on a [`Vec`].
#[derive(Default, Clone, Debug)]
pub struct PeakSetVec<P: IndexedCoordinate<C>, C> {
    pub peaks: Vec<P>,
    phantom: marker::PhantomData<C>,
}

impl<P: IndexedCoordinate<C>, C> PeakSetVec<P, C> {
    pub fn new(mut peaks: Vec<P>) -> Self {
        Self::_sort(&mut peaks);
        Self {
            peaks,
            phantom: marker::PhantomData,
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            peaks: Vec::with_capacity(capacity),
            phantom: marker::PhantomData,
        }
    }

    pub fn empty() -> Self {
        Self::with_capacity(0)
    }

    pub fn from<V: Iterator>(peaks: V, sort: bool) -> Self
    where
        V: Iterator<Item = P>,
    {
        let peaks: Vec<P> = peaks.collect();
        if sort {
            Self::new(peaks)
        } else {
            Self::wrap(peaks)
        }
    }

    pub fn wrap(peaks: Vec<P>) -> Self {
        Self {
            peaks,
            phantom: marker::PhantomData,
        }
    }

    fn _sort(peaks: &mut Vec<P>) {
        peaks.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (i, p) in peaks.iter_mut().enumerate() {
            p.set_index(i as u32);
        }
    }

    pub fn iter(&self) -> PeakSetIter<P, C> {
        PeakSetIter::new(self)
    }

    pub fn iter_mut(&mut self) -> PeakSetIterMut<P, C> {
        PeakSetIterMut::new(self)
    }
}

impl<P: IndexedCoordinate<C>, C> PeakCollection<P, C> for PeakSetVec<P, C> {
    fn sort(&mut self) {
        Self::_sort(&mut self.peaks);
    }

    fn push(&mut self, peak: P) {
        let n = self.len();
        match self.peaks.last() {
            Some(p) => {
                if p <= &peak {
                    self.peaks.push(peak);
                    let fin = &mut self.peaks[n];
                    fin.set_index(n as u32);
                } else {
                    self.peaks.push(peak);
                    self.sort();
                }
            }
            None => {
                self.peaks.push(peak);
                let fin = &mut self.peaks[n];
                fin.set_index(n as u32);
            }
        }
    }

    fn len(&self) -> usize {
        self.peaks.len()
    }

    fn get_item(&self, i: usize) -> &P {
        &self[i]
    }

    fn get_slice(&self, i: ops::Range<usize>) -> &[P] {
        &self.peaks[i]
    }

    fn _search_by(&self, query: f64) -> Result<usize, usize> {
        self.peaks
            .binary_search_by(|peak| peak.get_coordinate().partial_cmp(&query).unwrap())
    }
}

impl<P: IndexedCoordinate<C>, C> ops::Index<usize> for PeakSetVec<P, C> {
    type Output = P;

    fn index(&self, i: usize) -> &Self::Output {
        &(self.peaks[i])
    }
}

impl<P: IndexedCoordinate<C>, C> fmt::Display for PeakSetVec<P, C> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "PeakSetVec(<{} Peaks>)", self.len())?;
        Ok(())
    }
}

impl<P: IndexedCoordinate<C>, C> From<Vec<P>> for PeakSetVec<P, C> {
    fn from(v: Vec<P>) -> PeakSetVec<P, C> {
        PeakSetVec::wrap(v)
    }
}

impl<'a, P: IndexedCoordinate<C>, C> IntoIterator for &'a PeakSetVec<P, C> {
    type Item = &'a P;
    type IntoIter = PeakSetIter<'a, P, C>;

    fn into_iter(self) -> Self::IntoIter {
        PeakSetIter::new(self)
    }
}

impl<'a, P: IndexedCoordinate<C>, C> IntoIterator for &'a mut PeakSetVec<P, C> {
    type Item = &'a mut P;
    type IntoIter = PeakSetIterMut<'a, P, C>;

    fn into_iter(self) -> Self::IntoIter {
        PeakSetIterMut::new(self)
    }
}

/// Iterators

// Reference Iterator

pub struct PeakSetIter<'a, P, C> {
    iter: std::slice::Iter<'a, P>,
    phantom: marker::PhantomData<C>,
}

impl<'a, P: IndexedCoordinate<C>, C> PeakSetIter<'a, P, C> {
    fn new(peaks: &'a PeakSetVec<P, C>) -> PeakSetIter<'a, P, C> {
        PeakSetIter {
            iter: peaks.peaks.iter(),
            phantom: marker::PhantomData,
        }
    }
}

impl<'a, P: IndexedCoordinate<C>, C> Iterator for PeakSetIter<'a, P, C> {
    type Item = &'a P;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

// Mutable Reference Iterator

pub struct PeakSetIterMut<'a, P, C> {
    iter: std::slice::IterMut<'a, P>,
    phantom: marker::PhantomData<C>,
}

impl<'a, P: IndexedCoordinate<C>, C> PeakSetIterMut<'a, P, C> {
    fn new(peaks: &'a mut PeakSetVec<P, C>) -> PeakSetIterMut<'a, P, C> {
        PeakSetIterMut {
            iter: peaks.peaks.iter_mut(),
            phantom: marker::PhantomData,
        }
    }
}

impl<'a, P, C> Iterator for PeakSetIterMut<'a, P, C> {
    type Item = &'a mut P;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

pub type PeakSet = PeakSetVec<CentroidPeak, MZ>;
