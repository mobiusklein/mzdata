use std::{marker::PhantomData, mem};

use mzpeaks::{CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak};

use super::super::{MultiLayerSpectrum, SpectrumLike};

use super::util::GroupIterState;


/// An abstraction over [`SpectrumGroup`](crate::spectrum::SpectrumGroup)'s interface.
pub trait SpectrumGrouping<
    C: CentroidLike = CentroidPeak,
    D: DeconvolutedCentroidLike = DeconvolutedPeak,
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
        let product = self.products().iter().max_by(|a, b| {
            a.acquisition()
                .start_time()
                .total_cmp(&b.acquisition().start_time())
        });
        match (self.precursor(), product) {
            (None, None) => None,
            (None, Some(c)) => Some(c),
            (Some(p), None) => Some(p),
            (Some(p), Some(c)) => {
                if p.start_time() > c.start_time() {
                    Some(p)
                } else {
                    Some(c)
                }
            },
        }
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

/**
A pairing of an optional MS1 spectrum with all its associated MSn spectra.
*/
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SpectrumGroup<C = CentroidPeak, D = DeconvolutedPeak, S = MultiLayerSpectrum<C, D>>
where
    C: CentroidLike,
    D: DeconvolutedCentroidLike,
    S: SpectrumLike<C, D>,
{
    /// The MS1 spectrum of a group. This may be absent when the source does not contain any MS1 spectra
    pub precursor: Option<S>,
    /// The collection of related MSn spectra. If MSn for n > 2 is used, all levels are present in this
    /// collection, though there is no ordering guarantee.
    pub products: Vec<S>,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
}

impl<C, D, S> IntoIterator for SpectrumGroup<C, D, S>
where
    C: CentroidLike,
    D: DeconvolutedCentroidLike,
    S: SpectrumLike<C, D> + Default,
{
    type Item = S;

    type IntoIter = SpectrumGroupIntoIter<C, D, S, Self>;

    fn into_iter(self) -> Self::IntoIter {
        SpectrumGroupIntoIter::new(self)
    }
}

impl<'a, C, D, S> SpectrumGroup<C, D, S>
where
    C: CentroidLike,
    D: DeconvolutedCentroidLike,
    S: SpectrumLike<C, D>,
{
    pub fn new(precursor: Option<S>, products: Vec<S>) -> Self {
        Self { precursor, products, centroid_type: PhantomData, deconvoluted_type: PhantomData }
    }

    pub fn iter(&'a self) -> SpectrumGroupIter<'a, C, D, S> {
        SpectrumGroupIter::new(self)
    }
}

pub struct SpectrumGroupIntoIter<
    C: CentroidLike = CentroidPeak,
    D: DeconvolutedCentroidLike = DeconvolutedPeak,
    S: SpectrumLike<C, D> + Default = MultiLayerSpectrum<C, D>,
    G: SpectrumGrouping<C, D, S> = SpectrumGroup<C, D, S>,
> {
    group: G,
    state: GroupIterState,
    _c: PhantomData<C>,
    _d: PhantomData<D>,
    _s: PhantomData<S>,
}

impl<
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
        S: SpectrumLike<C, D> + Default,
        G: SpectrumGrouping<C, D, S>,
    > Iterator for SpectrumGroupIntoIter<C, D, S, G>
{
    type Item = S;

    fn next(&mut self) -> Option<Self::Item> {
        {
            let n = self.n_products();
            let emission = match self.state {
                GroupIterState::Precursor => match self.group.precursor_mut() {
                    Some(prec) => {
                        if n > 0 {
                            self.state = GroupIterState::Product(0);
                        } else {
                            self.state = GroupIterState::Done;
                        }
                        Some(mem::take(prec))
                    }
                    None => {
                        if n > 0 {
                            self.state = if n > 1 {
                                GroupIterState::Product(1)
                            } else {
                                GroupIterState::Done
                            };
                            Some(mem::take(&mut self.group.products_mut()[0]))
                        } else {
                            self.state = GroupIterState::Done;
                            None
                        }
                    }
                },
                GroupIterState::Product(i) => {
                    if i < n.saturating_sub(1) {
                        self.state = GroupIterState::Product(i + 1);
                        Some(mem::take(&mut self.group.products_mut()[i]))
                    } else {
                        self.state = GroupIterState::Done;
                        Some(mem::take(&mut self.group.products_mut()[i]))
                    }
                }
                GroupIterState::Done => None,
            };
            emission
        }
    }
}

impl<
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
        S: SpectrumLike<C, D> + Default,
        G: SpectrumGrouping<C, D, S>,
    > SpectrumGroupIntoIter<C, D, S, G>
{
    pub fn new(group: G) -> Self {
        Self {
            group,
            state: GroupIterState::Precursor,
            _c: PhantomData,
            _d: PhantomData,
            _s: PhantomData,
        }
    }

    fn n_products(&self) -> usize {
        self.group.products().len()
    }
}

/// Iterate over the spectra in [`SpectrumGroup`]
pub struct SpectrumGroupIter<
    'a,
    C: CentroidLike = CentroidPeak,
    D: DeconvolutedCentroidLike = DeconvolutedPeak,
    S: SpectrumLike<C, D> = MultiLayerSpectrum<C, D>,
    G: SpectrumGrouping<C, D, S> = SpectrumGroup<C, D, S>,
> {
    group: &'a G,
    state: GroupIterState,
    _c: PhantomData<C>,
    _d: PhantomData<D>,
    _s: PhantomData<S>,
}

impl<
        'a,
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
        S: SpectrumLike<C, D> + 'a,
        G: SpectrumGrouping<C, D, S>,
    > Iterator for SpectrumGroupIter<'a, C, D, S, G>
{
    type Item = &'a S;

    fn next(&mut self) -> Option<Self::Item> {
        {
            let n = self.n_products();
            let emission = match self.state {
                GroupIterState::Precursor => match self.group.precursor() {
                    Some(prec) => {
                        if n > 0 {
                            self.state = GroupIterState::Product(0);
                        } else {
                            self.state = GroupIterState::Done;
                        }
                        Some(prec)
                    }
                    None => {
                        if n > 0 {
                            self.state = if n > 1 {
                                GroupIterState::Product(1)
                            } else {
                                GroupIterState::Done
                            };
                            Some(&self.group.products()[0])
                        } else {
                            self.state = GroupIterState::Done;
                            None
                        }
                    }
                },
                GroupIterState::Product(i) => {
                    if i < n.saturating_sub(1) {
                        self.state = GroupIterState::Product(i + 1);
                        Some(&self.group.products()[i])
                    } else {
                        self.state = GroupIterState::Done;
                        Some(&self.group.products()[i])
                    }
                }
                GroupIterState::Done => None,
            };
            emission
        }
    }
}

impl<
        'a,
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
        S: SpectrumLike<C, D>,
        G: SpectrumGrouping<C, D, S>,
    > SpectrumGroupIter<'a, C, D, S, G>
{
    pub fn new(group: &'a G) -> Self {
        Self {
            group,
            state: GroupIterState::Precursor,
            _c: PhantomData,
            _d: PhantomData,
            _s: PhantomData,
        }
    }

    fn n_products(&self) -> usize {
        self.group.products().len()
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike, S: SpectrumLike<C, D>>
    Default for SpectrumGroup<C, D, S>
{
    fn default() -> Self {
        Self {
            precursor: None,
            products: Vec::new(),
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
        }
    }
}

impl<C, D, S> SpectrumGrouping<C, D, S> for SpectrumGroup<C, D, S>
where
    C: CentroidLike,
    D: DeconvolutedCentroidLike,
    S: SpectrumLike<C, D>,
{
    fn precursor(&self) -> Option<&S> {
        match &self.precursor {
            Some(prec) => Some(prec),
            None => None,
        }
    }

    fn precursor_mut(&mut self) -> Option<&mut S> {
        match &mut self.precursor {
            Some(prec) => Some(prec),
            None => None,
        }
    }

    fn set_precursor(&mut self, prec: S) {
        self.precursor = Some(prec)
    }

    fn products(&self) -> &[S] {
        &self.products
    }

    fn products_mut(&mut self) -> &mut Vec<S> {
        &mut self.products
    }

    fn into_parts(self) -> (Option<S>, Vec<S>) {
        (self.precursor, self.products)
    }
}


#[cfg(test)]
mod test {
    use crate::Spectrum;

    use super::*;

    fn make_group() -> SpectrumGroup {
        let mut spec1 = Spectrum::default();
        {
            let desc1 = spec1.description_mut();
            desc1.id = "index=0".into();
            desc1.ms_level = 1;
            desc1.index = 0;
            desc1.acquisition.first_scan_mut().unwrap().start_time = 100.0;
        }

        let mut spec2 = Spectrum::default();
        {
            let desc2 = spec2.description_mut();
            desc2.id = "index=1".into();
            desc2.ms_level = 2;
            desc2.index = 1;
            desc2.acquisition.first_scan_mut().unwrap().start_time = 101.0;
        }

        let mut spec3 = Spectrum::default();
        {
            let desc3 = spec3.description_mut();
            desc3.id = "index=2".into();
            desc3.ms_level = 2;
            desc3.index = 2;
            desc3.acquisition.first_scan_mut().unwrap().start_time = 101.5;
        }

        SpectrumGroup::new(Some(spec1), vec![
            spec2,
            spec3
        ])
    }

    #[test]
    fn test_construct() {
        let group = make_group();
        assert_eq!(group.lowest_ms_level().unwrap(), 1);
        assert_eq!(group.highest_ms_level().unwrap(), 2);

        assert_eq!(group.earliest_spectrum().unwrap().ms_level(), 1);
        assert_eq!(group.latest_spectrum().unwrap().ms_level(), 2);

        assert_eq!(group.iter().count(), 3);
        assert_eq!(group.into_iter().count(), 3);
    }
}