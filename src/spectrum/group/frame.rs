use std::{marker::PhantomData, mem};

use mzpeaks::{
    coordinate::{IonMobility, Mass, MZ},
    feature::FeatureLike, KnownCharge,
};

use crate::{io::IonMobilityFrameGrouping, spectrum::IonMobilityFrameLike};

use super::util::GroupIterState;

/**
A pairing of an optional MS1 ion mobility frame with all its associated MSn ion mobility frames.
*/
#[derive(Debug, Clone)]
pub struct IonMobilityFrameGroup<C, D, S>
where
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
    S: IonMobilityFrameLike<C, D>,
{
    /// The MS1 ion mobility frame of a group. This may be absent when the source does not contain any MS1 ion mobility frames
    pub precursor: Option<S>,
    /// The collection of related MSn ion mobility frames. If MSn for n > 2 is used, all levels are present in this
    /// collection, though there is no ordering guarantee.
    pub products: Vec<S>,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
}

impl<C, D, S> Default for IonMobilityFrameGroup<C, D, S>
where
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
    S: IonMobilityFrameLike<C, D>,
{
    fn default() -> Self {
        Self { precursor: None, products: Vec::new(), centroid_type: Default::default(), deconvoluted_type: Default::default() }
    }
}



impl<C, D, S> IonMobilityFrameGrouping<C, D, S> for IonMobilityFrameGroup<C, D, S>
where
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
    S: IonMobilityFrameLike<C, D>,
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


pub struct IonMobilityFrameGroupIntoIter<
    C: FeatureLike<MZ, IonMobility> + Default,
    D: FeatureLike<Mass, IonMobility> + KnownCharge + Default,
    S: IonMobilityFrameLike<C, D> + Default,
    G: IonMobilityFrameGrouping<C, D, S>,
> {
    group: G,
    state: GroupIterState,
    _c: PhantomData<C>,
    _d: PhantomData<D>,
    _s: PhantomData<S>,
}

impl<
        C: FeatureLike<MZ, IonMobility> + Default,
        D: FeatureLike<Mass, IonMobility> + KnownCharge + Default,
        S: IonMobilityFrameLike<C, D> + Default,
        G: IonMobilityFrameGrouping<C, D, S>,
    > Iterator for IonMobilityFrameGroupIntoIter<C, D, S, G>
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
        C: FeatureLike<MZ, IonMobility> + Default,
        D: FeatureLike<Mass, IonMobility> + KnownCharge + Default,
        S: IonMobilityFrameLike<C, D> + Default,
        G: IonMobilityFrameGrouping<C, D, S>,
    > IonMobilityFrameGroupIntoIter<C, D, S, G>
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

/// Iterate over the spectra in [`IonMobilityFrameGroup`]
pub struct IonMobilityFrameGroupIter<
    'a,
    C: FeatureLike<MZ, IonMobility> + Default,
    D: FeatureLike<Mass, IonMobility> + KnownCharge + Default,
    S: IonMobilityFrameLike<C, D> + Default + 'a,
    G: IonMobilityFrameGrouping<C, D, S>,
> {
    group: &'a G,
    state: GroupIterState,
    _c: PhantomData<C>,
    _d: PhantomData<D>,
    _s: PhantomData<S>,
}

impl<
        'a,
        C: FeatureLike<MZ, IonMobility> + Default,
        D: FeatureLike<Mass, IonMobility> + KnownCharge + Default,
        S: IonMobilityFrameLike<C, D> + Default + 'a,
        G: IonMobilityFrameGrouping<C, D, S>,
    > Iterator for IonMobilityFrameGroupIter<'a, C, D, S, G>
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
        C: FeatureLike<MZ, IonMobility> + Default,
        D: FeatureLike<Mass, IonMobility> + KnownCharge + Default,
        S: IonMobilityFrameLike<C, D> + Default,
        G: IonMobilityFrameGrouping<C, D, S>,
    > IonMobilityFrameGroupIter<'a, C, D, S, G>
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

impl<C, D, S> IntoIterator for IonMobilityFrameGroup<C, D, S>
where
    C: FeatureLike<MZ, IonMobility> + Default,
    D: FeatureLike<Mass, IonMobility> + KnownCharge + Default,
    S: IonMobilityFrameLike<C, D> + Default,
{
    type Item = S;

    type IntoIter = IonMobilityFrameGroupIntoIter<C, D, S, Self>;

    fn into_iter(self) -> Self::IntoIter {
        IonMobilityFrameGroupIntoIter::new(self)
    }
}

impl<'a, C, D, S> IonMobilityFrameGroup<C, D, S>
where
    C: FeatureLike<MZ, IonMobility> + Default,
    D: FeatureLike<Mass, IonMobility> + KnownCharge + Default,
    S: IonMobilityFrameLike<C, D> + Default,
{
    pub fn iter(&'a self) -> IonMobilityFrameGroupIter<'a, C, D, S, Self> {
        IonMobilityFrameGroupIter::new(self)
    }
}
