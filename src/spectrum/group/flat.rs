use std::marker::PhantomData;

use mzpeaks::{CentroidLike, DeconvolutedCentroidLike, IonMobility, Mass, MZ};

use crate::{io::IonMobilityFrameAccessError, prelude::*};

pub struct FlatSpectrumGroupingIterator<
    R: Iterator<Item = S>,
    C: CentroidLike,
    D: DeconvolutedCentroidLike,
    S: SpectrumLike<C, D>,
    G: SpectrumGrouping<C, D, S>,
> {
    source: R,
    _t: PhantomData<(C, D, S, G)>,
}

impl<
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
        S: SpectrumLike<C, D>,
        G: SpectrumGrouping<C, D, S>,
        I: Iterator<Item = S>,
    > Iterator for FlatSpectrumGroupingIterator<I, C, D, S, G>
{
    type Item = G;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_group()
    }
}

impl<
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
        S: SpectrumLike<C, D>,
        G: SpectrumGrouping<C, D, S>,
        I: Iterator<Item = S>,
    > FlatSpectrumGroupingIterator<I, C, D, S, G>
{
    pub fn new(source: I) -> Self {
        Self {
            source,
            _t: PhantomData,
        }
    }

    pub fn next_group(&mut self) -> Option<G> {
        self.source.next().map(|s| {
            let mut group = G::default();
            if s.ms_level() == 1 {
                group.set_precursor(s);
            } else {
                group.products_mut().push(s);
            }
            group
        })
    }
}


impl<
        R: RandomAccessSpectrumIterator<C, D, S>,
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
        S: SpectrumLike<C, D>,
        G: SpectrumGrouping<C, D, S>,
    > RandomAccessSpectrumGroupingIterator<C, D, S, G> for FlatSpectrumGroupingIterator<R, C, D, S, G>
{
    fn start_from_id(&mut self, id: &str) -> Result<&Self, SpectrumAccessError> {
        self.source.start_from_id(id)?;
        Ok(self)
    }

    fn start_from_index(&mut self, index: usize) -> Result<&Self, SpectrumAccessError> {
        self.source.start_from_index(index)?;
        Ok(self)
    }

    fn start_from_time(&mut self, time: f64) -> Result<&Self, SpectrumAccessError> {
        self.source.start_from_time(time)?;
        Ok(self)
    }

    fn reset_state(&mut self) {}
}

impl<
        R: Iterator<Item = S>,
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
        S: SpectrumLike<C, D>,
        G: SpectrumGrouping<C, D, S>,
    > MSDataFileMetadata for FlatSpectrumGroupingIterator<R, C, D, S, G>
where
    R: MSDataFileMetadata,
{
    crate::delegate_impl_metadata_trait!(source);
}


pub struct FlatIonMobilityGroupingIterator<
    R: Iterator<Item = S>,
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
    S: IonMobilityFrameLike<C, D>,
    G: IonMobilityFrameGrouping<C, D, S>,
> {
    source: R,
    _t: PhantomData<(C, D, S, G)>,
}

impl<
    R: Iterator<Item = S>,
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
    S: IonMobilityFrameLike<C, D>,
    G: IonMobilityFrameGrouping<C, D, S>
    > Iterator for FlatIonMobilityGroupingIterator<R, C, D, S, G>
{
    type Item = G;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_group()
    }
}

impl<
    R: Iterator<Item = S>,
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
    S: IonMobilityFrameLike<C, D>,
    G: IonMobilityFrameGrouping<C, D, S>
    > FlatIonMobilityGroupingIterator<R, C, D, S, G>
{
    pub fn new(source: R) -> Self {
        Self {
            source,
            _t: PhantomData,
        }
    }

    pub fn next_group(&mut self) -> Option<G> {
        self.source.next().map(|s| {
            let mut group = G::default();
            if s.ms_level() == 1 {
                group.set_precursor(s);
            } else {
                group.products_mut().push(s);
            }
            group
        })
    }
}


impl<
        R: RandomAccessIonMobilityFrameIterator<C, D, S>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        G: IonMobilityFrameGrouping<C, D, S>
    > RandomAccessIonMobilityFrameGroupingIterator<C, D, S, G> for FlatIonMobilityGroupingIterator<R, C, D, S, G>
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, IonMobilityFrameAccessError> {
        self.source.start_from_id(id)?;
        Ok(self)
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, IonMobilityFrameAccessError> {
        self.source.start_from_index(index)?;
        Ok(self)
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, IonMobilityFrameAccessError> {
        self.source.start_from_time(time)?;
        Ok(self)
    }

    fn reset_state(&mut self) {}
}

impl<
        R: Iterator<Item = S>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        G: IonMobilityFrameGrouping<C, D, S>
    > MSDataFileMetadata for FlatIonMobilityGroupingIterator<R, C, D, S, G>
where
    R: MSDataFileMetadata,
{
    crate::delegate_impl_metadata_trait!(source);
}