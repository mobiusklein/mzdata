#![allow(unused)]
use std::{collections::HashSet, marker::PhantomData};

use mzpeaks::{IonMobility, Mass, MZ};

use crate::{
    params::{Value, ValueRef},
    prelude::*,
    spectrum::{
        IonMobilityFrameDescription, MultiLayerIonMobilityFrame, MultiLayerSpectrum,
        SpectrumDescription,
    },
    RawSpectrum,
};

pub(crate) trait HasScanConfiguration {
    fn get_configuration(&self) -> Option<ValueRef>;
}

impl HasScanConfiguration for SpectrumDescription {
    fn get_configuration(&self) -> Option<ValueRef> {
        self.acquisition
            .first_scan()
            .and_then(|s| s.scan_configuration())
    }
}

impl HasScanConfiguration for IonMobilityFrameDescription {
    fn get_configuration(&self) -> Option<ValueRef> {
        self.acquisition
            .first_scan()
            .and_then(|s| s.scan_configuration())
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> HasScanConfiguration
    for MultiLayerSpectrum<C, D>
{
    fn get_configuration(&self) -> Option<ValueRef> {
        self.description.get_configuration()
    }
}

impl HasScanConfiguration for RawSpectrum {
    fn get_configuration(&self) -> Option<ValueRef> {
        self.description.get_configuration()
    }
}

impl<C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge>
    HasScanConfiguration for MultiLayerIonMobilityFrame<C, D>
{
    fn get_configuration(&self) -> Option<ValueRef> {
        self.description.get_configuration()
    }
}

pub struct MSEIterator<T: HasScanConfiguration, I: Iterator<Item = T>> {
    iterator: I,
    current_low: Option<T>,
    current_high: Vec<T>,
    pub low_configuration: Value,
    pub ignored_configurations: HashSet<Value>,
}

impl<T: HasScanConfiguration, I: Iterator<Item = T>> MSEIterator<T, I> {
    pub fn new(
        iterator: I,
        low_configuration: Value,
        ignored_configurations: HashSet<Value>,
    ) -> Self {
        Self {
            iterator,
            current_low: None,
            current_high: Vec::new(),
            low_configuration,
            ignored_configurations,
        }
    }

    fn make_group(&mut self) -> (Option<T>, Vec<T>) {
        let low = self.current_low.take();
        let high = core::mem::take(&mut self.current_high);
        (low, high)
    }

    fn feed_item(&mut self, item: T) -> Option<(Option<T>, Vec<T>)> {
        if let Some(val) = item.get_configuration() {
            if val == self.low_configuration {
                let low = self.current_low.take();
                let high = core::mem::take(&mut self.current_high);
                self.current_low = Some(item);
                Some((low, high))
            } else {
                if !self.ignored_configurations.contains(&(val.into())) {
                    self.current_high.push(item)
                }
                None
            }
        } else {
            None
        }
    }

    pub fn next_group(&mut self) -> Option<(Option<T>, Vec<T>)> {
        loop {
            if let Some(item) = self.iterator.next() {
                if let Some((low, highs)) = self.feed_item(item) {
                    if low.is_some() || !highs.is_empty() {
                        return Some((low, highs));
                    }
                }
            } else {
                let (low, highs) = self.make_group();
                if low.is_some() || !highs.is_empty() {
                    return Some((low, highs));
                } else {
                    return None;
                }
            }
        }
    }

    pub fn low_configuration(&self) -> &Value {
        &self.low_configuration
    }

    pub fn ignored_configurations(&self) -> &HashSet<Value> {
        &self.ignored_configurations
    }
}

pub struct SpectrumMSEIterator<
    C: CentroidLike,
    D: DeconvolutedCentroidLike,
    I: Iterator<Item = MultiLayerSpectrum<C, D>>,
    G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
> {
    inner: MSEIterator<MultiLayerSpectrum<C, D>, I>,
    _t: PhantomData<G>,
}

impl<
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
        I: Iterator<Item = MultiLayerSpectrum<C, D>>,
        G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
    > Iterator for SpectrumMSEIterator<C, D, I, G>
{
    type Item = G;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(parts) = self.inner.next_group() {
            let mut g = G::default();
            if let Some(prec) = parts.0 {
                g.set_precursor(prec);
            }
            *g.products_mut() = parts.1;
            Some(g)
        } else {
            None
        }
    }
}

impl<
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
        I: Iterator<Item = MultiLayerSpectrum<C, D>>,
        G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
    > SpectrumMSEIterator<C, D, I, G>
{
    pub fn new(
        iterator: I,
        low_configuration: Value,
        ignored_configurations: HashSet<Value>,
    ) -> Self {
        Self {
            inner: MSEIterator::new(iterator, low_configuration, ignored_configurations),
            _t: PhantomData,
        }
    }

    pub fn low_configuration(&self) -> &Value {
        self.inner.low_configuration()
    }

    pub fn ignored_configurations(&self) -> &HashSet<Value> {
        self.inner.ignored_configurations()
    }
}
