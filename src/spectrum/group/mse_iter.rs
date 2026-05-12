use std::marker::PhantomData;

use mzpeaks::{IonMobility, Mass, MZ};

use crate::{
    params::{Value, ValueRef},
    prelude::*,
    spectrum::{
        IonMobilityFrameDescription, IonMobilityFrameGroup, MultiLayerIonMobilityFrame,
        MultiLayerSpectrum, SpectrumDescription, SpectrumGroup,
    },
    RawSpectrum,
};

pub trait HasScanConfiguration {
    fn get_configuration(&'_ self) -> Option<ValueRef<'_>>;
}

impl HasScanConfiguration for SpectrumDescription {
    fn get_configuration(&'_ self) -> Option<ValueRef<'_>> {
        self.acquisition
            .first_scan()
            .and_then(|s| s.scan_configuration())
    }
}

impl HasScanConfiguration for IonMobilityFrameDescription {
    fn get_configuration(&'_ self) -> Option<ValueRef<'_>> {
        self.acquisition
            .first_scan()
            .and_then(|s| s.scan_configuration())
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> HasScanConfiguration
    for MultiLayerSpectrum<C, D>
{
    fn get_configuration(&self) -> Option<ValueRef<'_>> {
        self.description.get_configuration()
    }
}

impl HasScanConfiguration for RawSpectrum {
    fn get_configuration(&'_ self) -> Option<ValueRef<'_>> {
        self.description.get_configuration()
    }
}

impl<C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge>
    HasScanConfiguration for MultiLayerIonMobilityFrame<C, D>
{
    fn get_configuration(&'_ self) -> Option<ValueRef<'_>> {
        self.description.get_configuration()
    }
}

/// A low-level implementation detail of the MS^e iterator strategy that is independent of the
/// actual *thing* being looped over.
pub struct MSEIterator<T: HasScanConfiguration, I: Iterator<Item = T>> {
    iterator: I,
    current_low: Option<T>,
    current_high: Vec<T>,
    /// The low energy configuration to treat as the precursor frame
    pub low_configuration: Value,
    /// All configurations to skip and exclude from any grouping, discarding them
    pub ignored_configurations: Vec<Value>,
}

impl<T: HasScanConfiguration, I: Iterator<Item = T>> MSEIterator<T, I> {
    pub fn new(iterator: I, low_configuration: Value, ignored_configurations: Vec<Value>) -> Self {
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

    pub fn ignored_configurations(&self) -> &[Value] {
        &self.ignored_configurations
    }
}

macro_rules! impl_mse_iter {
    () => {
        pub fn new(
            iterator: I,
            low_configuration: Value,
            ignored_configurations: Vec<Value>,
        ) -> Self {
            Self {
                inner: MSEIterator::new(iterator, low_configuration, ignored_configurations),
                _t: PhantomData,
            }
        }

        /// The low energy configuration to treat as the precursor frame
        pub fn low_configuration(&self) -> &Value {
            self.inner.low_configuration()
        }

        /// All configurations to skip and exclude from any grouping, discarding them
        pub fn ignored_configurations(&self) -> &[Value] {
            self.inner.ignored_configurations()
        }
    };
}

macro_rules! impl_mse_next {
    () => {
        fn next(&mut self) -> Option<Self::Item> {
            let parts = self.inner.next_group()?;
            let mut g = G::default();
            if let Some(prec) = parts.0 {
                g.set_precursor(prec);
            }
            *g.products_mut() = parts.1;
            Some(g)
        }
    };
}

/// An MS^e iterator strategy that mirrors the pattern that is used for Waters instruments, applicable to spectra.
pub struct SpectrumMSEIterator<
    C: CentroidLike,
    D: DeconvolutedCentroidLike,
    I: Iterator<Item = MultiLayerSpectrum<C, D>>,
    G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>> = SpectrumGroup<
        C,
        D,
        MultiLayerSpectrum<C, D>,
    >,
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

    impl_mse_next!();
}

#[allow(unused)]
impl<
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
        I: Iterator<Item = MultiLayerSpectrum<C, D>>,
        G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
    > SpectrumMSEIterator<C, D, I, G>
{
    impl_mse_iter!();
}

pub trait SpectrumMSEIteratorExt<
    C: CentroidLike,
    D: DeconvolutedCentroidLike,
    I: Iterator<Item = MultiLayerSpectrum<C, D>>,
    G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>> = SpectrumGroup<
        C,
        D,
        MultiLayerSpectrum<C, D>,
    >,
>
{
    fn into_mse_iterator(
        self,
        low_configuration: Value,
        ignored_configurations: Vec<Value>,
    ) -> SpectrumMSEIterator<C, D, I, G>;
}

impl<
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
        I: Iterator<Item = MultiLayerSpectrum<C, D>>,
    > SpectrumMSEIteratorExt<C, D, I, SpectrumGroup<C, D, MultiLayerSpectrum<C, D>>> for I
{
    fn into_mse_iterator(
        self,
        low_configuration: Value,
        ignored_configurations: Vec<Value>,
    ) -> SpectrumMSEIterator<C, D, I, SpectrumGroup<C, D, MultiLayerSpectrum<C, D>>> {
        SpectrumMSEIterator::new(self, low_configuration, ignored_configurations)
    }
}

pub struct IonMobilityFrameMSEIterator<
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
    I: Iterator<Item = MultiLayerIonMobilityFrame<C, D>>,
    G: IonMobilityFrameGrouping<C, D, MultiLayerIonMobilityFrame<C, D>> = IonMobilityFrameGroup<
        C,
        D,
        MultiLayerIonMobilityFrame<C, D>,
    >,
> {
    inner: MSEIterator<MultiLayerIonMobilityFrame<C, D>, I>,
    _t: PhantomData<G>,
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        G: IonMobilityFrameGrouping<C, D, MultiLayerIonMobilityFrame<C, D>>,
        I: Iterator<Item = MultiLayerIonMobilityFrame<C, D>>,
    > Iterator for IonMobilityFrameMSEIterator<C, D, I, G>
{
    type Item = G;

    impl_mse_next!();
}

#[allow(unused)]
impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        G: IonMobilityFrameGrouping<C, D, MultiLayerIonMobilityFrame<C, D>>,
        I: Iterator<Item = MultiLayerIonMobilityFrame<C, D>>,
    > IonMobilityFrameMSEIterator<C, D, I, G>
{
    impl_mse_iter!();
}

/// An MS^e iterator strategy that mirrors the pattern that is used for Waters instruments, applicable to [`MultiLayerIonMobilityFrame`].
pub trait IonMobilityFrameMSEIteratorExt<
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
    I: Iterator<Item = MultiLayerIonMobilityFrame<C, D>>,
    G: IonMobilityFrameGrouping<C, D, MultiLayerIonMobilityFrame<C, D>>,
>
{
    fn into_mse_iterator(
        self,
        low_configuration: Value,
        ignored_configurations: Vec<Value>,
    ) -> IonMobilityFrameMSEIterator<C, D, I, G>;
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        I: Iterator<Item = MultiLayerIonMobilityFrame<C, D>>,
    >
    IonMobilityFrameMSEIteratorExt<
        C,
        D,
        I,
        IonMobilityFrameGroup<C, D, MultiLayerIonMobilityFrame<C, D>>,
    > for I
{
    fn into_mse_iterator(
        self,
        low_configuration: Value,
        ignored_configurations: Vec<Value>,
    ) -> IonMobilityFrameMSEIterator<
        C,
        D,
        I,
        IonMobilityFrameGroup<C, D, MultiLayerIonMobilityFrame<C, D>>,
    > {
        IonMobilityFrameMSEIterator::new(self, low_configuration, ignored_configurations)
    }
}
