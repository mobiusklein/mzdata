use std::{
    collections::{HashMap, VecDeque},
    marker::PhantomData,
};

use mzpeaks::{
    feature::FeatureLike, CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak,
    IonMobility, KnownCharge, Mass, MZ,
};

use crate::{
    io::{
        traits::{RandomAccessSpectrumGroupingIterator, SpectrumGrouping},
        IonMobilityFrameGrouping, RandomAccessSpectrumIterator, SpectrumAccessError,
    },
    prelude::{MSDataFileMetadata, SpectrumLike},
};

use super::{IonMobilityFrameLike, MultiLayerIonMobilityFrame, MultiLayerSpectrum};

mod frame;
mod spectrum;
mod util;

pub use frame::{IonMobilityFrameGroup, IonMobilityFrameGroupIntoIter, IonMobilityFrameGroupIter};
pub use spectrum::{SpectrumGroup, SpectrumGroupIntoIter, SpectrumGroupIter};
pub(crate) use util::GenerationTracker;

const MAX_GROUP_DEPTH: u32 = 512u32;

/**
A wrapper for [`Iterator`]-implementors that will batch together
all MSn spectra with their associated MS1 spectrum, producing [`SpectrumGroup`]
instances.

This type emulates the same interface that [`Iterator`] exposes, save that instead
of yield individual [`Spectrum`](crate::spectrum::Spectrum), it yields [`SpectrumGroup`] instead.
*/
#[derive(Debug)]
pub struct SpectrumGroupingIterator<
    R: Iterator<Item = S>,
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumLike<C, D> = MultiLayerSpectrum<C, D>,
    G: SpectrumGrouping<C, D, S> = SpectrumGroup<C, D, S>,
> {
    pub source: R,
    pub queue: VecDeque<S>,
    pub product_scan_mapping: HashMap<String, Vec<S>>,
    generation_tracker: GenerationTracker,
    buffering: usize,
    highest_ms_level: u8,
    generation: usize,
    depth: u32,
    passed_first_ms1: bool,
    phantom: PhantomData<S>,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
    grouping_type: PhantomData<G>,
}

const MISSING_SCAN_ID: &str = "___MISSING_PRECURSOR_ID___";

impl<
        R: Iterator<Item = S>,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        G: SpectrumGrouping<C, D, S>,
    > SpectrumGroupingIterator<R, C, D, S, G>
{
    /// Construct a new [`SpectrumGroupingIterator`] around a [`Iterator`] with a default
    /// buffering level of 3.
    pub fn new(source: R) -> SpectrumGroupingIterator<R, C, D, S, G> {
        SpectrumGroupingIterator::<R, C, D, S, G> {
            source,
            generation_tracker: GenerationTracker::default(),
            phantom: PhantomData,
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
            grouping_type: PhantomData,
            buffering: 3,
            product_scan_mapping: HashMap::new(),
            queue: VecDeque::new(),
            highest_ms_level: 0,
            generation: 0,
            depth: 0,
            passed_first_ms1: false,
        }
    }

    fn add_product(&mut self, scan: S) {
        self.depth += 1;
        if let Some(precursor) = scan.precursor() {
            match precursor.precursor_id.as_ref() {
                Some(prec_id) => {
                    let ent = self.product_scan_mapping.entry(prec_id.clone());
                    self.generation_tracker
                        .add(prec_id.clone(), self.generation);
                    ent.or_default().push(scan);
                }
                None => {
                    let buffer = self
                        .product_scan_mapping
                        .entry(MISSING_SCAN_ID.to_owned())
                        .or_default();
                    buffer.push(scan);
                    if buffer.len() % 1000 == 0 && !buffer.is_empty() {
                        log::warn!("Unassociated MSn scan buffer size is {}", buffer.len());
                    }
                }
            }
        } else if !self.queue.is_empty() {
            // Consider replacing with normal get_mut to avoid re-copying
            let last = self.queue.back().unwrap();
            self.generation_tracker
                .add(last.id().to_owned(), self.generation);
            self.product_scan_mapping
                .entry(last.id().to_owned())
                .or_default()
                .push(scan);
        } else {
            let ent = self.product_scan_mapping.entry(MISSING_SCAN_ID.to_owned());
            let buffer = ent.or_default();
            buffer.push(scan);
            if buffer.len() % 1000 == 0 && !buffer.is_empty() {
                log::warn!("Unassociated MSn scan buffer size is {}", buffer.len());
            }
        }
    }

    fn add_precursor(&mut self, scan: S) -> bool {
        self.queue.push_back(scan);
        self.queue.len() >= self.buffering
    }

    fn pop_precursor(&mut self, precursor_id: &str) -> Vec<S> {
        match self.product_scan_mapping.remove(precursor_id) {
            Some(v) => v,
            None => Vec::new(),
        }
    }

    fn flush_first_ms1(&mut self, group: &mut G) {
        let current_ms1_time = group.precursor().unwrap().start_time();
        let mut ids_to_remove = Vec::new();
        for (prec_id, prods) in self.product_scan_mapping.iter_mut() {
            let mut hold = vec![];
            for prod in prods.drain(..) {
                if prod.start_time() <= current_ms1_time {
                    group.products_mut().push(prod);
                } else {
                    hold.push(prod);
                }
            }
            if hold.is_empty() {
                ids_to_remove.push(prec_id.clone());
            } else {
                prods.extend(hold);
            }
        }
        for prec_id in ids_to_remove {
            self.product_scan_mapping.remove(&prec_id);
        }
    }

    fn flush_all_products(&mut self, group: &mut G) {
        for (_prec_id, prods) in self.product_scan_mapping.drain() {
            group.products_mut().extend(prods);
        }
    }

    fn include_higher_msn(&mut self, group: &mut G) {
        let mut blocks = Vec::new();
        let mut new_block = Vec::new();
        for msn_spec in group.products() {
            new_block.extend(self.pop_precursor(msn_spec.id()));
        }
        blocks.push(new_block);
        for _ in 0..self.highest_ms_level - 2 {
            if let Some(last_block) = blocks.last() {
                let mut new_block = Vec::new();
                for msn_spec in last_block {
                    new_block.extend(self.pop_precursor(msn_spec.id()));
                }
                blocks.push(new_block);
            }
        }
        if !blocks.is_empty() {
            for block in blocks {
                group.products_mut().extend(block);
            }
        }
    }

    fn flush_generations(&mut self, group: &mut G) {
        if self.buffering > self.generation {
            return;
        }
        for prec_id in self
            .generation_tracker
            .older_than(self.generation - self.buffering)
        {
            if let Some(prods) = self.product_scan_mapping.remove(&prec_id) {
                group.products_mut().extend(prods);
            }
        }
    }

    fn deque_group_without_precursor(&mut self) -> G {
        let mut group = G::default();
        self.flush_all_products(&mut group);
        self.generation += 1;
        self.depth = 0;
        group
    }

    fn deque_group(&mut self, flush_all: bool) -> Option<G> {
        let mut group = G::default();
        if let Some(precursor) = self.queue.pop_front() {
            group.set_precursor(precursor);
            let mut products = self.pop_precursor(group.precursor().unwrap().id());
            if self.product_scan_mapping.contains_key(MISSING_SCAN_ID) {
                products.extend(self.pop_precursor(MISSING_SCAN_ID));
            }
            group.products_mut().extend(products);

            self.flush_generations(&mut group);

            // Handle interleaving MS3 and up
            if self.highest_ms_level > 2 {
                self.include_higher_msn(&mut group);
            }

            if !self.passed_first_ms1 && !self.queue.is_empty() {
                self.flush_first_ms1(&mut group);
            }
            self.generation += 1;
            if flush_all {
                self.flush_all_products(&mut group);
            }
            self.depth = 0;
            Some(group)
        } else {
            None
        }
    }

    pub fn clear(&mut self) {
        self.product_scan_mapping.clear();
        self.queue.clear();
        self.generation_tracker.clear();
        self.generation = 0;
        self.depth = 0;
        self.passed_first_ms1 = false;
    }

    /**
    Retrieve the next group of spectra from the iterator, buffering all intermediate and
    interleaved spectra until the next complete group is available or the MS1 buffer is
    full.
    */
    pub fn next_group(&mut self) -> Option<G> {
        loop {
            if let Some(spectrum) = self.source.next() {
                let level = spectrum.ms_level();
                if level > self.highest_ms_level {
                    self.highest_ms_level = level;
                }
                if level > 1 {
                    self.add_product(spectrum);
                    if self.depth > MAX_GROUP_DEPTH {
                        return Some(self.deque_group_without_precursor());
                    }
                } else if self.add_precursor(spectrum) {
                    return self.deque_group(false);
                }
            } else {
                return match self.queue.len() {
                    d if d > 1 => self.deque_group(false),
                    1 => self.deque_group(true),
                    _ => {
                        if !self.product_scan_mapping.is_empty() {
                            Some(self.deque_group_without_precursor())
                        } else {
                            None
                        }
                    }
                };
            }
        }
    }
}

impl<
        R: Iterator<Item = S>,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        G: SpectrumGrouping<C, D, S> + Default,
    > Iterator for SpectrumGroupingIterator<R, C, D, S, G>
{
    type Item = G;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_group()
    }
}

impl<
        R: RandomAccessSpectrumIterator<C, D, S>,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        G: SpectrumGrouping<C, D, S> + Default,
    > SpectrumGroupingIterator<R, C, D, S, G>
{
    pub fn start_from_id(&mut self, id: &str) -> Result<&Self, SpectrumAccessError> {
        match self.source.start_from_id(id) {
            Ok(_) => {
                self.clear();
                Ok(self)
            }
            Err(err) => Err(err),
        }
    }

    pub fn start_from_index(&mut self, index: usize) -> Result<&Self, SpectrumAccessError> {
        match self.source.start_from_index(index) {
            Ok(_) => {
                self.clear();
                Ok(self)
            }
            Err(err) => Err(err),
        }
    }

    pub fn start_from_time(&mut self, time: f64) -> Result<&Self, SpectrumAccessError> {
        match self.source.start_from_time(time) {
            Ok(_) => {
                self.clear();
                Ok(self)
            }
            Err(err) => Err(err),
        }
    }
}

impl<
        R: RandomAccessSpectrumIterator<C, D, S>,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        G: SpectrumGrouping<C, D, S> + Default,
    > RandomAccessSpectrumGroupingIterator<C, D, S, G> for SpectrumGroupingIterator<R, C, D, S, G>
{
    fn start_from_id(&mut self, id: &str) -> Result<&Self, SpectrumAccessError> {
        self.start_from_id(id)
    }

    fn start_from_index(&mut self, index: usize) -> Result<&Self, SpectrumAccessError> {
        self.start_from_index(index)
    }

    fn start_from_time(&mut self, time: f64) -> Result<&Self, SpectrumAccessError> {
        self.start_from_time(time)
    }

    fn reset_state(&mut self) {
        self.clear()
    }
}

impl<
        R: Iterator<Item = S>,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumLike<C, D>,
        G: SpectrumGrouping<C, D, S> + Default,
    > MSDataFileMetadata for SpectrumGroupingIterator<R, C, D, S, G>
where
    R: MSDataFileMetadata,
{
    crate::delegate_impl_metadata_trait!(source);
}

#[derive(Debug)]
pub struct IonMobilityFrameGroupingIterator<
    R: Iterator<Item = S>,
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
    S: IonMobilityFrameLike<C, D> = MultiLayerIonMobilityFrame<C, D>,
    G: IonMobilityFrameGrouping<C, D, S> = IonMobilityFrameGroup<C, D, S>,
> {
    pub source: R,
    pub queue: VecDeque<S>,
    pub product_frame_mapping: HashMap<String, Vec<S>>,
    generation_tracker: GenerationTracker,
    buffering: usize,
    highest_ms_level: u8,
    generation: usize,
    depth: u32,
    passed_first_ms1: bool,
    phantom: PhantomData<S>,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
    grouping_type: PhantomData<G>,
}

impl<
        R: Iterator<Item = S>,
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        S: IonMobilityFrameLike<C, D>,
        G: IonMobilityFrameGrouping<C, D, S>,
    > Iterator for IonMobilityFrameGroupingIterator<R, C, D, S, G>
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
        G: IonMobilityFrameGrouping<C, D, S>,
    > IonMobilityFrameGroupingIterator<R, C, D, S, G>
{
    /// Construct a new [`IonMobilityFrameGroupingIterator`] around a [`Iterator`] with a default
    /// buffering level of 3.
    pub fn new(source: R) -> Self {
        Self {
            source,
            generation_tracker: GenerationTracker::default(),
            phantom: PhantomData,
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
            grouping_type: PhantomData,
            buffering: 3,
            product_frame_mapping: HashMap::new(),
            queue: VecDeque::new(),
            highest_ms_level: 0,
            generation: 0,
            depth: 0,
            passed_first_ms1: false,
        }
    }

    fn add_product(&mut self, scan: S) {
        self.depth += 1;
        if let Some(precursor) = scan.precursor() {
            match precursor.precursor_id.as_ref() {
                Some(prec_id) => {
                    let ent = self.product_frame_mapping.entry(prec_id.clone());
                    self.generation_tracker
                        .add(prec_id.clone(), self.generation);
                    ent.or_default().push(scan);
                }
                None => {
                    let buffer = self
                        .product_frame_mapping
                        .entry(MISSING_SCAN_ID.to_owned())
                        .or_default();
                    buffer.push(scan);
                    if buffer.len() % 1000 == 0 && !buffer.is_empty() {
                        log::warn!("Unassociated MSn frame buffer size is {}", buffer.len());
                    }
                }
            }
        } else if !self.queue.is_empty() {
            // Consider replacing with normal get_mut to avoid re-copying
            let last = self.queue.back().unwrap();
            self.generation_tracker
                .add(last.id().to_owned(), self.generation);
            self.product_frame_mapping
                .entry(last.id().to_owned())
                .or_default()
                .push(scan);
        } else {
            let ent = self.product_frame_mapping.entry(MISSING_SCAN_ID.to_owned());
            let buffer = ent.or_default();
            buffer.push(scan);
            if buffer.len() % 1000 == 0 && !buffer.is_empty() {
                log::warn!("Unassociated MSn frame buffer size is {}", buffer.len());
            }
        }
    }

    fn add_precursor(&mut self, scan: S) -> bool {
        self.queue.push_back(scan);
        self.queue.len() >= self.buffering
    }

    fn pop_precursor(&mut self, precursor_id: &str) -> Vec<S> {
        match self.product_frame_mapping.remove(precursor_id) {
            Some(v) => v,
            None => Vec::new(),
        }
    }

    fn flush_first_ms1(&mut self, group: &mut G) {
        let current_ms1_time = group.precursor().unwrap().start_time();
        let mut ids_to_remove = Vec::new();
        for (prec_id, prods) in self.product_frame_mapping.iter_mut() {
            let mut hold = vec![];
            for prod in prods.drain(..) {
                if prod.start_time() <= current_ms1_time {
                    group.products_mut().push(prod);
                } else {
                    hold.push(prod);
                }
            }
            if hold.is_empty() {
                ids_to_remove.push(prec_id.clone());
            } else {
                prods.extend(hold);
            }
        }
        for prec_id in ids_to_remove {
            self.product_frame_mapping.remove(&prec_id);
        }
    }

    fn flush_all_products(&mut self, group: &mut G) {
        for (_prec_id, prods) in self.product_frame_mapping.drain() {
            group.products_mut().extend(prods);
        }
    }

    fn include_higher_msn(&mut self, group: &mut G) {
        let mut blocks = Vec::new();
        let mut new_block = Vec::new();
        for msn_spec in group.products() {
            new_block.extend(self.pop_precursor(msn_spec.id()));
        }
        blocks.push(new_block);
        for _ in 0..self.highest_ms_level - 2 {
            if let Some(last_block) = blocks.last() {
                let mut new_block = Vec::new();
                for msn_spec in last_block {
                    new_block.extend(self.pop_precursor(msn_spec.id()));
                }
                blocks.push(new_block);
            }
        }
        if !blocks.is_empty() {
            for block in blocks {
                group.products_mut().extend(block);
            }
        }
    }

    fn flush_generations(&mut self, group: &mut G) {
        if self.buffering > self.generation {
            return;
        }
        for prec_id in self
            .generation_tracker
            .older_than(self.generation - self.buffering)
        {
            if let Some(prods) = self.product_frame_mapping.remove(&prec_id) {
                group.products_mut().extend(prods);
            }
        }
    }

    fn deque_group_without_precursor(&mut self) -> G {
        let mut group = G::default();
        self.flush_all_products(&mut group);
        self.generation += 1;
        self.depth = 0;
        group
    }

    fn deque_group(&mut self, flush_all: bool) -> Option<G> {
        let mut group = G::default();
        if let Some(precursor) = self.queue.pop_front() {
            group.set_precursor(precursor);
            let mut products = self.pop_precursor(group.precursor().unwrap().id());
            if self.product_frame_mapping.contains_key(MISSING_SCAN_ID) {
                products.extend(self.pop_precursor(MISSING_SCAN_ID));
            }
            group.products_mut().extend(products);

            self.flush_generations(&mut group);

            // Handle interleaving MS3 and up
            if self.highest_ms_level > 2 {
                self.include_higher_msn(&mut group);
            }

            if !self.passed_first_ms1 && !self.queue.is_empty() {
                self.flush_first_ms1(&mut group);
            }
            self.generation += 1;
            if flush_all {
                self.flush_all_products(&mut group);
            }
            self.depth = 0;
            Some(group)
        } else {
            None
        }
    }

    pub fn clear(&mut self) {
        self.product_frame_mapping.clear();
        self.queue.clear();
        self.generation_tracker.clear();
        self.generation = 0;
        self.depth = 0;
        self.passed_first_ms1 = false;
    }

    /**
    Retrieve the next group of spectra from the iterator, buffering all intermediate and
    interleaved spectra until the next complete group is available or the MS1 buffer is
    full.
    */
    pub fn next_group(&mut self) -> Option<G> {
        loop {
            if let Some(spectrum) = self.source.next() {
                let level = spectrum.ms_level();
                if level > self.highest_ms_level {
                    self.highest_ms_level = level;
                }
                if level > 1 {
                    self.add_product(spectrum);
                    if self.depth > MAX_GROUP_DEPTH {
                        return Some(self.deque_group_without_precursor());
                    }
                } else if self.add_precursor(spectrum) {
                    return self.deque_group(false);
                }
            } else {
                return match self.queue.len() {
                    d if d > 1 => self.deque_group(false),
                    1 => self.deque_group(true),
                    _ => {
                        if !self.product_frame_mapping.is_empty() {
                            Some(self.deque_group_without_precursor())
                        } else {
                            None
                        }
                    }
                };
            }
        }
    }
}

#[cfg(feature = "mzsignal")]
mod mzsignal_impl {
    use std::borrow::Cow;
    use std::sync::Arc;

    use crate::spectrum::bindata::{to_bytes, BuildArrayMapFrom, BuildFromArrayMap};
    use crate::spectrum::{
        ArrayType, BinaryArrayMap, BinaryDataArrayType, DataArray, RefPeakDataLevel,
        SignalContinuity,
    };

    use super::*;

    use log::warn;
    use mzpeaks::{MZLocated, PeakCollection};
    use mzsignal::average::{average_signal, SignalAverager};
    use mzsignal::reprofile::{reprofile, PeakSetReprofiler, PeakShape, PeakShapeModel};
    use mzsignal::{ArrayPair, FittedPeak};

    impl From<ArrayPair<'_>> for BinaryArrayMap {
        fn from(value: ArrayPair<'_>) -> Self {
            let mz_array = DataArray::wrap(
                &ArrayType::MZArray,
                BinaryDataArrayType::Float64,
                to_bytes(&value.mz_array),
            );

            let intensity_array = DataArray::wrap(
                &ArrayType::IntensityArray,
                BinaryDataArrayType::Float32,
                to_bytes(&value.intensity_array),
            );

            let mut array_map = BinaryArrayMap::new();
            array_map.add(mz_array);
            array_map.add(intensity_array);
            array_map
        }
    }

    /// Average a series of [`SpectrumLike`] together. The supplied dx will be used to [`reprofile`]
    /// centroid spectra. The resulting [`ArrayPair`] can be made into a [`BinaryArrayMap`] using
    /// `.into()`. Which in turn can be used to create a new [`MultiLayerSpectrum`] or
    /// [`RawSpectrum`](crate::spectrum::RawSpectrum). Internally it uses [`SpectrumLike::peaks`]
    /// to retrieve the spectrum data. Any deconvoluted spectra will be skipped.
    ///
    /// Note: only available with feature `mzsignal`.
    pub fn average_spectra<
        'lifetime,
        S: SpectrumLike<C, D> + 'lifetime,
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
    >(
        spectra: impl IntoIterator<Item = &'lifetime S>,
        dx: f64,
    ) -> ArrayPair<'static> {
        average_signal(
            &spectra
                .into_iter()
                .flat_map(|spectrum| match spectrum.peaks() {
                    RefPeakDataLevel::Missing | RefPeakDataLevel::Deconvoluted(_) => None,
                    RefPeakDataLevel::RawData(array_map) => {
                        let mzs = array_map.mzs().unwrap().to_vec();
                        let intensities = array_map.intensities().unwrap().to_vec();
                        Some(ArrayPair::from((mzs, intensities)))
                    }
                    RefPeakDataLevel::Centroid(peaks) => {
                        let fitted_peaks: Vec<_> = peaks
                            .iter()
                            .map(|p| FittedPeak::from(p.as_centroid()))
                            .collect();
                        Some(reprofile(fitted_peaks.iter(), dx))
                    }
                })
                .collect::<Vec<_>>(),
            dx,
        )
    }

    #[derive(Debug, Clone)]
    pub struct ArcArrays {
        pub mz_array: Arc<Vec<f64>>,
        pub intensity_array: Arc<Vec<f32>>,
        pub is_profile: bool,
    }

    impl ArcArrays {
        pub fn new(
            mz_array: Arc<Vec<f64>>,
            intensity_array: Arc<Vec<f32>>,
            is_profile: bool,
        ) -> Self {
            Self {
                mz_array,
                intensity_array,
                is_profile,
            }
        }

        pub fn reprofile_on(
            &self,
            reprofiler: &PeakSetReprofiler,
            fwhm: f32,
        ) -> ArrayPair<'static> {
            let models: Vec<_> = self
                .mz_array
                .iter()
                .zip(self.intensity_array.iter())
                .map(|(mz, inten)| {
                    PeakShapeModel::from_centroid(*mz, *inten, fwhm, PeakShape::Gaussian)
                })
                .collect();
            let profiles = reprofiler.reprofile_from_models(&models).to_owned();
            profiles
        }

        pub fn reprofile_with(mut self, reprofiler: &PeakSetReprofiler, fwhm: f32) -> Self {
            if self.is_profile {
                self
            } else {
                let profiles = self.reprofile_on(reprofiler, fwhm);
                self.mz_array = Arc::new(profiles.mz_array.into());
                self.intensity_array = Arc::new(profiles.intensity_array.into());
                self.is_profile = true;
                self
            }
        }
    }

    #[derive(Debug, Default)]
    pub struct SpectrumAveragingContext<
        C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
        D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
        G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
    > {
        pub group: G,
        ms1_context: Vec<ArcArrays>,
        _c: PhantomData<C>,
        _d: PhantomData<D>,
    }

    impl<
            C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
        > SpectrumAveragingContext<C, D, G>
    {
        pub fn new(group: G, ms1_context: Vec<ArcArrays>) -> Self {
            Self {
                group,
                ms1_context,
                _c: PhantomData,
                _d: PhantomData,
            }
        }

        pub fn iter(
            &self,
        ) -> SpectrumGroupIter<'_, C, D, MultiLayerSpectrum<C, D>, SpectrumAveragingContext<C, D, G>>
        {
            SpectrumGroupIter::new(self)
        }

        pub fn reprofile_with(mut self, reprofiler: &PeakSetReprofiler, fwhm: f32) -> Self {
            if self.ms1_context.iter().all(|ctx| ctx.is_profile) {
                return self;
            }
            let ms1_context: Vec<_> = self
                .ms1_context
                .into_iter()
                .map(|ctx| ctx.reprofile_with(reprofiler, fwhm))
                .collect();
            self.ms1_context = ms1_context;
            self
        }

        pub fn reprofile_with_average_with(
            self,
            averager: &mut SignalAverager,
            reprofiler: &PeakSetReprofiler,
        ) -> (G, ArrayPair<'static>) {
            averager.array_pairs.clear();

            for arrays in self.ms1_context.iter() {
                if arrays.is_profile {
                    let mp = arrays.mz_array.as_ptr();
                    let ip = arrays.intensity_array.as_ptr();
                    let unsafe_mzs =
                        unsafe { std::slice::from_raw_parts(mp, arrays.mz_array.len()) };
                    let unsafe_intens =
                        unsafe { std::slice::from_raw_parts(ip, arrays.intensity_array.len()) };
                    let pair = ArrayPair::from((unsafe_mzs, unsafe_intens));
                    averager.push(pair);
                } else {
                    let pair = arrays.reprofile_on(reprofiler, 0.01);
                    averager.push(pair);
                }
            }

            let avg_intens = averager.interpolate();
            let mz_grid = averager.mz_grid.clone();

            let pair_avgd = ArrayPair::from((mz_grid, avg_intens));

            averager.array_pairs.clear();
            (self.group, pair_avgd)
        }

        pub fn average_with(self, averager: &mut SignalAverager) -> (G, ArrayPair<'static>) {
            averager.array_pairs.clear();

            for arrays in self.ms1_context.iter() {
                if arrays.is_profile {
                    let mp = arrays.mz_array.as_ptr();
                    let ip = arrays.intensity_array.as_ptr();
                    let unsafe_mzs =
                        unsafe { std::slice::from_raw_parts(mp, arrays.mz_array.len()) };
                    let unsafe_intens =
                        unsafe { std::slice::from_raw_parts(ip, arrays.intensity_array.len()) };
                    let pair = ArrayPair::from((unsafe_mzs, unsafe_intens));
                    averager.push(pair);
                } else {
                    let peaks: Vec<_> = arrays
                        .mz_array
                        .iter()
                        .zip(arrays.intensity_array.iter())
                        .map(|(mz, intensity)| FittedPeak {
                            mz: *mz,
                            intensity: *intensity,
                            full_width_at_half_max: 0.01,
                            ..Default::default()
                        })
                        .collect();
                    let pair = reprofile(peaks.iter(), averager.dx).to_owned();
                    averager.push(pair);
                }
            }

            let avg_intens = averager.interpolate();
            let mz_grid = averager.mz_grid.clone();

            let pair_avgd = ArrayPair::from((mz_grid, avg_intens));

            averager.array_pairs.clear();
            (self.group, pair_avgd)
        }
    }

    impl<
            C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
        > SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>> for SpectrumAveragingContext<C, D, G>
    {
        fn precursor(&self) -> Option<&MultiLayerSpectrum<C, D>> {
            self.group.precursor()
        }

        fn precursor_mut(&mut self) -> Option<&mut MultiLayerSpectrum<C, D>> {
            self.group.precursor_mut()
        }

        fn set_precursor(&mut self, prec: MultiLayerSpectrum<C, D>) {
            self.group.set_precursor(prec);
        }

        fn products(&self) -> &[MultiLayerSpectrum<C, D>] {
            self.group.products()
        }

        fn products_mut(&mut self) -> &mut Vec<MultiLayerSpectrum<C, D>> {
            self.group.products_mut()
        }

        fn into_parts(
            self,
        ) -> (
            Option<MultiLayerSpectrum<C, D>>,
            Vec<MultiLayerSpectrum<C, D>>,
        ) {
            self.group.into_parts()
        }
    }

    pub struct DeferredSpectrumAveragingIterator<
        C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
        D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
        R: Iterator<Item = G>,
        G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
    > {
        source: R,
        averaging_width_index: usize,
        buffer: VecDeque<G>,
        context_buffer: VecDeque<ArcArrays>,
        output_buffer: VecDeque<SpectrumAveragingContext<C, D, G>>,
    }

    /// If the underlying iterator implements [`MSDataFileMetadata`] then [`DeferredSpectrumAveragingIterator`] will
    /// forward that implementation, assuming it is available.
    impl<
            C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            R: Iterator<Item = G>,
            G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
        > MSDataFileMetadata for DeferredSpectrumAveragingIterator<C, D, R, G>
    where
        R: MSDataFileMetadata,
    {
        crate::delegate_impl_metadata_trait!(source);
    }

    impl<
            C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            R: Iterator<Item = G>,
            G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
        > Iterator for DeferredSpectrumAveragingIterator<C, D, R, G>
    {
        type Item = SpectrumAveragingContext<C, D, G>;

        fn next(&mut self) -> Option<Self::Item> {
            self.next_group()
        }
    }

    impl<
            C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            R: Iterator<Item = G>,
            G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
        > DeferredSpectrumAveragingIterator<C, D, R, G>
    {
        pub fn new(source: R, averaging_width_index: usize) -> Self {
            let mut inst = Self {
                source,
                averaging_width_index,
                buffer: VecDeque::with_capacity(averaging_width_index * 2 + 1),
                context_buffer: VecDeque::with_capacity(averaging_width_index * 2 + 1),
                output_buffer: VecDeque::with_capacity(averaging_width_index * 2 + 1),
            };
            inst.initial_feed();
            inst
        }

        fn reset_buffers(&mut self) {
            self.buffer.clear();
            self.context_buffer.clear();
            self.output_buffer.clear();
            self.initial_feed();
        }

        pub fn into_inner(self) -> R {
            self.source
        }

        pub fn get_ref(&self) -> &R {
            &self.source
        }

        const fn half_capacity(&self) -> usize {
            self.averaging_width_index + 1
        }

        const fn capacity(&self) -> usize {
            self.averaging_width_index * 2 + 1
        }

        fn copy_out_arrays(&self, group: &G) -> Option<ArcArrays> {
            group.precursor().and_then(|scan| {
                if scan.signal_continuity() == SignalContinuity::Profile {
                    if let Some(array_map) = scan.raw_arrays() {
                        let mz = array_map.mzs().unwrap().to_vec();
                        let inten = array_map.intensities().unwrap().to_vec();
                        let mz_a = Arc::new(mz);
                        let inten_a = Arc::new(inten);
                        Some(ArcArrays::new(mz_a, inten_a, true))
                    } else {
                        warn!("{} did not have raw data arrays", scan.id());
                        None
                    }
                } else {
                    if let Some(peaks) = &scan.peaks {
                        let mut mz = Vec::with_capacity(peaks.len());
                        let mut inten = Vec::with_capacity(peaks.len());
                        for p in peaks {
                            mz.push(p.mz());
                            inten.push(p.intensity());
                        }
                        let mz_a = Arc::new(mz);
                        let inten_a = Arc::new(inten);
                        Some(ArcArrays::new(mz_a, inten_a, false))
                    } else if let Some(array_map) = scan.raw_arrays() {
                        let mz = array_map.mzs().unwrap().to_vec();
                        let inten = array_map.intensities().unwrap().to_vec();
                        let mz_a = Arc::new(mz);
                        let inten_a = Arc::new(inten);
                        Some(ArcArrays::new(mz_a, inten_a, false))
                    } else {
                        warn!("{} did not have raw data arrays", scan.id());
                        None
                    }
                }
            })
        }

        /// Populate the context tracking buffers and pre-fill the output queue
        /// with the first productions
        fn initial_feed(&mut self) {
            for _ in 0..self.capacity() {
                let group_opt = self.source.next();
                self.buffer.extend(group_opt.into_iter());
            }
            let mut arrays = VecDeque::with_capacity(self.capacity());
            self.buffer.iter().for_each(|group| {
                if let Some(pair) = self.copy_out_arrays(group) {
                    arrays.push_back(pair)
                }
            });

            for _ in 0..self.half_capacity() {
                if let Some(pair) = arrays.pop_front() {
                    self.context_buffer.push_back(pair);
                }
            }
            let next_group = self._next_group_no_pop();
            self.output_buffer.extend(next_group.into_iter());
            for _ in 0..self.averaging_width_index {
                if let Some(pair) = arrays.pop_front() {
                    self.context_buffer.push_back(pair);
                    let next_group = self._next_group_no_pop();
                    self.output_buffer.extend(next_group.into_iter());
                }
            }
        }

        fn _next_group_no_pop(&mut self) -> Option<SpectrumAveragingContext<C, D, G>> {
            let ms1_context = self.process_block();
            if let Some(group) = self.buffer.pop_front() {
                Some(SpectrumAveragingContext::new(group, ms1_context))
            } else {
                None
            }
        }

        fn update_group(&mut self) -> Option<SpectrumAveragingContext<C, D, G>> {
            let result = self._next_group_no_pop();
            if let Some(group) = self.source.next() {
                self.context_buffer.pop_front();
                if let Some(pair) = self.copy_out_arrays(&group) {
                    self.context_buffer.push_back(pair);
                }
                self.buffer.push_back(group);
            }
            result
        }

        fn process_block(&self) -> Vec<ArcArrays> {
            self.context_buffer.iter().cloned().collect()
        }

        fn next_group(&mut self) -> Option<SpectrumAveragingContext<C, D, G>> {
            if let Some(group) = self.output_buffer.pop_front() {
                Some(group)
            } else {
                self.update_group()
            }
        }
    }

    /// Wrap a `[SpectrumGroupingIterator]` to average MS1 spectra over
    pub struct SpectrumAveragingIterator<
        'lifespan,
        C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
        D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
        G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
        R: Iterator<Item = G>,
    > {
        source: R,
        averaging_width_index: usize,
        buffer: VecDeque<G>,
        output_buffer: VecDeque<G>,
        averager: SignalAverager<'lifespan>,
        reprofiler: PeakSetReprofiler,
        _c: PhantomData<C>,
        _d: PhantomData<D>,
    }

    /// If the underlying iterator implements [`MSDataFileMetadata`] then [`SpectrumAveragingIterator`] will
    /// forward that implementation, assuming it is available.
    impl<
            'lifespan,
            C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
            R: Iterator<Item = G>,
        > MSDataFileMetadata for SpectrumAveragingIterator<'lifespan, C, D, G, R>
    where
        R: MSDataFileMetadata,
    {
        crate::delegate_impl_metadata_trait!(source);
    }

    impl<
            'lifespan,
            C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
            R: Iterator<Item = G>,
        > SpectrumAveragingIterator<'lifespan, C, D, G, R>
    {
        pub fn new(
            source: R,
            averaging_width_index: usize,
            mz_start: f64,
            mz_end: f64,
            dx: f64,
        ) -> Self {
            let mut inst = Self {
                source,
                averaging_width_index,
                buffer: VecDeque::with_capacity(averaging_width_index * 2 + 1),
                output_buffer: VecDeque::with_capacity(averaging_width_index * 2 + 1),
                averager: SignalAverager::new(mz_start, mz_end, dx),
                reprofiler: PeakSetReprofiler::new(mz_start, mz_end, dx),
                _c: PhantomData,
                _d: PhantomData,
            };
            inst.initial_feed();
            inst
        }

        fn reset_buffers(&mut self) {
            self.buffer.clear();
            self.output_buffer.clear();
            self.initial_feed();
        }

        pub fn into_inner(self) -> R {
            self.source
        }

        pub const fn get_ref(&self) -> &R {
            &self.source
        }

        const fn half_capacity(&self) -> usize {
            self.averaging_width_index + 1
        }

        const fn capacity(&self) -> usize {
            self.averaging_width_index * 2 + 1
        }

        fn copy_out_arrays(&self, group: &G) -> Option<ArrayPair<'static>> {
            group.precursor().and_then(|scan| {
                if scan.signal_continuity() == SignalContinuity::Profile {
                    if let Some(array_map) = scan.raw_arrays() {
                        let mz = array_map.mzs().unwrap().to_vec();
                        let inten = array_map.intensities().unwrap().to_vec();
                        Some(ArrayPair::from((mz, inten)))
                    } else {
                        warn!("{} did not have raw data arrays", scan.id());
                        None
                    }
                } else {
                    if let Some(peaks) = scan.peaks.as_ref() {
                        let fpeaks: Vec<_> = peaks
                            .iter()
                            .map(|p| {
                                let fp = FittedPeak::from(p.as_centroid());
                                PeakShapeModel {
                                    peak: Cow::Owned(fp),
                                    shape: PeakShape::Gaussian,
                                }
                            })
                            .collect();
                        let signal = self.reprofiler.reprofile(&fpeaks);
                        Some(ArrayPair::from((
                            signal.mz_array.to_vec(),
                            signal.intensity_array.to_vec(),
                        )))
                    } else {
                        warn!(
                            "{} was not in profile mode but no centroids found",
                            scan.id()
                        );
                        None
                    }
                }
            })
        }

        /// Populate the context tracking buffers and pre-fill the output queue
        /// with the first productions
        fn initial_feed(&mut self) {
            for _ in 0..self.capacity() {
                let group_opt = self.source.next();
                self.buffer.extend(group_opt.into_iter());
            }
            let mut arrays = VecDeque::with_capacity(self.capacity());
            self.buffer.iter().for_each(|group| {
                if let Some(pair) = self.copy_out_arrays(group) {
                    arrays.push_back(pair)
                }
            });

            for _ in 0..self.half_capacity() {
                if let Some(pair) = arrays.pop_front() {
                    self.averager.push(pair);
                }
            }
            let next_group = self._next_group_no_pop();
            self.output_buffer.extend(next_group.into_iter());
            for _ in 0..self.averaging_width_index {
                if let Some(pair) = arrays.pop_front() {
                    self.averager.push(pair);
                    let next_group = self._next_group_no_pop();
                    self.output_buffer.extend(next_group.into_iter());
                }
            }
        }

        fn _next_group_no_pop(&mut self) -> Option<G> {
            let array_map = self.average_spectra();
            if let Some(mut group) = self.buffer.pop_front() {
                group.precursor_mut().and_then(|precursor| {
                    precursor.arrays = Some(array_map);
                    Some(())
                });
                Some(group)
            } else {
                None
            }
        }

        fn update_group(&mut self) -> Option<G> {
            let result = self._next_group_no_pop();
            if let Some(group) = self.source.next() {
                self.averager.pop();
                if let Some(pair) = self.copy_out_arrays(&group) {
                    self.averager.push(pair);
                }
                self.buffer.push_back(group);
            }
            result
        }

        fn average_spectra(&self) -> BinaryArrayMap {
            let first_intensity = self.averager.interpolate();
            let first_mz = self.averager.mz_grid.as_slice();

            let mz_array = DataArray::wrap(
                &ArrayType::MZArray,
                BinaryDataArrayType::Float64,
                to_bytes(first_mz),
            );
            let intensity_array = DataArray::wrap(
                &ArrayType::IntensityArray,
                BinaryDataArrayType::Float32,
                to_bytes(&first_intensity),
            );

            let mut array_map = BinaryArrayMap::new();
            array_map.add(mz_array);
            array_map.add(intensity_array);
            array_map
        }

        fn next_group(&mut self) -> Option<G> {
            if let Some(group) = self.output_buffer.pop_front() {
                Some(group)
            } else {
                self.update_group()
            }
        }
    }

    impl<
            'lifespan,
            C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
            R: Iterator<Item = G>,
        > Iterator for SpectrumAveragingIterator<'lifespan, C, D, G, R>
    {
        type Item = G;

        fn next(&mut self) -> Option<Self::Item> {
            self.next_group()
        }
    }

    impl<
            'lifespan,
            C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
            R: RandomAccessSpectrumGroupingIterator<C, D, MultiLayerSpectrum<C, D>, G>,
        > RandomAccessSpectrumGroupingIterator<C, D, MultiLayerSpectrum<C, D>, G>
        for SpectrumAveragingIterator<'lifespan, C, D, G, R>
    {
        fn reset_state(&mut self) {
            self.source.reset_state();
            self.reset_buffers()
        }

        fn start_from_id(&mut self, id: &str) -> Result<&Self, SpectrumAccessError> {
            match self.source.start_from_id(id) {
                Ok(_) => {
                    self.reset_state();
                    Ok(self)
                }
                Err(err) => Err(err),
            }
        }

        fn start_from_index(&mut self, index: usize) -> Result<&Self, SpectrumAccessError> {
            match self.source.start_from_index(index) {
                Ok(_) => {
                    self.reset_state();
                    Ok(self)
                }
                Err(err) => Err(err),
            }
        }

        fn start_from_time(&mut self, time: f64) -> Result<&Self, SpectrumAccessError> {
            match self.source.start_from_time(time) {
                Ok(_) => {
                    self.reset_state();
                    Ok(self)
                }
                Err(err) => Err(err),
            }
        }
    }

    impl<
            'lifespan,
            C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            R: RandomAccessSpectrumGroupingIterator<C, D, MultiLayerSpectrum<C, D>, G>,
            G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
        >
        RandomAccessSpectrumGroupingIterator<
            C,
            D,
            MultiLayerSpectrum<C, D>,
            SpectrumAveragingContext<C, D, G>,
        > for DeferredSpectrumAveragingIterator<C, D, R, G>
    {
        fn start_from_id(&mut self, id: &str) -> Result<&Self, SpectrumAccessError> {
            match self.source.start_from_id(id) {
                Ok(_) => {
                    self.reset_state();
                    Ok(self)
                }
                Err(err) => Err(err),
            }
        }

        fn start_from_index(&mut self, index: usize) -> Result<&Self, SpectrumAccessError> {
            match self.source.start_from_index(index) {
                Ok(_) => {
                    self.reset_state();
                    Ok(self)
                }
                Err(err) => Err(err),
            }
        }

        fn start_from_time(&mut self, time: f64) -> Result<&Self, SpectrumAccessError> {
            match self.source.start_from_time(time) {
                Ok(_) => {
                    self.reset_state();
                    Ok(self)
                }
                Err(err) => Err(err),
            }
        }

        fn reset_state(&mut self) {
            self.source.reset_state();
            self.reset_buffers();
        }
    }

    /// Adds signal averaging to an [`Iterator`] that produces [`SpectrumGrouping`] implementations
    /// of the appropriate type.
    pub trait SpectrumGroupAveraging<
        'lifespan,
        C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
        D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
        G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
    >: Iterator<Item = G> + Sized
    {
        /// Average MS1 spectra across [`SpectrumGrouping`](crate::io::traits::SpectrumGrouping) from this iterator
        ///
        /// # Arguments
        ///
        /// * `averaging_width_index` - The number of groups before and after the current group to average MS1 scans across
        /// * `mz_start` - The minimum m/z to average from
        /// * `mz_end` - The maximum m/z to average up to
        /// * `dx` - The m/z spacing in the averaged spectra
        ///
        fn averaging(
            self,
            averaging_width_index: usize,
            mz_start: f64,
            mz_end: f64,
            dx: f64,
        ) -> SpectrumAveragingIterator<'lifespan, C, D, G, Self> {
            SpectrumAveragingIterator::new(self, averaging_width_index, mz_start, mz_end, dx)
        }

        /// Create an iterator that defers averaging MS1 spectra across [`SpectrumGrouping`](crate::io::traits::SpectrumGrouping)
        /// to a future point by separately returning the machinery for consuming the prepared data, i.e. on another thread.
        ///
        /// # Arguments
        ///
        /// * `averaging_width_index` - The number of groups before and after the current group to average MS1 scans across
        /// * `mz_start` - The minimum m/z to average from
        /// * `mz_end` - The maximum m/z to average up to
        /// * `dx` - The m/z spacing in the averaged spectra
        ///
        fn averaging_deferred(
            self,
            averaging_width_index: usize,
            mz_start: f64,
            mz_end: f64,
            dx: f64,
        ) -> (
            DeferredSpectrumAveragingIterator<C, D, Self, G>,
            SignalAverager<'lifespan>,
            PeakSetReprofiler,
        ) {
            let iter = DeferredSpectrumAveragingIterator::new(self, averaging_width_index);
            let averager = SignalAverager::new(mz_start, mz_end, dx);
            let reprofiler = PeakSetReprofiler::new(mz_start, mz_end, dx);
            (iter, averager, reprofiler)
        }
    }

    impl<
            'lifespan,
            T,
            C: CentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + BuildFromArrayMap,
            G: SpectrumGrouping<C, D, MultiLayerSpectrum<C, D>>,
        > SpectrumGroupAveraging<'lifespan, C, D, G> for T
    where
        T: Iterator<Item = G>,
    {
    }
}

#[cfg(feature = "mzsignal")]
pub use mzsignal_impl::{
    average_spectra, DeferredSpectrumAveragingIterator, SpectrumAveragingIterator,
    SpectrumGroupAveraging,
};

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_group_iter() {
        let group: SpectrumGroup<
            CentroidPeak,
            DeconvolutedPeak,
            MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak>,
        > = SpectrumGroup::default();
        let entries: Vec<_> = group.iter().collect();
        assert_eq!(entries.len(), 0);
    }
}
