use log::warn;
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fs;
use std::io;
use std::marker::PhantomData;
use std::path::{self, PathBuf};

use mzpeaks::{CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak};

use crate::spectrum::spectrum::{MultiLayerSpectrum, SpectrumBehavior};

use super::utils::FileSource;
use super::OffsetIndex;

pub trait SeekRead: io::Read + io::Seek {}
impl<T: io::Read + io::Seek> SeekRead for T {}

/// A base trait defining the behaviors of a source of spectra.
pub trait ScanSource<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumBehavior<C, D> = MultiLayerSpectrum<C, D>,
>: Iterator<Item = S>
{
    fn reset(&mut self);

    /// Retrieve a spectrum by it's native ID
    fn get_spectrum_by_id(&mut self, id: &str) -> Option<S>;

    /// Retrieve a spectrum by it's integer index
    fn get_spectrum_by_index(&mut self, index: usize) -> Option<S>;

    /// Retrieve a spectrum by its scan start time
    /// Considerably more complex than seeking by ID or index.
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

    /// Retrieve the number of spectra in source file
    fn len(&self) -> usize {
        self.get_index().len()
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn get_index(&self) -> &OffsetIndex;

    fn set_index(&mut self, index: OffsetIndex);

    /// Helper method to support seeking to an ID
    fn _offset_of_id(&self, id: &str) -> Option<u64> {
        self.get_index().get(id)
    }

    /// Helper method to support seeking to an index
    fn _offset_of_index(&self, index: usize) -> Option<u64> {
        self.get_index()
            .get_index(index)
            .map(|(_id, offset)| offset)
    }

    /// Helper method to support seeking to a specific time.
    /// Considerably more complex than seeking by ID or index.
    fn _offset_of_time(&mut self, time: f64) -> Option<u64> {
        match self.get_spectrum_by_time(time) {
            Some(scan) => self._offset_of_index(scan.index()),
            None => None,
        }
    }

    /// Open a new iterator over this stream.
    fn iter(&mut self) -> SpectrumIterator<C, D, S, Self>
    where
        Self: Sized,
    {
        SpectrumIterator::new(self)
    }

    fn groups(&mut self) -> SpectrumGroupingIterator<Self, C, D, S>
    where
        Self: Sized,
    {
        SpectrumGroupingIterator::new(self.iter())
    }
}

/// A generic iterator over a [`ScanSource`] implementer that assumes the
/// source has already been indexed. Otherwise, the source's own iterator
/// behavior should be used.
pub struct SpectrumIterator<
    'lifespan,
    C: CentroidLike + Default,
    D: DeconvolutedCentroidLike + Default,
    S: SpectrumBehavior<C, D>,
    R: ScanSource<C, D, S>,
> {
    source: &'lifespan mut R,
    phantom: PhantomData<S>,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
    index: usize,
    back_index: usize,
}

impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        R: ScanSource<C, D, S>,
        S: SpectrumBehavior<C, D>,
    > SpectrumIterator<'lifespan, C, D, S, R>
{
    pub fn new(source: &mut R) -> SpectrumIterator<C, D, S, R> {
        SpectrumIterator::<C, D, S, R> {
            source,
            index: 0,
            back_index: 0,
            phantom: PhantomData,
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
        }
    }
}

impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        R: ScanSource<C, D, S>,
        S: SpectrumBehavior<C, D>,
    > Iterator for SpectrumIterator<'lifespan, C, D, S, R>
{
    type Item = S;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index + self.back_index >= self.len() {
            return None;
        }
        let result = self.source.get_spectrum_by_index(self.index);
        self.index += 1;
        result
    }
}

impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        R: ScanSource<C, D, S>,
        S: SpectrumBehavior<C, D>,
    > ExactSizeIterator for SpectrumIterator<'lifespan, C, D, S, R>
{
    fn len(&self) -> usize {
        self.source.len()
    }
}

impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        R: ScanSource<C, D, S>,
        S: SpectrumBehavior<C, D>,
    > DoubleEndedIterator for SpectrumIterator<'lifespan, C, D, S, R>
{
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.index + self.back_index >= self.len() {
            return None;
        };
        let i = self.len() - (self.back_index + 1);
        let result = self.source.get_spectrum_by_index(i);
        self.back_index += 1;
        result
    }
}

impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumBehavior<C, D>,
        R: ScanSource<C, D, S>,
    > ScanSource<C, D, S> for SpectrumIterator<'lifespan, C, D, S, R>
{
    fn reset(&mut self) {
        self.index = 0;
        self.back_index = 0;
    }

    fn get_spectrum_by_id(&mut self, id: &str) -> Option<S> {
        self.source.get_spectrum_by_id(id)
    }

    fn get_spectrum_by_index(&mut self, index: usize) -> Option<S> {
        self.source.get_spectrum_by_index(index)
    }

    fn get_spectrum_by_time(&mut self, time: f64) -> Option<S> {
        self.source.get_spectrum_by_time(time)
    }

    fn get_index(&self) -> &OffsetIndex {
        self.source.get_index()
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.source.set_index(index);
    }
}

/// A trait defining some helper methods to make efficient use of indices
/// automatic when opening a file from a path-like object.
pub trait MZFileReader<
    C: CentroidLike + Default,
    D: DeconvolutedCentroidLike + Default,
    S: SpectrumBehavior<C, D>,
>: ScanSource<C, D, S> + Sized
{
    /// An on-trait method of constructing an index. Assumed
    /// to be a trivial wrapper.
    fn construct_index_from_stream(&mut self) -> u64;

    /// Re-construct an offset index from this readable object, assuming
    /// it is a JSON stream over the serialized index.
    fn read_index(&mut self, reader: Box<dyn io::Read>) -> Result<&Self, serde_json::Error> {
        match OffsetIndex::from_reader(reader) {
            Ok(index) => {
                self.set_index(index);
                Ok(self)
            }
            Err(err) => Err(err),
        }
    }

    fn write_index(&self, writer: Box<dyn io::Write>) -> Result<&Self, serde_json::Error> {
        match self.get_index().to_writer(writer) {
            Ok(_) => Ok(self),
            Err(err) => Err(err),
        }
    }

    /// The preferred method of opening a file from a path-like object.
    /// This method will open the file at the provided path, test whether
    /// there is an accompanied index file next to it on the file system,
    /// and if not, build one and save it or otherwise read in the index.
    ///
    /// The index building process is usually neglible on "regular" IO file
    /// systems.
    fn open_path<P>(path: P) -> io::Result<Self>
    where
        P: Into<path::PathBuf> + Clone,
    {
        let source: FileSource<fs::File> = FileSource::from(path.clone());
        let index_file_name = source.index_file_name();

        match fs::File::open(path.into()) {
            Ok(file) => {
                let mut reader = Self::open_file(file);
                if let Some(index_path) = &index_file_name {
                    if index_path.exists() {
                        let index_stream = fs::File::open(index_path)?;
                        match reader.read_index(Box::new(io::BufReader::new(index_stream))) {
                            Ok(_) => {}
                            Err(_err) => {
                                reader.construct_index_from_stream();
                                _save_index(&index_path, &reader)?;
                            }
                        }
                    } else {
                        reader.construct_index_from_stream();
                        _save_index(&index_path, &reader)?;
                    }
                }
                Ok(reader)
            }
            Err(err) => Err(err),
        }
    }

    /// Given a regular file, construct a new instance without indexing.
    fn open_file(source: fs::File) -> Self;
}

fn _save_index<
    C: CentroidLike + Default,
    D: DeconvolutedCentroidLike + Default,
    S: SpectrumBehavior<C, D>,
>(index_path: &PathBuf, reader: &impl MZFileReader<C, D, S>) -> io::Result<()> {
    let index_stream = fs::File::create(index_path)?;
    match reader.write_index(Box::new(io::BufWriter::new(index_stream))) {
        Ok(_) => {}
        Err(err) => {
            warn!(
                "Failed to write index to {} because {:?}",
                index_path.display(),
                err
            );
        }
    }
    Ok(())
}

#[derive(Debug)]
pub enum ScanAccessError {
    ScanNotFound,
    IOError(Option<io::Error>),
}

impl std::fmt::Display for ScanAccessError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::error::Error for ScanAccessError {}

pub trait RandomAccessSpectrumIterator<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumBehavior<C, D> = MultiLayerSpectrum<C, D>,
>: ScanSource<C, D, S>
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, ScanAccessError>;
    fn start_from_index(&mut self, id: usize) -> Result<&mut Self, ScanAccessError>;
    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, ScanAccessError>;
}

impl<
        'lifespan,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumBehavior<C, D>,
        R: ScanSource<C, D, S>,
    > RandomAccessSpectrumIterator<C, D, S> for SpectrumIterator<'lifespan, C, D, S, R>
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, ScanAccessError> {
        if let Some(scan) = self.get_spectrum_by_id(id) {
            self.index = scan.index();
            self.back_index = 0;
            Ok(self)
        } else {
            if self.get_index().contains_key(id) {
                Err(ScanAccessError::IOError(None))
            } else {
                Err(ScanAccessError::ScanNotFound)
            }
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, ScanAccessError> {
        if index < self.len() {
            self.index = index;
            self.back_index = 0;
            Ok(self)
        } else {
            Err(ScanAccessError::ScanNotFound)
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, ScanAccessError> {
        if let Some(scan) = self.get_spectrum_by_time(time) {
            self.index = scan.index();
            self.back_index = 0;
            Ok(self)
        } else {
            if self
                .get_spectrum_by_index(self.len() - 1)
                .expect("Failed to fetch spectrum for boundary testing")
                .start_time()
                < time
            {
                Err(ScanAccessError::ScanNotFound)
            } else {
                Err(ScanAccessError::IOError(None))
            }
        }
    }
}

pub trait SpectrumGrouping<
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumBehavior<C, D> = MultiLayerSpectrum<C, D>,
>: Default
{
    /// Get the precursor spectrum, which may be absent
    fn precursor(&self) -> Option<&S>;
    /// Get a mutable reference to the precursor spectrum, which may be absent
    fn precursor_mut(&mut self) -> Option<&mut S>;
    /// Explicitly set the precursor spectrum directly.
    fn set_precursor(&mut self, prec: S);

    /// Get a reference to the collection of product spectra
    fn products(&self) -> &Vec<S>;

    /// Get a mutable reference to the collection of product spectra
    fn products_mut(&mut self) -> &mut Vec<S>;
}

/**
A pairing of an optional MS1 spectrum with all its associated MSn spectra.
*/
#[derive(Debug, Clone)]
pub struct SpectrumGroup<C = CentroidPeak, D = DeconvolutedPeak, S = MultiLayerSpectrum<C, D>>
where
    C: CentroidLike + Default,
    D: DeconvolutedCentroidLike + Default,
    S: SpectrumBehavior<C, D>,
{
    /// The MS1 spectrum of a group. This may be absent when the source does not contain any MS1 spectra
    pub precursor: Option<S>,
    /// The collection of related MSn spectra. If MSn for n > 2 is used, all levels are present in this
    /// collection, though there is no ordering guarantee.
    pub products: Vec<S>,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
}

impl<
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumBehavior<C, D>,
    > Default for SpectrumGroup<C, D, S>
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
    C: CentroidLike + Default,
    D: DeconvolutedCentroidLike + Default,
    S: SpectrumBehavior<C, D>,
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

    fn products(&self) -> &Vec<S> {
        &self.products
    }

    fn products_mut(&mut self) -> &mut Vec<S> {
        &mut self.products
    }
}

#[derive(Default)]
struct GenerationTracker {
    generation_to_id: HashMap<usize, HashSet<String>>,
    id_to_generation: HashMap<String, usize>,
    generations: VecDeque<usize>,
}

impl GenerationTracker {
    fn add_generation(&mut self, generation: usize) {
        match self.generations.binary_search(&generation) {
            Ok(i) => {
                self.generations.insert(i, generation);
            }
            Err(i) => {
                self.generations.insert(i, generation);
            }
        }
    }

    pub fn clear(&mut self) {
        self.generation_to_id.clear();
        self.id_to_generation.clear();
        self.generations.clear();
    }

    pub fn add(&mut self, identifier: String, generation: usize) {
        if !self.generation_to_id.contains_key(&generation) {
            self.add_generation(generation);
        }
        self.generation_to_id
            .entry(generation)
            .or_default()
            .insert(identifier.clone());
        match self.id_to_generation.entry(identifier) {
            Entry::Occupied(mut e) => {
                e.insert(generation);
            }
            Entry::Vacant(e) => {
                e.insert(generation);
            }
        }
    }

    #[allow(unused)]
    pub fn remove(&mut self, identifier: String) -> bool {
        match self.id_to_generation.entry(identifier.clone()) {
            Entry::Occupied(ent) => {
                let generation = ent.get();
                let id_set = self
                    .generation_to_id
                    .get_mut(generation)
                    .unwrap_or_else(|| {
                        panic!("Generation {} did not contain {}", generation, identifier)
                    });
                id_set.remove(&identifier);
                if id_set.is_empty() {
                    self.generations.remove(*generation);
                }
                true
            }
            Entry::Vacant(_ent) => false,
        }
    }

    pub fn older_than(&mut self, generation: usize) -> Vec<String> {
        let mut result = Vec::new();
        for gen in self.generations.iter() {
            if *gen < generation {
                if let Some(members) = self.generation_to_id.remove(gen) {
                    result.extend(members);
                }
            } else {
                break;
            }
        }
        for r in result.iter() {
            self.id_to_generation.remove(r);
        }
        result
    }
}

/**
A wrapper for [`SpectrumIterator`]-implementors that will batch together
all MSn spectra with their associated MS1 spectrum, producing [`SpectrumGroup`]
instances.

This type emulates the same interface that [`SpectrumIterator`] exposes, save that instead
of yield individual [`Spectrum`](crate::spectrum::spectrum::Spectrum), it yields [`SpectrumGroup`] instead. Naturally, it aslo
implements the [`Iterator`] trait.
*/
pub struct SpectrumGroupingIterator<
    'lifespan,
    R: ScanSource<C, D, S>,
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumBehavior<C, D> = MultiLayerSpectrum<C, D>,
    G: SpectrumGrouping<C, D, S> = SpectrumGroup<C, D, S>,
> {
    pub source: SpectrumIterator<'lifespan, C, D, S, R>,
    pub queue: VecDeque<S>,
    pub product_scan_mapping: HashMap<String, Vec<S>>,
    generation_tracker: GenerationTracker,
    buffering: usize,
    highest_ms_level: u8,
    generation: usize,
    passed_first_ms1: bool,
    phantom: PhantomData<S>,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
    grouping_type: PhantomData<G>,
}

const MISSING_SCAN_ID: &str = "___MISSING_PRECURSOR_ID___";

impl<
        'lifespan,
        R: ScanSource<C, D, S>,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumBehavior<C, D>,
        G: SpectrumGrouping<C, D, S>,
    > SpectrumGroupingIterator<'lifespan, R, C, D, S, G>
{
    /// Construct a new [`SpectrumGroupingIterator`] around a [`SpectrumIterator`] with a default
    /// buffering level of 3.
    pub fn new(
        source: SpectrumIterator<'lifespan, C, D, S, R>,
    ) -> SpectrumGroupingIterator<'lifespan, R, C, D, S, G> {
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
            passed_first_ms1: false,
        }
    }

    fn add_product(&mut self, scan: S) {
        if let Some(precursor) = scan.precursor() {
            match precursor.precursor_id.as_ref() {
                Some(prec_id) => {
                    let ent = self
                        .product_scan_mapping
                        .entry(prec_id.clone());
                    self.generation_tracker
                        .add(prec_id.clone(), self.generation);
                    ent.or_default().push(scan);
                },
                None => {
                    self.product_scan_mapping.entry(MISSING_SCAN_ID.to_owned()).or_default().push(scan);
                }
            }
        } else {
            if !self.queue.is_empty() {
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
                ent.or_default().push(scan);
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
            Some(group)
        } else {
            None
        }
    }

    fn clear(&mut self) {
        self.product_scan_mapping.clear();
        self.queue.clear();
        self.generation_tracker.clear();
    }

    /**
    Retrieve the next group of spectra from the iterator, buffering all intermediate and
    interleaved spectra until the next complete group is available or the MS1 buffer is
    full.
    */
    pub fn next_group(&mut self) -> Option<G> {
        if let Some(spectrum) = self.source.next() {
            let level = spectrum.ms_level();
            if level > self.highest_ms_level {
                self.highest_ms_level = level;
            }
            if level > 1 {
                self.add_product(spectrum);
                self.next_group()
            } else {
                if self.add_precursor(spectrum) {
                    self.deque_group(false)
                } else {
                    self.next_group()
                }
            }
        } else {
            match self.queue.len() {
                d if d > 1 => self.deque_group(false),
                1 => self.deque_group(true),
                _ => None,
            }
        }
    }
}

impl<
        'lifespan,
        R: ScanSource<C, D, S>,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumBehavior<C, D>,
        G: SpectrumGrouping<C, D, S> + Default,
    > Iterator for SpectrumGroupingIterator<'lifespan, R, C, D, S, G>
{
    type Item = G;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_group()
    }
}

impl<
        'lifespan,
        R: RandomAccessSpectrumIterator<C, D, S>,
        C: CentroidLike + Default,
        D: DeconvolutedCentroidLike + Default,
        S: SpectrumBehavior<C, D>,
        G: SpectrumGrouping<C, D, S> + Default,
    > SpectrumGroupingIterator<'lifespan, R, C, D, S, G>
{
    pub fn start_from_id(&mut self, id: &str) -> Result<&Self, ScanAccessError> {
        match self.source.start_from_id(id) {
            Ok(_) => {
                self.clear();
                Ok(self)
            }
            Err(err) => Err(err),
        }
    }

    pub fn start_from_index(&mut self, index: usize) -> Result<&Self, ScanAccessError> {
        match self.source.start_from_index(index) {
            Ok(_) => {
                self.clear();
                Ok(self)
            }
            Err(err) => Err(err),
        }
    }

    pub fn start_from_time(&mut self, time: f64) -> Result<&Self, ScanAccessError> {
        match self.source.start_from_time(time) {
            Ok(_) => {
                self.clear();
                Ok(self)
            }
            Err(err) => Err(err),
        }
    }
}

/// Common interface for spectrum writing
pub trait ScanWriter<
    'a,
    W: io::Write,
    C: CentroidLike + Default = CentroidPeak,
    D: DeconvolutedCentroidLike + Default = DeconvolutedPeak,
    S: SpectrumBehavior<C, D> + 'static = MultiLayerSpectrum<C, D>,
>
{
    /// Write out a single spectrum, returning the number of bytes written
    fn write(&mut self, spectrum: &'a S) -> io::Result<usize>;

    /// As [`std::io::Write::flush`]
    fn flush(&mut self) -> io::Result<()>;

    /// Consume an [`Iterator`] over [`Spectrum`](crate::spectrum::MultiLayerSpectrum) references
    fn write_all<T: Iterator<Item = &'a S>>(&mut self, iterator: T) -> io::Result<usize> {
        let mut n = 0;
        for spectrum in iterator {
            n += self.write(spectrum)?;
        }
        Ok(n)
    }

    /// Write a [`SpectrumGroup`] out in order, returning the number of bytes written
    fn write_group<G: SpectrumGrouping<C, D, S> + 'static>(
        &mut self,
        group: &'a G,
    ) -> io::Result<usize> {
        let mut n = 0;
        if let Some(precursor) = group.precursor() {
            n += self.write(precursor)?;
        }
        for product in group.products() {
            n += self.write(product)?;
        }
        Ok(n)
    }

    /// Consume an [`Iterator`] over [`SpectrumGroup`] references
    fn write_all_groups<G: SpectrumGrouping<C, D, S> + 'static, T: Iterator<Item = &'a G>>(
        &mut self,
        iterator: T,
    ) -> io::Result<usize> {
        let mut n = 0;
        for group in iterator {
            n += self.write_group(group)?;
        }
        Ok(n)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_object_safe() {
        // If `ScanSource` were not object safe, this code
        // couldn't compile.
        let _f = |_x: &dyn ScanSource| {};
    }
}
