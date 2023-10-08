#![allow(unused)]
use std::borrow::Cow;
use std::collections::HashMap;
use std::fmt::Display;
use std::io::{self, prelude::*, SeekFrom};
use std::ops::Range;
use std::path::Path;
use std::{fs, mem};

use filename;
use hdf5::{self, filters, Dataset, Selection};
use log::{debug, warn};
use ndarray::Ix1;

use mzpeaks::{CentroidPeak, DeconvolutedPeak};

use crate::io::mzml::{
    CVParamParse, IncrementingIdMap, MzMLParserError, MzMLParserState, MzMLReaderType, MzMLSAX,
    MzMLSpectrumBuilder, ParserResult, SpectrumBuilding,
};
use crate::io::prelude::MSDataFileMetadata;
use crate::io::traits::MZFileReader;
use crate::io::{OffsetIndex, RandomAccessSpectrumIterator, ScanAccessError, ScanSource};
// use crate::io::traits::{ScanSource, RandomAccessSpectrumIterator};

use crate::meta::{DataProcessing, FileDescription, InstrumentConfiguration, Software};
use crate::params::{ControlledVocabulary, Param};
use crate::spectrum::signal::{
    to_bytes, ArrayRetrievalError, ArrayType, BinaryCompressionType, BinaryDataArrayType,
    ByteArrayView, DataArray, as_bytes,
};
use crate::spectrum::spectrum::{
    CentroidPeakAdapting, DeconvolutedPeakAdapting, MultiLayerSpectrum,
};
use crate::spectrum::{IsolationWindow, ScanWindow, SelectedIon};

#[derive(Debug)]
pub enum MzMLbError {
    HDF5Error(hdf5::Error),
    MzMLError(MzMLParserError),
    ArrayRetrievalError(ArrayRetrievalError),
}

impl Display for MzMLbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(format!("{:?}", self).as_str())
    }
}

impl std::error::Error for MzMLbError {}

impl From<MzMLbError> for io::Error {
    fn from(value: MzMLbError) -> Self {
        Self::new(io::ErrorKind::Other, value)
    }
}

impl From<MzMLParserError> for MzMLbError {
    fn from(value: MzMLParserError) -> Self {
        Self::MzMLError(value)
    }
}

impl From<hdf5::Error> for MzMLbError {
    fn from(value: hdf5::Error) -> Self {
        Self::HDF5Error(value)
    }
}

impl From<ArrayRetrievalError> for MzMLbError {
    fn from(value: ArrayRetrievalError) -> Self {
        Self::ArrayRetrievalError(value)
    }
}

#[derive(Debug, Default, Clone)]
struct CacheInterval {
    pub start: usize,
    pub end: usize,
    pub data: DataArray,
}

impl<'transient, 'lifespan: 'transient> ByteArrayView<'transient, 'lifespan> for CacheInterval {
    fn view(&'lifespan self) -> Result<std::borrow::Cow<'lifespan, [u8]>, ArrayRetrievalError> {
        <DataArray as ByteArrayView<'transient, 'lifespan>>::view(&self.data)
    }

    fn dtype(&self) -> BinaryDataArrayType {
        <DataArray as ByteArrayView<'transient, 'lifespan>>::dtype(&self.data)
    }
}

impl PartialEq for CacheInterval {
    fn eq(&self, other: &Self) -> bool {
        self.start == other.start && self.end == other.end
    }
}

impl CacheInterval {
    pub fn new(start: usize, end: usize, data: DataArray) -> Self {
        Self { start, end, data }
    }

    #[inline]
    pub fn contains(&self, start: usize, end: usize) -> bool {
        if self.start <= start {
            if end < self.end {
                return true;
            }
        }
        return false;
    }

    pub fn copy_from_slice(&mut self, src: &[u8]) {
        let size = self.data.data.len();
        let source_size = src.len();
        if size < source_size {
            self.data.data.resize(source_size, 0)
        }
        self.data.data.copy_from_slice(src)
    }

    #[inline]
    pub fn get(&self, start: usize, end: usize) -> Result<Cow<'_, [u8]>, ArrayRetrievalError> {
        let buffer_start = start - self.start;
        let buffer_end = end - self.start;
        let buffer = self.data.slice_buffer(
            buffer_start * self.dtype().size_of(),
            buffer_end * self.dtype().size_of(),
        );
        buffer
    }
}

impl From<hdf5::Datatype> for BinaryDataArrayType {
    fn from(value: hdf5::Datatype) -> Self {
        match value.size() {
            1 => Self::ASCII,
            4 => {
                if value.is::<i32>() {
                    Self::Int32
                } else if value.is::<f32>() {
                    Self::Float32
                } else {
                    Self::Unknown
                }
            }
            8 => {
                if value.is::<i64>() {
                    Self::Int64
                } else if value.is::<f64>() {
                    Self::Float64
                } else {
                    Self::Unknown
                }
            }
            _ => Self::Unknown,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct DataRangeRequest {
    name: String,
    offset: usize,
    length: usize,
}

#[derive(Debug)]
pub struct ExternalDataRegistry {
    chunk_size: usize,
    registry: HashMap<String, hdf5::Dataset>,
    chunk_cache: HashMap<String, CacheInterval>,
}

impl Default for ExternalDataRegistry {
    fn default() -> Self {
        Self {
            chunk_size: 2usize.pow(22),
            registry: Default::default(),
            chunk_cache: Default::default(),
        }
    }
}

impl ExternalDataRegistry {
    pub fn new(chunk_size: usize) -> Self {
        let mut inst = Self::default();
        inst.chunk_size = chunk_size;
        inst
    }

    pub fn from_hdf5(handle: &hdf5::File) -> Result<Self, MzMLbError> {
        let mut storage = Self::default();
        match handle.datasets() {
            Ok(datasets) => {
                datasets.into_iter().for_each(|ds| {
                    let name = ds.name();
                    match name.as_str() {
                        "/mzML"
                        | "/mzML_spectrumIndex"
                        | "/mzML_spectrumIndex_idRef"
                        | "/mzML_chromatogramIndex"
                        | "/mzML_chromatogramIndex_idRef" => {}
                        _ => {
                            storage.registry.insert(name, ds);
                        }
                    };
                });
                Ok(storage)
            }
            Err(err) => return Err(err.into()),
        }
    }

    fn read_slice_into(
        dataset: &Dataset,
        start: usize,
        end: usize,
        buffer: &mut Vec<u8>
    ) -> Result<(), hdf5::Error> {
        let dtype = dataset.dtype()?;
        let sel: Selection = (start..end).into();
        let mztype: BinaryDataArrayType = dtype.into();
        match mztype {
            BinaryDataArrayType::Unknown => todo!(),
            BinaryDataArrayType::Float64 => {
                let block = dataset
                    .read_slice_1d::<f64, _>(sel)
                    .expect("Expected to read block from dataset");
                let view = as_bytes(block.as_slice().unwrap());
                buffer.resize(view.len(), 0);
                buffer.copy_from_slice(view);
            }
            BinaryDataArrayType::Float32 => {
                let block = dataset
                    .read_slice_1d::<f32, _>(sel)
                    .expect("Expected to read block from dataset");
                let view = as_bytes(block.as_slice().unwrap());
                buffer.resize(view.len(), 0);
                buffer.copy_from_slice(view);
            }
            BinaryDataArrayType::Int64 => {
                let block = dataset
                    .read_slice_1d::<i64, _>(sel)
                    .expect("Expected to read block from dataset");
                let view = as_bytes(block.as_slice().unwrap());
                buffer.resize(view.len(), 0);
                buffer.copy_from_slice(view);
            }
            BinaryDataArrayType::Int32 => {
                let block = dataset
                    .read_slice_1d::<i32, _>(sel)
                    .expect("Expected to read block from dataset");
                let view = as_bytes(block.as_slice().unwrap());
                buffer.resize(view.len(), 0);
                buffer.copy_from_slice(view);
            }
            BinaryDataArrayType::ASCII => {
                let block = dataset
                    .read_slice_1d::<u8, _>(sel)
                    .expect("Expected to read block from dataset");
                let view = as_bytes(block.as_slice().unwrap());
                buffer.resize(view.len(), 0);
                buffer.copy_from_slice(view);
            }
        };
        Ok(())
    }

    #[inline]
    fn read_slice_to_bytes(
        dataset: &Dataset,
        start: usize,
        end: usize,
    ) -> Result<Vec<u8>, hdf5::Error> {
        let mut block = Vec::new();
        Self::read_slice_into(dataset, start, end, &mut block)?;
        Ok(block)
    }

    pub fn get(
        &mut self,
        range_request: &DataRangeRequest,
        destination: &mut DataArray,
    ) -> Result<(), MzMLbError> {
        let z = destination.dtype().size_of();
        let start = range_request.offset;
        let end = range_request.offset + (range_request.length);
        destination.compression = BinaryCompressionType::Decoded;
        if let Some(chunk) = self.chunk_cache.get(&range_request.name) {
            if chunk.contains(start, end) {
                let block = chunk.get(start, end)?;
                destination.data.extend_from_slice(&block);
                assert_eq!(destination.data.len(),  range_request.length * z);
                return Ok(());
            }
        }
        if let Some(dset) = self.registry.get(&range_request.name).as_ref() {
            let dtype = dset.dtype()?;
            let block_end = (start + self.chunk_size).min(dset.size());
            let cache_block = Self::read_slice_to_bytes(&dset, start, block_end)?;
            let cache_block = DataArray::wrap(&destination.name, dtype.into(), cache_block);
            let mut cache_block = CacheInterval::new(start, block_end, cache_block);
            let block = if let Some(cache_block) = self
                .chunk_cache
                .get_mut(&range_request.name)
                .and_then(|prev| {
                    mem::swap(prev, &mut cache_block);
                    Some(prev)
                }) {
                debug!(
                    "Updated {} cache block: {:?}",
                    range_request.name, cache_block
                );
                let block = cache_block.get(start, end)?;
                block
            } else {
                self.chunk_cache
                    .insert(range_request.name.clone(), cache_block);
                let cache_block = self.chunk_cache.get_mut(&range_request.name).unwrap();
                let block = cache_block.get(start, end)?;
                block
            };
            destination.data.extend_from_slice(&block);
            assert_eq!(destination.data.len(),  range_request.length * z);
            Ok(())
        } else {
            Err(hdf5::Error::Internal(format!("Group {} not found", range_request.name)).into())
        }
    }
}

struct ByteReader {
    handle: hdf5::Dataset,
    position: usize,
}

impl ByteReader {
    fn new(handle: hdf5::Dataset, position: usize) -> Self {
        Self { handle, position }
    }
}

impl Read for ByteReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.handle.size();
        match self
            .handle
            .read_slice::<u8, _, Ix1>(self.position..(self.position + buf.len()).min(n))
        {
            Ok(slc) => {
                let slc_buf = slc.as_slice().unwrap();
                buf[..slc_buf.len()].copy_from_slice(slc_buf);
                self.position += slc.len();
                Ok(slc.len())
            }
            Err(e) => Err(MzMLbError::HDF5Error(e).into()),
        }
    }
}

impl Seek for ByteReader {
    fn seek(&mut self, pos: io::SeekFrom) -> io::Result<u64> {
        match pos {
            io::SeekFrom::Start(offset) => {
                self.position = (self.position + offset as usize).min(0);
            }
            io::SeekFrom::End(offset) => {
                let mut n = self.handle.size();
                self.position = (n + offset as usize).min(0);
            }
            io::SeekFrom::Current(offset) => {
                self.position = (self.position + offset as usize).max(0);
            }
        }
        Ok(self.position as u64)
    }
}

impl From<hdf5::Dataset> for ByteReader {
    fn from(value: hdf5::Dataset) -> Self {
        Self::new(value, 0)
    }
}

pub trait DataRegistryBorrower<'a> {
    fn borrow_data_registry(self, data_registry: &'a mut ExternalDataRegistry) -> Self;
}

#[derive(Default)]
pub struct MzMLbSpectrumBuilder<'a, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> {
    inner: MzMLSpectrumBuilder<'a, C, D>,
    data_registry: Option<&'a mut ExternalDataRegistry>,
    current_data_range_query: DataRangeRequest,
    // compression: Param,
}

impl<'a, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MzMLbSpectrumBuilder<'a, C, D> {
    pub fn new() -> Self {
        Self::default()
    }
}

impl<'a, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MzMLSAX
    for MzMLbSpectrumBuilder<'a, C, D>
{
    fn start_element(
        &mut self,
        event: &quick_xml::events::BytesStart,
        state: MzMLParserState,
    ) -> crate::io::mzml::ParserResult {
        self.inner.start_element(event, state)
    }

    fn empty_element(
        &mut self,
        event: &quick_xml::events::BytesStart,
        state: MzMLParserState,
        reader_position: usize,
    ) -> ParserResult {
        let elt_name = event.name();
        match elt_name.as_ref() {
            b"cvParam" | b"userParam" => {
                match self.inner.handle_param(event, reader_position, state) {
                    Ok(param) => {
                        match &state {
                            MzMLParserState::BinaryDataArray => {
                                if !param.is_controlled()
                                    || param.controlled_vocabulary.unwrap()
                                        != ControlledVocabulary::MS
                                {
                                    self.inner.fill_param_into(param, state)
                                } else {
                                    match param.accession.unwrap() {
                                        // external HDF5 dataset
                                        1002841 => {
                                            if self.current_data_range_query.name.len() == 0
                                                && !param.value.starts_with("/")
                                            {
                                                self.current_data_range_query
                                                    .name
                                                    .extend("/".chars());
                                            }
                                            self.current_data_range_query
                                                .name
                                                .extend(param.value.chars());
                                        }
                                        // external offset
                                        1002842 => {
                                            self.current_data_range_query.offset = param
                                                .value
                                                .parse()
                                                .expect("Failed to extract external offset")
                                        }
                                        // external array length
                                        1002843 => {
                                            self.current_data_range_query.length = param
                                                .value
                                                .parse()
                                                .expect("Failed to extract external array length")
                                        }
                                        _ => self.inner.fill_param_into(param, state),
                                    }
                                }
                            }
                            _ => self.inner.fill_param_into(param, state),
                        }
                        return Ok(state);
                    }
                    Err(err) => return Err(err),
                }
            }
            &_ => {}
        }
        Ok(state)
    }

    fn end_element(
        &mut self,
        event: &quick_xml::events::BytesEnd,
        state: MzMLParserState,
    ) -> ParserResult {
        let elt_name = event.name();
        match elt_name.as_ref() {
            b"binaryDataArray" => {
                if self.current_data_range_query.name.len() == 0 {
                    return Err(MzMLParserError::IncompleteElementError(
                        "The external data array name was missing or empty".to_owned(),
                        MzMLParserState::BinaryDataArray,
                    ));
                }
                let data_request = mem::take(&mut self.current_data_range_query);
                let array = self.inner.current_array_mut();
                self.data_registry
                    .as_mut()
                    .expect("Did not provide data registry")
                    .get(&data_request, array);
            }
            _ => {}
        };
        self.inner.end_element(event, state)
    }

    fn text(
        &mut self,
        event: &quick_xml::events::BytesText,
        state: MzMLParserState,
    ) -> ParserResult {
        self.inner.text(event, state)
    }
}

impl<'a, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> DataRegistryBorrower<'a>
    for MzMLbSpectrumBuilder<'a, C, D>
{
    fn borrow_data_registry(mut self, data_registry: &'a mut ExternalDataRegistry) -> Self {
        self.data_registry = Some(data_registry);
        self
    }
}

impl<'a, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    SpectrumBuilding<'a, C, D, MultiLayerSpectrum<C, D>> for MzMLbSpectrumBuilder<'a, C, D>
{
    fn isolation_window_mut(&mut self) -> &mut IsolationWindow {
        self.inner.isolation_window_mut()
    }

    fn scan_window_mut(&mut self) -> &mut ScanWindow {
        self.inner.scan_window_mut()
    }

    fn selected_ion_mut(&mut self) -> &mut SelectedIon {
        self.inner.selected_ion_mut()
    }

    fn current_array_mut(&mut self) -> &mut DataArray {
        self.inner.current_array_mut()
    }

    fn into_spectrum(self, spectrum: &mut MultiLayerSpectrum<C, D>) {
        self.inner.into_spectrum(spectrum)
    }

    fn fill_spectrum(&mut self, param: Param) {
        self.inner.fill_spectrum(param)
    }

    fn fill_binary_data_array(&mut self, param: Param) {
        self.inner.fill_binary_data_array(param)
    }

    fn borrow_instrument_configuration(
        mut self,
        instrument_configurations: &'a mut IncrementingIdMap,
    ) -> Self {
        self.inner = self
            .inner
            .borrow_instrument_configuration(instrument_configurations);
        self
    }
}

pub struct MzMLbReaderType<
    C: CentroidPeakAdapting = CentroidPeak,
    D: DeconvolutedPeakAdapting = DeconvolutedPeak,
> {
    pub handle: hdf5::File,
    /// A spectrum ID to byte offset for fast random access
    pub index: OffsetIndex,
    /// The description of the file's contents and the previous data files that were
    /// consumed to produce it.
    pub file_description: FileDescription,
    /// A mapping of different instrument configurations (source, analyzer, detector) components
    /// by ID string.
    pub instrument_configurations: HashMap<u32, InstrumentConfiguration>,
    /// The different software components that were involved in the processing and creation of this
    /// file.
    pub softwares: Vec<Software>,
    /// The data processing and signal transformation operations performed on the raw data in previous
    /// source files to produce this file's contents.
    pub data_processings: Vec<DataProcessing>,

    mzml_parser: MzMLReaderType<ByteReader, C, D>,
    data_buffers: ExternalDataRegistry,
}

impl<'a, 'b: 'a, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MzMLbReaderType<C, D> {
    pub fn new<P: AsRef<Path>>(path: &P) -> io::Result<Self> {
        let handle = match hdf5::File::open(path.as_ref()) {
            Ok(handle) => handle,
            Err(e) => Err(io::Error::new(io::ErrorKind::Other, e))?,
        };

        let index = Self::parse_index(&handle)?;

        let mzml_ds = match handle.dataset("mzML") {
            Ok(ds) => ds,
            Err(e) => Err(io::Error::new(io::ErrorKind::Other, e))?,
        };

        let meta_parser = MzMLReaderType::<ByteReader, C, D>::new(mzml_ds.into());

        let mzml_ds = match handle.dataset("mzML") {
            Ok(ds) => ds,
            Err(e) => Err(io::Error::new(io::ErrorKind::Other, e))?,
        };

        let data_buffers = ExternalDataRegistry::from_hdf5(&handle)?;

        let inst = Self {
            handle,
            index: index,
            file_description: meta_parser.file_description,
            instrument_configurations: meta_parser.instrument_configurations,
            softwares: meta_parser.softwares,
            data_processings: meta_parser.data_processings,
            mzml_parser: MzMLReaderType::<_, C, D>::new(mzml_ds.into()),
            data_buffers: data_buffers,
        };

        Ok(inst)
    }

    fn parse_index(handle: &hdf5::File) -> io::Result<OffsetIndex> {
        let mut index = OffsetIndex::new("spectrum".to_string());
        let mut spectrum_index_ids_ds: ByteReader = match handle.dataset("mzML_spectrumIndex_idRef")
        {
            Ok(ds) => ByteReader::from(ds),
            Err(e) => Err(io::Error::new(io::ErrorKind::Other, e))?,
        };

        let spectrum_index_offsets_ds = match handle.dataset("mzML_spectrumIndex") {
            Ok(ds) => ds,
            Err(e) => Err(io::Error::new(io::ErrorKind::Other, e))?,
        };

        let mut idx_buffer = Vec::new();
        spectrum_index_ids_ds.read_to_end(&mut idx_buffer)?;
        let idx_splits = idx_buffer.split(|b| *b == b'\x00');
        let offsets = match spectrum_index_offsets_ds.read_1d::<u64>() {
            Ok(series) => series,
            Err(e) => Err(io::Error::new(io::ErrorKind::Other, e))?,
        };

        idx_splits.zip(offsets).for_each(|(id, off)| {
            if id.len() == 0 {
                return;
            }
            index.insert(
                String::from_utf8(id.to_vec()).expect("Faild to decode spectrum id as UTF8"),
                off,
            );
        });
        index.init = true;
        Ok(index)
    }

    fn _parse_into<
        B: MzMLSAX
            + DataRegistryBorrower<'a>
            + SpectrumBuilding<'a, C, D, MultiLayerSpectrum<C, D>>
            + 'a,
    >(
        &'a mut self,
        mut accumulator: B,
    ) -> Result<(B, usize), MzMLbError> {
        accumulator = accumulator.borrow_data_registry(&mut self.data_buffers);
        match self.mzml_parser._parse_into(accumulator) {
            Ok(val) => Ok(val),
            Err(e) => Err(e.into()),
        }
    }

    /// Populate a new [`Spectrum`] in-place on the next available spectrum data.
    /// This allocates memory to build the spectrum's attributes but then moves it
    /// into `spectrum` rather than copying it.
    pub fn read_into(
        &mut self,
        spectrum: &mut MultiLayerSpectrum<C, D>,
    ) -> Result<usize, MzMLbError> {
        let accumulator = MzMLbSpectrumBuilder::<C, D>::new();
        match self._parse_into(accumulator) {
            Ok((accumulator, sz)) => {
                accumulator.into_spectrum(spectrum);
                Ok(sz)
            }
            Err(err) => Err(err),
        }
    }

    /// Read the next spectrum directly. Used to implement iteration.
    pub fn read_next(&mut self) -> Option<MultiLayerSpectrum<C, D>> {
        let mut spectrum = MultiLayerSpectrum::<C, D>::default();
        match self.read_into(&mut spectrum) {
            Ok(_sz) => Some(spectrum),
            Err(err) => {
                debug!("Failed to read next spectrum: {err}");
                None
            }
        }
    }

    pub fn get_blosc_available() -> bool {
        filters::blosc_available()
    }

    pub fn get_blosc_nthreads() -> u8 {
        filters::blosc_get_nthreads()
    }

    pub fn set_blosc_nthreads(num_threads: u8) -> u8 {
        filters::blosc_set_nthreads(num_threads)
    }
}

/// [`MzMLReaderType`] instances are [`Iterator`]s over [`Spectrum`]
impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> Iterator for MzMLbReaderType<C, D> {
    type Item = MultiLayerSpectrum<C, D>;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_next()
    }
}

/// They can also be used to fetch specific spectra by ID, index, or start
/// time when the underlying file stream supports [`io::Seek`].
impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    ScanSource<C, D, MultiLayerSpectrum<C, D>> for MzMLbReaderType<C, D>
{
    /// Retrieve a spectrum by it's native ID
    fn get_spectrum_by_id(&mut self, id: &str) -> Option<MultiLayerSpectrum<C, D>> {
        let offset_ref = self.index.get(id);
        let offset = offset_ref.expect("Failed to retrieve offset");
        let start = self
            .mzml_parser
            .stream_position()
            .expect("Failed to save checkpoint");
        self.mzml_parser
            .seek(SeekFrom::Start(offset))
            .expect("Failed to move seek to offset");
        let result = self.read_next();
        self.mzml_parser
            .seek(SeekFrom::Start(start))
            .expect("Failed to restore offset");
        result
    }

    /// Retrieve a spectrum by it's integer index
    fn get_spectrum_by_index(&mut self, index: usize) -> Option<MultiLayerSpectrum<C, D>> {
        let (_id, offset) = self.index.get_index(index)?;
        let byte_offset = offset;
        let start = self
            .mzml_parser
            .stream_position()
            .expect("Failed to save checkpoint");
        self.mzml_parser.seek(SeekFrom::Start(byte_offset)).ok()?;
        let result = self.read_next();
        self.mzml_parser
            .seek(SeekFrom::Start(start))
            .expect("Failed to restore offset");
        result
    }

    /// Return the data stream to the beginning
    fn reset(&mut self) {
        self.mzml_parser
            .seek(SeekFrom::Start(0))
            .expect("Failed to reset file stream");
    }

    fn get_index(&self) -> &OffsetIndex {
        if !self.index.init {
            warn!("Attempting to use an uninitialized offset index on MzMLbReaderType")
        }
        &self.index
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.index = index
    }
}

/// The iterator can also be updated to move to a different location in the
/// stream efficiently.
impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    RandomAccessSpectrumIterator<C, D, MultiLayerSpectrum<C, D>> for MzMLbReaderType<C, D>
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, ScanAccessError> {
        match self._offset_of_id(id) {
            Some(offset) => match self.mzml_parser.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(ScanAccessError::IOError(Some(err))),
            },
            None => Err(ScanAccessError::ScanNotFound),
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, ScanAccessError> {
        match self._offset_of_index(index) {
            Some(offset) => match self.mzml_parser.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(ScanAccessError::IOError(Some(err))),
            },
            None => Err(ScanAccessError::ScanNotFound),
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, ScanAccessError> {
        match self._offset_of_time(time) {
            Some(offset) => match self.mzml_parser.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(ScanAccessError::IOError(Some(err))),
            },
            None => Err(ScanAccessError::ScanNotFound),
        }
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    MZFileReader<C, D, MultiLayerSpectrum<C, D>> for MzMLbReaderType<C, D>
{
    fn open_file(source: fs::File) -> Self {
        let name = filename::file_name(&source).unwrap();
        Self::new(&name).unwrap()
    }

    fn construct_index_from_stream(&mut self) -> u64 {
        self.index.len() as u64
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MSDataFileMetadata
    for MzMLbReaderType<C, D>
{
    crate::impl_metadata_trait!();
}

pub type MzMLbReader = MzMLbReaderType<CentroidPeak, DeconvolutedPeak>;

#[cfg(test)]
mod test {
    use crate::{SpectrumBehavior, MzMLReader};

    use super::*;

    #[test]
    fn test_open() -> io::Result<()> {
        let mut reader = MzMLbReader::new(&"test/data/small.mzMLb")?;
        let mut ref_reader = MzMLReader::open_path("test/data/small.mzML")?;
        assert_eq!(reader.index.len(), 48);
        assert_eq!(reader.softwares.len(), 3);
        let scan = reader.next().unwrap();
        let ref_scan = ref_reader.next().unwrap();
        assert_eq!(scan.ms_level(), ref_scan.ms_level());
        assert_eq!(scan.id(), ref_scan.id());
        let arrays = scan.arrays.as_ref().unwrap();
        let ref_arrays = ref_scan.arrays.as_ref().unwrap();
        dbg!(arrays);
        assert!(arrays.mzs().len() == 19914);
        assert_eq!(arrays.mzs().len(), ref_arrays.mzs().len());
        assert_eq!(arrays.intensities().len(), ref_arrays.intensities().len());

        for (scan, ref_scan) in reader.zip(ref_reader) {
            dbg!(scan.id());
            assert_eq!(scan.ms_level(), ref_scan.ms_level());
            assert_eq!(scan.id(), ref_scan.id());
            let arrays = scan.arrays.as_ref().unwrap();
            let ref_arrays = ref_scan.arrays.as_ref().unwrap();
            dbg!(arrays, ref_arrays);
            assert_eq!(arrays.mzs().len(), ref_arrays.mzs().len());
            assert_eq!(arrays.intensities().len(), ref_arrays.intensities().len());
        }
        Ok(())
    }
}
