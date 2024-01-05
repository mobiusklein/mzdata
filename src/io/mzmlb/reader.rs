use std::borrow::Cow;
use std::collections::HashMap;
use std::io::{self, prelude::*, SeekFrom};
use std::path::Path;
use std::{fs, mem};

#[cfg(feature = "filename")]
use filename;
use hdf5::types::{FixedAscii, FixedUnicode, VarLenAscii, VarLenUnicode};
use hdf5::{self, filters, Dataset, Selection};
use log::{debug, warn};
use ndarray::Ix1;
use thiserror::Error;

use mzpeaks::{CentroidPeak, DeconvolutedPeak};

use crate::io::mzml::{
    CVParamParse, IncrementingIdMap, MzMLParserError, MzMLParserState, MzMLReaderType, MzMLSAX,
    MzMLSpectrumBuilder, ParserResult, SpectrumBuilding,
};
use crate::io::traits::MZFileReader;
use crate::io::utils::DetailLevel;
use crate::io::{OffsetIndex, RandomAccessSpectrumIterator, SpectrumAccessError, ScanSource};
use crate::prelude::{MSDataFileMetadata, ParamLike};

use crate::meta::{DataProcessing, FileDescription, InstrumentConfiguration, Software};
use crate::params::{ControlledVocabulary, Param};
use crate::spectrum::bindata::{
    as_bytes, delta_decoding, linear_prediction_decoding, ArrayRetrievalError,
    BinaryCompressionType, BinaryDataArrayType, ByteArrayView, ByteArrayViewMut, DataArray, BuildFromArrayMap,
};

#[cfg(feature = "numpress")]
use crate::spectrum::bindata::vec_as_bytes;
#[cfg(feature = "numpress")]
use numpress::numpress_decompress;

use crate::spectrum::spectrum::{
    CentroidPeakAdapting, DeconvolutedPeakAdapting, MultiLayerSpectrum,
};
use crate::spectrum::{IsolationWindow, ScanWindow, SelectedIon};

#[derive(Debug, Error)]
pub enum MzMLbError {
    #[error("An HDF5-related error occurred: {0}")]
    HDF5Error(#[from] hdf5::Error),
    #[error("An mzML-related error occurred: {0}")]
    MzMLError(#[from] MzMLParserError),
    #[error("An error occurred while decoding binary data: {0}")]
    ArrayRetrievalError(#[from] ArrayRetrievalError),
}

impl From<MzMLbError> for io::Error {
    fn from(value: MzMLbError) -> Self {
        Self::new(io::ErrorKind::Other, value)
    }
}

#[allow(unused)]
pub(crate) fn is_mzmlb(buf: &[u8]) -> bool {
    const MAGIC_NUMBER: &[u8] = br#"\211HDF\r\n\032\n"#;
    buf.starts_with(MAGIC_NUMBER)
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
        if self.start <= start && end < self.end {
            return true;
        }
        false
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
            chunk_size: Self::default_chunk_size(),
            registry: Default::default(),
            chunk_cache: Default::default(),
        }
    }
}

impl ExternalDataRegistry {
    pub fn new(chunk_size: usize) -> Self {
        Self {
            chunk_size,
            ..Default::default()
        }
    }

    fn default_chunk_size() -> usize {
        2usize.pow(20)
    }

    fn populate_registry(&mut self, handle: &hdf5::File) -> Result<(), MzMLbError> {
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
                            self.registry.insert(name, ds);
                        }
                    };
                });
            }
            Err(err) => return Err(err.into()),
        }
        Ok(())
    }

    #[allow(unused)]
    pub fn from_hdf5(handle: &hdf5::File) -> Result<Self, MzMLbError> {
        let mut storage = Self::default();
        storage.populate_registry(handle)?;
        Ok(storage)
    }

    pub fn from_hdf5_with_chunk_size(
        handle: &hdf5::File,
        chunk_size: usize,
    ) -> Result<Self, MzMLbError> {
        let mut storage = Self::new(chunk_size);
        storage.populate_registry(handle)?;
        Ok(storage)
    }

    fn read_slice_into(
        dataset: &Dataset,
        start: usize,
        end: usize,
        buffer: &mut Vec<u8>,
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

    #[allow(unused)]
    fn handle_encoding(data: &mut DataArray) -> Result<(), ArrayRetrievalError> {
        match data.compression {
            BinaryCompressionType::NoCompression => Ok(()),
            BinaryCompressionType::Zlib => Err(ArrayRetrievalError::DecompressionError(
                data.compression.unsupported_msg(None),
            )),
            #[cfg(feature = "numpress")]
            BinaryCompressionType::NumpressLinear => {
                match data.dtype {
                    BinaryDataArrayType::Float64 => {
                        let buffer = data.coerce::<u8>()?;
                        let buffer = numpress_decompress(&buffer)?;
                        data.data = vec_as_bytes(buffer);
                        data.compression = BinaryCompressionType::Decoded;
                    }
                    _ => {
                        return Err(ArrayRetrievalError::DecompressionError(
                            data.compression.unsupported_msg(Some(
                                format!("Not compatible with {:?}", data.dtype).as_str(),
                            )),
                        ))
                    }
                }
                Err(ArrayRetrievalError::DecompressionError(
                    data.compression.unsupported_msg(None),
                ))
            }
            #[cfg(not(feature = "numpress"))]
            BinaryCompressionType::NumpressLinear => Err(ArrayRetrievalError::DecompressionError(
                data.compression.unsupported_msg(None),
            )),
            BinaryCompressionType::NumpressSLOF => Err(ArrayRetrievalError::DecompressionError(
                data.compression.unsupported_msg(None),
            )),
            BinaryCompressionType::NumpressPIC => Err(ArrayRetrievalError::DecompressionError(
                data.compression.unsupported_msg(None),
            )),
            BinaryCompressionType::NumpressLinearZlib => Err(
                ArrayRetrievalError::DecompressionError(data.compression.unsupported_msg(None)),
            ),
            BinaryCompressionType::NumpressSLOFZlib => Err(
                ArrayRetrievalError::DecompressionError(data.compression.unsupported_msg(None)),
            ),
            BinaryCompressionType::NumpressPICZlib => Err(ArrayRetrievalError::DecompressionError(
                data.compression.unsupported_msg(None),
            )),
            BinaryCompressionType::LinearPrediction => {
                match data.dtype {
                    BinaryDataArrayType::Float64 => {
                        let buffer = data.coerce_mut::<f64>()?;
                        linear_prediction_decoding(buffer);
                        data.compression = BinaryCompressionType::Decoded;
                    }
                    BinaryDataArrayType::Float32 => {
                        let buffer = data.coerce_mut::<f32>()?;
                        linear_prediction_decoding(buffer);
                        data.compression = BinaryCompressionType::Decoded;
                    }
                    _ => {
                        return Err(ArrayRetrievalError::DecompressionError(
                            data.compression.unsupported_msg(Some(
                                format!("Not compatible with {:?}", data.dtype).as_str(),
                            )),
                        ))
                    }
                }
                Ok(())
            }
            BinaryCompressionType::DeltaPrediction => {
                match data.dtype {
                    BinaryDataArrayType::Float64 => {
                        let buffer = data.coerce_mut::<f64>()?;
                        delta_decoding(buffer);
                        data.compression = BinaryCompressionType::Decoded;
                    }
                    BinaryDataArrayType::Float32 => {
                        let buffer = data.coerce_mut::<f32>()?;
                        delta_decoding(buffer);
                        data.compression = BinaryCompressionType::Decoded;
                    }
                    _ => {
                        return Err(ArrayRetrievalError::DecompressionError(
                            data.compression.unsupported_msg(Some(
                                format!("Not compatible with {:?}", data.dtype).as_str(),
                            )),
                        ))
                    }
                }
                Ok(())
            }
            BinaryCompressionType::Decoded => Ok(()),
        }
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
                assert_eq!(destination.data.len(), range_request.length * z);
                return Ok(());
            }
        }
        if let Some(dset) = self.registry.get(&range_request.name).as_ref() {
            let dtype = dset.dtype()?;
            let block_end = (start + self.chunk_size).min(dset.size());
            let cache_block = Self::read_slice_to_bytes(dset, start, block_end)?;
            let cache_block = DataArray::wrap(&destination.name, dtype.into(), cache_block);
            let mut cache_block = CacheInterval::new(start, block_end, cache_block);
            let block = if let Some(cache_block) =
                self.chunk_cache.get_mut(&range_request.name).map(|prev| {
                    mem::swap(prev, &mut cache_block);
                    prev
                }) {
                debug!(
                    "Updated {} cache block: {:?}",
                    range_request.name, cache_block
                );
                cache_block.get(start, end)?
            } else {
                self.chunk_cache
                    .insert(range_request.name.clone(), cache_block);
                let cache_block = self.chunk_cache.get_mut(&range_request.name).unwrap();
                cache_block.get(start, end)?
            };
            destination.data.extend_from_slice(&block);
            assert_eq!(destination.data.len(), range_request.length * z);
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
                self.position = (offset as usize).max(0);
            }
            io::SeekFrom::End(offset) => {
                let n = self.handle.size();
                self.position = (n + offset as usize).max(n);
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
pub struct MzMLbSpectrumBuilder<'a, C: CentroidPeakAdapting + BuildFromArrayMap, D: DeconvolutedPeakAdapting + BuildFromArrayMap> {
    inner: MzMLSpectrumBuilder<'a, C, D>,
    data_registry: Option<&'a mut ExternalDataRegistry>,
    current_data_range_query: DataRangeRequest,
}

impl<'a, C: CentroidPeakAdapting + BuildFromArrayMap, D: DeconvolutedPeakAdapting + BuildFromArrayMap> MzMLbSpectrumBuilder<'a, C, D> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_detail_level(detail_level: DetailLevel) -> Self {
        Self { inner: MzMLSpectrumBuilder::with_detail_level(detail_level), ..Default::default() }
    }
}

impl<'a, C: CentroidPeakAdapting + BuildFromArrayMap, D: DeconvolutedPeakAdapting + BuildFromArrayMap> MzMLSAX
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
                match &state {
                    MzMLParserState::BinaryDataArray => {
                        match MzMLSpectrumBuilder::<'a, C, D>::handle_param(
                            event,
                            reader_position,
                            state,
                        ) {
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
                                                    if self.current_data_range_query.name.is_empty()
                                                        && !param.value.starts_with('/')
                                                    {
                                                        self.current_data_range_query
                                                            .name
                                                            .push('/');
                                                    }
                                                    self.current_data_range_query
                                                        .name
                                                        .push_str(&param.value);
                                                }
                                                // external offset
                                                1002842 => {
                                                    self.current_data_range_query.offset = param
                                                        .value
                                                        .parse()
                                                        .expect("Failed to extract external offset")
                                                }
                                                // external array length
                                                1002843 => self.current_data_range_query.length =
                                                    param.value.parse().expect(
                                                        "Failed to extract external array length",
                                                    ),
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
                    _ => return self.inner.empty_element(event, state, reader_position),
                }
            }
            _ => {}
        }
        Ok(state)
    }

    fn end_element(
        &mut self,
        event: &quick_xml::events::BytesEnd,
        state: MzMLParserState,
    ) -> ParserResult {
        let elt_name = event.name();
        let res = match elt_name.as_ref() {
            b"binaryDataArray" => {
                if self.current_data_range_query.name.is_empty() {
                    return Err(MzMLParserError::IncompleteElementError(
                        "The external data array name was missing or empty".to_owned(),
                        MzMLParserState::BinaryDataArray,
                    ));
                }
                let detail_level = self.inner.detail_level;
                let data_request = mem::take(&mut self.current_data_range_query);
                let array = self.inner.current_array_mut();
                if !matches!(detail_level, DetailLevel::MetadataOnly) {
                    self.data_registry
                        .as_mut()
                        .expect("Did not provide data registry")
                        .get(&data_request, array)
                } else {
                    Ok(())
                }
            }
            _ => Ok(()),
        };
        match res {
            Ok(()) => {}
            Err(e) => match e {
                MzMLbError::HDF5Error(e) => match e {
                    hdf5::Error::HDF5(_) => Err(MzMLParserError::IOError(
                        state,
                        io::Error::new(io::ErrorKind::Other, e),
                    ))?,
                    hdf5::Error::Internal(e) => Err(MzMLParserError::IOError(
                        state,
                        io::Error::new(io::ErrorKind::Other, e),
                    ))?,
                },
                MzMLbError::MzMLError(e) => Err(e)?,
                MzMLbError::ArrayRetrievalError(e) => {
                    Err(MzMLParserError::IOError(state, e.into()))?
                }
            },
        }
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

impl<'a, C: CentroidPeakAdapting + BuildFromArrayMap, D: DeconvolutedPeakAdapting + BuildFromArrayMap> DataRegistryBorrower<'a>
    for MzMLbSpectrumBuilder<'a, C, D>
{
    fn borrow_data_registry(mut self, data_registry: &'a mut ExternalDataRegistry) -> Self {
        self.data_registry = Some(data_registry);
        self
    }
}

impl<'a, C: CentroidPeakAdapting + BuildFromArrayMap, D: DeconvolutedPeakAdapting + BuildFromArrayMap>
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

    fn fill_spectrum<P: ParamLike + Into<Param>>(&mut self, param: P) {
        self.inner.fill_spectrum(param)
    }

    fn fill_binary_data_array<P: ParamLike + Into<Param>>(&mut self, param: P) {
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

    pub schema_version: String,

    pub detail_level: DetailLevel,

    mzml_parser: MzMLReaderType<ByteReader, C, D>,
    data_buffers: ExternalDataRegistry,
}

impl<'a, 'b: 'a, C: CentroidPeakAdapting + BuildFromArrayMap, D: DeconvolutedPeakAdapting + BuildFromArrayMap> MzMLbReaderType<C, D> {
    /// Create a new `[MzMLbReader]` with an internal cache size of `chunk_size` elements
    /// per data array to reduce the number of disk reads needed to populate spectra, and
    /// sets the `[DetailLevel]`.
    pub fn with_chunk_size_and_detail_level<P: AsRef<Path>>(
        path: &P,
        chunk_size: usize,
        detail_level: DetailLevel,
    ) -> io::Result<Self> {
        let handle = match hdf5::File::open(path.as_ref()) {
            Ok(handle) => handle,
            Err(e) => Err(io::Error::new(io::ErrorKind::Other, e))?,
        };

        let index = Self::parse_spectrum_index(&handle)?;

        let mzml_ds = match handle.dataset("mzML") {
            Ok(ds) => ds,
            Err(e) => Err(io::Error::new(io::ErrorKind::Other, e))?,
        };

        let schema_version = match Self::read_schema_version(&mzml_ds) {
            Ok(schema_version) => schema_version,
            Err(e) => Err(io::Error::new(io::ErrorKind::Other, e))?,
        };

        let mut mzml_parser = MzMLReaderType::<ByteReader, C, D>::new(mzml_ds.into());
        mzml_parser.seek(SeekFrom::Start(0))?;

        let data_buffers = ExternalDataRegistry::from_hdf5_with_chunk_size(&handle, chunk_size)?;

        let inst = Self {
            handle,
            index,
            file_description: mzml_parser.file_description.clone(),
            instrument_configurations: mzml_parser.instrument_configurations.clone(),
            softwares: mzml_parser.softwares.clone(),
            data_processings: mzml_parser.data_processings.clone(),
            detail_level,
            mzml_parser,
            schema_version,
            data_buffers,
        };

        Ok(inst)
    }

    /// Create a new `[MzMLbReader]` with an internal cache size of `chunk_size` elements
    /// per data array to reduce the number of disk reads needed to populate spectra.
    ///
    /// The default chunk size is 2**20 elements, which can use as much as 8.4 MB for 64-bit
    /// data types like those used to store m/z.
    pub fn with_chunk_size<P: AsRef<Path>>(path: &P, chunk_size: usize) -> io::Result<Self> {
        Self::with_chunk_size_and_detail_level(path, chunk_size, DetailLevel::Full)
    }

    fn read_schema_version(mzml_ds: &Dataset) -> Result<String, hdf5::Error> {
        let schema_version_attr = mzml_ds.attr("version")?;
        let dtype = schema_version_attr.dtype()?;
        let td = dtype.to_descriptor()?;
        let buf = match td {
            hdf5::types::TypeDescriptor::FixedAscii(_) => {
                let val: FixedAscii<9> = schema_version_attr.read_scalar()?;
                val.to_string()
            }
            hdf5::types::TypeDescriptor::FixedUnicode(_) => {
                let val: FixedUnicode<9> = schema_version_attr.read_scalar()?;
                val.to_string()
            }
            hdf5::types::TypeDescriptor::VarLenAscii => {
                let val: VarLenAscii = schema_version_attr.read_scalar()?;
                val.to_string()
            }
            hdf5::types::TypeDescriptor::VarLenUnicode => {
                let val: VarLenUnicode = schema_version_attr.read_scalar()?;
                val.to_string()
            }
            _ => {
                let val: [u8; 9] = schema_version_attr.as_reader().read_scalar()?;
                String::from_utf8_lossy(&val).to_string()
            }
        };
        debug!("Parsed mzMLb version: {}", buf);
        Ok(buf)
    }

    /// Create a new `[MzMLbReader]` with the default caching behavior.
    pub fn new<P: AsRef<Path>>(path: &P) -> io::Result<Self> {
        Self::with_chunk_size(path, ExternalDataRegistry::default_chunk_size())
    }

    /// Parses the regular spectrum index if it is present.
    fn parse_spectrum_index(handle: &hdf5::File) -> io::Result<OffsetIndex> {
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
            if id.is_empty() || off == 0 {
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

    /// Populate a new [`MultiLayerSpectrum`] in-place on the next available spectrum data.
    /// This allocates memory to build the spectrum's attributes but then moves it
    /// into `spectrum` rather than copying it.
    pub fn read_into(
        &mut self,
        spectrum: &mut MultiLayerSpectrum<C, D>,
    ) -> Result<usize, MzMLbError> {
        let accumulator = MzMLbSpectrumBuilder::<C, D>::with_detail_level(self.detail_level);
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

/// [`MzMLReaderType`] instances are [`Iterator`]s over [`MultiLayerSpectrum`]
impl<C: CentroidPeakAdapting + BuildFromArrayMap, D: DeconvolutedPeakAdapting + BuildFromArrayMap> Iterator for MzMLbReaderType<C, D> {
    type Item = MultiLayerSpectrum<C, D>;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_next()
    }
}

/// They can also be used to fetch specific spectra by ID, index, or start
/// time when the underlying file stream supports [`io::Seek`].
impl<C: CentroidPeakAdapting + BuildFromArrayMap, D: DeconvolutedPeakAdapting + BuildFromArrayMap>
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
            .expect("Failed to seek to offset");
        debug_assert!(
            self.mzml_parser.check_stream("spectrum").unwrap(),
            "The next XML tag was not `spectrum`"
        );
        self.mzml_parser.state = MzMLParserState::Resume;
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
        self.mzml_parser
            .seek(SeekFrom::Start(byte_offset))
            .expect("Failed to seek to offset");
        debug_assert!(
            self.mzml_parser.check_stream("spectrum").unwrap(),
            "The next XML tag was not `spectrum`"
        );
        self.mzml_parser.state = MzMLParserState::Resume;
        let result = self.read_next();
        self.mzml_parser
            .seek(SeekFrom::Start(start))
            .expect("Failed to restore offset");
        result
    }

    /// Return the data stream to the beginning
    fn reset(&mut self) {
        self.mzml_parser.state = MzMLParserState::Resume;
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
impl<C: CentroidPeakAdapting + BuildFromArrayMap, D: DeconvolutedPeakAdapting + BuildFromArrayMap>
    RandomAccessSpectrumIterator<C, D, MultiLayerSpectrum<C, D>> for MzMLbReaderType<C, D>
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, SpectrumAccessError> {
        match self._offset_of_id(id) {
            Some(offset) => match self.mzml_parser.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(SpectrumAccessError::IOError(Some(err))),
            },
            None => Err(SpectrumAccessError::SpectrumIdNotFound(id.to_string())),
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, SpectrumAccessError> {
        match self._offset_of_index(index) {
            Some(offset) => match self.mzml_parser.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(SpectrumAccessError::IOError(Some(err))),
            },
            None => Err(SpectrumAccessError::SpectrumIndexNotFound(index)),
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, SpectrumAccessError> {
        match self._offset_of_time(time) {
            Some(offset) => match self.mzml_parser.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(SpectrumAccessError::IOError(Some(err))),
            },
            None => Err(SpectrumAccessError::SpectrumNotFound),
        }
    }
}

impl<C: CentroidPeakAdapting + BuildFromArrayMap, D: DeconvolutedPeakAdapting + BuildFromArrayMap>
    MZFileReader<C, D, MultiLayerSpectrum<C, D>> for MzMLbReaderType<C, D>
{
    #[allow(unused)]
    /// The underlying HDF5 library for Rust, [`hdf5`](https://docs.rs/hdf5/latest/hdf5/) doesn't
    /// support reading directly from Rust [`io::Read`]-implementing objects yet. Without a means
    /// of retrieving a [`path::Path`]-like value from a file handle, with the [`filename`](https://docs.rs/filename/latest/filename/)
    /// library, this method **panics**. Enable this extra feature if you would like this method to
    /// work, but it is reported to have obscure compilation errors.
    fn open_file(source: fs::File) -> Self {
        #[cfg(feature = "filename")]
        {
            let name = filename::file_name(&source).unwrap();
            Self::new(&name).unwrap()
        }
        #[cfg(not(feature = "filename"))]
        panic!("Cannot read an mzMLb file from an open file handle without the `filename` crate")
    }

    fn construct_index_from_stream(&mut self) -> u64 {
        self.index.len() as u64
    }

    fn open_path<P>(path: P) -> io::Result<Self>
    where
        P: Into<std::path::PathBuf> + Clone,
    {
        Self::new(&path.into())
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MSDataFileMetadata
    for MzMLbReaderType<C, D>
{
    crate::impl_metadata_trait!();

    fn spectrum_count_hint(&self) -> Option<u64> {
        Some(self.index.len() as u64)
    }
}

pub type MzMLbReader = MzMLbReaderType<CentroidPeak, DeconvolutedPeak>;

#[cfg(test)]
mod test {
    use crate::{MzMLReader, SpectrumLike};

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
        assert!(arrays.mzs()?.len() == 19914);
        assert_eq!(arrays.mzs()?.len(), ref_arrays.mzs()?.len());
        assert_eq!(arrays.intensities()?.len(), ref_arrays.intensities()?.len());

        for (scan, ref_scan) in reader.zip(ref_reader) {
            assert_eq!(scan.ms_level(), ref_scan.ms_level());
            assert_eq!(scan.id(), ref_scan.id());
            let arrays = scan.arrays.as_ref().unwrap();
            let ref_arrays = ref_scan.arrays.as_ref().unwrap();
            assert_eq!(arrays.mzs()?.len(), ref_arrays.mzs()?.len());
            assert_eq!(arrays.intensities()?.len(), ref_arrays.intensities()?.len());
        }
        Ok(())
    }
}
