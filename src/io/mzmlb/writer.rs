use std::borrow::Cow;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::ffi::CString;
use std::io::{self, prelude::*, Cursor};
use std::marker::PhantomData;
use std::path::{Path, PathBuf};

use bytemuck;

use hdf5::from_id;
use hdf5::globals::H5T_C_S1;
use hdf5::types::{FixedAscii, IntSize, TypeDescriptor};
use hdf5::{self, filters, Attribute, Dataspace, Extents};
use hdf5::{h5lock, h5try};
use hdf5_sys;
use hdf5_sys::h5i::hid_t;
use hdf5_sys::h5t::{H5T_cset_t, H5T_str_t, H5Tcopy, H5Tset_cset, H5Tset_size, H5Tset_strpad};

use ndarray::Array1;

use mzpeaks::{CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak};
use quick_xml::events::{BytesStart, Event};
use thiserror::Error;

use crate::io::traits::ScanWriter;
use crate::meta::{
    DataProcessing, FileDescription, InstrumentConfiguration, MassSpectrometryRun, Software,
};
use crate::params::{ControlledVocabulary, ParamDescribed};
use crate::prelude::{MSDataFileMetadata, SpectrumLike};
use crate::spectrum::bindata::{
    ArrayRetrievalError, BinaryDataArrayType, BuildArrayMapFrom, ByteArrayView, DataArray,
};
use crate::spectrum::{ArrayType, BinaryArrayMap, Chromatogram, ChromatogramLike, PeakDataLevel};

use crate::io::mzml::{MzMLWriterError, MzMLWriterState, MzMLWriterType};

macro_rules! bstart {
    ($e:tt) => {
        BytesStart::from_content($e, $e.len())
    };
}

macro_rules! attrib {
    ($name:expr, $value:expr, $elt:ident) => {
        let key = $name.as_bytes();
        let value = $value.as_bytes();
        $elt.push_attribute((key, value));
    };
}

macro_rules! start_event {
    ($writer:ident, $target:ident) => {
        $writer
            .mzml_writer
            .write_event(Event::Start($target.borrow()))?;
    };
}

macro_rules! end_event {
    ($writer:ident, $target:ident) => {
        $writer
            .mzml_writer
            .write_event(Event::End($target.to_end()))?;
    };
}

const DEFAULT_CHUNK_SIZE: usize = 1000000;

#[derive(Debug, Error)]
pub enum MzMLbWriterError {
    #[error("An HDF5-related error occurred: {0}")]
    HDF5Error(#[from] hdf5::Error),
    #[error("An mzML-related error occurred: {0}")]
    MzMLError(#[from] MzMLWriterError),
    #[error("An error occurred while writing: {0}")]
    IOError(#[from] io::Error),
    #[error("An error occured while manipulating binary data: {0}")]
    ArrayRetrievalError(#[from] ArrayRetrievalError),
}

impl From<MzMLbWriterError> for io::Error {
    fn from(value: MzMLbWriterError) -> Self {
        Self::new(io::ErrorKind::Other, value)
    }
}

impl<C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default> MSDataFileMetadata
    for MzMLbWriterType<C, D>
where
    C: BuildArrayMapFrom,
    D: BuildArrayMapFrom,
{
    fn data_processings(&self) -> &Vec<DataProcessing> {
        self.mzml_writer.data_processings()
    }

    fn instrument_configurations(&self) -> &HashMap<u32, InstrumentConfiguration> {
        self.mzml_writer.instrument_configurations()
    }

    fn file_description(&self) -> &FileDescription {
        self.mzml_writer.file_description()
    }

    fn softwares(&self) -> &Vec<Software> {
        self.mzml_writer.softwares()
    }

    fn data_processings_mut(&mut self) -> &mut Vec<DataProcessing> {
        self.mzml_writer.data_processings_mut()
    }

    fn instrument_configurations_mut(&mut self) -> &mut HashMap<u32, InstrumentConfiguration> {
        self.mzml_writer.instrument_configurations_mut()
    }

    fn file_description_mut(&mut self) -> &mut FileDescription {
        self.mzml_writer.file_description_mut()
    }

    fn softwares_mut(&mut self) -> &mut Vec<Software> {
        self.mzml_writer.softwares_mut()
    }

    fn copy_metadata_from<T: MSDataFileMetadata>(&mut self, source: &T) {
        self.mzml_writer.copy_metadata_from(source)
    }

    fn run_description(&self) -> Option<&MassSpectrometryRun> {
        Some(&self.mzml_writer.run)
    }

    fn run_description_mut(&mut self) -> Option<&mut MassSpectrometryRun> {
        Some(&mut self.mzml_writer.run)
    }
}

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
pub enum BufferContext {
    Spectrum,
    Chromatogram,
}

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub struct BufferName {
    context: BufferContext,
    array_type: ArrayType,
    dtype: BinaryDataArrayType,
}

impl BufferName {
    pub fn new(context: BufferContext, array_type: ArrayType, dtype: BinaryDataArrayType) -> Self {
        Self {
            context,
            array_type,
            dtype,
        }
    }
}

impl ToString for BufferName {
    fn to_string(&self) -> String {
        let context = match self.context {
            BufferContext::Spectrum => "spectrum",
            BufferContext::Chromatogram => "chromatogram",
        };
        let tp_name = match &self.array_type {
            ArrayType::Unknown => Cow::Borrowed("unknown"),
            ArrayType::MZArray => Cow::Borrowed("mz"),
            ArrayType::IntensityArray => Cow::Borrowed("intensity"),
            ArrayType::ChargeArray => Cow::Borrowed("charge"),
            ArrayType::SignalToNoiseArray => Cow::Borrowed("snr"),
            ArrayType::TimeArray => Cow::Borrowed("time"),
            ArrayType::WavelengthArray => Cow::Borrowed("wavelength"),
            ArrayType::IonMobilityArray => Cow::Borrowed("ion_mobility"),
            ArrayType::MeanIonMobilityArray => Cow::Borrowed("mean_ion_mobility"),
            ArrayType::RawIonMobilityArray => Cow::Borrowed("raw_ion_mobility"),
            ArrayType::DeconvolutedIonMobilityArray => Cow::Borrowed("deconvoluted_ion_mobility"),
            ArrayType::NonStandardDataArray { name } => Cow::Owned(name.replace(['/', ' '], "_")),
        };
        let dtype = match self.dtype {
            BinaryDataArrayType::Unknown => "unknown",
            BinaryDataArrayType::Float64 => "f64",
            BinaryDataArrayType::Float32 => "f32",
            BinaryDataArrayType::Int64 => "i64",
            BinaryDataArrayType::Int32 => "i32",
            BinaryDataArrayType::ASCII => "u8",
        };
        format!("{}_{}_{}", context, tp_name, dtype)
    }
}

#[derive(Debug)]
struct BinaryDataArrayBuffer {
    pub buffer: Cursor<Vec<u8>>,
    pub dataset: hdf5::Dataset,
    pub dtype: BinaryDataArrayType,
    pub chunk_size: usize,
    pub offset: usize,
}

impl BinaryDataArrayBuffer {
    pub fn create_dataset(
        name: BufferName,
        builder: hdf5::DatasetBuilder,
        chunk_size: usize,
        filters: &[hdf5::filters::Filter],
        dtype: BinaryDataArrayType,
    ) -> Result<Self, MzMLbWriterError> {
        let td = TypeDescriptor::from(&dtype);
        let builder = builder
            .set_filters(filters)
            .empty_as(&td)
            .chunk((chunk_size,))
            .shape((1..,));
        let dataset: hdf5::Dataset = builder.create(name.to_string().as_str())?;
        let inst = Self {
            buffer: Cursor::new(Vec::with_capacity((chunk_size as f64 * 1.25) as usize)),
            dataset,
            dtype,
            chunk_size,
            offset: 0,
        };
        Ok(inst)
    }

    pub fn flush(&mut self, force: bool) -> WriterResult {
        let size = self.buffer.position();
        if size >= self.chunk_size as u64 || (force && size > 0) {
            self.buffer.seek(io::SeekFrom::Start(0))?;
            let chunk = &self.buffer.get_ref()[..size as usize];
            match self.dtype {
                BinaryDataArrayType::Unknown => {
                    let start = self.offset;
                    let end = start + chunk.len();
                    self.dataset.resize((end,))?;
                    self.dataset.as_writer().write_slice(chunk, start..end)?;
                    self.offset = end;
                }
                BinaryDataArrayType::Float64 => {
                    let chunk: &[f64] = bytemuck::cast_slice(chunk);
                    let start = self.offset;
                    let end = start + chunk.len();
                    self.dataset.resize((end,))?;
                    self.dataset.as_writer().write_slice(chunk, start..end)?;
                    self.offset = end;
                }
                BinaryDataArrayType::Float32 => {
                    let chunk: &[f32] = bytemuck::cast_slice(chunk);
                    let start = self.offset;
                    let end = start + chunk.len();
                    self.dataset.resize((end,))?;
                    self.dataset.as_writer().write_slice(chunk, start..end)?;
                    self.offset = end;
                }
                BinaryDataArrayType::Int64 => {
                    let chunk: &[i64] = bytemuck::cast_slice(chunk);
                    let start = self.offset;
                    let end = start + chunk.len();
                    self.dataset.resize((end,))?;
                    self.dataset.as_writer().write_slice(chunk, start..end)?;
                    self.offset = end;
                }
                BinaryDataArrayType::Int32 => {
                    let chunk: &[i32] = bytemuck::cast_slice(chunk);
                    let start = self.offset;
                    let end = start + chunk.len();
                    self.dataset.resize((end,))?;
                    self.dataset.as_writer().write_slice(chunk, start..end)?;
                    self.offset = end;
                }
                BinaryDataArrayType::ASCII => {
                    let start = self.offset;
                    let end = start + chunk.len();
                    self.dataset.resize((end,))?;
                    self.dataset.as_writer().write_slice(chunk, start..end)?;
                    self.offset = end;
                }
            };
            self.buffer.set_position(0);
            if log::log_enabled!(log::Level::Debug) {
                log::debug!(
                    "Flush {}:{}:{} {}/{} {:0.3}",
                    self.dataset.name(),
                    self.dtype,
                    self.offset,
                    self.buffer.position(),
                    self.buffer.get_ref().capacity(),
                    (self.buffer.position() as f64 / self.buffer.get_ref().capacity() as f64)
                        * 100.0
                );
            }
        };
        Ok(())
    }

    pub fn add(&mut self, data: &DataArray) -> Result<u64, MzMLbWriterError> {
        let offset = self.offset as u64 + self.buffer.position() / self.dtype.size_of() as u64;
        self.buffer.write_all(&data.data)?;
        self.flush(false)?;
        Ok(offset)
    }
}

#[derive(Debug)]
struct ByteWriter {
    dataset: hdf5::Dataset,
    buffer: io::Cursor<Vec<u8>>,
    chunk_size: usize,
    offset: usize,
}

impl ByteWriter {
    fn new(
        dataset: hdf5::Dataset,
        chunk_size: usize,
    ) -> Self {
        let offset = 0;
        let buffer = io::Cursor::new(Vec::with_capacity((chunk_size as f64 * 1.25) as usize));
        Self {
            dataset,
            buffer,
            chunk_size,
            offset,
        }
    }

    fn write_to_dataset(&mut self) -> Result<(), hdf5::Error> {
        let size = self.buffer.position();
        if size == 0 {
            return Ok(())
        }
        let start = self.offset;
        let chunk = &self.buffer.get_ref()[..size as usize];
        let end = start + chunk.len();
        self.dataset.resize((end,))?;
        self.dataset.as_writer().write_slice(chunk, start..end)?;
        self.offset = end;
        if log::log_enabled!(log::Level::Debug) {
            log::debug!(
                "Flush {}:{} {}/{} {:0.3}",
                self.dataset.name(),
                self.offset,
                self.buffer.position(),
                self.buffer.get_ref().capacity(),
                (self.buffer.position() as f64 / self.buffer.get_ref().capacity() as f64)
                    * 100.0
            );
        }
        self.buffer.set_position(0);
        Ok(())
    }

    fn resize_to_fit(&mut self) -> Result<(), hdf5::Error> {
        self.write_to_dataset()?;
        self.dataset.resize((self.offset, ))
    }
}


impl io::Write for ByteWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let z = self.buffer.write(buf)?;
        if self.buffer.position() > self.chunk_size as u64 {
            if let Err(e) = self.write_to_dataset() {
                return Err(io::Error::new(io::ErrorKind::Other, e))
            };
        }
        Ok(z)
    }

    fn flush(&mut self) -> io::Result<()> {
        if let Err(e) = self.write_to_dataset() {
            return Err(io::Error::new(io::ErrorKind::Other, e))
        };
        Ok(())
    }
}


pub type WriterResult = Result<(), MzMLbWriterError>;

impl<'a, C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default> ScanWriter<'a, C, D>
    for MzMLbWriterType<C, D>
where
    C: BuildArrayMapFrom,
    D: BuildArrayMapFrom,
{
    fn write<S: SpectrumLike<C, D> + 'static>(&mut self, spectrum: &S) -> io::Result<usize> {
        match self.write_spectrum(spectrum) {
            Ok(()) => {
                let pos = self.mzml_writer.stream_position()?;
                Ok(pos as usize)
            }
            Err(err) => Err(err.into()),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        self.mzml_writer.flush()?;
        Ok(())
    }

    fn close(&mut self) -> io::Result<()> {
        self.close()?;
        Ok(())
    }
}

#[derive(Debug, Default)]
pub struct MzMLbWriterBuilder<
    C: CentroidLike + Default + 'static,
    D: DeconvolutedCentroidLike + Default + 'static,
> where
    C: BuildArrayMapFrom,
    D: BuildArrayMapFrom,
{
    path: PathBuf,
    chunk_size: Option<usize>,
    filters: Option<Vec<filters::Filter>>,
    _c: PhantomData<C>,
    _d: PhantomData<D>,
}

impl<C: CentroidLike + Default + 'static, D: DeconvolutedCentroidLike + Default + 'static>
    MzMLbWriterBuilder<C, D>
where
    C: BuildArrayMapFrom,
    D: BuildArrayMapFrom,
{
    pub fn new<P: Into<PathBuf>>(path: P) -> Self {
        Self {
            path: path.into(),
            ..Default::default()
        }
    }

    pub fn with_chunk_size(mut self, chunk_size: usize) -> Self {
        self.chunk_size = Some(chunk_size);
        self
    }

    pub fn with_zlib_compression(mut self, compression_level: u8) -> Self {
        self.filters = Some(MzMLbWriterType::<C, D>::zlib_compression(compression_level));
        self
    }

    pub fn with_blosc_zstd_compression(mut self, compression_level: u8) -> Self {
        self.filters = Some(MzMLbWriterType::<C, D>::blosc_zstd(compression_level));
        self
    }

    pub fn create(self) -> io::Result<MzMLbWriterType<C, D>> {
        let mut inst = MzMLbWriterType::new_with_chunk_size(
            &self.path,
            self.chunk_size.unwrap_or(DEFAULT_CHUNK_SIZE),
        )?;
        inst.filters = self
            .filters
            .unwrap_or(MzMLbWriterType::<C, D>::zlib_compression(9));
        Ok(inst)
    }
}

impl<C: CentroidLike + Default + 'static, D: DeconvolutedCentroidLike + Default + 'static>
    From<MzMLbWriterBuilder<C, D>> for MzMLbWriterType<C, D>
where
    C: BuildArrayMapFrom,
    D: BuildArrayMapFrom,
{
    fn from(value: MzMLbWriterBuilder<C, D>) -> Self {
        value.create().unwrap()
    }
}

#[derive(Debug)]
pub struct MzMLbWriterType<
    C: CentroidLike + Default + 'static,
    D: DeconvolutedCentroidLike + Default + 'static,
> where
    C: BuildArrayMapFrom,
    D: BuildArrayMapFrom,
{
    handle: hdf5::File,
    mzml_writer: MzMLWriterType<ByteWriter, C, D>,
    buffers: HashMap<BufferName, BinaryDataArrayBuffer>,
    chunk_size: usize,
    filters: Vec<filters::Filter>,
}

impl<C: CentroidLike + Default + 'static, D: DeconvolutedCentroidLike + Default + 'static>
    MzMLbWriterType<C, D>
where
    C: BuildArrayMapFrom,
    D: BuildArrayMapFrom,
{
    pub fn new<P: AsRef<Path>>(path: &P) -> io::Result<Self> {
        Self::new_with_chunk_size(path, DEFAULT_CHUNK_SIZE)
    }

    pub fn new_with_chunk_size_and_filters<P: AsRef<Path>>(path: &P, chunk_size: usize, filters: Vec<hdf5::filters::Filter>) -> io::Result<Self> {
        let mut handle = match hdf5::File::create(path) {
            Ok(handle) => handle,
            Err(e) => return Err(MzMLbWriterError::HDF5Error(e).into()),
        };

        let buffer = match Self::make_mzml_buffer(&mut handle, chunk_size, &filters) {
            Ok(buffer) => {
                ByteWriter::new(buffer, chunk_size)
            },
            Err(e) => {
                return Err(e.into())
            },
        };

        let mzml_writer: MzMLWriterType<ByteWriter, C, D> =
            MzMLWriterType::new_with_index(buffer, false);

        Ok(Self {
            handle,
            mzml_writer,
            buffers: HashMap::new(),
            chunk_size,
            filters: Self::zlib_compression(9),
        })
    }

    pub fn new_with_chunk_size<P: AsRef<Path>>(path: &P, chunk_size: usize) -> io::Result<Self> {
        let filters = Self::zlib_compression(9);
        Self::new_with_chunk_size_and_filters(path, chunk_size, filters)
    }

    fn get_ms_cv(&self) -> &ControlledVocabulary {
        self.mzml_writer.get_ms_cv()
    }

    fn write_binary_data_array(
        &mut self,
        array: &DataArray,
        context: BufferContext,
        default_array_size: usize,
    ) -> WriterResult {
        let mut outer = bstart!("binaryDataArray");

        let size = array.data_len()?;
        let size_str = size.to_string();
        attrib!("encodedLength", "0", outer);
        if size != default_array_size {
            attrib!("arrayLength", size_str, outer);
        }

        start_event!(self, outer);
        match &array.dtype {
            BinaryDataArrayType::Float32 => self
                .mzml_writer
                .write_param(&self.get_ms_cv().param("MS:1000521", "32-bit float"))?,
            BinaryDataArrayType::Float64 => self
                .mzml_writer
                .write_param(&self.get_ms_cv().param("MS:1000523", "64-bit float"))?,
            BinaryDataArrayType::Int32 => self
                .mzml_writer
                .write_param(&self.get_ms_cv().param("MS:1000519", "32-bit integer"))?,
            BinaryDataArrayType::Int64 => self
                .mzml_writer
                .write_param(&self.get_ms_cv().param("MS:1000522", "64-bit integer"))?,
            BinaryDataArrayType::ASCII => self.mzml_writer.write_param(
                &self
                    .get_ms_cv()
                    .param("MS:1001479 ", "null-terminated ASCII string"),
            )?,
            _ => {
                panic!(
                    "Could not determine data type for binary data array. Found {:?}",
                    array.dtype
                )
            }
        }
        match &array.name {
            ArrayType::MZArray => self.mzml_writer.write_param(&array.name.as_param_const())?,
            ArrayType::IntensityArray => {
                self.mzml_writer.write_param(&array.name.as_param_const())?
            }
            ArrayType::ChargeArray => self.mzml_writer.write_param(&array.name.as_param_const())?,
            ArrayType::TimeArray => self.mzml_writer.write_param(&array.name.as_param_const())?,
            ArrayType::RawIonMobilityArray => self
                .mzml_writer
                .write_param(&array.name.as_param_with_unit_const(array.unit))?,
            ArrayType::MeanIonMobilityArray => self
                .mzml_writer
                .write_param(&array.name.as_param_with_unit_const(array.unit))?,
            ArrayType::DeconvolutedIonMobilityArray => self
                .mzml_writer
                .write_param(&array.name.as_param_with_unit_const(array.unit))?,
            ArrayType::SignalToNoiseArray => {
                self.mzml_writer.write_param(&array.name.as_param_const())?
            }
            ArrayType::NonStandardDataArray { name } => {
                let mut p =
                    self.get_ms_cv()
                        .param_val("MS:1000786", "non-standard data array", name);
                p = p.with_unit_t(&array.unit);
                self.mzml_writer.write_param(&p)?;
            }
            _ => {
                panic!("Could not determine how to name for {:?}", array.name);
            }
        }

        let key = BufferName::new(context, array.name.clone(), array.dtype);
        let dset_name = key.to_string();

        self.mzml_writer
            .write_param(&self.mzml_writer.get_ms_cv().param_val(
                "MS:1002841",
                "external HDF5 dataset",
                dset_name,
            ))?;
        self.mzml_writer
            .write_param(&self.mzml_writer.get_ms_cv().param_val(
                "MS:1002843",
                "external array length",
                size_str,
            ))?;
        let offset = match self.buffers.entry(key) {
            Entry::Occupied(mut buf) => buf.get_mut().add(array)?,
            Entry::Vacant(opening) => {
                let builder = self.handle.new_dataset_builder();
                let key2 = BufferName::new(context, array.name.clone(), array.dtype);
                let mut buf = BinaryDataArrayBuffer::create_dataset(
                    key2,
                    builder,
                    self.chunk_size,
                    &self.filters,
                    array.dtype,
                )?;
                let offset = buf.add(array)?;
                opening.insert(buf);
                offset
            }
        };
        self.mzml_writer
            .write_param(&self.mzml_writer.get_ms_cv().param_val(
                "MS:1002842",
                "external offset",
                offset.to_string(),
            ))?;
        let bin = bstart!("binary");
        start_event!(self, bin);
        end_event!(self, bin);
        end_event!(self, outer);
        Ok(())
    }

    fn write_binary_data_arrays(
        &mut self,
        arrays: &BinaryArrayMap,
        context: BufferContext,
        default_array_size: usize,
    ) -> WriterResult {
        let count = arrays.len().to_string();
        let mut outer = bstart!("binaryDataArrayList");
        attrib!("count", count, outer);
        start_event!(self, outer);
        let mut array_pairs: Vec<(&ArrayType, &DataArray)> = arrays.iter().collect();
        array_pairs.sort_by_key(|f| f.0);
        for (_tp, array) in array_pairs {
            self.write_binary_data_array(array, context, default_array_size)?
        }
        end_event!(self, outer);
        Ok(())
    }

    pub fn write_spectrum<S: SpectrumLike<C, D> + 'static>(
        &mut self,
        spectrum: &S,
    ) -> WriterResult {
        match self.mzml_writer.state {
            MzMLWriterState::SpectrumList => {}
            state if state < MzMLWriterState::SpectrumList => {
                self.mzml_writer.start_spectrum_list()?;
            }
            state if state > MzMLWriterState::SpectrumList => {
                // Cannot write spectrum, currently in state which happens
                // after spectra may be written
                return Err(MzMLWriterError::InvalidActionError(self.mzml_writer.state).into());
            }
            _ => {}
        }
        let pos = self.mzml_writer.stream_position()?;
        self.mzml_writer
            .spectrum_offset_index
            .insert(spectrum.id().to_string(), pos);

        let mut outer = bstart!("spectrum");
        let default_array_size = self.mzml_writer.start_spectrum(spectrum, &mut outer)?;

        self.mzml_writer.write_spectrum_descriptors(spectrum)?;

        self.mzml_writer
            .tic_collector
            .add(spectrum.start_time(), spectrum.peaks().tic());
        self.mzml_writer.bic_collector.add(
            spectrum.start_time(),
            spectrum.peaks().base_peak().intensity,
        );

        match spectrum.peaks() {
            PeakDataLevel::RawData(arrays) => {
                self.write_binary_data_arrays(arrays, BufferContext::Spectrum, default_array_size)?
            }
            PeakDataLevel::Centroid(arrays) => self.write_binary_data_arrays(
                &C::as_arrays(&arrays[0..]),
                BufferContext::Spectrum,
                default_array_size,
            )?,
            PeakDataLevel::Deconvoluted(arrays) => self.write_binary_data_arrays(
                &D::as_arrays(&arrays[0..]),
                BufferContext::Spectrum,
                default_array_size,
            )?,
            PeakDataLevel::Missing => todo!(),
        }

        end_event!(self, outer);
        Ok(())
    }

    pub fn write_chromatogram(&mut self, chromatogram: &Chromatogram) -> WriterResult {
        match self.mzml_writer.state {
            MzMLWriterState::ChromatogramList => {}
            state if state < MzMLWriterState::ChromatogramList => {
                self.mzml_writer.start_chromatogram_list()?;
            }
            state if state > MzMLWriterState::ChromatogramList => {
                // Cannot write chromatogram, currently in state which happens
                // after spectra may be written
                return Err(MzMLWriterError::InvalidActionError(self.mzml_writer.state).into());
            }
            _ => {}
        }
        let pos = self.mzml_writer.stream_position()?;
        self.mzml_writer
            .chromatogram_offset_index
            .insert(chromatogram.id().to_string(), pos);

        let mut outer = bstart!("chromatogram");
        let default_array_size = self
            .mzml_writer
            .start_chromatogram(chromatogram, &mut outer)?;
        self.mzml_writer
            .write_param_list(chromatogram.params().iter())?;

        if let Some(precursor) = chromatogram.precursor() {
            self.mzml_writer.write_precursor(precursor)?;
        }
        self.write_binary_data_arrays(
            &chromatogram.arrays,
            BufferContext::Chromatogram,
            default_array_size,
        )?;

        end_event!(self, outer);
        Ok(())
    }

    pub fn write_summary_chromatograms(&mut self) -> WriterResult {
        if !self.mzml_writer.wrote_summaries {
            self.write_chromatogram(&self.mzml_writer.tic_collector.to_chromatogram())?;
            self.write_chromatogram(&self.mzml_writer.bic_collector.to_chromatogram())?;
            self.mzml_writer.wrote_summaries = true;
        }
        Ok(())
    }

    /// Get a reference to the mzML writer's spectrum count.
    pub fn spectrum_count(&self) -> &u64 {
        &self.mzml_writer.spectrum_count
    }

    /// Set the mzML writer's spectrum count.
    pub fn set_spectrum_count(&mut self, spectrum_count: u64) {
        self.mzml_writer.spectrum_count = spectrum_count;
    }

    /// Get a mutable reference to the mzML writer's spectrum count to modify in-place.
    pub fn spectrum_count_mut(&mut self) -> &mut u64 {
        &mut self.mzml_writer.spectrum_count
    }

    /// Get a reference to the mzML writer's chromatogram count.
    pub fn chromatogram_count(&self) -> &u64 {
        &self.mzml_writer.chromatogram_count
    }

    /// Set the mzML writer's chromatogram count.
    pub fn set_chromatogram_count(&mut self, chromatogram_count: u64) {
        self.mzml_writer.chromatogram_count = chromatogram_count;
    }

    /// Get a mutable reference to the mzML writer's chromatogram count to modify in-place.
    pub fn chromatogram_count_mut(&mut self) -> &mut u64 {
        &mut self.mzml_writer.chromatogram_count
    }

    fn make_mzml_buffer(handle: &mut hdf5::File, chunk_size: usize, filters: &[hdf5::filters::Filter]) -> Result<hdf5::Dataset, MzMLbWriterError> {
        let mzml_buffer = handle
            .new_dataset_builder()
            .chunk(chunk_size)
            .empty_as(&TypeDescriptor::Unsigned(IntSize::U1))
            .shape((1..,))
            .set_filters(filters)
            .create("/mzML")?;
        Ok(mzml_buffer)
    }

    fn store_mzml_buffer(&mut self) -> WriterResult {
        self.mzml_writer.flush()?;
        let mzml_buffer = self.mzml_writer.get_mut()?;
        mzml_buffer.resize_to_fit()?;

        let version_str: FixedAscii<10> = match FixedAscii::from_ascii(b"mzMLb 1.0") {
            Ok(val) => val,
            Err(e) => {
                return Err(MzMLbWriterError::IOError(io::Error::new(
                    io::ErrorKind::InvalidData,
                    e,
                )))
            }
        };

        let version_attr = create_fixed_length_str_attribute(&mut mzml_buffer.dataset, "version", 10)?;
        version_attr.write_scalar(&version_str)?;
        Ok(())
    }

    fn finalize_data_buffers(&mut self) -> WriterResult {
        let result: WriterResult = self
            .buffers
            .iter_mut()
            .try_for_each(|(_, v)| -> WriterResult { v.flush(true) });
        result?;
        Ok(())
    }

    pub fn zlib_compression(level: u8) -> Vec<filters::Filter> {
        vec![
            filters::Filter::Shuffle,
            filters::Filter::Deflate(level),
            filters::Filter::Fletcher32,
        ]
    }

    pub fn blosc_zstd(level: u8) -> Vec<filters::Filter> {
        vec![filters::Filter::Blosc(
            filters::Blosc::ZStd,
            level,
            filters::BloscShuffle::Byte,
        )]
    }

    fn write_spectrum_offset_index(&mut self) -> WriterResult {
        let index_name = "mzML_spectrumIndex";
        let index_id_ref = "mzML_spectrumIndex_idRef";
        let mut id_arrays = Cursor::new(Vec::new());
        let mut offset_arrays: Vec<u64> = Vec::new();
        for (k, v) in self.mzml_writer.spectrum_offset_index.iter() {
            id_arrays.write_all(k.as_bytes())?;
            id_arrays.write_all(b"\x00")?;
            offset_arrays.push(*v);
        }
        offset_arrays.push(0);

        let offset_arrays: Array1<u64> = offset_arrays.into();
        let id_arrays: Array1<u8> = id_arrays.into_inner().into();

        self.handle
            .new_dataset_builder()
            .with_data_as(offset_arrays.view(), &TypeDescriptor::Unsigned(IntSize::U8))
            .set_filters(&self.filters)
            .chunk(self.chunk_size.min(offset_arrays.len()))
            .create(index_name)?;

        self.handle
            .new_dataset_builder()
            .with_data_as(id_arrays.view(), &TypeDescriptor::Unsigned(IntSize::U1))
            .set_filters(&self.filters)
            .chunk(self.chunk_size.min(id_arrays.len()))
            .create(index_id_ref)?;
        Ok(())
    }

    fn write_chromatogram_offset_index(&mut self) -> WriterResult {
        let index_name = "mzML_chromatogramIndex";
        let index_id_ref = "mzML_chromatogramIndex_idRef";
        let mut id_arrays = Cursor::new(Vec::new());
        let mut offset_arrays: Vec<u64> = Vec::new();
        for (k, v) in self.mzml_writer.chromatogram_offset_index.iter() {
            id_arrays.write_all(k.as_bytes())?;
            id_arrays.write_all(b"\x00")?;
            offset_arrays.push(*v);
        }
        offset_arrays.push(0);

        let offset_arrays: Array1<u64> = offset_arrays.into();
        let id_arrays: Array1<u8> = id_arrays.into_inner().into();

        self.handle
            .new_dataset_builder()
            .with_data_as(offset_arrays.view(), &TypeDescriptor::Unsigned(IntSize::U8))
            .set_filters(&self.filters)
            .chunk(self.chunk_size.min(offset_arrays.len()))
            .create(index_name)?;

        self.handle
            .new_dataset_builder()
            .with_data_as(id_arrays.view(), &TypeDescriptor::Unsigned(IntSize::U1))
            .set_filters(&self.filters)
            .chunk(self.chunk_size.min(id_arrays.len()))
            .create(index_id_ref)?;
        Ok(())
    }

    pub fn close(&mut self) -> WriterResult {
        if self.mzml_writer.state == MzMLWriterState::End {
            return Ok(());
        }
        if self.mzml_writer.state < MzMLWriterState::ChromatogramListClosed {
            self.write_summary_chromatograms()?;
        }
        self.mzml_writer.close()?;
        self.store_mzml_buffer()?;
        self.finalize_data_buffers()?;
        self.write_spectrum_offset_index()?;
        self.write_chromatogram_offset_index()?;
        Ok(())
    }
}

impl<
        C: CentroidLike + Default + BuildArrayMapFrom,
        D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom,
    > Drop for MzMLbWriterType<C, D>
{
    fn drop(&mut self) {
        MzMLbWriterType::close(self).unwrap();
    }
}

pub type MzMLbWriter = MzMLbWriterType<CentroidPeak, DeconvolutedPeak>;

fn build_nullterm_fixed_ascii_type_id(size: usize) -> Result<hid_t, hdf5::Error> {
    let string_id = hdf5::h5try!(H5Tcopy(*H5T_C_S1));
    h5try!(H5Tset_cset(string_id, H5T_cset_t::H5T_CSET_ASCII));
    h5try!(H5Tset_strpad(string_id, H5T_str_t::H5T_STR_NULLTERM));
    h5try!(H5Tset_size(string_id, size));
    Ok(string_id)
}

fn create_fixed_length_str_attribute(
    dset: &mut hdf5::Dataset,
    name: &str,
    size: usize,
) -> Result<Attribute, hdf5::Error> {
    let type_id = build_nullterm_fixed_ascii_type_id(size)?;
    let dataspace = Dataspace::try_new(Extents::Scalar)?;
    let name = CString::new(name).expect("Failed to convert name to a C string");
    let attr_id = h5try!(hdf5_sys::h5a::H5Acreate2(
        dset.id(),
        name.as_ptr(),
        type_id,
        dataspace.id(),
        hdf5_sys::h5p::H5P_DEFAULT,
        hdf5_sys::h5p::H5P_DEFAULT
    ));

    let attr = h5lock!(from_id::<Attribute>(attr_id))?;
    Ok(attr)
}

#[cfg(test)]
mod test {
    use tempfile;
    use test_log;

    use crate::prelude::*;

    use super::super::super::MzMLReader;
    use super::super::MzMLbReader;
    use super::*;

    #[test_log::test]
    fn test_writer() -> WriterResult {
        let tmpdir = tempfile::tempdir()?;
        let path = tmpdir.path().join("duplicate.mzMLb");
        let mut reader = MzMLReader::open_path("test/data/small.mzML")?;
        let mut writer = MzMLbWriterBuilder::new(&path)
            .with_zlib_compression(9)
            .with_chunk_size(DEFAULT_CHUNK_SIZE)
            .create()?;
        writer.copy_metadata_from(&reader);

        for s in reader.iter() {
            writer.write(&s)?;
        }
        writer.close()?;

        eprintln!("Path written to: {}", path.display());
        eprintln!("Exists? {}", path.exists());

        let mut reader2 = MzMLbReader::new(&path)?;
        assert_eq!(reader.file_description(), reader2.file_description());

        for (a, b) in reader.iter().zip(reader2.iter()) {
            assert_eq!(a.index(), b.index());
            assert_eq!(a.id(), b.id());
            assert_eq!(a.ms_level(), b.ms_level());
            for (x, y) in a
                .arrays
                .unwrap()
                .mzs()
                .unwrap()
                .iter()
                .zip(b.arrays.unwrap().mzs().unwrap().iter())
            {
                assert!((x - y).abs() < 1e-3)
            }
        }
        Ok(())
    }
}
