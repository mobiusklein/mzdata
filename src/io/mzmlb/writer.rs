#![allow(unused)]
use std::collections::HashMap;
use std::io::{self, prelude::*};
use std::mem;
use std::path::Path;

use filename;
use hdf5::types::{FixedUnicode, IntSize, TypeDescriptor};
use hdf5::{self, filters, Dataset};
use log::{debug, warn};
use mzpeaks::peak_set::PeakSetVec;
use mzpeaks::{CentroidLike, DeconvolutedCentroidLike, Mass, MZ};
use ndarray::{Array1, Ix1};
use thiserror::Error;

use crate::io::prelude::MSDataFileMetadata;
use crate::io::traits::ScanWriter;
use crate::meta::{DataProcessing, FileDescription, InstrumentConfiguration, Software};
use crate::params::{ControlledVocabulary, Param};
use crate::spectrum::MultiLayerSpectrum;
use crate::spectrum::signal::{
    delta_decoding, linear_prediction_decoding, BinaryCompressionType, BinaryDataArrayType,
    BuildArrayMapFrom, DataArray,
};

use crate::io::mzml::{MzMLSpectrumWriter, MzMLWriterError, MzMLWriterState, MzMLWriterType};

#[derive(Debug, Error)]
pub enum MzMLbWriterError {
    #[error("An HDF5-related error occurred: {0}")]
    HDF5Error(#[from] hdf5::Error),
    #[error("An mzML-related error occurred: {0}")]
    MzMLError(#[from] MzMLWriterError),
    #[error("An error occurred while decoding binary data: {0}")]
    IOError(#[from] io::Error),
}

impl From<MzMLbWriterError> for io::Error {
    fn from(value: MzMLbWriterError) -> Self {
        Self::new(io::ErrorKind::Other, value)
    }
}

impl<C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default> MSDataFileMetadata
    for MzMLbWriterType<C, D>
where
    PeakSetVec<C, MZ>: BuildArrayMapFrom,
    PeakSetVec<D, Mass>: BuildArrayMapFrom,
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
}

pub type WriterResult = Result<(), MzMLbWriterError>;


impl<'a, W: Write + Seek, C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default>
    ScanWriter<'a, W, C, D> for MzMLbWriterType<C, D>
where
    PeakSetVec<C, MZ>: BuildArrayMapFrom,
    PeakSetVec<D, Mass>: BuildArrayMapFrom,
{
    fn write(&mut self, spectrum: &'a MultiLayerSpectrum<C, D>) -> io::Result<usize> {
        match self.write_spectrum(spectrum) {
            Ok(()) => {
                let pos = self.mzml_writer.stream_position()?;
                Ok(pos as usize)
            }
            Err(err) => {
                let msg = err.to_string();
                Err(io::Error::new(io::ErrorKind::InvalidData, msg))
            }
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        self.mzml_writer.flush()?;
        Ok(())
    }
}


#[derive(Debug)]
pub struct MzMLbWriterType<
    C: CentroidLike + Default + 'static,
    D: DeconvolutedCentroidLike + Default + 'static,
> where
    PeakSetVec<C, MZ>: BuildArrayMapFrom,
    PeakSetVec<D, Mass>: BuildArrayMapFrom,
{
    handle: hdf5::File,
    mzml_writer: MzMLWriterType<io::Cursor<Vec<u8>>, C, D>,
    chunk_size: usize,
}

impl<C: CentroidLike + Default + 'static, D: DeconvolutedCentroidLike + Default + 'static>
    MzMLbWriterType<C, D>
where
    PeakSetVec<C, MZ>: BuildArrayMapFrom,
    PeakSetVec<D, Mass>: BuildArrayMapFrom,
{
    pub fn new<P: AsRef<Path>>(path: &P) -> io::Result<Self> {
        let handle = match hdf5::File::create(path) {
            Ok(handle) => handle,
            Err(e) => return Err(MzMLbWriterError::HDF5Error(e).into()),
        };

        let buffer = io::Cursor::new(Vec::new());
        let mzml_writer: MzMLWriterType<io::Cursor<Vec<u8>>, C, D> = MzMLWriterType::new(buffer);

        Ok(Self {
            handle,
            mzml_writer,
            chunk_size: 1000000,
        })
    }

    pub fn write_spectrum(&mut self, spectrum: &MultiLayerSpectrum<C, D>) -> WriterResult {

        Ok(())
    }

    /**
    Close the wrapping `<indexedmzML>` document, which will trigger writing
    out the offset indices and file checksum at the tail of the document.
    */
    pub fn close(&mut self) -> WriterResult {
        self.mzml_writer.close()?;
        let mut buffer = self.mzml_writer.get_mut()?;
        let raw_bytes = mem::take(buffer.get_mut());
        let mut mzml_buffer = self
            .handle
            .new_dataset_builder()
            .chunk(self.chunk_size)
            .with_data_as(&raw_bytes, &TypeDescriptor::Unsigned(IntSize::U1))
            .create("/mzML")?;

        mzml_buffer
            .new_attr_builder()
            .with_data_as("mzMLb 1.0", &TypeDescriptor::FixedUnicode(9))
            .create("version")?;

        Ok(())
    }
}
