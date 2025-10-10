use std::{
    collections::HashMap,
    convert::TryInto,
    fs,
    io::{self, BufRead, BufReader, Read, Seek},
    marker::PhantomData,
    num::ParseFloatError,
};

use indexmap::IndexMap;
use mzpeaks::{CentroidLike, DeconvolutedCentroidLike, DeconvolutedPeak};

use crate::{
    io::{utils::DetailLevel, MZFileReader, RandomAccessSpectrumIterator, SpectrumSource},
    meta::{
        DataProcessing, FileDescription, InstrumentConfiguration, MSDataFileMetadata, Sample,
        Software,
    },
    spectrum::{
        bindata::{
            ArrayType, BinaryArrayMap, BinaryDataArrayType, BuildArrayMapFrom, BuildFromArrayMap,
            DataArray,
        },
        scan_properties::*,
        spectrum_types::{CentroidSpectrumType, MultiLayerSpectrum, RawSpectrum},
    },
};

pub type Bytes = Vec<u8>;

#[doc(hidden)]
impl<R: std::io::Read, C: CentroidLike, D> From<XyReaderType<R, C, D>> for CentroidSpectrumType<C>
where
    C: BuildFromArrayMap + BuildArrayMapFrom,
{
    fn from(mut val: XyReaderType<R, C, D>) -> Self {
        let mut spec = MultiLayerSpectrum::<C, DeconvolutedPeak>::default();
        if let Some(raw) = val.read_next() {
            spec = raw.into();
        }
        spec.try_into().unwrap()
    }
}

#[doc(hidden)]
impl<
        R: std::io::Read,
        C: CentroidLike + BuildFromArrayMap,
        D: DeconvolutedCentroidLike + BuildFromArrayMap,
    > From<XyReaderType<R, C, D>> for MultiLayerSpectrum<C, D>
{
    fn from(mut val: XyReaderType<R, C, D>) -> Self {
        let mut spec = MultiLayerSpectrum::<C, D>::default();
        if let Some(raw) = val.read_next() {
            spec = raw.into();
        }
        spec
    }
}

#[doc(hidden)]
impl<R: std::io::Read, C, D> From<XyReaderType<R, C, D>> for RawSpectrum {
    fn from(mut val: XyReaderType<R, C, D>) -> Self {
        val.read_next().unwrap_or_default()
    }
}

#[derive(Debug)]
pub enum XyParserError {
    EOF,
    InvalidNumber(String, ParseFloatError),
    MissingColumns(String),
    IOError(std::io::Error),
}

#[derive(Debug, Default)]
pub enum XyParserState {
    #[default]
    Initial,
    Error(XyParserError),
    EOF,
}

/**
A parser that reads .xy files. These files solely contain a single profile mode MS spectrum. No
metadata whatshoever is saved. These files can only contain one spectrum. So when iterating over
the file either one spectrum or none are returned. None if there was an error during reading. When
using the idexing functions [`SpectrumSource::get_spectrum_by_index`] and
[`SpectrumSource::get_spectrum_by_id`] with every single call the file is parsed again and the same
spectrum is returned.

The format assumes `<mz> <intensity>` separated by a space or a tab. The parse is set up to handle
any amount of surrounding whitespace for robustness.
*/
pub struct XyReaderType<R: Read, C, D> {
    /// The raw reader
    handle: BufReader<R>,
    state: XyParserState,
    index: crate::io::OffsetIndex,
    data: PhantomData<(C, D)>,
    file_description: FileDescription,
    instrument_configurations: HashMap<u32, InstrumentConfiguration>,
    softwares: Vec<Software>,
    samples: Vec<Sample>,
    data_processings: Vec<DataProcessing>,
}

const BUFFER_SIZE: usize = 8192;

impl<R: Read, C, D> XyReaderType<R, C, D> {
    /// Create a new [`XyReaderType`] instance, wrapping the [`io::Read`] handle
    /// provided with an [`io::BufReader`].
    pub fn new(file: R) -> XyReaderType<R, C, D> {
        Self::with_buffer_capacity(file, BUFFER_SIZE)
    }

    /// Get the error if parsing failed
    pub fn error(&self) -> Option<&XyParserError> {
        match &self.state {
            XyParserState::Error(e) => Some(e),
            _ => None,
        }
    }

    pub fn with_buffer_capacity(file: R, capacity: usize) -> XyReaderType<R, C, D> {
        let handle = BufReader::with_capacity(capacity, file);
        XyReaderType {
            handle,
            state: XyParserState::Initial,
            index: crate::io::OffsetIndex {
                name: String::new(),
                offsets: IndexMap::from([("".to_string().into_boxed_str(), 0)]),
                init: true,
            },
            data: PhantomData,
            file_description: FileDescription::default(),
            instrument_configurations: Default::default(),
            softwares: Default::default(),
            samples: Default::default(),
            data_processings: Default::default(),
        }
    }

    /// Read the arrays.
    pub fn read_arrays(&mut self) -> Option<BinaryArrayMap> {
        if matches!(self.state, XyParserState::Initial) {
            match self.read_arrays_inner() {
                Ok(arrays) => Some(arrays),
                Err(e) => {
                    self.state = XyParserState::Error(e);
                    None
                }
            }
        } else {
            None
        }
    }

    /// Read the arrays.
    fn read_arrays_inner(&mut self) -> Result<BinaryArrayMap, XyParserError> {
        let mut intensity_array =
            DataArray::from_name_and_type(&ArrayType::IntensityArray, BinaryDataArrayType::Float32);
        let mut mz_array =
            DataArray::from_name_and_type(&ArrayType::MZArray, BinaryDataArrayType::Float64);
        let mut line = String::new();
        loop {
            line.clear();
            let z = self
                .handle
                .read_line(&mut line)
                .map_err(XyParserError::IOError)?;
            if z == 0 {
                break;
            }
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue; // Ignore empty lines
            }
            match trimmed.split_once(' ').or(trimmed.split_once('\t')) {
                None => return Err(XyParserError::MissingColumns(line)),
                Some((mz, intensity)) => {
                    let mz = mz
                        .trim()
                        .parse::<f64>()
                        .map_err(|e| XyParserError::InvalidNumber(line.clone(), e))?;
                    let intensity = intensity
                        .trim()
                        .parse::<f32>()
                        .map_err(|e| XyParserError::InvalidNumber(line.clone(), e))?;
                    mz_array.push(mz).unwrap();
                    intensity_array.push(intensity).unwrap();
                }
            }
        }

        let mut arrays = BinaryArrayMap::new();
        arrays.add(mz_array);
        arrays.add(intensity_array);
        Ok(arrays)
    }

    /// Read the next spectrum directly. Used to implement iteration.
    pub fn read_next(&mut self) -> Option<RawSpectrum> {
        let arrays = self.read_arrays()?;
        let description = SpectrumDescription {
            signal_continuity: SignalContinuity::Profile,
            ..Default::default()
        };

        Some(RawSpectrum::new(description, arrays))
    }
}

/// [`XyReaderType`] instances are [`Iterator`]s over [`Spectrum`]
impl<
        R: io::Read,
        C: CentroidLike + BuildFromArrayMap,
        D: DeconvolutedCentroidLike + BuildFromArrayMap,
    > Iterator for XyReaderType<R, C, D>
{
    type Item = MultiLayerSpectrum<C, D>;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_next().map(|r| r.into())
    }
}

impl<
        R: crate::prelude::SeekRead,
        C: CentroidLike + BuildFromArrayMap,
        D: DeconvolutedCentroidLike + BuildFromArrayMap,
    > SpectrumSource<C, D, MultiLayerSpectrum<C, D>> for XyReaderType<R, C, D>
{
    fn reset(&mut self) {
        self.state = XyParserState::Initial;
        let _ = self.handle.rewind();
    }

    fn detail_level(&self) -> &DetailLevel {
        &DetailLevel::Full
    }

    fn set_detail_level(&mut self, _detail_level: DetailLevel) {}

    fn get_spectrum_by_id(&mut self, _id: &str) -> Option<MultiLayerSpectrum<C, D>> {
        self.reset();
        self.next()
    }

    fn get_spectrum_by_index(&mut self, _index: usize) -> Option<MultiLayerSpectrum<C, D>> {
        self.reset();
        self.next()
    }

    fn get_index(&self) -> &crate::io::OffsetIndex {
        &self.index
    }

    fn set_index(&mut self, index: crate::io::OffsetIndex) {
        self.index = index;
    }
}

impl<
        C: CentroidLike + BuildFromArrayMap,
        D: DeconvolutedCentroidLike + BuildFromArrayMap,
        R: io::Read + io::Seek,
    > MSDataFileMetadata for XyReaderType<R, C, D>
{
    fn data_processings(&self) -> &Vec<crate::meta::DataProcessing> {
        &self.data_processings
    }

    fn instrument_configurations(
        &self,
    ) -> &std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
        &self.instrument_configurations
    }

    fn file_description(&self) -> &crate::meta::FileDescription {
        &self.file_description
    }

    fn softwares(&self) -> &Vec<crate::meta::Software> {
        &self.softwares
    }

    fn samples(&self) -> &Vec<crate::meta::Sample> {
        &self.samples
    }

    fn data_processings_mut(&mut self) -> &mut Vec<crate::meta::DataProcessing> {
        &mut self.data_processings
    }

    fn instrument_configurations_mut(
        &mut self,
    ) -> &mut std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
        &mut self.instrument_configurations
    }

    fn file_description_mut(&mut self) -> &mut crate::meta::FileDescription {
        &mut self.file_description
    }

    fn softwares_mut(&mut self) -> &mut Vec<crate::meta::Software> {
        &mut self.softwares
    }

    fn samples_mut(&mut self) -> &mut Vec<crate::meta::Sample> {
        &mut self.samples
    }
}

impl<
        C: CentroidLike + BuildFromArrayMap,
        D: DeconvolutedCentroidLike + BuildFromArrayMap,
        R: io::Read + io::Seek,
    > RandomAccessSpectrumIterator<C, D, MultiLayerSpectrum<C, D>> for XyReaderType<R, C, D>
{
    fn start_from_id(
        &mut self,
        _id: &str,
    ) -> Result<&mut Self, crate::prelude::SpectrumAccessError> {
        Ok(self)
    }

    fn start_from_index(
        &mut self,
        _index: usize,
    ) -> Result<&mut Self, crate::prelude::SpectrumAccessError> {
        Ok(self)
    }

    fn start_from_time(
        &mut self,
        _time: f64,
    ) -> Result<&mut Self, crate::prelude::SpectrumAccessError> {
        Ok(self)
    }
}

impl<C: CentroidLike + BuildFromArrayMap, D: DeconvolutedCentroidLike + BuildFromArrayMap>
    MZFileReader<C, D, MultiLayerSpectrum<C, D>> for XyReaderType<std::fs::File, C, D>
{
    fn open_file(source: fs::File) -> io::Result<Self> {
        Ok(Self::new(source))
    }

    fn construct_index_from_stream(&mut self) -> u64 {
        0
    }
}

#[cfg(test)]
mod test {

    use mzpeaks::CentroidPeak;

    use super::*;

    #[test]
    fn simple_xy() {
        let data = "60.406 140051.00
61.680 140877.00
 63.589\t 141602.00  
65.496 142758.00
67.403 138627.00
69.309 139428.00\t\t
\t\t71.217 140428.00
73.129 138416.00
  75.031 137233.00
76.937 141583.00
78.845\t \t134504.00
80.753 141383.00
82.661\t144141.00
84.568 141220.00
86.474 136657.00
88.383   142199.00
90.290 143039.00
92.197 142881.00
94.105 143715.00

96.012 139754.00
97.918 141515.00";
        let mut reader: XyReaderType<&[u8], CentroidPeak, DeconvolutedPeak> =
            XyReaderType::new(data.as_bytes());
        let spectrum = reader.next();
        println!("Error: {:?}", reader.error());
        let spectrum = spectrum.unwrap();
        assert_eq!(spectrum.arrays.as_ref().unwrap().mzs().unwrap().len(), 21);
        assert_eq!(
            spectrum
                .arrays
                .as_ref()
                .unwrap()
                .intensities()
                .unwrap()
                .len(),
            21
        );
    }
}
