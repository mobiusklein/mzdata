use std::{
    collections::HashMap,
    convert::TryInto,
    fs,
    io::{self, prelude::*, SeekFrom},
    marker::PhantomData,
    str,
};

use log::warn;
use thiserror::Error;

use mzpeaks::{CentroidPeak, DeconvolutedPeak};

use super::super::{
    offset_index::OffsetIndex,
    traits::{
        ChromatogramSource, MZFileReader, RandomAccessSpectrumIterator, SeekRead,
        SpectrumAccessError, SpectrumSource,
    },
    utils::DetailLevel,
};

use crate::meta::{
    DataProcessing, FileDescription, InstrumentConfiguration, MSDataFileMetadata,
    MassSpectrometryRun, Sample, Software,
};

use crate::params::{ControlledVocabulary, Param, ParamDescribed};

use crate::spectrum::{
    bindata::{
        vec_as_bytes, ArrayType, BinaryArrayMap, BinaryDataArrayType, BuildArrayMapFrom,
        BuildFromArrayMap, DataArray,
    },
    spectrum_types::{
        CentroidPeakAdapting, CentroidSpectrumType, DeconvolutedPeakAdapting, MultiLayerSpectrum,
    },
    Chromatogram, Precursor, PrecursorSelection, SelectedIon, SignalContinuity,
    SpectrumDescription,
};
use crate::utils::neutral_mass;

#[derive(PartialEq, Debug)]
pub enum MGFParserState {
    Start,
    FileHeader,
    ScanHeaders,
    Peaks,
    Between,
    Done,
    Error,
}

#[derive(Debug, Error)]
pub enum MGFError {
    #[error("No error occurred")]
    NoError,
    #[error("Encountered a malformed peak line")]
    MalformedPeakLine,
    #[error("Encountered a malformed header line: {0}")]
    MalformedHeaderLine(String),
    #[error("Too many columns for peak line encountered")]
    NotEnoughColumnsForPeakLine,
    #[error("Encountered an IO error: {0}")]
    IOError(
        #[from]
        #[source]
        io::Error,
    ),
}

#[derive(Debug)]
struct SpectrumBuilder<
    C: CentroidPeakAdapting = CentroidPeak,
    D: DeconvolutedPeakAdapting = DeconvolutedPeak,
> {
    pub description: SpectrumDescription,
    pub mz_array: Vec<f64>,
    pub intensity_array: Vec<f32>,
    pub charge_array: Vec<i32>,
    pub has_charge: u32,
    pub precursor_charge: Option<i32>,
    pub detail_level: DetailLevel,
    empty_metadata: bool,
    centroided_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> SpectrumBuilder<C, D> {
    pub fn is_empty(&self) -> bool {
        self.empty_metadata && self.mz_array.is_empty()
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> Default for SpectrumBuilder<C, D> {
    fn default() -> Self {
        let mut description = SpectrumDescription::default();
        description.signal_continuity = SignalContinuity::Centroid;
        description.ms_level = 2;
        Self {
            description,
            empty_metadata: true,
            mz_array: Default::default(),
            intensity_array: Default::default(),
            charge_array: Default::default(),
            has_charge: Default::default(),
            precursor_charge: Default::default(),
            detail_level: Default::default(),
            centroided_type: Default::default(),
            deconvoluted_type: Default::default(),
        }
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> From<SpectrumBuilder<C, D>>
    for BinaryArrayMap
{
    fn from(val: SpectrumBuilder<C, D>) -> Self {
        let mut arrays = BinaryArrayMap::new();
        arrays.add(DataArray::wrap(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            vec_as_bytes(val.mz_array),
        ));
        arrays.add(DataArray::wrap(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            vec_as_bytes(val.intensity_array),
        ));
        if val.has_charge > 0 {
            arrays.add(DataArray::wrap(
                &ArrayType::ChargeArray,
                BinaryDataArrayType::Int32,
                vec_as_bytes(val.charge_array),
            ));
        }
        arrays
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> SpectrumBuilder<C, D> {
    pub fn into_spectrum(self, spectrum: &mut MultiLayerSpectrum<C, D>) {
        if self.has_charge > 0 {
            spectrum.deconvoluted_peaks = Some(
                self.mz_array
                    .into_iter()
                    .zip(self.intensity_array)
                    .zip(self.charge_array)
                    .map(|((mz, inten), z)| {
                        DeconvolutedPeak {
                            neutral_mass: neutral_mass(mz, if z != 0 { z } else { 1 }),
                            intensity: inten,
                            charge: if z != 0 { z } else { 1 },
                            ..Default::default()
                        }
                        .into()
                    })
                    .collect(),
            )
        } else {
            spectrum.peaks = Some(
                self.mz_array
                    .into_iter()
                    .zip(self.intensity_array)
                    .map(|(mz, inten)| {
                        CentroidPeak {
                            mz,
                            intensity: inten,
                            ..Default::default()
                        }
                        .into()
                    })
                    .collect(),
            )
        }
        spectrum.description = self.description;
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> From<SpectrumBuilder<C, D>>
    for MultiLayerSpectrum<C, D>
{
    fn from(builder: SpectrumBuilder<C, D>) -> MultiLayerSpectrum<C, D> {
        let mut spec = MultiLayerSpectrum::default();
        builder.into_spectrum(&mut spec);
        spec
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> From<SpectrumBuilder<C, D>>
    for CentroidSpectrumType<C>
where
    C: BuildFromArrayMap + BuildArrayMapFrom,
    D: BuildFromArrayMap + BuildArrayMapFrom,
{
    fn from(builder: SpectrumBuilder<C, D>) -> CentroidSpectrumType<C> {
        let spec: MultiLayerSpectrum<C, D> = builder.into();
        spec.try_into().unwrap()
    }
}

/// An MGF (Mascot Generic Format) file parser that supports iteration and random access.
/// The parser produces [`Spectrum`](crate::spectrum::Spectrum) instances. These may be
/// converted directly into [`CentroidSpectrum`](crate::spectrum::CentroidSpectrum)
/// instances which better represent the nature of this preprocessed data type.
pub struct MGFReaderType<
    R: io::Read,
    C: CentroidPeakAdapting = CentroidPeak,
    D: DeconvolutedPeakAdapting = DeconvolutedPeak,
> {
    pub handle: io::BufReader<R>,
    pub state: MGFParserState,
    pub offset: usize,
    pub error: Option<MGFError>,
    index: OffsetIndex,
    file_description: FileDescription,
    instrument_configurations: HashMap<u32, InstrumentConfiguration>,
    softwares: Vec<Software>,
    samples: Vec<Sample>,
    data_processings: Vec<DataProcessing>,
    run: MassSpectrometryRun,
    pub detail_level: DetailLevel,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
}

impl<R: io::Read, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MGFReaderType<R, C, D> {
    fn parse_peak_from_line(
        &mut self,
        line: &str,
        builder: &mut SpectrumBuilder<C, D>,
    ) -> Option<bool> {
        let mut chars = line.chars();
        let first = chars.next().unwrap();
        if first.is_numeric() {
            let mut it = line.split_ascii_whitespace();
            let mz_token = it.next().unwrap();
            let mut intensity_token = "";
            let mut charge_token_opt = None;
            let nparts = if let Some(i) = it.next() {
                intensity_token = i;
                charge_token_opt = it.next();
                if charge_token_opt.is_some() {
                    3
                } else {
                    2
                }
            } else {
                1
            };

            if nparts < 2 {
                self.state = MGFParserState::Error;
                self.error = Some(MGFError::NotEnoughColumnsForPeakLine);
                return None;
            }
            if !matches!(builder.detail_level, DetailLevel::MetadataOnly) {
                let mz: f64 = mz_token.parse().unwrap();
                let intensity: f32 = intensity_token.parse().unwrap();
                builder.mz_array.push(mz);
                builder.intensity_array.push(intensity);

                if nparts == 3 {
                    let charge = charge_token_opt.unwrap().parse().unwrap();
                    builder.charge_array.push(charge);
                    builder.has_charge += 1;
                } else {
                    builder.charge_array.push(0);
                }
            }

            Some(true)
        } else {
            None
        }
    }

    fn handle_scan_header(&mut self, line: &str, builder: &mut SpectrumBuilder<C, D>) -> bool {
        let peak_line = self.parse_peak_from_line(line, builder).unwrap_or(false);
        if peak_line {
            self.state = MGFParserState::Peaks;
            true
        } else if line == "END IONS" {
            self.state = MGFParserState::Between;
            true
        } else if line.contains('=') {
            let (key, value) = line.split_once('=').unwrap();
            let value = value.trim();
            builder.empty_metadata = false;
            match key {
                "TITLE" => builder.description.id = value.to_string(),
                "RTINSECONDS" => {
                    let scan_ev = builder
                        .description
                        .acquisition
                        .first_scan_mut()
                        .expect("Automatically adds scan event");
                    scan_ev.start_time = value.parse::<f64>().unwrap() / 60.0
                }
                "PEPMASS" => {
                    let mut parts = value.split_ascii_whitespace();
                    let mz = match parts.next() {
                        Some(s) => s,
                        None => {
                            self.state = MGFParserState::Error;
                            self.error = Some(MGFError::MalformedHeaderLine(
                                "No m/z value in PEPMASS header".into(),
                            ));
                            return false;
                        }
                    };
                    let mz: f64 = match mz.parse() {
                        Ok(mz) => mz,
                        Err(e) => {
                            self.state = MGFParserState::Error;
                            self.error = Some(MGFError::MalformedHeaderLine(format!(
                                "Malformed m/z value in PEPMASS header {value}: {e}"
                            )));
                            return false;
                        }
                    };

                    let intensity: f32 = parts
                        .next()
                        .map(|v| v.parse())
                        .unwrap_or_else(|| Ok(0.0))
                        .map_err(|e| warn!("Failed to parse PEPMASS intensity {value}: {e}"))
                        .unwrap_or_default();
                    let charge: Option<i32> = match parts.next() {
                        Some(c) => self.parse_charge(c),
                        None => builder.precursor_charge,
                    };
                    builder.description.precursor = Some(Precursor {
                        ions: vec![SelectedIon {
                            mz,
                            intensity,
                            charge,
                            ..Default::default()
                        }],
                        ..Default::default()
                    });
                }
                "CHARGE" => {
                    builder.precursor_charge = self.parse_charge(value);
                    if let Some(ion) = builder
                        .description
                        .precursor
                        .get_or_insert_with(Precursor::default)
                        .iter_mut()
                        .last()
                    {
                        if ion.charge.is_none() {
                            ion.charge = builder.precursor_charge
                        }
                    }
                }
                &_ => {
                    builder
                        .description
                        .add_param(Param::new_key_value(key.to_lowercase(), value));
                }
            };

            true
        } else {
            self.state = MGFParserState::Error;
            self.error = Some(MGFError::MalformedHeaderLine(
                "No '=' in header line".into(),
            ));
            false
        }
    }

    fn parse_charge(&mut self, value: &str) -> Option<i32> {
        let (sign, value, tail_sign) = if let Some(stripped) = value.strip_suffix('+') {
            (1, stripped, true)
        } else if let Some(stripped) = value.strip_suffix('-') {
            (-1, stripped, true)
        } else {
            (1, value, false)
        };

        if tail_sign && (value.starts_with('-') || value.starts_with('+')) {
            self.state = MGFParserState::Error;
            self.error = Some(MGFError::MalformedHeaderLine(format!(
                "Could not parse charge value {value}"
            )));
            return None;
        }

        match value.parse::<i32>() {
            Ok(z) => Some(sign * z),
            Err(e) => {
                self.state = MGFParserState::Error;
                self.error = Some(MGFError::MalformedHeaderLine(format!(
                    "Could not parse charge value {value} : {e}"
                )));
                return None;
            }
        }
    }

    fn handle_peak(&mut self, line: &str, builder: &mut SpectrumBuilder<C, D>) -> bool {
        let peak_line = self.parse_peak_from_line(line, builder).unwrap_or(false);
        if peak_line {
            true
        } else if line == "END IONS" {
            self.state = MGFParserState::Between;
            false
        } else {
            self.state = MGFParserState::Error;
            self.error = Some(MGFError::MalformedPeakLine);
            false
        }
    }

    fn handle_start(&mut self, line: &str) -> bool {
        if line.contains('=') {
            match self.state {
                MGFParserState::Start => {
                    self.state = MGFParserState::FileHeader;
                    true
                }
                MGFParserState::FileHeader => true,
                _ => false,
            }
        } else if line == "BEGIN IONS" {
            self.state = MGFParserState::ScanHeaders;
            true
        } else {
            false
        }
    }

    fn handle_between(&mut self, line: &str) -> bool {
        if line == "BEGIN IONS" {
            self.state = MGFParserState::ScanHeaders;
        }
        true
    }

    fn read_line(&mut self, buffer: &mut String) -> io::Result<usize> {
        self.handle.read_line(buffer)
    }

    /// Read the next spectrum from the file, if there is one.
    pub fn read_next(&mut self) -> Option<MultiLayerSpectrum<C, D>> {
        let mut builder = SpectrumBuilder::<C, D>::default();
        self._parse_into(&mut builder)
            .ok()
            .and_then(|(_, started_spectrum)| {
                (started_spectrum && !builder.is_empty()).then(|| builder.into())
            })
    }

    /// Read the next spectrum's contents directly into the passed [`SpectrumBuilder`].
    fn _parse_into(
        &mut self,
        builder: &mut SpectrumBuilder<C, D>,
    ) -> Result<(usize, bool), MGFError> {
        let mut buffer = String::new();
        let mut work = true;
        let mut offset: usize = 0;
        let mut had_begin_ions = false;

        while work {
            buffer.clear();
            let b = match self.read_line(&mut buffer) {
                Ok(b) => b,
                Err(err) => {
                    self.state = MGFParserState::Error;
                    return Err(MGFError::IOError(err));
                }
            };

            // Count how many bytes we've read from the source
            offset += b;
            if b == 0 {
                self.state = MGFParserState::Done;
                break;
            }

            let line = buffer.trim();

            // Skip empty lines
            if line.is_empty() {
                continue;
            }

            work = match self.state {
                MGFParserState::Start | MGFParserState::FileHeader => self.handle_start(line),
                MGFParserState::ScanHeaders => {
                    had_begin_ions = true;
                    self.handle_scan_header(line, builder)
                }
                MGFParserState::Peaks => self.handle_peak(line, builder),
                MGFParserState::Between => self.handle_between(line),
                MGFParserState::Done => false,
                MGFParserState::Error => {
                    return Err(self.error.take().unwrap_or(MGFError::NoError));
                }
            };

            if self.state == MGFParserState::Error {
                return Err(self.error.take().unwrap_or(MGFError::NoError));
            }
        }
        Ok((offset, had_begin_ions))
    }

    pub fn read_into(
        &mut self,
        spectrum: &mut MultiLayerSpectrum<C, D>,
    ) -> Result<usize, MGFError> {
        let mut accumulator = SpectrumBuilder::default();
        match self._parse_into(&mut accumulator) {
            Ok((sz, started_spectrum)) => {
                if !started_spectrum {
                    Err(MGFError::IOError(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "EOF found before spectrum started",
                    )))
                } else {
                    accumulator.into_spectrum(spectrum);
                    Ok(sz)
                }
            }
            Err(err) => Err(err),
        }
    }

    fn default_file_description() -> FileDescription {
        let mut fd = FileDescription::default();
        let mut term = Param::new();
        term.name = "MSn spectrum".to_owned();
        term.accession = Some(1000580);
        term.controlled_vocabulary = Some(ControlledVocabulary::MS);
        fd.add_param(term);
        fd
    }

    /// Create a new, unindexed MGF parser
    pub fn new(file: R) -> MGFReaderType<R, C, D> {
        let handle = io::BufReader::with_capacity(500, file);
        MGFReaderType {
            handle,
            state: MGFParserState::Start,
            offset: 0,
            error: None,
            index: OffsetIndex::new("spectrum".to_owned()),
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
            instrument_configurations: HashMap::new(),
            data_processings: Vec::new(),
            softwares: Vec::new(),
            samples: Vec::new(),
            file_description: Self::default_file_description(),
            detail_level: DetailLevel::Full,
            run: MassSpectrometryRun::default(),
        }
    }
}

impl<R: io::Read, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> Iterator
    for MGFReaderType<R, C, D>
{
    type Item = MultiLayerSpectrum<C, D>;

    /// Read the next spectrum from the file.
    fn next(&mut self) -> Option<Self::Item> {
        self.read_next()
    }
}

impl<R: SeekRead, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MGFReaderType<R, C, D> {
    /// Construct a new MGFReaderType and build an offset index
    /// using [`Self::build_index`]
    pub fn new_indexed(file: R) -> MGFReaderType<R, C, D> {
        let mut reader = Self::new(file);
        reader.build_index();
        reader
    }

    pub fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.handle.seek(pos)
    }

    /// Builds an offset index to each `BEGIN IONS` line
    /// by doing a fast pre-scan of the text file.
    pub fn build_index(&mut self) -> u64 {
        let mut offset: u64 = 0;
        let mut last_start: u64 = 0;

        let mut found_start = false;

        let start = self
            .handle
            .stream_position()
            .expect("Failed to save restore location");
        self.seek(SeekFrom::Start(0))
            .expect("Failed to reset stream to beginning");

        let mut buffer: Vec<u8> = Vec::new();

        loop {
            buffer.clear();
            let b = match self.handle.read_until(b'\n', &mut buffer) {
                Ok(b) => b,
                Err(err) => {
                    panic!("Error while reading file: {}", err);
                }
            };
            if b == 0 {
                break;
            }
            if buffer.starts_with(b"BEGIN IONS") {
                found_start = true;
                last_start = offset;
            } else if found_start && buffer.starts_with(b"TITLE=") {
                match str::from_utf8(&buffer[6..]) {
                    Ok(string) => {
                        self.index.insert(string.to_owned(), last_start);
                    }
                    Err(_err) => {}
                };
                found_start = false;
                last_start = 0;
            }
            offset += b as u64;
        }
        self.seek(SeekFrom::Start(start))
            .expect("Failed to restore location");
        self.index.init = true;
        if self.index.is_empty() {
            warn!("An index was built but no entries were found")
        }
        offset
    }
}

impl<R: SeekRead, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    SpectrumSource<C, D, MultiLayerSpectrum<C, D>> for MGFReaderType<R, C, D>
{
    fn detail_level(&self) -> &DetailLevel {
        &self.detail_level
    }

    fn set_detail_level(&mut self, detail_level: DetailLevel) {
        self.detail_level = detail_level;
    }

    /// Retrieve a spectrum by it's native ID
    fn get_spectrum_by_id(&mut self, id: &str) -> Option<MultiLayerSpectrum<C, D>> {
        let offset = self.index.get(id)?;
        let index = self.index.index_of(id)?;
        let start = self
            .handle
            .stream_position()
            .expect("Failed to save checkpoint");
        self.seek(SeekFrom::Start(offset))
            .expect("Failed to move seek to offset");
        let result = self.read_next();
        self.seek(SeekFrom::Start(start))
            .expect("Failed to restore offset");
        result.map(|mut scan| {
            scan.description.index = index;
            scan
        })
    }

    /// Retrieve a spectrum by it's integer index
    fn get_spectrum_by_index(&mut self, index: usize) -> Option<MultiLayerSpectrum<C, D>> {
        let (_id, byte_offset) = self.index.get_index(index)?;
        let start = self
            .handle
            .stream_position()
            .expect("Failed to save checkpoint");
        self.seek(SeekFrom::Start(byte_offset)).ok()?;
        let result = self.read_next();
        self.seek(SeekFrom::Start(start))
            .expect("Failed to restore offset");
        result.map(|mut scan| {
            scan.description.index = index;
            scan
        })
    }

    /// Return the data stream to the beginning
    fn reset(&mut self) {
        self.seek(SeekFrom::Start(0))
            .expect("Failed to reset file stream");
    }

    fn get_index(&self) -> &OffsetIndex {
        if !self.index.init {
            warn!("Attempting to use an uninitialized offset index on MGFReaderType")
        }
        &self.index
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.index = index;
    }
}

impl<R: SeekRead, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    RandomAccessSpectrumIterator<C, D, MultiLayerSpectrum<C, D>> for MGFReaderType<R, C, D>
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, SpectrumAccessError> {
        match self._offset_of_id(id) {
            Some(offset) => match self.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(SpectrumAccessError::IOError(Some(err))),
            },
            None => Err(SpectrumAccessError::SpectrumIdNotFound(id.to_string())),
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, SpectrumAccessError> {
        match self._offset_of_index(index) {
            Some(offset) => match self.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(SpectrumAccessError::IOError(Some(err))),
            },
            None => Err(SpectrumAccessError::SpectrumIndexNotFound(index)),
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, SpectrumAccessError> {
        match self._offset_of_time(time) {
            Some(offset) => match self.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(SpectrumAccessError::IOError(Some(err))),
            },
            None => Err(SpectrumAccessError::SpectrumNotFound),
        }
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    MZFileReader<C, D, MultiLayerSpectrum<C, D>> for MGFReaderType<fs::File, C, D>
{
    fn open_file(source: fs::File) -> io::Result<Self> {
        Ok(Self::new(source))
    }

    fn construct_index_from_stream(&mut self) -> u64 {
        self.build_index()
    }
}

/// The MGF format does not contain any consistent metadata, but additional
/// information can be included after creation.
impl<R: io::Read, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MSDataFileMetadata
    for MGFReaderType<R, C, D>
{
    crate::impl_metadata_trait!();

    fn spectrum_count_hint(&self) -> Option<u64> {
        if self.index.init {
            Some(self.index.len() as u64)
        } else {
            None
        }
    }

    fn run_description(&self) -> Option<&MassSpectrometryRun> {
        Some(&self.run)
    }

    fn run_description_mut(&mut self) -> Option<&mut MassSpectrometryRun> {
        Some(&mut self.run)
    }
}

impl<R: Read, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> ChromatogramSource
    for MGFReaderType<R, C, D>
{
    fn get_chromatogram_by_id(&mut self, _: &str) -> Option<Chromatogram> {
        None
    }

    fn get_chromatogram_by_index(&mut self, _: usize) -> Option<Chromatogram> {
        None
    }
}

pub type MGFReader<R> = MGFReaderType<R, CentroidPeak, DeconvolutedPeak>;

pub fn is_mgf(buf: &[u8]) -> bool {
    let needle = b"BEGIN IONS";
    buf.windows(needle.len()).any(|window| window == needle)
}
