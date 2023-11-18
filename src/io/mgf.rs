use std::collections::HashMap;
use std::convert::TryInto;
use std::io::BufWriter;
use std::mem;
use std::fs;
use std::io;
use std::io::prelude::*;
use std::io::SeekFrom;
use std::marker::PhantomData;
use std::str;

use log::warn;
use thiserror::Error;

use lazy_static::lazy_static;
use mzpeaks::MZPeakSetType;
use mzpeaks::MassPeakSetType;
use regex::Regex;

use mzpeaks::peak::KnownCharge;
use mzpeaks::{CentroidPeak, DeconvolutedPeak, IntensityMeasurement, MZLocated};

use crate::meta::{
    DataProcessing, FileDescription, InstrumentConfiguration, MSDataFileMetadata, Software,
};
use crate::params::{ControlledVocabulary, Param, ParamDescribed, ParamLike};
use crate::spectrum::signal::{
    vec_as_bytes, ArrayType, BinaryArrayMap, BinaryDataArrayType, BuildArrayMapFrom,
    BuildFromArrayMap, DataArray,
};
use crate::spectrum::spectrum::{
    CentroidPeakAdapting, CentroidSpectrumType, DeconvolutedPeakAdapting, MultiLayerSpectrum,
};
use crate::spectrum::{
    Precursor, PrecursorSelection, SelectedIon, SpectrumBehavior, SpectrumDescription,
};
use crate::utils::neutral_mass;

use super::offset_index::OffsetIndex;
use super::traits::{
    MZFileReader, RandomAccessSpectrumIterator, ScanAccessError, ScanSource, ScanWriter, SeekRead,
};

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
    #[error("Encountered a malformed header line")]
    MalformedHeaderLine,
    #[error("Too many columns for peak line encountered")]
    TooManyColumnsForPeakLine,
    #[error("Encountered an IO error: {0}")]
    IOError(#[from] #[source] io::Error),
}

#[derive(Debug, Default)]
struct SpectrumBuilderFlex<
    C: CentroidPeakAdapting = CentroidPeak,
    D: DeconvolutedPeakAdapting = DeconvolutedPeak,
> {
    pub description: SpectrumDescription,
    pub mz_array: Vec<f64>,
    pub intensity_array: Vec<f32>,
    pub charge_array: Vec<i32>,
    pub has_charge: u32,
    centroided_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> From<SpectrumBuilderFlex<C, D>>
    for BinaryArrayMap
{
    fn from(val: SpectrumBuilderFlex<C, D>) -> Self {
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

impl<
        C: CentroidPeakAdapting + From<CentroidPeak>,
        D: DeconvolutedPeakAdapting + From<DeconvolutedPeak>,
    > SpectrumBuilderFlex<C, D>
{
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

impl<
        C: CentroidPeakAdapting + From<CentroidPeak>,
        D: DeconvolutedPeakAdapting + From<DeconvolutedPeak>,
    > From<SpectrumBuilderFlex<C, D>> for MultiLayerSpectrum<C, D>
{
    fn from(builder: SpectrumBuilderFlex<C, D>) -> MultiLayerSpectrum<C, D> {
        let mut spec = MultiLayerSpectrum::default();
        builder.into_spectrum(&mut spec);
        spec
    }
}

impl<
        C: CentroidPeakAdapting + From<CentroidPeak>,
        D: DeconvolutedPeakAdapting + From<DeconvolutedPeak>,
    > From<SpectrumBuilderFlex<C, D>> for CentroidSpectrumType<C>
where
    MZPeakSetType<C>: BuildFromArrayMap + BuildArrayMapFrom,
    MassPeakSetType<D>: BuildFromArrayMap + BuildArrayMapFrom,
{
    fn from(builder: SpectrumBuilderFlex<C, D>) -> CentroidSpectrumType<C> {
        let spec: MultiLayerSpectrum<C, D> = builder.into();
        spec.try_into().unwrap()
    }
}

/// An MGF (Mascot Generic Format) file parser that supports iteration and random access.
/// The parser produces [`crate::spectrum::Spectrum`] instances. These may be converted directly into [`crate::spectrum::CentroidSpectrum`]
/// instances using [`spectrum::Spectrum.into_centroid`] or the [`From`] trait which better represent
/// the nature of this preprocessed data type.
pub struct MGFReaderType<
    R: io::Read,
    C: CentroidPeakAdapting = CentroidPeak,
    D: DeconvolutedPeakAdapting = DeconvolutedPeak,
> {
    pub handle: io::BufReader<R>,
    pub state: MGFParserState,
    pub offset: usize,
    pub error: Option<MGFError>,
    pub index: OffsetIndex,
    file_description: FileDescription,
    instrument_configurations: HashMap<u32, InstrumentConfiguration>,
    softwares: Vec<Software>,
    data_processings: Vec<DataProcessing>,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
}

// A lazily created static regular expression to parse peak separators
lazy_static! {
    static ref PEAK_SEPERATOR: Regex = Regex::new(r"\t|\s+").unwrap();
}

impl<
        R: io::Read,
        C: CentroidPeakAdapting + From<CentroidPeak>,
        D: DeconvolutedPeakAdapting + From<DeconvolutedPeak>,
    > MGFReaderType<R, C, D>
{
    fn parse_peak_from_line_flex(
        &mut self,
        line: &str,
        builder: &mut SpectrumBuilderFlex<C, D>,
    ) -> Option<bool> {
        let mut chars = line.chars();
        let first = chars.next().unwrap();
        if first.is_numeric() {
            let parts: Vec<&str> = PEAK_SEPERATOR.split(line).collect();
            let nparts = parts.len();
            if !(2..=3).contains(&nparts) {
                self.state = MGFParserState::Error;
                self.error = Some(MGFError::TooManyColumnsForPeakLine);
                return None;
            }
            let mz: f64 = parts[0].parse().unwrap();
            let intensity: f32 = parts[1].parse().unwrap();
            builder.mz_array.push(mz);
            builder.intensity_array.push(intensity);

            if nparts == 3 {
                let charge = parts[2].parse().unwrap();
                builder.charge_array.push(charge);
                builder.has_charge += 1;
            } else {
                builder.charge_array.push(0);
            }

            Some(true)
        } else {
            None
        }
    }

    fn handle_scan_header_flex(
        &mut self,
        line: &str,
        builder: &mut SpectrumBuilderFlex<C, D>,
    ) -> bool {
        let peak_line = self.parse_peak_from_line_flex(line, builder).unwrap_or(false);
        if peak_line {
            self.state = MGFParserState::Peaks;
            true
        } else if line == "END IONS" {
            self.state = MGFParserState::Between;
            true
        } else if line.contains('=') {
            let (key, value) = line.split_once('=').unwrap();
            
            
            match key {
                "TITLE" => builder.description.id = String::from(value),
                "RTINSECONDS" => {
                    let scan_ev = builder
                        .description
                        .acquisition
                        .first_scan_mut()
                        .expect("Automatically adds scan event");
                    scan_ev.start_time = value.parse().unwrap()
                }
                "PEPMASS" => {
                    let mut parts = value.split_ascii_whitespace();
                    let mz: f64 = parts.next().unwrap().parse().unwrap();
                    let intensity: f32 = parts.next().unwrap().parse().unwrap();
                    let charge: Option<i32> = parts.next().map(|c| c.parse().unwrap());
                    builder.description.precursor = Some(Precursor {
                        ion: SelectedIon {
                            mz,
                            intensity,
                            charge,
                            ..Default::default()
                        },
                        ..Default::default()
                    });
                }
                &_ => {
                    builder.description.add_param(Param::new_key_value(
                        key.to_lowercase(),
                        String::from(value),
                    ));
                }
            };

            true
        } else {
            self.state = MGFParserState::Error;
            self.error = Some(MGFError::MalformedHeaderLine);
            false
        }
    }

    fn handle_peak_flex(&mut self, line: &str, builder: &mut SpectrumBuilderFlex<C, D>) -> bool {
        let peak_line = self.parse_peak_from_line_flex(line, builder).unwrap_or(false);
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
        } else if line == "BEGIN IONS" {
            self.state = MGFParserState::ScanHeaders;
        }
        true
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
        let mut builder = SpectrumBuilderFlex::<C, D>::default();
        match self._parse_into_flex(&mut builder) {
            Ok(offset) => {
                if offset > 0 {
                    Some(builder.into())
                } else {
                    None
                }
            }
            Err(err) => {
                println!("An error was encountered: {:?}", err);
                None
            }
        }
    }

    /// Read the next spectrum's contents directly into the passed struct.
    fn _parse_into_flex(
        &mut self,
        builder: &mut SpectrumBuilderFlex<C, D>,
    ) -> Result<usize, MGFError> {
        let mut buffer = String::new();
        let mut work = true;
        let mut offset: usize = 0;
        while work {
            buffer.clear();
            let b = match self.read_line(&mut buffer) {
                Ok(b) => {
                    if b == 0 {
                        work = false;
                    }
                    b
                }
                Err(err) => {
                    self.state = MGFParserState::Error;
                    return Err(MGFError::IOError(err));
                }
            };
            offset += b;
            if b == 0 {
                self.state = MGFParserState::Done;
                break;
            }
            let line = buffer.trim();
            let n = line.len();
            if n == 0 {
                continue;
            }
            if self.state == MGFParserState::Start {
                work = self.handle_start(line);
            } else if self.state == MGFParserState::Between {
                work = self.handle_between(line);
            } else if self.state == MGFParserState::ScanHeaders {
                work = self.handle_scan_header_flex(line, builder)
            } else if self.state == MGFParserState::Peaks {
                work = self.handle_peak_flex(line, builder);
            }
            if matches!(self.state, MGFParserState::Error) {
                let mut err = None;
                mem::swap(&mut self.error, &mut err);
                self.error = None;
                return Err(err.unwrap())
            }
        }
        Ok(offset)
    }

    pub fn read_into(
        &mut self,
        spectrum: &mut MultiLayerSpectrum<C, D>,
    ) -> Result<usize, MGFError> {
        let mut accumulator = SpectrumBuilderFlex::default();
        match self._parse_into_flex(&mut accumulator) {
            Ok(sz) => {
                accumulator.into_spectrum(spectrum);
                Ok(sz)
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
            file_description: Self::default_file_description(),
        }
    }
}

impl<
        R: io::Read,
        C: CentroidPeakAdapting + From<CentroidPeak>,
        D: DeconvolutedPeakAdapting + From<DeconvolutedPeak>,
    > Iterator for MGFReaderType<R, C, D>
{
    type Item = MultiLayerSpectrum<C, D>;

    /// Read the next spectrum from the file.
    fn next(&mut self) -> Option<Self::Item> {
        self.read_next()
    }
}

impl<
        R: SeekRead,
        C: CentroidPeakAdapting + From<CentroidPeak>,
        D: DeconvolutedPeakAdapting + From<DeconvolutedPeak>,
    > MGFReaderType<R, C, D>
{
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

impl<
        R: SeekRead,
        C: CentroidPeakAdapting + From<CentroidPeak>,
        D: DeconvolutedPeakAdapting + From<DeconvolutedPeak>,
    > ScanSource<C, D, MultiLayerSpectrum<C, D>> for MGFReaderType<R, C, D>
{
    /// Retrieve a spectrum by it's native ID
    fn get_spectrum_by_id(&mut self, id: &str) -> Option<MultiLayerSpectrum<C, D>> {
        let offset_ref = self.index.get(id);
        let offset = offset_ref.expect("Failed to retrieve offset");
        let index = self.index.index_of(id).unwrap();
        let start = self
            .handle
            .stream_position()
            .expect("Failed to save checkpoint");
        self.seek(SeekFrom::Start(offset))
            .expect("Failed to move seek to offset");
        let result = self.read_next();
        self.seek(SeekFrom::Start(start))
            .expect("Failed to restore offset");
        match result {
            Some(mut scan) => {
                scan.description.index = index;
                Some(scan)
            }
            None => None,
        }
    }

    /// Retrieve a spectrum by it's integer index
    fn get_spectrum_by_index(&mut self, index: usize) -> Option<MultiLayerSpectrum<C, D>> {
        let (_id, offset) = self.index.get_index(index)?;
        let byte_offset = offset;
        let start = self
            .handle
            .stream_position()
            .expect("Failed to save checkpoint");
        self.seek(SeekFrom::Start(byte_offset)).ok()?;
        let result = self.read_next();
        self.seek(SeekFrom::Start(start))
            .expect("Failed to restore offset");
        match result {
            Some(mut scan) => {
                scan.description.index = index;
                Some(scan)
            }
            None => None,
        }
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

impl<
        R: SeekRead,
        C: CentroidPeakAdapting + From<CentroidPeak>,
        D: DeconvolutedPeakAdapting + From<DeconvolutedPeak>,
    > RandomAccessSpectrumIterator<C, D, MultiLayerSpectrum<C, D>> for MGFReaderType<R, C, D>
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, ScanAccessError> {
        match self._offset_of_id(id) {
            Some(offset) => match self.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(ScanAccessError::IOError(Some(err))),
            },
            None => Err(ScanAccessError::ScanNotFound),
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, ScanAccessError> {
        match self._offset_of_index(index) {
            Some(offset) => match self.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(ScanAccessError::IOError(Some(err))),
            },
            None => Err(ScanAccessError::ScanNotFound),
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, ScanAccessError> {
        match self._offset_of_time(time) {
            Some(offset) => match self.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(ScanAccessError::IOError(Some(err))),
            },
            None => Err(ScanAccessError::ScanNotFound),
        }
    }
}

impl<
        C: CentroidPeakAdapting + From<CentroidPeak>,
        D: DeconvolutedPeakAdapting + From<DeconvolutedPeak>,
    > MZFileReader<C, D, MultiLayerSpectrum<C, D>> for MGFReaderType<fs::File, C, D>
{
    fn open_file(source: fs::File) -> Self {
        Self::new(source)
    }

    fn construct_index_from_stream(&mut self) -> u64 {
        self.build_index()
    }
}

impl<
        R: io::Read,
        C: CentroidPeakAdapting + From<CentroidPeak>,
        D: DeconvolutedPeakAdapting + From<DeconvolutedPeak>,
    > MSDataFileMetadata for MGFReaderType<R, C, D>
{
    crate::impl_metadata_trait!();
}

pub type MGFReader<R> = MGFReaderType<R, CentroidPeak, DeconvolutedPeak>;

pub(crate) fn is_mgf(buf: &[u8]) -> bool {
    let needle = b"BEGIN IONS";
    if let Some(_loc) = buf
        .windows(needle.len())
        .position(|window| window == needle)
    {
        true
    } else {
        false
    }
}

pub struct MGFWriterType<
    R: io::Write,
    C: CentroidPeakAdapting + From<CentroidPeak> = CentroidPeak,
    D: DeconvolutedPeakAdapting + From<DeconvolutedPeak> = DeconvolutedPeak,
> {
    pub handle: io::BufWriter<R>,
    pub offset: usize,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
}

impl<
        W: io::Write,
        C: CentroidPeakAdapting + From<CentroidPeak>,
        D: DeconvolutedPeakAdapting + From<DeconvolutedPeak>,
    > MGFWriterType<W, C, D>
{
    pub fn new(file: W) -> MGFWriterType<W, C, D> {
        let handle = io::BufWriter::with_capacity(500, file);
        MGFWriterType {
            handle,
            offset: 0,
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
        }
    }

    pub fn into_inner(self) ->BufWriter<W> {
        self.handle
    }

    fn write_param<P: ParamLike>(&mut self, param: &P) -> io::Result<()> {
        self.handle.write_all(param.name().to_uppercase().as_bytes())?;
        self.handle.write_all(param.value().as_bytes())?;
        self.handle.write_all(b"\n")?;
        Ok(())
    }

    pub fn write_header<T: SpectrumBehavior<C, D>>(&mut self, spectrum: &T) -> io::Result<()> {
        let desc = spectrum.description();
        if desc.ms_level == 1 {
            log::warn!(
                "Attempted to write an MS1 spectrum to MGF, {}, skipping.",
                desc.id
            );
            return Ok(());
        }
        self.handle.write_all(
            br#"BEGIN IONS
TITLE="#,
        )?;
        self.handle.write_all(desc.id.as_bytes())?;
        self.handle.write_all(b"\nRTINSECONDS=")?;
        self
            .handle
            .write_all(spectrum.start_time().to_string().as_bytes())?;
        self.handle.write_all(b"\n")?;
        match &desc.precursor {
            Some(precursor) => {
                let ion = precursor.ion();
                self.handle.write_all(b"PEPMASS=")?;
                self.handle.write_all(ion.mz.to_string().as_bytes())?;
                self.handle.write_all(b" ")?;
                self.handle.write_all(ion.intensity.to_string().as_bytes())?;
                if let Some(charge) = ion.charge {
                    self.handle.write_all(b" ")?;
                    self.handle.write_all(charge.to_string().as_bytes())?;
                }
                self.handle.write_all(b"\n")?;

                for param in precursor
                    .ion()
                    .params()
                    .iter()
                    .chain(precursor.activation.params())
                {
                    self.write_param(param)?;
                }
                if let Some(pid) = precursor.precursor_id() {
                    self.handle.write_all(b"PRECURSORSCAN=")?;
                    self.handle.write_all(pid.as_bytes())?;
                    self.handle.write_all(b"\n")?;
                }
            }
            None => {}
        }
        for param in desc.params() {
            self.handle.write_all(param.name.to_uppercase().as_bytes())?;
            self.handle.write_all(param.value.as_bytes())?;
        }

        Ok(())
    }

    pub fn write_deconvoluted_centroids(
        &mut self,
        spectrum: &MultiLayerSpectrum<C, D>,
    ) -> io::Result<()> {
        match &spectrum.deconvoluted_peaks {
            Some(centroids) => {
                for peak in centroids.iter().map(|p| p.as_centroid()) {
                    self.handle.write_all(peak.mz().to_string().as_bytes())?;
                    self.handle.write_all(b" ")?;
                    self.handle.write_all(peak.intensity().to_string().as_bytes())?;
                    self.handle.write_all(b" ")?;
                    self.handle.write_all(peak.charge().to_string().as_bytes())?;
                    self.handle.write_all(b"\n")?;
                }
            }
            None => {}
        }
        Ok(())
    }

    pub fn write_centroids(&mut self, spectrum: &MultiLayerSpectrum<C, D>) -> io::Result<()> {
        match &spectrum.peaks {
            Some(centroids) => {
                for peak in centroids {
                    self.handle.write_all(peak.mz().to_string().as_bytes())?;
                    self.handle.write_all(b" ")?;
                    self.handle.write_all(peak.intensity().to_string().as_bytes())?;
                    self.handle.write_all(b"\n")?;
                }
            }
            None => {}
        }
        Ok(())
    }
}

impl<
        'a,
        W: io::Write,
        C: CentroidPeakAdapting + From<CentroidPeak> + 'static,
        D: DeconvolutedPeakAdapting + From<DeconvolutedPeak> + 'static,
    > ScanWriter<'a, C, D> for MGFWriterType<W, C, D>
{
    fn write(&mut self, spectrum: &MultiLayerSpectrum<C, D>) -> io::Result<usize> {
        let desc = spectrum.description();
        if desc.ms_level == 1 {
            log::warn!(
                "Attempted to write an MS1 spectrum to MGF, {}, skipping.",
                desc.id
            );
            return Ok(0);
        }
        self.write_header(spectrum)?;
        if spectrum.deconvoluted_peaks.is_some() {
            self.write_deconvoluted_centroids(spectrum)?;
        } else if spectrum.peaks.is_some() {
            self.write_centroids(spectrum)?;
        } 
        self.handle.write_all(b"END IONS\n")?;
        Ok(0)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.handle.flush()
    }
}

pub type MGFWriter<W> = MGFWriterType<W, CentroidPeak, DeconvolutedPeak>;

#[cfg(test)]
mod test {
    use mzpeaks::IndexedCoordinate;

    use crate::CentroidSpectrum;

    use super::*;
    use std::fs;
    use std::path;

    #[test]
    fn test_reader() {
        let path = path::Path::new("./test/data/small.mgf");
        let file = fs::File::open(path).expect("Test file doesn't exist");
        let reader = MGFReaderType::<_, CentroidPeak, DeconvolutedPeak>::new(file);
        let mut ms1_count = 0;
        let mut msn_count = 0;
        for scan in reader {
            let level = scan.ms_level();
            if level == 1 {
                ms1_count += 1;
            } else {
                msn_count += 1;
            }
        }
        assert_eq!(ms1_count, 0);
        assert_eq!(msn_count, 34);
    }

    #[test]
    fn test_reader_indexed() {
        let path = path::Path::new("./test/data/small.mgf");
        let file = fs::File::open(path).expect("Test file doesn't exist");
        let mut reader = MGFReaderType::<_, CentroidPeak, DeconvolutedPeak>::new_indexed(file);

        let n = reader.len();
        let mut ms1_count = 0;
        let mut msn_count = 0;

        for i in (0..n).rev() {
            let scan = reader.get_spectrum_by_index(i).expect("Missing spectrum");
            let level = scan.ms_level();
            if level == 1 {
                ms1_count += 1;
            } else {
                msn_count += 1;
            }
            let centroided: CentroidSpectrum = scan.try_into().unwrap();
            centroided.peaks.iter().for_each(|p| {
                (centroided.peaks[p.get_index() as usize]).mz();
            })
        }
        assert_eq!(ms1_count, 0);
        assert_eq!(msn_count, 34);
    }

    #[test]
    fn test_writer() -> io::Result<()> {
        let buff: Vec<u8> = Vec::new();
        let inner_writer = io::Cursor::new(buff);
        let mut writer = MGFWriter::new(inner_writer);

        let path = path::Path::new("./test/data/small.mgf");
        let file = fs::File::open(path).expect("Test file doesn't exist");
        let mut reader = MGFReader::new(file);

        for scan in reader.iter() {
            writer.write(&scan)?;
        }
        writer.flush()?;
        let inner_writer = writer.handle.into_inner()?;
        let buffer = inner_writer.into_inner();
        let reader2 = MGFReader::new(io::Cursor::new(buffer));
        assert_eq!(reader2.len(), reader.len());


        // Not including platform-specific line endings
        Ok(())
    }
}
