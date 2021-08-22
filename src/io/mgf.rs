use std::collections::HashMap;
use std::convert::TryInto;
use std::fs;
use std::io;
use std::io::prelude::*;
use std::io::SeekFrom;
use std::marker::PhantomData;
use std::str;

use log::warn;

use lazy_static::lazy_static;
use regex::Regex;

use mzpeaks::peak::KnownCharge;
use mzpeaks::{
    CentroidPeak, DeconvolutedPeak, IntensityMeasurement, MZLocated, MZPeakSetType, PeakCollection,
};

use crate::meta::{
    DataProcessing, FileDescription, InstrumentConfiguration, MSDataFileMetadata, Software,
};
use crate::params::{Param, ParamDescribed};
use crate::spectrum::spectrum::{
    CentroidPeakAdapting, CentroidSpectrumType, DeconvolutedPeakAdapting, MultiLayerSpectrum,
};
use crate::spectrum::{
    scan_properties, Precursor, PrecursorSelection, SelectedIon, SpectrumBehavior,
    SpectrumDescription,
};

use super::offset_index::OffsetIndex;
use super::traits::{
    MZFileReader, RandomAccessScanIterator, ScanAccessError, ScanSource, ScanWriter, SeekRead,
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

#[derive(Debug, Clone, Copy)]
pub enum MGFError {
    NoError,
    MalformedPeakLine,
    MalformedHeaderLine,
    TooManyColumnsForPeakLine,
    IOError,
}

#[derive(Debug, Clone)]
struct SpectrumBuilder<
    C: CentroidPeakAdapting = CentroidPeak,
    D: DeconvolutedPeakAdapting = DeconvolutedPeak,
> {
    pub peaks: MZPeakSetType<C>,
    pub description: SpectrumDescription,
    deconvoluted_type: PhantomData<D>,
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> Default for SpectrumBuilder<C, D> {
    fn default() -> SpectrumBuilder<C, D> {
        SpectrumBuilder {
            peaks: MZPeakSetType::<C>::default(),
            description: SpectrumDescription {
                ms_level: 2,
                signal_continuity: scan_properties::SignalContinuity::Centroid,
                polarity: scan_properties::ScanPolarity::Unknown,
                ..Default::default()
            },
            deconvoluted_type: PhantomData,
        }
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> SpectrumBuilder<C, D> {
    pub fn into_spectrum(self, spectrum: &mut MultiLayerSpectrum<C, D>) {
        spectrum.description = self.description;
        spectrum.peaks = Some(self.peaks);
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> From<SpectrumBuilder<C, D>>
    for MultiLayerSpectrum<C, D>
{
    fn from(builder: SpectrumBuilder<C, D>) -> MultiLayerSpectrum<C, D> {
        let mut spectrum = MultiLayerSpectrum::<C, D>::default();
        builder.into_spectrum(&mut spectrum);
        spectrum
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> From<SpectrumBuilder<C, D>>
    for CentroidSpectrumType<C>
{
    fn from(builder: SpectrumBuilder<C, D>) -> CentroidSpectrumType<C> {
        let spec: MultiLayerSpectrum<C, D> = builder.into();
        spec.try_into().unwrap()
    }
}

/// An MGF (Mascot Generic Format) file parser that supports iteration and random access.
/// The parser produces [`Spectrum`] instances. These may be converted directly into [`CentroidSpectrum`]
/// instances using [`Spectrum.into_centroid`] or the [`From`] trait which better represent
/// the nature of this preprocessed data type.
pub struct MGFReaderType<
    R: io::Read,
    C: CentroidPeakAdapting = CentroidPeak,
    D: DeconvolutedPeakAdapting = DeconvolutedPeak,
> {
    pub handle: io::BufReader<R>,
    pub state: MGFParserState,
    pub offset: usize,
    pub error: MGFError,
    pub index: OffsetIndex,
    file_description: FileDescription,
    instrument_configurations: HashMap<String, InstrumentConfiguration>,
    softwares: Vec<Software>,
    data_processings: Vec<DataProcessing>,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
}

// A lazily created static regular expression to parse peak separators
lazy_static! {
    static ref PEAK_SEPERATOR: Regex = Regex::new(r"\t|\s+").unwrap();
}

impl<R: io::Read, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MGFReaderType<R, C, D> {
    fn parse_peak_from_line(&mut self, line: &str) -> Option<C> {
        let mut chars = line.chars();
        let first = chars.next().unwrap();
        if first.is_numeric() {
            let parts: Vec<&str> = PEAK_SEPERATOR.split(line).collect();
            let nparts = parts.len();
            if nparts < 2 {
                self.state = MGFParserState::Error;
                self.error = MGFError::TooManyColumnsForPeakLine;
            }
            let mz: f64 = parts[0].parse().unwrap();
            let intensity: f32 = parts[1].parse().unwrap();
            return Some(
                CentroidPeak {
                    mz,
                    intensity,
                    ..Default::default()
                }
                .into(),
            );
        }
        None
    }

    fn handle_scan_header(
        &mut self,
        line: &str,
        description: &mut SpectrumDescription,
        peaks: &mut MZPeakSetType<C>,
    ) -> bool {
        let peak_line = match self.parse_peak_from_line(line) {
            Some(peak) => {
                peaks.push(peak);
                true
            }
            None => false,
        };
        if peak_line {
            self.state = MGFParserState::Peaks;
            true
        } else if line == "END IONS" {
            self.state = MGFParserState::Between;
            true
        } else if line.contains('=') {
            let parts: Vec<&str> = line.splitn(2, '=').collect();
            let key = parts[0];
            let value = parts[1];
            match key {
                "TITLE" => description.id = String::from(value),
                "RTINSECONDS" => {
                    let scan_ev = description
                        .acquisition
                        .first_scan_mut()
                        .expect("Automatically adds scan event");
                    scan_ev.start_time = value.parse().unwrap()
                }
                "PEPMASS" => {
                    let parts: Vec<&str> = value.split_ascii_whitespace().collect();
                    let mz: f64 = parts[0].parse().unwrap();
                    let intensity: f32 = parts[1].parse().unwrap();
                    let mut charge: Option<i32> = None;

                    if parts.len() > 2 {
                        charge = Some(parts[2].parse().unwrap());
                    }
                    description.precursor = Some(Precursor {
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
                    description.add_param(Param::new_key_value(
                        String::from(key.to_lowercase()),
                        String::from(value),
                    ));
                }
            };

            true
        } else {
            self.state = MGFParserState::Error;
            self.error = MGFError::MalformedHeaderLine;
            false
        }
    }

    fn handle_peak(&mut self, line: &str, peaks: &mut MZPeakSetType<C>) -> bool {
        let peak_line = match self.parse_peak_from_line(line) {
            Some(peak) => {
                peaks.push(peak);
                return true;
            }
            None => false,
        };
        if peak_line {
            true
        } else if line == "END IONS" {
            self.state = MGFParserState::Between;
            false
        } else {
            self.state = MGFParserState::Error;
            self.error = MGFError::MalformedPeakLine;
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
        let mut builder = SpectrumBuilder::<C, D>::default();
        match self._parse_into(&mut builder) {
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
    fn _parse_into(&mut self, spectrum: &mut SpectrumBuilder<C, D>) -> Result<usize, MGFError> {
        let mut buffer = String::new();
        let mut work = true;
        let mut offset: usize = 0;
        let description = &mut spectrum.description;
        let peaks = &mut spectrum.peaks;

        while work {
            buffer.clear();
            let b = match self.read_line(&mut buffer) {
                Ok(b) => {
                    if b == 0 {
                        work = false;
                    }
                    b
                }
                Err(_err) => {
                    self.error = MGFError::IOError;
                    self.state = MGFParserState::Error;
                    return Err(self.error);
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
                work = self.handle_scan_header(line, description, peaks)
            } else if self.state == MGFParserState::Peaks {
                work = self.handle_peak(line, peaks);
            }
            if matches!(self.state, MGFParserState::Error) {
                panic!("MGF Parsing Error: {:?}", self.error);
            }
        }
        Ok(offset)
    }
    pub fn read_into(
        &mut self,
        spectrum: &mut MultiLayerSpectrum<C, D>,
    ) -> Result<usize, MGFError> {
        let mut accumulator = SpectrumBuilder::default();
        match self._parse_into(&mut accumulator) {
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
        term.accession = "MS:1000580".to_owned();
        term.controlled_vocabulary = Some("MS".to_owned());
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
            error: MGFError::NoError,
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
        if self.index.len() == 0 {
            warn!("An index was built but no entries were found")
        }
        offset
    }
}

impl<R: SeekRead, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    ScanSource<C, D, MultiLayerSpectrum<C, D>> for MGFReaderType<R, C, D>
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
    fn reset(&mut self) -> &Self {
        self.seek(SeekFrom::Start(0))
            .expect("Failed to reset file stream");
        self
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
    RandomAccessScanIterator<C, D, MultiLayerSpectrum<C, D>> for MGFReaderType<R, C, D>
{
    fn start_from_id(&mut self, id: &str) -> Result<&Self, ScanAccessError> {
        match self._offset_of_id(id) {
            Some(offset) => match self.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(ScanAccessError::IOError(Some(err))),
            },
            None => Err(ScanAccessError::ScanNotFound),
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&Self, ScanAccessError> {
        match self._offset_of_index(index) {
            Some(offset) => match self.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(ScanAccessError::IOError(Some(err))),
            },
            None => Err(ScanAccessError::ScanNotFound),
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&Self, ScanAccessError> {
        match self._offset_of_time(time) {
            Some(offset) => match self.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(ScanAccessError::IOError(Some(err))),
            },
            None => Err(ScanAccessError::ScanNotFound),
        }
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    MZFileReader<C, D, MultiLayerSpectrum<C, D>> for MGFReaderType<fs::File, C, D>
{
    fn open_file(source: fs::File) -> Self {
        Self::new(source)
    }

    fn construct_index_from_stream(&mut self) -> u64 {
        self.build_index()
    }
}

impl<R: io::Read, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MSDataFileMetadata
    for MGFReaderType<R, C, D>
{
    crate::impl_metadata_trait!();
}

pub type MGFReader<R> = MGFReaderType<R, CentroidPeak, DeconvolutedPeak>;

pub struct MGFWriterType<
    R: io::Write,
    C: CentroidPeakAdapting = CentroidPeak,
    D: DeconvolutedPeakAdapting = DeconvolutedPeak,
> {
    pub handle: io::BufWriter<R>,
    pub offset: usize,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
}

impl<W: io::Write, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MGFWriterType<W, C, D> {
    pub fn new(file: W) -> MGFWriterType<W, C, D> {
        let handle = io::BufWriter::with_capacity(500, file);
        MGFWriterType {
            handle,
            offset: 0,
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
        }
    }

    pub fn write_header(&mut self, spectrum: &MultiLayerSpectrum<C, D>) -> io::Result<usize> {
        let desc = spectrum.description();
        if desc.ms_level == 1 {
            log::warn!(
                "Attempted to write an MS1 spectrum to MGF, {}, skipping.",
                desc.id
            );
            return Ok(0);
        }
        let mut count = self.handle.write(
            br#"BEGIN IONS
TITLE="#,
        )?;
        count += self.handle.write(desc.id.as_bytes())?;
        count += self.handle.write(b"\nRTINSECONDS=")?;
        // let s = spectrum.start_time().to_string();
        count += self
            .handle
            .write(spectrum.start_time().to_string().as_bytes())?;
        count += self.handle.write(b"\n")?;
        match &desc.precursor {
            Some(precursor) => {
                let ion = precursor.ion();
                count += self.handle.write(b"PEPMASS=")?;
                count += self.handle.write(ion.mz.to_string().as_bytes())?;
                count += self.handle.write(b" ")?;
                count += self.handle.write(ion.intensity.to_string().as_bytes())?;
                match ion.charge {
                    Some(charge) => {
                        count += self.handle.write(b" ")?;
                        count += self.handle.write(charge.to_string().as_bytes())?;
                    }
                    None => {}
                }
                count += self.handle.write(b"\n")?;
                for param in precursor.params() {
                    count += self.handle.write(param.name.to_uppercase().as_bytes())?;
                    count += self.handle.write(param.value.as_bytes())?;
                }
            }
            None => {}
        }
        for param in desc.params() {
            count += self.handle.write(param.name.to_uppercase().as_bytes())?;
            count += self.handle.write(param.value.as_bytes())?;
        }

        Ok(count)
    }

    pub fn write_deconvoluted_centroids(
        &mut self,
        spectrum: &MultiLayerSpectrum<C, D>,
    ) -> io::Result<usize> {
        let mut count = 0;
        match &spectrum.deconvoluted_peaks {
            Some(centroids) => {
                for peak in centroids.iter().map(|p| p.as_centroid()) {
                    count += self.handle.write(peak.mz().to_string().as_bytes())?;
                    count += self.handle.write(b" ")?;
                    count += self.handle.write(peak.intensity().to_string().as_bytes())?;
                    count += self.handle.write(b" ")?;
                    count += self.handle.write(peak.charge().to_string().as_bytes())?;
                    count += self.handle.write(b"\n")?;
                }
            }
            None => {}
        }
        Ok(count)
    }

    pub fn write_centroids(&mut self, spectrum: &MultiLayerSpectrum<C, D>) -> io::Result<usize> {
        let mut count = 0;
        match &spectrum.peaks {
            Some(centroids) => {
                for peak in centroids {
                    count += self.handle.write(peak.mz().to_string().as_bytes())?;
                    count += self.handle.write(b" ")?;
                    count += self.handle.write(peak.intensity().to_string().as_bytes())?;
                    count += self.handle.write(b"\n")?;
                }
            }
            None => {}
        }
        Ok(count)
    }
}

impl<
        'a,
        W: io::Write,
        C: CentroidPeakAdapting + 'static,
        D: DeconvolutedPeakAdapting + 'static,
    > ScanWriter<'a, W, C, D> for MGFWriterType<W, C, D>
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
        let mut count = self.write_header(spectrum)?;
        if matches!(&spectrum.deconvoluted_peaks, Some(_)) {
            count += self.write_deconvoluted_centroids(spectrum)?;
            count += self.handle.write(b"END IONS\n")?;
        } else if matches!(&spectrum.peaks, Some(_)) {
            count += self.write_centroids(spectrum)?;
            count += self.handle.write(b"END IONS\n")?;
        } else {
            count += self.handle.write(b"END IONS\n")?;
        }
        self.offset += count;
        Ok(count)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.handle.flush()
    }
}

pub type MGFWriter<W> = MGFWriterType<W, CentroidPeak, DeconvolutedPeak>;

#[cfg(test)]
mod test {
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
        let reader = MGFReader::new(file);

        let mut count = 0;
        for scan in reader {
            count += writer.write(&scan)?;
        }
        writer.flush()?;
        let inner_writer = writer.handle.into_inner()?;
        let buffer = inner_writer.get_ref();
        let text = str::from_utf8(buffer).unwrap();
        // Not including platform-specific line endings
        assert_eq!(count, 548732);
        assert_eq!(count, text.len());
        Ok(())
    }
}
