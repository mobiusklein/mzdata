use std::fs;
use std::io;
use std::io::prelude::*;
use std::io::SeekFrom;
use std::str;

use log::warn;

use lazy_static::lazy_static;
use regex::Regex;

use crate::params::{Param, ParamDescribed};
use crate::peaks::{CentroidPeak, PeakCollection, PeakSet};
use crate::spectrum::{
    scan_properties, CentroidSpectrum, Precursor, SelectedIon, Spectrum,
    SpectrumDescription,
};

use super::offset_index::OffsetIndex;
use super::traits::{
    MZFileReader, RandomAccessScanIterator, ScanAccessError, ScanSource, SeekRead,
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
struct SpectrumBuilder {
    pub peaks: PeakSet,
    pub description: SpectrumDescription,
}

impl Default for SpectrumBuilder {
    fn default() -> SpectrumBuilder {
        SpectrumBuilder {
            peaks: PeakSet::default(),
            description: SpectrumDescription {
                ms_level: 2,
                signal_continuity: scan_properties::SignalContinuity::Centroid,
                polarity: scan_properties::ScanPolarity::Unknown,
                ..Default::default()
            },
        }
    }
}

impl SpectrumBuilder {
    pub fn into_spectrum(self, spectrum: &mut Spectrum) {
        spectrum.description = self.description;
        spectrum.peaks = Some(self.peaks);
    }
}

impl From<SpectrumBuilder> for Spectrum {
    fn from(builder: SpectrumBuilder) -> Spectrum {
        let mut spectrum = Spectrum::default();
        builder.into_spectrum(&mut spectrum);
        spectrum
    }
}

macro_rules! impl_from_spectrum_builder_for_spec {
    ($($t:ty), +) => {$(
        impl From<SpectrumBuilder> for $t {
            fn from(builder: SpectrumBuilder) -> $t {
                let spec: Spectrum = builder.into();
                spec.into()
            }
        }
    )+};
}

impl_from_spectrum_builder_for_spec!(CentroidSpectrum);

/// An MGF (Mascot Generic Format) file parser that supports iteration and random access.
/// The parser produces [`Spectrum`] instances. These may be converted directly into [`CentroidSpectrum`]
/// instances using [`Spectrum.into_centroid`] or the [`From`] trait which better represent
/// the nature of this preprocessed data type.
pub struct MGFReader<R: io::Read> {
    pub handle: io::BufReader<R>,
    pub state: MGFParserState,
    pub offset: usize,
    pub error: MGFError,
    pub index: OffsetIndex,
}

// A lazily created static regular expression to parse peak separators
lazy_static! {
    static ref PEAK_SEPERATOR: Regex = Regex::new(r"\t|\s+").unwrap();
}

impl<R: io::Read> MGFReader<R> {
    fn parse_peak_from_line(&mut self, line: &str) -> Option<CentroidPeak> {
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
            return Some(CentroidPeak {
                mz,
                intensity,
                ..Default::default()
            });
        }
        None
    }

    fn handle_scan_header(
        &mut self,
        line: &str,
        description: &mut SpectrumDescription,
        peaks: &mut PeakSet,
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
                    description.add_param(Param::key_value(String::from(key.to_lowercase()), String::from(value)));
                }
            };

            true
        } else {
            self.state = MGFParserState::Error;
            self.error = MGFError::MalformedHeaderLine;
            false
        }
    }

    fn handle_peak(&mut self, line: &str, peaks: &mut PeakSet) -> bool {
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
    pub fn read_next(&mut self) -> Option<Spectrum> {
        let mut builder = SpectrumBuilder::default();
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
    fn _parse_into(&mut self, spectrum: &mut SpectrumBuilder) -> Result<usize, MGFError> {
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
    pub fn read_into(&mut self, spectrum: &mut Spectrum) -> Result<usize, MGFError> {
        let mut accumulator = SpectrumBuilder::default();
        match self._parse_into(&mut accumulator) {
            Ok(sz) => {
                accumulator.into_spectrum(spectrum);
                Ok(sz)
            }
            Err(err) => Err(err),
        }
    }
    /// Create a new, unindexed MGF parser
    pub fn new(file: R) -> MGFReader<R> {
        let handle = io::BufReader::with_capacity(500, file);
        MGFReader {
            handle,
            state: MGFParserState::Start,
            offset: 0,
            error: MGFError::NoError,
            index: OffsetIndex::new("spectrum".to_owned()),
        }
    }
}

impl<R: io::Read> Iterator for MGFReader<R> {
    type Item = Spectrum;

    /// Read the next spectrum from the file.
    fn next(&mut self) -> Option<Self::Item> {
        self.read_next()
    }
}

impl<R: SeekRead> MGFReader<R> {
    /// Construct a new MGFReader and build an offset index
    /// using [`Self::build_index`]
    pub fn new_indexed(file: R) -> MGFReader<R> {
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

impl<R: SeekRead> ScanSource<Spectrum> for MGFReader<R> {
    /// Retrieve a spectrum by it's native ID
    fn get_spectrum_by_id(&mut self, id: &str) -> Option<Spectrum> {
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
            },
            None => None
        }
    }

    /// Retrieve a spectrum by it's integer index
    fn get_spectrum_by_index(&mut self, index: usize) -> Option<Spectrum> {
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
            },
            None => None
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
            warn!("Attempting to use an uninitialized offset index on MGFReader")
        }
        &self.index
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.index = index;
    }
}

impl<R: SeekRead> RandomAccessScanIterator<Spectrum> for MGFReader<R> {
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

impl MZFileReader<Spectrum> for MGFReader<fs::File> {
    fn open_file(source: fs::File) -> Self {
        Self::new(source)
    }

    fn construct_index_from_stream(&mut self) -> u64 {
        self.build_index()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::spectrum::spectrum::SpectrumBehavior;
    use std::fs;
    use std::path;

    #[test]
    fn test_reader() {
        let path = path::Path::new("./test/data/small.mgf");
        let file = fs::File::open(path).expect("Test file doesn't exist");
        let reader = MGFReader::new(file);
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
        let mut reader = MGFReader::new_indexed(file);

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
}
