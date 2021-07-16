use std::io::prelude::*;
use std::io;

use regex::Regex;
use lazy_static::lazy_static;

use crate::peak::{CentroidPeak};
use crate::peak_set::{PeakSet, PeakCollection};
use crate::spectrum::{Precursor, SpectrumDescription,
                      SelectedIon, CentroidSpectrum, scan_properties};


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


pub struct MGFReader<R: io::Read> {
    pub handle: io::BufReader<R>,
    pub state: MGFParserState,
    pub offset: usize,
    pub error: String,
}


impl<R: io::Read> MGFReader<R>{
    fn parse_peak_from_line(&mut self, line: &str) -> Option<CentroidPeak> {
        let mut chars = line.chars();
        let first = chars.next().unwrap();
        if first.is_numeric() {
            // A lazily created static regular expression to parse peak separators
            lazy_static! {
                static ref PEAK_SEPERATOR: Regex = Regex::new(r"\t|\s+").unwrap();
            }
            let parts: Vec<&str> = PEAK_SEPERATOR.split(line).collect();
            let nparts = parts.len();
            if nparts < 2 {
                self.state = MGFParserState::Error;
                self.error = String::from("Incorrect number of numerical columns");
            }
            let mz: f64 = parts[0].parse().unwrap();
            let intensity: f32 = parts[1].parse().unwrap();
            return Some(CentroidPeak {mz, intensity, ..Default::default()});
        }
        return None
    }

    fn handle_scan_header(&mut self, line: &str, description: &mut SpectrumDescription, peaks: &mut PeakSet) -> bool {
        let peak_line = match self.parse_peak_from_line(line) {
            Some(peak) => {peaks.push(peak); true}
            None => false
        };
        if peak_line {
            self.state = MGFParserState::Peaks;
            return true
        }
        else if line == "END IONS" {
            self.state = MGFParserState::Between;
            return true
        }
        else if line.contains('=') {
            let parts: Vec<&str> = line.splitn(2, '=').collect();
            let key = parts[0];
            let value = parts[1];
            match key {
                "TITLE" => {description.id = String::from(value)},
                "RTINSECONDS" => {
                    let scan_ev = description.acquisition.first_scan_mut().expect("Automatically adds scan event");
                    scan_ev.start_time = value.parse().unwrap()
                },
                "PEPMASS" => {
                    let parts: Vec<&str> = value.split_ascii_whitespace().collect();
                    let mz: f64 = parts[0].parse().unwrap();
                    let intensity: f32 = parts[1].parse().unwrap();
                    let mut charge: i32 = 0;

                    if parts.len() > 2 {
                        charge = parts[2].parse().unwrap();
                    }
                    description.precursor = Some(Precursor {
                        ion: SelectedIon {mz, intensity, charge, ..Default::default()},
                        ..Default::default()});
                }
                &_ => {
                    description.annotations.insert(String::from(key.to_lowercase()), String::from(value));
                }
            };

            return true;
        } else {
            self.state = MGFParserState::Error;
            self.error = String::from(format!("Unexpected content {} in scan header", line));
            return false;
        }
    }

    fn handle_peak(&mut self, line: &str, peaks: &mut PeakSet) -> bool{
        let peak_line = match self.parse_peak_from_line(line) {
            Some(peak) => {peaks.push(peak); return true}
            None => false
        };
        if peak_line {
            return true
        }
        else if line == "END IONS" {
            self.state = MGFParserState::Between;
            return false;
        }
        else {
            self.state = MGFParserState::Error;
            self.error = String::from(format!("Unexpected content {} in peak list", line));
            return false;
        }
    }

    fn handle_start(&mut self, line: &str) -> bool {
        if line.contains("=") {}
        else if line == "BEGIN IONS" {
            self.state = MGFParserState::ScanHeaders;
        }
        return true;
    }

    fn handle_between(&mut self, line: &str) -> bool {
        if line == "BEGIN IONS" {
            self.state = MGFParserState::ScanHeaders;
        }
        return true;
    }

    pub fn new_scan(&self) -> CentroidSpectrum {
        let description: SpectrumDescription = SpectrumDescription {
            ms_level: 2,
            is_profile: scan_properties::ScanSiganlContinuity::Centroid,
            polarity: scan_properties::ScanPolarity::Unknown,
            ..Default::default() };

        let peaks: PeakSet = PeakSet::empty();
        let scan = CentroidSpectrum {description, peaks};
        return scan;
    }

    fn read_line(&mut self, buffer: &mut String) -> io::Result<usize> {
        return self.handle.read_line(buffer);
    }

    pub fn read_next(&mut self) -> Option<CentroidSpectrum> {
        let mut scan = self.new_scan();
        match self.read_into(&mut scan) {
            Ok(offset) => if offset > 0 { Some(scan) } else { None },
            Err(message) => {
                panic!("{}", message);
            }
        }
    }

    pub fn read_into(&mut self, spectrum: &mut CentroidSpectrum) -> Result<usize, String>{
        let mut buffer = String::new();
        let mut work = true;
        let mut offset: usize = 0;
        let description = &mut spectrum.description;
        let peaks = &mut spectrum.peaks;

        while work {
            buffer.truncate(0);
            let b = match self.read_line(&mut buffer) {
                Ok(b) => {
                    if b == 0 {
                        work = false;
                    }
                    b
                },
                Err(err) => {
                    return Err(String::from(format!("Error while reading file: {}", err)));
                }
            };
            offset += b;
            if b == 0 {
                self.state = MGFParserState::Done;
                break
            }
            let line = buffer.trim();
            let n = line.len();
            if n == 0 {
                continue;
            }
            if self.state == MGFParserState::Start {
                work = self.handle_start(line);
            }
            else if self.state == MGFParserState::Between {
                work = self.handle_between(line);
            }
            else if self.state == MGFParserState::ScanHeaders {
                work = self.handle_scan_header(line, description, peaks)
            }
            else if self.state == MGFParserState::Peaks {
                work = self.handle_peak(line, peaks);
            }
            if matches!(self.state, MGFParserState::Error) {
                panic!("MGF Parsing Error: {}", self.error);
            }
        }
        return Ok(offset);
    }

    pub fn new(file: R) -> MGFReader<R> {
        let handle = io::BufReader::with_capacity(500, file);
        return MGFReader {handle, state: MGFParserState::Start, offset: 0, error: String::from("")}
    }
}

impl<R: io::Read> Iterator for MGFReader<R> {
    type Item = CentroidSpectrum;

    fn next(&mut self) -> Option<Self::Item> {
        return self.read_next();
    }
}