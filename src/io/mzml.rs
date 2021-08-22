//! Implements a parser for the PSI-MS mzML and indexedmzML XML file formats
//! for representing raw and processed mass spectra.

use std::collections::HashMap;
use std::convert::TryInto;
use std::fmt::Debug;
use std::fs;
use std::io;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::marker::PhantomData;

use log::warn;

use quick_xml::events::BytesDecl;
use quick_xml::events::{BytesEnd, BytesStart, BytesText, Event};
use quick_xml::Error as XMLError;
use quick_xml::{Reader, Writer};

use super::offset_index::OffsetIndex;
use super::traits::{
    MZFileReader, RandomAccessScanIterator, ScanAccessError, ScanSource, ScanWriter, SeekRead,
};
use super::utils::MD5HashingStream;

use mzpeaks::{peak_set::PeakSetVec, CentroidPeak, DeconvolutedPeak, Mass, MZ};

use crate::meta::file_description::{FileDescription, SourceFile};
use crate::meta::instrument::{Component, ComponentType, InstrumentConfiguration};
use crate::meta::{DataProcessing, MSDataFileMetadata, ProcessingMethod, Software};
use crate::params::{ControlledVocabulary, Param, ParamList};
use crate::spectrum::scan_properties::*;
use crate::spectrum::signal::{
    ArrayType, BinaryArrayMap, BinaryCompressionType, BinaryDataArrayType, DataArray,
};
use crate::spectrum::spectrum::{
    CentroidPeakAdapting, CentroidSpectrumType, DeconvolutedPeakAdapting, MultiLayerSpectrum,
    RawSpectrum, Spectrum,
};
use crate::ParamDescribed;
use crate::SpectrumBehavior;

pub type Bytes = Vec<u8>;

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum MzMLParserState {
    Start = 0,
    Resume,

    // Top-level metadata
    CVList,
    FileDescription,
    FileContents,
    SourceFileList,
    SourceFile,

    ReferenceParamGroupList,
    ReferenceParamGroup,

    SoftwareList,
    Software,

    InstrumentConfigurationList,
    InstrumentConfiguration,
    ComponentList,
    Source,
    Analyzer,
    Detector,

    DataProcessingList,
    DataProcessing,
    ProcessingMethod,

    Run,

    // Spectrum and Chromatogram List Elements
    Spectrum,
    SpectrumDone,

    SpectrumList,
    SpectrumListDone,

    BinaryDataArrayList,
    BinaryDataArray,
    Binary,

    ScanList,
    Scan,
    ScanWindowList,
    ScanWindow,

    PrecursorList,
    Precursor,
    IsolationWindow,
    SelectedIonList,
    SelectedIon,
    Activation,

    Chromatogram,
    ChromatogramDone,

    ParserError,
}

pub trait XMLParseBase {
    fn handle_xml_error(&self, error: quick_xml::Error, state: MzMLParserState) -> MzMLParserError {
        MzMLParserError::XMLError(state, format!("{:?}", error))
    }
}

pub trait CVParamParse: XMLParseBase {
    fn handle_param<B: io::BufRead>(
        &self,
        event: &BytesStart,
        reader: &Reader<B>,
        state: MzMLParserState,
    ) -> Result<Param, MzMLParserError> {
        let mut param = Param::new();
        for attr_parsed in event.attributes() {
            match attr_parsed {
                Ok(attr) => match attr.key {
                    b"name" => {
                        param.name = attr.unescape_and_decode_value(reader).expect(&format!(
                            "Error decoding CV param name at {}",
                            reader.buffer_position()
                        ));
                    }
                    b"value" => {
                        param.value = attr.unescape_and_decode_value(reader).expect(&format!(
                            "Error decoding CV param value at {}",
                            reader.buffer_position()
                        ));
                    }
                    b"cvRef" => {
                        param.controlled_vocabulary =
                            Some(attr.unescape_and_decode_value(reader).expect(&format!(
                                "Error decoding CV param reference at {}",
                                reader.buffer_position()
                            )));
                    }
                    b"accession" => {
                        param.accession = attr.unescape_and_decode_value(reader).expect(&format!(
                            "Error decoding CV param accession at {}",
                            reader.buffer_position()
                        ));
                    }
                    b"unitName" => {
                        param.unit_name =
                            Some(attr.unescape_and_decode_value(reader).expect(&format!(
                                "Error decoding CV param unit name at {}",
                                reader.buffer_position()
                            )));
                    }
                    b"unitAccession" => {
                        param.unit_accession =
                            Some(attr.unescape_and_decode_value(reader).expect(&format!(
                                "Error decoding CV param unit name at {}",
                                reader.buffer_position()
                            )));
                    }
                    b"unitCvRef" => {}
                    _ => {}
                },
                Err(msg) => return Err(self.handle_xml_error(msg, state)),
            }
        }
        Ok(param)
    }
}

#[derive(Debug, Clone)]
pub enum MzMLParserError {
    NoError,
    UnknownError(MzMLParserState),
    IncompleteSpectrum,
    IncompleteElementError(String, MzMLParserState),
    XMLError(MzMLParserState, String),
}

impl Default for MzMLParserError {
    fn default() -> MzMLParserError {
        MzMLParserError::NoError
    }
}

const BUFFER_SIZE: usize = 10000;

#[derive(Default)]
struct MzMLSpectrumBuilder<
    C: CentroidPeakAdapting = CentroidPeak,
    D: DeconvolutedPeakAdapting = DeconvolutedPeak,
> {
    pub params: ParamList,
    pub acquisition: Acquisition,
    pub precursor: Precursor,

    pub arrays: BinaryArrayMap,
    pub current_array: DataArray,

    pub index: usize,
    pub scan_id: String,
    pub ms_level: u8,
    pub polarity: ScanPolarity,
    pub signal_continuity: SignalContinuity,
    pub has_precursor: bool,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> XMLParseBase
    for MzMLSpectrumBuilder<C, D>
{
}
impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> CVParamParse
    for MzMLSpectrumBuilder<C, D>
{
}

pub trait SpectrumBuilding<
    C: CentroidPeakAdapting,
    D: DeconvolutedPeakAdapting,
    S: SpectrumBehavior<C, D>,
>
{
    fn isolation_window_mut(&mut self) -> &mut IsolationWindow;
    fn scan_window_mut(&mut self) -> &mut ScanWindow;
    fn selected_ion_mut(&mut self) -> &mut SelectedIon;
    fn current_array_mut(&mut self) -> &mut DataArray;
    fn into_spectrum(self, spectrum: &mut S);

    fn fill_binary_data_array(&mut self, param: Param) {
        match param.name.as_ref() {
            // Compression types
            "zlib compression" => {
                self.current_array_mut().compression = BinaryCompressionType::Zlib;
            }
            "no compression" => {
                self.current_array_mut().compression = BinaryCompressionType::NoCompression;
            }

            // Array data types
            "64-bit float" => {
                self.current_array_mut().dtype = BinaryDataArrayType::Float64;
            }
            "32-bit float" => {
                self.current_array_mut().dtype = BinaryDataArrayType::Float32;
            }
            "64-bit integer" => {
                self.current_array_mut().dtype = BinaryDataArrayType::Int64;
            }
            "32-bit integer" => {
                self.current_array_mut().dtype = BinaryDataArrayType::Int32;
            }
            "null-terminated ASCII string" => {
                self.current_array_mut().dtype = BinaryDataArrayType::ASCII;
            }

            // Array types
            "m/z array" => self.current_array_mut().name = ArrayType::MZArray,
            "intensity array" => self.current_array_mut().name = ArrayType::IntensityArray,
            "charge array" => self.current_array_mut().name = ArrayType::ChargeArray,
            "non-standard data array" => {
                self.current_array_mut().name =
                    ArrayType::NonStandardDataArray { name: param.value };
            }
            "mean ion mobility array"
            | "mean drift time array"
            | "mean inverse reduced ion mobility array" => {
                self.current_array_mut().name = ArrayType::MeanIonMobilityArray
            }
            "ion mobility array" | "drift time array" | "inverse reduced ion mobility array" => {
                self.current_array_mut().name = ArrayType::IonMobilityArray
            }
            "deconvoluted ion mobility array"
            | "deconvoluted drift time array"
            | "deconvoluted inverse reduced ion mobility array" => {
                self.current_array_mut().name = ArrayType::DeconvolutedIonMobilityArray
            }

            &_ => {
                self.current_array_mut().params.push(param);
            }
        }
    }

    fn fill_selected_ion(&mut self, param: Param) {
        match param.name.as_ref() {
            "selected ion m/z" => {
                self.selected_ion_mut().mz = param.coerce().expect("Failed to parse ion m/z");
            }
            "peak intensity" => {
                self.selected_ion_mut().intensity =
                    param.coerce().expect("Failed to parse peak intensity");
            }
            "charge state" => {
                self.selected_ion_mut().charge =
                    Some(param.coerce().expect("Failed to parse ion charge"));
            }
            &_ => {
                self.selected_ion_mut().params.push(param);
            }
        };
    }

    fn fill_isolation_window(&mut self, param: Param) {
        let window = self.isolation_window_mut();
        match param.name.as_ref() {
            "isolation window target m/z" => {
                window.target = param
                    .coerce()
                    .expect("Failed to parse isolation window target");
                window.flags = match window.flags {
                    IsolationWindowState::Unknown => IsolationWindowState::Complete,
                    IsolationWindowState::Explicit => IsolationWindowState::Complete,
                    IsolationWindowState::Offset => {
                        window.lower_bound = window.target - window.lower_bound;
                        window.upper_bound += window.target;
                        IsolationWindowState::Complete
                    }
                    IsolationWindowState::Complete => IsolationWindowState::Complete,
                };
            }
            "isolation window lower offset" => {
                let lower_bound: f64 = param
                    .coerce()
                    .expect("Failed to parse isolation window limit");
                match window.flags {
                    IsolationWindowState::Unknown => {
                        window.flags = IsolationWindowState::Offset;
                        window.lower_bound = lower_bound;
                    }
                    IsolationWindowState::Complete => {
                        window.lower_bound = window.target - lower_bound;
                    }
                    _ => {}
                }
            }
            "isolation window upper offset" => {
                let upper_bound: f64 = param
                    .coerce()
                    .expect("Failed to parse isolation window limit");
                match window.flags {
                    IsolationWindowState::Unknown => {
                        window.flags = IsolationWindowState::Offset;
                        window.upper_bound = upper_bound;
                    }
                    IsolationWindowState::Complete => {
                        window.upper_bound = window.target + upper_bound;
                    }
                    _ => {}
                }
            }
            "isolation window lower limit" => {
                let lower_bound: f64 = param
                    .coerce()
                    .expect("Failed to parse isolation window limit");
                match window.flags {
                    IsolationWindowState::Unknown => {
                        window.flags = IsolationWindowState::Explicit;
                        window.lower_bound = lower_bound;
                    }
                    _ => {}
                }
            }
            "isolation window upper limit" => {
                let upper_bound: f64 = param
                    .coerce()
                    .expect("Failed to parse isolation window limit");
                match window.flags {
                    IsolationWindowState::Unknown => {
                        window.flags = IsolationWindowState::Explicit;
                        window.upper_bound = upper_bound;
                    }
                    _ => {}
                }
            }
            &_ => {}
        }
    }

    fn fill_scan_window(&mut self, param: Param) {
        let window = self.scan_window_mut();
        match param.name.as_ref() {
            "scan window lower limit" => {
                window.lower_bound = param.coerce().expect("Failed to parse scan window limit");
            }
            "scan window upper limit" => {
                window.upper_bound = param.coerce().expect("Failed to parse scan window limit");
            }
            &_ => {}
        }
    }
}

pub type ParserResult = Result<MzMLParserState, MzMLParserError>;

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    SpectrumBuilding<C, D, MultiLayerSpectrum<C, D>> for MzMLSpectrumBuilder<C, D>
{
    fn isolation_window_mut(&mut self) -> &mut IsolationWindow {
        &mut self.precursor.isolation_window
    }

    fn scan_window_mut(&mut self) -> &mut ScanWindow {
        self.acquisition
            .scans
            .last_mut()
            .unwrap()
            .scan_windows
            .last_mut()
            .unwrap()
    }

    fn selected_ion_mut(&mut self) -> &mut SelectedIon {
        &mut self.precursor.ion
    }

    fn current_array_mut(&mut self) -> &mut DataArray {
        &mut self.current_array
    }

    fn into_spectrum(self, spectrum: &mut MultiLayerSpectrum<C, D>) {
        let description = &mut spectrum.description;

        description.id = self.scan_id;
        description.index = self.index;
        description.signal_continuity = self.signal_continuity;
        description.ms_level = self.ms_level;
        description.polarity = self.polarity;

        description.params = self.params;
        description.acquisition = self.acquisition;
        if self.has_precursor {
            description.precursor = Some(self.precursor);
        } else {
            description.precursor = None;
        }

        spectrum.arrays = Some(self.arrays);
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MzMLSpectrumBuilder<C, D> {
    pub fn new() -> MzMLSpectrumBuilder<C, D> {
        MzMLSpectrumBuilder {
            ..Default::default()
        }
    }

    pub fn _reset(&mut self) {
        self.params.clear();
        self.acquisition = Acquisition::default();
        self.arrays.clear();
        self.current_array.clear();
        self.scan_id.clear();

        self.precursor = Precursor::default();
        self.index = 0;
        self.has_precursor = false;
        self.signal_continuity = SignalContinuity::Unknown;
        self.polarity = ScanPolarity::Unknown;
    }

    pub fn _to_spectrum(&self, spectrum: &mut MultiLayerSpectrum<C, D>) {
        let description = &mut spectrum.description;

        description.id = self.scan_id.clone();
        description.index = self.index;
        description.signal_continuity = self.signal_continuity;
        description.ms_level = self.ms_level;
        description.polarity = self.polarity;

        description.params = self.params.clone();
        description.acquisition = self.acquisition.clone();
        if self.has_precursor {
            description.precursor = Some(self.precursor.clone());
        } else {
            description.precursor = None;
        }

        spectrum.arrays = Some(self.arrays.clone());
    }

    pub fn fill_spectrum(&mut self, param: Param) {
        match param.name.as_ref() {
            "ms level" => {
                self.ms_level = param.coerce().expect("Failed to parse ms level");
            }
            "positive scan" => {
                self.polarity = ScanPolarity::Positive;
            }
            "negative scan" => {
                self.polarity = ScanPolarity::Negative;
            }
            "profile spectrum" => {
                self.signal_continuity = SignalContinuity::Profile;
            }
            "centroid spectrum" => {
                self.signal_continuity = SignalContinuity::Centroid;
            }
            &_ => {
                self.params.push(param);
            }
        };
    }

    pub fn fill_param_into(&mut self, param: Param, state: MzMLParserState) {
        match state {
            MzMLParserState::Spectrum => {
                self.fill_spectrum(param);
            }
            MzMLParserState::ScanList => self.acquisition.params.push(param),
            MzMLParserState::Scan => {
                let event = self.acquisition.scans.last_mut().unwrap();
                match param.name.as_bytes() {
                    b"scan start time" => {
                        let value: f64 = param
                            .coerce()
                            .expect("Expected floating point number for scan time");
                        let value = if let Some(unit) = &param.unit_name {
                            match unit.as_bytes() {
                                b"minute" => value,
                                b"second" => 60.0 * value,
                                _ => {
                                    warn!("Could not infer unit for {:?}", param);
                                    value
                                }
                            }
                        } else {
                            value
                        };
                        event.start_time = value;
                    }
                    b"ion injection time" => {
                        event.injection_time = param
                            .coerce()
                            .expect("Expected floating point number for injection time");
                    }
                    _ => event.params.push(param),
                }
            }
            MzMLParserState::ScanWindowList => self
                .acquisition
                .scans
                .last_mut()
                .unwrap()
                .params
                .push(param),
            MzMLParserState::ScanWindow => {
                self.fill_scan_window(param);
            }
            MzMLParserState::IsolationWindow => {
                self.fill_isolation_window(param);
            }
            MzMLParserState::SelectedIon | MzMLParserState::SelectedIonList => {
                self.fill_selected_ion(param);
            }
            MzMLParserState::Activation => match param.name.as_ref() {
                "collision energy" | "activation energy" => {
                    self.precursor.activation.energy =
                        param.coerce().expect("Failed to parse collision energy");
                }
                &_ => {
                    self.precursor.activation.params.push(param);
                }
            },
            MzMLParserState::BinaryDataArrayList => {}
            MzMLParserState::BinaryDataArray => {
                self.fill_binary_data_array(param);
            }
            MzMLParserState::Precursor | MzMLParserState::PrecursorList => {
                self.precursor.params.push(param);
            }
            _ => {}
        };
    }

    pub fn start_element<B: io::BufRead>(
        &mut self,
        event: &BytesStart,
        state: MzMLParserState,
        reader: &Reader<B>,
    ) -> ParserResult {
        let elt_name = event.name();
        match elt_name {
            b"spectrum" => {
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => match attr.key {
                            b"id" => {
                                self.scan_id = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding id");
                            }
                            b"index" => {
                                self.index = (&String::from_utf8_lossy(&attr.value))
                                    .parse::<usize>()
                                    .expect("Failed to parse index");
                            }
                            _ => {}
                        },
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg, state));
                        }
                    }
                }
                return Ok(MzMLParserState::Spectrum);
            }
            b"spectrumList" => {
                return Ok(MzMLParserState::SpectrumList);
            }
            b"scanList" => {
                return Ok(MzMLParserState::ScanList);
            }
            b"scan" => {
                let mut scan_event = ScanEvent::default();
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key == b"instrumentConfigurationRef" {
                                scan_event.instrument_configuration_id = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding id");
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg, state));
                        }
                    }
                }
                self.acquisition.scans.push(scan_event);
                return Ok(MzMLParserState::Scan);
            }
            b"scanWindow" => {
                let window = ScanWindow::default();
                self.acquisition
                    .scans
                    .last_mut()
                    .expect("Scan window without scan")
                    .scan_windows
                    .push(window);
                return Ok(MzMLParserState::ScanWindow);
            }
            b"scanWindowList" => {
                return Ok(MzMLParserState::ScanWindowList);
            }
            b"precursorList" => {
                return Ok(MzMLParserState::PrecursorList);
            }
            b"precursor" => {
                self.has_precursor = true;
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key == b"spectrumRef" {
                                self.precursor.precursor_id = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding id");
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg, state));
                        }
                    }
                }
                return Ok(MzMLParserState::Precursor);
            }
            b"isolationWindow" => {
                return Ok(MzMLParserState::IsolationWindow);
            }
            b"selectedIonList" => {
                return Ok(MzMLParserState::SelectedIonList);
            }
            b"selectedIon" => {
                return Ok(MzMLParserState::SelectedIon);
            }
            b"activation" => {
                return Ok(MzMLParserState::Activation);
            }
            b"binaryDataArrayList" => {
                return Ok(MzMLParserState::BinaryDataArrayList);
            }
            b"binaryDataArray" => {
                return Ok(MzMLParserState::BinaryDataArray);
            }
            b"binary" => {
                return Ok(MzMLParserState::Binary);
            }
            _ => {}
        };
        Ok(state)
    }

    pub fn empty_element<B: io::BufRead>(
        &mut self,
        event: &BytesStart,
        state: MzMLParserState,
        reader: &Reader<B>,
    ) -> ParserResult {
        let elt_name = event.name();
        match elt_name {
            b"cvParam" | b"userParam" => match self.handle_param(event, reader, state) {
                Ok(param) => {
                    self.fill_param_into(param, state);
                    return Ok(state);
                }
                Err(err) => return Err(err),
            },
            &_ => {}
        }
        Ok(state)
    }

    pub fn end_element(&mut self, event: &BytesEnd, state: MzMLParserState) -> ParserResult {
        let elt_name = event.name();
        match elt_name {
            b"spectrum" => return Ok(MzMLParserState::SpectrumDone),
            b"scanList" => return Ok(MzMLParserState::Spectrum),
            b"scan" => return Ok(MzMLParserState::ScanList),
            b"scanWindow" => return Ok(MzMLParserState::ScanWindowList),
            b"scanWindowList" => return Ok(MzMLParserState::Scan),
            b"precursorList" => return Ok(MzMLParserState::Spectrum),
            b"precursor" => return Ok(MzMLParserState::PrecursorList),
            b"isolationWindow" => return Ok(MzMLParserState::Precursor),
            b"selectedIonList" => return Ok(MzMLParserState::Precursor),
            b"selectedIon" => return Ok(MzMLParserState::SelectedIonList),
            b"activation" => return Ok(MzMLParserState::Precursor),
            b"binaryDataArrayList" => {
                return Ok(MzMLParserState::Spectrum);
            }
            b"binaryDataArray" => {
                let mut array = self.current_array.clone();
                array
                    .decode_and_store()
                    .expect("Error during decoding and storing of array data");
                self.arrays.add(array);
                self.current_array.clear();
                return Ok(MzMLParserState::BinaryDataArrayList);
            }
            b"binary" => return Ok(MzMLParserState::BinaryDataArray),
            b"spectrumList" => return Ok(MzMLParserState::SpectrumListDone),
            _ => {}
        };
        Ok(state)
    }

    pub fn text(&mut self, event: &BytesText, state: MzMLParserState) -> ParserResult {
        if state == MzMLParserState::Binary {
            let bin = event
                .unescaped()
                .expect("Failed to unescape binary data array content");
            self.current_array.data = Bytes::from(&*bin);
        }
        Ok(state)
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> Into<CentroidSpectrumType<C>>
    for MzMLSpectrumBuilder<C, D>
{
    fn into(self) -> CentroidSpectrumType<C> {
        let mut spec = MultiLayerSpectrum::<C, D>::default();
        self.into_spectrum(&mut spec);
        spec.try_into().unwrap()
    }
}

impl Into<Spectrum> for MzMLSpectrumBuilder {
    fn into(self) -> Spectrum {
        let mut spec = Spectrum::default();
        self.into_spectrum(&mut spec);
        spec
    }
}

impl Into<RawSpectrum> for MzMLSpectrumBuilder {
    fn into(self) -> RawSpectrum {
        let mut spec = Spectrum::default();
        self.into_spectrum(&mut spec);
        spec.into()
    }
}

#[derive(Debug, Default, Clone)]
pub struct FileMetadataBuilder {
    pub file_description: FileDescription,
    pub instrument_configurations: Vec<InstrumentConfiguration>,
    pub softwares: Vec<Software>,
    pub data_processings: Vec<DataProcessing>,
    pub reference_param_groups: HashMap<String, Vec<Param>>,
    pub last_group: String,
}

impl XMLParseBase for FileMetadataBuilder {}
impl CVParamParse for FileMetadataBuilder {}

impl FileMetadataBuilder {
    pub fn start_element<B: io::BufRead>(
        &mut self,
        event: &BytesStart,
        state: MzMLParserState,
        reader: &Reader<B>,
    ) -> ParserResult {
        let elt_name = event.name();
        match elt_name {
            b"fileDescription" => return Ok(MzMLParserState::FileDescription),
            b"fileContent" => return Ok(MzMLParserState::FileContents),
            b"sourceFileList" => return Ok(MzMLParserState::SourceFileList),
            b"sourceFile" => {
                let mut source_file = SourceFile::default();
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key == b"id" {
                                source_file.id = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding id");
                            } else if attr.key == b"name" {
                                source_file.name = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding name");
                            } else if attr.key == b"location" {
                                source_file.location = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding location");
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg, state));
                        }
                    }
                }
                self.file_description.source_files.push(source_file);
                return Ok(MzMLParserState::SourceFile);
            }
            b"softwareList" => return Ok(MzMLParserState::SoftwareList),
            b"software" => {
                let mut software = Software::default();
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key == b"id" {
                                software.id = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding id");
                            } else if attr.key == b"version" {
                                software.version = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding version");
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg, state));
                        }
                    }
                }
                self.softwares.push(software);
                return Ok(MzMLParserState::Software);
            }
            b"referenceableParamGroupList" => {
                return Ok(MzMLParserState::ReferenceParamGroupList);
            }
            b"referenceableParamGroup" => {
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key == b"id" {
                                let key = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding id");
                                self.reference_param_groups.entry(key.clone()).or_default();
                                self.last_group = key;
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg, state));
                        }
                    }
                }
                return Ok(MzMLParserState::ReferenceParamGroup);
            }
            b"instrumentConfigurationList" => {
                return Ok(MzMLParserState::InstrumentConfigurationList)
            }
            b"instrumentConfiguration" => {
                let mut ic = InstrumentConfiguration::default();
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key == b"id" {
                                ic.id = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding id");
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg, state));
                        }
                    }
                }
                self.instrument_configurations.push(ic);
                return Ok(MzMLParserState::InstrumentConfiguration);
            }
            b"componentList" => return Ok(MzMLParserState::ComponentList),
            b"source" => {
                let mut source = Component::default();
                source.component_type = ComponentType::IonSource;
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key == b"order" {
                                source.order = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding order")
                                    .parse()
                                    .expect("Failed to parse integer from `order`");
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg, state));
                        }
                    }
                }
                self.instrument_configurations
                    .last_mut()
                    .unwrap()
                    .components
                    .push(source);
                return Ok(MzMLParserState::Source);
            }
            b"analyzer" => {
                let mut analyzer = Component::default();
                analyzer.component_type = ComponentType::Analyzer;
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key == b"order" {
                                analyzer.order = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding order")
                                    .parse()
                                    .expect("Failed to parse integer from `order`");
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg, state));
                        }
                    }
                }
                self.instrument_configurations
                    .last_mut()
                    .unwrap()
                    .components
                    .push(analyzer);
                return Ok(MzMLParserState::Analyzer);
            }
            b"detector" => {
                let mut detector = Component::default();
                detector.component_type = ComponentType::Detector;
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key == b"order" {
                                detector.order = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding order")
                                    .parse()
                                    .expect("Failed to parse integer from `order`");
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg, state));
                        }
                    }
                }
                self.instrument_configurations
                    .last_mut()
                    .unwrap()
                    .components
                    .push(detector);
                return Ok(MzMLParserState::Detector);
            }
            b"dataProcessingList" => return Ok(MzMLParserState::DataProcessingList),
            b"dataProcessing" => {
                let mut dp = DataProcessing::default();
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key == b"id" {
                                dp.id = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding id");
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg, state));
                        }
                    }
                }
                self.data_processings.push(dp);
                return Ok(MzMLParserState::DataProcessing);
            }
            b"processingMethod" => {
                let mut method = ProcessingMethod::default();
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key == b"order" {
                                method.order = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding order")
                                    .parse()
                                    .expect("Failed to parse order");
                            } else if attr.key == b"softwareRef" {
                                method.software_reference = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding softwareRef");
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg, state));
                        }
                    }
                }
                self.data_processings
                    .last_mut()
                    .unwrap()
                    .methods
                    .push(method);
                return Ok(MzMLParserState::ProcessingMethod);
            }
            b"run" => return Ok(MzMLParserState::Run),
            _ => {}
        }

        Ok(state)
    }

    pub fn fill_param_into(&mut self, param: Param, state: MzMLParserState) {
        match state {
            MzMLParserState::SourceFile => {
                let sf = self.file_description.source_files.last_mut().unwrap();
                sf.add_param(param)
            }
            MzMLParserState::FileContents => {
                self.file_description.add_param(param);
            }
            MzMLParserState::InstrumentConfiguration => {
                self.instrument_configurations
                    .last_mut()
                    .unwrap()
                    .add_param(param);
            }
            MzMLParserState::Software => self.softwares.last_mut().unwrap().add_param(param),
            MzMLParserState::ProcessingMethod => self
                .data_processings
                .last_mut()
                .unwrap()
                .methods
                .last_mut()
                .unwrap()
                .add_param(param),
            MzMLParserState::Detector | MzMLParserState::Analyzer | MzMLParserState::Source => self
                .instrument_configurations
                .last_mut()
                .unwrap()
                .components
                .last_mut()
                .unwrap()
                .add_param(param),
            MzMLParserState::ReferenceParamGroup => {
                self.reference_param_groups
                    .get_mut(&self.last_group)
                    .unwrap()
                    .push(param);
            }
            _ => {}
        }
    }

    pub fn empty_element<B: io::BufRead>(
        &mut self,
        event: &BytesStart,
        state: MzMLParserState,
        reader: &Reader<B>,
    ) -> ParserResult {
        let elt_name = event.name();
        match elt_name {
            b"cvParam" | b"userParam" => match self.handle_param(event, reader, state) {
                Ok(param) => {
                    self.fill_param_into(param, state);
                    return Ok(state);
                }
                Err(err) => return Err(err),
            },
            b"softwareRef" => match state {
                MzMLParserState::InstrumentConfiguration => {
                    let ic = self.instrument_configurations.last_mut().unwrap();
                    for attr_parsed in event.attributes() {
                        match attr_parsed {
                            Ok(attr) => {
                                if attr.key == b"ref" {
                                    ic.software_reference = attr
                                        .unescape_and_decode_value(reader)
                                        .expect("Error decoding software reference");
                                }
                            }
                            Err(msg) => {
                                return Err(self.handle_xml_error(msg, state));
                            }
                        }
                    }
                }
                _ => {}
            },
            b"referenceableParamGroupRef" => {
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key == b"ref" {
                                let group_id = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding reference group");

                                let param_group = match self.reference_param_groups.get(&group_id) {
                                    Some(params) => params.clone(),
                                    None => {
                                        panic!("Encountered a referenceableParamGroupRef without a group definition")
                                    }
                                };

                                for param in param_group {
                                    self.fill_param_into(param, state)
                                }
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg, state));
                        }
                    }
                }
            }
            &_ => {}
        }
        Ok(state)
    }

    pub fn end_element(&mut self, event: &BytesEnd, state: MzMLParserState) -> ParserResult {
        let elt_name = event.name();
        match elt_name {
            b"fileDescription" => return Ok(MzMLParserState::FileDescription),
            b"fileContent" => return Ok(MzMLParserState::FileDescription),
            b"sourceFile" => return Ok(MzMLParserState::SourceFileList),
            b"softwareList" => return Ok(MzMLParserState::SoftwareList),
            b"software" => return Ok(MzMLParserState::SoftwareList),
            b"referenceableParamGroupList" => {
                return Ok(MzMLParserState::ReferenceParamGroupList);
            }
            b"instrumentConfigurationList" => {
                return Ok(MzMLParserState::InstrumentConfigurationList)
            }
            b"instrumentConfiguration" => return Ok(MzMLParserState::InstrumentConfigurationList),
            b"componentList" => return Ok(MzMLParserState::InstrumentConfiguration),
            b"source" => return Ok(MzMLParserState::ComponentList),
            b"analyzer" => return Ok(MzMLParserState::ComponentList),
            b"detector" => return Ok(MzMLParserState::ComponentList),
            b"dataProcessingList" => return Ok(MzMLParserState::DataProcessingList),
            b"dataProcessing" => return Ok(MzMLParserState::DataProcessingList),
            b"processingMethod" => return Ok(MzMLParserState::DataProcessing),
            b"run" => {
                // TODO
            }
            _ => {}
        }
        Ok(state)
    }

    pub fn text(&mut self, _event: &BytesText, state: MzMLParserState) -> ParserResult {
        Ok(state)
    }
}

/// An mzML parser that supports iteration and random access. The parser produces
/// [`Spectrum`] instances, which may be converted to [`RawSpectrum`](crate::spectrum::spectrum::RawSpectrum)
/// or [`CentroidSpectrum`](crate::spectrum::CentroidSpectrum) as is appropriate to the data.
///
/// When the readable stream the parser is wrapped around supports [`io::Seek`],
/// additional random access operations are available.
pub struct MzMLReaderType<
    R: Read,
    C: CentroidPeakAdapting = CentroidPeak,
    D: DeconvolutedPeakAdapting = DeconvolutedPeak,
> {
    /// The state the parser was in last.
    pub state: MzMLParserState,
    pub handle: BufReader<R>,
    /// A place to store the last error the parser encountered
    pub error: MzMLParserError,
    pub index: OffsetIndex,
    /// The description of the file's contents and the previous data files that were
    /// consumed to produce it.
    pub file_description: FileDescription,
    /// A mapping of different instrument configurations (source, analyzer, detector) components
    /// by ID string.
    pub instrument_configurations: HashMap<String, InstrumentConfiguration>,
    /// The different software components that were involved in the processing and creation of this
    /// file.
    pub softwares: Vec<Software>,
    /// The data processing and signal transformation operations performed on the raw data in previous
    /// source files to produce this file's contents.
    pub data_processings: Vec<DataProcessing>,
    /// A cache of repeated paramters
    pub reference_param_groups: HashMap<String, Vec<Param>>,
    buffer: Bytes,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
}

impl<R: Read, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MzMLReaderType<R, C, D> {
    pub fn new(file: R) -> MzMLReaderType<R, C, D> {
        let handle = BufReader::with_capacity(BUFFER_SIZE, file);
        let mut inst = MzMLReaderType {
            handle,
            state: MzMLParserState::Start,
            error: MzMLParserError::default(),
            buffer: Bytes::new(),
            index: OffsetIndex::new("spectrum".to_owned()),

            file_description: FileDescription::default(),
            instrument_configurations: HashMap::new(),
            softwares: Vec::new(),
            data_processings: Vec::new(),
            reference_param_groups: HashMap::new(),

            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
        };
        match inst.parse_metadata() {
            Ok(()) => {}
            Err(err) => {
                warn!(
                    "Encountered error {:?} while parsing mzML file metadata",
                    err
                );
            }
        }
        inst
    }

    fn parse_metadata(&mut self) -> Result<(), MzMLParserError> {
        let mut reader = Reader::from_reader(&mut self.handle);
        reader.trim_text(true);
        let mut accumulator = FileMetadataBuilder::default();
        loop {
            match reader.read_event(&mut self.buffer) {
                Ok(Event::Start(ref e)) => {
                    match accumulator.start_element(e, self.state, &reader) {
                        Ok(state) => {
                            self.state = state;
                            match &self.state {
                                MzMLParserState::Run
                                | MzMLParserState::SpectrumList
                                | MzMLParserState::Spectrum => break,
                                _ => {}
                            }
                        }
                        Err(message) => {
                            self.state = MzMLParserState::ParserError;
                            self.error = message;
                        }
                    };
                }
                Ok(Event::End(ref e)) => {
                    match accumulator.end_element(e, self.state) {
                        Ok(state) => {
                            self.state = state;
                        }
                        Err(message) => {
                            self.state = MzMLParserState::ParserError;
                            self.error = message;
                        }
                    };
                }
                Ok(Event::Text(ref e)) => {
                    match accumulator.text(e, self.state) {
                        Ok(state) => {
                            self.state = state;
                        }
                        Err(message) => {
                            self.state = MzMLParserState::ParserError;
                            self.error = message;
                        }
                    };
                }
                Ok(Event::Empty(ref e)) => {
                    match accumulator.empty_element(e, self.state, &reader) {
                        Ok(state) => {
                            self.state = state;
                        }
                        Err(message) => {
                            self.state = MzMLParserState::ParserError;
                            self.error = message;
                        }
                    }
                }
                Ok(Event::Eof) => {
                    break;
                }
                Err(err) => match &err {
                    XMLError::EndEventMismatch {
                        expected,
                        found: _found,
                    } => {
                        if expected.is_empty() && self.state == MzMLParserState::Resume {
                            continue;
                        } else {
                            self.error = MzMLParserError::IncompleteElementError(
                                String::from_utf8_lossy(&self.buffer).to_owned().to_string(),
                                self.state,
                            );
                            self.state = MzMLParserState::ParserError;
                        }
                    }
                    _ => {
                        self.error = MzMLParserError::IncompleteElementError(
                            String::from_utf8_lossy(&self.buffer).to_owned().to_string(),
                            self.state,
                        );
                        self.state = MzMLParserState::ParserError;
                    }
                },
                _ => {}
            };
            self.buffer.clear();
            match self.state {
                MzMLParserState::Run | MzMLParserState::ParserError => {
                    break;
                }
                _ => {}
            };
        }
        self.file_description = accumulator.file_description;
        self.instrument_configurations = accumulator
            .instrument_configurations
            .into_iter()
            .map(|ic| (ic.id.clone(), ic))
            .collect();
        self.softwares = accumulator.softwares;
        self.data_processings = accumulator.data_processings;
        self.reference_param_groups = accumulator.reference_param_groups;

        match self.state {
            MzMLParserState::SpectrumDone => Ok(()),
            MzMLParserState::ParserError => Err(self.error.clone()),
            _ => Err(MzMLParserError::IncompleteSpectrum),
        }
    }

    fn _parse_into(
        &mut self,
        accumulator: &mut MzMLSpectrumBuilder<C, D>,
    ) -> Result<usize, MzMLParserError> {
        let mut reader = Reader::from_reader(&mut self.handle);
        reader.trim_text(true);
        let mut offset: usize = 0;
        loop {
            match reader.read_event(&mut self.buffer) {
                Ok(Event::Start(ref e)) => {
                    match accumulator.start_element(e, self.state, &reader) {
                        Ok(state) => {
                            self.state = state;
                        }
                        Err(message) => {
                            self.state = MzMLParserState::ParserError;
                            self.error = message;
                        }
                    };
                }
                Ok(Event::End(ref e)) => {
                    match accumulator.end_element(e, self.state) {
                        Ok(state) => {
                            self.state = state;
                        }
                        Err(message) => {
                            self.state = MzMLParserState::ParserError;
                            self.error = message;
                        }
                    };
                }
                Ok(Event::Text(ref e)) => {
                    match accumulator.text(e, self.state) {
                        Ok(state) => {
                            self.state = state;
                        }
                        Err(message) => {
                            self.state = MzMLParserState::ParserError;
                            self.error = message;
                        }
                    };
                }
                Ok(Event::Empty(ref e)) => {
                    match accumulator.empty_element(e, self.state, &reader) {
                        Ok(state) => {
                            self.state = state;
                        }
                        Err(message) => {
                            self.state = MzMLParserState::ParserError;
                            self.error = message;
                        }
                    }
                }
                Ok(Event::Eof) => {
                    break;
                }
                Err(err) => match &err {
                    XMLError::EndEventMismatch {
                        expected,
                        found: _found,
                    } => {
                        if expected.is_empty() && self.state == MzMLParserState::Resume {
                            continue;
                        } else {
                            self.error = MzMLParserError::IncompleteElementError(
                                String::from_utf8_lossy(&self.buffer).to_owned().to_string(),
                                self.state,
                            );
                            self.state = MzMLParserState::ParserError;
                        }
                    }
                    _ => {
                        self.error = MzMLParserError::IncompleteElementError(
                            String::from_utf8_lossy(&self.buffer).to_owned().to_string(),
                            self.state,
                        );
                        self.state = MzMLParserState::ParserError;
                    }
                },
                _ => {}
            };
            offset += self.buffer.len();
            self.buffer.clear();
            match self.state {
                MzMLParserState::SpectrumDone | MzMLParserState::ParserError => {
                    break;
                }
                _ => {}
            };
        }
        match self.state {
            MzMLParserState::SpectrumDone => Ok(offset),
            MzMLParserState::ParserError => Err(self.error.clone()),
            _ => Err(MzMLParserError::IncompleteSpectrum),
        }
    }

    /// Populate a new [`Spectrum`] in-place on the next available spectrum data.
    /// This allocates memory to build the spectrum's attributes but then moves it
    /// into `spectrum` rather than copying it.
    pub fn read_into(
        &mut self,
        spectrum: &mut MultiLayerSpectrum<C, D>,
    ) -> Result<usize, MzMLParserError> {
        let mut accumulator = MzMLSpectrumBuilder::<C, D>::new();
        if self.state == MzMLParserState::SpectrumDone {
            self.state = MzMLParserState::Resume;
        }
        match self._parse_into(&mut accumulator) {
            Ok(sz) => {
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
            Ok(_sz) => Some(spectrum.into()),
            Err(_err) => None,
        }
    }
}

/// [`MzMLReaderType`] instances are [`Iterator`]s over [`Spectrum`]
impl<R: io::Read, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> Iterator
    for MzMLReaderType<R, C, D>
{
    type Item = MultiLayerSpectrum<C, D>;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_next()
    }
}

/// They can also be used to fetch specific spectra by ID, index, or start
/// time when the underlying file stream supports [`io::Seek`].
impl<R: io::Read + io::Seek, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    ScanSource<C, D, MultiLayerSpectrum<C, D>> for MzMLReaderType<R, C, D>
{
    /// Retrieve a spectrum by it's native ID
    fn get_spectrum_by_id(&mut self, id: &str) -> Option<MultiLayerSpectrum<C, D>> {
        let offset_ref = self.index.get(id);
        let offset = offset_ref.expect("Failed to retrieve offset");
        let start = self
            .handle
            .stream_position()
            .expect("Failed to save checkpoint");
        self.seek(SeekFrom::Start(offset))
            .expect("Failed to move seek to offset");
        let result = self.read_next();
        self.seek(SeekFrom::Start(start))
            .expect("Failed to restore offset");
        result
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
        result
    }

    /// Return the data stream to the beginning
    fn reset(&mut self) -> &Self {
        self.seek(SeekFrom::Start(0))
            .expect("Failed to reset file stream");
        self
    }

    fn get_index(&self) -> &OffsetIndex {
        if !self.index.init {
            warn!("Attempting to use an uninitialized offset index on MzMLReaderType")
        }
        &self.index
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.index = index
    }
}

/// The iterator can also be updated to move to a different location in the
/// stream efficiently.
impl<R: SeekRead, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    RandomAccessScanIterator<C, D, MultiLayerSpectrum<C, D>> for MzMLReaderType<R, C, D>
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

impl<R: SeekRead, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MzMLReaderType<R, C, D> {
    /// Construct a new MzMLReaderType and build an offset index
    /// using [`Self::build_index`]
    pub fn new_indexed(file: R) -> MzMLReaderType<R, C, D> {
        let mut reader = Self::new(file);
        reader.build_index();
        reader
    }

    pub fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.handle.seek(pos)
    }

    /// Builds an offset index to each `<spectrum>` XML element
    /// by doing a fast pre-scan of the XML file.
    pub fn build_index(&mut self) -> u64 {
        let start = self
            .handle
            .stream_position()
            .expect("Failed to save restore location");
        self.seek(SeekFrom::Start(0))
            .expect("Failed to reset stream to beginning");
        let mut reader = Reader::from_reader(&mut self.handle);
        reader.trim_text(true);
        loop {
            match reader.read_event(&mut self.buffer) {
                Ok(Event::Start(ref e)) => {
                    let element_name = e.name();
                    if element_name == b"spectrum" {
                        // Hit a spectrum, extract ID and save current offset

                        for attr_parsed in e.attributes() {
                            match attr_parsed {
                                Ok(attr) => {
                                    match attr.key {
                                        b"id" => {
                                            let scan_id = attr
                                                .unescape_and_decode_value(&reader)
                                                .expect("Error decoding id");
                                            // This count is off by 2 because somehow the < and > bytes are removed?
                                            self.index.insert(
                                                scan_id,
                                                (reader.buffer_position() - e.len() - 2) as u64,
                                            );
                                            break;
                                        }
                                        &_ => {}
                                    };
                                }
                                Err(_msg) => {}
                            }
                        }
                    }
                }
                Ok(Event::End(ref e)) => {
                    let element_name = e.name();
                    if element_name == b"spectrumList" {
                        break;
                    }
                }
                Ok(Event::Eof) => {
                    break;
                }
                _ => {}
            };
            self.buffer.clear();
        }
        let offset = reader.buffer_position() as u64;
        self.handle
            .seek(SeekFrom::Start(start))
            .expect("Failed to restore location");
        self.index.init = true;
        if self.index.len() == 0 {
            warn!("An index was built but no entries were found")
        }
        offset
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    MZFileReader<C, D, MultiLayerSpectrum<C, D>> for MzMLReaderType<fs::File, C, D>
{
    fn open_file(source: fs::File) -> Self {
        Self::new(source)
    }

    fn construct_index_from_stream(&mut self) -> u64 {
        self.build_index()
    }
}

impl<R: Read, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MSDataFileMetadata
    for MzMLReaderType<R, C, D>
{
    crate::impl_metadata_trait!();
}

pub type MzMLReader<R> = MzMLReaderType<R, CentroidPeak, DeconvolutedPeak>;

macro_rules! bstart {
    ($e:tt) => {
        BytesStart::owned($e.as_bytes().to_vec(), $e.len())
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
            .handle
            .write_event(Event::Start($target.to_borrowed()))?;
    };
}

macro_rules! end_event {
    ($writer:ident, $target:ident) => {
        $writer.handle.write_event(Event::End($target.to_end()))?;
    };
}

pub type XMLResult = Result<(), XMLError>;

struct InnerXMLWriter<W: io::Write> {
    pub handle: Writer<BufWriter<MD5HashingStream<W>>>,
}

impl<W: Write> Debug for InnerXMLWriter<W> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("InnerXMLWriter")
            .field("handle", &"...")
            .finish()
    }
}

impl<W: io::Write> InnerXMLWriter<W> {
    const INDENT_SIZE: u64 = 2;

    pub fn new(file: W) -> InnerXMLWriter<W> {
        let handle = BufWriter::with_capacity(BUFFER_SIZE, MD5HashingStream::new(file));
        Self {
            handle: Writer::new_with_indent(handle, b' ', 2),
        }
    }

    pub fn digest(&mut self) -> String {
        let digest = self.handle.inner().get_ref().compute();
        format!("{:x}", digest)
    }

    pub fn flush(&mut self) -> io::Result<()> {
        self.handle.inner().flush()
    }

    pub fn write_param(&mut self, param: &Param) -> XMLResult {
        let mut elt = if param.accession.is_empty() {
            bstart!("userParam")
        } else {
            let mut elt = bstart!("cvParam");
            attrib!("accession", param.accession, elt);
            if let Some(cv_ref) = &param.controlled_vocabulary {
                attrib!("cvRef", cv_ref, elt);
            }
            elt
        };

        attrib!("name", param.name, elt);
        if !param.value.is_empty() {
            attrib!("value", param.value, elt);
        }

        if param.unit_accession.is_some() || param.unit_name.is_some() {
            if let Some(unit_acc) = &param.unit_accession {
                let mut split = unit_acc.split(":");
                if let Some(prefix) = split.next() {
                    attrib!("unitCvRef", prefix, elt);
                } else {
                    attrib!("unitCvRef", "UO", elt);
                }
                attrib!("unitAccession", unit_acc, elt);
            } else {
                attrib!("unitCvRef", "UO", elt);
            }
            if let Some(unit_name) = &param.unit_name {
                attrib!("unitName", unit_name, elt);
            }
        }
        self.handle.write_event(Event::Empty(elt))
    }

    pub fn write_param_list<'a, T: Iterator<Item = &'a Param>>(&mut self, params: T) -> XMLResult {
        for param in params {
            self.write_param(param)?
        }
        Ok(())
    }

    pub fn write_event(&mut self, event: Event) -> XMLResult {
        self.handle.write_event(event)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord)]
pub enum MzMLWriterState {
    Start,
    DocumentOpen,
    Header,
    Run,
    SpectrumList,
    SpectrumListClosed,
    ChromatogramList,
    ChromatogramListClosed,
    RunClosed,
    MzMLClosed,
    IndexList,
    IndexListClosed,
    End,
}

#[derive(Debug)]
pub struct MzMLWriterType<
    W: Write + Seek,
    C: CentroidPeakAdapting + 'static = CentroidPeak,
    D: DeconvolutedPeakAdapting + 'static = DeconvolutedPeak,
> {
    pub offset: usize,

    pub spectrum_count: u64,
    pub spectrum_counter: u64,

    pub chromatogram_count: u64,
    pub chromatogram_counter: u64,

    pub data_array_compression: BinaryCompressionType,

    pub file_description: FileDescription,
    pub softwares: Vec<Software>,
    pub data_processings: Vec<DataProcessing>,
    pub instrument_configurations: HashMap<String, InstrumentConfiguration>,

    pub state: MzMLWriterState,
    pub offset_index: OffsetIndex,

    handle: InnerXMLWriter<W>,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
    ms_cv: ControlledVocabulary,
}

impl<
        'a,
        W: Write + Seek,
        C: CentroidPeakAdapting + 'static,
        D: DeconvolutedPeakAdapting + 'static,
    > ScanWriter<'a, W, C, D> for MzMLWriterType<W, C, D>
where
    &'a PeakSetVec<C, MZ>: Into<BinaryArrayMap>,
    &'a PeakSetVec<D, Mass>: Into<BinaryArrayMap>,
{
    fn write(&mut self, spectrum: &'a MultiLayerSpectrum<C, D>) -> io::Result<usize> {
        match self.write_spectrum(spectrum) {
            Ok(()) => {
                let pos = self.stream_position()?;
                Ok(pos as usize)
            }
            Err(err) => {
                let msg = err.to_string();
                Err(io::Error::new(io::ErrorKind::InvalidData, msg))
            }
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        self.handle.flush()
    }
}

impl<W: Write + Seek, C: CentroidPeakAdapting + 'static, D: DeconvolutedPeakAdapting + 'static>
    MSDataFileMetadata for MzMLWriterType<W, C, D>
{
    crate::impl_metadata_trait!();
}

impl<'a, W: Write + Seek, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>
    MzMLWriterType<W, C, D>
where
    &'a PeakSetVec<C, MZ>: Into<BinaryArrayMap>,
    &'a PeakSetVec<D, Mass>: Into<BinaryArrayMap>,
{
    const PSIMS_VERSION: &'static str = "4.1.57";
    const UNIT_VERSION: &'static str = "releases/2020-03-10";

    pub fn new(file: W) -> MzMLWriterType<W, C, D> {
        let handle = InnerXMLWriter::new(file);
        let inst = MzMLWriterType {
            handle,
            file_description: FileDescription::default(),
            instrument_configurations: HashMap::new(),
            softwares: Vec::new(),
            data_processings: Vec::new(),
            offset: 0,
            offset_index: OffsetIndex::new("spectrum".into()),
            state: MzMLWriterState::Start,
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
            spectrum_count: 0,
            spectrum_counter: 0,
            chromatogram_count: 0,
            chromatogram_counter: 0,
            ms_cv: ControlledVocabulary::new("MS".to_string()),
            data_array_compression: BinaryCompressionType::Zlib,
        };
        inst
    }

    fn stream_position(&mut self) -> io::Result<u64> {
        self.handle.handle.inner().stream_position()
    }

    fn write_cv_list(&mut self) -> XMLResult {
        let mut cv_list = BytesStart::owned(b"cvList".to_vec(), 6);
        cv_list.push_attribute(("count", "2"));
        self.handle.write_event(Event::Start(cv_list))?;

        let mut cv = BytesStart::owned(b"cv".to_vec(), 2);
        cv.push_attribute(("id", "MS"));
        cv.push_attribute(("fullName", "PSI-MS"));
        cv.push_attribute((
            "URI",
            "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo",
        ));
        cv.push_attribute(("version", Self::PSIMS_VERSION));
        self.handle.write_event(Event::Empty(cv))?;

        let mut cv = BytesStart::owned(b"cv".to_vec(), 2);
        cv.push_attribute(("id", "UO"));
        cv.push_attribute(("fullName", "UNIT-ONTOLOGY"));
        cv.push_attribute(("URI", "http://ontologies.berkeleybop.org/uo.obo"));
        cv.push_attribute(("version", Self::UNIT_VERSION));
        self.handle.write_event(Event::Empty(cv))?;

        self.handle
            .write_event(Event::End(BytesEnd::owned(b"cvList".to_vec())))?;
        Ok(())
    }

    fn start_document(&mut self) -> XMLResult {
        self.handle
            .write_event(Event::Decl(BytesDecl::new(b"1.0", Some(b"utf-8"), None)))?;
        let mut indexed = BytesStart::owned(b"indexedmzML".to_vec(), 11);
        indexed.push_attribute(("xmlns", "http://psi.hupo.org/ms/mzml"));
        indexed.push_attribute(("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance"));
        indexed.push_attribute((
            "xsi:schemaLocation",
            "http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd",
        ));
        self.handle.write_event(Event::Start(indexed))?;

        let mut mzml = BytesStart::owned(b"mzML".to_vec(), 4);
        mzml.push_attribute(("xmlns", "http://psi.hupo.org/ms/mzml"));
        mzml.push_attribute(("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance"));
        mzml.push_attribute((
            "xsi:schemaLocation",
            "http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd",
        ));
        mzml.push_attribute(("version", "1.1.0"));
        self.handle.write_event(Event::Start(mzml))?;

        self.state = MzMLWriterState::DocumentOpen;
        Ok(())
    }

    fn writer_header(&mut self) -> XMLResult {
        if self.state < MzMLWriterState::DocumentOpen {
            self.start_document()?;
        } else {
            panic!(
                "Cannot start writing the header of mzML, currently in state {:?} which happens after.", self.state)
        }
        self.write_cv_list()?;
        self.write_file_description()?;
        self.write_software_list()?;
        self.write_instrument_configuration()?;
        self.write_data_processing()?;

        self.state = MzMLWriterState::Header;
        Ok(())
    }

    fn write_file_description(&mut self) -> XMLResult {
        let fd = bstart!("fileDescription");
        start_event!(self, fd);

        let fc_tag = bstart!("fileContents");
        start_event!(self, fc_tag);
        for param in self.file_description.params() {
            self.handle.write_param(param)?
        }
        end_event!(self, fc_tag);

        let mut outer = bstart!("sourceFileList");
        let count = self.file_description.source_files.len().to_string();
        attrib!("count", count, outer);
        self.handle.write_event(Event::Start(outer.to_borrowed()))?;
        for sf in self.file_description.source_files.iter() {
            let mut tag = bstart!("sourceFile");
            attrib!("id", sf.id, tag);
            attrib!("name", sf.name, tag);
            attrib!("location", sf.location, tag);
            self.handle.write_event(Event::Start(tag.to_borrowed()))?;
            for param in sf.params() {
                self.handle.write_param(param)?
            }
            self.handle.write_event(Event::End(tag.to_end()))?;
        }
        self.handle.write_event(Event::End(outer.to_end()))?;

        end_event!(self, fd);
        Ok(())
    }

    fn write_software_list(&mut self) -> XMLResult {
        let mut outer = bstart!("softwareList");
        let count = self.softwares.len().to_string();
        attrib!("count", count, outer);
        self.handle.write_event(Event::Start(outer.to_borrowed()))?;
        for soft in self.softwares.iter() {
            let mut tag = bstart!("software");
            attrib!("id", soft.id, tag);
            attrib!("version", soft.version, tag);
            self.handle.write_event(Event::Start(tag.to_borrowed()))?;
            for param in soft.params() {
                self.handle.write_param(param)?
            }
            self.handle.write_event(Event::End(tag.to_end()))?;
        }
        self.handle.write_event(Event::End(outer.to_end()))?;
        Ok(())
    }

    fn write_instrument_configuration(&mut self) -> XMLResult {
        let mut outer = bstart!("instrumentConfigurationList");
        let count = self.instrument_configurations.len().to_string();
        attrib!("count", count, outer);
        self.handle.write_event(Event::Start(outer.to_borrowed()))?;
        for (_key, ic) in self.instrument_configurations.iter() {
            let mut tag = bstart!("instrumentConfiguration");
            attrib!("id", ic.id, tag);
            self.handle.write_event(Event::Start(tag.to_borrowed()))?;
            for param in ic.params() {
                self.handle.write_param(param)?
            }
            for comp in ic.components.iter() {
                let mut cmp_tag = match comp.component_type {
                    ComponentType::Analyzer => bstart!("analyzer"),
                    ComponentType::Detector => bstart!("detector"),
                    ComponentType::IonSource => bstart!("source"),
                    ComponentType::Unknown => {
                        panic!("Could not identify component tag for {:?}", comp)
                    }
                };
                let order = comp.order.to_string();
                attrib!("order", order, cmp_tag);
                self.handle
                    .write_event(Event::Start(cmp_tag.to_borrowed()))?;
                for param in comp.params() {
                    self.handle.write_param(param)?
                }
                self.handle.write_event(Event::End(cmp_tag.to_end()))?;
            }
            let mut sw = bstart!("sofwareRef");
            attrib!("ref", ic.software_reference, sw);
            self.handle.write_event(Event::Empty(sw))?;
            self.handle.write_event(Event::End(tag.to_end()))?;
        }
        self.handle.write_event(Event::End(outer.to_end()))?;
        Ok(())
    }

    fn write_data_processing(&mut self) -> XMLResult {
        let mut outer = bstart!("dataProcessingList");
        let count = self.data_processings.len().to_string();
        attrib!("count", count, outer);
        self.handle.write_event(Event::Start(outer.to_borrowed()))?;
        for dp in self.data_processings.iter() {
            let mut tag = bstart!("dataProcessing");
            attrib!("id", dp.id, tag);
            self.handle.write_event(Event::Start(tag.to_borrowed()))?;
            for proc in dp.methods.iter() {
                let mut mtag = bstart!("processingMethod");
                let order = proc.order.to_string();
                attrib!("order", order, mtag);
                attrib!("softwareRef", proc.software_reference, mtag);
                self.handle.write_event(Event::Start(mtag.to_borrowed()))?;
                for param in proc.params() {
                    self.handle.write_param(param)?
                }
                self.handle.write_event(Event::End(mtag.to_end()))?;
            }
            self.handle.write_event(Event::End(tag.to_end()))?;
        }
        self.handle.write_event(Event::End(outer.to_end()))?;
        Ok(())
    }

    fn start_run(&mut self) -> XMLResult {
        if self.state < MzMLWriterState::Run {
            self.writer_header()?;
        } else {
            panic!(
                "Cannot start writing the run of mzML, currently in state {:?} which happens after.", self.state)
        }
        let mut run = bstart!("run");
        attrib!("id", "1", run);
        if let Some(ic_ref) = self.instrument_configurations.keys().next() {
            attrib!("defaultInstrumentConfigurationRef", ic_ref, run);
        }
        if let Some(sf_ref) = self.file_description.source_files.first() {
            attrib!("defaultSourceFileRef", sf_ref.id, run);
        };
        self.handle.write_event(Event::Start(run))?;
        self.state = MzMLWriterState::Run;
        Ok(())
    }

    fn start_spectrum_list(&mut self) -> XMLResult {
        if self.state < MzMLWriterState::SpectrumList {
            self.start_run()?;
        } else if MzMLWriterState::SpectrumList > self.state {
            panic!("Cannot start writing the run of mzML, currently in state {:?} which happens after the run has already begin", self.state)
        }
        let mut list = bstart!("spectrumList");
        let count = self.spectrum_count.to_string();
        attrib!("count", count, list);
        if let Some(dp) = self.data_processings.first() {
            attrib!("defaultDataProcessingRef", dp.id, list);
        }
        self.handle.write_event(Event::Start(list))?;
        self.state = MzMLWriterState::SpectrumList;
        Ok(())
    }

    fn close_spectrum_list(&mut self) -> XMLResult {
        let tag = bstart!("spectrumList");
        end_event!(self, tag);
        self.state = MzMLWriterState::SpectrumListClosed;
        Ok(())
    }

    fn close_run(&mut self) -> XMLResult {
        if self.state < MzMLWriterState::Run {
            self.start_run()?;
        } else if self.state == MzMLWriterState::SpectrumList {
            self.close_spectrum_list()?;
        } else if self.state == MzMLWriterState::ChromatogramList {
            // TODO
        } else if self.state > MzMLWriterState::RunClosed {
            panic!("Cannot close the run of mzML, currently in state {:?} which happens after the run has already ended", self.state)
        }
        let tag = bstart!("run");
        end_event!(self, tag);
        self.state = MzMLWriterState::RunClosed;
        Ok(())
    }

    fn close_mzml(&mut self) -> XMLResult {
        if self.state < MzMLWriterState::RunClosed {
            self.close_run()?;
        }
        let tag = bstart!("mzML");
        self.state = MzMLWriterState::MzMLClosed;
        end_event!(self, tag);
        Ok(())
    }

    fn close_indexed_mzml(&mut self) -> XMLResult {
        if self.state < MzMLWriterState::MzMLClosed {
            self.close_mzml()?;
        }
        self.write_index_list()?;
        let tag = bstart!("indexedmzML");
        end_event!(self, tag);
        self.state = MzMLWriterState::End;
        Ok(())
    }

    pub fn close(&mut self) -> XMLResult {
        if self.state < MzMLWriterState::End {
            self.close_indexed_mzml()
        } else {
            Ok(())
        }
    }

    fn write_scan_list(&mut self, acq: &Acquisition) -> XMLResult {
        let mut scan_list_tag = bstart!("scanList");
        let count = acq.scans.len().to_string();
        attrib!("count", count, scan_list_tag);
        start_event!(self, scan_list_tag);
        self.handle
            .write_param(&self.ms_cv.param("MS:1000016", "no combination"))?;

        for scan in acq.scans.iter() {
            let mut scan_tag = bstart!("scan");
            attrib!(
                "instrumentConfigurationRef",
                scan.instrument_configuration_id,
                scan_tag
            );
            self.handle
                .write_event(Event::Start(scan_tag.to_borrowed()))?;

            self.handle.write_param(
                &self
                    .ms_cv
                    .param_val("MS:1000016", "scan start time", scan.start_time)
                    .with_unit("UO:0000031", "minute"),
            )?;

            self.handle.write_param(
                &self
                    .ms_cv
                    .param_val("MS:1000927", "ion injection time", scan.injection_time)
                    .with_unit("UO:0000028", "millisecond"),
            )?;

            for param in scan.params() {
                self.handle.write_param(param)?
            }

            let mut scan_window_list_tag = bstart!("scanWindowList");
            let scan_window_list_count = scan.scan_windows.len().to_string();

            attrib!("count", scan_window_list_count, scan_window_list_tag);
            self.handle
                .write_event(Event::Start(scan_window_list_tag.to_borrowed()))?;
            for window in scan.scan_windows.iter() {
                let window_tag = bstart!("scanWindow");
                self.handle
                    .write_event(Event::Start(window_tag.to_borrowed()))?;
                self.handle.write_param(
                    &self
                        .ms_cv
                        .param_val(
                            "MS:1000501",
                            "scan window lower limit",
                            window.lower_bound.to_string(),
                        )
                        .with_unit("MS:1000040", "m/z"),
                )?;
                self.handle.write_param(
                    &self
                        .ms_cv
                        .param_val(
                            "MS:1000500",
                            "scan window upper limit",
                            window.upper_bound.to_string(),
                        )
                        .with_unit("MS:1000040", "m/z"),
                )?;
                self.handle.write_event(Event::End(window_tag.to_end()))?;
            }
            self.handle
                .write_event(Event::End(scan_window_list_tag.to_end()))?;
            self.handle.write_event(Event::End(scan_tag.to_end()))?;
        }
        end_event!(self, scan_list_tag);
        Ok(())
    }

    fn write_isolation_window(&mut self, iw: &IsolationWindow) -> XMLResult {
        let iw_tag = bstart!("isolationWindow");
        self.handle
            .write_event(Event::Start(iw_tag.to_borrowed()))?;
        self.handle.write_param(
            &self
                .ms_cv
                .param_val(
                    "MS:1000827",
                    "isolation window target m/z",
                    iw.target.to_string(),
                )
                .with_unit("MS:1000040", "m/z"),
        )?;
        self.handle.write_param(
            &self
                .ms_cv
                .param_val(
                    "MS:1000828",
                    "isolation window lower offset",
                    (iw.target - iw.lower_bound).to_string(),
                )
                .with_unit("MS:1000040", "m/z"),
        )?;
        self.handle.write_param(
            &self
                .ms_cv
                .param_val(
                    "MS:1000829",
                    "isolation window upper offset",
                    (iw.upper_bound - iw.target).to_string(),
                )
                .with_unit("MS:1000040", "m/z"),
        )?;
        self.handle.write_event(Event::End(iw_tag.to_end()))
    }

    fn write_selected_ions(&mut self, precursor: &Precursor) -> XMLResult {
        let mut outer = bstart!("selectedIonList");
        attrib!("count", "1", outer);
        start_event!(self, outer);
        let tag = bstart!("selectedIon");
        start_event!(self, tag);

        let ion = precursor.ion();
        self.handle.write_param(
            &self
                .ms_cv
                .param_val("MS:1000744", "selected ion m/z", ion.mz)
                .with_unit("MS:1000040", "m/z"),
        )?;
        self.handle.write_param(
            &self
                .ms_cv
                .param_val("MS:1000042", "peak intensity", ion.intensity)
                .with_unit("MS:1000131", "number of detector counts"),
        )?;
        if let Some(charge) = &ion.charge {
            self.handle
                .write_param(&self.ms_cv.param_val("MS:1000041", "charge state", charge))?;
        }
        for param in ion.params() {
            self.handle.write_param(param)?
        }
        end_event!(self, tag);
        end_event!(self, outer);
        Ok(())
    }

    fn write_activation(&mut self, precursor: &Precursor) -> XMLResult {
        let act = precursor.activation();
        let tag = bstart!("activation");
        start_event!(self, tag);
        self.handle.write_param_list(act.params().iter())?;
        self.handle.write_param(
            &self
                .ms_cv
                .param_val("MS:1000045", "collision energy", act.energy)
                .with_unit("UO:0000266", "electronvolt"),
        )?;
        end_event!(self, tag);
        Ok(())
    }

    fn write_precursor(&mut self, precursor: &Precursor) -> XMLResult {
        let mut precursor_list_tag = bstart!("precursorList");
        attrib!("count", "1", precursor_list_tag);
        start_event!(self, precursor_list_tag);

        let mut precursor_tag = bstart!("precursor");
        attrib!("spectrumRef", precursor.precursor_id, precursor_tag);
        self.handle
            .write_event(Event::Start(precursor_tag.to_borrowed()))?;

        let iw = precursor.isolation_window();
        self.write_isolation_window(iw)?;
        self.write_selected_ions(precursor)?;
        self.write_activation(precursor)?;
        end_event!(self, precursor_tag);
        end_event!(self, precursor_list_tag);
        Ok(())
    }

    fn write_binary_data_array(&mut self, array: &DataArray) -> XMLResult {
        let mut outer = bstart!("binaryDataArray");

        let encoded = array.encode_bytestring(self.data_array_compression);
        let encoded_len = encoded.len().to_string();
        attrib!("encodedLength", encoded_len, outer);

        start_event!(self, outer);
        match &array.dtype {
            BinaryDataArrayType::Float32 => self
                .handle
                .write_param(&self.ms_cv.param("MS:1000521", "32-bit float"))?,
            BinaryDataArrayType::Float64 => self
                .handle
                .write_param(&self.ms_cv.param("MS:1000523", "64-bit float"))?,
            BinaryDataArrayType::Int32 => self
                .handle
                .write_param(&self.ms_cv.param("MS:1000519", "32-bit integer"))?,
            BinaryDataArrayType::Int64 => self
                .handle
                .write_param(&self.ms_cv.param("MS:1000522", "64-bit integer"))?,
            BinaryDataArrayType::ASCII => self.handle.write_param(
                &self
                    .ms_cv
                    .param("MS:1001479 ", "null-terminated ASCII string"),
            )?,
            _ => {
                panic!(
                    "Could not determine data type for binary data array. Found {:?}",
                    array.dtype
                )
            }
        }
        if self.data_array_compression == BinaryCompressionType::NoCompression {
            self.handle
                .write_param(&self.ms_cv.param("MS:1000576", "no compression"))?;
        } else {
            self.handle
                .write_param(&self.ms_cv.param("MS:1000574", "zlib compression"))?;
        }
        match &array.name {
            ArrayType::MZArray => self.handle.write_param(
                &self
                    .ms_cv
                    .param("MS:1000514", "m/z array")
                    .with_unit("MS:1000040", "m/z"),
            )?,
            ArrayType::IntensityArray => self.handle.write_param(
                &self
                    .ms_cv
                    .param("MS:1000515", "intensity array")
                    .with_unit("MS:1000131", "number of detector counts"),
            )?,
            ArrayType::ChargeArray => self
                .handle
                .write_param(&self.ms_cv.param("MS:1000516", "charge array"))?,
            ArrayType::TimeArray => self.handle.write_param(
                &self
                    .ms_cv
                    .param("MS:1000595", "time array")
                    .with_unit("UO:0000031", "minute"),
            )?,
            ArrayType::RawIonMobilityArray => self.handle.write_param(
                &self
                    .ms_cv
                    .param("MS:1003007", "raw ion mobility array")
                    .with_unit("UO:0000028", "millisecond"),
            )?,
            ArrayType::MeanIonMobilityArray => self.handle.write_param(
                &self
                    .ms_cv
                    .param("MS:1002816", "mean ion mobility array")
                    .with_unit("UO:0000028", "millisecond"),
            )?,
            ArrayType::DeconvolutedIonMobilityArray => self.handle.write_param(
                &self
                    .ms_cv
                    .param("MS:1003154", "deconvoluted ion mobility array")
                    .with_unit("UO:0000028", "millisecond"),
            )?,
            ArrayType::NonStandardDataArray { name } => self.handle.write_param(
                &self
                    .ms_cv
                    .param_val("MS:1000786", "non-standard data array", name)
                    .with_unit("UO:0000028", "millisecond"),
            )?,
            _ => {
                panic!("Could not determine how to name for {:?}", array.name);
            }
        }

        let bin = bstart!("binary");
        start_event!(self, bin);
        self.handle
            .write_event(Event::Text(BytesText::from_plain(&encoded)))?;
        end_event!(self, bin);
        end_event!(self, outer);
        Ok(())
    }

    fn write_binary_data_arrays(&mut self, arrays: &BinaryArrayMap) -> XMLResult {
        let count = arrays.len().to_string();
        let mut outer = bstart!("binaryDataArrayList");
        attrib!("count", count, outer);
        start_event!(self, outer);
        let mut array_pairs: Vec<(&ArrayType, &DataArray)> = arrays.iter().collect();
        array_pairs.sort_by_key(|f| f.0);
        for (_tp, array) in array_pairs {
            self.write_binary_data_array(array)?
        }
        end_event!(self, outer);
        Ok(())
    }

    pub fn write_spectrum(&mut self, spectrum: &'a MultiLayerSpectrum<C, D>) -> XMLResult {
        if self.state < MzMLWriterState::SpectrumList {
            self.start_spectrum_list()?;
        } else if self.state > MzMLWriterState::SpectrumList {
            panic!("Cannot write spectrum, currently in state {:?} which happens after spectra may be written", self.state)
        }
        let pos = self.stream_position()? - (1 + (4 * InnerXMLWriter::<W>::INDENT_SIZE));
        self.offset_index.insert(spectrum.id().to_string(), pos);
        let mut outer = bstart!("spectrum");
        attrib!("id", spectrum.id(), outer);
        let count = self.spectrum_counter.to_string();
        attrib!("index", count, outer);
        self.handle.write_event(Event::Start(outer.to_borrowed()))?;
        self.spectrum_counter += 1;

        let ms_level = spectrum.ms_level();
        if ms_level == 1 {
            self.handle
                .write_param(&self.ms_cv.param("MS:1000579", "MS1 spectrum"))?;
        } else {
            self.handle
                .write_param(&self.ms_cv.param("MS:1000580", "MSn spectrum"))?;
        }
        self.handle.write_param(&self.ms_cv.param_val(
            "MS:1000511",
            "ms level",
            ms_level.to_string(),
        ))?;

        match spectrum.polarity() {
            ScanPolarity::Negative => self
                .handle
                .write_param(&self.ms_cv.param("MS:1000129", "negative scan"))?,
            ScanPolarity::Positive => self
                .handle
                .write_param(&self.ms_cv.param("MS:1000130", "positive scan"))?,
            ScanPolarity::Unknown => {
                warn!(
                    "Could not determine scan polarity for {}, assuming positive",
                    spectrum.id()
                );
                self.handle
                    .write_param(&self.ms_cv.param("MS:1000130", "positive scan"))?
            }
        }

        match spectrum.signal_continuity() {
            SignalContinuity::Profile => self
                .handle
                .write_param(&self.ms_cv.param("MS:1000128", "profile spectrum"))?,
            SignalContinuity::Unknown => {
                warn!(
                    "Could not determine scan polarity for {}, assuming centroid",
                    spectrum.id()
                );
                self.handle
                    .write_param(&self.ms_cv.param("MS:1000127", "centroid spectrum"))?;
            }
            _ => {
                self.handle
                    .write_param(&self.ms_cv.param("MS:1000127", "centroid spectrum"))?;
            }
        }

        let acq = spectrum.acquisition();

        self.write_scan_list(acq)?;

        if let Some(precursor) = spectrum.precursor() {
            self.write_precursor(precursor)?;
        }

        if let Some(mass_peaks) = &spectrum.deconvoluted_peaks {
            let arrays: BinaryArrayMap = mass_peaks.into();
            self.write_binary_data_arrays(&arrays)?
        } else if let Some(mz_peaks) = &spectrum.peaks {
            let arrays = mz_peaks.into();
            self.write_binary_data_arrays(&arrays)?
        } else if let Some(arrays) = &spectrum.arrays {
            self.write_binary_data_arrays(arrays)?
        }

        end_event!(self, outer);
        Ok(())
    }

    fn write_index(&mut self, index: &OffsetIndex) -> XMLResult {
        let mut outer = bstart!("index");
        attrib!("name", index.name, outer);
        start_event!(self, outer);
        for (id, offset) in index.iter() {
            let mut tag = bstart!("offset");
            attrib!("idRef", id, tag);
            start_event!(self, tag);
            let content = offset.to_string();
            let text = BytesText::from_plain_str(&content);
            self.handle.write_event(Event::Text(text))?;
            end_event!(self, tag);
        }
        end_event!(self, outer);
        Ok(())
    }

    fn write_index_list(&mut self) -> XMLResult {
        if self.state < MzMLWriterState::IndexList {
            self.close_mzml()?;
        }
        let offset = self.stream_position()?;
        let mut outer = bstart!("indexList");
        attrib!("count", "1", outer);
        start_event!(self, outer);
        self.write_index(&self.offset_index.clone())?;
        end_event!(self, outer);

        let tag = bstart!("indexListOffset");
        start_event!(self, tag);
        let content = offset.to_string();
        let text = BytesText::from_plain_str(&content);
        self.handle.write_event(Event::Text(text))?;
        end_event!(self, tag);

        let tag = bstart!("fileChecksum");
        start_event!(self, tag);
        let content = self.handle.digest();
        let text = BytesText::from_plain_str(&content);
        self.handle.write_event(Event::Text(text))?;
        end_event!(self, tag);
        Ok(())
    }

    /// Get a reference to the mz m l writer type's spectrum count.
    pub fn spectrum_count(&self) -> &u64 {
        &self.spectrum_count
    }

    /// Set the mz m l writer type's spectrum count.
    pub fn set_spectrum_count(&mut self, spectrum_count: u64) {
        self.spectrum_count = spectrum_count;
    }

    /// Get a mutable reference to the mz m l writer type's spectrum count.
    pub fn spectrum_count_mut(&mut self) -> &mut u64 {
        &mut self.spectrum_count
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::spectrum::spectrum::SpectrumBehavior;
    use std::fs;
    use std::path;

    #[test]
    fn reader_from_file() {
        let path = path::Path::new("./test/data/small.mzML");
        let file = fs::File::open(path).expect("Test file doesn't exist");
        let reader = MzMLReaderType::<_, CentroidPeak, DeconvolutedPeak>::new(file);
        let mut ms1_count = 0;
        let mut msn_count = 0;

        assert_eq!(reader.data_processings().len(), 1);
        assert_eq!(reader.instrument_configurations().len(), 2);
        assert_eq!(reader.softwares().len(), 2);
        assert_eq!(reader.file_description().source_files.len(), 1);
        assert_eq!(reader.file_description().contents.len(), 2);

        assert!(reader
            .file_description()
            .get_param_by_accession("MS:1000579")
            .is_some());
        assert_eq!(reader.file_description().source_files[0].name, "small.RAW");

        for scan in reader {
            let level = scan.ms_level();
            if level == 1 {
                ms1_count += 1;
            } else {
                msn_count += 1;
            }
        }
        assert_eq!(ms1_count, 14);
        assert_eq!(msn_count, 34);
    }

    #[test]
    fn reader_from_file_indexed() {
        let path = path::Path::new("./test/data/small.mzML");
        let file = fs::File::open(path).expect("Test file doesn't exist");
        let mut reader = MzMLReaderType::<_, CentroidPeak, DeconvolutedPeak>::new_indexed(file);

        let n = reader.len();
        assert_eq!(n, 48);

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
        assert_eq!(ms1_count, 14);
        assert_eq!(msn_count, 34);

        assert_eq!(reader.data_processings().len(), 1);
        assert_eq!(reader.instrument_configurations().len(), 2);
        assert_eq!(reader.softwares().len(), 2);
        assert_eq!(reader.file_description().source_files.len(), 1);
        assert_eq!(reader.file_description().contents.len(), 2);
    }

    #[test]
    fn reader_from_path() {
        let path = path::Path::new("./test/data/small.mzML");
        let mut reader = MzMLReaderType::<_, CentroidPeak, DeconvolutedPeak>::open_path(path)
            .expect("Test file doesn't exist?");

        let n = reader.len();
        assert_eq!(n, 48);

        let mut ms1_count = 0;
        let mut msn_count = 0;

        for i in (0..n).rev() {
            let scan = match reader.get_spectrum_by_index(i) {
                Some(scan) => scan,
                None => {
                    if let Some(offset) = reader._offset_of_index(i) {
                        panic!(
                            "Failed to locate spectrum {} at offset {}, parser state {:?}",
                            i, offset, reader.state,
                        );
                    } else {
                        panic!("Failed to locate spectrum or offset {}", i);
                    }
                }
            };
            let level = scan.ms_level();
            if level == 1 {
                ms1_count += 1;
            } else {
                msn_count += 1;
            }
        }
        assert_eq!(ms1_count, 14);
        assert_eq!(msn_count, 34);
    }

    #[test]
    fn grouped_iteration() {
        let path = path::Path::new("./test/data/small.mzML");
        let mut reader = MzMLReaderType::<_, CentroidPeak, DeconvolutedPeak>::open_path(path)
            .expect("Test file doesn't exist?");

        let n = reader.len();
        assert_eq!(n, 48);

        let mut ms1_count = 0;
        let mut msn_count = 0;

        for group in reader.groups() {
            ms1_count += group.precursor.is_some() as usize;
            msn_count += group.products.len();
        }
        assert_eq!(ms1_count, 14);
        assert_eq!(msn_count, 34);
    }

    #[test]
    fn write_test() -> XMLResult {
        let path = path::Path::new("./test/data/small.mzML");
        let mut reader = MzMLReaderType::<_, CentroidPeak, DeconvolutedPeak>::open_path(path)
            .expect("Test file doesn't exist?");

        let n = reader.len();
        assert_eq!(n, 48);

        let dest = fs::File::create("./test/data/duplicate.mzML")?;
        let mut writer = MzMLWriterType::new(dest);
        writer.copy_metadata_from(&reader);
        *writer.spectrum_count_mut() = reader.len() as u64;
        for group in reader.groups() {
            for prec in group.precursor.iter() {
                writer.write_spectrum(prec)?
            }
            for prod in group.products.iter() {
                writer.write_spectrum(prod)?
            }
        }
        writer.close()?;

        let mut reader2 = MzMLReader::open_path(path)?;
        assert_eq!(reader.file_description(), reader2.file_description());

        for (a, b) in reader.iter().zip(reader2.iter()) {
            assert_eq!(a.id(), b.id());
            assert_eq!(a.ms_level(), b.ms_level());
            assert_eq!(a.index(), b.index());
            for (x, y) in a
                .arrays
                .unwrap()
                .mzs()
                .iter()
                .zip(b.arrays.unwrap().mzs().iter())
            {
                assert!((x - y).abs() < 1e-3)
            }
        }

        Ok(())
    }
}
