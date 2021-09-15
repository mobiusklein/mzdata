use std::collections::HashMap;
use std::convert::TryInto;
use std::fmt::Debug;
use std::fs;
use std::io;
use std::mem;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::marker::PhantomData;

use log::warn;

use regex;

use quick_xml::events::{BytesEnd, BytesStart, BytesText, Event};
use quick_xml::Error as XMLError;
use quick_xml::Reader;

use super::super::offset_index::OffsetIndex;
use super::super::traits::{
    MZFileReader, RandomAccessScanIterator, ScanAccessError, ScanSource, SeekRead,
};

use mzpeaks::{CentroidPeak, DeconvolutedPeak};

use crate::meta::file_description::{FileDescription, SourceFile};
use crate::meta::instrument::{Component, ComponentType, InstrumentConfiguration};
use crate::meta::{DataProcessing, MSDataFileMetadata, ProcessingMethod, Software};
use crate::params::{Param, ParamList, Unit};
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
        MzMLParserError::XMLError(state, error)
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
                        param.name = attr.unescape_and_decode_value(reader).unwrap_or_else(|e| {
                            panic!(
                                "Error decoding CV param name at {}: {}",
                                reader.buffer_position(),
                                e
                            )
                        });
                    }
                    b"value" => {
                        param.value = attr.unescape_and_decode_value(reader).unwrap_or_else(|e| {
                            panic!(
                                "Error decoding CV param value at {}: {}",
                                reader.buffer_position(),
                                e
                            )
                        });
                    }
                    b"cvRef" => {
                        param.controlled_vocabulary =
                            Some(attr.unescape_and_decode_value(reader).unwrap_or_else(|e| {
                                panic!(
                                    "Error decoding CV param reference at {}: {}",
                                    reader.buffer_position(),
                                    e
                                )
                            }));
                    }
                    b"accession" => {
                        param.accession =
                            attr.unescape_and_decode_value(reader).unwrap_or_else(|e| {
                                panic!(
                                    "Error decoding CV param accession at {}: {}",
                                    reader.buffer_position(),
                                    e
                                )
                            });
                    }
                    b"unitName" => {
                        param.unit_name =
                            Some(attr.unescape_and_decode_value(reader).unwrap_or_else(|e| {
                                panic!(
                                    "Error decoding CV param unit name at {}: {}",
                                    reader.buffer_position(),
                                    e
                                )
                            }));
                    }
                    b"unitAccession" => {
                        param.unit_accession =
                            Some(attr.unescape_and_decode_value(reader).unwrap_or_else(|e| {
                                panic!(
                                    "Error decoding CV param unit name at {}: {}",
                                    reader.buffer_position(),
                                    e
                                )
                            }));
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

#[derive(Debug)]
pub enum MzMLParserError {
    NoError,
    UnknownError(MzMLParserState),
    IncompleteSpectrum,
    IncompleteElementError(String, MzMLParserState),
    XMLError(MzMLParserState, XMLError),
    IOError(MzMLParserState, io::Error)
}

impl std::fmt::Display for MzMLParserError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::error::Error for MzMLParserError {}


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
            "mean ion mobility array" => {
                self.current_array_mut().name = ArrayType::MeanIonMobilityArray(Unit::Unknown)
            }
            "mean drift time array" => {
                self.current_array_mut().name = ArrayType::MeanIonMobilityArray(Unit::Millisecond)
            }
            "mean inverse reduced ion mobility array" => {
                self.current_array_mut().name =
                    ArrayType::MeanIonMobilityArray(Unit::VoltSecondPerSquareCentimeter)
            }
            "ion mobility array" => {
                self.current_array_mut().name = ArrayType::IonMobilityArray(Unit::Unknown)
            }
            "drift time array" => {
                self.current_array_mut().name = ArrayType::IonMobilityArray(Unit::Millisecond)
            }
            "inverse reduced ion mobility array" => {
                self.current_array_mut().name =
                    ArrayType::IonMobilityArray(Unit::VoltSecondPerSquareCentimeter)
            }
            "deconvoluted ion mobility array" => {
                self.current_array_mut().name =
                    ArrayType::DeconvolutedIonMobilityArray(Unit::Unknown)
            }
            "deconvoluted drift time array" => {
                self.current_array_mut().name =
                    ArrayType::DeconvolutedIonMobilityArray(Unit::Millisecond)
            }
            "deconvoluted inverse reduced ion mobility array" => {
                self.current_array_mut().name =
                    ArrayType::DeconvolutedIonMobilityArray(Unit::VoltSecondPerSquareCentimeter)
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
                if let IsolationWindowState::Unknown = window.flags {
                    window.flags = IsolationWindowState::Explicit;
                    window.lower_bound = lower_bound;
                }
            }
            "isolation window upper limit" => {
                let upper_bound: f64 = param
                    .coerce()
                    .expect("Failed to parse isolation window limit");
                if let IsolationWindowState::Unknown = window.flags {
                    window.flags = IsolationWindowState::Explicit;
                    window.upper_bound = upper_bound;
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

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> Into<MultiLayerSpectrum<C, D>>
    for MzMLSpectrumBuilder<C, D>
{
    fn into(self) -> MultiLayerSpectrum<C, D> {
        let mut spec = MultiLayerSpectrum::<C, D>::default();
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

/**A SAX-style parser for building up the metadata section prior to the `<run>` element
of an mzML file.*/
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
                let mut source = Component {
                    component_type: ComponentType::IonSource,
                    ..Default::default()
                };
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
                let mut analyzer = Component {
                    component_type: ComponentType::Analyzer,
                    ..Default::default()
                };
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
                let mut detector = Component {
                    component_type: ComponentType::Detector,
                    ..Default::default()
                };
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

#[derive(Debug, Default, Clone)]
pub struct IndexedMzMLIndexExtractor {
    spectrum_index: OffsetIndex,
    chromatogram_index: OffsetIndex,
    index_offset: Option<u64>,
    last_id: String
}

#[derive(Debug, Clone, Copy)]
pub enum IndexParserState {
    Start,
    SpectrumIndexList,
    ChromatogramIndexList,
    Done
}

impl XMLParseBase for IndexedMzMLIndexExtractor {}

impl IndexedMzMLIndexExtractor {
    pub fn new() -> IndexedMzMLIndexExtractor {
        IndexedMzMLIndexExtractor {
            spectrum_index: OffsetIndex::new("spectrum".into()),
            chromatogram_index: OffsetIndex::new("chromatogram".into()),
            index_offset: None,
            last_id: String::new(),
        }
    }

    pub fn find_offset_from_reader<R: SeekRead>(&self, reader: &mut R) -> io::Result<Option<u64>>{
        reader.seek(SeekFrom::End(-200))?;
        let mut buf = String::new();
        reader.read_to_string(&mut buf)?;
        let pattern = regex::Regex::new("<indexListOffset>(\\d+)</indexListOffset>").unwrap();
        if let Some(captures) = pattern.captures(&buf) {
            if let Some(offset) = captures.get(1) {
                if let Ok(offset) = offset.as_str().parse::<u64>() {
                    return Ok(Some(offset));
                }
            }
        }
        Ok(None)
    }

    pub fn start_element<B: io::BufRead>(
        &mut self,
        event: &BytesStart,
        state: IndexParserState,
        reader: &Reader<B>,
    ) -> Result<IndexParserState, XMLError> {
        let elt_name = event.name();
        match elt_name {
            b"offset" => {
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key == b"idRef" {
                                self.last_id = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding idRef");
                            }
                        },
                        Err(err) => {
                            return Err(err);
                        }
                    }
                }
            }
            b"index" => {
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key == b"name" {
                                let index_name = attr
                                    .unescape_and_decode_value(reader)
                                    .expect("Error decoding idRef");
                                match index_name.as_ref() {
                                    "spectrum" => {
                                        return Ok(IndexParserState::SpectrumIndexList)
                                    },
                                    "chromatogram" => {
                                        return Ok(IndexParserState::ChromatogramIndexList)
                                    },
                                    _ => {},
                                }
                            }
                        },
                        Err(err) => {
                            return Err(err);
                        }
                    }
                }
            }
            b"indexList" => {}
            _ => {}
        }

        Ok(state)
    }

    pub fn end_element<B: io::BufRead>(
        &mut self,
        event: &BytesEnd,
        state: IndexParserState,
        _reader: &Reader<B>,
    ) -> Result<IndexParserState, XMLError> {
        let elt_name = event.name();
        match elt_name {
            b"offset" => {}
            b"index" => {}
            b"indexList" => {
                return Ok(IndexParserState::Done)
            }
            _ => {}
        }
        Ok(state)
    }

    pub fn text(&mut self, event: &BytesText, state: IndexParserState) -> Result<IndexParserState, XMLError> {
        match state {
            IndexParserState::SpectrumIndexList => {
                let bin = event
                    .unescaped()
                    .expect("Failed to unescape spectrum offset");
                if let Ok(offset) = String::from_utf8_lossy(&*bin).parse::<u64>() {
                    if self.last_id != "" {
                        let key = mem::take(&mut self.last_id);
                        self.spectrum_index.insert(key, offset);
                    } else {
                        warn!("Out of order text in index")
                    }
                }
            }
            IndexParserState::ChromatogramIndexList => {
                let bin = event
                    .unescaped()
                    .expect("Failed to unescape chromatogram offset");
                if let Ok(offset) = String::from_utf8_lossy(&*bin).parse::<u64>() {
                    if self.last_id != "" {
                        let key = mem::take(&mut self.last_id);
                        self.chromatogram_index.insert(key, offset);
                    } else {
                        warn!("Out of order text in index")
                    }
                }
            }
            _ => {}
        }
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
    error: MzMLParserError,
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
    /// Create a new [`MzMLReaderType`] instance, wrapping the [`io::Read`] handle
    /// provided with an `[io::BufReader`] and parses the metadata section of the file.
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

    /**Parse the metadata section of the file using [`FileMetadataBuilder`]
     */
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
            MzMLParserState::ParserError => {
                let mut error = MzMLParserError::NoError;
                mem::swap(&mut error, &mut self.error);
                Err(error)
            },
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
            MzMLParserState::ParserError => {
                let mut error = MzMLParserError::NoError;
                mem::swap(&mut error, &mut self.error);
                Err(error)
            },
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
            Ok(_sz) => Some(spectrum),
            Err(_err) => None,
        }
    }
}


#[derive(Debug)]
pub enum MzMLIndexingError {
    OffsetNotFound,
    XMLError(XMLError),
    IOError(io::Error)
}

impl<R: io::Read + io::Seek, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MzMLReaderType<R, C, D> {}


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
        match reader.read_index_from_end() {
            Ok(()) => {}
            Err(err) => {
                match reader.seek(SeekFrom::Start(0)) {
                    Ok(_) => {
                        reader.build_index();
                    },
                    Err(error) => {
                        panic!("Unrecoverable IO Error during file pointer reset {} while handling {:?}", error, err);
                    }
                }
            }
        }
        reader
    }

    pub fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.handle.seek(pos)
    }

    /// Read the offset index at the end of an `<indexedmzML>` document,
    /// though this index may be malformed in some older files.
    pub fn read_index_from_end(&mut self) -> Result<(), MzMLIndexingError> {
        let mut indexer = IndexedMzMLIndexExtractor::new();

        let offset = match indexer.find_offset_from_reader(&mut self.handle) {
            Ok(offset) => {
                if let Some(offset) = offset {
                    offset
                } else {
                    return Err(MzMLIndexingError::OffsetNotFound)
                }
            },
            Err(err) => {
                return Err(MzMLIndexingError::IOError(err))
            }
        };

        let current_position = match self.handle.stream_position() {
            Ok(position) => {
                position
            },
            Err(err) => {
                return Err(MzMLIndexingError::IOError(err))
            }
        };

        self.handle.seek(SeekFrom::Start(offset)).unwrap();
        let mut indexer_state = IndexParserState::Start;
        let mut reader = Reader::from_reader(&mut self.handle);
        reader.trim_text(true);

        loop {
            match reader.read_event(&mut self.buffer) {
                Ok(Event::Start(ref e)) => {
                    match indexer.start_element(e, indexer_state, &reader) {
                        Ok(state) => {
                            indexer_state = state;
                            match &indexer_state {
                                IndexParserState::Done => break,
                                _ => {}
                            }
                        }
                        Err(message) => {
                            return Err(MzMLIndexingError::XMLError(message))
                        }
                    };
                }
                Ok(Event::End(ref e)) => {
                    match indexer.end_element(e, indexer_state, &reader) {
                        Ok(state) => {
                            indexer_state = state;
                        }
                        Err(message) => {
                            return Err(MzMLIndexingError::XMLError(message))
                        }
                    };
                }
                Ok(Event::Text(ref e)) => {
                    match indexer.text(e, indexer_state) {
                        Ok(state) => {
                            indexer_state = state;
                        }
                        Err(message) => {
                            return Err(MzMLIndexingError::XMLError(message))
                        }
                    };
                }
                Ok(Event::Eof) => {
                    break;
                }
                Err(err) => {
                    return Err(MzMLIndexingError::XMLError(err))
                }
                _ => {}
            }
        }
        self.buffer.clear();
        self.index = indexer.spectrum_index;
        self.index.init = true;
        self.handle.seek(SeekFrom::Start(current_position)).unwrap();
        Ok(())
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
        if self.index.is_empty() {
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
    fn find_offset() -> io::Result<()> {
        let path = path::Path::new("./test/data/small.mzML");
        let mut f = fs::File::open(path)?;

        let index = IndexedMzMLIndexExtractor::new();
        if let Ok(offset) = index.find_offset_from_reader(&mut f) {
            if let Some(offset) = offset {
                assert_eq!(offset, 3183414);
            } else {
                panic!("Failed to parse offset from element")
            }
        } else {
            panic!("Failed to find offset element")
        }

        Ok(())
    }

    #[test]
    fn read_index() -> io::Result<()> {
        let path = path::Path::new("./test/data/duplicate.mzML");
        let file = fs::File::open(path)?;
        let mut reader = MzMLReader::new(file);
        match reader.read_index_from_end() {
            Ok(()) => {
            },
            Err(err) => {
                panic!("Failed to parse out index {:?}", err);
            }
        };
        assert!(reader.index.len() > 0);
        Ok(())
    }
}
