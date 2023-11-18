use std::collections::HashMap;
use std::io::SeekFrom;
use std::{io, mem};

use log::warn;
use quick_xml::events::{BytesEnd, BytesStart, BytesText};
use quick_xml::Error as XMLError;

use thiserror::Error;

use crate::ParamDescribed;
use crate::io::OffsetIndex;
use crate::io::traits::SeekRead;
use crate::meta::{FileDescription, InstrumentConfiguration, Software, DataProcessing, SourceFile, ComponentType, Component, ProcessingMethod};
use crate::params::{ParamCow, Param, Unit, ControlledVocabulary, curie_to_num};

use super::reader::Bytes;

/**
The different states the [`MzMLReaderType`] can enter while parsing
different phases of the document. This information is really only
needed by the module consumer to determine where in the document an
error occurred.
*/
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

/**
All the ways that mzML parsing can go wrong
*/
#[derive(Debug, Error)]
pub enum MzMLParserError {
    #[error("An error occurred outside of normal conditions {0:?}")]
    UnknownError(MzMLParserState),
    #[error("An incomplete spectrum was parsed")]
    IncompleteSpectrum,
    #[error("An incomplete element {0} was encountered in {1:?}")]
    IncompleteElementError(String, MzMLParserState),
    #[error("An XML error {1:?} was encountered in {0:?}")]
    XMLError(MzMLParserState, #[source] XMLError),
    #[error("An IO error {1} was encountered in {0:?}")]
    IOError(MzMLParserState, #[source] io::Error),
}


impl From<MzMLParserError> for io::Error {
    fn from(value: MzMLParserError) -> Self {
        match value {
            MzMLParserError::IOError(_, ref e) => io::Error::new(e.kind(), value),
            _ => io::Error::new(io::ErrorKind::InvalidData, value),
        }
    }
}


pub type ParserResult = Result<MzMLParserState, MzMLParserError>;

/**
Common XML error handling behaviors
*/
pub trait XMLParseBase {
    fn handle_xml_error(&self, error: quick_xml::Error, state: MzMLParserState) -> MzMLParserError {
        MzMLParserError::XMLError(state, error)
    }
}

/**
Common `CVParam` parsing behaviors
*/
pub trait CVParamParse: XMLParseBase {
    fn handle_param_borrowed<'inner, 'outer: 'inner + 'event, 'event: 'inner>(
        event: &'event BytesStart<'event>,
        reader_position: usize,
        state: MzMLParserState,
    ) -> Result<ParamCow<'inner>, MzMLParserError> {
        let mut name = None;
        let mut value = None;
        let mut accession = None;
        let mut controlled_vocabulary = None;
        let mut unit = Unit::Unknown;

        for attr_parsed in event.attributes() {
            match attr_parsed {
                Ok(attr) => match attr.key.as_ref() {
                    b"name" => {
                        name = Some(attr.unescape_value().unwrap_or_else(|e| {
                            panic!("Error decoding CV param name at {}: {}", reader_position, e)
                        }));
                    }
                    b"value" => {
                        value = Some(attr.unescape_value().unwrap_or_else(|e| {
                            panic!(
                                "Error decoding CV param value at {}: {}",
                                reader_position, e
                            )
                        }));
                    }
                    b"cvRef" => {
                        let cv_id = attr.unescape_value().unwrap_or_else(|e| {
                            panic!(
                                "Error decoding CV param reference at {}: {}",
                                reader_position, e
                            )
                        });
                        controlled_vocabulary = cv_id
                            .parse::<ControlledVocabulary>()
                            .unwrap_or_else(|_|
                                panic!("Failed to parse controlled vocabulary ID {}", cv_id)
                            )
                            .as_option();
                    }
                    b"accession" => {
                        let v = attr.unescape_value().unwrap_or_else(|e| {
                            panic!(
                                "Error decoding CV param accession at {}: {}",
                                reader_position, e
                            )
                        });
                        let (_, acc) = curie_to_num(&v);
                        accession = acc;
                    }
                    b"unitName" => {
                        let v = attr.unescape_value().unwrap_or_else(|e| {
                            panic!(
                                "Error decoding CV param unit name at {}: {}",
                                reader_position, e
                            )
                        });
                        unit = Unit::from_name(&v);
                    }
                    b"unitAccession" => {
                        let v = attr.unescape_value().unwrap_or_else(|e| {
                            panic!(
                                "Error decoding CV param unit name at {}: {}",
                                reader_position, e
                            )
                        });
                        unit = Unit::from_accession(&v);
                    }
                    b"unitCvRef" => {}
                    _ => {}
                },
                Err(msg) => return Err(MzMLParserError::XMLError(state, msg.into())),
            }
        }
        let param = ParamCow::new(
            name.unwrap(),
            value.unwrap_or_default(),
            accession,
            controlled_vocabulary,
            unit,
        );
        Ok(param)
    }

    fn handle_param(
        event: &BytesStart,
        reader_position: usize,
        state: MzMLParserState,
    ) -> Result<Param, MzMLParserError> {
        let mut param = Param::new();
        let mut unit_name = None;
        let mut unit_accession = None;
        for attr_parsed in event.attributes() {
            match attr_parsed {
                Ok(attr) => match attr.key.as_ref() {
                    b"name" => {
                        param.name = attr
                            .unescape_value()
                            .unwrap_or_else(|e| {
                                panic!("Error decoding CV param name at {}: {}", reader_position, e)
                            })
                            .to_string();
                    }
                    b"value" => {
                        param.value = attr
                            .unescape_value()
                            .unwrap_or_else(|e| {
                                panic!(
                                    "Error decoding CV param value at {}: {}",
                                    reader_position, e
                                )
                            })
                            .to_string();
                    }
                    b"cvRef" => {
                        let cv_id = attr.unescape_value().unwrap_or_else(|e| {
                            panic!(
                                "Error decoding CV param reference at {}: {}",
                                reader_position, e
                            )
                        });
                        param.controlled_vocabulary = cv_id
                            .parse::<ControlledVocabulary>()
                            .unwrap_or_else(|_|
                                panic!("Failed to parse controlled vocabulary ID {}", cv_id)
                            )
                            .as_option();
                    }
                    b"accession" => {
                        let v = attr.unescape_value().unwrap_or_else(|e| {
                            panic!(
                                "Error decoding CV param accession at {}: {}",
                                reader_position, e
                            )
                        });
                        let (_, acc) = curie_to_num(&v);
                        param.accession = acc;
                    }
                    b"unitName" => {
                        let v = attr.unescape_value().unwrap_or_else(|e| {
                            panic!(
                                "Error decoding CV param unit name at {}: {}",
                                reader_position, e
                            )
                        });
                        unit_name = Some(Unit::from_name(&v));
                    }
                    b"unitAccession" => {
                        let v = attr.unescape_value().unwrap_or_else(|e| {
                            panic!(
                                "Error decoding CV param unit name at {}: {}",
                                reader_position, e
                            )
                        });
                        unit_accession = Some(Unit::from_accession(&v));
                    }
                    b"unitCvRef" => {}
                    _ => {}
                },
                Err(msg) => return Err(MzMLParserError::XMLError(state, msg.into())),
            }
        }
        if let Some(unit_acc) = unit_accession {
            match unit_acc {
                Unit::Unknown => {}
                _ => {
                    param.unit = unit_acc;
                }
            }
        }
        if let Some(unit_name) = unit_name {
            match unit_name {
                Unit::Unknown => {}
                _ => {
                    param.unit = unit_name;
                }
            }
        }
        Ok(param)
    }
}

pub trait MzMLSAX {
    fn start_element(&mut self, event: &BytesStart, state: MzMLParserState) -> ParserResult;

    fn empty_element(
        &mut self,
        event: &BytesStart,
        state: MzMLParserState,
        reader_position: usize,
    ) -> ParserResult;

    fn end_element(&mut self, event: &BytesEnd, state: MzMLParserState) -> ParserResult;

    fn text(&mut self, event: &BytesText, state: MzMLParserState) -> ParserResult;
}


#[derive(Debug, Error)]
pub enum MzMLIndexingError {
    #[error("Offset index not found")]
    OffsetNotFound,
    #[error("XML error {0} occurred while reading out mzML index")]
    XMLError(#[from] #[source] XMLError),
    #[error("IO error {0} occurred while reading out mzML index")]
    IOError(#[from] #[source] io::Error),
}


#[derive(Debug, Default, Clone)]
pub struct IndexedMzMLIndexExtractor {
    pub spectrum_index: OffsetIndex,
    pub chromatogram_index: OffsetIndex,
    last_id: String,
}

#[derive(Debug, Clone, Copy)]
pub enum IndexParserState {
    Start,
    SpectrumIndexList,
    ChromatogramIndexList,
    Done,
}

impl XMLParseBase for IndexedMzMLIndexExtractor {}

impl IndexedMzMLIndexExtractor {
    pub fn new() -> IndexedMzMLIndexExtractor {
        IndexedMzMLIndexExtractor {
            spectrum_index: OffsetIndex::new("spectrum".into()),
            chromatogram_index: OffsetIndex::new("chromatogram".into()),
            last_id: String::new(),
        }
    }

    pub fn find_offset_from_reader<R: SeekRead>(&self, reader: &mut R) -> io::Result<Option<u64>> {
        reader.seek(SeekFrom::End(-200))?;
        let mut buf = Bytes::new();
        reader.read_to_end(&mut buf)?;
        let pattern = regex::Regex::new("<indexListOffset>(\\d+)</indexListOffset>").unwrap();
        if let Some(captures) = pattern.captures(&String::from_utf8_lossy(&buf)) {
            if let Some(offset) = captures.get(1) {
                if let Ok(offset) = offset.as_str().parse::<u64>() {
                    return Ok(Some(offset));
                }
            }
        }
        Ok(None)
    }

    pub fn start_element(
        &mut self,
        event: &BytesStart,
        state: IndexParserState,
    ) -> Result<IndexParserState, XMLError> {
        let elt_name = event.name();
        match elt_name.as_ref() {
            b"offset" => {
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key.as_ref() == b"idRef" {
                                self.last_id = attr
                                    .unescape_value()
                                    .expect("Error decoding idRef")
                                    .to_string();
                            }
                        }
                        Err(err) => {
                            return Err(err.into());
                        }
                    }
                }
            }
            b"index" => {
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key.as_ref() == b"name" {
                                let index_name = attr
                                    .unescape_value()
                                    .expect("Error decoding idRef")
                                    .to_string();
                                match index_name.as_ref() {
                                    "spectrum" => return Ok(IndexParserState::SpectrumIndexList),
                                    "chromatogram" => {
                                        return Ok(IndexParserState::ChromatogramIndexList)
                                    }
                                    _ => {}
                                }
                            }
                        }
                        Err(err) => {
                            return Err(err.into());
                        }
                    }
                }
            }
            b"indexList" => {}
            _ => {}
        }

        Ok(state)
    }

    pub fn end_element(
        &mut self,
        event: &BytesEnd,
        state: IndexParserState,
    ) -> Result<IndexParserState, XMLError> {
        let elt_name = event.name();
        match elt_name.as_ref() {
            b"offset" => {}
            b"index" => {}
            b"indexList" => return Ok(IndexParserState::Done),
            _ => {}
        }
        Ok(state)
    }

    pub fn text(
        &mut self,
        event: &BytesText,
        state: IndexParserState,
    ) -> Result<IndexParserState, XMLError> {
        match state {
            IndexParserState::SpectrumIndexList => {
                let bin = event
                    .unescape()
                    .expect("Failed to unescape spectrum offset");
                if let Ok(offset) = bin.parse::<u64>() {
                    if !self.last_id.is_empty() {
                        let key = mem::take(&mut self.last_id);
                        self.spectrum_index.insert(key, offset);
                    } else {
                        warn!("Out of order text in index")
                    }
                }
            }
            IndexParserState::ChromatogramIndexList => {
                let bin = event
                    .unescape()
                    .expect("Failed to unescape chromatogram offset");
                if let Ok(offset) = bin.parse::<u64>() {
                    if !self.last_id.is_empty() {
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


/**A SAX-style parser for building up the metadata section prior to the `<run>` element
of an mzML file.*/
#[derive(Debug, Default)]
pub struct FileMetadataBuilder<'a> {
    pub file_description: FileDescription,
    pub instrument_configurations: Vec<InstrumentConfiguration>,
    pub softwares: Vec<Software>,
    pub data_processings: Vec<DataProcessing>,
    pub reference_param_groups: HashMap<String, Vec<Param>>,
    pub last_group: String,
    pub(crate) instrument_id_map: Option<&'a mut IncrementingIdMap>,
}

impl<'a> XMLParseBase for FileMetadataBuilder<'a> {}
impl<'a> CVParamParse for FileMetadataBuilder<'a> {}

impl<'a> FileMetadataBuilder<'a> {
    pub fn start_element(&mut self, event: &BytesStart, state: MzMLParserState) -> ParserResult {
        let elt_name = event.name();
        match elt_name.as_ref() {
            b"fileDescription" => return Ok(MzMLParserState::FileDescription),
            b"fileContent" => return Ok(MzMLParserState::FileContents),
            b"sourceFileList" => return Ok(MzMLParserState::SourceFileList),
            b"sourceFile" => {
                let mut source_file = SourceFile::default();
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key.as_ref() == b"id" {
                                source_file.id = attr
                                    .unescape_value()
                                    .expect("Error decoding id")
                                    .to_string();
                            } else if attr.key.as_ref() == b"name" {
                                source_file.name = attr
                                    .unescape_value()
                                    .expect("Error decoding name")
                                    .to_string();
                            } else if attr.key.as_ref() == b"location" {
                                source_file.location = attr
                                    .unescape_value()
                                    .expect("Error decoding location")
                                    .to_string();
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg.into(), state));
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
                            if attr.key.as_ref() == b"id" {
                                software.id = attr
                                    .unescape_value()
                                    .expect("Error decoding id")
                                    .to_string();
                            } else if attr.key.as_ref() == b"version" {
                                software.version = attr
                                    .unescape_value()
                                    .expect("Error decoding version")
                                    .to_string();
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg.into(), state));
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
                            if attr.key.as_ref() == b"id" {
                                let key = attr
                                    .unescape_value()
                                    .expect("Error decoding id")
                                    .to_string();
                                self.reference_param_groups.entry(key.clone()).or_default();
                                self.last_group = key;
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg.into(), state));
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
                            if attr.key.as_ref() == b"id" {
                                ic.id = self
                                    .instrument_id_map
                                    .as_mut()
                                    .expect("An instrument ID map was not provided")
                                    .get(
                                        attr.unescape_value().expect("Error decoding id").as_ref(),
                                    );
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg.into(), state));
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
                            if attr.key.as_ref() == b"order" {
                                source.order = attr
                                    .unescape_value()
                                    .expect("Error decoding order")
                                    .parse()
                                    .expect("Failed to parse integer from `order`");
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg.into(), state));
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
                            if attr.key.as_ref() == b"order" {
                                analyzer.order = attr
                                    .unescape_value()
                                    .expect("Error decoding order")
                                    .parse()
                                    .expect("Failed to parse integer from `order`");
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg.into(), state));
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
                            if attr.key.as_ref() == b"order" {
                                detector.order = attr
                                    .unescape_value()
                                    .expect("Error decoding order")
                                    .parse()
                                    .expect("Failed to parse integer from `order`");
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg.into(), state));
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
                            if attr.key.as_ref() == b"id" {
                                dp.id = attr
                                    .unescape_value()
                                    .expect("Error decoding id")
                                    .to_string();
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg.into(), state));
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
                            if attr.key.as_ref() == b"order" {
                                method.order = attr
                                    .unescape_value()
                                    .expect("Error decoding order")
                                    .parse()
                                    .expect("Failed to parse order");
                            } else if attr.key.as_ref() == b"softwareRef" {
                                method.software_reference = attr
                                    .unescape_value()
                                    .expect("Error decoding softwareRef")
                                    .to_string();
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg.into(), state));
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

    pub fn empty_element(
        &mut self,
        event: &BytesStart,
        state: MzMLParserState,
        reader_position: usize,
    ) -> ParserResult {
        let elt_name = event.name();
        match elt_name.as_ref() {
            b"cvParam" | b"userParam" => match Self::handle_param(event, reader_position, state) {
                Ok(param) => {
                    self.fill_param_into(param, state);
                    return Ok(state);
                }
                Err(err) => return Err(err),
            },
            b"softwareRef" => if state == MzMLParserState::InstrumentConfiguration {
                    let ic = self.instrument_configurations.last_mut().unwrap();
                    for attr_parsed in event.attributes() {
                        match attr_parsed {
                            Ok(attr) => {
                                if attr.key.as_ref() == b"ref" {
                                    ic.software_reference = attr
                                        .unescape_value()
                                        .expect("Error decoding software reference")
                                        .to_string();
                                }
                            }
                            Err(msg) => {
                                return Err(self.handle_xml_error(msg.into(), state));
                            }
                        }
                    }
                },
            b"referenceableParamGroupRef" => {
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key.as_ref() == b"ref" {
                                let group_id = attr
                                    .unescape_value()
                                    .expect("Error decoding reference group")
                                    .to_string();

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
                            return Err(self.handle_xml_error(msg.into(), state));
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
        match elt_name.as_ref() {
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
pub struct IncrementingIdMap {
    id_map: HashMap<String, u32>,
    next_id: u32,
}

impl IncrementingIdMap {
    pub fn get(&mut self, key: &str) -> u32 {
        if let Some(value) = self.id_map.get(key) {
            *value
        } else {
            let value = self.next_id;
            self.id_map.insert(key.to_string(), self.next_id);
            self.next_id += 1;
            value
        }
    }
}
