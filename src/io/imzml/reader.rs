use std::{
    fs,
    io::{self, BufReader},
    mem,
    path::{Path, PathBuf},
};

use log::warn;
use quick_xml::{events::Event, se, Reader};
use thiserror::Error;
use uuid::Uuid;

use crate::io::{
    mzml::{
        CVParamParse, EntryType, FileMetadataBuilder, MzMLParserError, MzMLParserState, MzMLReader,
        MzMLSAX, MzMLSpectrumBuilder, ParserResult, SpectrumBuilding, IncrementingIdMap
    },
    utils::DetailLevel,
    SpectrumSource,
};
use mzpeaks::prelude::*;

use crate::params::{ControlledVocabulary, ParamValue};
use crate::prelude::*;
use crate::spectrum::bindata::{ArrayRetrievalError, BinaryDataArrayType, BuildFromArrayMap};

#[derive(Debug, Default)]
pub struct ImzMLFileMetadata {
    pub uuid: Option<Uuid>,
    pub data_mode: Option<IbdDataMode>,
    pub ibd_checksum: Option<String>,
    pub ibd_checksum_type: Option<String>,
}

use super::{
    ibd::{IbdDataMode, IbdFile},
    IbdError,
};

use crate::spectrum::spectrum_types::MultiLayerSpectrum;

#[derive(Debug, Default, Clone)]
pub struct DataRangeQuery {
    pub name: String,
    pub offset: usize,
    pub length: usize,
    pub datatype: BinaryDataArrayType,
}

#[derive(Debug, Error)]
pub enum ImzMLError {
    #[error("An IBD-related error occurred: {0}")]
    IbdError(#[from] IbdError),
    #[error("An mzML-related error occurred: {0}")]
    MzMLError(#[from] MzMLParserError),
    #[error("An error occurred while decoding binary data: {0}")]
    ArrayRetrievalError(#[from] ArrayRetrievalError),
}

impl From<ImzMLError> for io::Error {
    fn from(value: ImzMLError) -> Self {
        Self::new(io::ErrorKind::Other, Box::new(value))
    }
}

/// Check if the buffer contains an imzML file by looking for the IMS controlled vocabulary
/// There isn't AFAIK a formal mechanism to identify imzML files other than the presence of
/// the IMS controlled vocabulary (CV) in the cvList section of the mzML file.
pub fn is_imzml(buffer: &[u8]) -> bool {
    let mut bufread = BufReader::new(io::Cursor::new(buffer));
    let mut reader = Reader::from_reader(&mut bufread);
    let mut buf = Vec::new();
    let mut in_cv_list = false;

    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Start(ref e)) => {
                match e.name().as_ref() {
                    b"cvList" => {
                        in_cv_list = true;
                    }
                    b"cv" if in_cv_list => {
                        // Check if this cv element has id="IMS"
                        for attr in e.attributes() {
                            if let Ok(attr) = attr {
                                if attr.key.as_ref() == b"id" {
                                    if let Ok(value) = attr.unescape_value() {
                                        if value.as_ref() == "IMS" {
                                            return true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    _ => {}
                }
            }
            Ok(Event::End(ref e)) => {
                if e.name().as_ref() == b"cvList" {
                    in_cv_list = false;
                    // If we've finished processing cvList and haven't found IMS, it's not imzML
                    return false;
                }
            }
            Ok(Event::Eof) => return false,
            Ok(_) => {}
            Err(_) => return false,
        }
        buf.clear();
    }
}

pub struct ImzMLSpectrumBuilder<
    'a,
    C: CentroidLike + BuildFromArrayMap,
    D: DeconvolutedCentroidLike + BuildFromArrayMap,
> {
    inner: MzMLSpectrumBuilder<'a, C, D>,
    data_registry: Option<&'a mut IbdFile>,
    current_data_range_query: DataRangeQuery,
    metadata: ImzMLFileMetadata,
}

impl<C: CentroidLike + BuildFromArrayMap, D: DeconvolutedCentroidLike + BuildFromArrayMap> Default
    for ImzMLSpectrumBuilder<'_, C, D>
{
    fn default() -> Self {
        Self {
            inner: Default::default(),
            data_registry: None,
            current_data_range_query: DataRangeQuery::default(),
            metadata: ImzMLFileMetadata::default(),
        }
    }
}

/**A SAX-style parser for building up the metadata section prior to the `<run>` element
of an mzML file.*/
#[derive(Debug, Default)]
pub struct ImzmlMetadataBuilder<'a> {
    pub mzml_metadata_builder: FileMetadataBuilder<'a>,
    pub imzml_metadata: ImzMLFileMetadata,
}



impl XMLParseBase for FileMetadataBuilder<'_> {}
impl CVParamParse for FileMetadataBuilder<'_> {}

impl ImzmlMetadataBuilder<'_> {
    pub fn start_element(&mut self, event: &BytesStart, state: MzMLParserState) -> ParserResult {
        // Always delegate to the underlying builder
        self.mzml_metadata_builder.start_element(event, state)
    }

    pub fn empty_element(
        &mut self,
        event: &BytesStart,
        state: MzMLParserState,
        reader_position: usize,
    ) -> ParserResult {
        let elt_name = event.name();
        match elt_name.as_ref() {
            b"cvParam" | b"userParam" => {
                match &state {
                    MzMLParserState::FileContents => {
                        // Try to handle as IMS parameter first
                        match MzMLSpectrumBuilder::<C, D>::handle_param(event, reader_position, state) {
                            Ok(param) => {
                                if param.is_controlled() 
                                    && param.controlled_vocabulary == Some(ControlledVocabulary::IMS) 
                                {
                                    // Handle IMS-specific parameters
                                    match param.accession.unwrap() {
                                        1000080 => {
                                            // Handle UUID...
                                            let uuid_str = param.value.as_str().trim_matches(|c| c == '{' || c == '}');
                                            match Uuid::parse_str(uuid_str) {
                                                Ok(uuid) => self.imzml_metadata.uuid = Some(uuid),
                                                Err(e) => warn!("Failed to parse UUID '{}': {}", uuid_str, e),
                                            }
                                        }
                                        1000031 => {
                                            self.imzml_metadata.data_mode = Some(IbdDataMode::Processed);
                                        }
                                        1000030 => {
                                            self.imzml_metadata.data_mode = Some(IbdDataMode::Continuous);
                                        }
                                        1000090 => {
                                            self.imzml_metadata.ibd_checksum = Some(param.value.to_string());
                                            self.imzml_metadata.ibd_checksum_type = Some("MD5".to_string());
                                        }
                                        _ => {
                                            // Other IMS parameters - delegate to base builder
                                            return self.mzml_metadata_builder.empty_element(event, state, reader_position);
                                        }
                                    }
                                    return Ok(state);
                                } else {
                                    // Non-IMS parameters - delegate to base builder
                                    return self.mzml_metadata_builder.empty_element(event, state, reader_position);
                                }
                            }
                            Err(err) => return Err(err),
                        }
                    }
                    _ => {
                        return self.mzml_metadata_builder.empty_element(event, state, reader_position);
                    }
                }
            }
            _ => {
                return self.mzml_metadata_builder.empty_element(event, state, reader_position);
            }
        }
    }

    // Add other missing methods for complete passthrough
    pub fn end_element(&mut self, event: &BytesEnd, state: MzMLParserState) -> ParserResult {
        self.mzml_metadata_builder.end_element(event, state)
    }

    pub fn text(&mut self, event: &BytesText, state: MzMLParserState) -> ParserResult {
        self.mzml_metadata_builder.text(event, state)
    }
}

impl<C: CentroidLike + BuildFromArrayMap, D: DeconvolutedCentroidLike + BuildFromArrayMap>
    ImzMLSpectrumBuilder<'_, C, D>
{
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_detail_level(detail_level: DetailLevel) -> Self {
        Self {
            inner: MzMLSpectrumBuilder::with_detail_level(detail_level),
            data_registry: None,
            current_data_range_query: DataRangeQuery::default(),
            metadata: ImzMLFileMetadata::default(),
        }
    }

    pub fn is_chromatogram_entry(&self) -> bool {
        self.inner.is_chromatogram_entry()
    }

    pub fn is_spectrum_entry(&self) -> bool {
        self.inner.is_spectrum_entry()
    }

    pub fn entry_type(&self) -> EntryType {
        self.inner.entry_type()
    }

    pub fn set_entry_type(&mut self, entry_type: EntryType) {
        self.inner.set_entry_type(entry_type)
    }
}

impl<'a, C: CentroidLike + BuildFromArrayMap, D: DeconvolutedCentroidLike + BuildFromArrayMap>
    MzMLSAX for ImzMLSpectrumBuilder<'a, C, D>
{
    fn start_element(
        &mut self,
        event: &quick_xml::events::BytesStart,
        state: MzMLParserState,
    ) -> crate::io::mzml::ParserResult {
        self.inner.start_element(event, state)
    }

    fn empty_element(
        &mut self,
        event: &quick_xml::events::BytesStart,
        state: MzMLParserState,
        reader_position: usize,
    ) -> ParserResult {
        let elt_name = event.name();
        match elt_name.as_ref() {
            b"cvParam" | b"userParam" => {
                match &state {
                    MzMLParserState::BinaryDataArray => {
                        match MzMLSpectrumBuilder::<'a, C, D>::handle_param(
                            event,
                            reader_position,
                            state,
                        ) {
                            Ok(param) => {
                                match &state {
                                    MzMLParserState::BinaryDataArray => {
                                        // Check if this is an IMS namespace parameter
                                        if param.is_controlled()
                                            && param.controlled_vocabulary
                                                == Some(ControlledVocabulary::IMS)
                                        {
                                            match param.accession.unwrap() {
                                                // IMS:1000102 - external offset
                                                1000102 => {
                                                    self.current_data_range_query.offset = param
                                                        .value
                                                        .to_u64()
                                                        .expect("Failed to extract external offset")
                                                        as usize
                                                }
                                                // IMS:1000103 - external array length
                                                1000103 => self.current_data_range_query.length =
                                                    param.value.to_u64().expect(
                                                        "Failed to extract external array length",
                                                    )
                                                        as usize,
                                                // IMS:1000104 - external encoded length (not used for IBD)
                                                1000104 => {
                                                    // This is for encoded length, we can ignore it for IBD files
                                                    // as we use the raw array length (IMS:1000103)
                                                    // TODO maybe store this to sanity check that the raw length
                                                    // divided by the encoded length is matches the size of the data type?
                                                }
                                                _ => self.inner.fill_param_into(param, state),
                                            }
                                        } else {
                                            // Delegate all non-IMS parameters to the inner mzML parser
                                            self.inner.fill_param_into(param, state)
                                        }
                                    }
                                    _ => self.inner.fill_param_into(param, state),
                                }
                                return Ok(state);
                            }
                            Err(err) => return Err(err),
                        }
                    }
                    _ => return self.inner.empty_element(event, state, reader_position),
                }
            }
            _ => {}
        }
        Ok(state)
    }

    fn end_element(
        &mut self,
        event: &quick_xml::events::BytesEnd,
        state: MzMLParserState,
    ) -> ParserResult {
        let elt_name = event.name();
        let res = match elt_name.as_ref() {
            b"binaryDataArray" => {
                // For imzML/IBD, we don't need a dataset name like mzMLb/HDF5
                // We just need offset and length to read from the IBD file
                if self.current_data_range_query.offset == 0
                    && self.current_data_range_query.length == 0
                {
                    return Err(MzMLParserError::IncompleteElementError(
                        "The external data offset and length were missing".to_owned(),
                        MzMLParserState::BinaryDataArray,
                    ));
                }
                let detail_level = self.inner.detail_level;
                let data_request = mem::take(&mut self.current_data_range_query);
                let array = self.inner.current_array_mut();
                if !matches!(detail_level, DetailLevel::MetadataOnly) {
                    self.data_registry
                        .as_mut()
                        .expect("Did not provide data registry")
                        .get(&data_request, array)
                } else {
                    Ok(())
                }
            }
            _ => Ok(()),
        };
        match res {
            Ok(()) => {}
            Err(e) => Err(MzMLParserError::IOError(state, io::Error::from(e)))?,
        }
        self.inner.end_element(event, state)
    }

    fn text(
        &mut self,
        event: &quick_xml::events::BytesText,
        state: MzMLParserState,
    ) -> ParserResult {
        self.inner.text(event, state)
    }
}


/*
An mzML parser that supports iteration and random access. The parser produces
[`Spectrum`] instances, which may be converted to [`RawSpectrum`](crate::spectrum::RawSpectrum)
or [`CentroidSpectrum`](crate::spectrum::CentroidSpectrum) as is appropriate to the data.

When the readable stream the parser is wrapped around supports [`io::Seek`],
additional random access operations are available.
*/
pub struct ImzMLReaderType<
    R: Read,
    C: CentroidLike = CentroidPeak,
    D: DeconvolutedCentroidLike = DeconvolutedPeak,
> {
    /// The state the parser was in last.
    pub state: MzMLParserState,
    /// The raw reader
    handle: BufReader<R>,
    /// A place to store the last error the parser encountered
    error: Option<Box<MzMLParserError>>,
    /// A spectrum ID to byte offset for fast random access
    pub spectrum_index: OffsetIndex,
    pub chromatogram_index: Box<OffsetIndex>,
    /// The description of the file's contents and the previous data files that were
    /// consumed to produce it.
    pub(crate) file_description: FileDescription,
    /// A mapping of different instrument configurations (source, analyzer, detector) components
    /// by ID string.
    pub(crate) instrument_configurations: HashMap<u32, InstrumentConfiguration>,
    /// The different software components that were involved in the processing and creation of this
    /// file.
    pub(crate) softwares: Vec<Software>,
    pub(crate) samples: Vec<Sample>,
    /// The data processing and signal transformation operations performed on the raw data in previous
    /// source files to produce this file's contents.
    pub(crate) data_processings: Vec<DataProcessing>,
    pub(crate) scan_settings: Vec<ScanSettings>,
    /// A cache of repeated paramters
    pub reference_param_groups: HashMap<String, Vec<Param>>,
    pub detail_level: DetailLevel,

    // SpectrumList attributes
    pub run: MassSpectrometryRun,
    num_spectra: Option<u64>,

    buffer: Bytes,
    instrument_id_map: Box<IncrementingIdMap>,

    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,

    // ✅ ADD THESE FIELDS:
    pub imzml_metadata: ImzMLFileMetadata,
    pub ibd_file: Option<IbdFile>,  // This will be opened later when needed
}

impl<
        R: Read + Seek,
        C: CentroidLike + BuildFromArrayMap,
        D: DeconvolutedCentroidLike + BuildFromArrayMap,
    > IntoIonMobilityFrameSource<C, D> for ImzMLReaderType<R, C, D>
{
    type IonMobilityFrameSource<
        CF: FeatureLike<mzpeaks::MZ, mzpeaks::IonMobility>,
        DF: FeatureLike<mzpeaks::Mass, mzpeaks::IonMobility> + KnownCharge,
    > = Generic3DIonMobilityFrameSource<C, D, Self, CF, DF>;

    fn try_into_frame_source<
        CF: FeatureLike<mzpeaks::MZ, mzpeaks::IonMobility>,
        DF: FeatureLike<mzpeaks::Mass, mzpeaks::IonMobility> + KnownCharge,
    >(
        mut self,
    ) -> Result<Self::IonMobilityFrameSource<CF, DF>, crate::io::IntoIonMobilityFrameSourceError>
    {
        if let Some(state) = self.has_ion_mobility() {
            if matches!(state, HasIonMobility::Dimension) {
                Ok(Self::IonMobilityFrameSource::new(self))
            } else {
                Err(crate::io::IntoIonMobilityFrameSourceError::ConversionNotPossible)
            }
        } else {
            Err(crate::io::IntoIonMobilityFrameSourceError::NoIonMobilityFramesFound)
        }
    }
}

impl<
        'a,
        'b: 'a,
        R: Read,
        C: CentroidLike + BuildFromArrayMap,
        D: DeconvolutedCentroidLike + BuildFromArrayMap,
    > ImzMLReaderType<R, C, D>
{
    /// Create a new [`ImzMLReaderType`] instance, wrapping the [`io::Read`] handle
    /// provided with an [`io::BufReader`] and parses the metadata section of the file.
    pub fn new(file: R) -> ImzMLReaderType<R, C, D> {
        Self::with_buffer_capacity_and_detail_level(file, BUFFER_SIZE, DetailLevel::Full)
    }

    pub fn with_buffer_capacity_and_detail_level(
        file: R,
        capacity: usize,
        detail_level: DetailLevel,
    ) -> ImzMLReaderType<R, C, D> {
        let handle = BufReader::with_capacity(capacity, file);
        let mut inst = ImzMLReaderType {
            handle,
            state: MzMLParserState::Start,
            error: None,
            buffer: Bytes::new(),
            spectrum_index: OffsetIndex::new("spectrum".to_owned()),
            chromatogram_index: Box::new(OffsetIndex::new("chromatogram".to_owned())),

            file_description: FileDescription::default(),
            scan_settings: Default::default(),
            instrument_configurations: HashMap::new(),
            softwares: Vec::new(),
            samples: Vec::new(),
            data_processings: Vec::new(),
            reference_param_groups: HashMap::new(),
            detail_level,

            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
            instrument_id_map: Box::new(IncrementingIdMap::default()),
            num_spectra: None,
            run: MassSpectrometryRun::default(),

            imzml_metadata: ImzMLFileMetadata::default(),
            ibd_file: None,
        };
        match inst.parse_metadata() {
            Ok(()) => {}
            Err(_err) => {}
        }
        inst
    }

    /**Parse the metadata section of the file using [`FileMetadataBuilder`]
     */
    fn parse_metadata(&mut self) -> Result<(), MzMLParserError> {
        let mut reader = Reader::from_reader(&mut self.handle);
        reader.trim_text(true);
        let mut accumulator: ImzmlMetadataBuilder<'_> = ImzmlMetadataBuilder {
            mzml_metadata_builder: FileMetadataBuilder {
                instrument_id_map: Some(&mut self.instrument_id_map),
                ..Default::default()
            },
            imzml_metadata: ImzMLFileMetadata::default(),
        };
        loop {
            match reader.read_event_into(&mut self.buffer) {
                Ok(Event::Start(ref e)) => {
                    match accumulator.start_element(e, self.state) {
                        Ok(state) => {
                            self.state = state;
                            match &self.state {
                                MzMLParserState::SpectrumList | MzMLParserState::Spectrum => break,
                                MzMLParserState::ParserError => {
                                    log::error!(
                                        "Encountered an error while starting {:?}",
                                        String::from_utf8_lossy(&self.buffer)
                                    );
                                }
                                _ => {}
                            }
                        }
                        Err(message) => {
                            self.state = MzMLParserState::ParserError;
                            self.error = Some(Box::new(message));
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
                            self.error = Some(Box::new(message));
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
                            self.error = Some(Box::new(message));
                        }
                    };
                }
                Ok(Event::Empty(ref e)) => {
                    match accumulator.empty_element(e, self.state, reader.buffer_position()) {
                        Ok(state) => {
                            self.state = state;
                        }
                        Err(message) => {
                            self.state = MzMLParserState::ParserError;
                            self.error = Some(Box::new(message));
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
                            self.error = Some(Box::new(MzMLParserError::IncompleteElementError(
                                String::from_utf8_lossy(&self.buffer).to_string(),
                                self.state,
                            )));
                            self.state = MzMLParserState::ParserError;
                        }
                    }
                    _ => {
                        self.error = Some(Box::new(MzMLParserError::IncompleteElementError(
                            String::from_utf8_lossy(&self.buffer).to_string(),
                            self.state,
                        )));
                        self.state = MzMLParserState::ParserError;
                    }
                },
                _ => {}
            };
            self.buffer.clear();
            match self.state {
                MzMLParserState::SpectrumList | MzMLParserState::ParserError => {
                    break;
                }
                _ => {}
            };
        }
        // Extract standard mzML metadata (already present)
        self.file_description = accumulator.mzml_metadata_builder.file_description;
        self.instrument_configurations = accumulator
            .mzml_metadata_builder
            .instrument_configurations
            .into_iter()
            .map(|ic| (ic.id, ic))
            .collect();
        self.softwares = accumulator.mzml_metadata_builder.softwares;
        self.samples = accumulator.mzml_metadata_builder.samples;
        self.data_processings = accumulator.mzml_metadata_builder.data_processings;
        self.reference_param_groups = accumulator.mzml_metadata_builder.reference_param_groups;
        self.scan_settings = accumulator.mzml_metadata_builder.scan_settings;
        self.run.id = accumulator.mzml_metadata_builder.run_id;
        self.run.default_instrument_id = accumulator.mzml_metadata_builder.default_instrument_config;
        self.run.default_source_file_id = accumulator.mzml_metadata_builder.default_source_file;
        self.run.start_time = accumulator.mzml_metadata_builder.start_timestamp;
        self.run.default_data_processing_id = accumulator.mzml_metadata_builder.default_data_processing;
        self.num_spectra = accumulator.mzml_metadata_builder.num_spectra;


        let imzml_metadata = accumulator.imzml_metadata;

        // Validate required metadata
        let data_mode = imzml_metadata.data_mode.ok_or_else(|| {
            MzMLParserError::IncompleteElementError(
                "Missing required imzML data mode (IMS:1000030 or IMS:1000031)".to_string(),
                MzMLParserState::FileContents,
            )
        })?;
        
        let uuid = imzml_metadata.uuid.ok_or_else(|| {
            MzMLParserError::IncompleteElementError(
                "Missing required imzML UUID (IMS:1000080)".to_string(),
                MzMLParserState::FileContents,
            )
        })?;
        
        let checksum = imzml_metadata.ibd_checksum.ok_or_else(|| {
            MzMLParserError::IncompleteElementError(
                "Missing required imzML IBD checksum".to_string(),
                MzMLParserState::FileContents,
            )
        })?;
        
        // TODO: Store these in the ImzMLReaderType struct for later use
        log::debug!("Parsed imzML metadata - Mode: {:?}, UUID: {}, Checksum: {}", 
                   data_mode, uuid, checksum);

        match self.state {
            MzMLParserState::SpectrumDone | MzMLParserState::ChromatogramDone => Ok(()),
            MzMLParserState::ParserError => {
                Err(*self
                    .error
                    .take()
                    .unwrap_or(Box::new(MzMLParserError::UnknownError(
                        MzMLParserState::ParserError,
                    ))))
            }
            _ => Err(MzMLParserError::IncompleteSpectrum),
        }
    }

    pub(crate) fn _parse_into<
        B: MzMLSAX + SpectrumBuilding<'a, C, D, MultiLayerSpectrum<C, D>> + 'a,
    >(
        &'b mut self,
        mut accumulator: B,
    ) -> Result<(B, usize), MzMLParserError> {
        if self.state == MzMLParserState::EOF {
            return Err(MzMLParserError::SectionOver("spectrum"));
        }

        let mut reader = Reader::from_reader(&mut self.handle);
        reader.trim_text(true);
        accumulator = accumulator.borrow_instrument_configuration(&mut self.instrument_id_map);
        accumulator.set_run_data_processing(self.run.default_data_processing_id.clone().map(|v| v.into_boxed_str()));
        let mut offset: usize = 0;

        macro_rules! err_state {
            ($message:ident) => {{
                self.state = MzMLParserState::ParserError;
                self.error = Some(Box::new($message));
            }};
        }

        loop {
            match reader.read_event_into(&mut self.buffer) {
                Ok(Event::Start(ref e)) => {
                    if log::log_enabled!(log::Level::Trace) {
                        log::trace!(
                            "Starting mzML element: {}",
                            String::from_utf8_lossy(e.name().as_ref())
                        );
                    }
                    match accumulator.start_element(e, self.state) {
                        Ok(state) => {
                            self.state = state;
                            if state == MzMLParserState::ParserError {
                                warn!(
                                    "Encountered an error while starting {:?}",
                                    String::from_utf8_lossy(&self.buffer)
                                );
                            }
                        }
                        Err(message) => err_state!(message),
                    };
                }
                Ok(Event::End(ref e)) => {
                    log::trace!(
                        "Ending mzML element: {}",
                        String::from_utf8_lossy(e.name().as_ref())
                    );
                    match accumulator.end_element(e, self.state) {
                        Ok(state) => {
                            self.state = state;
                        }
                        Err(message) => err_state!(message),
                    };
                }
                Ok(Event::Text(ref e)) => {
                    match accumulator.text(e, self.state) {
                        Ok(state) => {
                            self.state = state;
                        }
                        Err(message) => err_state!(message),
                    };
                }
                Ok(Event::Empty(ref e)) => {
                    match accumulator.empty_element(e, self.state, reader.buffer_position()) {
                        Ok(state) => {
                            self.state = state;
                        }
                        Err(message) => err_state!(message),
                    }
                }
                Ok(Event::Eof) => {
                    log::trace!("Reached EOF");
                    self.state = MzMLParserState::EOF;
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
                            self.error = Some(Box::new(MzMLParserError::IncompleteElementError(
                                String::from_utf8_lossy(&self.buffer).to_string(),
                                self.state,
                            )));
                            self.state = MzMLParserState::ParserError;
                            log::trace!("Expected element {expected}, found {_found}");
                        }
                    }
                    e => {
                        self.error = Some(Box::new(MzMLParserError::IncompleteElementError(
                            e.to_string(),
                            self.state,
                        )));
                        self.state = MzMLParserState::ParserError;
                    }
                },
                _ => {}
            };
            offset += self.buffer.len();
            self.buffer.clear();
            match self.state {
                MzMLParserState::SpectrumDone
                | MzMLParserState::ChromatogramDone
                | MzMLParserState::ParserError => {
                    break;
                }
                _ => {}
            };
        }
        match self.state {
            MzMLParserState::SpectrumDone | MzMLParserState::ChromatogramDone => {
                Ok((accumulator, offset))
            }
            MzMLParserState::ParserError if self.error.is_some() => {
                Err(*self.error.take().unwrap())
            }
            MzMLParserState::ParserError if self.error.is_none() => {
                warn!(
                    "Terminated with ParserError but no error set: {:?}",
                    self.error
                );
                Ok((accumulator, offset))
            }
            MzMLParserState::EOF => {
                Err(MzMLParserError::EOF)
            }
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
        let accumulator = MzMLSpectrumBuilder::<C, D>::with_detail_level(self.detail_level);
        match self.state {
            MzMLParserState::SpectrumDone => {
                self.state = MzMLParserState::Resume;
            }
            MzMLParserState::ParserError => {
                log::error!("Starting parsing from error: {:?}", self.error);
            }
            state if state > MzMLParserState::SpectrumDone => {
                log::error!(
                    "Attempting to start parsing a spectrum in state {}",
                    self.state
                );
            }
            _ => {}
        }
        match self._parse_into(accumulator) {
            Ok((accumulator, sz)) => {
                accumulator.into_spectrum(spectrum);
                if self.detail_level == DetailLevel::Full {
                    if let Err(e) = spectrum.try_build_peaks() {
                        log::debug!("Failed to eagerly load peaks from centroid spectrum: {e}");
                    }
                }
                Ok(sz)
            }
            Err(err) => {
                match &err {
                    MzMLParserError::EOF => {},
                    err => log::error!("Error while reading mzML spectrum: {err}")
                };
                Err(err)
            },
        }
    }

    /// Read the next spectrum directly. Used to implement iteration.
    pub fn read_next(&mut self) -> Option<MultiLayerSpectrum<C, D>> {
        if self.state == MzMLParserState::EOF {
            return None;
        }
        let mut spectrum = MultiLayerSpectrum::<C, D>::default();
        match self.read_into(&mut spectrum) {
            Ok(_sz) => Some(spectrum),
            Err(err) => {
                match err {
                    MzMLParserError::EOF => {},
                    err => {
                        trace!("Failed to read next spectrum: {err}");
                    }
                }
                None
            }
        }
    }

    fn _read_next_chromatogram(&mut self) -> Result<Chromatogram, MzMLParserError> {
        let accumulator = MzMLSpectrumBuilder::<C, D>::with_detail_level(self.detail_level);

        match self.state {
            MzMLParserState::ChromatogramDone => {
                self.state = MzMLParserState::Resume;
            }
            MzMLParserState::ParserError => {
                warn!("Starting parsing from error: {:?}", self.error);
            }
            state
                if state > MzMLParserState::ChromatogramDone
                    && state < MzMLParserState::Chromatogram =>
            {
                warn!(
                    "Attempting to start parsing a spectrum in state {}",
                    self.state
                );
            }
            _ => {}
        }
        match self._parse_into(accumulator) {
            Ok((accumulator, _sz)) => {
                if accumulator.is_chromatogram_entry() {
                    let mut chrom = Chromatogram::default();
                    accumulator.into_chromatogram(&mut chrom);
                    Ok(chrom)
                } else {
                    Err(MzMLParserError::UnknownError(self.state))
                }
            }
            Err(err) => {
                log::error!("Error while reading mzML chromatogram: {err}");
                Err(err)
            },
        }
    }
}