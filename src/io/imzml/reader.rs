use std::{
    fs,
    io::{self, BufReader},
    mem,
    path::{Path, PathBuf},
};

use log::warn;
use quick_xml::{Reader, events::Event};
use thiserror::Error;
use uuid::Uuid;

use mzpeaks::prelude::*;
use crate::io::{
    mzml::{
        CVParamParse, EntryType, MzMLParserError, MzMLParserState, MzMLReader, MzMLSAX,
        MzMLSpectrumBuilder, ParserResult, SpectrumBuilding, IncrementingIdMap
    },
    utils::DetailLevel,
    SpectrumSource,
};

use crate::prelude::*;
use crate::params::{ControlledVocabulary, Param, ParamValue};
use crate::spectrum::{IsolationWindow, ScanWindow, SelectedIon, Precursor};
use crate::spectrum::bindata::{
    ArrayRetrievalError,
    BinaryDataArrayType, BuildFromArrayMap, DataArray
};

#[derive(Debug, Default)]
pub struct ImzMLFileMetadata {
    pub uuid: Option<Uuid>,
    pub data_mode: Option<IbdDataMode>,
    pub ibd_checksum: Option<String>,
    pub ibd_checksum_type: Option<String>,
    pub ibd_file_path: Option<String>,
}

use super::{ibd::{IbdFile, IbdDataMode}, IbdError};

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
                    MzMLParserState::FileContents => {
                        match MzMLSpectrumBuilder::<'a, C, D>::handle_param(
                            event,
                            reader_position,
                            state,
                        ) {
                            Ok(param) => {
                                // Check if this is an IMS parameter for file-level metadata
                                if param.is_controlled() && param.controlled_vocabulary == Some(ControlledVocabulary::IMS) {
                                    match param.accession.unwrap() {
                                        // IMS:1000080 - universally unique identifier
                                        1000080 => {
                                            // Parse UUID from parameter value
                                            match Uuid::parse_str(param.value.as_str().trim_matches(|c| c == '{' || c == '}')) {
                                                Ok(uuid) => {
                                                    self.metadata.uuid = Some(uuid);
                                                }
                                                Err(e) => {
                                                    warn!("Failed to parse UUID '{}': {}", param.value.as_str(), e);
                                                }
                                            }
                                        }
                                        // IMS:1000031 - processed mode
                                        1000031 => {
                                            if self.metadata.data_mode.is_none() {
                                                self.metadata.data_mode = Some(IbdDataMode::Processed);
                                            }
                                            else {
                                                warn!("Multiple data modes found");
                                            }
                                        }
                                        // IMS:1000030 - continuous mode (if exists)
                                        1000030 => {
                                            if self.metadata.data_mode.is_none() {
                                                self.metadata.data_mode = Some(IbdDataMode::Continuous);
                                            }
                                            else {
                                                warn!("Multiple data modes found");
                                            }
                                        }
                                        // IMS:1000090 - ibd MD5
                                        1000090 => {
                                            if self.metadata.ibd_checksum.is_none() {
                                                self.metadata.ibd_checksum = Some(param.value.to_string());
                                                self.metadata.ibd_checksum_type = Some("MD5".to_string());
                                            } else {
                                                warn!("Multiple IBD checksums found");
                                            }
                                        }
                                        // IMS:1000091 - ibd SHA-1
                                        1000091 => {
                                            if self.metadata.ibd_checksum.is_none() {
                                                self.metadata.ibd_checksum = Some(param.value.to_string());
                                                self.metadata.ibd_checksum_type = Some("SHA-1".to_string());
                                            } else {
                                                warn!("Multiple IBD checksums found");
                                            }
                                        }
                                        // IMS:1000092 - ibd SHA-256
                                        1000092 => {
                                            if self.metadata.ibd_checksum.is_none() {
                                                self.metadata.ibd_checksum = Some(param.value.to_string());
                                                self.metadata.ibd_checksum_type = Some("SHA-256".to_string());
                                            } else {
                                                warn!("Multiple IBD checksums found");
                                            }
                                        }
                                        // IMS:1000100 - ibd file path
                                        1000070 => {
                                            if self.metadata.ibd_file_path.is_none() {
                                                self.metadata.ibd_file_path = Some(param.value.to_string());
                                            } else {
                                                warn!("Multiple IBD file paths found");
                                            }
                                        }
                                        _ => {
                                            // Other IMS parameters, delegate to inner parser
                                            self.inner.fill_param_into(param, state);
                                        }
                                    }
                                } else {
                                    // Non-IMS parameters, delegate to inner parser
                                    self.inner.fill_param_into(param, state);
                                }
                                return Ok(state);
                            }
                            Err(err) => return Err(err),
                        }
                    }
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
                                        if param.is_controlled() && param.controlled_vocabulary == Some(ControlledVocabulary::IMS) {
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
                if self.current_data_range_query.offset == 0 && self.current_data_range_query.length == 0 {
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
            b"fileContent" => {
                // Validate that we have all required metadata
                if self.metadata.data_mode.is_none() {
                    return Err(MzMLParserError::IncompleteElementError(
                        "Data mode (IMS:1000030 or IMS:1000031) not found in fileContent".to_owned(),
                        state,
                    ));
                }
                
                if self.metadata.uuid.is_none() {
                    return Err(MzMLParserError::IncompleteElementError(
                        "UUID (IMS:1000080) not found in fileContent".to_owned(),
                        state,
                    ));
                }
                
                if self.metadata.ibd_checksum.is_none() || self.metadata.ibd_checksum_type.is_none() {
                    return Err(MzMLParserError::IncompleteElementError(
                        "IBD checksum (IMS:1000090, IMS:1000091, or IMS:1000092) not found in fileContent".to_owned(),
                        state,
                    ));
                }
                
                // If ibd_file_path is not set, derive it from the imzML path
                if self.metadata.ibd_file_path.is_none() {
                    // We need access to the imzML file path here
                    // This is a limitation of the current design - we don't have access to the file path
                    // in the spectrum builder. For now, we'll set it to None and handle it in the reader
                    warn!("IBD file path not specified in metadata, will derive from imzML path");
                }
                
                // TODO: Try to open IBD file here once we have access to the file path
                // For now, we'll defer this to the reader level
                
                Ok(())
            }

            _ => Ok(()),
        };
        match res {
            Ok(()) => {}
            Err(e) => Err(MzMLParserError::IOError(
                state,
                io::Error::from(e),
            ))?,
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


impl<'a, C: CentroidLike + BuildFromArrayMap, D: DeconvolutedCentroidLike + BuildFromArrayMap>
    SpectrumBuilding<'a, C, D, MultiLayerSpectrum<C, D>> for ImzMLSpectrumBuilder<'a, C, D>
{
    fn isolation_window_mut(&mut self) -> &mut IsolationWindow {
        self.inner.isolation_window_mut()
    }

    fn scan_window_mut(&mut self) -> &mut ScanWindow {
        self.inner.scan_window_mut()
    }

    fn selected_ion_mut(&mut self) -> &mut SelectedIon {
        self.inner.selected_ion_mut()
    }

    fn current_array_mut(&mut self) -> &mut DataArray {
        self.inner.current_array_mut()
    }

    fn into_spectrum(self, spectrum: &mut MultiLayerSpectrum<C, D>) {
        self.inner.into_spectrum(spectrum)
    }

    fn fill_spectrum<P: ParamLike + Into<Param> + ParamValue>(&mut self, param: P) {
        self.inner.fill_spectrum(param)
    }

    fn fill_binary_data_array<P: ParamLike + Into<Param> + ParamValue>(&mut self, param: P) {
        self.inner.fill_binary_data_array(param)
    }

    fn borrow_instrument_configuration(
        mut self,
        instrument_configurations: &'a mut IncrementingIdMap,
    ) -> Self {
        self.inner = self
            .inner
            .borrow_instrument_configuration(instrument_configurations);
        self
    }

    fn new_selected_ion(&mut self) -> &mut SelectedIon {
        self.inner.new_selected_ion()
    }

    fn into_chromatogram(self, chromatogram: &mut crate::spectrum::Chromatogram) {
        self.inner.into_chromatogram(chromatogram)
    }

    fn new_precursor_mut(&mut self) -> &mut Precursor {
        self.inner.new_precursor_mut()
    }

    fn precursor_mut(&mut self) -> &mut Precursor {
        self.inner.precursor_mut()
    }
}





/// imzML file reader that extends mzML functionality with external binary data support
///
/// This is a simple wrapper around [`MzMLReader`] that adds support for reading
/// binary data from an external `.ibd` file as specified in the imzML format.
pub struct ImzMLReader {
    /// The underlying mzML reader
    mzml_reader: MzMLReader<fs::File>,
    /// Handle to the IBD binary data file  
    ibd_file: Option<IbdFile>,
    /// The imzML file path
    imzml_path: PathBuf,
    /// Data storage mode (continuous vs processed)
    data_mode: IbdDataMode,
    /// File-level metadata extracted from imzML
    file_metadata: ImzMLFileMetadata,
}

pub type ImzMLReaderType = ImzMLReader;

impl ImzMLReader {
    /// Extract file-level metadata from imzML file
    fn extract_file_metadata<P: AsRef<Path>>(path: P) -> io::Result<ImzMLFileMetadata> {
        let file = fs::File::open(&path)?;
        let mut reader = Reader::from_reader(BufReader::new(file));
        let mut buf = Vec::new();
        let mut metadata = ImzMLFileMetadata::default();
        let mut in_file_content = false;
        
        loop {
            match reader.read_event_into(&mut buf) {
                Ok(Event::Start(ref e)) => {
                    if e.name().as_ref() == b"fileContent" {
                        in_file_content = true;
                    }
                }
                Ok(Event::Empty(ref e)) if in_file_content => {
                    if e.name().as_ref() == b"cvParam" {
                        let mut accession = None;
                        let mut value = None;
                        
                        // Parse attributes
                        for attr in e.attributes() {
                            if let Ok(attr) = attr {
                                match attr.key.as_ref() {
                                    b"accession" => {
                                        if let Ok(acc) = attr.unescape_value() {
                                            accession = Some(acc.to_string());
                                        }
                                    }
                                    b"value" => {
                                        if let Ok(val) = attr.unescape_value() {
                                            value = Some(val.to_string());
                                        }
                                    }
                                    _ => {}
                                }
                            }
                        }
                        
                        // Process IMS parameters
                        if let Some(acc) = accession {
                            match acc.as_str() {
                                "IMS:1000080" => {
                                    // UUID - strip braces and parse
                                    if let Some(uuid_str) = value {
                                        let cleaned_uuid = uuid_str.trim_matches(|c| c == '{' || c == '}');
                                        match Uuid::parse_str(cleaned_uuid) {
                                            Ok(uuid) => metadata.uuid = Some(uuid),
                                            Err(e) => warn!("Failed to parse UUID '{}': {}", cleaned_uuid, e),
                                        }
                                    }
                                }
                                "IMS:1000031" => {
                                    metadata.data_mode = Some(IbdDataMode::Processed);
                                }
                                "IMS:1000030" => {
                                    metadata.data_mode = Some(IbdDataMode::Continuous);
                                }
                                "IMS:1000090" => {
                                    metadata.ibd_checksum = value.clone();
                                    metadata.ibd_checksum_type = Some("MD5".to_string());
                                }
                                "IMS:1000091" => {
                                    metadata.ibd_checksum = value.clone();
                                    metadata.ibd_checksum_type = Some("SHA-1".to_string());
                                }
                                "IMS:1000092" => {
                                    metadata.ibd_checksum = value.clone();
                                    metadata.ibd_checksum_type = Some("SHA-256".to_string());
                                }
                                "IMS:1000070" => {
                                    metadata.ibd_file_path = value;
                                }
                                _ => {}
                            }
                        }
                    }
                }
                Ok(Event::End(ref e)) => {
                    if e.name().as_ref() == b"fileContent" {
                        in_file_content = false;
                        break; // We've processed all we need
                    }
                }
                Ok(Event::Eof) => break,
                Ok(_) => {}
                Err(_) => break,
            }
            buf.clear();
        }
        
        // If ibd_file_path is not set, derive it from imzML path
        if metadata.ibd_file_path.is_none() {
            let ibd_path = IbdFile::derive_ibd_path(&path);
            metadata.ibd_file_path = Some(ibd_path.to_string_lossy().to_string());
        }
        
        Ok(metadata)
    }

    /// Create a new ImzMLReader from a file path
    pub fn open_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let imzml_path = path.as_ref().to_path_buf();
        
        // First, extract file-level metadata
        let file_metadata = Self::extract_file_metadata(&imzml_path)?;
        
        // Create mzML reader
        let mzml_reader = MzMLReader::open_path(&imzml_path)?;
        
        // Determine data mode from metadata
        let data_mode = file_metadata.data_mode.unwrap_or(IbdDataMode::Unknown);
        
        // Validate required metadata
        if file_metadata.data_mode.is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Data mode (IMS:1000030 or IMS:1000031) not found in imzML file"
            ));
        }
        
        if file_metadata.uuid.is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "UUID (IMS:1000080) not found in imzML file"
            ));
        }
        
        if file_metadata.ibd_checksum.is_none() || file_metadata.ibd_checksum_type.is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "IBD checksum not found in imzML file"
            ));
        }
        
        // Determine IBD file path
        let ibd_path = if let Some(ref custom_path) = file_metadata.ibd_file_path {
            PathBuf::from(custom_path)
        } else {
            // Derive IBD path from imzML path
            IbdFile::derive_ibd_path(&imzml_path)
        };
        
        // Try to open the associated IBD file
        let ibd_file = if ibd_path.exists() {
            let ibd = IbdFile::open(&ibd_path, data_mode)?;
            
            // Validate UUID
            let expected_uuid = file_metadata.uuid.as_ref().unwrap();
            let ibd_uuid_bytes = ibd.uuid(); // This returns &[u8; 16]
            let ibd_uuid = Uuid::from_bytes(*ibd_uuid_bytes);
            if *expected_uuid != ibd_uuid {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("UUID mismatch between imzML ({}) and IBD ({})", expected_uuid, ibd_uuid)
                ));
            }
            
            // TODO: Validate checksum if needed
            
            Some(ibd)
        } else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("IBD file not found at expected location: {:?}", ibd_path)
            ));
        };

        Ok(Self {
            mzml_reader,
            ibd_file,
            imzml_path,
            data_mode,
            file_metadata,
        })
    }

    /// Get the next spectrum from the imzML file
    pub fn next_spectrum(&mut self) -> Option<MultiLayerSpectrum> {
        // TODO: Handle IBD data loading here
        // For now, just delegate to the mzML reader
        self.mzml_reader.next()
    }

    /// Get a spectrum by its ID
    pub fn get_spectrum_by_id(&mut self, id: &str) -> Option<MultiLayerSpectrum> {
        // TODO: Handle IBD data loading here
        self.mzml_reader.get_spectrum_by_id(id)
    }

    /// Get a spectrum by its index
    pub fn get_spectrum_by_index(&mut self, index: usize) -> Option<MultiLayerSpectrum> {
        // TODO: Handle IBD data loading here
        self.mzml_reader.get_spectrum_by_index(index)
    }

    /// Reset the reader to the beginning
    pub fn reset(&mut self) {
        self.mzml_reader.reset()
    }

    /// Set the IBD file for external binary data
    pub fn set_ibd_file(&mut self, ibd_file: IbdFile) {
        self.ibd_file = Some(ibd_file);
    }

    /// Get a reference to the IBD file
    pub fn ibd_file(&self) -> Option<&IbdFile> {
        self.ibd_file.as_ref()
    }

    /// Get a mutable reference to the IBD file
    pub fn ibd_file_mut(&mut self) -> Option<&mut IbdFile> {
        self.ibd_file.as_mut()
    }

    /// Get a reference to the underlying mzML reader
    pub fn mzml_reader(&self) -> &MzMLReader<fs::File> {
        &self.mzml_reader
    }

    /// Get a mutable reference to the underlying mzML reader
    pub fn mzml_reader_mut(&mut self) -> &mut MzMLReader<fs::File> {
        &mut self.mzml_reader
    }

    /// Get the file metadata extracted from the imzML
    pub fn file_metadata(&self) -> &ImzMLFileMetadata {
        &self.file_metadata
    }
}

impl Iterator for ImzMLReader {
    type Item = MultiLayerSpectrum;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_spectrum()
    }
}
