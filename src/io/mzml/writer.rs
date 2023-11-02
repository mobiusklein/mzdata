use std::collections::HashMap;
use std::fmt::{Debug, Display};
use std::io;
use std::io::{BufWriter, Seek, Write};
use std::marker::PhantomData;

use log::warn;

use mzpeaks::PeakCollection;
use quick_xml::events::BytesDecl;
use quick_xml::events::{BytesEnd, BytesStart, BytesText, Event};
use quick_xml::Error as XMLError;
use quick_xml::Writer;

use super::super::offset_index::OffsetIndex;
use super::super::traits::ScanWriter;
use super::super::utils::MD5HashingStream;

use mzpeaks::{peak_set::PeakSetVec, CentroidPeak, DeconvolutedPeak, Mass, MZ};

use crate::meta::file_description::FileDescription;
use crate::meta::instrument::{ComponentType, InstrumentConfiguration};
use crate::meta::{DataProcessing, MSDataFileMetadata, Software};
use crate::params::{ControlledVocabulary, Param, ParamLike, Unit, ParamCow};
use crate::spectrum::scan_properties::*;
use crate::spectrum::signal::{
    ArrayType, BinaryArrayMap, BinaryCompressionType, BinaryDataArrayType, DataArray,
};
use crate::spectrum::spectrum::{
    CentroidPeakAdapting, DeconvolutedPeakAdapting, MultiLayerSpectrum,
};
use crate::ParamDescribed;
use crate::SpectrumBehavior;

const BUFFER_SIZE: usize = 10000;

macro_rules! bstart {
    ($e:tt) => {
        BytesStart::from_content($e, $e.len())
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
        $writer.handle.write_event(Event::Start($target.borrow()))?;
    };
}

macro_rules! end_event {
    ($writer:ident, $target:ident) => {
        $writer.handle.write_event(Event::End($target.to_end()))?;
    };
}


fn instrument_id(id: &u32) -> String {
    format!("IC{}", *id + 1)
}


const MS1_SPECTRUM: ParamCow = ControlledVocabulary::MS.const_param_ident("MS1 spectrum", 1000579);
const MSN_SPECTRUM: ParamCow = ControlledVocabulary::MS.const_param_ident("MSn spectrum", 1000580);
const NEGATIVE_SCAN: ParamCow = ControlledVocabulary::MS.const_param_ident("negative scan", 1000129);
const POSITIVE_SCAN: ParamCow = ControlledVocabulary::MS.const_param_ident("positive scan", 1000130);
const PROFILE_SPECTRUM: ParamCow = ControlledVocabulary::MS.const_param_ident("profile spectrum", 1000128);
const CENTROID_SPECTRUM: ParamCow = ControlledVocabulary::MS.const_param_ident("centroid spectrum", 1000127);


#[derive(Debug)]
pub enum MzMLWriterError {
    XMLError(XMLError),
    StateTransitionError {from_state: MzMLWriterState, to_state: MzMLWriterState},
    IOError(io::Error),
    InvalidActionError(MzMLWriterState),
}

impl Display for MzMLWriterError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::error::Error for MzMLWriterError {}

impl From<XMLError> for MzMLWriterError {
    fn from(value: XMLError) -> Self {
        Self::XMLError(value)
    }
}

impl From<io::Error> for MzMLWriterError {
    fn from(value: io::Error) -> Self {
        Self::IOError(value)
    }
}


pub type WriterResult = Result<(), MzMLWriterError>;

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
        let digest = self.handle.get_mut().get_ref().compute();
        format!("{:x}", digest)
    }

    pub fn flush(&mut self) -> io::Result<()> {
        self.handle.get_mut().flush()
    }

    pub fn write_param<P: ParamLike + Debug>(&mut self, param: &P) -> WriterResult {
        let mut elt = if !param.is_controlled() {
            bstart!("userParam")
        } else {
            let mut elt = bstart!("cvParam");
            let accession_str = param.curie().unwrap();
            attrib!("accession", accession_str, elt);
            if let Some(cv_ref) = &param.controlled_vocabulary() {
                attrib!("cvRef", cv_ref, elt);
            }
            elt
        };

        attrib!("name", param.name(), elt);
        if !param.value().is_empty() {
            attrib!("value", param.value(), elt);
        }
        match param.unit() {
            Unit::Unknown => {}
            unit => {
                let (unit_acc, unit_name) = unit.for_param();
                let mut split = unit_acc.split(':');
                if let Some(prefix) = split.next() {
                    attrib!("unitCvRef", prefix, elt);
                } else {
                    attrib!("unitCvRef", "UO", elt);
                }
                attrib!("unitAccession", unit_acc, elt);
                attrib!("unitName", unit_name, elt);
            }
        }
        self.handle.write_event(Event::Empty(elt))?;
        Ok(())
    }

    pub fn write_param_list<'a, T: Iterator<Item = &'a Param>>(
        &mut self,
        params: T,
    ) -> WriterResult {
        for param in params {
            self.write_param(param)?
        }
        Ok(())
    }

    pub fn write_event(&mut self, event: Event) -> WriterResult {
        self.handle.write_event(event)?;
        Ok(())
    }
}

/**
The different states that [`MzMLWriterType`] can enter while
writing an mzML document. This is only necessary for the module
consumer when determining where something may have gone wrong.
*/
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

pub trait MzMLSpectrumWriter<'a,
    C: CentroidPeakAdapting + 'static = CentroidPeak,
    D: DeconvolutedPeakAdapting + 'static = DeconvolutedPeak,
>
{
    fn write_scan_list(&mut self, acq: &Acquisition) -> WriterResult;
    fn write_isolation_window(&mut self, iw: &IsolationWindow) -> WriterResult;
    fn write_activation(&mut self, precursor: &Precursor) -> WriterResult;
    fn write_precursor(&mut self, precursor: &Precursor) -> WriterResult;

    fn write_binary_data_array(&mut self, array: &DataArray) -> WriterResult;
    fn write_binary_data_arrays(&mut self, arrays: &BinaryArrayMap) -> WriterResult;

    fn write_spectrum(&'a mut self, spectrum: &'a MultiLayerSpectrum<C, D>) -> WriterResult;
}

/**
An indexed mzML writer that handles [`Spectrum`](crate::spectrum::Spectrum)
by reference.

Does not buffer spectra in-memory, writing them out immediately.

Currently, no chromatograms are written (including no TIC or base peak chromatograms).
*/
#[derive(Debug)]
pub struct MzMLWriterType<
    W: Write + Seek,
    C: CentroidPeakAdapting + 'static = CentroidPeak,
    D: DeconvolutedPeakAdapting + 'static = DeconvolutedPeak,
> {
    /// The current offset from the stream start
    pub offset: usize,

    /// The total number of spectra this mzML document will contain.
    /// This value will appear in the `spectrumList` element's count attribute
    pub spectrum_count: u64,
    /// The number of `spectrum` elements written so far.
    pub spectrum_counter: u64,

    /// The total number of chromatograms this mzML document will contain.
    /// This value will appear in the `chromatogramList` element's count attribute
    pub chromatogram_count: u64,
    /// The number of chromatograms written so far
    pub chromatogram_counter: u64,

    /// The compression type to use when generating binary data arrays.
    pub data_array_compression: BinaryCompressionType,

    /// The file-level metadata describing the provenance of the original data
    pub file_description: FileDescription,
    /// The list of software components that were used to process the data into
    /// its current state
    pub softwares: Vec<Software>,
    /// The types of data transformations applied to (parts of) the data
    pub data_processings: Vec<DataProcessing>,
    /// The different instrument configurations that were in use during the
    /// data acquisition.
    pub instrument_configurations: HashMap<u32, InstrumentConfiguration>,

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

    /// Wrap a new [`std::io::Write`]-able type, constructing a new [`MzMLWriterType`]
    pub fn new(file: W) -> MzMLWriterType<W, C, D> {
        let handle = InnerXMLWriter::new(file);
        MzMLWriterType {
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
            ms_cv: ControlledVocabulary::MS,
            data_array_compression: BinaryCompressionType::Zlib,
        }
    }

    fn transition_err(&self, to_state: MzMLWriterState) -> WriterResult {
        Err(MzMLWriterError::StateTransitionError { from_state: self.state, to_state })
    }

    fn stream_position(&mut self) -> io::Result<u64> {
        self.handle.handle.get_mut().stream_position()
    }

    fn make_psi_ms_cv(&self) -> BytesStart<'static> {
        let mut cv = BytesStart::from_content("cv", 2);
        cv.push_attribute(("id", "MS"));
        cv.push_attribute(("fullName", "PSI-MS"));
        cv.push_attribute((
            "URI",
            "http://purl.obolibrary.org/obo/ms.obo",
        ));
        cv.push_attribute(("version", Self::PSIMS_VERSION));
        cv
    }

    fn make_unit_cv(&self) -> BytesStart<'static> {
        let mut cv = BytesStart::from_content("cv", 2);
        cv.push_attribute(("id", "UO"));
        cv.push_attribute(("fullName", "UNIT-ONTOLOGY"));
        cv.push_attribute(("URI", "http://ontologies.berkeleybop.org/uo.obo"));
        cv.push_attribute(("version", Self::UNIT_VERSION));
        cv
    }

    fn write_cv_list(&mut self) -> WriterResult {
        let mut cv_list = BytesStart::from_content("cvList", 6);
        cv_list.push_attribute(("count", "2"));
        self.handle.write_event(Event::Start(cv_list))?;

        let cv = self.make_psi_ms_cv();
        self.handle.write_event(Event::Empty(cv))?;

        let cv = self.make_unit_cv();
        self.handle.write_event(Event::Empty(cv))?;

        self.handle
            .write_event(Event::End(BytesEnd::new("cvList")))?;
        Ok(())
    }

    fn start_document(&mut self) -> WriterResult {
        self.handle
            .write_event(Event::Decl(BytesDecl::new("1.0", Some("utf-8"), None)))?;
        let mut indexed = BytesStart::from_content("indexedmzML", 11);
        indexed.push_attribute(("xmlns", "http://psi.hupo.org/ms/mzml"));
        indexed.push_attribute(("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance"));
        indexed.push_attribute((
            "xsi:schemaLocation",
            "http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.3_idx.xsd",
        ));
        self.handle.write_event(Event::Start(indexed))?;

        let mut mzml = BytesStart::from_content("mzML", 4);
        mzml.push_attribute(("xmlns", "http://psi.hupo.org/ms/mzml"));
        mzml.push_attribute(("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance"));
        mzml.push_attribute((
            "xsi:schemaLocation",
            "http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.1.xsd",
        ));
        mzml.push_attribute(("version", "1.1.1"));
        self.handle.write_event(Event::Start(mzml))?;

        self.state = MzMLWriterState::DocumentOpen;
        Ok(())
    }

    fn write_header(&mut self) -> WriterResult {
        if self.state < MzMLWriterState::DocumentOpen {
            self.start_document()?;
        } else {
            return self.transition_err(MzMLWriterState::Header)
        }
        self.write_cv_list()?;
        self.write_file_description()?;
        self.write_software_list()?;
        self.write_instrument_configuration()?;
        self.write_data_processing()?;

        self.state = MzMLWriterState::Header;
        Ok(())
    }

    fn write_file_description(&mut self) -> WriterResult {
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
        self.handle.write_event(Event::Start(outer.borrow()))?;
        for sf in self.file_description.source_files.iter() {
            let mut tag = bstart!("sourceFile");
            attrib!("id", sf.id, tag);
            attrib!("name", sf.name, tag);
            attrib!("location", sf.location, tag);
            self.handle.write_event(Event::Start(tag.borrow()))?;
            for param in sf.params() {
                self.handle.write_param(param)?
            }
            self.handle.write_event(Event::End(tag.to_end()))?;
        }
        self.handle.write_event(Event::End(outer.to_end()))?;

        end_event!(self, fd);
        Ok(())
    }

    fn write_software_list(&mut self) -> WriterResult {
        let mut outer = bstart!("softwareList");
        let count = self.softwares.len().to_string();
        attrib!("count", count, outer);
        self.handle.write_event(Event::Start(outer.borrow()))?;
        for soft in self.softwares.iter() {
            let mut tag = bstart!("software");
            attrib!("id", soft.id, tag);
            attrib!("version", soft.version, tag);
            self.handle.write_event(Event::Start(tag.borrow()))?;
            for param in soft.params() {
                self.handle.write_param(param)?
            }
            self.handle.write_event(Event::End(tag.to_end()))?;
        }
        self.handle.write_event(Event::End(outer.to_end()))?;
        Ok(())
    }

    fn write_instrument_configuration(&mut self) -> WriterResult {
        let mut outer = bstart!("instrumentConfigurationList");
        let count = self.instrument_configurations.len().to_string();
        attrib!("count", count, outer);
        self.handle.write_event(Event::Start(outer.borrow()))?;

        // Sort the keys so the ordering is consistent on every run
        let mut configs: Vec<_> = self.instrument_configurations.keys().collect();
        configs.sort();
        for key in configs {
            let ic = self.instrument_configurations.get(key).unwrap();
            let mut tag = bstart!("instrumentConfiguration");
            let inst_id = instrument_id(&ic.id);
            attrib!("id", inst_id, tag);
            self.handle.write_event(Event::Start(tag.borrow()))?;
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
                self.handle.write_event(Event::Start(cmp_tag.borrow()))?;
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

    fn write_data_processing(&mut self) -> WriterResult {
        let mut outer = bstart!("dataProcessingList");
        let count = self.data_processings.len().to_string();
        attrib!("count", count, outer);
        self.handle.write_event(Event::Start(outer.borrow()))?;
        for dp in self.data_processings.iter() {
            let mut tag = bstart!("dataProcessing");
            attrib!("id", dp.id, tag);
            self.handle.write_event(Event::Start(tag.borrow()))?;
            for proc in dp.methods.iter() {
                let mut mtag = bstart!("processingMethod");
                let order = proc.order.to_string();
                attrib!("order", order, mtag);
                attrib!("softwareRef", proc.software_reference, mtag);
                self.handle.write_event(Event::Start(mtag.borrow()))?;
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

    fn start_run(&mut self) -> WriterResult {
        if self.state < MzMLWriterState::Run {
            self.write_header()?;
        } else {
            return self.transition_err(MzMLWriterState::Run);
        }
        let mut run = bstart!("run");
        attrib!("id", "1", run);
        let mut keys: Vec<_> = self.instrument_configurations.keys().collect();
        keys.sort();
        if let Some(ic_ref) = keys.first() {
            let inst_id = instrument_id(&ic_ref);
            attrib!("defaultInstrumentConfigurationRef", inst_id, run);
        }
        if let Some(sf_ref) = self.file_description.source_files.first() {
            attrib!("defaultSourceFileRef", sf_ref.id, run);
        };
        self.handle.write_event(Event::Start(run))?;
        self.state = MzMLWriterState::Run;
        Ok(())
    }

    fn start_spectrum_list(&mut self) -> WriterResult {
        match self.state {
            MzMLWriterState::SpectrumList => {}
            state if state < MzMLWriterState::SpectrumList => {
                self.start_run()?;
            }
            state if state > MzMLWriterState::SpectrumList => {
                return self.transition_err(MzMLWriterState::SpectrumList);
            }
            _ => {}
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

    fn close_spectrum_list(&mut self) -> WriterResult {
        let tag = bstart!("spectrumList");
        end_event!(self, tag);
        self.state = MzMLWriterState::SpectrumListClosed;
        Ok(())
    }

    fn close_run(&mut self) -> WriterResult {
        if self.state < MzMLWriterState::Run {
            self.start_run()?;
        } else if self.state == MzMLWriterState::SpectrumList {
            self.close_spectrum_list()?;
        } else if self.state == MzMLWriterState::ChromatogramList {
            // TODO
        } else if self.state > MzMLWriterState::RunClosed {
            // Cannot close the run of mzML, currently in state which happens after the run has already ended
            return self.transition_err(MzMLWriterState::RunClosed);
        }
        let tag = bstart!("run");
        end_event!(self, tag);
        self.state = MzMLWriterState::RunClosed;
        Ok(())
    }

    fn close_mzml(&mut self) -> WriterResult {
        if self.state < MzMLWriterState::RunClosed {
            self.close_run()?;
        }
        let tag = bstart!("mzML");
        self.state = MzMLWriterState::MzMLClosed;
        end_event!(self, tag);
        Ok(())
    }

    fn close_indexed_mzml(&mut self) -> WriterResult {
        if self.state < MzMLWriterState::MzMLClosed {
            self.close_mzml()?;
        }
        self.write_index_list()?;
        let tag = bstart!("indexedmzML");
        end_event!(self, tag);
        self.state = MzMLWriterState::End;
        Ok(())
    }

    /**
    Close the wrapping `<indexedmzML>` document, which will trigger writing
    out the offset indices and file checksum at the tail of the document.
    */
    pub fn close(&mut self) -> WriterResult {
        if self.state < MzMLWriterState::End {
            self.close_indexed_mzml()
        } else {
            Ok(())
        }
    }

    fn write_scan_list(&mut self, acq: &Acquisition) -> WriterResult {
        let mut scan_list_tag = bstart!("scanList");
        let count = acq.scans.len().to_string();
        attrib!("count", count, scan_list_tag);
        start_event!(self, scan_list_tag);
        self.handle.write_param(&acq.combination.to_param())?;

        for scan in acq.scans.iter() {
            let mut scan_tag = bstart!("scan");
            let id = instrument_id(&scan.instrument_configuration_id);
            attrib!("instrumentConfigurationRef", id, scan_tag);
            self.handle.write_event(Event::Start(scan_tag.borrow()))?;

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
                .write_event(Event::Start(scan_window_list_tag.borrow()))?;
            for window in scan.scan_windows.iter() {
                let window_tag = bstart!("scanWindow");
                self.handle.write_event(Event::Start(window_tag.borrow()))?;
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

    fn write_isolation_window(&mut self, iw: &IsolationWindow) -> WriterResult {
        let iw_tag = bstart!("isolationWindow");
        self.handle.write_event(Event::Start(iw_tag.borrow()))?;
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

    fn write_selected_ions(&mut self, precursor: &Precursor) -> WriterResult {
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

    fn write_activation(&mut self, precursor: &Precursor) -> WriterResult {
        let act = precursor.activation();
        let tag = bstart!("activation");
        start_event!(self, tag);
        match act.method() {
            Some(meth) => {
                self.handle.write_param(meth)?;
            }
            None => {}
        }
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

    fn write_precursor(&mut self, precursor: &Precursor) -> WriterResult {
        let mut precursor_list_tag = bstart!("precursorList");
        attrib!("count", "1", precursor_list_tag);
        start_event!(self, precursor_list_tag);

        let mut precursor_tag = bstart!("precursor");
        if let Some(prec_id) = precursor.precursor_id() {
            attrib!("spectrumRef", prec_id, precursor_tag);
        }
        self.handle
            .write_event(Event::Start(precursor_tag.borrow()))?;

        let iw = precursor.isolation_window();
        self.write_isolation_window(iw)?;
        self.write_selected_ions(precursor)?;
        self.write_activation(precursor)?;
        end_event!(self, precursor_tag);
        end_event!(self, precursor_list_tag);
        Ok(())
    }

    fn write_binary_data_array(&mut self, array: &DataArray) -> WriterResult {
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
            ArrayType::MZArray => self.handle.write_param(&array.name.as_param_const())?,
            ArrayType::IntensityArray => self.handle.write_param(&array.name.as_param_const())?,
            ArrayType::ChargeArray => self.handle.write_param(&array.name.as_param_const())?,
            ArrayType::TimeArray => self.handle.write_param(&array.name.as_param_const())?,
            ArrayType::RawIonMobilityArray => self.handle.write_param(&array.name.as_param_with_unit_const(array.unit))?,
            ArrayType::MeanIonMobilityArray => self.handle.write_param(&array.name.as_param_with_unit_const(array.unit))?,
            ArrayType::DeconvolutedIonMobilityArray => self.handle.write_param(&array.name.as_param_with_unit_const(array.unit))?,
            ArrayType::NonStandardDataArray { name } => {
                // self.handle.write_param()?;
                let mut p = self
                    .ms_cv
                    .param_val("MS:1000786", "non-standard data array", name);
                p = p.with_unit_t(&array.unit);
                self.handle.write_param(&p)?;
            }
            _ => {
                panic!("Could not determine how to name for {:?}", array.name);
            }
        }

        let bin = bstart!("binary");
        start_event!(self, bin);
        self.handle.write_event(Event::Text(BytesText::new(
            String::from_utf8_lossy(&encoded).as_ref(),
        )))?;
        end_event!(self, bin);
        end_event!(self, outer);
        Ok(())
    }

    fn write_binary_data_arrays(&mut self, arrays: &BinaryArrayMap) -> WriterResult {
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

    /**
    Write a [`MultiLayerSpectrum`] out to the mzML file, encoding the highest procressing
    degree peak data present.

    ## Side-Effects
    If the writer has not already started writing the spectra, this will cause all the metadata
    to be written out and the `<spectrumList>` element will be opened, preventing no new metadata
    from being written to this stream. Furthermore, this writes the spectrum count out, so the value
    may no longer be changed.
    */
    pub fn write_spectrum(&mut self, spectrum: &'a MultiLayerSpectrum<C, D>) -> WriterResult {
        match self.state {
            MzMLWriterState::SpectrumList => {}
            state if state < MzMLWriterState::SpectrumList => {
                self.start_spectrum_list()?;
            }
            state if state > MzMLWriterState::SpectrumList => {
                // Cannot write spectrum, currently in state which happens
                // after spectra may be written
                return Err(MzMLWriterError::InvalidActionError(self.state))
            }
            _ => {}
        }
        let pos = self.stream_position()? - (1 + (4 * InnerXMLWriter::<W>::INDENT_SIZE));
        self.offset_index.insert(spectrum.id().to_string(), pos);
        let mut outer = bstart!("spectrum");
        attrib!("id", spectrum.id(), outer);
        let count = self.spectrum_counter.to_string();
        attrib!("index", count, outer);
        let default_array_len = if let Some(mass_peaks) = &spectrum.deconvoluted_peaks {
            mass_peaks.len()
        } else if let Some(mz_peaks) = &spectrum.peaks {
            mz_peaks.len()
        } else if let Some(arrays) = &spectrum.arrays {
            arrays.mzs().len()
        } else {
            0
        }.to_string();

        attrib!("defaultArrayLength", default_array_len, outer);

        self.handle.write_event(Event::Start(outer.borrow()))?;
        self.spectrum_counter += 1;

        let ms_level = spectrum.ms_level();
        if ms_level == 1 {
            self.handle
                .write_param(&MS1_SPECTRUM)?;
        } else {
            self.handle
                .write_param(&MSN_SPECTRUM)?;
        }
        self.handle.write_param(&self.ms_cv.param_val(
            "MS:1000511",
            "ms level",
            ms_level.to_string(),
        ))?;

        match spectrum.polarity() {
            ScanPolarity::Negative => self
                .handle
                .write_param(&NEGATIVE_SCAN)?,
            ScanPolarity::Positive => self
                .handle
                .write_param(&POSITIVE_SCAN)?,
            ScanPolarity::Unknown => {
                warn!(
                    "Could not determine scan polarity for {}, assuming positive",
                    spectrum.id()
                );
                self.handle
                    .write_param(&POSITIVE_SCAN)?
            }
        }

        match spectrum.signal_continuity() {
            SignalContinuity::Profile => self
                .handle
                .write_param(&PROFILE_SPECTRUM)?,
            SignalContinuity::Unknown => {
                warn!(
                    "Could not determine scan polarity for {}, assuming centroid",
                    spectrum.id()
                );
                self.handle
                    .write_param(&CENTROID_SPECTRUM)?;
            }
            _ => {
                self.handle
                    .write_param(&CENTROID_SPECTRUM)?;
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

    fn write_index(&mut self, index: &OffsetIndex) -> WriterResult {
        let mut outer = bstart!("index");
        attrib!("name", index.name, outer);
        start_event!(self, outer);
        for (id, offset) in index.iter() {
            let mut tag = bstart!("offset");
            attrib!("idRef", id, tag);
            start_event!(self, tag);
            let content = offset.to_string();
            let text = BytesText::new(&content);
            self.handle.write_event(Event::Text(text))?;
            end_event!(self, tag);
        }
        end_event!(self, outer);
        Ok(())
    }

    fn write_index_list(&mut self) -> WriterResult {
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
        let text = BytesText::new(&content);
        self.handle.write_event(Event::Text(text))?;
        end_event!(self, tag);

        let tag = bstart!("fileChecksum");
        start_event!(self, tag);
        let content = self.handle.digest();
        let text = BytesText::new(&content);
        self.handle.write_event(Event::Text(text))?;
        end_event!(self, tag);
        Ok(())
    }

    /// Get a reference to the mzML writer's spectrum count.
    pub fn spectrum_count(&self) -> &u64 {
        &self.spectrum_count
    }

    /// Set the mzML writer's spectrum count.
    pub fn set_spectrum_count(&mut self, spectrum_count: u64) {
        self.spectrum_count = spectrum_count;
    }

    /// Get a mutable reference to the mzML writer's spectrum count to modify in-place.
    pub fn spectrum_count_mut(&mut self) -> &mut u64 {
        &mut self.spectrum_count
    }
}


impl<'a, C: CentroidPeakAdapting + 'static, D: DeconvolutedPeakAdapting + 'static, W: Write + Seek> MzMLSpectrumWriter<'a, C, D> for MzMLWriterType<W, C, D> where
    &'a PeakSetVec<C, MZ>: Into<BinaryArrayMap>,
    &'a PeakSetVec<D, Mass>: Into<BinaryArrayMap> {

    fn write_scan_list(&mut self, acq: &Acquisition) -> WriterResult {
        self.write_scan_list(acq)
    }

    fn write_isolation_window(&mut self, iw: &IsolationWindow) -> WriterResult {
        self.write_isolation_window(iw)
    }

    fn write_activation(&mut self, precursor: &Precursor) -> WriterResult {
        self.write_activation(precursor)
    }

    fn write_precursor(&mut self, precursor: &Precursor) -> WriterResult {
        self.write_precursor(precursor)
    }

    fn write_binary_data_array(&mut self, array: &DataArray) -> WriterResult {
        self.write_binary_data_array(array)
    }

    fn write_binary_data_arrays(&mut self, arrays: &BinaryArrayMap) -> WriterResult {
        self.write_binary_data_arrays(arrays)
    }

    fn write_spectrum(&'a mut self, spectrum: &'a MultiLayerSpectrum<C, D>) -> WriterResult {
        self.write_spectrum(spectrum)
    }
}

/// A specialization of [`MzMLWriterType`] for the default peak types, for common use.
pub type MzMLWriter<W> = MzMLWriterType<CentroidPeak, DeconvolutedPeak, W>;

#[cfg(test)]
mod test {
    use super::super::reader::{MzMLReader, MzMLReaderType};
    use super::*;
    use crate::io::prelude::*;
    use std::fs;
    use std::path;

    #[test]
    fn write_test() -> WriterResult {
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
