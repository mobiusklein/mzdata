use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::io::{BufWriter, Write};
use std::marker::PhantomData;
use std::{borrow::Cow, io, mem};

use log::warn;
use mzpeaks::feature::FeatureLike;
#[cfg(feature = "parallelism")]
use rayon::prelude::*;
use thiserror::Error;

use mzpeaks::{CentroidLike, DeconvolutedCentroidLike, IonMobility, KnownCharge, Mass, MZ};
use quick_xml::escape;
use quick_xml::events::{BytesDecl, BytesEnd, BytesStart, BytesText, Event};
use quick_xml::{Error as XMLError, Writer};

use super::super::offset_index::OffsetIndex;
use super::super::traits::SpectrumWriter;
use super::super::utils::MD5HashingStream;

use mzpeaks::{CentroidPeak, DeconvolutedPeak};

use crate::io::traits::IonMobilityFrameWriter;
use crate::meta::{
    ComponentType, DataProcessing, FileDescription, InstrumentConfiguration, MSDataFileMetadata, MassSpectrometryRun, Sample, Software
};
use crate::params::{
    ControlledVocabulary, Param, ParamCow, ParamDescribed, ParamLike, ParamValue, Unit, ValueRef,
};
use crate::spectrum::bindata::{
    to_bytes, ArrayRetrievalError, ArrayType, BinaryArrayMap, BinaryCompressionType,
    BinaryDataArrayType, BuildArrayMap3DFrom, BuildArrayMapFrom, ByteArrayView, DataArray,
};
use crate::spectrum::spectrum_types::SpectrumLike;
use crate::spectrum::{scan_properties::*, Chromatogram, ChromatogramLike, RefPeakDataLevel};
use crate::{curie, impl_param_described, RawSpectrum};

const BUFFER_SIZE: usize = 10000;

#[cfg(feature = "parallelism")]
const PARALLEL_COMPRESSION_FAN: usize = 3;

macro_rules! bstart {
    ($e:tt) => {
        BytesStart::from_content($e, $e.len())
    };
}

macro_rules! attrib {
    ($name:expr, $value:expr, $elt:ident) => {
        let key = $name.as_bytes();
        let value = $value.as_bytes();
        // Because quick_xml::escape does not escape newlines
        let decoded = unsafe { std::str::from_utf8_unchecked(&value) };
        let mut escaped_value = escape::escape(&decoded);
        if escaped_value.contains(['\n', '\r']) {
            escaped_value = Cow::Owned(escaped_value.replace('\n', "&#10;"));
            if escaped_value.contains('\r') {
                escaped_value = Cow::Owned(escaped_value.replace('\r', "&#13;"));
            }
        }
        $elt.push_attribute((key, escaped_value.as_bytes()));
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

fn param_val_xsd<V: ParamValue>(param: &V) -> Option<&[u8]> {
    if param.is_f64() {
        Some(b"xsd:double")
    } else if param.is_i64() {
        Some(b"xsd:integer")
    } else if param.is_str() {
        Some(b"xsd:string")
    } else if param.is_boolean() {
        Some(b"xsd:boolean")
    } else {
        None
    }
}

const MS1_SPECTRUM: ParamCow = ControlledVocabulary::MS.const_param_ident("MS1 spectrum", 1000579);
const MSN_SPECTRUM: ParamCow = ControlledVocabulary::MS.const_param_ident("MSn spectrum", 1000580);
const NEGATIVE_SCAN: ParamCow =
    ControlledVocabulary::MS.const_param_ident("negative scan", 1000129);
const POSITIVE_SCAN: ParamCow =
    ControlledVocabulary::MS.const_param_ident("positive scan", 1000130);
const PROFILE_SPECTRUM: ParamCow =
    ControlledVocabulary::MS.const_param_ident("profile spectrum", 1000128);
const CENTROID_SPECTRUM: ParamCow =
    ControlledVocabulary::MS.const_param_ident("centroid spectrum", 1000127);

#[allow(unused)]
struct ParamPack {
    name: &'static str,
    accession: u32,
    unit: Unit,
}

#[allow(unused)]
impl ParamPack {
    const fn new(name: &'static str, accession: u32, unit: Unit) -> Self {
        Self {
            name,
            accession,
            unit,
        }
    }

    fn pack<V: Into<ValueRef<'static>>>(&self, value: V) -> ParamCow<'static> {
        let value: ValueRef<'static> = value.into();
        ControlledVocabulary::MS.const_param(self.name, value, self.accession, self.unit)
    }
}

const TIC_TERM: ParamPack = ParamPack::new("total ion current", 1000285, Unit::DetectorCounts);
#[allow(unused)]
const BPI_TERM: ParamPack = ParamPack::new("base peak intensity", 1000505, Unit::DetectorCounts);
const BPMZ_TERM: ParamPack = ParamPack::new("base peak m/z", 1000504, Unit::MZ);

/// All the ways that mzML writing can go wrong
#[derive(Debug, Error)]
pub enum MzMLWriterError {
    #[error("An XML error occurred: {0}")]
    XMLError(
        #[from]
        #[source]
        XMLError,
    ),

    #[error("Attempted to transition from {from_state:?} to {to_state:?}")]
    StateTransitionError {
        from_state: MzMLWriterState,
        to_state: MzMLWriterState,
    },
    #[error("An IO error occurred: {0}")]
    IOError(
        #[from]
        #[source]
        io::Error,
    ),

    #[error("An error occurred while retrieving array: {0}")]
    ArrayRetrievalError(
        #[from]
        #[source]
        ArrayRetrievalError,
    ),

    #[error("Attempted to perform an invalid action {0:?}")]
    InvalidActionError(MzMLWriterState),
}

impl From<MzMLWriterError> for io::Error {
    fn from(value: MzMLWriterError) -> Self {
        match value {
            MzMLWriterError::XMLError(e) => match e {
                XMLError::Io(o) => io::Error::new(o.kind(), o),
                _ => io::Error::new(io::ErrorKind::InvalidData, e),
            },
            MzMLWriterError::StateTransitionError {
                from_state: _,
                to_state: _,
            } => io::Error::new(io::ErrorKind::InvalidData, value),
            MzMLWriterError::IOError(o) => o,
            MzMLWriterError::ArrayRetrievalError(e) => {
                io::Error::new(io::ErrorKind::InvalidData, e)
            }
            MzMLWriterError::InvalidActionError(_) => {
                io::Error::new(io::ErrorKind::InvalidData, value)
            }
        }
    }
}

pub type WriterResult = Result<(), MzMLWriterError>;

struct ByteCountingStream<W: io::Write> {
    stream: BufWriter<MD5HashingStream<W>>,
    bytes_written: u64,
}

impl<W: io::Write> ByteCountingStream<W> {
    pub fn new(stream: BufWriter<MD5HashingStream<W>>) -> Self {
        Self {
            stream,
            bytes_written: 0,
        }
    }

    pub fn bytes_written(&self) -> u64 {
        self.bytes_written
    }

    pub fn checksum(&self) -> md5::Digest {
        self.stream.get_ref().compute()
    }

    pub fn get_mut(&mut self) -> &mut W {
        self.stream.get_mut().get_mut()
    }
}

impl<W: Write> Write for ByteCountingStream<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let wrote = self.stream.write(buf)?;
        self.bytes_written += wrote as u64;
        Ok(wrote)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.stream.flush()
    }
}

struct InnerXMLWriter<W: io::Write> {
    pub handle: Writer<ByteCountingStream<W>>,
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
        let handle = ByteCountingStream::new(BufWriter::with_capacity(
            BUFFER_SIZE,
            MD5HashingStream::new(file),
        ));
        Self {
            handle: Writer::new_with_indent(handle, b' ', 2),
        }
    }

    pub fn digest(&mut self) -> String {
        let digest = self.handle.get_ref().checksum();
        format!("{:x}", digest)
    }

    pub fn flush(&mut self) -> io::Result<()> {
        self.handle.get_mut().flush()
    }

    pub fn write_param<P: ParamLike>(&mut self, param: &P) -> WriterResult {
        let mut elt = if !param.is_controlled() {
            let mut elt = bstart!("userParam");
            if let Some(tp) = param_val_xsd(&param.value()) {
                elt.push_attribute(("type".as_bytes(), tp));
            }
            elt
        } else {
            let mut elt = bstart!("cvParam");
            let accession_str = param.curie().unwrap().to_string();
            attrib!("accession", accession_str, elt);
            if let Some(cv_ref) = &param.controlled_vocabulary() {
                attrib!("cvRef", cv_ref, elt);
            }
            elt
        };

        attrib!("name", param.name(), elt);
        let param_val = param.value();
        if !param_val.is_empty() {
            attrib!("value", param_val, elt);
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

    pub fn write_param_list<'a, P: ParamLike + 'a, T: Iterator<Item = &'a P>>(
        &mut self,
        params: T,
    ) -> WriterResult {
        for param in params {
            self.write_param(param)?
        }
        Ok(())
    }

    pub fn write_reference_group(&mut self, param_group: &ParamGroup) -> WriterResult {
        let mut tag_ref = bstart!("referenceableParamGroupRef");
        attrib!("ref", param_group.id, tag_ref);
        self.handle.write_event(Event::Empty(tag_ref))?;
        Ok(())
    }

    pub fn write_param_list_with_refs<'a, T: Iterator<Item = &'a Param>>(
        &mut self,
        params: T,
        param_groups: &[ParamGroup],
    ) -> WriterResult {
        let mut params: Vec<_> = params.collect();
        let mut ref_tags = Vec::new();
        for param_group in param_groups.iter() {
            if param_group.covers(params.as_ref()) {
                param_group.filter(&mut params);
                ref_tags.push(param_group);
            }
        }

        for p in params {
            self.write_param(p)?;
        }

        for pg in ref_tags {
            self.write_reference_group(pg)?;
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

#[derive(Debug, Clone, Default)]
pub struct ChromatogramCollector {
    name: ChromatogramType,
    time: Vec<u8>,
    intensity: Vec<u8>,
}

impl ChromatogramCollector {
    pub fn of(chromatogram_type: ChromatogramType) -> Self {
        Self {
            name: chromatogram_type,
            ..Default::default()
        }
    }

    pub fn add(&mut self, time: f64, intensity: f32) {
        self.time.extend(time.to_le_bytes());
        self.intensity.extend(intensity.to_le_bytes());
    }

    pub fn to_chromatogram(&self) -> Chromatogram {
        let mut descr = ChromatogramDescription::default();
        let mut arrays = BinaryArrayMap::default();
        match self.name {
            ChromatogramType::TotalIonCurrentChromatogram => {
                descr.id = "TIC".to_string();
                descr.add_param(
                    ControlledVocabulary::MS
                        .const_param_ident("total ion current chromatogram", 1000235)
                        .into(),
                );
                let mut time_array = DataArray::wrap(
                    &ArrayType::TimeArray,
                    BinaryDataArrayType::Float64,
                    to_bytes(&self.time),
                );
                time_array.unit = Unit::Minute;

                let intensity_array = DataArray::wrap(
                    &ArrayType::IntensityArray,
                    BinaryDataArrayType::Float32,
                    to_bytes(&self.intensity),
                );
                arrays.add(time_array);
                arrays.add(intensity_array);
            }
            ChromatogramType::BasePeakChromatogram => {
                descr.id = "BIC".to_string();
                descr.add_param(
                    ControlledVocabulary::MS
                        .const_param_ident("basepeak chromatogram", 1000628)
                        .into(),
                );
                let mut time_array = DataArray::wrap(
                    &ArrayType::TimeArray,
                    BinaryDataArrayType::Float64,
                    to_bytes(&self.time),
                );
                time_array.unit = Unit::Minute;

                let intensity_array = DataArray::wrap(
                    &ArrayType::IntensityArray,
                    BinaryDataArrayType::Float32,
                    to_bytes(&self.intensity),
                );
                arrays.add(time_array);
                arrays.add(intensity_array);
            }
            _ => panic!("Don't know how to construct {:?}", &self.name),
        };
        Chromatogram::new(descr, arrays)
    }
}

impl From<ChromatogramCollector> for Chromatogram {
    fn from(value: ChromatogramCollector) -> Self {
        value.to_chromatogram()
    }
}

/**
An indexed mzML writer that writes [`MultiLayerSpectrum`](crate::spectrum::MultiLayerSpectrum).

Does not buffer spectra in-memory, writing them out immediately but summary chromatogram information
is accumulated.
*/
#[derive(Debug)]
pub struct MzMLWriterType<
    W: Write,
    C: CentroidLike + Default + BuildArrayMapFrom + 'static = CentroidPeak,
    D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom + 'static = DeconvolutedPeak,
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
    pub samples: Vec<Sample>,
    /// The types of data transformations applied to (parts of) the data
    pub data_processings: Vec<DataProcessing>,
    /// The different instrument configurations that were in use during the
    /// data acquisition.
    pub instrument_configurations: HashMap<u32, InstrumentConfiguration>,

    pub state: MzMLWriterState,
    pub write_index: bool,

    pub spectrum_offset_index: OffsetIndex,
    pub chromatogram_offset_index: OffsetIndex,

    pub tic_collector: ChromatogramCollector,
    pub bic_collector: ChromatogramCollector,
    pub wrote_summaries: bool,

    pub run: MassSpectrometryRun,

    handle: InnerXMLWriter<W>,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
    ms_cv: ControlledVocabulary,

    param_groups: Vec<ParamGroup>,
}

impl<
        W: Write,
        C: CentroidLike + Default + BuildArrayMapFrom,
        D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom,
    > SpectrumWriter<C, D> for MzMLWriterType<W, C, D>
{
    fn write<S: SpectrumLike<C, D> + 'static>(&mut self, spectrum: &S) -> io::Result<usize> {
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

    fn close(&mut self) -> io::Result<()> {
        self.close()?;
        Ok(())
    }
}

impl<
        W: Write,
        C: CentroidLike + Default + BuildArrayMapFrom,
        D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom,
        CF: FeatureLike<MZ, IonMobility> + BuildArrayMap3DFrom,
        DF: FeatureLike<Mass, IonMobility> + KnownCharge + BuildArrayMap3DFrom,
    > IonMobilityFrameWriter<CF, DF> for MzMLWriterType<W, C, D>
{
    fn write_frame<S: crate::spectrum::IonMobilityFrameLike<CF, DF> + 'static>(
        &mut self,
        frame: &S,
    ) -> io::Result<usize> {
        let state = frame.description().clone().into();
        let peak_data = match frame.features() {
            crate::spectrum::frame::RefFeatureDataLevel::Missing => BinaryArrayMap::default(),
            crate::spectrum::frame::RefFeatureDataLevel::RawData(a) => a.unstack()?,
            crate::spectrum::frame::RefFeatureDataLevel::Centroid(c) => CF::as_arrays(&c[..]),
            crate::spectrum::frame::RefFeatureDataLevel::Deconvoluted(d) => DF::as_arrays(&d[..]),
        };
        let spectrum = RawSpectrum::new(state, peak_data);
        self.write_owned(spectrum)
    }

    fn write_frame_owned<S: crate::spectrum::IonMobilityFrameLike<CF, DF> + 'static>(
        &mut self,
        frame: S,
    ) -> io::Result<usize> {
        let (features, state) = frame.into_features_and_parts();
        let peak_data = match features {
            crate::spectrum::frame::FeatureDataLevel::Missing => BinaryArrayMap::default(),
            crate::spectrum::frame::FeatureDataLevel::RawData(a) => a.unstack()?,
            crate::spectrum::frame::FeatureDataLevel::Centroid(c) => CF::as_arrays(&c[..]),
            crate::spectrum::frame::FeatureDataLevel::Deconvoluted(d) => DF::as_arrays(&d[..]),
        };
        let spectrum = RawSpectrum::new(state.into(), peak_data);
        self.write_owned(spectrum)
    }

    fn flush_frame(&mut self) -> io::Result<()> {
        self.flush()
    }

    fn close_frames(&mut self) -> io::Result<()> {
        if let Err(e) = self.close() {
            return Err(e.into());
        } else {
            Ok(())
        }
    }
}

impl<
        W: Write,
        C: CentroidLike + Default + BuildArrayMapFrom,
        D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom,
    > MSDataFileMetadata for MzMLWriterType<W, C, D>
{
    crate::impl_metadata_trait!();

    fn spectrum_count_hint(&self) -> Option<u64> {
        Some(self.spectrum_count)
    }

    fn set_spectrum_count_hint(&mut self, value: Option<u64>) {
        match value {
            Some(value) => {
                self.spectrum_count = value;
            },
            None => {},
        };
    }

    fn run_description(&self) -> Option<&MassSpectrometryRun> {
        Some(&self.run)
    }

    fn run_description_mut(&mut self) -> Option<&mut MassSpectrometryRun> {
        Some(&mut self.run)
    }
}

#[derive(Debug, Default, Clone)]
pub struct SpectrumHasSummary {
    pub has_tic: bool,
    pub has_bp: bool,
    pub has_mz_range: bool,
    pub count: usize,
}

impl SpectrumHasSummary {
    pub fn len(&self) -> usize {
        self.count
    }

    pub fn is_empty(&self) -> bool {
        self.count == 0
    }
}

#[derive(Debug, Default, Clone)]
pub struct ParamGroup {
    pub id: String,
    pub params: Vec<Param>,
}

impl ParamGroup {
    pub fn new(id: String, params: Vec<Param>) -> Self {
        Self { id, params }
    }

    pub fn iter(&self) -> impl Iterator<Item = &Param> {
        self.params.iter()
    }

    #[allow(unused)]
    pub fn contains(&self, item: &Param) -> bool {
        self.params.contains(item)
    }

    pub fn covers(&self, params: &[&Param]) -> bool {
        let mut covered = self.params.clone();
        for param in params {
            if let Some(i) = covered.iter().position(|p| **param == *p) {
                covered.remove(i);
            }
        }
        covered.is_empty()
    }

    pub fn filter(&self, params: &mut Vec<&Param>) {
        for p in self.iter() {
            if let Some(i) = params.iter().position(|p2| *p == **p2) {
                params.remove(i);
            }
        }
    }
}

impl IntoIterator for ParamGroup {
    type Item = Param;

    type IntoIter = <Vec<Param> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.params.into_iter()
    }
}

impl_param_described!(ParamGroup);

impl<W: Write, C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default>
    MzMLWriterType<W, C, D>
where
    C: BuildArrayMapFrom,
    D: BuildArrayMapFrom,
{
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_cv_metadata.py', 'data-version']).decode('utf8').strip()
    cog.outl(f'const PSIMS_VERSION: &\'static str = "{buf}";')
    ]]]*/
    const PSIMS_VERSION: &'static str = "4.1.174";
    //[[[end]]] (checksum: c20f130bf2cb029cf820fb56ecf3075c)
    const UNIT_VERSION: &'static str = "releases/2020-03-10";

    pub const fn get_indent_size() -> u64 {
        InnerXMLWriter::<W>::INDENT_SIZE
    }

    pub fn new_with_index_and_compression(
        file: W,
        write_index: bool,
        data_array_compression: BinaryCompressionType,
    ) -> MzMLWriterType<W, C, D> {
        let handle = InnerXMLWriter::new(file);
        let data_array_compression = match data_array_compression {
            BinaryCompressionType::Decoded => {
                warn!("The mzML writer was asked to use the `Decoded` array compression, using `Zlib` instead");
                BinaryCompressionType::Zlib
            }
            _ => data_array_compression,
        };
        MzMLWriterType {
            handle,
            file_description: FileDescription::default(),
            instrument_configurations: HashMap::new(),
            softwares: Vec::new(),
            samples: Vec::new(),
            data_processings: Vec::new(),
            offset: 0,
            spectrum_offset_index: OffsetIndex::new("spectrum".into()),
            chromatogram_offset_index: OffsetIndex::new("chromatogram".into()),
            state: MzMLWriterState::Start,
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
            write_index,
            spectrum_count: 0,
            spectrum_counter: 0,
            chromatogram_count: 2,
            chromatogram_counter: 0,
            tic_collector: ChromatogramCollector::of(ChromatogramType::TotalIonCurrentChromatogram),
            bic_collector: ChromatogramCollector::of(ChromatogramType::BasePeakChromatogram),
            ms_cv: ControlledVocabulary::MS,
            data_array_compression,
            wrote_summaries: false,
            run: MassSpectrometryRun::default(),
            param_groups: Vec::default(),
        }
    }

    pub fn new_with_index(file: W, write_index: bool) -> MzMLWriterType<W, C, D> {
        Self::new_with_index_and_compression(file, write_index, BinaryCompressionType::Zlib)
    }

    /// Wrap a new [`std::io::Write`]-able type, constructing a new [`MzMLWriterType`]
    pub fn new(file: W) -> MzMLWriterType<W, C, D> {
        Self::new_with_index(file, true)
    }

    fn transition_err(&self, to_state: MzMLWriterState) -> WriterResult {
        Err(MzMLWriterError::StateTransitionError {
            from_state: self.state,
            to_state,
        })
    }

    /// Imitate the [`io::Seek`] method using an internal byte counter
    pub fn stream_position(&mut self) -> io::Result<u64> {
        Ok(self.handle.handle.get_mut().bytes_written())
    }

    fn count_params<'a>(&self, entry: &'a impl ParamDescribed) -> HashMap<&'a Param, usize> {
        let mut acc = HashMap::new();
        for p in entry.params().iter() {
            *acc.entry(p).or_default() += 1;
        }
        acc
    }

    fn find_repeated_params<'a, P: ParamDescribed + 'a>(
        &self,
        entries: impl Iterator<Item = &'a P>,
        tag: &str,
    ) -> Option<ParamGroup> {
        let mut limit_counters: HashMap<_, (usize, usize)> = HashMap::new();
        let mut total_counters: HashMap<_, usize> = HashMap::new();

        let blocks: Vec<_> = entries.map(|entry| self.count_params(entry)).collect();
        for block in blocks.iter() {
            for (k, v) in block.iter() {
                let b = limit_counters.entry(*k).or_insert((usize::MAX, usize::MIN));
                *b = (*v.min(&b.0), *v.max(&b.1));
                *total_counters.entry(*k).or_default() += *v;
            }
        }

        let mut candidates = Vec::new();
        for (param, (min_count, _max_count)) in limit_counters.iter() {
            if *total_counters.get(param).unwrap() > *min_count {
                candidates.push((*param, *min_count))
            }
        }

        let spanned_blocks: HashSet<_> = blocks
            .iter()
            .map(|block| {
                let members: Vec<_> = candidates
                    .iter()
                    .filter(|(p, count)| *block.get(p).unwrap_or(&0) >= *count)
                    .collect();
                members
            })
            .collect();

        if let Some(block) = spanned_blocks.iter().min_by_key(|block| block.len()) {
            let members: Vec<_> = block
                .iter()
                .map(|(p, s)| {
                    let mut v = Vec::with_capacity(*s);
                    for _ in 0..*s {
                        v.push((*p).clone());
                    }
                    v.into_iter()
                })
                .flatten()
                .collect();
            if members.len() < 2 {
                None
            } else {
                let id = format!("{}_params_shared", tag);
                Some(ParamGroup::new(id, members))
            }
        } else {
            None
        }
    }

    fn analyze_header_params(&mut self) {
        self.param_groups.extend(self.find_repeated_params(
            self.instrument_configurations.values(),
            "instrument_configuration_",
        ));
    }

    fn make_psi_ms_cv(&self) -> BytesStart<'static> {
        let mut cv = BytesStart::from_content("cv", 2);
        cv.push_attribute(("id", "MS"));
        cv.push_attribute(("fullName", "PSI-MS"));
        cv.push_attribute(("URI", "http://purl.obolibrary.org/obo/ms.obo"));
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

        if self.write_index {
            let mut indexed = BytesStart::from_content("indexedmzML", 11);
            indexed.push_attribute(("xmlns", "http://psi.hupo.org/ms/mzml"));
            indexed.push_attribute(("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance"));
            indexed.push_attribute((
                "xsi:schemaLocation",
                "http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.3_idx.xsd",
            ));
            self.handle.write_event(Event::Start(indexed))?;
        }

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
            return self.transition_err(MzMLWriterState::Header);
        }
        self.analyze_header_params();
        self.write_cv_list()?;
        self.write_file_description()?;
        self.write_referenceable_param_group_list()?;
        self.write_sample_list()?;
        self.write_software_list()?;
        self.write_instrument_configuration()?;
        self.write_data_processing()?;

        self.state = MzMLWriterState::Header;
        Ok(())
    }

    fn write_file_description(&mut self) -> WriterResult {
        let fd = bstart!("fileDescription");
        start_event!(self, fd);

        let fc_tag = bstart!("fileContent");
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
            if let Some(param) = sf.file_format.as_ref() {
                self.handle.write_param(param)?
            }
            if let Some(param) = sf.id_format.as_ref() {
                self.handle.write_param(param)?
            }
            for param in sf.params() {
                self.handle.write_param(param)?
            }
            self.handle.write_event(Event::End(tag.to_end()))?;
        }
        self.handle.write_event(Event::End(outer.to_end()))?;

        end_event!(self, fd);
        Ok(())
    }

    fn write_referenceable_param_group_list(&mut self) -> WriterResult {
        if !self.param_groups.is_empty() {
            let mut list = bstart!("referenceableParamGroupList");
            let count = self.param_groups.len().to_string();
            attrib!("count", count, list);
            start_event!(self, list);
            for param_group in self.param_groups.iter() {
                let mut group_tag = bstart!("referenceableParamGroup");
                attrib!("id", param_group.id, group_tag);
                start_event!(self, group_tag);
                for param in param_group.iter() {
                    self.handle.write_param(param)?;
                }
                end_event!(self, group_tag);
            }
            end_event!(self, list);
        }
        Ok(())
    }

    fn write_sample_list(&mut self) -> WriterResult {
        if self.samples.is_empty() {
            return Ok(())
        }
        let mut outer = bstart!("sampleList");
        let count = self.samples.len().to_string();
        attrib!("count", count, outer);
        self.handle.write_event(Event::Start(outer.borrow()))?;
        for sample in self.samples.iter() {
            let mut tag = bstart!("sample");
            attrib!("id", sample.id, tag);
            if let Some(name) = sample.name.as_ref() {
                attrib!("name", name, tag);
            }
            self.handle.write_event(Event::Start(tag.borrow()))?;
            for param in sample.params() {
                self.handle.write_param(param)?
            }
            self.handle.write_event(Event::End(tag.to_end()))?;
        }
        self.handle.write_event(Event::End(outer.to_end()))?;
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
            self.handle
                .write_param_list_with_refs(ic.params().iter(), &self.param_groups)?;
            let mut comp_list_tag = bstart!("componentList");
            let comp_count = ic.components.len().to_string();
            attrib!("count", comp_count, comp_list_tag);
            self.handle
                .write_event(Event::Start(comp_list_tag.borrow()))?;
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
            self.handle
                .write_event(Event::End(comp_list_tag.to_end()))?;
            let mut sw = bstart!("softwareRef");
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
            let inst_id = instrument_id(ic_ref);
            attrib!("defaultInstrumentConfigurationRef", inst_id, run);
        }
        if let Some(sf_ref) = self.file_description.source_files.first() {
            attrib!("defaultSourceFileRef", sf_ref.id, run);
        };
        self.handle.write_event(Event::Start(run))?;
        self.state = MzMLWriterState::Run;
        Ok(())
    }

    pub fn write_param<P: ParamLike>(&mut self, param: &P) -> WriterResult {
        self.handle.write_param(param)
    }

    pub fn write_param_list<'a, P: ParamLike + 'a, T: Iterator<Item = &'a P>>(
        &mut self,
        params: T,
    ) -> WriterResult {
        self.handle.write_param_list(params)
    }

    pub fn write_param_list_with_refs<'a, T: Iterator<Item = &'a Param>>(
        &mut self,
        params: T,
    ) -> WriterResult {
        let mut params: Vec<_> = params.collect();
        let mut ref_tags = Vec::new();
        for param_group in self.param_groups.iter() {
            if param_group.covers(params.as_ref()) {
                param_group.filter(&mut params);
                ref_tags.push(param_group);
            }
        }

        for p in params {
            self.handle.write_param(p)?;
        }

        for pg in ref_tags {
            self.handle.write_reference_group(pg)?;
        }
        Ok(())
    }

    pub const fn get_ms_cv(&self) -> &ControlledVocabulary {
        &self.ms_cv
    }
}

impl<W: Write, C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default>
    MzMLWriterType<W, C, D>
where
    C: BuildArrayMapFrom,
    D: BuildArrayMapFrom,
{
    pub fn start_spectrum_list(&mut self) -> WriterResult {
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
        if self.state < MzMLWriterState::SpectrumList {
            self.start_spectrum_list()?;
        }
        let tag = bstart!("spectrumList");
        end_event!(self, tag);
        self.state = MzMLWriterState::SpectrumListClosed;
        Ok(())
    }

    pub fn write_summary_chromatograms(&mut self) -> WriterResult {
        if !self.wrote_summaries {
            self.write_chromatogram(&self.tic_collector.to_chromatogram())?;
            self.write_chromatogram(&self.bic_collector.to_chromatogram())?;
            self.wrote_summaries = true;
        }
        Ok(())
    }

    pub fn start_chromatogram_list(&mut self) -> WriterResult {
        match self.state {
            MzMLWriterState::ChromatogramList => {}
            MzMLWriterState::SpectrumList => {
                self.close_spectrum_list()?;
            }
            state if state < MzMLWriterState::SpectrumList => {
                self.start_run()?;
                self.start_spectrum_list()?;
                self.close_spectrum_list()?;
            }
            state if state > MzMLWriterState::ChromatogramList => {
                return self.transition_err(MzMLWriterState::ChromatogramList);
            }
            _ => {}
        }
        let mut list = bstart!("chromatogramList");

        let count = self.chromatogram_count.to_string();
        attrib!("count", count, list);
        if let Some(dp) = self.data_processings.first() {
            attrib!("defaultDataProcessingRef", dp.id, list);
        }
        self.handle.write_event(Event::Start(list))?;
        self.state = MzMLWriterState::ChromatogramList;
        Ok(())
    }

    fn close_chromatogram_list(&mut self) -> WriterResult {
        if self.state < MzMLWriterState::ChromatogramList {
            self.start_chromatogram_list()?;
        }
        self.write_summary_chromatograms()?;
        let tag = bstart!("chromatogramList");
        end_event!(self, tag);
        self.state = MzMLWriterState::SpectrumListClosed;
        Ok(())
    }

    fn close_run(&mut self) -> WriterResult {
        if self.state < MzMLWriterState::Run {
            self.start_run()?;
        } else if self.state == MzMLWriterState::SpectrumList {
            self.close_spectrum_list()?;
            self.start_chromatogram_list()?;
            self.close_chromatogram_list()?;
        } else if self.state == MzMLWriterState::ChromatogramList {
            self.close_chromatogram_list()?;
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
        let tag = BytesEnd::new("mzML");
        self.write_event(Event::End(tag))?;
        self.state = MzMLWriterState::MzMLClosed;
        Ok(())
    }

    fn close_indexed_mzml(&mut self) -> WriterResult {
        if self.state < MzMLWriterState::MzMLClosed {
            self.close_mzml()?;
        }
        if self.write_index {
            self.state = MzMLWriterState::IndexList;
            self.write_index_list()?;
            let tag = bstart!("indexedmzML");
            end_event!(self, tag);
        }
        self.state = MzMLWriterState::End;
        Ok(())
    }

    /**
    Close the wrapping `<indexedmzML>` document, which will trigger writing
    out the offset indices and file checksum at the tail of the document.
    */
    pub fn close(&mut self) -> WriterResult {
        if self.state < MzMLWriterState::End {
            self.close_indexed_mzml()?;
            self.handle.flush()?;
            Ok(())
        } else {
            Ok(())
        }
    }
}

impl<W: Write, C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default>
    MzMLWriterType<W, C, D>
where
    C: BuildArrayMapFrom,
    D: BuildArrayMapFrom,
{
    pub fn write_scan_list(&mut self, acq: &Acquisition) -> WriterResult {
        let mut scan_list_tag = bstart!("scanList");
        let count = acq.scans.len().to_string();
        attrib!("count", count, scan_list_tag);
        start_event!(self, scan_list_tag);
        self.handle.write_param(&acq.combination.to_param())?;

        for scan in acq.scans.iter() {
            let mut scan_tag = bstart!("scan");
            let had_ref = if let Some(sref) = scan.spectrum_reference.as_ref() {
                attrib!("spectrumRef", sref, scan_tag);
                true
            } else {
                false
            };

            let id = instrument_id(&scan.instrument_configuration_id);
            attrib!("instrumentConfigurationRef", id, scan_tag);
            self.handle.write_event(Event::Start(scan_tag.borrow()))?;

            // If these are zero and we had a spectrumRef attribute, this must be an empty reference
            if !(had_ref && scan.start_time == 0.0 && scan.injection_time == 0.0) {
                self.handle.write_param(&self.ms_cv.const_param(
                    "scan start time",
                    ValueRef::Float(scan.start_time),
                    1000016,
                    Unit::Minute,
                ))?;

                self.handle.write_param(&self.ms_cv.const_param(
                    "ion injection time",
                    ValueRef::Float(scan.injection_time as f64),
                    1000927,
                    Unit::Millisecond,
                ))?;
            }

            for param in scan.params() {
                self.handle.write_param(param)?
            }

            if !scan.scan_windows.is_empty() {
                let mut scan_window_list_tag = bstart!("scanWindowList");
                let scan_window_list_count = scan.scan_windows.len().to_string();

                attrib!("count", scan_window_list_count, scan_window_list_tag);
                self.handle
                    .write_event(Event::Start(scan_window_list_tag.borrow()))?;
                for window in scan.scan_windows.iter() {
                    let window_tag = bstart!("scanWindow");
                    self.handle.write_event(Event::Start(window_tag.borrow()))?;

                    self.handle.write_param(&self.ms_cv.const_param(
                        "scan window lower limit",
                        ValueRef::Float(window.lower_bound as f64),
                        1000501,
                        Unit::MZ,
                    ))?;

                    self.handle.write_param(&self.ms_cv.const_param(
                        "scan window upper limit",
                        ValueRef::Float(window.upper_bound as f64),
                        1000500,
                        Unit::MZ,
                    ))?;

                    self.handle.write_event(Event::End(window_tag.to_end()))?;
                }
                self.handle
                    .write_event(Event::End(scan_window_list_tag.to_end()))?;
            }
            self.handle.write_event(Event::End(scan_tag.to_end()))?;
        }
        end_event!(self, scan_list_tag);
        Ok(())
    }

    pub fn write_isolation_window(&mut self, iw: &IsolationWindow) -> WriterResult {
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

    pub fn write_selected_ions(&mut self, precursor: &impl PrecursorSelection) -> WriterResult {
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
        self.handle.write_param_list(ion.params().iter())?;
        end_event!(self, tag);
        end_event!(self, outer);
        Ok(())
    }

    pub fn write_activation(&mut self, precursor: &impl PrecursorSelection) -> WriterResult {
        let act = precursor.activation();
        let tag = bstart!("activation");
        start_event!(self, tag);
        if let Some(meth) = act.method() {
            let meth_param: Param = meth.clone().into();
            self.handle.write_param(&meth_param)?;
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

    pub fn write_precursor(&mut self, precursor: &impl PrecursorSelection) -> WriterResult {
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

    fn write_ms_level<
        C1: CentroidLike + Default + BuildArrayMapFrom,
        D1: DeconvolutedCentroidLike + Default + BuildArrayMapFrom,
        S: SpectrumLike<C1, D1> + 'static
    >(&mut self, spectrum: &S) -> WriterResult {
        let ms_level = spectrum.ms_level();
        if ms_level == 1 {
            self.handle.write_param(&MS1_SPECTRUM)?;
        } else if ms_level > 1 {
            self.handle.write_param(&MSN_SPECTRUM)?;
        }
        self.write_param(&self.ms_cv.const_param(
            "ms level",
            ValueRef::Int(ms_level as i64),
            1000511,
            Unit::Unknown,
        ))?;
        Ok(())
    }

    fn write_polarity<
        C1: CentroidLike + Default + BuildArrayMapFrom,
        D1: DeconvolutedCentroidLike + Default + BuildArrayMapFrom,
        S: SpectrumLike<C1, D1> + 'static
    >(&mut self, spectrum: &S) -> WriterResult {
        match spectrum.polarity() {
            ScanPolarity::Negative => self.handle.write_param(&NEGATIVE_SCAN),
            ScanPolarity::Positive => self.handle.write_param(&POSITIVE_SCAN),
            ScanPolarity::Unknown => {
                warn!(
                    "Could not determine scan polarity for {}, assuming positive",
                    spectrum.id()
                );
                self.handle.write_param(&POSITIVE_SCAN)
            }
        }
    }

    fn write_continuity<
        C1: CentroidLike + Default + BuildArrayMapFrom,
        D1: DeconvolutedCentroidLike + Default + BuildArrayMapFrom,
        S: SpectrumLike<C1, D1> + 'static
    >(&mut self, spectrum: &S) -> WriterResult {
        match spectrum.signal_continuity() {
            SignalContinuity::Profile => self.handle.write_param(&PROFILE_SPECTRUM),
            SignalContinuity::Unknown => {
                warn!(
                    "Could not determine scan polarity for {}, assuming centroid",
                    spectrum.id()
                );
                self.handle.write_param(&CENTROID_SPECTRUM)
            }
            _ => self.handle.write_param(&CENTROID_SPECTRUM),
        }
    }

    fn write_signal_properties<
        C1: CentroidLike + Default + BuildArrayMapFrom,
        D1: DeconvolutedCentroidLike + Default + BuildArrayMapFrom,
        S: SpectrumLike<C1, D1> + 'static
    >(
        &mut self,
        spectrum: &S,
    ) -> WriterResult {
        if spectrum.ms_level() > 0 {
            self.handle.write_param_list(
                spectrum
                    .params()
                    .iter()
                    .filter(|p| **p != MS1_SPECTRUM && **p != MSN_SPECTRUM),
            )?
        } else {
            self.handle.write_param_list(spectrum.params().iter())?
        }
        Ok(())
    }
}

impl<W: Write, C: CentroidLike + Default, D: DeconvolutedCentroidLike + Default>
    MzMLWriterType<W, C, D>
where
    C: BuildArrayMapFrom,
    D: BuildArrayMapFrom,
{
    /// Write a `binaryDataArray` from a [`DataArray`] whose contents have been translated into a base64
    /// encoded string.
    ///
    /// Unless the `array` has already been encoded, use [`Self::write_binary_data_array`] instead.
    ///
    /// # Panics
    ///
    /// Panics if `array.dtype()` is [`BinaryDataArrayType::Unknown`] as these cannot be
    /// encoded correctly, or if `array.name` cannot be converted to a `cvParam` or `userParam`.
    ///
    /// # Errors
    /// This function will return an error if a [`MzMLWriterError`] error occurs during
    /// writing any underlying data occurs.
    /// .
    pub fn write_binary_data_array_pre_encoded(
        &mut self,
        array: &DataArray,
        default_array_len: usize,
        encoded_array: &[u8],
    ) -> WriterResult {
        let mut outer = bstart!("binaryDataArray");

        let encoded_len = encoded_array.len().to_string();
        attrib!("encodedLength", encoded_len, outer);
        let array_len = array.data_len()?;
        if array_len != default_array_len {
            let array_len = array_len.to_string();
            attrib!("arrayLength", array_len, outer);
        }

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

        self.handle.write_param(
            self.data_array_compression
                .clone()
                .as_param()
                .as_ref()
                .unwrap(),
        )?;

        match &array.name {
            ArrayType::MZArray | ArrayType::IntensityArray | ArrayType::ChargeArray => {
                self.handle.write_param(&array.name.as_param_const())?
            }
            ArrayType::TimeArray
            | ArrayType::RawIonMobilityArray
            | ArrayType::MeanIonMobilityArray
            | ArrayType::DeconvolutedIonMobilityArray => self
                .handle
                .write_param(&array.name.as_param_with_unit_const(array.unit))?,
            ArrayType::NonStandardDataArray { name } => {
                let mut p =
                    self.ms_cv
                        .param_val("MS:1000786", "non-standard data array", name.as_str());
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
            String::from_utf8_lossy(encoded_array).as_ref(),
        )))?;
        end_event!(self, bin);
        end_event!(self, outer);
        Ok(())
    }

    /// Write a `binaryDataArray` from a [`DataArray`].
    ///
    /// # Errors
    ///
    /// This function will return an error if a [`MzMLWriterError`] error occurs during
    /// writing any underlying data occurs.
    pub fn write_binary_data_array(
        &mut self,
        array: &DataArray,
        default_array_len: usize,
    ) -> WriterResult {
        let encoded_array = array.encode_bytestring(self.data_array_compression);
        self.write_binary_data_array_pre_encoded(array, default_array_len, &encoded_array)
    }

    pub fn write_binary_data_arrays(
        &mut self,
        arrays: &BinaryArrayMap,
        default_array_len: usize,
    ) -> WriterResult {
        let count = arrays.len().to_string();
        let mut outer = bstart!("binaryDataArrayList");
        attrib!("count", count, outer);
        start_event!(self, outer);
        #[cfg(feature = "parallelism")]
        {
            let compression = self.data_array_compression;
            let mut array_pairs: Vec<(&ArrayType, &DataArray, Vec<u8>)> =
                if arrays.len() < PARALLEL_COMPRESSION_FAN {
                    arrays
                        .iter()
                        .map(|(t, d)| {
                            let encoded = d.encode_bytestring(compression);
                            (t, d, encoded)
                        })
                        .collect()
                } else {
                    arrays
                        .par_iter()
                        .map(|(t, d)| {
                            let encoded = d.encode_bytestring(compression);
                            (t, d, encoded)
                        })
                        .collect()
                };
            array_pairs.sort_by_key(|f| f.0);
            for (_tp, array, encoded) in array_pairs {
                self.write_binary_data_array_pre_encoded(array, default_array_len, &encoded)?
            }
        }

        #[cfg(not(feature = "parallelism"))]
        {
            let mut array_pairs: Vec<(&ArrayType, &DataArray)> = arrays.iter().collect();
            array_pairs.sort_by_key(|f| f.0);
            for (_tp, array) in array_pairs {
                self.write_binary_data_array(array, default_array_len)?
            }
        }
        end_event!(self, outer);
        Ok(())
    }

    /// Write the opening tag for a `spectrum` tag and write it out.
    ///
    /// *NOTE*: This function isn't useful unless you are modifying the writing of
    /// the spectrum content but not the index management and tag level
    /// information. Instead use [`Self::write_spectrum`].
    ///
    /// # See also
    /// [`Self::write_spectrum`]
    ///
    /// # Errors
    ///
    /// This function will return an error if a [`MzMLWriterError`] error occurs during
    /// writing any underlying data occurs.
    pub fn start_spectrum<
        C1: CentroidLike + Default + BuildArrayMapFrom,
        D1: DeconvolutedCentroidLike + Default + BuildArrayMapFrom,
        S: SpectrumLike<C1, D1> + 'static,
    >(
        &mut self,
        spectrum: &S,
        outer: &mut BytesStart,
        summary_metrics: &SpectrumHasSummary,
    ) -> Result<usize, MzMLWriterError> {
        attrib!("id", spectrum.id(), outer);
        let count = self.spectrum_counter.to_string();
        attrib!("index", count, outer);
        let default_array_len_u = summary_metrics.count;
        let default_array_len = default_array_len_u.to_string();

        attrib!("defaultArrayLength", default_array_len, outer);

        self.handle.write_event(Event::Start(outer.borrow()))?;
        self.spectrum_counter += 1;
        Ok(default_array_len_u)
    }

    /// Checks if spectrum-level summaries are already calculated for
    /// `spectrum`.
    pub fn spectrum_has_summaries<
        C1: CentroidLike + Default + BuildArrayMapFrom,
        D1: DeconvolutedCentroidLike + Default + BuildArrayMapFrom,
        S: SpectrumLike<C1, D1> + 'static
    >(
        &self,
        spectrum: &S,
    ) -> SpectrumHasSummary {
        let peaks = spectrum.peaks();
        let peak_count = match peaks {
            RefPeakDataLevel::RawData(arrays) => {
                if let Some(arr) = arrays.get(&ArrayType::MZArray) {
                    if let Ok(count) = arr.data_len() {
                        count
                    } else {
                        0
                    }
                } else {
                    0
                }
            }
            _ => peaks.len(),
        };

        let mut summary = SpectrumHasSummary {
            count: peak_count,
            ..Default::default()
        };
        for p in spectrum.params().iter().filter(|p| p.is_ms()) {
            let a = p.accession.unwrap();
            if a == TIC_TERM.accession {
                summary.has_tic = true;
            } else if a == BPMZ_TERM.accession {
                summary.has_bp = true;
            }
        }
        summary
    }

    /// Write spectrum-level descriptive metadata, acquisition scan metadata,
    /// and precursors if any are present.
    pub fn write_spectrum_descriptors<
        C1: CentroidLike + Default + BuildArrayMapFrom,
        D1: DeconvolutedCentroidLike + Default + BuildArrayMapFrom,
        S: SpectrumLike<C1, D1> + 'static
    >(
        &mut self,
        spectrum: &S,
        _summary_metrics: &SpectrumHasSummary,
    ) -> WriterResult {
        self.write_ms_level(spectrum)?;
        self.write_polarity(spectrum)?;
        self.write_continuity(spectrum)?;
        self.write_signal_properties(spectrum)?;

        self.write_scan_list(spectrum.acquisition())?;
        for precursor in spectrum.precursor_iter() {
            self.write_precursor(precursor)?;
        }
        Ok(())
    }

    /**
    Write a spectrum  out to the mzML file, encoding the highest procressing degree peak data present.

    ## Side-Effects
    If the writer has not already started writing the spectra, this will cause all the metadata
    to be written out and the `<spectrumList>` element will be opened, preventing new metadata
    from being written to this stream. Furthermore, this writes the spectrum count out, so the value
    may no longer be changed.

    # Errors
    This function will return an error if a [`MzMLWriterError`] error occurs during
    writing any underlying data occurs.
    */
    pub fn write_spectrum<
        C1: CentroidLike + Default + BuildArrayMapFrom,
        D1: DeconvolutedCentroidLike + Default + BuildArrayMapFrom,
        S: SpectrumLike<C1, D1> + 'static
    >(
        &mut self,
        spectrum: &S,
    ) -> WriterResult {
        match self.state {
            MzMLWriterState::SpectrumList => {}
            state if state < MzMLWriterState::SpectrumList => {
                self.start_spectrum_list()?;
            }
            state if state > MzMLWriterState::SpectrumList => {
                // Cannot write spectrum, currently in state which happens
                // after spectra may be written
                return Err(MzMLWriterError::InvalidActionError(self.state));
            }
            _ => {}
        }
        let pos = self.stream_position()?;
        self.spectrum_offset_index
            .insert(spectrum.id().to_string(), pos);

        let summary_metrics = self.spectrum_has_summaries(spectrum);

        let mut outer = bstart!("spectrum");
        let default_array_len = self.start_spectrum(spectrum, &mut outer, &summary_metrics)?;

        self.write_spectrum_descriptors(spectrum, &summary_metrics)?;

        let tic = spectrum
            .params()
            .get_param_by_curie(&curie!(MS:1000285))
            .map(|p| p.to_f32().unwrap())
            .unwrap_or_else(|| spectrum.peaks().tic());

        let bpi = spectrum
            .params()
            .get_param_by_curie(&curie!(MS:1000505))
            .map(|p| p.to_f32().unwrap())
            .unwrap_or_else(|| spectrum.peaks().base_peak().intensity);

        let time = spectrum.start_time();

        self.tic_collector.add(time, tic);
        self.bic_collector.add(time, bpi);

        match spectrum.peaks() {
            RefPeakDataLevel::RawData(arrays) => {
                self.write_binary_data_arrays(arrays, default_array_len)?;
            }
            RefPeakDataLevel::Centroid(arrays) => {
                self.write_binary_data_arrays(&C1::as_arrays(&arrays[0..]), default_array_len)?
            }
            RefPeakDataLevel::Deconvoluted(arrays) => {
                self.write_binary_data_arrays(&D1::as_arrays(&arrays[0..]), default_array_len)?
            }
            RefPeakDataLevel::Missing => {
                self.write_binary_data_arrays(&BinaryArrayMap::new(), default_array_len)?
            }
        }

        end_event!(self, outer);
        Ok(())
    }

    /// Write the opening tag for a `chromatogram` tag and write it out.
    ///
    /// *NOTE*: This function isn't useful unless you are modifying the writing of
    /// the chromatogram content but not the index management and tag level
    /// information. Instead use [`Self::write_chromatogram`].
    ///
    /// # See also
    /// [`Self::write_chromatogram`]
    ///
    /// # Errors
    ///
    /// This function will return an error if a [`MzMLWriterError`] error occurs during
    /// writing any underlying data occurs.
    pub fn start_chromatogram(
        &mut self,
        chromatogram: &Chromatogram,
        outer: &mut BytesStart,
    ) -> Result<usize, MzMLWriterError> {
        attrib!("id", chromatogram.id(), outer);
        let count = self.chromatogram_counter.to_string();
        attrib!("index", count, outer);
        let default_array_len_u = chromatogram
            .arrays
            .get(&ArrayType::TimeArray)
            .unwrap()
            .data_len()?;
        let default_array_len = default_array_len_u.to_string();

        attrib!("defaultArrayLength", default_array_len, outer);

        self.handle.write_event(Event::Start(outer.borrow()))?;
        self.chromatogram_counter += 1;
        Ok(default_array_len_u)
    }

    /**
    Write a chromatogram  out to the mzML file.

    ## Side-Effects
    If the writer has not already started writing the spectra, this will cause all the metadata
    to be written out and the `<chromatogramList>` element will be opened, preventing new metadata
    or spectra from being written to this stream. Furthermore, this writes the chromatogram count out,
    so the value may no longer be changed.

    # Errors
    This function will return an error if a [`MzMLWriterError`] error occurs during
    writing any underlying data occurs.
    */
    pub fn write_chromatogram(&mut self, chromatogram: &Chromatogram) -> WriterResult {
        match self.state {
            MzMLWriterState::ChromatogramList => {}
            state if state < MzMLWriterState::ChromatogramList => {
                self.start_chromatogram_list()?;
            }
            state if state > MzMLWriterState::ChromatogramList => {
                // Cannot write chromatogram, currently in state which happens
                // after spectra may be written
                return Err(MzMLWriterError::InvalidActionError(self.state));
            }
            _ => {}
        }
        let pos = self.stream_position()?;
        self.chromatogram_offset_index
            .insert(chromatogram.id().to_string(), pos);

        let mut outer = bstart!("chromatogram");
        let default_array_len = self.start_chromatogram(chromatogram, &mut outer)?;
        self.write_param_list(chromatogram.params().iter())?;

        if let Some(precursor) = chromatogram.precursor() {
            self.write_precursor(precursor)?;
        }

        self.write_binary_data_arrays(&chromatogram.arrays, default_array_len)?;
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
        if !self.write_index {
            return Ok(());
        }
        if self.state < MzMLWriterState::IndexList {
            self.close_mzml()?;
        }
        let offset = self.stream_position()?;
        let mut outer = bstart!("indexList");
        attrib!("count", "2", outer);
        start_event!(self, outer);
        let mut tmp_index = OffsetIndex::new("tmp".to_string());
        mem::swap(&mut tmp_index, &mut self.spectrum_offset_index);
        self.write_index(&tmp_index)?;
        mem::swap(&mut tmp_index, &mut self.spectrum_offset_index);

        mem::swap(&mut tmp_index, &mut self.chromatogram_offset_index);
        self.write_index(&tmp_index)?;
        mem::swap(&mut tmp_index, &mut self.chromatogram_offset_index);

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

    pub fn get_mut(&mut self) -> io::Result<&mut W> {
        let inner = self.handle.handle.get_mut();
        let inner = inner.get_mut();
        Ok(inner)
    }

    pub fn write_event(&mut self, event: Event) -> WriterResult {
        self.handle.write_event(event)
    }
}

impl<
        W: io::Write,
        C: CentroidLike + Default + BuildArrayMapFrom,
        D: DeconvolutedCentroidLike + Default + BuildArrayMapFrom,
    > Drop for MzMLWriterType<W, C, D>
{
    fn drop(&mut self) {
        MzMLWriterType::close(self).unwrap();
    }
}

/// A specialization of [`MzMLWriterType`] for the default peak types, for common use.
pub type MzMLWriter<W> = MzMLWriterType<W, CentroidPeak, DeconvolutedPeak>;

#[cfg(test)]
mod test {
    use super::super::reader::MzMLReader;
    use super::*;
    use crate::prelude::*;
    use std::fs;
    use std::path;
    use tempfile;

    #[test_log::test]
    fn write_test() -> WriterResult {
        let tmpdir = tempfile::tempdir()?;
        let dest_path = tmpdir.path().join("duplicate.mzML");
        let path = path::Path::new("./test/data/small.mzML");
        let mut reader = MzMLReader::<_>::open_path(path).expect("Test file doesn't exist?");

        let n = reader.len();
        assert_eq!(n, 48);

        let dest = fs::File::create(dest_path.clone())?;
        let mut writer = MzMLWriterType::new(dest);
        writer.copy_metadata_from(&reader);
        assert_eq!(reader.samples(), writer.samples());
        *writer.spectrum_count_mut() = reader.len() as u64;
        for mut group in reader.groups() {
            if let Some(prec) = group.precursor_mut() {
                if let Some(arrays) = prec.arrays.as_mut() {
                    arrays.iter_mut().for_each(|(_, arr)| {
                        arr.store_compressed(BinaryCompressionType::Zlib).unwrap();
                    });
                }
                writer.write(prec)?;
            }
            for prod in group.products.iter() {
                writer.write(prod)?;
            }
        }
        writer.close()?;

        let mut reader2 = MzMLReader::open_path(dest_path)?;
        assert_eq!(reader.file_description(), reader2.file_description());

        for (a, b) in reader.iter().zip(reader2.iter()) {
            assert_eq!(a.id(), b.id());
            assert_eq!(a.ms_level(), b.ms_level());
            assert_eq!(a.index(), b.index());
            for (x, y) in a
                .arrays
                .as_ref()
                .unwrap()
                .mzs()
                .unwrap()
                .iter()
                .zip(b.arrays.unwrap().mzs().unwrap().iter())
            {
                assert!(
                    (x - y).abs() < 1e-3,
                    "{}: {} - {} = {}",
                    a.id(),
                    x,
                    y,
                    x - y
                )
            }
        }

        Ok(())
    }
}
