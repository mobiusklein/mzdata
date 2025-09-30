use std::{
    collections::HashMap,
    convert::TryInto,
    fs,
    io::{self, BufReader, Read, Seek, SeekFrom},
    marker::PhantomData,
    mem,
};

use log::{debug, trace, warn};

use mzpeaks::{CentroidLike, CentroidPeak, DeconvolutedPeak};
use quick_xml::{
    Error as XMLError, Reader,
    events::{BytesEnd, BytesStart, BytesText, Event},
};

use crate::{
    io::{Generic3DIonMobilityFrameSource, IntoIonMobilityFrameSource, utils::DetailLevel},
    meta::{
        DataProcessing, DissociationEnergyTerm, FileDescription, InstrumentConfiguration,
        MSDataFileMetadata, MassSpectrometryRun, Sample, ScanSettings, Software,
    },
    params::{Param, ParamList, Unit},
    prelude::{ParamLike, *},
    spectrum::{
        HasIonMobility,
        bindata::{
            ArrayType, BinaryArrayMap, BinaryCompressionType, BinaryDataArrayType,
            BuildArrayMapFrom, BuildFromArrayMap, DataArray,
        },
        chromatogram::{Chromatogram, ChromatogramLike},
        scan_properties::*,
        spectrum_types::{CentroidSpectrumType, MultiLayerSpectrum, RawSpectrum, Spectrum},
    },
};

use super::{
    super::{
        offset_index::OffsetIndex,
        traits::{
            ChromatogramSource, MZFileReader, RandomAccessSpectrumIterator, SeekRead,
            SpectrumAccessError, SpectrumSource,
        },
    },
    reading_shared::{
        CVParamParse, EntryType, FileMetadataBuilder, IncrementingIdMap, IndexParserState,
        IndexedMzMLIndexExtractor, MzMLIndexingError, MzMLParserError, MzMLParserState, MzMLSAX,
        ParserResult, XMLParseBase,
    },
};

pub type Bytes = Vec<u8>;

/// Convert mzML spectrum XML into [`MultiLayerSpectrum`](crate::spectrum::MultiLayerSpectrum)
pub trait SpectrumBuilding<'a, C: CentroidLike, D: DeconvolutedCentroidLike, S: SpectrumLike<C, D>>
{
    /// Get the last isolation window being constructed
    fn isolation_window_mut(&mut self) -> &mut IsolationWindow;
    /// Get the last scan window being constructed.
    fn scan_window_mut(&mut self) -> &mut ScanWindow;

    /// Get the current [`SelectedIon`] being built.
    fn selected_ion_mut(&mut self) -> &mut SelectedIon;

    /// Add a new [`SelectedIon`] to the stack for the current [`Precursor`].
    /// This will be the current [`SelectedIon`].
    fn new_selected_ion(&mut self) -> &mut SelectedIon;

    /// Add a new [`Precursor`] to the stack. This will be the current [`Precursor`].
    fn new_precursor_mut(&mut self) -> &mut Precursor;

    /// Get the current [`Precursor`] being built.
    fn precursor_mut(&mut self) -> &mut Precursor;

    /// Get the current [`DataArray`] being built. This may be an empty instance if
    /// if there is no array being built currently.
    fn current_array_mut(&mut self) -> &mut DataArray;

    /// Move all the data into the provided `spectrum` reference
    fn into_spectrum(self, spectrum: &mut S);

    /// Optionally set the global data processing identifier for the run to be used if
    /// a data processing reference isn't specified locally.
    ///
    /// This is a no-op if not explicitly implemented.
    fn set_run_data_processing(&mut self, _identifier: Option<Box<str>>) {}

    /// Move all the data into the provided `chromatogram` reference
    fn into_chromatogram(self, chromatogram: &mut Chromatogram);

    /// Put a parameter-like instance into the current top-level instance
    fn fill_spectrum<P: ParamLike + Into<Param> + ParamValue>(&mut self, param: P);

    /// Set the compression method for the current [`DataArray`]
    fn set_current_compressiion(&mut self, compression: BinaryCompressionType) {
        trace!(
            "Setting current compression method for {:?} to {compression:?}",
            self.current_array_mut().name()
        );
        self.current_array_mut().compression = compression;
    }

    /// Put a parameter-like instance into the current [`DataArray`]
    fn fill_binary_data_array<P: ParamLike + Into<Param> + ParamValue>(&mut self, param: P) {
        if param.is_ms() {
            match param.accession().unwrap() {
                // Compression types
                x if x == unsafe { BinaryCompressionType::Zlib.accession().unwrap_unchecked() } => {
                    self.set_current_compressiion(BinaryCompressionType::Zlib);
                }
                x if x
                    == unsafe {
                        BinaryCompressionType::NoCompression
                            .accession()
                            .unwrap_unchecked()
                    } =>
                {
                    self.set_current_compressiion(BinaryCompressionType::NoCompression);
                }
                x if x
                    == unsafe {
                        BinaryCompressionType::NumpressLinear
                            .accession()
                            .unwrap_unchecked()
                    } =>
                {
                    self.set_current_compressiion(BinaryCompressionType::NumpressLinear);
                }
                x if x
                    == unsafe {
                        BinaryCompressionType::NumpressPIC
                            .accession()
                            .unwrap_unchecked()
                    } =>
                {
                    self.set_current_compressiion(BinaryCompressionType::NumpressPIC);
                }
                x if x
                    == unsafe {
                        BinaryCompressionType::NumpressSLOF
                            .accession()
                            .unwrap_unchecked()
                    } =>
                {
                    self.set_current_compressiion(BinaryCompressionType::NumpressSLOF);
                }
                x if x
                    == unsafe {
                        BinaryCompressionType::NumpressLinearZlib
                            .accession()
                            .unwrap_unchecked()
                    } =>
                {
                    self.set_current_compressiion(BinaryCompressionType::NumpressLinearZlib);
                }
                x if x
                    == unsafe {
                        BinaryCompressionType::NumpressPICZlib
                            .accession()
                            .unwrap_unchecked()
                    } =>
                {
                    self.set_current_compressiion(BinaryCompressionType::NumpressPICZlib);
                }
                x if x
                    == unsafe {
                        BinaryCompressionType::NumpressSLOFZlib
                            .accession()
                            .unwrap_unchecked()
                    } =>
                {
                    self.set_current_compressiion(BinaryCompressionType::NumpressSLOFZlib);
                }
                x if x
                    == unsafe {
                        BinaryCompressionType::DeltaPrediction
                            .accession()
                            .unwrap_unchecked()
                    } =>
                {
                    self.set_current_compressiion(BinaryCompressionType::DeltaPrediction);
                }
                x if x
                    == unsafe {
                        BinaryCompressionType::LinearPrediction
                            .accession()
                            .unwrap_unchecked()
                    } =>
                {
                    self.set_current_compressiion(BinaryCompressionType::LinearPrediction);
                }
                x if x
                    == unsafe {
                        BinaryCompressionType::ShuffleZstd
                            .accession()
                            .unwrap_unchecked()
                    } =>
                {
                    self.set_current_compressiion(BinaryCompressionType::ShuffleZstd);
                }
                x if x
                    == unsafe {
                        BinaryCompressionType::DeltaShuffleZstd
                            .accession()
                            .unwrap_unchecked()
                    } =>
                {
                    self.set_current_compressiion(BinaryCompressionType::DeltaShuffleZstd);
                }
                x if x == unsafe { BinaryCompressionType::Zstd.accession().unwrap_unchecked() } => {
                    self.set_current_compressiion(BinaryCompressionType::Zstd);
                }
                x if x
                    == unsafe {
                        BinaryCompressionType::ZstdDict
                            .accession()
                            .unwrap_unchecked()
                    } =>
                {
                    self.set_current_compressiion(BinaryCompressionType::ZstdDict);
                }
                x if x
                    == unsafe {
                        BinaryCompressionType::NumpressLinearZstd
                            .accession()
                            .unwrap_unchecked()
                    } =>
                {
                    self.set_current_compressiion(BinaryCompressionType::NumpressLinearZstd);
                }
                x if x
                    == unsafe {
                        BinaryCompressionType::NumpressSLOFZstd
                            .accession()
                            .unwrap_unchecked()
                    } =>
                {
                    self.set_current_compressiion(BinaryCompressionType::NumpressSLOFZstd);
                }
                // Array data types
                1000523 => {
                    self.current_array_mut().dtype = BinaryDataArrayType::Float64;
                }
                1000521 => {
                    self.current_array_mut().dtype = BinaryDataArrayType::Float32;
                }
                1000522 => {
                    self.current_array_mut().dtype = BinaryDataArrayType::Int64;
                }
                1000519 => {
                    self.current_array_mut().dtype = BinaryDataArrayType::Int32;
                }
                1001479 => {
                    self.current_array_mut().dtype = BinaryDataArrayType::ASCII;
                }

                // Array types
                1000514 => {
                    self.current_array_mut().name = ArrayType::MZArray;
                    *self.current_array_mut().unit_mut() = param.unit();
                }
                1000515 => {
                    self.current_array_mut().name = ArrayType::IntensityArray;
                    *self.current_array_mut().unit_mut() = param.unit();
                }
                1000516 => {
                    self.current_array_mut().name = ArrayType::ChargeArray;
                    *self.current_array_mut().unit_mut() = param.unit();
                }
                1000517 => {
                    self.current_array_mut().name = ArrayType::SignalToNoiseArray;
                    *self.current_array_mut().unit_mut() = param.unit();
                }
                1000595 => {
                    self.current_array_mut().name = ArrayType::TimeArray;
                    let unit = param.unit();
                    match unit {
                        Unit::Minute | Unit::Second | Unit::Millisecond => {
                            self.current_array_mut().unit = unit
                        }
                        _ => {
                            warn!("Invalid unit {} found for time array", unit)
                        }
                    }
                }
                1000617 => {
                    self.current_array_mut().name = ArrayType::WavelengthArray;
                    self.current_array_mut().unit = param.unit();
                }
                1000786 => {
                    self.current_array_mut().name = ArrayType::NonStandardDataArray {
                        name: Box::new(param.value().to_string()),
                    };
                    *self.current_array_mut().unit_mut() = param.unit();
                }
                1002477 => {
                    self.current_array_mut().name = ArrayType::MeanDriftTimeArray;
                    self.current_array_mut().unit = param.unit();
                }
                1002816 => {
                    self.current_array_mut().name = ArrayType::MeanIonMobilityArray;
                    self.current_array_mut().unit = param.unit();
                }
                1003006 => {
                    self.current_array_mut().name = ArrayType::MeanInverseReducedIonMobilityArray;
                    self.current_array_mut().unit = param.unit();
                }
                1003007 => {
                    self.current_array_mut().name = ArrayType::RawIonMobilityArray;
                    self.current_array_mut().unit = param.unit();
                }
                1003153 => {
                    self.current_array_mut().name = ArrayType::RawDriftTimeArray;
                    self.current_array_mut().unit = param.unit();
                }
                1003156 => {
                    self.current_array_mut().name = ArrayType::DeconvolutedDriftTimeArray;
                    self.current_array_mut().unit = param.unit();
                }
                1003008 => {
                    self.current_array_mut().name = ArrayType::RawInverseReducedIonMobilityArray;
                    self.current_array_mut().unit = param.unit();
                }
                1003154 => {
                    self.current_array_mut().name = ArrayType::DeconvolutedDriftTimeArray;
                    self.current_array_mut().unit = param.unit();
                }
                1003155 => {
                    self.current_array_mut().name =
                        ArrayType::DeconvolutedInverseReducedIonMobilityArray;
                    self.current_array_mut().unit = param.unit();
                }
                _ => {
                    self.current_array_mut().add_param(param.into());
                }
            }
        } else {
            self.current_array_mut().add_param(param.into());
        }
    }

    /// Put a parameter-like instance into the current [`SelectedIon`]
    fn fill_selected_ion(&mut self, param: Param) {
        match param.name.as_ref() {
            "selected ion m/z" => {
                self.selected_ion_mut().mz = param.to_f64().expect("Failed to parse ion m/z");
            }
            "peak intensity" => {
                self.selected_ion_mut().intensity =
                    param.to_f32().expect("Failed to parse peak intensity");
            }
            "charge state" => {
                self.selected_ion_mut().charge =
                    Some(param.to_i32().expect("Failed to parse ion charge"));
            }
            &_ => {
                self.selected_ion_mut().add_param(param);
            }
        };
    }

    /// Put a parameter-like instance into the current [`IsolationWindow`]
    fn fill_isolation_window(&mut self, param: Param) {
        let window = self.isolation_window_mut();
        match param.name.as_ref() {
            "isolation window target m/z" => {
                window.target = param
                    .to_f32()
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
                let lower_bound = param
                    .to_f32()
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
                let upper_bound = param
                    .to_f32()
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
                let lower_bound = param
                    .to_f32()
                    .expect("Failed to parse isolation window limit");
                if matches!(
                    window.flags,
                    IsolationWindowState::Unknown | IsolationWindowState::Explicit
                ) {
                    window.flags = IsolationWindowState::Explicit;
                    window.lower_bound = lower_bound;
                }
            }
            "isolation window upper limit" => {
                let upper_bound = param
                    .to_f32()
                    .expect("Failed to parse isolation window limit");
                if matches!(
                    window.flags,
                    IsolationWindowState::Unknown | IsolationWindowState::Explicit
                ) {
                    window.flags = IsolationWindowState::Explicit;
                    window.upper_bound = upper_bound;
                }
            }
            &_ => {}
        }
    }

    /// Put a parameter-like instance into the current [`ScanWindow`]
    fn fill_scan_window(&mut self, param: Param) {
        let window = self.scan_window_mut();
        match param.name.as_ref() {
            "scan window lower limit" => {
                window.lower_bound = param.to_f32().expect("Failed to parse scan window limit");
            }
            "scan window upper limit" => {
                window.upper_bound = param.to_f32().expect("Failed to parse scan window limit");
            }
            &_ => {}
        }
    }

    fn borrow_instrument_configuration(
        self,
        instrument_configurations: &'a mut IncrementingIdMap,
    ) -> Self;
}

macro_rules! xml_error {
    ($state:ident, $xml_err:ident) => {
        MzMLParserError::XMLError($state, $xml_err)
    };
    ($state:ident, $xml_err:ident, $ctx:expr) => {
        MzMLParserError::XMLErrorContext($state, $xml_err, $ctx)
    };
}

const BUFFER_SIZE: usize = 10000;

/// An accumulator for the attributes of a spectrum as it is read from an
/// mzML document.
///
/// While this type is public, it is unnecessary for most users. Instead
/// just use [`MzMLReaderType::read_next`].
pub struct MzMLSpectrumBuilder<
    'a,
    C: CentroidLike = CentroidPeak,
    D: DeconvolutedCentroidLike = DeconvolutedPeak,
> {
    pub params: ParamList,
    pub acquisition: Acquisition,
    pub precursor: Vec<Precursor>,

    pub arrays: BinaryArrayMap,
    pub current_array: DataArray,

    pub index: usize,
    pub entry_id: String,
    pub ms_level: u8,
    pub polarity: ScanPolarity,
    pub signal_continuity: SignalContinuity,
    pub has_precursor: bool,
    pub detail_level: DetailLevel,
    pub instrument_id_map: Option<&'a mut IncrementingIdMap>,
    pub run_level_data_processing: Option<Box<str>>,
    pub spectrum_data_processing_ref: Option<Box<str>>,
    entry_type: EntryType,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> Default for MzMLSpectrumBuilder<'_, C, D> {
    fn default() -> Self {
        Self {
            params: Default::default(),
            acquisition: Default::default(),
            precursor: Default::default(),
            arrays: Default::default(),
            current_array: Default::default(),
            index: Default::default(),
            entry_id: Default::default(),
            ms_level: Default::default(),
            polarity: Default::default(),
            signal_continuity: Default::default(),
            has_precursor: Default::default(),
            detail_level: Default::default(),
            instrument_id_map: Default::default(),
            run_level_data_processing: None,
            spectrum_data_processing_ref: None,
            entry_type: Default::default(),
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
        }
    }
}

impl<C: CentroidLike, D: DeconvolutedCentroidLike> XMLParseBase for MzMLSpectrumBuilder<'_, C, D> {}
impl<C: CentroidLike, D: DeconvolutedCentroidLike> CVParamParse for MzMLSpectrumBuilder<'_, C, D> {}

#[doc(hidden)]
impl<'inner, C: CentroidLike, D: DeconvolutedCentroidLike>
    SpectrumBuilding<'inner, C, D, MultiLayerSpectrum<C, D>> for MzMLSpectrumBuilder<'inner, C, D>
{
    fn set_run_data_processing(&mut self, identifier: Option<Box<str>>) {
        self.run_level_data_processing = identifier;
    }

    fn isolation_window_mut(&mut self) -> &mut IsolationWindow {
        &mut self.precursor_mut().isolation_window
    }

    fn scan_window_mut(&mut self) -> &mut ScanWindow {
        let event = self.acquisition.last_scan_mut().unwrap();
        if event.scan_windows.is_empty() {
            event.scan_windows.push(ScanWindow::default());
        }
        event.scan_windows.last_mut().unwrap()
    }

    fn selected_ion_mut(&mut self) -> &mut SelectedIon {
        self.precursor_mut().ion_mut()
    }

    fn current_array_mut(&mut self) -> &mut DataArray {
        &mut self.current_array
    }

    fn into_spectrum(self, spectrum: &mut MultiLayerSpectrum<C, D>) {
        let description = &mut spectrum.description;

        description.id = self.entry_id;
        description.index = self.index;
        description.signal_continuity = self.signal_continuity;
        description.ms_level = self.ms_level;
        description.polarity = self.polarity;

        description.params = self.params;
        description.acquisition = self.acquisition;
        if self.has_precursor {
            description.precursor = self.precursor;
        }

        spectrum.arrays = Some(self.arrays);
    }

    fn fill_spectrum<P: ParamLike + Into<Param> + ParamValue>(&mut self, param: P) {
        match param.name() {
            "ms level" => {
                self.ms_level = param.to_i32().expect("Failed to parse ms level") as u8;
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
                self.params.push(param.into());
            }
        };
    }

    fn borrow_instrument_configuration(
        mut self,
        instrument_configurations: &'inner mut IncrementingIdMap,
    ) -> Self {
        self.instrument_id_map = Some(instrument_configurations);
        self
    }

    fn new_selected_ion(&mut self) -> &mut SelectedIon {
        let prec = self.precursor_mut();
        prec.last_ion_mut()
    }

    fn into_chromatogram(self, chromatogram: &mut Chromatogram) {
        let description = chromatogram.description_mut();
        description.id = self.entry_id;
        description.index = self.index;
        description.polarity = self.polarity;
        description.ms_level = Some(self.ms_level);

        let mut params = Vec::with_capacity(self.params.len());
        for param in self.params.into_iter() {
            if param.is_controlled() {
                if param.is_ms() {
                    if let Some(chrom_type) =
                        ChromatogramType::from_accession(param.accession.unwrap())
                    {
                        description.chromatogram_type = chrom_type;
                    } else {
                        params.push(param)
                    }
                } else {
                    params.push(param)
                }
            } else {
                params.push(param)
            }
        }

        description.params = params;
        if self.has_precursor {
            description.precursor = self.precursor;
        }

        chromatogram.arrays = self.arrays;
    }

    fn new_precursor_mut(&mut self) -> &mut Precursor {
        self.precursor.push(Default::default());
        self.precursor.last_mut().unwrap()
    }

    fn precursor_mut(&mut self) -> &mut Precursor {
        if self.precursor.is_empty() {
            self.new_precursor_mut()
        } else {
            self.precursor.last_mut().unwrap()
        }
    }
}

impl<'inner, C: CentroidLike + BuildFromArrayMap, D: DeconvolutedCentroidLike + BuildFromArrayMap>
    MzMLSpectrumBuilder<'inner, C, D>
{
    pub fn new() -> MzMLSpectrumBuilder<'inner, C, D> {
        Self::with_detail_level(DetailLevel::Full)
    }

    pub fn set_run_data_processing(&mut self, identifier: Option<Box<str>>) {
        self.run_level_data_processing = identifier;
    }

    pub fn with_detail_level(detail_level: DetailLevel) -> MzMLSpectrumBuilder<'inner, C, D> {
        Self {
            detail_level,
            ..Default::default()
        }
    }

    fn warning_context(&self) -> String {
        if self.is_spectrum_entry() {
            format!("spectrum entry {} ({})", self.index, self.entry_id)
        } else if self.is_chromatogram_entry() {
            format!("chromatogram entry {} ({})", self.index, self.entry_id)
        } else {
            format!("unknown entry {} ({})", self.index, self.entry_id)
        }
    }

    pub fn _reset(&mut self) {
        self.params.clear();
        self.acquisition = Acquisition::default();
        self.arrays.clear();
        self.current_array.clear();
        self.entry_id.clear();
        self.entry_type = EntryType::Spectrum;

        self.precursor = Vec::new();
        self.index = 0;
        self.has_precursor = false;
        self.signal_continuity = SignalContinuity::Unknown;
        self.polarity = ScanPolarity::Unknown;
    }

    pub fn set_entry_type(&mut self, entry_type: EntryType) {
        self.entry_type = entry_type;
    }

    pub fn entry_type(&self) -> EntryType {
        self.entry_type
    }

    pub fn is_spectrum_entry(&self) -> bool {
        matches!(self.entry_type, EntryType::Spectrum)
    }

    pub fn is_chromatogram_entry(&self) -> bool {
        matches!(self.entry_type, EntryType::Chromatogram)
    }

    pub fn fill_param_into(&mut self, param: Param, state: MzMLParserState) {
        match state {
            MzMLParserState::Spectrum => {
                self.fill_spectrum(param);
            }
            MzMLParserState::ScanList => {
                if param.is_controlled() {
                    if let Some(comb) = ScanCombination::from_accession(
                        param.controlled_vocabulary.unwrap(),
                        param.accession.unwrap(),
                    ) {
                        self.acquisition.combination = comb
                    } else {
                        self.acquisition.add_param(param)
                    }
                } else {
                    self.acquisition.add_param(param)
                }
            }
            MzMLParserState::Scan => {
                let event = self.acquisition.last_scan_mut().unwrap();
                match param.name.as_bytes() {
                    b"scan start time" => {
                        let value: f64 = param
                            .to_f64()
                            .expect("Expected floating point number for scan time");
                        let value = match &param.unit {
                            Unit::Minute => value,
                            Unit::Second => value / 60.0,
                            Unit::Millisecond => value / 60000.0,
                            _ => {
                                warn!("Could not infer unit for {:?}", param);
                                value
                            }
                        };
                        event.start_time = value;
                    }
                    b"ion injection time" => {
                        event.injection_time = param
                            .to_f64()
                            .expect("Expected floating point number for injection time")
                            as f32;
                    }
                    _ => event.add_param(param),
                }
            }
            MzMLParserState::ScanWindowList => {
                self.acquisition.last_scan_mut().unwrap().add_param(param)
            }
            MzMLParserState::ScanWindow => {
                self.fill_scan_window(param);
            }
            MzMLParserState::IsolationWindow => {
                self.fill_isolation_window(param);
            }
            MzMLParserState::SelectedIon | MzMLParserState::SelectedIonList => {
                self.fill_selected_ion(param);
            }
            MzMLParserState::Activation => {
                if Activation::is_param_activation(&param)
                    && self.precursor_mut().activation.method().is_none()
                {
                    self.precursor_mut()
                        .activation
                        .methods_mut()
                        .push(param.into());
                } else {
                    match param.name.as_ref() {
                        "collision energy" | "activation energy" => {
                            self.precursor_mut().activation.energy =
                                param.to_f32().expect("Failed to parse collision energy");
                        }
                        &_ => {
                            self.precursor_mut().activation.add_param(param);
                        }
                    }
                }
            }
            MzMLParserState::BinaryDataArrayList => {}
            MzMLParserState::BinaryDataArray => {
                self.fill_binary_data_array(param);
            }
            MzMLParserState::Precursor | MzMLParserState::PrecursorList => {
                warn!("cvParam found for {:?} where none are allowed", &state);
            }
            _ => {}
        };
    }
}

impl<C: CentroidLike + BuildFromArrayMap, D: DeconvolutedCentroidLike + BuildFromArrayMap> MzMLSAX
    for MzMLSpectrumBuilder<'_, C, D>
{
    fn start_element(&mut self, event: &BytesStart, state: MzMLParserState) -> ParserResult {
        let elt_name = event.name();
        match elt_name.as_ref() {
            b"spectrum" => {
                self.set_entry_type(EntryType::Spectrum);
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => match attr.key.as_ref() {
                            b"id" => {
                                self.entry_id = match attr.unescape_value() {
                                    Ok(value) => value.to_string(),
                                    Err(e) => {
                                        return Err(xml_error!(
                                            state,
                                            e,
                                            "Failed to decode spectrum id".into()
                                        ));
                                    }
                                };
                                trace!("Stored spectrum id = {}", self.entry_id);
                            }
                            b"index" => {
                                self.index = String::from_utf8_lossy(&attr.value)
                                    .parse::<usize>()
                                    .expect("Failed to parse index");
                                trace!("Stored spectrum index = {}", self.index);
                            }
                            b"dataProcessingRef" => {
                                let ident: Box<str> = String::from_utf8_lossy(&attr.value).into();
                                self.spectrum_data_processing_ref = Some(ident);
                            }
                            _ => {}
                        },
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg.into(), state));
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
                            if attr.key.as_ref() == b"instrumentConfigurationRef" {
                                scan_event.instrument_configuration_id = self
                                    .instrument_id_map
                                    .as_mut()
                                    .expect("An instrument ID map was not provided")
                                    .get(&attr.unescape_value().expect("Error decoding id"));
                            } else if attr.key.as_ref() == b"spectrumRef" {
                                let sref =
                                    attr.unescape_value().expect("Error decoding spectrumRef");
                                scan_event.spectrum_reference = Some(sref.into());
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg.into(), state));
                        }
                    }
                }
                self.acquisition.scans.push(scan_event);
                return Ok(MzMLParserState::Scan);
            }
            b"scanWindow" => {
                let window = ScanWindow::default();
                self.acquisition
                    .last_scan_mut()
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
                self.new_precursor_mut();
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key.as_ref() == b"spectrumRef" {
                                self.precursor_mut().precursor_id = Some(
                                    attr.unescape_value()
                                        .expect("Error decoding id")
                                        .to_string(),
                                );
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg.into(), state));
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
                let mut dp_set = false;
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            if attr.key.as_ref() == b"dataProcessingRef" {
                                match attr.unescape_value() {
                                    Ok(v) => {
                                        self.current_array
                                            .set_data_processing_reference(Some(v.into()));
                                        dp_set = true;
                                        break;
                                    }
                                    Err(msg) => return Err(self.handle_xml_error(msg, state)),
                                }
                            }
                        }
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg.into(), state));
                        }
                    }
                }
                if !dp_set {
                    if let Some(dp_ref) = self.spectrum_data_processing_ref.as_ref() {
                        self.current_array
                            .set_data_processing_reference(Some(dp_ref.clone()));
                    } else if let Some(dp_ref) = self.run_level_data_processing.as_ref() {
                        self.current_array
                            .set_data_processing_reference(Some(dp_ref.clone()));
                    }
                }

                return Ok(MzMLParserState::BinaryDataArray);
            }
            b"binary" => {
                return Ok(MzMLParserState::Binary);
            }
            b"chromatogramList" => return Ok(MzMLParserState::ChromatogramList),
            b"chromatogram" => {
                self.set_entry_type(EntryType::Chromatogram);
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => match attr.key.as_ref() {
                            b"id" => {
                                self.entry_id = attr
                                    .unescape_value()
                                    .expect("Error decoding id")
                                    .to_string();
                                trace!("Stored chromatogram id = {}", self.entry_id);
                            }
                            b"index" => {
                                self.index = String::from_utf8_lossy(&attr.value)
                                    .parse::<usize>()
                                    .expect("Failed to parse index");
                                trace!("Stored chromatogram index = {}", self.index);
                            }
                            _ => {}
                        },
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg.into(), state));
                        }
                    }
                }
                return Ok(MzMLParserState::Chromatogram);
            }
            _ => {}
        };
        Ok(state)
    }

    fn empty_element(
        &mut self,
        event: &BytesStart,
        state: MzMLParserState,
        reader_position: usize,
    ) -> ParserResult {
        let elt_name = event.name();
        match elt_name.as_ref() {
            // Inline the `fill_param_into` to avoid excessive copies.
            b"cvParam" | b"userParam" => {
                match Self::handle_param_borrowed(event, reader_position, state) {
                    Ok(param) => match state {
                        MzMLParserState::Spectrum | MzMLParserState::Chromatogram => {
                            self.fill_spectrum(param)
                        }
                        MzMLParserState::ScanList => {
                            if param.is_controlled() {
                                if let Some(comb) = ScanCombination::from_accession(
                                    param.controlled_vocabulary.unwrap(),
                                    param.accession.unwrap(),
                                ) {
                                    self.acquisition.combination = comb
                                } else {
                                    self.acquisition.add_param(param.into())
                                }
                            } else {
                                self.acquisition.add_param(param.into())
                            }
                        }
                        MzMLParserState::Scan => {
                            match param.name.as_bytes() {
                                b"scan start time" => {
                                    let value: f64 = param
                                    .to_f64()
                                    .unwrap_or_else(|e| panic!("Expected floating point number for scan time: {e} for {}", self.warning_context()));
                                    let value = match &param.unit {
                                        Unit::Minute => value,
                                        Unit::Second => value / 60.0,
                                        Unit::Millisecond => value / 60000.0,
                                        _ => {
                                            warn!(
                                                "Could not infer unit for {:?} for {}",
                                                param,
                                                self.warning_context()
                                            );
                                            value
                                        }
                                    };
                                    self.acquisition.last_scan_mut().unwrap().start_time = value;
                                }
                                b"ion injection time" => {
                                    self.acquisition.last_scan_mut().unwrap().injection_time = param.to_f32().unwrap_or_else(
                                            |e| panic!("Expected floating point number for injection time: {e} for {}", self.warning_context())
                                        );
                                }
                                _ => self
                                    .acquisition
                                    .last_scan_mut()
                                    .unwrap()
                                    .add_param(param.into()),
                            }
                        }
                        MzMLParserState::ScanWindowList => self
                            .acquisition
                            .last_scan_mut()
                            .unwrap()
                            .add_param(param.into()),
                        MzMLParserState::ScanWindow => {
                            self.fill_scan_window(param.into());
                        }
                        MzMLParserState::IsolationWindow => {
                            self.fill_isolation_window(param.into());
                        }
                        MzMLParserState::SelectedIon | MzMLParserState::SelectedIonList => {
                            self.fill_selected_ion(param.into());
                        }
                        MzMLParserState::Activation => {
                            if Activation::is_param_activation(&param) {
                                self.precursor_mut()
                                    .activation
                                    .methods_mut()
                                    .push(param.into());
                            } else {
                                let dissociation_energy = param.curie().and_then(|c| {
                                        DissociationEnergyTerm::from_curie(&c, param.value().to_f32().unwrap_or_else(|e| {
                                            warn!("Failed to convert dissociation energy: {e} for {} for {}", param.name(), self.warning_context());
                                            0.0
                                        }))
                                    });
                                match dissociation_energy {
                                    Some(t) => {
                                        if t.is_supplemental() {
                                            self.precursor_mut().activation.add_param(param.into())
                                        } else {
                                            if self.precursor_mut().activation.energy != 0.0 {
                                                warn!(
                                                    "Multiple dissociation energies detected. Saw {t} after already setting dissociation energy for {}",
                                                    self.warning_context()
                                                );
                                            }
                                            self.precursor_mut().activation.energy = t.energy();
                                        }
                                    }
                                    None => {
                                        self.precursor_mut().activation.add_param(param.into());
                                    }
                                }
                            }
                        }
                        MzMLParserState::BinaryDataArrayList => {}
                        MzMLParserState::BinaryDataArray => {
                            self.fill_binary_data_array(param);
                        }
                        MzMLParserState::Precursor | MzMLParserState::PrecursorList => {
                            warn!("cvParam found for {:?} where none are allowed", &state);
                        }
                        _ => {}
                    },
                    Err(err) => return Err(err),
                }
            }
            &_ => {}
        }
        Ok(state)
    }

    fn end_element(&mut self, event: &BytesEnd, state: MzMLParserState) -> ParserResult {
        let elt_name = event.name();
        match elt_name.as_ref() {
            b"spectrum" => return Ok(MzMLParserState::SpectrumDone),
            b"chromatogram" => return Ok(MzMLParserState::ChromatogramDone),
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
                let mut array = mem::take(&mut self.current_array);
                if self.detail_level == DetailLevel::Full {
                    array.decode_and_store().map_err(|e| {
                        let new_err =
                            MzMLParserError::ArrayDecodingError(state, array.name.clone(), e);
                        log::error!("Failed to decode mzML array: {new_err}");
                        new_err
                    })?;
                }
                self.arrays.add(array);
                return Ok(MzMLParserState::BinaryDataArrayList);
            }
            b"binary" => return Ok(MzMLParserState::BinaryDataArray),
            b"spectrumList" => return Ok(MzMLParserState::SpectrumListDone),
            b"chromatogramList" => return Ok(MzMLParserState::ChromatogramListDone),
            _ => {}
        };
        Ok(state)
    }

    fn text(&mut self, event: &BytesText, state: MzMLParserState) -> ParserResult {
        if state == MzMLParserState::Binary && self.detail_level != DetailLevel::MetadataOnly {
            let bin = event
                .unescape()
                .map_err(|e| MzMLParserError::XMLError(state, e))?;
            self.current_array.data = Bytes::from(bin.as_bytes());
        }
        Ok(state)
    }
}

#[doc(hidden)]
impl<'a, C: CentroidLike, D: DeconvolutedCentroidLike> From<MzMLSpectrumBuilder<'a, C, D>>
    for CentroidSpectrumType<C>
where
    C: BuildFromArrayMap + BuildArrayMapFrom,
    D: BuildFromArrayMap + BuildArrayMapFrom,
{
    fn from(val: MzMLSpectrumBuilder<'a, C, D>) -> Self {
        let mut spec = MultiLayerSpectrum::<C, D>::default();
        val.into_spectrum(&mut spec);
        spec.try_into().unwrap()
    }
}

#[doc(hidden)]
impl<'a, C: CentroidLike, D: DeconvolutedCentroidLike> From<MzMLSpectrumBuilder<'a, C, D>>
    for MultiLayerSpectrum<C, D>
{
    fn from(val: MzMLSpectrumBuilder<'a, C, D>) -> Self {
        let mut spec = MultiLayerSpectrum::<C, D>::default();
        val.into_spectrum(&mut spec);
        spec
    }
}

#[doc(hidden)]
impl<'a> From<MzMLSpectrumBuilder<'a>> for RawSpectrum {
    fn from(val: MzMLSpectrumBuilder<'a>) -> Self {
        let mut spec = Spectrum::default();
        val.into_spectrum(&mut spec);
        spec.into()
    }
}

/**
An mzML parser that supports iteration and random access. The parser produces
[`Spectrum`] instances, which may be converted to [`RawSpectrum`](crate::spectrum::RawSpectrum)
or [`CentroidSpectrum`](crate::spectrum::CentroidSpectrum) as is appropriate to the data.

When the readable stream the parser is wrapped around supports [`io::Seek`],
additional random access operations are available.
*/
pub struct MzMLReaderType<
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
}

impl<
    R: Read + Seek,
    C: CentroidLike + BuildFromArrayMap,
    D: DeconvolutedCentroidLike + BuildFromArrayMap,
> IntoIonMobilityFrameSource<C, D> for MzMLReaderType<R, C, D>
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
> MzMLReaderType<R, C, D>
{
    /// Create a new [`MzMLReaderType`] instance, wrapping the [`io::Read`] handle
    /// provided with an [`io::BufReader`] and parses the metadata section of the file.
    pub fn new(file: R) -> MzMLReaderType<R, C, D> {
        Self::with_buffer_capacity_and_detail_level(file, BUFFER_SIZE, DetailLevel::Full)
    }

    pub fn with_buffer_capacity_and_detail_level(
        file: R,
        capacity: usize,
        detail_level: DetailLevel,
    ) -> MzMLReaderType<R, C, D> {
        let handle = BufReader::with_capacity(capacity, file);
        let mut inst = MzMLReaderType {
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
        let mut accumulator = FileMetadataBuilder {
            instrument_id_map: Some(&mut self.instrument_id_map),
            ..Default::default()
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
        self.file_description = accumulator.file_description;
        self.instrument_configurations = accumulator
            .instrument_configurations
            .into_iter()
            .map(|ic| (ic.id, ic))
            .collect();
        self.softwares = accumulator.softwares;
        self.samples = accumulator.samples;
        self.data_processings = accumulator.data_processings;
        self.reference_param_groups = accumulator.reference_param_groups;
        self.scan_settings = accumulator.scan_settings;
        self.run.id = accumulator.run_id;
        self.run.default_instrument_id = accumulator.default_instrument_config;
        self.run.default_source_file_id = accumulator.default_source_file;
        self.run.start_time = accumulator.start_timestamp;
        self.run.default_data_processing_id = accumulator.default_data_processing;
        self.num_spectra = accumulator.num_spectra;

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
        accumulator.set_run_data_processing(
            self.run
                .default_data_processing_id
                .clone()
                .map(|v| v.into_boxed_str()),
        );
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
            MzMLParserState::EOF => Err(MzMLParserError::EOF),
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
                    MzMLParserError::EOF => {}
                    err => log::error!("Error while reading mzML spectrum: {err}"),
                };
                Err(err)
            }
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
                    MzMLParserError::EOF => {}
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
            }
        }
    }
}

/// When the underlying stream supports random access, this type can read the index at the end of
/// an `indexedmzML` document and use the offset map to jump to immediately jump to a specific spectrum
impl<
    R: SeekRead,
    C: CentroidLike + BuildFromArrayMap,
    D: DeconvolutedCentroidLike + BuildFromArrayMap,
> MzMLReaderType<R, C, D>
{
    pub fn check_stream(&mut self, next_tag: &str) -> Result<bool, MzMLParserError> {
        let position = match self.stream_position() {
            Ok(pos) => pos,
            Err(err) => return Err(MzMLParserError::IOError(self.state, err)),
        };
        let mut reader = Reader::from_reader(&mut self.handle);
        reader.trim_text(true);
        let matched_tag = match reader.read_event_into(&mut self.buffer) {
            Ok(event) => match event {
                Event::Start(ref e) => {
                    trace!(
                        "From {}, the next element started was {}",
                        position,
                        String::from_utf8_lossy(e.name().0)
                    );
                    e.name().0 == next_tag.as_bytes()
                }
                Event::End(ref e) => {
                    trace!(
                        "From {}, the next element ended was {}",
                        position,
                        String::from_utf8_lossy(e.name().0)
                    );
                    false
                }
                Event::Empty(ref e) => {
                    trace!(
                        "From {}, the next empty element was {}",
                        position,
                        String::from_utf8_lossy(e.name().0)
                    );
                    e.name().0 == next_tag.as_bytes()
                }
                Event::Text(ref e) => {
                    trace!(
                        "From {}, the next was a text node of {} bytes",
                        position,
                        e.len()
                    );
                    false
                }
                Event::Eof => {
                    trace!("From {}, the next was EOF", position);
                    false
                }
                e => {
                    trace!("From {}, the next was {:?}", position, e);
                    false
                }
            },
            Err(err) => {
                self.buffer.clear();
                return Err(MzMLParserError::XMLError(self.state, err));
            }
        };
        self.buffer.clear();
        match self.seek(SeekFrom::Start(position)) {
            Ok(_) => {}
            Err(err) => return Err(MzMLParserError::IOError(self.state, err)),
        }
        Ok(matched_tag)
    }

    pub fn get_chromatogram_by_id(&mut self, id: &str) -> Option<Chromatogram> {
        let offset = self.chromatogram_index.get(id)?;
        let start = self
            .handle
            .stream_position()
            .expect("Failed to save checkpoint");
        self.seek(SeekFrom::Start(offset))
            .expect("Failed to move seek to offset");
        debug_assert!(
            self.check_stream("chromatogram").unwrap(),
            "The next XML tag was not `chromatogram`"
        );
        self.state = MzMLParserState::Resume;
        let result = self._read_next_chromatogram();
        self.seek(SeekFrom::Start(start))
            .expect("Failed to restore offset");
        result.ok()
    }

    pub fn get_chromatogram_by_index(&mut self, index: usize) -> Option<Chromatogram> {
        let (_key, offset) = self.chromatogram_index.get_index(index)?;
        let start = self
            .handle
            .stream_position()
            .expect("Failed to save checkpoint");
        self.seek(SeekFrom::Start(offset))
            .expect("Failed to move seek to offset");
        debug_assert!(
            self.check_stream("chromatogram").unwrap(),
            "The next XML tag was not `chromatogram`"
        );
        self.state = MzMLParserState::Resume;
        let result = self._read_next_chromatogram();
        self.seek(SeekFrom::Start(start))
            .expect("Failed to restore offset");
        result.ok()
    }

    pub fn iter_chromatograms(&mut self) -> ChromatogramIter<'_, R, C, D> {
        ChromatogramIter::new(self)
    }
}

impl<
    R: SeekRead,
    C: CentroidLike + BuildFromArrayMap,
    D: DeconvolutedCentroidLike + BuildFromArrayMap,
> ChromatogramSource for MzMLReaderType<R, C, D>
{
    fn get_chromatogram_by_id(&mut self, id: &str) -> Option<Chromatogram> {
        self.get_chromatogram_by_id(id)
    }

    fn get_chromatogram_by_index(&mut self, index: usize) -> Option<Chromatogram> {
        self.get_chromatogram_by_index(index)
    }
}

/// [`MzMLReaderType`] instances are [`Iterator`]s over [`Spectrum`]
impl<
    R: io::Read,
    C: CentroidLike + BuildFromArrayMap,
    D: DeconvolutedCentroidLike + BuildFromArrayMap,
> Iterator for MzMLReaderType<R, C, D>
{
    type Item = MultiLayerSpectrum<C, D>;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_next()
    }
}

/// They can also be used to fetch specific spectra by ID, index, or start
/// time when the underlying file stream supports [`io::Seek`].
impl<
    R: SeekRead,
    C: CentroidLike + BuildFromArrayMap,
    D: DeconvolutedCentroidLike + BuildFromArrayMap,
> SpectrumSource<C, D, MultiLayerSpectrum<C, D>> for MzMLReaderType<R, C, D>
{
    /// Retrieve a spectrum by it's native ID
    fn get_spectrum_by_id(&mut self, id: &str) -> Option<MultiLayerSpectrum<C, D>> {
        let offset = self.spectrum_index.get(id)?;
        let start = self
            .handle
            .stream_position()
            .expect("Failed to save checkpoint");
        self.seek(SeekFrom::Start(offset))
            .expect("Failed to move seek to offset");
        debug_assert!(
            self.check_stream("spectrum").unwrap(),
            "The next XML tag was not `spectrum`"
        );
        self.state = MzMLParserState::Resume;
        let result = self.read_next();
        self.seek(SeekFrom::Start(start))
            .expect("Failed to restore offset");
        result
    }

    /// Retrieve a spectrum by it's integer index
    fn get_spectrum_by_index(&mut self, index: usize) -> Option<MultiLayerSpectrum<C, D>> {
        let (_id, offset) = self.spectrum_index.get_index(index)?;
        let byte_offset = offset;
        let start = self
            .handle
            .stream_position()
            .expect("Failed to save checkpoint");
        self.seek(SeekFrom::Start(byte_offset)).ok()?;
        debug_assert!(
            self.check_stream("spectrum").unwrap(),
            "The next XML tag was not `spectrum`"
        );
        self.state = MzMLParserState::Resume;
        let result = self.read_next();
        self.seek(SeekFrom::Start(start))
            .expect("Failed to restore offset");
        result
    }

    /// Return the data stream to the beginning
    fn reset(&mut self) {
        self.state = MzMLParserState::Resume;
        self.seek(SeekFrom::Start(0))
            .expect("Failed to reset file stream");
    }

    fn get_index(&self) -> &OffsetIndex {
        if !self.spectrum_index.init {
            warn!("Attempting to use an uninitialized offset index on MzMLReaderType")
        }
        &self.spectrum_index
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.spectrum_index = index
    }

    fn detail_level(&self) -> &DetailLevel {
        &self.detail_level
    }

    fn set_detail_level(&mut self, detail_level: DetailLevel) {
        self.detail_level = detail_level;
    }
}

/// The iterator can also be updated to move to a different location in the
/// stream efficiently.
impl<
    R: SeekRead,
    C: CentroidLike + BuildFromArrayMap,
    D: DeconvolutedCentroidLike + BuildFromArrayMap,
> RandomAccessSpectrumIterator<C, D, MultiLayerSpectrum<C, D>> for MzMLReaderType<R, C, D>
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, SpectrumAccessError> {
        match self._offset_of_id(id) {
            Some(offset) => match self.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(SpectrumAccessError::IOError(Some(err))),
            },
            None => Err(SpectrumAccessError::SpectrumIdNotFound(id.to_string())),
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, SpectrumAccessError> {
        match self._offset_of_index(index) {
            Some(offset) => match self.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(SpectrumAccessError::IOError(Some(err))),
            },
            None => Err(SpectrumAccessError::SpectrumIndexNotFound(index)),
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, SpectrumAccessError> {
        match self._offset_of_time(time) {
            Some(offset) => match self.seek(SeekFrom::Start(offset)) {
                Ok(_) => Ok(self),
                Err(err) => Err(SpectrumAccessError::IOError(Some(err))),
            },
            None => Err(SpectrumAccessError::SpectrumNotFound),
        }
    }
}

enum IndexRecoveryOperation {
    EOLMismatchSuspected,
    IOFailure(io::Error),
}

impl<
    R: SeekRead,
    C: CentroidLike + BuildFromArrayMap,
    D: DeconvolutedCentroidLike + BuildFromArrayMap,
> MzMLReaderType<R, C, D>
{
    /// Construct a new MzMLReaderType and build an offset index
    /// using [`Self::build_index`]
    pub fn new_indexed(file: R) -> MzMLReaderType<R, C, D> {
        Self::with_buffer_capacity_and_detail_level_indexed(file, BUFFER_SIZE, DetailLevel::Full)
    }

    fn _read_index(&mut self) {
        if let Err(err) = self.read_index_from_end() {
            debug!("Failed to read index from the end of the file: {}", err);
            match self.seek(SeekFrom::Start(0)) {
                Ok(_) => {
                    self.build_index();
                    return;
                }
                Err(error) => {
                    panic!(
                        "Unrecoverable IO Error during file pointer reset {} while handling {:?}",
                        error, err
                    );
                }
            }
        }
        match self.verify_index() {
            Ok(()) => {
                // No extra work required, the index checks out.
            }
            Err(e) => match e {
                IndexRecoveryOperation::EOLMismatchSuspected => {
                    warn!("Rebuilding index, EOL mismatch suspected");
                    self.seek(SeekFrom::Start(0)).unwrap_or_else(|e| {
                        panic!("An IO error occurred while trying to recover the index: {e}")
                    });
                    self.build_index();
                }
                IndexRecoveryOperation::IOFailure(err) => {
                    panic!("An IO error occurred while validating the index: {err}")
                }
            },
        }
    }

    fn verify_index(&mut self) -> Result<(), IndexRecoveryOperation> {
        let n = self.spectrum_index.len();
        trace!("Verifying offset index of length {n}");
        let position = self
            .handle
            .stream_position()
            .map_err(IndexRecoveryOperation::IOFailure)?;
        if n > 0 {
            // Try to pick a spectrum that's close to the beginning of the file to avoid large
            // amounts of wasted scanning for non-linear files, but pick one far enough in it would
            // be affected by byte drift.
            let center = (n / 2).min(100);
            trace!("Checking offset of {center}th spectrum");
            let dl = self.detail_level;
            self.set_detail_level(DetailLevel::MetadataOnly);
            let s = self.get_spectrum_by_index(center);
            self.set_detail_level(dl);
            let s_found = s.is_some_and(|s| s.index() == center);
            if s_found {
                self.seek(SeekFrom::Start(position))
                    .map_err(IndexRecoveryOperation::IOFailure)?;
                return Ok(());
            } else {
                match self.handle.fill_buf() {
                    Ok(buf) => {
                        if let Some(i) = buf.iter().position(|b| *b == b'\r') {
                            if let Some(b2) = buf.get(i + 1) {
                                let has_windows_eol = *b2 == b'\n';
                                if has_windows_eol {
                                    warn!(
                                        "Carriage return line endings detected and offset index is not valid"
                                    );
                                    self.seek(SeekFrom::Start(position))
                                        .map_err(IndexRecoveryOperation::IOFailure)?;
                                    return Err(IndexRecoveryOperation::EOLMismatchSuspected);
                                }
                            }
                        }
                    }
                    Err(e) => return Err(IndexRecoveryOperation::IOFailure(e)),
                }
            }
        }
        self.seek(SeekFrom::Start(position))
            .map_err(IndexRecoveryOperation::IOFailure)?;
        Ok(())
    }

    pub fn with_buffer_capacity_and_detail_level_indexed(
        file: R,
        capacity: usize,
        detail_level: DetailLevel,
    ) -> MzMLReaderType<R, C, D> {
        let mut reader = Self::with_buffer_capacity_and_detail_level(file, capacity, detail_level);
        reader._read_index();
        reader
    }

    pub fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.handle.seek(pos)
    }

    pub fn stream_position(&mut self) -> io::Result<u64> {
        self.handle.stream_position()
    }

    /// Read the checksum from the end of an `indexedmzML` document
    pub fn read_checksum(&mut self) -> io::Result<Option<String>> {
        let current_position = self.handle.stream_position()?;

        self.handle.seek(SeekFrom::End(-200))?;
        let mut buf = Bytes::new();
        self.handle.read_to_end(&mut buf)?;

        let pattern = regex::Regex::new("<fileChecksum>([0-9a-zA-Z]+)</fileChecksum>").unwrap();
        if let Some(captures) = pattern.captures(&String::from_utf8_lossy(&buf)) {
            if let Some(hit) = captures.get(1) {
                self.handle.seek(SeekFrom::Start(current_position))?;
                return Ok(Some(hit.as_str().to_string()));
            }
        }

        self.handle.seek(SeekFrom::Start(current_position))?;
        Ok(None)
    }

    /// Read the offset index at the end of an `<indexedmzML>` document,
    /// though this index may be malformed in some older files.
    pub fn read_index_from_end(&mut self) -> Result<u64, MzMLIndexingError> {
        let mut indexer = IndexedMzMLIndexExtractor::new();
        let current_position = match self.handle.stream_position() {
            Ok(position) => position,
            Err(err) => return Err(MzMLIndexingError::IOError(err)),
        };
        let offset = match indexer.find_offset_from_reader(&mut self.handle) {
            Ok(offset) => {
                if let Some(offset) = offset {
                    offset
                } else {
                    return Err(MzMLIndexingError::OffsetNotFound);
                }
            }
            Err(err) => return Err(err),
        };
        let mut indexer_state = IndexParserState::Start;
        self.handle
            .seek(SeekFrom::Start(offset))
            .expect("Failed to seek to the index offset");

        let mut reader = Reader::from_reader(&mut self.handle);
        reader.trim_text(true);

        loop {
            match reader.read_event_into(&mut self.buffer) {
                Ok(Event::Start(ref e)) => {
                    match indexer.start_element(e, indexer_state) {
                        Ok(state) => {
                            indexer_state = state;
                            if matches!(indexer_state, IndexParserState::Done) {
                                break;
                            }
                        }
                        Err(err) => return Err(err),
                    };
                }
                Ok(Event::End(ref e)) => {
                    match indexer.end_element(e, indexer_state) {
                        Ok(state) => {
                            indexer_state = state;
                        }
                        Err(err) => return Err(err),
                    };
                }
                Ok(Event::Text(ref e)) => {
                    match indexer.text(e, indexer_state) {
                        Ok(state) => {
                            indexer_state = state;
                        }
                        Err(message) => return Err(MzMLIndexingError::XMLError(message)),
                    };
                }
                Ok(Event::Eof) => {
                    break;
                }
                Err(err) => return Err(MzMLIndexingError::XMLError(err)),
                _ => {}
            }
        }
        self.buffer.clear();
        self.spectrum_index = indexer.spectrum_index;
        self.spectrum_index.init = true;
        *self.chromatogram_index = indexer.chromatogram_index;
        self.chromatogram_index.init = true;
        self.handle.seek(SeekFrom::Start(current_position))?;
        Ok(self.spectrum_index.len() as u64)
    }

    /// Builds an offset index to each `<spectrum>` XML element
    /// by doing a fast pre-scan of the XML file.
    pub fn build_index(&mut self) -> u64 {
        let start = self
            .handle
            .stream_position()
            .expect("Failed to save restore location");
        trace!(
            "Starting to build offset index by traversing the file, storing last position as {start}"
        );
        self.seek(SeekFrom::Start(0))
            .expect("Failed to reset stream to beginning");
        let mut reader = Reader::from_reader(&mut self.handle);
        reader.trim_text(true);
        loop {
            match reader.read_event_into(&mut self.buffer) {
                Ok(Event::Start(ref e)) => {
                    let element_name = e.name();
                    if element_name.as_ref() == b"spectrum" {
                        // Hit a spectrum, extract ID and save current offset

                        for attr_parsed in e.attributes() {
                            match attr_parsed {
                                Ok(attr) => {
                                    if attr.key.as_ref() == b"id" {
                                        let scan_id = attr
                                            .unescape_value()
                                            .expect("Error decoding spectrum id in streaming mzML index")
                                            .to_string();
                                        // This count is off by 2 because somehow the < and > bytes are removed?
                                        self.spectrum_index.insert(
                                            scan_id,
                                            (reader.buffer_position() - e.len() - 2) as u64,
                                        );
                                        break;
                                    };
                                }
                                Err(_msg) => {}
                            }
                        }
                    }
                }
                Ok(Event::End(ref e)) => {
                    let element_name = e.name();
                    if element_name.as_ref() == b"spectrumList" {
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
        trace!("Ended indexing scan at offset {offset}. Restoring starting position {start}");
        self.handle
            .seek(SeekFrom::Start(start))
            .expect("Failed to restore location");
        self.spectrum_index.init = true;
        if self.spectrum_index.is_empty() {
            warn!("An index was built but no entries were found")
        }
        offset
    }
}

impl<C: CentroidLike + BuildFromArrayMap, D: DeconvolutedCentroidLike + BuildFromArrayMap>
    MZFileReader<C, D, MultiLayerSpectrum<C, D>> for MzMLReaderType<fs::File, C, D>
{
    fn open_file(source: fs::File) -> io::Result<Self> {
        Ok(Self::new_indexed(source))
    }

    fn construct_index_from_stream(&mut self) -> u64 {
        trace!("Constructing index from stream");
        if let Ok(count) = self.read_index_from_end() {
            count
        } else {
            self.seek(SeekFrom::Start(0)).unwrap();
            if self.spectrum_index.is_empty() && !self.spectrum_index.init {
                self.build_index()
            } else {
                trace!("Index already constructed, skipping full scan");
                self.spectrum_index.len() as u64
            }
        }
    }
}

impl<R: Read, C: CentroidLike, D: DeconvolutedCentroidLike> MSDataFileMetadata
    for MzMLReaderType<R, C, D>
{
    crate::impl_metadata_trait!();

    fn scan_settings(&self) -> Option<&Vec<ScanSettings>> {
        Some(&self.scan_settings)
    }

    fn scan_settings_mut(&mut self) -> Option<&mut Vec<ScanSettings>> {
        Some(&mut self.scan_settings)
    }

    fn spectrum_count_hint(&self) -> Option<u64> {
        self.num_spectra
    }

    fn set_spectrum_count_hint(&mut self, _value: Option<u64>) {
        self.num_spectra = _value;
    }

    fn run_description(&self) -> Option<&MassSpectrometryRun> {
        Some(&self.run)
    }

    fn run_description_mut(&mut self) -> Option<&mut MassSpectrometryRun> {
        Some(&mut self.run)
    }
}

/// A specialization of [`MzMLReaderType`] for the default peak types, for common use.
pub type MzMLReader<R> = MzMLReaderType<R, CentroidPeak, DeconvolutedPeak>;

pub(crate) fn is_mzml(buf: &[u8]) -> bool {
    let mut bufread = BufReader::new(io::Cursor::new(buf));
    let mut reader = Reader::from_reader(&mut bufread);
    let mut buffer = Vec::new();
    loop {
        match reader.read_event_into(&mut buffer) {
            Ok(Event::Start(ref e)) => {
                let elt_name = e.name();
                match elt_name.as_ref() {
                    b"mzML" => return true,
                    b"indexedmzML" => return true,
                    _ => {}
                }
            }
            Ok(Event::Eof) => return false,
            Ok(_) => {}
            Err(_) => return false,
        }
    }
}

pub struct ChromatogramIter<
    'a,
    R: SeekRead,
    C: CentroidLike + BuildFromArrayMap,
    D: DeconvolutedCentroidLike + BuildFromArrayMap,
> {
    reader: &'a mut MzMLReaderType<R, C, D>,
    index: usize,
}

impl<
    'a,
    R: SeekRead,
    C: CentroidLike + BuildFromArrayMap,
    D: DeconvolutedCentroidLike + BuildFromArrayMap,
> ChromatogramIter<'a, R, C, D>
{
    pub fn new(reader: &'a mut MzMLReaderType<R, C, D>) -> Self {
        Self { reader, index: 0 }
    }
}

impl<
    R: SeekRead,
    C: CentroidLike + BuildFromArrayMap,
    D: DeconvolutedCentroidLike + BuildFromArrayMap,
> Iterator for ChromatogramIter<'_, R, C, D>
{
    type Item = Chromatogram;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.reader.chromatogram_index.len() {
            let result = self.reader.get_chromatogram_by_index(self.index);
            self.index += 1;
            result
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::io::traits::SpectrumGrouping;
    use crate::params::ControlledVocabulary;
    use crate::spectrum::spectrum_types::SpectrumLike;
    use std::fs;
    use std::path;

    fn test_metadata<T: MSDataFileMetadata>(reader: &T) {
        assert_eq!(reader.data_processings().len(), 1);
        assert_eq!(reader.instrument_configurations().len(), 2);
        assert_eq!(reader.softwares().len(), 2);
        assert_eq!(reader.file_description().source_files.len(), 1);
        assert_eq!(reader.file_description().contents.len(), 2);

        assert!(
            reader
                .file_description()
                .get_param_by_accession("MS:1000579")
                .is_some()
        );
        assert_eq!(reader.file_description().source_files[0].name, "small.RAW");

        assert!(reader.file_description().has_ms1_spectra());
        assert!(reader.file_description().has_msn_spectra());
        assert!(reader.file_description().has_contents());
        assert!(
            reader
                .file_description()
                .source_files
                .first()
                .unwrap()
                .native_id_format()
                .is_some()
        );

        let config = reader.instrument_configurations().get(&0).unwrap();
        let comp = config.iter().find_map(|c| c.mass_analyzer()).unwrap();
        assert_eq!(
            comp.name(),
            "fourier transform ion cyclotron resonance mass spectrometer"
        );
        assert_eq!(
            config.components.get(1).unwrap().name(),
            Some("fourier transform ion cyclotron resonance mass spectrometer")
        );

        let comp = config.iter().find_map(|c| c.detector()).unwrap();
        assert_eq!(comp.name(), "inductive detector");
        assert_eq!(
            config.components.get(2).unwrap().name(),
            Some("inductive detector")
        );

        let comp = config.iter().find_map(|c| c.ionization_type()).unwrap();
        assert_eq!(comp.name(), "electrospray ionization");
        assert_eq!(
            config.components.first().unwrap().name(),
            Some("electrospray ionization")
        );

        for comp in config.iter() {
            assert!(!comp.parent_types().is_empty());
        }

        assert_eq!(config.len(), 3);
        assert!(!config.is_empty());
        assert_eq!(config.last().unwrap().order, 3);

        assert_eq!(reader.samples().first().unwrap().number(), None);
    }

    #[test_log::test]
    fn reader_from_file() {
        let path = path::Path::new("./test/data/small.mzML");
        let file = fs::File::open(path).expect("Test file doesn't exist");
        let mut reader = MzMLReaderType::<_, CentroidPeak, DeconvolutedPeak>::new(file);
        let mut ms1_count = 0;
        let mut msn_count = 0;

        test_metadata(&reader);
        test_metadata(&reader.iter());

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

    #[test_log::test]
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

        test_metadata(&reader);
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
            let filter_string = scan
                .acquisition()
                .first_scan()
                .unwrap()
                .get_param_by_accession("MS:1000512")
                .unwrap();
            let configs = scan.acquisition().instrument_configuration_ids();
            let conf = configs[0];
            if filter_string.value.to_string().contains("ITMS") {
                assert_eq!(conf, 1);
            } else {
                assert_eq!(conf, 0);
            }
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

        test_metadata(&reader.groups());
    }

    #[test_log::test]
    fn find_offset() -> io::Result<()> {
        let path = path::Path::new("./test/data/small.mzML");
        let mut f = fs::File::open(path)?;

        let mut index = IndexedMzMLIndexExtractor::new();
        if let Ok(offset) = index.find_offset_from_reader(&mut f) {
            if let Some(offset) = offset {
                assert_eq!(offset, 5116653);
            } else {
                panic!("Failed to parse offset from element")
            }
        } else {
            panic!("Failed to find offset element")
        }

        Ok(())
    }

    #[test_log::test]
    fn read_index() -> io::Result<()> {
        let path = path::Path::new("./test/data/read_index_of.mzML");
        let file = fs::File::open(path)?;
        let mut reader = MzMLReader::new(file);
        match reader.read_index_from_end() {
            Ok(_count) => {}
            Err(err) => {
                panic!("Failed to parse out index {:?}", err);
            }
        };
        assert!(!reader.spectrum_index.is_empty());
        let mut reader2 = MzMLReader::new(fs::File::open(path)?);
        reader2.build_index();
        for (k, v) in reader2.spectrum_index.iter() {
            let v2 = reader.spectrum_index.get(k);
            assert_eq!(v2, Some(*v));
        }
        let checksum = reader2.read_checksum()?;
        assert_eq!(
            checksum,
            Some("148ffca890b2bc1701be942a91d7d8aad56c9557".to_string())
        );
        Ok(())
    }

    #[test_log::test]
    fn read_index_raw() -> io::Result<()> {
        let path = path::Path::new("./test/data/read_index_of.mzML");
        let file = fs::File::open(path)?;
        let mut reader = io::BufReader::new(file);

        let mut indexer = IndexedMzMLIndexExtractor::new();
        let offset = indexer
            .find_offset_from_reader(&mut reader)?
            .expect("Failed to find index offset");

        assert_eq!(3548711, offset);

        let stream_pos = reader.seek(SeekFrom::Start(offset))?;

        assert_eq!(offset, stream_pos);

        let mut buf = [0u8; 500];
        reader.read_exact(&mut buf)?;

        let decoded = String::from_utf8_lossy(&buf);
        let needle = r#"<indexList count="2">
<index name="spectrum">
<offset idRef="controllerType=0 controllerNumber=1 scan=1">5003</offset>
<offset idRef="controllerType=0 controllerNumber=1 scan=2">260791</offset>
<offset idRef="controllerType=0 controllerNumber=1 scan=3">451106</offset>
<offset idRef="controllerType=0 controllerNumber=1 scan=4">461270</offset>
<offset idRef="controllerType=0 controllerNumber=1 scan=5">476248</offset>
"#;
        for line in needle.split('\n') {
            assert!(
                decoded.contains(line),
                "Failed to find {} in {}",
                line,
                decoded
            );
        }

        Ok(())
    }

    #[test]
    fn test_random_access_iterator() -> io::Result<()> {
        let path = path::Path::new("./test/data/small.mzML");
        let mut reader = MzMLReaderType::<_, CentroidPeak, DeconvolutedPeak>::open_path(path)?;
        let mut counter = 0;
        for scan in reader.iter().start_from_index(30).expect("Seek failed") {
            counter += 1;
            assert!(scan.index() >= 30);
        }
        let n = reader.len();
        assert_eq!(n, counter + 30);
        Ok(())
    }

    #[test]
    fn test_random_access_failure() -> io::Result<()> {
        let path = path::Path::new("./test/data/small.mzML");
        let mut reader = MzMLReaderType::<_, CentroidPeak, DeconvolutedPeak>::open_path(path)?;

        let offset = reader._offset_of_index(20).unwrap().saturating_sub(200);
        reader.seek(SeekFrom::Start(offset))?;

        let found = reader.check_stream("spectrum")?;
        assert!(!found);

        let offset = reader._offset_of_index(20).unwrap();
        reader.seek(SeekFrom::Start(offset))?;
        let found = reader.check_stream("spectrum")?;
        assert!(found);
        Ok(())
    }

    #[test]
    fn test_get_by_index() -> io::Result<()> {
        let path = path::Path::new("./test/data/small.mzML");
        let mut reader = MzMLReader::open_path(path)?;
        let scan = reader.get_spectrum_by_index(0).unwrap();
        assert_eq!(scan.index(), 0);

        let scan = reader.get_spectrum_by_index(10).unwrap();
        assert_eq!(scan.index(), 10);

        let scan = reader.get_spectrum_by_index(30).unwrap();
        assert_eq!(scan.index(), 30);
        Ok(())
    }

    #[test]
    fn test_with_detail_level() -> io::Result<()> {
        let path = path::Path::new("./test/data/small.mzML");
        let mut reader = MzMLReader::open_path(path)?;
        assert_eq!(*reader.detail_level(), DetailLevel::Full);

        let scan_full = reader.get_spectrum_by_index(0).unwrap();
        scan_full
            .arrays
            .as_ref()
            .unwrap()
            .iter()
            .for_each(|(_, v)| {
                assert!(matches!(v.compression, BinaryCompressionType::Decoded));
            });

        reader.set_detail_level(DetailLevel::Lazy);
        let scan_lazy = reader.get_spectrum_by_index(0).unwrap();
        scan_lazy
            .arrays
            .as_ref()
            .unwrap()
            .iter()
            .for_each(|(_, v)| {
                assert!(matches!(
                    v.compression,
                    BinaryCompressionType::NoCompression
                ));
            });

        reader.detail_level = DetailLevel::MetadataOnly;
        let scan_lazy = reader.get_spectrum_by_index(0).unwrap();
        scan_lazy
            .arrays
            .as_ref()
            .unwrap()
            .iter()
            .for_each(|(_, v)| {
                assert!(matches!(
                    v.compression,
                    BinaryCompressionType::NoCompression
                ));
                assert!(v.data.is_empty());
            });

        Ok(())
    }

    #[test]
    fn test_random_start() -> io::Result<()> {
        let path = path::Path::new("./test/data/batching_test.mzML");
        let mut reader = MzMLReader::open_path(path)?;

        let scan = reader
            .start_from_id("controllerType=0 controllerNumber=1 scan=25869")?
            .next()
            .unwrap();
        assert_eq!(scan.id(), "controllerType=0 controllerNumber=1 scan=25869");

        let scan2 = reader.start_from_index(scan.index())?.next().unwrap();
        assert_eq!(scan.index(), scan2.index());

        let scan2 = reader.start_from_time(scan.start_time())?.next().unwrap();
        assert_eq!(scan.start_time(), scan2.start_time());
        Ok(())
    }

    #[test]
    fn test_interleaved_groups() -> io::Result<()> {
        let path = path::Path::new("./test/data/batching_test.mzML");
        let mut reader = MzMLReader::open_path(path)?;
        let groups: Vec<_> = reader.groups().collect();
        let grp = &groups[10];
        let prec = grp.precursor().unwrap();
        assert_eq!(prec.id(), "controllerType=0 controllerNumber=1 scan=25869");
        assert_eq!(prec.index(), 129);
        let min_msn_idx = grp
            .products()
            .iter()
            .map(|scan| scan.index())
            .min()
            .unwrap();
        assert_eq!(min_msn_idx, 142);
        assert_eq!(groups.len(), 188);

        reader.reset();
        let mut it = reader.groups();

        RandomAccessSpectrumGroupingIterator::start_from_id(
            &mut it,
            "controllerType=0 controllerNumber=1 scan=25869",
        )?;
        let grp2 = it.next().unwrap();
        assert_eq!(
            grp2.precursor().map(|s| s.description()),
            grp.precursor().map(|s| s.description())
        );
        for (a, b) in grp2
            .products()
            .iter()
            .map(|s| s.description())
            .zip(grp.products().iter().map(|s| s.description()))
        {
            assert_eq!(*a, *b)
        }

        RandomAccessSpectrumGroupingIterator::start_from_index(
            &mut it,
            grp.precursor().unwrap().index(),
        )?;
        let grp2 = it.next().unwrap();
        assert_eq!(
            grp2.precursor().map(|s| s.description()),
            grp.precursor().map(|s| s.description())
        );
        for (a, b) in grp2
            .products()
            .iter()
            .map(|s| s.description())
            .zip(grp.products().iter().map(|s| s.description()))
        {
            assert_eq!(*a, *b)
        }

        RandomAccessSpectrumGroupingIterator::start_from_time(
            &mut it,
            grp.precursor().unwrap().start_time(),
        )?;
        let grp2 = it.next().unwrap();
        assert_eq!(
            grp2.precursor().map(|s| s.description()),
            grp.precursor().map(|s| s.description())
        );
        for (a, b) in grp2
            .products()
            .iter()
            .map(|s| s.description())
            .zip(grp.products().iter().map(|s| s.description()))
        {
            assert_eq!(*a, *b)
        }
        Ok(())
    }

    #[test_log::test]
    fn test_interleaved_into_groups() -> io::Result<()> {
        let path = path::Path::new("./test/data/batching_test.mzML");
        let reader = MzMLReader::open_path(path)?;
        let groups: Vec<_> = reader.into_groups().collect();
        let grp = &groups[10];
        let prec = grp.precursor().unwrap();
        assert_eq!(prec.id(), "controllerType=0 controllerNumber=1 scan=25869");
        assert_eq!(prec.index(), 129);
        let min_msn_idx = grp
            .products()
            .iter()
            .map(|scan| scan.index())
            .min()
            .unwrap();
        assert_eq!(min_msn_idx, 142);
        assert_eq!(groups.len(), 188);
        Ok(())
    }

    #[test_log::test]
    fn test_get_chromatogram() -> io::Result<()> {
        let path = path::Path::new("./test/data/batching_test.mzML");
        let mut reader = MzMLReader::open_path(path)?;
        let chrom = reader.get_chromatogram_by_id("TIC").unwrap();
        assert_eq!(chrom.id(), "TIC");
        assert_eq!(chrom.index(), 0);
        assert_eq!(chrom.arrays.len(), 3);

        let key = ArrayType::NonStandardDataArray {
            name: Box::new("ms level".to_string()),
        };
        let arr = chrom.arrays.get(&key).unwrap();
        let view = arr.to_i64().unwrap();
        assert_eq!(view.len(), 73368);
        assert_eq!(chrom.time().unwrap().len(), 73368);

        let chrom2 = reader.get_chromatogram_by_index(0).unwrap();
        assert_eq!(chrom2.id(), "TIC");
        assert_eq!(chrom2.time()?.len(), 73368);

        reader.reset();
        assert_eq!(reader.chromatogram_index.len(), 1);
        let chrom3 = reader.iter_chromatograms().next().unwrap();
        assert_eq!(chrom3.id(), "TIC");
        assert_eq!(chrom3.time()?.len(), 73368);

        let chrom4 = ChromatogramSource::get_chromatogram_by_index(&mut reader, 0).unwrap();
        assert_eq!(chrom4.id(), "TIC");
        assert_eq!(chrom4.time()?.len(), 73368);

        let chrom5 = ChromatogramSource::get_chromatogram_by_id(&mut reader, "TIC").unwrap();
        assert_eq!(chrom5.id(), "TIC");
        assert_eq!(chrom5.time()?.len(), 73368);

        Ok(())
    }

    #[cfg(feature = "mzsignal")]
    #[test_log::test]
    fn test_averaging() -> io::Result<()> {
        use crate::spectrum::group::SpectrumAveragingIterator;

        let path = path::Path::new("./test/data/small.mzML");
        let mut reader = MzMLReader::open_path(path)?;
        let groups: Vec<_> = reader.groups().collect();

        reader.reset();
        let averaging_iter =
            SpectrumAveragingIterator::new(reader.groups(), 1, 100.0, 2200.0, 0.005);
        let avg_groups: Vec<_> = averaging_iter.collect();
        assert_eq!(groups.len(), avg_groups.len());

        let raw_tic = &groups[0].precursor.as_ref().unwrap().peaks().tic();

        let avg_group = &avg_groups[0];
        let prec = avg_group.precursor().unwrap();
        let prec_arrays = prec.arrays.as_ref().unwrap();

        let avg_mzs = prec_arrays.mzs()?;
        let avg_intens = prec_arrays.intensities()?;

        assert_eq!(avg_mzs.len(), avg_intens.len());

        let low = avg_mzs.first().unwrap();
        let high = avg_mzs.last().unwrap();
        eprintln!("{low} {high}");
        assert!((low - 100.0).abs() < 1e-3);
        assert!((high - (2200.0 - 0.005)).abs() < 1e-3);

        let tic = prec.peaks().tic();
        let base_peak = prec.peaks().base_peak();
        eprintln!("{base_peak:?} {tic} {raw_tic}");

        Ok(())
    }

    #[cfg(feature = "mzsignal")]
    #[test_log::test]
    fn test_averaging_deferred() -> io::Result<()> {
        let path = path::Path::new("./test/data/small.mzML");
        let mut reader = MzMLReader::open_path(path)?;
        let groups: Vec<_> = reader.groups().collect();

        reader.reset();
        // let averaging_iter =
        //     SpectrumAveragingIterator::new(reader.groups(), 1, 100.0, 2200.0, 0.005);
        let (averaging_iter, mut averager, reprofiler) =
            reader.groups().averaging_deferred(1, 100.0, 2200.0, 0.005);
        let avg_groups: Vec<_> = averaging_iter
            .map(|grp| {
                let (mut grp, arrays) = grp.reprofile_with_average_with(&mut averager, &reprofiler);
                *grp.precursor_mut().unwrap().arrays.as_mut().unwrap() = arrays.into();
                grp
            })
            .take(2)
            .collect();

        assert_eq!(2, avg_groups.len());

        let raw_tic = &groups[0].precursor.as_ref().unwrap().peaks().tic();
        let avg_group = &avg_groups[0];
        let prec = avg_group.precursor().unwrap();
        let prec_arrays = prec.arrays.as_ref().unwrap();

        let avg_mzs = prec_arrays.mzs()?;
        let avg_intens = prec_arrays.intensities()?;

        assert_eq!(avg_mzs.len(), avg_intens.len());

        let low = avg_mzs.first().unwrap();
        let high = avg_mzs.last().unwrap();
        eprintln!("{low} {high}");
        assert!((low - 100.0).abs() < 1e-3);
        assert!((high - (2200.0 - 0.005)).abs() < 1e-3);

        let tic = prec.peaks().tic();
        let base_peak = prec.peaks().base_peak();
        eprintln!("{base_peak:?} {tic} {raw_tic}");

        Ok(())
    }

    #[cfg(feature = "mzsignal")]
    #[test_log::test]
    fn test_averaging_deferred_partial() -> io::Result<()> {
        let path = path::Path::new("./test/data/small.mzML");
        let mut reader = MzMLReader::open_path(path)?;
        let groups: Vec<_> = reader.groups().collect();

        reader.reset();
        // let averaging_iter =
        //     SpectrumAveragingIterator::new(reader.groups(), 1, 100.0, 2200.0, 0.005);
        let (averaging_iter, mut averager, reprofiler) = reader
            .groups()
            .map(|mut grp| {
                grp.precursor_mut().unwrap().pick_peaks(1.0).unwrap();
                grp.precursor_mut()
                    .unwrap()
                    .description_mut()
                    .signal_continuity = SignalContinuity::Centroid;
                grp
            })
            .averaging_deferred(1, 100.0, 2200.0, 0.005);
        let avg_groups: Vec<_> = averaging_iter
            .map(|grp| {
                let grp = grp.reprofile_with(&reprofiler, 0.002);
                let (mut grp, ctx) = grp.average_with(&mut averager);
                *grp.precursor_mut().unwrap().arrays.as_mut().unwrap() = ctx.into();
                grp
            })
            .take(2)
            .collect();

        assert_eq!(2, avg_groups.len());

        let raw_tic = &groups[0].precursor.as_ref().unwrap().peaks().tic();
        let avg_group = &avg_groups[0];
        let prec = avg_group.precursor().unwrap();
        let prec_arrays = prec.arrays.as_ref().unwrap();

        let avg_mzs = prec_arrays.mzs()?;
        let avg_intens = prec_arrays.intensities()?;

        assert_eq!(avg_mzs.len(), avg_intens.len());

        let low = avg_mzs.first().unwrap();
        let high = avg_mzs.last().unwrap();
        eprintln!("{low} {high}");
        assert!((low - 100.0).abs() < 1e-3);
        assert!((high - (2200.0 - 0.005)).abs() < 1e-3);

        let tic = prec.peaks().tic();
        let base_peak = prec.peaks().base_peak();
        eprintln!("{base_peak:?} {tic} {raw_tic}");

        Ok(())
    }

    #[test_log::test]
    fn test_iterator_specialization() -> io::Result<()> {
        let path = path::Path::new("./test/data/small.mzML");
        let mut reader = MzMLReader::open_path(path)?;

        let spec = reader.nth(10).unwrap();
        reader.reset();
        let spec2 = reader.iter().nth(10).unwrap();
        assert_eq!(spec.id(), spec2.id());
        Ok(())
    }

    #[test]
    fn test_spectrum_builder() {
        let mut builder: MzMLSpectrumBuilder<'_, CentroidPeak, DeconvolutedPeak> =
            MzMLSpectrumBuilder::new();
        assert!(builder.is_spectrum_entry());
        assert!(!builder.is_chromatogram_entry());
        assert_eq!(builder.entry_type(), EntryType::Spectrum);
        assert_eq!(builder.warning_context(), "spectrum entry 0 ()");
        builder.set_entry_type(EntryType::Chromatogram);
        assert_eq!(builder.warning_context(), "chromatogram entry 0 ()");
        builder._reset();

        builder.fill_binary_data_array(ControlledVocabulary::MS.param(1002312, "numpress linear"));
        assert_eq!(
            builder.current_array.compression,
            BinaryCompressionType::NumpressLinear
        );

        let pairs = [
            (1000574, BinaryCompressionType::Zlib),
            (1000576, BinaryCompressionType::NoCompression),
            (1002312, BinaryCompressionType::NumpressLinear),
            (1002313, BinaryCompressionType::NumpressPIC),
            (1002314, BinaryCompressionType::NumpressSLOF),
            (1002746, BinaryCompressionType::NumpressLinearZlib),
            (1002747, BinaryCompressionType::NumpressPICZlib),
            (1002748, BinaryCompressionType::NumpressSLOFZlib),
            (1003089, BinaryCompressionType::DeltaPrediction),
            (1003090, BinaryCompressionType::LinearPrediction),
        ];

        for (acc, term) in pairs {
            builder.fill_binary_data_array(ControlledVocabulary::MS.param(acc, term.to_string()));
            assert_eq!(builder.current_array.compression, term);
            builder._reset();
        }

        let pairs = [
            (1000523, BinaryDataArrayType::Float64),
            (1000521, BinaryDataArrayType::Float32),
            (1000522, BinaryDataArrayType::Int64),
            (1000519, BinaryDataArrayType::Int32),
            (1001479, BinaryDataArrayType::ASCII),
        ];

        for (acc, term) in pairs {
            builder.fill_binary_data_array(ControlledVocabulary::MS.param(acc, term.to_string()));
            assert_eq!(builder.current_array.dtype, term);
            builder._reset();
        }
        let pairs = [
            (1002477, ArrayType::MeanDriftTimeArray, Unit::Millisecond),
            (
                1003006,
                ArrayType::MeanInverseReducedIonMobilityArray,
                Unit::VoltSecondPerSquareCentimeter,
            ),
            (1003007, ArrayType::RawIonMobilityArray, Unit::Unknown),
            (1003153, ArrayType::RawDriftTimeArray, Unit::Unknown),
            (
                1003156,
                ArrayType::DeconvolutedDriftTimeArray,
                Unit::Millisecond,
            ),
            (
                1003008,
                ArrayType::RawInverseReducedIonMobilityArray,
                Unit::VoltSecondPerSquareCentimeter,
            ),
            (
                1003154,
                ArrayType::DeconvolutedDriftTimeArray,
                Unit::Millisecond,
            ),
            (
                1003155,
                ArrayType::DeconvolutedInverseReducedIonMobilityArray,
                Unit::VoltSecondPerSquareCentimeter,
            ),
        ];

        for (acc, term, unit) in pairs {
            builder.fill_binary_data_array(
                ControlledVocabulary::MS
                    .param(acc, term.to_string())
                    .with_unit_t(&unit),
            );
            assert_eq!(builder.current_array.name, term);
            assert_eq!(builder.current_array.unit, unit);
            builder._reset();
        }

        let param = ControlledVocabulary::MS.const_param(
            "ion injection time",
            crate::params::ValueRef::Float(50.0),
            0,
            Unit::Millisecond,
        );
        builder.fill_param_into(param.into(), MzMLParserState::Scan);
        assert_eq!(
            builder.acquisition.first_scan().unwrap().injection_time,
            50.0
        );

        let param = ControlledVocabulary::MS.const_param(
            "scan start time",
            crate::params::ValueRef::Float(50.0),
            0,
            Unit::Minute,
        );
        builder.fill_param_into(param.into(), MzMLParserState::Scan);
        assert_eq!(builder.acquisition.first_scan().unwrap().start_time, 50.0);

        let param = ScanCombination::NoCombination.to_param();
        builder.fill_param_into(param, MzMLParserState::ScanList);
        assert_eq!(
            builder.acquisition.combination,
            ScanCombination::NoCombination
        );

        builder._reset();

        let param = ControlledVocabulary::MS.const_param(
            "isolation window target m/z",
            crate::params::ValueRef::Float(50.0),
            0,
            Unit::MZ,
        );
        builder.fill_param_into(param.into(), MzMLParserState::IsolationWindow);
        assert_eq!(builder.isolation_window_mut().target, 50.0);

        builder._reset();
        let param = ControlledVocabulary::MS.const_param(
            "isolation window lower limit",
            crate::params::ValueRef::Float(48.0),
            0,
            Unit::MZ,
        );
        builder.fill_param_into(param.into(), MzMLParserState::IsolationWindow);
        assert_eq!(builder.isolation_window_mut().lower_bound, 48.0);
        let param = ControlledVocabulary::MS.const_param(
            "isolation window upper limit",
            crate::params::ValueRef::Float(52.0),
            0,
            Unit::MZ,
        );
        builder.fill_param_into(param.into(), MzMLParserState::IsolationWindow);
        assert_eq!(builder.isolation_window_mut().upper_bound, 52.0);
    }
}
