use std::{collections::HashMap, io, marker::PhantomData, mem, path::PathBuf};

use chrono::DateTime;
use log::{debug, warn};

use crate::{
    io::{
        traits::ChromatogramSource, utils::checksum_file, DetailLevel,
        Generic3DIonMobilityFrameSource, OffsetIndex,
    },
    meta::{
        Component, ComponentType, DataProcessing, DetectorTypeTerm, DissociationMethodTerm,
        FileDescription, InstrumentConfiguration, IonizationTypeTerm, MassAnalyzerTerm,
        MassSpectrometryRun, Sample, Software, SourceFile,
    },
    params::{ControlledVocabulary, Unit, Value},
    prelude::*,
    spectrum::{
        ArrayType, BinaryArrayMap, BinaryDataArrayType, CentroidPeakAdapting, Chromatogram,
        ChromatogramDescription, ChromatogramType, DataArray, DeconvolutedPeakAdapting,
        MultiLayerSpectrum, Precursor, ScanEvent, ScanPolarity, ScanWindow, SelectedIon,
        SignalContinuity,
    },
    Param,
};

use mzpeaks::{peak_set::PeakSetVec, prelude::*, CentroidPeak, DeconvolutedPeak, MZ};

use super::instruments::{
    create_instrument_configurations, instrument_model_to_ion_sources,
    instrument_model_to_mass_analyzers,
};
#[allow(unused)]
use super::instruments::{
    instrument_model_to_detector, parse_instrument_model, InstrumentModelType,
};

macro_rules! param {
    ($name:expr, $acc:expr) => {
        ControlledVocabulary::MS.const_param_ident($name, $acc)
    };
}

/// Check to see if a buffer contains the header of a Thermo RAW file
///
/// Thermo RAW files start with a UTF-16 header with "Finnigan" at
/// codepoints 1-9.
pub fn is_thermo_raw_prefix(buffer: &[u8]) -> bool {
    if buffer.len() < 18 {
        debug!("Attempted to test a byte buffer for the Thermo prefix, buffer was less than 18 bytes long");
        return false;
    }
    let view: &[u16] = unsafe { mem::transmute(&buffer[2..18]) };
    let prefix = String::from_utf16_lossy(view);
    prefix.starts_with("Finnigan")
}

impl<C: CentroidLike + From<CentroidPeak>, D: DeconvolutedCentroidLike>
    MZFileReader<C, D, MultiLayerSpectrum<C, D>> for ThermoRawReaderType<C, D>
{
    /// The format is always indexed and requires no formal build step.
    fn construct_index_from_stream(&mut self) -> u64 {
        self.len() as u64
    }

    /// The underlying Thermo library requires an explicit file system path to open
    /// the file, and as such this method always fails.
    ///
    /// [`open_path`](Self::open_path) works as normal.
    #[allow(unused)]
    fn open_file(source: std::fs::File) -> io::Result<Self> {
        Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "Cannot read a Thermo RAW file from an open file handle, only directly from a path",
        ))
    }

    fn open_path<P>(path: P) -> io::Result<Self>
    where
        P: Into<std::path::PathBuf> + Clone,
    {
        Self::new(path.into())
    }
}

#[inline(always)]
fn make_native_id(index: i32) -> String {
    format!(
        "controllerType=0 controllerNumber=1 scan={}",
        index as usize + 1
    )
}

const SOURCE_FILE_ID: & str = "RAW1";

#[cfg(not(feature = "doc-only"))]
pub(crate) mod sealed {
    use std::path::Path;

    use crate::spectrum::bindata::to_bytes;

    use super::*;
    use thermorawfilereader::{
        schema::{
            AcquisitionT, DissociationMethod, Polarity, PrecursorT, SpectrumData, SpectrumMode,
        },
        ExtendedSpectrumData, FileDescription as ThermoFileDescription, IonizationMode,
        MassAnalyzer, RawFileReader,
    };

    /**
        A Thermo Fisher RAW file reader that supports iteration and random access.
    */
    pub struct ThermoRawReaderType<
        C: CentroidLike + From<CentroidPeak> = CentroidPeak,
        D: DeconvolutedCentroidLike = DeconvolutedPeak,
    > {
        pub path: PathBuf,
        pub detail_level: DetailLevel,
        pub(crate) handle: RawFileReader,
        pub(crate) index: usize,
        pub(crate) spectrum_index: OffsetIndex,
        pub(crate) file_description: FileDescription,
        pub(crate) instrument_configurations: HashMap<u32, InstrumentConfiguration>,
        pub(crate) components_to_instrument_id: HashMap<MassAnalyzer, u32>,
        pub(crate) softwares: Vec<Software>,
        pub(crate) samples: Vec<Sample>,
        pub(crate) data_processings: Vec<DataProcessing>,
        pub(crate) ms_run: MassSpectrometryRun,
        _c: PhantomData<C>,
        _d: PhantomData<D>,
        load_extended_spectrum_data: bool,
    }

    // The public API
    impl<C: CentroidLike + From<CentroidPeak>, D: DeconvolutedCentroidLike>
        ThermoRawReaderType<C, D>
    {
        /// Create a new [`ThermoRawReaderType`] from a path.
        /// This may trigger an expensive I/O operation to checksum the file
        pub fn new_with_detail_level_and_centroiding<P: Into<PathBuf>>(
            path: P,
            mut detail_level: DetailLevel,
            centroiding: bool,
        ) -> io::Result<Self> {
            let path: PathBuf = path.into();
            let mut handle = RawFileReader::open(&path)?;
            handle.set_centroid_spectra(centroiding);

            if matches!(detail_level, DetailLevel::Lazy) {
                log::warn!(
                    "ThermoRawReader does not support lazy loading. Using {:?}",
                    DetailLevel::Full
                );
                detail_level = DetailLevel::Full
            }

            let spectrum_index = Self::build_index(&handle);

            let thermo_file_description: ThermoFileDescription = handle.file_description();

            let file_description = Self::make_file_description(&path, &thermo_file_description)?;

            let (sw, instrument_configurations, components_to_instrument_id) =
                Self::make_instrument_configuration(&handle);

            let ms_run = Self::make_ms_run(&path, &thermo_file_description);
            let samples = Self::make_sample(&thermo_file_description)
                .into_iter()
                .collect();

            Ok(Self {
                path,
                detail_level,
                handle,
                index: 0,
                spectrum_index,
                file_description,
                instrument_configurations,
                samples,
                components_to_instrument_id,
                softwares: vec![sw],
                data_processings: vec![],
                ms_run,
                _c: PhantomData,
                _d: PhantomData,
                load_extended_spectrum_data: false,
            })
        }

        /// Create a new [`ThermoRawReaderType`] from a path.
        /// This may trigger an expensive I/O operation to checksum the file
        pub fn new<P: Into<PathBuf>>(path: P) -> io::Result<Self> {
            Self::new_with_detail_level_and_centroiding(path, DetailLevel::Full, false)
        }

        pub fn len(&self) -> usize {
            self.handle.len()
        }

        pub fn is_empty(&self) -> bool {
            self.handle.is_empty()
        }

        /// Get whether or not to centroid spectra on read using the vendor algorithm
        pub fn get_centroiding(&self) -> bool {
            self.handle.get_centroid_spectra()
        }

        /// Set whether or not to centroid spectra on read using the vendor algorithm
        pub fn set_centroiding(&mut self, value: bool) {
            self.handle.set_centroid_spectra(value)
        }

        fn unpack_chromatogram_signal(
            &self,
            descr: thermorawfilereader::ChromatogramDescription,
        ) -> BinaryArrayMap {
            let mut array_map = BinaryArrayMap::default();
            if let Some(data) = descr.data() {
                let time_array = data.time();
                let intensity_array = data.intensity();

                let mut time_array_in = DataArray::from_name_type_size(
                    &ArrayType::TimeArray,
                    BinaryDataArrayType::Float64,
                    time_array.len() * BinaryDataArrayType::Float64.size_of(),
                );
                time_array_in.unit = Unit::Minute;
                time_array_in.extend(&time_array).unwrap();
                array_map.add(time_array_in);

                let mut intensity_array_in = DataArray::from_name_type_size(
                    &ArrayType::IntensityArray,
                    BinaryDataArrayType::Float32,
                    time_array.len() * BinaryDataArrayType::Float32.size_of(),
                );
                intensity_array_in.unit = Unit::DetectorCounts;
                intensity_array_in.extend(&intensity_array).unwrap();
                array_map.add(intensity_array_in);
            } else {
                let mut time_array_in = DataArray::from_name_and_type(
                    &ArrayType::TimeArray,
                    BinaryDataArrayType::Float64,
                );
                time_array_in.unit = Unit::Minute;
                array_map.add(time_array_in);

                let mut intensity_array_in = DataArray::from_name_and_type(
                    &ArrayType::IntensityArray,
                    BinaryDataArrayType::Float32,
                );
                intensity_array_in.unit = Unit::DetectorCounts;
                array_map.add(intensity_array_in);
            }

            array_map
        }

        pub fn get_tic(&mut self) -> Chromatogram {
            let tic = self.handle.tic();
            let array_map = self.unpack_chromatogram_signal(tic);

            let descr = ChromatogramDescription {
                chromatogram_type: ChromatogramType::TotalIonCurrentChromatogram,
                id: "TIC".to_string(),
                ..Default::default()
            };
            Chromatogram::new(descr, array_map)
        }

        pub fn get_bpc(&mut self) -> Chromatogram {
            let bpc = self.handle.bpc();

            let array_map = self.unpack_chromatogram_signal(bpc);
            let descr = ChromatogramDescription {
                id: "BPC".to_string(),
                chromatogram_type: ChromatogramType::BasePeakChromatogram,
                ..Default::default()
            };
            Chromatogram::new(descr, array_map)
        }

        /// Get whether or not to load extended spectrum signal information for the spectrum.
        ///
        /// The loaded data isn't incorporated into a peak list, instead access them under
        /// the binary data arrays.
        pub fn get_load_extended_spectrum_data(&self) -> bool {
            self.load_extended_spectrum_data
        }

        /// Set whether or not to load extended spectrum signal information for the spectrum.
        ///
        /// The loaded data isn't incorporated into a peak list, instead access them under
        /// the binary data arrays.
        pub fn set_load_extended_spectrum_data(&mut self, load_extended_spectrum_data: bool) {
            self.load_extended_spectrum_data = load_extended_spectrum_data;
        }

        /// Enumerate the available status logs' names
        pub fn status_log_names(&self) -> impl Iterator<Item = String> {
            let it = self
                .handle
                .get_status_logs()
                .into_iter()
                .next()
                .map(|logs| {
                    let names: Vec<_> = logs
                        .str_logs()
                        .map(|log| log.name.clone())
                        .chain(logs.float_logs().map(|log| log.name.clone()))
                        .chain(logs.int_logs().map(|log| log.name.clone()))
                        .chain(logs.bool_logs().map(|log| log.name.clone()))
                        .collect();
                    names.into_iter()
                })
                .into_iter()
                .flatten();
            it
        }

        /// Get a status log as a [`Chromatogram`]
        pub fn get_status_log_trace_by_name(&self, trace_id: &str) -> Option<Chromatogram> {
            let temperature_pattern = regex::Regex::new(r#"Temp\S?\s|Temperature|°C"#).unwrap();
            if let Some(logs) = self.handle.get_status_logs() {
                macro_rules! make_description {
                    ($descr:ident, $log:ident) => {
                        $descr.id = $log.name.clone();
                        if temperature_pattern.is_match(&($log.name)) {
                            $descr.chromatogram_type = ChromatogramType::TemperatureChromatogram;
                        }
                    };
                }
                macro_rules! make_arrays {
                    ($log:ident, $arrays:ident) => {
                        let times = $log.times();

                        let mut time_array = DataArray::from_name_type_size(
                            &ArrayType::TimeArray,
                            BinaryDataArrayType::Float64,
                            times.len(),
                        );
                        time_array.unit = Unit::Minute;
                        time_array.extend(&times).unwrap();

                        $arrays.add(time_array);
                    };
                }
                let chrom = logs
                    .str_logs()
                    .filter_map(|log| -> Option<Chromatogram> {
                        if log.name == trace_id {
                            let mut descr = ChromatogramDescription::default();
                            let mut arrays = BinaryArrayMap::new();
                            make_description!(descr, log);
                            make_arrays!(log, arrays);

                            let mut str_array = DataArray::from_name_and_type(
                                &ArrayType::nonstandard(&log.name),
                                BinaryDataArrayType::ASCII,
                            );
                            for s in log.strings() {
                                str_array.extend(s.as_bytes()).unwrap();
                                str_array.push(b'\0').unwrap();
                            }
                            arrays.add(str_array);
                            Some(Chromatogram::new(descr, arrays))
                        } else {
                            None
                        }
                    })
                    .chain(logs.float_logs().filter_map(|log| {
                        if log.name == trace_id {
                            let mut descr = ChromatogramDescription::default();
                            let mut arrays = BinaryArrayMap::new();
                            make_description!(descr, log);
                            make_arrays!(log, arrays);

                            let data = log.values();
                            let mut array = DataArray::from_name_type_size(
                                &ArrayType::nonstandard(&log.name),
                                BinaryDataArrayType::Float64,
                                data.len(),
                            );

                            array.extend(&data).unwrap();
                            arrays.add(array);
                            Some(Chromatogram::new(descr, arrays))
                        } else {
                            None
                        }
                    }))
                    .chain(logs.int_logs().flat_map(|log| {
                        if log.name == trace_id {
                            let mut descr = ChromatogramDescription::default();
                            let mut arrays = BinaryArrayMap::new();
                            make_description!(descr, log);
                            make_arrays!(log, arrays);

                            let data = log.values();
                            let mut array = DataArray::from_name_type_size(
                                &ArrayType::nonstandard(&log.name),
                                BinaryDataArrayType::Int64,
                                data.len(),
                            );

                            array.extend(&data).unwrap();
                            arrays.add(array);
                            Some(Chromatogram::new(descr, arrays))
                        } else {
                            None
                        }
                    }))
                    .next();
                chrom
            } else {
                None
            }
        }
    }

    fn ionization_mode_to_inlet(value: IonizationMode) -> Option<Param> {
        match value {
            IonizationMode::CardNanoSprayIonization | IonizationMode::NanoSpray => {
                Some(param!("nanospray inlet", 1000485).into())
            }
            IonizationMode::ElectroSpray => Some(param!("electrospray inlet", 1000057).into()),
            IonizationMode::ThermoSpray => Some(param!("thermospray inlet", 1000069).into()),
            IonizationMode::FastAtomBombardment => {
                Some(param!("continuous flow fast atom bombardment", 1000055).into())
            }
            _ => None,
        }
    }

    fn translate_ionization_mode(value: IonizationMode) -> Param {
        match value {
            IonizationMode::CardNanoSprayIonization | IonizationMode::NanoSpray => {
                IonizationTypeTerm::Nanoelectrospray
            }
            IonizationMode::ElectroSpray => IonizationTypeTerm::ElectrosprayIonization,
            IonizationMode::AtmosphericPressureChemicalIonization => {
                IonizationTypeTerm::AtmosphericPressureChemicalIonization
            }
            IonizationMode::FastAtomBombardment => {
                IonizationTypeTerm::FastAtomBombardmentIonization
            }
            IonizationMode::GlowDischarge => IonizationTypeTerm::GlowDischargeIonization,
            IonizationMode::ElectronImpact => IonizationTypeTerm::ElectronIonization,
            IonizationMode::MatrixAssistedLaserDesorptionIonization => {
                IonizationTypeTerm::MatrixAssistedLaserDesorptionIonization
            }
            IonizationMode::ChemicalIonization => IonizationTypeTerm::ChemicalIonization,
            _ => IonizationTypeTerm::IonizationType,
        }
        .to_param()
        .into()
    }

    fn translate_mass_analyzer(value: MassAnalyzer) -> Param {
        match value {
            MassAnalyzer::ITMS => MassAnalyzerTerm::RadialEjectionLinearIonTrap,
            MassAnalyzer::FTMS => MassAnalyzerTerm::Orbitrap,
            MassAnalyzer::ASTMS => MassAnalyzerTerm::AsymmetricTrackLosslessTimeOfFlightAnalyzer,
            MassAnalyzer::TOFMS => MassAnalyzerTerm::TimeOfFlight,
            // MassAnalyzer::TQMS => {}
            // MassAnalyzer::SQMS => {}
            MassAnalyzer::Sector => MassAnalyzerTerm::MagneticSector,
            _ => MassAnalyzerTerm::MassAnalyzerType,
        }
        .to_param()
        .into()
    }

    fn translate_mass_analyzer_detector(value: MassAnalyzer) -> Param {
        match value {
            MassAnalyzer::ITMS | MassAnalyzer::ASTMS => DetectorTypeTerm::ElectronMultiplier,
            MassAnalyzer::FTMS => DetectorTypeTerm::InductiveDetector,
            _ => DetectorTypeTerm::DetectorType,
            // MassAnalyzer::TOFMS => {}
            // MassAnalyzer::TQMS => {}
            // MassAnalyzer::SQMS => {}
            // MassAnalyzer::Sector => {}
        }
        .to_param()
        .into()
    }

    fn translate_mass_analyzer_reverse(value: &MassAnalyzerTerm) -> MassAnalyzer {
        match value {
            MassAnalyzerTerm::RadialEjectionLinearIonTrap => MassAnalyzer::ITMS,
            MassAnalyzerTerm::Orbitrap => MassAnalyzer::FTMS,
            MassAnalyzerTerm::AsymmetricTrackLosslessTimeOfFlightAnalyzer => MassAnalyzer::ASTMS,
            MassAnalyzerTerm::TimeOfFlight => MassAnalyzer::TOFMS,
            // { => MassAnalyzer::TQMS}
            // { => MassAnalyzer::SQMS}
            MassAnalyzerTerm::MagneticSector => MassAnalyzer::Sector,
            _ => MassAnalyzer::Unknown,
        }
    }

    impl<C: CentroidLike + From<CentroidPeak>, D: DeconvolutedCentroidLike>
        ThermoRawReaderType<C, D>
    {
        pub(crate) fn make_ms_run(
            path: &Path,
            thermo_file_description: &ThermoFileDescription,
        ) -> MassSpectrometryRun {
            let run = MassSpectrometryRun {
                default_instrument_id: Some(0),
                default_source_file_id: Some(SOURCE_FILE_ID.to_string()),
                id: path
                .file_name()
                .map(|s| s.to_string_lossy().split(".").next().unwrap().to_string()),
                start_time: thermo_file_description.creation_date().map(|s| DateTime::parse_from_rfc3339(s).unwrap()),
                ..Default::default()
            };
            run
        }

        pub(crate) fn make_file_description(
            path: &PathBuf,
            thermo_file_description: &ThermoFileDescription,
        ) -> io::Result<FileDescription> {
            let mut sf = SourceFile::default();
            let description = thermo_file_description;
            sf.name = path.file_name().unwrap().to_string_lossy().to_string();
            sf.location = format!(
                "file:///{}",
                path.canonicalize()?.parent().unwrap().display()
            );
            sf.id = SOURCE_FILE_ID.to_string();
            sf.file_format = Some(
                ControlledVocabulary::MS
                    .const_param_ident("Thermo RAW format", 1000563)
                    .into(),
            );
            sf.id_format = Some(
                ControlledVocabulary::MS
                    .const_param_ident("Thermo nativeID format", 1000768)
                    .into(),
            );
            sf.add_param(ControlledVocabulary::MS.param_val(
                1000569,
                "SHA-1",
                checksum_file(path)?,
            ));

            let levels: Vec<_> = description
                .spectra_per_ms_level()
                .into_iter()
                .flatten()
                .collect();

            let mut contents = Vec::new();
            if levels.first().copied().unwrap_or_default() > 0 {
                contents.push(param!("MS1 spectrum", 1000579).into())
            }
            if levels[1..].iter().copied().sum::<u32>() > 0 {
                contents.push(param!("MSn spectrum", 1000580).into())
            }

            let file_description = FileDescription::new(contents, vec![sf]);
            Ok(file_description)
        }

        pub(crate) fn make_instrument_configuration(
            handle: &RawFileReader,
        ) -> (
            Software,
            HashMap<u32, InstrumentConfiguration>,
            HashMap<MassAnalyzer, u32>,
        ) {
            let descr = handle.instrument_model();

            let mut sw = Software {
                id: "thermo_xcalibur".to_string(),
                version: descr.software_version().unwrap().to_string(),
                ..Default::default()
            };
            sw.add_param(
                ControlledVocabulary::MS
                    .const_param_ident("Xcalibur", 1000532)
                    .into(),
            );

            let model_type = if let Some(model_name) = descr.model() {
                parse_instrument_model(model_name)
            } else {
                InstrumentModelType::Unknown
            };

            let mut configs = HashMap::new();
            let mut components_to_instrument_id = HashMap::new();

            let method_texts: Vec<Param> = (0..handle.instrument_method_count())
                .flat_map(|i| handle.instrument_method(i as u8))
                .flat_map(|m| {
                    m.text().map(|s| {
                        ControlledVocabulary::MS.param_val(1000032, "customization", s.to_string())
                    })
                })
                .collect();

            let serial_number_param = descr.serial_number().map(|serial| ControlledVocabulary::MS.param_val(
                    1000529,
                    "instrument serial number",
                    serial,
                ));

            // Try to build the instrument configuration from the metadata
            for (i, vconf) in descr.configurations().enumerate() {
                let mut config = InstrumentConfiguration::default();
                if vconf.mass_analyzer == MassAnalyzer::Unknown {
                    continue;
                }
                config.extend_params(
                    method_texts
                        .iter()
                        .cloned()
                        .chain(serial_number_param.clone())
                        .chain([model_type.to_param()]),
                );

                let ion_source = config.new_component(ComponentType::IonSource);
                ion_source.extend_params(ionization_mode_to_inlet(vconf.ionization_mode));
                ion_source.add_param(translate_ionization_mode(vconf.ionization_mode));

                let analyzer = config.new_component(ComponentType::Analyzer);
                analyzer.add_param(translate_mass_analyzer(vconf.mass_analyzer));

                let detector = config.new_component(ComponentType::Detector);
                detector.add_param(translate_mass_analyzer_detector(vconf.mass_analyzer));

                config.software_reference = sw.id.clone();
                config.id = i as u32;

                components_to_instrument_id.insert(vconf.mass_analyzer, i as u32);

                configs.insert(i as u32, config);
            }

            // If the configurations weren't detectable from the top level metadata,
            // try to guess them from the instrument model.
            if configs.is_empty() && !handle.is_empty() {
                debug!("Using instrument mode {model_type} to infer configurations");
                if let Some(first_ionization) = handle.get(0).and_then(|s| {
                    s.acquisition()
                        .map(|acq| -> thermorawfilereader::IonizationMode {
                            acq.ionization_mode().0.into()
                        })
                }) {
                    let mut source = Component {
                        component_type: ComponentType::IonSource,
                        ..Default::default()
                    };
                    source.extend_params(ionization_mode_to_inlet(first_ionization));
                    source.add_param(translate_ionization_mode(first_ionization));

                    for (i, mut config) in create_instrument_configurations(model_type, source)
                        .into_iter()
                        .enumerate()
                    {
                        config.extend_params(
                            method_texts
                                .iter()
                                .cloned()
                                .chain([model_type.to_param()])
                                .chain(serial_number_param.clone()),
                        );
                        config.id = i as u32;
                        config.software_reference = sw.id.clone();

                        if let Some(mass_analyzer) =
                            config.iter().rev().flat_map(|c| c.mass_analyzer()).next()
                        {
                            let vconf_mass_analyzer =
                                translate_mass_analyzer_reverse(&mass_analyzer);
                            components_to_instrument_id.insert(vconf_mass_analyzer, i as u32);
                        } else {
                            warn!("Failed to locate mass analyzer from {config:?}")
                        }

                        configs.insert(i as u32, config);
                    }
                }
            }

            // If the whole configurations weren't specified by the instrument model,
            // try to guess them piece-meal.
            if configs.is_empty() {
                debug!("Using instrument mode {model_type} to infer configurations by parts");
                let mass_analyzers = instrument_model_to_mass_analyzers(model_type);
                let ionization_types = instrument_model_to_ion_sources(model_type);
                let detectors = instrument_model_to_detector(model_type);
                let mut i = 0;
                for ionization in ionization_types.iter() {
                    for (mass_analyzer, detector_type) in
                        mass_analyzers.iter().zip(detectors.iter())
                    {
                        let mut config = InstrumentConfiguration::default();

                        config.extend_params(
                            method_texts
                                .iter()
                                .cloned()
                                .chain(serial_number_param.clone())
                                .chain([model_type.to_param()]),
                        );

                        let ion_source = config.new_component(ComponentType::IonSource);
                        ion_source.add_param(ionization.into());

                        let analyzer = config.new_component(ComponentType::Analyzer);
                        analyzer.add_param(mass_analyzer.into());

                        let detector = config.new_component(ComponentType::Detector);
                        detector.add_param(detector_type.into());

                        config.software_reference = sw.id.clone();
                        config.id = i;

                        let vconf_mass_analyzer = translate_mass_analyzer_reverse(mass_analyzer);
                        components_to_instrument_id.insert(vconf_mass_analyzer, i);

                        configs.insert(i, config);
                        i += 1;
                    }
                }
            }
            if configs.is_empty() {
                log::warn!("No instrument configurations were found in Thermo RAW file?")
            }
            (sw, configs, components_to_instrument_id)
        }

        pub(crate) fn build_index(handle: &RawFileReader) -> OffsetIndex {
            let mut spectrum_index: OffsetIndex = OffsetIndex::new("spectrum".to_string());
            (0..handle.len()).for_each(|i| {
                spectrum_index.insert(make_native_id(i as i32), i as u64);
            });
            spectrum_index
        }

        pub(crate) fn make_sample(
            thermo_file_description: &ThermoFileDescription,
        ) -> Option<Sample> {
            thermo_file_description.sample_id().map(|name| Sample::new(
                    name.to_string(),
                    Some(name.to_string()),
                    Vec::new(),
                ))
        }

        pub(crate) fn populate_precursor(&self, vprec: &PrecursorT, precursor: &mut Precursor) {
            let mut ion = SelectedIon {
                mz: vprec.mz(),
                intensity: vprec.intensity(),
                ..Default::default()
            };
            ion.intensity = vprec.intensity();
            ion.charge = match vprec.charge() {
                0 => None,
                z => Some(z),
            };
            *precursor.ion_mut() = ion;

            let activation = &mut precursor.activation;
            let vact = vprec.activation();
            activation.energy = vact.collision_energy() as f32;
            match vact.dissociation_method() {
                DissociationMethod::CID => {
                    activation
                        .methods_mut()
                        .push(DissociationMethodTerm::CollisionInducedDissociation);
                }
                DissociationMethod::HCD => {
                    activation
                        .methods_mut()
                        .push(DissociationMethodTerm::BeamTypeCollisionInducedDissociation);
                }
                DissociationMethod::ECD => {
                    activation
                        .methods_mut()
                        .push(DissociationMethodTerm::ElectronCaptureDissociation);
                }
                DissociationMethod::ETD => {
                    activation
                        .methods_mut()
                        .push(DissociationMethodTerm::ElectronTransferDissociation);
                }
                DissociationMethod::ETHCD => {
                    activation
                        .methods_mut()
                        .push(DissociationMethodTerm::ElectronTransferDissociation);
                    activation.add_param(
                        DissociationMethodTerm::SupplementalBeamTypeCollisionInducedDissociation
                            .into(),
                    );
                }
                DissociationMethod::ETCID => {
                    activation
                        .methods_mut()
                        .push(DissociationMethodTerm::ElectronTransferDissociation);
                    activation.methods_mut().push(
                        DissociationMethodTerm::SupplementalCollisionInducedDissociation,
                    );
                }
                DissociationMethod::NETD => {
                    activation
                        .methods_mut()
                        .push(DissociationMethodTerm::NegativeElectronTransferDissociation);
                }
                DissociationMethod::MPD => {
                    todo!("Need to define MPD")
                }
                DissociationMethod::PTD => {
                    todo!("Need to define PTD")
                }
                DissociationMethod::ECCID => {
                    activation
                        .methods_mut()
                        .push(DissociationMethodTerm::ElectronCaptureDissociation);
                    activation.add_param(
                        DissociationMethodTerm::SupplementalCollisionInducedDissociation.into(),
                    );
                }
                DissociationMethod::ECHCD => {
                    activation
                        .methods_mut()
                        .push(DissociationMethodTerm::ElectronCaptureDissociation);
                    activation.add_param(
                        DissociationMethodTerm::SupplementalBeamTypeCollisionInducedDissociation
                            .into(),
                    )
                }
                _ => {
                    warn!("No activation translation found for {:?}", vact);
                    activation
                        .methods_mut()
                        .push(DissociationMethodTerm::CollisionInducedDissociation);
                }
            }

            let iso_window = &mut precursor.isolation_window;
            let vwin = vprec.isolation_window();
            iso_window.lower_bound = vwin.lower() as f32;
            iso_window.target = vwin.target() as f32;
            iso_window.upper_bound = vwin.upper() as f32;

            precursor.precursor_id = Some(make_native_id(vprec.parent_index()));
        }

        fn populate_scan_event(&self, vevent: &AcquisitionT, event: &mut ScanEvent) {
            event.injection_time = vevent.injection_time();
            let window = ScanWindow::new(vevent.low_mz() as f32, vevent.high_mz() as f32);
            event.scan_windows.push(window);
            if let Some(cv) = vevent.compensation_voltage() {
                let mut param =
                    ControlledVocabulary::MS.param_val(1001581, "FAIMS compensation voltage", cv);
                param.unit = Unit::Volt;
                event.add_param(param);
            }
            if let Some(resolution) = vevent.resolution() {
                event.add_param(ControlledVocabulary::MS.param_val(
                    1000011,
                    "mass resolution",
                    resolution,
                ));
            }
            event.add_param(ControlledVocabulary::MS.param_val(
                1000616,
                "preset scan configuration",
                vevent.scan_event(),
            ));
            let ic_key = vevent.mass_analyzer().0.into();
            event.instrument_configuration_id =
                self.get_instrument_configuration_by_mass_analyzer(ic_key);
        }

        fn get_instrument_configuration_by_mass_analyzer(
            &self,
            mass_analyzer: MassAnalyzer,
        ) -> u32 {
            *self
                .components_to_instrument_id
                .get(&mass_analyzer)
                .unwrap_or_else(|| {
                    panic!(
                        "Failed to map instrument configuration for {:?} from among {:?}",
                        mass_analyzer, self.components_to_instrument_id,
                    )
                })
        }

        fn populate_raw_signal(&self, data: &SpectrumData) -> BinaryArrayMap {
            let mut arrays = BinaryArrayMap::default();

            if let Some(mz) = data.mz() {
                let buffer = mz.bytes();
                let mz_array = DataArray::wrap(
                    &ArrayType::MZArray,
                    BinaryDataArrayType::Float64,
                    buffer.to_vec(),
                );
                arrays.add(mz_array)
            }

            if let Some(intensity) = data.intensity() {
                let buffer = intensity.bytes();
                let intensity_array = DataArray::wrap(
                    &ArrayType::IntensityArray,
                    BinaryDataArrayType::Float32,
                    buffer.to_vec(),
                );
                arrays.add(intensity_array);
            }
            arrays
        }

        fn populate_extended_data(&self, arrays: &mut BinaryArrayMap, data: &ExtendedSpectrumData) {
            if let Some(charge) = data.charge() {
                let mut array = DataArray::from_name_type_size(
                    &ArrayType::ChargeArray,
                    BinaryDataArrayType::Int32,
                    charge.len(),
                );
                for z in charge.iter() {
                    array.push(*z as i32).unwrap();
                }
                arrays.add(array)
            }

            if let Some(baseline) = data.baseline() {
                let buffer = to_bytes(baseline.as_ref());
                let intensity_array = DataArray::wrap(
                    &ArrayType::BaselineArray,
                    BinaryDataArrayType::Float32,
                    buffer,
                );
                arrays.add(intensity_array);
            }

            if let Some(noise) = data.noise() {
                let buffer = to_bytes(noise.as_ref());
                let intensity_array = DataArray::wrap(
                    &ArrayType::SignalToNoiseArray,
                    BinaryDataArrayType::Float32,
                    buffer,
                );
                arrays.add(intensity_array);
            }
        }

        fn populate_peaks(&self, data: &SpectrumData) -> PeakSetVec<C, MZ> {
            let mut peaks = PeakSetVec::empty();
            if let (Some(mz), Some(intensity)) = (data.mz(), data.intensity()) {
                for (mz_i, intensity_i) in mz.iter().zip(intensity) {
                    let peak = C::from(CentroidPeak::new(mz_i, intensity_i, 0));
                    peaks.push(peak);
                }
            }
            peaks
        }

        pub(crate) fn get_spectrum(&mut self, index: usize) -> Option<MultiLayerSpectrum<C, D>> {
            if matches!(self.detail_level, DetailLevel::MetadataOnly) {
                self.handle.set_signal_loading(false);
            } else {
                self.handle.set_signal_loading(true);
            }
            let raw = self.handle.get(index)?;
            let view = raw.view();

            let mut spec = MultiLayerSpectrum::<C, D>::default();

            spec.description.index = view.index() as usize;
            spec.description.id = make_native_id(view.index());
            spec.description.polarity = match view.polarity() {
                Polarity::Negative => ScanPolarity::Negative,
                Polarity::Positive => ScanPolarity::Positive,
                _ => ScanPolarity::Unknown,
            };
            spec.description.ms_level = view.ms_level();
            spec.description.signal_continuity = match view.mode() {
                SpectrumMode::Centroid => SignalContinuity::Centroid,
                SpectrumMode::Profile => SignalContinuity::Profile,
                _ => SignalContinuity::Unknown,
            };

            if let Some(vprec) = view.precursor() {
                let mut prec = Precursor::default();
                self.populate_precursor(vprec, &mut prec);
                spec.description.precursor = Some(prec);
            }

            let event = spec.description.acquisition.first_scan_mut().unwrap();
            event.start_time = view.time();
            if let Some(vacq) = view.acquisition() {
                self.populate_scan_event(&vacq, event);
                if let Some(filter) = view.filter_string() {
                    let mut p: Param = param!("filter string", 1000512).into();
                    p.value = filter.into();
                    event.add_param(p);
                }
            }

            if let Some(data) = view.data() {
                let extra = self.handle.get_extended_spectrum_data(index, false);
                if spec.signal_continuity() == SignalContinuity::Centroid {
                    spec.peaks = Some(self.populate_peaks(&data));
                    if let Some(extra) = extra {
                        spec.arrays = Some(self.populate_raw_signal(&data));
                        self.populate_extended_data(spec.arrays.as_mut().unwrap(), &extra);
                    }
                } else {
                    spec.arrays = Some(self.populate_raw_signal(&data));
                    if let Some(extra) = extra {
                        self.populate_extended_data(spec.arrays.as_mut().unwrap(), &extra);
                    }
                }
            }

            macro_rules! trailer {
                ($kv:ident) => {
                    let name = format!("[Thermo Trailer Extra]{}", $kv.label);
                    if !$kv.value.is_empty() {
                        let param = Param::new_key_value(name, $kv.value.parse::<Value>().unwrap());
                        spec.params_mut().push(param);
                    }
                };
            }

            spec.update_summaries();
            let ms_level = spec.ms_level();

            let trailers = self.handle.get_raw_trailers_for(index)?;
            for kv in trailers.iter() {
                match kv.label {
                    "Micro Scan Count" => {
                        trailer!(kv);
                    }
                    "Scan Segment" => {
                        trailer!(kv);
                    }
                    "Scan Event" => {
                        trailer!(kv);
                    }
                    "Monoisotopic M/Z" if ms_level > 1 && kv.value != "0" => {
                        trailer!(kv);
                    }
                    "HCD Energy eV" if ms_level > 1 && kv.value != "0" => {
                        trailer!(kv);
                    }
                    "HCD Energy" if ms_level > 1 && kv.value != "0" => {
                        trailer!(kv);
                    }
                    _ => {}
                }
            }
            Some(spec)
        }

        pub(crate) fn read_next_spectrum(&mut self) -> Option<MultiLayerSpectrum<C, D>> {
            let i = self.index;
            if i < self.len() {
                let s = self.get_spectrum(i);
                self.index += 1;
                s
            } else {
                None
            }
        }
    }
}

#[allow(unused)]
#[cfg(feature = "doc-only")]
pub(crate) mod stub {
    use super::*;

    /**
        A Thermo Fisher RAW file reader that supports iteration and random access.
    */
    // This is a stub for documentation compilation when the dotnet runtime isn't available.
    // See the the [`sealed`](super::sealed) module for the real implementation.
    pub struct ThermoRawReaderType<
        C: CentroidLike + From<CentroidPeak> = CentroidPeak,
        D: DeconvolutedCentroidLike = DeconvolutedPeak,
    > {
        pub path: PathBuf,
        pub detail_level: DetailLevel,
        pub(crate) index: usize,
        pub(crate) spectrum_index: OffsetIndex,
        pub(crate) file_description: FileDescription,
        pub(crate) instrument_configurations: HashMap<u32, InstrumentConfiguration>,
        pub(crate) softwares: Vec<Software>,
        pub(crate) samples: Vec<Sample>,
        pub(crate) data_processings: Vec<DataProcessing>,
        pub(crate) ms_run: MassSpectrometryRun,
        _c: PhantomData<C>,
        _d: PhantomData<D>,
    }

    // The public API
    // This is a stub for documentation compilation when the dotnet runtime isn't available.
    // See the the [`sealed`](super::sealed) module for the real implementation.
    impl<C: CentroidLike + From<CentroidPeak>, D: DeconvolutedCentroidLike>
        ThermoRawReaderType<C, D>
    {
        /// Get whether or not to load extended spectrum signal information for the spectrum.
        ///
        /// The loaded data isn't incorporated into a peak list, instead access them under
        /// the binary data arrays.
        pub fn get_load_extended_spectrum_data(&self) -> bool {
            false
        }

        /// Set whether or not to load extended spectrum signal information for the spectrum.
        ///
        /// The loaded data isn't incorporated into a peak list, instead access them under
        /// the binary data arrays.
        pub fn set_load_extended_spectrum_data(&mut self, load_extended_spectrum_data: bool) {}

        /// Create a new [`ThermoRawReaderType`] from a path.
        /// This may trigger an expensive I/O operation to checksum the file
        pub fn new_with_detail_level_and_centroiding<P: Into<PathBuf>>(
            path: P,
            mut detail_level: DetailLevel,
            centroiding: bool,
        ) -> io::Result<Self> {
            let path: PathBuf = path.into();

            if matches!(detail_level, DetailLevel::Lazy) {
                log::warn!(
                    "ThermoRawReader does not support lazy loading. Using {:?}",
                    DetailLevel::Full
                );
                detail_level = DetailLevel::Full
            }

            Ok(Self {
                path,
                detail_level,
                index: 0,
                spectrum_index: OffsetIndex::new("spectrum".into()),
                file_description: FileDescription::default(),
                instrument_configurations: HashMap::default(),
                samples: Vec::new(),
                softwares: Vec::new(),
                data_processings: vec![],
                ms_run: MassSpectrometryRun::default(),
                _c: PhantomData,
                _d: PhantomData,
            })
        }

        /// Create a new [`ThermoRawReaderType`] from a path.
        /// This may trigger an expensive I/O operation to checksum the file
        pub fn new<P: Into<PathBuf>>(path: P) -> io::Result<Self> {
            Self::new_with_detail_level_and_centroiding(path, DetailLevel::Full, false)
        }

        pub fn len(&self) -> usize {
            1
        }

        pub fn is_empty(&self) -> bool {
            false
        }

        /// Get whether or not to centroid spectra on read using the vendor algorithm
        pub fn get_centroiding(&self) -> bool {
            false
        }

        /// Set whether or not to centroid spectra on read using the vendor algorithm
        pub fn set_centroiding(&mut self, value: bool) {}

        pub fn get_tic(&mut self) -> Chromatogram {
            Chromatogram::default()
        }

        pub fn get_bpc(&mut self) -> Chromatogram {
            Chromatogram::default()
        }

        pub(crate) fn get_spectrum(&mut self, index: usize) -> Option<MultiLayerSpectrum<C, D>> {
            None
        }

        pub(crate) fn read_next_spectrum(&mut self) -> Option<MultiLayerSpectrum<C, D>> {
            let i = self.index;
            if i < self.len() {
                let s = self.get_spectrum(i);
                self.index += 1;
                s
            } else {
                None
            }
        }

        /// Get a status log as a [`Chromatogram`]
        pub fn get_status_log_trace_by_name(&self, trace_id: &str) -> Option<Chromatogram> {
            None
        }

        /// Enumerate the available status logs' names
        pub fn status_log_names(&self) -> impl Iterator<Item = String> {
            Option::<String>::None.into_iter()
        }
    }
}

#[cfg(feature = "doc-only")]
pub use stub::ThermoRawReaderType;

#[cfg(not(feature = "doc-only"))]
pub use sealed::ThermoRawReaderType;

impl<C: CentroidLike + From<CentroidPeak>, D: DeconvolutedCentroidLike> Iterator
    for ThermoRawReaderType<C, D>
{
    type Item = MultiLayerSpectrum<C, D>;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_next_spectrum()
    }
}

impl<C: CentroidLike + From<CentroidPeak>, D: DeconvolutedCentroidLike>
    SpectrumSource<C, D, MultiLayerSpectrum<C, D>> for ThermoRawReaderType<C, D>
{
    fn detail_level(&self) -> &DetailLevel {
        &self.detail_level
    }

    fn set_detail_level(&mut self, detail_level: DetailLevel) {
        self.detail_level = detail_level;
    }

    fn reset(&mut self) {
        self.index = 0;
    }

    fn get_spectrum_by_id(&mut self, id: &str) -> Option<MultiLayerSpectrum<C, D>> {
        let offset = self.spectrum_index.get(id)?;
        self.get_spectrum(offset as usize)
    }

    fn get_spectrum_by_index(&mut self, index: usize) -> Option<MultiLayerSpectrum<C, D>> {
        self.get_spectrum(index)
    }

    fn get_index(&self) -> &OffsetIndex {
        &self.spectrum_index
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.spectrum_index = index;
    }
}

impl<C: CentroidLike + From<CentroidPeak>, D: DeconvolutedCentroidLike>
    RandomAccessSpectrumIterator<C, D, MultiLayerSpectrum<C, D>> for ThermoRawReaderType<C, D>
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, SpectrumAccessError> {
        if let Some(i) = self.spectrum_index.get(id) {
            self.index = i as usize;
            Ok(self)
        } else {
            Err(SpectrumAccessError::SpectrumIdNotFound(id.to_string()))
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, SpectrumAccessError> {
        if index < self.len() {
            self.index = index;
            Ok(self)
        } else {
            Err(SpectrumAccessError::SpectrumIndexNotFound(index))
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, SpectrumAccessError> {
        if let Some(i) = self._offset_of_time(time) {
            self.index = i as usize;
            Ok(self)
        } else {
            Err(SpectrumAccessError::SpectrumNotFound)
        }
    }
}

impl<C: CentroidLike + From<CentroidPeak>, D: DeconvolutedCentroidLike>
    MSDataFileMetadata for ThermoRawReaderType<C, D>
{
    fn data_processings(&self) -> &Vec<DataProcessing> {
        &self.data_processings
    }
    fn instrument_configurations(&self) -> &HashMap<u32, InstrumentConfiguration> {
        &self.instrument_configurations
    }
    fn file_description(&self) -> &FileDescription {
        &self.file_description
    }
    fn softwares(&self) -> &Vec<Software> {
        &self.softwares
    }
    fn data_processings_mut(&mut self) -> &mut Vec<DataProcessing> {
        &mut self.data_processings
    }
    fn instrument_configurations_mut(&mut self) -> &mut HashMap<u32, InstrumentConfiguration> {
        &mut self.instrument_configurations
    }
    fn file_description_mut(&mut self) -> &mut FileDescription {
        &mut self.file_description
    }
    fn softwares_mut(&mut self) -> &mut Vec<Software> {
        &mut self.softwares
    }

    fn spectrum_count_hint(&self) -> Option<u64> {
        Some(self.spectrum_index.len() as u64)
    }

    fn run_description(&self) -> Option<&MassSpectrometryRun> {
        Some(&self.ms_run)
    }

    fn samples(&self) -> &Vec<Sample> {
        &self.samples
    }

    fn samples_mut(&mut self) -> &mut Vec<Sample> {
        &mut self.samples
    }

    fn run_description_mut(&mut self) -> Option<&mut MassSpectrometryRun> {
        Some(&mut self.ms_run)
    }

    fn source_file_name(&self) -> Option<&str> {
        self.path.file_name().and_then(|s| s.to_str())
    }
}

impl<C: CentroidLike + From<CentroidPeak>, D: DeconvolutedCentroidLike> ChromatogramSource
    for ThermoRawReaderType<C, D>
{
    fn get_chromatogram_by_id(&mut self, id: &str) -> Option<Chromatogram> {
        match id {
            "TIC" => Some(self.get_tic()),
            "BPC" => Some(self.get_bpc()),
            _ => None,
        }
    }

    fn get_chromatogram_by_index(&mut self, index: usize) -> Option<Chromatogram> {
        match index {
            0 => Some(self.get_tic()),
            1 => Some(self.get_bpc()),
            _ => None,
        }
    }
}

/// A convenience alias for [`ThermoRawReaderType`] with peak types specified.
pub type ThermoRawReader = ThermoRawReaderType<CentroidPeak, DeconvolutedPeak>;

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> IntoIonMobilityFrameSource<C, D>
    for ThermoRawReaderType<C, D>
{
    type IonMobilityFrameSource<
        CF: FeatureLike<MZ, mzpeaks::IonMobility>,
        DF: FeatureLike<mzpeaks::Mass, mzpeaks::IonMobility> + KnownCharge,
    > = Generic3DIonMobilityFrameSource<C, D, Self, CF, DF>;

    fn try_into_frame_source<
        CF: FeatureLike<MZ, mzpeaks::IonMobility>,
        DF: FeatureLike<mzpeaks::Mass, mzpeaks::IonMobility> + KnownCharge,
    >(
        mut self,
    ) -> Result<Self::IonMobilityFrameSource<CF, DF>, crate::io::IntoIonMobilityFrameSourceError> {
        if let Some(state) = self.has_ion_mobility() {
            if matches!(state, crate::spectrum::HasIonMobility::Dimension) {
                Ok(Self::IonMobilityFrameSource::new(self))
            } else {
                Err(crate::io::IntoIonMobilityFrameSourceError::ConversionNotPossible)
            }
        } else {
            Err(crate::io::IntoIonMobilityFrameSourceError::NoIonMobilityFramesFound)
        }
    }
}

#[cfg(test)]
mod test {
    use std::fs;

    use super::*;
    use crate::MzMLReader;

    #[test]
    fn test_read_metadata() -> io::Result<()> {
        let reader = ThermoRawReader::open_path("./test/data/small.RAW")?;
        let sf = &reader.file_description().source_files[0];
        assert_eq!(sf.id, "RAW1");
        assert_eq!(sf.name, "small.RAW");
        assert_eq!(reader.source_file_name().unwrap(), sf.name);
        assert_eq!(
            sf.get_param_by_name("SHA-1").unwrap().value,
            "b43e9286b40e8b5dbc0dfa2e428495769ca96a96"
        );

        assert!(reader
            .file_description()
            .get_param_by_accession("MS:1000579")
            .is_some());
        assert!(reader
            .file_description()
            .get_param_by_accession("MS:1000580")
            .is_some());

        let confs = reader.instrument_configurations();
        assert_eq!(confs.len(), 2);

        let conf = confs.get(&0).unwrap();
        assert_eq!(conf.id, 0);
        assert_eq!(conf.components.len(), 3);
        assert_eq!(conf.software_reference, "thermo_xcalibur");
        assert_eq!(
            conf.get_param_by_accession("MS:1000448").unwrap().name(),
            "LTQ FT"
        );

        let sw = reader
            .softwares()
            .iter()
            .find(|s| s.id == "thermo_xcalibur")
            .unwrap();
        assert!(sw.is_acquisition());
        assert!(sw.is_data_processing());

        assert_eq!(reader.samples().len(), 1);

        assert_eq!(reader.spectrum_count_hint(), Some(reader.len() as u64));

        Ok(())
    }

    #[test]
    fn test_vs_mzml() -> io::Result<()> {
        let reader = ThermoRawReader::open_path("./test/data/small.RAW")?;
        let ref_reader = MzMLReader::open_path("./test/data/small.mzML")?;

        let spectra: Vec<_> = reader.collect();
        let ref_spectra: Vec<_> = ref_reader.collect();

        assert_eq!(spectra.len(), ref_spectra.len());

        spectra.iter().zip(ref_spectra.iter()).for_each(|(s, r)| {
            assert_eq!(s.id(), r.id());
            assert_eq!(s.index(), r.index());
            assert_eq!(s.ms_level(), r.ms_level());
            assert!(
                (s.start_time() - r.start_time()).abs() < 1e-3,
                "{} - {} = {}",
                s.start_time(),
                r.start_time(),
                s.start_time() - r.start_time()
            );
            if s.ms_level() == 2 {
                let ps = s.precursor().unwrap();
                let pr = r.precursor().unwrap();
                assert!(
                    (ps.ion().mz() - pr.ion().mz()).abs() < 1e-3,
                    "{} - {} = {}",
                    ps.ion().mz(),
                    pr.ion().mz(),
                    ps.ion().mz() - pr.ion().mz()
                );
                assert_eq!(ps.ion().charge, pr.ion().charge);
            }
        });

        let s1 = &spectra[0];
        let r1 = &ref_spectra[0];

        let as1 = s1.raw_arrays().unwrap();
        let ar1 = r1.raw_arrays().unwrap();

        assert_eq!(as1.mzs().unwrap().len(), ar1.mzs().unwrap().len());
        as1.mzs()
            .unwrap()
            .iter()
            .enumerate()
            .zip(ar1.mzs().unwrap().iter())
            .for_each(|((i, s), r)| {
                assert!((s - r).abs() < 1e-3, "[{i}]{s} - {r} = {}", s - r);
            });
        assert_eq!(as1.intensities().unwrap().len(), as1.mzs().unwrap().len());
        as1.intensities()
            .unwrap()
            .iter()
            .enumerate()
            .zip(ar1.intensities().unwrap().iter())
            .for_each(|((i, s), r)| {
                assert!((s - r).abs() < 1e-3, "[{i}]{s} - {r} = {}", s - r);
            });
        Ok(())
    }

    #[test]
    fn test_read_spectra() -> io::Result<()> {
        let mut reader = ThermoRawReader::new_with_detail_level_and_centroiding(
            "./test/data/small.RAW",
            DetailLevel::Lazy,
            false,
        )?;
        assert_eq!(reader.len(), 48);

        let groups: Vec<_> = reader.groups().collect();
        assert_eq!(groups.len(), 14);

        let spec = reader.get_spectrum_by_index(0).unwrap();
        let first_sid = spec.id();
        assert_eq!(spec.peaks().len(), 19913);
        assert!(
            (spec.peaks().tic() - 71263170.0).abs() < 1.0,
            "TIC = {}, diff {}",
            spec.peaks().tic(),
            spec.peaks().tic() - 71263170.0
        );

        let bp = spec.peaks().base_peak();
        assert!(
            (bp.mz - 810.41547).abs() < 1e-3,
            "Base Peak m/z = {}, Diff = {}",
            bp.mz,
            bp.mz - 810.41547
        );

        assert_eq!(spec.ms_level(), 1);
        assert_eq!(spec.signal_continuity(), SignalContinuity::Profile);
        assert_eq!(spec.polarity(), ScanPolarity::Positive);

        let event = spec.acquisition().first_scan().unwrap();
        assert_eq!(1, event.instrument_configuration_id);

        assert!((event.injection_time - 68.227486).abs() < 1e-3);
        assert!((event.start_time - 0.004935).abs() < 1e-3);

        for _ in reader.instrument_configurations().iter() {}

        assert!(!reader.get_centroiding());
        reader.set_centroiding(true);
        assert!(reader.get_centroiding());
        let spec_centr = reader.get_spectrum_by_index(spec.index()).unwrap();
        assert_eq!(spec_centr.signal_continuity(), SignalContinuity::Centroid);
        assert!(spec_centr.peaks.is_some());

        let scan = reader.start_from_index(20).unwrap().next().unwrap();
        assert_eq!(scan.index(), 20);
        let time = scan.start_time();
        let scan = reader
            .start_from_id(first_sid)
            .ok()
            .and_then(|it| it.next())
            .unwrap();
        assert_eq!(scan.index(), 0);
        let scan = reader
            .start_from_time(time)
            .ok()
            .and_then(|it| it.next())
            .unwrap();
        assert_eq!(scan.index(), 20);

        reader.reset();
        Ok(())
    }

    #[test]
    fn test_sniff() -> io::Result<()> {
        let mut handle = fs::File::open("./test/data/small.RAW")?;
        let mut buf: [u8; 128] = [0; 128];
        handle.read_exact(&mut buf)?;
        assert!(is_thermo_raw_prefix(&buf));
        buf.reverse();
        assert!(!is_thermo_raw_prefix(&buf));
        Ok(())
    }

    #[test]
    fn test_read_chromatograms() -> io::Result<()> {
        let mut reader = ThermoRawReader::open_path("./test/data/small.RAW")?;
        let tic = reader.get_chromatogram_by_id("TIC").unwrap();
        let exp_n = reader.len();
        let obs_n = tic.time()?.len();
        assert_eq!(obs_n, exp_n);

        let tic = reader.get_chromatogram_by_id("BPC").unwrap();
        let exp_n = reader.len();
        let obs_n = tic.time()?.len();
        assert_eq!(obs_n, exp_n);

        let tic = reader.get_chromatogram_by_index(0).unwrap();
        let exp_n = reader.len();
        let obs_n = tic.time()?.len();
        assert_eq!(obs_n, exp_n);

        let tic = reader.get_chromatogram_by_index(1).unwrap();
        let exp_n = reader.len();
        let obs_n = tic.time()?.len();
        assert_eq!(obs_n, exp_n);

        assert_eq!(reader.iter_chromatograms().count(), 2);
        Ok(())
    }

    #[test]
    fn test_extended_data_model() -> io::Result<()> {
        let mut reader = ThermoRawReader::open_path("./test/data/small.RAW")?;
        reader.set_centroiding(true);
        reader.set_load_extended_spectrum_data(true);

        let spec = reader.get_spectrum_by_index(0).unwrap();
        let arrays = spec.arrays.as_ref().unwrap();
        arrays.iter().find(|(k, _)| **k == ArrayType::ChargeArray);
        Ok(())
    }
}
