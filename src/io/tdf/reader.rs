use std::{
    collections::HashMap, io, iter::FusedIterator, marker::PhantomData, ops::Range, path::Path,
    sync::Arc,
};

use chrono::DateTime;

#[allow(unused)]
use crate::io::checksum_file;

use crate::{
    curie,
    io::{DetailLevel, IntoIonMobilityFrameSource, IonMobilityFrameAccessError, OffsetIndex},
    meta::{
        Component, ComponentType, DataProcessing, DetectorTypeTerm,
        DissociationMethodTerm::CollisionInducedDissociation, FileDescription,
        InstrumentConfiguration, MassAnalyzerTerm, MassSpectrometerFileFormatTerm,
        MassSpectrometryRun, NativeSpectrumIdentifierFormatTerm, Sample, Software, SoftwareTerm,
        SourceFile,
    },
    mzpeaks::{
        feature::{ChargedFeature, Feature},
        CentroidPeak, DeconvolutedPeak, IonMobility, Mass, MZ,
    },
    params::{ControlledVocabulary, Unit, Value},
    prelude::*,
    spectrum::{
        Activation, ArrayType, BinaryArrayMap, BinaryDataArrayType, Chromatogram,
        ChromatogramDescription, ChromatogramType, DataArray, IonMobilityFrameDescription,
        IsolationWindow, IsolationWindowState, MultiLayerIonMobilityFrame, MultiLayerSpectrum,
        Precursor, ScanCombination, ScanEvent, ScanWindow, SelectedIon, SignalContinuity,
    },
    Param,
};
use identity_hash::BuildIdentityHasher;
use rusqlite::Error;

use timsrust::{
    converters::ConvertableDomain,
    readers::{FrameReader, FrameReaderError, MetadataReader},
    Metadata, TimsRustError,
};

pub use super::arrays::FrameToArraysMapper;
use super::{
    arrays::consolidate_peaks,
    constants::{InstrumentSource, MsMsType},
    sql::{
        ChromatographyData, FromSQL, PasefPrecursor, RawTDFSQLReader, SQLDIAFrameMsMsWindow,
        SQLFrame, SQLPasefFrameMsMs, SQLPrecursor, TDFMSnFacet,
    },
};

const PEAK_MERGE_TOLERANCE: Tolerance = Tolerance::PPM(10.0);

fn inverse_reduce_ion_mobility_param(value: f64) -> Param {
    ControlledVocabulary::MS
        .param_val(1002815, "inverse reduced ion mobility", value)
        .with_unit_t(&Unit::VoltSecondPerSquareCentimeter)
}

#[derive(Debug, Clone)]
pub struct IndexExtry {
    pub frame: Arc<SQLFrame>,
    pub index: usize,
    pub parent_index: Option<usize>,
    pub tdf_facet: TDFMSnFacet,
}

impl IndexExtry {
    pub fn new(
        frame: Arc<SQLFrame>,
        index: usize,
        parent_index: Option<usize>,
        tdf_facet: TDFMSnFacet,
    ) -> Self {
        Self {
            frame,
            index,
            parent_index,
            tdf_facet,
        }
    }

    pub fn scan_range(&self) -> (usize, usize) {
        match self.ms_type() {
            MsMsType::MS1 => (0, self.frame.num_scans),
            MsMsType::DDAPASEF => {
                if let Some(pasef) = self.tdf_facet.pasef_msms() {
                    (pasef.scan_start, pasef.scan_end)
                } else {
                    (0, self.frame.num_scans)
                }
            }
            MsMsType::DIAPASEF => {
                if let Some(pasef) = self.tdf_facet.dia_window() {
                    (pasef.scan_start, pasef.scan_end)
                } else {
                    (0, self.frame.num_scans)
                }
            }
            MsMsType::MRM | MsMsType::PRMPASEF | MsMsType::Unknown => todo!(),
        }
    }

    pub fn format_native_id(&self) -> String {
        let (start_scan, end_scan) = self.scan_range();
        let start_scan = start_scan + 1;
        format!(
            "merged={} frame={} startScan={start_scan} endScan={end_scan}",
            self.index, self.frame.id,
        )
    }

    pub fn ms_type(&self) -> MsMsType {
        MsMsType::from(self.frame.msms_type)
    }

    pub fn ms_level(&self) -> u8 {
        self.ms_type().ms_level()
    }

    pub fn pasef_msms(&self) -> Option<&SQLPasefFrameMsMs> {
        self.tdf_facet.pasef_msms()
    }

    pub fn precursor(&self) -> Option<&SQLPrecursor> {
        self.tdf_facet.precursor()
    }

    pub fn dia_window(&self) -> Option<&SQLDIAFrameMsMsWindow> {
        self.tdf_facet.dia_window()
    }
}

/// An ion mobility frame reader for Bruker TDF file format. It supports the same
/// kind of PASEF frame traversal that is supported by [ProteoWizard](https://github.com/ProteoWizard/pwiz).
///
/// It implements the full range of [`IonMobilityFrameSource`]-derived traits. To view
/// these frames as flat peak lists, this type can be wrapped in a [`TDFSpectrumReaderType`]
/// using [`TDFFrameReaderType::into_spectrum_reader`].
#[derive(Debug)]
pub struct TDFFrameReaderType<
    C: FeatureLike<MZ, IonMobility> = Feature<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge = ChargedFeature<Mass, IonMobility>,
> {
    metadata: timsrust::Metadata,
    frame_reader: timsrust::readers::FrameReader,
    tdf_reader: RawTDFSQLReader,
    entry_index: Vec<IndexExtry>,
    index: usize,
    offset_index: OffsetIndex,
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
    pub detail_level: DetailLevel,

    // SpectrumList attributes
    pub run: MassSpectrometryRun,
    _c: PhantomData<C>,
    _d: PhantomData<D>,
}

type DIAFrameWindowMap = HashMap<u32, Vec<Arc<SQLDIAFrameMsMsWindow>>, BuildIdentityHasher<u32>>;

impl<C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge>
    TDFFrameReaderType<C, D>
{
    /// Construct a new reader from the specified file system path.
    ///
    /// # Errors
    /// This may fail if any of the component files fails to match the expected schema
    /// or layout.
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, timsrust::TimsRustError> {
        Self::new_with_detail_level(path, DetailLevel::Full)
    }

    /// Construct a new reader from the specified file system path with the specified
    /// [`DetailLevel`].
    ///
    /// # Errors
    /// This may fail if any of the component files fails to match the expected schema
    /// or layout.
    pub fn new_with_detail_level<P: AsRef<Path>>(
        path: P,
        detail_level: DetailLevel,
    ) -> Result<Self, timsrust::TimsRustError> {
        let path = path.as_ref();
        let tdf_path = path.join("analysis.tdf");
        if !tdf_path.exists() {
            return Err(timsrust::TimsRustError::FrameReaderError(
                FrameReaderError::FileNotFound(tdf_path.display().to_string()),
            ));
        }

        let metadata = MetadataReader::new(&tdf_path)?;
        let frame_reader = FrameReader::new(path)?;
        let tdf_reader = RawTDFSQLReader::new(&tdf_path)
            .map_err(|e| TimsRustError::FrameReaderError(FrameReaderError::SqlError(e.into())))?;

        let mut this = Self {
            metadata,
            frame_reader,
            tdf_reader,
            entry_index: Vec::new(),
            index: 0,
            offset_index: OffsetIndex::new("spectrum".into()),

            instrument_configurations: HashMap::default(),
            file_description: FileDescription::default(),
            softwares: Vec::new(),
            data_processings: Vec::new(),
            samples: Vec::new(),

            run: MassSpectrometryRun::default(),
            detail_level,

            _c: PhantomData,
            _d: PhantomData,
        };

        this.build_index()
            .map_err(|e| TimsRustError::FrameReaderError(FrameReaderError::SqlError(e.into())))?;
        this.build_metadata().unwrap();

        Ok(this)
    }

    /// Consume this reader, wrapping it in a [`TDFSpectrumReaderType`] with the default
    /// peak merging tolerance.
    pub fn into_spectrum_reader<CP: CentroidLike, DP: DeconvolutedCentroidLike>(
        self,
    ) -> TDFSpectrumReaderType<C, D, CP, DP> {
        self.into_spectrum_reader_with_peak_merging_tolerance(PEAK_MERGE_TOLERANCE)
    }

    /// Consume this reader, wrapping it in a [`TDFSpectrumReaderType`] with the
    /// specified peak merging error tolerance.
    pub fn into_spectrum_reader_with_peak_merging_tolerance<
        CP: CentroidLike,
        DP: DeconvolutedCentroidLike,
    >(
        self,
        peak_merging_tolerance: Tolerance,
    ) -> TDFSpectrumReaderType<C, D, CP, DP> {
        TDFSpectrumReaderType {
            frame_reader: self,
            peak_merging_tolerance,
            do_consolidate_peaks: true,
            _cp: PhantomData,
            _dp: PhantomData,
        }
    }

    fn guess_number_of_entries(&self) -> Result<(usize, Vec<MsMsType>), Error> {
        let conn = self.tdf_reader.connection();

        let mut stmt = conn
            .prepare("SELECT MsMsType, COUNT(1) FROM Frames GROUP BY MsMsType")
            .unwrap_or_else(|e| panic!("Failed to count scan types: {}", e));

        let counts_per_msms_type: Vec<(u8, usize)> = stmt
            .query([])?
            .mapped(|row| -> Result<(u8, usize), Error> { Ok((row.get(0)?, row.get(1)?)) })
            .flatten()
            .collect();

        let mut entry_count_guess = 0;
        let mut frame_types = Vec::new();
        for (msms_type, count_of) in counts_per_msms_type.iter().copied() {
            let ft = MsMsType::from(msms_type);
            frame_types.push(ft);
            let entry_count = match ft {
                MsMsType::MS1 | MsMsType::MRM | MsMsType::Unknown => count_of,
                MsMsType::DDAPASEF => conn.query_row(
                    "SELECT COUNT(1) FROM PasefFrameMsMsInfo pasef JOIN Precursors prec ON pasef.Precursor = prec.Id",
                    [], |row| {
                    row.get::<usize, usize>(0)
                }).unwrap_or_else(|e| panic!("Failed to count precursors: {}", e)),
                MsMsType::DIAPASEF => {
                    let mut q = conn.prepare(r#"SELECT frames_of.WindowGroup, FrameCount * GroupSize FROM
                    (SELECT info.WindowGroup, COUNT(1) AS FrameCount FROM Frames f JOIN DiaFrameMsMsInfo info ON f.Id = info.Frame GROUP BY info.WindowGroup) frames_of
                    JOIN
                    (SELECT WindowGroup AS WindowGroup, COUNT(1) AS GroupSize FROM DiaFrameMsMsWindows GROUP BY WindowGroup) as window_group_size
                    ON frames_of.WindowGroup = window_group_size.WindowGroup GROUP BY frames_of.WindowGroup"#)?;

                    let total = q.query_map([], |row| {
                        let count: usize = row.get(1)?;
                        Ok(count)
                    })?.flatten().sum::<usize>();
                    #[allow(clippy::let_and_return)]
                    total
                },
                MsMsType::PRMPASEF => todo!(),
            };
            entry_count_guess += entry_count;
        }

        Ok((entry_count_guess, frame_types))
    }

    fn build_window_index(&self) -> Result<DIAFrameWindowMap, Error> {
        let conn = self.tdf_reader.connection();
        let dia_windows = SQLDIAFrameMsMsWindow::read_from(&conn, []).unwrap_or_default();
        let mut window_groups: HashMap<u32, Vec<Arc<SQLDIAFrameMsMsWindow>>, _> =
            HashMap::default();
        for window in dia_windows {
            window_groups
                .entry(window.window_group)
                .or_default()
                .push(Arc::new(window));
        }
        Ok(window_groups)
    }

    fn build_index(&mut self) -> Result<(), Error> {
        let pasef_window_groups = self.build_window_index()?;

        let (entry_guess, _frame_types) = self
            .guess_number_of_entries()
            .unwrap_or_else(|e| panic!("Failed to estimate index size: {}", e));
        let mut index_entries = Vec::with_capacity(entry_guess);

        let conn = self.tdf_reader.connection();

        let precursors: Vec<_> = SQLPrecursor::read_from(&conn, [])
            .unwrap_or_default()
            .into_iter()
            .map(Arc::new)
            .collect();

        let pasef_msms_info: HashMap<u32, Vec<Arc<SQLPasefFrameMsMs>>, BuildIdentityHasher<u32>> =
            SQLPasefFrameMsMs::read_from(&conn, [])
                .unwrap_or_default()
                .into_iter()
                .map(Arc::new)
                .fold(HashMap::default(), |mut index, p| {
                    index.entry(p.frame as u32).or_default().push(p);
                    index
                });

        let q = format!("{} ORDER BY Id", SQLFrame::get_sql());
        let frames = self
            .tdf_reader
            .query::<SQLFrame>(&q, [])
            .unwrap_or_default();
        let mut frame_index = Vec::with_capacity(frames.len());
        let mut parent_index: HashMap<usize, usize> = HashMap::new();
        let mut last_parent = 0usize;
        for frame in frames.into_iter().map(Arc::new) {
            frame_index.push(frame.clone());
            if frame.num_peaks == 0 {
                continue;
            }
            let msms_type = MsMsType::from(frame.msms_type);
            match msms_type {
                MsMsType::MS1 => {
                    let i = index_entries.len();
                    last_parent = i;
                    parent_index.insert(frame.id, i);
                    index_entries.push(IndexExtry::new(frame, i, None, TDFMSnFacet::None));
                }
                MsMsType::DDAPASEF => {
                    let pasef_infos = pasef_msms_info.get(&(frame.id as u32)).unwrap();
                    for pasef_info in pasef_infos {
                        let precursors = precursors
                            .get(pasef_info.precursor.saturating_sub(1))
                            .cloned();
                        if let Some(precursor) = precursors {
                            let parent = parent_index.get(&precursor.precursor_frame).copied();
                            let entry = IndexExtry::new(
                                frame.clone(),
                                index_entries.len(),
                                parent,
                                PasefPrecursor::new(precursor, pasef_info.clone()).into(),
                            );
                            index_entries.push(entry);
                        }
                    }
                }
                MsMsType::DIAPASEF => {
                    let group_id = self.tdf_reader.dia_window_group_for(frame.id)?;
                    for window in pasef_window_groups.get(&group_id).into_iter().flatten() {
                        let entry = IndexExtry::new(
                            frame.clone(),
                            index_entries.len(),
                            Some(last_parent),
                            window.clone().into(),
                        );
                        index_entries.push(entry);
                    }
                }
                MsMsType::PRMPASEF => todo!(),
                MsMsType::MRM => todo!(),
                MsMsType::Unknown => todo!(),
            }
        }
        self.entry_index = index_entries;
        let mut offset_index = OffsetIndex::new("spectrum".into());
        self.entry_index.iter().enumerate().for_each(|(i, e)| {
            offset_index.insert(e.format_native_id().into_boxed_str(), i as u64);
        });
        self.offset_index = offset_index;
        Ok(())
    }

    fn build_chromatogram(&self, chromatogram_type: ChromatogramType) -> Option<Chromatogram> {
        let mut descr = ChromatogramDescription::default();
        match chromatogram_type {
            ChromatogramType::TotalIonCurrentChromatogram => {
                descr.chromatogram_type = ChromatogramType::TotalIonCurrentChromatogram;
                descr.id = "TIC".into();
                descr.index = 0;
            }
            ChromatogramType::BasePeakChromatogram => {
                descr.chromatogram_type = ChromatogramType::BasePeakChromatogram;
                descr.id = "BPC".into();
                descr.index = 1;
            }

            _ => return None,
        }

        let n = self.entry_index.len();
        let mut time_array: Vec<u8> =
            Vec::with_capacity(n * BinaryDataArrayType::Float64.size_of());
        let mut intensity_array: Vec<u8> =
            Vec::with_capacity(n * BinaryDataArrayType::Float32.size_of());

        for entry in self.entry_index.iter() {
            time_array.extend_from_slice(&(entry.frame.time / 60.0).to_le_bytes());
            match chromatogram_type {
                ChromatogramType::TotalIonCurrentChromatogram => {
                    intensity_array
                        .extend_from_slice(&entry.frame.summed_intensities.to_le_bytes());
                }
                ChromatogramType::BasePeakChromatogram => {
                    intensity_array.extend_from_slice(&entry.frame.max_intensity.to_le_bytes());
                }
                _ => {
                    unimplemented!()
                }
            }
        }

        let mut arrays = BinaryArrayMap::default();
        let mut time_array = DataArray::wrap(
            &ArrayType::TimeArray,
            BinaryDataArrayType::Float64,
            time_array,
        );
        time_array.unit = Unit::Minute;
        arrays.add(time_array);

        let intensity_array = DataArray::wrap(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            intensity_array,
        );
        arrays.add(intensity_array);

        Some(Chromatogram::new(descr, arrays))
    }

    pub(crate) fn tic(&self) -> Chromatogram {
        self.build_chromatogram(ChromatogramType::TotalIonCurrentChromatogram)
            .unwrap()
    }

    pub(crate) fn bpc(&self) -> Chromatogram {
        self.build_chromatogram(ChromatogramType::BasePeakChromatogram)
            .unwrap()
    }

    fn index_by_time(&self, time: f64) -> Option<&IndexExtry> {
        let n = self.entry_index.len();
        let mut lo: usize = 0;
        let mut hi: usize = n;

        let mut best_error: f64 = f64::INFINITY;
        let mut best_match: Option<&IndexExtry> = None;

        if lo == hi {
            return None;
        }
        while hi != lo {
            let mid = (hi + lo) / 2;
            let scan = self.entry_index.get(mid).unwrap();
            let scan_time = scan.frame.time / 60.0;
            let err = (scan_time - time).abs();

            if err < best_error {
                best_error = err;
                best_match = Some(scan);
            }
            if hi.saturating_sub(1) == lo {
                return best_match;
            } else if scan_time > time {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        best_match
    }

    fn frame_to_spectrum<CP: CentroidLike + From<CentroidPeak>, DP: DeconvolutedCentroidLike>(
        &self,
        frame: MultiLayerIonMobilityFrame<C, D>,
        error_tolerance: Tolerance,
        do_consolidate: bool,
    ) -> MultiLayerSpectrum<CP, DP> {
        let (feature_d, mut descr) = frame.into_features_and_parts();
        match feature_d {
            crate::spectrum::FeatureDataLevel::Missing => MultiLayerSpectrum {
                description: descr.into(),
                arrays: None,
                peaks: None,
                deconvoluted_peaks: None,
            },
            crate::spectrum::FeatureDataLevel::RawData(arrays) => {
                let peaks = if do_consolidate {
                    consolidate_peaks(
                        &arrays,
                        &(0..arrays.ion_mobility_dimension.len() as u32),
                        &self.metadata,
                        error_tolerance,
                    )
                    .ok()
                } else {
                    None
                };
                descr.signal_continuity = SignalContinuity::Centroid;
                let arrays = arrays.unstack().unwrap();
                let spec = MultiLayerSpectrum {
                    description: descr.into(),
                    arrays: Some(arrays),
                    peaks,
                    deconvoluted_peaks: None,
                };
                spec
            }
            _ => panic!("Failed to extract array data"),
        }
    }

    pub(crate) fn get(
        &self,
        index: usize,
    ) -> Result<Option<MultiLayerIonMobilityFrame<C, D>>, TimsRustError> {
        if let Some(entry) = self.entry_index.get(index) {
            // `timsrust` uses base-zero indexing, but frame IDs start at 1
            let frame = self
                .frame_reader
                .get(entry.frame.id.saturating_sub(1))
                .inspect_err(|e| {
                    log::error!("Failed to read frame {index}: {e}");
                })?;

            let mut descr = frame_to_description(&self.metadata, entry, None);

            if let Some(parent_entry) = entry.parent_index.and_then(|i| self.entry_index.get(i)) {
                descr.precursor = index_to_precursor(entry, &self.metadata, parent_entry);
            }

            let arrays = if !matches!(self.detail_level, DetailLevel::MetadataOnly) {
                if let Some(pasef) = entry.pasef_msms() {
                    let arrays = FrameToArraysMapper::new(&frame, &self.metadata)
                        .process_3d_slice(pasef.scan_start..pasef.scan_end);
                    Some(arrays)
                } else if let Some(dia_pasef) = entry.dia_window() {
                    let arrays = FrameToArraysMapper::new(&frame, &self.metadata)
                        .process_3d_slice(dia_pasef.scan_start..dia_pasef.scan_end);
                    Some(arrays)
                } else {
                    Some(FrameToArraysMapper::new(&frame, &self.metadata).process_3d_slice(..))
                }
            } else {
                None
            };

            let frame = MultiLayerIonMobilityFrame::new(arrays, None, None, descr);
            Ok(Some(frame))
        } else {
            Ok(None)
        }
    }

    pub fn get_trace_reader(&self) -> Result<ChromatographyData, TimsRustError> {
        let path = self
            .metadata
            .path
            .parent()
            .expect(".tdf file did not have an enclosing directory")
            .join(super::sql::ChromatographyData::FILE_NAME);
        let handle = super::sql::ChromatographyData::new(&path).map_err(|e| {
            TimsRustError::MetadataReaderError(timsrust::readers::MetadataReaderError::SqlError(
                e.into(),
            ))
        })?;
        Ok(handle)
    }
}

// Metadata construction routine
impl<C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge>
    TDFFrameReaderType<C, D>
{
    fn build_file_description(&self) -> io::Result<FileDescription> {
        let mut descr = FileDescription::default();

        let (ms1_count, msn_count) = self.entry_index.iter().fold((0, 0), |(ms1, msn), entry| {
            if entry.ms_level() > 1 {
                (ms1, msn + 1)
            } else {
                (ms1 + 1, msn)
            }
        });

        if ms1_count > 0 {
            descr.contents.add_param(
                Param::builder()
                    .name("MS1 spectrum")
                    .curie(curie!(MS:1000579))
                    .build(),
            );
        }

        if msn_count > 0 {
            descr.contents.add_param(
                Param::builder()
                    .name("MSn spectrum")
                    .curie(curie!(MS:1000580))
                    .build(),
            );
        }

        let mut sf = SourceFile::from_path(&self.metadata.path)?;
        sf.file_format = Some(MassSpectrometerFileFormatTerm::BrukerTDF.into());
        sf.id_format = Some(NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat.into());
        #[cfg(not(debug_assertions))]
        sf.add_param(ControlledVocabulary::MS.param_val(
            1000569u32,
            "SHA-1",
            checksum_file(&self.metadata.path)?,
        ));
        sf.id = self
            .metadata
            .path
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| "analysis_tdf".into());
        descr.source_files.push(sf);

        let tdf_bin = self.metadata.path.with_extension("tdf_bin");
        if tdf_bin.exists() {
            let mut sf = SourceFile::from_path(&tdf_bin)?;
            sf.file_format = Some(MassSpectrometerFileFormatTerm::BrukerTDF.into());
            sf.id_format = Some(NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat.into());
            #[cfg(not(debug_assertions))]
            sf.add_param(ControlledVocabulary::MS.param_val(
                1000569u32,
                "SHA-1",
                checksum_file(&tdf_bin)?,
            ));
            sf.id = tdf_bin
                .file_name()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| "analysis_tdf_bin".into());
            descr.source_files.push(sf);
        }

        Ok(descr)
    }

    fn build_software_from_metadata(&self, metadata: &HashMap<String, Value>) -> Vec<Software> {
        let mut sw_entries = Vec::new();

        let sdk = Software::new(
            "TIMS_SDK".into(),
            "*".into(),
            vec![
                SoftwareTerm::BrukerSoftware.into(),
                Param::new_key_value("software name", "timsrust"),
            ],
        );
        sw_entries.push(sdk);

        let acq_sw_name = metadata
            .get("AcquisitionSoftware")
            .map(|s| {
                let s = s.as_str();

                if s.contains("oTOFcontrol") {
                    SoftwareTerm::MicrOTOFcontrol
                } else if s.contains("HCT") {
                    SoftwareTerm::HCTcontrol
                } else {
                    SoftwareTerm::MicrOTOFcontrol
                }
                .to_param()
            })
            .unwrap();

        let acq_sw_vers = metadata
            .get("AcquisitionSoftwareVersion")
            .map(|v| v.to_string())
            .unwrap_or_else(|| "*".to_string());

        sw_entries.push(Software::new(
            "ACQ_SW".to_string(),
            acq_sw_vers,
            vec![acq_sw_name.into()],
        ));

        sw_entries
    }

    fn build_instrument_configuration(
        &self,
        metadata: &HashMap<String, Value>,
    ) -> InstrumentConfiguration {
        let components = {
            let source = metadata
                .get("InstrumentSourceType")
                .into_iter()
                .flat_map(|s| s.to_i32())
                .map(|s| InstrumentSource::from(s as u8).to_component())
                .next()
                .unwrap_or_default();

            let mut components = Vec::new();
            components.push(source);

            let mut analyzer = Component {
                order: 2,
                component_type: ComponentType::Analyzer,
                ..Default::default()
            };
            analyzer.add_param(MassAnalyzerTerm::Quadrupole.into());
            components.push(analyzer);

            let mut analyzer = Component {
                component_type: ComponentType::Analyzer,
                order: 3,
                ..Default::default()
            };
            analyzer.add_param(MassAnalyzerTerm::TimeOfFlight.into());
            components.push(analyzer);

            let mut detector = Component {
                component_type: ComponentType::Detector,
                order: 4,
                ..Default::default()
            };
            detector.add_param(DetectorTypeTerm::MicrochannelPlateDetector.into());
            components.push(detector);

            detector = Component {
                component_type: ComponentType::Detector,
                order: 5,
                ..Default::default()
            };
            detector.add_param(DetectorTypeTerm::Photomultiplier.into());
            components.push(detector);

            components
        };
        let mut config = InstrumentConfiguration {
            components,
            id: 0,
            ..Default::default()
        };
        if let Some(serial) = metadata.get("InstrumentSerialNumber").cloned() {
            config.add_param(
                Param::builder()
                    .curie(curie!(MS:1000529))
                    .name("instrument serial number")
                    .value(serial)
                    .build(),
            );
        }
        if let Some(model) = metadata.get("InstrumentName").cloned() {
            let (name, curie) = match model.as_str().as_ref().trim() {
                "timsTOF Pro" => ("timsTOF Pro", curie!(MS:1003005)),
                "timsTOF fleX" => ("timsTOF fleX", curie!(MS:1003124)),
                "timsTOF Pro 2" => ("timsTOF Pro 2", curie!(MS:1003230)),
                "timsTOF SCP" => ("timsTOF SCP", curie!(MS:1003231)),
                "timsTOF Ultra" => ("timsTOF Ultra", curie!(MS:1003383)),
                "timsTOF fleX MALDI-2" => ("timsTOF fleX MALDI-2", curie!(MS:1003397)),
                "timsTOF HT" => ("timsTOF HT", curie!(MS:1003404)),
                "timsTOF Ultra 2" => ("timsTOF Ultra 2", curie!(MS:1003412)),
                "timsTOF" => ("timsTOF", curie!(MS:1003229)),
                _ => ("Bruker Daltonics timsTOF series", curie!(MS:1003123)),
            };
            config.add_param(ControlledVocabulary::MS.param(curie.accession, name));
        }
        config.software_reference = "ACQ_SW".into();
        config
    }

    fn build_ms_run(
        &self,
        metadata: &HashMap<String, Value>,
        file_description: &FileDescription,
    ) -> MassSpectrometryRun {
        let sf = file_description.source_files.first().unwrap();
        let run_id = self
            .metadata
            .path
            .parent()
            .map(|s| s.as_os_str().to_string_lossy().to_string());

        let start_time = metadata.get("AcquisitionDateTime").and_then(|s| {
            let s = s.as_str();
            DateTime::parse_from_rfc3339(&s).ok()
        });

        MassSpectrometryRun::new(run_id, None, Some(1), Some(sf.id.clone()), start_time)
    }

    fn build_sample(&self, metadata: &HashMap<String, Value>) -> Option<Sample> {
        metadata.get("SampleName").map(|name| {
            let name = name.to_string();
            let mut params = Vec::new();
            if let Some(analysis_id) = metadata.get("AnalysisId").cloned() {
                params.push(Param::new_key_value("TDF:AnalysisId", analysis_id))
            }
            Sample::new("SAMPLE_1".into(), Some(name), params)
        })
    }

    fn build_metadata(&mut self) -> io::Result<()> {
        let metadata = self.tdf_reader.metadata().unwrap();
        let file_description = self.build_file_description()?;

        self.softwares = self.build_software_from_metadata(&metadata);
        self.instrument_configurations
            .insert(0, self.build_instrument_configuration(&metadata));
        self.samples.extend(self.build_sample(&metadata));
        self.run = self.build_ms_run(&metadata, &file_description);
        self.file_description = file_description;
        Ok(())
    }
}

impl<C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge>
    MSDataFileMetadata for TDFFrameReaderType<C, D>
{
    crate::impl_metadata_trait!();

    fn run_description(&self) -> Option<&MassSpectrometryRun> {
        Some(&self.run)
    }

    fn run_description_mut(&mut self) -> Option<&mut MassSpectrometryRun> {
        Some(&mut self.run)
    }

    fn spectrum_count_hint(&self) -> Option<u64> {
        Some(self.entry_index.len() as u64)
    }

    fn source_file_name(&self) -> Option<&str> {
        self.metadata.path.to_str()
    }
}

impl<C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge>
    ChromatogramSource for TDFFrameReaderType<C, D>
{
    fn get_chromatogram_by_id(&mut self, id: &str) -> Option<Chromatogram> {
        match id {
            "TIC" => Some(self.tic()),
            "BPC" => Some(self.bpc()),
            _ => None,
        }
    }

    fn get_chromatogram_by_index(&mut self, index: usize) -> Option<Chromatogram> {
        match index {
            0 => Some(self.tic()),
            1 => Some(self.bpc()),
            _ => None,
        }
    }
}

impl<C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge> Iterator
    for TDFFrameReaderType<C, D>
{
    type Item = MultiLayerIonMobilityFrame<C, D>;

    fn next(&mut self) -> Option<Self::Item> {
        let frame = self.get(self.index);
        self.index += 1;
        frame.ok().unwrap_or_default()
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.index += n;
        self.next()
    }
}

impl<C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge>
    IonMobilityFrameSource<C, D, MultiLayerIonMobilityFrame<C, D>> for TDFFrameReaderType<C, D>
{
    fn reset(&mut self) {
        self.index = 0;
    }

    fn get_frame_by_id(&mut self, id: &str) -> Option<MultiLayerIonMobilityFrame<C, D>> {
        let offset = self.offset_index.get(id)?;
        self.get(offset as usize).ok().unwrap_or_default()
    }

    fn get_frame_by_index(&mut self, index: usize) -> Option<MultiLayerIonMobilityFrame<C, D>> {
        self.get(index).ok().unwrap_or_default()
    }

    fn get_frame_by_time(&mut self, time: f64) -> Option<MultiLayerIonMobilityFrame<C, D>> {
        self.index_by_time(time)
            .and_then(|e| self.get(e.index).ok().unwrap_or_default())
    }

    fn get_index(&self) -> &OffsetIndex {
        &self.offset_index
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.offset_index = index
    }

    fn detail_level(&self) -> &DetailLevel {
        &self.detail_level
    }

    fn set_detail_level(&mut self, detail_level: DetailLevel) {
        self.detail_level = detail_level;
    }
}

impl<C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge> FusedIterator
    for TDFFrameReaderType<C, D>
{
}

impl<C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge>
    RandomAccessIonMobilityFrameIterator<C, D, MultiLayerIonMobilityFrame<C, D>>
    for TDFFrameReaderType<C, D>
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, IonMobilityFrameAccessError> {
        match self.offset_index.index_of(id) {
            Some(i) => {
                self.index = i;
                Ok(self)
            }
            None => Err(IonMobilityFrameAccessError::FrameIdNotFound(id.into())),
        }
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, IonMobilityFrameAccessError> {
        if index < self.offset_index.len() {
            self.index = index;
            Ok(self)
        } else {
            Err(IonMobilityFrameAccessError::FrameIndexNotFound(index))
        }
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, IonMobilityFrameAccessError> {
        match self.index_by_time(time) {
            Some(entry) => {
                self.index = entry.index;
                Ok(self)
            }
            None => Err(IonMobilityFrameAccessError::FrameNotFound),
        }
    }
}

pub type TDFFrameReader =
    TDFFrameReaderType<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>>;

/// A flat spectrum reader for Bruker TDF file format. It sums over ion mobility
/// spectra, consolidating features into peaks.
#[derive(Debug)]
pub struct TDFSpectrumReaderType<
    C: FeatureLike<MZ, IonMobility> = Feature<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge = ChargedFeature<Mass, IonMobility>,
    CP: CentroidLike = CentroidPeak,
    DP: DeconvolutedCentroidLike = DeconvolutedPeak,
> {
    frame_reader: TDFFrameReaderType<C, D>,
    peak_merging_tolerance: Tolerance,
    do_consolidate_peaks: bool,
    _cp: PhantomData<CP>,
    _dp: PhantomData<DP>,
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        CP: CentroidLike + From<CentroidPeak>,
        DP: DeconvolutedCentroidLike + From<DeconvolutedPeak>,
    > IntoIonMobilityFrameSource<CP, DP> for TDFSpectrumReaderType<C, D, CP, DP>
{
    type IonMobilityFrameSource<
        CF: FeatureLike<MZ, IonMobility>,
        DF: FeatureLike<Mass, IonMobility> + KnownCharge,
    > = TDFFrameReaderType<CF, DF>;

    fn try_into_frame_source<
        CF: FeatureLike<MZ, IonMobility>,
        DF: FeatureLike<Mass, IonMobility> + KnownCharge,
    >(
        self,
    ) -> Result<Self::IonMobilityFrameSource<CF, DF>, crate::io::IntoIonMobilityFrameSourceError>
    {
        let view = self.into_frame_reader();

        Ok(TDFFrameReaderType {
            tdf_reader: view.tdf_reader,
            metadata: view.metadata,
            frame_reader: view.frame_reader,
            entry_index: view.entry_index,
            index: view.index,
            offset_index: view.offset_index,
            file_description: view.file_description,
            instrument_configurations: view.instrument_configurations,
            softwares: view.softwares,
            samples: view.samples,
            data_processings: view.data_processings,
            detail_level: view.detail_level,
            run: view.run,
            _c: PhantomData,
            _d: PhantomData,
        })
    }
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        CP: CentroidLike + From<CentroidPeak>,
        DP: DeconvolutedCentroidLike,
    > RandomAccessSpectrumIterator<CP, DP, MultiLayerSpectrum<CP, DP>>
    for TDFSpectrumReaderType<C, D, CP, DP>
{
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, SpectrumAccessError> {
        self.frame_reader
            .start_from_id(id)
            .map_err(SpectrumAccessError::from)?;
        Ok(self)
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, SpectrumAccessError> {
        self.frame_reader
            .start_from_index(index)
            .map_err(SpectrumAccessError::from)?;
        Ok(self)
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, SpectrumAccessError> {
        self.frame_reader
            .start_from_time(time)
            .map_err(SpectrumAccessError::from)?;
        Ok(self)
    }
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        CP: CentroidLike + From<CentroidPeak>,
        DP: DeconvolutedCentroidLike,
    > SpectrumSource<CP, DP> for TDFSpectrumReaderType<C, D, CP, DP>
{
    fn reset(&mut self) {
        self.frame_reader.reset();
    }

    fn get_spectrum_by_id(&mut self, id: &str) -> Option<MultiLayerSpectrum<CP, DP>> {
        self.frame_reader.get_frame_by_id(id).map(|f| {
            self.frame_reader.frame_to_spectrum(
                f,
                self.peak_merging_tolerance,
                self.do_consolidate_peaks,
            )
        })
    }

    fn get_spectrum_by_index(&mut self, index: usize) -> Option<MultiLayerSpectrum<CP, DP>> {
        self.frame_reader.get_frame_by_index(index).map(|f| {
            self.frame_reader.frame_to_spectrum(
                f,
                self.peak_merging_tolerance,
                self.do_consolidate_peaks,
            )
        })
    }

    fn get_spectrum_by_time(&mut self, time: f64) -> Option<MultiLayerSpectrum<CP, DP>> {
        self.frame_reader.get_frame_by_time(time).map(|f| {
            self.frame_reader.frame_to_spectrum(
                f,
                self.peak_merging_tolerance,
                self.do_consolidate_peaks,
            )
        })
    }

    fn get_index(&self) -> &OffsetIndex {
        self.frame_reader.get_index()
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.frame_reader.set_index(index);
    }

    fn detail_level(&self) -> &DetailLevel {
        &self.frame_reader.detail_level
    }

    fn set_detail_level(&mut self, detail_level: DetailLevel) {
        self.frame_reader.set_detail_level(detail_level);
    }
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        CP: CentroidLike + From<CentroidPeak>,
        DP: DeconvolutedCentroidLike,
    > Iterator for TDFSpectrumReaderType<C, D, CP, DP>
{
    type Item = MultiLayerSpectrum<CP, DP>;

    fn next(&mut self) -> Option<Self::Item> {
        self.frame_reader.next().map(|f| {
            self.frame_reader.frame_to_spectrum(
                f,
                self.peak_merging_tolerance,
                self.do_consolidate_peaks,
            )
        })
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.frame_reader.nth(n).map(|f| {
            self.frame_reader.frame_to_spectrum(
                f,
                self.peak_merging_tolerance,
                self.do_consolidate_peaks,
            )
        })
    }
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        CP: CentroidLike + From<CentroidPeak>,
        DP: DeconvolutedCentroidLike,
    > FusedIterator for TDFSpectrumReaderType<C, D, CP, DP>
{
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        CP: CentroidLike,
        DP: DeconvolutedCentroidLike,
    > MSDataFileMetadata for TDFSpectrumReaderType<C, D, CP, DP>
{
    crate::delegate_impl_metadata_trait!(frame_reader);
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        CP: CentroidLike,
        DP: DeconvolutedCentroidLike,
    > ChromatogramSource for TDFSpectrumReaderType<C, D, CP, DP>
{
    fn get_chromatogram_by_id(&mut self, id: &str) -> Option<Chromatogram> {
        self.frame_reader.get_chromatogram_by_id(id)
    }

    fn get_chromatogram_by_index(&mut self, index: usize) -> Option<Chromatogram> {
        self.frame_reader.get_chromatogram_by_index(index)
    }
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        CP: CentroidLike,
        DP: DeconvolutedCentroidLike,
    > TDFSpectrumReaderType<C, D, CP, DP>
{
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, TimsRustError> {
        Self::new_with_peak_merging_tolerance(path, PEAK_MERGE_TOLERANCE, true)
    }

    pub fn new_with_peak_merging_tolerance<P: AsRef<Path>>(
        path: P,
        peak_merging_tolerance: Tolerance,
        do_consolidate_peaks: bool,
    ) -> Result<Self, TimsRustError> {
        TDFFrameReaderType::<C, D>::new(path).map(|s| Self {
            frame_reader: s,
            peak_merging_tolerance,
            do_consolidate_peaks,
            _cp: PhantomData,
            _dp: PhantomData,
        })
    }

    pub fn new_with_detail_level<P: AsRef<Path>>(
        path: P,
        peak_merging_tolerance: Option<Tolerance>,
        do_consolidate_peaks: bool,
        detail_level: DetailLevel,
    ) -> Result<Self, TimsRustError> {
        TDFFrameReaderType::<C, D>::new(path).map(|mut s| {
            s.detail_level = detail_level;
            Self {
                frame_reader: s,
                peak_merging_tolerance: peak_merging_tolerance.unwrap_or(PEAK_MERGE_TOLERANCE),
                do_consolidate_peaks,
                _cp: PhantomData,
                _dp: PhantomData,
            }
        })
    }

    /// The number of spectra available
    pub fn len(&self) -> usize {
        self.frame_reader.len()
    }

    pub fn is_empty(&self) -> bool {
        self.frame_reader.is_empty()
    }

    /// Retrieve a spectrum by index
    pub fn get(&self, index: usize) -> Result<Option<MultiLayerSpectrum>, TimsRustError> {
        self.frame_reader.get(index).map(|f| {
            f.map(|f| {
                self.frame_reader.frame_to_spectrum(
                    f,
                    self.peak_merging_tolerance,
                    self.do_consolidate_peaks,
                )
            })
        })
    }

    /// Consume the spectrum reader, retrieving the underlying [`TDFFrameReaderType`]
    pub fn into_inner(self) -> TDFFrameReaderType<C, D> {
        self.frame_reader
    }

    /// Retrieving a mutable reference the underlying [`TDFFrameReaderType`]
    pub fn frame_reader(&mut self) -> &mut TDFFrameReaderType<C, D> {
        &mut self.frame_reader
    }

    /// Consume the spectrum reader, retrieving the underlying [`TDFFrameReaderType`]
    ///
    /// # See also
    /// An alias of [`TDFSpectrumReaderType::into_inner`]
    pub fn into_frame_reader(self) -> TDFFrameReaderType<C, D> {
        self.into_inner()
    }

    /// Get the error tolerance margin for consolidating peaks over the
    /// ion mobility dimension.
    pub fn peak_merging_tolerance(&self) -> Tolerance {
        self.peak_merging_tolerance
    }

    /// Get a mutable reference to the error tolerance margin for consolidating
    /// peaks over the ion mobility dimension.
    pub fn peak_merging_tolerance_mut(&mut self) -> &mut Tolerance {
        &mut self.peak_merging_tolerance
    }

    pub fn get_trace_reader(&self) -> Result<ChromatographyData, TimsRustError> {
        self.frame_reader.get_trace_reader()
    }

    pub fn will_consolidate_peaks(&self) -> bool {
        self.do_consolidate_peaks
    }

    pub fn set_consolidate_peaks(&mut self, do_consolidate_peaks: bool) {
        self.do_consolidate_peaks = do_consolidate_peaks;
    }
}

pub type TDFSpectrumReader = TDFSpectrumReaderType<
    Feature<MZ, IonMobility>,
    ChargedFeature<Mass, IonMobility>,
    CentroidPeak,
    DeconvolutedPeak,
>;

fn index_to_precursor(
    index_entry: &IndexExtry,
    metadata: &Metadata,
    parent_entry: &IndexExtry,
) -> Option<Precursor> {
    if let Some(prec) = index_entry.precursor() {
        let mut ion = SelectedIon {
            mz: prec.mz,
            intensity: prec.intensity as f32,
            charge: if prec.charge != 0 {
                Some(prec.charge)
            } else {
                None
            },
            ..Default::default()
        };

        let im = metadata.im_converter.convert(prec.scan_average);

        let p = inverse_reduce_ion_mobility_param(im);
        ion.add_param(p);
        let mut act = Activation::default();
        let mut isolation = IsolationWindow::default();
        if let Some(pasef) = index_entry.pasef_msms() {
            act.energy = pasef.collision_energy as f32;
            act.methods_mut().push(CollisionInducedDissociation);

            let iso_width = pasef.isolation_width / 2.0;
            isolation.target = pasef.isolation_mz as f32;
            isolation.lower_bound = (pasef.isolation_mz - iso_width) as f32;
            isolation.upper_bound = (pasef.isolation_mz + iso_width) as f32;
            isolation.flags = IsolationWindowState::Complete;
        }
        if let Some(pasef) = index_entry.dia_window() {
            act.energy = pasef.collision_energy;
            act.methods_mut().push(CollisionInducedDissociation);

            let iso_width = pasef.isolation_width / 2.0;
            isolation.target = pasef.isolation_mz as f32;
            isolation.lower_bound = (pasef.isolation_mz - iso_width) as f32;
            isolation.upper_bound = (pasef.isolation_mz + iso_width) as f32;
            isolation.flags = IsolationWindowState::Complete;
        }

        let mut mz_prec = Precursor::default();
        mz_prec.add_ion(ion);
        mz_prec.activation = act;
        mz_prec.isolation_window = isolation;
        mz_prec.precursor_id = Some(parent_entry.format_native_id());
        Some(mz_prec)
    } else if let Some(dia_window) = index_entry.dia_window() {
        let mut ion = SelectedIon {
            mz: dia_window.isolation_mz,
            intensity: 0.0,
            charge: None,
            ..Default::default()
        };
        let im = metadata
            .im_converter
            .convert((dia_window.scan_end + dia_window.scan_start) as f64 / 2.0);

        let p = inverse_reduce_ion_mobility_param(im);
        ion.add_param(p);
        let mut act = Activation::default();
        let mut isolation = IsolationWindow::default();
        if let Some(pasef) = index_entry.pasef_msms() {
            act.energy = pasef.collision_energy as f32;
            act.methods_mut().push(CollisionInducedDissociation);

            let iso_width = pasef.isolation_width / 2.0;
            isolation.target = pasef.isolation_mz as f32;
            isolation.lower_bound = (pasef.isolation_mz - iso_width) as f32;
            isolation.upper_bound = (pasef.isolation_mz + iso_width) as f32;
            isolation.flags = IsolationWindowState::Complete;
        }
        if let Some(pasef) = index_entry.dia_window() {
            act.energy = pasef.collision_energy;
            act.methods_mut().push(CollisionInducedDissociation);

            let iso_width = pasef.isolation_width / 2.0;
            isolation.target = pasef.isolation_mz as f32;
            isolation.lower_bound = (pasef.isolation_mz - iso_width) as f32;
            isolation.upper_bound = (pasef.isolation_mz + iso_width) as f32;
            isolation.flags = IsolationWindowState::Complete;
        }
        let mut mz_prec = Precursor::default();
        mz_prec.add_ion(ion);
        mz_prec.activation = act;
        mz_prec.isolation_window = isolation;
        mz_prec.precursor_id = Some(parent_entry.format_native_id());
        Some(mz_prec)
    } else {
        None
    }
}

fn frame_to_description(
    metadata: &Metadata,
    index_entry: &IndexExtry,
    frame_slice: Option<Range<u32>>,
) -> IonMobilityFrameDescription {
    let mut descr = IonMobilityFrameDescription {
        id: index_entry.format_native_id(),
        index: index_entry.index,
        ms_level: index_entry.ms_level(),
        signal_continuity: SignalContinuity::Centroid,
        polarity: index_entry.frame.polarity,
        ..Default::default()
    };

    descr.add_param(
        Param::builder()
            .name("total ion current")
            .curie(curie!(MS:1000285))
            .value(index_entry.frame.summed_intensities)
            .build(),
    );
    descr.add_param(
        Param::builder()
            .name("base peak intensity")
            .curie(curie!(MS:1000505))
            .value(index_entry.frame.max_intensity)
            .build(),
    );

    let mut scan = ScanEvent::new(
        index_entry.frame.time / 60.0,
        index_entry.frame.accumulation_time,
        vec![ScanWindow {
            lower_bound: metadata.lower_mz as f32,
            upper_bound: metadata.upper_mz as f32,
        }],
        0,
        None,
    );

    if let Some(scan_range) = frame_slice.as_ref() {
        let im = (metadata.im_converter.convert(scan_range.start)
            + metadata.im_converter.convert(scan_range.end))
            / 2.0;
        let p = inverse_reduce_ion_mobility_param(im);
        scan.add_param(p);
    }

    if let Some(pasef) = index_entry.pasef_msms() {
        let im_low = metadata.im_converter.convert(pasef.scan_start as u32);
        let im_high = metadata.im_converter.convert(pasef.scan_end as u32);

        descr.add_param(
            Param::new_key_value("ion mobility lower limit", im_low)
                .with_unit_t(&Unit::VoltSecondPerSquareCentimeter),
        );
        descr.add_param(
            Param::new_key_value("ion mobility upper limit", im_high)
                .with_unit_t(&Unit::VoltSecondPerSquareCentimeter),
        );

        if frame_slice.is_none() {
            scan.add_param(inverse_reduce_ion_mobility_param((im_low + im_high) / 2.0));
        }
    }

    if let Some(pasef) = index_entry.dia_window() {
        let im_low = metadata.im_converter.convert(pasef.scan_start as u32);
        let im_high = metadata.im_converter.convert(pasef.scan_end as u32);

        descr.add_param(
            Param::new_key_value("ion mobility lower limit", im_low)
                .with_unit_t(&Unit::VoltSecondPerSquareCentimeter),
        );
        descr.add_param(
            Param::new_key_value("ion mobility upper limit", im_high)
                .with_unit_t(&Unit::VoltSecondPerSquareCentimeter),
        );

        scan.add_param(
            Param::builder()
                .name("window group")
                .value(Value::Int(pasef.window_group as i64))
                .build(),
        );

        if frame_slice.is_none() {
            scan.add_param(inverse_reduce_ion_mobility_param((im_low + im_high) / 2.0));
        }
    }

    *descr.acquisition.first_scan_mut().unwrap() = scan;
    descr.acquisition.combination = ScanCombination::Sum;
    descr
}

impl<
        C: FeatureLike<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge,
        CP: CentroidLike + From<CentroidPeak>,
        DP: DeconvolutedCentroidLike,
    > MZFileReader<CP, DP, MultiLayerSpectrum<CP, DP>> for TDFSpectrumReaderType<C, D, CP, DP>
{
    fn construct_index_from_stream(&mut self) -> u64 {
        self.frame_reader.entry_index.len() as u64
    }

    /// The underlying Bruker library requires an explicit file system path to the `.d` directory
    /// and as such this method always fails.
    ///
    /// [`open_path`](Self::open_path) works as normal.
    #[allow(unused)]
    fn open_file(source: std::fs::File) -> io::Result<Self> {
        Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "Cannot read a TDF dataset from an open file handle, only from a directory path",
        ))
    }

    fn open_path<P>(path: P) -> io::Result<Self>
    where
        P: Into<std::path::PathBuf> + Clone,
    {
        let frame_reader = TDFFrameReaderType::new(path.into())
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

        let spectrum_reader = frame_reader.into_spectrum_reader();
        Ok(spectrum_reader)
    }
}

/// Test if a file system path is a TDF directory
///
/// This will first probe the directory for the required
/// files, and only after verifying the files exist, attempt
/// to test the SQL schema.
pub fn is_tdf<P: AsRef<Path>>(path: P) -> bool {
    let path = path.as_ref();

    if !path.extension().map(|e| e == "d").unwrap_or_default() {
        return false;
    }

    if !path.exists() {
        return false;
    }

    let tdf_path = path.join("analysis.tdf");

    if !tdf_path.exists() {
        return false;
    }

    if !path.join("analysis.tdf_bin").exists() {
        return false;
    }

    if MetadataReader::new(tdf_path).is_err() {
        return false;
    }
    true
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_tdf_spectrum() -> io::Result<()> {
        let mut reader = TDFSpectrumReader::open_path("test/data/diaPASEF.d")?;
        eprintln!("{}", reader.len());
        let s = reader.get_spectrum_by_index(0).unwrap();
        assert!(s.peaks.is_some());
        assert_eq!(s.signal_continuity(), SignalContinuity::Centroid);
        assert_eq!(s.ms_level(), 1);
        Ok(())
    }

    #[test]
    fn test_tdf_frame() -> io::Result<()> {
        let mut reader = TDFFrameReader::new("test/data/diaPASEF.d")
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        eprintln!("{}", reader.len());
        let s = reader.get_frame_by_index(0).unwrap();
        assert!(s.features.is_none());
        assert_eq!(s.signal_continuity(), SignalContinuity::Centroid);
        assert_eq!(s.ms_level(), 1);

        let s = reader.get_frame_by_index(1).unwrap();
        assert!(s.features.is_none());
        assert_eq!(s.signal_continuity(), SignalContinuity::Centroid);
        assert_eq!(s.ms_level(), 2);
        Ok(())
    }
}
