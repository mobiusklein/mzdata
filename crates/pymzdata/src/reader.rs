use std::fs;

use mzdata::io::{DetailLevel, IMMZReaderType, MZReader, MZFileReader};
use mzdata::meta::MSDataFileMetadata;
use mzdata::prelude::*;
use mzpeaks::{
    feature::{ChargedFeature, Feature},
    IonMobility, Mass, MZ,
};
use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;

use crate::meta::{
    PyDataProcessing, PyFileDescription, PyInstrumentConfiguration, PyMSRun, PySample, PySoftware,
};
use crate::spectrum::{PyIonMobilityFrame, PySpectrum, RawIonMobilityFrame};

// ---------------------------------------------------------------------------
// Concrete type aliases for the file-based readers
// ---------------------------------------------------------------------------

type RawMZReader = MZReader<fs::File>;
type RawIMMZReader = IMMZReaderType<
    fs::File,
    Feature<MZ, IonMobility>,
    ChargedFeature<Mass, IonMobility>,
>;

// ---------------------------------------------------------------------------
// Helper: parse detail level string
// ---------------------------------------------------------------------------

fn parse_detail_level(s: &str) -> PyResult<DetailLevel> {
    match s {
        "full" => Ok(DetailLevel::Full),
        "lazy" => Ok(DetailLevel::Lazy),
        "metadata_only" => Ok(DetailLevel::MetadataOnly),
        other => Err(PyValueError::new_err(format!(
            "Unknown detail level {:?}. Choose from 'full', 'lazy', 'metadata_only'.",
            other
        ))),
    }
}

fn detail_level_str(dl: &DetailLevel) -> &'static str {
    match dl {
        DetailLevel::Full => "full",
        DetailLevel::Lazy => "lazy",
        DetailLevel::MetadataOnly => "metadata_only",
    }
}


#[pyclass(name = "DetailLevel", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PyDetailLevel(pub DetailLevel);

impl From<DetailLevel> for PyDetailLevel {
    fn from(value: DetailLevel) -> Self {
        Self(value)
    }
}

impl From<PyDetailLevel> for DetailLevel {
    fn from(value: PyDetailLevel) -> Self {
        value.0
    }
}

#[allow(non_snake_case)]
#[pymethods]
impl PyDetailLevel {
    #[classattr]
    fn Full() -> Self {
        Self(DetailLevel::Full)
    }

    #[classattr]
    fn Lazy() -> Self { Self(DetailLevel::Lazy) }

    #[classattr]
    fn MetadataOnly() -> Self { Self(DetailLevel::MetadataOnly) }

    fn __repr__(&self) -> String {
        format!("DetailLevel.{:?}", self.0)
    }
}

// ---------------------------------------------------------------------------
// PyMZReader
// ---------------------------------------------------------------------------

/// A file-based mass spectrometry reader.
///
/// Supports mzML, MGF, Thermo RAW, and Bruker TDF formats. Use as a context
/// manager (``with`` block) for automatic resource cleanup.
///
/// Example::
///
///     with pymzdata.MZReader("data.mzML") as reader:
///         for spectrum in reader:
///             print(spectrum.id, spectrum.ms_level)
///
#[pyclass(name = "MZReader", module = "pymzdata")]
pub struct PyMZReader {
    inner: Option<RawMZReader>,
}

#[pymethods]
impl PyMZReader {
    /// Open a mass spectrometry file for reading.
    #[new]
    fn open(path: &str) -> PyResult<Self> {
        let reader =
            RawMZReader::open_path(path).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(PyMZReader {
            inner: Some(reader),
        })
    }

    /// Check if the underlying reader is closed
    fn closed(&self) -> bool {
        self.inner.is_none()
    }

    /// Close the underlying reader
    fn close(&mut self) {
        self.inner = None;
    }

    // ---- Context manager --------------------------------------------------

    fn __enter__(slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf
    }

    fn __exit__(
        &mut self,
        _exc_type: Option<Bound<'_, PyAny>>,
        _exc_value: Option<Bound<'_, PyAny>>,
        _traceback: Option<Bound<'_, PyAny>>,
    ) -> bool {
        self.inner = None;
        false
    }

    // ---- Iterator protocol ------------------------------------------------

    fn __iter__(slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<PySpectrum> {
        self.inner.as_mut()?.next().map(|s| s.into())
    }

    // ---- Length ------------------------------------------------------------

    fn __len__(&self) -> usize {
        self.inner.as_ref().map(|r| r.len()).unwrap_or(0)
    }

    // ---- Random access ----------------------------------------------------

    /// Retrieve a spectrum by its native ID string.
    fn get_by_id(&mut self, id: &str) -> PyResult<Option<PySpectrum>> {
        self.require_open()?;
        Ok(self
            .inner
            .as_mut()
            .unwrap()
            .get_spectrum_by_id(id)
            .map(|s| s.into()))
    }

    /// Retrieve a spectrum by its zero-based ordinal index.
    fn get_by_index(&mut self, index: usize) -> PyResult<Option<PySpectrum>> {
        self.require_open()?;
        Ok(self
            .inner
            .as_mut()
            .unwrap()
            .get_spectrum_by_index(index)
            .map(|s| s.into()))
    }

    /// Retrieve the spectrum whose retention time is closest to `rt` (in minutes).
    fn get_by_time(&mut self, rt: f64) -> PyResult<Option<PySpectrum>> {
        self.require_open()?;
        Ok(self
            .inner
            .as_mut()
            .unwrap()
            .get_spectrum_by_time(rt)
            .map(|s| s.into()))
    }

    // ---- Detail level property --------------------------------------------

    /// The detail level used when loading spectra.
    ///
    /// One of ``"full"``, ``"lazy"``, or ``"metadata_only"``.
    #[getter]
    fn detail_level(&self) -> PyResult<PyDetailLevel> {
        let inner = self.inner.as_ref().ok_or_else(|| {
            PyRuntimeError::new_err("Reader is closed")
        })?;
        Ok((*inner.detail_level()).into())
    }

    #[setter]
    fn set_detail_level(&mut self, level: &PyDetailLevel) -> PyResult<()> {
        let inner = self.inner.as_mut().ok_or_else(|| {
            PyRuntimeError::new_err("Reader is closed")
        })?;
        inner.set_detail_level(level.0);
        Ok(())
    }

    /// Convert this reader to an IMMZReader for ion mobility frame access.
    ///
    /// Raises ``RuntimeError`` if the file does not contain ion mobility data.
    fn into_frame_reader(&mut self) -> PyResult<PyIMMZReader> {
        let inner = self.inner.take().ok_or_else(|| {
            PyRuntimeError::new_err("Reader is already closed")
        })?;
        let frame_source = inner
            .try_into_frame_source::<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>>()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(PyIMMZReader {
            inner: Some(frame_source),
        })
    }

    // ---- File metadata (MSDataFileMetadata) -----------------------------------

    /// The file-level description: content types and source files.
    fn file_description(&self) -> PyResult<PyFileDescription> {
        self.require_open()?;
        Ok(PyFileDescription(
            self.inner.as_ref().unwrap().file_description().clone(),
        ))
    }

    /// All instrument configurations keyed by their id, returned as a list sorted by id.
    fn instrument_configurations(&self) -> PyResult<Vec<PyInstrumentConfiguration>> {
        self.require_open()?;
        let mut configs: Vec<_> = self
            .inner
            .as_ref()
            .unwrap()
            .instrument_configurations()
            .values()
            .cloned()
            .map(PyInstrumentConfiguration)
            .collect();
        configs.sort_by_key(|c| c.0.id);
        Ok(configs)
    }

    /// All data processing pipelines defined in the file.
    fn data_processings(&self) -> PyResult<Vec<PyDataProcessing>> {
        self.require_open()?;
        Ok(self
            .inner
            .as_ref()
            .unwrap()
            .data_processings()
            .iter()
            .cloned()
            .map(PyDataProcessing)
            .collect())
    }

    /// All software entries defined in the file.
    fn softwares(&self) -> PyResult<Vec<PySoftware>> {
        self.require_open()?;
        Ok(self
            .inner
            .as_ref()
            .unwrap()
            .softwares()
            .iter()
            .cloned()
            .map(PySoftware)
            .collect())
    }

    /// All sample entries defined in the file.
    fn samples(&self) -> PyResult<Vec<PySample>> {
        self.require_open()?;
        Ok(self
            .inner
            .as_ref()
            .unwrap()
            .samples()
            .iter()
            .cloned()
            .map(PySample)
            .collect())
    }

    /// The run-level metadata, if present.
    fn run_description(&self) -> PyResult<Option<PyMSRun>> {
        self.require_open()?;
        Ok(self
            .inner
            .as_ref()
            .unwrap()
            .run_description()
            .cloned()
            .map(PyMSRun))
    }

    fn __repr__(&self) -> String {
        if self.inner.is_some() {
            "MZReader(<open>)".to_string()
        } else {
            "MZReader(<closed>)".to_string()
        }
    }
}

impl PyMZReader {
    fn require_open(&self) -> PyResult<()> {
        if self.inner.is_none() {
            Err(PyRuntimeError::new_err("Reader is closed"))
        } else {
            Ok(())
        }
    }
}

// ---------------------------------------------------------------------------
// PyIMMZReader
// ---------------------------------------------------------------------------

/// A file-based ion mobility frame reader.
///
/// Opens data files containing ion mobility data, such as Bruker TDF files or mzML files.
///
/// Use as a context manager (``with`` block) for automatic
/// resource cleanup.
///
/// Can also be obtained from :meth:`MZReader.into_frame_reader`.
///
/// Example::
///
///     with pymzdata.IMMZReader("data.d") as reader:
///         for frame in reader:
///             print(frame.id, frame.ms_level)
///
#[pyclass(name = "IMMZReader", module = "pymzdata")]
pub struct PyIMMZReader {
    pub inner: Option<RawIMMZReader>,
}

#[pymethods]
impl PyIMMZReader {
    /// Open an ion mobility data file for reading.
    #[new]
    fn open(path: &str) -> PyResult<Self> {
        // Open via MZReader and convert to IM frame source.
        let reader =
            RawMZReader::open_path(path).map_err(|e| PyIOError::new_err(e.to_string()))?;
        let frame_source = reader
            .try_into_frame_source::<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>>()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(PyIMMZReader {
            inner: Some(frame_source),
        })
    }

    /// Check if the underlying reader is closed
    fn closed(&self) -> bool {
        self.inner.is_none()
    }

    /// Close the underlying reader
    fn close(&mut self) {
        self.inner = None;
    }

    // ---- Context manager --------------------------------------------------

    fn __enter__(slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf
    }

    fn __exit__(
        &mut self,
        _exc_type: Option<Bound<'_, PyAny>>,
        _exc_value: Option<Bound<'_, PyAny>>,
        _traceback: Option<Bound<'_, PyAny>>,
    ) -> bool {
        self.inner = None;
        false
    }

    // ---- Iterator protocol ------------------------------------------------

    fn __iter__(slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<PyIonMobilityFrame> {
        let frame: RawIonMobilityFrame = self.inner.as_mut()?.next()?;
        Some(frame.into())
    }

    // ---- Length ------------------------------------------------------------

    fn __len__(&self) -> usize {
        self.inner
            .as_ref()
            .map(|r| r.get_index().len())
            .unwrap_or(0)
    }

    // ---- Random access ----------------------------------------------------

    /// Retrieve a frame by its native ID string.
    fn get_by_id(&mut self, id: &str) -> PyResult<Option<PyIonMobilityFrame>> {
        self.require_open()?;
        Ok(self
            .inner
            .as_mut()
            .unwrap()
            .get_frame_by_id(id)
            .map(|f| f.into()))
    }

    /// Retrieve a frame by its zero-based ordinal index.
    fn get_by_index(&mut self, index: usize) -> PyResult<Option<PyIonMobilityFrame>> {
        self.require_open()?;
        Ok(self
            .inner
            .as_mut()
            .unwrap()
            .get_frame_by_index(index)
            .map(|f| f.into()))
    }

    // ---- Detail level property --------------------------------------------

    #[getter]
    fn detail_level(&self) -> PyResult<&'static str> {
        let inner = self.inner.as_ref().ok_or_else(|| {
            PyRuntimeError::new_err("Reader is closed")
        })?;
        Ok(detail_level_str(inner.detail_level()))
    }

    #[setter]
    fn set_detail_level(&mut self, level: &str) -> PyResult<()> {
        let inner = self.inner.as_mut().ok_or_else(|| {
            PyRuntimeError::new_err("Reader is closed")
        })?;
        inner.set_detail_level(parse_detail_level(level)?);
        Ok(())
    }

    // ---- File metadata (MSDataFileMetadata) -----------------------------------

    /// The file-level description: content types and source files.
    fn file_description(&self) -> PyResult<PyFileDescription> {
        self.require_open()?;
        Ok(PyFileDescription(
            self.inner.as_ref().unwrap().file_description().clone(),
        ))
    }

    /// All instrument configurations keyed by their id, returned as a list sorted by id.
    fn instrument_configurations(&self) -> PyResult<Vec<PyInstrumentConfiguration>> {
        self.require_open()?;
        let mut configs: Vec<_> = self
            .inner
            .as_ref()
            .unwrap()
            .instrument_configurations()
            .values()
            .cloned()
            .map(PyInstrumentConfiguration)
            .collect();
        configs.sort_by_key(|c| c.0.id);
        Ok(configs)
    }

    /// All data processing pipelines defined in the file.
    fn data_processings(&self) -> PyResult<Vec<PyDataProcessing>> {
        self.require_open()?;
        Ok(self
            .inner
            .as_ref()
            .unwrap()
            .data_processings()
            .iter()
            .cloned()
            .map(PyDataProcessing)
            .collect())
    }

    /// All software entries defined in the file.
    fn softwares(&self) -> PyResult<Vec<PySoftware>> {
        self.require_open()?;
        Ok(self
            .inner
            .as_ref()
            .unwrap()
            .softwares()
            .iter()
            .cloned()
            .map(PySoftware)
            .collect())
    }

    /// All sample entries defined in the file.
    fn samples(&self) -> PyResult<Vec<PySample>> {
        self.require_open()?;
        Ok(self
            .inner
            .as_ref()
            .unwrap()
            .samples()
            .iter()
            .cloned()
            .map(PySample)
            .collect())
    }

    /// The run-level metadata, if present.
    fn run_description(&self) -> PyResult<Option<PyMSRun>> {
        self.require_open()?;
        Ok(self
            .inner
            .as_ref()
            .unwrap()
            .run_description()
            .cloned()
            .map(PyMSRun))
    }

    fn __repr__(&self) -> String {
        if self.inner.is_some() {
            "IMMZReader(<open>)".to_string()
        } else {
            "IMMZReader(<closed>)".to_string()
        }
    }
}

impl PyIMMZReader {
    fn require_open(&self) -> PyResult<()> {
        if self.inner.is_none() {
            Err(PyRuntimeError::new_err("Reader is closed"))
        } else {
            Ok(())
        }
    }
}
