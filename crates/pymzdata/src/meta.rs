use mzdata::meta::{
    Component, ComponentType, DataProcessing, FileDescription, InstrumentConfiguration,
    MassSpectrometryRun, ProcessingMethod, Sample, Software, SourceFile,
};
use mzdata::prelude::*;
use mzdata::params::{CURIEParsingError, ParamDescribed, CURIE};
use pyo3::prelude::*;

use crate::params::{PyParam, find_param};

// ---------------------------------------------------------------------------
// PyComponentType
// ---------------------------------------------------------------------------

#[pyclass(name = "ComponentType", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PyComponentType(pub ComponentType);

impl From<ComponentType> for PyComponentType {
    fn from(v: ComponentType) -> Self {
        Self(v)
    }
}

impl From<PyComponentType> for ComponentType {
    fn from(v: PyComponentType) -> Self {
        v.0
    }
}

#[allow(non_snake_case)]
#[pymethods]
impl PyComponentType {
    #[classattr]
    fn Analyzer() -> Self {
        Self(ComponentType::Analyzer)
    }

    #[classattr]
    fn IonSource() -> Self {
        Self(ComponentType::IonSource)
    }

    #[classattr]
    fn Detector() -> Self {
        Self(ComponentType::Detector)
    }

    #[classattr]
    fn Unknown() -> Self {
        Self(ComponentType::Unknown)
    }

    fn __repr__(&self) -> String {
        format!("ComponentType.{}", self.0)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __eq__(&self, other: &PyComponentType) -> bool {
        self.0 == other.0
    }
}

// ---------------------------------------------------------------------------
// PyComponent
// ---------------------------------------------------------------------------

#[pyclass(name = "Component", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyComponent(pub Component);

#[pymethods]
impl PyComponent {
    #[getter]
    fn component_type(&self) -> PyComponentType {
        PyComponentType(self.0.component_type)
    }

    #[getter]
    fn order(&self) -> u8 {
        self.0.order
    }

    /// The name of the mass analyzer term, if this is an analyzer component.
    #[getter]
    fn mass_analyzer(&self) -> Option<String> {
        self.0.mass_analyzer().map(|t| t.name().to_string())
    }

    /// The name of the detector type term, if this is a detector component.
    #[getter]
    fn detector(&self) -> Option<String> {
        self.0.detector().map(|t| t.name().to_string())
    }

    /// The name of the ionization type term, if this is an ion source component.
    #[getter]
    fn ionization_type(&self) -> Option<String> {
        self.0.ionization_type().map(|t| t.name().to_string())
    }

    /// The human-readable name of this component's primary term.
    #[getter]
    fn name(&self) -> Option<String> {
        self.0.name().map(|s| s.to_string())
    }

    fn params(&self) -> Vec<PyParam> {
        self.0.params().iter().cloned().map(PyParam).collect()
    }

    #[pyo3(signature = (name = None, accession = None))]
    fn find_param(
        &self,
        name: Option<&str>,
        accession: Option<&str>,
    ) -> PyResult<Option<PyParam>> {
        find_param(&self.0, name, accession)
    }

    fn __repr__(&self) -> String {
        format!(
            "Component(type={}, order={})",
            self.0.component_type, self.0.order
        )
    }
}

// ---------------------------------------------------------------------------
// PyInstrumentConfiguration
// ---------------------------------------------------------------------------

#[pyclass(name = "InstrumentConfiguration", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyInstrumentConfiguration(pub InstrumentConfiguration);

#[pymethods]
impl PyInstrumentConfiguration {
    #[getter]
    fn id(&self) -> u32 {
        self.0.id
    }

    #[getter]
    fn software_reference(&self) -> String {
        self.0.software_reference.clone()
    }

    fn components(&self) -> Vec<PyComponent> {
        self.0.components.iter().cloned().map(PyComponent).collect()
    }

    fn params(&self) -> Vec<PyParam> {
        self.0.params().iter().cloned().map(PyParam).collect()
    }

    #[pyo3(signature = (name = None, accession = None))]
    fn find_param(
        &self,
        name: Option<&str>,
        accession: Option<&str>,
    ) -> PyResult<Option<PyParam>> {
        find_param(&self.0, name, accession)
    }

    fn __repr__(&self) -> String {
        format!(
            "InstrumentConfiguration(id={}, components={})",
            self.0.id,
            self.0.components.len()
        )
    }
}

// ---------------------------------------------------------------------------
// PySourceFile
// ---------------------------------------------------------------------------

#[pyclass(name = "SourceFile", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PySourceFile(pub SourceFile);

#[pymethods]
impl PySourceFile {
    #[getter]
    fn name(&self) -> &str {
        &self.0.name
    }

    #[getter]
    fn location(&self) -> &str {
        &self.0.location
    }

    #[getter]
    fn id(&self) -> &str {
        &self.0.id
    }

    /// The file format name (e.g. `"Thermo RAW format"`), if known.
    #[getter]
    fn file_format(&self) -> Option<String> {
        self.0.file_format.as_ref().map(|p| p.name().to_string())
    }

    /// The native spectrum ID format name (e.g. `"Thermo nativeID format"`), if known.
    #[getter]
    fn id_format(&self) -> Option<String> {
        self.0
            .id_format
            .as_ref()
            .map(|p| p.name().to_string())
    }

    fn params(&self) -> Vec<PyParam> {
        self.0.params.iter().cloned().map(PyParam).collect()
    }

    #[pyo3(signature = (name = None, accession = None))]
    fn find_param(
        &self,
        name: Option<&str>,
        accession: Option<&str>,
    ) -> PyResult<Option<PyParam>> {
        if name.is_none() && accession.is_none() {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Must provide one of `name` or `accession`",
            ));
        }
        if let Some(n) = name {
            Ok(self
                .0
                .params
                .iter()
                .find(|p| p.name() == n)
                .cloned()
                .map(PyParam))
        } else if let Some(acc) = accession {
            let curie: Option<CURIE> = Some(
                acc.parse()
                    .map_err(|e: CURIEParsingError| {
                        pyo3::exceptions::PyValueError::new_err(e.to_string())
                    })?,
            );
            Ok(self
                .0
                .params
                .iter()
                .find(|p| p.curie() == curie)
                .cloned()
                .map(PyParam))
        } else {
            Ok(None)
        }
    }

    fn __repr__(&self) -> String {
        format!("SourceFile(name={:?}, location={:?})", self.0.name, self.0.location)
    }
}

// ---------------------------------------------------------------------------
// PyFileDescription
// ---------------------------------------------------------------------------

#[pyclass(name = "FileDescription", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyFileDescription(pub FileDescription);

#[pymethods]
impl PyFileDescription {
    /// Whether the file's metadata indicates it contains MS1 spectra.
    #[getter]
    fn has_ms1_spectra(&self) -> bool {
        self.0.has_ms1_spectra()
    }

    /// Whether the file's metadata indicates it contains MSn spectra.
    #[getter]
    fn has_msn_spectra(&self) -> bool {
        self.0.has_msn_spectra()
    }

    fn source_files(&self) -> Vec<PySourceFile> {
        self.0
            .source_files
            .iter()
            .cloned()
            .map(PySourceFile)
            .collect()
    }

    /// The file content descriptors (CV terms describing what kinds of spectra are present).
    fn params(&self) -> Vec<PyParam> {
        self.0.params().iter().cloned().map(PyParam).collect()
    }

    #[pyo3(signature = (name = None, accession = None))]
    fn find_param(
        &self,
        name: Option<&str>,
        accession: Option<&str>,
    ) -> PyResult<Option<PyParam>> {
        if name.is_none() && accession.is_none() {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Must provide one of `name` or `accession`",
            ));
        }
        if let Some(n) = name {
            Ok(self
                .0
                .params()
                .iter()
                .find(|p| p.name() == n)
                .cloned()
                .map(PyParam))
        } else if let Some(acc) = accession {
            let curie: Option<CURIE> = Some(
                acc.parse()
                    .map_err(|e: CURIEParsingError| {
                        pyo3::exceptions::PyValueError::new_err(e.to_string())
                    })?,
            );
            Ok(self
                .0
                .params()
                .iter()
                .find(|p| p.curie() == curie)
                .cloned()
                .map(PyParam))
        } else {
            Ok(None)
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "FileDescription(ms1={}, msn={}, source_files={})",
            self.0.has_ms1_spectra(),
            self.0.has_msn_spectra(),
            self.0.source_files.len()
        )
    }
}

// ---------------------------------------------------------------------------
// PyProcessingMethod
// ---------------------------------------------------------------------------

#[pyclass(name = "ProcessingMethod", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyProcessingMethod(pub ProcessingMethod);

#[pymethods]
impl PyProcessingMethod {
    #[getter]
    fn order(&self) -> i8 {
        self.0.order
    }

    #[getter]
    fn software_reference(&self) -> &str {
        &self.0.software_reference
    }

    fn params(&self) -> Vec<PyParam> {
        self.0.params().iter().cloned().map(PyParam).collect()
    }

    #[pyo3(signature = (name = None, accession = None))]
    fn find_param(
        &self,
        name: Option<&str>,
        accession: Option<&str>,
    ) -> PyResult<Option<PyParam>> {
        if name.is_none() && accession.is_none() {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Must provide one of `name` or `accession`",
            ));
        }
        if let Some(n) = name {
            Ok(self
                .0
                .params()
                .iter()
                .find(|p| p.name() == n)
                .cloned()
                .map(PyParam))
        } else if let Some(acc) = accession {
            let curie: Option<CURIE> = Some(
                acc.parse()
                    .map_err(|e: CURIEParsingError| {
                        pyo3::exceptions::PyValueError::new_err(e.to_string())
                    })?,
            );
            Ok(self
                .0
                .params()
                .iter()
                .find(|p| p.curie() == curie)
                .cloned()
                .map(PyParam))
        } else {
            Ok(None)
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "ProcessingMethod(order={}, software={:?})",
            self.0.order, self.0.software_reference
        )
    }
}

// ---------------------------------------------------------------------------
// PyDataProcessing
// ---------------------------------------------------------------------------

#[pyclass(name = "DataProcessing", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyDataProcessing(pub DataProcessing);

#[pymethods]
impl PyDataProcessing {
    #[getter]
    fn id(&self) -> &str {
        &self.0.id
    }

    fn methods(&self) -> Vec<PyProcessingMethod> {
        self.0
            .methods
            .iter()
            .cloned()
            .map(PyProcessingMethod)
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "DataProcessing(id={:?}, methods={})",
            self.0.id,
            self.0.methods.len()
        )
    }
}

// ---------------------------------------------------------------------------
// PySoftware
// ---------------------------------------------------------------------------

#[pyclass(name = "Software", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PySoftware(pub Software);

#[pymethods]
impl PySoftware {
    #[getter]
    fn id(&self) -> &str {
        &self.0.id
    }

    #[getter]
    fn version(&self) -> &str {
        &self.0.version
    }

    /// Whether this software is categorised as analysis software.
    #[getter]
    fn is_analysis(&self) -> bool {
        self.0.is_analysis()
    }

    /// Whether this software is categorised as data-processing software.
    #[getter]
    fn is_data_processing(&self) -> bool {
        self.0.is_data_processing()
    }

    /// Whether this software is categorised as acquisition software.
    #[getter]
    fn is_acquisition(&self) -> bool {
        self.0.is_acquisition()
    }

    /// The controlled-vocabulary name of the software (e.g. `"Xcalibur"`), if known.
    #[getter]
    fn software_name(&self) -> Option<String> {
        self.0
            .find_software_term()
            .map(|t| t.name().to_string())
    }

    fn params(&self) -> Vec<PyParam> {
        self.0.params().iter().cloned().map(PyParam).collect()
    }

    #[pyo3(signature = (name = None, accession = None))]
    fn find_param(
        &self,
        name: Option<&str>,
        accession: Option<&str>,
    ) -> PyResult<Option<PyParam>> {
        if name.is_none() && accession.is_none() {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Must provide one of `name` or `accession`",
            ));
        }
        if let Some(n) = name {
            Ok(self
                .0
                .params()
                .iter()
                .find(|p| p.name() == n)
                .cloned()
                .map(PyParam))
        } else if let Some(acc) = accession {
            let curie: Option<CURIE> = Some(
                acc.parse()
                    .map_err(|e: CURIEParsingError| {
                        pyo3::exceptions::PyValueError::new_err(e.to_string())
                    })?,
            );
            Ok(self
                .0
                .params()
                .iter()
                .find(|p| p.curie() == curie)
                .cloned()
                .map(PyParam))
        } else {
            Ok(None)
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Software(id={:?}, version={:?})",
            self.0.id, self.0.version
        )
    }
}

// ---------------------------------------------------------------------------
// PySample
// ---------------------------------------------------------------------------

#[pyclass(name = "Sample", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PySample(pub Sample);

#[pymethods]
impl PySample {
    #[getter]
    fn id(&self) -> &str {
        &self.0.id
    }

    #[getter]
    fn name(&self) -> Option<&str> {
        self.0.name.as_deref()
    }

    fn params(&self) -> Vec<PyParam> {
        self.0.params().iter().cloned().map(PyParam).collect()
    }

    #[pyo3(signature = (name = None, accession = None))]
    fn find_param(
        &self,
        name: Option<&str>,
        accession: Option<&str>,
    ) -> PyResult<Option<PyParam>> {
        if name.is_none() && accession.is_none() {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Must provide one of `name` or `accession`",
            ));
        }
        if let Some(n) = name {
            Ok(self
                .0
                .params()
                .iter()
                .find(|p| p.name() == n)
                .cloned()
                .map(PyParam))
        } else if let Some(acc) = accession {
            let curie: Option<CURIE> = Some(
                acc.parse()
                    .map_err(|e: CURIEParsingError| {
                        pyo3::exceptions::PyValueError::new_err(e.to_string())
                    })?,
            );
            Ok(self
                .0
                .params()
                .iter()
                .find(|p| p.curie() == curie)
                .cloned()
                .map(PyParam))
        } else {
            Ok(None)
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Sample(id={:?}, name={:?})",
            self.0.id, self.0.name
        )
    }
}

// ---------------------------------------------------------------------------
// PyMSRun
// ---------------------------------------------------------------------------

#[pyclass(name = "MSRun", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyMSRun(pub MassSpectrometryRun);

#[pymethods]
impl PyMSRun {
    #[getter]
    fn id(&self) -> Option<&str> {
        self.0.id.as_deref()
    }

    #[getter]
    fn default_data_processing_id(&self) -> Option<&str> {
        self.0.default_data_processing_id.as_deref()
    }

    #[getter]
    fn default_instrument_id(&self) -> Option<u32> {
        self.0.default_instrument_id
    }

    #[getter]
    fn default_source_file_id(&self) -> Option<&str> {
        self.0.default_source_file_id.as_deref()
    }

    /// The run start time as an RFC 3339 string, if present.
    #[getter]
    fn start_time(&self) -> Option<String> {
        self.0.start_time.map(|t| t.to_rfc3339())
    }

    fn __repr__(&self) -> String {
        format!(
            "MSRun(id={:?}, start_time={:?})",
            self.0.id,
            self.0.start_time.map(|t| t.to_rfc3339())
        )
    }
}
