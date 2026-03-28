use mzdata::params::{Param, ParamDescribed, ParamLike, ValueRef, CURIE, CURIEParsingError, Unit};
use mzdata::spectrum::{Activation, IsolationWindow, Precursor, ScanEvent, ScanWindow, SelectedIon};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyFloat, PyList, PyString, PyInt};

// ---------------------------------------------------------------------------
// PyParam
// ---------------------------------------------------------------------------

#[pyclass(name = "Param", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyParam(pub Param);

#[pymethods]
impl PyParam {
    #[getter]
    fn name(&self) -> &str {
        self.0.name()
    }

    /// The CURIE-formatted accession string (e.g. `"MS:1000579"`), if any.
    #[getter]
    fn id(&self) -> Option<String> {
        self.0.curie_str()
    }

    /// The parameter value as a Python primitive (str, float, int, bool, list, or None).
    #[getter]
    fn value(&self, py: Python<'_>) -> Py<PyAny> {
        match self.0.value() {
            ValueRef::String(s) => PyString::new(py, s.as_ref()).into_any().unbind(),
            ValueRef::Float(f) => PyFloat::new(py, f).into_any().unbind(),
            ValueRef::Int(i) => PyInt::new(py, i).into_any().unbind(),
            ValueRef::Boolean(b) => PyBool::new(py, b).to_owned().into_any().unbind(),
            ValueRef::Buffer(buf) => {
                let list = PyList::new(py, buf.iter().copied().map(|b| b as u32))
                    .expect("Failed to create list");
                list.into_any().unbind()
            }
            ValueRef::Empty => py.None(),
        }
    }

    /// The unit name (e.g. `"minute"`, `"dalton"`, or `"none"`).
    #[getter]
    fn unit(&self) -> String {
        self.0.unit.to_string()
    }

    fn __repr__(&self) -> String {
        let id_part = self
            .0
            .curie_str()
            .map(|s| format!("{s}|"))
            .unwrap_or_default();
        let unit_part = if !matches!(self.0.unit, Unit::Unknown) {
            format!(", unit={}", self.0.unit)
        } else {
            "".to_string()
        };
        format!(
            "Param({}{}, value={}{})",
            id_part,
            self.0.name(),
            self.0.value,
            unit_part
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

impl From<Param> for PyParam {
    fn from(p: Param) -> Self {
        PyParam(p)
    }
}


pub fn find_param(source: &impl ParamDescribed,  name: Option<&str>, accession: Option<&str>) -> PyResult<Option<PyParam>> {
    if name.is_none() && accession.is_none() {
        return Err(pyo3::exceptions::PyTypeError::new_err("Must provide one of `name` or `accession`"));
    }
    if let Some(name) = name {
        Ok(source.params().iter().find(|p| p.name() == name).cloned().map(PyParam))
    } else if let Some(accession) = accession {
        let acc: Option<CURIE> = Some(accession.parse().map_err(|e: CURIEParsingError| pyo3::exceptions::PyValueError::new_err(e.to_string()))?);
        Ok(source.params().iter().find(|p| p.curie() == acc).cloned().map(PyParam))
    } else {
        Ok(None)
    }
}


// ---------------------------------------------------------------------------
// PyIsolationWindow
// ---------------------------------------------------------------------------

#[pyclass(name = "IsolationWindow", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyIsolationWindow(pub IsolationWindow);

#[pymethods]
impl PyIsolationWindow {
    #[getter]
    fn target(&self) -> f32 {
        self.0.target
    }

    #[getter]
    fn lower_bound(&self) -> f32 {
        self.0.lower_bound
    }

    #[getter]
    fn upper_bound(&self) -> f32 {
        self.0.upper_bound
    }

    fn contains(&self, mz: f32) -> bool {
        self.0.contains(mz)
    }

    fn __repr__(&self) -> String {
        format!(
            "IsolationWindow(target={}, lower={}, upper={})",
            self.0.target, self.0.lower_bound, self.0.upper_bound
        )
    }
}

// ---------------------------------------------------------------------------
// PyScanWindow
// ---------------------------------------------------------------------------

#[pyclass(name = "ScanWindow", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyScanWindow(pub ScanWindow);

#[pymethods]
impl PyScanWindow {
    #[getter]
    fn lower_bound(&self) -> f32 {
        self.0.lower_bound
    }

    #[getter]
    fn upper_bound(&self) -> f32 {
        self.0.upper_bound
    }

    fn contains(&self, mz: f32) -> bool {
        self.0.contains(mz)
    }

    fn __repr__(&self) -> String {
        format!(
            "ScanWindow(lower={}, upper={})",
            self.0.lower_bound, self.0.upper_bound
        )
    }
}

// ---------------------------------------------------------------------------
// PySelectedIon
// ---------------------------------------------------------------------------

#[pyclass(name = "SelectedIon", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PySelectedIon(pub SelectedIon);

#[pymethods]
impl PySelectedIon {
    #[getter]
    fn mz(&self) -> f64 {
        self.0.mz
    }

    #[getter]
    fn intensity(&self) -> f32 {
        self.0.intensity
    }

    #[getter]
    fn charge(&self) -> Option<i32> {
        self.0.charge
    }

    fn params(&self) -> Vec<PyParam> {
        self.0.params().iter().cloned().map(PyParam).collect()
    }

    #[pyo3(signature = (name = None, accession = None))]
    fn find_param(&self, name: Option<&str>, accession: Option<&str>) -> PyResult<Option<PyParam>> {
        find_param(&self.0, name, accession)
    }

    fn __repr__(&self) -> String {
        format!(
            "SelectedIon(mz={:.4}, intensity={}, charge={:?})",
            self.0.mz, self.0.intensity, self.0.charge
        )
    }
}

// ---------------------------------------------------------------------------
// PyActivation
// ---------------------------------------------------------------------------

#[pyclass(name = "Activation", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyActivation(pub Activation);

#[pymethods]
impl PyActivation {
    /// The primary activation method name, if present.
    #[getter]
    fn method(&self) -> Option<String> {
        self.0.method().map(|m| m.name().to_string())
    }

    #[getter]
    fn energy(&self) -> f32 {
        self.0.energy
    }

    /// All activation methods as Param objects.
    fn methods(&self) -> Vec<PyParam> {
        self.0
            .methods()
            .iter()
            .map(|m| PyParam(m.to_param().into()))
            .collect()
    }

    fn params(&self) -> Vec<PyParam> {
        self.0.params.iter().cloned().map(PyParam).collect()
    }

    #[pyo3(signature = (name = None, accession = None))]
    fn find_param(&self, name: Option<&str>, accession: Option<&str>) -> PyResult<Option<PyParam>> {
        find_param(&self.0, name, accession)
    }

    fn __repr__(&self) -> String {
        format!(
            "Activation(method={:?}, energy={})",
            self.method(),
            self.0.energy
        )
    }
}

// ---------------------------------------------------------------------------
// PyPrecursor
// ---------------------------------------------------------------------------

#[pyclass(name = "Precursor", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyPrecursor(pub Precursor);

#[pymethods]
impl PyPrecursor {
    /// The native ID of the precursor spectrum, if known.
    #[getter]
    fn precursor_id(&self) -> Option<String> {
        self.0.precursor_id.clone()
    }

    #[getter]
    fn isolation_window(&self) -> PyIsolationWindow {
        PyIsolationWindow(self.0.isolation_window.clone())
    }

    #[getter]
    fn activation(&self) -> PyActivation {
        PyActivation(self.0.activation.clone())
    }

    fn ions(&self) -> Vec<PySelectedIon> {
        self.0.ions.iter().cloned().map(PySelectedIon).collect()
    }

    fn __repr__(&self) -> String {
        let mz: Vec<_> = self.0.ions.iter().map(|i| format!("{:.4}", i.mz)).collect();
        format!(
            "Precursor(ions=[{}], precursor_id={:?})",
            mz.join(", "),
            self.0.precursor_id
        )
    }
}

// ---------------------------------------------------------------------------
// PyScanEvent
// ---------------------------------------------------------------------------

#[pyclass(name = "ScanEvent", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyScanEvent(pub ScanEvent);

#[pymethods]
impl PyScanEvent {
    /// Scan start time in minutes relative to run start.
    #[getter]
    fn start_time(&self) -> f64 {
        self.0.start_time
    }

    /// Ion injection/trapping time in milliseconds.
    #[getter]
    fn injection_time(&self) -> f32 {
        self.0.injection_time
    }

    #[getter]
    fn instrument_configuration_id(&self) -> u32 {
        self.0.instrument_configuration_id
    }

    fn scan_windows(&self) -> Vec<PyScanWindow> {
        self.0.scan_windows.iter().cloned().map(PyScanWindow).collect()
    }

    fn filter_string(&self) -> Option<String> {
        self.0.filter_string().map(|s| s.to_string())
    }

    fn params(&self) -> Vec<PyParam> {
        self.0.params().iter().cloned().map(PyParam).collect()
    }

    #[pyo3(signature = (name = None, accession = None))]
    fn find_param(&self, name: Option<&str>, accession: Option<&str>) -> PyResult<Option<PyParam>> {
        find_param(&self.0, name, accession)
    }

    fn __repr__(&self) -> String {
        format!(
            "ScanEvent(start_time={:.4}, injection_time={})",
            self.0.start_time, self.0.injection_time
        )
    }
}
