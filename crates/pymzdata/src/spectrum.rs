use mzdata::prelude::{IonMobilityFrameLike, SpectrumLike, ParamDescribed, ParamLike};
use mzdata::params::{CURIE, CURIEParsingError};
use mzdata::spectrum::{
    MultiLayerIonMobilityFrame, MultiLayerSpectrum, SignalContinuity,
};
use mzdata::io::proxi;

use mzpeaks::{
    feature::{ChargedFeature, Feature},
    prelude::*,
    CentroidPeak, DeconvolutedPeak, IonMobility, Mass, MZ,
};

use pyo3::prelude::*;
use numpy::{IntoPyArray, PyArray1, self};

use crate::params::{PyParam, PyPrecursor, PyScanEvent};

#[pyclass(name = "SignalContinuity", module = "pymzdata", skip_from_py_object, eq, )]
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PySignalContinuity {
    #[default]
    Unknown = 0,
    /// The spectrum is centroided, indicating that its primary representation is that of a
    /// discrete peak list. There may be multiple peak lists and a profile spectrum may still
    /// be present on the same spectrum.
    Centroid = 3,
    /// The spectrum is profile, indicating that its primary representation is a continuous
    /// profile.
    Profile = 5,
}

impl std::fmt::Display for PySignalContinuity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SignalContinuity.{:?}", self)
    }
}

#[pymethods]
impl PySignalContinuity {
    fn __repr__(&self) -> String {
        format!("{}", SignalContinuity::from(*self))
    }
}

impl From<SignalContinuity> for PySignalContinuity {
    fn from(value: SignalContinuity) -> Self {
        match value {
            SignalContinuity::Unknown => Self::Unknown,
            SignalContinuity::Centroid => Self::Centroid,
            SignalContinuity::Profile => Self::Profile,
        }
    }
}

impl From<PySignalContinuity> for SignalContinuity {
    fn from(value: PySignalContinuity) -> Self {
        match value {
            PySignalContinuity::Unknown => Self::Unknown,
            PySignalContinuity::Centroid => Self::Centroid,
            PySignalContinuity::Profile => Self::Profile,
        }
    }
}

// ---------------------------------------------------------------------------
// PySpectrum
// ---------------------------------------------------------------------------

#[pyclass(name = "Spectrum", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PySpectrum {
    pub inner: MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak>,
}

impl From<MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak>> for PySpectrum {
    fn from(inner: MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak>) -> Self {
        PySpectrum { inner }
    }
}

#[pymethods]
impl PySpectrum {
    #[getter]
    fn id(&self) -> &str {
        self.inner.id()
    }

    #[getter]
    fn index(&self) -> usize {
        self.inner.index()
    }

    #[getter]
    fn ms_level(&self) -> u8 {
        self.inner.ms_level()
    }

    #[getter]
    fn start_time(&self) -> f64 {
        self.inner.start_time()
    }

    #[getter]
    fn signal_continuity(&self) -> PySignalContinuity {
        self.inner.signal_continuity().into()
    }

    #[getter]
    fn is_profile(&self) -> bool {
        matches!(self.inner.signal_continuity(), SignalContinuity::Profile)
    }

    #[getter]
    fn is_centroid(&self) -> bool {
        matches!(self.inner.signal_continuity(), SignalContinuity::Centroid)
    }

    /// The first precursor, if any.
    #[getter]
    fn precursor(&self) -> Option<PyPrecursor> {
        self.inner.precursor().cloned().map(PyPrecursor)
    }

    fn params(&self) -> Vec<PyParam> {
        self.inner
            .description()
            .params()
            .iter()
            .cloned()
            .map(PyParam)
            .collect()
    }

    #[pyo3(signature = (name = None, accession = None))]
    fn find_param(&self, name: Option<&str>, accession: Option<&str>) -> PyResult<Option<PyParam>> {
        use mzdata::params::ParamDescribed;
        if name.is_none() && accession.is_none() {
            return Err(pyo3::exceptions::PyTypeError::new_err("Must provide one of `name` or `accession`"));
        }
        if let Some(name) = name {
            Ok(self.inner.description().params().iter().find(|p| p.name() == name).cloned().map(PyParam))
        } else if let Some(accession) = accession {
            let acc: Option<CURIE> = Some(accession.parse().map_err(|e: CURIEParsingError| pyo3::exceptions::PyValueError::new_err(e.to_string()))?);
            Ok(self.inner.description().params().iter().find(|p| p.curie() == acc).cloned().map(PyParam))
        } else {
            Ok(None)
        }
    }

    fn scan_events(&self) -> Vec<PyScanEvent> {
        self.inner
            .description()
            .acquisition
            .iter()
            .cloned()
            .map(PyScanEvent)
            .collect()
    }

    /// Decoded m/z array from the raw binary data, if present.
    fn mz_array<'py>(slf: &'py Bound<'py, Self>) -> Option<Bound<'py, PyArray1<f64>>> {
        slf.borrow().inner.arrays.as_ref().and_then(|a| a.mzs().ok()).map(|cow| cow.to_vec().into_pyarray(slf.py()))
    }

    /// Decoded intensity array from the raw binary data, if present.
    fn intensity_array<'py>(slf: &'py Bound<'py, Self>) -> Option<Bound<'py, PyArray1<f32>>> {
        slf.borrow().inner
            .arrays
            .as_ref()
            .and_then(|a| a.intensities().ok())
            .map(|cow| cow.to_vec().into_pyarray(slf.py()))
    }

    /// Centroid peaks as a list of `(mz, intensity)` tuples.
    fn centroid_peaks<'py>(slf: &'py Bound<'py, Self>) -> Option<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f32>>)> {
        if let Some(peaks) = slf.borrow().inner.peaks.as_ref() {
            let mut mzs = Vec::with_capacity(peaks.len());
            let mut ints = Vec::with_capacity(peaks.len());
            for p in peaks.iter() {
                mzs.push(p.mz());
                ints.push(p.intensity());
            }
            let mzs = mzs.into_pyarray(slf.py());
            let ints = ints.into_pyarray(slf.py());
            Some((mzs, ints))
        } else {
            None
        }
    }

    /// Deconvoluted peaks as a list of `(mz, intensity, charge)` tuples.
    fn deconvoluted_peaks<'py>(slf: &'py Bound<'py, Self>) -> Option<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f32>>, Bound<'py, PyArray1<i32>>)> {
        if let Some(peaks) = slf.borrow().inner.deconvoluted_peaks.as_ref() {
            let mut mzs = Vec::with_capacity(peaks.len());
            let mut ints = Vec::with_capacity(peaks.len());
            let mut zs = Vec::with_capacity(peaks.len());
            for p in peaks.iter() {
                mzs.push(p.mz());
                ints.push(p.intensity());
                zs.push(p.charge())
            }
            let mzs = mzs.into_pyarray(slf.py());
            let ints = ints.into_pyarray(slf.py());
            let zs = zs.into_pyarray(slf.py());
            Some((mzs, ints, zs))
        } else {
            None
        }
    }

    /// Peak-pick the profile signal with a given signal-to-noise threshold.
    ///
    /// Requires the `nalgebra` feature (mzsignal backend).
    #[pyo3(signature = (snr_threshold=1.0))]
    fn pick_peaks(&mut self, snr_threshold: f32) -> PyResult<()> {
        self.inner
            .pick_peaks(snr_threshold)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    }

    /// Denoise the profile signal with a local baseline subtraction.
    fn denoise(&mut self, scale: f32) -> PyResult<()> {
        self.inner
            .denoise(scale)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    }

    /// Re-profile centroid peaks onto a uniform m/z grid.
    ///
    /// `dx` is the grid spacing; `fwhm` is the assumed peak width.
    fn reprofile(&mut self, dx: f64, fwhm: f32) -> PyResult<()> {
        self.inner
            .reprofile_with_shape(dx, fwhm)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        self.inner.description_mut().signal_continuity = SignalContinuity::Profile;
        Ok(())
    }

    fn __len__(&self) -> usize {
        use mzdata::spectrum::RefPeakDataLevel;
        match self.inner.peaks() {
            RefPeakDataLevel::Missing => 0,
            RefPeakDataLevel::RawData(_) => {
                self.inner
                    .arrays
                    .as_ref()
                    .and_then(|a| a.mzs().ok())
                    .map(|c| c.len())
                    .unwrap_or(0)
            }
            RefPeakDataLevel::Centroid(p) => p.len(),
            RefPeakDataLevel::Deconvoluted(p) => p.len(),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Spectrum(id={:?}, ms_level={}, rt={:.4}, continuity={})",
            self.inner.id(),
            self.inner.ms_level(),
            self.inner.start_time(),
            self.signal_continuity(),
        )
    }

    fn to_proxi<'py>(slf: &'py Bound<'py, Self>) -> Result<Bound<'py, PyAny>, serde_pyobject::Error> {
        let spec = proxi::PROXISpectrum::from(&slf.borrow().inner);
        serde_pyobject::to_pyobject(slf.py(), &spec)
    }
}

// ---------------------------------------------------------------------------
// PyIonMobilityFrame
// ---------------------------------------------------------------------------

pub type RawIonMobilityFrame =
    MultiLayerIonMobilityFrame<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>>;

#[pyclass(name = "IonMobilityFrame", module = "pymzdata", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyIonMobilityFrame {
    pub inner: RawIonMobilityFrame,
}

impl From<RawIonMobilityFrame> for PyIonMobilityFrame {
    fn from(inner: RawIonMobilityFrame) -> Self {
        PyIonMobilityFrame { inner }
    }
}

#[pymethods]
impl PyIonMobilityFrame {
    #[getter]
    fn id(&self) -> &str {
        self.inner.id()
    }

    #[getter]
    fn index(&self) -> usize {
        self.inner.index()
    }

    #[getter]
    fn ms_level(&self) -> u8 {
        self.inner.ms_level()
    }

    #[getter]
    fn start_time(&self) -> f64 {
        self.inner.start_time()
    }

    #[getter]
    fn signal_continuity(&self) -> PySignalContinuity {
        self.inner.description().signal_continuity.into()
    }

    #[getter]
    fn is_profile(&self) -> bool {
        matches!(
            self.inner.description().signal_continuity,
            SignalContinuity::Profile
        )
    }

    #[getter]
    fn is_centroid(&self) -> bool {
        matches!(
            self.inner.description().signal_continuity,
            SignalContinuity::Centroid
        )
    }

    #[getter]
    fn precursor(&self) -> Option<PyPrecursor> {
        self.inner.precursor().cloned().map(PyPrecursor)
    }

    fn params(&self) -> Vec<PyParam> {

        self.inner
            .description()
            .params()
            .iter()
            .cloned()
            .map(PyParam)
            .collect()
    }

    #[pyo3(signature = (name = None, accession = None))]
    fn find_param(&self, name: Option<&str>, accession: Option<&str>) -> PyResult<Option<PyParam>> {
        use mzdata::params::ParamDescribed;
        if name.is_none() && accession.is_none() {
            return Err(pyo3::exceptions::PyTypeError::new_err("Must provide one of `name` or `accession`"));
        }
        if let Some(name) = name {
            Ok(self.inner.description().params().iter().find(|p| p.name() == name).cloned().map(PyParam))
        } else if let Some(accession) = accession {
            let acc: Option<CURIE> = Some(accession.parse().map_err(|e: CURIEParsingError| pyo3::exceptions::PyValueError::new_err(e.to_string()))?);
            Ok(self.inner.description().params().iter().find(|p| p.curie() == acc).cloned().map(PyParam))
        } else {
            Ok(None)
        }
    }

    fn scan_events(&self) -> Vec<PyScanEvent> {
        self.inner
            .description()
            .acquisition
            .iter()
            .cloned()
            .map(PyScanEvent)
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "IonMobilityFrame(id={:?}, ms_level={}, rt={:.4})",
            self.inner.id(),
            self.inner.ms_level(),
            self.inner.start_time(),
        )
    }
}
