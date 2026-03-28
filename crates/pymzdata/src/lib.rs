mod meta;
mod params;
mod reader;
mod spectrum;

use pyo3::prelude::*;

/// Python bindings for the mzdata mass spectrometry library.
///
/// Provides :class:`MZReader` and :class:`IMMZReader` for reading spectra and
/// ion mobility frames from mzML, MGF, Thermo RAW, and Bruker TDF files.
#[pymodule]
fn pymzdata(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<meta::PyComponentType>()?;
    m.add_class::<meta::PyComponent>()?;
    m.add_class::<meta::PyInstrumentConfiguration>()?;
    m.add_class::<meta::PySourceFile>()?;
    m.add_class::<meta::PyFileDescription>()?;
    m.add_class::<meta::PyProcessingMethod>()?;
    m.add_class::<meta::PyDataProcessing>()?;
    m.add_class::<meta::PySoftware>()?;
    m.add_class::<meta::PySample>()?;
    m.add_class::<meta::PyMSRun>()?;
    m.add_class::<params::PyParam>()?;
    m.add_class::<params::PyIsolationWindow>()?;
    m.add_class::<params::PyScanWindow>()?;
    m.add_class::<params::PySelectedIon>()?;
    m.add_class::<params::PyActivation>()?;
    m.add_class::<params::PyPrecursor>()?;
    m.add_class::<params::PyScanEvent>()?;
    m.add_class::<spectrum::PySpectrum>()?;
    m.add_class::<spectrum::PyIonMobilityFrame>()?;
    m.add_class::<spectrum::PySignalContinuity>()?;
    m.add_class::<reader::PyMZReader>()?;
    m.add_class::<reader::PyIMMZReader>()?;
    Ok(())
}
