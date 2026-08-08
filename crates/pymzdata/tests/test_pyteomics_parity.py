"""Cross-check: pymzdata and pyteomics must parse the same mzML identically.

pyteomics is part of the ``[test]`` extra, so this runs as a normal part of the suite. It validates
pymzdata's decoded output against an independent, widely-used reader: peak arrays, MS level, retention
time (unit-normalised to minutes), and the precursor selected ion. The ``importorskip`` below is a
defensive fallback so the rest of the suite still runs on a platform where pyteomics (or its lxml /
psims dependencies) fails to install.
"""

from __future__ import annotations

import numpy as np
import pytest

import pymzdata

pyteomics_mzml = pytest.importorskip("pyteomics.mzml", reason="pyteomics not installed")


def _rt_minutes(scan_start_time) -> float:
    """Normalise a pyteomics 'scan start time' to minutes (pymzdata reports minutes)."""
    unit = getattr(scan_start_time, "unit_info", None)
    value = float(scan_start_time)
    return value / 60.0 if unit == "second" else value


@pytest.fixture(scope="module")
def pyteomics_by_id(mzml_path: str) -> dict:
    """All pyteomics spectra keyed by native id."""
    out = {}
    with pyteomics_mzml.MzML(mzml_path) as handle:
        for spec in handle:
            out[spec["id"]] = spec
    return out


def test_same_spectrum_count(pyteomics_by_id: dict, spectra: list) -> None:
    assert len(spectra) == len(pyteomics_by_id)


def test_ms_level_and_peaks_match(pyteomics_by_id: dict, spectra: list) -> None:
    compared = 0
    for spectrum in spectra:
        reference = pyteomics_by_id.get(spectrum.id)
        assert reference is not None, f"pyteomics has no spectrum with id {spectrum.id!r}"

        assert spectrum.ms_level == int(reference["ms level"])

        mz = spectrum.mz_array()
        if mz is None or "m/z array" not in reference:
            continue
        ref_mz = np.asarray(reference["m/z array"], dtype=np.float64)
        ref_intensity = np.asarray(reference["intensity array"], dtype=np.float64)
        intensity = spectrum.intensity_array()

        assert len(mz) == len(ref_mz), f"m/z length mismatch for {spectrum.id!r}"
        np.testing.assert_allclose(mz, ref_mz, rtol=1e-6, atol=1e-4)
        # intensities round-trip through f32 in pymzdata, so allow f32-scale tolerance.
        np.testing.assert_allclose(intensity, ref_intensity, rtol=1e-3, atol=1e-2)
        compared += 1
    assert compared > 0, "no spectra with peak arrays were compared"


def test_retention_time_matches(pyteomics_by_id: dict, spectra: list) -> None:
    for spectrum in spectra:
        reference = pyteomics_by_id[spectrum.id]
        scans = reference.get("scanList", {}).get("scan", [])
        if not scans or "scan start time" not in scans[0]:
            continue
        assert spectrum.start_time == pytest.approx(_rt_minutes(scans[0]["scan start time"]), abs=1e-3)


def test_precursor_selected_ion_matches(pyteomics_by_id: dict, spectra: list) -> None:
    compared = 0
    for spectrum in spectra:
        if spectrum.ms_level < 2 or spectrum.precursor is None:
            continue
        reference = pyteomics_by_id[spectrum.id]
        precursors = reference.get("precursorList", {}).get("precursor", [])
        if not precursors:
            continue
        ref_ions = precursors[0].get("selectedIonList", {}).get("selectedIon", [])
        if not ref_ions:
            continue
        ref_mz = float(ref_ions[0]["selected ion m/z"])
        ion = spectrum.precursor.ions()[0]
        assert ion.mz == pytest.approx(ref_mz, rel=1e-6, abs=1e-4)
        if "charge state" in ref_ions[0] and ion.charge is not None:
            assert ion.charge == int(ref_ions[0]["charge state"])
        compared += 1
    if compared == 0:
        pytest.skip("no MS2 precursors to compare")
