"""Shared fixtures for the pymzdata Python test suite.

The suite reads the small sample files under the repository's ``test/data`` directory and asserts on
invariants of the parsed output (types, array shapes, sort order, m/z-in-isolation-window, etc.) rather
than hard-coded numbers, so it is robust to the exact sample files. A cross-check against ``pyteomics``
(a ``[test]`` dependency) validates the numeric values against an independent reader.
"""

from __future__ import annotations

import pathlib

import pytest

import pymzdata

# This file: <repo-root>/crates/pymzdata/tests/conftest.py -> repo root is three levels up.
_REPO_ROOT = pathlib.Path(__file__).resolve().parents[3]
_DATA_DIR = _REPO_ROOT / "test" / "data"


@pytest.fixture(scope="session")
def data_dir() -> pathlib.Path:
    """The repository's sample-data directory."""
    if not _DATA_DIR.is_dir():
        pytest.skip(f"sample data directory not found: {_DATA_DIR}")
    return _DATA_DIR


@pytest.fixture(scope="session")
def mzml_path(data_dir: pathlib.Path) -> str:
    """Path to the small mzML sample."""
    path = data_dir / "small.mzML"
    if not path.exists():
        pytest.skip("small.mzML not available")
    return str(path)


@pytest.fixture
def reader(mzml_path: str):
    """An open MZReader over the small mzML sample (closed on teardown)."""
    with pymzdata.MZReader(mzml_path) as r:
        yield r


@pytest.fixture(scope="session")
def spectra(mzml_path: str) -> list:
    """Every spectrum in the small mzML sample, materialised once."""
    with pymzdata.MZReader(mzml_path) as r:
        return list(r)


@pytest.fixture(scope="session")
def ms1_spectrum(spectra: list):
    """The first MS1 spectrum."""
    for spectrum in spectra:
        if spectrum.ms_level == 1:
            return spectrum
    pytest.skip("no MS1 spectrum in the sample file")


@pytest.fixture(scope="session")
def ms2_spectrum(spectra: list):
    """The first MS2 spectrum that carries a precursor."""
    for spectrum in spectra:
        if spectrum.ms_level == 2 and spectrum.precursor is not None:
            return spectrum
    pytest.skip("no MS2 spectrum with a precursor in the sample file")
