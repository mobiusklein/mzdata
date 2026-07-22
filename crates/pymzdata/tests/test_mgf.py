"""Reading MGF files through the same MZReader interface."""

from __future__ import annotations

import pathlib

import pytest

import pymzdata


@pytest.fixture(scope="module")
def mgf_spectra(data_dir: pathlib.Path) -> list:
    path = data_dir / "small.mgf"
    if not path.exists():
        pytest.skip("small.mgf not available")
    with pymzdata.MZReader(str(path)) as reader:
        return list(reader)


def test_mgf_reads_spectra(mgf_spectra: list) -> None:
    assert len(mgf_spectra) > 0


def test_mgf_spectra_are_msn_with_peaks(mgf_spectra: list) -> None:
    for spectrum in mgf_spectra[:25]:
        # MGF holds fragmentation spectra: MS level > 1, a precursor, and a non-empty peak list.
        assert spectrum.ms_level >= 2
        assert spectrum.precursor is not None
        # MGF peaks are decoded into the centroid peak list (not the raw binary arrays), so
        # mz_array() is None here; the peaks come through centroid_peaks().
        peaks = spectrum.centroid_peaks()
        assert peaks is not None
        mz, intensity = peaks
        assert len(mz) == len(intensity) > 0
        assert len(spectrum) == len(mz)


def test_mgf_precursor(mgf_spectra: list) -> None:
    ions = mgf_spectra[0].precursor.ions()
    assert len(ions) >= 1
    assert ions[0].mz > 0
