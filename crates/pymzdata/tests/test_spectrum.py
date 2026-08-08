"""Spectrum fields, peak arrays, params, and scan events."""

from __future__ import annotations

import numpy as np
import pytest

import pymzdata


def test_basic_fields(spectra: list) -> None:
    for spectrum in spectra:
        assert isinstance(spectrum.id, str) and spectrum.id
        assert isinstance(spectrum.index, int) and spectrum.index >= 0
        assert isinstance(spectrum.ms_level, int) and spectrum.ms_level >= 1
        assert isinstance(spectrum.start_time, float)


def test_index_matches_iteration_order(spectra: list) -> None:
    for i, spectrum in enumerate(spectra):
        assert spectrum.index == i


def test_signal_continuity_flags(spectra: list) -> None:
    profile = pymzdata.SignalContinuity.Profile
    centroid = pymzdata.SignalContinuity.Centroid
    for spectrum in spectra:
        sc = spectrum.signal_continuity
        assert sc in (pymzdata.SignalContinuity.Unknown, centroid, profile)
        assert spectrum.is_profile == (sc == profile)
        assert spectrum.is_centroid == (sc == centroid)


def test_mz_and_intensity_arrays(spectra: list) -> None:
    saw_arrays = False
    for spectrum in spectra:
        mz = spectrum.mz_array()
        intensity = spectrum.intensity_array()
        if mz is None:
            continue
        saw_arrays = True
        assert isinstance(mz, np.ndarray) and mz.dtype == np.float64
        assert isinstance(intensity, np.ndarray) and intensity.dtype == np.float32
        assert mz.shape == intensity.shape
        assert len(mz) == len(spectrum)
        assert np.all(np.diff(mz) >= 0), "m/z array must be sorted ascending"
        assert np.all(intensity >= 0), "intensities must be non-negative"
    assert saw_arrays, "no spectrum exposed an m/z array"


def test_len_matches_peak_count(spectra: list) -> None:
    for spectrum in spectra:
        mz = spectrum.mz_array()
        assert len(spectrum) == (0 if mz is None else len(mz))


def test_raw_arrays_contains_mz(ms1_spectrum) -> None:
    raw = ms1_spectrum.raw_arrays()
    assert isinstance(raw, dict)
    joined = " ".join(raw.keys()).lower()
    assert "m/z" in joined or "mz" in joined
    for value in raw.values():
        assert isinstance(value, np.ndarray)


def test_centroid_peaks_shape(spectra: list) -> None:
    for spectrum in spectra:
        peaks = spectrum.centroid_peaks()
        if peaks is None:
            continue
        mz, intensity = peaks
        assert isinstance(mz, np.ndarray) and isinstance(intensity, np.ndarray)
        assert mz.shape == intensity.shape


def test_scan_events(ms1_spectrum) -> None:
    events = ms1_spectrum.scan_events()
    assert isinstance(events, list)
    for event in events:
        assert isinstance(event.start_time, float)
        assert event.injection_time >= 0
        filter_string = event.filter_string()
        assert filter_string is None or isinstance(filter_string, str)
        for window in event.scan_windows():
            assert window.lower_bound <= window.upper_bound
            assert window.contains((window.lower_bound + window.upper_bound) / 2)


def test_repr(ms1_spectrum) -> None:
    assert "Spectrum(" in repr(ms1_spectrum)


def test_to_proxi(ms1_spectrum) -> None:
    proxi = ms1_spectrum.to_proxi()
    assert isinstance(proxi, dict)
