"""Precursor, isolation window, selected ions, and activation for MS2 spectra."""

from __future__ import annotations


def test_ms2_has_precursor(ms2_spectrum) -> None:
    assert ms2_spectrum.precursor is not None


def test_isolation_window(ms2_spectrum) -> None:
    window = ms2_spectrum.precursor.isolation_window
    assert window.lower_bound <= window.target <= window.upper_bound
    assert window.contains(window.target)


def test_selected_ions(ms2_spectrum) -> None:
    ions = ms2_spectrum.precursor.ions()
    assert len(ions) >= 1
    ion = ions[0]
    assert ion.mz > 0
    assert ion.intensity >= 0
    assert ion.charge is None or (isinstance(ion.charge, int) and ion.charge != 0)


def test_selected_ion_near_isolation_window(ms2_spectrum) -> None:
    precursor = ms2_spectrum.precursor
    ion = precursor.ions()[0]
    window = precursor.isolation_window
    # The selected ion should lie within the isolation window (allow a small margin for windows
    # recorded only as a target width).
    assert window.contains(float(ion.mz)) or (window.lower_bound - 5.0 <= ion.mz <= window.upper_bound + 5.0)


def test_activation(ms2_spectrum) -> None:
    activation = ms2_spectrum.precursor.activation
    assert activation.energy >= 0
    assert activation.method is None or isinstance(activation.method, str)
    assert isinstance(activation.methods(), list)


def test_precursor_id_type(ms2_spectrum) -> None:
    precursor_id = ms2_spectrum.precursor.precursor_id
    assert precursor_id is None or isinstance(precursor_id, str)
