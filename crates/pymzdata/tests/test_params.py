"""Controlled-vocabulary Param access and lookup."""

from __future__ import annotations

import pytest


def test_param_fields(ms1_spectrum) -> None:
    for param in ms1_spectrum.params():
        assert isinstance(param.name, str)
        assert param.id is None or isinstance(param.id, str)
        assert isinstance(param.unit, str)
        assert param.value is None or isinstance(param.value, (str, float, int, bool, list))


def test_find_param_requires_an_argument(ms1_spectrum) -> None:
    with pytest.raises(TypeError):
        ms1_spectrum.find_param()


def test_find_param_by_name(ms1_spectrum) -> None:
    params = ms1_spectrum.params()
    if not params:
        pytest.skip("spectrum exposes no params")
    name = params[0].name
    found = ms1_spectrum.find_param(name=name)
    assert found is not None
    assert found.name == name


def test_find_param_by_accession(ms1_spectrum) -> None:
    with_accession = [p for p in ms1_spectrum.params() if p.id]
    if not with_accession:
        pytest.skip("no param with an accession")
    accession = with_accession[0].id
    found = ms1_spectrum.find_param(accession=accession)
    assert found is not None
    assert found.id == accession


def test_find_param_missing_name_returns_none(ms1_spectrum) -> None:
    assert ms1_spectrum.find_param(name="this parameter does not exist") is None


def test_find_param_invalid_accession_raises(ms1_spectrum) -> None:
    with pytest.raises(ValueError):
        ms1_spectrum.find_param(accession="not-a-valid-curie")
