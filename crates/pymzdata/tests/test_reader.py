"""Reader lifecycle, iteration, random access, detail level, and error handling."""

from __future__ import annotations

import pytest

import pymzdata


def test_open_and_close(mzml_path: str) -> None:
    r = pymzdata.MZReader(mzml_path)
    assert not r.closed()
    assert len(r) > 0
    r.close()
    assert r.closed()


def test_context_manager_closes(mzml_path: str) -> None:
    with pymzdata.MZReader(mzml_path) as r:
        assert not r.closed()
        n = len(r)
    assert r.closed()
    assert n > 0


def test_iteration_count_matches_len(mzml_path: str) -> None:
    with pymzdata.MZReader(mzml_path) as r:
        expected = len(r)
        count = sum(1 for _ in r)
    assert count == expected


def test_repr(reader) -> None:
    assert "MZReader" in repr(reader)
    assert "open" in repr(reader)


def test_get_by_index(mzml_path: str) -> None:
    with pymzdata.MZReader(mzml_path) as r:
        first = r.get_by_index(0)
        assert first is not None
        assert first.index == 0
        # Out-of-range index returns None rather than raising.
        assert r.get_by_index(10**9) is None


def test_get_by_id_roundtrip(mzml_path: str) -> None:
    with pymzdata.MZReader(mzml_path) as r:
        first = r.get_by_index(0)
        again = r.get_by_id(first.id)
        assert again is not None
        assert again.id == first.id
        assert again.index == first.index


def test_get_by_id_missing(mzml_path: str) -> None:
    with pymzdata.MZReader(mzml_path) as r:
        assert r.get_by_id("controllerType=0 controllerNumber=1 scan=999999999") is None


def test_get_by_time(mzml_path: str, spectra: list) -> None:
    target = spectra[len(spectra) // 2].start_time
    with pymzdata.MZReader(mzml_path) as r:
        found = r.get_by_time(target)
    assert found is not None
    # The returned spectrum's time is at least as close as the extremes of the run.
    span = max(abs(spectra[0].start_time - target), abs(spectra[-1].start_time - target))
    assert abs(found.start_time - target) <= span + 1e-6


def test_file_not_found() -> None:
    with pytest.raises((IOError, OSError)):
        pymzdata.MZReader("/no/such/file/definitely_missing.mzML")


def test_operations_on_closed_reader(mzml_path: str) -> None:
    r = pymzdata.MZReader(mzml_path)
    r.close()
    with pytest.raises(RuntimeError):
        r.get_by_index(0)
    with pytest.raises(RuntimeError):
        r.file_description()


def test_detail_level_getter_and_roundtrip(reader) -> None:
    level = reader.detail_level
    assert "DetailLevel" in repr(level)
    # Round-tripping the current level must not raise.
    reader.detail_level = level
    assert "DetailLevel" in repr(reader.detail_level)
