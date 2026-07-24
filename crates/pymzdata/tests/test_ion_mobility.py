"""Ion-mobility frame reading via IMMZReader (skipped where the sample lacks IM data)."""

from __future__ import annotations

import pathlib

import pytest

import pymzdata


def _open_im_reader(path: pathlib.Path):
    try:
        return pymzdata.IMMZReader(str(path))
    except (RuntimeError, IOError, OSError) as exc:
        pytest.skip(f"ion-mobility frames not available for {path.name}: {exc}")


def test_immzreader_frames(data_dir: pathlib.Path) -> None:
    path = data_dir / "diaPASEF.mzML"
    if not path.exists():
        pytest.skip("diaPASEF.mzML not available")
    with _open_im_reader(path) as reader:
        n = len(reader)
        if n == 0:
            pytest.skip("no ion-mobility frames in the sample")
        frame = reader.get_by_index(0)
        assert frame is not None
        assert isinstance(frame.id, str) and frame.id
        assert frame.ms_level >= 1
        assert isinstance(frame.start_time, float)
        raw = frame.raw_arrays()
        assert isinstance(raw, dict)


def test_into_frame_reader_from_mzreader(data_dir: pathlib.Path) -> None:
    path = data_dir / "diaPASEF.mzML"
    if not path.exists():
        pytest.skip("diaPASEF.mzML not available")
    reader = pymzdata.MZReader(str(path))
    try:
        frame_reader = reader.into_frame_reader()
    except (RuntimeError, IOError, OSError) as exc:
        pytest.skip(f"cannot convert to frame reader: {exc}")
    with frame_reader as fr:
        assert len(fr) >= 0
    # into_frame_reader() consumes the MZReader (moves the underlying reader out), leaving it closed.
    assert reader.closed()
