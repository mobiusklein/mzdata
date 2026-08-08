# pymzdata

`pymzdata` provides Python bindings for [mzdata](https://github.com/mobiusklein/mzdata),
a fast Rust library for reading mass spectrometry data files.

## Supported formats

`MZReader` supports mzML and indexed mzML, MGF, Bruker TDF, and imzML files.
`IMMZReader` provides ion-mobility frame access for compatible mzML and Bruker TDF
files. Thermo RAW is also supported when its native runtime dependencies are
available on the host system.

## Installation

```console
pip install pymzdata
```

Python 3.9 or later is required. NumPy is installed automatically and is used for
the peak arrays returned by the reader.

## Quickstart

```python
import pymzdata

with pymzdata.MZReader("example.mzML") as reader:
    print(f"{len(reader)} spectra")

    for spectrum in reader:
        print(spectrum.id, spectrum.ms_level, spectrum.start_time)
        mz = spectrum.mz_array()
        intensity = spectrum.intensity_array()
        if mz is not None:
            print(f"  {len(mz)} peaks")
```

`MZReader` supports iteration and random access with `get_by_index`, `get_by_id`,
and `get_by_time` (retention time in minutes). Spectrum metadata, precursor
information, and decoded NumPy peak arrays are available from each `Spectrum`.

For Bruker ion-mobility data, use `IMMZReader` directly or convert an open
`MZReader` with `into_frame_reader()` to access ion-mobility frames.

## Development

To build the bindings from a checkout, install [maturin](https://www.maturin.rs/)
and run:

```console
cd crates/pymzdata
maturin develop
```
