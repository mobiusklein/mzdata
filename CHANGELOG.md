# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog],
and this project adheres to [Semantic Versioning].

## [0.7.0] - 2023-12-25

### Added

- Limited async mzML reader support
- mzMLb read and write support
- Reading mzML and MGF from `STDIN`. HDF5, and ergo mzMLb is not supported on non-seek-able I/O devices. See the `from_stdin` example.
- Parsing of mzML `run` and `spectrumList` metadata, although they are still not a part of the common data model
- Spectrum averaging now has eager `averaging` and `averaging_deferred` adapter implementations as iterator adapters on `SpectrumGroupIterator`.
  The deferred adapter is preferred for distributing the process with `rayon`. See the `averaging_writer` example.
- Added ordered parallel iteration collation with `mzdata::spectrum::Collator` to make consuming `rayon` iterators easier while preserving the
  original order. See `averaging_writer` example.
- The mzML and mzMLb writers now write the total ion chromatogram and base peak chromatogram

### Changed

- `RandomAccessScanIterator` methods now return mutable references to self, making them actually useful in a chain.
- Make some window size attributes smaller as they do not require double precision.
- Clean up the internal implementation of the various internal `SpectrumBuilder` types.
- Factor up `mzdata::spectrum::signal` to be less monolithic and a complete redesign of the traits used to convert `mzpeaks` to and from binary arrays.
- Massive refactoring of `mzdata::io::traits` to make more traits depend upon `ScanSource` instead of `SpectrumIterator` and to make things slightly less verbose.
- Switched the default `mzsignal` backend to `nalgebra` instead of `intel-mkl` for simplicity.

## [0.5.0] - 2021-09-22

### Added

- MzML writing via `mzdata::io::mzml::MzMLWriter`
- Added feature flags to allow the user to choose amongst more `flate2` backends (zlib _default_, zlib-ng-compat, miniz_oxide)
- Grouped iteration mode for connecting precursor and product spectra over an iterator stream using the `groups` method of `ScanSource`.

### Changed

- Re-structuring and renaming of the various iterator mechanisms for more
  consistency. `ScanIterator` -> `SpectrumIterator`, et cetera. Minor refactoring
  of this sort expected to come for `ScanSource` as responsibilities are worked out.

### Deprecated

### Removed

### Fixed

- Fixed documentation in several places, particularly where it was substantially out of date.

### Security

<!-- Links -->

[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html

<!-- Versions -->

[unreleased]: https://github.com/mobiusklein/mzdata/compare/v0.7.0...HEAD
[0.7.0]: https://github.com/mobiusklein/mzdata/compare/v0.5.0...v0.7.0
[0.5.0]: https://github.com/mobiusklein/mzdata/compare/v0.1.0...v0.5.0
[0.1.0]: https://github.com/mobiusklein/mzdata/releases/tag/v0.1.0
