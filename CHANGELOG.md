# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog],
and this project adheres to [Semantic Versioning].

## [Unreleased]

### Added
- Limited async MzML reader support
- mzMLb reader support

### Changed
- `RandomAccessScanIterator` methods now return mutable references to self, making them actually useful in a chain.
- Make some window size attributes smaller as they do not require double precision.
- Clean up the internal implementation of the various internal `SpectrumBuilder` types.
- Factor up `mzdata::spectrum::signal` to be less monolithic.

### Deprecated

### Removed

### Fixed

### Security


## [0.5.0] - 2021-09-22

### Added
- MzML writing via `mzdata::io::mzml::MzMLWriter`
- Added feature flags to allow the user to choose amongst more `flate2` backends (zlib *default*, zlib-ng-compat, miniz_oxide)
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
[unreleased]: https://github.com/mobiusklein/mzdata/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/mobiusklein/mzdata/compare/v0.1.0...v0.5.0
[0.1.0]: https://github.com/mobiusklein/mzdata/releases/tag/v0.1.0