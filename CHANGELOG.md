# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog],
and this project adheres to [Semantic Versioning].

## [0.28.1] - 2024-09-01

### Fixed

- Revert `mzsignal` change to fix downstream crashes while averaging

## [0.28.0] - 2024-08-30

### Added

- Add `PeakDataIterDispatch` to the public API
- Add more descriptive documentation for `meta`-metadata types.

### Changed

- Upgraded `mzpeaks` v0.20.0, `mzsignal` v0.21.0, `thermorawfilereader` v0.2.7

## [0.27.0] - 2024-08-25

### Added

- Add dispatching implementation to `PeakDataIter` and `RefPeakDataIter`
- Add basic `USI` parsing

### Fixed

- Do not `panic` when requested MGF  ID not in the index (#5)
- Port #5 to all other reader types

## [0.26.0] - 2024-08-10

### Added

- Added `Unit::VoltSecondPerSquareCentimeter`
- Added `BinaryArrayMap3D::ion_mobility_unit`
- Add `get` and `iter` methods for `PeakDataLevel` and `RefPeakDataLevel`
- Add tutorial for spectrum types
- Add `Sample` to the `MSDataFileMetadata`
- Add basic `usi` parsing.

### Changed

- Change `MzMLWriterType` methods to be generic to peaks.
  All the writing methods are generic over peak types
  so long as they satisfy the `BuildArrayMapFrom` constraint.
  The type still needs type parameters for default behavior.
  This may cause parameter inference issues for `RawSpectrum`
  in which case use `SpectrumWriter::write`.

### Documentation

- More docs

### Fixed

- Fix mzpeaks API breakage

## [0.25.0] - 2024-07-26

### Added

- Add `NativeSpectrumIdentifierFormatTerm` and `MassSpectrometerFileFormatTerm` enums for CV subsets
	- `NativeSpectrumIdentifierFormatTerm` corresponds to "native spectrum identifier format" and
	  includes a regular expression from the controlled vocabulary
	  for parsing ids.
	- `MassSpectrometerFileFormatTerm` corresponds to "mass spectrometer file format". No extra
	  behavior.
- Add `IonMobilityFrameGroup` plus associated iterators and traits.
- Add `BuildArrayMap3DFrom` and `BuildFromArrayMap3D` traits
- Add `IonMobilityFrameWriter` trait
- Add `read_checksum` to `MzMLReaderType`
- Add `BuildArrayMap3DFrom` and `BuildFromArrayMap3D` to the public API
- Add `IonMobilityFrameWriter` implementation for mzML and mzMLb
    - `FeatureDataLevel` and `RefFeatureDataLevel` to mirror the peak set equivalents for `IonMobilityFrameLike`
    - Implement `From<MultiLayerIonMobilityFrame>` for `RawSpectrum` to let `IonMobilityFrameWriter` re-use most of `SpectrumWriter`.

### Changed

- Change `RawSpectrum` to be `SpectrumLike` over all peak types.
- Upgrade controlled vocabulary
- Rename methods of `IonMobilityFrameWriter` to not clash with `SpectrumWriter`

### Fixed

- Fixed error in setting highest observed m/z in `SpectrumLike::update_summaries`

## [0.24.0] - 2024-07-18

### Changed

- Upgrade `mzsignal` to v0.18.0

## [0.23.0] - 2024-07-17

### Changed

- Upgraded `serde`, `rayon`

### Fixed

- Don't write duplicate CV params about MS level
- Write the correct size of the chromatogram list in mzML
- Fix up type names in the spectrum documentation

## [0.22.0] - 2024-07-16

### Added

- Add owning variants for `write_all` `SpectrumWriter`
- Add `ChromatogramSource` trait for producing `Chromatogram`s
- Add `name` and `unit` to `ByteArrayView`

### Changed

- Change `OffsetIndex` key to `Box<str>`
	This reduces memory consumption.
- Upgraded `mzsignal` minimum version

### Fixed

- Fix error in search-by-time methods

## [0.21.0] - 2024-07-01

### Added

- Add Dimensionless to `Unit` for arrays without units
- Add ion mobility stacked spectra support, initial build out of trait system.
- Add Boolean variant to `ParamValue` types

### Fixed

- Fixed #3 empty spectrum from trailing newlines in MGF

## [0.20.0] - 2024-05-30

### Fixed

- Handle unknown mass analyzers slightly more gracefully.
- Actually use the fast path when pre-encoding arrays

## [0.19.0] - 2024-05-26

### Added

- Add more convenience methods for populating `InstrumentConfiguration` and `ParamDescribed` types
- Support writing to STDOUT properly in `mzconvert` example

### Fixed

- Fixed the creation of `InstrumentConfiguration` from Thermo RAW files

## [0.18.0] - 2024-05-25

### Changed

- Upgraded `mzpeaks` minimum version to 0.12.0, `mzsignal` minimum version to 0.13.0
- Change how spectrum summary descriptions are calculated
	- Do not pre-emptively re-calculate summaries when writing to avoid decoding overhead.
	- Added `update_summaries` to `SpectrumLike` and `fetch_summaries` to `RefPeakDataLevel` to do
	  this work.
- Upgrade mzsignal minimum version to 0.14.0

### Fixed

- Fix changelog

## [0.17.0] - 2024-05-18

### Changed

- Upgraded `mzpeaks` minimum version to 0.12.0, `mzsignal` minimum version to 0.13.0
- Change how spectrum summary descriptions are calculated
	- Do not pre-emptively re-calculate summaries when writing to avoid decoding overhead.
	- Added `update_summaries` to `SpectrumLike` and `fetch_summaries` to `RefPeakDataLevel` to do
	  this work.

### Fixed

- Fix changelog

## [0.16.0] - 2024-05-07

### Added

- Added Thermo instrument methods to Thermo instrumentConfiguration
- Add `Hash` implementation for `Param`, `Value`, and associated ref types
- Add `referenceParamGroup` inference for `instrumentConfiguration` when writing mzML, and escape newlines in `Param` values

### Fixed

- Properly set compression in store_compressed
- Empty buffer coercion failure on Linux (#2)
- Prevent slice coercion from failing in coerce_mut
- Prevent potential unwrap failures in PeakData.base_peak when dealing with unorderable floats

## [0.15.0] - 2024-05-02

### Added

- Added `store_compressed` to `DataArray` to support compression at-will
- Added `encode_array` to specify compression state of a `DataArray` in `BinaryArrayMap`

### Changed

- Export the utils module for convenience

### Documentation

- Add documentation for more of the param module
- Simplify and document peak type blanket traits
- Document new compression methods

### Fixed

- MGF spectra are properly classified as MS2 and centroid
- PeakDataLevel is now exported as was originally intended

### Maintenance

- Use git-cliff

## [0.14.0] - 2024-04-25

### Added
- Added dependency on `chrono` to parse the `startTimeStamp` in mzML(b), and the equivalent Thermo property.
- Added more conversion convenience methods for parameter-like enums:
  - `MassSpectrometryFormat`
  - `FormatConversion`
- Added `mzdata::meta::custom_software_name` convenience for creating a CV parameter describing an unreleased tool.
  Helpers to deal with *published* tools will be worked on soon.
- Added `MassSpectrometryRun` to `ThermoRawReaderType`.
- Added `Value` and `ValueRef` types which are used to hold `cvParam`/`userParam` values with eager parsing. While
  not an actual space savings, this removes the need to repeatedly parse values if they are accessed more than once.
  They support coercion via the `ParamValue` trait in the prelude. `Param` and `ParamCow` also implement this trait
  and delegate it to their `value` field.
- Added `ion_mobility`, `filter_string`, and `mass_resolution` methods to `ScanEvent` that perform `cvParam` look ups.
- Added `MGFStyle` traits and marker types that control how MGF spectrum headers are populated. The default behavior
  is now encapsulated in the default style, `MZDataMGFStyle`.
- Added `SoftwareTerm` enumeration generated from the [`psi-ms.obo`](https://github.com/HUPO-PSI/psi-ms-CV).

### Changed
- The `MassSpectrometryRun.start_time` field is now a `chrono::DateTime` instead of `String`.
- `Param` and `ParamCow` now store their values as `Value`/`ValueRef`, which eagerly parse their values
  to elide the parsing cost on repeated access.
- `CURIE` can now be constructed with the `mzdata::curie` macro
- `MGFWriterType` now takes a fourth template parameter, a type implementing the `MGFStyle` trait. This
  is parameterized by default.
- The `MGFWriterType`'s public API has been expanded and the internal responsibilities refactored so that
  the `write_header` and `write_peaks` methods are more succinct.

## [0.13.0] - 2024-03-28

### Added
- `MGFReaderType` and `MGFWriterType` implement `MSDataFileMetadata`
- `ThermoRawFileReader` has been added to read Thermo RAW files when the .NET 8 runtime is available, using [`thermorawfilereader`](https://crates.io/crates/thermorawfilereader/0.2.1)
- `Source` and `Sink` algebraic types to represent things that spectra can be read from or written to.
- `mz_read` and `mz_write` are macros to open files for reading and writing in unboxed context, but which
  only live within a scoped closure.
- `MassSpectrometryReadWriteProcess` trait for orchestrating reading from a `Source`, writing to a `Sink`, and transforming the
  data through an arbitrary function specified as part of the trait implementation. Like `mz_read`/`mz_write`, the scope enclosed
  by the trait method.

### Changed
- `MGFWriterType` now generates a spectrum title when one is absent, rather than defaulting to
  the spectrum's native ID.
- `CURIE` can now be compared to `Param`
- Renamed `ScanWriter` to `SpectrumWriter` and `ScanSource` to `SpectrumSource` for consistency with other trait naming conventions.
- `MZFileReader::open_file` now returns an `io::Result` in keeping with the idea that reading a `File` might fail as well,
  even if it is already open, because it is the wrong type of file. This also allows file formats that cannot be read from
  arbitrary `io::Read` objects to signal failure without crashing the whole system.
- `Collator`, `std::sync::mpsc::{Sender, SyncSender}` now implement `SpectrumWriter` when properly parameterized.
- `PeakDataLevel` has been refactored into two types, `PeakDataLevel` is an owning type and `RefPeakDataLevel`
  is a borrowing type.

### Deprecated

### Removed

### Fixed

### Security


## [0.12.0] - 2024-01-29

### Changed
- Require a newer version of `mzsignal`, fixing the rather embarrassing error of swapping FWHM
  and SNR during peak picking.
- Thicken the use of internal abstraction around `PrecursorSelection` for the future of allowing
  more than one `SelectedIon` per `Precursor`.

## [0.11.0] - 2024-01-24

### Added
- Added an `hdf5_static` flag to build the HDF5 library from source at build time.

### Removed
- Removed `SpectrumLike` from the crate root, it is more appropriately imported with
  the rest of the prelude.


## [0.9.0] - 2024-01-20

### Added
- `RestartableGzDecoder` added to `mzdata::io` to provide an `io::Seek` compatible GZIP decompressor,
  with the caveat that seeking backwards requires re-reading everything read thus far and that seeks
  cannot be relative to the EOF.
- Added documentation to more examples.

### Fixed
- `infer_from_stream` handles compressed streams correctly without triggering an unexpected EOF error.
- `MGFWriterType` properly formats additional annotations.

## [0.8.0] - 2024-01-10

### Added
- Added `close` to the `SpectrumWriter` trait which "closes" the formatted structure of the file. As Rust lacks a notion of a "closed"
  `io::Write`, the underlying writer isn't actually "closed" until the whole struct is dropped.
- Added `Drop` implementation for `MzMLWriterType` and `MzMLbWriterType` which ensures that the `close` method is called to make the
  resulting file well-formed.
- Added new peak picking methods `MultiLayerSpectrum::pick_peaks_in_intervals` and `RawSpectrum::pick_peaks_in_intervals_into` that
  pick peaks only in selected m/z regions.
- Added `MassSpectrometryRun` to record information found on or near the `<run>` element in mzML which includes several default values
  otherwise absent from the data model. It is accessible via `MSDataFileMetadata::run_description` trait method, included in the prelude.
- Added the `MSDataFileMetadata::spectrum_count_hint` which returns the total number of spectra in the file, if available.
- Added `MSDataFileMetadata` implementations for `SpectrumIterator`, `StreamingSpectrumIterator`, and `SpectrumGroupingIterator` where
  their sources do.
- `SpectrumGroupingIterator` and other such iterator support `RandomAccessSpectrumGroupingIterator`.

### Changed
- `SpectrumWriter` no longer applies a lifespan requirement on individual writing operations.
- `filename` is no longer a required dependency, it is only needed to use `MzMLbReaderType::from_file` which otherwise
  panics. It introduces unpredictable and difficult to diagnose compilation errors.
- `MGFWriterType` skips MS1 spectra automatically.

### Fixed
- Properly track whether the `spectrumList` has been finished by mzML parsers.
- Iteration over `SpectrumGrouping` when there are missing components has been fixed.
- Binary data arrays may panic if empty on certain platforms with certain compression schemes.


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
- Massive refactoring of `mzdata::io::traits` to make more traits depend upon `SpectrumSource` instead of `SpectrumIterator` and to make things slightly less verbose.
- Switched the default `mzsignal` backend to `nalgebra` instead of `intel-mkl` for simplicity.

## [0.5.0] - 2021-09-22

### Added

- MzML writing via `mzdata::io::mzml::MzMLWriter`
- Added feature flags to allow the user to choose amongst more `flate2` backends (zlib _default_, zlib-ng-compat, miniz_oxide)
- Grouped iteration mode for connecting precursor and product spectra over an iterator stream using the `groups` method of `SpectrumSource`.

### Changed

- Re-structuring and renaming of the various iterator mechanisms for more
  consistency. `ScanIterator` -> `SpectrumIterator`, et cetera. Minor refactoring
  of this sort expected to come for `SpectrumSource` as responsibilities are worked out.

### Deprecated

### Removed

### Fixed

- Fixed documentation in several places, particularly where it was substantially out of date.

### Security

<!-- Links -->

[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html

<!-- Versions -->

[unreleased]: https://github.com/mobiusklein/mzdata/compare/v0.28.1...HEAD
[0.28.1]: https://github.com/mobiusklein/mzdata/compare/v0.28.0...v0.28.1
[0.28.0]: https://github.com/mobiusklein/mzdata/compare/v0.27.0...v0.28.0
[0.27.0]: https://github.com/mobiusklein/mzdata/compare/v0.26.0...v0.27.0
[0.26.0]: https://github.com/mobiusklein/mzdata/compare/v0.25.0...v0.26.0
[0.25.0]: https://github.com/mobiusklein/mzdata/compare/v0.24.0...v0.25.0
[0.24.0]: https://github.com/mobiusklein/mzdata/compare/v0.23.0...v0.24.0
[0.23.0]: https://github.com/mobiusklein/mzdata/compare/v0.22.0...v0.23.0
[0.22.0]: https://github.com/mobiusklein/mzdata/compare/v0.21.0...v0.22.0
[0.21.0]: https://github.com/mobiusklein/mzdata/compare/v0.20.0...v0.21.0
[0.20.0]: https://github.com/mobiusklein/mzdata/compare/v0.19.0...v0.20.0
[0.19.0]: https://github.com/mobiusklein/mzdata/compare/v0.18.0...v0.19.0
[0.18.0]: https://github.com/mobiusklein/mzdata/compare/v0.17.0...v0.18.0
[0.17.0]: https://github.com/mobiusklein/mzdata/compare/v0.16.0...v0.17.0
[0.16.0]: https://github.com/mobiusklein/mzdata/compare/v0.15.0...v0.16.0
[0.15.0]: https://github.com/mobiusklein/mzdata/compare/v0.14.0...v0.15.0
[0.14.0]: https://github.com/mobiusklein/mzdata/compare/v0.13.0...v0.14.0
[0.13.0]: https://github.com/mobiusklein/mzdata/compare/v0.12.0...v0.13.0
[0.12.0]: https://github.com/mobiusklein/mzdata/compare/v0.11.0...v0.12.0
[0.11.0]: https://github.com/mobiusklein/mzdata/compare/v0.10.0...v0.11.0
[0.10.0]: https://github.com/mobiusklein/mzdata/compare/v0.9.0...v0.10.0
[0.9.0]: https://github.com/mobiusklein/mzdata/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/mobiusklein/mzdata/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/mobiusklein/mzdata/compare/v0.5.0...v0.7.0
[0.5.0]: https://github.com/mobiusklein/mzdata/compare/v0.1.0...v0.5.0
[0.1.0]: https://github.com/mobiusklein/mzdata/releases/tag/v0.1.0