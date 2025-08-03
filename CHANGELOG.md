# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog],
and this project adheres to [Semantic Versioning].

## [0.57.0] - 2025-08-03

### Added

- Add `get_data_arrays_for` to `ThermoRawReaderType` for loading data arrays without the additional metadata
- Add `iter_buffer` to `StreamingSpectrumIterator` to borrow the spectrum buffer

## [0.56.0] - 2025-07-20

### Fixed

- Fix bug in `BinaryArrayMap3D::unstack`

## [0.55.0] - 2025-07-17

### Added

- Add `ChromatogramLike` to the prelude
- Add `ArrayType::from_accession` and `BinaryDataArrayType::from_accession` to convert from CURIEs
- Add ISQ7000, Orbitrap Exploris GC 240, support multiple FAIMS voltages (#26)
* Add Orbitrap Exploris GC 240, to supported instruments. Added some debugging statements to make diagnosing missing instrument types easier in the future

* Add ISQ 7000, fix typo in reverse_transpose

* remove unknown detector type, bump number of instruments

* support multiple CompensationVoltageValues

* Add handling of FAIMS start and end ramp, add warning if there are more than 2 faims voltages

* Update Cargo.toml

---------
- Add support for `ScanSettings` section; fix bug in decompression of Numpress

### Fixed

- Fix recognition of wavelength array in mzML

## [0.53.0] - 2025-05-26

### Added

- Add prototype Zstandard, Numpress, and dictionary-based binary data array compression methods (#24)
* add more `MS-Numpress` encodings
* feature: add units to common arrays
* change: change TDF peak collapsing to be off by default
* doc: document the behavior of `TDFSpectrumReaderType`
* feature: add byte shuffling to the dictionary encoding
* chore: pin Rust toolchain version in test to 1.86 due to hostfxr-sys + syn issue

## [0.52.0] - 2025-04-08

### Added

- Add temperature and pressure units
- Avoid attempting to serialize a parameter value that is a buffer.

### Changed

- Change the order of merges in TDF peak flattening

### Fixed

- Change attribute newline detection to use `memchr`
- Do not try to marshal peaks if peaks already exist with `MultiLayerSpectrum::try_build_peaks`
- Do not include a tag when serializing `PROXIAccession`
- Fix comparison in `HasIonMobility`

### Removed

- Remove newline escaping in mzML attribute values, trust the XML library to escape what needs to be

## [0.51.0] - 2025-03-27

### Added

- Add temperature and pressure units
- Avoid attempting to serialize a parameter value that is a buffer.

### Fixed

- Change attribute newline detection to use `memchr`

### Removed

- Remove newline escaping in mzML attribute values, trust the XML library to escape what needs to be

## [0.50.0] - 2025-03-20

### Added

- Add a more efficient implementation of `nth` to some `SpectrumSource`-based iterators and add `get_group_by_index` to `SpectrumSource`
- Add `MZReaderType::open_gzipped_read_seek` to conveniently handle gzipped files
- Add `MZReaderType::open_gzipped_read` to read non-seek-able gzipped streams. Less efficient buffering but does the job.

### Changed

- Change `BinaryArrayMap3D`'s serialization to serialize IM-array map pairs in an array instead of a map, requiring `serde_with`

## [0.49.0] - 2025-03-12

### Added

- Add basic tests for Bruker TDF file reading

### Changed

- Changed `MZReaderType::Unknown` to require `Sync`

## [0.48.3] - 2025-03-12

### Added

- Add BrukerTDF support to `MZReaderType::open_path`

## [0.48.2] - 2025-02-20

### Added

- Add proper handling of TDF DIA windows

## [0.48.1] - 2025-02-19

### Fixed

- Fix reading of DIA TDF files

## [0.48.0] - 2025-02-18

### Added

- Add ion mobility unit propagation to `IonMobilityFrameLike` types

### Fixed

- Upgrade minimum `mzsignal` version to fix malformed mz peak picking numerical errors
- Formally use the expanded ion mobility `ArrayType`s, can't solve units

### Removed

- Remove the requirement for `Default` on peak types throughout the library

## [0.47.0] - 2025-02-15

### Added

- Add chromatogram data reader to `bruker_tdf`

### Changed

- `Generic3DIonMobilityFrameSource::next` will skip spectra without ion mobility dimension

## [0.46.1] - 2025-02-11

### Added

- Add automatic item count caching prior to compressing a `DataArray`

### Documentation

- Improve documentation about file format-related features

## [0.46.0] - 2025-02-09

### Fixed

- Fix `mz_read` and `mz_write` macro feature leakage

## [0.45.0] - 2025-02-08

### Changed

- Refactor to allow completely disabling parser machinery and create proxy crate `mzdata-spectra` to re-export the remaining symbols (#21)

## [0.44.0] - 2025-02-02

### Added

- Add `BuildArrayMap3DFrom` and `BuildFromArrayMap3D` to the prelude

### Fixed

- Fix first instrument configuration to in TDF to count ID from 0
- Prevent TDF reader index-out-of-bounds slice in `FrameToArrayMapper`

## [0.43.0] - 2025-02-01

### Added

- Add `IntoIonMobilityFrameSource` for `StreamingSpectrumIterator`

### Fixed

- Fix error in `is_tdf`

## [0.42.0] - 2025-01-26

### Added

- Add `AsyncThermoRawReaderType` and `AsyncMZReaderType`, with associated supporting types
- Add support for more `ArrayType`s
- Add optional arrays when available to Thermo data
- Add `MultiLayerIonMobilityFrame::try_build_features`
- Add `FileMetadataConfig` to consolidate `MSDataFileMetadata` details into a single type that's easier to compose.
- Async thermo import in doc-only mode
- Add ion mobility frame dispatch via `IntoIonMobilityFrameSource` and `IMMZReaderType`

### Fixed

- Fix off-by-one error in feature conversion

## [0.41.0] - 2025-01-02

### Added

- Add `AsyncSpectrumSource` trait and add `async_partial` for WASM-compatibility
- Add feature-gated `serde` support for `mzdata` types

### Fixed

- Make USI provenance identifer parsing more robust and abstract over repository codes
* Fixed parsing bug for USI continaining a 'shielded' colon.

* Changed USI parsing of provinances to use hardcoded provinance identifiers

* Moved to rsplit_once

* Handled provenance repositories better

## [0.40.0] - 2024-12-15

### Added

- Add `AsyncMGFReaderType`

### Fixed

- Make mzML parser tolerate some misformatted XML attributes better

## [0.39.0] - 2024-11-30

### Added

- Add EFO, OBI, HANCESTRO and BFO controlled vocabularies
- Add NCIT, BTO, and PRIDE to `ControlledVocabulary`

### Fixed

- Incrementally move towards decoupling `CURIE` from `u32`
- Fix `FileMetadataBuilder` consumption to set mzML spectrum count hint

## [0.38.0] - 2024-11-30

### Added

- Add `set_spectrum_count_hint` to `MSDataFileMetadata`

## [0.37.0] - 2024-11-29

### Added

- Add `detail_level` manipulation to `SpectrumSource` and equivalents
- Add `has_ion_mobility_dimension` to `SpectrumLike`
- MSV collections appear to be unsupported by PROXI backends, and the MassiVE server is broken
- Add Bruker TDF support
- Add more documentation and tutorial

### Changed

- MzML(b) readers attempt to eagerly instantiate peak lists with `try_build_peaks`

### Fixed

- Include `source_file_name` to `delegate_impl_metadata_trait`
- Thermo and MZReaderType time random access fixes
- Fix parameter conversion for mzML in `MassSpectrometryFormat::as_param`
- Fix typo in `ChromatogramLike::chromatogram_type`
- Fix FWHM assignment during reprofiling in `average_spectra` (GH #17)

## [0.36.0] - 2024-11-11

### Added

- Adjust test using `extract_features_simple` for new `mzsignal` behavior

### Fixed

- Upgrade to `mzpeaks v0.23.0`  and `mzsignal v0.26.0` and organize dependencies
- Split the `thermo` implementation so it does not require a .NET library during documentation building
- Refactor scan unit grouping trait locations, preserve import paths
- Upgrade to `mzsignal v0.27.0`

## [0.35.0] - 2024-11-01

### Added

- Add `ParamBuilder` to make building `Param` easier
- Add `unit_mut` to `ByteArrayViewMut`
- Add error handling logic to nativeID format parsing logic
- Add `MZReaderBuilder` type to the public API
- Add more metadata types to the public API
- Add `spectrum_reference` to `ScanEvent`
- Add source-manipulation functions to `Generic3DIonMobilityFrameSource`
- Added method to retrieve the raw spectrum from any USI (#15)
    * Added `USI::download_spectrum_blocking` and `USI::download_spectrum_async` to retrieve the designated spectrum from a PROXI server

### Changed

- Change the hashing function for `mzdata::spectrum::Collator` to `IdentityHasher`

### Fixed

- Fix multiple errors in `BinaryArrayMap3D`

## [0.35.0] - 2024-11-01

### Added

- Add `ParamBuilder` to make building `Param` easier
- Add `unit_mut` to `ByteArrayViewMut`
- Add error handling logic to nativeID format parsing logic
- Add `MZReaderBuilder` type to the public API
- Add more metadata types to the public API
- Add `spectrum_reference` to `ScanEvent`
- Add source-manipulation functions to `Generic3DIonMobilityFrameSource`
- Added method to retrieve the raw spectrum from any USI (#15)
    * Added `USI::download_spectrum_blocking` and `USI::download_spectrum_async` to retrieve the designated spectrum from a PROXI server

### Changed

- Change the hashing function for `mzdata::spectrum::Collator` to `IdentityHasher`

### Fixed

- Fix multiple errors in `BinaryArrayMap3D`

## [0.34.0] - 2024-10-31

### Fixed

- Fix parsing precursor charge in MGF files (#13)
* fix parsing precursor charge in MGF files

  This fixes two issues:
  - the charge value in the PEPMASS header can have the same suffixed
    format as the CHARGE header
  - the order of the headers isn't determined, so the charge value should
    still be saved correctly if the CHARGE header is before the PEPMASS
    header

## [0.33.0] - 2024-10-17

### Fixed

- Fixed crash while mgf parsing (#12)
* Fixed rare crash while mgf parsing
* Fixed silent error that returns an unfinished spectrum
* Handle charge when the sign is at the tail

## [0.32.0] - 2024-10-14

### Added

- Add some Thermo scan trailer values to `ThermoRawReaderType`
- Add `NativeSpectrumIdentifierFormatTerm::format`

### Changed

- Upgraded `mzpeaks` and `mzsignal` dependencies

## [0.31.0] - 2024-10-03

### Changed

- Upgrade to `thermorawfilereader` v0.3.0

### Fixed

- Test on thermo feature instead of dependency (#10)
- Fix `mass_charge_ratio` order of operations error

## [0.30.0] - 2024-09-23

### Added

- Add `MZReaderBuilder` to make configuring `MZReader`s easier
- Add CV-derived dissociation methods and energies, update vendored PSI-MS CV
- Add more descriptive failure messages when parsing MGF files
- Add more descriptive errors and failure messages when encountering errors while parsing mzML files

### Fixed

- Publicly re-export `mzpeaks` dependency to guarantee `mz_read` works (#9)
Prevent downstream users needing to declare a dependency on mzpeaks when
using mz_read macro. This also prevents potential version mismatches.

## [0.29.0] - 2024-09-07

### Added

- Add PROXI format parser/serializer in `mzdata:io:proxi`, no network request feature.

### Changed

- Upgrade `thermorawfilereader` to v0.2.9, theoretically making an absent .NET framework to be non-fatal
- Upgrade `mzsignal` and enable new AVX feature

### Fixed

- Change `average_spectra` to take an iterator over `SpectrumLike` (#6)
* Allowed average from an iterator and added more documentation

---------

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

[unreleased]: https://github.com/mobiusklein/mzdata/compare/v0.57.0...HEAD
[0.57.0]: https://github.com/mobiusklein/mzdata/compare/v0.56.0...v0.57.0
[0.56.0]: https://github.com/mobiusklein/mzdata/compare/v0.55.0...v0.56.0
[0.55.0]: https://github.com/mobiusklein/mzdata/compare/v0.54.0...v0.55.0
[0.54.0]: https://github.com/mobiusklein/mzdata/compare/v0.53.0...v0.54.0
[0.53.0]: https://github.com/mobiusklein/mzdata/compare/v0.52.0...v0.53.0
[0.52.0]: https://github.com/mobiusklein/mzdata/compare/v0.51.0...v0.52.0
[0.51.0]: https://github.com/mobiusklein/mzdata/compare/v0.50.0...v0.51.0
[0.50.0]: https://github.com/mobiusklein/mzdata/compare/v0.49.0...v0.50.0
[0.49.0]: https://github.com/mobiusklein/mzdata/compare/v0.48.3...v0.49.0
[0.48.3]: https://github.com/mobiusklein/mzdata/compare/v0.48.2...v0.48.3
[0.48.2]: https://github.com/mobiusklein/mzdata/compare/v0.48.1...v0.48.2
[0.48.1]: https://github.com/mobiusklein/mzdata/compare/v0.48.0...v0.48.1
[0.48.0]: https://github.com/mobiusklein/mzdata/compare/v0.47.0...v0.48.0
[0.47.0]: https://github.com/mobiusklein/mzdata/compare/v0.46.1...v0.47.0
[0.46.1]: https://github.com/mobiusklein/mzdata/compare/v0.46.0...v0.46.1
[0.46.0]: https://github.com/mobiusklein/mzdata/compare/v0.45.0...v0.46.0
[0.45.0]: https://github.com/mobiusklein/mzdata/compare/v0.44.0...v0.45.0
[0.44.0]: https://github.com/mobiusklein/mzdata/compare/v0.43.0...v0.44.0
[0.43.0]: https://github.com/mobiusklein/mzdata/compare/v0.42.0...v0.43.0
[0.42.0]: https://github.com/mobiusklein/mzdata/compare/v0.41.0...v0.42.0
[0.41.0]: https://github.com/mobiusklein/mzdata/compare/v0.40.0...v0.41.0
[0.40.0]: https://github.com/mobiusklein/mzdata/compare/v0.39.0...v0.40.0
[0.39.0]: https://github.com/mobiusklein/mzdata/compare/v0.38.0...v0.39.0
[0.38.0]: https://github.com/mobiusklein/mzdata/compare/v0.37.0...v0.38.0
[0.37.0]: https://github.com/mobiusklein/mzdata/compare/v0.36.0...v0.37.0
[0.36.0]: https://github.com/mobiusklein/mzdata/compare/v0.35.0...v0.36.0
[0.35.0]: https://github.com/mobiusklein/mzdata/compare/v0.35.0...v0.35.0
[0.35.0]: https://github.com/mobiusklein/mzdata/compare/v0.34.0...v0.35.0
[0.34.0]: https://github.com/mobiusklein/mzdata/compare/v0.33.0...v0.34.0
[0.33.0]: https://github.com/mobiusklein/mzdata/compare/v0.32.0...v0.33.0
[0.32.0]: https://github.com/mobiusklein/mzdata/compare/v0.31.0...v0.32.0
[0.31.0]: https://github.com/mobiusklein/mzdata/compare/v0.30.0...v0.31.0
[0.30.0]: https://github.com/mobiusklein/mzdata/compare/v0.29.0...v0.30.0
[0.29.0]: https://github.com/mobiusklein/mzdata/compare/v0.28.1...v0.29.0
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