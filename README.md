# mzdata

A Rust library for reading mass spectrometry data file formats.

### Supported Formats
1. mzML and indexedmzML
2. MGF

## TODO
    1. mzML Implementation
        1. Interpret file-level metadata (file description, software, instrument configuration,
           data processing).
        2. Improve indexing performance, with the ability to store the computed index.
        3. Improve unit handling.
    2. Iteration behavior
       1. Implement a batch iterator to group together related MS1 and MSN spectra.
       2. Implement a caching strategy for re-using recently yielded spectra.
    3. Other formats?

### Disclaimer
This library was made in part to learn Rust, so it may not use the preferred idioms,
patterns, or libraries. Any recommendations are welcome.