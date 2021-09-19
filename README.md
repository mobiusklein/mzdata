# mzdata
[![Latest Version](https://img.shields.io/crates/v/mzdata.svg)](https://crates.io/crates/mzdata)

A Rust library for reading mass spectrometry data file formats.

## Quickstart
```rust
use std::fs::File;
use mzdata::io::prelude::*;
use mzpeaks::{MassErrorType, prelude::*};
use mzdata::io::MzMLReader;
use mzdata::spectrum::{SignalContinuity};

let reader = MzMLReader::new(File::open("./test/data/small.mzML").unwrap());
for spectrum in reader {
    println!("Scan {} => BP {}", spectrum.id(), spectrum.peaks().base_peak().mz);
    if spectrum.signal_continuity() < SignalContinuity::Profile {
        let peak_picked = spectrum.into_centroid().unwrap();
        println!("Matches for 579.155: {:?}", peak_picked.peaks.all_peaks_for(579.155, 0.02, MassErrorType::Absolute));
    }
}
```

### Supported Formats
1. `mzML` and `indexedmzML`
2. MGF

## TODO
 1. mzML Implementation
    1. Improve unit handling.
 2. Iteration behavior
    1. Implement a caching strategy for re-using recently yielded spectra.
 3. Other formats?

### Disclaimer
This library was made in part to learn Rust, so it may not use the preferred idioms,
patterns, or libraries. Any recommendations are welcome.