# mzdata
[![Latest Version](https://img.shields.io/crates/v/mzdata.svg)](https://crates.io/crates/mzdata)

A Rust library for reading mass spectrometry data file formats.

## Quickstart
```rust
use std::fs::File;
use mzdata::io::prelude::*;
use mzdata::MassErrorType;
use mzdata::io::MzMLReader;

let reader = MzMLReader::new(File::open("./test/data/small.mzML").unwrap());
for spectrum in reader {
    println!("Scan ID: {} => {} data points", spectrum.id(), spectrum.mzs().len());
    for scan in iter {
        println!("Scan {} => BP {}", scan.id(), scan.peaks().base_peak().1);
        if scan.signal_continuity() < SignalContinuity::Profile {
            let peak_picked = scan.into_centroid().unwrap();
            println!("Matches for 579.155: {:?}", peak_picked.peaks.all_peaks_for(579.155, 0.02, MassErrorType::Exact));
        }
    }
}
```

### Supported Formats
1. `mzML` and `indexedmzML`
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