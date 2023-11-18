# mzdata
[![Latest Version](https://img.shields.io/crates/v/mzdata.svg)](https://crates.io/crates/mzdata)

A Rust library for reading mass spectrometry data file formats.

## Quickstart
```rust
use std::fs;
use mzdata::io::prelude::*;
use mzpeaks::{Tolerance, prelude::*};
use mzdata::io::MzMLReader;
use mzdata::spectrum::{SignalContinuity};

let reader = MzMLReader::new(File::open("./test/data/small.mzML").unwrap());
for spectrum in reader {
    println!("Scan {} => BP {}", spectrum.id(), spectrum.peaks().base_peak().mz);
    if spectrum.signal_continuity() < SignalContinuity::Profile {
        let peak_picked = spectrum.into_centroid().unwrap();
        println!("Matches for 579.155: {:?}", peak_picked.peaks.all_peaks_for(579.155, Tolerance::Da(0.02)));
    }
}
println!("MS1 Count: {}\nMSn Count: {}", ms1_count, msn_count);
assert_eq!(ms1_count, 14);
assert_eq!(msn_count, 34);
```

### Supported Formats
1. `mzML` and `indexedmzML`
2. `MGF`
3. `mzMLb`

### Disclaimer
This library was made in part to learn Rust, so it may not use the preferred idioms,
patterns, or libraries. Any recommendations are welcome.