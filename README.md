# mzdata
[![Latest Version](https://img.shields.io/crates/v/mzdata?style=for-the-badge&color=mediumpurple&logo=rust)](https://crates.io/crates/mzdata)
[![docs.rs](https://img.shields.io/docsrs/mzdata?style=for-the-badge&logo=docs.rs&color=mediumseagreen)](https://docs.rs/mzdata/latest/mzdata/)

A Rust library for reading mass spectrometry data file formats.

## Quickstart
```rust
use std::fs;
use mzdata::prelude::*;
use mzpeaks::Tolerance;
use mzdata::MzMLReader;
use mzdata::spectrum::SignalContinuity;

fn main() {
    let mut ms1_count = 0;
    let mut msn_count = 0;
    let reader = MzMLReader::open_path("./test/data/small.mzML").unwrap();
    for spectrum in reader {
        if spectrum.ms_level() == 1 {
            ms1_count += 1;
        } else {
            msn_count += 1;
        }
        println!("Scan {} => BP {}", spectrum.id(), spectrum.peaks().base_peak().mz);
        if spectrum.signal_continuity() == SignalContinuity::Centroid {
            let peak_picked = spectrum.into_centroid().unwrap();
            println!("Matches for 579.155: {:?}", peak_picked.peaks.all_peaks_for(579.155, Tolerance::Da(0.02)));
        }
    }
    println!("MS1 Count: {}\nMSn Count: {}", ms1_count, msn_count);
    assert_eq!(ms1_count, 14);
    assert_eq!(msn_count, 34);
}


```

## Supported Formats
1. `mzML` and `indexedmzML`
2. `MGF`
3. `mzMLb`
4. Thermo RAW
5. Bruker TDF
6. `imzML`
7. PROXI

## Disclaimer
This library was made in part to learn Rust, so it may not use the preferred idioms,
patterns, or libraries. Any recommendations are welcome.