#![allow(unused)]
use std::env;
use std::fs;
use std::path;

use mzdata::io::prelude::*;
use mzdata::io::{mgf, mzml, ScanSource};
use mzdata::peaks::PeakCollection;
use mzdata::spectrum::{CentroidSpectrum, SignalContinuity, SpectrumBehavior};
use mzdata::MassErrorType;

fn main() {
    let args: Vec<String> = env::args().collect();
    let path: &path::Path;
    if args.len() > 1 {
        path = path::Path::new(&args[1]);
    } else {
        path = path::Path::new("./test/data/small.mzML");
    }
    println!("Path: {}", path.to_str().unwrap());
    let file = fs::File::open(path).unwrap();
    let mut parser = mzml::MzMLReader::new_indexed(file);
    let iter = parser.iter();
    for scan in iter {
        println!("Scan {} => BP {}", scan.id(), scan.peaks().base_peak().1);
        if scan.signal_continuity() < SignalContinuity::Profile {
            let peak_picked = scan.into_centroid().unwrap();
            println!(
                "Matches for 579.155: {:?}",
                peak_picked
                    .peaks
                    .all_peaks_for(579.155, 0.02, MassErrorType::Exact)
            );
        }
    }
}
