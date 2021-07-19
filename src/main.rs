#![allow(unused)]
use std::fs;
use std::path;
use std::env;

use rayon::prelude::*;

use mzdata::io::{mgf, mzml, ScanSource};
use mzdata::io::prelude::*;
use mzdata::peaks::PeakCollection;
use mzdata::spectrum::{SignalContinuity, SpectrumBehavior, CentroidSpectrum};

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
        println!("Scan {} => TIC {}", scan.id(), scan.intensities().iter().sum::<f32>())
    }
}
