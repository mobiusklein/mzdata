#![allow(unused)]
use std::fs;
use std::path;
// use std::env;

use rayon::prelude::*;

use mzdata::io::{mgf, mzml, ScanSource};
use mzdata::io::prelude::*;
use mzdata::peaks::PeakCollection;
use mzdata::spectrum::{SignalContinuity, SpectrumBehavior, CentroidSpectrum};

fn main() {
    // let args: Vec<String> = env::args().collect();
    let path: &path::Path;
    path = path::Path::new("./test/data/small.mzML");
    println!("Path: {}", path.to_str().unwrap());
    let file = fs::File::open(path).unwrap();
    let mut parser = mzml::MzMLReader::new_indexed(file);
    let iter = parser.iter();
    let tot: f32 = iter.map(|s|s.intensities().iter().sum::<f32>()).sum();
    let iter = parser.iter();
    let par_tot: f32 = iter.par_bridge().map(|s|s.intensities().iter().sum::<f32>()).sum();
    println!("Linear: {}\nParallel: {}", tot, par_tot);
}
