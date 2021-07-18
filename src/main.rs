#![allow(unused)]
use std::fs;
use std::path;
// use std::env;

use mzdata::io::{mzml, mgf, ScanSource};
use mzdata::peaks::PeakCollection;
use mzdata::spectrum::{SignalContinuity, SpectrumBehavior};

fn main() {
    // let args: Vec<String> = env::args().collect();
    let path: &path::Path;
    path = path::Path::new(
        "./test/data/small.mgf",
    );
    println!("Path: {}", path.to_str().unwrap());
    let file = fs::File::open(path).unwrap();
    let mut parser = mgf::MGFReader::new_indexed(file);
    for scan in &mut parser {
        println!(
            "Scan ID {}, MS Level {}",
            scan.id(),
            scan.description().ms_level
        );
        scan.peaks.len();
    }
}
