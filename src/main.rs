use std::fs;
use std::path;
// use std::env;

use std::io::prelude::{Seek};
use std::io::SeekFrom;

use structured::spectrum::{SpectrumBehavior, ScanSiganlContinuity};
use structured::peak_set::{PeakCollection};
use structured::io::mzml;

fn main() {
    // let args: Vec<String> = env::args().collect();
    let path: &path::Path;
    // if args.len() > 1 {
    //     path = path::Path::new(&args[1]);
    // } else {
    //     path = path::Path::new("C:\\Users\\Joshua\\Dev\\ms_deisotope\\ms_deisotope\\test\\test_data\\small.mgf");
    // }
    // println!("Path: {}", path.to_str().unwrap());
    // let file = fs::File::open(path).unwrap();
    // let reader = MGFReader::new(file);
    // let mut counter = 0;
    // for _scan in reader {
    //     // println!("Scan ID: {}", _scan.description.id);
    //     counter += 1;
    // }
    // println!("{} scans in {}", counter, path.to_str().unwrap());
    path = path::Path::new("C:\\Users\\Joshua\\Dev\\ms_deisotope\\ms_deisotope\\test\\test_data\\three_test_scans.mzML");
    println!("Path: {}", path.to_str().unwrap());
    let file = fs::File::open(path).unwrap();
    let mut parser = mzml::MzMLReader::new(file);
    println!("{}", parser.seek(SeekFrom::Start(158797)).expect("what?"));
    let scan = parser.read_next();
    println!("{}", parser.read_next().expect("Read next spectrum failed").get_id());
    for scan in parser {
        println!("Scan ID {}, MS Level {}", scan.get_id(), scan.get_description().ms_level);
        if scan.get_description().is_profile != ScanSiganlContinuity::Centroid {
            println!("Profile spectrum");
        } else {
            let cscan = scan.into_centroid().expect("Coercion to centroid failed");
            println!("{} Peaks", cscan.peaks.len());
        }
    }
}