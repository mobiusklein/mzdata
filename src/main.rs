use std::fs;
use std::path;
// use std::env;

use structured::spectrum::{SpectrumBehavior, ScanSiganlContinuity};
use structured::peak_set::{PeakCollection};
use structured::io::mzml;

fn main() {
    // let args: Vec<String> = env::args().collect();
    let path: &path::Path;
    path = path::Path::new("C:\\Users\\Joshua\\Dev\\ms_deisotope\\ms_deisotope\\test\\test_data\\small.mzML");
    println!("Path: {}", path.to_str().unwrap());
    let file = fs::File::open(path).unwrap();
    let mut parser = mzml::MzMLReader::new_indexed(file);
    for scan in &mut parser {
        println!("Scan ID {}, MS Level {}", scan.get_id(), scan.get_description().ms_level);
        if scan.get_description().is_profile != ScanSiganlContinuity::Centroid {
            println!("Profile spectrum");
        } else {
            let cscan = scan.into_centroid().expect("Coercion to centroid failed");
            println!("{} Peaks", cscan.peaks.len());
        }
    }
    let scan = parser.get_spectrum_by_index(0).expect("Failed to get scan by index");
    println!("Scan ID {}, MS Level {}", scan.get_id(), scan.get_description().ms_level);

    let scan = parser.get_spectrum_by_id(
        "controllerType=0 controllerNumber=1 scan=44").expect("Failed to get scan by id");
    println!("Scan ID {}, MS Level {}", scan.get_id(), scan.get_description().ms_level);
}