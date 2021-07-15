use std::fs;
use std::path;
// use std::env;

// use structured::MGFReader;
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
    path = path::Path::new("C:\\Users\\Joshua\\Dev\\ms_deisotope\\ms_deisotope\\test\\test_data\\small.mzML");
    println!("Path: {}", path.to_str().unwrap());
    let file = fs::File::open(path).unwrap();
    let mut parser = mzml::MzMLReader::new(file);
    let desc = parser.read_next();
    println!("Result: {:?}", desc);
    let desc = parser.read_next();
    println!("Result: {:?}", desc);
}