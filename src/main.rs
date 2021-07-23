#![allow(unused)]
use std::env;
use std::io;
use std::path;
use std::time;
use std::fs;

use mzdata::io::prelude::*;
use mzdata::io::{mgf, mzml, offset_index, ScanSource};
use mzdata::peaks::PeakCollection;
use mzdata::spectrum::{CentroidSpectrum, SignalContinuity, SpectrumBehavior};

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let path: &path::Path;
    if args.len() > 1 {
        path = path::Path::new(&args[1]);
    } else {
        path = path::Path::new("./test/data/small.mzML");
    }
    let start = time::Instant::now();
    let fh = fs::File::open(path)?;
    // let mut reader = mzml::MzMLReader::open_path(path)?;
    let mut reader = mzml::MzMLReader::new(fh);
    let index_path = path.with_extension("index.json");
    println!("Index Path: {}", index_path.display());
    let index_fh = fs::File::open(index_path)?;
    reader.read_index(index_fh);
    let end = time::Instant::now();
    println!(
        "Loaded in {} seconds",
        (end - start).as_secs()
    );
    println!("{} spectra", reader.len());

    Ok(())
}
