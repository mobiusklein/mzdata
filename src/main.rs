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
    let mut reader = mzml::MzMLReader::open_path(path)?;
    let end = time::Instant::now();
    println!(
        "Loaded in {} seconds",
        (end - start).as_secs()
    );
    println!("{} spectra", reader.len());

    Ok(())
}
