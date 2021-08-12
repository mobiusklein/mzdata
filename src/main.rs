#![allow(unused)]
use std::collections::HashMap;
use std::env;
use std::fs;
use std::io;
use std::path;
use std::time;

use mzdata::io::prelude::*;
use mzdata::io::{mgf, mzml, offset_index, ScanSource};
use mzdata::spectrum::SpectrumBehavior;
use mzpeaks::PeakCollection;

fn load_file<P: Into<path::PathBuf> + Clone>(path: P) -> io::Result<mzml::MzMLReader<fs::File>> {
    let start = time::Instant::now();
    let reader = mzml::MzMLReader::open_path(path)?;
    let end = time::Instant::now();
    println!("Loaded in {} seconds", (end - start).as_secs());
    println!("{} spectra", reader.len());
    Ok(reader)
}

fn scan_file<R: SeekRead>(
    reader: &mut mzml::MzMLReader<R>,
) -> (HashMap<u8, usize>, HashMap<i32, usize>, usize) {
    let start = time::Instant::now();
    let mut level_table: HashMap<u8, usize> = HashMap::new();
    let mut charge_table: HashMap<i32, usize> = HashMap::new();
    let mut peak_count: usize = 0;
    for (i, scan) in reader.enumerate() {
        if i % 10000 == 0 {
            println!(
                "\tScan {}: {} ({} seconds)",
                i,
                scan.id(),
                (time::Instant::now() - start).as_secs()
            );
        }
        let level = scan.ms_level();
        *level_table.entry(level).or_default() += 1;
        if level > 1 {
            if let Some(charge) = scan.precursor().unwrap().ion.charge {
                *charge_table.entry(charge).or_default() += 1;
            } else {
                *charge_table.entry(0).or_default() += 1;
            }
        }
        peak_count += scan.arrays.unwrap().mzs().len();
    }
    let end = time::Instant::now();
    println!("Loaded in {} seconds", (end - start).as_secs());
    (level_table, charge_table, peak_count)
}

fn main() -> io::Result<()> {
    let path = path::PathBuf::from(
        env::args()
            .skip(1)
            .next()
            .unwrap_or("./test/data/small.mzML".to_owned()),
    );

    let mut reader = load_file(path)?;
    let (level_table, charge_table, peak_count) = scan_file(&mut reader);

    println!("MS Levels:");
    let mut level_set: Vec<(&u8, &usize)> = level_table.iter().collect();
    level_set.sort_by_key(|(a, b)| *a);
    for (level, count) in level_set.iter() {
        println!("\t{}: {}", level, count);
    }

    println!("Charge States:");
    let mut charge_set: Vec<(&i32, &usize)> = charge_table.iter().collect();
    charge_set.sort_by_key(|(a, b)| *a);
    for (charge, count) in charge_set.iter() {
        println!("\t{}: {}", charge, count);
    }

    println!("{} peaks read", peak_count);

    Ok(())
}
