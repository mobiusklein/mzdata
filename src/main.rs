#![allow(unused)]
use std::collections::HashMap;
use std::env;
use std::fs;
use std::io;
use std::path;
use std::time;

use mzdata::io::prelude::*;
use mzdata::io::{mgf, mzml, offset_index, ScanSource};
#[cfg(feature = "mzmlb")]
use mzdata::io::mzmlb;
use mzdata::spectrum::{DeconvolutedSpectrum, MultiLayerSpectrum, SpectrumBehavior};
use mzpeaks::PeakCollection;
use mzpeaks::{CentroidPeak, DeconvolutedPeak};

fn load_file<P: Into<path::PathBuf> + Clone>(path: P) -> io::Result<mzml::MzMLReader<fs::File>> {
    let start = time::Instant::now();
    let reader = mzml::MzMLReader::open_path(path)?;
    let end = time::Instant::now();
    println!("Index Loaded/Built in {} seconds", (end - start).as_secs());
    println!("{} spectra", reader.len());
    Ok(reader)
}

#[cfg(feature = "mzmlb")]
fn load_mzmlb_file<P: Into<path::PathBuf> + Clone>(path: P) -> io::Result<mzmlb::MzMLbReader> {
    let start = time::Instant::now();
    let reader = mzmlb::MzMLbReader::open_path(&path.into())?;
    let end = time::Instant::now();
    println!("Index Loaded/Built in {} seconds", (end - start).as_secs());
    let blosc_threads = match std::env::var("BLOSC_NUM_THREADS") {
        Ok(val) => {
            match val.parse() {
                Ok(nt) => nt,
                Err(e) => {
                    eprintln!("Failed to parse BLOSC_NUM_THREADS env var: {}", e);
                    4
                },
            }
        },
        Err(e) => 4,
    };
    mzmlb::MzMLbReader::set_blosc_nthreads(blosc_threads);
    println!("{} spectra", reader.len());
    Ok(reader)
}

fn scan_file<
    R: MZFileReader<
        CentroidPeak,
        DeconvolutedPeak,
        MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak>,
    >,
>(
    reader: &mut R,
) -> (
    HashMap<u8, usize>,
    HashMap<i32, usize>,
    usize,
    HashMap<u8, HashMap<i32, usize>>,
) {
    let start = time::Instant::now();
    let mut level_table: HashMap<u8, usize> = HashMap::new();
    let mut charge_table: HashMap<i32, usize> = HashMap::new();
    let mut peak_charge_table: HashMap<u8, HashMap<i32, usize>> = HashMap::new();
    let mut peak_count: usize = 0;
    for (i, scan) in reader.enumerate() {
        if i % 10000 == 0 {
            println!(
                "\tScan {}: {} ({} seconds, {} peaks)",
                i,
                scan.id(),
                (time::Instant::now() - start).as_secs_f64(),
                peak_count
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
        let has_charge = if let Some(ref arrays) = scan.arrays {
            peak_count += arrays.mzs().unwrap().len();
            arrays.charges().is_ok()
        } else {
            false
        };
        if has_charge {
            let deconv_scan: DeconvolutedSpectrum = scan.try_into().unwrap();
            deconv_scan.deconvoluted_peaks.iter().for_each(|p| {
                *(*peak_charge_table.entry(deconv_scan.ms_level()).or_default())
                    .entry(p.charge)
                    .or_default() += 1;
                assert!((p.index as usize) < deconv_scan.deconvoluted_peaks.len())
            })
        }
    }
    let end = time::Instant::now();
    println!("Loaded in {} seconds", (end - start).as_secs_f64());
    (level_table, charge_table, peak_count, peak_charge_table)
}

fn main() -> io::Result<()> {
    let path = path::PathBuf::from(
        env::args()
            .skip(1)
            .next()
            .unwrap_or("./test/data/small.mzML".to_owned()),
    );

    let (level_table, charge_table, peak_count, peak_charge_table_by_level) =
        if let Some(ext) = path.extension() {
            if ext.to_string_lossy().to_lowercase() == "mzmlb" {
                #[cfg(feature = "mzmlb")]
                {
                    let mut reader = load_mzmlb_file(path)?;
                    scan_file(&mut reader)
                }
                #[cfg(not(feature = "mzmlb"))]
                {
                    panic!("Cannot read mzMLb file. Recompile enabling the `mzmlb` feature")
                }
            } else {
                let mut reader = load_file(path)?;
                scan_file(&mut reader)
            }
        } else {
            let mut reader = load_file(path)?;
            scan_file(&mut reader)
        };

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
        if **charge == 0 {
            println!("\tCharge Not Reported: {}", count);
        } else {
            println!("\t{}: {}", charge, count);
        }
    }

    let mut peak_charge_levels: Vec<_> = peak_charge_table_by_level.iter().collect();

    peak_charge_levels.sort_by(|(level_a, _), (level_b, _)| level_a.cmp(level_b));

    for (level, peak_charge_table) in peak_charge_levels {
        if peak_charge_table.len() > 0 {
            println!("Peak Charge States for MS level {}:", level);
            let mut peak_charge_set: Vec<(&i32, &usize)> = peak_charge_table.iter().collect();
            peak_charge_set.sort_by_key(|(a, b)| *a);
            for (charge, count) in peak_charge_set.iter() {
                if **charge == 0 {
                    println!("\tCharge Not Reported: {}", count);
                } else {
                    println!("\t{}: {}", charge, count);
                }
            }
        }
    }

    println!("{} peaks read", peak_count);

    Ok(())
}
