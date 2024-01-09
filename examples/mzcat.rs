use std::{io, path, env, fs};
use std::time;

use rayon::prelude::*;

use mzdata::io::{mzml, mzmlb};
use mzdata::prelude::*;
use mzdata::spectrum::MultiLayerSpectrum;


fn load_file<P: Into<path::PathBuf> + Clone>(path: P) -> io::Result<mzml::MzMLReader<fs::File>> {
    let reader = mzml::MzMLReader::open_path(path)?;
    Ok(reader)
}

#[cfg(feature = "mzmlb")]
fn load_mzmlb_file<P: Into<path::PathBuf> + Clone>(path: P) -> io::Result<mzmlb::MzMLbReader> {
    let reader = mzmlb::MzMLbReader::open_path(&path.into())?;
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
        Err(_) => 4,
    };
    mzmlb::MzMLbReader::set_blosc_nthreads(blosc_threads);
    Ok(reader)
}

fn scan_file<
    R: MZFileReader + Iterator<Item=MultiLayerSpectrum> + Send,
>(
    reader: &mut R,
) {
    let start = time::Instant::now();
    reader.enumerate().par_bridge().for_each(|(i, scan)| {
        if i % 10000 == 0 {
            println!(
                "\tScan {}: {}|{} ({} seconds)",
                i,
                scan.id(),
                scan.index(),
                (time::Instant::now() - start).as_secs_f64(),
            );
        }
    });
    let end = time::Instant::now();
    println!("Loaded in {} seconds", (end - start).as_secs_f64());
}


fn main() -> io::Result<()> {
    let path = path::PathBuf::from(
        env::args().nth(1)
            .expect("Please pass an MS data file path"),
    );
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
        } else if ext.to_string_lossy().to_lowercase() == "mzml" {
            let mut reader = load_file(path)?;
            scan_file(&mut reader)
        } else {
            panic!("Could not infer the file format")
        }
    } else {
        let mut reader = load_file(path)?;
        scan_file(&mut reader)
    };
    Ok(())
}