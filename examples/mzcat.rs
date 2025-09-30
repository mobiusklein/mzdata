use std::time;
use std::{env, io, path};

use mzdata::spectrum::MultiLayerSpectrum;
use mzdata::{MZReader, prelude::*};
use rayon::prelude::*;

fn scan_file<R: MZFileReader + Iterator<Item = MultiLayerSpectrum> + Send>(reader: &mut R) {
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
    env_logger::init();
    let path = path::PathBuf::from(
        env::args()
            .nth(1)
            .expect("Please pass an MS data file path"),
    );

    let mut reader = MZReader::open_path(&path)?;
    scan_file(&mut reader);
    Ok(())
}
