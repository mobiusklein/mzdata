/*!
 * Demo how to read an mzML file that is compressed
 */
use std::env;
use std::fs;
use std::io;
use std::process::exit;
use std::time::Instant;

use mzdata::io::PreBufferedStream;
use mzdata::io::{MzMLReader, RestartableGzDecoder};
use mzdata::prelude::*;

fn main() -> io::Result<()> {
    let input = env::args().nth(1).unwrap_or_else(|| {
        eprintln!("Please provide a file path or '-' for STDIN");
        exit(1)
    });
    let start = Instant::now();
    let groups = if input == "-" {
        let stream =
            RestartableGzDecoder::new(io::BufReader::new(PreBufferedStream::new(io::stdin())?));
        let reader = MzMLReader::new(stream);
        let groups: Vec<_> = reader.into_groups().collect();
        groups
    } else {
        let stream = RestartableGzDecoder::new(io::BufReader::new(fs::File::open(input)?));
        let reader = MzMLReader::new(stream);
        let groups: Vec<_> = reader.into_groups().collect();
        groups
    };
    let spectra: Vec<_> = groups
        .iter()
        .flat_map(|g| g.precursor.iter().chain(g.products.iter()))
        .collect();
    let end = Instant::now();
    eprintln!(
        "Read {} groups with {} spectra in {:0.3?}",
        groups.len(),
        spectra.len(),
        end - start
    );
    assert!(!spectra.is_empty());

    Ok(())
}
