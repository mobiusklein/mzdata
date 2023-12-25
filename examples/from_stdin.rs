use std::io;
use std::time::Instant;

use mzdata::MzMLReader;
use mzdata::io::{PreBufferedStream, ScanSource};



fn main() -> io::Result<()> {
    let start = Instant::now();
    let stream = io::stdin();
    let stream = PreBufferedStream::new(stream)?;

    let reader = MzMLReader::new(stream);
    let groups: Vec<_> = reader.into_groups().collect();
    let spectra: Vec<_> = groups.iter().flat_map(|g| {
        g.precursor.iter().chain(g.products.iter())
    }).collect();
    let end = Instant::now();
    eprintln!("Read {} groups with {} spectra in {:0.3?}", groups.len(), spectra.len(), end - start);
    assert!(spectra.len() > 0);
    Ok(())
}