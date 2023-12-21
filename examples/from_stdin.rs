use std::io;
use std::time::Instant;

use mzdata::MzMLReader;
use mzdata::io::PreBufferedStream;



fn main() -> io::Result<()> {
    let start = Instant::now();
    let stream = io::stdin();
    let stream = PreBufferedStream::new(stream)?;

    let reader = MzMLReader::new(stream);
    let spectra: Vec<_> = reader.collect();
    let end = Instant::now();
    eprintln!("Read {} spectra in {:0.3?}", spectra.len(), end - start);
    assert!(spectra.len() > 0);
    Ok(())
}