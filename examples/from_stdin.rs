use std::io;

use mzdata::MzMLReader;
use mzdata::io::PreBufferedStream;



fn main() -> io::Result<()> {
    let stream = io::stdin();
    let stream = PreBufferedStream::new(stream)?;

    let reader = MzMLReader::new(stream);
    let spectra: Vec<_> = reader.collect();
    eprintln!("Read {} spectra", spectra.len());
    assert!(spectra.len() > 0);

    Ok(())
}