use std::io;
use std::time::Instant;

use mzdata::{MzMLReader, MGFReader};
use mzdata::io::{PreBufferedStream, ScanSource, infer_from_stream, MassSpectrometryFormat};

fn main() -> io::Result<()> {
    let start = Instant::now();
    let stream = io::stdin();
    let mut stream = PreBufferedStream::new(stream)?;
    let (fmt, compressed) = infer_from_stream(&mut stream)?;
    if compressed {
        panic!("Compression not supported!")
    }
    let groups: Vec<_> = match fmt {
        MassSpectrometryFormat::MGF => {
            MGFReader::new(stream).into_groups().collect()
        },
        MassSpectrometryFormat::MzML => {
            MzMLReader::new(stream).into_groups().collect()
        },
        _ => {
            panic!("Cannot identify file format")
        }
    };
    let spectra: Vec<_> = groups.iter().flat_map(|g| {
        g.precursor.iter().chain(g.products.iter())
    }).collect();
    let end = Instant::now();
    eprintln!("Read {} groups with {} spectra in {:0.3?}", groups.len(), spectra.len(), end - start);
    assert!(spectra.len() > 0);
    Ok(())
}