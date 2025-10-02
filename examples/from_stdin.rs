/*!
 * Demo how to read a data file from STDIN, and then batch collect its spectra
 */
use std::io::{self, Seek};
use std::time::Instant;

use mzdata::io::{
    infer_from_stream, MassSpectrometryFormat, PreBufferedStream, RestartableGzDecoder,
    SpectrumSource,
};
use mzdata::{MGFReader, MzMLReader};

fn main() -> io::Result<()> {
    env_logger::init();
    let start = Instant::now();
    let stream = io::stdin();
    let mut stream = PreBufferedStream::new(stream)?;
    let (fmt, compressed) = infer_from_stream(&mut stream)?;
    stream.seek(io::SeekFrom::Start(0))?;
    let groups: Vec<_> = match fmt {
        MassSpectrometryFormat::MGF => {
            if compressed {
                MGFReader::new(RestartableGzDecoder::new(io::BufReader::new(stream)))
                    .into_groups()
                    .collect()
            } else {
                MGFReader::new(stream).into_groups().collect()
            }
        }
        MassSpectrometryFormat::MzML => {
            if compressed {
                MzMLReader::new(RestartableGzDecoder::new(io::BufReader::new(stream)))
                    .into_groups()
                    .collect()
            } else {
                MzMLReader::new(stream).into_groups().collect()
            }
        }
        x => {
            panic!("Cannot identify file format ({:?})", x)
        }
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
