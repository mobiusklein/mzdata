use std::{env, io, path};

use mzdata::io::mzml;
use mzdata::prelude::*;

fn main() -> io::Result<()> {
    let path = path::PathBuf::from(
        env::args()
            .nth(1)
            .expect("Please pass an MS data file path"),
        // "test/data/batching_test.mzML"
    );

    let mut reader = mzml::MzMLReader::open_path(path)?;

    let n_spectra = reader.len();

    // Find the spectrum at the midpoint of the run
    let spec = reader.get_spectrum_by_index(n_spectra / 2).unwrap();
    eprintln!(
        "Midpoint spectrum {} (level {}) at time {}",
        spec.id(),
        spec.ms_level(),
        spec.start_time()
    );

    // Jump the iterator to that point in time
    reader.start_from_time(spec.start_time())?;
    let s = reader.next().unwrap();
    eprintln!(
        "Resuming at {} (level {}) at time {}",
        s.id(),
        s.ms_level(),
        s.start_time()
    );

    // Convert the iterator into a group iterator
    let mut group_iter = reader.into_groups();
    // Jump the group iterator to that point in time (If an MSn spectrum was found, the next MS1 may be shown instead)
    group_iter.start_from_time(spec.start_time())?;
    let g = group_iter.next().unwrap();
    eprintln!(
        "Resuming at group having {:?} at time {:?}",
        g.earliest_spectrum().and_then(|s| Some(s.id())),
        g.earliest_spectrum().and_then(|s| Some(s.start_time()))
    );

    // Convert the group iterator into an averaging group iterator
    let mut avg_iter = group_iter.averaging(1, 200.0, 2200.0, 0.005);
    // Jump the group iterator to that point in time (If an MSn spectrum was found, the next MS1 may be shown instead)
    avg_iter.start_from_time(spec.start_time())?;
    let g = avg_iter.next().unwrap();
    eprintln!(
        "Resuming at group having {:?} at time {:?}",
        g.earliest_spectrum().and_then(|s| Some(s.id())),
        g.earliest_spectrum().and_then(|s| Some(s.start_time()))
    );

    Ok(())
}
