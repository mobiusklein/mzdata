use std::{env, io, path};

use log::info;
use mzdata::io::MZReader;
use mzdata::prelude::*;

fn main() -> io::Result<()> {
    env_logger::init();
    let mut args = env::args().skip(1);

    let path = path::PathBuf::from(args.next().expect("Please pass an MS data file path"));

    let key = args
        .next()
        .expect("Please provide a key type, \"id\", \"index\" or \"time\"");
    let key_value = args
        .next()
        .expect("Please provide a key value matching the key type");

    info!("Opening {}", path.display());
    let mut reader = MZReader::open_path(path)?;

    let spectrum = match key.as_str() {
        "id" => reader.get_spectrum_by_id(&key_value).unwrap(),
        "index" => reader
            .get_spectrum_by_index(key_value.parse().unwrap())
            .unwrap(),
        "time" => reader
            .get_spectrum_by_time(key_value.parse().unwrap())
            .unwrap(),
        _ => {
            panic!("Unknown key type {}", key);
        }
    };

    // dbg!(spectrum);
    println!("ID: {}; Index: {}; Time: {}", spectrum.id(), spectrum.index(), spectrum.start_time());
    println!("Num data points: {}", spectrum.raw_arrays().unwrap().mzs().unwrap().len());

    Ok(())
}
