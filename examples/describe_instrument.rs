use std::env;
use std::io;

use mzdata::mz_read;
use mzdata::prelude::*;

fn main() -> io::Result<()> {
    env_logger::init();
    let inpath = env::args().nth(1).unwrap_or_else(|| {
        eprintln!("Please provide a path to read an MS data file from, or '-'");
        std::process::exit(1)
    });

    let configs = mz_read!(inpath, reader => {
        reader.instrument_configurations().clone()
    })?;

    for (k, config) in configs.iter() {
        println!("Configuration ID: {k}");
        for component in config.components.iter() {
            println!("\t{:?} -> {}\n", component.component_type, component.order);
            for p in component.params() {
                println!("\t{p}");
            }
        }
        for p in config.params() {
            println!("{p}");
        }
    }

    Ok(())
}
