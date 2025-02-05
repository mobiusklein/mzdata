use std::time;
use std::{env, io, path};

use futures::StreamExt;
use tokio::fs;

use mzdata::io::mzml;
use mzdata::prelude::*;

async fn load_file<P: Into<path::PathBuf> + Clone>(
    path: P,
) -> io::Result<mzml::AsyncMzMLReader<fs::File>> {
    let fh = fs::File::open(path.into()).await?;
    let mut reader = mzml::AsyncMzMLReader::new(fh).await;
    reader
        .read_index_from_end()
        .await
        .expect("Failed to read index from the file");
    Ok(reader)
}

async fn scan_file(reader: &mut mzml::AsyncMzMLReader<fs::File>) {
    let start = time::Instant::now();
    let mut i = 0;

    let mut stream = reader.as_stream();

    while let Some(scan) = stream.next().await {
        if i % 10000 == 0 {
            println!(
                "\tScan {}: {}|{} ({} seconds)",
                i,
                scan.id(),
                scan.index(),
                (time::Instant::now() - start).as_secs_f64(),
            );
        }
        i += 1;
    }
    let end = time::Instant::now();
    println!(
        "Loaded in {} spectra {} seconds",
        i,
        (end - start).as_secs_f64()
    );
}

#[tokio::main(flavor = "multi_thread", worker_threads = 10)]
async fn main() -> io::Result<()> {
    let path = path::PathBuf::from(
        env::args().nth(1)
            .expect("Please pass an MS data file path"),
    );
    if let Some(ext) = path.extension() {
        if ext.to_string_lossy().to_lowercase() == "mzml" {
            let mut reader = load_file(path).await?;
            scan_file(&mut reader).await;
        } else {
            panic!("Could not infer the file format")
        }
    } else {
        let mut reader = load_file(path).await?;
        scan_file(&mut reader).await;
    };
    Ok(())
}
