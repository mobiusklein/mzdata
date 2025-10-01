/*!
 * Demo the minimum code needed to infer the input file format from a path
 * or STDIN using `infer_format`, `infer_from_stream` and
 */
use std::env;
use std::io;
use std::process::exit;

use mzdata::io::{PreBufferedStream, infer_format, infer_from_stream};

fn main() -> io::Result<()> {
    let input = env::args().nth(1).unwrap_or_else(|| {
        eprintln!("Please provide a file path or '-' for STDIN");
        exit(1)
    });

    let (inferred, gzipped) = if input == "-" {
        let mut stream = PreBufferedStream::new(io::stdin())?;
        infer_from_stream(&mut stream)?
    } else {
        infer_format(input)?
    };
    println!("{:?} (gzip: {})", inferred, gzipped);
    Ok(())
}
