use std::{io, fs, path, env};

use mzdata_param::MSVocabulary;
use mzcv::{CVIndex, CVSource, CVStructure, HashBufReader};


fn main() -> io::Result<()> {
    pretty_env_logger::init_timed();

    let mut args = env::args().skip(1);

    let infile: path::PathBuf = args.next().unwrap_or_else(|| panic!("Please provide a path to an OBO file for PSI-MS")).into();
    let outfile: path::PathBuf = args.next().unwrap_or_else(|| panic!("Please provide a path to write the static data to")).into();

    if !infile.exists() {
        panic!("Input file {} not found", infile.display())
    }

    let cv_reader: HashBufReader<Box<dyn std::io::Read>, sha1::Sha1> = if infile.extension().is_some_and(|s| s == "gz") {
        HashBufReader::boxed(flate2::read::GzDecoder::new(fs::File::open(&infile)?))
    } else {
        HashBufReader::boxed(fs::File::open(&infile)?)
    };

    let (ver, data, errs) = MSVocabulary::parse([cv_reader].into_iter()).unwrap();
    log::info!("Read {:?}/{} with {} items", ver.version, ver.hash_hex(), data.len());
    if !errs.is_empty() {
        for e in errs {
            log::error!("Failed to parse: {e}")
        }
        panic!("One or more errors occurred while parsing CV")
    }

    let mut cv = CVIndex::<MSVocabulary>::empty();

    cv.update_from_structure(ver, data).unwrap();

    log::info!("Writing {} items to {}", cv.len(), outfile.display());
    cv.save_to_cache_at(&outfile).unwrap();
    let buf = fs::read(&outfile)?;
    log::info!("{} bytes uncompressed", buf.len());
    let level = flate2::Compression::best();
    let mut compress = flate2::Compress::new(level, true);
    let mut outbuf = Vec::new();
    outbuf.reserve(buf.len() * 2);
    compress.compress_vec(&buf, &mut outbuf, flate2::FlushCompress::Sync)?;
    log::info!("{} bytes compressed", outbuf.len());
    fs::write(&outfile, outbuf)?;
    Ok(())
}