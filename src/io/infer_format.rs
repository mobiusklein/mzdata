use std::fs;
use std::io;
use std::io::BufReader;
use std::path;

use std::io::prelude::*;

use flate2::bufread::GzDecoder;

use crate::MGFReader;
use crate::MzMLReader;

#[cfg(feature = "mzmlb")]
pub use crate::MzMLbReader;

use crate::io::traits::ScanSource;
use crate::io::mzml::is_mzml;
use crate::io::mgf::is_mgf;
use crate::io::compression::{is_gzipped, is_gzipped_extension};

/// Mass spectrometry file formats that [`mzdata`](crate)
/// supports
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MassSpectrometryFormat {
    MGF,
    MzML,
    #[cfg(feature = "mzmlb")]
    MzMLb,
    Unknown
}


/// Given a path, infer the file format and whether or not the file at that path is
/// GZIP compressed
pub fn infer_from_path<P: Into<path::PathBuf>,>(path: P) -> (MassSpectrometryFormat, bool) {
    let path: path::PathBuf = path.into();
    let (is_gzipped, path) = is_gzipped_extension(path);
    if let Some(ext) = path.extension() {
        if let Some(ext) = ext.to_ascii_lowercase().to_str() {
            let form = match ext {
                "mzml" => MassSpectrometryFormat::MzML,
                "mgf" => MassSpectrometryFormat::MGF,
                #[cfg(feature = "mzmlb")]
                "mzmlb" => MassSpectrometryFormat::MzMLb,
                _ => MassSpectrometryFormat::Unknown
            };
            (form, is_gzipped)
        } else {
            (MassSpectrometryFormat::Unknown, is_gzipped)
        }
    } else {
        (MassSpectrometryFormat::Unknown, is_gzipped)
    }
}


/// Given a stream of bytes, infer the file format and whether or not the
/// stream is GZIP compressed. This assumes the stream is seekable.
pub fn infer_from_stream<R: Read + Seek>(stream: &mut R) -> io::Result<(MassSpectrometryFormat, bool)> {
    let mut buf = Vec::with_capacity(100);
    buf.resize(100, b'\0');
    let current_pos = stream.stream_position()?;
    stream.read_exact(buf.as_mut_slice())?;
    let is_stream_gzipped = is_gzipped(buf.as_slice());
    if is_stream_gzipped {
        let mut decoder = GzDecoder::new(buf.as_slice());
        let mut decompressed_buf = Vec::new();
        decoder.read_to_end(&mut decompressed_buf)?;
        buf = decompressed_buf;
    }
    stream.seek(io::SeekFrom::Start(current_pos))?;
    if is_mzml(&buf) {
        Ok((MassSpectrometryFormat::MzML, is_stream_gzipped))
    }
    else if is_mgf(&buf) {
        Ok((MassSpectrometryFormat::MGF, is_stream_gzipped))
    } else {
        Ok((MassSpectrometryFormat::Unknown, is_stream_gzipped))
    }
}


/// Given a path, infer the file format and whether or not the file at that path is
/// GZIP compressed, using both the file name and by trying to open and read the file
/// header
pub fn infer_format<P: Into<path::PathBuf>>(path: P) -> io::Result<(MassSpectrometryFormat, bool)> {
    let path: path::PathBuf = path.into();

    let (format, is_gzipped) = infer_from_path(&path);
    match format {
        MassSpectrometryFormat::Unknown => {
            let handle = fs::File::open(path.clone())?;
            let mut stream = BufReader::new(handle);
            let (format, is_gzipped) = infer_from_stream(&mut stream)?;
            Ok((format, is_gzipped))
        },
        _ => {
            Ok((format, is_gzipped))
        }
    }
}


/// Given a local file system path, infer the file format, and attempt to open it
/// for reading.
pub fn open_file<P: Into<path::PathBuf>>(path: P) -> io::Result<Box<dyn ScanSource>>{
    let path = path.into();
    let (format, is_gzipped) = infer_format(path.clone())?;

    if is_gzipped {
        Err(io::Error::new(io::ErrorKind::Unsupported, "Gzipped files are not supported"))
    } else {
        match format {
            MassSpectrometryFormat::MGF => {
                let handle = fs::File::open(path)?;
                let reader = MGFReader::new_indexed(handle);
                Ok(Box::new(reader))
            },
            MassSpectrometryFormat::MzML => {
                let handle = fs::File::open(path)?;
                let reader = MzMLReader::new_indexed(handle);
                Ok(Box::new(reader))
            },
            #[cfg(feature = "mzmlb")]
            MassSpectrometryFormat::MzMLb => {
                let reader = MzMLbReader::new(&path)?;
                Ok(Box::new(reader))
            }
            _ => {
                Err(io::Error::new(io::ErrorKind::Unsupported, "File format not supported"))
            }
        }
    }
}


#[cfg(test)]
mod test {
    use crate::{SpectrumLike, spectrum::{Spectrum, ArrayType}};

    use super::*;

    #[test]
    fn infer_mzml() {
        let path = path::Path::new("./test/data/small.mzML");
        assert!(path.exists());
        let (fmt, zipped) = infer_from_path(path);
        assert_eq!(fmt, MassSpectrometryFormat::MzML);
        assert!(!zipped);
    }

    #[test]
    fn infer_mgf() {
        let path = path::Path::new("./test/data/small.mgf");
        assert!(path.exists());
        let (fmt, zipped) = infer_from_path(path);
        assert_eq!(fmt, MassSpectrometryFormat::MGF);
        assert!(!zipped);
    }

    #[test]
    fn infer_open() {
        let path = path::Path::new("./test/data/small.mzML");
        assert!(path.exists());
        if let Ok(mut reader) = open_file(path) {
            assert_eq!(reader.len(), 48);

            if let Some(spec) = reader.get_spectrum_by_index(10) {
                let spec: Spectrum = spec;
                assert!(spec.index() == 10);
                assert!(spec.id() == "controllerType=0 controllerNumber=1 scan=11");
                if let Some(data_arrays) = &spec.arrays {
                    assert!(data_arrays.has_array(&ArrayType::MZArray));
                    assert!(data_arrays.has_array(&ArrayType::IntensityArray));
                    let mzs = data_arrays.mzs().unwrap();
                    assert!(mzs.len() == 941);
                }
            }
        } else {
            panic!("Failed to open file")
        }
    }
}