use std::{
    convert::TryFrom,
    fmt::Display,
    fs,
    io::{self, prelude::*, BufReader},
    path,
};

use flate2::bufread::GzDecoder;

use crate::{
    io::compression::{is_gzipped, is_gzipped_extension},
    meta::FormatConversion,
    params::ControlledVocabulary,
    Param,
};

#[cfg(feature = "mgf")]
use crate::io::mgf::is_mgf;

#[cfg(feature = "mzml")]
use crate::io::mzml::is_mzml;

#[cfg(feature = "thermo")]
use crate::io::thermo::is_thermo_raw_prefix;

#[cfg(feature = "bruker_tdf")]
use crate::io::tdf::is_tdf;

/// Mass spectrometry file formats that [`mzdata`](crate)
/// supports
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum MassSpectrometryFormat {
    MGF,
    MzML,
    MzMLb,
    ThermoRaw,
    BrukerTDF,
    Unknown,
}

impl MassSpectrometryFormat {
    pub fn as_conversion(&self) -> Option<FormatConversion> {
        match self {
            MassSpectrometryFormat::MzML => Some(FormatConversion::ConversionToMzML),
            MassSpectrometryFormat::MzMLb => Some(FormatConversion::ConversionToMzMLb),
            _ => None,
        }
    }

    pub fn as_param(&self) -> Option<Param> {
        let p = match self {
            MassSpectrometryFormat::MGF => {
                ControlledVocabulary::MS.const_param_ident("Mascot MGF format", 1001062)
            }
            MassSpectrometryFormat::MzML => {
                ControlledVocabulary::MS.const_param_ident("mzML format", 1000584)
            }
            MassSpectrometryFormat::MzMLb => {
                ControlledVocabulary::MS.const_param_ident("mzMLb format", 1002838)
            }
            MassSpectrometryFormat::ThermoRaw => {
                ControlledVocabulary::MS.const_param_ident("Thermo RAW format", 1000563)
            }
            MassSpectrometryFormat::BrukerTDF => {
                ControlledVocabulary::MS.const_param_ident("Bruker TDF format", 1002817)
            }
            MassSpectrometryFormat::Unknown => return None,
        };
        Some(p.into())
    }
}

impl TryFrom<MassSpectrometryFormat> for Param {
    type Error = &'static str;

    fn try_from(value: MassSpectrometryFormat) -> Result<Self, Self::Error> {
        if let Some(p) = value.as_param() {
            Ok(p)
        } else {
            Err("No conversion")
        }
    }
}

impl TryFrom<MassSpectrometryFormat> for FormatConversion {
    type Error = &'static str;

    fn try_from(value: MassSpectrometryFormat) -> Result<Self, Self::Error> {
        if let Some(p) = value.as_conversion() {
            Ok(p)
        } else {
            Err("No conversion")
        }
    }
}

impl Display for MassSpectrometryFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// Given a path, infer the file format and whether or not the file at that path is
/// GZIP compressed
pub fn infer_from_path<P: Into<path::PathBuf>>(path: P) -> (MassSpectrometryFormat, bool) {
    let path: path::PathBuf = path.into();
    if path.is_dir() {
        #[cfg(feature = "bruker_tdf")]
        if is_tdf(path) {
            return (MassSpectrometryFormat::BrukerTDF, false);
        } else {
            return (MassSpectrometryFormat::Unknown, false);
        }
    }
    let (is_gzipped, path) = is_gzipped_extension(path);
    if let Some(ext) = path.extension() {
        if let Some(ext) = ext.to_ascii_lowercase().to_str() {
            let form = match ext {
                "mzml" => MassSpectrometryFormat::MzML,
                "mgf" => MassSpectrometryFormat::MGF,
                #[cfg(feature = "mzmlb")]
                "mzmlb" => MassSpectrometryFormat::MzMLb,
                #[cfg(feature = "thermo")]
                "raw" => MassSpectrometryFormat::ThermoRaw,
                _ => MassSpectrometryFormat::Unknown,
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
pub fn infer_from_stream<R: Read + Seek>(
    stream: &mut R,
) -> io::Result<(MassSpectrometryFormat, bool)> {
    // We need to read in at least enough bytes to span a complete XML head plus the
    // end of an opening tag
    let mut buf = Vec::with_capacity(500);
    buf.resize(500, b'\0');
    let current_pos = stream.stream_position()?;
    // record how many bytes were actually read so we know the upper bound
    let bytes_read = stream.read(buf.as_mut_slice())?;
    buf.shrink_to(bytes_read);
    let is_stream_gzipped = is_gzipped(buf.as_slice());
    if is_stream_gzipped {
        let mut decompressed_buf = Vec::new();
        // In the worst case, we can't have fewer bytes than those that were read in (minus the size of the gzip header)
        // and we assume the compression ratio means we have recouped that. We read in only that many bytes
        // decompressed because the decompressor treats an incomplete segment as an error and thus using
        // io::Read::read_to_end is not an option.
        decompressed_buf.resize(bytes_read, b'\0');
        let mut decoder = GzDecoder::new(io::Cursor::new(buf));
        decoder.read_exact(&mut decompressed_buf)?;
        buf = decompressed_buf;
    }
    stream.seek(io::SeekFrom::Start(current_pos))?;

    match &buf {
        #[cfg(feature = "mzml")]
        _ if is_mzml(&buf) => Ok((MassSpectrometryFormat::MzML, is_stream_gzipped)),
        #[cfg(feature = "mgf")]
        _ if is_mgf(&buf) => Ok((MassSpectrometryFormat::MGF, is_stream_gzipped)),
        #[cfg(feature = "thermo")]
        _ if is_thermo_raw_prefix(&buf) => {
            Ok((MassSpectrometryFormat::ThermoRaw, is_stream_gzipped))
        }
        _ => Ok((MassSpectrometryFormat::Unknown, is_stream_gzipped)),
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
        }
        _ => Ok((format, is_gzipped)),
    }
}
