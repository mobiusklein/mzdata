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

#[cfg(feature = "imzml")]
use crate::io::imzml::is_imzml;

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
    IMzML,
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
            MassSpectrometryFormat::IMzML => {
                ControlledVocabulary::MS.const_param_ident("imzML format", 1003577)
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
                #[cfg(feature = "mzml")]
                "mzml" => MassSpectrometryFormat::MzML,
                #[cfg(feature = "mgf")]
                "mgf" => MassSpectrometryFormat::MGF,
                #[cfg(feature = "mzmlb")]
                "mzmlb" => MassSpectrometryFormat::MzMLb,
                #[cfg(feature = "thermo")]
                "raw" => MassSpectrometryFormat::ThermoRaw,
                #[cfg(feature = "imzml")]
                "imzml" => MassSpectrometryFormat::IMzML,
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
/// The stream position is restored even if inference fails. If restoring it
/// fails, the seek error is returned. Only the prefix is inspected, not the
/// complete file or gzip checksum. Short reads may require further input before
/// inference returns.
/// Inference uses up to 500 decoded bytes. Gzip decoding reads at most 64 KiB
/// of compressed input.
pub fn infer_from_stream<R: Read + Seek>(
    stream: &mut R,
) -> io::Result<(MassSpectrometryFormat, bool)> {
    let current_pos = stream.stream_position()?;
    let result = (|| {
        let mut buf = vec![0u8; 500];
        let mut n = 0;
        while n < 2 {
            match stream.read(&mut buf[n..]) {
                Ok(0) => return Ok((MassSpectrometryFormat::Unknown, false)),
                Ok(read) => n += read,
                Err(err) if err.kind() == io::ErrorKind::Interrupted => continue,
                Err(err) => return Err(err),
            }
        }
        let gzipped = is_gzipped(&buf[..n]);
        let format = if gzipped {
            stream.seek(io::SeekFrom::Start(current_pos))?;
            let source = (&mut *stream).take(2u64.pow(16));
            let mut decoder = GzDecoder::new(BufReader::with_capacity(buf.len(), source));
            infer_from_prefix(&mut decoder, &mut buf, 0)?
        } else {
            infer_from_prefix(stream, &mut buf, n)?
        };
        Ok((format, gzipped))
    })();
    stream.seek(io::SeekFrom::Start(current_pos))?;
    result
}

fn infer_from_prefix<R: Read>(
    stream: &mut R,
    buf: &mut [u8],
    mut n: usize,
) -> io::Result<MassSpectrometryFormat> {
    loop {
        let format = match &buf[..n] {
            [] => MassSpectrometryFormat::Unknown,
            #[cfg(feature = "imzml")]
            _ if is_imzml(&buf[..n]) => MassSpectrometryFormat::IMzML,
            #[cfg(feature = "mzml")]
            _ if is_mzml(&buf[..n]) => MassSpectrometryFormat::MzML,
            #[cfg(feature = "mgf")]
            _ if is_mgf(&buf[..n]) => MassSpectrometryFormat::MGF,
            #[cfg(feature = "thermo")]
            _ if is_thermo_raw_prefix(&buf[..n]) => MassSpectrometryFormat::ThermoRaw,
            _ => MassSpectrometryFormat::Unknown,
        };
        // An XML comment can contain the MGF signature before the mzML opening tag.
        let ambiguous_mgf = cfg!(feature = "mzml")
            && format == MassSpectrometryFormat::MGF
            && buf[..n]
                .strip_prefix(b"\xef\xbb\xbf")
                .unwrap_or(&buf[..n])
                .iter()
                .find(|byte| !byte.is_ascii_whitespace())
                == Some(&b'<');
        // imzML declares its IMS vocabulary after the mzML opening tag.
        let identified = format != MassSpectrometryFormat::Unknown
            && !(cfg!(feature = "imzml") && format == MassSpectrometryFormat::MzML)
            && !ambiguous_mgf;
        if identified || n == buf.len() {
            return Ok(format);
        }
        match stream.read(&mut buf[n..]) {
            Ok(0) => return Ok(format),
            Ok(read) => n += read,
            Err(err) if err.kind() == io::ErrorKind::Interrupted => continue,
            Err(err) => return Err(err),
        }
    }
}

/// Given a path, infer the file format and whether or not the file at that path is
/// GZIP compressed, using both the file name and by trying to open and read the file
/// header
pub fn infer_format<P: Into<path::PathBuf>>(path: P) -> io::Result<(MassSpectrometryFormat, bool)> {
    let path: path::PathBuf = path.into();

    let (format, is_gzipped) = infer_from_path(&path);
    log::debug!("Inferred format from path: {:?} (gzip: {})", format, is_gzipped);
    match format {
        MassSpectrometryFormat::Unknown => {
            // If the path is a directory, don't try to open it as a file
            if path.is_dir() {
                Ok((MassSpectrometryFormat::Unknown, false))
            } else {
                let handle = fs::File::open(path.clone())?;
                let mut stream = BufReader::new(handle);
                let (format, is_gzipped) = infer_from_stream(&mut stream)?;
                log::debug!(
                    "Inferred format from stream: {:?} (gzip: {})",
                    format,
                    is_gzipped
                );
                Ok((format, is_gzipped))
            }
        }
        _ => Ok((format, is_gzipped)),
    }
}


pub trait _SourceFileExt: Sized {

    fn _new_path(name: String, location: String, file_format: Option<Param>) -> Self;

    /// Create a new [`SourceFile`] from a path.
    ///
    /// This function makes a minimal effort to infer information about the file,
    /// using [`infer_format`] to populate [`SourceFile::file_format`]
    fn from_path<P: AsRef<path::Path>>(path: P) -> io::Result<Self> {
        let path = path.as_ref();
        let format = infer_format(path)
            .ok()
            .and_then(|(format, _)| format.as_param());
        let inst = Self::_new_path(
            path
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_default(),
            path
                .canonicalize()?
                .parent()
                .map(|s| format!("file://{}", s.to_string_lossy()))
                .unwrap_or_else(|| "file://".to_string()),
            format,
        );
        Ok(inst)
    }
}

impl _SourceFileExt for crate::meta::SourceFile {
    fn _new_path(name: String, location: String, file_format: Option<Param>) -> Self {
        Self {
            name,
            location,
            file_format,
            ..Default::default()
        }
    }
}
