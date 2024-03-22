use std::fs;
use std::io;
use std::io::BufReader;
use std::path;

use std::io::prelude::*;
use std::path::PathBuf;

use crate::spectrum::bindata::{BuildArrayMapFrom, BuildFromArrayMap};
use crate::MGFReader;
use crate::MzMLReader;
use flate2::{bufread::GzDecoder, write::GzEncoder};
use mzpeaks::CentroidLike;
use mzpeaks::CentroidPeak;
use mzpeaks::DeconvolutedCentroidLike;
use mzpeaks::DeconvolutedPeak;

#[cfg(feature = "mzmlb")]
pub use crate::{
    io::mzmlb::{MzMLbReaderType, MzMLbWriterBuilder},
    MzMLbReader,
};

use crate::io::compression::{is_gzipped, is_gzipped_extension};
use crate::io::mgf::{is_mgf, MGFReaderType, MGFWriterType};
use crate::io::mzml::{is_mzml, MzMLReaderType, MzMLWriterType};
use crate::io::traits::{RandomAccessSpectrumIterator, ScanSource, ScanWriter};
use crate::meta::MSDataFileMetadata;

#[cfg(feature = "thermorawfilereader")]
use super::thermo::{ThermoRawReader, ThermoRawReaderType, is_thermo_raw_prefix};

use super::RestartableGzDecoder;
use super::StreamingSpectrumIterator;

/// Mass spectrometry file formats that [`mzdata`](crate)
/// supports
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MassSpectrometryFormat {
    MGF,
    MzML,
    #[cfg(feature = "mzmlb")]
    MzMLb,
    #[cfg(feature = "thermorawfilereader")]
    ThermoRaw,
    Unknown,
}

/// Given a path, infer the file format and whether or not the file at that path is
/// GZIP compressed
pub fn infer_from_path<P: Into<path::PathBuf>>(path: P) -> (MassSpectrometryFormat, bool) {
    let path: path::PathBuf = path.into();
    let (is_gzipped, path) = is_gzipped_extension(path);
    if let Some(ext) = path.extension() {
        if let Some(ext) = ext.to_ascii_lowercase().to_str() {
            let form = match ext {
                "mzml" => MassSpectrometryFormat::MzML,
                "mgf" => MassSpectrometryFormat::MGF,
                #[cfg(feature = "mzmlb")]
                "mzmlb" => MassSpectrometryFormat::MzMLb,
                #[cfg(feature = "thermorawfilereader")]
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
        decoder.read(&mut decompressed_buf)?;
        buf = decompressed_buf;
    }
    stream.seek(io::SeekFrom::Start(current_pos))?;

    match &buf {
        _ if is_mzml(&buf) => Ok((MassSpectrometryFormat::MzML, is_stream_gzipped)),
        _ if is_mgf(&buf) => Ok((MassSpectrometryFormat::MGF, is_stream_gzipped)),
        #[cfg(feature = "thermorawfilereader")]
        _ if is_thermo_raw_prefix(&buf) => Ok((MassSpectrometryFormat::ThermoRaw, is_stream_gzipped)),
        _ => Ok((MassSpectrometryFormat::Unknown, is_stream_gzipped))
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

/// Given a local file system path, infer the file format, and attempt to open it
/// for reading.
pub fn open_file<P: Into<path::PathBuf>>(path: P) -> io::Result<Box<dyn ScanSource>> {
    let path = path.into();
    let (format, is_gzipped) = infer_format(path.clone())?;

    if is_gzipped {
        Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "Gzipped files are not supported",
        ))
    } else {
        match format {
            MassSpectrometryFormat::MGF => {
                let handle = fs::File::open(path)?;
                let reader = MGFReader::new_indexed(handle);
                Ok(Box::new(reader))
            }
            MassSpectrometryFormat::MzML => {
                let handle = fs::File::open(path)?;
                let reader = MzMLReader::new_indexed(handle);
                Ok(Box::new(reader))
            }
            #[cfg(feature = "thermorawfilereader")]
            MassSpectrometryFormat::ThermoRaw => {
                let reader = ThermoRawReader::new(&path)?;
                Ok(Box::new(reader))
            }
            #[cfg(feature = "mzmlb")]
            MassSpectrometryFormat::MzMLb => {
                let reader = MzMLbReader::new(&path)?;
                Ok(Box::new(reader))
            }
            _ => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "File format not supported",
            )),
        }
    }
}

pub trait MassSpectrometryReadWriteProcess<
    C: CentroidLike
        + Default
        + From<CentroidPeak>
        + BuildArrayMapFrom
        + BuildFromArrayMap
        + 'static
        + Send=CentroidPeak,
    D: DeconvolutedCentroidLike
        + Default
        + From<DeconvolutedPeak>
        + BuildArrayMapFrom
        + BuildFromArrayMap
        + 'static
        + Send=DeconvolutedPeak,
>
{
    type ErrorType: From<io::Error>;

    fn open_reader<P: Into<path::PathBuf>, Q: Into<path::PathBuf>>(
        &self,
        read_path: P,
        write_path: Q,
    ) -> Result<(), Self::ErrorType> {
        let read_path = read_path.into();
        let (format, is_gzipped) = infer_format(read_path.clone())?;
        {
            match format {
                MassSpectrometryFormat::MGF => {
                    let handle = fs::File::open(read_path)?;
                    if is_gzipped {
                        let fh = RestartableGzDecoder::new(io::BufReader::new(handle));
                        let reader = StreamingSpectrumIterator::new(MGFReaderType::new(fh));
                        let reader = self.transform_reader(reader, format)?;
                        self.open_writer(reader, format, write_path)?;
                    } else {
                        let reader = MGFReaderType::new_indexed(handle);
                        let reader = self.transform_reader(reader, format)?;
                        self.open_writer(reader, format, write_path)?;
                    };
                    Ok(())
                }
                MassSpectrometryFormat::MzML => {
                    let handle = fs::File::open(read_path)?;

                    if is_gzipped {
                        let fh = RestartableGzDecoder::new(io::BufReader::new(handle));
                        let reader = StreamingSpectrumIterator::new(MzMLReaderType::new(fh));
                        let reader = self.transform_reader(reader, format)?;
                        self.open_writer(reader, format, write_path)?;
                    } else {
                        let reader = MzMLReaderType::new_indexed(handle);
                        let reader = self.transform_reader(reader, format)?;
                        self.open_writer(reader, format, write_path)?;
                    };
                    Ok(())
                }
                #[cfg(feature = "mzmlb")]
                MassSpectrometryFormat::MzMLb => {
                    let reader = MzMLbReaderType::new(&read_path)?;
                    let reader = self.transform_reader(reader, format)?;
                    self.open_writer(reader, format, write_path)?;
                    Ok(())
                },
                #[cfg(feature = "thermorawfilereader")]
                MassSpectrometryFormat::ThermoRaw => {
                    let reader = ThermoRawReaderType::new(&read_path)?;
                    let reader = self.transform_reader(reader, format)?;
                    self.open_writer(reader, format, write_path)?;
                    Ok(())
                },
                _ => Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    format!(
                        "Input file format for {} not supported",
                        read_path.display()
                    ),
                )
                .into()),
            }
        }
    }

    fn open_writer<
        Q: Into<path::PathBuf>,
        R: RandomAccessSpectrumIterator<C, D> + MSDataFileMetadata + ScanSource<C, D> + Send + 'static,
    >(
        &self,
        reader: R,
        reader_format: MassSpectrometryFormat,
        write_path: Q,
    ) -> Result<(), Self::ErrorType> {
        let write_path: PathBuf = write_path.into();
        let (writer_format, is_gzip) = infer_from_path(&write_path);
        match writer_format {
            MassSpectrometryFormat::MGF => {
                let handle = io::BufWriter::new(fs::File::create(&write_path)?);
                if is_gzip {
                    let handle = GzEncoder::new(handle, flate2::Compression::best());
                    let mut writer = MGFWriterType::new(
                        handle,
                    );
                    writer.copy_metadata_from(&reader);
                    let (reader, writer) =
                        self.transform_writer(reader, reader_format, writer, writer_format)?;
                    self.task(reader, writer)?;
                } else {
                    let mut writer = MGFWriterType::new(
                        handle,
                    );
                    writer.copy_metadata_from(&reader);
                    let (reader, writer) =
                        self.transform_writer(reader, reader_format, writer, writer_format)?;
                    self.task(reader, writer)?;
                }
                Ok(())
            }
            MassSpectrometryFormat::MzML => {
                let handle = io::BufWriter::new(fs::File::create(&write_path)?);
                if is_gzip {
                    let handle = GzEncoder::new(handle, flate2::Compression::best());
                    let mut writer = MzMLWriterType::new(
                        handle,
                    );
                    writer.copy_metadata_from(&reader);
                    let (reader, writer) =
                        self.transform_writer(reader, reader_format, writer, writer_format)?;
                    self.task(reader, writer)?;
                } else {
                    let mut writer = MzMLWriterType::new(
                        handle,
                    );
                    writer.copy_metadata_from(&reader);
                    let (reader, writer) =
                        self.transform_writer(reader, reader_format, writer, writer_format)?;
                    self.task(reader, writer)?;
                }
                Ok(())
            }
            #[cfg(feature = "mzmlb")]
            MassSpectrometryFormat::MzMLb => {
                let mut writer = MzMLbWriterBuilder::<C, D>::new(&write_path)
                    .with_zlib_compression(9)
                    .create()?;
                writer.copy_metadata_from(&reader);
                let (reader, writer) =
                    self.transform_writer(reader, reader_format, writer, writer_format)?;
                self.task(reader, writer)?;
                Ok(())
            }
            _ => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                format!(
                    "Output file format for {} not supported",
                    write_path.display()
                ),
            )
            .into()),
        }
    }

    #[allow(unused)]
    fn transform_reader<
        R: RandomAccessSpectrumIterator<C, D> + ScanSource<C, D> + Send + 'static,
    >(
        &self,
        reader: R,
        format: MassSpectrometryFormat,
    ) -> Result<R, Self::ErrorType> {
        Ok(reader)
    }

    #[allow(unused)]
    fn transform_writer<
        R: RandomAccessSpectrumIterator<C, D> + ScanSource<C, D> + Send + 'static,
        W: ScanWriter<C, D> + Send + 'static,
    >(
        &self,
        reader: R,
        reader_format: MassSpectrometryFormat,
        writer: W,
        writer_format: MassSpectrometryFormat,
    ) -> Result<(R, W), Self::ErrorType> {
        Ok((reader, writer))
    }

    fn task<
        R: RandomAccessSpectrumIterator<C, D> + ScanSource<C, D> + Send + 'static,
        W: ScanWriter<C, D> + Send + 'static,
    >(
        &self,
        reader: R,
        writer: W,
    ) -> Result<(), Self::ErrorType>;
}

#[cfg(test)]
mod test {
    use crate::{
        prelude::*,
        spectrum::{ArrayType, Spectrum},
    };

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
