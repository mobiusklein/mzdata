use std::fmt::Display;
use std::fs;
use std::io::{self, prelude::*, BufReader};
use std::path::{self, Path, PathBuf};
use std::sync::mpsc::{Receiver, Sender, SyncSender};


use flate2::{bufread::GzDecoder, write::GzEncoder};
use mzpeaks::{CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak};

use crate::io::PreBufferedStream;
use crate::params::ControlledVocabulary;
#[cfg(feature = "mzmlb")]
pub use crate::{
    io::mzmlb::{MzMLbReaderType, MzMLbWriterBuilder},
    MzMLbReader,
};

use crate::io::compression::{is_gzipped, is_gzipped_extension, RestartableGzDecoder};
use crate::io::mgf::{is_mgf, MGFReaderType, MGFWriterType};
use crate::io::mzml::{is_mzml, MzMLReaderType, MzMLWriterType};
use crate::io::traits::{RandomAccessSpectrumIterator, SpectrumSource, SpectrumWriter};
use crate::meta::MSDataFileMetadata;
use crate::spectrum::bindata::{BuildArrayMapFrom, BuildFromArrayMap};
use crate::spectrum::MultiLayerSpectrum;
use crate::{MGFReader, MzMLReader, Param};

#[cfg(feature = "thermorawfilereader")]
use super::thermo::{ThermoRawReader, ThermoRawReaderType, is_thermo_raw_prefix};

use super::traits::{SeekRead, SpectrumReceiver, StreamingSpectrumIterator};

/// Mass spectrometry file formats that [`mzdata`](crate)
/// supports
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MassSpectrometryFormat {
    MGF,
    MzML,
    MzMLb,
    ThermoRaw,
    Unknown,
}

impl MassSpectrometryFormat {

    pub fn as_param(&self) -> Option<Param> {
        let p = match self {
            MassSpectrometryFormat::MGF => ControlledVocabulary::MS.const_param_ident("Mascot MGF format", 1001062),
            MassSpectrometryFormat::MzML => ControlledVocabulary::MS.const_param_ident("MzML format", 1000584),
            MassSpectrometryFormat::MzMLb => ControlledVocabulary::MS.const_param_ident("mzMLb format", 1002838),
            MassSpectrometryFormat::ThermoRaw => ControlledVocabulary::MS.const_param_ident("Thermo RAW format", 1000563),
            MassSpectrometryFormat::Unknown => return None,
        };
        Some(p.into())
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
pub fn open_file<P: Into<path::PathBuf>>(path: P) -> io::Result<Box<dyn SpectrumSource>> {
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


/// An abstraction over different ways to get a [`SpectrumSource`] from a file path,
/// buffer, or pipe.
pub enum Source<C: CentroidLike
        + Default
        + From<CentroidPeak>
        + BuildArrayMapFrom
        + BuildFromArrayMap
        + Clone
        + 'static
        + Sync
        + Send=CentroidPeak,
        D: DeconvolutedCentroidLike
        + Default
        + From<DeconvolutedPeak>
        + BuildArrayMapFrom
        + BuildFromArrayMap
        + Clone
        + Sync
        + 'static
        + Send=DeconvolutedPeak> {
    /// A concrete path on the file system
    PathLike(PathBuf),
    /// An in-memory channel of existing spectra
    Receiver(SpectrumReceiver<C, D, MultiLayerSpectrum<C, D>>),
    /// Read from Stdin
    Stdin,
    /// A thing implementing [`std::io::Read `] and [`std::io::Seek`], along with an expected format
    Reader(Box<dyn SeekRead + Send>, Option<MassSpectrometryFormat>)
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + 'static + Sync + Send, D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + Sync + 'static + Send> Source<C, D> {

    pub fn index_file_name(&self) -> Option<PathBuf> {
        match &self {
            Self::PathLike(path) => {
                if let Some(stem) = path.file_name() {
                    if let Some(parent) = path.parent() {
                        let base = parent.join(stem);
                        let name = base.with_extension("index.json");
                        return Some(name);
                    }
                }
                None
            }
            _ => None
        }
    }

    pub fn has_index_file(&self) -> bool {
        match self.index_file_name() {
            Some(path) => path.exists(),
            None => false,
        }
    }
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + 'static + Sync + Send,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + Sync + 'static + Send> From<&Path> for Source<C, D> {
    fn from(value: &Path) -> Self {
        Self::PathLike(value.into())
    }
}


impl<C: CentroidLike + Default + From<CentroidPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + 'static + Sync + Send,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + Sync + 'static + Send> From<String> for Source<C, D> {
    fn from(value: String) -> Self {
        Self::PathLike(value.into())
    }
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + 'static + Sync + Send,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + Sync + 'static + Send> From<SpectrumReceiver<C, D, MultiLayerSpectrum<C, D>>> for Source<C, D> {
    fn from(value: SpectrumReceiver<C, D, MultiLayerSpectrum<C, D>>) -> Self {
        Self::Receiver(value.into())
    }
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + 'static + Sync + Send,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + Sync + 'static + Send> From<Receiver<MultiLayerSpectrum<C, D>>> for Source<C, D> {
    fn from(value: Receiver<MultiLayerSpectrum<C, D>>) -> Self {
        Self::Receiver(value.into())
    }
}

impl<
     C: CentroidLike + Default + From<CentroidPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + 'static + Sync + Send,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + Sync + 'static + Send> From<(Box<dyn SeekRead + Send>, MassSpectrometryFormat)> for Source<C, D> {
    fn from(value: (Box<dyn SeekRead + Send>, MassSpectrometryFormat)) -> Self {
        Self::Reader(value.0, Some(value.1))
    }
}

impl<
     C: CentroidLike + Default + From<CentroidPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + 'static + Sync + Send,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + Sync + 'static + Send> From<Box<dyn SeekRead + Send>> for Source<C, D> {
    fn from(value: Box<dyn SeekRead + Send>) -> Self {
        Self::Reader(value, None)
    }
}

/// An abstraction over places to write spectra
pub enum Sink<C: CentroidLike
        + Default
        + From<CentroidPeak>
        + BuildArrayMapFrom
        + BuildFromArrayMap
        + Clone
        + 'static
        + Sync
        + Send=CentroidPeak,
        D: DeconvolutedCentroidLike
        + Default
        + From<DeconvolutedPeak>
        + BuildArrayMapFrom
        + BuildFromArrayMap
        + Clone
        + Sync
        + 'static
        + Send=DeconvolutedPeak> {
    /// A concrete path on the file system
    PathLike(PathBuf),
    /// An in-memory channel for spectra
    Sender(Sender<MultiLayerSpectrum<C, D>>),
    /// An in-memory channel for spectra
    SyncSender(SyncSender<MultiLayerSpectrum<C, D>>),
    /// A thing implementing [`std::io::Write `], along with an expected format
    Writer(Box<dyn io::Write + Send>, MassSpectrometryFormat)
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + 'static + Sync + Send,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + Sync + 'static + Send>
     From<(Box<dyn io::Write + Send>, MassSpectrometryFormat)> for Sink<C, D> {
    fn from(value: (Box<dyn io::Write + Send>, MassSpectrometryFormat)) -> Self {
        Self::Writer(value.0, value.1)
    }
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + 'static + Sync + Send,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + Sync + 'static + Send> From<&Path> for Sink<C, D> {
    fn from(value: &Path) -> Self {
        Self::PathLike(value.into())
    }
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + 'static + Sync + Send,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + Sync + 'static + Send> From<String> for Sink<C, D> {
    fn from(value: String) -> Self {
        Self::PathLike(value.into())
    }
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + 'static + Sync + Send,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + Sync + 'static + Send> From<Sender<MultiLayerSpectrum<C, D>>> for Sink<C, D> {
    fn from(value: Sender<MultiLayerSpectrum<C, D>>) -> Self {
        Self::Sender(value.into())
    }
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + 'static + Sync + Send,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + Sync + 'static + Send> From<SyncSender<MultiLayerSpectrum<C, D>>> for Sink<C, D> {
    fn from(value: SyncSender<MultiLayerSpectrum<C, D>>) -> Self {
        Self::SyncSender(value.into())
    }
}


/// Encapsulate the read-transform-write process for mass spectrometry data sources.
///
/// This trait handles all the gory details of file format inference with [`open_reader`](MassSpectrometryReadWriteProcess::open_reader)
/// and [`open_writer`](MassSpectrometryReadWriteProcess::open_writer), leaving open the chance to customize those objects after their
/// creation in [`transform_reader`](MassSpectrometryReadWriteProcess::transform_reader) and [`transform_writer`](MassSpectrometryReadWriteProcess::transform_writer) respectively.
///
/// The only function that must be implemented explicitly is [`task`](MassSpectrometryReadWriteProcess::task) which receives
/// the reader and writer, and must contain the logic to transmit one from the other
/// with whatever transformations you wish to apply between them.
pub trait MassSpectrometryReadWriteProcess<
    C: CentroidLike
        + Default
        + From<CentroidPeak>
        + BuildArrayMapFrom
        + BuildFromArrayMap
        + Clone
        + 'static
        + Sync
        + Send=CentroidPeak,
    D: DeconvolutedCentroidLike
        + Default
        + From<DeconvolutedPeak>
        + BuildArrayMapFrom
        + BuildFromArrayMap
        + Clone
        + Sync
        + 'static
        + Send=DeconvolutedPeak,
>
{
    type ErrorType: From<io::Error>;

    /// The main entry point that starts the whole system running on a reader [`Source`]
    /// and a writer [`Sink`], or equivalent objects.
    ///
    /// By default this just invokes [`MassSpectrometryReadWriteProcess::open_reader`], but if any additional
    /// configuration needs to be done before that happens, it can be done here.
    /// Examples include creating a thread pool, temporary files or directories,
    /// or some other scoped activity.
    fn main<P: Into<Source<C, D>>, Q: Into<Sink<C, D>>>(
        &self,
        read_path: P,
        write_path: Q,
    ) -> Result<(), Self::ErrorType> {
        self.open_reader(read_path, write_path)
    }

    /// Opens the reader, transforms it with [`MassSpectrometryReadWriteProcess::transform_reader`], and then passes control to [`MassSpectrometryReadWriteProcess::open_writer`]
    fn open_reader<P: Into<Source<C, D>>, Q: Into<Sink<C, D>>>(
        &self,
        read_path: P,
        write_path: Q,
    ) -> Result<(), Self::ErrorType> {
        let read_path = read_path.into();
        match read_path {
            Source::PathLike(read_path) => {
                let (format, is_gzipped) = infer_format(&read_path)?;
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
            },
            Source::Reader(mut handle, format) => {
                let (format, _is_gzipped) = if let Some(format) = format {
                    (format, false)
                } else {
                    infer_from_stream(&mut handle)?
                };
                match format {
                    MassSpectrometryFormat::MGF => {
                        let handle = io::BufReader::new(handle);
                        let reader = MGFReaderType::new_indexed(handle);
                        let reader = self.transform_reader(reader, format)?;
                        self.open_writer(reader, format, write_path)?;
                        Ok(())
                    },
                    MassSpectrometryFormat::MzML => {
                        let handle = io::BufReader::new(handle);

                        let reader = MzMLReaderType::new_indexed(handle);
                        let reader = self.transform_reader(reader, format)?;
                        self.open_writer(reader, format, write_path)?;

                        Ok(())
                    },
                    _ => Err(io::Error::new(
                        io::ErrorKind::Unsupported,
                        format!(
                            "Input file format for {:?} not supported",
                            format
                        ),
                    )
                    .into()),
                }
            },
            Source::Receiver(receiver) => {
                let reader = StreamingSpectrumIterator::new(receiver);
                let reader = self.transform_reader(reader, MassSpectrometryFormat::Unknown)?;
                self.open_writer(reader, MassSpectrometryFormat::Unknown, write_path)?;
                Ok(())
            },
            Source::Stdin => {
                let mut buffered =
                    PreBufferedStream::new_with_buffer_size(io::stdin(), 2usize.pow(20))?;
                let (ms_format, compressed) = infer_from_stream(&mut buffered)?;
                log::debug!("Detected {ms_format:?} from STDIN (compressed? {compressed})");
                match ms_format {
                    MassSpectrometryFormat::MGF => {
                        if compressed {
                            let reader = StreamingSpectrumIterator::new(MGFReaderType::new(
                                RestartableGzDecoder::new(io::BufReader::new(buffered)),
                            ));
                            let reader = self.transform_reader(reader, ms_format)?;
                            self.open_writer(reader, ms_format, write_path)?;
                        } else {
                            let reader = StreamingSpectrumIterator::new(MGFReaderType::new(buffered));
                            let reader = self.transform_reader(reader, ms_format)?;
                            self.open_writer(reader, ms_format, write_path)?;
                        }
                        Ok(())
                    }
                    MassSpectrometryFormat::MzML => {
                        if compressed {
                            let reader = StreamingSpectrumIterator::new(MzMLReaderType::new(
                                RestartableGzDecoder::new(io::BufReader::new(buffered)),
                            ));
                            let reader = self.transform_reader(reader, ms_format)?;
                            self.open_writer(reader, ms_format, write_path)?;
                        } else {
                            let reader = StreamingSpectrumIterator::new(MzMLReaderType::new(buffered));
                            let reader = self.transform_reader(reader, ms_format)?;
                            self.open_writer(reader, ms_format, write_path)?;
                        }
                        Ok(())
                    }
                    _ => {
                        return Err(io::Error::new(
                            io::ErrorKind::Unsupported,
                            "{ms_format:?} format is not supported over Stdin",
                        ).into())
                    }
                }
            },
        }
    }

    /// Opens the writer, transforms it with [`MassSpectrometryReadWriteProcess::transform_writer`], and then passes control to [`MassSpectrometryReadWriteProcess::task`]
    fn open_writer<
        Q: Into<Sink<C, D>>,
        R: RandomAccessSpectrumIterator<C, D> + MSDataFileMetadata + SpectrumSource<C, D> + Send + 'static,
    >(
        &self,
        reader: R,
        reader_format: MassSpectrometryFormat,
        write_path: Q,
    ) -> Result<(), Self::ErrorType> {
        let write_path = write_path.into();

        match write_path {
            Sink::Sender(writer) =>  {
                return self.task(reader, writer)
            },
            Sink::SyncSender(writer) => {
                return self.task(reader, writer)
            },
            Sink::PathLike(write_path) => {
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
            },
            Sink::Writer(handle, writer_format) => {
                match writer_format {
                    MassSpectrometryFormat::MGF => {
                        let handle = io::BufWriter::new(handle);
                        let mut writer = MGFWriterType::new(
                            handle,
                        );
                        writer.copy_metadata_from(&reader);
                        let (reader, writer) =
                            self.transform_writer(reader, reader_format, writer, writer_format)?;
                        self.task(reader, writer)?;

                        Ok(())
                    }
                    MassSpectrometryFormat::MzML => {
                        let handle = io::BufWriter::new(handle);
                        let mut writer = MzMLWriterType::new(
                            handle,
                        );
                        writer.copy_metadata_from(&reader);
                        let (reader, writer) =
                            self.transform_writer(reader, reader_format, writer, writer_format)?;
                        self.task(reader, writer)?;
                        Ok(())
                    }
                    _ => {
                        Err(io::Error::new(
                                io::ErrorKind::Unsupported,
                                format!(
                                    "Output file format for {:?} not supported",
                                    writer_format
                                ),
                            )
                            .into())
                    }
                }
            }
        }
    }

    /// Customize the reader in some way. The format is passed along to allow each format
    /// to be customized explicitly.
    ///
    /// A no-op by default.
    #[allow(unused)]
    fn transform_reader<
        R: RandomAccessSpectrumIterator<C, D> + MSDataFileMetadata + SpectrumSource<C, D> + Send + 'static,
    >(
        &self,
        reader: R,
        format: MassSpectrometryFormat,
    ) -> Result<R, Self::ErrorType> {
        Ok(reader)
    }

    /// Customize the writer in some way. The format is passed along to allow each format
    /// to be customized explicitly, and the reader is provided side-by-side to permit additional
    /// information to be used.
    ///
    /// A no-op by default.
    ///
    /// # Note
    /// The caller already invokes [`MSDataFileMetadata::copy_metadata_from`]
    #[allow(unused)]
    fn transform_writer<
        R: RandomAccessSpectrumIterator<C, D> + MSDataFileMetadata + SpectrumSource<C, D> + Send + 'static,
        W: SpectrumWriter<C, D> + MSDataFileMetadata + Send + 'static,
    >(
        &self,
        reader: R,
        reader_format: MassSpectrometryFormat,
        writer: W,
        writer_format: MassSpectrometryFormat,
    ) -> Result<(R, W), Self::ErrorType> {
        Ok((reader, writer))
    }

    /// The place where the work happens to transmit data from `reader` to `writer` with whatever transformations
    /// need to take place.
    fn task<
        R: RandomAccessSpectrumIterator<C, D> + MSDataFileMetadata + SpectrumSource<C, D> + Send + 'static,
        W: SpectrumWriter<C, D> + Send + 'static,
    >(
        &self,
        reader: R,
        writer: W,
    ) -> Result<(), Self::ErrorType>;
}


/// A macro that dynamically works out how to get a [`SpectrumSource`]-derived object
/// from a path or [`io::Read`](std::io::Read) + [`io::Seek`](std::io::Seek) boxed object.
/// This is meant to be a convenience for working with a scoped file reader
/// without penalty.
///
/// `$source` is coerced into a [`Source`] which the macro in turn probes to determine
/// the appropriate file reading type. Unlike [`open_file`], this macro does not actually
/// return the reading type behind an opaque `Box<dyn SpectrumSource>`, but lets you interact
/// with the concrete type intersection without concern with object safety in an anonymous closure:
///
/// ```
/// # use std::io;
/// # use mzdata::prelude::*;
/// # use mzdata::Spectrum;
/// # fn main() -> io::Result<()> {
/// let spectra: Vec<Spectrum> = mzdata::mz_read!("./test/data/small.mzML".as_ref(), reader => { reader.collect() })?;
/// # Ok(())
/// # }
/// ```
/// The closure will return a `std::io::Result` whose success value is inferred from context. The
/// reader's lifetime is bound to the closure, and cannot be extracted without substantial type system
/// torture.
///
/// If you want to use peak types *other than* the simple defaults, pass them as additional parameters after
/// the closure.
#[macro_export]
macro_rules! mz_read {
    ($source:expr, $reader:ident => $impl:tt) => {
        $crate::mz_read!($source, $reader => $impl, mzpeaks::CentroidPeak, mzpeaks::DeconvolutedPeak)
    };
    ($source:expr, $reader:ident => $impl:tt, $C:ty, $D:ty) => {{
        let source = $crate::io::Source::<_, _>::from($source);
        match source {
            $crate::io::Source::PathLike(read_path) => {
                let (format, is_gzipped) = $crate::io::infer_format(&read_path)?;
                match format {
                    $crate::io::MassSpectrometryFormat::MGF => {
                        let handle = std::fs::File::open(read_path)?;
                        if is_gzipped {
                            let fh = $crate::io::RestartableGzDecoder::new(std::io::BufReader::new(handle));
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new($crate::io::mgf::MGFReaderType::<_, $C, $D>::new(fh));
                            Ok($impl)
                        } else {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::mgf::MGFReaderType<_, $C, $D> = $crate::io::mgf::MGFReaderType::new_indexed(handle);
                            Ok($impl)
                        }
                    }
                    $crate::io::MassSpectrometryFormat::MzML => {
                        let handle = std::fs::File::open(read_path)?;

                        if is_gzipped {
                            let fh = $crate::io::RestartableGzDecoder::new(std::io::BufReader::new(handle));
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new($crate::io::mzml::MzMLReaderType::<_, $C, $D>::new(fh));
                            Ok($impl)
                        } else {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::mzml::MzMLReaderType<_, $C, $D> = $crate::io::mzml::MzMLReaderType::<_, $C, $D>::new_indexed(handle);
                            Ok($impl)
                        }
                    }
                    #[cfg(feature = "mzmlb")]
                    $crate::io::MassSpectrometryFormat::MzMLb => {
                        #[allow(unused_mut)]
                        let mut $reader: $crate::io::mzmlb::MzMLbReaderType<$C, $D> = $crate::io::mzmlb::MzMLbReaderType::<$C, $D>::new(&read_path)?;
                        Ok($impl)
                    },
                    #[cfg(feature = "thermorawfilereader")]
                    $crate::io::MassSpectrometryFormat::ThermoRaw => {
                        #[allow(unused_mut)]
                        let mut $reader: $crate::io::thermo::ThermoRawReaderType<$C, $D> = $crate::io::thermo::ThermoRawReaderType::<$C, $D>::new(&read_path)?;
                        Ok($impl)
                    },
                    _ => Err(std::io::Error::new(
                        std::io::ErrorKind::Unsupported,
                        format!(
                            "Input file format for {} not supported",
                            read_path.display()
                        ),
                    )),
                }
            },
            $crate::io::Source::Reader(mut handle, format) => {
                let (format, is_gzipped) = if let Some(format) = format { (format, false) } else { $crate::io::infer_from_stream(&mut handle)? };
                match format {
                    $crate::io::MassSpectrometryFormat::MGF => {
                        let handle = std::io::BufReader::new(handle);
                        #[allow(unused_mut)]
                        if is_gzipped {
                            let fh = $crate::io::RestartableGzDecoder::new(std::io::BufReader::new(handle));
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new($crate::io::mgf::MGFReaderType::<_, $C, $D>::new(fh));
                            Ok($impl)
                        } else {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::mgf::MGFReaderType<_, $C, $D> = $crate::io::mgf::MGFReaderType::new_indexed(handle);
                            Ok($impl)
                        }
                    },
                    $crate::io::MassSpectrometryFormat::MzML => {
                        let handle = std::io::BufReader::new(handle);
                        #[allow(unused_mut)]
                        if is_gzipped {
                            let fh = $crate::io::RestartableGzDecoder::new(std::io::BufReader::new(handle));
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new($crate::io::mzml::MzMLReaderType::<_, $C, $D>::new(fh));
                            Ok($impl)
                        } else {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::mzml::MzMLReaderType<_, $C, $D> = $crate::io::mzml::MzMLReaderType::<_, $C, $D>::new_indexed(handle);
                            Ok($impl)
                        }
                    },
                    _ => Err(std::io::Error::new(
                        std::io::ErrorKind::Unsupported,
                        format!(
                            "Input file format for {:?} not supported from an io::Read",
                            format
                        ),
                    )),
                }
            },
            $crate::io::Source::Receiver(receiver) => {
                #[allow(unused_mut)]
                let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new(receiver);
                Ok($impl)
            },
            $crate::io::Source::Stdin => {
                let mut buffered =
                    $crate::io::PreBufferedStream::new_with_buffer_size(std::io::stdin(), 2usize.pow(20))?;
                let (ms_format, compressed) = $crate::io::infer_from_stream(&mut buffered)?;
                match ms_format {
                    $crate::io::MassSpectrometryFormat::MGF => {
                        if compressed {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new(
                                $crate::io::mgf::MGFReaderType::new(
                                    $crate::io::RestartableGzDecoder::new(std::io::BufReader::new(buffered)),
                            ));
                            Ok($impl)
                        } else {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new(
                                $crate::io::mgf::MGFReaderType::new(buffered));
                            Ok($impl)
                        }
                    }
                    $crate::io::MassSpectrometryFormat::MzML => {
                        if compressed {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new(
                                $crate::io::mzml::MzMLReaderType::new($crate::io::RestartableGzDecoder::new(std::io::BufReader::new(buffered)),
                            ));
                            Ok($impl)
                        } else {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new(
                                $crate::io::mzml::MzMLReaderType::new(buffered));
                            Ok($impl)
                        }
                    }
                    _ => {
                        return Err(std::io::Error::new(
                            std::io::ErrorKind::Unsupported,
                            "{ms_format:?} format is not supported over Stdin",
                        ).into())
                    }
                }
            },
        }
    }};
}

/// A macro that dynamically works out how to get a [`SpectrumWriter`](crate::io::SpectrumWriter) from a path
/// or [`io::Write`] boxed object.
///
/// `$sink` is coerced to a [`Sink`] which in turn the macro probes in order to determine how
/// to create the appropriate writer type. Unlike other uses of [`Sink`], `Sender` and `SyncSender`
/// are not supported.  It lets you interact with the concrete type intersection in an anonymous closure:
///
/// ```
/// # use std::io;
/// # use mzdata::prelude::*;
/// # fn main() -> io::Result<()> {
///     use mzdata::{mz_read, mz_write};
///     mzdata::mz_read!("./test/data/small.mzML".as_ref(), reader => {
///         mzdata::mz_write!("./tmp/test.mzML".as_ref(), writer => {
///             writer.copy_metadata_from(&reader);
///             for s in reader {
///                 writer.write_owned(s)?;
///             }
///         })?;
///     })?;
/// #   Ok(())
/// # }
/// ```
///
/// The closure will return a `std::io::Result` whose success value is inferred from context. The
/// writer's lifetime is bound to the closure, and cannot be extracted without substantial type system
/// torture.
///
/// If you want to use peak types *other than* the simple defaults, pass them as additional parameters after
/// the closure
#[macro_export]
macro_rules! mz_write {
    ($sink:expr, $writer:ident => $impl:tt) => {
        mz_write!($sink, $writer => $impl, mzpeaks::CentroidPeak, mzpeaks::DeconvolutedPeak)
    };
    ($sink:expr, $writer:ident => $impl:tt, $C:ty, $D:ty) => {{
        let sink = $crate::io::Sink::<$C, $D>::from($sink);
        match sink {
           $crate::io:: Sink::Sender(_) | $crate::io::Sink::SyncSender(_) =>  {
                Err(std::io::Error::new(std::io::ErrorKind::Unsupported, "Sender writers aren't supported by `mz_write`"))
            }
            $crate::io::Sink::PathLike(write_path) => {
                let (writer_format, is_gzip) = $crate::io::infer_from_path(&write_path);
                match writer_format {
                    $crate::io::MassSpectrometryFormat::MGF => {
                        let handle = std::io::BufWriter::new(std::fs::File::create(&write_path)?);
                        if is_gzip {
                            let handle = flate2::write::GzEncoder::new(handle, flate2::Compression::best());
                            let mut $writer: $crate::io::mgf::MGFWriterType<_, $C, $D> = $crate::io::mgf::MGFWriterType::new(
                                handle,
                            );
                            Ok($impl)
                        } else {
                            let mut $writer: $crate::io::mgf::MGFWriterType<_, $C, $D> = $crate::io::mgf::MGFWriterType::new(
                                handle,
                            );
                            Ok($impl)

                        }
                    }
                    $crate::io::MassSpectrometryFormat::MzML => {
                        let handle = std::io::BufWriter::new(std::fs::File::create(&write_path)?);
                        if is_gzip {
                            let handle = flate2::write::GzEncoder::new(handle, flate2::Compression::best());
                            let mut $writer: $crate::io::mzml::MzMLWriterType<_, $C, $D> = $crate::io::mzml::MzMLWriterType::new(
                                handle,
                            );
                            Ok($impl)
                        } else {
                            let mut $writer: $crate::io::mzml::MzMLWriterType<_, $C, $D> = $crate::io::mzml::MzMLWriterType::new(
                                handle,
                            );
                            Ok($impl)
                        }
                    }
                    #[cfg(feature = "mzmlb")]
                    $crate::io::MassSpectrometryFormat::MzMLb => {
                        let mut $writer = $crate::io::mzmlb::MzMLbWriterBuilder::<$C, $D>::new(&write_path)
                            .with_zlib_compression(9)
                            .create()?;
                        Ok($impl)
                    }
                    _ => Err(std::io::Error::new(
                        std::io::ErrorKind::Unsupported,
                        format!(
                            "Output file format {:?} for {} not supported",
                            writer_format,
                            write_path.display()
                        ),
                    )),
                }
            },
            $crate::io::Sink::Writer(handle, writer_format) => {
                match writer_format {
                    $crate::io::MassSpectrometryFormat::MGF => {
                        let handle = std::io::BufWriter::new(handle);
                        let mut $writer: $crate::io::mgf::MGFWriterType<_, $C, $D> = $crate::io::mgf::MGFWriterType::new(
                            handle,
                        );
                        Ok($impl)
                    }
                    $crate::io::MassSpectrometryFormat::MzML => {
                        let handle = std::io::BufWriter::new(handle);
                        let mut $writer: $crate::io::mzml::MzMLWriterType<_, $C, $D> = $crate::io::mzml::MzMLWriterType::new(
                            handle,
                        );
                        Ok($impl)
                    }
                    _ => {
                        Err(std::io::Error::new(
                                std::io::ErrorKind::Unsupported,
                                format!(
                                    "Output file format for {:?} not supported",
                                    writer_format
                                ),
                            ))
                    }
                }
            }
        }
    }};
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

    #[test]
    fn test_source_conv() -> io::Result<()> {
        let s = Source::<CentroidPeak, DeconvolutedPeak>::from("text/path".as_ref());
        assert!(matches!(s, Source::PathLike(_)));

        let fh = Box::new(io::BufReader::new(fs::File::open("./test/data/small.mgf")?)) as Box<dyn SeekRead + Send>;
        let rs: Source<CentroidPeak, DeconvolutedPeak> = (fh, MassSpectrometryFormat::MGF).into();
        assert!(matches!(rs, Source::Reader(_, _)));

        Ok(())
    }

    #[test]
    fn test_mz_read() -> io::Result<()> {
        let val: Vec<_> = mz_read!("./test/data/small.mzML".as_ref(), reader => { reader.collect() })?;
        assert_eq!(val.len(), 48);
        let val: Vec<_> = mz_read!("./test/data/small.mgf".as_ref(), reader => { reader.collect() })?;
        assert_eq!(val.len(), 34);
        let val = mz_read!("./test/data/small.mzML".as_ref(), reader => { reader.file_description().clone() })?;
        assert_eq!(val.source_files.len(), 1);
        Ok(())
    }

    #[test]
    fn test_mz_read_nested() -> io::Result<()> {
        mz_read!("./test/data/small.mzML".as_ref(), reader => {
            mz_read!("./test/data/small.mzML".as_ref(), reader2 => {
                assert_eq!(reader.len(), reader2.len());
            })?;
        })?;

        Ok(())
    }

    #[test]
    fn test_mz_write() -> io::Result<()> {
        let tmpdir = tempfile::tempdir()?;
        let path = tmpdir.path().join("test.mzML");
        mz_read!("./test/data/small.mzML".as_ref(), reader => {
            mz_write!(path.as_ref(), writer => {
                writer.copy_metadata_from(&reader);
                for s in reader {
                    writer.write_owned(s)?;
                }
            })?;
        })?;
        Ok(())
    }
}
