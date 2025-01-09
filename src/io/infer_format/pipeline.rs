use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{Receiver, Sender, SyncSender};
use std::any::Any;

use flate2::write::GzEncoder;
use mzpeaks::{CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak};
#[cfg(feature = "bruker_tdf")]
use mzpeaks::{feature::{ChargedFeature, Feature}, IonMobility, Mass, MZ};

use crate::io::PreBufferedStream;
#[cfg(feature = "mzmlb")]
pub use crate::io::mzmlb::{MzMLbReaderType, MzMLbWriterBuilder};

use crate::io::compression::RestartableGzDecoder;
use crate::io::mgf::{MGFReaderType, MGFWriterType};
use crate::io::mzml::{MzMLReaderType, MzMLWriterType};
use crate::io::traits::{RandomAccessSpectrumIterator, SpectrumSource, SpectrumWriter, MZFileReader};
use crate::io::SpectrumReceiver;
use crate::io::StreamingSpectrumIterator;
use crate::meta::MSDataFileMetadata;
use crate::prelude::*;
use crate::spectrum::bindata::{BuildArrayMapFrom, BuildFromArrayMap};
use crate::spectrum::MultiLayerSpectrum;

#[cfg(feature = "thermo")]
use crate::io::thermo::ThermoRawReaderType;

#[cfg(feature = "bruker_tdf")]
use crate::io::tdf::TDFSpectrumReaderType;

use super::infer_format;
use super::infer_from_path;
use super::infer_from_stream;
use super::MassSpectrometryFormat;



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
        Self::Receiver(value)
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
        Self::Sender(value)
    }
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + 'static + Sync + Send,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildArrayMapFrom + BuildFromArrayMap + Clone + Sync + 'static + Send> From<SyncSender<MultiLayerSpectrum<C, D>>> for Sink<C, D> {
    fn from(value: SyncSender<MultiLayerSpectrum<C, D>>) -> Self {
        Self::SyncSender(value)
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
                    #[cfg(feature = "thermo")]
                    MassSpectrometryFormat::ThermoRaw => {
                        let reader = ThermoRawReaderType::new(&read_path)?;
                        let reader = self.transform_reader(reader, format)?;
                        self.open_writer(reader, format, write_path)?;
                        Ok(())
                    },
                    #[cfg(feature = "bruker_tdf")]
                    MassSpectrometryFormat::BrukerTDF => {
                        let reader: TDFSpectrumReaderType<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>, C, D> = TDFSpectrumReaderType::open_path(read_path)?;
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
                        Err(io::Error::new(
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
        R: RandomAccessSpectrumIterator<C, D> + MSDataFileMetadata + SpectrumSource<C, D> + Send + Any + 'static,
    >(
        &self,
        reader: R,
        reader_format: MassSpectrometryFormat,
        write_path: Q,
    ) -> Result<(), Self::ErrorType> {
        let write_path = write_path.into();

        match write_path {
            Sink::Sender(writer) =>  {
                self.task(reader, writer)
            },
            Sink::SyncSender(writer) => {
                self.task(reader, writer)
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
        R: RandomAccessSpectrumIterator<C, D> + MSDataFileMetadata + SpectrumSource<C, D> + Send + Any + 'static,
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
        R: RandomAccessSpectrumIterator<C, D> + MSDataFileMetadata + SpectrumSource<C, D> + Any + Send + 'static,
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
        R: RandomAccessSpectrumIterator<C, D> + MSDataFileMetadata + SpectrumSource<C, D> + Send + Any + 'static,
        W: SpectrumWriter<C, D> + Send + Any + 'static,
    >(
        &self,
        reader: R,
        writer: W,
    ) -> Result<(), Self::ErrorType>;
}

