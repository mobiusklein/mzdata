#[allow(unused)]
use super::{Sink, Source, SpectrumSource};

/// A macro that dynamically works out how to get a [`SpectrumSource`]-derived object
/// from a path or [`io::Read`](std::io::Read) + [`io::Seek`](std::io::Seek) boxed object.
/// This is meant to be a convenience for working with a scoped file reader
/// without penalty.
///
/// `$source` is coerced into a [`Source`] which the macro in turn probes to determine
/// the appropriate file reading type. This lets you interact with the concrete type intersection
/// without concern with object safety in an anonymous closure:
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
        $crate::mz_read!($source, $reader => $impl, $crate::mzpeaks::CentroidPeak, $crate::mzpeaks::DeconvolutedPeak)
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
                    #[cfg(feature = "thermo")]
                    $crate::io::MassSpectrometryFormat::ThermoRaw => {
                        #[allow(unused_mut)]
                        let mut $reader: $crate::io::thermo::ThermoRawReaderType<$C, $D> = $crate::io::thermo::ThermoRawReaderType::<$C, $D>::new(&read_path)?;
                        Ok($impl)
                    },
                    #[cfg(feature = "bruker_tdf")]
                    $crate::io::MassSpectrometryFormat::BrukerTDF => {
                        #[allow(unused_mut)]
                        let mut $reader: $crate::io::tdf::TDFSpectrumReaderType<$crate::mzpeaks::feature::Feature<$crate::mzpeaks::MZ, $crate::mzpeaks::IonMobility>, $crate::mzpeaks::feature::ChargedFeature<$crate::mzpeaks::Mass, $crate::mzpeaks::IonMobility>, $C, $D> = $crate::io::tdf::TDFSpectrumReaderType::open_path(read_path)?;
                        Ok($impl)
                    }
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
/// or [`io::Write`](std::io::Write) boxed object.
///
/// `$sink` is coerced to a [`Sink`] which in turn the macro probes in order to determine how
/// to create the appropriate writer type. Unlike other uses of [`Sink`], `Sender` and `SyncSender`
/// are not supported.  It lets you interact with the concrete type intersection in an anonymous closure:
///
/// ```ignore
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
        mz_write!($sink, $writer => $impl, $crate::mzpeaks::CentroidPeak, $crate::mzpeaks::DeconvolutedPeak)
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
    use crate::prelude::*;
    use std::io;

    #[test]
    fn test_mz_read() -> io::Result<()> {
        let val: Vec<_> =
            mz_read!("./test/data/small.mzML".as_ref(), reader => { reader.collect() })?;
        assert_eq!(val.len(), 48);
        let val: Vec<_> =
            mz_read!("./test/data/small.mgf".as_ref(), reader => { reader.collect() })?;
        assert_eq!(val.len(), 35);
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
