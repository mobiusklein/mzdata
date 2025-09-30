use std::{
    io,
    marker::PhantomData,
    path::{Path, PathBuf},
};

#[allow(unused)]
use mzpeaks::{
    feature::{ChargedFeature, Feature},
    IonMobility, Mass, MZ,
};

#[allow(unused)]
use crate::{
    prelude::*,
    spectrum::{CentroidPeakAdapting, DeconvolutedPeakAdapting},
};

#[allow(unused)]
use super::{infer_format::MZReaderType, Sink, Source, SpectrumSource};

#[doc(hidden)]
#[cfg(not(feature = "mgf"))]
pub use super::infer_format::MZReaderType as MGFReaderType;
#[doc(hidden)]
#[cfg(feature = "mgf")]
pub use super::mgf::{MGFReaderType, MGFWriterType};

#[doc(hidden)]
#[cfg(not(feature = "mgf"))]
pub type MGFWriterType<W, C, D> = DummyWriter<W, C, D>;

#[doc(hidden)]
#[cfg(feature = "mzml")]
pub use super::mzml::{MzMLReaderType, MzMLWriterType};

#[doc(hidden)]
#[cfg(not(feature = "mzml"))]
pub use super::infer_format::MZReaderType as MzMLReaderType;

#[cfg(not(feature = "mzml"))]
pub type MzMLWriterType<W, C, D> = DummyWriter<W, C, D>;

#[doc(hidden)]
#[cfg(feature = "thermo")]
pub use super::thermo::ThermoRawReaderType;

#[doc(hidden)]
#[cfg(not(feature = "thermo"))]
pub type ThermoRawReaderType<C, D> = MZReaderType<std::fs::File, C, D>;

#[doc(hidden)]
#[cfg(feature = "mzmlb")]
pub use super::mzmlb::{MzMLbReaderType, MzMLbWriterType};

#[doc(hidden)]
#[cfg(not(feature = "mzmlb"))]
pub type MzMLbReaderType<C, D> = MZReaderType<std::fs::File, C, D>;

#[doc(hidden)]
#[cfg(not(feature = "mzmlb"))]
pub type MzMLbWriterType<C, D> = DummyWriter2<C, D>;

#[doc(hidden)]
#[cfg(feature = "bruker_tdf")]
use super::tdf::TDFSpectrumReaderType as _TDFSpectrumReaderType;

#[doc(hidden)]
#[cfg(feature = "bruker_tdf")]
pub type TDFSpectrumReaderType<C, D> =
    _TDFSpectrumReaderType<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>, C, D>;

#[doc(hidden)]
#[cfg(not(feature = "bruker_tdf"))]
pub type TDFSpectrumReaderType<C, D> = MZReaderType<std::fs::File, C, D>;

pub const fn mzml_support() -> bool {
    #[cfg(feature = "mzml")]
    return true;
    #[cfg(not(feature = "mzml"))]
    return false;
}

pub const fn mgf_support() -> bool {
    #[cfg(feature = "mgf")]
    return true;
    #[cfg(not(feature = "mgf"))]
    return false;
}

pub const fn mzmlb_support() -> bool {
    #[cfg(feature = "mzmlb")]
    return true;
    #[cfg(not(feature = "mzmlb"))]
    return false;
}

pub const fn thermo_support() -> bool {
    #[cfg(feature = "thermo")]
    return true;
    #[cfg(not(feature = "thermo"))]
    return false;
}

pub const fn bruker_tdf_support() -> bool {
    #[cfg(feature = "bruker_tdf")]
    return true;
    #[cfg(not(feature = "bruker_tdf"))]
    return false;
}

#[doc(hidden)]
#[allow(unused)]
pub fn mgf_new<
    R: io::Read + io::Seek,
    C: CentroidPeakAdapting + BuildFromArrayMap,
    D: DeconvolutedPeakAdapting + BuildFromArrayMap,
>(
    handle: R,
) -> MGFReaderType<R, C, D> {
    #[cfg(feature = "mgf")]
    return MGFReaderType::new(handle);
    #[cfg(not(feature = "mgf"))]
    panic!("MGF reading is not enabled")
}

#[doc(hidden)]
#[allow(unused)]
pub fn mgf_new_indexed<
    R: io::Read + io::Seek,
    C: CentroidPeakAdapting + BuildFromArrayMap,
    D: DeconvolutedPeakAdapting + BuildFromArrayMap,
>(
    handle: R,
) -> MGFReaderType<R, C, D> {
    #[cfg(not(feature = "mgf"))]
    panic!("MGF reading is not enabled");
    #[cfg(feature = "mgf")]
    return MGFReaderType::new_indexed(handle);
}

#[doc(hidden)]
#[allow(unused)]
pub fn mzml_new<
    R: io::Read + io::Seek,
    C: CentroidPeakAdapting + BuildFromArrayMap,
    D: DeconvolutedPeakAdapting + BuildFromArrayMap,
>(
    handle: R,
) -> MzMLReaderType<R, C, D> {
    #[cfg(feature = "mzml")]
    return MzMLReaderType::new(handle);
    #[cfg(not(feature = "mzml"))]
    panic!("mzML reading is not enabled")
}

#[doc(hidden)]
#[allow(unused)]
pub fn mzml_new_indexed<
    R: io::Read + io::Seek,
    C: CentroidPeakAdapting + BuildFromArrayMap,
    D: DeconvolutedPeakAdapting + BuildFromArrayMap,
>(
    handle: R,
) -> MzMLReaderType<R, C, D> {
    #[cfg(feature = "mzml")]
    return MzMLReaderType::new_indexed(handle);
    #[cfg(not(feature = "mzml"))]
    panic!("mzML reading is not enabled")
}

#[doc(hidden)]
#[allow(unused)]
pub fn thermo_new<
    R: Into<PathBuf>,
    C: CentroidPeakAdapting + BuildFromArrayMap,
    D: DeconvolutedPeakAdapting + BuildFromArrayMap,
>(
    handle: R,
) -> io::Result<ThermoRawReaderType<C, D>> {
    #[cfg(feature = "thermo")]
    return ThermoRawReaderType::new(handle);
    #[cfg(not(feature = "thermo"))]
    panic!("Thermo RAW file reading not enabled")
}

#[doc(hidden)]
#[allow(unused)]
pub fn mzmlb_new<
    R: AsRef<Path>,
    C: CentroidPeakAdapting + BuildFromArrayMap,
    D: DeconvolutedPeakAdapting + BuildFromArrayMap,
>(
    handle: R,
) -> io::Result<MzMLbReaderType<C, D>> {
    #[cfg(feature = "mzmlb")]
    return MzMLbReaderType::new(&handle);
    #[cfg(not(feature = "mzmlb"))]
    panic!("mzMLb file reading not enabled")
}

#[doc(hidden)]
#[allow(unused)]
pub fn bruker_tdf_new<
    R: AsRef<Path>,
    C: CentroidPeakAdapting + BuildFromArrayMap,
    D: DeconvolutedPeakAdapting + BuildFromArrayMap,
>(
    handle: R,
) -> io::Result<TDFSpectrumReaderType<C, D>> {
    #[cfg(feature = "bruker_tdf")]
    return TDFSpectrumReaderType::new(handle).map_err(|e| io::Error::new(io::ErrorKind::Other, e));
    #[cfg(not(feature = "bruker_tdf"))]
    panic!("Bruker TDF file reading not enabled")
}

#[doc(hidden)]
pub struct DummyWriter<W, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>(
    PhantomData<(W, C, D)>,
);

impl<W, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MSDataFileMetadata
    for DummyWriter<W, C, D>
{
    fn data_processings(&self) -> &Vec<crate::meta::DataProcessing> {
        panic!("stub");
    }

    fn instrument_configurations(
        &self,
    ) -> &std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
        panic!("stub");
    }

    fn file_description(&self) -> &crate::meta::FileDescription {
        panic!("stub");
    }

    fn softwares(&self) -> &Vec<crate::meta::Software> {
        panic!("stub");
    }

    fn samples(&self) -> &Vec<crate::meta::Sample> {
        panic!("stub");
    }

    fn data_processings_mut(&mut self) -> &mut Vec<crate::meta::DataProcessing> {
        panic!("stub");
    }

    fn instrument_configurations_mut(
        &mut self,
    ) -> &mut std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
        panic!("stub");
    }

    fn file_description_mut(&mut self) -> &mut crate::meta::FileDescription {
        panic!("stub");
    }

    fn softwares_mut(&mut self) -> &mut Vec<crate::meta::Software> {
        panic!("stub");
    }

    fn samples_mut(&mut self) -> &mut Vec<crate::meta::Sample> {
        panic!("stub");
    }
}

impl<W, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> DummyWriter<W, C, D> {
    pub fn new(_: W) -> Self {
        Self(PhantomData)
    }
}

impl<W, C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> SpectrumWriter<C, D>
    for DummyWriter<W, C, D>
{
    #[allow(unused)]
    fn write<S: SpectrumLike<C, D> + 'static>(&mut self, spectrum: &S) -> io::Result<usize> {
        Err(io::Error::new(io::ErrorKind::Unsupported, "Dummy Writer"))
    }

    fn flush(&mut self) -> io::Result<()> {
        Err(io::Error::new(io::ErrorKind::Unsupported, "Dummy Writer"))
    }

    fn close(&mut self) -> io::Result<()> {
        Err(io::Error::new(io::ErrorKind::Unsupported, "Dummy Writer"))
    }
}

#[doc(hidden)]
pub struct DummyWriter2<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting>(PhantomData<(C, D)>);

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> MSDataFileMetadata
    for DummyWriter2<C, D>
{
    fn data_processings(&self) -> &Vec<crate::meta::DataProcessing> {
        todo!()
    }

    fn instrument_configurations(
        &self,
    ) -> &std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
        todo!()
    }

    fn file_description(&self) -> &crate::meta::FileDescription {
        todo!()
    }

    fn softwares(&self) -> &Vec<crate::meta::Software> {
        todo!()
    }

    fn samples(&self) -> &Vec<crate::meta::Sample> {
        todo!()
    }

    fn data_processings_mut(&mut self) -> &mut Vec<crate::meta::DataProcessing> {
        todo!()
    }

    fn instrument_configurations_mut(
        &mut self,
    ) -> &mut std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
        todo!()
    }

    fn file_description_mut(&mut self) -> &mut crate::meta::FileDescription {
        todo!()
    }

    fn softwares_mut(&mut self) -> &mut Vec<crate::meta::Software> {
        todo!()
    }

    fn samples_mut(&mut self) -> &mut Vec<crate::meta::Sample> {
        todo!()
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> DummyWriter2<C, D> {
    #[allow(unused)]
    pub fn new<P>(path: P) -> io::Result<Self> {
        Ok(Self(PhantomData))
    }
}

impl<C: CentroidPeakAdapting, D: DeconvolutedPeakAdapting> SpectrumWriter<C, D>
    for DummyWriter2<C, D>
{
    #[allow(unused)]
    fn write<S: SpectrumLike<C, D> + 'static>(&mut self, spectrum: &S) -> io::Result<usize> {
        Err(io::Error::new(io::ErrorKind::Unsupported, "Dummy Writer"))
    }

    fn flush(&mut self) -> io::Result<()> {
        Err(io::Error::new(io::ErrorKind::Unsupported, "Dummy Writer"))
    }

    fn close(&mut self) -> io::Result<()> {
        Err(io::Error::new(io::ErrorKind::Unsupported, "Dummy Writer"))
    }
}

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
    ($source:expr_2021, $reader:ident => $impl:tt) => {
        $crate::mz_read!($source, $reader => $impl, $crate::mzpeaks::CentroidPeak, $crate::mzpeaks::DeconvolutedPeak)
    };
    ($source:expr_2021, $reader:ident => $impl:tt, $C:ty, $D:ty) => {{
        let source = $crate::io::Source::<_, _>::from($source);
        match source {
            $crate::io::Source::PathLike(read_path) => {
                let (format, is_gzipped) = $crate::io::infer_format(&read_path)?;
                match format {
                    $crate::io::MassSpectrometryFormat::MGF if $crate::io::_impl::mgf_support() => {
                        let handle = std::fs::File::open(read_path)?;
                        if is_gzipped {
                            let fh = $crate::io::RestartableGzDecoder::new(std::io::BufReader::new(handle));
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new($crate::io::_impl::mgf_new::<_, $C, $D>(fh));
                            Ok($impl)
                        } else {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::_impl::MGFReaderType<_, $C, $D> = $crate::io::_impl::mgf_new_indexed(handle);
                            Ok($impl)
                        }
                    }
                    $crate::io::MassSpectrometryFormat::MzML if $crate::io::_impl::mzml_support() => {
                        let handle = std::fs::File::open(read_path)?;

                        if is_gzipped {
                            let fh = $crate::io::RestartableGzDecoder::new(std::io::BufReader::new(handle));
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new($crate::io::_impl::mzml_new::<_, $C, $D>(fh));
                            Ok($impl)
                        } else {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::mzml::MzMLReaderType<_, $C, $D> = $crate::io::_impl::mzml_new_indexed::<_, $C, $D>(handle);
                            Ok($impl)
                        }
                    }
                    $crate::io::MassSpectrometryFormat::MzMLb if $crate::io::_impl::mzmlb_support() => {
                        #[allow(unused_mut)]
                        let mut $reader: $crate::io::_impl::MzMLbReaderType<$C, $D> = $crate::io::_impl::mzmlb_new(&read_path)?;
                        Ok($impl)
                    },
                    $crate::io::MassSpectrometryFormat::ThermoRaw if $crate::io::_impl::thermo_support() => {
                        #[allow(unused_mut)]
                        let mut $reader: $crate::io::_impl::ThermoRawReaderType<$C, $D> = $crate::io::_impl::thermo_new::<_, $C, $D>(&read_path)?;
                        Ok($impl)
                    },
                    $crate::io::MassSpectrometryFormat::BrukerTDF if $crate::io::_impl::bruker_tdf_support() => {
                        #[allow(unused_mut)]
                        let mut $reader: $crate::io::_impl::TDFSpectrumReaderType<$C, $D> = $crate::io::_impl::bruker_tdf_new(read_path)?;
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
                    $crate::io::MassSpectrometryFormat::MGF if $crate::io::_impl::mgf_support() => {
                        let handle = std::io::BufReader::new(handle);
                        #[allow(unused_mut)]
                        if is_gzipped {
                            let fh = $crate::io::RestartableGzDecoder::new(std::io::BufReader::new(handle));
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new($crate::io::_impl::mgf_new::<_, $C, $D>(fh));
                            Ok($impl)
                        } else {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::_impl::MGFReaderType<_, $C, $D> = $crate::io::_impl::mgf_new_indexed(handle);
                            Ok($impl)
                        }
                    },
                    $crate::io::MassSpectrometryFormat::MzML if $crate::io::_impl::mzml_support() => {
                        let handle = std::io::BufReader::new(handle);
                        #[allow(unused_mut)]
                        if is_gzipped {
                            let fh = $crate::io::RestartableGzDecoder::new(std::io::BufReader::new(handle));
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new($crate::io::_impl::mzml_new::<_, $C, $D>(fh));
                            Ok($impl)
                        } else {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::mzml::MzMLReaderType<_, $C, $D> = $crate::io::_impl::mzml_new_indexed::<_, $C, $D>(handle);
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
                    $crate::io::MassSpectrometryFormat::MGF if $crate::io::_impl::mgf_support() => {
                        if compressed {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new(
                                $crate::io::_impl::mgf_new(
                                    $crate::io::RestartableGzDecoder::new(std::io::BufReader::new(buffered)),
                            ));
                            Ok($impl)
                        } else {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new(
                                $crate::io::_impl::mgf_new(buffered));
                            Ok($impl)
                        }
                    }
                    $crate::io::MassSpectrometryFormat::MzML if $crate::io::_impl::mzml_support() => {
                        if compressed {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new(
                                $crate::io::_impl::mzml_new($crate::io::RestartableGzDecoder::new(std::io::BufReader::new(buffered)),
                            ));
                            Ok($impl)
                        } else {
                            #[allow(unused_mut)]
                            let mut $reader: $crate::io::StreamingSpectrumIterator<$C, $D, _, _> = $crate::io::StreamingSpectrumIterator::new(
                                $crate::io::_impl::mzml_new(buffered));
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
    ($sink:expr_2021, $writer:ident => $impl:tt) => {
        mz_write!($sink, $writer => $impl, $crate::mzpeaks::CentroidPeak, $crate::mzpeaks::DeconvolutedPeak)
    };
    ($sink:expr_2021, $writer:ident => $impl:tt, $C:ty, $D:ty) => {{
        let sink = $crate::io::Sink::<$C, $D>::from($sink);
        match sink {
           $crate::io:: Sink::Sender(_) | $crate::io::Sink::SyncSender(_) =>  {
                Err(std::io::Error::new(std::io::ErrorKind::Unsupported, "Sender writers aren't supported by `mz_write`"))
            }
            $crate::io::Sink::PathLike(write_path) => {
                let (writer_format, is_gzip) = $crate::io::infer_from_path(&write_path);
                match writer_format {
                    $crate::io::MassSpectrometryFormat::MGF if $crate::io::_impl::mgf_support() => {
                        let handle = std::io::BufWriter::new(std::fs::File::create(&write_path)?);
                        if is_gzip {
                            let handle = flate2::write::GzEncoder::new(handle, flate2::Compression::best());
                            let mut $writer: $crate::io::_impl::MGFWriterType<_, $C, $D> = $crate::io::_impl::MGFWriterType::new(
                                handle,
                            );
                            Ok($impl)
                        } else {
                            let mut $writer: $crate::io::_impl::MGFWriterType<_, $C, $D> = $crate::io::_impl::MGFWriterType::new(
                                handle,
                            );
                            Ok($impl)

                        }
                    }
                    $crate::io::MassSpectrometryFormat::MzML if $crate::io::_impl::mzml_support() => {
                        let handle = std::io::BufWriter::new(std::fs::File::create(&write_path)?);
                        if is_gzip {
                            let handle = flate2::write::GzEncoder::new(handle, flate2::Compression::best());
                            let mut $writer: $crate::io::_impl::MzMLWriterType<_, $C, $D> = $crate::io::_impl::MzMLWriterType::new(
                                handle,
                            );
                            Ok($impl)
                        } else {
                            let mut $writer: $crate::io::_impl::MzMLWriterType<_, $C, $D> = $crate::io::_impl::MzMLWriterType::new(
                                handle,
                            );
                            Ok($impl)
                        }
                    }

                    $crate::io::MassSpectrometryFormat::MzMLb if $crate::io::_impl::mzmlb_support() => {
                        let mut $writer = $crate::io::_impl::MzMLbWriterType::<$C, $D>::new(&write_path)?;
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
                    $crate::io::MassSpectrometryFormat::MGF if $crate::io::_impl::mgf_support() => {
                        let handle = std::io::BufWriter::new(handle);
                        let mut $writer: $crate::io::_impl::MGFWriterType<_, $C, $D> = $crate::io::_impl::MGFWriterType::new(
                            handle,
                        );
                        Ok($impl)
                    }
                    $crate::io::MassSpectrometryFormat::MzML if $crate::io::_impl::mzml_support() => {
                        let handle = std::io::BufWriter::new(handle);
                        let mut $writer: $crate::io::_impl::MzMLWriterType<_, $C, $D> = $crate::io::_impl::MzMLWriterType::new(
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
