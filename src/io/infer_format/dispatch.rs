#![allow(clippy::type_complexity, clippy::large_enum_variant)]
use std::{fmt::Debug, fs, io, marker::PhantomData, path::{self, Path}};

use mzpeaks::{prelude::FeatureLike, CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak, KnownCharge};
use mzpeaks::{feature::{ChargedFeature, Feature}, IonMobility, Mass, MZ};

use crate::{io::{Generic3DIonMobilityFrameSource, IonMobilityFrameSource, IntoIonMobilityFrameSource, PreBufferedStream, RandomAccessIonMobilityFrameIterator}, spectrum::MultiLayerIonMobilityFrame};
#[cfg(feature = "mzmlb")]
pub use crate::io::mzmlb::MzMLbReaderType;

use crate::io::mgf::MGFReaderType;
use crate::io::mzml::MzMLReaderType;
use crate::io::traits::{RandomAccessSpectrumIterator, SpectrumSource, MZFileReader};
use crate::meta::MSDataFileMetadata;
use crate::spectrum::bindata::BuildFromArrayMap;
use crate::spectrum::MultiLayerSpectrum;

#[cfg(feature = "thermo")]
use crate::io::thermo::ThermoRawReaderType;

#[cfg(feature = "bruker_tdf")]
use crate::io::tdf::{TDFSpectrumReaderType, TDFFrameReaderType};

use crate::io::traits::{ChromatogramSource, StreamingSpectrumIterator};
use crate::io::{DetailLevel, SpectrumSourceWithMetadata};

use super::{infer_format, infer_from_stream, MassSpectrometryFormat};


/// An explicit file format dispatching ADT that provides the complete [`SpectrumSource`],
/// [`RandomAccessSpectrumIterator`], [`MZFileReader`] and [`MSDataFileMetadata`] APIs.
/// The preferred means of creating an instance is through the [`MZReaderType::open_path`]
/// function.
///
/// This type internally wraps a concrete type and each operation has to perform a very
/// fast test to decide which implementation to invoke. This overhead should be trivial
/// for the vast majority of operations, as compared to the I/O being performed otherwise.
///
/// Please refer to each concrete implementation type for details about those types. Keep
/// in mind that formats that require specific features to be enabled to be available won't
/// even be visible if their required features aren't set since their variants would have an
/// unknown size at compile time.
#[non_exhaustive]
pub enum MZReaderType<
        R: io::Read + io::Seek,
        C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap=CentroidPeak,
        D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap=DeconvolutedPeak> {
    MzML(MzMLReaderType<R, C, D>),
    MGF(MGFReaderType<R, C, D>),
    #[cfg(feature = "thermo")]
    ThermoRaw(ThermoRawReaderType<C, D>),
    #[cfg(feature = "mzmlb")]
    MzMLb(Box<MzMLbReaderType<C, D>>),
    #[cfg(feature = "bruker_tdf")]
    BrukerTDF(TDFSpectrumReaderType<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>, C, D>),
    Unknown(Box<dyn SpectrumSourceWithMetadata<C, D, MultiLayerSpectrum<C, D>> + Send>),
}

impl<
        R: io::Read + io::Seek,
        C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap,
        D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap> Debug for MZReaderType<R, C, D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        match self {
            Self::MzML(arg0) => f.debug_tuple("MzML").field(&arg0.source_file_name()).finish(),
            Self::MGF(arg0) => f.debug_tuple("MGF").field(&arg0.source_file_name()).finish(),
            #[cfg(feature = "thermo")]
            Self::ThermoRaw(arg0) => f.debug_tuple("ThermoRaw").field(&arg0.source_file_name()).finish(),
            #[cfg(feature = "mzmlb")]
            Self::MzMLb(arg0) => f.debug_tuple("MzMLb").field(&arg0.source_file_name()).finish(),
            #[cfg(feature = "bruker_tdf")]
            Self::BrukerTDF(arg0) => f.debug_tuple("BrukerTDF").field(&arg0.source_file_name()).finish(),
            Self::Unknown(arg0) => f.debug_tuple("Unknown").field(&arg0.source_file_name()).finish(),
        }
    }
}



/// A builder type for [`MZReaderType`].
///
/// To create an instance, see [`MZReaderType::builder`]
#[derive(Debug)]
pub struct MZReaderBuilder<
        C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap=CentroidPeak,
        D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap=DeconvolutedPeak> {
    buffer_size: Option<usize>,
    detail_level: DetailLevel,
    _c: PhantomData<C>,
    _d: PhantomData<D>,
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap, D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap> Default for MZReaderBuilder<C, D> {
    fn default() -> Self {
        Self { buffer_size: None, detail_level: Default::default(), _c: Default::default(), _d: Default::default() }
    }
}

#[allow(unused)]
impl<C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap, D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap> MZReaderBuilder<C, D> {

    /// Set the buffer capacity for a streaming reader.
    pub fn buffer_size(mut self, capacity: usize) -> Self {
        self.buffer_size = Some(capacity);
        self
    }

    /// Set the detail level for controlling how much work the reader
    /// will do to load peak information from spectra.
    pub fn detail_level(mut self, detail_level: DetailLevel) -> Self {
        self.detail_level = detail_level;
        self
    }

    /// Create a reader from a file on the local file system denoted by `path`.
    pub fn from_path<P: AsRef<Path>>(self, path: P) -> io::Result<MZReaderType<fs::File, C, D>> {
        let mut reader = MZReaderType::open_path(path.as_ref())?;
        reader.set_detail_level(self.detail_level);
        Ok(reader)
    }

    /// Create a reader from a type that supports [`io::Read`] and
    /// [`io::Seek`].
    ///
    /// # Note
    /// Not all formats can be read from an `io` type, these will
    /// fail to open and an error will be returned
    pub fn from_read_seek<R: io::Read + io::Seek>(self, source: R) -> io::Result<MZReaderType<R, C, D>> {
        let mut reader = MZReaderType::open_read_seek(source)?;
        reader.set_detail_level(self.detail_level);
        Ok(reader)
    }

    /// Create a reader from a type that supports [`io::Read`].
    ///
    /// This will internally wrap the file in a [`PreBufferedStream`] for metadata
    /// reading, but does not construct an index for full random access. Attempting
    /// to use the reader to access spectra may move the reader forwards, but it can
    /// never go backwards.
    ///
    /// # Note
    /// Not all formats can be read from an `io` type, these will
    /// fail to open and an error will be returned
    pub fn from_read<R: io::Read>(self, source: R) -> io::Result<StreamingSpectrumIterator<C, D, MultiLayerSpectrum<C, D>, MZReaderType<PreBufferedStream<R>, C, D>>> {
        let mut reader = if let Some(buffer_size) = self.buffer_size {
            MZReaderType::open_read_with_buffer_size(source, buffer_size)
        } else {
            MZReaderType::open_read(source)
        }?;
        reader.get_mut().set_detail_level(self.detail_level);
        Ok(reader)
    }
}



macro_rules! msfmt_dispatch {
    ($d:ident, $r:ident, $e:expr) => {
        match $d {
            MZReaderType::MzML($r) => $e,
            MZReaderType::MGF($r) => $e,
            #[cfg(feature = "thermo")]
            MZReaderType::ThermoRaw($r) => $e,
            #[cfg(feature = "mzmlb")]
            MZReaderType::MzMLb($r) => $e,
            #[cfg(feature = "bruker_tdf")]
            MZReaderType::BrukerTDF($r) => $e,
            MZReaderType::Unknown($r) => $e,
        }
    };
}


impl<R: io::Read + io::Seek,
     C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap> MZReaderType<R, C, D> {

    /// Create a [`MZReaderBuilder`] which can be used to configure the
    /// created reader, setting [`DetailLevel`] and buffer capacity
    /// The builder can create a reader from any type that the
    /// direct creation functions support.
    pub fn builder() -> MZReaderBuilder<C, D> {
        MZReaderBuilder::default()
    }

    /// Get the file format for this reader
    pub fn as_format(&self) -> MassSpectrometryFormat {
        match &self {
            MZReaderType::MzML(_) => MassSpectrometryFormat::MzML,
            MZReaderType::MGF(_) => MassSpectrometryFormat::MGF,
            #[cfg(feature = "thermo")]
            MZReaderType::ThermoRaw(_) => MassSpectrometryFormat::ThermoRaw,
            #[cfg(feature = "mzmlb")]
            MZReaderType::MzMLb(_) => MassSpectrometryFormat::MzMLb,
            #[cfg(feature = "bruker_tdf")]
            MZReaderType::BrukerTDF(_) => MassSpectrometryFormat::BrukerTDF,
            _ => MassSpectrometryFormat::Unknown
        }
    }

    /// Get the [`DetailLevel`] the reader currently uses
    pub fn detail_level(&self) -> &DetailLevel {
        msfmt_dispatch!(self, reader, {
            &reader.detail_level()
        })
    }

    /// Set the [`DetailLevel`] for the reader, changing
    /// the amount of work done immediately on loading a
    /// spectrum.
    ///
    /// # Note
    /// Not all readers support all detail levels, and the
    /// behavior when requesting one of those levels will
    /// depend upon the underlying reader.
    pub fn set_detail_level(&mut self, detail_level: DetailLevel) {
        msfmt_dispatch!(self, reader, {
            reader.set_detail_level(detail_level);
        });
    }

    /// Create a reader from a type that supports [`io::Read`] and
    /// [`io::Seek`].
    ///
    /// # Note
    /// Not all formats can be read from an `io` type, these will
    /// fail to open and an error will be returned
    pub fn open_read_seek(mut stream: R) -> io::Result<Self> {
        let (fmt, gzipped) = infer_from_stream(&mut stream)?;
        if gzipped {
            return Err(io::Error::new(io::ErrorKind::Unsupported, "This method does not support gzipped streams"))
        }
        match fmt {
            MassSpectrometryFormat::MGF => Ok(Self::MGF(MGFReaderType::new_indexed(stream))),
            MassSpectrometryFormat::MzML => Ok(Self::MzML(MzMLReaderType::new_indexed(stream))),
            _ => {
                Err(io::Error::new(io::ErrorKind::Unsupported, format!("This method does not support {fmt}")))
            }
        }
    }
}

impl<R: io::Read + io::Seek,
     C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap> ChromatogramSource for MZReaderType<R, C, D> {

    fn get_chromatogram_by_id(&mut self, id: &str) -> Option<crate::spectrum::Chromatogram> {
        match self {
            MZReaderType::MzML(r) => r.get_chromatogram_by_id(id),
            MZReaderType::MGF(r) => r.get_chromatogram_by_id(id),
            #[cfg(feature = "thermo")]
            MZReaderType::ThermoRaw(r) => r.get_chromatogram_by_id(id),
            #[cfg(feature = "mzmlb")]
            MZReaderType::MzMLb(r) => r.get_chromatogram_by_id(id),
            #[cfg(feature = "bruker_tdf")]
            MZReaderType::BrukerTDF(r) => r.get_chromatogram_by_id(id),
            _ => None
        }
    }

    fn get_chromatogram_by_index(&mut self, index: usize) -> Option<crate::spectrum::Chromatogram> {
        match self {
            MZReaderType::MzML(r) => r.get_chromatogram_by_index(index),
            MZReaderType::MGF(r) => r.get_chromatogram_by_index(index),
            #[cfg(feature = "thermo")]
            MZReaderType::ThermoRaw(r) => r.get_chromatogram_by_index(index),
            #[cfg(feature = "mzmlb")]
            MZReaderType::MzMLb(r) => r.get_chromatogram_by_index(index),
            #[cfg(feature = "bruker_tdf")]
            MZReaderType::BrukerTDF(r) => r.get_chromatogram_by_index(index),
            _ => None
        }
    }
}

impl<R: io::Read,
     C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap> MZReaderType<PreBufferedStream<R>, C, D> {

    /// Create a reader from a type that supports [`io::Read`].
    ///
    /// This will internally wrap the file in a [`PreBufferedStream`] for metadata
    /// reading, but does not construct an index for full random access. Attempting
    /// to use the reader to access spectra may move the reader forwards, but it can
    /// never go backwards.
    ///
    /// # Note
    /// Not all formats can be read from an `io` type, these will
    /// fail to open and an error will be returned
    pub fn open_read(stream: R) -> io::Result<StreamingSpectrumIterator<C, D, MultiLayerSpectrum<C, D>, Self>> {
        let mut stream = PreBufferedStream::new(stream)?;
        let (fmt, gzipped) = infer_from_stream(&mut stream)?;

        if gzipped {
            return Err(io::Error::new(io::ErrorKind::Unsupported, "This method does not support gzipped streams"))
        }

        let reader = match fmt {
            MassSpectrometryFormat::MGF => Self::MGF(MGFReaderType::new(stream)),
            MassSpectrometryFormat::MzML => Self::MzML(MzMLReaderType::new(stream)),
            _ => {
                return Err(io::Error::new(io::ErrorKind::Unsupported, format!("This method does not support {fmt}")))
            }
        };
        Ok(StreamingSpectrumIterator::new(reader))
    }

    /// See [`MZReaderType::open_read`].
    ///
    /// This function lets the caller specify the prebuffering size for files with large
    /// headers that exceed the default buffer size.
    pub fn open_read_with_buffer_size(stream: R, buffer_size: usize) -> io::Result<StreamingSpectrumIterator<C, D, MultiLayerSpectrum<C, D>, Self>> {
        let mut stream = PreBufferedStream::new_with_buffer_size(stream, buffer_size)?;
        let (fmt, gzipped) = infer_from_stream(&mut stream)?;

        if gzipped {
            return Err(io::Error::new(io::ErrorKind::Unsupported, "This method does not support gzipped streams"))
        }

        let reader = match fmt {
            MassSpectrometryFormat::MGF => Self::MGF(MGFReaderType::new(stream)),
            MassSpectrometryFormat::MzML => Self::MzML(MzMLReaderType::new(stream)),
            _ => {
                return Err(io::Error::new(io::ErrorKind::Unsupported, format!("This method does not support {fmt}")))
            }
        };
        Ok(StreamingSpectrumIterator::new(reader))
    }
}

/// A specialization of [`MZReaderType`] for the default peak types, for common use. The preferred means
/// of creating an instance is using the [`MZReader::open_path`] function.
pub type MZReader<R> = MZReaderType<R, CentroidPeak, DeconvolutedPeak>;

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap> MZFileReader<C, D, MultiLayerSpectrum<C, D>> for MZReaderType<fs::File, C, D> {

    fn construct_index_from_stream(&mut self) -> u64 {
        match self {
            MZReaderType::MzML(reader) => reader.construct_index_from_stream(),
            MZReaderType::MGF(reader) => reader.construct_index_from_stream(),
            #[cfg(feature = "thermo")]
            MZReaderType::ThermoRaw(reader) => reader.construct_index_from_stream(),
            #[cfg(feature = "mzmlb")]
            MZReaderType::MzMLb(reader) => reader.construct_index_from_stream(),
            #[cfg(feature = "bruker_tdf")]
            MZReaderType::BrukerTDF(reader) => reader.construct_index_from_stream(),
            MZReaderType::Unknown(reader) => reader.get_index().len() as u64,
        }
    }

    fn open_path<P>(path: P) -> io::Result<Self>
        where
            P: Into<path::PathBuf> + Clone, {
        let (format, is_gzipped) = infer_format(path.clone())?;
        if is_gzipped {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "Gzipped files are not supported",
            ))
        }
        match format {
            MassSpectrometryFormat::MGF => {
                let reader = MGFReaderType::open_path(path)?;
                Ok(Self::MGF(reader))
            }
            MassSpectrometryFormat::MzML => {
                let reader = MzMLReaderType::open_path(path)?;
                Ok(Self::MzML(reader))
            }
            #[cfg(feature = "thermo")]
            MassSpectrometryFormat::ThermoRaw => {
                let reader = ThermoRawReaderType::open_path(path)?;
                Ok(Self::ThermoRaw(reader))
            }
            #[cfg(feature = "mzmlb")]
            MassSpectrometryFormat::MzMLb => {
                let reader = MzMLbReaderType::open_path(path)?;
                Ok(Self::MzMLb(Box::new(reader)))
            }
            _ => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "File format not supported",
            )),
        }

    }

    fn open_file(mut source: fs::File) -> io::Result<Self> {
        let (format, is_gzipped) = infer_from_stream(&mut source)?;

        if is_gzipped {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "Gzipped files are not supported",
            ))
        }
        match format {
            MassSpectrometryFormat::MGF => {
                let reader = MGFReaderType::open_file(source)?;
                Ok(Self::MGF(reader))
            }
            MassSpectrometryFormat::MzML => {
                let reader = MzMLReaderType::open_file(source)?;
                Ok(Self::MzML(reader))
            }
            #[cfg(feature = "thermo")]
            MassSpectrometryFormat::ThermoRaw => {
                let reader = ThermoRawReaderType::open_file(source)?;
                Ok(Self::ThermoRaw(reader))
            }
            #[cfg(feature = "mzmlb")]
            MassSpectrometryFormat::MzMLb => {
                let reader = MzMLbReaderType::open_file(source)?;
                Ok(Self::MzMLb(Box::new(reader)))
            }
            _ => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "File format not supported",
            )),
        }
    }
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap,
     R: io::Read + io::Seek> Iterator for MZReaderType<R, C, D> {

    type Item = MultiLayerSpectrum<C, D>;

    fn next(&mut self) -> Option<Self::Item> {
        msfmt_dispatch!(self, reader, reader.next())
    }
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap,
     R: io::Read + io::Seek> SpectrumSource<C, D, MultiLayerSpectrum<C, D>> for MZReaderType<R, C, D> {

    fn reset(&mut self) {
        msfmt_dispatch!(self, reader, reader.reset())
    }

    fn get_spectrum_by_id(&mut self, id: &str) -> Option<MultiLayerSpectrum<C, D>> {
        msfmt_dispatch!(self, reader, reader.get_spectrum_by_id(id))
    }

    fn get_spectrum_by_index(&mut self, index: usize) -> Option<MultiLayerSpectrum<C, D>> {
        msfmt_dispatch!(self, reader, reader.get_spectrum_by_index(index))
    }

    fn get_spectrum_by_time(&mut self, time: f64) -> Option<MultiLayerSpectrum<C, D>> {
        match self {
            MZReaderType::MzML(reader) => reader.get_spectrum_by_time(time),
            MZReaderType::MGF(reader) => reader.get_spectrum_by_time(time),
            #[cfg(feature = "thermo")]
            MZReaderType::ThermoRaw(reader) => reader.get_spectrum_by_time(time),
            #[cfg(feature = "mzmlb")]
            MZReaderType::MzMLb(reader) => reader.get_spectrum_by_time(time),
            #[cfg(feature = "bruker_tdf")]
            MZReaderType::BrukerTDF(r) => r.get_spectrum_by_time(time),
            MZReaderType::Unknown(r) => r.get_spectrum_by_time(time),
        }
    }

    fn get_index(&self) -> &crate::io::OffsetIndex {
        msfmt_dispatch!(self, reader, reader.get_index())
    }

    fn set_index(&mut self, index: crate::io::OffsetIndex) {
        msfmt_dispatch!(self, reader, reader.set_index(index))
    }

    fn detail_level(&self) -> &DetailLevel {
        self.detail_level()
    }

    fn set_detail_level(&mut self, detail_level: DetailLevel) {
        self.set_detail_level(detail_level);
    }
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap,
     R: io::Read + io::Seek> MSDataFileMetadata for MZReaderType<R, C, D> {

    fn data_processings(&self) -> &Vec<crate::meta::DataProcessing> {
        msfmt_dispatch!(self, reader, reader.data_processings())
    }

    fn instrument_configurations(&self) -> &std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
        msfmt_dispatch!(self, reader, reader.instrument_configurations())
    }

    fn file_description(&self) -> &crate::meta::FileDescription {
        msfmt_dispatch!(self, reader, reader.file_description())
    }

    fn softwares(&self) -> &Vec<crate::meta::Software> {
        msfmt_dispatch!(self, reader, reader.softwares())
    }

    fn samples(&self) -> &Vec<crate::meta::Sample> {
        msfmt_dispatch!(self, reader, reader.samples())
    }

    fn data_processings_mut(&mut self) -> &mut Vec<crate::meta::DataProcessing> {
        msfmt_dispatch!(self, reader, reader.data_processings_mut())
    }

    fn instrument_configurations_mut(&mut self) -> &mut std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
        msfmt_dispatch!(self, reader, reader.instrument_configurations_mut())
    }

    fn file_description_mut(&mut self) -> &mut crate::meta::FileDescription {
        msfmt_dispatch!(self, reader, reader.file_description_mut())
    }

    fn softwares_mut(&mut self) -> &mut Vec<crate::meta::Software> {
        msfmt_dispatch!(self, reader, reader.softwares_mut())
    }

    fn samples_mut(&mut self) -> &mut Vec<crate::meta::Sample> {
        msfmt_dispatch!(self, reader, reader.samples_mut())
    }

    fn source_file_name(&self) -> Option<&str> {
        msfmt_dispatch!(self, reader, reader.source_file_name())
    }

    fn run_description(&self) -> Option<&crate::meta::MassSpectrometryRun> {
        msfmt_dispatch!(self, reader, reader.run_description())
    }

    fn spectrum_count_hint(&self) -> Option<u64> {
        msfmt_dispatch!(self, reader, reader.spectrum_count_hint())
    }
}

macro_rules! msfmt_dispatch_cap {
    ($d:ident, $r:ident, $e:expr) => {
        match $d {
            MZReaderType::MzML($r) => {
                $e?;
            },
            MZReaderType::MGF($r) => {
                $e?;
            },
            #[cfg(feature = "thermo")]
            MZReaderType::ThermoRaw($r) => {
                $e?;
            },
            #[cfg(feature = "mzmlb")]
            MZReaderType::MzMLb($r) => {
                $e?;
            },
            #[cfg(feature = "bruker_tdf")]
            MZReaderType::BrukerTDF($r) => {
                $e?;
            }
            MZReaderType::Unknown(_) => {
                Err(crate::io::SpectrumAccessError::IOError(Some(io::Error::new(io::ErrorKind::Unsupported, "Dynamic adaptor doesn't know how to do random access iterators"))))?
            }
        };
    };
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap,
     R: io::Read + io::Seek> RandomAccessSpectrumIterator<C, D, MultiLayerSpectrum<C, D>> for MZReaderType<R, C, D> {
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, crate::io::SpectrumAccessError> {
        msfmt_dispatch_cap!(self, reader, reader.start_from_id(id));
        Ok(self)
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, crate::io::SpectrumAccessError> {
        msfmt_dispatch_cap!(self, reader, reader.start_from_index(index));
        Ok(self)
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, crate::io::SpectrumAccessError> {
        msfmt_dispatch_cap!(self, reader, reader.start_from_time(time));
        Ok(self)
    }
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap,
     R: io::Read + io::Seek> IntoIonMobilityFrameSource<C, D> for MZReaderType<R, C, D> {

    type IonMobilityFrameSource<CF: FeatureLike<MZ, IonMobility>, DF: FeatureLike<Mass, IonMobility> + KnownCharge> = IMMZReaderType<R, CF, DF, C, D>;

    fn try_into_frame_source<CF: FeatureLike<MZ, IonMobility>, DF: FeatureLike<Mass, IonMobility> + KnownCharge>(self) -> Result<Self::IonMobilityFrameSource<CF, DF>, crate::io::IntoIonMobilityFrameSourceError> {
        let view = match self {
            MZReaderType::MzML(reader) => IMMZReaderType::MzML(reader.try_into_frame_source()?),
            MZReaderType::MGF(_) => return Err(crate::io::IntoIonMobilityFrameSourceError::NoIonMobilityFramesFound),
            #[cfg(feature="thermo")]
            MZReaderType::ThermoRaw(_) => return Err(crate::io::IntoIonMobilityFrameSourceError::NoIonMobilityFramesFound),
            #[cfg(feature="mzmlb")]
            MZReaderType::MzMLb(reader) => IMMZReaderType::MzMLb(Box::new(reader.try_into_frame_source()?)),
            #[cfg(feature="bruker_tdf")]
            MZReaderType::BrukerTDF(reader) => IMMZReaderType::BrukerTDF(reader.try_into_frame_source()?),
            MZReaderType::Unknown(_) => todo!(),
        };
        Ok(view)
    }
}

#[cfg(feature = "async_partial")]
mod async_impl {
    use super::*;
    use io::SeekFrom;
    use tokio::io::{AsyncRead, AsyncSeek, AsyncSeekExt, AsyncReadExt};
    use crate::io::{
        mgf::AsyncMGFReaderType,
        mzml::AsyncMzMLReaderType,
        traits::{AsyncMZFileReader, AsyncRandomAccessSpectrumIterator, AsyncSpectrumSource}};

    #[cfg(feature = "thermo")]
    use crate::io::thermo::AsyncThermoRawReaderType;

    #[non_exhaustive]
    pub enum AsyncMZReaderType<
        R: AsyncRead + AsyncSeek + Unpin + Send,
        C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap + Send + Sync=CentroidPeak,
        D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap + Send + Sync=DeconvolutedPeak,
    > {
        MzML(AsyncMzMLReaderType<R, C, D>),
        MGF(AsyncMGFReaderType<R, C, D>),
        #[cfg(feature = "thermo")]
        ThermoRaw(AsyncThermoRawReaderType<C, D>),
    }

    impl<R: AsyncRead + AsyncSeek + Unpin + Send,
         C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap + Send + Sync,
         D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap + Send + Sync> AsyncMZReaderType<R, C, D> {

        pub fn builder() -> AsyncMZReaderBuilder<C, D> {
            AsyncMZReaderBuilder::default()
        }

        pub async fn open_read_seek(mut source: R) -> io::Result<AsyncMZReaderType<R, C, D>> {
            let mut buffer: Vec<u8> = vec![0; 250];
            let current = source.stream_position().await?;
            source.read_exact(&mut buffer).await?;
            source.seek(SeekFrom::Start(current)).await?;
            let mut stream = io::Cursor::new(buffer);
            let (ms_format, gzipped) = infer_from_stream(&mut stream)?;

            if gzipped {
                return Err(io::Error::new(io::ErrorKind::Unsupported, "Cannot read compressed streams with this interface"))
            }

            match ms_format {
                MassSpectrometryFormat::MGF => {
                    Ok(Self::MGF(AsyncMGFReaderType::new_indexed(source).await))
                },
                MassSpectrometryFormat::MzML => {
                    Ok(Self::MzML(AsyncMzMLReaderType::new_indexed(source).await))
                },
                _ => Err(io::Error::new(io::ErrorKind::Unsupported, format!("Cannot read {ms_format} files from streams")))
            }
        }

    }

    pub type AsyncMZReader<R> = AsyncMZReaderType<R, CentroidPeak, DeconvolutedPeak>;

    macro_rules! amsfmt_dispatch {
        ($d:ident, $r:ident, $e:expr) => {
            match $d {
                AsyncMZReaderType::MzML($r) => $e,
                AsyncMZReaderType::MGF($r) => $e,
                #[cfg(feature = "thermo")]
                AsyncMZReaderType::ThermoRaw($r) => $e,
            }
        };
    }

    impl<R: AsyncRead + AsyncSeek + Unpin + Send,
         C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap + Send + Sync + 'static,
         D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap + Send + Sync + 'static> MSDataFileMetadata for AsyncMZReaderType<R, C, D> {

        fn data_processings(&self) -> &Vec<crate::meta::DataProcessing> {
            amsfmt_dispatch!(self, reader, reader.data_processings())
        }

        fn instrument_configurations(&self) -> &std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
            amsfmt_dispatch!(self, reader, reader.instrument_configurations())
        }

        fn file_description(&self) -> &crate::meta::FileDescription {
            amsfmt_dispatch!(self, reader, reader.file_description())
        }

        fn softwares(&self) -> &Vec<crate::meta::Software> {
            amsfmt_dispatch!(self, reader, reader.softwares())
        }

        fn samples(&self) -> &Vec<crate::meta::Sample> {
            amsfmt_dispatch!(self, reader, reader.samples())
        }

        fn data_processings_mut(&mut self) -> &mut Vec<crate::meta::DataProcessing> {
            amsfmt_dispatch!(self, reader, reader.data_processings_mut())
        }

        fn instrument_configurations_mut(&mut self) -> &mut std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
            amsfmt_dispatch!(self, reader, reader.instrument_configurations_mut())
        }

        fn file_description_mut(&mut self) -> &mut crate::meta::FileDescription {
            amsfmt_dispatch!(self, reader, reader.file_description_mut())
        }

        fn softwares_mut(&mut self) -> &mut Vec<crate::meta::Software> {
            amsfmt_dispatch!(self, reader, reader.softwares_mut())
        }

        fn samples_mut(&mut self) -> &mut Vec<crate::meta::Sample> {
            amsfmt_dispatch!(self, reader, reader.samples_mut())
        }

        fn source_file_name(&self) -> Option<&str> {
            amsfmt_dispatch!(self, reader, reader.source_file_name())
        }

        fn run_description(&self) -> Option<&crate::meta::MassSpectrometryRun> {
            amsfmt_dispatch!(self, reader, reader.run_description())
        }

        fn spectrum_count_hint(&self) -> Option<u64> {
            amsfmt_dispatch!(self, reader, reader.spectrum_count_hint())
        }
    }


    impl<R: AsyncRead + AsyncSeek + Unpin + Send,
         C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap + Send + Sync + 'static,
         D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap + Send + Sync + 'static> AsyncSpectrumSource<C, D, MultiLayerSpectrum<C, D>> for AsyncMZReaderType<R, C, D> {
        async fn reset(&mut self) {
            amsfmt_dispatch!(self, reader, reader.reset().await)
        }

        fn detail_level(&self) -> &DetailLevel {
            amsfmt_dispatch!(self, reader, reader.detail_level())
        }

        fn set_detail_level(&mut self, detail_level: DetailLevel) {
            amsfmt_dispatch!(self, reader, reader.set_detail_level(detail_level))
        }

        async fn get_spectrum_by_id(&mut self, id: &str) -> Option<MultiLayerSpectrum<C, D>> {
            amsfmt_dispatch!(self, reader, reader.get_spectrum_by_id(id).await)
        }

        async fn get_spectrum_by_index(&mut self, index: usize) -> Option<MultiLayerSpectrum<C, D>> {
            amsfmt_dispatch!(self, reader, reader.get_spectrum_by_index(index).await)
        }

        async fn get_spectrum_by_time(&mut self, time: f64) -> Option<MultiLayerSpectrum<C, D>> {
            amsfmt_dispatch!(self, reader, reader.get_spectrum_by_time(time).await)
        }

        fn get_index(&self) -> &crate::io::OffsetIndex {
            amsfmt_dispatch!(self, reader, reader.get_index())
        }

        fn set_index(&mut self, index: crate::io::OffsetIndex) {
            amsfmt_dispatch!(self, reader, reader.set_index(index))
        }

        async fn read_next(&mut self) -> Option<MultiLayerSpectrum<C, D>> {
            amsfmt_dispatch!(self, reader, reader.read_next().await)
        }
    }

    impl<R: AsyncRead + AsyncSeek + Unpin + Send,
         C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap + Send + Sync + 'static,
         D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap + Send + Sync + 'static> AsyncRandomAccessSpectrumIterator<C, D, MultiLayerSpectrum<C, D>> for AsyncMZReaderType<R, C, D> {

        async fn start_from_id(&mut self, id: &str) -> Result<&mut Self, crate::prelude::SpectrumAccessError> {
            amsfmt_dispatch!(self, reader,{ reader.start_from_id(id).await?;});
            Ok(self)
        }

        async fn start_from_index(&mut self, index: usize) -> Result<&mut Self, crate::prelude::SpectrumAccessError> {
            amsfmt_dispatch!(self, reader,{ reader.start_from_index(index).await?;});
            Ok(self)
        }

        async fn start_from_time(&mut self, time: f64) -> Result<&mut Self, crate::prelude::SpectrumAccessError> {
            amsfmt_dispatch!(self, reader,{ reader.start_from_time(time).await?;});
            Ok(self)
        }
    }

    #[cfg(feature = "async")]
    impl<C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap + Send + Sync + 'static,
         D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap + Send + Sync + 'static> AsyncMZFileReader<C, D, MultiLayerSpectrum<C,D>> for AsyncMZReaderType<tokio::fs::File, C, D> {

        async fn construct_index_from_stream(&mut self) -> u64 {
            amsfmt_dispatch!(self, reader, reader.construct_index_from_stream().await)
        }

        async fn open_path<P>(path: P) -> io::Result<Self>
                where
                    P: Into<path::PathBuf>, {
            let path: path::PathBuf = path.into();
            let (ms_format, gzipped) = infer_format(path.clone())?;
            if gzipped {
                return Err(io::Error::new(io::ErrorKind::Unsupported, "Compressed files are not supported"))
            }
            match ms_format {
                MassSpectrometryFormat::MGF => {
                    Ok(Self::MGF(AsyncMGFReaderType::open_path(path).await?))
                },
                MassSpectrometryFormat::MzML => {
                    Ok(Self::MzML(AsyncMzMLReaderType::open_path(path).await?))
                },
                MassSpectrometryFormat::MzMLb => {
                    Err(io::Error::new(io::ErrorKind::Unsupported, "MzMLb files are not supported in async mode"))
                },
                MassSpectrometryFormat::ThermoRaw => {
                    #[cfg(feature = "thermo")]
                    {Ok(Self::ThermoRaw(AsyncThermoRawReaderType::open_path(path).await?))}
                    #[cfg(not(feature = "thermo"))]
                    {Err(io::Error::new(io::ErrorKind::Unsupported, "Thermo RAW files are not supported. Enable the 'thermo' feature"))}
                },
                MassSpectrometryFormat::BrukerTDF => Err(io::Error::new(io::ErrorKind::Unsupported, "Bruker TDF files are not supported in async mode")),
                MassSpectrometryFormat::Unknown => Err(io::Error::new(io::ErrorKind::Unsupported, "Unknown file format unsupported")),
            }
        }

        async fn open_file(mut source: tokio::fs::File) -> io::Result<Self> {
            let mut sync_source = source.into_std().await;
            let (ms_format, compressed) = infer_from_stream(&mut sync_source)?;
            if compressed {
                return Err(io::Error::new(io::ErrorKind::Unsupported, "Compressed files are not supported"))
            }
            source = tokio::fs::File::from_std(sync_source);
            match ms_format {
                MassSpectrometryFormat::MGF => {
                    Ok(Self::MGF(AsyncMGFReaderType::open_file(source).await?))
                },
                MassSpectrometryFormat::MzML => {
                    Ok(Self::MzML(AsyncMzMLReaderType::open_file(source).await?))
                },
                _ => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "File format not supported",
                )),
            }
        }
    }


    /// A builder type for [`AsyncMZReaderType`].
    ///
    /// To create an instance, see [`AsyncMZReaderType::builder`]
    #[derive(Debug)]
    pub struct AsyncMZReaderBuilder<
            C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap + Send + Sync + 'static,
            D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap + Send + Sync + 'static> {
        buffer_size: Option<usize>,
        detail_level: DetailLevel,
        _c: PhantomData<C>,
        _d: PhantomData<D>,
    }

    impl<C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap + Send + Sync + 'static,
         D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap + Send + Sync + 'static> Default for AsyncMZReaderBuilder<C, D> {
        fn default() -> Self {
            Self { buffer_size: None, detail_level: Default::default(), _c: Default::default(), _d: Default::default() }
        }
    }

    #[allow(unused)]
    impl<C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap + Send + Sync + 'static,
         D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap + Send + Sync + 'static> AsyncMZReaderBuilder<C, D> {

        /// Set the buffer capacity for a streaming reader.
        pub fn buffer_size(mut self, capacity: usize) -> Self {
            self.buffer_size = Some(capacity);
            self
        }

        /// Set the detail level for controlling how much work the reader
        /// will do to load peak information from spectra.
        pub fn detail_level(mut self, detail_level: DetailLevel) -> Self {
            self.detail_level = detail_level;
            self
        }

        #[cfg(feature = "async")]
        /// Create a reader from a file on the local file system denoted by `path`.
        pub async fn from_path<P: AsRef<Path>>(self, path: P) -> io::Result<AsyncMZReaderType<tokio::fs::File, C, D>> {
            let mut reader = AsyncMZReaderType::open_path(path::PathBuf::from(path.as_ref())).await?;
            reader.set_detail_level(self.detail_level);
            Ok(reader)
        }

        /// Create a reader from a type that supports [`tokio::io::Read`] and
        /// [`tokio::io::Seek`].
        ///
        /// # Note
        /// Not all formats can be read from an `io` type, these will
        /// fail to open and an error will be returned
        pub async fn from_read_seek<R: tokio::io::AsyncRead + tokio::io::AsyncSeek + Unpin + Send>(self, source: R) -> io::Result<AsyncMZReaderType<R, C, D>> {
            let mut reader = AsyncMZReaderType::open_read_seek(source).await?;
            reader.set_detail_level(self.detail_level);
            Ok(reader)
        }

        // /// Create a reader from a type that supports [`io::Read`].
        // ///
        // /// This will internally wrap the file in a [`PreBufferedStream`] for metadata
        // /// reading, but does not construct an index for full random access. Attempting
        // /// to use the reader to access spectra may move the reader forwards, but it can
        // /// never go backwards.
        // ///
        // /// # Note
        // /// Not all formats can be read from an `io` type, these will
        // /// fail to open and an error will be returned
        // pub fn from_read<R: io::Read>(self, source: R) -> io::Result<StreamingSpectrumIterator<C, D, MultiLayerSpectrum<C, D>, MZReaderType<PreBufferedStream<R>, C, D>>> {
        //     let mut reader = if let Some(buffer_size) = self.buffer_size {
        //         MZReaderType::open_read_with_buffer_size(source, buffer_size)
        //     } else {
        //         MZReaderType::open_read(source)
        //     }?;
        //     reader.get_mut().set_detail_level(self.detail_level);
        //     Ok(reader)
        // }
    }


}

#[cfg(feature = "async_partial")]
pub use async_impl::{AsyncMZReaderType, AsyncMZReader, AsyncMZReaderBuilder};

pub enum IMMZReaderType<
        R: io::Read + io::Seek,
        C: FeatureLike<MZ, IonMobility> = Feature<MZ, IonMobility>,
        D: FeatureLike<Mass, IonMobility> + KnownCharge = ChargedFeature<Mass, IonMobility>,
        CP: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap=CentroidPeak,
        DP: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap=DeconvolutedPeak,
    > {
    MzML(Generic3DIonMobilityFrameSource<CP, DP, MzMLReaderType<R, CP, DP>, C, D>),
    #[cfg(feature = "mzmlb")]
    MzMLb(Box<Generic3DIonMobilityFrameSource<CP, DP, MzMLbReaderType<CP, DP>, C, D>>),
    #[cfg(feature = "bruker_tdf")]
    BrukerTDF(TDFFrameReaderType<C, D>),
}


macro_rules! immsfmt_dispatch {
    ($d:ident, $r:ident, $e:expr) => {
        match $d {
            IMMZReaderType::MzML($r) => $e,
            #[cfg(feature = "mzmlb")]
            IMMZReaderType::MzMLb($r) => $e,
            #[cfg(feature = "bruker_tdf")]
            IMMZReaderType::BrukerTDF($r) => $e,
        }
    };
}


impl<R: io::Read + io::Seek, C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge, CP: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap, DP: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap> MSDataFileMetadata for IMMZReaderType<R, C, D, CP, DP> {
    fn data_processings(&self) -> &Vec<crate::meta::DataProcessing> {
        immsfmt_dispatch!(self, reader, reader.data_processings())
    }

    fn instrument_configurations(&self) -> &std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
        immsfmt_dispatch!(self, reader, reader.instrument_configurations())
    }

    fn file_description(&self) -> &crate::meta::FileDescription {
        immsfmt_dispatch!(self, reader, reader.file_description())
    }

    fn softwares(&self) -> &Vec<crate::meta::Software> {
        immsfmt_dispatch!(self, reader, reader.softwares())
    }

    fn samples(&self) -> &Vec<crate::meta::Sample> {
        immsfmt_dispatch!(self, reader, reader.samples())
    }

    fn data_processings_mut(&mut self) -> &mut Vec<crate::meta::DataProcessing> {
        immsfmt_dispatch!(self, reader, reader.data_processings_mut())
    }

    fn instrument_configurations_mut(&mut self) -> &mut std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
        immsfmt_dispatch!(self, reader, reader.instrument_configurations_mut())
    }

    fn file_description_mut(&mut self) -> &mut crate::meta::FileDescription {
        immsfmt_dispatch!(self, reader, reader.file_description_mut())
    }

    fn softwares_mut(&mut self) -> &mut Vec<crate::meta::Software> {
        immsfmt_dispatch!(self, reader, reader.softwares_mut())
    }

    fn samples_mut(&mut self) -> &mut Vec<crate::meta::Sample> {
        immsfmt_dispatch!(self, reader, reader.samples_mut())
    }

    fn source_file_name(&self) -> Option<&str> {
        immsfmt_dispatch!(self, reader, reader.source_file_name())
    }

    fn run_description(&self) -> Option<&crate::meta::MassSpectrometryRun> {
        immsfmt_dispatch!(self, reader, reader.run_description())
    }

    fn spectrum_count_hint(&self) -> Option<u64> {
        immsfmt_dispatch!(self, reader, reader.spectrum_count_hint())
    }
}

impl<R: io::Read + io::Seek, C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge, CP: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap, DP: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap> Iterator for IMMZReaderType<R, C, D, CP, DP> {
    type Item = MultiLayerIonMobilityFrame<C, D>;

    fn next(&mut self) -> Option<Self::Item> {
        immsfmt_dispatch!(self, reader, reader.next())
    }
}

impl<R: io::Read + io::Seek, C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge, CP: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap, DP: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap> IonMobilityFrameSource<C, D, MultiLayerIonMobilityFrame<C, D>> for IMMZReaderType<R, C, D, CP, DP> {
    fn reset(&mut self) {
        immsfmt_dispatch!(self, reader, reader.reset())
    }

    fn detail_level(&self) -> &DetailLevel {
        immsfmt_dispatch!(self, reader, reader.detail_level())
    }

    fn set_detail_level(&mut self, detail_level: DetailLevel) {
        immsfmt_dispatch!(self, reader, reader.set_detail_level(detail_level))
    }

    fn get_frame_by_id(&mut self, id: &str) -> Option<MultiLayerIonMobilityFrame<C, D>> {
        immsfmt_dispatch!(self, reader, reader.get_frame_by_id(id))
    }

    fn get_frame_by_index(&mut self, index: usize) -> Option<MultiLayerIonMobilityFrame<C, D>> {
        immsfmt_dispatch!(self, reader, reader.get_frame_by_index(index))
    }

    fn get_index(&self) -> &crate::io::OffsetIndex {
        immsfmt_dispatch!(self, reader, reader.get_index())
    }

    fn set_index(&mut self, index: crate::io::OffsetIndex) {
        immsfmt_dispatch!(self, reader, reader.set_index(index))
    }
}

impl<R: io::Read + io::Seek, C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge, CP: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap, DP: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap> RandomAccessIonMobilityFrameIterator<C, D, MultiLayerIonMobilityFrame<C, D>> for IMMZReaderType<R, C, D, CP, DP> {
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, crate::io::IonMobilityFrameAccessError> {
        immsfmt_dispatch!(self, reader, {reader.start_from_id(id)?;});
        Ok(self)
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, crate::io::IonMobilityFrameAccessError> {
        immsfmt_dispatch!(self, reader, {reader.start_from_index(index)?;});
        Ok(self)
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, crate::io::IonMobilityFrameAccessError> {
        immsfmt_dispatch!(self, reader, {reader.start_from_time(time)?;});
        Ok(self)
    }
}

