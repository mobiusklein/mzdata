use std::convert::TryFrom;
use std::fmt::Display;
use std::fs;
use std::io::{self, prelude::*, BufReader};
use std::marker::PhantomData;
use std::path::{self, Path, PathBuf};
use std::sync::mpsc::{Receiver, Sender, SyncSender};


use flate2::{bufread::GzDecoder, write::GzEncoder};
use mzpeaks::{CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak};
#[cfg(feature = "bruker_tdf")]
use mzpeaks::{feature::{ChargedFeature, Feature}, IonMobility, Mass, MZ};

use crate::io::PreBufferedStream;
use crate::params::ControlledVocabulary;
#[cfg(feature = "mzmlb")]
pub use crate::io::mzmlb::{MzMLbReaderType, MzMLbWriterBuilder};

use crate::io::compression::{is_gzipped, is_gzipped_extension, RestartableGzDecoder};
use crate::io::mgf::{is_mgf, MGFReaderType, MGFWriterType};
use crate::io::mzml::{is_mzml, MzMLReaderType, MzMLWriterType};
use crate::io::traits::{RandomAccessSpectrumIterator, SpectrumSource, SpectrumWriter, MZFileReader};
use crate::meta::{FormatConversion, MSDataFileMetadata};
use crate::spectrum::bindata::{BuildArrayMapFrom, BuildFromArrayMap};
use crate::spectrum::MultiLayerSpectrum;
use crate::Param;

#[cfg(feature = "thermo")]
use super::thermo::{ThermoRawReaderType, is_thermo_raw_prefix};

#[cfg(feature = "bruker_tdf")]
use super::tdf::{is_tdf, TDFSpectrumReaderType};

use super::traits::{ChromatogramSource, SeekRead, SpectrumReceiver, StreamingSpectrumIterator};
use super::DetailLevel;

/// Mass spectrometry file formats that [`mzdata`](crate)
/// supports
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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
            _ => None
        }
    }

    pub fn as_param(&self) -> Option<Param> {
        let p = match self {
            MassSpectrometryFormat::MGF => ControlledVocabulary::MS.const_param_ident("Mascot MGF format", 1001062),
            MassSpectrometryFormat::MzML => ControlledVocabulary::MS.const_param_ident("MzML format", 1000584),
            MassSpectrometryFormat::MzMLb => ControlledVocabulary::MS.const_param_ident("mzMLb format", 1002838),
            MassSpectrometryFormat::ThermoRaw => ControlledVocabulary::MS.const_param_ident("Thermo RAW format", 1000563),
            MassSpectrometryFormat::BrukerTDF => ControlledVocabulary::MS.const_param_ident("Bruker TDF format", 1002817),
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
    MzMLb(MzMLbReaderType<C, D>),
    #[cfg(feature = "bruker_tdf")]
    BrukerTDF(TDFSpectrumReaderType<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>, C, D>)
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
            MZReaderType::BrukerTDF(_) => MassSpectrometryFormat::BrukerTDF
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
            MZReaderType::BrukerTDF(r) => r.get_chromatogram_by_index(index)
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
        msfmt_dispatch!(self, reader, reader.construct_index_from_stream())
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
                Ok(Self::MzMLb(reader))
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
                Ok(Self::MzMLb(reader))
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
            MZReaderType::BrukerTDF(r) => r.get_spectrum_by_time(time)
        }
    }

    fn get_index(&self) -> &super::OffsetIndex {
        msfmt_dispatch!(self, reader, reader.get_index())
    }

    fn set_index(&mut self, index: super::OffsetIndex) {
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
        };
    };
}

impl<C: CentroidLike + Default + From<CentroidPeak> + BuildFromArrayMap,
     D: DeconvolutedCentroidLike + Default + From<DeconvolutedPeak> + BuildFromArrayMap,
     R: io::Read + io::Seek> RandomAccessSpectrumIterator<C, D, MultiLayerSpectrum<C, D>> for MZReaderType<R, C, D> {
    fn start_from_id(&mut self, id: &str) -> Result<&mut Self, super::SpectrumAccessError> {
        msfmt_dispatch_cap!(self, reader, reader.start_from_id(id));
        Ok(self)
    }

    fn start_from_index(&mut self, index: usize) -> Result<&mut Self, super::SpectrumAccessError> {
        msfmt_dispatch_cap!(self, reader, reader.start_from_index(index));
        Ok(self)
    }

    fn start_from_time(&mut self, time: f64) -> Result<&mut Self, super::SpectrumAccessError> {
        msfmt_dispatch_cap!(self, reader, reader.start_from_time(time));
        Ok(self)
    }
}

/// Given a path, infer the file format and whether or not the file at that path is
/// GZIP compressed
pub fn infer_from_path<P: Into<path::PathBuf>>(path: P) -> (MassSpectrometryFormat, bool) {
    let path: path::PathBuf = path.into();
    if path.is_dir() {
        #[cfg(feature = "bruker_tdf")]
        if is_tdf(path) {
            return (MassSpectrometryFormat::BrukerTDF, false)
        } else {
            return (MassSpectrometryFormat::Unknown, false)
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
        _ if is_mzml(&buf) => Ok((MassSpectrometryFormat::MzML, is_stream_gzipped)),
        _ if is_mgf(&buf) => Ok((MassSpectrometryFormat::MGF, is_stream_gzipped)),
        #[cfg(feature = "thermo")]
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

    #[cfg(feature = "thermo")]
    #[test]
    fn infer_thermo() {
        let path = path::Path::new("./test/data/small.RAW");
        let (fmt, zipped) = infer_from_path(path);
        assert_eq!(fmt, MassSpectrometryFormat::ThermoRaw);
        assert!(!zipped);
    }

    #[test]
    fn infer_open() {
        let path = path::Path::new("./test/data/small.mzML");
        assert!(path.exists());
        if let Ok(mut reader) = MZReader::open_path(path) {
            assert_eq!(reader.len(), 48);
            assert_eq!(*reader.detail_level(), DetailLevel::Full);
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
            } else {
                panic!("Failed to retrieve spectrum by index")
            }

            assert_eq!(reader.get_spectrum_by_id("controllerType=0 controllerNumber=1 scan=11").unwrap().index(), 10);

            if let Some(spec) = reader.get_spectrum_by_time(0.358558333333) {
                assert_eq!(spec.index(), 34);
            } else {
                panic!("Failed to retrieve spectrum by time")
            }

        } else {
            panic!("Failed to open file")
        }
    }

    #[cfg(feature = "thermo")]
    #[test]
    fn infer_open_thermo() {
        let path = path::Path::new("./test/data/small.RAW");
        assert!(path.exists());
        if let Ok(mut reader) = MZReader::open_path(path) {
            assert_eq!(reader.len(), 48);
            assert_eq!(*reader.detail_level(), DetailLevel::Full);
            if let Some(spec) = reader.get_spectrum_by_index(10) {
                let spec: Spectrum = spec;
                assert_eq!(spec.index(), 10);
                assert_eq!(spec.id(), "controllerType=0 controllerNumber=1 scan=11");
                assert_eq!(spec.peaks().len(), 941);
            } else {
                panic!("Failed to retrieve spectrum by index")
            }

            assert_eq!(reader.get_spectrum_by_id("controllerType=0 controllerNumber=1 scan=11").unwrap().index(), 10);

            if let Some(spec) = reader.get_spectrum_by_time(0.358558333333) {
                assert_eq!(spec.index(), 34);
            } else {
                panic!("Failed to retrieve spectrum by time")
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
    fn test_dispatch_mzreader() -> io::Result<()> {
        let mut reader = MZReader::open_path("./test/data/small.mzML")?;

        let n = reader.len();
        let n_ms1 = reader.iter().filter(|s| s.ms_level() == 1).count();
        let n_msn = reader.iter().filter(|s| s.ms_level() == 2).count();

        assert_eq!(n, 48);
        assert_eq!(n, n_ms1 + n_msn);
        Ok(())
    }
}
