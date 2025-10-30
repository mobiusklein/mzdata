
use std::collections::HashMap;
use std::io::{self, prelude::*, BufWriter};
use std::marker::PhantomData;
use std::str;


use mzpeaks::{
    peak::KnownCharge, CentroidPeak, DeconvolutedPeak, IntensityMeasurement, MZLocated,
    PeakCollection,
};

use crate::prelude::*;

use crate::meta::{
    DataProcessing, FileDescription, InstrumentConfiguration, MSDataFileMetadata,
    MassSpectrometryRun, Sample, Software,
};
use crate::params::{
    ControlledVocabulary, ParamDescribed, ParamLike, ParamValue as _, CURIE,
};

use crate::spectrum::{
    bindata::BinaryArrayMap,
    IonProperties, Precursor, PrecursorSelection, RefPeakDataLevel, SignalContinuity,
    SpectrumDescription, SpectrumLike,
};


const TITLE_CV: CURIE = ControlledVocabulary::MS.curie(1000796);
const MS_LEVEL_CV: CURIE = ControlledVocabulary::MS.curie(1000511);
const MSN_SPECTRUM_CV: CURIE = ControlledVocabulary::MS.curie(1000580);

/// A trait that controls what additional descriptive entries
/// are written in the spectrum header of an MGF file, not including
/// the essential items like `RTINSECONDS` and `PEPMASS`
pub trait MGFHeaderStyle: Sized {
    fn write_header<
        W: io::Write,
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
        S: SpectrumLike<C, D>,
    >(
        writer: &mut MGFWriterType<W, C, D, Self>,
        spectrum: &S,
    ) -> io::Result<()> {
        let desc = spectrum.description();
        writer.write_kv("SCANS", &desc.index.to_string())?;
        Ok(())
    }

    fn write_precursor<W: io::Write, C: CentroidLike, D: DeconvolutedCentroidLike>(
        writer: &mut MGFWriterType<W, C, D, Self>,
        precursor: &Precursor,
    ) -> io::Result<()> {
        if let Some(ion) = precursor.ion() {
            writer.handle.write_all(b"PEPMASS=")?;
            writer.handle.write_all(ion.mz.to_string().as_bytes())?;
            writer.handle.write_all(b" ")?;
            writer
                .handle
                .write_all(ion.intensity.to_string().as_bytes())?;
            if let Some(charge) = ion.charge {
                writer.handle.write_all(b" ")?;
                writer.handle.write_all(charge.to_string().as_bytes())?;
            }
            writer.handle.write_all(b"\n")?;
            for param in ion
                .params()
                .iter()
                .chain(precursor.activation.params())
            {
                writer.write_param(param)?;
            }

        }

        if let Some(pid) = precursor.precursor_id() {
            writer.handle.write_all(b"PRECURSORSCAN=")?;
            writer.handle.write_all(pid.as_bytes())?;
            writer.handle.write_all(b"\n")?;
        }
        Ok(())
    }
}

/// An MGF style that writes the `SCANS` entry, but no additional
/// descriptions beyond the minimum.
#[derive(Debug, Clone, Copy)]
pub struct SimpleMGFStyle();

impl MGFHeaderStyle for SimpleMGFStyle {}

/// An MGF style that writes the contents of [`SpectrumLike`]'s
/// [`Param`](crate::params::Param) as spectrum header entries. This is the default style.
#[derive(Debug, Clone, Copy)]
pub struct MZDataMGFStyle();

impl MGFHeaderStyle for MZDataMGFStyle {
    fn write_header<
        W: io::Write,
        C: CentroidLike,
        D: DeconvolutedCentroidLike,
        S: SpectrumLike<C, D>,
    >(
        writer: &mut MGFWriterType<W, C, D, Self>,
        spectrum: &S,
    ) -> io::Result<()> {
        let desc = spectrum.description();
        writer.write_kv("NATIVEID", spectrum.id())?;
        writer.write_kv("SCANS", &desc.index.to_string())?;
        for param in desc
            .params()
            .iter()
            .filter(|p| TITLE_CV != **p && MSN_SPECTRUM_CV != **p && MS_LEVEL_CV != **p)
        {
            writer.write_param(param)?;
        }
        Ok(())
    }
}

/// An MGF writer type that only writes centroided MSn spectra.
///
/// To customize the way that spectrum metadata is written, provide
/// a type implementing [`MGFHeaderStyle`]. The default style, [`MZDataMGFStyle`]
/// writes all parameters it can find.
pub struct MGFWriterType<
    W: io::Write,
    C: CentroidLike = CentroidPeak,
    D: DeconvolutedCentroidLike = DeconvolutedPeak,
    Y: MGFHeaderStyle = MZDataMGFStyle,
> {
    pub handle: io::BufWriter<W>,
    pub offset: usize,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
    file_description: FileDescription,
    instrument_configurations: HashMap<u32, InstrumentConfiguration>,
    softwares: Vec<Software>,
    samples: Vec<Sample>,
    data_processings: Vec<DataProcessing>,
    style_type: PhantomData<Y>,
    run: MassSpectrometryRun,
}

impl<W: io::Write, C: CentroidLike, D: DeconvolutedCentroidLike, Y: MGFHeaderStyle>
    MGFWriterType<W, C, D, Y>
{
    pub fn new(file: W) -> MGFWriterType<W, C, D, Y> {
        let handle = io::BufWriter::with_capacity(500, file);
        MGFWriterType {
            handle,
            offset: 0,
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
            file_description: Default::default(),
            instrument_configurations: Default::default(),
            softwares: Default::default(),
            samples: Default::default(),
            data_processings: Default::default(),
            run: Default::default(),
            style_type: PhantomData,
        }
    }

    /// Format a spectrum title similarly to the [Trans-Proteomic Pipeline](https://tools.proteomecenter.org/software.php)
    /// compatibility.
    pub fn make_title<S: SpectrumLike<C, D>>(&self, spectrum: &S) -> String {
        let idx = spectrum.index();
        let charge = spectrum
            .precursor()
            .and_then(|prec| prec.ion().and_then(|i| i.charge()))
            .unwrap_or_default();
        let id = spectrum.id();
        let run_id = self.run_description().and_then(|d| d.id.as_ref());
        let source_file = self.source_file_name();
        match (run_id, source_file) {
            (None, None) => format!("run.{idx}.{idx}.{charge} NativeID:\"{id}\""),
            (None, Some(source_name)) => {
                format!("run.{idx}.{idx}.{charge} SourceFile:\"{source_name}\"")
            }
            (Some(run_id), None) => format!("{run_id}.{idx}.{idx}.{charge} NativeID:\"{id}\""),
            (Some(run_id), Some(source_name)) => format!(
                "{run_id}.{idx}.{idx}.{charge} SourceFile:\"{source_name}\", NativeID:\"{id}\""
            ),
        }
    }

    pub fn into_inner(self) -> BufWriter<W> {
        self.handle
    }

    /// Convert a [`ParamLike`] value into a spectrum header `key=value`
    /// pair.
    ///
    /// Calles [`MGFWriterType::write_kv`].
    pub fn write_param<P: ParamLike>(&mut self, param: &P) -> io::Result<()> {
        self.handle
            .write_all(param.name().to_uppercase().replace(' ', "_").as_bytes())?;
        self.handle.write_all(b"=")?;
        self.handle.write_all(&param.value().as_bytes())?;
        self.handle.write_all(b"\n")?;
        Ok(())
    }

    /// Write a spectrum header `KEY=value`
    pub fn write_kv(&mut self, key: &str, value: &str) -> io::Result<()> {
        self.handle.write_all(key.as_bytes())?;
        self.handle.write_all(b"=")?;
        self.handle.write_all(value.as_bytes())?;
        self.handle.write_all(b"\n")?;
        Ok(())
    }

    fn write_precursor(&mut self, precursor: &Precursor) -> io::Result<()> {
        Y::write_precursor(self, precursor)?;
        Ok(())
    }

    /// Write the header of a spectrum, everything after `BEGIN IONS`, before writing
    /// the peak list.
    pub fn write_header<T: SpectrumLike<C, D>>(&mut self, spectrum: &T) -> io::Result<()> {
        let desc = spectrum.description();
        let (title, _had_title) = desc
            .get_param_by_curie(&TITLE_CV)
            .map(|p| (p.value.clone(), true))
            .unwrap_or_else(|| (self.make_title(spectrum).into(), false));
        self.handle.write_all(&title.as_bytes())?;
        self.handle.write_all(b"\nRTINSECONDS=")?;
        self.handle
            .write_all((spectrum.start_time() * 60.0).to_string().as_bytes())?;
        self.handle.write_all(b"\n")?;
        if let Some(precursor) = &desc.precursor.first() {
            self.write_precursor(precursor)?;
        }

        Y::write_header(self, spectrum)?;
        Ok(())
    }

    fn write_deconvoluted_centroids(&mut self, centroids: &[D]) -> io::Result<()> {
        let mut centroids: Vec<DeconvolutedPeak> =
            centroids.iter().map(|p| p.as_centroid()).collect();
        centroids.sort_by(|a, b| a.mz().total_cmp(&b.mz()));
        for peak in centroids.into_iter() {
            self.handle.write_all(peak.mz().to_string().as_bytes())?;
            self.handle.write_all(b" ")?;
            self.handle
                .write_all(peak.intensity().to_string().as_bytes())?;
            self.handle.write_all(b" ")?;
            self.handle
                .write_all(peak.charge().to_string().as_bytes())?;
            self.handle.write_all(b"\n")?;
        }
        Ok(())
    }

    fn write_centroids(&mut self, centroids: &[C]) -> io::Result<()> {
        for peak in centroids {
            self.handle.write_all(peak.mz().to_string().as_bytes())?;
            self.handle.write_all(b" ")?;
            self.handle
                .write_all(peak.intensity().to_string().as_bytes())?;
            self.handle.write_all(b"\n")?;
        }
        Ok(())
    }

    fn write_arrays(
        &mut self,
        description: &SpectrumDescription,
        arrays: &BinaryArrayMap,
    ) -> io::Result<()> {
        match description.signal_continuity {
            SignalContinuity::Centroid => {
                for (mz, inten) in arrays.mzs()?.iter().zip(arrays.intensities()?.iter()) {
                    self.handle.write_all(mz.to_string().as_bytes())?;
                    self.handle.write_all(b" ")?;
                    self.handle.write_all(inten.to_string().as_bytes())?;
                    self.handle.write_all(b"\n")?;
                }
            }
            SignalContinuity::Profile | SignalContinuity::Unknown => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "MGF spectrum must be centroided",
                ))
            }
        }
        Ok(())
    }

    /// Write the peak list of a spectrum, everything unitl the `END IONS`
    pub fn write_peaks<S: SpectrumLike<C, D> + 'static>(&mut self, spectrum: &S) -> io::Result<()> {
        let description = spectrum.description();
        match spectrum.peaks() {
            RefPeakDataLevel::Missing => {
                log::warn!(
                    "Attempting to write a spectrum without any peak data, {}",
                    description.id
                )
            }
            RefPeakDataLevel::RawData(arrays) => {
                if description.signal_continuity == SignalContinuity::Profile {
                    return Err(io::Error::new(
                        io::ErrorKind::Unsupported,
                        "Cannot write profile spectrum to MGF",
                    ));
                }
                self.write_arrays(description, arrays)?
            }
            RefPeakDataLevel::Centroid(centroids) => {
                self.write_centroids(&centroids[0..centroids.len()])?
            }
            RefPeakDataLevel::Deconvoluted(deconvoluted) => {
                self.write_deconvoluted_centroids(&deconvoluted[0..deconvoluted.len()])?
            }
        }
        Ok(())
    }

    /// Write a spectrum from start to finish. It will skip spectra where `ms_level() == 1`
    pub fn write<S: SpectrumLike<C, D> + 'static>(&mut self, spectrum: &S) -> io::Result<usize> {
        let description = spectrum.description();
        if description.ms_level == 1 {
            log::warn!(
                "Attempted to write an MS1 spectrum to MGF, {}, skipping.",
                description.id
            );
            return Ok(0);
        }
        // spectrum
        self.handle.write_all(
            br#"BEGIN IONS
TITLE="#,
        )?;
        self.write_header(spectrum)?;
        self.write_peaks(spectrum)?;
        self.handle.write_all(b"END IONS\n")?;
        Ok(0)
    }
}

impl<W: io::Write, C: CentroidLike, D: DeconvolutedCentroidLike, Y: MGFHeaderStyle>
    MSDataFileMetadata for MGFWriterType<W, C, D, Y>
{
    crate::impl_metadata_trait!();

    fn run_description(&self) -> Option<&MassSpectrometryRun> {
        Some(&self.run)
    }

    fn run_description_mut(&mut self) -> Option<&mut MassSpectrometryRun> {
        Some(&mut self.run)
    }
}

impl<W: io::Write, C: CentroidLike + 'static, D: DeconvolutedCentroidLike + 'static>
    SpectrumWriter<C, D> for MGFWriterType<W, C, D>
{
    fn write<S: SpectrumLike<C, D> + 'static>(&mut self, spectrum: &S) -> io::Result<usize> {
        if spectrum.ms_level() != 1 {
            self.write(spectrum)
        } else {
            log::trace!("Skipping writing MS1 spectrum {} to MGF", spectrum.id());
            Ok(0)
        }
    }

    fn write_group<
        S: SpectrumLike<C, D> + 'static,
        G: super::super::SpectrumGrouping<C, D, S> + 'static,
    >(
        &mut self,
        group: &G,
    ) -> io::Result<usize> {
        let mut c = 0;
        for s in group.products() {
            c += self.write(s)?;
        }
        Ok(c)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.handle.flush()
    }

    fn close(&mut self) -> io::Result<()> {
        self.handle.flush()
    }
}

/// A convenient alias for [`MGFWriterType`] with the peak types specified
pub type MGFWriter<W> = MGFWriterType<W, CentroidPeak, DeconvolutedPeak, MZDataMGFStyle>;
pub type SimpleMGFWriter<W> = MGFWriterType<W, CentroidPeak, DeconvolutedPeak, SimpleMGFStyle>;
pub type SimpleMGFWriterType<W, C, D> = MGFWriterType<W, C, D, SimpleMGFStyle>;