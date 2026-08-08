use std::any::Any;
use std::fs::File;
use std::io;
use std::path::PathBuf;
use std::sync::{
    atomic::{AtomicU64, Ordering as AtomicOrdering},
    mpsc::sync_channel,
    Arc,
};
use std::thread;
use std::time;

use clap::Parser;

use log::info;
use mzdata::io::MassSpectrometryFormat;
use mzdata::io::{checksum_file, MassSpectrometryReadWriteProcess, Sink, Source};
use mzdata::meta::custom_software_name;
use mzdata::meta::Software;
use mzdata::meta::{DataProcessing, ProcessingMethod, SourceFile};
use mzdata::params::ControlledVocabulary;
use mzdata::prelude::*;

use mzdata::spectrum::bindata::BinaryCompressionType;
use mzdata::spectrum::ArrayType;
use mzdata::spectrum::ArrayType::IntensityArray;
use mzdata::spectrum::ArrayType::MZArray;
use mzdata::spectrum::BinaryDataArrayType;
use mzdata::MzMLWriter;
use mzpeaks::{CentroidPeak, DeconvolutedPeak};

fn compression_parser(compression: &str) -> Result<BinaryCompressionType, String> {
    let compression = if !compression.ends_with(" compression") {
        format!("{compression} compression")
    } else {
        compression.to_string()
    };

    BinaryCompressionType::COMPRESSION_METHODS
        .iter()
        .find(|x| x.as_param().unwrap().name() == compression)
        .copied()
        .ok_or_else(|| compression.to_string())
}

#[derive(Debug, Clone, Parser)]
pub struct MZConvert {
    #[arg()]
    pub inpath: String,
    #[arg()]
    pub outpath: String,

    #[arg(short = 'b', long, default_value_t = 8192)]
    pub buffer_size: usize,

    #[arg(long, value_parser=compression_parser, default_value="zlib compression")]
    pub mz_compression: BinaryCompressionType,

    #[arg(long, value_parser=compression_parser, default_value="zlib compression")]
    pub intensity_compression: BinaryCompressionType,

    #[arg(long, value_parser=compression_parser, default_value="zlib compression")]
    pub ion_mobility_compression: BinaryCompressionType,
}

impl MZConvert {
    pub fn new(inpath: String, outpath: String) -> Self {
        Self {
            inpath,
            outpath,
            buffer_size: 8192,
            mz_compression: BinaryCompressionType::Zlib,
            intensity_compression: BinaryCompressionType::Zlib,
            ion_mobility_compression: BinaryCompressionType::Zlib,
        }
    }

    pub fn main(&self) -> io::Result<()> {
        let source = if self.inpath == "-" {
            Source::Stdin
        } else {
            Source::<_, _>::from(self.inpath.as_ref())
        };

        let sink = if self.outpath == "-" {
            Sink::Writer(Box::new(io::stdout()), MassSpectrometryFormat::MzML)
        } else {
            Sink::<CentroidPeak, DeconvolutedPeak>::from(self.outpath.as_ref())
        };

        self.open_reader(source, sink)
    }

    fn task<R: SpectrumSource + Send + 'static, W: SpectrumWriter + Send + 'static>(
        &self,
        reader: R,
        mut writer: W,
    ) -> io::Result<()> {
        let (send, recv) = sync_channel(self.buffer_size);
        let buffered = Arc::new(AtomicU64::default());
        let buffered_w = Arc::clone(&buffered);
        let reader_handle = thread::spawn(move || {
            reader.enumerate().for_each(|(i, s)| {
                let waiting_cnt = buffered.fetch_add(1, AtomicOrdering::SeqCst);
                if i % 5000 == 0 && i > 0 {
                    log::info!("Reading {} {} ({waiting_cnt} to write)", i, s.id());
                }
                if i % 100 == 0 && i > 0 {
                    log::debug!("Reading {} {} ({waiting_cnt} to write)", i, s.id());
                }

                send.send(s).unwrap()
            });
        });

        let writer_handle = thread::spawn(move || {
            for s in recv.iter() {
                let i = s.index();
                buffered_w.fetch_sub(1, AtomicOrdering::SeqCst);
                writer
                    .write_owned(s)
                    .inspect_err(|e| log::error!("Failed to write spectrum {i}: {e}"))
                    .unwrap();
            }
            writer.close().unwrap();
        });

        reader_handle.join().unwrap();
        writer_handle.join().unwrap();
        Ok(())
    }

    #[allow(unused)]
    fn configure_reader(&self, reader: &mut dyn Any) {
        #[cfg(feature = "thermo")]
        if let Some(reader) = reader.downcast_mut::<mzdata::io::thermo::ThermoRawReader>() {
            reader.set_load_extended_spectrum_data(false);
        }
        #[cfg(feature = "bruker_tdf")]
        if let Some(reader) = reader.downcast_mut::<mzdata::io::tdf::TDFSpectrumReader>() {
            reader.set_consolidate_peaks(false);
        }
    }

    fn configure_writer(&self, writer: &mut dyn Any) {
        info!("Configuring writer...");
        if let Some(writer) = writer.downcast_mut::<MzMLWriter<io::BufWriter<File>>>() {
            log::debug!("Configuring compression methods: {:?}", self.mz_compression);
            writer.set_compression_method(
                MZArray,
                BinaryDataArrayType::Float32,
                self.mz_compression,
            );
            writer.set_compression_method(
                MZArray,
                BinaryDataArrayType::Float64,
                self.mz_compression,
            );

            log::debug!(
                "Configuring compression methods: {:?}",
                self.intensity_compression
            );
            writer.set_compression_method(
                IntensityArray,
                BinaryDataArrayType::Float32,
                self.intensity_compression,
            );
            writer.set_compression_method(
                IntensityArray,
                BinaryDataArrayType::Float64,
                self.intensity_compression,
            );

            log::debug!(
                "Configuring compression methods: {:?}",
                self.ion_mobility_compression
            );

            let im_arrays = [
                ArrayType::MeanInverseReducedIonMobilityArray,
                ArrayType::RawDriftTimeArray,
                ArrayType::MeanDriftTimeArray,
                ArrayType::RawIonMobilityArray,
                ArrayType::MeanIonMobilityArray,
            ];

            for im_array in im_arrays {
                writer.set_compression_method(
                    im_array.clone(),
                    BinaryDataArrayType::Float64,
                    self.ion_mobility_compression,
                );
                writer.set_compression_method(
                    im_array,
                    BinaryDataArrayType::Float32,
                    self.ion_mobility_compression,
                );
            }
        }
    }
}

impl MassSpectrometryReadWriteProcess<CentroidPeak, DeconvolutedPeak> for MZConvert {
    type ErrorType = io::Error;

    fn task<
        R: RandomAccessSpectrumIterator<CentroidPeak, DeconvolutedPeak>
            + SpectrumSource<CentroidPeak, DeconvolutedPeak>
            + Send
            + Any
            + 'static,
        W: SpectrumWriter<CentroidPeak, DeconvolutedPeak> + Send + Any + 'static,
    >(
        &self,
        reader: R,
        writer: W,
    ) -> Result<(), Self::ErrorType> {
        self.task(reader, writer)
    }

    fn transform_reader<
        R: RandomAccessSpectrumIterator<CentroidPeak, DeconvolutedPeak>
            + MSDataFileMetadata
            + SpectrumSource<CentroidPeak, DeconvolutedPeak>
            + Send
            + Any
            + 'static,
    >(
        &self,
        mut reader: R,
        _format: MassSpectrometryFormat,
    ) -> Result<R, Self::ErrorType> {
        self.configure_reader(&mut reader);
        Ok(reader)
    }

    fn transform_writer<
        R: RandomAccessSpectrumIterator<CentroidPeak, DeconvolutedPeak>
            + MSDataFileMetadata
            + SpectrumSource<CentroidPeak, DeconvolutedPeak>
            + Send
            + Any
            + 'static,
        W: SpectrumWriter<CentroidPeak, DeconvolutedPeak> + MSDataFileMetadata + Send + Any + 'static,
    >(
        &self,
        reader: R,
        reader_format: MassSpectrometryFormat,
        mut writer: W,
        writer_format: MassSpectrometryFormat,
    ) -> Result<(R, W), Self::ErrorType> {
        self.configure_writer(&mut writer);
        if self.inpath != "-" {
            let pb: PathBuf = self.inpath.clone().into();
            let checksum = if pb.is_dir() {
                Default::default()
            } else {
                info!("Computing checksum for {}", pb.display());

                checksum_file(&pb)?
            };
            let has_already = reader
                .file_description()
                .source_files
                .iter()
                .flat_map(|f| {
                    f.get_param_by_name("SHA-1")
                        .map(|c| c.value.as_str() == checksum)
                })
                .all(|a| a);

            if !has_already {
                let mut sf = SourceFile {
                    location: pb
                        .parent()
                        .map(|p| format!("file:///{}", p.to_string_lossy()))
                        .unwrap_or("file:///".to_string()),
                    name: pb
                        .file_name()
                        .map(|p| p.to_string_lossy().to_string())
                        .unwrap_or("".to_string()),
                    ..Default::default()
                };
                let par = ControlledVocabulary::MS.param_val(1000569, "SHA-1", checksum);
                sf.add_param(par);
                sf.file_format = reader_format.as_param();

                if let Some(ref_sf) = reader.file_description().source_files.last() {
                    sf.id_format = ref_sf.id_format.clone()
                }
                writer.file_description_mut().source_files.push(sf);
            }
        };

        let mut sw = Software {
            version: format!(
                "v{}",
                option_env!("CARGO_PKG_VERSION").unwrap_or_else(|| "?")
            ),
            id: Software::find_unique_id("mzconvert", writer.softwares()),
            ..Default::default()
        };
        sw.add_param(custom_software_name("mzconvert"));

        let mut method = ProcessingMethod {
            software_reference: sw.id.clone(),
            ..Default::default()
        };
        writer.softwares_mut().push(sw);

        if let Some(conv) = writer_format.as_conversion() {
            method.add_param(conv.into())
        }

        if writer.data_processings().is_empty() {
            let mut dp = DataProcessing::default();
            method.order = 0;
            dp.push(method.clone());
            dp.id = "DP1".to_string();
            writer.data_processings_mut().push(dp)
        } else {
            for dp in writer.data_processings_mut() {
                let mut next_step = method.clone();
                next_step.order = dp.highest_order() + 1;
                dp.push(next_step);
            }
        }
        Ok((reader, writer))
    }
}

fn main() -> io::Result<()> {
    env_logger::init();

    let job = MZConvert::parse();
    let start = time::Instant::now();
    job.main()?;
    let end = time::Instant::now();
    let elapsed = end - start;
    eprintln!("Conversion finished: {:0.2} seconds", elapsed.as_secs_f64());
    Ok(())
}
