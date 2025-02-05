use std::any::Any;
use std::env;
use std::io;
use std::path::PathBuf;
use std::process::exit;
use std::sync::mpsc::sync_channel;
use std::thread;
use std::time;

use log::info;
use mzdata::io::MassSpectrometryFormat;
use mzdata::io::{checksum_file, MassSpectrometryReadWriteProcess, Sink, Source};
use mzdata::meta::custom_software_name;
use mzdata::meta::Software;
use mzdata::meta::{DataProcessing, ProcessingMethod, SourceFile};
use mzdata::params::ControlledVocabulary;
use mzdata::prelude::*;

use mzpeaks::{CentroidPeak, DeconvolutedPeak};

#[derive(Debug, Clone)]
pub struct MZConvert {
    pub inpath: String,
    pub outpath: String,
}

impl MZConvert {
    pub fn new(inpath: String, outpath: String) -> Self {
        Self { inpath, outpath }
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
        let (send, recv) = sync_channel(2usize.pow(14));

        let reader_handle = thread::spawn(move || {
            reader.enumerate().for_each(|(i, s)| {
                if i % 10000 == 0 && i > 0 {
                    log::info!("Reading {} {}", i, s.id());
                }
                send.send(s).unwrap()
            });
        });

        let writer_handle = thread::spawn(move || {
            for s in recv.iter() {
                writer.write_owned(s).unwrap();
            }
            writer.close().unwrap();
        });

        reader_handle.join().unwrap();
        writer_handle.join().unwrap();
        Ok(())
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
        if self.inpath != "-" {
            let pb: PathBuf = self.inpath.clone().into();
            info!("Computing checksum for {}", pb.display());
            let checksum = checksum_file(&pb)?;
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
    let inpath = env::args().nth(1).unwrap_or_else(|| {
        eprintln!("Please provide a path to read an MS data file from, or '-'");
        exit(1)
    });

    let outpath = env::args().nth(2).unwrap_or_else(|| {
        eprintln!("Please provide a path to write an MS file to, or '-'");
        exit(1)
    });
    let start = time::Instant::now();
    let job = MZConvert::new(inpath, outpath);
    job.main()?;
    let end = time::Instant::now();
    let elapsed = end - start;
    eprintln!("Conversion finished: {:0.2} seconds", elapsed.as_secs_f64());
    Ok(())
}
