use std::env;
use std::io;
use std::path::PathBuf;
use std::process::exit;
use std::thread;
use std::time;
use std::sync::mpsc::sync_channel;

use mzdata::io::{
    Sink, Source, MassSpectrometryReadWriteProcess,
    checksum_file
};
use mzdata::meta::SourceFile;
use mzdata::params::ControlledVocabulary;
use mzdata::prelude::*;

use env_logger;
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
        } else {Source::<_, _>::from(self.inpath.as_ref())};
        let sink = Sink::<CentroidPeak, DeconvolutedPeak>::from(self.outpath.as_ref());
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
            + 'static,
        W: SpectrumWriter<CentroidPeak, DeconvolutedPeak> + Send + 'static,
    >(
        &self,
        reader: R,
        writer: W,
    ) -> Result<(), Self::ErrorType> {
        self.task(reader, writer)
    }

    #[allow(unused)]
    fn transform_writer<
            R: RandomAccessSpectrumIterator<CentroidPeak, DeconvolutedPeak> + MSDataFileMetadata + SpectrumSource<CentroidPeak, DeconvolutedPeak> + Send + 'static,
            W: SpectrumWriter<CentroidPeak, DeconvolutedPeak> + MSDataFileMetadata + Send + 'static,
        >(
            &self,
            reader: R,
            reader_format: mzdata::io::MassSpectrometryFormat,
            mut writer: W,
            writer_format: mzdata::io::MassSpectrometryFormat,
        ) -> Result<(R, W), Self::ErrorType> {
        if self.inpath != "-" {
            let pb: PathBuf = self.inpath.clone().into();
            let checksum = checksum_file(&pb)?;
            let has_already = reader.file_description().source_files.iter().flat_map(|f| f.get_param_by_name("SHA-1").map(|c| c.value == checksum)).all(|a| a);
            if !has_already {
                let mut sf = SourceFile::default();
                sf.location = pb.parent().map(|p| format!("file://{}", p.to_string_lossy())).unwrap_or("file://".to_string());
                sf.name = pb.file_name().map(|p| p.to_string_lossy().to_string()).unwrap_or("".to_string());
                let par = ControlledVocabulary::MS.param_val(1000569u32, "SHA-1", checksum);
                sf.add_param(par);
                sf.file_format = reader_format.as_param();

                if let Some(ref_sf) = reader.file_description().source_files.last() {
                    sf.id_format = ref_sf.id_format.clone()
                }
                writer.file_description_mut().source_files.push(sf);
            }
        };
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
