use std::env;
use std::fs;
use std::io;
use std::process::exit;
use std::thread;
use std::time;

use std::sync::mpsc::sync_channel;

#[cfg(feature = "mzmlb")]
use mzdata::io::mzmlb;

#[cfg(feature = "thermorawfilereader")]
use mzdata::io::ThermoRawReader;

use mzdata::io::MassSpectrometryReadWriteProcess;
use mzdata::io::{
    infer_format, infer_from_path, infer_from_stream, MassSpectrometryFormat, PreBufferedStream,
};
use mzdata::prelude::*;
use mzdata::{MGFReader, MGFWriter, MzMLReader, MzMLWriter};

use env_logger;
use mzpeaks::CentroidPeak;
use mzpeaks::DeconvolutedPeak;

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
        self.reader_then()
    }

    fn reader_then(&self) -> io::Result<()> {
        if self.inpath == "-" {
            let mut stream = PreBufferedStream::new(io::stdin())?;
            let (ms_format, _compressed) = infer_from_stream(&mut stream)?;
            match ms_format {
                MassSpectrometryFormat::MGF => self.writer_then(MGFReader::new(stream))?,
                MassSpectrometryFormat::MzML => {
                    self.writer_then(MzMLReader::new(io::BufReader::new(stream)))?
                }
                _ => {
                    eprintln!("Could not infer input format from STDIN");
                    exit(1)
                }
            }
        } else {
            let (ms_format, _compressed) = infer_format(&self.inpath)?;
            match ms_format {
                MassSpectrometryFormat::MGF => {
                    let reader = MGFReader::open_path(&self.inpath)?;
                    self.writer_then(reader)?;
                }
                MassSpectrometryFormat::MzML => {
                    let reader = MzMLReader::open_path(&self.inpath)?;
                    self.writer_then(reader)?;
                }
                #[cfg(feature = "mzmlb")]
                MassSpectrometryFormat::MzMLb => {
                    let reader = mzmlb::MzMLbReader::open_path(&self.inpath)?;
                    self.writer_then(reader)?;
                }
                #[cfg(feature = "thermorawfilereader")]
                MassSpectrometryFormat::ThermoRaw => {
                    let reader = ThermoRawReader::open_path(&self.inpath)?;
                    self.writer_then(reader)?;
                }
                _ => {
                    eprintln!("Could not infer input format from {}", self.inpath);
                    exit(1)
                }
            }
        };
        Ok(())
    }

    fn writer_then<R: ScanSource + MSDataFileMetadata + Send + 'static>(
        &self,
        reader: R,
    ) -> io::Result<()> {
        match infer_from_path(&self.outpath).0 {
            MassSpectrometryFormat::MGF => {
                let mut writer =
                    MGFWriter::new(io::BufWriter::new(fs::File::create(&self.outpath)?));
                writer.copy_metadata_from(&reader);
                self.task(reader, writer)?;
            }
            MassSpectrometryFormat::MzML => {
                let mut writer =
                    MzMLWriter::new(io::BufWriter::new(fs::File::create(&self.outpath)?));
                writer.copy_metadata_from(&reader);
                self.task(reader, writer)?;
            }
            #[cfg(feature = "mzmlb")]
            MassSpectrometryFormat::MzMLb => {
                let mut writer = mzmlb::MzMLbWriterBuilder::new(&self.outpath)
                    .with_zlib_compression(9)
                    .create()?;
                writer.copy_metadata_from(&reader);
                self.task(reader, writer)?;
            }
            _ => {
                eprintln!("Could not infer output format from {}", self.outpath);
                exit(1)
            }
        }
        Ok(())
    }

    fn task<R: ScanSource + Send + 'static, W: ScanWriter + Send + 'static>(
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
            + ScanSource<CentroidPeak, DeconvolutedPeak>
            + Send
            + 'static,
        W: ScanWriter<CentroidPeak, DeconvolutedPeak> + Send + 'static,
    >(
        &self,
        reader: R,
        writer: W,
    ) -> Result<(), Self::ErrorType> {
        self.task(reader, writer)
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
