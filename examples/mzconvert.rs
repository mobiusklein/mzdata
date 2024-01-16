use std::env;
use std::fs;
use std::io;
use std::process::exit;
use std::thread;

use std::time;

use std::sync::mpsc::sync_channel;

#[cfg(feature = "mzmlb")]
use mzdata::io::mzmlb;

use mzdata::io::{
    infer_format, infer_from_path, infer_from_stream,
    MassSpectrometryFormat, PreBufferedStream,
};
use mzdata::prelude::*;
use mzdata::{MGFReader, MGFWriter, MzMLReader, MzMLWriter};

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
                MassSpectrometryFormat::MzML => self.writer_then(MzMLReader::new(io::BufReader::new(stream)))?,
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
                let writer = MGFWriter::new(io::BufWriter::new(fs::File::create(&self.outpath)?));
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

    fn task<R: ScanSource + Send + 'static, W: ScanWriter<'static> + Send + 'static>(
        &self,
        reader: R,
        mut writer: W,
    ) -> io::Result<()> {
        let (send, recv) = sync_channel(2usize.pow(16));

        let reader_handle = thread::spawn(move || {
            reader.for_each(|s| send.send(s).unwrap());
        });

        let writer_handle = thread::spawn(move || {
            for s in recv.iter() {
                writer.write(&s).unwrap();
            }
            writer.close().unwrap();
        });

        reader_handle.join().unwrap();
        writer_handle.join().unwrap();
        Ok(())
    }
}

fn main() -> io::Result<()> {
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
