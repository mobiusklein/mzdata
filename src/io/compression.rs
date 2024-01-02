use std::mem::swap;
use std::{io, path};

use flate2::bufread::MultiGzDecoder;
use std::io::prelude::*;

pub fn is_gzipped(header: &[u8]) -> bool {
    header.starts_with(b"\x1f\x8b")
}

pub fn is_gzipped_extension(path: path::PathBuf) -> (bool, path::PathBuf) {
    if let Some(ext) = path.extension() {
        if ext.to_ascii_lowercase() == "gz" {
            (true, path.with_extension(""))
        } else {
            (false, path)
        }
    } else {
        (false, path)
    }
}

#[allow(unused)]
pub struct RestartableGzDecoder<R: BufRead + Seek> {
    handle: Option<MultiGzDecoder<R>>,
    offset: u64,
}

#[allow(unused)]
impl<R: BufRead + Seek> RestartableGzDecoder<R> {
    pub fn new(handle: R) -> Self {
        Self {
            handle: Some(MultiGzDecoder::new(handle)),
            offset: 0,
        }
    }

    fn reset(&mut self) -> io::Result<u64> {
        let mut handle = None;
        swap(&mut self.handle, &mut handle);
        let handle = handle.unwrap();
        let mut inner = handle.into_inner();
        let res = inner.seek(io::SeekFrom::Start(0));
        self.handle = Some(MultiGzDecoder::new(inner));
        self.offset = 0;
        res
    }
}

impl<R: BufRead + Seek> Read for RestartableGzDecoder<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let handle = self.handle.as_mut().unwrap();
        match handle.read(buf) {
            Ok(b) => {
                self.offset += b as u64;
                Ok(b)
            }
            Err(e) => Err(e),
        }
    }
}

impl<R: BufRead + Seek> Seek for RestartableGzDecoder<R> {
    fn seek(&mut self, pos: io::SeekFrom) -> io::Result<u64> {
        match pos {
            io::SeekFrom::Start(o) => {
                self.reset()?;
                let mut buf = Vec::new();
                buf.resize(o as usize, 0);
                self.read_exact(&mut buf)?;
                Ok(o)
            }
            io::SeekFrom::End(_) => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "Cannot seek relative to end of a gzip stream",
            )),
            io::SeekFrom::Current(o) => {
                if o == 0 {
                    return Ok(self.offset);
                } else if o < 0 {
                    if o.abs() as u64 > self.offset {
                        Err(io::Error::new(
                            io::ErrorKind::Unsupported,
                            "Cannot earlier than the start of the stream",
                        ))
                    } else {
                        self.seek(io::SeekFrom::Start((self.offset as i64 + o) as u64))
                    }
                } else {
                    let mut buf = Vec::new();
                    buf.resize(o as usize, 0);
                    self.read_exact(&mut buf)?;
                    Ok(self.offset)
                }
            }
        }
    }
}
