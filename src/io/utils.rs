#![allow(dead_code)]

use std::fs;
use std::io;
use std::path;

use md5::Context as MD5Context;
use md5::Digest;

type ByteBuffer = io::Cursor<Vec<u8>>;

#[derive(Debug, Clone, Default)]
pub enum FileWrapper<T: io::Read> {
    FileSystem(path::PathBuf),
    Stream(T),
    #[default]
    Empty,
}

/// Controls the level of spectral detail read from an MS data file
#[derive(Debug, Default, Clone, Copy, Hash, PartialEq, Eq)]
pub enum DetailLevel {
    #[default]
    /// Read all spectral data, including peak data, eagerly decoding it. This is the default
    Full,
    /// Read all spectral data, including peak data but defer decoding until later if possible
    Lazy,
    /// Read only the metadata of spectra, ignoring peak data entirely
    MetadataOnly,
}

#[derive(Debug, Clone, Default)]
pub struct FileSource<T: io::Read> {
    pub source: FileWrapper<T>,
}

// This really should be a full file-like object abstraction, but that
// feels like it is beyond the scope of this crate. Something like
// https://github.com/bnjjj/chicon-rs
impl<'lifespan, T: io::Read> FileSource<T> {
    pub fn from_path<P>(path: P) -> FileSource<T>
    where
        P: Into<path::PathBuf>,
    {
        FileSource {
            source: FileWrapper::FileSystem(path.into()),
        }
    }

    pub fn from_stream(stream: T) -> FileSource<T> {
        FileSource {
            source: FileWrapper::Stream(stream),
        }
    }

    pub fn file_name(&self) -> Option<&path::Path> {
        match &self.source {
            FileWrapper::FileSystem(path) => Some(path),
            FileWrapper::Stream(_stream) => None,
            FileWrapper::Empty => None,
        }
    }

    pub fn index_file_name(&self) -> Option<path::PathBuf> {
        match &self.source {
            FileWrapper::Empty => None,
            FileWrapper::Stream(_stream) => None,
            FileWrapper::FileSystem(path) => {
                if let Some(stem) = path.file_name() {
                    if let Some(parent) = path.parent() {
                        let base = parent.join(stem);
                        let name = base.with_extension("index.json");
                        return Some(name);
                    }
                }
                None
            }
        }
    }

    pub fn has_index_file(&self) -> bool {
        match self.index_file_name() {
            Some(path) => path.exists(),
            None => false,
        }
    }
}

pub fn from_path<P>(path: P) -> FileSource<fs::File>
where
    P: Into<path::PathBuf>,
{
    FileSource::from_path(path)
}

impl<T, P> From<P> for FileSource<T>
where
    P: Into<path::PathBuf>,
    T: io::Read,
{
    fn from(path: P) -> FileSource<T> {
        FileSource::from_path(path)
    }
}

/// A writable stream that keeps a running MD5 checksum of all bytes
#[derive(Clone)]
pub(crate) struct MD5HashingStream<T: io::Write> {
    pub stream: T,
    pub context: MD5Context,
}

impl<T: io::Write> MD5HashingStream<T> {
    pub fn new(file: T) -> MD5HashingStream<T> {
        Self {
            stream: file,
            context: MD5Context::new(),
        }
    }

    pub fn compute(&self) -> Digest {
        self.context.clone().compute()
    }

    pub fn get_mut(&mut self) -> &mut T {
        &mut self.stream
    }

    pub fn into_inner(self) -> T {
        self.stream
    }
}

impl<T: io::Write> io::Write for MD5HashingStream<T> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.context.consume(buf);
        self.stream.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.stream.flush()
    }
}

impl<T: io::Seek + io::Write> io::Seek for MD5HashingStream<T> {
    fn seek(&mut self, pos: io::SeekFrom) -> io::Result<u64> {
        self.stream.seek(pos)
    }
}

pub struct PreBufferedStream<R: io::Read> {
    stream: R,
    buffer: io::Cursor<Vec<u8>>,
    buffer_size: usize,
    position: usize,
}

impl<R: io::Read> io::Seek for PreBufferedStream<R> {
    fn seek(&mut self, pos: io::SeekFrom) -> io::Result<u64> {
        match pos {
            io::SeekFrom::Start(offset) => {
                if self.position > self.buffer_size {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Seeking after leaving buffered prefix",
                    ));
                } else if self.position + offset as usize > self.buffer_size {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Cannot seeking beyond buffered prefix",
                    ));
                }
                self.position = offset as usize;
                let r = self.buffer.seek(pos);
                log::trace!("{pos:?} Position {0} -> {1}: {r:?}", self.position, self.buffer.stream_position().unwrap());
                r
            }
            io::SeekFrom::End(_) => {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "Cannot seek relative the end of PreBufferedStream",
                ))
            }
            io::SeekFrom::Current(offset) => {
                if self.position > self.buffer_size {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Seeking after leaving buffered prefix",
                    ));
                }
                if offset < 0 {
                    if offset.abs() as usize > self.position {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "Cannot seek to negative position",
                        ));
                    } else {
                        self.position = self.position.saturating_sub(offset.abs() as usize);
                        let r = self.buffer.seek(io::SeekFrom::Start(self.position as u64));
                        log::trace!("{pos:?} Position {0} -> {1}: {r:?}", self.position, self.buffer.stream_position().unwrap());
                        r
                    }
                } else {
                    if offset as usize + self.position > self.buffer_size {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "Cannot seeking beyond buffered prefix",
                        ));
                    } else {
                        let r = self.buffer.seek(io::SeekFrom::Current(offset));
                        log::trace!("{pos:?} Position {0} -> {1}: {r:?}", self.position, self.buffer.stream_position().unwrap());
                        r
                    }
                }
            }
        }
    }

    fn stream_position(&mut self) -> io::Result<u64> {
        Ok(self.position as u64)
    }
}

impl<R: io::Read> io::Read for PreBufferedStream<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n_total = buf.len();
        let before = self.position;
        let (n_remaining, n_from_buffer) = if self.position < self.buffer_size {
            let n_from_buffer = self.buffer.read(buf)?;
            self.position += n_from_buffer;
            (n_total.saturating_sub(n_from_buffer), n_from_buffer)
        } else {
            (n_total, 0)
        };
        if n_remaining > 0 {
            let n_rest = self.stream.read(&mut buf[n_from_buffer..])?;
            self.position += n_rest;
        }
        let total_read = self.position - before;
        Ok(total_read)
    }
}

const BUFFER_SIZE: usize = 2usize.pow(16);

impl<R: io::Read> PreBufferedStream<R> {
    pub fn new(stream: R) -> io::Result<Self> {
        Self::new_with_buffer_size(stream, BUFFER_SIZE)
    }

    pub fn new_with_buffer_size(stream: R, buffer_size: usize) -> io::Result<Self> {
        let buffer = io::Cursor::new(Vec::with_capacity(buffer_size));
        let mut inst = Self {
            stream,
            buffer_size,
            buffer,
            position: 0,
        };
        inst.prefill_buffer()?;
        Ok(inst)
    }

    fn prefill_buffer(&mut self) -> io::Result<usize> {
        let mut buffer = self.buffer.get_mut();
        buffer.resize(self.buffer_size, 0);
        let bytes_read = self.stream.read(&mut buffer)?;
        buffer.shrink_to(bytes_read);
        self.buffer_size = bytes_read;
        Ok(bytes_read)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::io::prelude::*;

    #[test]
    fn test_from_buffer() {
        let mut buff: Vec<u8> = Vec::new();
        buff.extend(b"foobar");
        let stream = ByteBuffer::new(buff);
        let mut out: Vec<u8> = Vec::new();
        let desc = FileSource::<ByteBuffer>::from_stream(stream);
        assert!(desc.file_name().is_none());
        if let FileWrapper::Stream(mut buff) = desc.source {
            buff.read_to_end(&mut out).unwrap();
            assert_eq!(out, b"foobar");
        }
    }

    #[test]
    fn test_prebuffering() -> io::Result<()> {
        let mut fh = fs::File::open("./test/data/batching_test.mzML")?;
        let mut data = Vec::new();
        fh.read_to_end(&mut data)?;
        let content = io::Cursor::new(data);
        let mut stream = PreBufferedStream::new_with_buffer_size(content, 512)?;

        assert_eq!(stream.buffer_size, 512);

        let mut buffer = [0u8; 128];
        stream.read(&mut buffer)?;
        assert_eq!(buffer.len(), 128);
        assert!(buffer.starts_with(b"<?xml version=\"1.0\" encoding=\"utf-8\"?>"));

        let mut buffer2 = [0u8; 128];
        stream.seek(io::SeekFrom::Start(0))?;
        stream.read(&mut buffer2)?;

        assert_eq!(buffer, buffer2);

        assert!(stream.seek(io::SeekFrom::Start(556)).is_err());

        Ok(())
    }
}
