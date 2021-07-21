#![allow(dead_code)]

use std::io;
use std::path;

type ByteBuffer = io::Cursor<Vec<u8>>;

#[derive(Debug, Clone)]
pub enum FileSource<T: io::Read> {
    FileSystem(path::PathBuf),
    Stream(T),
    Empty,
}

impl<T: io::Read> Default for FileSource<T> {
    fn default() -> FileSource<T> {
        FileSource::Empty
    }
}

#[derive(Debug, Clone, Default)]
pub struct FileDescription<T: io::Read> {
    pub source: FileSource<T>,
}

impl<'lifespan, T: io::Read> FileDescription<T> {
    pub fn from_path<P>(path: P) -> FileDescription<T>
    where
        P: Into<path::PathBuf>,
    {
        FileDescription {
            source: FileSource::FileSystem(path.into()),
        }
    }

    pub fn from_stream(stream: T) -> FileDescription<T> {
        FileDescription {
            source: FileSource::Stream(stream),
        }
    }

    pub fn file_name(&self) -> Option<&path::Path> {
        match &self.source {
            FileSource::FileSystem(path) => Some(path),
            FileSource::Stream(_stream) => None,
            FileSource::Empty => None,
        }
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
        let desc = FileDescription::<ByteBuffer>::from_stream(stream);
        assert!(matches!(desc.file_name(), None));
        if let FileSource::Stream(mut buff) = desc.source {
            buff.read_to_end(&mut out).unwrap();
            assert_eq!(out, b"foobar");
        }
    }
}
