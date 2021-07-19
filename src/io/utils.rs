#![allow(dead_code)]

use std::path;
use std::io;

type ByteBuffer = io::Cursor<Vec<u8>>;


#[derive(Debug, Clone)]
pub enum FileSource {
    FileSystem(path::PathBuf),
    InMemory(ByteBuffer),
    Stream
}

impl Default for FileSource {
    fn default() -> FileSource {
        FileSource::InMemory(io::Cursor::new(Vec::new()))
    }
}


#[derive(Debug, Clone, Default)]
pub struct FileDescription {
    pub source: FileSource
}


impl<'lifespan> FileDescription {

    pub fn from_path<T>(path: T) -> FileDescription where T: Into<path::PathBuf> {
        FileDescription {
            source: FileSource::FileSystem(path.into()),
        }
    }

    pub fn from_buffer(buffer: ByteBuffer) -> FileDescription {
        FileDescription {
            source: FileSource::InMemory(buffer)
        }
    }

    pub fn file_name(&self) -> Option<&path::Path> {
        match &self.source {
            FileSource::FileSystem(path) => Some(path),
            FileSource::InMemory(_stream) => None,
            FileSource::Stream => None
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::io::Read;

    #[test]
    fn test_from_buffer() {
        let mut buff: Vec<u8> = Vec::new();
        buff.extend(b"foobar");
        let stream = ByteBuffer::new(buff);
        let mut out: Vec<u8> = Vec::new();
        let desc = FileDescription::from_buffer(stream);
        assert!(matches!(desc.file_name(), None));
        if let FileSource::InMemory(mut buff) = desc.source {
            buff.read_to_end(&mut out).unwrap();
            assert_eq!(out, b"foobar");
        }
    }
}