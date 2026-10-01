#![allow(dead_code)]

use std::fs;
use std::io;
#[allow(unused)]
use std::io::prelude::*;
#[allow(unused)]
use std::path::PathBuf;

#[cfg(feature = "checksum")]
use sha1::{self, Digest as _};

type ByteBuffer = io::Cursor<Vec<u8>>;

/// Controls the level of spectral detail read from an MS data file
#[derive(Debug, Default, Clone, Copy, Hash, PartialEq, Eq)]
pub enum DetailLevel {
    #[default]
    /// Read all spectral data, including peak data, eagerly decoding it. This is the default
    Full,
    /// Read all spectral data, including peak data but defer decoding until later if possible.
    /// Check a format reader's documentation to see if it supports lazy loading. Lazy loading
    /// is only really of value for dense profile mode data or very, very long peak lists
    /// that are large and expensive to decode.
    Lazy,
    /// Read only the metadata of spectra, ignoring peak data entirely
    MetadataOnly,
}

/// A wrapper around an [`io::Read`] to provide limited [`io::Seek`] access even if the
/// underlying stream does not support it. It retains the first *n* bytes as they are
/// read and permits seeking within the retained prefix until reading beyond that range.
///
/// This is useful for working with [`io::stdin`] or a network stream.
pub struct PreBufferedStream<R: io::Read> {
    stream: R,
    buffer: io::Cursor<Vec<u8>>,
    buffer_size: usize,
    position: u64,
}

impl<R: io::Read> io::Seek for PreBufferedStream<R> {
    fn seek(&mut self, pos: io::SeekFrom) -> io::Result<u64> {
        let target = match pos {
            io::SeekFrom::Start(offset) => offset,
            io::SeekFrom::End(_) => {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "Cannot seek relative the end of PreBufferedStream",
                ))
            }
            io::SeekFrom::Current(offset) => {
                self.position.checked_add_signed(offset).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        if offset < 0 {
                            "Cannot seek to negative position"
                        } else {
                            "Cannot seeking beyond buffered prefix"
                        },
                    )
                })?
            }
        };
        let retained = self.buffer.get_ref().len() as u64;
        if self.position > retained {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Seeking after leaving buffered prefix",
            ));
        }
        if target > retained {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Cannot seeking beyond buffered prefix",
            ));
        }
        let before = self.position;
        self.buffer.seek(io::SeekFrom::Start(target))?;
        self.position = target;
        log::trace!("{pos:?} Position {before} -> {target}");
        Ok(target)
    }

    fn stream_position(&mut self) -> io::Result<u64> {
        Ok(self.position)
    }
}

impl<R: io::Read> io::Read for PreBufferedStream<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        let buffered = self.buffer.read(buf)?;
        if buffered > 0 {
            self.position += buffered as u64;
            return Ok(buffered);
        }
        let n = self.stream.read(buf)?;
        if self.position < self.buffer_size as u64 {
            let retained = n.min(self.buffer_size - self.position as usize);
            self.buffer.get_mut().extend_from_slice(&buf[..retained]);
            self.buffer.set_position(self.buffer.get_ref().len() as u64);
        }
        self.position += n as u64;
        Ok(n)
    }
}

const BUFFER_SIZE: usize = 2usize.pow(16);

impl<R: io::Read> PreBufferedStream<R> {
    /// Create a new pre-buffered stream wrapping `stream` with a buffer size of 2<sup>16</sup> bytes.
    ///
    /// This method fails if the initial read fails.
    pub fn new(stream: R) -> io::Result<Self> {
        Self::new_with_buffer_size(stream, BUFFER_SIZE)
    }

    /// Create a new pre-buffered stream wrapping `stream` with a buffer size of `buffer_size` bytes.
    ///
    /// This method fails if the initial read fails.
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
        if self.buffer_size == 0 {
            return Ok(0);
        }
        let buffer = self.buffer.get_mut();
        buffer.resize(self.buffer_size, 0);
        let bytes_read = loop {
            match self.stream.read(buffer) {
                Err(err) if err.kind() == io::ErrorKind::Interrupted => continue,
                result => break result?,
            }
        };
        buffer.truncate(bytes_read);
        Ok(bytes_read)
    }
}

#[cfg(feature = "checksum")]
/// Compute a SHA-1 digest of a file path
pub fn checksum_file(path: &PathBuf) -> io::Result<String> {
    let mut checksum = sha1::Sha1::new();
    let mut reader = io::BufReader::new(fs::File::open(path)?);
    let mut buf = vec![0; 2usize.pow(20)];
    while let Ok(i) = reader.read(&mut buf) {
        if i == 0 {
            break;
        }
        checksum.update(&buf[..i]);
    }
    Ok(hex::encode(checksum.finalize()))
}

#[cfg(feature = "checksum")]
/// A writable stream that keeps a running SHA-1 checksum of all bytes
#[derive(Clone)]
pub(crate) struct SHA1HashingStream<T: io::Write> {
    pub stream: T,
    pub context: sha1::Sha1,
}

#[cfg(feature = "checksum")]
impl<T: io::Write> SHA1HashingStream<T> {
    pub fn new(file: T) -> SHA1HashingStream<T> {
        Self {
            stream: file,
            context: sha1::Sha1::new(),
        }
    }

    pub fn compute(&self) -> sha1::Sha1 {
        self.context.clone()
    }

    pub fn get_mut(&mut self) -> &mut T {
        &mut self.stream
    }

    pub fn into_inner(self) -> T {
        self.stream
    }
}


#[cfg(feature = "checksum")]
impl<T: io::Write> io::Write for SHA1HashingStream<T> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.context.update(buf);
        self.stream.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.stream.flush()
    }
}

#[cfg(feature = "checksum")]
impl<T: io::Seek + io::Write> io::Seek for SHA1HashingStream<T> {
    fn seek(&mut self, pos: io::SeekFrom) -> io::Result<u64> {
        self.stream.seek(pos)
    }
}


#[cfg(feature = "parallelism")]
mod parallelism {
    use rayon::prelude::*;
    use crate::prelude::*;

    use super::*;

    /// A helper type to load spectra concurrently across multiple threads.
    /// Requires the reader type implement [`MZFileReader`].
    ///
    /// # Note
    /// This helper is still too low level. Expect a higher level API to eventually become available.
    pub struct ConcurrentLoader {
        path: PathBuf,
        num_threads: Option<usize>,
    }

    impl ConcurrentLoader {
        pub fn new(path: PathBuf, num_threads: Option<usize>) -> Self {
            Self { path, num_threads }
        }

        /// Do the actual concurrent loading
        pub fn load<F: MZFileReader<C, D, S>, C: CentroidLike, D: DeconvolutedCentroidLike, S: SpectrumLike<C, D> + Send>(self) -> io::Result<Vec<S>> {
            let guide = F::open_path(&self.path)?;
            let n = guide.len();

            let num_threads = self.num_threads.unwrap_or_else(|| rayon::max_num_threads());

            let task = || -> Vec<S>{
                let mut chunks: Vec<_> = (0..n).into_par_iter().chunks((n / num_threads / 3).max(10)).map(|ii| {
                    let start = ii[0];
                    let mut local_reader = F::open_path(&self.path).unwrap();
                    let spectra: Vec<_> = ii.into_iter().flat_map(|i| local_reader.get_spectrum_by_index(i)).collect();
                    (start ,spectra)
                }).collect();
                chunks.par_sort_by(|a, b| a.0.cmp(&b.0));
                chunks.into_iter().map(|(_, chunk)| chunk).flatten().collect()
            };

            let out = if let Some(num_threads) = self.num_threads {
                let pool = rayon::ThreadPoolBuilder::new().num_threads(num_threads).thread_name(|i| format!("mzdata-concurrent-loader-{i}")).build().unwrap();
                pool.install(|| task())
            } else {
                task()
            };

            Ok(out)
        }
    }
}


#[cfg(feature = "parallelism")]
pub use parallelism::ConcurrentLoader;


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_prebuffering() -> io::Result<()> {
        let mut fh = fs::File::open("./test/data/batching_test.mzML")?;
        let mut data = Vec::new();
        fh.read_to_end(&mut data)?;
        let content = io::Cursor::new(data);
        let mut stream = PreBufferedStream::new_with_buffer_size(content, 512)?;

        assert_eq!(stream.buffer_size, 512);

        let mut buffer = [0u8; 128];
        stream.read_exact(&mut buffer)?;
        assert_eq!(buffer.len(), 128);
        assert!(buffer.starts_with(b"<?xml version=\"1.0\" encoding=\"utf-8\"?>"));

        let mut buffer2 = [0u8; 128];
        stream.seek(io::SeekFrom::Start(0))?;
        stream.read_exact(&mut buffer2)?;

        assert_eq!(buffer, buffer2);

        assert!(stream.seek(io::SeekFrom::Start(556)).is_err());

        Ok(())
    }

    #[test]
    fn test_prebuffer_short_reads() -> io::Result<()> {
        let input = b"abcdefghijklmnop";
        for size in [0, 2, 4, 8, 64] {
            let source = io::Cursor::new(&input[..2]).chain(io::Cursor::new(&input[2..]));
            let mut stream = PreBufferedStream::new_with_buffer_size(source, size)?;
            let mut output = Vec::new();
            stream.read_to_end(&mut output)?;
            assert_eq!(output, input);
            assert!(stream.buffer.get_ref().len() <= size);
        }

        let source = io::Cursor::new(&input[..2]).chain(io::Cursor::new(&input[2..]));
        let mut stream = PreBufferedStream::new_with_buffer_size(source, 8)?;
        assert!(stream.seek(io::SeekFrom::Start(4)).is_err());
        let mut prefix = [0; 4];
        stream.read_exact(&mut prefix)?;
        assert_eq!(stream.seek(io::SeekFrom::Start(0))?, 0);
        stream.read_exact(&mut prefix)?;
        assert_eq!(&prefix, b"abcd");
        Ok(())
    }

    #[test]
    fn test_prebuffer_seek() -> io::Result<()> {
        let mut stream =
            PreBufferedStream::new_with_buffer_size(io::Cursor::new(b"abcdefghijkl"), 8)?;
        let mut prefix = [0; 2];
        stream.read_exact(&mut prefix)?;
        assert_eq!(stream.seek(io::SeekFrom::Current(2))?, 4);
        assert_eq!(stream.stream_position()?, 4);
        assert_eq!(stream.seek(io::SeekFrom::Start(6))?, 6);
        assert_eq!(stream.seek(io::SeekFrom::Current(-1))?, 5);
        for pos in [
            io::SeekFrom::Start(u64::MAX),
            io::SeekFrom::Current(i64::MAX),
            io::SeekFrom::Current(i64::MIN),
            io::SeekFrom::End(0),
        ] {
            assert!(stream.seek(pos).is_err());
            assert_eq!(stream.stream_position()?, 5);
        }
        let mut next = [0];
        stream.read_exact(&mut next)?;
        assert_eq!(&next, b"f");
        assert_eq!(stream.seek(io::SeekFrom::Start(8))?, 8);
        stream.read_exact(&mut next)?;
        assert_eq!(&next, b"i");
        assert!(stream.seek(io::SeekFrom::Start(0)).is_err());
        assert_eq!(stream.stream_position()?, 9);
        Ok(())
    }

    #[test]
    fn test_prebuffer_pending_source() -> io::Result<()> {
        struct PendingSource {
            prefix: io::Cursor<&'static [u8]>,
            interrupted: bool,
        }
        impl io::Read for PendingSource {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                if self.interrupted {
                    self.interrupted = false;
                    return Err(io::ErrorKind::Interrupted.into());
                }
                match self.prefix.read(buf)? {
                    0 => Err(io::ErrorKind::WouldBlock.into()),
                    n => Ok(n),
                }
            }
        }
        let source = PendingSource {
            prefix: io::Cursor::new(b"abcd"),
            interrupted: true,
        };
        let mut stream = PreBufferedStream::new_with_buffer_size(source, 8)?;
        assert_eq!(stream.read(&mut [])?, 0);
        let mut bytes = [0; 8];
        assert_eq!(stream.read(&mut bytes)?, 4);
        assert_eq!(&bytes[..4], b"abcd");
        assert_eq!(
            stream.read(&mut bytes).unwrap_err().kind(),
            io::ErrorKind::WouldBlock
        );
        assert_eq!(stream.stream_position()?, 4);
        stream.seek(io::SeekFrom::Start(0))?;
        assert_eq!(stream.read(&mut bytes)?, 4);
        assert_eq!(&bytes[..4], b"abcd");
        Ok(())
    }

    #[cfg(feature = "parallelism")]
    #[test]
    fn test_parallel_load() -> io::Result<()> {
        use crate::prelude::*;

        let loader= ConcurrentLoader::new("./test/data/batching_test.mzML".into(), Some(4));
        let spectra = loader.load::<crate::MzMLReader<fs::File>, _, _, _>()?;
        assert_eq!(spectra.len(), 2232);

        let _ = spectra.iter().fold(0, |last, spec| {
            assert_eq!(last, spec.index());
            spec.index() + 1
        });

        Ok(())
    }
}
