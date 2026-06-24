use std::{future::Future, io, pin::Pin, task::Poll};

use futures::{io::AsyncRead, lock::Mutex, AsyncReadExt, FutureExt};
use tokio::io::AsyncRead as TokioAsyncRead;
use gloo_file::futures::read_as_bytes;
use js_sys::Number;
use wasm_bindgen::{prelude::wasm_bindgen, JsValue};
use web_sys::{Blob, File};

pub enum BufferHandle {
    Blob(Blob),
    File(File),
}

impl From<Blob> for BufferHandle {
    fn from(value: Blob) -> Self {
        Self::Blob(value)
    }
}

impl From<File> for BufferHandle {
    fn from(value: File) -> Self {
        Self::File(value)
    }
}

fn f64_to_u64_safe(val: f64) -> Option<u64> {
    if 0.0 <= val && val <= Number::MAX_SAFE_INTEGER {
        Some(val as u64)
    } else {
        None
    }
}

// fn u64_to_f64(val: u64) -> Option<f64> {
//     let val = val as f64;
//     if val <= Number::MAX_SAFE_INTEGER {
//         Some(val)
//     } else {
//         None
//     }
// }

impl BufferHandle {
    pub fn size(&self) -> u64 {
        let val = match self {
            BufferHandle::Blob(handle) => handle.size(),
            BufferHandle::File(handle) => handle.size(),
        };
        f64_to_u64_safe(val).expect("Could not convert buffer size to valid integer")
    }

    pub fn slice_with_f64_and_f64(&self, start: f64, end: f64) -> Result<Blob, JsValue> {
        match self {
            BufferHandle::Blob(handle) => handle.slice_with_f64_and_f64(start, end),
            BufferHandle::File(handle) => handle.slice_with_f64_and_f64(start, end),
        }
    }
}

#[wasm_bindgen]
pub struct WebIO {
    handle: BufferHandle,
    position: u64,
}

impl WebIO {
    pub fn new<B: Into<BufferHandle>>(handle: B) -> Self {
        Self {
            handle: handle.into(),
            position: 0,
        }
    }

    pub fn size(&self) -> u64 {
        self.handle.size()
    }

    pub fn slice(&self, length: usize) -> Result<Blob, JsValue> {
        let n_bytes_total = self.handle.size();
        let n_bytes_remaining = n_bytes_total - self.position;
        let bytes_to_read = (length as u64).min(n_bytes_remaining);

        let blob = self
            .handle
            .slice_with_f64_and_f64(self.position as f64, bytes_to_read as f64)?
            .into();
        Ok(blob)
    }

    pub async fn read(&mut self, buffer: &mut [u8]) -> Result<u64, JsValue> {
        let n_bytes_total = self.handle.size();
        let n_bytes_remaining = n_bytes_total - self.position;
        let bytes_to_read = (buffer.len() as u64).min(n_bytes_remaining);

        let blob = self
            .handle
            .slice_with_f64_and_f64(self.position as f64, bytes_to_read as f64)?
            .into();

        let buf = gloo_file::futures::read_as_bytes(&blob)
            .await
            .map_err(|e| JsValue::from_str(&e.to_string()))?;

        let z = buf.len() as u64;
        self.position += z;
        buffer.copy_from_slice(&buf);
        Ok(z)
    }
}

enum State {
    Idle(Option<Vec<u8>>),
    Busy(Pin<Box<dyn Future<Output = Result<Vec<u8>, gloo_file::FileReadError>>>>),
    Complete,
}

struct Inner {
    state: Option<State>,
}

impl Inner {}

pub struct WebReaderAsyncRead {
    stream_reader: WebIO,
    inner: Mutex<Inner>,
}

impl WebReaderAsyncRead {
    pub fn new(stream_reader: WebIO) -> Self {
        Self { stream_reader, inner: Mutex::new(Inner { state: Some(State::Idle(None))}) }
    }
}

impl AsyncRead for WebReaderAsyncRead {
    fn poll_read(
        self: std::pin::Pin<&mut Self>,
        cx: &mut std::task::Context<'_>,
        buf: &mut [u8],
    ) -> std::task::Poll<std::io::Result<usize>> {
        let me = self.get_mut();
        let inner = me.inner.get_mut();

        loop {
            let state = inner.state.take().unwrap();

            match state {
                State::Idle(bytes_mut) => {
                    if let Some(bytes_read) = bytes_mut {
                        buf.copy_from_slice(&bytes_read);
                        inner.state = Some(State::Idle(None));
                        return Poll::Ready(Ok(bytes_read.len()));
                    } else {
                        let blob = match me.stream_reader.slice(buf.len()) {
                            Ok(blob) => gloo_file::Blob::from(blob),
                            Err(e) => {
                                return Poll::Ready(Err(io::Error::new(
                                    io::ErrorKind::Other,
                                    e.as_string().unwrap(),
                                )))
                            }
                        };
                        let t = read_as_bytes(&blob).boxed_local();
                        inner.state = Some(State::Busy(t))
                    }
                }
                State::Busy(mut js_future) => {
                    log::info!("Busy, polling future");
                    let j = js_future.as_mut().poll(cx);
                    match j {
                        Poll::Ready(read) => match read {
                            Ok(data) => {
                                log::info!("Read successful: {}", data.len());
                                me.stream_reader.position += data.len() as u64;
                                inner.state = Some(State::Idle(Some(data)));
                            }
                            Err(e) => {
                                log::info!("Read Failed: {e}");
                                let err = match e {
                                    gloo_file::FileReadError::AbortedEarly => {
                                        io::Error::new(io::ErrorKind::ConnectionAborted, e)
                                    }
                                    gloo_file::FileReadError::NotFound(_) => {
                                        io::Error::new(io::ErrorKind::NotFound, e)
                                    }
                                    gloo_file::FileReadError::NotReadable(_) => {
                                        io::Error::new(io::ErrorKind::Unsupported, e)
                                    }
                                    gloo_file::FileReadError::Security(_) => {
                                        io::Error::new(io::ErrorKind::ConnectionRefused, e)
                                    }
                                };
                                inner.state = Some(State::Complete);
                                return Poll::Ready(Err(err));
                            }
                        },
                        Poll::Pending => {
                            inner.state = Some(State::Busy(js_future));
                        }
                    }
                }
                State::Complete => {
                    inner.state = Some(State::Complete);
                    return Poll::Ready(Ok(0));
                }
            }
        }
    }
}


impl TokioAsyncRead for WebReaderAsyncRead {
    fn poll_read(
        self: Pin<&mut Self>,
        cx: &mut std::task::Context<'_>,
        buf: &mut tokio::io::ReadBuf<'_>,
    ) -> Poll<io::Result<()>> {
        let me = self.get_mut();
        let inner = me.inner.get_mut();

        loop {
            let state = inner.state.take().unwrap();

            match state {
                State::Idle(bytes_mut) => {
                    if let Some(bytes_read) = bytes_mut {
                        buf.initialize_unfilled().copy_from_slice(&bytes_read);
                        inner.state = Some(State::Idle(None));
                        return Poll::Ready(Ok(()));
                    } else {
                        let blob = match me.stream_reader.slice(buf.capacity()) {
                            Ok(blob) => gloo_file::Blob::from(blob),
                            Err(e) => {
                                return Poll::Ready(Err(io::Error::new(
                                    io::ErrorKind::Other,
                                    e.as_string().unwrap(),
                                )))
                            }
                        };
                        let t = read_as_bytes(&blob).boxed_local();
                        inner.state = Some(State::Busy(t))
                    }
                }
                State::Busy(mut js_future) => {
                    log::info!("Busy, polling future");
                    let j = js_future.as_mut().poll(cx);
                    match j {
                        Poll::Ready(read) => match read {
                            Ok(data) => {
                                log::info!("Read successful: {}", data.len());
                                me.stream_reader.position += data.len() as u64;
                                inner.state = Some(State::Idle(Some(data)));
                            }
                            Err(e) => {
                                log::info!("Read unsuccessful: {}", e);
                                let err = match e {
                                    gloo_file::FileReadError::AbortedEarly => {
                                        io::Error::new(io::ErrorKind::ConnectionAborted, e)
                                    }
                                    gloo_file::FileReadError::NotFound(_) => {
                                        io::Error::new(io::ErrorKind::NotFound, e)
                                    }
                                    gloo_file::FileReadError::NotReadable(_) => {
                                        io::Error::new(io::ErrorKind::Unsupported, e)
                                    }
                                    gloo_file::FileReadError::Security(_) => {
                                        io::Error::new(io::ErrorKind::ConnectionRefused, e)
                                    }
                                };
                                inner.state = Some(State::Complete);
                                return Poll::Ready(Err(err));
                            }
                        }
                        Poll::Pending => {
                            inner.state = Some(State::Busy(js_future));
                        }
                    }
                }
                State::Complete => {
                    inner.state = Some(State::Complete);
                    return Poll::Ready(Ok(()));
                }
            }
        }
    }
}


#[wasm_bindgen]
pub async fn test_reader_blob(stream_reader: File) -> String {
    let mut handle = WebReaderAsyncRead::new(WebIO::new(stream_reader));
    let mut buf = String::new();
    handle.read_to_string(&mut buf).await.unwrap();
    buf
}