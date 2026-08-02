use std::{fmt::Debug, future::Future, io, pin::Pin, task::Poll};

use futures::{AsyncReadExt, AsyncSeek, FutureExt, io::AsyncRead, lock::Mutex};
use gloo_file::{
    FileReadError,
    callbacks::{FileReader, read_as_bytes as read_as_bytes_cb},
    futures::read_as_bytes,
};
use js_sys::Uint8Array;
use tokio::io::{AsyncRead as TokioAsyncRead, AsyncSeek as TokioAsyncSeek};
use wasm_bindgen::{JsValue, prelude::*};
use web_sys::{Blob, File};

#[derive(Debug, Clone)]
pub enum BufferHandle {
    Blob(gloo_file::Blob),
    File(gloo_file::File),
}

impl From<Blob> for BufferHandle {
    fn from(value: Blob) -> Self {
        Self::Blob(gloo_file::Blob::from(value))
    }
}

impl From<File> for BufferHandle {
    fn from(value: File) -> Self {
        Self::File(gloo_file::File::from(value))
    }
}

impl BufferHandle {
    pub fn size(&self) -> u64 {
        let val = match self {
            BufferHandle::Blob(handle) => handle.size(),
            BufferHandle::File(handle) => handle.size(),
        };
        val
    }

    pub fn slice(&self, start: u64, end: u64) -> Self {
        match self {
            Self::Blob(handle) => Self::Blob(handle.slice(start, end)),
            Self::File(handle) => Self::File(handle.slice(start, end)),
        }
    }

    pub fn slice_with_f64_and_f64(&self, start: f64, end: f64) -> Result<Blob, JsValue> {
        match self {
            BufferHandle::Blob(handle) => {
                let inner: &web_sys::Blob = handle.as_ref();
                inner.slice_with_f64_and_f64(start, end)
            }
            BufferHandle::File(handle) => {
                let inner: &web_sys::File = handle.as_ref();
                inner.slice_with_f64_and_f64(start, end)
            }
        }
    }

    pub fn read(self) -> impl Future<Output = Result<Vec<u8>, FileReadError>> {
        match self {
            BufferHandle::Blob(blob) => read_as_bytes(&blob),
            BufferHandle::File(file) => read_as_bytes(&file),
        }
    }

    pub fn read_callback<F>(&self, callback: F) -> FileReader
    where
        F: FnOnce(Result<Vec<u8>, FileReadError>) + 'static,
    {
        match self {
            BufferHandle::Blob(blob) => read_as_bytes_cb(blob, callback),
            BufferHandle::File(file) => read_as_bytes_cb(file, callback),
        }
    }
}

enum State {
    Idle(Option<Vec<u8>>),
    Busy(
        (
            Operation,
            Vec<u8>,
            Option<Pin<Box<dyn Future<Output = Result<Result<Vec<u8>, FileReadError>, futures::channel::oneshot::Canceled>>>>>,
        ),
    ),
}

impl Debug for State {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Idle(arg0) => f.debug_tuple("Idle").field(arg0).finish(),
            Self::Busy(arg0) => f
                .debug_tuple("Busy")
                .field(&arg0.0)
                .field(&arg0.1)
                .field(&if arg0.2.is_some() { "<callback>" } else { "-" })
                .finish(),
        }
    }
}

#[derive(Debug)]
enum Operation {
    Read,
}

#[derive(Debug)]
struct Inner {
    state: State,
    pos: u64,
}

#[derive(Debug)]
#[wasm_bindgen]
pub struct WebIO {
    handle: BufferHandle,
    position: u64,
    state: Mutex<Inner>,
}

impl Clone for WebIO {
    fn clone(&self) -> Self {
        let state = Inner {
            pos: self.position,
            state: State::Idle(None),
        };
        Self {
            handle: self.handle.clone(),
            position: self.position.clone(),
            state: Mutex::new(state),
        }
    }
}

#[wasm_bindgen]
impl WebIO {
    pub fn blob(blob: web_sys::Blob) -> Self {
        Self::new(blob)
    }

    pub fn file(file: web_sys::File) -> Self {
        Self::new(file)
    }

    #[wasm_bindgen(getter, js_name = "size")]
    pub fn js_size(&self) -> u64 {
        self.size()
    }

    pub async fn read(&mut self, n: usize) -> Result<Uint8Array, JsError> {
        let mut data = vec![0u8; n];
        let result = <Self as AsyncReadExt>::read(self, &mut data).await;
        let v = match result {
            Ok(_z) => Ok(Uint8Array::new_from_slice(&data)),
            Err(e) => Err(JsError::from(e)),
        };
        v
    }
}

impl WebIO {
    pub fn new<B: Into<BufferHandle>>(handle: B) -> Self {
        Self {
            handle: handle.into(),
            position: 0,
            state: Mutex::new(Inner {
                state: State::Idle(None),
                pos: 0,
            }),
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

    pub fn set_position(&mut self, pos: io::SeekFrom) {
        let size = self.handle.size();
        match pos {
            io::SeekFrom::Start(offset) => {
                if offset > size {
                    self.position = size;
                } else {
                    self.position = offset;
                }
            }
            io::SeekFrom::End(offset) => {
                let offset_u = offset.unsigned_abs();
                if offset < 0 {
                    if offset_u > size {
                        self.position = 0;
                    } else {
                        self.position = size - offset_u
                    }
                } else {
                    self.position = size;
                }
            }
            io::SeekFrom::Current(offset) => {
                let offset_u = offset.unsigned_abs();
                if offset < 0 {
                    if offset_u > self.position {
                        self.position = 0;
                    } else {
                        self.position = self.position - offset_u
                    }
                } else {
                    self.position = size;
                }
            }
        }
        self.state.get_mut().pos = self.position;
    }

    pub async fn read_direct(&mut self, buffer: &mut [u8]) -> Result<u64, JsValue> {
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

impl AsyncRead for WebIO {
    fn poll_read(
        self: Pin<&mut Self>,
        cx: &mut std::task::Context<'_>,
        buf: &mut [u8],
    ) -> Poll<io::Result<usize>> {
        let me = self.get_mut();
        if let Some(mut inner) = me.state.try_lock() {
            log::info!("Polling read, state: {:?}", inner.state);
            match &mut inner.state {
                State::Idle(items) => {
                    if let Some(readbuf) = items.take() {
                        log::info!("Found buffer");
                        if !readbuf.is_empty() || buf.len() == 0 {
                            buf.copy_from_slice(&readbuf);
                            return Poll::Ready(Ok(readbuf.len()));
                        }
                    }
                    log::info!("In idle branch, no buffer");
                    let v = me.handle.slice(me.position, me.position + buf.len() as u64);
                    log::info!("Slicing {}-{}", me.position, me.position + buf.len() as u64);
                    let (send, recv) = futures::channel::oneshot::channel();
                    let task = async move || {
                        log::info!("In callback");
                        let buf = v.read().await;
                        log::info!("Read completed in callback");
                        send.send(buf).unwrap();
                    };
                    log::info!("Spawning the task locally");
                    wasm_bindgen_futures::spawn_local(task());
                    log::info!("Updating state to busy");
                    inner.state = State::Busy((Operation::Read, Vec::new(), Some(Box::pin(recv))));
                    log::info!("Returning pending");
                    Poll::Pending
                }
                State::Busy((_op, _opbuf, task)) => match Pin::new(task.as_mut().unwrap()).poll(cx) {
                    Poll::Ready(value) => {
                        match value {
                            Ok(value) => {
                                log::info!("In busy branch, buffer read: {value:?}");
                                let value = value.map_err(|e| io::Error::other(e))?;
                                buf.copy_from_slice(&value);
                                me.position += value.len() as u64;
                                inner.pos += value.len() as u64;
                                inner.state = State::Idle(None);
                                Poll::Ready(Ok(value.len()))
                            },
                            Err(futures::channel::oneshot::Canceled) => {
                                inner.state = State::Idle(None);
                                Poll::Ready(Err(io::Error::new(io::ErrorKind::ConnectionAborted, "The read operation was cancelled")))
                            }
                        }
                    }
                    Poll::Pending => Poll::Pending,
                },
            }
        } else {
            log::info!("State locked, pending");
            Poll::Pending
        }
    }
}

impl AsyncSeek for WebIO {
    fn poll_seek(
        mut self: Pin<&mut Self>,
        _cx: &mut std::task::Context<'_>,
        pos: io::SeekFrom,
    ) -> Poll<io::Result<u64>> {
        self.set_position(pos);
        Poll::Ready(Ok(self.position))
    }
}

impl TokioAsyncSeek for WebIO {
    fn start_seek(self: Pin<&mut Self>, position: io::SeekFrom) -> io::Result<()> {
        let me = self.get_mut();
        let inner = me.state.get_mut();
        match &mut inner.state {
            State::Idle(_) => {
                me.set_position(position);
                Ok(())
            }
            State::Busy((op, _buf, _task)) => Err(io::Error::new(
                io::ErrorKind::ResourceBusy,
                format!("WebIO is currently executing another operation {op:?}"),
            )),
        }
    }

    fn poll_complete(
        self: Pin<&mut Self>,
        _cx: &mut std::task::Context<'_>,
    ) -> Poll<io::Result<u64>> {
        Poll::Ready(Ok(self.position))
    }
}

#[wasm_bindgen]
pub async fn read_all(handle: &mut WebIO) -> Vec<u8> {
    let mut buf = Vec::new();
    handle.read_to_end(&mut buf).await.unwrap();
    buf
}

// impl TokioAsyncRead for WebIO {
//     fn poll_read(
//         self: Pin<&mut Self>,
//         cx: &mut std::task::Context<'_>,
//         buf: &mut tokio::io::ReadBuf<'_>,
//     ) -> Poll<io::Result<()>> {
//         let me = self.get_mut();
//         let inner = me.state.get_mut();
//         match &mut inner.state {
//             State::Idle(items) => {
//                 if let Some(readbuf) = items.take() {
//                     if !readbuf.is_empty() || buf.capacity() == 0 {
//                         buf.put_slice(&readbuf);
//                         return Poll::Ready(Ok(()));
//                     }
//                 }
//                 let v = me
//                     .handle
//                     .slice(me.position, me.position + buf.capacity() as u64);
//                 let task = v.read().boxed_local();
//                 inner.state = State::Busy((Operation::Read, Vec::new(), Some(task)));
//                 Poll::Pending
//             }
//             State::Busy((_op, _opbuf, task)) => match Pin::new(task.as_mut().unwrap()).poll(cx) {
//                 Poll::Ready(value) => {
//                     let value = value.map_err(|e| io::Error::other(e))?;
//                     buf.put_slice(&value);
//                     me.position += value.len() as u64;
//                     inner.pos += value.len() as u64;
//                     inner.state = State::Idle(None);
//                     Poll::Ready(Ok(()))
//                 }
//                 Poll::Pending => Poll::Pending,
//             },
//         }
//     }
// }
