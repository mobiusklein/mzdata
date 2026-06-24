use std::future::Future;
use std::io::{Cursor, self};
use std::task::Poll;
use std::pin::Pin;

use bytes::BytesMut;
use futures::Stream;
use js_sys::{Reflect, Uint8Array};
use log::info;
use tokio::io::AsyncRead;
use wasm_bindgen::prelude::*;
use wasm_bindgen_futures::JsFuture;

use tokio::{io::{simplex, AsyncReadExt, AsyncWriteExt, ReadHalf, SimplexStream, WriteHalf}, sync::Mutex};

use web_sys::{ReadableStream, ReadableStreamDefaultReader};

#[pin_project::pin_project]
pub struct WebReaderPipe {
    stream_reader: ReadableStreamDefaultReader,
    write_half: WriteHalf<SimplexStream>,
    #[pin]
    read_half: ReadHalf<SimplexStream>,
    stream_done: bool,
    buffer_done: bool,
}

impl WebReaderPipe {
    fn new(
        stream_reader: ReadableStreamDefaultReader,
        write_half: WriteHalf<SimplexStream>,
        read_half: ReadHalf<SimplexStream>,
    ) -> Self {
        Self {
            stream_reader,
            write_half,
            read_half,
            stream_done: false,
            buffer_done: false,
        }
    }

    pub fn from_stream(stream: ReadableStream) -> WebReaderPipe {
        let reader = stream.get_reader();
        info!("Retrieving stream reader");
        let reader: ReadableStreamDefaultReader = reader.dyn_into().unwrap();
        info!("Allocating pipe");
        let (read_half, write_half) = simplex(2usize.pow(20u32));
        info!("Pipe Ready");
        Self::new(reader, write_half, read_half)
    }

    pub fn into_stream_reader(
        self,
    ) -> tokio_util::io::StreamReader<
        impl Stream<Item = std::io::Result<bytes::BytesMut>>,
        bytes::BytesMut,
    > {
        let chunk_stream = futures::stream::unfold(self, |mut state| async move {
            if state.buffer_done {
                return None;
            }
            state.pump().await;
            let mut buf = bytes::BytesMut::with_capacity(65536);
            let yielded = state.read_half.read(&mut buf).await;
            match yielded {
                Ok(z) => {
                    if z == 0 && state.stream_done {
                        state.buffer_done = true;
                    }
                    Some((Ok(buf), state))
                }
                Err(e) => Some((Err(e), state)),
            }
        });
        tokio_util::io::StreamReader::new(chunk_stream)
    }

    async fn pump(&mut self) {
        let pinned = Pin::new(self);
        pinned.pump_pinned().await
    }

    async fn pump_pinned(self: std::pin::Pin<&mut Self>) {
        if self.stream_done {
            return;
        }
        let next_chunk = JsFuture::from(self.stream_reader.read());
        let this = self.project();
        if let Some(chunk) = next_chunk.await.ok() {
            let done_key = JsValue::from_str("done");
            let value_key = JsValue::from_str("value");

            info!("Checking if chunk is good");
            if !(Reflect::has(&chunk, &done_key).unwrap()
                && Reflect::has(&chunk, &value_key).unwrap())
            {
                return;
            }

            info!("Updating stream state");
            let value = Reflect::get(&chunk, &done_key).unwrap();
            *this.stream_done = value.is_falsy();

            info!("Fetching new buffer");
            let value = Reflect::get(&chunk, &value_key).unwrap();
            if value.is_falsy() {
                return;
            }

            info!("Coercing buffer");
            let value = Uint8Array::new(&value);

            info!("Writing {} bytes to simplex", value.length());
            this.write_half.write_all(&value.to_vec()).await.unwrap();
        }
    }
}

enum State {
    Idle(Option<BytesMut>),
    Busy(JsFuture, BytesMut),
    Complete
}


struct Inner {
    state: Option<State>,
    buffer: Cursor<Vec<u8>>,
    stream_done: bool,
}

impl Inner {

}

pub struct WebReaderAsyncRead {
    stream_reader: ReadableStreamDefaultReader,
    inner: Mutex<Inner>,
}

impl WebReaderAsyncRead {
    pub fn new(stream_reader: ReadableStreamDefaultReader) -> Self {
        let inner = Mutex::new(Inner {
            state: Some(State::Idle(Some(BytesMut::new()))),
            buffer: Cursor::default(),
            stream_done: false,
        });

        Self { stream_reader, inner }
    }
}

impl AsyncRead for WebReaderAsyncRead {
    fn poll_read(
        self: Pin<&mut Self>,
        cx: &mut std::task::Context<'_>,
        dst: &mut tokio::io::ReadBuf<'_>,
    ) -> std::task::Poll<io::Result<()>> {

        let me = self.get_mut();
        let inner = me.inner.get_mut();
        match io::Read::read(&mut inner.buffer, &mut dst.initialize_unfilled()) {
            Ok(z) => {
                info!("Wrote {z} bytes from buffer and inner stream is done? {}", inner.stream_done);
                if z > 0 {
                    return std::task::Poll::Ready(Ok(()))
                } else if inner.stream_done {
                    return std::task::Poll::Ready(Ok(()))
                }
            },
            Err(e) => {
                info!("Buffer read error: {e}");
                return std::task::Poll::Ready(Err(e))
            }
        }
        if inner.stream_done {
            info!("Stream is done after buffer check, we're done!");
            return std::task::Poll::Ready(Ok(()))
        }
        loop {
            let mut state = inner.state.take().unwrap();
            match state {
                State::Idle(ref mut buf_cell) => {
                    log::info!("Idled, destination buffer: {}", dst.filled().len());
                    let mut buf = buf_cell.take().unwrap();

                    if !buf.is_empty() {
                        log::info!("Buffer is not empty! {}", buf.len());
                        inner.buffer = io::Cursor::new(buf.to_vec());
                        buf.clear();
                    }

                    if !inner.stream_done {
                        info!("scheduling next stream read");
                        let f = JsFuture::from(me.stream_reader.read());
                        inner.state = Some(State::Busy(f, buf));
                    } else {
                        info!("Inner state is complete");
                        inner.state = Some(State::Complete);
                    }

                    match io::Read::read(&mut inner.buffer, &mut dst.initialize_unfilled()) {
                        Ok(z) => {
                            info!("Idle: Wrote {z} bytes from buffer and inner stream is done? {}", inner.stream_done);
                            if z > 0 {
                                return std::task::Poll::Ready(Ok(()))
                            } else if inner.stream_done {
                                return std::task::Poll::Ready(Ok(()))
                            }
                        },
                        Err(e) => {
                            info!("Idle: Read buffer with error: {e}");
                            return std::task::Poll::Ready(Err(e))
                        }
                    }

                },
                State::Busy(mut fut, mut buf) => {
                    log::info!("Busy, polling future");
                    let j = Pin::new(&mut fut).poll(cx);

                    match j {
                        Poll::Ready(jres) => {
                            match jres {
                                Ok(res) => {
                                    let value = Reflect::get(&res, &JsValue::from_str("value")).unwrap();
                                    let done = Reflect::get(&res, &JsValue::from_str("done")).unwrap().as_bool().unwrap();
                                    if done || value.is_null() || value.is_undefined() {
                                        log::info!("Exiting, done? {done}, value bad? {}", value.is_null() || value.is_undefined());
                                        inner.state = Some(State::Idle(Some(buf)));
                                        inner.stream_done = true;
                                    } else {
                                        let value_buf = Uint8Array::from(value).to_vec();
                                        log::info!("Received chunk, {}", value_buf.len());
                                        buf.extend_from_slice(&value_buf);
                                        inner.state = Some(State::Idle(Some(buf)));
                                        inner.stream_done = done;
                                    }
                                },
                                Err(err) => {
                                    log::error!("Error occurred while polling stream: {:?}", err.as_string());
                                    inner.state = Some(State::Idle(Some(BytesMut::new())));
                                    return Poll::Ready(Err(io::Error::new(io::ErrorKind::Other, format!("{err:?}"))))
                                },
                            }
                        },
                        Poll::Pending => {
                            log::info!("Not ready yet");
                            inner.state = Some(State::Busy(fut, buf));
                            return Poll::Pending
                        },
                    }
                },
                State::Complete => {
                    inner.state = Some(State::Complete);
                    return std::task::Poll::Ready(Ok(()))
                }
            }
            assert!(inner.state.is_some());
        }
    }
}


#[wasm_bindgen]
pub async fn test_reader(stream_reader: ReadableStreamDefaultReader) -> String {
    let mut handle = WebReaderAsyncRead::new(stream_reader);
    let mut buf = String::new();
    handle.read_to_string(&mut buf).await.unwrap();
    buf
}