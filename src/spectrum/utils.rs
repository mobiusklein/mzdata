use std::{
    collections::HashMap,
    io,
    sync::mpsc::{Receiver, Sender, SyncSender, TryRecvError},
    time::Duration,
};

use mzpeaks::{CentroidLike, DeconvolutedCentroidLike};

use crate::prelude::*;

use super::MultiLayerSpectrum;

/// A helper for consuming parallel iteration in the original ordering sequentially later.
/// Useful for things like splitting work up with `rayon` and then merging it back together
/// for writing out at the end.
///
/// **Warning**: In the worst-case scenario, this may lead to buffering the entire data stream
/// in memory waiting for a slow-to-arrive item. No effort is made to mitigate this, but this is
/// expected to be quite irregular.
#[derive(Debug)]
pub struct Collator<T: Send> {
    pub waiting: HashMap<usize, T>,
    pub next_key: usize,
    pub ticks: usize,
    pub done: bool,
}

impl<T: Send> Default for Collator<T> {
    fn default() -> Self {
        Self {
            waiting: Default::default(),
            next_key: Default::default(),
            ticks: Default::default(),
            done: Default::default(),
        }
    }
}

impl<T: Send> Collator<T> {
    /// Take the next item `group` with ordering key `group_idx` and add it to the waiting
    /// item queue
    pub fn receive(&mut self, group_idx: usize, group: T) {
        self.waiting.insert(group_idx, group);
    }

    /// Block on reading from `receiver` until it takes more than a microsecond to
    /// retrieve the next item from it, or until `batch_size` items have been read
    pub fn receive_from(&mut self, receiver: &Receiver<(usize, T)>, batch_size: usize) {
        self.receive_from_timeout(receiver, batch_size, Duration::from_micros(1))
    }

    /// Block on reading from `receiver` until it takes more than `timeout` time
    /// to retrieve the next item from it, or until `batch_size` items have been read
    pub fn receive_from_timeout(
        &mut self,
        receiver: &Receiver<(usize, T)>,
        batch_size: usize,
        timeout: Duration,
    ) {
        let mut counter = 0usize;
        while let Ok((group_idx, group)) = receiver.recv_timeout(timeout) {
            self.receive(group_idx, group);
            counter += 1;
            if counter > batch_size {
                break;
            }
        }
    }

    pub fn receive_from_map_timeout<U, F: Fn(usize, U) -> (usize, T)>(
        &mut self,
        receiver: &Receiver<(usize, U)>,
        batch_size: usize,
        timeout: Duration,
        cb: F,
    ) {
        let mut counter = 0usize;
        while let Ok((group_idx, group)) = receiver.recv_timeout(timeout) {
            let (group_idx, group) = cb(group_idx, group);
            self.receive(group_idx, group);
            counter += 1;
            if counter > batch_size {
                break;
            }
        }
    }

    pub fn receive_from_map_iter_timeout<
        U,
        I: Iterator<Item = (usize, T)>,
        F: Fn(usize, U) -> I,
    >(
        &mut self,
        receiver: &Receiver<(usize, U)>,
        batch_size: usize,
        timeout: Duration,
        cb: F,
    ) {
        let mut counter = 0usize;
        while let Ok((group_idx, group)) = receiver.recv_timeout(timeout) {
            self.receive_map_iter(group_idx, group, &cb);
            counter += 1;
            if counter > batch_size {
                break;
            }
        }
    }

    pub fn receive_map<U, F: Fn(usize, U) -> (usize, T)>(
        &mut self,
        group_idx: usize,
        group: U,
        cb: F,
    ) {
        let (group_idx, group) = cb(group_idx, group);
        self.receive(group_idx, group);
    }

    pub fn receive_map_iter<U, I: Iterator<Item = (usize, T)>, F: Fn(usize, U) -> I>(
        &mut self,
        group_idx: usize,
        group: U,
        cb: F,
    ) {
        cb(group_idx, group).for_each(|(i, x)| {
            self.receive(i, x);
        })
    }

    /// Check if the next item is already available
    pub fn has_next(&self) -> bool {
        self.waiting.contains_key(&self.next_key)
    }

    /// Try to get the next item according to the collation, or None if it hasn't arrived yet
    pub fn try_next(&mut self) -> Option<(usize, T)> {
        self.waiting.remove_entry(&self.next_key).map(|op| {
            self.next_key += 1;
            op
        })
    }

    /// Explicitly set the next key waiting key
    pub fn set_next_key(&mut self, key: usize) {
        self.next_key = key
    }

    /// Given a channel `receiver` to read from and a `sender` to send to another channel
    /// in collation order, while there is still something to receive, repeatedly read from the `receiver`
    /// and write as many items to the `sender` in collation order every cycle.
    ///
    /// **Note**: This function blocks, so it should be run in a separate thread.
    pub fn collate_sync(receiver: Receiver<(usize, T)>, sender: SyncSender<(usize, T)>) {
        let mut collator = Self::default();
        loop {
            match receiver.try_recv() {
                Ok((group_idx, group)) => {
                    collator.receive(group_idx, group);
                    collator.receive_from(&receiver, 100);
                }
                Err(e) => match e {
                    TryRecvError::Empty => {}
                    TryRecvError::Disconnected => {
                        collator.done = true;
                        break;
                    }
                },
            }

            while let Some((group_idx, group)) = collator.try_next() {
                match sender.send((group_idx, group)) {
                    Ok(()) => {}
                    Err(e) => {
                        log::error!("Failed to send {group_idx} for writing: {e}")
                    }
                }
            }
        }
    }

    /// As [`collate_sync`](Collator::collate_sync), but with an unbounded channel
    pub fn collate(receiver: Receiver<(usize, T)>, sender: Sender<(usize, T)>) {
        let mut collator = Self::default();
        loop {
            match receiver.try_recv() {
                Ok((group_idx, group)) => {
                    collator.receive(group_idx, group);
                    collator.receive_from(&receiver, 100);
                }
                Err(e) => match e {
                    TryRecvError::Empty => {}
                    TryRecvError::Disconnected => {
                        collator.done = true;
                        break;
                    }
                },
            }

            while let Some((group_idx, group)) = collator.try_next() {
                match sender.send((group_idx, group)) {
                    Ok(()) => {}
                    Err(e) => {
                        log::error!("Failed to send {group_idx} for writing: {e}")
                    }
                }
            }
        }
    }
}

impl<
        C: CentroidLike + Default + Send + BuildArrayMapFrom + BuildFromArrayMap + Clone,
        D: DeconvolutedCentroidLike + Default + Send + BuildArrayMapFrom + BuildFromArrayMap + Clone,
    > ScanWriter<C, D> for Collator<MultiLayerSpectrum<C, D>>
{
    fn write<S: SpectrumLike<C, D> + 'static>(&mut self, spectrum: &S) -> std::io::Result<usize> {
        let k = spectrum.index();
        let peaks = spectrum.peaks().cloned();
        let descr = spectrum.description().clone();
        let t = MultiLayerSpectrum::from_peaks_data_levels_and_description(peaks, descr);
        self.receive(k, t);
        Ok(k)
    }

    fn write_owned<S: SpectrumLike<C, D> + 'static>(
        &mut self,
        spectrum: S,
    ) -> std::io::Result<usize> {
        let k = spectrum.index();
        let (peaks, description) = spectrum.into_peaks_and_description();
        let t = MultiLayerSpectrum::from_peaks_data_levels_and_description(peaks, description);
        self.receive(k, t);
        Ok(k)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }

    fn close(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl<
        C: CentroidLike + Default + Send + BuildArrayMapFrom + BuildFromArrayMap + Clone,
        D: DeconvolutedCentroidLike + Default + Send + BuildArrayMapFrom + BuildFromArrayMap + Clone,
    > ScanWriter<C, D> for Sender<MultiLayerSpectrum<C, D>> {
    fn write<S: SpectrumLike<C, D> + 'static>(&mut self, spectrum: &S) -> std::io::Result<usize> {
        let k = spectrum.index();
        let peaks = spectrum.peaks().cloned();
        let descr = spectrum.description().clone();
        let t = MultiLayerSpectrum::from_peaks_data_levels_and_description(peaks, descr);
        match self.send(t) {
            Ok(_) => {Ok(k)},
            Err(e) => {
                Err(
                    io::Error::new(io::ErrorKind::BrokenPipe, e.to_string())
                )
            },
        }
    }

    fn write_owned<S: SpectrumLike<C, D> + 'static>(
        &mut self,
        spectrum: S,
    ) -> std::io::Result<usize> {
        let k = spectrum.index();
        let (peaks, description) = spectrum.into_peaks_and_description();
        let t = MultiLayerSpectrum::from_peaks_data_levels_and_description(peaks, description);
        match self.send(t) {
            Ok(_) => {Ok(k)},
            Err(e) => {
                Err(
                    io::Error::new(io::ErrorKind::BrokenPipe, e.to_string())
                )
            },
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }

    fn close(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}



impl<
        C: CentroidLike + Default + Send + BuildArrayMapFrom + BuildFromArrayMap + Clone,
        D: DeconvolutedCentroidLike + Default + Send + BuildArrayMapFrom + BuildFromArrayMap + Clone,
    > ScanWriter<C, D> for SyncSender<MultiLayerSpectrum<C, D>> {
    fn write<S: SpectrumLike<C, D> + 'static>(&mut self, spectrum: &S) -> std::io::Result<usize> {
        let k = spectrum.index();
        let peaks = spectrum.peaks().cloned();
        let descr = spectrum.description().clone();
        let t = MultiLayerSpectrum::from_peaks_data_levels_and_description(peaks, descr);
        match self.send(t) {
            Ok(_) => {Ok(k)},
            Err(e) => {
                Err(
                    io::Error::new(io::ErrorKind::BrokenPipe, e.to_string())
                )
            },
        }
    }

    fn write_owned<S: SpectrumLike<C, D> + 'static>(
        &mut self,
        spectrum: S,
    ) -> std::io::Result<usize> {
        let k = spectrum.index();
        let (peaks, description) = spectrum.into_peaks_and_description();
        let t = MultiLayerSpectrum::from_peaks_data_levels_and_description(peaks, description);
        match self.send(t) {
            Ok(_) => {Ok(k)},
            Err(e) => {
                Err(
                    io::Error::new(io::ErrorKind::BrokenPipe, e.to_string())
                )
            },
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }

    fn close(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}