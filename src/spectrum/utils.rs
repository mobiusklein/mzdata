use std::{
    collections::HashMap,
    sync::mpsc::{Receiver, Sender, TryRecvError, SyncSender},
    time::Duration,
};


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
    pub fn receive_from_timeout(&mut self, receiver: &Receiver<(usize, T)>, batch_size: usize, timeout: Duration) {
        let mut counter = 0usize;
        while let Ok((group_idx, group)) = receiver.recv_timeout(timeout) {
            self.receive(group_idx, group);
            counter += 1;
            if counter > batch_size {
                break;
            }
        }
    }

    /// Check if the next item is already available
    pub fn has_next(&self) -> bool {
        self.waiting.contains_key(&self.next_key)
    }

    /// Try to get the next item according to the collation, or None if it hasn't arrived yet
    pub fn try_next(&mut self) -> Option<(usize, T)> {
        self.waiting.remove_entry(&self.next_key).and_then(|op| {
            self.next_key += 1;
            Some(op)
        })
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