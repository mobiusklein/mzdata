use std::{
    collections::HashMap,
    sync::mpsc::{Receiver, Sender, TryRecvError},
    time::Duration,
};

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
    pub fn receive(&mut self, group_idx: usize, group: T) {
        self.waiting.insert(group_idx, group);
    }

    pub fn receive_from(&mut self, receiver: &Receiver<(usize, T)>, batch_size: usize) {
        let mut counter = 0usize;
        while let Ok((group_idx, group)) = receiver.recv_timeout(Duration::from_micros(1)) {
            self.receive(group_idx, group);
            counter += 1;
            if counter > batch_size {
                break;
            }
        }
    }

    pub fn has_next(&self) -> bool {
        self.waiting.contains_key(&self.next_key)
    }

    pub fn try_next(&mut self) -> Option<(usize, T)> {
        self.waiting.remove_entry(&self.next_key).and_then(|op| {
            self.next_key += 1;
            Some(op)
        })
    }

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
