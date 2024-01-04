use std::cmp::Ordering;
use std::collections::VecDeque;
use std::io;
use std::ops::Range;
use std::{env, fs, path};

use mzdata::prelude::*;
use mzdata::spectrum::{MultiLayerSpectrum, SpectrumGroup, SpectrumGroupingIterator};
use mzpeaks::{CentroidPeak, DeconvolutedPeak, Tolerance};

#[derive(Default, Debug, Clone, PartialEq)]
pub struct SelectionTargetSpecification {
    pub mz: f64,
    pub charge: Option<i32>,
    pub time_range: Range<f64>,
}

impl SelectionTargetSpecification {
    pub fn new(mz: f64, charge: Option<i32>, time_range: Range<f64>) -> Self {
        Self {
            mz,
            charge,
            time_range,
        }
    }

    pub fn spans(&self, time: f64) -> bool {
        self.time_range.start <= time && self.time_range.end > time
    }

    pub fn from_spectrum(spectrum: &MultiLayerSpectrum) -> Self {
        let prec = spectrum.precursor().unwrap();
        Self {
            mz: prec.mz(),
            charge: prec.charge(),
            time_range: (spectrum.start_time()..spectrum.start_time()),
        }
    }

    pub fn observe(&mut self, spectrum: &MultiLayerSpectrum) {
        let t = spectrum.start_time();
        if t < self.time_range.start {
            self.time_range.start = t;
        }
        if t > self.time_range.end {
            self.time_range.end = t;
        }
    }
}

impl PartialOrd for SelectionTargetSpecification {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.mz.partial_cmp(&other.mz) {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }
        match self.charge.partial_cmp(&other.charge) {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }
        match self.time_range.start.partial_cmp(&other.time_range.start) {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }
        self.time_range.end.partial_cmp(&other.time_range.end)
    }
}

#[derive(Default, Debug, Clone, Copy, PartialEq)]
pub struct SelectedTarget {
    pub mz: f64,
    pub charge: Option<i32>,
}

impl SelectedTarget {
    pub fn new(mz: f64, charge: Option<i32>) -> Self {
        Self { mz, charge }
    }
}


pub struct MSnTargetTrackingIterator<R: ScanSource> {
    source: SpectrumGroupingIterator<R, CentroidPeak, DeconvolutedPeak, MultiLayerSpectrum>,
    time_width: f64,
    error_tolerance: Tolerance,
    buffer: VecDeque<(SpectrumGroup, f64)>,
    pushback_buffer: Option<(SpectrumGroup, f64)>,
    targets: VecDeque<SelectionTargetSpecification>,
}

impl<R: ScanSource> MSnTargetTrackingIterator<R> {
    pub fn new(
        source: SpectrumGroupingIterator<R, CentroidPeak, DeconvolutedPeak, MultiLayerSpectrum>,
        time_width: f64,
        error_tolerance: Tolerance,
    ) -> Self {
        let mut inst = Self {
            source,
            time_width,
            error_tolerance,
            buffer: Default::default(),
            pushback_buffer: Default::default(),
            targets: Default::default(),
        };
        inst.initial_feed();
        inst
    }

    fn observe(&mut self, group: &SpectrumGroup) {
        let error_tolerance = self.error_tolerance;
        let time_width_half = self.time_width / 2.0;
        let new_targets: usize = group
            .products
            .iter()
            .map(|s| {
                let prec = s.precursor().unwrap();
                let mz = prec.mz();
                let t = s.start_time();
                let hits: usize = self
                    .targets
                    .iter_mut()
                    .map(|p| {
                        if error_tolerance.test(p.mz, mz) && (t > p.time_range.end) && (t - p.time_range.end) < time_width_half {
                            p.time_range.end = t;
                            1
                        } else {
                            0
                        }
                    })
                    .sum();
                if hits == 0 {
                    let p = SelectionTargetSpecification::new(mz, prec.charge(), t-0.5..t+0.5);
                    // eprintln!("Added {p:?}");
                    self.targets.push_back(p);
                    1
                } else {
                    0
                }
            })
            .sum();
        if new_targets > 0 {
            log::debug!("Added {new_targets} new targets");
        }
    }

    fn initial_feed(&mut self) {
        let group = self.source.next();
        let start_time = group
            .as_ref()
            .and_then(|s| s.earliest_spectrum().map(|s| s.start_time()))
            .unwrap();
        let end_time = start_time + self.time_width;
        eprintln!("Initial time window {start_time} to {end_time}");
        if let Some(g) = group {
            self.observe(&g);
            self.buffer.push_back((g, start_time));
        }
        while let Some(group) = self.source.next() {
            let t = group.earliest_spectrum().map(|s| s.start_time()).unwrap();
            if t < end_time {
                self.observe(&group);
                self.buffer.push_back((group, t));
            } else {
                self.pushback_buffer = Some((group, t));
                break;
            }
        }
        eprintln!("{} targets extracted from buffer size {}", self.targets.len(), self.buffer.len());
    }


    fn get_current_window_end(&self) -> f64 {
        if self.buffer.len() == 0 {
            return f64::NEG_INFINITY;
        }
        let start = self.buffer.front().and_then(|(_, t)| Some(*t)).unwrap();
        let end = self.buffer.back().and_then(|(_, t)| Some(*t)).unwrap();
        let mid = start + (end - start) / 2.0;
        mid + self.time_width
    }

    fn feed_next(&mut self) {
        let threshold = self.get_current_window_end();
        let (use_pushback, pushback_populated) = if let Some((_, t)) = self.pushback_buffer.as_ref()
        {
            (*t < threshold, true)
        } else {
            (false, false)
        };
        if use_pushback {
            let (group, t) = self.pushback_buffer.take().unwrap();
            self.observe(&group);
            self.buffer.push_back((group, t));
        }
        if !pushback_populated {
            if let Some(group) = self.source.next() {
                let t = group.earliest_spectrum().map(|s| s.start_time()).unwrap();
                if t < threshold {
                    self.observe(&group);
                    self.buffer.push_back((group, t));
                } else {
                    self.pushback_buffer = Some((group, t));
                }
            }
        }
    }

    fn step(&mut self) -> Option<(SpectrumGroup, Vec<SelectedTarget>)> {
        if let Some((group, t)) = self.buffer.pop_front() {
            let targets: Vec<_> = self
                .targets
                .iter()
                .filter(|p| p.spans(t))
                .map(|p| SelectedTarget::new(p.mz, p.charge))
                .collect();
            self.targets = self
                .targets
                .drain(..)
                .filter(|target| {
                    // Keep targets which have not ended by this time point
                    let cond = target.time_range.end >= t;
                    if !cond {
                        log::debug!("Dropping {target:?} at {t}")
                    }
                    cond
                })
                .collect();
            Some((group, targets))
        } else {
            None
        }
    }
}

impl<R: ScanSource> Iterator for MSnTargetTrackingIterator<R> {
    type Item = (SpectrumGroup, Vec<SelectedTarget>);

    fn next(&mut self) -> Option<Self::Item> {
        let value = self.step();
        self.feed_next();
        value
    }
}

fn main() -> io::Result<()> {
    let path = path::PathBuf::from(
        env::args()
            .nth(1)
            .expect("Please pass an MS data file path"),
    );
    env_logger::init();

    let reader = mzdata::MzMLReader::new(fs::File::open(path)?);
    let group_iter = reader.into_groups();
    let tracking_iter = MSnTargetTrackingIterator::new(group_iter, 2.0, Tolerance::PPM(5.0));
    for (group, targets) in tracking_iter {
        if group.precursor.is_some() {
            let prec = group.precursor.as_ref().unwrap();
            eprintln!("{} {}: {}", prec.id(), prec.start_time(), targets.len());
        }
    }
    Ok(())
}
