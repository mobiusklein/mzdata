use std::collections::HashMap;
use std::env;
use std::io;
use std::path;
use std::process;
use std::thread::spawn;
use std::time;

use std::sync::mpsc::sync_channel;

use mzdata::io::Source;
use mzdata::meta::SpectrumType;
use mzdata::prelude::*;
use mzdata::spectrum::ChromatogramType;
use mzdata::spectrum::{
    DeconvolutedSpectrum, MultiLayerSpectrum, RefPeakDataLevel, SignalContinuity, SpectrumLike,
};
use mzdata::MZReader;

struct MSDataFileSummary {
    pub start_time: f64,
    pub end_time: f64,
    pub level_table: HashMap<u8, usize>,
    pub charge_table: HashMap<i32, usize>,
    pub peak_charge_table: HashMap<u8, HashMap<i32, usize>>,
    pub peak_mode_table: HashMap<SignalContinuity, usize>,
    pub has_ion_mobility: bool,
    pub spectrum_types: HashMap<Option<SpectrumType>, usize>,
    pub chromatogram_types: HashMap<ChromatogramType, usize>
}

impl Default for MSDataFileSummary {
    fn default() -> Self {
        Self {
            start_time: f64::INFINITY,
            end_time: f64::NEG_INFINITY,
            level_table: Default::default(),
            charge_table: Default::default(),
            peak_charge_table: Default::default(),
            peak_mode_table: Default::default(),
            has_ion_mobility: false,
            spectrum_types: Default::default(),
            chromatogram_types: Default::default(),
        }
    }
}

impl MSDataFileSummary {
    pub fn handle_scan(&mut self, scan: MultiLayerSpectrum) {
        let time = scan.start_time();
        self.start_time = self.start_time.min(time);
        self.end_time = self.end_time.max(time);
        let level = scan.ms_level();
        *self.level_table.entry(level).or_default() += 1;
        *self.spectrum_types.entry(scan.spectrum_type()).or_default() += 1;
        if level > 1 {
            if let Some(charge) = scan.precursor().unwrap().ion().and_then(|i| i.charge()) {
                *self.charge_table.entry(charge).or_default() += 1;
            } else {
                *self.charge_table.entry(0).or_default() += 1;
            }
        }
        *self
            .peak_mode_table
            .entry(scan.signal_continuity())
            .or_default() += scan.peaks().len();

        let has_charge = match scan.peaks() {
            RefPeakDataLevel::Missing => false,
            RefPeakDataLevel::RawData(arrays) => arrays.charges().is_ok(),
            RefPeakDataLevel::Centroid(_) => false,
            RefPeakDataLevel::Deconvoluted(_) => true,
        };

        let has_ion_mobility = match scan.peaks() {
            RefPeakDataLevel::RawData(arrays) => arrays.has_ion_mobility(),
            _ => false,
        } || scan.has_ion_mobility();
        self.has_ion_mobility |= has_ion_mobility;

        if has_charge {
            let deconv_scan: DeconvolutedSpectrum = scan.try_into().unwrap();
            deconv_scan.deconvoluted_peaks.iter().for_each(|p| {
                *(*self
                    .peak_charge_table
                    .entry(deconv_scan.ms_level())
                    .or_default())
                .entry(p.charge)
                .or_default() += 1;
                assert!((p.index as usize) < deconv_scan.deconvoluted_peaks.len())
            })
        }
    }

    pub fn _scan_file<R: SpectrumSource>(&mut self, reader: &mut R) {
        let start = time::Instant::now();
        reader.enumerate().for_each(|(i, scan)| {
            if i % 10000 == 0 && i > 0 {
                println!(
                    "\tScan {}: {} ({:0.3} seconds, {} peaks|points)",
                    i,
                    scan.id(),
                    (time::Instant::now() - start).as_secs_f64(),
                    self.peak_mode_table.values().sum::<usize>()
                );
            }
            self.handle_scan(scan);
        });
        let end = time::Instant::now();
        let elapsed = end - start;
        println!("{:0.3} seconds elapsed", elapsed.as_secs_f64());
    }

    pub fn scan_file<R: SpectrumSource + Send + 'static>(&mut self, reader: R) {
        self.scan_file_threaded(reader)
    }

    pub fn scan_file_threaded<R: SpectrumSource + Send + 'static>(&mut self, reader: R) {
        let start = time::Instant::now();
        let (sender, receiver) = sync_channel(2usize.pow(12));
        let read_handle = spawn(move || {
            reader.into_iter()
                .enumerate()
                .for_each(|(i, scan)| {
                    sender.send((i, scan)).unwrap()
                });
        });
        let i = receiver.iter().fold(0, |_, (i, scan)| {
            if i % 10000 == 0 && i > 0 {
                println!(
                    "\tScan {}: {} ({:0.3} seconds, {} peaks|points)",
                    i,
                    scan.id(),
                    (time::Instant::now() - start).as_secs_f64(),
                    self.peak_mode_table.values().sum::<usize>()
                );
            }
            self.handle_scan(scan);
            i
        });
        read_handle.join().unwrap();
        let end = time::Instant::now();
        let elapsed = end - start;
        println!("{:0.3} seconds elapsed, handled {i} spectra", elapsed.as_secs_f64());
    }

    pub fn scan_chromatograms<R: ChromatogramSource + Send + 'static>(&mut self, reader: &mut R) {
        for chrom in reader.iter_chromatograms() {
            *self.chromatogram_types.entry(chrom.chromatogram_type()).or_default() += 1;
        }
    }

    pub fn write_out(&self) {
        println!("Start Time: {:0.2}", self.start_time);
        println!("End Time: {:0.2}", self.end_time);
        println!("Has Ion Mobility: {}", self.has_ion_mobility);
        println!("MS Levels:");
        let mut level_set: Vec<(&u8, &usize)> = self.level_table.iter().collect();
        level_set.sort_by_key(|(a, _)| *a);
        for (level, count) in level_set.iter() {
            println!("\t{}: {}", level, count);
        }

        let mut spec_types: Vec<(_, _)> = self.spectrum_types.iter().collect();
        spec_types.sort_by(|a, b| a.0.cmp(b.0));

        if spec_types.iter().any(|v| v.0.is_some_and(|v| v != SpectrumType::MS1Spectrum && v != SpectrumType::MSnSpectrum)) {
            println!("Spectrum types:");
            for (level, count) in spec_types.iter() {
                if let Some(level) = level {
                    println!("\t{}: {}", level.name(), count);
                } else {
                    println!("\tUnknown: {}", count);
                }
            }
        }

        println!("Precursor Charge States:");
        let mut charge_set: Vec<(&i32, &usize)> = self.charge_table.iter().collect();
        charge_set.sort_by_key(|(a, _)| *a);
        for (charge, count) in charge_set.iter() {
            if **charge == 0 {
                println!("\tCharge Not Reported: {}", count);
            } else {
                println!("\t{}: {}", charge, count);
            }
        }

        let mut peak_charge_levels: Vec<_> = self.peak_charge_table.iter().collect();

        peak_charge_levels.sort_by(|(level_a, _), (level_b, _)| level_a.cmp(level_b));

        for (level, peak_charge_table) in peak_charge_levels {
            if !peak_charge_table.is_empty() {
                println!("Peak Charge States for MS level {}:", level);
                let mut peak_charge_set: Vec<(&i32, &usize)> = peak_charge_table.iter().collect();
                peak_charge_set.sort_by_key(|(a, _)| *a);
                for (charge, count) in peak_charge_set.iter() {
                    if **charge == 0 {
                        println!("\tCharge Not Reported: {}", count);
                    } else {
                        println!("\t{}: {}", charge, count);
                    }
                }
            }
        }
        self.peak_mode_table
            .iter()
            .for_each(|(mode, count)| match mode {
                SignalContinuity::Unknown => println!("Unknown continuity: {}", count),
                SignalContinuity::Centroid => println!("Peaks: {}", count),
                SignalContinuity::Profile => println!("Points: {}", count),
            });

        if !self.chromatogram_types.is_empty() {
            let mut spec_types: Vec<(_, _)> = self.chromatogram_types.iter().collect();
            spec_types.sort_by(|a, b| a.0.cmp(b.0));
            println!("Chromatogram types:");
            for (level, count) in spec_types.iter() {
                println!("\t{:?}: {}", level, count);
            }
        }
    }
}

fn main() -> io::Result<()> {
    env_logger::init();
    let path = path::PathBuf::from(env::args().nth(1).unwrap_or_else(|| {
        eprintln!("Please provide a path to an MS data file");
        process::exit(1)
    }));
    let mut summarizer = MSDataFileSummary::default();

    if path.as_os_str() == "-" {
        mzdata::mz_read!(Source::Stdin, reader => {
            summarizer.scan_file(reader)
        })?;
    } else {
        let mut reader = MZReader::open_path(path)?;
        eprintln!("Format: {}", reader.as_format());
        if reader.count_chromatograms() > 0 {
            summarizer.scan_chromatograms(&mut reader);
        }
        summarizer.scan_file(reader)
    };

    summarizer.write_out();
    Ok(())
}
