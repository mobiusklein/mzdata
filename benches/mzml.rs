use criterion::{black_box, criterion_group, criterion_main, Criterion};

use std::fs;

use mzdata::io::mzml::MzMLReader;
use mzdata::io::prelude::*;

fn serial(file_path: &str) {
    let file = fs::File::open(file_path).unwrap();
    let reader = MzMLReader::new(file);
    let total: usize = reader.map(|s| s.mzs().len()).sum();
    assert_eq!(total, 305214);
}

fn serial_with_index(file_path: &str) {
    let file = fs::File::open(file_path).unwrap();
    let reader = MzMLReader::new_indexed(file);
    let total: usize = reader.map(|s| s.mzs().len()).sum();
    assert_eq!(total, 305214);
}

fn serial_with_external_iterator(file_path: &str) {
    let file = fs::File::open(file_path).unwrap();
    let mut reader = MzMLReader::new_indexed(file);
    let total: usize = reader.iter().map(|s| s.mzs().len()).sum();
    assert_eq!(total, 305214);
}

fn mzml_totaling(c: &mut Criterion) {
    c.bench_function("serial_execution", |b| {
        b.iter(|| serial(black_box("./test/data/small.mzML")))
    });
    c.bench_function("serial_execution_with_index", |b| {
        b.iter(|| serial_with_index(black_box("./test/data/small.mzML")))
    });
    c.bench_function("serial_with_external_iterator", |b| {
        b.iter(|| serial_with_external_iterator(black_box("./test/data/small.mzML")))
    });
}

criterion_group!(benches, mzml_totaling);
criterion_main!(benches);
