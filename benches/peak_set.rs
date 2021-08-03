use criterion::{black_box, criterion_group, criterion_main, Criterion};

use std::fs;

use mzdata::io::MGFReader;
use mzdata::MassErrorType;
use mzpeaks::{PeakCollection, PeakSet};

fn find_with_search(peak_set: &PeakSet, queries: &Vec<f64>) {
    for q in queries.iter() {
        peak_set.search(*q, 10.0, MassErrorType::PPM);
    }
}

fn search(
    peak_set: &PeakSet,
    query: f64,
    error_tolerance: f64,
    error_type: MassErrorType,
) -> Option<usize> {
    let lower_bound = error_type.lower_bound(query, error_tolerance);
    match _search_by(peak_set, lower_bound) {
        Ok(j) => peak_set._closest_peak(query, error_tolerance, j, error_type),
        Err(j) => peak_set._closest_peak(query, error_tolerance, j, error_type),
    }
}

fn _search_by(peak_set: &PeakSet, query: f64) -> Result<usize, usize> {
    peak_set
        .peaks
        .binary_search_by(|peak| peak.mz.partial_cmp(&query).unwrap())
}

fn peak_set_find_inline(peak_set: &PeakSet, queries: &Vec<f64>) {
    for q in queries.iter() {
        search(peak_set, *q, 10.0, MassErrorType::PPM);
    }
}

fn peak_search(c: &mut Criterion) {
    let peaks = MGFReader::new(fs::File::open("./test/data/small.mgf").unwrap())
        .next()
        .unwrap()
        .into_centroid()
        .unwrap()
        .peaks;
    let queries: Vec<f64> = peaks.iter().map(|p| p.mz).collect();
    c.bench_function("peak_set_search", |b| {
        b.iter(|| find_with_search(black_box(&peaks), black_box(&queries)))
    });
    c.bench_function("peak_set_find_inline", |b| {
        b.iter(|| peak_set_find_inline(black_box(&peaks), black_box(&queries)))
    });
}

criterion_group!(benches, peak_search);
criterion_main!(benches);
