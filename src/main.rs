use structured::{Peak, PeakSet, MassErrorType, PeakCollection};


fn main() {
    let peak = Peak {mz: 100.0, .. Default::default()};
    let peak2 = Peak {mz: 101.0, .. Default::default()};
    let peak3 = Peak {mz: 102.0, .. Default::default()};
    println!("{}, {}", peak, peak2);
    println!("{}", peak <= peak2);

    let peaks = PeakSet::new(vec![peak3, peak2, peak]);
    println!("{}", peaks);
    println!("{}, {}", peaks[0], peaks[2]);
    match peaks.has_peak(101.01, 0.02, MassErrorType::Exact) {
        Some(p) => println!("{}", p),
        None => println!("No Match")
    }
    for (i, p) in (&peaks).into_iter().enumerate() {
        println!("{} => {}", i, p);
    }
    for (i, p) in (&peaks).into_iter().enumerate() {
        println!("{} => {}", i, p);
    }
}
