pub use mzdata_spectrum::{
    CentroidPeakAdapting, CentroidSpectrum, CentroidSpectrumType, DeconvolutedPeakAdapting,
    DeconvolutedSpectrum, DeconvolutedSpectrumType, MultiLayerSpectrum, RawSpectrum, Spectrum,
    SpectrumConversionError, SpectrumLike, SpectrumProcessingError,
};

#[cfg(all(test, feature = "mzml"))]
mod test {
    use std::collections::HashMap;
    use std::{fs, io};

    use crate::io::mzml::MzMLReader;
    use crate::io::DetailLevel;
    use crate::prelude::*;
    use crate::spectrum::BinaryArrayMap3D;
    use crate::spectrum::*;

    #[test_log::test]
    fn test_peakdata_lazy() -> io::Result<()> {
        let mut reader = MzMLReader::open_path("./test/data/small.mzML")?;
        reader.detail_level = DetailLevel::Lazy;
        let spec = reader.get_spectrum_by_index(0).unwrap();

        let peaks = spec.peaks();
        let n = peaks.len();
        assert_eq!(n, 19913);

        let iter = peaks.iter();
        let data: Vec<_> = iter.collect();

        let n2 = data.len();
        assert_eq!(n2, 19913);

        let p1 = peaks.get(5000);
        let p2 = data.get(5000).cloned();
        assert_eq!(p1, p2);
        Ok(())
    }

    macro_rules! behaviors {
        ($spec:ident) => {
            assert_eq!($spec.id(), "controllerType=0 controllerNumber=1 scan=10014");
            assert_eq!($spec.start_time(), 22.12829);
            assert_eq!($spec.ms_level(), 1);
            assert_eq!($spec.polarity(), ScanPolarity::Positive);
            assert!($spec.precursor().is_none());
            assert_eq!($spec.precursor_iter().count(), 0);
            assert!(!$spec.has_ion_mobility());
            assert!(!$spec.has_ion_mobility_dimension());
            assert!(!$spec.params().is_empty());
        };
    }

    #[allow(unused)]
    fn test_spectrum_behavior<T: SpectrumLike>(spec: &T) {
        assert_eq!(
            spec.spectrum_type(),
            Some(crate::meta::SpectrumType::MS1Spectrum)
        );
        behaviors!(spec);
    }

    #[test_log::test]
    fn test_peakdata() -> io::Result<()> {
        let mut reader = MzMLReader::open_path("./test/data/small.mzML")?;
        let spec = reader.get_spectrum_by_index(0).unwrap();

        let peaks = spec.peaks();
        let n = peaks.len();
        assert_eq!(n, 19913);

        let iter = peaks.iter();
        let data: Vec<_> = iter.collect();

        let n2 = data.len();
        assert_eq!(n2, 19913);

        let p1 = peaks.get(5000);
        let p2 = data.get(5000).cloned();
        assert_eq!(p1, p2);
        Ok(())
    }

    #[cfg(feature = "mzsignal")]
    #[test_log::test]
    fn test_profile_read() {
        let mut reader = MzMLReader::open_path("./test/data/three_test_scans.mzML")
            .expect("Failed to open test file");
        reader.reset();
        let mut scan = reader.next().expect("Failed to read spectrum");
        assert_eq!(scan.signal_continuity(), SignalContinuity::Profile);
        test_spectrum_behavior(&scan);
        if let Err(err) = scan.pick_peaks(1.0) {
            panic!("Should not have an error! {}", err);
        }

        if let Some(peaks) = &scan.peaks {
            assert_eq!(peaks.len(), 2107);
            let n = peaks.len();
            for i in 1..n {
                let diff = peaks[i].get_index() - peaks[i - 1].get_index();
                assert_eq!(diff, 1);
            }

            peaks
                .has_peak(562.741, Tolerance::PPM(3f64))
                .expect("Expected to find peak");
            peaks
                .has_peak(563.240, Tolerance::PPM(3f64))
                .expect("Expected to find peak");
            let p = peaks
                .has_peak(563.739, Tolerance::PPM(1f64))
                .expect("Expected to find peak");
            assert!((p.mz() - 563.739).abs() < 1e-3)
        }

        scan.pick_peaks_with(&mzsignal::PeakPicker::default())
            .unwrap();
        if let Some(peaks) = &scan.peaks {
            assert_eq!(peaks.len(), 2107);
            let n = peaks.len();
            for i in 1..n {
                let diff = peaks[i].get_index() - peaks[i - 1].get_index();
                assert_eq!(diff, 1);
            }

            peaks
                .has_peak(562.741, Tolerance::PPM(3f64))
                .expect("Expected to find peak");
            peaks
                .has_peak(563.240, Tolerance::PPM(3f64))
                .expect("Expected to find peak");
            let p = peaks
                .has_peak(563.739, Tolerance::PPM(1f64))
                .expect("Expected to find peak");
            assert!((p.mz() - 563.739).abs() < 1e-3)
        }

        let mut cent_scan = scan.clone().into_centroid().unwrap();
        assert_eq!(cent_scan.signal_continuity(), SignalContinuity::Centroid);
        assert_eq!(SpectrumLike::index(&cent_scan), 0);
        behaviors!(cent_scan);
        cent_scan.update_summaries();

        let tmp: Spectrum = cent_scan.into_spectrum().unwrap();
        let raw_scan = tmp.into_raw().unwrap();
        behaviors!(raw_scan);
        assert_eq!(raw_scan.index(), 0);
        assert_eq!(raw_scan.signal_continuity(), SignalContinuity::Centroid);
    }

    #[cfg(feature = "mzsignal")]
    #[test]
    fn test_reprofile() {
        let mut reader = MzMLReader::open_path("./test/data/three_test_scans.mzML")
            .expect("Failed to open test file");
        reader.reset();
        let mut scan = reader.next().expect("Failed to read spectrum");
        assert_eq!(scan.signal_continuity(), SignalContinuity::Profile);
        assert_eq!(scan.ms_level(), 1);
        assert_eq!(scan.polarity(), ScanPolarity::Positive);
        assert!(scan.precursor().is_none());

        if let Err(err) = scan.pick_peaks(1.0) {
            panic!("Should not have an error! {}", err);
        }

        let mut duplicate = scan.clone();
        duplicate.reprofile_with_shape(0.001, 0.01).unwrap();
        duplicate.peaks = None;
        duplicate.pick_peaks(1.0).unwrap();
        let peak = duplicate.peaks.as_ref().unwrap().base_peak().unwrap();
        eprintln!("{}", peak);
    }

    #[test]
    fn test_3d_stack_unstack() -> io::Result<()> {
        let mut reader = crate::MZReader::open_gzipped_read(io::BufReader::new(fs::File::open(
            "test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz",
        )?))?;

        let spec = reader
            .get_spectrum_by_id("merged=42869 frame=9717 scanStart=1 scanEnd=705")
            .unwrap();
        let mut arrays = spec.arrays.unwrap();
        let units_map: HashMap<_, _> = arrays.iter().map(|(k, v)| (k.clone(), v.unit)).collect();
        let mzs = arrays.mzs()?;
        assert!(!mzs.is_sorted());
        drop(mzs);
        arrays.sort_by_array(&ArrayType::MZArray)?;
        let mzs = arrays.mzs()?;
        assert!(mzs.is_sorted());
        let n = mzs.len();

        let arrays_3d = BinaryArrayMap3D::stack(&arrays)?;
        let stacked_n: usize = arrays_3d
            .iter()
            .map(|(_, va)| va.mzs().unwrap().len())
            .sum();

        assert_eq!(
            units_map[&arrays_3d.ion_mobility_type],
            arrays_3d.ion_mobility_unit
        );

        assert_eq!(n, stacked_n);

        let unstacked = arrays_3d.unstack()?;
        let unstacked_n = unstacked.mzs()?.len();

        assert_eq!(unstacked_n, n);
        for (k, v) in unstacked.iter() {
            assert_eq!(units_map[k], v.unit);
        }
        Ok(())
    }
}
