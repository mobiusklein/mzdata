pub use mzdata_spectrum::{
    MultiLayerIonMobilityFrame, IonMobilityFrameDescription, IonMobilityFrameLike,
    RefFeatureDataLevel, FeatureDataLevel,
};

#[cfg(test)]
mod test {
    use std::{fs, io};

    use super::*;

    use crate::io::Generic3DIonMobilityFrameSource;
    use crate::prelude::*;
    use super::super::*;
    use mzpeaks::{feature::{Feature, ChargedFeature}, MZ, Mass, IonMobility};

    macro_rules! assert_is_close {
        ($t1:expr, $t2:expr, $tol:expr, $label:literal) => {
            assert!(
                ($t1 - $t2).abs() < $tol,
                "Observed {} {}, expected {}, difference {}",
                $label,
                $t1,
                $t2,
                $t1 - $t2,
            );
        };
        ($t1:expr, $t2:expr, $tol:expr, $label:literal, $obj:ident) => {
            assert!(
                ($t1 - $t2).abs() < $tol,
                "Observed {} {}, expected {}, difference {} from {:?}",
                $label,
                $t1,
                $t2,
                $t1 - $t2,
                $obj
            );
        };
    }

    #[cfg(feature = "mzml")]
    #[test]
    fn test_loader_conversion() -> io::Result<()> {
        let fh = fs::File::open(
            "./test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz",
        )?;
        let handle = crate::io::compression::RestartableGzDecoder::new(io::BufReader::new(fh));
        let mzml_reader = crate::MzMLReader::new_indexed(handle);
        let mut frame_reader = mzml_reader
            .try_into_frame_source::<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>>()
            .unwrap();
        let frame = frame_reader
            .get_frame_by_id("merged=42926 frame=9728 scanStart=1 scanEnd=705")
            .unwrap();
        assert_eq!(frame.ms_level(), 1);
        Ok(())
    }

    #[cfg(feature = "mzml")]
    #[test]
    fn test_reader_wrapper_iter() -> io::Result<()> {
        use mzpeaks::feature::ChargedFeature;

        let group = crate::mz_read!("./test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz".as_ref(), reader => {
            let mut wrapper: Generic3DIonMobilityFrameSource<_, _, _, Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>> = Generic3DIonMobilityFrameSource::new(reader);
            let mut group_iter = wrapper.into_groups();
            group_iter.start_from_id("merged=42926 frame=9728 scanStart=1 scanEnd=705")?;
            group_iter.next().unwrap()
        })?;

        assert_eq!(group.lowest_ms_level().unwrap(), 1);
        assert_eq!(group.highest_ms_level().unwrap(), 2);

        Ok(())
    }

    #[test]
    fn test_reader_wrapper_extract() -> io::Result<()> {
        crate::mz_read!("./test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz".as_ref(), reader => {
            let mut wrapper: Generic3DIonMobilityFrameSource<_, _, _, Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>> = Generic3DIonMobilityFrameSource::new(reader);
            let mut frame = wrapper.get_frame_by_id("merged=42926 frame=9728 scanStart=1 scanEnd=705").unwrap();
            assert_eq!(frame.ms_level(), 1);
            assert_eq!(frame.polarity(), ScanPolarity::Positive);

            #[cfg(feature = "mzsignal")]
            {
                frame.extract_features_simple(Tolerance::PPM(15.0), 2, 0.1, None)?;

                let fm = frame.features.as_ref().unwrap();
                let query = 1456.95;
                let hits = fm.all_features_for(query, Tolerance::PPM(15.0));
                assert_eq!(hits.len(), 2);

                let feat = &hits[0];

                if false {
                    let mut writer = io::BufWriter::new(std::fs::File::create("tmp/peaks_over_time.txt")?);
                    let arrays = frame.arrays.as_ref().unwrap();
                    for (im, arrs) in arrays.ion_mobility_dimension.iter().zip(arrays.arrays.iter()) {
                        let mzs = arrs.mzs()?;
                        let intensities = arrs.intensities()?;
                        for (mz, int) in mzs.iter().zip(intensities.iter()) {
                            writeln!(writer, "{mz}\t{int}\t{im}")?;
                        }
                    }
                }

                if false {
                    let mut writer = io::BufWriter::new(std::fs::File::create("tmp/features_graph.txt")?);
                    writer.write_all(b"feature_id\tmz\trt\tintensity\n")?;
                    for (i, f) in fm.iter().enumerate() {
                        for (mz, rt, inten) in f.iter() {
                            writer.write_all(format!("{i}\t{mz}\t{rt}\t{inten}\n").as_bytes())?;
                        }
                    }
                }

                eprintln!("{:?}\t{:?}\t{:?}\t{}", feat.start_time(), feat.end_time(), feat.apex_time(), feat.len());

                assert_is_close!(feat.start_time().unwrap(), 0.9547949817459156, 1e-3, "start_time");
                assert_is_close!(feat.end_time().unwrap(), 1.2638564827665548, 1e-3, "end_time");
                assert_is_close!(feat.apex_time().unwrap(), 1.212666, 1e-3, "apex_time");
                assert_eq!(feat.len(), 83);
            }
        })?;
        Ok(())
    }
}
