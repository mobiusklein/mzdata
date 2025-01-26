mod inference;
mod dispatch;
mod pipeline;

pub use dispatch::{MZReader, MZReaderType, MZReaderBuilder, IMMZReaderType};
#[cfg(feature = "async_partial")]
pub use dispatch::{AsyncMZReaderType, AsyncMZReader, AsyncMZReaderBuilder};



pub use inference::{infer_from_path, infer_from_stream, infer_format, MassSpectrometryFormat};

pub use pipeline::{MassSpectrometryReadWriteProcess, Source, Sink};

#[cfg(test)]
mod test {
    use std::{fs, io, path};

    use mzpeaks::{CentroidPeak, DeconvolutedPeak};

    use crate::{
        prelude::*,
        spectrum::{ArrayType, Spectrum},
        io::DetailLevel
    };

    use super::*;

    #[test]
    fn infer_mzml() {
        let path = path::Path::new("./test/data/small.mzML");
        assert!(path.exists());
        let (fmt, zipped) = infer_from_path(path);
        assert_eq!(fmt, MassSpectrometryFormat::MzML);
        assert!(!zipped);
    }

    #[test]
    fn infer_mgf() {
        let path = path::Path::new("./test/data/small.mgf");
        assert!(path.exists());
        let (fmt, zipped) = infer_from_path(path);
        assert_eq!(fmt, MassSpectrometryFormat::MGF);
        assert!(!zipped);
    }

    #[cfg(feature = "thermo")]
    #[test]
    fn infer_thermo() {
        let path = path::Path::new("./test/data/small.RAW");
        let (fmt, zipped) = infer_from_path(path);
        assert_eq!(fmt, MassSpectrometryFormat::ThermoRaw);
        assert!(!zipped);
    }

    #[test]
    fn infer_open() {
        let path = path::Path::new("./test/data/small.mzML");
        assert!(path.exists());
        if let Ok(mut reader) = MZReader::open_path(path) {
            assert_eq!(reader.len(), 48);
            assert_eq!(*reader.detail_level(), DetailLevel::Full);
            if let Some(spec) = reader.get_spectrum_by_index(10) {
                let spec: Spectrum = spec;
                assert!(spec.index() == 10);
                assert!(spec.id() == "controllerType=0 controllerNumber=1 scan=11");
                if let Some(data_arrays) = &spec.arrays {
                    assert!(data_arrays.has_array(&ArrayType::MZArray));
                    assert!(data_arrays.has_array(&ArrayType::IntensityArray));
                    let mzs = data_arrays.mzs().unwrap();
                    assert!(mzs.len() == 941);
                }
            } else {
                panic!("Failed to retrieve spectrum by index")
            }

            assert_eq!(reader.get_spectrum_by_id("controllerType=0 controllerNumber=1 scan=11").unwrap().index(), 10);

            if let Some(spec) = reader.get_spectrum_by_time(0.358558333333) {
                assert_eq!(spec.index(), 34);
            } else {
                panic!("Failed to retrieve spectrum by time")
            }

        } else {
            panic!("Failed to open file")
        }
    }

    #[cfg(feature = "thermo")]
    #[test]
    fn infer_open_thermo() {
        let path = path::Path::new("./test/data/small.RAW");
        assert!(path.exists());
        if let Ok(mut reader) = MZReader::open_path(path) {
            assert_eq!(reader.len(), 48);
            assert_eq!(*reader.detail_level(), DetailLevel::Full);
            if let Some(spec) = reader.get_spectrum_by_index(10) {
                let spec: Spectrum = spec;
                assert_eq!(spec.index(), 10);
                assert_eq!(spec.id(), "controllerType=0 controllerNumber=1 scan=11");
                assert_eq!(spec.peaks().len(), 941);
            } else {
                panic!("Failed to retrieve spectrum by index")
            }

            assert_eq!(reader.get_spectrum_by_id("controllerType=0 controllerNumber=1 scan=11").unwrap().index(), 10);

            if let Some(spec) = reader.get_spectrum_by_time(0.358558333333) {
                assert_eq!(spec.index(), 34);
            } else {
                panic!("Failed to retrieve spectrum by time")
            }

        } else {
            panic!("Failed to open file")
        }
    }

    #[test]
    fn test_source_conv() -> io::Result<()> {
        let s = Source::<CentroidPeak, DeconvolutedPeak>::from("text/path".as_ref());
        assert!(matches!(s, Source::PathLike(_)));

        let fh = Box::new(io::BufReader::new(fs::File::open("./test/data/small.mgf")?)) as Box<dyn SeekRead + Send>;
        let rs: Source<CentroidPeak, DeconvolutedPeak> = (fh, MassSpectrometryFormat::MGF).into();
        assert!(matches!(rs, Source::Reader(_, _)));

        Ok(())
    }

    #[test]
    fn test_dispatch_mzreader() -> io::Result<()> {
        let mut reader = MZReader::open_path("./test/data/small.mzML")?;

        let n = reader.len();
        let n_ms1 = reader.iter().filter(|s| s.ms_level() == 1).count();
        let n_msn = reader.iter().filter(|s| s.ms_level() == 2).count();

        assert_eq!(n, 48);
        assert_eq!(n, n_ms1 + n_msn);
        Ok(())
    }

    #[test]
    fn test_infer_stream() -> io::Result<()> {
        let mut mzml_file = fs::File::open("./test/data/small.mzML")?;
        let (form, gzip) = infer_from_stream(&mut mzml_file)?;
        assert_eq!(form, MassSpectrometryFormat::MzML);
        assert!(!gzip);

        mzml_file = fs::File::open("./test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz")?;
        let (form, gzip) = infer_from_stream(&mut mzml_file)?;
        assert_eq!(form, MassSpectrometryFormat::MzML);
        assert!(gzip);

        let mut mgf_file = fs::File::open("./test/data/small.mgf")?;
        let (form, gzip) = infer_from_stream(&mut mgf_file)?;
        assert_eq!(form, MassSpectrometryFormat::MGF);
        assert!(!gzip);
        Ok(())
    }
}
