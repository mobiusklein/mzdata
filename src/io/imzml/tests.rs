//! Tests for imzML reader functionality

#[cfg(test)]
mod test {
    use crate::io::imzml::{is_imzml, ImzMLReader};
    #[cfg(feature = "imzml")]
    use crate::io::mzml::MzMLReader;
    use crate::meta::MSDataFileMetadata;
    use crate::prelude::*;
    use std::io;

    #[test]
    fn test_open_path_reports_errors() -> io::Result<()> {
        let dir = std::env::temp_dir().join(format!("mzdata-imzml-open-{}", std::process::id()));
        std::fs::create_dir(&dir)?;
        let path = dir.join("input.imzML");
        let xml = std::fs::read("test/data/imaging/Example_Processed.imzML")?;
        std::fs::write(&path, &xml)?;
        std::fs::write(path.with_extension("ibd"), [])?;
        let truncated = ImzMLReader::open_path(&path).err().map(|e| e.kind());
        let dispatched = crate::MZReader::open_path(&path).err().map(|e| e.kind());
        let mut xml = xml;
        let offset = xml.windows(11).position(|s| s == b"IMS:1000031").unwrap();
        xml[offset..offset + 11].copy_from_slice(b"IMS:1009999");
        std::fs::write(&path, xml)?;
        std::fs::write(
            path.with_extension("ibd"),
            std::fs::read("test/data/imaging/Example_Processed.ibd")?,
        )?;
        let metadata = ImzMLReader::open_path(&path).err().map(|e| e.kind());
        std::fs::remove_dir_all(dir)?;
        assert_eq!(truncated, Some(io::ErrorKind::UnexpectedEof));
        assert_eq!(dispatched, Some(io::ErrorKind::UnexpectedEof));
        assert_eq!(metadata, Some(io::ErrorKind::InvalidData));
        Ok(())
    }

    #[test]
    fn test_fallible_construction() -> io::Result<()> {
        use crate::io::DetailLevel;
        use std::io::Cursor;

        for name in ["Example_Continuous", "Example_Processed"] {
            let xml = std::fs::read(format!("test/data/imaging/{name}.imzML"))?;
            let ibd = std::fs::read(format!("test/data/imaging/{name}.ibd"))?;
            for level in [
                DetailLevel::Full,
                DetailLevel::Lazy,
                DetailLevel::MetadataOnly,
            ] {
                let mut reader = ImzMLReader::try_with_buffer_capacity_and_detail_level(
                    Cursor::new(&xml),
                    Cursor::new(&ibd),
                    512,
                    level,
                )?;
                let mut reference = ImzMLReader::new(Cursor::new(&xml), Cursor::new(&ibd));
                reference.set_detail_level(level);
                assert_eq!(reader.len(), 9);
                assert_eq!(reader.imzml_metadata.uuid, reference.imzml_metadata.uuid);
                assert_eq!(
                    reader.imzml_metadata.data_mode,
                    reference.imzml_metadata.data_mode
                );
                for index in 0..9 {
                    let expected = reference.read_next().unwrap();
                    let actual = reader.get_spectrum_by_index(index).unwrap();
                    assert_eq!(actual.description, expected.description);
                    if level == DetailLevel::MetadataOnly {
                        assert_eq!(actual.arrays.is_some(), expected.arrays.is_some());
                        for (_, array) in actual.raw_arrays().unwrap().iter() {
                            assert!(array.data.is_empty());
                        }
                    } else {
                        let actual = actual.raw_arrays().unwrap();
                        let expected = expected.raw_arrays().unwrap();
                        assert_eq!(actual.mzs()?.len(), 8399);
                        assert_eq!(actual.mzs()?, expected.mzs()?);
                        assert_eq!(actual.intensities()?, expected.intensities()?);
                    }
                }
            }
            for size in 0..16 {
                let error = ImzMLReader::try_new(Cursor::new(&xml), Cursor::new(&ibd[..size]))
                    .err()
                    .unwrap();
                assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof);
            }
            let mut mismatched = ibd.clone();
            mismatched[0] ^= 1;
            let mut reader = ImzMLReader::try_new(Cursor::new(&xml), Cursor::new(&mismatched))?;
            assert_eq!(
                reader
                    .read_next()
                    .unwrap()
                    .raw_arrays()
                    .unwrap()
                    .mzs()?
                    .len(),
                8399
            );
        }
        assert_eq!(
            ImzMLReader::try_new(Cursor::new([]), Cursor::new([]))
                .err()
                .unwrap()
                .kind(),
            io::ErrorKind::InvalidData
        );
        Ok(())
    }

    #[test]
    fn test_is_imzml_detection() {
        // Test with proper imzML cvList containing IMS
        let imzml_content = br#"
            <mzML xmlns="http://psi.hupo.org/ms/mzml">
                <cvList count="3">
                    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology"/>
                    <cv id="UO" fullName="Unit Ontology"/>
                    <cv id="IMS" fullName="Imaging MS Ontology"/>
                </cvList>
            </mzML>
        "#;
        assert!(is_imzml(imzml_content));

        // Test with regular mzML content (no IMS in cvList)
        let mzml_content = br#"
            <mzML xmlns="http://psi.hupo.org/ms/mzml">
                <cvList count="2">
                    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology"/>
                    <cv id="UO" fullName="Unit Ontology"/>
                </cvList>
            </mzML>
        "#;
        assert!(!is_imzml(mzml_content));

        // Test with empty/minimal content
        let minimal_content = b"<mzML xmlns=\"http://psi.hupo.org/ms/mzml\"";
        assert!(!is_imzml(minimal_content));
    }

    #[test]
    fn test_imzml_reader_type() {
        // This is a basic compilation test to ensure the type is correctly defined
        use std::fs::File;

        // The reader should be creatable with the correct type parameters
        let _reader_type = std::marker::PhantomData::<ImzMLReader<File, File>>;

        // Test IBD path derivation logic
        let imzml_path = std::path::Path::new("test.imzML");
        let ibd_path_lower = imzml_path.with_extension("ibd");
        let ibd_path_upper = imzml_path.with_extension("IBD");

        assert_eq!(ibd_path_lower, std::path::Path::new("test.ibd"));
        assert_eq!(ibd_path_upper, std::path::Path::new("test.IBD"));
    }

    #[test]
    fn test_imzml_read_operation() -> io::Result<()> {
        let mut reader = ImzMLReader::open_path("test/data/imaging/Example_Continuous.imzML")?;
        let spec = reader.get_spectrum_by_index(0).unwrap();
        let acq = spec.acquisition();
        let event = &acq.scans[0];
        let x = event
            .get_param_by_curie(&crate::curie!(IMS:1000050))
            .unwrap();
        assert_eq!(x.to_i64(), Ok(1));
        let y = event
            .get_param_by_curie(&crate::curie!(IMS:1000051))
            .unwrap();
        assert_eq!(y.to_i64(), Ok(1));

        let arrays = spec.raw_arrays().unwrap();
        let arr = arrays.mzs()?;
        assert_eq!(arr.len(), 8399);

        reader = ImzMLReader::open_path("test/data/imaging/Example_Processed.imzML")?;
        let spec = reader.get_spectrum_by_index(0).unwrap();
        let acq = spec.acquisition();
        let event = &acq.scans[0];
        let x = event
            .get_param_by_curie(&crate::curie!(IMS:1000050))
            .unwrap();
        assert_eq!(x.to_i64(), Ok(1));
        let y = event
            .get_param_by_curie(&crate::curie!(IMS:1000051))
            .unwrap();
        assert_eq!(y.to_i64(), Ok(1));

        let arrays = spec.raw_arrays().unwrap();
        let arr = arrays.mzs()?;
        assert_eq!(arr.len(), 8399);
        Ok(())
    }

    #[test]
    fn test_imzml_scan_settings_processed() -> io::Result<()> {
        let reader = ImzMLReader::open_path("test/data/imaging/Example_Processed.imzML")?;
        let settings_list = reader
            .scan_settings()
            .expect("ImzMLReader should expose scan_settings");
        assert_eq!(settings_list.len(), 1, "expected one scanSettings entry");

        let settings = &settings_list[0];
        assert!(
            !settings.id.is_empty(),
            "ScanSettings id should be non-empty"
        );
        assert!(
            !settings.params.is_empty(),
            "ScanSettings params should be non-empty"
        );
        assert!(settings.source_file_refs.is_empty());
        assert!(settings.targets.is_empty());
        Ok(())
    }

    #[test]
    fn test_mzml_scan_settings_empty() -> io::Result<()> {
        let reader = MzMLReader::open_path("test/data/small.mzML")?;
        let settings_list = reader
            .scan_settings()
            .expect("MzMLReader should expose scan_settings (even if empty)");
        assert!(
            settings_list.is_empty(),
            "plain mzML should have no scanSettings entries"
        );
        Ok(())
    }
}
