//! Tests for imzML reader functionality

#[cfg(test)]
mod test {
    use std::io;
    use crate::prelude::*;
    use crate::io::imzml::{is_imzml, ImzMLReader};
    #[cfg(feature = "imzml")]
    use crate::io::mzml::MzMLReader;
    use crate::meta::MSDataFileMetadata;

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
        let x = event.get_param_by_curie(&crate::curie!(IMS:1000050)).unwrap();
        assert_eq!(x.to_i64(), Ok(1));
        let y = event.get_param_by_curie(&crate::curie!(IMS:1000051)).unwrap();
        assert_eq!(y.to_i64(), Ok(1));

        let arrays = spec.raw_arrays().unwrap();
        let arr = arrays.mzs()?;
        assert_eq!(arr.len(), 8399);

        reader = ImzMLReader::open_path("test/data/imaging/Example_Processed.imzML")?;
        let spec = reader.get_spectrum_by_index(0).unwrap();
        let acq = spec.acquisition();
        let event = &acq.scans[0];
        let x = event.get_param_by_curie(&crate::curie!(IMS:1000050)).unwrap();
        assert_eq!(x.to_i64(), Ok(1));
        let y = event.get_param_by_curie(&crate::curie!(IMS:1000051)).unwrap();
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
        assert!(!settings.id.is_empty(), "ScanSettings id should be non-empty");
        assert!(!settings.params.is_empty(), "ScanSettings params should be non-empty");
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
