//! Tests for imzML reader functionality

#[cfg(test)]
mod test {
    use crate::io::imzml::{is_imzml, ImzMLReader};

    #[test]
    fn test_is_imzml_detection() {
        // Test with imagingMS.obo reference
        let imaging_content = b"<cv URI=\"http://www.maldi-msi.org/download/imzml/imagingMS.obo\"";
        assert!(is_imzml(imaging_content));

        // Test with IMS term
        let ims_content = b"<cvParam cvRef=\"IMS\" accession=\"IMS:1000080\"";
        assert!(is_imzml(ims_content));

        // Test with regular mzML content
        let mzml_content = b"<mzML xmlns=\"http://psi.hupo.org/ms/mzml\"";
        assert!(!is_imzml(mzml_content));
    }

    #[test] 
    fn test_imzml_reader_creation() {
        // This is a basic compilation test since we don't have test files
        // In a real scenario, you would test with actual imzML files
        
        // The reader should be creatable
        let _reader_type = std::marker::PhantomData::<ImzMLReader>;
        
        // Test IBD path derivation
        use crate::io::imzml::IbdFile;
        let imzml_path = std::path::Path::new("test.imzML");
        let ibd_path = IbdFile::derive_ibd_path(imzml_path);
        assert_eq!(ibd_path, std::path::Path::new("test.ibd"));
    }
}
