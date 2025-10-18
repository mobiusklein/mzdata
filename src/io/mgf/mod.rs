//! Read and write [MGF](https://www.matrixscience.com/help/data_file_help.html#GEN) files.
//! Supports random access when reading from a source that supports [`io::Seek`](std::io::Seek).
//!
//! **Requires the `mgf` feature, enabled by default**
#![cfg(feature = "mgf")]
mod reader;
mod writer;

pub use reader::{
    DefaultTitleIndexing, DefaultTitleParser, MGFError, MGFParserState, MGFReader, MGFReaderType,
    ScansIndexing, TPPTitleParser, TPPTitleParsingNativeIDIndexing,
    TPPTitleParsingScanNumberIndexing, MGFIndexing, MGFTitle,
};
pub use writer::{
    MGFHeaderStyle, MGFWriter, MGFWriterType, MZDataMGFStyle, SimpleMGFStyle, SimpleMGFWriter,
    SimpleMGFWriterType,
};

#[cfg(feature = "async_partial")]
mod async_reader;

#[cfg(feature = "async_partial")]
pub use crate::io::mgf::async_reader::{
    MGFReader as AsyncMGFReader, MGFReaderType as AsyncMGFReaderType,
};

pub fn is_mgf(buf: &[u8]) -> bool {
    let needle = b"BEGIN IONS";
    buf.windows(needle.len()).any(|window| window == needle)
}

#[cfg(test)]
mod test {
    use crate::io::mgf::reader::{
        TPPTitleParsingNativeIDIndexing, TPPTitleParsingScanNumberIndexing,
    };
    use crate::io::DetailLevel;
    use crate::spectrum::RefPeakDataLevel;
    use crate::CentroidSpectrum;
    use crate::{io::RestartableGzDecoder, prelude::*};
    use mzpeaks::{CentroidPeak, DeconvolutedPeak, IndexedCoordinate};

    use super::*;
    use std::{fs, io, path};

    #[test]
    fn test_reader() {
        let path = path::Path::new("./test/data/small.mgf");
        let file = fs::File::open(path).expect("Test file doesn't exist");
        let reader = MGFReaderType::<_, CentroidPeak, DeconvolutedPeak>::new(file);
        let mut ms1_count = 0;
        let mut msn_count = 0;
        for scan in reader {
            let level = scan.ms_level();
            if level == 1 {
                ms1_count += 1;
            } else {
                msn_count += 1;
            }
        }
        assert_eq!(ms1_count, 0);
        assert_eq!(msn_count, 35);
    }

    #[test]
    fn test_reader_indexed() {
        let path = path::Path::new("./test/data/small.mgf");
        let file = fs::File::open(path).expect("Test file doesn't exist");
        let mut reader = MGFReaderType::<_, CentroidPeak, DeconvolutedPeak>::new_indexed(file);

        let n = reader.len();
        let mut ms1_count = 0;
        let mut msn_count = 0;

        for i in (0..n).rev() {
            let scan = reader.get_spectrum_by_index(i).expect("Missing spectrum");
            let level = scan.ms_level();
            if level == 1 {
                ms1_count += 1;
            } else {
                msn_count += 1;
            }
            let centroided: CentroidSpectrum = scan.try_into().unwrap();
            centroided.peaks.iter().for_each(|p| {
                (centroided.peaks[p.get_index() as usize]).mz();
            })
        }
        assert_eq!(ms1_count, 0);
        assert_eq!(msn_count, 35);
    }

    #[test]
    fn test_writer() -> io::Result<()> {
        let buff: Vec<u8> = Vec::new();
        let inner_writer = io::Cursor::new(buff);
        let mut writer = MGFWriter::new(inner_writer);

        let path = path::Path::new("./test/data/small.mgf");
        let file = fs::File::open(path).expect("Test file doesn't exist");
        let mut reader = MGFReader::new(file);

        for scan in reader.iter() {
            writer.write(&scan)?;
        }
        writer.flush()?;
        let inner_writer = writer.handle.into_inner()?;
        let buffer = inner_writer.into_inner();
        let reader2 = MGFReader::new(io::Cursor::new(buffer));
        assert_eq!(reader2.len(), reader.len());

        // Not including platform-specific line endings
        Ok(())
    }

    #[test]
    fn test_read_unsupported() -> io::Result<()> {
        let path = path::Path::new("./test/data/small.mgf");
        let file = fs::File::open(path).expect("Test file doesn't exist");
        let mut reader = MGFReader::new(file);

        assert!(reader.get_chromatogram_by_id("not real").is_none());
        assert!(reader.get_chromatogram_by_index(0).is_none());
        assert!(reader.spectrum_count_hint().is_none());

        reader = MGFReaderType::open_path(path)?;
        assert_eq!(reader.spectrum_count_hint().unwrap() as usize, reader.len());

        assert!(reader.run_description().is_some());
        Ok(())
    }

    #[test_log::test]
    fn test_reader_tpp_native_id() {
        let path = path::Path::new("./test/data/tpp/small.mgf");
        let file = fs::File::open(path).expect("Test file doesn't exist");
        let mut reader = MGFReaderType::<_, CentroidPeak, DeconvolutedPeak>::new_indexed_with(
            file,
            Box::new(TPPTitleParsingNativeIDIndexing::default()),
        );

        let mut ms1_count = 0;
        let mut msn_count = 0;

        let scan = reader
            .get_spectrum_by_id("controllerType=0 controllerNumber=1 scan=48")
            .expect("Failed to retrieve by ID");
        assert_eq!(scan.index(), 33);

        for (k, _) in reader.get_index().clone().iter().rev() {
            let scan = reader.get_spectrum_by_id(k).expect("Missing spectrum");
            let level = scan.ms_level();
            if level == 1 {
                ms1_count += 1;
            } else {
                msn_count += 1;
            }
            let centroided: CentroidSpectrum = scan.try_into().unwrap();
            centroided.peaks.iter().for_each(|p| {
                (centroided.peaks[p.get_index() as usize]).mz();
            })
        }
        assert_eq!(ms1_count, 0);
        assert_eq!(msn_count, 34);
    }

    #[test_log::test]
    fn test_reader_tpp_scan_num() {
        let path = path::Path::new("./test/data/tpp/small.mgf");
        let file = fs::File::open(path).expect("Test file doesn't exist");
        let mut reader = MGFReaderType::<_, CentroidPeak, DeconvolutedPeak>::new_indexed_with(
            file,
            Box::new(TPPTitleParsingScanNumberIndexing::default()),
        );

        let mut ms1_count = 0;
        let mut msn_count = 0;

        let scan = reader
            .get_spectrum_by_id("48")
            .expect("Failed to retrieve by ID");
        assert_eq!(scan.index(), 33);

        for (k, _) in reader.get_index().clone().iter().rev() {
            let scan = reader.get_spectrum_by_id(k).expect("Missing spectrum");
            let level = scan.ms_level();
            if level == 1 {
                ms1_count += 1;
            } else {
                msn_count += 1;
            }
            let centroided: CentroidSpectrum = scan.try_into().unwrap();
            centroided.peaks.iter().for_each(|p| {
                (centroided.peaks[p.get_index() as usize]).mz();
            })
        }
        assert_eq!(ms1_count, 0);
        assert_eq!(msn_count, 34);
    }

    #[test]
    fn test_read_charged_complex() -> io::Result<()> {
        let fh = io::BufReader::new(fs::File::open("./test/data/processed_batch.mgf.gz")?);
        let fh = RestartableGzDecoder::new(fh);
        let mut reader = MGFReader::new_indexed(fh);

        let mut scan = reader.next().unwrap();
        scan.try_build_deconvoluted_centroids().unwrap();
        assert!(matches!(scan.peaks(), RefPeakDataLevel::Deconvoluted(_)));
        scan.try_build_centroids()
            .expect_err("Expected not to find");
        assert!(matches!(scan.peaks(), RefPeakDataLevel::Deconvoluted(_)));
        scan.update_summaries();

        assert_eq!(scan.index(), 0);
        let summaries = scan.peaks().fetch_summaries();
        assert!(
            (summaries.tic - 3758148.3).abs() < 1e-3,
            "TIC error {}",
            summaries.tic - 3758148.3
        );
        assert!(
            (summaries.base_peak.mz - 443.2600402832031).abs() < 1e-3,
            "BP m/z {}",
            (summaries.base_peak.mz - 443.2600402832031)
        );
        assert!(
            (summaries.mz_range.0 > 120.0 && summaries.mz_range.0 < 120.1),
            "{:?}",
            summaries.mz_range
        );
        assert_eq!(summaries.count, 492);

        reader.read_into(&mut scan).unwrap();
        assert_eq!(scan.index(), 1);

        let sid = "MouseBrain-Z-T-1.25740.25740.2 File:\"MouseBrain-Z-T-1.raw\", NativeID:\"controllerType=0 controllerNumber=1 scan=25740\"";
        scan = reader.get_spectrum_by_id(sid).unwrap();
        assert_eq!(scan.id(), sid);
        assert_eq!(scan.index(), 0);
        reader.reset();

        assert_eq!(*reader.detail_level(), DetailLevel::Full);
        scan = reader.start_from_index(30).unwrap().next().unwrap();
        assert!(!scan.peaks().is_empty());
        assert_eq!(scan.index(), 30);
        let time = scan.start_time();

        reader.set_detail_level(DetailLevel::MetadataOnly);
        scan = reader.start_from_id(sid).unwrap().next().unwrap();
        assert_eq!(scan.index(), 0);
        assert!(scan.peaks().is_empty());

        scan = reader.start_from_time(time).unwrap().next().unwrap();
        assert_eq!(scan.index(), 30);
        assert!(scan.peaks().is_empty());
        Ok(())
    }

    #[cfg(feature = "async")]
    mod async_tests {
        use super::*;
        use futures::StreamExt;
        use tokio::fs;

        #[tokio::test]
        async fn test_reader() {
            let path = path::Path::new("./test/data/small.mgf");
            let file = fs::File::open(path).await.expect("Test file doesn't exist");
            let mut reader = AsyncMGFReaderType::<_>::new(file).await;
            let mut ms1_count = 0;
            let mut msn_count = 0;
            while let Some(scan) = reader.read_next().await {
                let level = scan.ms_level();
                if level == 1 {
                    ms1_count += 1;
                } else {
                    msn_count += 1;
                }
            }
            assert_eq!(ms1_count, 0);
            assert_eq!(msn_count, 35);
        }

        #[tokio::test]
        async fn test_reader_indexed() {
            let path = path::Path::new("./test/data/small.mgf");
            let file = fs::File::open(path).await.expect("Test file doesn't exist");
            let mut reader =
                AsyncMGFReaderType::<_, CentroidPeak, DeconvolutedPeak>::new_indexed(file).await;

            let n = reader.len();
            let mut ms1_count = 0;
            let mut msn_count = 0;

            for i in (0..n).rev() {
                let scan = reader.get_spectrum_by_index(i).await.unwrap();
                let level = scan.ms_level();
                if level == 1 {
                    ms1_count += 1;
                } else {
                    msn_count += 1;
                }
                let centroided: CentroidSpectrum = scan.try_into().unwrap();
                centroided.peaks.iter().for_each(|p| {
                    (centroided.peaks[p.get_index() as usize]).mz();
                })
            }
            assert_eq!(ms1_count, 0);
            assert_eq!(msn_count, 35);

            ms1_count = 0;
            msn_count = 0;
            reader.reset().await;

            let mut stream = reader.as_stream();
            while let Some(scan) = stream.next().await {
                let level = scan.ms_level();
                if level == 1 {
                    ms1_count += 1;
                } else {
                    msn_count += 1;
                }
                let centroided: CentroidSpectrum = scan.try_into().unwrap();
                centroided.peaks.iter().for_each(|p| {
                    (centroided.peaks[p.get_index() as usize]).mz();
                })
            }

            assert_eq!(ms1_count, 0);
            assert_eq!(msn_count, 35);
        }
    }
}
