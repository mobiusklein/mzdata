/*!
Read and write [MGF](https://www.matrixscience.com/help/data_file_help.html#GEN) files.
Supports random access when reading from a source that supports [`io::Seek`].
*/
mod reader;
mod writer;

pub use reader::{MGFReaderType, MGFReader, is_mgf, MGFError, MGFParserState};
pub use writer::{MGFHeaderStyle, MGFWriter, MGFWriterType, SimpleMGFStyle, MZDataMGFStyle};

#[cfg(feature = "async")]
mod async_reader;

#[cfg(feature = "async")]
pub use crate::io::mgf::async_reader::{
    MGFReaderType as AsyncMGFReaderType,
    MGFReader as AsyncMGFReader
};


#[cfg(test)]
mod test {
    use mzpeaks::{IndexedCoordinate, CentroidPeak, DeconvolutedPeak};
    use crate::prelude::*;
    use crate::CentroidSpectrum;

    use super::*;
    use std::{fs, io, path};

    #[test]
    fn test_reader() {
        let path = path::Path::new("./test/data/small.mgf");
        let file = fs::File::open(path).expect("Test file doesn't exist");
        let reader = MGFReaderType::<_>::new(file);
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


    #[cfg(feature = "async")]
    mod async_tests {
        use super::*;
        use tokio::fs;

        #[tokio::test]
        async fn test_reader() {
            let path = path::Path::new("./test/data/small.mgf");
            let file = fs::File::open(path).await.expect("Test file doesn't exist");
            let mut reader = AsyncMGFReaderType::<_>::new(file);
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
            let mut reader = AsyncMGFReaderType::<_, CentroidPeak, DeconvolutedPeak>::new_indexed(file).await;

            let n = reader.len();
            let mut ms1_count = 0;
            let mut msn_count = 0;

            for i in (0..n).rev() {
                let scan = reader.get_spectrum_by_index(i).await.expect("Missing spectrum");
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
