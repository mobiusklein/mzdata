use std::{
    fs::File,
    io::{self, BufReader, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use thiserror::Error;

use crate::spectrum::bindata::{BinaryDataArrayType, DataArray};

#[derive(Debug, Error)]
pub enum IbdError {
    #[error("An IO error occurred while reading IBD file: {0}")]
    IoError(#[from] io::Error),
    #[error("Invalid data type for IBD array: {0:?}")]
    InvalidDataType(BinaryDataArrayType),
    #[error("IBD file UUID mismatch: expected {expected:?}, found {found:?}")]
    UuidMismatch { expected: [u8; 16], found: [u8; 16] },
    #[error("Invalid offset or length for IBD data: offset={offset}, length={length}")]
    InvalidRange { offset: u64, length: u64 },
}

impl From<IbdError> for io::Error {
    fn from(value: IbdError) -> Self {
        match value {
            IbdError::IoError(e) => e,
            _ => Self::new(io::ErrorKind::Other, value),
        }
    }
}

/// UUID is 16 bytes stored in big-endian format at the start of IBD files
const UUID_SIZE: usize = 16;

/// Represents the two data storage modes in imzML
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IbdDataMode {
    /// All spectra share the same m/z values
    Continuous,
    /// Each spectrum has its own m/z and intensity arrays
    Processed,
    Unknown,
}

/// Handle for reading binary data from an imzML .ibd file
#[derive(Debug)]
pub struct IbdFile {
    reader: BufReader<File>,
    uuid: [u8; UUID_SIZE],
    data_mode: IbdDataMode,
    shared_mz_array: Option<Vec<f64>>,
}

impl IbdFile {
    /// Open an IBD file from a path
    pub fn open<P: AsRef<Path>>(path: P, data_mode: IbdDataMode) -> io::Result<Self> {
        let file = File::open(&path)?;
        let mut reader = BufReader::new(file);
        
        // Read the UUID from the first 16 bytes
        let mut uuid = [0u8; UUID_SIZE];
        reader.read_exact(&mut uuid)?;
        
        Ok(Self {
            reader,
            uuid,
            data_mode,
            shared_mz_array: None,
        })
    }

    /// Get the UUID of this IBD file
    pub fn uuid(&self) -> &[u8; UUID_SIZE] {
        &self.uuid
    }

    /// Set the data mode (determined from imzML metadata)
    pub fn set_data_mode(&mut self, mode: IbdDataMode) {
        self.data_mode = mode;
    }

    /// Set the shared m/z array for continuous mode
    pub fn set_shared_mz_array(&mut self, mz_array: Vec<f64>) {
        self.shared_mz_array = Some(mz_array);
    }

    /// Read binary data array from the IBD file at the specified offset and length
    pub fn read_array(
        &mut self,
        offset: usize,
        array_length: usize,
        data_type: BinaryDataArrayType,
    ) -> Result<Vec<u8>, IbdError> {
        // Seek to the specified offset (skip UUID)
        self.reader.seek(SeekFrom::Start(offset as u64))?;
        
        // Calculate the number of bytes to read based on data type and array length
        let element_size = data_type.size_of();
        
        let total_bytes = array_length * element_size;
        let mut buffer = vec![0u8; total_bytes];
        self.reader.read_exact(&mut buffer)?;
        
        Ok(buffer)
    }

    /// Read and decode a data array from the IBD file
    pub fn read_data_array(
        &mut self,
        offset: u64,
        array_length: u64,
        data_type: BinaryDataArrayType,
    ) -> Result<DataArray, IbdError> {
        let raw_data = self.read_array(offset as usize, array_length as usize, data_type)?;
        
        let mut data_array = DataArray::new();
        data_array.dtype = data_type;
        data_array.data = raw_data; // DataArray.data is Vec<u8>, so store the raw bytes
        
        Ok(data_array)
    }

    /// Read m/z data for a spectrum in continuous mode
    pub fn read_mz_array_continuous(&self) -> Option<&Vec<f64>> {
        self.shared_mz_array.as_ref()
    }

    /// Derive the IBD file path from an imzML file path
    /// TODO: Add support for custom IBD paths from imzML metadata
    pub fn derive_ibd_path<P: AsRef<Path>>(imzml_path: P) -> PathBuf {
        let path = imzml_path.as_ref();
        path.with_extension("ibd")
    }

    /// Get binary data from IBD file based on a data range query and populate the provided array
    pub fn get(&mut self, query: &crate::io::imzml::DataRangeQuery, array: &mut DataArray) -> Result<(), IbdError> {
        // Read the raw data from the IBD file
        let raw_data = self.read_array(
            query.offset,
            query.length,
            array.dtype
        )?;
        
        // Copy the data into the provided ByteArrayView
        array.data = raw_data;
        
        Ok(())
    }
}


