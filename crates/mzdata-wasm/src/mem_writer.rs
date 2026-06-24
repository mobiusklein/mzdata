use std::io;

use mzdata::meta::MSDataFileMetadata;
use wasm_bindgen::prelude::*;

use mzdata::prelude::*;
use mzdata::io::mgf::{MGFWriterType, MZDataMGFStyle};
use mzdata::io::mzml::MzMLWriterType;


use crate::WebSpectrum;

#[wasm_bindgen]
pub fn write_spectra_to_mgf(spectra: Vec<WebSpectrum>) -> String {
    let mut buffer = Vec::new();
    let writer = io::Cursor::new(&mut buffer);
    let mut writer: MGFWriterType<
        io::Cursor<&mut Vec<u8>>,
        mzpeaks::CentroidPeak,
        mzpeaks::DeconvolutedPeak,
        MZDataMGFStyle,
    > = MGFWriterType::new(writer);

    for spec in spectra.iter() {
        writer.write(spec.as_ref()).unwrap();
    }
    drop(writer);
    let content = String::from_utf8_lossy(&buffer);
    content.to_string()
}


#[wasm_bindgen]
pub fn _write_spectra_to_mzml(spectra: Vec<WebSpectrum>, reader: &crate::mem_reader::MemWebMZReader) -> String {
    let mut buffer = Vec::new();
    let writer = io::Cursor::new(&mut buffer);
    let mut writer: MzMLWriterType<
        io::Cursor<&mut Vec<u8>>,
        mzpeaks::CentroidPeak,
        mzpeaks::DeconvolutedPeak,
    > = MzMLWriterType::new(writer);
    writer.copy_metadata_from(reader.get_ref());
    for spec in spectra.iter() {
        writer.write(spec.as_ref()).unwrap();
    }
    drop(writer);
    let content = String::from_utf8_lossy(&buffer);
    content.to_string()
}