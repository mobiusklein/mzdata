use std::io;
use std::sync::Arc;

use js_sys::{Array, Object, Reflect};
use mzdata::spectrum::utils::HasIonMobility;
use mzpeaks::feature::{Feature, ChargedFeature};
use mzpeaks::{CentroidPeak, IonMobility, MZ, Mass, DeconvolutedPeak};
use wasm_bindgen::prelude::*;

use mzdata::io::{IMMZReaderType, MZReaderType};
use mzdata::prelude::*;
use mzdata::spectrum::{MultiLayerIonMobilityFrame, MultiLayerSpectrum, SignalContinuity};

use crate::binds::{WebIonMobilityFrame, WebSpectrum};

#[derive(Debug)]
pub struct SharedBuffer {
    buffer: Arc<Vec<u8>>,
}

impl SharedBuffer {
    fn new(buffer: Arc<Vec<u8>>) -> Self {
        Self { buffer }
    }
}

impl Clone for SharedBuffer {
    fn clone(&self) -> Self {
        Self {
            buffer: self.buffer.clone(),
        }
    }
}

impl AsRef<[u8]> for SharedBuffer {
    fn as_ref(&self) -> &[u8] {
        &self.buffer
    }
}

type BufferType = SharedBuffer;

type ReaderType = MZReaderType<io::Cursor<BufferType>, CentroidPeak, mzpeaks::DeconvolutedPeak>;
type IMReaderType = IMMZReaderType<
    io::Cursor<BufferType>,
    Feature<MZ, IonMobility>,
    ChargedFeature<Mass, IonMobility>,
    CentroidPeak,
    mzpeaks::DeconvolutedPeak,
>;

#[wasm_bindgen]
pub struct MemWebMZReader {
    handle: ReaderType,
    peak_picking: bool,
    buffer_handle: Option<SharedBuffer>,
}

impl MemWebMZReader {
    pub fn get_mut(&mut self) -> &mut ReaderType {
        &mut self.handle
    }

    pub fn get_ref(&self) -> &ReaderType {
        &self.handle
    }

    pub fn buffer(&self) -> Option<&SharedBuffer> {
        self.buffer_handle.as_ref()
    }
}

#[wasm_bindgen]
impl MemWebMZReader {
    pub fn from_buffer(handle: js_sys::Uint8Array) -> Self {
        let n = handle.length() as usize;
        let mut buf = Vec::with_capacity(n);
        buf.resize(n, 0);
        handle.copy_to(&mut buf);

        let buf = SharedBuffer::new(Arc::new(buf));
        Self {
            handle: MZReaderType::open_read_seek(io::Cursor::new(buf.clone())).unwrap(),
            peak_picking: false,
            buffer_handle: Some(buf),
        }
    }

    pub fn set_data_loading(&mut self, load_data: bool) {
        if load_data {
            self.handle.set_detail_level(mzdata::io::DetailLevel::Full);
        } else {
            self.handle
                .set_detail_level(mzdata::io::DetailLevel::MetadataOnly);
        }
    }

    pub fn set_peak_picking(&mut self, pick_peaks: bool) {
        self.peak_picking = pick_peaks;
    }

    #[wasm_bindgen(getter)]
    pub fn file_format(&self) -> String {
        self.handle.as_format().to_string()
    }

    #[wasm_bindgen(getter)]
    pub fn length(&self) -> usize {
        self.handle.len()
    }

    fn convert_spectrum(
        &self,
        mut spectrum: MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak>,
    ) -> WebSpectrum {
        // if spectrum.peaks.is_none() && spectrum.deconvoluted_peaks.is_none() {
        //     spectrum.pick_peaks(1.0).unwrap();
        //     spectrum.description_mut().signal_continuity = SignalContinuity::Centroid;
        // }
        if self.peak_picking && spectrum.signal_continuity() == SignalContinuity::Profile {
            spectrum.pick_peaks(1.0).unwrap();
            spectrum.description_mut().signal_continuity = SignalContinuity::Centroid;
        }
        WebSpectrum::from(spectrum)
    }

    pub fn get_spectrum_by_id(&mut self, id: &str) -> Option<WebSpectrum> {
        let spectrum = self.handle.get_spectrum_by_id(id)?;
        Some(self.convert_spectrum(spectrum))
    }

    pub fn get_spectrum_by_index(&mut self, index: usize) -> Option<WebSpectrum> {
        let spectrum = self.handle.get_spectrum_by_index(index)?;
        Some(self.convert_spectrum(spectrum))
    }

    pub fn get_spectrum_by_time(&mut self, time: f64) -> Option<WebSpectrum> {
        let spectrum = self.handle.get_spectrum_by_time(time)?;
        Some(self.convert_spectrum(spectrum))
    }

    pub fn next(&mut self) -> Option<WebSpectrum> {
        self.handle.next().map(|s| self.convert_spectrum(s))
    }

    pub fn start_from_index(&mut self, index: usize) {
        self.handle.start_from_index(index).unwrap();
    }

    pub fn start_from_time(&mut self, time: f64) {
        self.handle.start_from_time(time).unwrap();
    }

    pub fn group_at(&mut self, index: usize) -> Option<Object> {
        let mut it = self.handle.iter();
        let mut it = it.groups();
        it.start_from_index(index).unwrap();
        let group = it.next();
        drop(it);

        let obj = Object::new();
        let group = match group {
            Some(group) => group,
            None => return None,
        };

        let (prec, products) = group.into_parts();
        if let Some(prec) = prec {
            let spec = self.convert_spectrum(prec);
            let val = JsValue::from(spec);
            Reflect::set(&obj, &JsValue::from_str(&"precursor"), &val).unwrap();
        } else {
            Reflect::set(&obj, &JsValue::from_str(&"precursor"), &JsValue::null()).unwrap();
        }

        let products: Array = products
            .into_iter()
            .map(|spec| JsValue::from(self.convert_spectrum(spec)))
            .collect();
        Reflect::set(&obj, &JsValue::from_str(&"products"), &products).unwrap();
        Some(obj)
    }

    pub fn to_frame_reader(&mut self) -> Result<MemWebIMMZReader, String> {
        if let Some(im) = self.handle.has_ion_mobility() {
            if matches!(im, HasIonMobility::Dimension) {
                let buffer = self.buffer().ok_or("No shared buffer found")?.clone();
                let reader = ReaderType::open_read_seek(io::Cursor::new(buffer))
                    .map_err(|e| format!("Failed to reuse shared buffer: {e}"))?;
                let handle = reader.into_frame_source();
                let this = MemWebIMMZReader {
                    handle,
                    feature_extraction: self.peak_picking,
                };
                return Ok(this);
            } else {
                return Err("No ion mobility data found".to_string());
            }
        }
        return Err("No ion mobility data found".to_string());
    }
}

#[wasm_bindgen]
pub struct MemWebIMMZReader {
    handle: IMReaderType,
    feature_extraction: bool,
}

impl MemWebIMMZReader {
    pub fn get_mut(&mut self) -> &mut IMReaderType {
        &mut self.handle
    }

    pub fn get_ref(&self) -> &IMReaderType {
        &self.handle
    }
}

#[wasm_bindgen]
impl MemWebIMMZReader {
    pub fn from_buffer(handle: js_sys::Uint8Array) -> Self {
        let n = handle.length() as usize;
        let mut buf = Vec::with_capacity(n);
        buf.resize(n, 0);
        handle.copy_to(&mut buf);
        let buf = SharedBuffer::new(Arc::new(buf));
        let inner_reader = MZReaderType::open_read_seek(io::Cursor::new(buf)).unwrap();

        let handle = inner_reader.into_frame_source();

        Self {
            handle: handle,
            feature_extraction: false,
        }
    }

    pub fn set_data_loading(&mut self, load_data: bool) {
        if load_data {
            self.handle.set_detail_level(mzdata::io::DetailLevel::Full);
        } else {
            self.handle
                .set_detail_level(mzdata::io::DetailLevel::MetadataOnly);
        }
    }

    pub fn set_feature_extraction(&mut self, feature_extraction: bool) {
        self.feature_extraction = feature_extraction;
    }

    #[wasm_bindgen(getter)]
    pub fn length(&self) -> usize {
        self.handle.len()
    }

    fn convert_frame(
        &self,
        mut frame: MultiLayerIonMobilityFrame<
            Feature<MZ, IonMobility>,
            ChargedFeature<Mass, IonMobility>,
        >,
    ) -> WebIonMobilityFrame {
        if self.feature_extraction && frame.features.is_none() {
            frame
                .extract_features_simple(Tolerance::PPM(15.0), 2, 0.1, None)
                .unwrap();
            frame.description_mut().signal_continuity = SignalContinuity::Centroid;
        }
        WebIonMobilityFrame::from(frame)
    }

    pub fn get_frame_by_id(&mut self, id: &str) -> Option<WebIonMobilityFrame> {
        let frame = self.handle.get_frame_by_id(id)?;
        Some(self.convert_frame(frame))
    }

    pub fn get_frame_by_index(&mut self, index: usize) -> Option<WebIonMobilityFrame> {
        let frame = self.handle.get_frame_by_index(index)?;
        Some(self.convert_frame(frame))
    }

    pub fn get_frame_by_time(&mut self, time: f64) -> Option<WebIonMobilityFrame> {
        let frame = self.handle.get_frame_by_time(time)?;
        Some(self.convert_frame(frame))
    }

    pub fn next(&mut self) -> Option<WebIonMobilityFrame> {
        self.handle.next().map(|s| self.convert_frame(s))
    }

    pub fn start_from_index(&mut self, index: usize) {
        self.handle.start_from_index(index).unwrap();
    }

    pub fn start_from_time(&mut self, time: f64) {
        self.handle.start_from_time(time).unwrap();
    }

    pub fn group_at(&mut self, index: usize) -> Option<Object> {
        let mut it = self.handle.iter();
        let mut it = it.groups();
        it.start_from_index(index).unwrap();
        let group = it.next();
        drop(it);

        let obj = Object::new();
        let group = match group {
            Some(group) => group,
            None => return None,
        };

        let (prec, products) = group.into_parts();
        if let Some(prec) = prec {
            let spec = self.convert_frame(prec);
            let val = JsValue::from(spec);
            Reflect::set(&obj, &JsValue::from_str(&"precursor"), &val).unwrap();
        } else {
            Reflect::set(&obj, &JsValue::from_str(&"precursor"), &JsValue::null()).unwrap();
        }

        let products: Array = products
            .into_iter()
            .map(|spec| JsValue::from(self.convert_frame(spec)))
            .collect();
        Reflect::set(&obj, &JsValue::from_str(&"products"), &products).unwrap();
        Some(obj)
    }
}
