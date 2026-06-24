use js_sys::{Array, Float32Array, Float64Array, Int32Array, Object, Reflect, Uint8Array};
use log::{self, info};

use mzdata::{
    prelude::*,
    spectrum::{
        bindata::ArrayRetrievalError, Activation, ArrayType, BinaryArrayMap,
        MultiLayerIonMobilityFrame, MultiLayerSpectrum, PeakDataLevel, RefPeakDataLevel, ScanEvent,
    },
    utils::mass_charge_ratio,
};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use mzpeaks::{
    feature::{Feature, SimpleFeature, ChargedFeature},
    peak::MZPoint,
    prelude::*,
    CentroidPeak, DeconvolutedPeak, IonMobility, Time, Tolerance, MZ, Mass
};

use mzsignal::feature_statistics::{FeatureTransform, FitPeaksOn, MultiPeakShapeFit};

use mzdata::{
    params::{Param, ParamDescribed, ValueRef},
    spectrum::{
        IsolationWindow, Precursor, ScanWindow, SelectedIon, SignalContinuity, SpectrumDescription,
    },
};

pub fn array_map_to_js(arrays: &BinaryArrayMap) -> Result<Object, ArrayRetrievalError> {
    let mut acc = Vec::new();
    for (name, data) in arrays.iter() {
        let name = if let ArrayType::NonStandardDataArray { name } = name {
            name.to_string()
        } else {
            name.as_param_const().name.to_string()
        };
        match data.dtype() {
            mzdata::spectrum::BinaryDataArrayType::Unknown
            | mzdata::spectrum::BinaryDataArrayType::ASCII => {
                let view = data.decode()?;
                let buf = Uint8Array::new_with_length(view.len() as u32);
                buf.copy_from(&view);
                acc.push(Array::of2(&JsValue::from_str(&name), &JsValue::from(buf)))
            }
            mzdata::spectrum::BinaryDataArrayType::Float64 => {
                let view = data.to_f64()?;
                let buf = Float64Array::new_with_length(view.len() as u32);
                buf.copy_from(&view);
                acc.push(Array::of2(&JsValue::from_str(&name), &JsValue::from(buf)))
            }
            mzdata::spectrum::BinaryDataArrayType::Float32 => {
                let view = data.to_f32()?;
                let buf = Float32Array::new_with_length(view.len() as u32);
                buf.copy_from(&view);
                acc.push(Array::of2(&JsValue::from_str(&name), &JsValue::from(buf)))
            }
            mzdata::spectrum::BinaryDataArrayType::Int64 => todo!(),
            mzdata::spectrum::BinaryDataArrayType::Int32 => {
                let view = data.to_i32()?;
                let buf = Int32Array::new_with_length(view.len() as u32);
                buf.copy_from(&view);
                acc.push(Array::of2(&JsValue::from_str(&name), &JsValue::from(buf)))
            }
        }
    }

    let res = Object::from_entries(&JsValue::from(acc)).expect("Failed to build array map");
    Ok(res)
}

#[wasm_bindgen(js_name = "Tolerance")]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct WebTolerance(Tolerance);

#[wasm_bindgen(js_class = "Tolerance")]
impl WebTolerance {
    pub fn ppm(value: f64) -> WebTolerance {
        WebTolerance(Tolerance::PPM(value))
    }

    pub fn da(value: f64) -> WebTolerance {
        WebTolerance(Tolerance::Da(value))
    }

    pub fn tol(&self) -> f64 {
        self.0.tol()
    }

    pub fn parse(text: &str) -> Result<Self, JsError> {
        let tol = text.parse().map_err(|e| JsError::new(&format!("{e}")))?;

        Ok(WebTolerance(tol))
    }

    #[wasm_bindgen(js_name = "toString")]
    pub fn to_string(&self) -> String {
        self.0.to_string()
    }

    #[wasm_bindgen(js_name = "toJSON")]
    pub fn to_json(&self) -> Result<JsValue, JsError> {
        serde_wasm_bindgen::to_value(&self.0).map_err(|e| JsError::new(&e.to_string()))
    }

    #[wasm_bindgen(js_name = "fromJSON")]
    pub fn from_json(val: JsValue) -> Result<Self, JsError> {
        let inner = serde_wasm_bindgen::from_value(val).map_err(|e| JsError::new(&e.to_string()))?;
        Ok(Self(inner))
    }

    pub fn copy(&self) -> Self {
        *self
    }
}

#[wasm_bindgen(inspectable, js_name = "SimplePeak")]
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize)]
pub struct SimpleWebPeak {
    pub mz: f64,
    pub intensity: f32,
}

impl From<MZPoint> for SimpleWebPeak {
    fn from(value: MZPoint) -> Self {
        (&value).into()
    }
}

impl From<&MZPoint> for SimpleWebPeak {
    fn from(value: &MZPoint) -> Self {
        SimpleWebPeak {
            mz: value.mz,
            intensity: value.intensity,
            ..Default::default()
        }
    }
}

impl From<&CentroidPeak> for SimpleWebPeak {
    fn from(value: &CentroidPeak) -> Self {
        SimpleWebPeak {
            mz: value.mz,
            intensity: value.intensity as f32,
            ..Default::default()
        }
    }
}
impl From<&DeconvolutedPeak> for SimpleWebPeak {
    fn from(value: &DeconvolutedPeak) -> Self {
        SimpleWebPeak {
            mz: value.mz(),
            intensity: value.intensity as f32,
            ..Default::default()
        }
    }
}

#[wasm_bindgen(inspectable, js_name = "SimpleChargedPeak")]
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize)]
pub struct SimpleWebChargedPeak {
    pub mz: f64,
    pub intensity: f32,
    pub charge: i32,
    pub score: f32,
}

impl From<&DeconvolutedPeak> for SimpleWebChargedPeak {
    fn from(value: &DeconvolutedPeak) -> Self {
        SimpleWebChargedPeak {
            mz: value.mz(),
            intensity: value.intensity as f32,
            charge: value.charge,
            ..Default::default()
        }
    }
}

#[wasm_bindgen(js_name = "Param")]
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize)]
pub struct WebParam(Param);


fn value_to_js(val: &ValueRef<'_>) -> JsValue {
    match val {
        ValueRef::String(cow) => JsValue::from_str(&cow),
        ValueRef::Float(x) => JsValue::from_f64(*x),
        ValueRef::Int(x) => JsValue::from_f64(*x as f64),
        ValueRef::Buffer(cow) => {
            let val: Array = cow
                .to_vec()
                .into_iter()
                .map(|x| JsValue::from_f64(x as f64))
                .collect();
            val.into()
        }
        ValueRef::Empty => JsValue::null(),
        ValueRef::Boolean(x) => JsValue::from_bool(*x),
        ValueRef::List(v) => {
            let jsv: Vec<_> = v.iter().map(|vi| value_to_js(&vi.as_ref())).collect();
            jsv.into()
        }
    }
}


#[wasm_bindgen(js_class = "Param")]
impl WebParam {
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.0.name().to_string()
    }

    #[wasm_bindgen(getter)]
    pub fn id(&self) -> Option<String> {
        self.0.curie_str()
    }

    #[wasm_bindgen(getter)]
    pub fn value(&self) -> JsValue {
        let val = self.0.value();
        value_to_js(&val)
    }

    #[wasm_bindgen(getter)]
    pub fn unit(&self) -> String {
        self.0.unit.to_string()
    }

    #[wasm_bindgen(js_name = "toJSON")]
    pub fn to_json(&self) -> Result<JsValue, JsError> {
        serde_wasm_bindgen::to_value(&self.0).map_err(|e| JsError::new(&e.to_string()))
    }

    #[wasm_bindgen(js_name = "fromJSON")]
    pub fn from_json(val: JsValue) -> Result<Self, JsError> {
        let inner = serde_wasm_bindgen::from_value(val).map_err(|e| JsError::new(&e.to_string()))?;
        Ok(Self(inner))
    }
    #[wasm_bindgen(js_name = "toString")]
    pub fn to_string(&self) -> String {
        format!(
            "WebParam {{ name: {}, value: {}, id: {:?} }} ",
            self.name(),
            self.0.value.to_string(),
            self.id()
        )
    }
}

#[wasm_bindgen(js_name = "IsolationWindow")]
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize)]
pub struct WebIsolationWindow(IsolationWindow);

#[wasm_bindgen(js_class = "IsolationWindow")]
impl WebIsolationWindow {
    #[wasm_bindgen(getter, js_name = "lowerBound")]
    pub fn lower_bound(&self) -> f32 {
        self.0.lower_bound
    }

    #[wasm_bindgen(getter, js_name = "upperBound")]
    pub fn upper_bound(&self) -> f32 {
        self.0.upper_bound
    }

    pub fn contains(&self, x: f32) -> bool {
        self.0.contains(x)
    }

    #[wasm_bindgen(js_name = "toJSON")]
    pub fn to_json(&self) -> Object {
        let inst = Object::new();
        Reflect::set(
            &inst,
            &JsValue::from_str("lower_bound"),
            &JsValue::from_f64(self.0.lower_bound as f64),
        )
        .unwrap();
        Reflect::set(
            &inst,
            &JsValue::from_str("upper_bound"),
            &JsValue::from_f64(self.0.upper_bound as f64),
        )
        .unwrap();
        inst
    }
}

#[wasm_bindgen(js_name = "ScanWindow")]
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize)]
pub struct WebScanWindow(ScanWindow);

#[wasm_bindgen(js_class = "ScanWindow")]
impl WebScanWindow {
    #[wasm_bindgen(getter, js_name = "lowerBound")]
    pub fn lower_bound(&self) -> f32 {
        self.0.lower_bound
    }

    #[wasm_bindgen(getter, js_name = "upperBound")]
    pub fn upper_bound(&self) -> f32 {
        self.0.upper_bound
    }

    pub fn contains(&self, x: f32) -> bool {
        self.0.contains(x)
    }

    #[wasm_bindgen(js_name = "toJSON")]
    pub fn to_json(&self) -> Object {
        let inst = Object::new();
        Reflect::set(
            &inst,
            &JsValue::from_str("lower_bound"),
            &JsValue::from_f64(self.0.lower_bound as f64),
        )
        .unwrap();
        Reflect::set(
            &inst,
            &JsValue::from_str("upper_bound"),
            &JsValue::from_f64(self.0.upper_bound as f64),
        )
        .unwrap();
        inst
    }
}

#[wasm_bindgen(js_name = "ScanEvent")]
pub struct WebScanEvent(ScanEvent);

#[wasm_bindgen(js_class = "ScanEvent")]
impl WebScanEvent {
    #[wasm_bindgen(getter, js_name = "startTime")]
    pub fn start_time(&self) -> f64 {
        self.0.start_time
    }

    #[wasm_bindgen(getter, js_name = "injectionTime")]
    pub fn injection_time(&self) -> f32 {
        self.0.injection_time
    }

    #[wasm_bindgen(getter, js_name = "scanWindows")]
    pub fn scan_windows(&self) -> Vec<WebScanWindow> {
        self.0
            .scan_windows
            .iter()
            .map(|w| WebScanWindow(w.clone()))
            .collect()
    }

    #[wasm_bindgen(getter, js_name = "instrumentConfigurationID")]
    pub fn instrument_configuration_id(&self) -> u32 {
        self.0.instrument_configuration_id
    }

    #[wasm_bindgen(js_name = "filterString")]
    pub fn filter_string(&self) -> Option<String> {
        self.0.filter_string().map(|s| s.to_string())
    }

    pub fn params(&self) -> Vec<WebParam> {
        self.0
            .params()
            .iter()
            .map(|p| WebParam(p.clone()))
            .collect()
    }
}

#[wasm_bindgen(js_name = "SelectedIon")]
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct WebSelectedIon(SelectedIon);

#[wasm_bindgen(js_class = "SelectedIon")]
impl WebSelectedIon {
    #[wasm_bindgen(getter)]
    pub fn mz(&self) -> f64 {
        self.0.mz
    }

    #[wasm_bindgen(getter)]
    pub fn intensity(&self) -> f32 {
        self.0.intensity
    }

    #[wasm_bindgen(getter)]
    pub fn charge(&self) -> Option<i32> {
        self.0.charge
    }

    pub fn params(&self) -> Vec<WebParam> {
        self.0
            .params()
            .iter()
            .cloned()
            .map(|p| WebParam(p))
            .collect()
    }
}

#[wasm_bindgen(js_name = "Activation")]
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct WebActivation(Activation);

#[wasm_bindgen(js_class = "Activation")]
impl WebActivation {
    #[wasm_bindgen(getter)]
    pub fn method(&self) -> Option<String> {
        self.0.method().map(|m| m.name().to_string())
    }

    #[wasm_bindgen(getter)]
    pub fn energy(&self) -> f32 {
        self.0.energy
    }

    pub fn methods(&self) -> Vec<WebParam> {
        self.0
            .methods()
            .iter()
            .map(|m| WebParam(m.to_param().into()))
            .collect()
    }

    pub fn params(&self) -> Vec<WebParam> {
        self.0.params.iter().cloned().map(|p| WebParam(p)).collect()
    }

    #[wasm_bindgen(js_name = "toJSON")]
    pub fn to_json(&self) -> Result<Object, JsValue> {
        let entries = Array::of3(
            &Array::of2(&JsValue::from_str("method"), &self.method().into()),
            &Array::of2(&JsValue::from_str("energy"), &self.energy().into()),
            &Array::of2(&JsValue::from_str("params"), &self.params().into()),
        );
        Object::from_entries(&entries)
    }
}

#[wasm_bindgen(getter_with_clone, inspectable, js_name = "Precursor")]
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct WebPrecursor {
    _inner: Precursor,
    pub ions: Vec<WebSelectedIon>,
    #[wasm_bindgen(js_name = "isolationWindow")]
    pub isolation_window: WebIsolationWindow,
    pub activation: WebActivation,
}

#[wasm_bindgen(js_class = "Precursor")]
impl WebPrecursor {
    fn new(precursor: Precursor) -> Self {
        let ions = precursor
            .ions
            .iter()
            .cloned()
            .map(|i| WebSelectedIon(i))
            .collect();
        let isolation_window = WebIsolationWindow(precursor.isolation_window.clone());
        let activation = WebActivation(precursor.activation.clone());
        Self {
            _inner: precursor,
            ions,
            isolation_window,
            activation,
        }
    }

    #[wasm_bindgen(getter, js_name = "precursorScanID")]
    pub fn precursor_id(&self) -> Option<String> {
        self._inner.precursor_id.clone()
    }
}

#[wasm_bindgen(js_name = "SignalContinuity")]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default, Hash, Eq)]
pub enum WebSignalContinuity {
    #[default]
    Unknown,
    Centroid,
    Profile,
}

impl From<SignalContinuity> for WebSignalContinuity {
    fn from(value: SignalContinuity) -> Self {
        match value {
            SignalContinuity::Unknown => Self::Unknown,
            SignalContinuity::Centroid => Self::Centroid,
            SignalContinuity::Profile => Self::Profile,
        }
    }
}

#[wasm_bindgen(js_name = "Spectrum")]
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct WebSpectrum {
    inner: MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak>,
}

impl From<MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak>> for WebSpectrum {
    fn from(mut value: MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak>) -> Self {
        value.arrays.as_ref().inspect(|a| {
            for (_k, v) in a.iter() {
                if let Some(n) = v.data_len().ok() {
                    if n > 0 {
                        // info!("Converting {} with {} items", k, n);
                        break;
                    }
                }
            }
        });
        match value.try_build_peaks() {
            Ok(_x) => {
                // if x.len() > 0 {
                //     match x {
                //         RefPeakDataLevel::Missing => info!("Missing"),
                //         RefPeakDataLevel::RawData(_) => info!("RawData"),
                //         RefPeakDataLevel::Centroid(_) => info!("Centroid"),
                //         RefPeakDataLevel::Deconvoluted(_) => info!("Deconvoluted"),
                //     }
                // }
            }
            Err(e) => {
                info!(
                    "Failed to locate data arrays for spectrum {}: {e}",
                    value.id()
                )
            }
        }

        Self { inner: value }
    }
}

impl WebSpectrum {
    pub fn as_ref(&self) -> &MultiLayerSpectrum<CentroidPeak, DeconvolutedPeak> {
        &self.inner
    }

    pub fn new(description: SpectrumDescription) -> Self {
        Self {
            inner: MultiLayerSpectrum::from_peaks_data_levels_and_description(
                PeakDataLevel::Missing,
                description,
            ),
        }
    }

    pub fn new_with_peaks(
        description: SpectrumDescription,
        peaks: PeakDataLevel<CentroidPeak, DeconvolutedPeak>,
    ) -> Self {
        Self {
            inner: MultiLayerSpectrum::from_peaks_data_levels_and_description(peaks, description),
        }
    }

    pub fn description(&self) -> &SpectrumDescription {
        self.inner.description()
    }

    pub fn description_mut(&mut self) -> &mut SpectrumDescription {
        self.inner.description_mut()
    }

    pub fn peaks(
        &self,
    ) -> mzdata::spectrum::RefPeakDataLevel<'_, CentroidPeak, DeconvolutedPeak> {
        self.inner.peaks()
    }
}

#[wasm_bindgen(js_class = "Spectrum")]
impl WebSpectrum {
    #[wasm_bindgen(js_name = "toJSON")]
    pub fn to_json(&self) -> Result<JsValue, JsError> {
        serde_wasm_bindgen::to_value(&self.inner).map_err(|e| JsError::new(&e.to_string()))
    }

    #[wasm_bindgen(js_name = "fromJSON")]
    pub fn from_json(val: JsValue) -> Result<Self, JsError> {
        let inner = serde_wasm_bindgen::from_value(val).map_err(|e| JsError::new(&e.to_string()))?;
        Ok(Self { inner })
    }

    #[wasm_bindgen(js_name = "fromBatchJSON")]
    pub fn batch_from_json(text: String) -> Result<Vec<Self>, JsError> {
        let inner = serde_json::from_str::<Vec<_>>(&text)
            .map_err(|e| JsError::new(&e.to_string()))?
            .into_iter()
            .map(|inner| Self { inner })
            .collect();
        Ok(inner)
    }

    #[wasm_bindgen(getter)]
    pub fn id(&self) -> String {
        self.description().id.clone()
    }

    #[wasm_bindgen(getter, js_name = "startTime")]
    pub fn start_time(&self) -> f64 {
        self.description().acquisition.start_time()
    }

    #[wasm_bindgen(getter, js_name = "signalContinuity")]
    pub fn signal_continuity(&self) -> WebSignalContinuity {
        self.description().signal_continuity.into()
    }

    #[wasm_bindgen(getter, js_name = "isProfile")]
    pub fn is_profile(&self) -> bool {
        matches!(
            self.description().signal_continuity,
            SignalContinuity::Profile
        )
    }

    #[wasm_bindgen(getter)]
    pub fn index(&self) -> usize {
        self.description().index
    }

    #[wasm_bindgen(getter, js_name = "msLevel")]
    pub fn ms_level(&self) -> u8 {
        self.description().ms_level
    }

    #[wasm_bindgen(getter)]
    pub fn precursor(&self) -> Option<WebPrecursor> {
        self.inner.precursor()
            .cloned()
            .map(|p| WebPrecursor::new(p))
    }

    pub fn params(&self) -> Vec<WebParam> {
        self.description()
            .params
            .iter()
            .cloned()
            .map(|p| WebParam(p))
            .collect()
    }

    #[wasm_bindgen(getter, js_name = "scanEvents")]
    pub fn scan_events(&self) -> Vec<WebScanEvent> {
        self.description()
            .acquisition
            .iter()
            .map(|e| WebScanEvent(e.clone()))
            .collect()
    }

    pub fn copy(&self) -> Self {
        Self {
            inner: self.inner.clone(),
        }
    }

    #[wasm_bindgen(js_name = "pickPeaks")]
    pub fn pick_peaks(&mut self, signal_to_noise_threshold: f32) {
        if !self.is_profile() && self.inner.raw_arrays().is_none() {
            return;
        }
        self.inner.pick_peaks(signal_to_noise_threshold).unwrap();
    }

    pub fn reprofile(&mut self, dx: f64, fwhm: f32) {
        if self.is_profile() {
            return;
        }
        self.inner.reprofile_with_shape(dx, fwhm).unwrap();
        self.description_mut().signal_continuity = SignalContinuity::Profile;
    }

    pub fn denoise(&mut self, scale: f32) {
        if self.inner.raw_arrays().is_none() {
            log::warn!("Cannot denoise a spectrum that has no profile signal")
        }
        self.inner.denoise(scale).unwrap();
    }

    #[wasm_bindgen(getter)]
    pub fn length(&self) -> usize {
        self.inner.peaks().len()
    }

    #[wasm_bindgen(js_name = "hasPeak")]
    pub fn has_peak(&self, query: f64, error_tolerance: &WebTolerance) -> Option<SimpleWebPeak> {
        self.peaks()
            .search(query, error_tolerance.0)
            .and_then(|i| self.peaks().get(i))
            .map(|p| SimpleWebPeak {
                mz: p.mz,
                intensity: p.intensity,
            })
    }

    #[wasm_bindgen(js_name = "allPeaksFor")]
    pub fn all_peaks_for(&self, query: f64, error_tolerance: &WebTolerance) -> Vec<SimpleWebPeak> {
        match self.peaks() {
            RefPeakDataLevel::Missing => Vec::new(),
            RefPeakDataLevel::RawData(_) => Vec::new(),
            RefPeakDataLevel::Centroid(peak_set_vec) => peak_set_vec
                .all_peaks_for(query, error_tolerance.0)
                .into_iter()
                .map(|p| p.into())
                .collect(),
            RefPeakDataLevel::Deconvoluted(peak_set_vec) => peak_set_vec
                .all_peaks_for(query, error_tolerance.0)
                .into_iter()
                .map(|p| p.into())
                .collect(),
        }
    }

    pub fn at(&self, index: usize) -> Option<SimpleWebPeak> {
        self.peaks().get(index).map(|p| SimpleWebPeak {
            mz: p.mz,
            intensity: p.intensity,
        })
    }

    #[wasm_bindgen(js_name = "basePeak")]
    pub fn base_peak(&self) -> SimpleWebPeak {
        let p = self.peaks().base_peak();
        SimpleWebPeak::from(&p)
    }

    pub fn tic(&self) -> f32 {
        self.peaks().tic()
    }

    pub fn between(&self, low: f64, high: f64) -> Vec<SimpleWebPeak> {
        match &self.peaks() {
            RefPeakDataLevel::Missing => Vec::new(),
            RefPeakDataLevel::RawData(_) => Vec::new(),
            RefPeakDataLevel::Centroid(peak_set_vec) => peak_set_vec
                .between(low, high, Tolerance::PPM(20.0))
                .iter()
                .map(|p| p.into())
                .collect(),
            RefPeakDataLevel::Deconvoluted(peak_set_vec) => peak_set_vec
                .between(low, high, Tolerance::PPM(20.0))
                .iter()
                .map(|p| p.into())
                .collect(),
        }
    }

    /// Convert whatever the most processed data are into [`SimplePeak`]s
    #[wasm_bindgen(js_name = "toArray")]
    pub fn to_array(&self) -> Vec<SimpleWebPeak> {
        self.peaks().iter().map(|p| p.into()).collect()
    }

    #[wasm_bindgen(js_name = "rawArrays")]
    pub fn raw_arrays(&self) -> Option<Object> {
        self.inner
            .arrays
            .as_ref()
            .map(array_map_to_js)
            .into_iter()
            .flatten()
            .next()
    }

    #[wasm_bindgen(js_name = "centroidPeaks")]
    pub fn centroid_peaks(&self) -> Option<Vec<SimpleWebPeak>> {
        self.inner
            .peaks
            .as_ref()
            .map(|peaks| peaks.iter().map(|p| p.into()).collect())
    }

    #[wasm_bindgen(js_name = "deconvolutedPeaks")]
    pub fn deconvoluted_peaks(&self) -> Option<Vec<SimpleWebChargedPeak>> {
        self.inner
            .deconvoluted_peaks
            .as_ref()
            .map(|peaks| peaks.iter().map(|p| p.into()).collect())
    }

    #[wasm_bindgen(js_name = "hasIonMobilityDimension")]
    pub fn has_ion_mobility_dimension(&self) -> bool {
        self.inner.has_ion_mobility_dimension()
    }

    #[wasm_bindgen(js_name = "asIonMobilityFrame")]
    pub fn as_ion_mobility_frame(&self) -> Result<WebIonMobilityFrame, String> {
        if self.has_ion_mobility_dimension() {
            let arrays = self
                .inner
                .arrays
                .as_ref()
                .ok_or_else(|| "No data arrays found".to_string())?;
            let arrays = arrays
                .clone()
                .stack_ion_mobility()
                .map_err(|e| e.to_string())?;

            let frame = MultiLayerIonMobilityFrame::new(
                Some(arrays),
                None,
                None,
                self.description().clone().into(),
            );

            let this = WebIonMobilityFrame { inner: frame };

            Ok(this)
        } else {
            Err("No ion mobility dimension found".into())
        }
    }
}

#[wasm_bindgen(inspectable, js_name = "FeaturePoint")]
pub struct WebFeaturePoint {
    pub mz: f64,
    pub time: f64,
    pub intensity: f32,
}

impl From<(f64, f64, f32)> for WebFeaturePoint {
    fn from(value: (f64, f64, f32)) -> Self {
        Self {
            mz: value.0,
            time: value.1,
            intensity: value.2,
        }
    }
}

#[wasm_bindgen(js_name = "Feature")]
#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct WebFeature(Feature<MZ, IonMobility>);

impl From<Feature<MZ, Time>> for WebFeature {
    fn from(value: Feature<MZ, Time>) -> Self {
        let (x, y, z) = value.into_inner();
        Self(Feature::new(x, y, z))
    }
}

impl From<Feature<MZ, IonMobility>> for WebFeature {
    fn from(value: Feature<MZ, IonMobility>) -> Self {
        let (x, y, z) = value.into_inner();
        Self(Feature::new(x, y, z))
    }
}

impl From<SimpleFeature<MZ, IonMobility>> for WebFeature {
    fn from(value: SimpleFeature<MZ, IonMobility>) -> Self {
        let (x, y, z) = value.into_inner();
        let x: Vec<_> = std::iter::repeat_n(x, y.len()).collect();
        Self(Feature::new(x, y, z))
    }
}

#[wasm_bindgen(js_class = "Feature")]
impl WebFeature {
    #[wasm_bindgen(js_name = "toJSON")]
    pub fn to_json(&self) -> Result<JsValue, JsError> {
        serde_wasm_bindgen::to_value(&self.0).map_err(|e| JsError::new(&e.to_string()))
    }

    #[wasm_bindgen(js_name = "fromJSON")]
    pub fn from_json(text: JsValue) -> Result<Self, JsError> {
        let inner = serde_wasm_bindgen::from_value(text).map_err(|e| JsError::new(&e.to_string()))?;
        Ok(Self(inner))
    }
    #[wasm_bindgen(constructor)]
    pub fn new(x: Vec<f64>, y: Vec<f64>, z: Vec<f32>) -> Self {
        Self(Feature::new(x, y, z))
    }

    pub fn clone(&self) -> Self {
        Self(self.0.clone())
    }

    #[wasm_bindgen(getter, js_name = "startTime")]
    pub fn start_time(&self) -> Option<f64> {
        self.0.start_time()
    }

    #[wasm_bindgen(getter, js_name = "endTime")]
    pub fn end_time(&self) -> Option<f64> {
        self.0.end_time()
    }

    #[wasm_bindgen(getter, js_name = "apexTime")]
    pub fn apex_time(&self) -> Option<f64> {
        self.0.apex_time()
    }

    #[wasm_bindgen(getter)]
    pub fn times(&self) -> Box<[f64]> {
        self.0.as_view().into_inner().1.into()
    }

    #[wasm_bindgen(getter)]
    pub fn intensities(&self) -> Box<[f32]> {
        self.0.as_view().into_inner().2.into()
    }

    #[wasm_bindgen(getter)]
    pub fn mzs(&self) -> Box<[f64]> {
        self.0.as_view().into_inner().0.into()
    }

    pub fn at(&self, index: usize) -> Option<WebFeaturePoint> {
        self.0.at(index).map(|pt| pt.into())
    }

    #[wasm_bindgen(getter)]
    pub fn length(&self) -> usize {
        self.0.len()
    }

    pub fn smooth(&mut self, size: usize) {
        self.0.smooth(size);
    }

    #[wasm_bindgen(js_name = "totalIonCurrent", getter)]
    pub fn tic(&self) -> f32 {
        self.0.total_intensity()
    }

    #[wasm_bindgen(js_name = "averageMz", getter)]
    pub fn average_mz(&self) -> f64 {
        self.0.mz()
    }

    #[wasm_bindgen(js_name = "fitPeaks")]
    pub fn fit_peaks(&self) -> FeatureFit {
        FeatureFit(self.0.fit_peaks_with(Default::default()).peak_fits)
    }

    pub fn area(&self) -> f32 {
        self.0.area()
    }
}

#[wasm_bindgen(js_name = "DeconvolvedFeature")]
#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct WebDeconvolvedFeature(ChargedFeature<Mass, IonMobility>);

impl From<ChargedFeature<Mass, IonMobility>> for WebDeconvolvedFeature {
    fn from(value: ChargedFeature<Mass, IonMobility>) -> Self {
        Self(value)
    }
}

#[wasm_bindgen(js_class = "DeconvolvedFeature")]
impl WebDeconvolvedFeature {
    pub fn clone(&self) -> Self {
        Self(self.0.clone())
    }

    #[wasm_bindgen(js_name = "toJSON")]
    pub fn to_json(&self) -> Result<JsValue, JsError> {
        serde_wasm_bindgen::to_value(&self.0).map_err(|e| JsError::new(&e.to_string()))
    }

    #[wasm_bindgen(js_name = "fromJSON")]
    pub fn from_json(text: JsValue) -> Result<Self, JsError> {
        let inner = serde_wasm_bindgen::from_value(text).map_err(|e| JsError::new(&e.to_string()))?;
        Ok(Self(inner))
    }

    #[wasm_bindgen(getter, js_name = "startTime")]
    pub fn start_time(&self) -> Option<f64> {
        self.0.start_time()
    }

    #[wasm_bindgen(getter, js_name = "endTime")]
    pub fn end_time(&self) -> Option<f64> {
        self.0.end_time()
    }

    #[wasm_bindgen(getter, js_name = "apexTime")]
    pub fn apex_time(&self) -> Option<f64> {
        self.0.apex_time()
    }

    #[wasm_bindgen(getter)]
    pub fn charge(&self) -> i32 {
        self.0.charge()
    }

    #[wasm_bindgen(getter)]
    pub fn times(&self) -> Box<[f64]> {
        self.0
            .as_inner().0
            .as_view()
            .into_inner()
            .1
            .into()
    }

    #[wasm_bindgen(getter)]
    pub fn intensities(&self) -> Box<[f32]> {
        self.0
            .as_inner().0
            .as_view()
            .into_inner()
            .2
            .into()
    }

    #[wasm_bindgen(getter)]
    pub fn mzs(&self) -> Box<[f64]> {
        let z = self.charge();
        let view = self.0.as_inner().0.as_view().into_inner().0;
        view.iter()
            .copied()
            .map(|mass| mass_charge_ratio(mass, z))
            .collect()
    }

    pub fn at(&self, index: usize) -> Option<WebFeaturePoint> {
        self.0.at(index).map(|pt| pt.into())
    }

    #[wasm_bindgen(getter)]
    pub fn length(&self) -> usize {
        self.0.len()
    }

    pub fn smooth(&mut self, size: usize) {
        self.0.smooth(size);
    }

    #[wasm_bindgen(js_name = "totalIonCurrent", getter)]
    pub fn tic(&self) -> f32 {
        self.0.intensity()
    }

    #[wasm_bindgen(js_name = "weightedNeutralMass", getter)]
    pub fn neutral_mass(&self) -> f64 {
        self.0.neutral_mass()
    }

    #[wasm_bindgen(js_name = "weightedMZ", getter)]
    pub fn mz(&self) -> f64 {
        self.0.mz()
    }

    #[wasm_bindgen(js_name = "fitPeaks")]
    pub fn fit_peaks(&self) -> FeatureFit {
        FeatureFit(
            self.0
                .as_inner().0
                .fit_peaks_with(Default::default())
                .peak_fits,
        )
    }

    pub fn area(&self) -> f32 {
        self.0.area()
    }
}

#[wasm_bindgen]
pub struct FeatureFit(MultiPeakShapeFit);

#[wasm_bindgen]
impl FeatureFit {
    pub fn models(&self) -> Vec<JsValue> {
        self.0
            .iter()
            .map(|f| serde_wasm_bindgen::to_value(f).unwrap())
            .collect()
    }

    pub fn predict(&self, times: &[f64]) -> Vec<f64> {
        self.0.predict(times).to_vec()
    }

    pub fn density(&self, time: f64) -> f64 {
        self.0.density(time)
    }

    #[wasm_bindgen(js_name = "toJSON")]
    pub fn to_json(&self) -> JsValue {
        self.models().into()
    }
}

#[wasm_bindgen(js_name = "IonMobilityFrame")]
#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct WebIonMobilityFrame {
    inner: MultiLayerIonMobilityFrame<
        Feature<MZ, IonMobility>,
        ChargedFeature<Mass, IonMobility>,
    >,
}

impl WebIonMobilityFrame {
    fn description(&self) -> &mzdata::spectrum::IonMobilityFrameDescription {
        self.inner.description()
    }
}

impl
    From<
        MultiLayerIonMobilityFrame<
            Feature<MZ, IonMobility>,
            ChargedFeature<Mass, IonMobility>,
        >,
    > for WebIonMobilityFrame
{
    fn from(
        mut value: MultiLayerIonMobilityFrame<
            Feature<MZ, IonMobility>,
            ChargedFeature<Mass, IonMobility>,
        >,
    ) -> Self {
        value.try_build_features().ok();
        Self { inner: value }
    }
}

#[wasm_bindgen(js_class = "IonMobilityFrame")]
impl WebIonMobilityFrame {
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> String {
        self.description().id.clone()
    }

    #[wasm_bindgen(js_name = "toJSON")]
    pub fn to_json(&self) -> Result<JsValue, JsError> {
        serde_wasm_bindgen::to_value(&self.inner).map_err(|e| JsError::new(&e.to_string()))
    }

    #[wasm_bindgen(js_name = "fromJSON")]
    pub fn from_json(val: JsValue) -> Result<Self, JsError> {
        let inner = serde_wasm_bindgen::from_value(val).map_err(|e| JsError::new(&e.to_string()))?;
        Ok(Self { inner })
    }

    #[wasm_bindgen(js_name = "fromBatchJSON")]
    pub fn batch_from_json(text: String) -> Result<Vec<Self>, JsError> {
        let inner = serde_json::from_str::<Vec<_>>(&text)
            .map_err(|e| JsError::new(&e.to_string()))?
            .into_iter()
            .map(|inner| Self { inner })
            .collect();
        Ok(inner)
    }

    #[wasm_bindgen(getter, js_name = "startTime")]
    pub fn start_time(&self) -> f64 {
        self.description().acquisition.start_time()
    }

    #[wasm_bindgen(getter, js_name = "signalContinuity")]
    pub fn signal_continuity(&self) -> WebSignalContinuity {
        self.description().signal_continuity.into()
    }

    #[wasm_bindgen(getter, js_name = "isProfile")]
    pub fn is_profile(&self) -> bool {
        matches!(
            self.description().signal_continuity,
            SignalContinuity::Profile
        )
    }

    #[wasm_bindgen(getter)]
    pub fn index(&self) -> usize {
        self.description().index
    }

    #[wasm_bindgen(getter, js_name = "msLevel")]
    pub fn ms_level(&self) -> u8 {
        self.description().ms_level
    }

    #[wasm_bindgen(getter)]
    pub fn precursor(&self) -> Option<WebPrecursor> {
        self.inner.precursor()
            .cloned()
            .map(|p| WebPrecursor::new(p))
    }

    pub fn params(&self) -> Vec<WebParam> {
        self.description()
            .params
            .iter()
            .cloned()
            .map(|p| WebParam(p))
            .collect()
    }

    #[wasm_bindgen(getter, js_name = "scanEvents")]
    pub fn scan_events(&self) -> Vec<WebScanEvent> {
        self.description()
            .acquisition
            .iter()
            .map(|e| WebScanEvent(e.clone()))
            .collect()
    }

    pub fn copy(&self) -> Self {
        Self {
            inner: self.inner.clone(),
        }
    }

    #[wasm_bindgen(getter)]
    pub fn length(&self) -> usize {
        match self.inner.features() {
            mzdata::spectrum::RefFeatureDataLevel::Missing => 0,
            mzdata::spectrum::RefFeatureDataLevel::RawData(binary_array_map3_d) => {
                binary_array_map3_d.iter().map(|(_t, v)| v.len()).sum()
            }
            mzdata::spectrum::RefFeatureDataLevel::Centroid(feature_map) => feature_map.len(),
            mzdata::spectrum::RefFeatureDataLevel::Deconvoluted(feature_map) => feature_map.len(),
        }
    }

    pub fn at(&self, index: usize) -> Option<WebFeature> {
        match self.inner.features() {
            mzdata::spectrum::RefFeatureDataLevel::Missing => None,
            mzdata::spectrum::RefFeatureDataLevel::RawData(_) => None,
            mzdata::spectrum::RefFeatureDataLevel::Centroid(feature_map) => feature_map
                .as_slice()
                .get(index)
                .map(|f| WebFeature(f.clone())),
            mzdata::spectrum::RefFeatureDataLevel::Deconvoluted(feature_map) => {
                feature_map.as_slice().get(index).map(|f| {
                    let (inner, _z) = f.clone().into_inner();
                    let (x, y, z) = inner.into_inner();
                    WebFeature::new(x, y, z)
                })
            }
        }
    }

    #[wasm_bindgen(js_name = "extractFeatures")]
    pub fn extract_features(
        &mut self,
        min_length: usize,
        maximum_gap_size: f64,
        error_tolerance: Option<WebTolerance>,
    ) {
        let error_tolerance = error_tolerance.map(|e| e.0).unwrap_or(Tolerance::PPM(15.0));
        self.inner
            .extract_features_simple(error_tolerance, min_length, maximum_gap_size, None)
            .unwrap();
    }

    pub fn features(&self) -> Option<Vec<WebFeature>> {
        self.inner
            .features
            .as_ref()
            .map(|fs| fs.iter().map(|f| WebFeature(f.clone())).collect())
    }

    #[wasm_bindgen(js_name = "deconvolutedFeatures")]
    pub fn deconvoluted_features(&self) -> Option<Vec<WebDeconvolvedFeature>> {
        self.inner.deconvoluted_features.as_ref().map(|fs| {
            fs.iter()
                .map(|f| WebDeconvolvedFeature(f.clone()))
                .collect()
        })
    }

    #[wasm_bindgen(js_name = "rawArrays")]
    pub fn raw_arrays(&self) -> Result<Vec<Object>, JsError> {
        let mut points = Vec::new();
        if let Some(arrays) = self.inner.arrays.as_ref() {
            for (im, arrays) in arrays.iter() {
                let mzs = arrays.mzs()?;
                let intens = arrays.intensities()?;
                let mzs_js = Float64Array::new_with_length(mzs.len() as u32);
                mzs_js.copy_from(&mzs);
                let intens_js = Float32Array::new_with_length(mzs.len() as u32);
                intens_js.copy_from(&intens);
                let arrays = Object::new();
                Reflect::set(&arrays, &JsValue::from_str("m/z array"), &mzs_js)
                    .map_err(|e| JsError::new(&e.as_string().unwrap_or_default()))?;
                Reflect::set(&arrays, &JsValue::from_str("intensity array"), &intens_js)
                    .map_err(|e| JsError::new(&e.as_string().unwrap_or_default()))?;
                Reflect::set(
                    &arrays,
                    &JsValue::from_str("ion mobility"),
                    &JsValue::from_f64(im),
                )
                .map_err(|e| JsError::new(&e.as_string().unwrap_or_default()))?;
                points.push(arrays);
            }
        }

        Ok(points)
    }

    #[wasm_bindgen(js_name = "clearArrays")]
    pub fn clear_arrays(&mut self) {
        self.inner.arrays = None;
    }

    #[wasm_bindgen(js_name = "clearFeatures")]
    pub fn clear_features(&mut self) {
        self.inner.features = None;
    }

    #[wasm_bindgen(js_name = "clearDeconvolutedFeatures")]
    pub fn clear_deconvoluted_features(&mut self) {
        self.inner.deconvoluted_features = None;
    }

    #[wasm_bindgen(js_name = "hasFeature")]
    pub fn has_peak(&self, query: f64, error_tolerance: WebTolerance) -> Option<WebFeature> {
        self.inner.features.as_ref().and_then(|fm| {
            fm.has_feature(query, error_tolerance.0)
                .map(|f| WebFeature(f.clone()))
        })
    }

    #[wasm_bindgen(js_name = "allPeaksFor")]
    pub fn all_peaks_for(&self, query: f64, error_tolerance: WebTolerance) -> Vec<WebFeature> {
        self.inner
            .features
            .as_ref()
            .map(|fm| {
                fm.all_features_for(query, error_tolerance.0)
                    .into_iter()
                    .map(|f| WebFeature(f.clone()))
                    .collect()
            })
            .unwrap_or_default()
    }
}
