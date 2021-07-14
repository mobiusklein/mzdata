use std::collections::{HashMap};

use super::params;

#[derive(Debug, Clone, Copy)]
pub enum IsolationWindowState {
    Unknown = 0,
    Offset,
    Explicit,
    Complete
}

impl Default for IsolationWindowState {
    fn default() -> IsolationWindowState {
        IsolationWindowState::Unknown
    }
}


#[derive(Default, Debug, Clone)]
pub struct IsolationWindow {
    pub target: f64,
    pub lower_bound: f64,
    pub upper_bound: f64,
    pub flags: IsolationWindowState
}

#[derive(Default, Debug, Clone)]
pub struct ScanWindow {
    pub lower_bound: f64,
    pub upper_bound: f64
}


pub type ScanWindowList = Vec<ScanWindow>;


#[derive(Default, Debug, Clone)]
pub struct ScanEventInformation {
    pub start_time: f64,
    pub injection_time: f32,
    pub scan_windows: ScanWindowList,
    pub params: params::ParamList
}


pub type ScanEventList = Vec<ScanEventInformation>;


#[derive(Default, Debug, Clone)]
pub struct AcquisitionInformation {
    pub scans: ScanEventList,
    pub params: params::ParamList
}


#[derive(Debug, Clone, Default)]
pub struct SelectedIon {
    pub mz: f64,
    pub intensity: f32,
    pub charge: i32,
    pub params: params::ParamList
}

#[derive(Debug, Default, Clone)]
pub struct Activation {
    pub method: String,
    pub energy: f32,
    pub params: params::ParamList
}

#[derive(Debug, Default, Clone)]
pub struct Precursor {
    pub ion: SelectedIon,
    pub isolation_window: IsolationWindow,
    pub precursor_id: String,
    pub product_id: String,
    pub params: params::ParamList,
    pub activation: Activation
}

#[derive(Debug, Default, Clone)]
pub struct SpectrumDescription {
    pub id: String,
    pub index: usize,
    pub ms_level: i8,
    pub time: f32,
    pub polarity: i8,
    pub is_profile: i8,
    pub params: params::ParamList,
    pub acquisition: Option<AcquisitionInformation>,
    pub precursor_information: Option<Precursor>,
    pub annotations: HashMap<String, String>,
}
