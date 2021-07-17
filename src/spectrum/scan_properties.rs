use std::collections::HashMap;

use super::params;

#[derive(Debug, Clone, Copy)]
pub enum IsolationWindowState {
    Unknown = 0,
    Offset,
    Explicit,
    Complete,
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
    pub flags: IsolationWindowState,
}

#[derive(Default, Debug, Clone)]
pub struct ScanWindow {
    pub lower_bound: f64,
    pub upper_bound: f64,
}

pub type ScanWindowList = Vec<ScanWindow>;

#[derive(Default, Debug, Clone)]
pub struct ScanEvent {
    pub start_time: f64,
    pub injection_time: f32,
    pub scan_windows: ScanWindowList,
    pub instrument_configuration_id: String,
    pub params: params::ParamList,
}

pub type ScanEventList = Vec<ScanEvent>;

#[derive(Default, Debug, Clone)]
pub struct Acquisition {
    pub scans: ScanEventList,
    pub params: params::ParamList,
}

impl Acquisition {
    pub fn first_scan(&self) -> Option<&ScanEvent> {
        self.scans.first()
    }

    pub fn first_scan_mut(&mut self) -> Option<&mut ScanEvent> {
        if self.scans.is_empty() {
            self.scans.push(ScanEvent::default());
        }
        self.scans.first_mut()
    }
}

#[derive(Debug, Clone, Default)]
pub struct SelectedIon {
    pub mz: f64,
    pub intensity: f32,
    pub charge: i32,
    pub params: params::ParamList,
}

#[derive(Debug, Default, Clone)]
pub struct Activation {
    pub method: String,
    pub energy: f32,
    pub params: params::ParamList,
}

#[derive(Debug, Default, Clone)]
pub struct Precursor {
    pub ion: SelectedIon,
    pub isolation_window: IsolationWindow,
    pub precursor_id: String,
    pub product_id: String,
    pub params: params::ParamList,
    pub activation: Activation,
}

#[repr(i8)]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum ScanPolarity {
    Unknown = 0,
    Positive = 1,
    Negative = -1,
}

impl Default for ScanPolarity {
    fn default() -> ScanPolarity {
        ScanPolarity::Unknown
    }
}

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum ScanSiganlContinuity {
    Unknown = 0,
    Centroid = 1,
    Profile = 5,
}

impl Default for ScanSiganlContinuity {
    fn default() -> ScanSiganlContinuity {
        ScanSiganlContinuity::Unknown
    }
}

#[derive(Debug, Default, Clone)]
pub struct SpectrumDescription {
    pub id: String,
    pub index: usize,
    pub ms_level: u8,

    pub polarity: ScanPolarity,
    pub is_profile: ScanSiganlContinuity,

    pub params: params::ParamList,
    pub acquisition: Acquisition,
    pub precursor: Option<Precursor>,
    pub annotations: HashMap<String, String>,
}
