use std::collections::HashMap;

use crate::params;

#[derive(Debug, Clone, Copy)]
#[repr(i8)]
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
/// The interval around the precursor ion that was isolated in the precursor scan.
/// Although an isolation window may be specified either with explicit bounds or
/// offsets from the target, this data structure always uses explicit bounds.
pub struct IsolationWindow {
    pub target: f64,
    pub lower_bound: f64,
    pub upper_bound: f64,
    /// Describes the decision making process used to establish the bounds of the
    /// window from the source file.
    pub flags: IsolationWindowState,
}

#[derive(Default, Debug, Clone)]
pub struct ScanWindow {
    pub lower_bound: f64,
    pub upper_bound: f64,
}

pub type ScanWindowList = Vec<ScanWindow>;

#[derive(Default, Debug, Clone)]
/// Describes a single scan event. Unless additional post-processing is done,
/// there is usually only one event per spectrum.
pub struct ScanEvent {
    pub start_time: f64,
    pub injection_time: f32,
    pub scan_windows: ScanWindowList,
    pub instrument_configuration_id: String,
    pub params: params::ParamList,
}

pub type ScanEventList = Vec<ScanEvent>;

#[derive(Default, Debug, Clone)]
/// Describe the series of acquisition events that constructed the spectrum
/// being described.
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
/// Describes a single selected ion from a precursor isolation
pub struct SelectedIon {
    /// The selected ion's m/z as reported, may not be the monoisotopic peak.
    pub mz: f64,
    pub intensity: f32,
    /// The reported precursor ion's charge state. May be absent in
    /// some source files.
    pub charge: Option<i32>,
    pub params: params::ParamList,
}

#[derive(Debug, Default, Clone)]
/// Describes the activation method used to dissociate the precursor ion
pub struct Activation {
    pub method: String,
    pub energy: f32,
    pub params: params::ParamList,
}

#[derive(Debug, Default, Clone)]
/// Describes the precursor ion of the owning spectrum.
pub struct Precursor {
    /// Describes the selected ion's properties
    pub ion: SelectedIon,
    /// Describes the isolation window around the selected ion
    pub isolation_window: IsolationWindow,
    /// The precursor scan ID, if given
    pub precursor_id: String,
    /// The product scan ID, if given
    pub product_id: String,
    /// The activation process applied to the precursor ion
    pub activation: Activation,
    /// Additional parameters describing this precursor ion
    pub params: params::ParamList,
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
pub enum SignalContinuity {
    Unknown = 0,
    Centroid = 1,
    Profile = 5,
}

impl Default for SignalContinuity {
    fn default() -> SignalContinuity {
        SignalContinuity::Unknown
    }
}

#[derive(Debug, Default, Clone)]
pub struct SpectrumDescription {
    pub id: String,
    pub index: usize,
    pub ms_level: u8,

    pub polarity: ScanPolarity,
    pub signal_continuity: SignalContinuity,

    pub params: params::ParamList,
    pub acquisition: Acquisition,
    pub precursor: Option<Precursor>,
    pub annotations: HashMap<String, String>,
}
