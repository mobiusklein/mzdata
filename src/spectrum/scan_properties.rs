use super::spectrum::{CentroidPeakAdapting, DeconvolutedPeakAdapting, SpectrumBehavior};
use crate::io::traits::ScanSource;
use crate::params::ControlledVocabulary;
use crate::{impl_param_described, ParamList};
use crate::{params, Param};

/**
Describe the initialization stage of an isolation window
*/
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
    pub target: f32,
    pub lower_bound: f32,
    pub upper_bound: f32,
    /// Describes the decision making process used to establish the bounds of the
    /// window from the source file.
    pub flags: IsolationWindowState,
}

impl PartialEq for IsolationWindow {
    fn eq(&self, other: &Self) -> bool {
        self.target == other.target
            && self.lower_bound == other.lower_bound
            && self.upper_bound == other.upper_bound
    }
}

#[derive(Default, Debug, Clone, PartialEq)]
pub struct ScanWindow {
    pub lower_bound: f32,
    pub upper_bound: f32,
}

pub type ScanWindowList = Vec<ScanWindow>;

#[derive(Default, Debug, Clone, PartialEq)]
/// Describes a single scan event. Unless additional post-processing is done,
/// there is usually only one event per spectrum.
pub struct ScanEvent {
    pub start_time: f64,
    pub injection_time: f32,
    pub scan_windows: ScanWindowList,
    pub instrument_configuration_id: u32,
    pub params: ParamList,
}

pub type ScanEventList = Vec<ScanEvent>;

#[derive(Default, Debug, Clone, PartialEq)]
/// Describe the series of acquisition events that constructed the spectrum
/// being described.
pub struct Acquisition {
    pub scans: ScanEventList,
    pub params: ParamList,
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

    pub fn instrument_configuration_ids(&self) -> Vec<u32> {
        self.scans.iter().map(|s| {
            s.instrument_configuration_id
        }).collect()
    }
}

pub trait IonProperties {
    fn neutral_mass(&self) -> f64;
    fn charge(&self) -> Option<i32>;
    fn has_charge(&self) -> bool {
        self.charge().is_some()
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
/// Describes a single selected ion from a precursor isolation
pub struct SelectedIon {
    /// The selected ion's m/z as reported, may not be the monoisotopic peak.
    pub mz: f64,
    pub intensity: f32,
    /// The reported precursor ion's charge state. May be absent in
    /// some source files.
    pub charge: Option<i32>,
    pub params: ParamList,
}

impl IonProperties for SelectedIon {
    fn neutral_mass(&self) -> f64 {
        let charge = match self.charge {
            Some(z) => z,
            None => 1,
        };
        crate::utils::neutral_mass(self.mz, charge)
    }

    fn charge(&self) -> Option<i32> {
        self.charge
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
/// Describes the activation method used to dissociate the precursor ion
pub struct Activation {
    _method: Option<Param>,
    pub energy: f32,
    pub params: ParamList,
}

impl Activation {
    pub fn method(&self) -> Option<&Param> {
        if self._method.is_some() {
            self._method.as_ref()
        } else {
            None
        }
    }

    pub fn method_mut(&mut self) -> &mut Option<Param> {
        &mut self._method
    }

    pub fn is_param_activation(p: &Param) -> bool {
        if p.is_controlled() && p.controlled_vocabulary.unwrap() == ControlledVocabulary::MS {
            Self::accession_to_activation(p.accession.unwrap())
        } else {
            false
        }
    }

    pub fn accession_to_activation(accession: u32) -> bool {
        match accession {
            1000133 => true,
            1000134 => true,
            1000135 => true,
            1000136 => true,
            1000242 => true,
            1000250 => true,
            1000282 => true,
            1000433 => true,
            1000435 => true,
            1000598 => true,
            1000599 => true,
            1001880 => true,
            1002000 => true,
            1003181 => true,
            1003247 => true,
            1000422 => true,
            1002472 => true,
            1002679 => true,
            1003294 => true,
            1000262 => true,
            1003246 => true,
            1002631 => true,
            1003182 => true,
            1002481 => true,
            1002678 => true,
            _ => false,
        }
    }

    pub fn _set_method(&mut self) {
        let found = self
            .params
            .iter()
            .enumerate()
            .filter(|(_, p)| {
                if p.is_controlled() && p.controlled_vocabulary.unwrap() == ControlledVocabulary::MS
                {
                    Self::accession_to_activation(p.accession.unwrap())
                } else {
                    false
                }
            })
            .map(|(i, _)| i)
            .next();
        if let Some(hit) = found {
            self._method = Some(self.params.remove(hit));
        }
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
/// Describes the precursor ion of the owning spectrum.
pub struct Precursor {
    /// Describes the selected ion's properties
    pub ion: SelectedIon,
    /// Describes the isolation window around the selected ion
    pub isolation_window: IsolationWindow,
    /// The precursor scan ID, if given
    pub precursor_id: Option<String>,
    /// The product scan ID, if given
    pub product_id: Option<String>,
    /// The activation process applied to the precursor ion
    pub activation: Activation,
    /// Additional parameters describing this precursor ion
    pub params: ParamList,
}

impl Precursor {
    /// Given a ScanSource object, look up the precursor scan in it.
    /// This is useful when examining the area *around* where the precursor
    /// ion was or to obtain a snapshot of the retention time when the spectrum
    /// was scheduled.
    pub fn precursor_spectrum<C, D, S, R>(&self, source: &mut R) -> Option<S>
    where
        C: CentroidPeakAdapting,
        D: DeconvolutedPeakAdapting,
        S: SpectrumBehavior<C, D>,
        R: ScanSource<C, D, S>,
    {
        match self.precursor_id.as_ref() {
            Some(id) => source.get_spectrum_by_id(id),
            None => None,
        }
    }

    /// Given a ScanSource object, look up the product scan in it.
    /// This is rarely needed unless you have manually separated [`Precursor`]
    /// objects from their spectra.
    pub fn product_spectrum<C, D, S, R>(&self, source: &mut R) -> Option<S>
    where
        C: CentroidPeakAdapting,
        D: DeconvolutedPeakAdapting,
        S: SpectrumBehavior<C, D>,
        R: ScanSource<C, D, S>,
    {
        match self.product_id.as_ref() {
            Some(id) => source.get_spectrum_by_id(id),
            None => None,
        }
    }
}

/**
A trait for abstracting over how a precursor ion is described, immutably.
*/
pub trait PrecursorSelection: params::ParamDescribed {
    /// Describes the selected ion's properties
    fn ion(&self) -> &SelectedIon;
    /// Describes the isolation window around the selected ion
    fn isolation_window(&self) -> &IsolationWindow;
    /// The precursor scan ID, if given
    fn precursor_id(&self) -> Option<&String>;
    /// The product scan ID, if given
    fn product_id(&self) -> Option<&String>;
    /// The activation process applied to the precursor ion
    fn activation(&self) -> &Activation;
}

impl PrecursorSelection for Precursor {
    fn ion(&self) -> &SelectedIon {
        &self.ion
    }

    fn isolation_window(&self) -> &IsolationWindow {
        &self.isolation_window
    }

    fn precursor_id(&self) -> Option<&String> {
        self.precursor_id.as_ref()
    }

    fn product_id(&self) -> Option<&String> {
        self.product_id.as_ref()
    }

    fn activation(&self) -> &Activation {
        &self.activation
    }
}

impl<T> IonProperties for T
where
    T: PrecursorSelection,
{
    fn neutral_mass(&self) -> f64 {
        self.ion().neutral_mass()
    }

    fn charge(&self) -> Option<i32> {
        self.ion().charge()
    }
}

/**
Describes the polarity of a mass spectrum. A spectrum is either `Positive` (1+), `Negative` (-1)
or `Unknown` (0). The `Unknown` state is the default.
*/
#[repr(i8)]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq)]
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

/**
Describes the initial representation of the signal of a spectrum.

Though most formats explicitly have a method of either conveying a processing level
or an assumed level, the `Unknown` option is retained for partial initialization.
*/
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum SignalContinuity {
    Unknown = 0,
    Centroid = 3,
    Profile = 5,
}

impl Default for SignalContinuity {
    fn default() -> SignalContinuity {
        SignalContinuity::Unknown
    }
}

/**
The set of descriptive metadata that give context for how a mass spectrum was acquired
within a particular run. This forms the basis for a large portion of the [`SpectrumBehavior`]
trait.
*/
#[derive(Debug, Default, Clone, PartialEq)]
pub struct SpectrumDescription {
    pub id: String,
    pub index: usize,
    pub ms_level: u8,

    pub polarity: ScanPolarity,
    pub signal_continuity: SignalContinuity,

    pub params: ParamList,
    pub acquisition: Acquisition,
    pub precursor: Option<Precursor>,
}

impl_param_described!(
    Acquisition,
    Activation,
    Precursor,
    SelectedIon,
    ScanEvent,
    SpectrumDescription
);
