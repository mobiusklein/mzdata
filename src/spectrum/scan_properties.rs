use super::spectrum::{CentroidPeakAdapting, DeconvolutedPeakAdapting, SpectrumBehavior};
use crate::io::traits::ScanSource;
use crate::params::{ControlledVocabulary, Param, ParamLike, Unit};
use crate::{impl_param_described, ParamList};

/**
Describe the initialization stage of an isolation window
*/
#[derive(Debug, Clone, Copy, Default)]
#[repr(i8)]
pub enum IsolationWindowState {
    #[default]
    Unknown = 0,
    Offset,
    Explicit,
    Complete,
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
    pub params: Option<Box<ParamList>>,
}

pub type ScanEventList = Vec<ScanEvent>;

#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq, Default)]
pub enum ScanCombination {
    // MS:1000795
    #[default]
    NoCombination,
    // MS:1000571
    Sum,
    // MS:1000573
    Median,
}


impl ScanCombination {
    pub fn from_accession(
        controlled_vocabulary: ControlledVocabulary,
        accession: u32,
    ) -> Option<ScanCombination> {
        match controlled_vocabulary {
            ControlledVocabulary::MS => match accession {
                1000795 => Some(Self::NoCombination),
                1000571 => Some(Self::Sum),
                1000573 => Some(Self::Median),
                _ => None,
            },
            _ => None,
        }
    }

    pub fn name(&self) -> &str {
        match self {
            ScanCombination::NoCombination => "no combination",
            ScanCombination::Sum => "sum of spectra",
            ScanCombination::Median => "median of spectra",
        }
    }

    pub fn accession(&self) -> u32 {
        match self {
            ScanCombination::NoCombination => 1000795,
            ScanCombination::Sum => 1000571,
            ScanCombination::Median => 1000573,
        }
    }

    pub fn to_param(&self) -> Param {
        Param {
            name: self.name().to_string(),
            value: String::default(),
            accession: Some(self.accession()),
            controlled_vocabulary: Some(ControlledVocabulary::MS),
            unit: Unit::Unknown,
        }
    }
}

#[derive(Default, Debug, Clone, PartialEq)]
/// Describe the series of acquisition events that constructed the spectrum
/// being described.
pub struct Acquisition {
    pub scans: ScanEventList,
    pub combination: ScanCombination,
    pub params: Option<Box<ParamList>>,
}

impl Acquisition {

    pub fn start_time(&self) -> f64 {
        self.first_scan().unwrap().start_time
    }

    pub fn first_scan(&self) -> Option<&ScanEvent> {
        self.scans.first()
    }

    pub fn first_scan_mut(&mut self) -> Option<&mut ScanEvent> {
        if self.scans.is_empty() {
            self.scans.push(ScanEvent::default());
        }
        self.scans.first_mut()
    }

    pub fn last_scan(&self) -> Option<&ScanEvent> {
        self.scans.last()
    }

    pub fn last_scan_mut(&mut self) -> Option<&mut ScanEvent> {
        if self.scans.is_empty() {
            self.scans.push(ScanEvent::default());
        }
        self.scans.last_mut()
    }

    pub fn instrument_configuration_ids(&self) -> Vec<u32> {
        self.scans
            .iter()
            .map(|s| s.instrument_configuration_id)
            .collect()
    }
}

pub trait IonProperties {
    fn mz(&self) -> f64;
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
    pub params: Option<Box<ParamList>>,
}

impl IonProperties for SelectedIon {
    #[inline]
    fn mz(&self) -> f64 {
        self.mz
    }

    #[inline]
    fn neutral_mass(&self) -> f64 {
        let charge = match self.charge {
            Some(z) => z,
            None => 1,
        };
        crate::utils::neutral_mass(self.mz, charge)
    }

    #[inline]
    fn charge(&self) -> Option<i32> {
        self.charge
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum ActivationMethod {
    CollisionInducedDissociation,
    HighEnergyCollisionInducedDissociation,
    LowEnergyCollisionInducedDissociation,
    InSourceCollisionInducedDissociation,
    BeamTypeCollisionInducedDissociation,
    TrapTypeCollisionInducedDissociation,

    SupplementalCollisionInducedDissociation,
    SupplementalBeamTypeCollisionInducedDissociation,

    ElectronTransferDissociation,
    ElectronCaptureDissociation,
    ElectronActivationDissociation,
    NegativeElectronTransferDissociation,

    Photodissociation,
    UltravioletPhotodissociation,
    Other(Box<Param>),
}

impl ActivationMethod {
    pub fn accession(&self) -> Option<u32> {
        match self {
            ActivationMethod::CollisionInducedDissociation => Some(1000133),
            ActivationMethod::LowEnergyCollisionInducedDissociation => Some(1000433),
            ActivationMethod::BeamTypeCollisionInducedDissociation => Some(1000422),
            ActivationMethod::TrapTypeCollisionInducedDissociation => Some(1002472),
            ActivationMethod::SupplementalCollisionInducedDissociation => Some(1002678),
            ActivationMethod::SupplementalBeamTypeCollisionInducedDissociation => Some(1002679),
            ActivationMethod::HighEnergyCollisionInducedDissociation => Some(1002481),
            ActivationMethod::InSourceCollisionInducedDissociation => Some(1001880),

            ActivationMethod::ElectronCaptureDissociation => Some(1000250),
            ActivationMethod::ElectronTransferDissociation => Some(1000598),
            ActivationMethod::NegativeElectronTransferDissociation => Some(1003247),
            ActivationMethod::ElectronActivationDissociation => Some(1003294),

            ActivationMethod::Photodissociation => Some(1000435),
            ActivationMethod::UltravioletPhotodissociation => Some(1003246),

            ActivationMethod::Other(param) => param.accession(),
        }
    }

    pub fn is_collisional(&self) -> Option<bool> {
        match self {
            ActivationMethod::CollisionInducedDissociation | ActivationMethod::LowEnergyCollisionInducedDissociation |
            ActivationMethod::BeamTypeCollisionInducedDissociation | Self::TrapTypeCollisionInducedDissociation |
            Self::InSourceCollisionInducedDissociation | Self::SupplementalBeamTypeCollisionInducedDissociation |
            Self::SupplementalCollisionInducedDissociation | Self::HighEnergyCollisionInducedDissociation => Some(true),
            Self::Other(_) => None,
            _ => Some(false)
        }
    }

    pub fn controlled_vocabulary(&self) -> Option<ControlledVocabulary> {
        match self {
            ActivationMethod::Other(param) => param.controlled_vocabulary(),
            _ => Some(ControlledVocabulary::MS),
        }
    }
}

impl From<ActivationMethod> for Param {
    fn from(value: ActivationMethod) -> Self {
        match value {
            ActivationMethod::Other(val) => *val,
            ActivationMethod::CollisionInducedDissociation => ControlledVocabulary::MS
                .const_param_ident("collision-induced dissociation", value.accession().unwrap())
                .into(),
            ActivationMethod::LowEnergyCollisionInducedDissociation => ControlledVocabulary::MS
                .const_param_ident(
                    "low-energy collision-induced dissociation",
                    value.accession().unwrap(),
                )
                .into(),
            ActivationMethod::BeamTypeCollisionInducedDissociation => ControlledVocabulary::MS
                .const_param_ident(
                    "beam-type collision-induced dissociation",
                    value.accession().unwrap(),
                )
                .into(),
            ActivationMethod::TrapTypeCollisionInducedDissociation => ControlledVocabulary::MS
                .const_param_ident(
                    "trap-type collision-induced dissociation",
                    value.accession().unwrap(),
                )
                .into(),
            ActivationMethod::HighEnergyCollisionInducedDissociation => ControlledVocabulary::MS
                .const_param_ident(
                    "higher energy collision-induced dissociation",
                    value.accession().unwrap(),
                )
                .into(),
            ActivationMethod::InSourceCollisionInducedDissociation => ControlledVocabulary::MS
                .const_param_ident(
                    "in-source collision-induced dissociation",
                    value.accession().unwrap(),
                )
                .into(),
            ActivationMethod::SupplementalCollisionInducedDissociation => ControlledVocabulary::MS
                .const_param_ident(
                    "supplemental collision-induced dissociation",
                    value.accession().unwrap(),
                )
                .into(),
            ActivationMethod::SupplementalBeamTypeCollisionInducedDissociation => {
                ControlledVocabulary::MS
                    .const_param_ident(
                        "supplemental beam-type collision-induced dissociation",
                        value.accession().unwrap(),
                    )
                    .into()
            }
            ActivationMethod::ElectronTransferDissociation => ControlledVocabulary::MS
                .const_param_ident("electron transfer dissociation", value.accession().unwrap())
                .into(),
            ActivationMethod::ElectronCaptureDissociation => ControlledVocabulary::MS
                .const_param_ident("electron capture dissociation", value.accession().unwrap())
                .into(),
            ActivationMethod::ElectronActivationDissociation => ControlledVocabulary::MS
                .const_param_ident(
                    "electron activation dissociation",
                    value.accession().unwrap(),
                )
                .into(),
            ActivationMethod::NegativeElectronTransferDissociation => ControlledVocabulary::MS
                .const_param_ident(
                    "negative electron transfer dissociation",
                    value.accession().unwrap(),
                )
                .into(),
            ActivationMethod::Photodissociation => ControlledVocabulary::MS
                .const_param_ident("photodissociation", value.accession().unwrap())
                .into(),
            ActivationMethod::UltravioletPhotodissociation => ControlledVocabulary::MS
                .const_param_ident("ultraviolet photodissociation", value.accession().unwrap())
                .into(),
        }
    }
}

impl<P: ParamLike + Into<Param>> From<P> for ActivationMethod {
    fn from(value: P) -> Self {
        if value.is_ms() {
            match value.accession().unwrap() {
                1000133 => Self::CollisionInducedDissociation,
                1000433 => Self::LowEnergyCollisionInducedDissociation,
                1002481 => Self::HighEnergyCollisionInducedDissociation,
                1000422 => Self::BeamTypeCollisionInducedDissociation,
                1002472 => Self::TrapTypeCollisionInducedDissociation,
                1002678 => Self::SupplementalCollisionInducedDissociation,
                1002679 => Self::SupplementalBeamTypeCollisionInducedDissociation,

                1001880 => Self::InSourceCollisionInducedDissociation,

                1000250 => Self::ElectronCaptureDissociation,
                1000598 => Self::ElectronTransferDissociation,
                1003247 => Self::NegativeElectronTransferDissociation,
                1003294 => Self::ElectronActivationDissociation,

                1000435 => Self::Photodissociation,
                1003246 => Self::UltravioletPhotodissociation,

                _ => Self::Other(Box::new(value.into())),
            }
        } else {
            Self::Other(Box::new(value.into()))
        }
    }
}


#[derive(Debug, Default, Clone, PartialEq)]
/// Describes the activation method used to dissociate the precursor ion
pub struct Activation {
    _method: Option<ActivationMethod>,
    pub energy: f32,
    pub params: ParamList,
}

impl Activation {
    pub fn method(&self) -> Option<&ActivationMethod> {
        if self._method.is_some() {
            self._method.as_ref()
        } else {
            None
        }
    }

    pub fn method_mut(&mut self) -> &mut Option<ActivationMethod> {
        &mut self._method
    }

    pub fn is_param_activation<P: ParamLike>(p: &P) -> bool {
        if p.is_controlled() && p.controlled_vocabulary().unwrap() == ControlledVocabulary::MS {
            Self::accession_to_activation(p.accession().unwrap())
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
            self._method = Some(self.params.remove(hit).into());
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
pub trait PrecursorSelection {
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

    #[inline]
    fn mz(&self) -> f64 {
        self.ion().mz()
    }

    #[inline]
    fn neutral_mass(&self) -> f64 {
        self.ion().neutral_mass()
    }

    #[inline]
    fn charge(&self) -> Option<i32> {
        self.ion().charge()
    }
}

/**
Describes the polarity of a mass spectrum. A spectrum is either `Positive` (1+), `Negative` (-1)
or `Unknown` (0). The `Unknown` state is the default.
*/
#[repr(i8)]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Default)]
pub enum ScanPolarity {
    #[default]
    Unknown = 0,
    Positive = 1,
    Negative = -1,
}


/**
Describes the initial representation of the signal of a spectrum.

Though most formats explicitly have a method of either conveying a processing level
or an assumed level, the `Unknown` option is retained for partial initialization.
*/
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
pub enum SignalContinuity {
    #[default]
    Unknown = 0,
    Centroid = 3,
    Profile = 5,
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

impl_param_described!(Activation, SpectrumDescription);
impl_param_described_deferred!(SelectedIon, Acquisition, ScanEvent);


#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum ChromatogramType {
    #[default]
    Unknown,
    TotalIonCurrentChromatogram,
    BasePeakChromatogram,
    SelectedIonCurrentChromatogram,
    SelectedReactionMonitoringChromatogram,
}


#[derive(Debug, Default, Clone, PartialEq)]
pub struct ChromatogramDescription {
    pub id: String,
    pub index: usize,
    pub ms_level: Option<u8>,
    pub polarity: ScanPolarity,
    pub chromatogram_type: ChromatogramType,

    pub params: ParamList,
    pub precursor: Option<Precursor>,
}

impl_param_described!(ChromatogramDescription);