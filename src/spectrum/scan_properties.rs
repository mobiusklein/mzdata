use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fmt::Display;

use log::warn;
use num_traits::Float;

use super::spectrum_types::{CentroidPeakAdapting, DeconvolutedPeakAdapting, SpectrumLike};
use crate::io::traits::SpectrumSource;
use crate::meta::DissociationMethodTerm;
use crate::params::{
    AccessionIntCode, ControlledVocabulary, Param, ParamDescribed, ParamLike, ParamValue, Unit,
    CURIE,
};
use crate::{curie, impl_param_described, ParamList};

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

impl Display for IsolationWindowState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[derive(Default, Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
/// The interval around the precursor ion that was isolated in the precursor scan.
/// Although an isolation window may be specified either with explicit bounds or
/// offsets from the target, this data structure always uses explicit bounds once
/// it is in a [`IsolationWindowState::Complete`] .
pub struct IsolationWindow {
    /// The recorded isolation window target m/z, which may actually be outside the window
    pub target: f32,
    /// The lower m/z boundary of the isolation window if `flags` is
    /// [`IsolationWindowState::Explicit`], or the offset from `target`
    /// if `flags` is [`IsolationWindowState::Offset`]
    pub lower_bound: f32,
    /// The upper m/z boundary of the isolation window if `flags` is
    /// [`IsolationWindowState::Explicit`], or the offset from `target`
    /// if `flags` is [`IsolationWindowState::Offset`]
    pub upper_bound: f32,
    /// Describes the decision making process used to establish the bounds of the
    /// window from the source file.
    #[cfg_attr(feature = "serde", serde(skip))]
    pub flags: IsolationWindowState,
}

impl IsolationWindow {
    pub fn new(
        target: f32,
        lower_bound: f32,
        upper_bound: f32,
        flags: IsolationWindowState,
    ) -> Self {
        Self {
            target,
            lower_bound,
            upper_bound,
            flags,
        }
    }

    pub fn around(target: f32, width: f32) -> Self {
        let lower_bound = target - width;
        let upper_bound = target + width;
        Self::new(
            target,
            lower_bound,
            upper_bound,
            IsolationWindowState::Complete,
        )
    }

    pub fn contains<F: Float>(&self, point: F) -> bool {
        let point = point.to_f32().unwrap();
        self.lower_bound <= point && self.upper_bound <= point
    }

    pub fn is_empty(&self) -> bool {
        self.lower_bound == 0.0 && self.upper_bound == 0.0
    }
}

impl PartialEq for IsolationWindow {
    fn eq(&self, other: &Self) -> bool {
        self.target == other.target
            && self.lower_bound == other.lower_bound
            && self.upper_bound == other.upper_bound
    }
}

/// The m/z range which was scanned
#[derive(Default, Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ScanWindow {
    /// The minimum m/z scanned
    pub lower_bound: f32,
    /// The maximum m/z scanned
    pub upper_bound: f32,
}

impl ScanWindow {
    pub fn new(lower_bound: f32, upper_bound: f32) -> Self {
        Self {
            lower_bound,
            upper_bound,
        }
    }

    pub fn contains<F: Float>(&self, point: F) -> bool {
        let point = point.to_f32().unwrap();
        self.lower_bound <= point && self.upper_bound <= point
    }

    pub fn is_empty(&self) -> bool {
        self.lower_bound == 0.0 && self.upper_bound == 0.0
    }
}

type ScanWindowList = Vec<ScanWindow>;

#[derive(Default, Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
/// Describes a single scan event. Unless additional post-processing is done,
/// there is usually only one event per spectrum.
pub struct ScanEvent {
    /// When did the scan start acquisition relative to the run start, in minutes
    pub start_time: f64,
    /// The time spent filling the ion trapping device, in milliseconds
    pub injection_time: f32,
    /// The windows of the primary axis covered by this event
    pub scan_windows: ScanWindowList,
    /// The identifier of the [`InstrumentConfiguration`](crate::meta::InstrumentConfiguration) associated
    /// with the event
    pub instrument_configuration_id: u32,
    pub params: Option<Box<ParamList>>,
    pub spectrum_reference: Option<Box<str>>,
}

pub(crate) const ION_MOBILITY_SCAN_TERMS: [CURIE; 4] = [
    // ion mobility drift time
    curie!(MS:1002476),
    // inverse reduced ion mobility drift time
    curie!(MS:1002815),
    // FAIMS compensation voltage
    curie!(MS:1001581),
    // SELEXION compensation voltage
    curie!(MS:1003371),
];

pub trait IonMobilityMeasure: ParamDescribed {
    fn ion_mobility(&'_ self) -> Option<f64> {
        for u in ION_MOBILITY_SCAN_TERMS {
            if let Some(v) = self.get_param_by_curie(&u).map(|p| p.value()) {
                return v.to_f64().map(Some).unwrap_or_else(|e| {
                    warn!("Failed to parse ion mobility {u} value {v}: {e}");
                    None
                });
            }
        }
        None
    }

    fn has_ion_mobility(&self) -> bool {
        self.ion_mobility().is_some()
    }

    fn ion_mobility_type(&self) -> Option<&Param> {
        for u in ION_MOBILITY_SCAN_TERMS {
            if let Some(v) = self.get_param_by_curie(&u) {
                return Some(v);
            }
        }
        None
    }
}

pub(crate) const PRESET_SCAN_CONFIGURATION: CURIE = curie!(MS:1000616);
pub(crate) const MASS_RESOLUTION: CURIE = curie!(MS:1000011);
pub(crate) const FILTER_STRING: CURIE = curie!(MS:1000512);
pub(crate) const SCAN_TITLE: CURIE = curie!(MS:1000499);

impl ScanEvent {
    pub fn new(
        start_time: f64,
        injection_time: f32,
        scan_windows: ScanWindowList,
        instrument_configuration_id: u32,
        params: Option<Box<ParamList>>,
    ) -> Self {
        Self {
            start_time,
            injection_time,
            scan_windows,
            instrument_configuration_id,
            params,
            spectrum_reference: None,
        }
    }

    crate::find_param_method!(
        filter_string,
        &FILTER_STRING,
        |p| { p.as_str() },
        Option<Cow<'_, str>>
    );
    crate::find_param_method!(resolution, &MASS_RESOLUTION);
    crate::find_param_method!(scan_configuration, &PRESET_SCAN_CONFIGURATION);
}

impl IonMobilityMeasure for ScanEvent {}

type ScanEventList = Vec<ScanEvent>;

/// Represents means by which a spectrum is generated using
/// one or more instrument analyzers
#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ScanCombination {
    // MS:1000795
    #[default]
    NoCombination,
    // MS:1000571
    Sum,
    // MS:1000573
    Median,
}

impl Display for ScanCombination {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl ScanCombination {
    pub fn from_accession(
        controlled_vocabulary: ControlledVocabulary,
        accession: AccessionIntCode,
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

    pub const fn name(&self) -> &str {
        match self {
            ScanCombination::NoCombination => "no combination",
            ScanCombination::Sum => "sum of spectra",
            ScanCombination::Median => "median of spectra",
        }
    }

    pub const fn accession(&self) -> AccessionIntCode {
        match self {
            ScanCombination::NoCombination => 1000795,
            ScanCombination::Sum => 1000571,
            ScanCombination::Median => 1000573,
        }
    }

    pub fn to_param(&self) -> Param {
        Param {
            name: self.name().to_string(),
            value: Default::default(),
            accession: Some(self.accession()),
            controlled_vocabulary: Some(ControlledVocabulary::MS),
            unit: Unit::Unknown,
        }
    }
}

#[derive(Default, Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
/// Describe the series of acquisition events that constructed the spectrum
/// being described.
pub struct Acquisition {
    pub scans: ScanEventList,
    pub combination: ScanCombination,
    pub params: Option<Box<ParamList>>,
}

impl Acquisition {
    pub fn start_time(&self) -> f64 {
        self.first_scan().map(|v| v.start_time).unwrap_or_default()
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

    pub fn len(&self) -> usize {
        self.scans.len()
    }

    pub fn is_empty(&self) -> bool {
        self.scans.is_empty()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, ScanEvent> {
        self.scans.iter()
    }

    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, ScanEvent> {
        self.scans.iter_mut()
    }
}

/// Describe the precursor ion that this entity represents
pub trait IonProperties {
    /// The selected ion m/z
    fn mz(&self) -> f64;
    /// The selected ion's estimated neutral mass, given the m/z and charge
    fn neutral_mass(&self) -> f64;
    /// The expected charge state of the selected ion, if known
    fn charge(&self) -> Option<i32>;
    fn has_charge(&self) -> bool {
        self.charge().is_some()
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
        crate::utils::neutral_mass(self.mz, self.charge.unwrap_or(1))
    }

    #[inline]
    fn charge(&self) -> Option<i32> {
        self.charge
    }
}

impl IonMobilityMeasure for SelectedIon {}

#[derive(Debug, Default, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
/// Describes the activation method used to dissociate the precursor ion
pub struct Activation {
    _methods: Vec<DissociationMethodTerm>,
    pub energy: f32,
    pub params: ParamList,
}

impl Activation {
    /// Get a reference to the first activation method, if it exists
    pub fn method(&self) -> Option<&DissociationMethodTerm> {
        self._methods.first()
    }

    /// Get a mutable reference to the first activation method, if it exists
    pub fn method_mut(&mut self) -> Option<&mut DissociationMethodTerm> {
        self._methods.first_mut()
    }

    /// Get a slice over all of the activation methods used.
    pub fn methods(&self) -> &[DissociationMethodTerm] {
        &self._methods
    }

    /// Get a mutable reference to a [`Vec`] of activation methods, which may be useful
    /// for adding or removing methods.
    pub fn methods_mut(&mut self) -> &mut Vec<DissociationMethodTerm> {
        &mut self._methods
    }

    /// Check if multiple dissociation methods were used
    pub fn is_combined(&self) -> bool {
        self._methods.len() > 1
    }

    /// Check if a [`ParamLike`] type references an activation method
    pub fn is_param_activation<P: ParamLike>(p: &P) -> bool {
        if p.is_controlled() && p.controlled_vocabulary().unwrap() == ControlledVocabulary::MS {
            Self::accession_to_activation(p.accession().unwrap())
        } else {
            false
        }
    }

    pub fn accession_to_activation(accession: AccessionIntCode) -> bool {
        DissociationMethodTerm::from_accession(accession).is_some()
    }

    pub fn _extract_methods_from_params(&mut self) {
        let mut methods = Vec::with_capacity(1);
        let mut rest = Vec::with_capacity(self.params.len());
        for p in self.params.drain(..) {
            if Self::is_param_activation(&p) {
                methods.push(p.into())
            } else {
                rest.push(p)
            }
        }
        self.params = rest;
        self._methods = methods;
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
/// Describes the precursor ion of the owning spectrum.
pub struct Precursor {
    /// Describes the selected ion's properties
    pub ions: Vec<SelectedIon>,
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
    /// Given a [`SpectrumSource`] object, look up the precursor scan in it.
    /// This is useful when examining the area *around* where the precursor
    /// ion was or to obtain a snapshot of the retention time when the spectrum
    /// was scheduled.
    pub fn precursor_spectrum<C, D, S, R>(&self, source: &mut R) -> Option<S>
    where
        C: CentroidPeakAdapting,
        D: DeconvolutedPeakAdapting,
        S: SpectrumLike<C, D>,
        R: SpectrumSource<C, D, S>,
    {
        match self.precursor_id.as_ref() {
            Some(id) => source.get_spectrum_by_id(id),
            None => None,
        }
    }

    /// Check if this precursor selection selected any ions
    pub fn has_ions(&self) -> bool {
        self.ions.is_empty()
    }

    /// Given a [`SpectrumSource`] object, look up the product scan in it.
    /// This is rarely needed unless you have manually separated [`Precursor`]
    /// objects from their spectra.
    pub fn product_spectrum<C, D, S, R>(&self, source: &mut R) -> Option<S>
    where
        C: CentroidPeakAdapting,
        D: DeconvolutedPeakAdapting,
        S: SpectrumLike<C, D>,
        R: SpectrumSource<C, D, S>,
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
    fn ion(&self) -> Option<&SelectedIon>;
    /// Describes the isolation window around the selected ion
    fn isolation_window(&self) -> &IsolationWindow;
    /// The precursor scan ID, if given
    fn precursor_id(&self) -> Option<&String>;
    /// The product scan ID, if given
    fn product_id(&self) -> Option<&String>;
    /// The activation process applied to the precursor ion
    fn activation(&self) -> &Activation;

    fn iter(&self) -> impl Iterator<Item = &SelectedIon>;

    /// Check if this precursor selection selected any ions
    fn has_ions(&self) -> bool {
        self.iter().count() > 0
    }

    fn last_ion(&self) -> Option<&SelectedIon> {
        self.iter().last()
    }

    fn iter_mut(&mut self) -> impl Iterator<Item = &mut SelectedIon>;
    fn ion_mut(&mut self) -> Option<&mut SelectedIon>;
    fn activation_mut(&mut self) -> &mut Activation;
    fn isolation_window_mut(&mut self) -> &mut IsolationWindow;
    fn add_ion(&mut self, ion: SelectedIon);
    fn last_ion_mut(&mut self) -> Option<&mut SelectedIon> {
        self.iter_mut().last()
    }
}

impl PrecursorSelection for Precursor {
    fn ion(&self) -> Option<&SelectedIon> {
        self.ions.first()
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

    fn ion_mut(&mut self) -> Option<&mut SelectedIon> {
        if self.ions.is_empty() {
            self.ions.push(SelectedIon::default())
        }
        self.ions.first_mut()
    }

    fn activation_mut(&mut self) -> &mut Activation {
        &mut self.activation
    }

    fn isolation_window_mut(&mut self) -> &mut IsolationWindow {
        &mut self.isolation_window
    }

    fn iter(&self) -> impl Iterator<Item = &SelectedIon> {
        self.ions.iter()
    }

    fn iter_mut(&mut self) -> impl Iterator<Item = &mut SelectedIon> {
        self.ions.iter_mut()
    }

    fn add_ion(&mut self, ion: SelectedIon) {
        self.ions.push(ion);
    }

    fn last_ion_mut(&mut self) -> Option<&mut SelectedIon> {
        if self.ions.is_empty() {
            self.ions.push(SelectedIon::default())
        }
        self.iter_mut().last()
    }
}


/**
Describes the polarity of a mass spectrum. A spectrum is either `Positive` (1+), `Negative` (-1)
or `Unknown` (0). The `Unknown` state is the default.
*/
#[repr(i8)]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ScanPolarity {
    #[default]
    /// The polarity of the spectrum is unknown
    Unknown = 0,
    /// The polarity is positive, where charge states represent the addition of positive charge
    Positive = 1,
    /// The polarity is negative, where the charge states reprsent the addition of electrons or other
    /// negatively charged adduction
    Negative = -1,
}

impl ScanPolarity {
    /// Get a signed integer representing the polarity of the spectrum. This
    /// assumes that unknown spectra are more likely to be positive so it
    /// is distinct from the raw integer value of the enum where [`Unknown`](Self::Unknown)
    /// is 0.
    pub fn sign(&self) -> i32 {
        match self {
            ScanPolarity::Unknown => 1,
            ScanPolarity::Positive => 1,
            ScanPolarity::Negative => -1,
        }
    }
}

impl Display for ScanPolarity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

/**
Describes the initial representation of the signal of a spectrum.

Though most formats explicitly have a method of either conveying a processing level
or an assumed level, the `Unknown` option is retained for partial initialization.
*/
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default, Hash, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SignalContinuity {
    #[default]
    Unknown = 0,
    /// The spectrum is centroided, indicating that its primary representation is that of a
    /// discrete peak list. There may be multiple peak lists and a profile spectrum may still
    /// be present on the same spectrum.
    Centroid = 3,
    /// The spectrum is profile, indicating that its primary representation is a continuous
    /// profile.
    Profile = 5,
}

impl SignalContinuity {

    pub const fn is_profile(&self) -> bool {
        matches!(self, Self::Profile)
    }

    pub const fn is_centroid(&self) -> bool {
        matches!(self, Self::Centroid)
    }
}

impl Display for SignalContinuity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}


/// An adapter type to make it possible to pass an `Option<Precursor>`, `Vec<Precursor>`,
/// or [`Precursor`] in an argument context.
#[derive(Debug)]
pub enum AsPrecursorCollection {
    Single(Option<Precursor>),
    Multiple(Vec<Precursor>)
}

impl From<AsPrecursorCollection> for Vec<Precursor> {
    fn from(value: AsPrecursorCollection) -> Self {
        match value {
            AsPrecursorCollection::Single(precursor) => precursor.map(|v| vec![v]).unwrap_or_default(),
            AsPrecursorCollection::Multiple(precursors) => precursors,
        }
    }
}

impl From<Precursor> for AsPrecursorCollection {
    fn from(value: Precursor) -> Self {
        Self::Single(Some(value))
    }
}

impl From<Option<Precursor>> for AsPrecursorCollection {
    fn from(value: Option<Precursor>) -> Self {
        Self::Single(value)
    }
}

impl From<Vec<Precursor>> for AsPrecursorCollection {
    fn from(value: Vec<Precursor>) -> Self {
        Self::Multiple(value)
    }
}


/**
The set of descriptive metadata that give context for how a mass spectrum was acquired
within a particular run. This forms the basis for a large portion of the [`SpectrumLike`]
trait.
*/
#[derive(Debug, Default, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SpectrumDescription {
    /// The spectrum's native identifier
    pub id: String,
    /// The ordinal sequence number for the spectrum
    pub index: usize,
    /// The degree of exponentiation of the spectrum, e.g MS1, MS2, MS3, etc
    pub ms_level: u8,

    /// The spectrum is in positive or negative mode
    pub polarity: ScanPolarity,

    /// The spectrum's main representation is as a peak list or a continuous
    /// profile
    pub signal_continuity: SignalContinuity,

    /// A set of controlled or uncontrolled descriptors of the spectrum not already
    /// covered by fields
    pub params: ParamList,

    /// A description of how the spectrum was acquired including time, scan windows, and more
    pub acquisition: Acquisition,

    /// The parent ion or ions and their isolation and activation description
    pub precursor: Vec<Precursor>,
}

impl SpectrumDescription {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        id: String,
        index: usize,
        ms_level: u8,
        polarity: ScanPolarity,
        signal_continuity: SignalContinuity,
        params: ParamList,
        acquisition: Acquisition,
        precursor: impl Into<AsPrecursorCollection>,
    ) -> Self {
        Self {
            id,
            index,
            ms_level,
            polarity,
            signal_continuity,
            params,
            acquisition,
            precursor: precursor.into().into(),
        }
    }

    crate::find_param_method!(title, &SCAN_TITLE, |p| p.as_str(), Option<Cow<'_, str>>);

    /// Find the type of spectrum described.
    ///
    /// A spectrum in `mzdata` is *usually* a mass spectrum of some sort, but that's not guaranteed
    /// to be the case. `mzdata` can handle non-MS spectra, but little of the signal processing
    /// machinery it provides currently supports those other kinds of data.
    pub fn spectrum_type(&self) -> Option<crate::meta::SpectrumType> {
        const SPECTRUM_TYPES: &[(crate::meta::SpectrumType, crate::params::ParamCow<'static>)] = crate::meta::SpectrumType::all_types();

        let conv_table: HashMap<CURIE, crate::meta::SpectrumType> = SPECTRUM_TYPES
            .iter()
            .map(|(t, v)| (v.curie().unwrap(), *t))
            .collect();

        for param in self.params().iter() {
            if let Some(c) = param.curie() {
                if let Some(t) = conv_table.get(&c) {
                    return Some(*t);
                }
            }
        }
        None
    }

    /// Set the kind of spectrum represented.
    pub fn set_spectrum_type(&mut self, spectrum_type: crate::meta::SpectrumType) {
        const SPECTRUM_TYPES: &[(crate::meta::SpectrumType, crate::params::ParamCow<'static>)] = crate::meta::SpectrumType::all_types();

        let to_insert: crate::params::ParamCow<'_> = spectrum_type.to_param();

        let conv_table: HashSet<CURIE> = SPECTRUM_TYPES
            .iter()
            .map(|(_, v)| v.curie().unwrap())
            .collect();

        for param in self.params_mut().iter_mut() {
            if let Some(c) = param.curie() {
                if conv_table.contains(&c) {
                    *param = to_insert.into();
                    return;
                }
            }
        }

        self.add_param(to_insert.into());
    }


}

impl_param_described!(Activation, SpectrumDescription);
impl_param_described_deferred!(SelectedIon, Acquisition, ScanEvent);

/// Types of chromatograms enumerated in the PSI-MS controlled vocabulary
#[derive(Debug, Clone, Copy, PartialEq, Default, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ChromatogramType {
    #[default]
    Unknown,
    TotalIonCurrentChromatogram,
    BasePeakChromatogram,
    SelectedIonCurrentChromatogram,
    SelectedIonMonitoringChromatogram,
    SelectedReactionMonitoringChromatogram,
    AbsorptionChromatogram,
    EmissionChromatogram,
    FlowRateChromatogram,
    PressureChromatogram,
    TemperatureChromatogram,
    ElectromagneticRadiationChromatogram,
}

impl ChromatogramType {
    pub fn from_accession(accession: AccessionIntCode) -> Option<Self> {
        let tp = match accession {
            1000235 => Self::TotalIonCurrentChromatogram,
            1000628 => Self::BasePeakChromatogram,
            1000627 => Self::SelectedIonCurrentChromatogram,
            1000472 => Self::SelectedIonMonitoringChromatogram,
            1000473 => Self::SelectedReactionMonitoringChromatogram,
            1000812 => Self::AbsorptionChromatogram,
            1000813 => Self::EmissionChromatogram,
            1003020 => Self::FlowRateChromatogram,
            1003019 => Self::PressureChromatogram,
            1000811 => Self::ElectromagneticRadiationChromatogram,
            1000626 => Self::Unknown,
            1002715 => Self::TemperatureChromatogram,
            _ => return None,
        };
        Some(tp)
    }

    pub fn is_electromagnetic_radiation(&self) -> bool {
        matches!(
            self,
            Self::AbsorptionChromatogram | Self::EmissionChromatogram | Self::ElectromagneticRadiationChromatogram
        )
    }

    pub fn is_aggregate(&self) -> bool {
        matches!(
            self,
            Self::TotalIonCurrentChromatogram
                | Self::BasePeakChromatogram
                | Self::PressureChromatogram
                | Self::FlowRateChromatogram
        )
    }

    pub fn is_ion_current(&self) -> bool {
        matches!(
            self,
            Self::SelectedIonCurrentChromatogram
                | Self::SelectedIonMonitoringChromatogram
                | Self::SelectedReactionMonitoringChromatogram
                | Self::TotalIonCurrentChromatogram
                | Self::BasePeakChromatogram
        )
    }

    pub fn to_curie(&self) -> CURIE {
        match self {
            Self::TotalIonCurrentChromatogram => CURIE::new(ControlledVocabulary::MS, 1000235),
            Self::BasePeakChromatogram => CURIE::new(ControlledVocabulary::MS, 1000628),
            Self::SelectedIonCurrentChromatogram => CURIE::new(ControlledVocabulary::MS, 1000627),
            Self::SelectedIonMonitoringChromatogram => {
                CURIE::new(ControlledVocabulary::MS, 1000472)
            }
            Self::SelectedReactionMonitoringChromatogram => {
                CURIE::new(ControlledVocabulary::MS, 1000473)
            }
            Self::AbsorptionChromatogram => CURIE::new(ControlledVocabulary::MS, 1000812),
            Self::ElectromagneticRadiationChromatogram => CURIE::new(ControlledVocabulary::MS, 1000811),
            Self::EmissionChromatogram => CURIE::new(ControlledVocabulary::MS, 1000813),
            Self::FlowRateChromatogram => CURIE::new(ControlledVocabulary::MS, 1003020),
            Self::PressureChromatogram => CURIE::new(ControlledVocabulary::MS, 1003019),
            Self::TemperatureChromatogram => CURIE::new(ControlledVocabulary::MS, 1002715),
            Self::Unknown => CURIE::new(ControlledVocabulary::MS, 1000626),
        }
    }
}

/// The set of descriptive metadata that give context for how a chromatogram was
/// recorded.
#[derive(Debug, Default, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ChromatogramDescription {
    pub id: String,
    pub index: usize,
    pub ms_level: Option<u8>,
    pub polarity: ScanPolarity,
    pub chromatogram_type: ChromatogramType,

    pub params: ParamList,
    pub precursor: Vec<Precursor>,
}

impl ChromatogramDescription {
    pub fn is_aggregate(&self) -> bool {
        self.chromatogram_type.is_aggregate()
    }

    pub fn is_electromagnetic_radiation(&self) -> bool {
        self.chromatogram_type.is_electromagnetic_radiation()
    }

    pub fn is_ion_current(&self) -> bool {
        self.chromatogram_type.is_ion_current()
    }
}

impl_param_described!(ChromatogramDescription);
