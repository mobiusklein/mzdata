//! Elements of controlled vocabularies used to describe mass spectra and their components.
//!
//! Directly maps to the usage of the PSI-MS controlled vocabulary in mzML
use std::borrow::Cow;
use std::collections::HashMap;
use std::fmt::Display;
use std::num;
use std::str::{self, FromStr};

use thiserror::Error;

#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum ValueType {
    String(Box<String>),
    Integer(i64),
    Float(f64),
    Other(Box<Vec<u8>>),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CURIE {
    controlled_vocabulary: ControlledVocabulary,
    accession: u32,
}

impl CURIE {
    pub const fn new(cv_id: ControlledVocabulary, accession: u32) -> Self {
        Self {
            controlled_vocabulary: cv_id,
            accession,
        }
    }
}

impl<T: ParamLike> PartialEq<T> for CURIE {
    fn eq(&self, other: &T) -> bool {
        if other.is_controlled() {
            false
        } else {
            if other.controlled_vocabulary().unwrap() != self.controlled_vocabulary {
                false
            } else if other.accession().unwrap() != self.accession {
                false
            } else {
                true
            }
        }
    }
}

#[derive(Debug, Error)]
pub enum CURIEParsingError {
    #[error("{0} is not a recognized controlled vocabulary")]
    UnknownControlledVocabulary(
        #[from]
        #[source]
        ControlledVocabularyResolutionError,
    ),
    #[error("Failed to parse accession number {0}")]
    AccessionParsingError(
        #[from]
        #[source]
        num::ParseIntError,
    ),
    #[error("Did not detect a namespace separator ':' token")]
    MissingNamespaceSeparator,
}

impl FromStr for CURIE {
    type Err = CURIEParsingError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut tokens = s.split(":");
        let cv = tokens.next().unwrap();
        let accession = tokens.next();
        if accession.is_none() {
            Err(CURIEParsingError::MissingNamespaceSeparator)
        } else {
            let cv: ControlledVocabulary = cv.parse::<ControlledVocabulary>()?;

            let accession = accession.unwrap().parse()?;
            Ok(CURIE::new(cv, accession))
        }
    }
}

pub fn curie_to_num(curie: &str) -> (Option<ControlledVocabulary>, Option<u32>) {
    let mut parts = curie.split(':');
    let prefix = match parts.next() {
        Some(v) => v.parse::<ControlledVocabulary>().unwrap().as_option(),
        None => None,
    };
    if let Some(k) = curie.split(':').nth(1) {
        match k.parse() {
            Ok(v) => (prefix, Some(v)),
            Err(_) => (prefix, None),
        }
    } else {
        (prefix, None)
    }
}

/// Describe a controlled vocabulary parameter or a user-defined parameter
pub trait ParamLike {
    fn name(&self) -> &str;
    fn value(&self) -> &str;
    fn accession(&self) -> Option<u32>;
    fn controlled_vocabulary(&self) -> Option<ControlledVocabulary>;
    fn unit(&self) -> Unit;
    fn is_ms(&self) -> bool {
        if let Some(cv) = self.controlled_vocabulary() {
            cv == ControlledVocabulary::MS
        } else {
            false
        }
    }

    fn parse<T: str::FromStr>(&self) -> Result<T, T::Err> {
        self.value().parse::<T>()
    }

    fn is_controlled(&self) -> bool {
        self.accession().is_some()
    }

    fn curie(&self) -> Option<String> {
        if !self.is_controlled() {
            None
        } else {
            let cv = &self.controlled_vocabulary().unwrap();
            let acc = self.accession().unwrap();
            let accession_str = format!("{}:{:07}", cv.prefix(), acc);
            Some(accession_str)
        }
    }
}

/// A statically allocate-able or non-owned data version of [`Param`]
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ParamCow<'a> {
    pub name: Cow<'a, str>,
    pub value: Cow<'a, str>,
    pub accession: Option<u32>,
    pub controlled_vocabulary: Option<ControlledVocabulary>,
    pub unit: Unit,
}

impl ParamCow<'static> {
    pub const fn const_new(
        name: &'static str,
        value: &'static str,
        accession: Option<u32>,
        controlled_vocabulary: Option<ControlledVocabulary>,
        unit: Unit,
    ) -> Self {
        Self {
            name: Cow::Borrowed(name),
            value: Cow::Borrowed(value),
            accession,
            controlled_vocabulary,
            unit,
        }
    }
}

impl<'a> ParamCow<'a> {
    pub fn new(
        name: Cow<'a, str>,
        value: Cow<'a, str>,
        accession: Option<u32>,
        controlled_vocabulary: Option<ControlledVocabulary>,
        unit: Unit,
    ) -> Self {
        Self {
            name,
            value,
            accession,
            controlled_vocabulary,
            unit,
        }
    }

    pub fn parse<T: str::FromStr>(&self) -> Result<T, T::Err> {
        self.value.parse::<T>()
    }

    pub fn is_controlled(&self) -> bool {
        self.accession.is_some()
    }
}

impl<'a> ParamLike for ParamCow<'a> {
    fn name(&self) -> &str {
        &self.name
    }

    fn value(&self) -> &str {
        &self.value
    }

    fn accession(&self) -> Option<u32> {
        self.accession
    }

    fn controlled_vocabulary(&self) -> Option<ControlledVocabulary> {
        self.controlled_vocabulary
    }

    fn unit(&self) -> Unit {
        self.unit
    }
}

impl<'a> From<ParamCow<'a>> for Param {
    fn from(value: ParamCow<'a>) -> Self {
        Param {
            name: value.name.into_owned(),
            value: value.value.into_owned(),
            accession: value.accession,
            controlled_vocabulary: value.controlled_vocabulary,
            unit: value.unit,
        }
    }
}

/// A controlled vocabulary or user parameter
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Param {
    pub name: String,
    pub value: String,
    pub accession: Option<u32>,
    pub controlled_vocabulary: Option<ControlledVocabulary>,
    pub unit: Unit,
}

impl Display for Param {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut body = if self.is_controlled() {
            format!(
                "{}:{}|{}={}",
                String::from_utf8_lossy(self.controlled_vocabulary.unwrap().as_bytes()),
                self.accession.unwrap(),
                self.name,
                self.value
            )
        } else {
            format!("{}={}", self.name, self.value)
        };
        if self.unit != Unit::Unknown {
            body.extend(format!(" {}", self.unit).chars());
        };
        f.write_str(body.as_str())
    }
}

impl Param {
    pub fn new() -> Param {
        Param {
            ..Default::default()
        }
    }

    pub fn new_key_value(name: String, value: String) -> Param {
        let mut inst = Self::new();
        inst.name = name;
        inst.value = value;
        inst
    }

    pub fn parse<T: str::FromStr>(&self) -> Result<T, T::Err> {
        self.value.parse::<T>()
    }

    pub fn is_controlled(&self) -> bool {
        self.accession.is_some()
    }

    pub fn curie(&self) -> Option<String> {
        if !self.is_controlled() {
            None
        } else {
            let cv = &self.controlled_vocabulary.unwrap();
            let acc = self.accession.unwrap();
            let accession_str = format!("{}:{:07}", cv.prefix(), acc);
            Some(accession_str)
        }
    }

    pub fn with_unit<S: AsRef<str>, A: AsRef<str>>(mut self, accession: S, name: A) -> Param {
        self.unit = Unit::from_accession(accession.as_ref());
        if matches!(self.unit, Unit::Unknown) {
            self.unit = Unit::from_name(name.as_ref());
        }
        self
    }

    pub fn with_unit_t(mut self, unit: &Unit) -> Param {
        self.unit = *unit;
        self
    }
}

impl ParamLike for Param {
    fn name(&self) -> &str {
        &self.name
    }

    fn value(&self) -> &str {
        &self.value
    }

    fn accession(&self) -> Option<u32> {
        self.accession
    }

    fn controlled_vocabulary(&self) -> Option<ControlledVocabulary> {
        self.controlled_vocabulary
    }

    fn unit(&self) -> Unit {
        self.unit
    }
}

/// Controlled vocabularies used in mass spectrometry data files
#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
pub enum ControlledVocabulary {
    MS,
    UO,
    Unknown,
}

const MS_CV: &str = "MS";
const UO_CV: &str = "UO";
const MS_CV_BYTES: &[u8] = MS_CV.as_bytes();
const UO_CV_BYTES: &[u8] = UO_CV.as_bytes();

impl ControlledVocabulary {
    pub fn prefix(&self) -> Cow<'static, str> {
        match &self {
            Self::MS => Cow::Borrowed(MS_CV),
            Self::UO => Cow::Borrowed(UO_CV),
            Self::Unknown => panic!("Cannot encode unknown CV"),
        }
    }

    pub fn as_bytes(&self) -> &'static [u8] {
        match &self {
            Self::MS => MS_CV_BYTES,
            Self::UO => UO_CV_BYTES,
            Self::Unknown => panic!("Cannot encode unknown CV"),
        }
    }

    pub fn as_option(&self) -> Option<Self> {
        match self {
            Self::Unknown => None,
            _ => Some(*self),
        }
    }

    pub fn param<A: AsRef<str>, S: Into<String>>(&self, accession: A, name: S) -> Param {
        let mut param = Param::new();
        param.controlled_vocabulary = Some(*self);
        param.name = name.into();
        if let Some(nb) = accession.as_ref().split(':').nth(1) {
            param.accession = Some(nb.parse().unwrap_or_else(|_| {
                panic!(
                    "Expected accession to be numeric, got {}",
                    accession.as_ref()
                )
            }));
        }
        param
    }

    pub const fn const_param(
        &self,
        name: &'static str,
        value: &'static str,
        accession: u32,
        unit: Unit,
    ) -> ParamCow<'static> {
        ParamCow {
            name: Cow::Borrowed(name),
            value: Cow::Borrowed(value),
            accession: Some(accession),
            controlled_vocabulary: Some(*self),
            unit,
        }
    }

    pub const fn const_param_ident(&self, name: &'static str, accession: u32) -> ParamCow<'static> {
        self.const_param(name, "", accession, Unit::Unknown)
    }

    pub const fn const_param_ident_unit(
        &self,
        name: &'static str,
        accession: u32,
        unit: Unit,
    ) -> ParamCow<'static> {
        self.const_param(name, "", accession, unit)
    }

    pub fn param_val<S: Into<String>, A: AsRef<str>, V: ToString>(
        &self,
        accession: A,
        name: S,
        value: V,
    ) -> Param {
        let mut param = self.param(accession, name);
        param.value = value.to_string();
        param
    }
}

#[doc(hidden)]
#[derive(Debug, Clone, Error)]
pub enum ControlledVocabularyResolutionError {
    #[error("Unrecognized controlled vocabulary {0}")]
    UnknownControlledVocabulary(String),
}

impl FromStr for ControlledVocabulary {
    type Err = ControlledVocabularyResolutionError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "MS" | "PSI-MS" => Ok(Self::MS),
            "UO" => Ok(Self::UO),
            _ => Ok(Self::Unknown),
        }
    }
}

pub type ParamList = Vec<Param>;

pub type ParamMap = HashMap<String, Param>;

pub trait ParamDescribed {
    fn params(&self) -> &[Param];
    fn params_mut(&mut self) -> &mut ParamList;

    fn add_param(&mut self, param: Param) {
        self.params_mut().push(param);
    }

    fn remove_param(&mut self, index: usize) -> Param {
        self.params_mut().remove(index)
    }

    fn get_param_by_name(&self, name: &str) -> Option<&Param> {
        self.params().iter().find(|&param| param.name == name)
    }

    fn get_param_by_curie(&self, curie: &CURIE) -> Option<&Param> {
        self.params().iter().find(|&param| curie == param)
    }

    fn get_param_by_accession(&self, accession: &str) -> Option<&Param> {
        let (cv, acc_num) = curie_to_num(accession);
        return self
            .params()
            .iter()
            .find(|&param| param.accession == acc_num && param.controlled_vocabulary == cv);
    }
}

#[macro_export]
macro_rules! impl_param_described {
    ($($t:ty), +) => {$(

        impl $crate::params::ParamDescribed for $t {
            fn params(&self) -> &[$crate::params::Param] {
                return &self.params
            }

            fn params_mut(&mut self) -> &mut $crate::params::ParamList {
                return &mut self.params
            }
        }
    )+};
}

#[doc(hidden)]
pub const _EMPTY_PARAM: &[Param] = &[];

#[macro_export]
macro_rules! impl_param_described_deferred {
    ($($t:ty), +) => {$(
        impl $crate::params::ParamDescribed for $t {
            fn params(&self) -> &[$crate::params::Param] {
                match &self.params {
                    Some(val) => &val,
                    None => {
                        $crate::params::_EMPTY_PARAM
                    }
                }
            }

            fn params_mut(&mut self) -> &mut $crate::params::ParamList {
                let val = &mut self.params;
                if val.is_some() {
                    return val.as_deref_mut().unwrap()
                } else {
                    *val = Some(Box::default());
                    return val.as_deref_mut().unwrap()
                }
            }
        }
    )+};
}

/// Units that a term's value might have
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Unit {
    // Mass
    MZ,
    Mass,
    PartsPerMillion,
    Nanometer,

    // Time
    Minute,
    Second,
    Millisecond,
    VoltSecondPerSquareCentimeter,

    // Intensity
    DetectorCounts,
    PercentBasePeak,
    PercentBasePeakTimes100,
    AbsorbanceUnit,
    CountsPerSecond,

    // Collision Energy
    Electronvolt,
    PercentElectronVolt,
    Volt,

    Unknown,
}

impl Unit {
    pub const fn for_param(&self) -> (&'static str, &'static str) {
        match self {
            Self::Millisecond => ("UO:0000028", "millisecond"),
            Self::Second => ("UO:0000010", "second"),
            Self::Minute => ("UO:0000031", "minute"),

            Self::MZ => ("MS:1000040", "m/z"),
            Self::Mass => ("UO:000221", "dalton"),

            Self::DetectorCounts => ("MS:1000131", "number of detector counts"),
            Self::PercentBasePeak => ("MS:1000132", "percent of base peak"),
            Self::PercentBasePeakTimes100 => ("MS:1000905", "percent of base peak times 100"),
            Self::AbsorbanceUnit => ("UO:0000269", "absorbance unit"),
            Self::CountsPerSecond => ("MS:1000814", "counts per second"),

            Self::Electronvolt => ("UO:0000266", "electronvolt"),
            Self::PercentElectronVolt => ("UO:0000187", "percent"),

            _ => ("", ""),
        }
    }

    pub fn from_name(name: &str) -> Unit {
        match name {
            "millisecond" => Self::Millisecond,
            "second" => Self::Second,
            "minute" => Self::Minute,

            "m/z" => Self::MZ,
            "dalton" => Self::Mass,

            "number of detector counts" => Self::DetectorCounts,
            "percent of base peak" => Self::PercentBasePeak,
            "percent of base peak times 100" => Self::PercentBasePeakTimes100,
            "absorbance unit" => Self::AbsorbanceUnit,
            "counts per second" => Self::CountsPerSecond,

            "electronvolt" => Self::Electronvolt,
            "percent" => Self::PercentElectronVolt,
            _ => Unit::Unknown,
        }
    }

    pub fn from_accession(acc: &str) -> Unit {
        match acc {
            "UO:0000028" => Self::Millisecond,
            "UO:0000010" => Self::Second,
            "UO:0000031" => Self::Minute,

            "MS:1000040" => Self::MZ,
            "UO:000221" => Self::Mass,

            "MS:1000131" => Self::DetectorCounts,
            "MS:1000132" => Self::PercentBasePeak,
            "MS:1000905" => Self::PercentBasePeakTimes100,
            "UO:0000269" => Self::AbsorbanceUnit,
            "MS:1000814" => Self::CountsPerSecond,

            "UO:0000266" => Self::Electronvolt,
            "UO:0000187" => Self::PercentElectronVolt,
            _ => Unit::Unknown,
        }
    }

    pub const fn from_param(param: &Param) -> Unit {
        param.unit
    }
}

impl Default for Unit {
    fn default() -> Self {
        Self::Unknown
    }
}

impl Display for Unit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(format!("{:?}", self).as_str())
    }
}
