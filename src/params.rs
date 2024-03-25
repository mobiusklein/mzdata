//! Elements of controlled vocabularies used to describe mass spectra and their components.
//!
//! Directly maps to the usage of the PSI-MS controlled vocabulary in mzML
use std::borrow::Cow;
use std::collections::HashMap;
use std::fmt::Display;
use std::{io, num};
use std::str::{self, FromStr};

use thiserror::Error;

#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, PartialOrd, Default)]
pub enum Value {
    String(String),
    Float(f64),
    Int(i64),
    Buffer(Box<[u8]>),
    #[default]
    Empty
}

impl Value {}

#[doc(hidden)]
#[derive(Debug, Clone, Error, PartialEq)]
pub enum ParamValueParseError {
    #[error("Failed to extract a float from {0:?}")]
    FailedToExtractFloat(Option<String>),
    #[error("Failed to extract a int from {0:?}")]
    FailedToExtractInt(Option<String>),
    #[error("Failed to extract a string")]
    FailedToExtractString,
    #[error("Failed to extract a buffer")]
    FailedToExtractBuffer
}

impl FromStr for Value {
    type Err = ParamValueParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Ok(Self::Empty)
        }
        if let Ok(value) = s.parse::<i64>() {
            Ok(Self::Int(value))
        } else if let Ok(value) = s.parse::<f64>() {
            Ok(Self::Float(value))
        } else {
            Ok(Self::String(s.to_string()))
        }
    }
}

impl ToString for Value {
    fn to_string(&self) -> String {
        match self {
            Value::String(v) => v.to_string(),
            Value::Float(v) => v.to_string(),
            Value::Int(v) => v.to_string(),
            Value::Buffer(v) => String::from_utf8_lossy(v).to_string(),
            Value::Empty => "".to_string(),
        }
    }
}

impl From<ParamValueParseError> for io::Error {
    fn from(value: ParamValueParseError) -> Self {
        Self::new(io::ErrorKind::InvalidData, value)
    }
}

impl Value {

    pub fn wrap(s: &str) -> Self {
        if s.is_empty() {
            Self::Empty
        }
        else if let Ok(value) = s.parse::<i64>() {
            Self::Int(value)
        } else if let Ok(value) = s.parse::<f64>() {
            Self::Float(value)
        } else {
            Self::String(s.to_string())
        }
    }

    pub fn is_empty(&self) -> bool {
        matches!(self, Self::Empty)
    }

    pub fn is_int(&self) -> bool {
        matches!(self, Self::Int(_))
    }

    pub fn is_float(&self) -> bool {
        matches!(self, Self::Float(_))
    }

    pub fn is_buffer(&self) -> bool {
        matches!(self, Self::Buffer(_))
    }

    pub fn is_string(&self) -> bool {
        matches!(self, Self::String(_))
    }

    pub fn to_f64(&self) -> Result<f64, ParamValueParseError> {
        if let Self::Float(val) = self {
            Ok(*val)
        } else {
            Err(ParamValueParseError::FailedToExtractFloat(Some(self.to_string())))
        }
    }

    pub fn to_i64(&self) -> Result<i64, ParamValueParseError> {
        if let Self::Int(val) = self {
            Ok(*val)
        } else {
            Err(ParamValueParseError::FailedToExtractInt(Some(self.to_string())))
        }
    }

    pub fn to_str(&self) -> Result<Cow<'_, str>, ParamValueParseError> {
        if let Self::String(val) = self {
            Ok(Cow::Borrowed(val))
        } else {
            Ok(Cow::Owned(self.to_string()))
        }
    }

    pub fn to_buffer(&self) -> Result<Cow<'_, [u8]>, ParamValueParseError> {
        if let  Self::Buffer(val) = self {
            Ok(Cow::Borrowed(val))
        } else if let Self::String(val) = self {
            Ok(Cow::Borrowed(val.as_bytes()))
        } else {
            Err(ParamValueParseError::FailedToExtractBuffer)
        }
    }
}


/// A CURIE is a namespace + accession identifier
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

    pub fn as_param(&self) -> Param {
        let mut param = Param::new();
        param.controlled_vocabulary = Some(self.controlled_vocabulary);
        param.accession = Some(self.accession);
        param
    }
}

impl Display for CURIE {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}",
            self.controlled_vocabulary.prefix(),
            self.accession
        )
    }
}

impl<T: ParamLike> PartialEq<T> for CURIE {
    fn eq(&self, other: &T) -> bool {
        if !other.is_controlled()
            || other.controlled_vocabulary().unwrap() != self.controlled_vocabulary
        {
            return false;
        } else {
            other.accession().unwrap() == self.accession
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
        let mut tokens = s.split(':');
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

    fn curie(&self) -> Option<CURIE> {
        if !self.is_controlled() {
            None
        } else {
            let cv = self.controlled_vocabulary().unwrap();
            let acc = self.accession().unwrap();
            // let accession_str = format!("{}:{:07}", cv.prefix(), acc);
            Some(CURIE::new(cv, acc))
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

impl<'a> PartialEq<CURIE> for ParamCow<'a> {
    fn eq(&self, other: &CURIE) -> bool {
        other.eq(self)
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

    pub fn new_key_value<K: Into<String>, V: Into<String>>(name: K, value: V) -> Param {
        let mut inst = Self::new();
        inst.name = name.into();
        inst.value = value.into();
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

impl PartialEq<CURIE> for Param {
    fn eq(&self, other: &CURIE) -> bool {
        other.eq(self)
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


/// Anything that can be converted into an accession code portion of a [`CURIE`]
#[derive(Debug, Clone)]
pub enum AccessionLike<'a> {
    Text(Cow<'a, str>),
    Number(u32)
}

impl<'a> From<u32> for AccessionLike<'a> {
    fn from(value: u32) -> Self {
        Self::Number(value)
    }
}

impl<'a> From<&'a str> for AccessionLike<'a> {
    fn from(value: &'a str) -> Self {
        Self::Text(Cow::Borrowed(value))
    }
}

impl<'a> From<String> for AccessionLike<'a> {
    fn from(value: String) -> Self {
        Self::Text(Cow::Owned(value))
    }
}

impl<'a> ControlledVocabulary {
    pub const fn prefix(&self) -> Cow<'static, str> {
        match &self {
            Self::MS => Cow::Borrowed(MS_CV),
            Self::UO => Cow::Borrowed(UO_CV),
            Self::Unknown => panic!("Cannot encode unknown CV"),
        }
    }

    pub const fn as_bytes(&self) -> &'static [u8] {
        match &self {
            Self::MS => MS_CV_BYTES,
            Self::UO => UO_CV_BYTES,
            Self::Unknown => panic!("Cannot encode unknown CV"),
        }
    }

    pub const fn as_option(&self) -> Option<Self> {
        match self {
            Self::Unknown => None,
            _ => Some(*self),
        }
    }

    pub fn param<A: Into<AccessionLike<'a>>, S: Into<String>>(&self, accession: A, name: S) -> Param {
        let mut param = Param::new();
        param.controlled_vocabulary = Some(*self);
        param.name = name.into();

        let accession: AccessionLike = accession.into();

        match accession {
            AccessionLike::Text(s) => {
                if let Some(nb) = s.split(":").last() {
                    param.accession = Some(nb.parse().unwrap_or_else(|_| {
                        panic!(
                            "Expected accession to be numeric, got {}",
                            s
                        )
                    }))
                }
            },
            AccessionLike::Number(n) => {
                param.accession = Some(n)
            },
        }
        param
    }

    pub const fn curie(&self, accession: u32) -> CURIE {
        CURIE::new(*self, accession)
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

    pub fn param_val<S: Into<String>, A: Into<AccessionLike<'a>>, V: ToString>(
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

impl ParamDescribed for ParamList {
    fn params(&self) -> &[Param] {
        self
    }

    fn params_mut(&mut self) -> &mut ParamList {
        self
    }
}

/// Implement the [`ParamDescribed`] trait for type `$t`, referencing a `params` member
/// of type `Vec<`[`Param`]`>`.
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

/// Implement the [`ParamDescribed`] trait for type `$t`, referencing a `params` member
/// that is an `Option<Vec<`[`Param`]`>>` that will lazily be initialized automatically
/// when it is accessed mutably.
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
    Unknown,

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

    pub const fn from_curie(acc: &CURIE) -> Unit {
        match acc {
            CURIE { controlled_vocabulary: ControlledVocabulary::UO, accession: 28 } => Self::Millisecond,
            CURIE { controlled_vocabulary: ControlledVocabulary::UO, accession: 10 } => Self::Second,
            CURIE { controlled_vocabulary: ControlledVocabulary::UO, accession: 31 } => Self::Minute,

            CURIE { controlled_vocabulary: ControlledVocabulary::MS, accession: 1000040 } => Self::MZ,
            CURIE { controlled_vocabulary: ControlledVocabulary::UO, accession: 221 } => Self::Mass,

            CURIE { controlled_vocabulary: ControlledVocabulary::MS, accession: 1000131 } => Self::DetectorCounts,
            CURIE { controlled_vocabulary: ControlledVocabulary::MS, accession: 1000132 } => Self::PercentBasePeak,
            CURIE { controlled_vocabulary: ControlledVocabulary::UO, accession: 269 } => Self::AbsorbanceUnit,
            CURIE { controlled_vocabulary: ControlledVocabulary::MS, accession: 1000814 } => Self::CountsPerSecond,

            CURIE { controlled_vocabulary: ControlledVocabulary::UO, accession: 266 } => Self::Electronvolt,
            CURIE { controlled_vocabulary: ControlledVocabulary::UO, accession: 187 } => Self::PercentElectronVolt,

            _ => Unit::Unknown
        }
    }

    pub const fn to_curie(&self) -> Option<CURIE> {
        match self {
            Self::Millisecond => Some(CURIE { controlled_vocabulary: ControlledVocabulary::UO, accession: 28 }),
            Self::Second => Some(CURIE { controlled_vocabulary: ControlledVocabulary::UO, accession: 10 }),
            Self::Minute => Some(CURIE { controlled_vocabulary: ControlledVocabulary::UO, accession: 31 }),

            Self::MZ => Some(CURIE { controlled_vocabulary: ControlledVocabulary::MS, accession: 1000040 }),
            Self::Mass => Some(CURIE { controlled_vocabulary: ControlledVocabulary::UO, accession: 221 }),

            Self::DetectorCounts => Some(CURIE { controlled_vocabulary: ControlledVocabulary::MS, accession: 1000131 }),
            Self::PercentBasePeak => Some(CURIE { controlled_vocabulary: ControlledVocabulary::MS, accession: 1000132 }),
            Self::AbsorbanceUnit => Some(CURIE { controlled_vocabulary: ControlledVocabulary::UO, accession: 269 }),
            Self::CountsPerSecond => Some(CURIE { controlled_vocabulary: ControlledVocabulary::MS, accession: 1000814 }),

            Self::Electronvolt => Some(CURIE { controlled_vocabulary: ControlledVocabulary::UO, accession: 266 }),
            Self::PercentElectronVolt => Some(CURIE { controlled_vocabulary: ControlledVocabulary::UO, accession: 187 }),

            _ => None
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
