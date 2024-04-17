//! Elements of controlled vocabularies used to describe mass spectra and their components.
//!
//! Directly maps to the usage of the PSI-MS controlled vocabulary in mzML
use std::borrow::Cow;
use std::convert::TryFrom;
use std::fmt::Display;
use std::str::{self, FromStr};
use std::{io, num};

use thiserror::Error;

/// An owned parameter value that may be a string, a number, or empty. It is intended to
/// be paired with the [`ParamValue`] trait.
///
/// The borrowed equivalent of this type is [`ValueRef`].
#[derive(Debug, Clone, PartialEq, PartialOrd, Default)]
pub enum Value {
    /// A text value of arbitrary length
    String(String),
    /// A floating point number
    Float(f64),
    /// A integral number
    Int(i64),
    /// Arbitrary binary data
    Buffer(Box<[u8]>),
    /// No value specified
    #[default]
    Empty,
}

impl Eq for Value {}

impl From<String> for Value {
    fn from(value: String) -> Self {
        Value::new(value)
    }
}

impl From<&str> for Value {
    fn from(value: &str) -> Self {
        Value::wrap(value)
    }
}

impl From<Cow<'_, str>> for Value {
    fn from(value: Cow<'_, str>) -> Self {
        Value::wrap(&value)
    }
}

/// Access a parameter's value, with specific coercion rules
/// and eager type conversion.
pub trait ParamValue {
    /// Check if the value is empty
    fn is_empty(&self) -> bool;

    /// Check if the value is an integer
    fn is_i64(&self) -> bool;

    /// Check if the value is a floating point
    /// number explicitly. An integral number might
    /// still be usable as a floating point number
    fn is_f64(&self) -> bool;

    /// Check if the value is an arbitrary buffer
    fn is_buffer(&self) -> bool;

    /// Check if the value is stored as an explicit string.
    /// All variants can be coerced to a string.
    fn is_str(&self) -> bool;

    /// Check if the value is of either numeric type.
    fn is_numeric(&self) -> bool {
        self.is_i64() | self.is_f64()
    }

    /// Get the value as an `f64`, if possible
    fn to_f64(&self) -> Result<f64, ParamValueParseError>;

    /// Get the value as an `f32`, if possible
    fn to_f32(&self) -> Result<f32, ParamValueParseError> {
        let v = self.to_f64()?;
        Ok(v as f32)
    }

    /// Get the value as an `i64`, if possible
    fn to_i64(&self) -> Result<i64, ParamValueParseError>;

    /// Get the value as an `i32`, if possible
    fn to_i32(&self) -> Result<i32, ParamValueParseError> {
        let v = self.to_i64()?;
        Ok(v as i32)
    }

    /// Get the value as an `u64`, if possible
    fn to_u64(&self) -> Result<u64, ParamValueParseError> {
        let v = self.to_i64()?;
        Ok(v as u64)
    }

    /// Get the value as a string
    fn to_str(&self) -> Cow<'_, str>;

    /// Get the value as a string, possibly borrowed
    fn as_str(&self) -> Cow<'_, str> {
        self.to_str()
    }

    /// Get the value as a byte buffer, if possible.
    ///
    /// The intent here is distinct from [`ParamValue::as_bytes`]. The byte buffer
    /// represents the byte representation of the native value, while
    /// [`ParamValue::as_bytes`] is a byte string of the string representation.
    fn to_buffer(&self) -> Result<Cow<'_, [u8]>, ParamValueParseError>;

    /// Convert the value's string representation to `T` if possible
    fn parse<T: FromStr>(&self) -> Result<T, T::Err>;

    /// Convert the value to a byte string, the bytes
    /// of the string representation.
    fn as_bytes(&self) -> Cow<'_, [u8]>;

    /// Get a reference to the stored value
    fn as_ref(&self) -> ValueRef<'_>;
}

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
    FailedToExtractBuffer,
}

impl FromStr for Value {
    type Err = ParamValueParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Ok(Self::Empty);
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

impl Display for Value {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Value::String(v) => f.write_str(&v),
            Value::Float(v) => v.fmt(f),
            Value::Int(v) => v.fmt(f),
            Value::Buffer(v) => f.write_str(&String::from_utf8_lossy(v)),
            Value::Empty => f.write_str(""),
        }
    }
}

impl From<ParamValueParseError> for io::Error {
    fn from(value: ParamValueParseError) -> Self {
        Self::new(io::ErrorKind::InvalidData, value)
    }
}

impl Value {
    /// Convert a string value into a precise value type by trying
    /// successive types to parse, defaulting to storing the string
    /// as-is.
    ///
    /// This takes ownership of the string. To coerce from a borrowed
    /// string see [`Value::wrap`].
    pub fn new(s: String) -> Self {
        if s.is_empty() {
            Self::Empty
        } else if let Ok(value) = s.parse::<i64>() {
            Self::Int(value)
        } else if let Ok(value) = s.parse::<f64>() {
            Self::Float(value)
        } else {
            Self::String(s)
        }
    }

    /// Convert a borrowed string value into a precise value type by trying
    /// successive types to parse, defaulting to storing the string
    /// as-is.
    ///
    /// This only makes a copy of the string if it cannot be parsed into a
    /// numeric type.
    pub fn wrap(s: &str) -> Self {
        if s.is_empty() {
            Self::Empty
        } else if let Ok(value) = s.parse::<i64>() {
            Self::Int(value)
        } else if let Ok(value) = s.parse::<f64>() {
            Self::Float(value)
        } else {
            Self::String(s.to_string())
        }
    }

    fn is_empty(&self) -> bool {
        matches!(self, Self::Empty)
    }

    fn is_i64(&self) -> bool {
        matches!(self, Self::Int(_))
    }

    fn is_f64(&self) -> bool {
        matches!(self, Self::Float(_))
    }

    fn is_buffer(&self) -> bool {
        matches!(self, Self::Buffer(_))
    }

    fn is_str(&self) -> bool {
        matches!(self, Self::String(_))
    }

    /// Store the value as a floating point number
    pub fn coerce_f64(&mut self) -> Result<(), ParamValueParseError> {
        let value = self.to_f64()?;
        *self = Self::Float(value);
        Ok(())
    }

    /// Store the value as an integer
    pub fn coerce_i64(&mut self) -> Result<(), ParamValueParseError> {
        let value = self.to_i64()?;
        *self = Self::Int(value);
        Ok(())
    }

    /// Store the value as a string
    pub fn coerce_str(&mut self) -> Result<(), ParamValueParseError> {
        let value = self.to_string();
        *self = Self::String(value);
        Ok(())
    }

    /// Discard the value, leaving this value [`Value::Empty`]
    pub fn coerce_empty(&mut self) {
        *self = Self::Empty;
    }

    /// Store the value as a byte buffer
    pub fn coerce_buffer(&mut self) -> Result<(), ParamValueParseError> {
        let buffer = self.to_buffer()?;
        *self = Self::Buffer(buffer.into());
        Ok(())
    }

    fn parse<T: FromStr>(&self) -> Result<T, T::Err> {
        match self {
            Value::String(s) => s.parse(),
            Value::Float(v) => v.to_string().parse(),
            Value::Int(i) => i.to_string().parse(),
            Value::Buffer(b) => String::from_utf8_lossy(b).parse(),
            Value::Empty => "".parse(),
        }
    }

    fn to_f64(&self) -> Result<f64, ParamValueParseError> {
        if let Self::Float(val) = self {
            return Ok(*val);
        } else if let Self::Int(val) = self {
            return Ok(*val as f64);
        } else if let Self::String(val) = self {
            if let Ok(v) = val.parse() {
                return Ok(v);
            }
        }
        return Err(ParamValueParseError::FailedToExtractFloat(Some(
            self.to_string(),
        )));
    }

    fn to_i64(&self) -> Result<i64, ParamValueParseError> {
        if let Self::Int(val) = self {
            return Ok(*val);
        } else if let Self::Float(val) = self {
            return Ok(*val as i64);
        } else if let Self::String(val) = self {
            if let Ok(v) = val.parse() {
                return Ok(v);
            }
        }
        return Err(ParamValueParseError::FailedToExtractInt(Some(
            self.to_string(),
        )));
    }

    fn to_str(&self) -> Cow<'_, str> {
        if let Self::String(val) = self {
            Cow::Borrowed(val)
        } else {
            Cow::Owned(self.to_string())
        }
    }

    fn to_buffer(&self) -> Result<Cow<'_, [u8]>, ParamValueParseError> {
        if let Self::Buffer(val) = self {
            Ok(Cow::Borrowed(val))
        } else if let Self::String(val) = self {
            Ok(Cow::Borrowed(val.as_bytes()))
        } else {
            Err(ParamValueParseError::FailedToExtractBuffer)
        }
    }

    fn as_ref(&self) -> ValueRef<'_> {
        self.into()
    }
}

impl ParamValue for Value {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn is_i64(&self) -> bool {
        self.is_i64()
    }

    fn is_f64(&self) -> bool {
        self.is_f64()
    }

    fn is_buffer(&self) -> bool {
        self.is_buffer()
    }

    fn is_str(&self) -> bool {
        self.is_str()
    }

    fn to_f64(&self) -> Result<f64, ParamValueParseError> {
        self.to_f64()
    }

    fn to_i64(&self) -> Result<i64, ParamValueParseError> {
        self.to_i64()
    }

    fn to_str(&self) -> Cow<'_, str> {
        self.to_str()
    }

    fn to_buffer(&self) -> Result<Cow<'_, [u8]>, ParamValueParseError> {
        self.to_buffer()
    }

    fn parse<T: FromStr>(&self) -> Result<T, T::Err> {
        self.parse()
    }

    fn as_bytes(&self) -> Cow<'_, [u8]> {
        match self {
            Self::String(v) => Cow::Borrowed(v.as_bytes()),
            Self::Buffer(v) => Cow::Borrowed(v.as_ref()),
            Self::Float(v) => Cow::Owned(v.to_string().into_bytes()),
            Self::Int(v) => Cow::Owned(v.to_string().into_bytes()),
            Self::Empty => Cow::Borrowed(b""),
        }
    }

    fn as_ref(&self) -> ValueRef<'_> {
        self.into()
    }
}

impl PartialEq<String> for Value {
    fn eq(&self, other: &String) -> bool {
        self.as_str() == other.as_str()
    }
}

impl PartialEq<str> for Value {
    fn eq(&self, other: &str) -> bool {
        self.as_str() == other
    }
}

impl PartialEq<&str> for Value {
    fn eq(&self, other: &&str) -> bool {
        self.as_str() == *other
    }
}

impl PartialEq<i64> for Value {
    fn eq(&self, other: &i64) -> bool {
        if let Self::Int(val) = self {
            val == other
        } else {
            false
        }
    }
}

impl PartialEq<f64> for Value {
    fn eq(&self, other: &f64) -> bool {
        if let Self::Float(val) = self {
            val == other
        } else {
            false
        }
    }
}

macro_rules! param_value_int {
    ($val:ty) => {
        impl From<$val> for Value {
            fn from(value: $val) -> Self {
                Self::Int(value as i64)
            }
        }

        impl From<&$val> for Value {
            fn from(value: &$val) -> Self {
                Self::Int(*value as i64)
            }
        }

        impl From<Option<$val>> for Value {
            fn from(value: Option<$val>) -> Self {
                if let Some(v) = value {
                    Self::Int(v as i64)
                } else {
                    Self::Empty
                }
            }
        }
    };
}

macro_rules! param_value_float {
    ($val:ty) => {
        impl From<$val> for Value {
            fn from(value: $val) -> Self {
                Self::Float(value as f64)
            }
        }

        impl From<&$val> for Value {
            fn from(value: &$val) -> Self {
                Self::Float(*value as f64)
            }
        }

        impl From<Option<$val>> for Value {
            fn from(value: Option<$val>) -> Self {
                if let Some(v) = value {
                    Self::Float(v as f64)
                } else {
                    Self::Empty
                }
            }
        }
    };
}

param_value_int!(i8);
param_value_int!(i16);
param_value_int!(i32);
param_value_int!(i64);

param_value_int!(u8);
param_value_int!(u16);
param_value_int!(u32);
param_value_int!(u64);
param_value_int!(usize);

param_value_float!(f32);
param_value_float!(f64);

/// A borrowed parameter value that may be a string, a number, or empty. It is intended to
/// be paired with the [`ParamValue`] trait.
///
/// The owned equivalent of this type is [`Value`].
#[derive(Debug, Clone, PartialEq, PartialOrd, Default)]
pub enum ValueRef<'a> {
    /// A text value of arbitrary length
    String(Cow<'a, str>),
    /// A floating point number
    Float(f64),
    /// A integral number
    Int(i64),
    /// Arbitrary binary data
    Buffer(Cow<'a, [u8]>),
    /// No value specified
    #[default]
    Empty,
}

impl<'a> Eq for ValueRef<'a> {}

impl<'a> From<String> for ValueRef<'a> {
    fn from(value: String) -> Self {
        value.parse().unwrap()
    }
}

impl<'a> From<&'a str> for ValueRef<'a> {
    fn from(value: &'a str) -> Self {
        ValueRef::new(value)
    }
}

impl<'a> From<Cow<'a, str>> for ValueRef<'a> {
    fn from(value: Cow<'a, str>) -> Self {
        match value {
            Cow::Borrowed(s) => Self::new(s),
            Cow::Owned(s) => s.parse().unwrap(),
        }
    }
}

impl<'a> PartialEq<String> for ValueRef<'a> {
    fn eq(&self, other: &String) -> bool {
        self.as_str() == other.as_str()
    }
}

impl<'a> PartialEq<str> for ValueRef<'a> {
    fn eq(&self, other: &str) -> bool {
        self.as_str() == other
    }
}

impl<'a> PartialEq<&str> for ValueRef<'a> {
    fn eq(&self, other: &&str) -> bool {
        self.as_str() == *other
    }
}

impl<'a> PartialEq<i64> for ValueRef<'a> {
    fn eq(&self, other: &i64) -> bool {
        if let Self::Int(val) = self {
            val == other
        } else {
            false
        }
    }
}

impl<'a> PartialEq<f64> for ValueRef<'a> {
    fn eq(&self, other: &f64) -> bool {
        if let Self::Float(val) = self {
            val == other
        } else {
            false
        }
    }
}

impl FromStr for ValueRef<'_> {
    type Err = ParamValueParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Ok(Self::Empty);
        }
        if let Ok(value) = s.parse::<i64>() {
            Ok(Self::Int(value))
        } else if let Ok(value) = s.parse::<f64>() {
            Ok(Self::Float(value))
        } else {
            Ok(Self::String(Cow::Owned(s.to_string())))
        }
    }
}

impl<'a> Display for ValueRef<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::String(v) => f.write_str(&v),
            Self::Float(v) => v.fmt(f),
            Self::Int(v) => v.fmt(f),
            Self::Buffer(v) => f.write_str(&String::from_utf8_lossy(v)),
            Self::Empty => f.write_str(""),
        }
    }
}

impl<'a> ValueRef<'a> {
    /// Convert a string value into a precise value type by trying
    /// successive types to parse, defaulting to storing the string
    /// as-is.
    pub fn new(s: &'a str) -> Self {
        if s.is_empty() {
            return Self::Empty;
        }
        if let Ok(value) = s.parse::<i64>() {
            Self::Int(value)
        } else if let Ok(value) = s.parse::<f64>() {
            Self::Float(value)
        } else {
            Self::String(Cow::Borrowed(s))
        }
    }

    /// Create a string [`ValueRef`]
    pub const fn wrap(s: &'a str) -> Self {
        Self::String(Cow::Borrowed(s))
    }

    fn is_empty(&self) -> bool {
        matches!(self, Self::Empty)
    }

    fn is_i64(&self) -> bool {
        matches!(self, Self::Int(_))
    }

    fn is_f64(&self) -> bool {
        matches!(self, Self::Float(_))
    }

    fn is_buffer(&self) -> bool {
        matches!(self, Self::Buffer(_))
    }

    fn is_str(&self) -> bool {
        matches!(self, Self::String(_))
    }

    /// Store the value as a floating point number
    pub fn coerce_f64(&mut self) -> Result<(), ParamValueParseError> {
        let value = self.to_f64()?;
        *self = Self::Float(value);
        Ok(())
    }

    /// Store the value as an integer
    pub fn coerce_i64(&mut self) -> Result<(), ParamValueParseError> {
        let value = self.to_i64()?;
        *self = Self::Int(value);
        Ok(())
    }

    /// Store the value as a string
    pub fn coerce_str(&mut self) -> Result<(), ParamValueParseError> {
        if self.is_str() {
        } else {
            let value = self.to_string();
            *self = Self::String(Cow::Owned(value));
        }
        Ok(())
    }

    /// Discard the value, leaving this value [`ValueRef::Empty`]
    pub fn coerce_empty(&mut self) {
        *self = Self::Empty;
    }

    /// Store the value as a byte buffer
    pub fn coerce_buffer(&mut self) -> Result<(), ParamValueParseError> {
        if self.is_buffer() {
            Ok(())
        } else {
            let buffer = Cow::Owned(self.to_buffer()?.to_vec());
            *self = Self::Buffer(buffer);
            Ok(())
        }
    }

    fn parse<T: FromStr>(&self) -> Result<T, T::Err> {
        match self {
            Self::String(s) => s.parse(),
            Self::Float(v) => v.to_string().parse(),
            Self::Int(i) => i.to_string().parse(),
            Self::Buffer(b) => String::from_utf8_lossy(b).parse(),
            Self::Empty => "".parse(),
        }
    }

    fn to_f64(&self) -> Result<f64, ParamValueParseError> {
        if let Self::Float(val) = self {
            return Ok(*val);
        } else if let Self::Int(val) = self {
            return Ok(*val as f64);
        } else if let Self::String(val) = self {
            if let Ok(v) = val.parse() {
                return Ok(v);
            }
        }
        return Err(ParamValueParseError::FailedToExtractFloat(Some(
            self.to_string(),
        )));
    }

    fn to_i64(&self) -> Result<i64, ParamValueParseError> {
        if let Self::Int(val) = self {
            return Ok(*val);
        } else if let Self::Float(val) = self {
            return Ok(*val as i64);
        } else if let Self::String(val) = self {
            if let Ok(v) = val.parse() {
                return Ok(v);
            }
        }
        return Err(ParamValueParseError::FailedToExtractInt(Some(
            self.to_string(),
        )));
    }

    fn to_str(&self) -> Cow<'_, str> {
        if let Self::String(val) = self {
            Cow::Borrowed(val)
        } else {
            Cow::Owned(self.to_string())
        }
    }

    fn to_buffer(&self) -> Result<Cow<'_, [u8]>, ParamValueParseError> {
        if let Self::Buffer(val) = self {
            match val {
                Cow::Borrowed(v) => Ok(Cow::Borrowed(*v)),
                Cow::Owned(v) => Ok(Cow::Borrowed(v)),
            }
        } else if let Self::String(val) = self {
            Ok(Cow::Borrowed(val.as_bytes()))
        } else {
            Err(ParamValueParseError::FailedToExtractBuffer)
        }
    }
}

impl<'a> ParamValue for ValueRef<'a> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn is_i64(&self) -> bool {
        self.is_i64()
    }

    fn is_f64(&self) -> bool {
        self.is_f64()
    }

    fn is_buffer(&self) -> bool {
        self.is_buffer()
    }

    fn is_str(&self) -> bool {
        self.is_str()
    }

    fn to_f64(&self) -> Result<f64, ParamValueParseError> {
        self.to_f64()
    }

    fn to_i64(&self) -> Result<i64, ParamValueParseError> {
        self.to_i64()
    }

    fn to_str(&self) -> Cow<'_, str> {
        self.to_str()
    }

    fn to_buffer(&self) -> Result<Cow<'_, [u8]>, ParamValueParseError> {
        self.to_buffer()
    }

    fn parse<T: FromStr>(&self) -> Result<T, T::Err> {
        self.parse()
    }

    fn as_bytes(&self) -> Cow<'_, [u8]> {
        match self {
            Self::String(v) => Cow::Borrowed(v.as_bytes()),
            Self::Buffer(v) => Cow::Borrowed(v.as_ref()),
            Self::Float(v) => Cow::Owned(v.to_string().into_bytes()),
            Self::Int(v) => Cow::Owned(v.to_string().into_bytes()),
            Self::Empty => Cow::Borrowed(b""),
        }
    }

    fn as_ref(&self) -> ValueRef<'_> {
        self.clone()
    }
}

impl<'a> From<&'a Value> for ValueRef<'a> {
    fn from(value: &'a Value) -> Self {
        match value {
            Value::String(s) => Self::String(Cow::Borrowed(s)),
            Value::Float(v) => Self::Float(*v),
            Value::Int(v) => Self::Int(*v),
            Value::Buffer(v) => Self::Buffer(Cow::Borrowed(v)),
            Value::Empty => Self::Empty,
        }
    }
}

impl<'a> From<ValueRef<'a>> for Value {
    fn from(value: ValueRef<'a>) -> Self {
        match value {
            ValueRef::String(s) => match s {
                Cow::Borrowed(s) => Self::String(s.to_string()),
                Cow::Owned(s) => Self::String(s),
            },
            ValueRef::Float(v) => Self::Float(v),
            ValueRef::Int(v) => Self::Int(v),
            ValueRef::Buffer(v) => Self::Buffer(v.to_vec().into_boxed_slice()),
            ValueRef::Empty => Self::Empty,
        }
    }
}

macro_rules! param_value_ref_int {
    ($val:ty) => {
        impl<'a> From<$val> for ValueRef<'a> {
            fn from(value: $val) -> Self {
                Self::Int(value as i64)
            }
        }

        impl<'a> From<&$val> for ValueRef<'a> {
            fn from(value: &$val) -> Self {
                Self::Int(*value as i64)
            }
        }

        impl<'a> From<Option<$val>> for ValueRef<'a> {
            fn from(value: Option<$val>) -> Self {
                if let Some(v) = value {
                    Self::Int(v as i64)
                } else {
                    Self::Empty
                }
            }
        }
    };
}

macro_rules! param_value_ref_float {
    ($val:ty) => {
        impl<'a> From<$val> for ValueRef<'a> {
            fn from(value: $val) -> Self {
                Self::Float(value as f64)
            }
        }

        impl<'a> From<&$val> for ValueRef<'a> {
            fn from(value: &$val) -> Self {
                Self::Float(*value as f64)
            }
        }

        impl<'a> From<Option<$val>> for ValueRef<'a> {
            fn from(value: Option<$val>) -> Self {
                if let Some(v) = value {
                    Self::Float(v as f64)
                } else {
                    Self::Empty
                }
            }
        }
    };
}

param_value_ref_int!(i8);
param_value_ref_int!(i16);
param_value_ref_int!(i32);
param_value_ref_int!(i64);

param_value_ref_int!(u8);
param_value_ref_int!(u16);
param_value_ref_int!(u32);
param_value_ref_int!(u64);
param_value_ref_int!(usize);

param_value_ref_float!(f32);
param_value_ref_float!(f64);

/// A CURIE is a namespace + accession identifier
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CURIE {
    pub controlled_vocabulary: ControlledVocabulary,
    pub accession: u32,
}

#[macro_export]
macro_rules! curie {
    (MS:$acc:literal) => {
        $crate::params::CURIE::new($crate::params::ControlledVocabulary::MS, $acc)
    };
    (UO:$acc:literal) => {
        $crate::params::CURIE::new($crate::params::ControlledVocabulary::UO, $acc)
    };
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

impl TryFrom<&Param> for CURIE {
    type Error = String;

    fn try_from(value: &Param) -> Result<Self, Self::Error> {
        if value.is_controlled() {
            Ok(CURIE::new(
                value.controlled_vocabulary.unwrap(),
                value.accession.unwrap(),
            ))
        } else {
            Err(format!(
                "{} does is not a controlled vocabulary term",
                value.name()
            ))
        }
    }
}

impl<'a> TryFrom<&ParamCow<'a>> for CURIE {
    type Error = String;

    fn try_from(value: &ParamCow<'a>) -> Result<Self, Self::Error> {
        if value.is_controlled() {
            Ok(CURIE::new(
                value.controlled_vocabulary.unwrap(),
                value.accession.unwrap(),
            ))
        } else {
            Err(format!(
                "{} does is not a controlled vocabulary term",
                value.name()
            ))
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
    fn value(&self) -> ValueRef;
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
    pub value: ValueRef<'a>,
    pub accession: Option<u32>,
    pub controlled_vocabulary: Option<ControlledVocabulary>,
    pub unit: Unit,
}

impl<'a> ParamValue for ParamCow<'a> {
    fn is_empty(&self) -> bool {
        <ValueRef<'a> as ParamValue>::is_empty(&self.value)
    }

    fn is_i64(&self) -> bool {
        <ValueRef<'a> as ParamValue>::is_i64(&self.value)
    }

    fn is_f64(&self) -> bool {
        <ValueRef<'a> as ParamValue>::is_f64(&self.value)
    }

    fn is_buffer(&self) -> bool {
        <ValueRef<'a> as ParamValue>::is_buffer(&self.value)
    }

    fn is_str(&self) -> bool {
        <ValueRef<'a> as ParamValue>::is_str(&self.value)
    }

    fn to_f64(&self) -> Result<f64, ParamValueParseError> {
        <ValueRef<'a> as ParamValue>::to_f64(&self.value)
    }

    fn to_i64(&self) -> Result<i64, ParamValueParseError> {
        <ValueRef<'a> as ParamValue>::to_i64(&self.value)
    }

    fn to_str(&self) -> Cow<'_, str> {
        <ValueRef<'a> as ParamValue>::to_str(&self.value)
    }

    fn to_buffer(&self) -> Result<Cow<'_, [u8]>, ParamValueParseError> {
        <ValueRef<'a> as ParamValue>::to_buffer(&self.value)
    }

    fn parse<T: FromStr>(&self) -> Result<T, T::Err> {
        <ValueRef<'a> as ParamValue>::parse(&self.value)
    }

    fn as_bytes(&self) -> Cow<'_, [u8]> {
        <ValueRef<'a> as ParamValue>::as_bytes(&self.value)
    }

    fn as_ref(&self) -> ValueRef<'_> {
        <ValueRef<'a> as ParamValue>::as_ref(&self.value)
    }
}

impl ParamCow<'static> {
    pub const fn const_new(
        name: &'static str,
        value: ValueRef<'static>,
        accession: Option<u32>,
        controlled_vocabulary: Option<ControlledVocabulary>,
        unit: Unit,
    ) -> Self {
        Self {
            name: Cow::Borrowed(name),
            value: value,
            accession,
            controlled_vocabulary,
            unit,
        }
    }
}

impl<'a> ParamCow<'a> {
    pub fn new(
        name: Cow<'a, str>,
        value: ValueRef<'a>,
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

    fn value(&self) -> ValueRef<'a> {
        self.value.clone()
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
            value: value.value.into(),
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

impl<'a> AsRef<ValueRef<'a>> for ParamCow<'a> {
    fn as_ref(&self) -> &ValueRef<'a> {
        &self.value
    }
}

/// A controlled vocabulary or user parameter
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Param {
    pub name: String,
    pub value: Value,
    pub accession: Option<u32>,
    pub controlled_vocabulary: Option<ControlledVocabulary>,
    pub unit: Unit,
}

impl AsRef<Value> for Param {
    fn as_ref(&self) -> &Value {
        &self.value
    }
}

impl ParamValue for Param {
    fn is_empty(&self) -> bool {
        <Value as ParamValue>::is_empty(&self.value)
    }

    fn is_i64(&self) -> bool {
        <Value as ParamValue>::is_i64(&self.value)
    }

    fn is_f64(&self) -> bool {
        <Value as ParamValue>::is_f64(&self.value)
    }

    fn is_buffer(&self) -> bool {
        <Value as ParamValue>::is_buffer(&self.value)
    }

    fn is_str(&self) -> bool {
        <Value as ParamValue>::is_str(&self.value)
    }

    fn to_f64(&self) -> Result<f64, ParamValueParseError> {
        <Value as ParamValue>::to_f64(&self.value)
    }

    fn to_i64(&self) -> Result<i64, ParamValueParseError> {
        <Value as ParamValue>::to_i64(&self.value)
    }

    fn to_str(&self) -> Cow<'_, str> {
        <Value as ParamValue>::to_str(&self.value)
    }

    fn to_buffer(&self) -> Result<Cow<'_, [u8]>, ParamValueParseError> {
        <Value as ParamValue>::to_buffer(&self.value)
    }

    fn parse<T: FromStr>(&self) -> Result<T, T::Err> {
        <Value as ParamValue>::parse(&self.value)
    }

    fn as_bytes(&self) -> Cow<'_, [u8]> {
        <Value as ParamValue>::as_bytes(&self.value)
    }

    fn as_ref(&self) -> ValueRef<'_> {
        <Value as ParamValue>::as_ref(&self.value)
    }
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
        inst.value = value.into().into();
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

    fn value(&self) -> ValueRef {
        self.value.as_ref()
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
    Number(u32),
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

    pub fn param<A: Into<AccessionLike<'a>>, S: Into<String>>(
        &self,
        accession: A,
        name: S,
    ) -> Param {
        let mut param = Param::new();
        param.controlled_vocabulary = Some(*self);
        param.name = name.into();

        let accession: AccessionLike = accession.into();

        match accession {
            AccessionLike::Text(s) => {
                if let Some(nb) = s.split(":").last() {
                    param.accession =
                        Some(nb.parse().unwrap_or_else(|_| {
                            panic!("Expected accession to be numeric, got {}", s)
                        }))
                }
            }
            AccessionLike::Number(n) => param.accession = Some(n),
        }
        param
    }

    pub const fn curie(&self, accession: u32) -> CURIE {
        CURIE::new(*self, accession)
    }

    pub const fn const_param(
        &self,
        name: &'static str,
        value: ValueRef<'static>,
        accession: u32,
        unit: Unit,
    ) -> ParamCow<'static> {
        ParamCow {
            name: Cow::Borrowed(name),
            value: value,
            accession: Some(accession),
            controlled_vocabulary: Some(*self),
            unit,
        }
    }

    pub const fn const_param_ident(&self, name: &'static str, accession: u32) -> ParamCow<'static> {
        self.const_param(name, ValueRef::Empty, accession, Unit::Unknown)
    }

    pub const fn const_param_ident_unit(
        &self,
        name: &'static str,
        accession: u32,
        unit: Unit,
    ) -> ParamCow<'static> {
        self.const_param(name, ValueRef::Empty, accession, unit)
    }

    pub fn param_val<S: Into<String>, A: Into<AccessionLike<'a>>, V: Into<Value>>(
        &self,
        accession: A,
        name: S,
        value: V,
    ) -> Param {
        let mut param = self.param(accession, name);
        param.value = value.into();
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

    pub const fn from_name(name: &str) -> Unit {
        let bytes = name.as_bytes();
        match bytes {
            b"millisecond" => Self::Millisecond,
            b"second" => Self::Second,
            b"minute" => Self::Minute,

            b"m/z" => Self::MZ,
            b"dalton" => Self::Mass,

            b"number of detector counts" => Self::DetectorCounts,
            b"percent of base peak" => Self::PercentBasePeak,
            b"percent of base peak times 100" => Self::PercentBasePeakTimes100,
            b"absorbance unit" => Self::AbsorbanceUnit,
            b"counts per second" => Self::CountsPerSecond,

            b"electronvolt" => Self::Electronvolt,
            b"percent" => Self::PercentElectronVolt,
            _ => Unit::Unknown,
        }
    }

    pub const fn from_accession(acc: &str) -> Unit {
        let bytes = acc.as_bytes();
        match bytes {
            b"UO:0000028" => Self::Millisecond,
            b"UO:0000010" => Self::Second,
            b"UO:0000031" => Self::Minute,

            b"MS:1000040" => Self::MZ,
            b"UO:000221" => Self::Mass,

            b"MS:1000131" => Self::DetectorCounts,
            b"MS:1000132" => Self::PercentBasePeak,
            b"MS:1000905" => Self::PercentBasePeakTimes100,
            b"UO:0000269" => Self::AbsorbanceUnit,
            b"MS:1000814" => Self::CountsPerSecond,

            b"UO:0000266" => Self::Electronvolt,
            b"UO:0000187" => Self::PercentElectronVolt,
            _ => Unit::Unknown,
        }
    }

    pub const fn from_curie(acc: &CURIE) -> Unit {
        match acc {
            CURIE {
                controlled_vocabulary: ControlledVocabulary::UO,
                accession: 28,
            } => Self::Millisecond,
            CURIE {
                controlled_vocabulary: ControlledVocabulary::UO,
                accession: 10,
            } => Self::Second,
            CURIE {
                controlled_vocabulary: ControlledVocabulary::UO,
                accession: 31,
            } => Self::Minute,

            CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1000040,
            } => Self::MZ,
            CURIE {
                controlled_vocabulary: ControlledVocabulary::UO,
                accession: 221,
            } => Self::Mass,

            CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1000131,
            } => Self::DetectorCounts,
            CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1000132,
            } => Self::PercentBasePeak,
            CURIE {
                controlled_vocabulary: ControlledVocabulary::UO,
                accession: 269,
            } => Self::AbsorbanceUnit,
            CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1000814,
            } => Self::CountsPerSecond,

            CURIE {
                controlled_vocabulary: ControlledVocabulary::UO,
                accession: 266,
            } => Self::Electronvolt,
            CURIE {
                controlled_vocabulary: ControlledVocabulary::UO,
                accession: 187,
            } => Self::PercentElectronVolt,

            _ => Unit::Unknown,
        }
    }

    pub const fn to_curie(&self) -> Option<CURIE> {
        match self {
            Self::Millisecond => Some(CURIE {
                controlled_vocabulary: ControlledVocabulary::UO,
                accession: 28,
            }),
            Self::Second => Some(CURIE {
                controlled_vocabulary: ControlledVocabulary::UO,
                accession: 10,
            }),
            Self::Minute => Some(CURIE {
                controlled_vocabulary: ControlledVocabulary::UO,
                accession: 31,
            }),

            Self::MZ => Some(CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1000040,
            }),
            Self::Mass => Some(CURIE {
                controlled_vocabulary: ControlledVocabulary::UO,
                accession: 221,
            }),

            Self::DetectorCounts => Some(CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1000131,
            }),
            Self::PercentBasePeak => Some(CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1000132,
            }),
            Self::AbsorbanceUnit => Some(CURIE {
                controlled_vocabulary: ControlledVocabulary::UO,
                accession: 269,
            }),
            Self::CountsPerSecond => Some(CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1000814,
            }),

            Self::Electronvolt => Some(CURIE {
                controlled_vocabulary: ControlledVocabulary::UO,
                accession: 266,
            }),
            Self::PercentElectronVolt => Some(CURIE {
                controlled_vocabulary: ControlledVocabulary::UO,
                accession: 187,
            }),

            _ => None,
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
