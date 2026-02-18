//! Elements of controlled vocabularies used to describe mass spectra and their components.
//!
//! Directly maps to the usage of the PSI-MS controlled vocabulary in mzML
use std::borrow::Cow;
use std::convert::TryFrom;
use std::fmt::Display;
use std::hash::Hash;
use std::str::{self, FromStr};
use std::{io, mem, num};

use thiserror::Error;

/// An owned parameter value that may be a string, a number, or empty. It is intended to
/// be paired with the [`ParamValue`] trait.
///
/// The borrowed equivalent of this type is [`ValueRef`].
#[derive(Debug, Clone, PartialEq, PartialOrd, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Value {
    /// A text value of arbitrary length
    String(String),
    /// A floating point number
    Float(f64),
    /// A integral number
    Int(i64),
    /// Arbitrary binary data
    Buffer(Box<[u8]>),
    /// true/false value
    Boolean(bool),
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

    /// Check if the value is a boolean
    fn is_boolean(&self) -> bool;

    /// Get the value as an `f64`, if possible
    fn to_f64(&self) -> Result<f64, ParamValueParseError>;

    /// Get the value as an `f32`, if possible
    fn to_f32(&self) -> Result<f32, ParamValueParseError> {
        let v = self.to_f64()?;
        Ok(v as f32)
    }

    /// Get the value as a `bool`, if possible
    fn to_bool(&self) -> Result<bool, ParamValueParseError>;

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

    /// Get the size of the stored data type
    fn data_len(&self) -> usize;
}

/// Errors that might occur while trying to convert to a particular value type
/// from whatever type might be stored in a [`ParamValue`]-like object.
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

/// A [`Value`] can be parsed from a string infallibly, even though an error type
/// is given, but unless one of the numerical parsers succeeds, the stored value will
/// just be a string.
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
        } else if let Ok(value) = s.parse::<bool>() {
            Ok(Self::Boolean(value))
        } else {
            Ok(Self::String(s.to_string()))
        }
    }
}

impl Display for Value {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Value::String(v) => f.write_str(v),
            Value::Float(v) => v.fmt(f),
            Value::Int(v) => v.fmt(f),
            Value::Buffer(v) => f.write_str(&String::from_utf8_lossy(v)),
            Value::Empty => f.write_str(""),
            Value::Boolean(v) => v.fmt(f),
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
        } else if let Ok(value) = s.parse::<bool>() {
            Self::Boolean(value)
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

    pub fn coerce_bool(&mut self) -> Result<(), ParamValueParseError> {
        let value = self.to_bool()?;
        *self = Self::Boolean(value);
        Ok(())
    }

    /// Convert the value to a text string and then try to parse it into `T`.
    /// If the type is one of the common numeric types, prefer one of the provided
    /// methods with a `to_` prefix as they avoid the string conversions.
    pub fn parse<T: FromStr>(&self) -> Result<T, T::Err> {
        match self {
            Value::String(s) => s.parse(),
            Value::Float(v) => v.to_string().parse(),
            Value::Int(i) => i.to_string().parse(),
            Value::Buffer(b) => String::from_utf8_lossy(b).parse(),
            Value::Empty => "".parse(),
            Value::Boolean(b) => b.to_string().parse(),
        }
    }

    fn to_bool(&self) -> Result<bool, ParamValueParseError> {
        if let Self::Boolean(val) = self {
            Ok(*val)
        } else if self.is_numeric() {
            Ok(self.to_i64()? != 0)
        } else if let Self::Empty = self {
            Ok(false)
        } else if let Ok(v) = self.parse() {
            Ok(v)
        } else {
            Err(ParamValueParseError::FailedToExtractInt(Some(
                self.to_string(),
            )))
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
        Err(ParamValueParseError::FailedToExtractFloat(Some(
            self.to_string(),
        )))
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
        Err(ParamValueParseError::FailedToExtractInt(Some(
            self.to_string(),
        )))
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
            Self::Boolean(v) => Cow::Owned(v.to_string().into_bytes()),
        }
    }

    fn as_ref(&self) -> ValueRef<'_> {
        self.into()
    }

    fn data_len(&self) -> usize {
        match self {
            Self::String(v) => v.len(),
            Self::Buffer(v) => v.len(),
            Self::Float(_) => 8,
            Self::Int(_) => 8,
            Self::Empty => 0,
            Self::Boolean(_) => mem::size_of::<bool>(),
        }
    }

    fn is_boolean(&self) -> bool {
        matches!(self, Self::Boolean(_))
    }

    fn to_bool(&self) -> Result<bool, ParamValueParseError> {
        self.to_bool()
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

impl PartialEq<bool> for Value {
    fn eq(&self, other: &bool) -> bool {
        if let Self::Boolean(val) = self {
            val == other
        } else {
            false
        }
    }
}

impl Hash for Value {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        core::mem::discriminant(self).hash(state);
        match self {
            Self::String(s) => s.hash(state),
            Self::Float(v) => v.to_bits().hash(state),
            Self::Int(v) => (*v).hash(state),
            Self::Buffer(v) => v.hash(state),
            Self::Empty => 0u8.hash(state),
            Self::Boolean(v) => v.hash(state),
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
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
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
    /// A true/false value
    Boolean(bool),
}

impl Eq for ValueRef<'_> {}

impl From<String> for ValueRef<'_> {
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

impl PartialEq<String> for ValueRef<'_> {
    fn eq(&self, other: &String) -> bool {
        self.as_str() == other.as_str()
    }
}

impl PartialEq<str> for ValueRef<'_> {
    fn eq(&self, other: &str) -> bool {
        self.as_str() == other
    }
}

impl PartialEq<&str> for ValueRef<'_> {
    fn eq(&self, other: &&str) -> bool {
        self.as_str() == *other
    }
}

impl PartialEq<i64> for ValueRef<'_> {
    fn eq(&self, other: &i64) -> bool {
        if let Self::Int(val) = self {
            val == other
        } else {
            false
        }
    }
}

impl PartialEq<f64> for ValueRef<'_> {
    fn eq(&self, other: &f64) -> bool {
        if let Self::Float(val) = self {
            val == other
        } else {
            false
        }
    }
}

impl PartialEq<bool> for ValueRef<'_> {
    fn eq(&self, other: &bool) -> bool {
        if let Self::Boolean(val) = self {
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
        if let Ok(value) = s.parse() {
            Ok(Self::Int(value))
        } else if let Ok(value) = s.parse() {
            Ok(Self::Float(value))
        } else if let Ok(value) = s.parse() {
            Ok(Self::Boolean(value))
        } else {
            Ok(Self::String(Cow::Owned(s.to_string())))
        }
    }
}

impl Display for ValueRef<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::String(v) => f.write_str(v),
            Self::Float(v) => v.fmt(f),
            Self::Int(v) => v.fmt(f),
            Self::Buffer(v) => f.write_str(&String::from_utf8_lossy(v)),
            Self::Empty => f.write_str(""),
            Self::Boolean(v) => v.fmt(f),
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
        } else if let Ok(value) = s.parse() {
            Self::Boolean(value)
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

    fn to_bool(&self) -> Result<bool, ParamValueParseError> {
        if let Self::Boolean(val) = self {
            Ok(*val)
        } else if self.is_numeric() {
            Ok(self.to_i64()? != 0)
        } else if let Self::Empty = self {
            Ok(false)
        } else if let Ok(v) = self.parse() {
            Ok(v)
        } else {
            Err(ParamValueParseError::FailedToExtractInt(Some(
                self.to_string(),
            )))
        }
    }

    pub fn coerce_bool(&mut self) -> Result<(), ParamValueParseError> {
        let value = self.to_bool()?;
        *self = Self::Boolean(value);
        Ok(())
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
            Self::Boolean(v) => v.to_string().parse(),
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
        Err(ParamValueParseError::FailedToExtractFloat(Some(
            self.to_string(),
        )))
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
        Err(ParamValueParseError::FailedToExtractInt(Some(
            self.to_string(),
        )))
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

impl ParamValue for ValueRef<'_> {
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

    fn is_boolean(&self) -> bool {
        matches!(self, Self::Boolean(_))
    }

    fn to_bool(&self) -> Result<bool, ParamValueParseError> {
        self.to_bool()
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
            Self::Boolean(v) => Cow::Owned(v.to_string().into_bytes()),
        }
    }

    fn as_ref(&self) -> ValueRef<'_> {
        self.clone()
    }

    fn data_len(&self) -> usize {
        match self {
            Self::String(v) => v.len(),
            Self::Buffer(v) => v.len(),
            Self::Float(_) => 8,
            Self::Int(_) => 8,
            Self::Empty => 0,
            Self::Boolean(_) => mem::size_of::<bool>(),
        }
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
            Value::Boolean(v) => Self::Boolean(*v),
        }
    }
}

impl PartialEq<Value> for ValueRef<'_> {
    fn eq(&self, other: &Value) -> bool {
        *self == other.as_ref()
    }
}

impl<'a> PartialEq<ValueRef<'a>> for Value {
    fn eq(&self, other: &ValueRef<'a>) -> bool {
        self.as_ref() == *other
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
            ValueRef::Boolean(v) => Self::Boolean(v),
        }
    }
}

impl Hash for ValueRef<'_> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        core::mem::discriminant(self).hash(state);
        match self {
            Self::String(s) => s.hash(state),
            Self::Float(v) => v.to_bits().hash(state),
            Self::Int(v) => (*v).hash(state),
            Self::Buffer(v) => v.hash(state),
            Self::Empty => 0u8.hash(state),
            Self::Boolean(v) => (*v).hash(state),
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

impl From<ValueRef<'_>> for f32 {
    fn from(value: ValueRef<'_>) -> Self {
        value.to_f32().unwrap()
    }
}

impl From<ValueRef<'_>> for f64 {
    fn from(value: ValueRef<'_>) -> Self {
        value.to_f64().unwrap()
    }
}

impl From<ValueRef<'_>> for i32 {
    fn from(value: ValueRef<'_>) -> Self {
        value.to_i32().unwrap()
    }
}

impl From<ValueRef<'_>> for i64 {
    fn from(value: ValueRef<'_>) -> Self {
        value.to_i64().unwrap()
    }
}

impl From<bool> for Value {
    fn from(value: bool) -> Self {
        Self::Boolean(value)
    }
}

impl From<bool> for ValueRef<'_> {
    fn from(value: bool) -> Self {
        Self::Boolean(value)
    }
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

#[cfg(feature = "serde")]
impl From<Value> for serde_json::Value {
    fn from(value: Value) -> Self {
        match value {
            Value::Boolean(val) => serde_json::Value::Bool(val),
            Value::Float(val) => {
                serde_json::Value::Number(serde_json::Number::from_f64(val).unwrap())
            }
            Value::Int(val) => {
                serde_json::Value::Number(serde_json::Number::from_i128(val as i128).unwrap())
            }
            Value::String(val) => serde_json::Value::String(val),
            Value::Buffer(val) => serde_json::to_value(&val).unwrap(),
            Value::Empty => serde_json::Value::Null,
        }
    }
}

pub type AccessionIntCode = u32;
pub type AccessionByteCode7 = [u8; 7];

#[allow(unused)]
#[derive(Debug)]
#[repr(u8)]
enum AccessionCode {
    Int(AccessionIntCode),
    Byte7(AccessionByteCode7),
}

impl Display for AccessionCode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AccessionCode::Int(v) => write!(f, "{v:07}"),
            AccessionCode::Byte7(v) => write!(
                f,
                "{}",
                core::str::from_utf8(v)
                    .map_err(|e| format!("ERROR:{e}"))
                    .unwrap()
            ),
        }
    }
}

#[derive(Debug, thiserror::Error, Clone, PartialEq)]
pub enum AccessionCodeParseError {
    #[error("The acccession code was too long: {0}")]
    AccessionCodeTooLong(String),
    #[error("The acccession code was not in range: {0}")]
    AccessionCodeNotInRange(String),
}

impl FromStr for AccessionCode {
    type Err = AccessionCodeParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() > 7 {
            return Err(AccessionCodeParseError::AccessionCodeTooLong(s.to_string()));
        }
        if !s.is_ascii() {
            return Err(AccessionCodeParseError::AccessionCodeNotInRange(
                s.to_string(),
            ));
        }
        if let Ok(u) = s.parse::<AccessionIntCode>() {
            Ok(Self::Int(u))
        } else {
            let mut bytes = AccessionByteCode7::default();
            for (byte_from, byte_to) in s.as_bytes().iter().rev().zip(bytes.iter_mut().rev()) {
                *byte_to = *byte_from;
            }
            Ok(Self::Byte7(bytes))
        }
    }
}

/// A CURIE is a namespace + accession identifier
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CURIE {
    pub controlled_vocabulary: ControlledVocabulary,
    pub accession: AccessionIntCode,
}

#[allow(unused)]
macro_rules! accessioncode {
    ($acc:literal) => {
        AccessionCode::Int($acc)
    };
    ($acc:tt) => {
        match stringify!($acc).as_bytes() {
            [a, b, c, d, e, f, g] => AccessionCode::Byte7([*a, *b, *c, *d, *e, *f, *g]),
            _ => panic!(concat!(
                "Cannot convert ",
                stringify!($acc),
                " to accession code. Expected exactly 7 bytes"
            )),
        }
    };
}

#[macro_export]
macro_rules! find_param_method {
    ($meth:ident, $curie:expr) => {
        $crate::find_param_method!($meth, $curie, "Find a parameter by its CURIE");
    };
    ($meth:ident, $curie:expr, $desc:literal) => {
        #[doc=$desc]
        pub fn $meth(&self) -> Option<$crate::params::ValueRef<'_>> {
            self.get_param_by_curie($curie)
                .map(|p| $crate::params::ParamLike::value(p))
        }
    };
    ($meth:ident, $curie:expr, $conv:expr, $result:ty) => {
        $crate::find_param_method!(
            $meth,
            $curie,
            $conv,
            $result,
            "Find a parameter by its CURIE"
        );
    };
    ($meth:ident, $curie:expr, $conv:expr, $result:ty, $desc:literal) => {
        #[doc=$desc]
        pub fn $meth(&self) -> $result {
            self.get_param_by_curie($curie).map($conv)
        }
    };
}

#[macro_export]
macro_rules! curie {
    ($ns:ident:$acc:literal) => {
        $crate::params::CURIE {
            controlled_vocabulary: $crate::params::ControlledVocabulary::$ns,
            accession: $acc,
        }
    };
    (IMZML:$acc:literal) => {
        $crate::params::CURIE { controlled_vocabulary: $crate::params::ControlledVocabulary::IMZML, accession: $acc }
    };

}

impl CURIE {
    pub const fn new(cv_id: ControlledVocabulary, accession: AccessionIntCode) -> Self {
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

    #[inline(always)]
    pub fn accession_int(&self) -> u32 {
        self.accession
    }

    #[inline(always)]
    pub fn controlled_vocabulary(&self) -> ControlledVocabulary {
        self.controlled_vocabulary
    }
}

/// A CURIE is a namespace + accession identifier
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct PackedCURIE(u64);

#[allow(unused)]
impl PackedCURIE {
    pub fn from_curie(curie: CURIE) -> Self {
        let cv = curie.controlled_vocabulary();
        let cv_id = cv as u64;
        let acc_code = curie.accession_int();
        let code = acc_code as u64;
        Self(cv_id << 56 | code)
    }

    pub fn accession(&self) -> u64 {
        self.0 & 0x00ffffffffffffff
    }

    pub fn controlled_vocabulary(&self) -> ControlledVocabulary {
        ((self.0 & 0xff00000000000000) as u8).try_into().unwrap()
    }
}

impl Display for CURIE {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{:07}",
            self.controlled_vocabulary.prefix(),
            self.accession
        )
    }
}

impl<T: ParamLike> PartialEq<T> for CURIE {
    fn eq(&self, other: &T) -> bool {
        if !other.is_controlled()
            || other
                .controlled_vocabulary()
                .map(|c| c != self.controlled_vocabulary)
                .unwrap_or_default()
        {
            false
        } else {
            other
                .accession()
                .map(|a| a == self.accession)
                .unwrap_or_default()
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
        let cv = tokens
            .next()
            .ok_or(CURIEParsingError::MissingNamespaceSeparator)?;
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
        match (value.controlled_vocabulary, value.accession) {
            (Some(cv), Some(acc)) => Ok(CURIE::new(cv, acc)),
            _ => Err(format!(
                "{} is missing controlled vocabulary or accession",
                value.name()
            )),
        }
    }
}

impl<'a> TryFrom<&ParamCow<'a>> for CURIE {
    type Error = String;

    fn try_from(value: &ParamCow<'a>) -> Result<Self, Self::Error> {
        match (value.controlled_vocabulary, value.accession) {
            (Some(cv), Some(acc)) => Ok(CURIE::new(cv, acc)),
            _ => Err(format!(
                "{} is missing controlled vocabulary or accession",
                value.name()
            )),
        }
    }
}

pub fn curie_to_num(curie: &str) -> (Option<ControlledVocabulary>, Option<AccessionIntCode>) {
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
    fn value(&self) -> ValueRef<'_>;
    fn accession(&self) -> Option<AccessionIntCode>;
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
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct ParamCow<'a> {
    pub name: Cow<'a, str>,
    pub value: ValueRef<'a>,
    pub accession: Option<AccessionIntCode>,
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

    fn data_len(&self) -> usize {
        <ValueRef<'a> as ParamValue>::data_len(&self.value)
    }

    fn is_boolean(&self) -> bool {
        <ValueRef<'a> as ParamValue>::is_boolean(&self.value)
    }

    fn to_bool(&self) -> Result<bool, ParamValueParseError> {
        <ValueRef<'a> as ParamValue>::to_bool(&self.value)
    }
}

impl ParamCow<'static> {
    pub const fn const_new(
        name: &'static str,
        value: ValueRef<'static>,
        accession: Option<AccessionIntCode>,
        controlled_vocabulary: Option<ControlledVocabulary>,
        unit: Unit,
    ) -> Self {
        Self {
            name: Cow::Borrowed(name),
            value,
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
        accession: Option<AccessionIntCode>,
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

    pub const fn is_controlled(&self) -> bool {
        self.accession.is_some()
    }

    pub const fn curie(&self) -> Option<CURIE> {
        match (self.controlled_vocabulary, self.accession) {
            (Some(cv), Some(acc)) => Some(CURIE::new(cv, acc)),
            _ => None,
        }
    }
}

impl<'a> ParamLike for ParamCow<'a> {
    fn name(&self) -> &str {
        &self.name
    }

    fn value(&self) -> ValueRef<'a> {
        self.value.clone()
    }

    fn accession(&self) -> Option<AccessionIntCode> {
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

impl PartialEq<CURIE> for ParamCow<'_> {
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Param {
    pub name: String,
    pub value: Value,
    pub accession: Option<AccessionIntCode>,
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

    fn data_len(&self) -> usize {
        <Value as ParamValue>::data_len(&self.value)
    }

    fn is_boolean(&self) -> bool {
        <Value as ParamValue>::is_boolean(&self.value)
    }

    fn to_bool(&self) -> Result<bool, ParamValueParseError> {
        <Value as ParamValue>::to_bool(&self.value)
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

/// Incrementally build up a [`Param`]
#[derive(Default, Debug, Clone)]
pub struct ParamBuilder {
    name: String,
    value: Value,
    accession: Option<AccessionIntCode>,
    controlled_vocabulary: Option<ControlledVocabulary>,
    unit: Unit,
}

impl ParamBuilder {
    pub fn name<S: ToString>(mut self, name: S) -> Self {
        self.name = name.to_string();
        self
    }

    pub fn value<V: Into<Value>>(mut self, value: V) -> Self {
        self.value = value.into();
        self
    }

    pub fn controlled_vocabulary(mut self, cv: ControlledVocabulary) -> Self {
        self.controlled_vocabulary = Some(cv);
        self
    }

    pub fn accession(mut self, accession: AccessionIntCode) -> Self {
        self.accession = Some(accession);
        self
    }

    /// A convenience method that configures the controlled vocabulary and accession number
    /// from a [`CURIE`]
    pub fn curie(mut self, curie: CURIE) -> Self {
        self.controlled_vocabulary = Some(curie.controlled_vocabulary);
        self.accession = Some(curie.accession);
        self
    }

    pub fn unit(mut self, unit: Unit) -> Self {
        self.unit = unit;
        self
    }

    /// Consume the builder to produce a [`Param`]
    pub fn build(self) -> Param {
        let mut this = Param::new();
        this.name = self.name;
        this.value = self.value;
        this.controlled_vocabulary = self.controlled_vocabulary;
        this.accession = self.accession;
        this.unit = self.unit;
        this
    }
}

impl Param {
    /// Create a new, empty [`Param`].
    ///
    /// # See also
    /// - [`Param::new_key_value`]
    /// - [`Param::builder`]
    /// - [`ControlledVocabulary::param`]
    /// - [`ControlledVocabulary::param_val`]
    pub fn new() -> Param {
        Param {
            ..Default::default()
        }
    }

    /// Create a new [`ParamBuilder`] to make creating a new [`Param`] more convenient.
    pub fn builder() -> ParamBuilder {
        ParamBuilder::default()
    }

    /// A construction method for [`Param`] that sets [`Param::name`] and [`Param::value`]
    /// but leaves all other attributes as default.
    pub fn new_key_value<K: Into<String>, V: Into<Value>>(name: K, value: V) -> Param {
        let mut inst = Self::new();
        inst.name = name.into();
        inst.value = value.into();
        inst
    }

    /// Attempt to parse the value of this parameter into `T`.
    ///
    /// See [`Value::parse`]
    pub fn parse<T: str::FromStr>(&self) -> Result<T, T::Err> {
        self.value.parse::<T>()
    }

    /// Check if this parameter is defined in a controlled vocabulary
    pub const fn is_controlled(&self) -> bool {
        self.accession.is_some()
    }

    /// Create a [`CURIE`] from [`Param::controlled_vocabulary`] and [`Param::accession`]
    pub const fn curie(&self) -> Option<CURIE> {
        match (self.controlled_vocabulary, self.accession) {
            (Some(cv), Some(acc)) => Some(CURIE::new(cv, acc)),
            _ => None,
        }
    }

    /// Format the [`Param::curie`] as a string, if it exists
    pub fn curie_str(&self) -> Option<String> {
        self.curie().map(|c| c.to_string())
    }

    /// Update [`Param::unit`] inferred from `accession`, failing that, `name`
    pub fn with_unit<S: AsRef<str>, A: AsRef<str>>(mut self, accession: S, name: A) -> Param {
        self.unit = Unit::from_accession(accession.as_ref());
        if matches!(self.unit, Unit::Unknown) {
            self.unit = Unit::from_name(name.as_ref());
        }
        self
    }

    /// Update [`Param::unit`] from `unit`
    pub fn with_unit_t(mut self, unit: &Unit) -> Param {
        self.unit = *unit;
        self
    }
}

impl ParamLike for Param {
    fn name(&self) -> &str {
        &self.name
    }

    fn value(&self) -> ValueRef<'_> {
        self.value.as_ref()
    }

    fn accession(&self) -> Option<AccessionIntCode> {
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

impl<'a> PartialEq<ParamCow<'a>> for Param {
    fn eq(&self, other: &ParamCow<'a>) -> bool {
        self.controlled_vocabulary == other.controlled_vocabulary
            && self.accession == other.accession
            && self.name == other.name
            && self.value == other.value
            && self.unit == other.unit
    }
}

impl PartialEq<Param> for ParamCow<'_> {
    fn eq(&self, other: &Param) -> bool {
        self.controlled_vocabulary == other.controlled_vocabulary
            && self.accession == other.accession
            && self.name == other.name
            && self.value == other.value
            && self.unit == other.unit
    }
}

impl Hash for Param {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.name.hash(state);
        self.value.hash(state);
        self.accession.hash(state);
        self.controlled_vocabulary.hash(state);
        self.unit.hash(state);
    }
}

/// Controlled vocabularies used in mass spectrometry data files
#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(u8)]
pub enum ControlledVocabulary {
    /// The PSI-MS Controlled Vocabulary [https://www.ebi.ac.uk/ols4/ontologies/ms](https://www.ebi.ac.uk/ols4/ontologies/ms)
    MS = 1,
    /// The Unit Ontology [https://www.ebi.ac.uk/ols4/ontologies/uo](https://www.ebi.ac.uk/ols4/ontologies/uo)
    UO,
    /// The Experimental Factor Ontology <https://www.ebi.ac.uk/ols4/ontologies/efo>
    EFO,
    /// The Ontology for Biomedical Investigations <https://www.ebi.ac.uk/ols4/ontologies/obi>
    OBI,
    /// The Human Ancestry Ontology <https://www.ebi.ac.uk/ols4/ontologies/hancestro>
    HANCESTRO,
    /// The Basic Formal Ontology <https://www.ebi.ac.uk/ols4/ontologies/bfo>
    BFO,
    /// The NCI Thesaurus OBO Edition <https://www.ebi.ac.uk/ols4/ontologies/ncit>
    NCIT,
    /// The BRENDA Tissue Ontology <https://www.ebi.ac.uk/ols4/ontologies/bto>
    BTO,
    /// The PRIDE Controlled Vocabulary <https://www.ebi.ac.uk/ols4/ontologies/pride>
    PRIDE,
    /// The Imaging MS Controlled Vocabulary <https://www.ms-imaging.org/imzml/>
    #[cfg(feature = "imzml")]
    IMS,
    Unknown,
}

const MS_CV: &str = "MS";
const UO_CV: &str = "UO";
const EFO_CV: &str = "EFO";
const OBI_CV: &str = "OBI";
const HANCESTRO_CV: &str = "HANCESTRO";
const BFO_CV: &str = "BFO";
const BTO_CV: &str = "BTO";
const NCIT_CV: &str = "NCIT";
const PRIDE_CV: &str = "PRIDE";
#[cfg(feature = "imzml")]
const IMS_CV: &str = "IMS";

const MS_CV_BYTES: &[u8] = MS_CV.as_bytes();
const UO_CV_BYTES: &[u8] = UO_CV.as_bytes();
const EFO_CV_BYTES: &[u8] = EFO_CV.as_bytes();
const OBI_CV_BYTES: &[u8] = OBI_CV.as_bytes();
const HANCESTRO_CV_BYTES: &[u8] = HANCESTRO_CV.as_bytes();
const BFO_CV_BYTES: &[u8] = BFO_CV.as_bytes();
const BTO_CV_BYTES: &[u8] = BTO_CV.as_bytes();
const NCIT_CV_BYTES: &[u8] = NCIT_CV.as_bytes();
const PRIDE_CV_BYTES: &[u8] = PRIDE_CV.as_bytes();
#[cfg(feature = "imzml")]
const IMS_CV_BYTES: &[u8] = IMS_CV.as_bytes();

impl TryFrom<u8> for ControlledVocabulary {
    type Error = ControlledVocabularyResolutionError;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            1 => Ok(Self::MS),
            2 => Ok(Self::UO),
            3 => Ok(Self::EFO),
            4 => Ok(Self::OBI),
            5 => Ok(Self::HANCESTRO),
            6 => Ok(Self::BFO),
            7 => Ok(Self::NCIT),
            8 => Ok(Self::BTO),
            9 => Ok(Self::PRIDE),
            #[cfg(feature = "imzml")]
            10 => Ok(Self::IMS),
            _ => Err(ControlledVocabularyResolutionError::UnknownControlledVocabularyCode(value)),
        }
    }
}

/// Anything that can be converted into an accession code portion of a [`CURIE`]
#[derive(Debug, Clone)]
pub enum AccessionLike<'a> {
    Text(Cow<'a, str>),
    Number(AccessionIntCode),
    CURIE(CURIE),
}

impl From<AccessionIntCode> for AccessionLike<'_> {
    fn from(value: AccessionIntCode) -> Self {
        Self::Number(value)
    }
}

impl<'a> From<&'a str> for AccessionLike<'a> {
    fn from(value: &'a str) -> Self {
        Self::Text(Cow::Borrowed(value))
    }
}

impl From<String> for AccessionLike<'_> {
    fn from(value: String) -> Self {
        Self::Text(Cow::Owned(value))
    }
}

impl<'a> ControlledVocabulary {
    /// Get the CURIE namespace prefix for this controlled vocabulary
    pub const fn prefix(&self) -> Cow<'static, str> {
        match &self {
            Self::MS => Cow::Borrowed(MS_CV),
            Self::UO => Cow::Borrowed(UO_CV),
            Self::EFO => Cow::Borrowed(EFO_CV),
            Self::OBI => Cow::Borrowed(OBI_CV),
            Self::HANCESTRO => Cow::Borrowed(HANCESTRO_CV),
            Self::BFO => Cow::Borrowed(BFO_CV),
            Self::NCIT => Cow::Borrowed(NCIT_CV),
            Self::BTO => Cow::Borrowed(BTO_CV),
            Self::PRIDE => Cow::Borrowed(PRIDE_CV),
            #[cfg(feature = "imzml")]
            Self::IMS => Cow::Borrowed(IMS_CV),
            Self::Unknown => panic!("Cannot encode unknown CV"),
        }
    }

    /// Like [`ControlledVocabulary::prefix`], but obtain a byte string instead
    pub const fn as_bytes(&self) -> &'static [u8] {
        match &self {
            Self::MS => MS_CV_BYTES,
            Self::UO => UO_CV_BYTES,
            Self::EFO => EFO_CV_BYTES,
            Self::OBI => OBI_CV_BYTES,
            Self::HANCESTRO => HANCESTRO_CV_BYTES,
            Self::BFO => BFO_CV_BYTES,
            Self::NCIT => NCIT_CV_BYTES,
            Self::BTO => BTO_CV_BYTES,
            Self::PRIDE => PRIDE_CV_BYTES,
            #[cfg(feature = "imzml")]
            Self::IMS => IMS_CV_BYTES,
            Self::Unknown => panic!("Cannot encode unknown CV"),
        }
    }

    pub const fn as_option(&self) -> Option<Self> {
        match self {
            Self::Unknown => None,
            _ => Some(*self),
        }
    }

    /// Create a [`Param`] whose accession comes from this controlled vocabulary namespace with
    /// an empty value.
    ///
    /// # Arguments
    /// - `accession`: The accession code for the [`Param`]. If specified as a [`CURIE`] or a string-like type,
    ///   any namespace is ignored.
    /// - `name`: The name of the parameter
    /// # See Also
    /// - [`ControlledVocabulary::param_val`]
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
                if let Some(nb) = s.split(':').next_back() {
                    param.accession =
                        Some(nb.parse().unwrap_or_else(|_| {
                            panic!("Expected accession to be numeric, got {}", s)
                        }))
                }
            }
            AccessionLike::Number(n) => param.accession = Some(n),
            AccessionLike::CURIE(c) => param.accession = Some(c.accession),
        }
        param
    }

    pub const fn curie(&self, accession: AccessionIntCode) -> CURIE {
        CURIE::new(*self, accession)
    }

    /// Create a [`ParamCow`] from this namespace in a `const` context, useful for preparing
    /// global constants or inlined variables.
    ///
    /// All parameters must have a `'static` lifetime.
    ///
    /// # Arguments
    /// - `name`: The name of the controlled vocabulary term.
    /// - `value`: The wrapped value as a constant.
    /// - `accession`: The a priori determined accession code for the term
    /// - `unit`: The unit associated with the value
    pub const fn const_param(
        &self,
        name: &'static str,
        value: ValueRef<'static>,
        accession: AccessionIntCode,
        unit: Unit,
    ) -> ParamCow<'static> {
        ParamCow {
            name: Cow::Borrowed(name),
            value,
            accession: Some(accession),
            controlled_vocabulary: Some(*self),
            unit,
        }
    }

    /// Create a [`ParamCow`] from this namespace in a `const` context with an empty
    /// value and no unit.
    ///
    /// See [`ControlledVocabulary::const_param`] for more details.
    pub const fn const_param_ident(
        &self,
        name: &'static str,
        accession: AccessionIntCode,
    ) -> ParamCow<'static> {
        self.const_param(name, ValueRef::Empty, accession, Unit::Unknown)
    }

    /// Create a [`ParamCow`] from this namespace in a `const` context with an empty
    /// value but a specified unit.
    ///
    /// This is intended to create a "template" that will be copied and have a value specified.
    ///
    /// See [`ControlledVocabulary::const_param`] for more details.
    pub const fn const_param_ident_unit(
        &self,
        name: &'static str,
        accession: AccessionIntCode,
        unit: Unit,
    ) -> ParamCow<'static> {
        self.const_param(name, ValueRef::Empty, accession, unit)
    }

    /// Create a [`Param`] whose accession comes from this controlled vocabulary namespace with
    /// the given value.
    ///
    /// # Arguments
    /// - `accession`: The accession code for the [`Param`]. If specified as a [`CURIE`] or a string-like type,
    ///   any namespace is ignored.
    /// - `name`: The name of the parameter
    /// - `value`: The value of the parameter
    ///
    /// # See Also
    /// - [`ControlledVocabulary::param`]
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

/// An error describing a failure to map a controlled vocabulary identifier
/// to a known namespace
#[derive(Debug, Clone, Error)]
pub enum ControlledVocabularyResolutionError {
    #[error("Unrecognized controlled vocabulary {0}")]
    UnknownControlledVocabulary(String),
    #[error("Unrecognized controlled vocabulary code {0}")]
    UnknownControlledVocabularyCode(u8),
}

impl FromStr for ControlledVocabulary {
    type Err = ControlledVocabularyResolutionError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "MS" | "PSI-MS" => Ok(Self::MS),
            "UO" => Ok(Self::UO),
            EFO_CV => Ok(Self::EFO),
            OBI_CV => Ok(Self::OBI),
            BFO_CV => Ok(Self::BFO),
            HANCESTRO_CV => Ok(Self::HANCESTRO),
            #[cfg(feature = "imzml")]
            IMS_CV => Ok(Self::IMS),
            _ => Ok(Self::Unknown),
        }
    }
}

pub type ParamList = Vec<Param>;

pub trait ParamDescribedRead {
    /// Obtain an immutable slice over the encapsulated [`Param`] list
    fn params(&self) -> &[Param];

    /// Find the first [`Param`] whose name matches `name`
    fn get_param_by_name(&self, name: &str) -> Option<&Param> {
        self.params().iter().find(|&param| param.name == name)
    }

    /// Find the first [`Param`] whose [`CURIE`] matches `curie`
    fn get_param_by_curie(&self, curie: &CURIE) -> Option<&Param> {
        self.params().iter().find(|&param| curie == param)
    }

    /// Find the first [`Param`] whose [`Param::accession`] matches `accession`
    ///
    /// This is equivalent to [`ParamDescribed::get_param_by_curie`] on `accession.parse::<CURIE>().unwrap()`
    fn get_param_by_accession(&self, accession: &str) -> Option<&Param> {
        let (cv, acc_num) = curie_to_num(accession);
        self.params()
            .iter()
            .find(|&param| param.accession == acc_num && param.controlled_vocabulary == cv)
    }

    /// Iterate over the encapsulated parameter list
    fn iter_params(&self) -> std::slice::Iter<'_, Param> {
        self.params().iter()
    }
}

pub trait ParamDescribedMut {
    /// Obtain an mutable slice over the encapsulated [`Param`] list
    fn params_mut(&mut self) -> &mut ParamList;

    /// Add a new [`Param`] to the entity
    fn add_param(&mut self, param: Param) {
        self.params_mut().push(param);
    }

    /// Add all parameters from an iterator of [`Param`] to the entity
    fn extend_params(&mut self, it: impl IntoIterator<Item = Param>) {
        self.params_mut().extend(it)
    }

    /// Remove the `i`th [`Param`] from the entity.
    fn remove_param(&mut self, index: usize) -> Param {
        self.params_mut().remove(index)
    }

    /// Iterate mutably over the encapsulated parameter list
    fn iter_params_mut(&mut self) -> std::slice::IterMut<'_, Param> {
        self.params_mut().iter_mut()
    }
}

impl ParamDescribedRead for &[Param] {
    fn params(&self) -> &[Param] {
        self
    }
}

impl ParamDescribedRead for Vec<Param> {
    fn params(&self) -> &[Param] {
        self.as_ref()
    }
}

pub trait ParamDescribed {
    /// Obtain an immutable slice over the encapsulated [`Param`] list
    fn params(&self) -> &[Param];

    /// Obtain an mutable slice over the encapsulated [`Param`] list
    fn params_mut(&mut self) -> &mut ParamList;

    /// Add a new [`Param`] to the entity
    fn add_param(&mut self, param: Param) {
        self.params_mut().push(param);
    }

    /// Add all parameters from an iterator of [`Param`] to the entity
    fn extend_params(&mut self, it: impl IntoIterator<Item = Param>) {
        self.params_mut().extend(it)
    }

    /// Remove the `i`th [`Param`] from the entity.
    fn remove_param(&mut self, index: usize) -> Param {
        self.params_mut().remove(index)
    }

    /// Find the first [`Param`] whose name matches `name`
    fn get_param_by_name(&self, name: &str) -> Option<&Param> {
        self.params().iter().find(|&param| param.name == name)
    }

    /// Find the first [`Param`] whose [`CURIE`] matches `curie`
    fn get_param_by_curie(&self, curie: &CURIE) -> Option<&Param> {
        self.params().iter().find(|&param| curie == param)
    }

    /// Find the first [`Param`] whose [`Param::accession`] matches `accession`
    ///
    /// This is equivalent to [`ParamDescribed::get_param_by_curie`] on `accession.parse::<CURIE>().unwrap()`
    fn get_param_by_accession(&self, accession: &str) -> Option<&Param> {
        let (cv, acc_num) = curie_to_num(accession);
        self.params()
            .iter()
            .find(|&param| param.accession == acc_num && param.controlled_vocabulary == cv)
    }

    /// Iterate over the encapsulated parameter list
    fn iter_params(&self) -> std::slice::Iter<'_, Param> {
        self.params().iter()
    }

    /// Iterate mutably over the encapsulated parameter list
    fn iter_params_mut(&mut self) -> std::slice::IterMut<'_, Param> {
        self.params_mut().iter_mut()
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

/// Easily define all units. See the usage for more details. The accession and name both have to be
/// included in `&'static str` and `byte slice` to allow for constant time matching.
macro_rules! units {
    [$($unit:ident, $accession:literal, $baccession:literal, $name:literal, $bname:literal, $cv:ident, $id:literal);*;] => {
        /// Units that a term's value might have
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
        #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
        pub enum Unit {
            Unknown,
            $($unit,)*
        }

        impl Unit {
            pub const fn for_param(&self) -> (&'static str, &'static str) {
                match self {
                    $(Self::$unit => ($accession, $name),)*
                    Self::Unknown => ("","")
                }
            }

            pub const fn from_name(name: &str) -> Unit {
                match name.as_bytes() {
                    $($bname => Self::$unit,)*
                    _ => Self::Unknown
                }
            }

            pub const fn from_accession(acc: &str) -> Unit {
                match acc.as_bytes() {
                    $($baccession => Self::$unit,)*
                    _ => Self::Unknown
                }
            }

            pub const fn from_curie(acc: &CURIE) -> Unit {
                match acc {
                    $(CURIE {
                        controlled_vocabulary: ControlledVocabulary::$cv,
                        accession: $id,
                    } => Self::$unit,)*
                    _ => Self::Unknown
                }
            }

            pub const fn to_curie(&self) -> Option<CURIE> {
                match self {
                    $(Self::$unit => Some(CURIE {
                        controlled_vocabulary: ControlledVocabulary::$cv,
                        accession: $id,
                    }),)*
                    Self::Unknown => None
                }
            }

            pub const fn from_param(param: &Param) -> Unit {
                param.unit
            }

            pub const fn is_unknown(&self) -> bool {
                matches!(self, Self::Unknown)
            }
        }
    };
}

units![
    AbsorbanceUnit, "UO:0000269", b"UO:0000269", "absorbance unit", b"absorbance unit", UO, 269;
    Celsius, "UO:0000027", b"UO:0000027", "degree Celsius", b"degree Celsius", UO, 27;
    CountsPerSecond, "MS:1000814", b"MS:1000814", "counts per second", b"counts per second", MS, 1000814;
    DetectorCounts, "MS:1000131", b"MS:1000131", "number of detector counts", b"number of detector counts", MS, 1000131;
    Dimensionless, "UO:0000186", b"UO:0000186", "dimensionless unit", b"dimensionless unit", UO, 186;
    Electronvolt, "UO:0000266", b"UO:0000266", "electronvolt", b"electronvolt", UO, 266;
    Kelvin, "UO:0000012", b"UO:0000012", "kelvin", b"kelvin", UO, 12;
    Mass, "UO:000221", b"UO:000221", "dalton", b"dalton", UO, 221;
    MicrolitersPerMinute, "UO:0000271", b"UO:0000271", "microliters per minute", b"microliters per minute", UO, 271;
    Millisecond, "UO:0000028", b"UO:0000028", "millisecond", b"millisecond", UO, 28;
    Minute, "UO:0000031", b"UO:0000031", "minute", b"minute", UO, 31;
    MZ, "MS:1000040", b"MS:1000040", "m/z", b"m/z", MS, 1000040;
    Nanometer, "UO:0000018", b"UO:0000018", "nanometer", b"nanometer", UO, 18;
    Micrometer, "UO:0000017", b"UO:0000017", "micrometer", b"micrometer", UO, 17;
    Millimeter, "UO:0000016", b"UO:0000016", "millimeter", b"millimeter", UO, 16;
    Centimeter, "UO:0000015", b"UO:0000015", "centimeter", b"centimeter", UO, 15;
    PartsPerMillion, "UO:0000169", b"UO:0000169", "parts per million ", b"parts per million ", UO, 169;
    Pascal, "UO:0000110", b"UO:0000110", "pascal", b"pascal", UO, 110;
    Percent, "UO:0000187", b"UO:0000187", "percent", b"percent", UO, 187;
    PercentBasePeak, "MS:1000132", b"MS:1000132", "percent of base peak", b"percent of base peak", MS, 1000132;
    PercentBasePeakTimes100, "MS:1000905", b"MS:1000905", "percent of base peak times 100", b"percent of base peak times 100", MS, 1000905;
    Psi, "UO:0010052", b"UO:0010052", "pounds per square inch", b"pounds per square inch", UO, 10052;
    Second, "UO:0000010", b"UO:0000010", "second", b"second", UO, 10;
    Volt, "UO:0000218", b"UO:0000218", "volt", b"volt", UO, 218;
    VoltSecondPerSquareCentimeter, "MS:1002814", b"MS:1002814", "volt-second per square centimeter", b"volt-second per square centimeter", MS, 1002814;
    Hertz, "UO:000106", b"UO:000106", "hertz", b"hertz", UO, 106;
    Liter, "UO:0000099", b"UO:0000099", "liter", b"liter", UO, 99;
    Milliliter, "UO:0000098", b"UO:0000098", "milliliter", b"milliliter", UO, 98;
    Microliter, "UO:0000101", b"UO:0000101", "microliter", b"microliter", UO, 101;
];

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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_build_param() {
        assert_eq!(
            ParamBuilder::default()
                .name("dalton")
                .curie(curie!(UO:221))
                .build(),
            ControlledVocabulary::UO.param("UO:000221", "dalton")
        );
        // <cvParam cvRef="MS" accession="MS:1000529" name="instrument serial number" value="FSN10375"/>
        let p = ParamBuilder::default()
            .controlled_vocabulary(ControlledVocabulary::MS)
            .accession(1000529)
            .name("instrument serial number")
            .value("FSN10375")
            .unit(Unit::Unknown)
            .build();
        assert_eq!(p.value(), "FSN10375");
        assert_eq!(p.unit(), Unit::Unknown);
    }

    #[test]
    fn test_value() {
        let x = 42;
        let mut val: Value = x.into();
        let mut val_ref: ValueRef = x.into();
        let mut val_ref2: ValueRef = (&x).into();
        assert_eq!(val, x);
        assert_eq!(val_ref, x);
        assert_eq!(val_ref, val);
        assert_eq!(val_ref2, val_ref);
        assert!(val.to_bool().unwrap());
        assert!(val_ref.to_bool().unwrap());
        assert_eq!(val.to_str(), x.to_string());
        assert_eq!(val_ref.to_str(), x.to_string());
        val = x.to_string().parse().unwrap();
        val_ref = x.to_string().parse().unwrap();
        assert_eq!(val, x);
        assert_eq!(val_ref, x);
        assert_eq!(val_ref, val);

        let x2 = Some(x);
        val = x2.into();
        val_ref = x2.into();
        val_ref2 = (&x).into();
        assert_eq!(val, x);
        assert_eq!(val_ref, x);
        assert_eq!(val_ref, val);
        assert_eq!(val_ref2, val_ref);
        assert!(val.to_bool().unwrap());
        assert!(val_ref.to_bool().unwrap());
        assert_eq!(val.to_str(), x.to_string());
        assert_eq!(val_ref.to_str(), x.to_string());
        val = x.to_string().parse().unwrap();
        val_ref = x.to_string().parse().unwrap();
        assert_eq!(val, x);
        assert_eq!(val_ref, x);
        assert_eq!(val_ref, val);

        let x = 42.01;
        val = x.into();
        val_ref = x.into();
        val_ref2 = (&x).into();
        assert_eq!(val, x);
        assert_eq!(val_ref, x);
        assert_eq!(val_ref, val);
        assert_eq!(val_ref2, val_ref);
        assert!(val.to_bool().unwrap());
        assert!(val_ref.to_bool().unwrap());
        assert_eq!(val.to_str(), x.to_string());
        assert_eq!(val_ref.to_str(), x.to_string());
        val = x.to_string().parse().unwrap();
        val_ref = x.to_string().parse().unwrap();
        assert_eq!(val, x);
        assert_eq!(val_ref, x);
        assert_eq!(val_ref, val);

        let x2 = Some(x);
        val = x2.into();
        val_ref = x2.into();
        assert_eq!(val, x);
        assert_eq!(val_ref, x);
        assert_eq!(val_ref, val);
        assert!(val.to_bool().unwrap());
        assert!(val_ref.to_bool().unwrap());
        assert_eq!(val.to_str(), x.to_string());
        assert_eq!(val_ref.to_str(), x.to_string());
        val = x.to_string().parse().unwrap();
        val_ref = x.to_string().parse().unwrap();
        assert_eq!(val, x);
        assert_eq!(val_ref, x);
        assert_eq!(val_ref, val);

        let x = true;
        val = x.into();
        val_ref = x.into();
        assert_eq!(val, x);
        assert_eq!(val_ref, x);
        assert_eq!(val_ref, val);
        assert!(val.to_bool().unwrap());
        assert!(val_ref.to_bool().unwrap());
        assert_eq!(val.to_str(), x.to_string());
        assert_eq!(val_ref.to_str(), x.to_string());
        val = x.to_string().parse().unwrap();
        val_ref = x.to_string().parse().unwrap();
        assert_eq!(val, x);
        assert_eq!(val_ref, x);
        assert_eq!(val_ref, val);

        let x = "Foobar".to_string();
        val = x.clone().into();
        val_ref = x.clone().into();
        val_ref2 = x.as_str().into();
        assert_eq!(val, x);
        assert_eq!(val_ref, x);
        assert_eq!(val_ref, val);
        assert_eq!(val_ref2, val_ref);
        assert_eq!(val.to_str(), x.to_string());
        assert_eq!(val_ref.to_str(), x.to_string());
        val = x.to_string().parse().unwrap();
        val_ref = x.to_string().parse().unwrap();
        assert_eq!(val, x);
        assert_eq!(val_ref, x);
        assert_eq!(val_ref, val);
        assert_eq!(val.to_buffer().unwrap(), x.as_bytes());
        assert_eq!(val_ref.to_buffer().unwrap(), x.as_bytes());
        assert_eq!(val_ref.to_buffer().unwrap(), val.to_buffer().unwrap());
    }
}
