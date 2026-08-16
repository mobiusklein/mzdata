use std::borrow::Cow;
use std::fmt::Display;
use std::hash::Hash;
use std::str::{self, FromStr};
use std::{io, mem};

use thiserror::Error;

use crate::ValueRef;


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
    /// A list of heterogenous [`Value`]
    List(Box<[Value]>),
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

    /// Check if the value is a list
    fn is_list(&self) -> bool;

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

    /// Get the value as a slice
    fn as_slice(&self) -> Cow<'_, [Value]>;

    /// Convert the value's string representation to `T` if possible
    fn parse<T: FromStr>(&self) -> Result<T, T::Err>;

    /// Convert the value to a byte string, the bytes
    /// of the string representation.
    fn as_bytes(&self) -> Cow<'_, [u8]>;

    /// Get a reference to the stored value
    fn as_ref(&self) -> crate::ValueRef<'_>;

    /// Get the size of the stored data type
    fn data_len(&self) -> usize;
}

/// Errors that might occur while trying to convert to a particular value type
/// from whatever type might be stored in a [`ParamValue`]-like object.
#[derive(Debug, Clone, Error, PartialEq)]
pub enum ParamValueParseError {
    /// The value could not be interpreted as a floating point number.
    #[error("Failed to extract a float from {0:?}")]
    FailedToExtractFloat(Option<String>),
    /// The value could not be interpreted as an integer.
    #[error("Failed to extract a int from {0:?}")]
    FailedToExtractInt(Option<String>),
    /// The value could not be interpreted as a string.
    #[error("Failed to extract a string")]
    FailedToExtractString,
    /// The value could not be interpreted as a byte buffer.
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
            Value::List(v) => {
                f.write_str("[ ")?;
                if let Some(vi) = v.first() {
                    vi.fmt(f)?;
                }
                for vi in v.iter().skip(1) {
                    f.write_str(", ")?;
                    vi.fmt(f)?;
                }
                f.write_str(" ]")
            }
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

    /// See [`ParamValue::is_empty`].
    pub fn is_empty(&self) -> bool {
        matches!(self, Self::Empty)
    }

    /// See [`ParamValue::is_i64`].
    pub fn is_i64(&self) -> bool {
        matches!(self, Self::Int(_))
    }

    /// See [`ParamValue::is_f64`].
    pub fn is_f64(&self) -> bool {
        matches!(self, Self::Float(_))
    }

    /// See [`ParamValue::is_buffer`].
    pub fn is_buffer(&self) -> bool {
        matches!(self, Self::Buffer(_))
    }

    /// See [`ParamValue::is_str`].
    pub fn is_str(&self) -> bool {
        matches!(self, Self::String(_))
    }

    /// See [`ParamValue::is_list`].
    pub fn is_list(&self) -> bool {
        matches!(self, Self::List(_))
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

    /// Store the value as a boolean
    pub fn coerce_bool(&mut self) -> Result<(), ParamValueParseError> {
        let value = self.to_bool()?;
        *self = Self::Boolean(value);
        Ok(())
    }

    /// Store the value as a list
    pub fn coerce_list(&mut self) -> Result<(), ParamValueParseError> {
        if !self.is_list() {
            let mut tmp = Self::Empty;
            core::mem::swap(&mut tmp, self);
            *self = Self::List([tmp].into());
        }
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
            Value::List(_) => self.to_string().parse(),
        }
    }

    /// See [`ParamValue::to_bool`].
    pub fn to_bool(&self) -> Result<bool, ParamValueParseError> {
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

    /// See [`ParamValue::to_f64`].
    pub fn to_f64(&self) -> Result<f64, ParamValueParseError> {
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

    /// See [`ParamValue::to_i64`].
    pub fn to_i64(&self) -> Result<i64, ParamValueParseError> {
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

    /// See [`ParamValue::to_str`].
    pub fn to_str(&self) -> Cow<'_, str> {
        if let Self::String(val) = self {
            Cow::Borrowed(val)
        } else {
            Cow::Owned(self.to_string())
        }
    }

    /// See [`ParamValue::to_buffer`].
    pub fn to_buffer(&self) -> Result<Cow<'_, [u8]>, ParamValueParseError> {
        if let Self::Buffer(val) = self {
            Ok(Cow::Borrowed(val))
        } else if let Self::String(val) = self {
            Ok(Cow::Borrowed(val.as_bytes()))
        } else {
            Err(ParamValueParseError::FailedToExtractBuffer)
        }
    }

    /// See [`ParamValue::as_ref`].
    pub fn as_ref(&self) -> ValueRef<'_> {
        self.into()
    }

    /// View this [`Value`] as a [`slice`]
    pub fn as_slice(&self) -> &[Self] {
        if let Self::List(val) = self {
            val.as_ref()
        } else {
            core::slice::from_ref(self)
        }
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
            Self::List(_) => Cow::Owned(self.to_string().into_bytes()),
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
            Self::List(v) => v.iter().map(|vi| vi.data_len()).sum(),
        }
    }

    fn is_boolean(&self) -> bool {
        matches!(self, Self::Boolean(_))
    }

    fn to_bool(&self) -> Result<bool, ParamValueParseError> {
        self.to_bool()
    }

    fn is_list(&self) -> bool {
        self.is_list()
    }

    fn as_slice(&self) -> Cow<'_, [Value]> {
        Cow::Borrowed(self.as_slice())
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
            Self::List(v) => {
                v.iter().for_each(|vi| vi.hash(state));
            }
        }
    }
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


/// When the `serde` feature is enabled, [`mzdata_param::Value`] can be converted from
/// `serde_json::Value` without going through the generic deserialization process.
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
            Value::List(val) => {
                let mut ve = Vec::new();
                for vi in val {
                    ve.push(vi.into());
                }
                serde_json::Value::Array(ve)
            }
        }
    }
}
