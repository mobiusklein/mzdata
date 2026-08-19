use std::borrow::Cow;
use std::fmt::Display;
use std::hash::Hash;
use std::str::{self, FromStr};
use std::mem;


use crate::{ParamValue, ParamValueParseError, Value};

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
    /// A collection of heterogenous [`Value`]
    List(Cow<'a, [Value]>),
}

impl Eq for ValueRef<'_> {}

impl<U: Into<Value>> FromIterator<U> for ValueRef<'static> {
    fn from_iter<T: IntoIterator<Item = U>>(iter: T) -> Self {
        let values: Vec<Value> = iter.into_iter().map(|v| v.into()).collect();
        Self::List(Cow::Owned(values))
    }
}


impl From<Vec<Value>> for ValueRef<'static> {
    fn from(value: Vec<Value>) -> Self {
        Self::List(Cow::Owned(value))
    }
}

impl<'a> From<Cow<'a, [Value]>> for ValueRef<'a> {
    fn from(value: Cow<'a, [Value]>) -> Self {
        Self::List(value)
    }
}

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
            Self::List(v) => {
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

    fn is_list(&self) -> bool {
        matches!(self, Self::List(_))
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

    /// Store the value as a boolean
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

    /// Store the value as a list
    pub fn coerce_list(&mut self) -> Result<(), ParamValueParseError> {
        if !self.is_list() {
            let dup = match self {
                Self::Boolean(v) => Value::Boolean(*v),
                Self::Empty => Value::Empty,
                Self::Float(v) => Value::Float(*v),
                Self::Int(v) => Value::Int(*v),
                Self::Buffer(v) => Value::Buffer(v.to_vec().into()),
                Self::String(v) => Value::String(v.to_string()),
                Self::List(_) => unimplemented!(),
            };
            *self = Self::List(Cow::Owned([dup].into()));
        }
        Ok(())
    }

    fn parse<T: FromStr>(&self) -> Result<T, T::Err> {
        match self {
            Self::String(s) => s.parse(),
            Self::Float(v) => v.to_string().parse(),
            Self::Int(i) => i.to_string().parse(),
            Self::Buffer(b) => String::from_utf8_lossy(b).parse(),
            Self::Empty => "".parse(),
            Self::Boolean(v) => v.to_string().parse(),
            Self::List(_) => self.to_string().parse(),
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
            Self::List(_) => Cow::Owned(self.to_string().into_bytes()),
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
            Self::List(v) => v.iter().map(|vi| vi.data_len()).sum(),
        }
    }

    fn is_list(&self) -> bool {
        self.is_list()
    }

    fn as_slice(&self) -> Cow<'_, [Value]> {
        match self {
            Self::List(v) => Cow::Borrowed(v),
            _ => {
                let dup = match self {
                    Self::Boolean(v) => Value::Boolean(*v),
                    Self::Empty => Value::Empty,
                    Self::Float(v) => Value::Float(*v),
                    Self::Int(v) => Value::Int(*v),
                    Self::Buffer(v) => Value::Buffer(v.to_vec().into()),
                    Self::String(v) => Value::String(v.to_string()),
                    Self::List(_) => unimplemented!(),
                };
                Cow::Owned([dup].into())
            }
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
            Value::List(v) => Self::List(Cow::Borrowed(v)),
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
            ValueRef::List(v) => {
                let mut ve = Vec::with_capacity(v.len());
                for vi in v.iter() {
                    ve.push(vi.clone())
                }
                Self::List(ve.into_boxed_slice())
            }
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
            Self::List(v) => v.iter().for_each(|vi| vi.hash(state)),
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
