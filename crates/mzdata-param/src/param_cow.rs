use std::borrow::Cow;
use std::str::{self, FromStr};


use crate::{
    AccessionIntCode, CURIE, ControlledVocabulary, Param, ParamLike, ParamValue, ParamValueParseError, Unit, Value, ValueRef
};


/// A statically allocate-able or non-owned data version of [`Param`]
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct ParamCow<'a> {
    /// The human-readable name of the parameter.
    pub name: Cow<'a, str>,
    /// The parameter's value.
    pub value: ValueRef<'a>,
    /// The numeric accession code within `controlled_vocabulary`, if this parameter is
    /// drawn from a controlled vocabulary.
    pub accession: Option<AccessionIntCode>,
    /// The controlled vocabulary this parameter's term belongs to. `None` for a
    /// free-text, user-defined parameter.
    pub controlled_vocabulary: Option<ControlledVocabulary>,
    /// The unit the value is expressed in, if known.
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

    fn is_list(&self) -> bool {
        <ValueRef<'a> as ParamValue>::is_list(&self.value)
    }

    fn as_slice(&self) -> Cow<'_, [Value]> {
        <ValueRef<'a> as ParamValue>::as_slice(&self.value)
    }
}

impl ParamCow<'static> {
    /// Build a [`ParamCow<'static>`] in a `const` context. Prefer
    /// [`ControlledVocabulary::const_param`] and its siblings, which fill in the namespace
    /// for you.
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
    /// Build a [`ParamCow`] from its parts directly.
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

    /// Create a [`CURIE`] from [`ParamCow::controlled_vocabulary`] and [`ParamCow::accession`]
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
