use std::borrow::Cow;
use std::fmt::Display;
use std::hash::Hash;
use std::str::{self, FromStr};


use crate::{
    AccessionIntCode, CURIE, ControlledVocabulary, ParamCow, ParamLike, ParamValue, ParamValueParseError, Unit, Value, ValueRef
};


/// A controlled vocabulary or user parameter
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Param {
    /// The human-readable name of the parameter.
    pub name: String,
    /// A type generic value associated with the parameter
    pub value: Value,
    /// The numeric component of the accession code in the identifier for the parameter's
    /// definition in [`Self::controlled_vocabulary`], if it is defined by one.
    pub accession: Option<AccessionIntCode>,
    /// The controlled vocabulary this parameter is defined in, if any.
    pub controlled_vocabulary: Option<ControlledVocabulary>,
    /// The unit the parameter's value is expressed in, if any.
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

    fn is_list(&self) -> bool {
        <Value as ParamValue>::is_list(&self.value)
    }

    fn as_slice(&self) -> Cow<'_, [Value]> {
        <Value as ParamValue>::as_slice(&self.value)
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
    /// Set the parameter's name.
    pub fn name<S: ToString>(mut self, name: S) -> Self {
        self.name = name.to_string();
        self
    }

    /// Set the parameter's value.
    pub fn value<V: Into<Value>>(mut self, value: V) -> Self {
        self.value = value.into();
        self
    }

    /// Set the parameter's controlled vocabulary namespace.
    ///
    /// # See also
    /// - [`ParamBuilder::curie`], to set the namespace and accession together
    pub fn controlled_vocabulary(mut self, cv: ControlledVocabulary) -> Self {
        self.controlled_vocabulary = Some(cv);
        self
    }

    /// Set the parameter's accession code.
    ///
    /// # See also
    /// - [`ParamBuilder::curie`], to set the namespace and accession together
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

    /// Set the parameter's unit of measure.
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

    /// Create a new [`ParamBuilder`] to make creating a new [`Param`] more convenient than
    /// setting fields on a bare [`Param`] directly.
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
    ///
    /// # See also
    /// - [`Param::with_unit_t`], to infer the unit from a resolved [`Unit`]
    pub fn with_unit<S: AsRef<str>, A: AsRef<str>>(mut self, accession: S, name: A) -> Param {
        self.unit = Unit::from_accession(accession.as_ref());
        if matches!(self.unit, Unit::Unknown) {
            self.unit = Unit::from_name(name.as_ref());
        }
        self
    }

    /// Set [`Param::unit`] directly from `unit`.
    ///
    /// # See also
    /// - [`Param::with_unit`], to infer the unit from an accession or name instead
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
