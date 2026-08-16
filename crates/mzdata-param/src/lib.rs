//! Elements of controlled vocabularies used to describe mass spectra and their components.
//!
//! This crate implements the "CV param" model used throughout the PSI-MS family of mass
//! spectrometry file formats (mzML, imzML, mzMLb, ...): every piece of metadata is either
//! a term drawn from a controlled vocabulary (like the PSI-MS ontology) or a free-text,
//! user-defined key/value pair. Both are represented uniformly by [`Param`].
//!
//! ## Core types
//!
//! - [`Param`] / [`ParamCow`]: a single CV or user parameter. [`Param`] owns its data;
//!   [`ParamCow`] borrows where possible (useful for `const` tables of well-known terms).
//! - [`Value`] / [`ValueRef`]: the value half of a [`Param`], which may be a string,
//!   integer, float, boolean, byte buffer, list, or empty. [`Value`] owns its data;
//!   [`ValueRef`] borrows where possible. Both implement [`ParamValue`] for reading and
//!   coercing the stored value.
//! - [`CURIE`]: a `namespace:accession` identifier (e.g. `MS:1000041`) that names a term
//!   within a [`ControlledVocabulary`].
//! - [`ControlledVocabulary`]: the vocabulary a term belongs to (`MS`, `UO`, ...), and a
//!   factory for building [`Param`]s within that namespace.
//! - [`Unit`]: a closed set of units of measure (from the Unit Ontology and PSI-MS) that
//!   a [`Param`]'s value may be expressed in.
//! - [`ParamDescribed`] / [`ParamDescribedRead`]: traits for types that carry a list of
//!   [`Param`]s (spectra, scans, instrument components, ...), with helpers to look a
//!   parameter up by name, [`CURIE`], or accession string.
//!
//! ## Building a parameter
//!
//! Construct a controlled-vocabulary parameter from a [`ControlledVocabulary`] namespace:
//!
//! ```rust
//! use mzdata_param::{ControlledVocabulary, ParamLike, ParamValue};
//!
//! let param = ControlledVocabulary::MS.param_val("MS:1000041", "charge state", 2i32);
//! assert_eq!(param.name(), "charge state");
//! assert_eq!(param.curie().unwrap().to_string(), "MS:1000041");
//! assert_eq!(param.value().to_i64().unwrap(), 2);
//! ```
//!
//! Or use [`ParamBuilder`] / [`Param::builder`] for more incremental construction,
//! including plain user-defined parameters that have no controlled vocabulary:
//!
//! ```rust
//! use mzdata_param::{Param, ParamValue, Unit, curie};
//!
//! let p = Param::builder()
//!     .name("scan start time")
//!     .curie(curie!(MS:1000016))
//!     .value(12.34)
//!     .unit(Unit::Minute)
//!     .build();
//! assert_eq!(p.to_f64().unwrap(), 12.34);
//! assert_eq!(p.unit, Unit::Minute);
//!
//! // A user-defined parameter has no CURIE at all.
//! let custom = Param::new_key_value("my custom field", "some value");
//! assert!(!custom.is_controlled());
//! ```
//!
//! ## Reading values generically
//!
//! [`ParamValue`] lets you read a [`Param`]'s value without caring whether it was parsed
//! as a string, integer, or float - coercions are attempted on demand:
//!
//! ```rust
//! use mzdata_param::{ControlledVocabulary, ParamValue};
//!
//! let p = ControlledVocabulary::MS.param_val("MS:1000827", "isolation window target m/z", "500.25");
//! assert_eq!(p.to_f64().unwrap(), 500.25);
//! assert_eq!(p.to_str(), "500.25");
//! ```
//!
//! ## Attaching parameters to types with [`ParamDescribed`]
//!
//! The [`ParamDescribed`] trait provides many methods that make it easier to operate
//! on instances of types which are described by a list of [`Param`]s.
//!
//! ### Simple implementation helpers
//!
//! Use [`impl_param_described!`] (for a plain `Vec<Param>` field) or
//! [`impl_param_described_deferred!`] (for an `Option<Vec<Param>>` field that is
//! lazily allocated on first write) to implement [`ParamDescribed`] for your type:
//!
//! ```rust
//! use mzdata_param::{impl_param_described, Param, ParamDescribed, ParamList, ParamValue};
//!
//! #[derive(Default)]
//! struct MyComponent {
//!     params: ParamList,
//! }
//!
//! impl_param_described!(MyComponent);
//!
//! let mut c = MyComponent::default();
//! c.add_param(Param::new_key_value("vendor", "Acme"));
//! assert_eq!(c.get_param_by_name("vendor").unwrap().to_str(), "Acme");
//! ```
//!
//! ## Matching against known terms with `CURIE`
//!
//! [`CURIE`] implements `PartialEq` against anything implementing [`ParamLike`], so you can
//! compare a parameter directly against a well-known accession, and the [`curie!`] macro
//! gives a compact way to write one inline:
//!
//! ```rust
//! use mzdata_param::{ControlledVocabulary, ParamLike, curie};
//!
//! let p = ControlledVocabulary::MS.param("MS:1000016", "scan start time");
//! assert!(curie!(MS:1000016) == p);
//! ```
//!
//! ## The PSI-MS ontology (`cv` feature)
//!
//! With the `cv` feature enabled, [`MSVocabulary`] gives access to the full PSI-MS
//! ontology (term names, synonyms, and the `is_a` parent/child hierarchy), backed by an
//! embedded static snapshot and an optional on-disk cache that can be refreshed from a
//! `.obo` file via [`MSVocabulary::update_from_obo`]:
//!
//! ```rust
//! # #[cfg(feature = "cv")]
//! # {
//! use mzdata_param::{MSVocabulary, curie};
//!
//! let term = MSVocabulary::get(curie!(MS:1000044)).unwrap();
//! assert_eq!(term.name.as_ref(), "dissociation method");
//! assert!(MSVocabulary::is_child_of(curie!(MS:1000133), curie!(MS:1000044)));
//! # }
//! ```
use std::borrow::Cow;
use std::convert::TryFrom;
use std::fmt::Display;
use std::hash::Hash;
use std::num;
use std::str::{self, FromStr};

use thiserror::Error;

pub(crate) mod value;

pub use value::{ParamValue, ParamValueParseError, Value};

pub(crate) mod value_ref;

pub use value_ref::ValueRef;

pub(crate) mod curie_;

pub use curie_::{
    AccessionCodeParseError, AccessionIntCode, CURIE, CURIEParsingError, ControlledVocabulary,
    ControlledVocabularyResolutionError, curie_to_num,
};

#[cfg(feature = "cv")]
pub(crate) mod cv;

#[cfg(feature = "cv")]
pub use cv::{CVTraversal, MSTerm, MSVocabulary, VocabularyData};

/// A helper to generate methods that find a value by a [`CURIE`]
#[macro_export]
macro_rules! find_param_method {
    ($meth:ident, $curie:expr) => {
        $crate::find_param_method!($meth, $curie, "Find a parameter by its CURIE");
    };
    ($meth:ident, $curie:expr, $desc:literal) => {
        #[doc=$desc]
        pub fn $meth(&self) -> Option<$crate::ValueRef<'_>> {
            self.get_param_by_curie($curie)
                .map(|p| $crate::ParamLike::value(p))
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

/// A syntactic shortcut for creating [`CURIE`] instances using compact notation.
///
/// The following are identical.
/// ```rust
/// # use mzdata_param::{CURIE, ControlledVocabulary, curie};
/// assert_eq!(CURIE::new(ControlledVocabulary::MS, 1000016), curie!(MS:1000016));
/// ```
///
/// This macro has two advantages for writing [`CURIE`] "constants", it does not
/// require a function call (so may be used in `match` without extra overhead) and
/// is far more succinct.
#[macro_export]
macro_rules! curie {
    ($ns:ident:$acc:literal) => {
        $crate::CURIE {
            controlled_vocabulary: $crate::ControlledVocabulary::$ns,
            accession: $acc,
        }
    };
}

/// A minimal, read-only view over a single CV or user-defined parameter.
///
/// Both [`Param`] and [`ParamCow`] implement this trait; prefer it when writing code that
/// should work with either the owned or borrowed representation.
pub trait ParamLike {
    /// The human-readable name of the parameter.
    fn name(&self) -> &str;
    /// A borrowed view of the parameter's value.
    fn value(&self) -> ValueRef<'_>;
    /// The numeric accession code within [`ParamLike::controlled_vocabulary`], if this
    /// parameter is drawn from a controlled vocabulary.
    fn accession(&self) -> Option<AccessionIntCode>;
    /// The controlled vocabulary this parameter's term belongs to, if any.
    fn controlled_vocabulary(&self) -> Option<ControlledVocabulary>;
    /// The unit the parameter's value is expressed in, if known.
    fn unit(&self) -> Unit;

    /// Check whether this parameter's term belongs to the PSI-MS controlled vocabulary.
    fn is_ms(&self) -> bool {
        if let Some(cv) = self.controlled_vocabulary() {
            cv == ControlledVocabulary::MS
        } else {
            false
        }
    }

    /// Parse the parameter's value as `T`. See [`ParamValue::parse`].
    fn parse<T: str::FromStr>(&self) -> Result<T, T::Err> {
        self.value().parse::<T>()
    }

    /// Check whether this parameter is drawn from a controlled vocabulary, as opposed to
    /// being a free-text user-defined parameter.
    fn is_controlled(&self) -> bool {
        self.accession().is_some()
    }

    /// Build a [`CURIE`] from [`ParamLike::controlled_vocabulary`] and
    /// [`ParamLike::accession`], if both are present.
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

pub(crate) mod param_cow;
pub use param_cow::ParamCow;

pub(crate) mod param;

pub use param::{Param, ParamBuilder};

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

/// The concrete container type used to hold a [`Param`] collection, as referenced by
/// [`ParamDescribed`] and the [`impl_param_described!`]/[`impl_param_described_deferred!`]
/// macros.
pub type ParamList = Vec<Param>;

/// A read-only form of [`ParamDescribed`], implemented directly for `&[Param]` so that a
/// borrowed slice of parameters can be queried the same way as a type that owns its
/// parameter list.
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

impl ParamDescribedRead for &[Param] {
    fn params(&self) -> &[Param] {
        self
    }
}


/// A type that has a [`ParamList`] that uses [`Param`] instances to describe an entity
/// with key-value pairs.
///
/// Most implementors get this via [`impl_param_described!`] or
/// [`impl_param_described_deferred!`] rather than implementing it by hand.
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

        impl $crate::ParamDescribed for $t {
            fn params(&self) -> &[$crate::Param] {
                return &self.params
            }

            fn params_mut(&mut self) -> &mut $crate::ParamList {
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
        impl $crate::ParamDescribed for $t {
            fn params(&self) -> &[$crate::Param] {
                match &self.params {
                    Some(val) => &val,
                    None => {
                        $crate::_EMPTY_PARAM
                    }
                }
            }

            fn params_mut(&mut self) -> &mut $crate::ParamList {
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

pub(crate) mod units;

pub use units::Unit;

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

    #[cfg(feature = "static_data")]
    #[test]
    fn test_mzcv_ms() {
        let cv = MSVocabulary::init();
        let term = cv.get_by_name("ms level").unwrap();
        assert_eq!(term.name.as_ref(), "ms level");
        let parents_of = CVTraversal::parents_of(cv, term.curie()).unwrap();
        let parent = cv
            .get_by_index(&parents_of.iter().next().unwrap().0)
            .unwrap();
        assert_eq!(parent.name.as_ref(), "spectrum attribute");
    }
}
