//! Elements of controlled vocabularies used to describe mass spectra and their components.
//!
//! Directly maps to the usage of the PSI-MS controlled vocabulary in mzML
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

/// A syntactic shortcut for creating [`CURIE`] instances using compact notation
#[macro_export]
macro_rules! curie {
    ($ns:ident:$acc:literal) => {
        $crate::CURIE {
            controlled_vocabulary: $crate::ControlledVocabulary::$ns,
            accession: $acc,
        }
    };
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

pub type ParamList = Vec<Param>;

/// A read-only form of [`ParamDescribed`]
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
