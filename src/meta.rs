/*! Metadata describing mass spectrometry data files and their contents.
 */
#[macro_use]
mod file_description;
mod activation;
mod data_processing;
mod instrument;
mod run;
mod sample;
mod software;
#[macro_use]
mod traits;

use std::borrow::Cow;

pub use data_processing::{
    DataProcessing, DataProcessingAction, DataTransformationAction, FormatConversion,
    ProcessingMethod,
};
pub use software::{
    custom_software_name, Software, SoftwareTerm, SoftwareType as SoftwareTypeFlags
};

pub use file_description::{
    FileDescription, MassSpectrometerFileFormatTerm, NativeIDFormatError, NativeSpectrumIDFormat,
    NativeSpectrumIdentifierFormatTerm, SourceFile
};

pub use instrument::{
    Component, ComponentType, DetectorTypeTerm, InletTypeTerm, InstrumentConfiguration,
    IonizationTypeTerm, MassAnalyzerTerm
};

pub use activation::{DissociationEnergy, DissociationEnergyTerm, DissociationMethodTerm};
pub use run::MassSpectrometryRun;
pub use sample::Sample;
pub use traits::MSDataFileMetadata;

use crate::params::{ParamValueParseError, Value, ValueRef};

#[macro_export]
macro_rules! cvmap {
    (
        #[flag_type=$flag_type:ty]
        $(#[$enum_attrs:meta])*
        $vis:vis enum $enum_name:ident {
            $(
                #[term(cv=$cv:ident, accession=$accession:literal, name=$term_name:literal, flags=$flags:tt, parents=$parents:tt)]
                $(#[$variant_attrs:meta])*
                $variant:ident
            ),*
            $(,)?
        }
    ) => {
        $(#[$enum_attrs])*
        $vis enum $enum_name {
            $(
                $(#[$variant_attrs])*
                $variant,
            )*
        }

        /// These methods are part of the controlled vocabulary mapping
        impl $enum_name {

            /// Retrieve the accession number for this term, independent of its controlled vocabulary
            pub const fn accession(&self) -> $crate::params::AccessionIntCode {
                match self {
                    $(Self::$variant => $accession,)*
                }
            }

            /// Retrieve the controlled vocabulary this term belongs to
            pub const fn controlled_vocabulary(&self) -> $crate::params::ControlledVocabulary {
                match self {
                    $(Self::$variant => $crate::params::ControlledVocabulary::$cv,)*
                }
            }

            /// Retrieve the plain text human readable name for this term
            pub const fn name(&self) -> &'static str {
                match self {
                    $(Self::$variant => $term_name,)*
                }
            }

            /// Attempt to map a string by name to retrieve one of the terms from this
            /// set.
            ///
            /// If no match is found, [`None`] is returned.
            pub fn from_name(name: &str) -> Option<Self> {
                match name {
                    $($term_name => Some(Self::$variant),)*
                    _ => None
                }
            }


            /// Attempt to map the numeric accession number to retrieve one of the terms from this
            /// set.
            ///
            /// If no match is found, [`None`] is returned.
            pub const fn from_accession(accession: $crate::params::AccessionIntCode) -> Option<Self> {
                match accession {
                    $($accession => Some(Self::$variant),)*
                    _ => None
                }
            }

            /// Convert this term into a [`ParamCow`](crate::params::ParamCow) without a value.
            pub const fn to_param(self) -> $crate::params::ParamCow<'static> {
                $crate::params::ParamCow::const_new(
                    self.name(),
                    $crate::params::ValueRef::Empty,
                    Some(self.accession()),
                    Some(self.controlled_vocabulary()),
                    $crate::params::Unit::Unknown
                )
            }

            /// Convert a [`CURIE`]($crate::params::CURIE) by accession.
            ///
            /// If no match is found, [`None`] is returned.
            pub const fn from_curie(curie: &$crate::params::CURIE) -> Option<Self> {
                if matches!(curie.controlled_vocabulary, $crate::params::ControlledVocabulary::MS) {
                    Self::from_accession(curie.accession)
                } else {
                    None
                }
            }

            /// Attempt to convert a [`ParamCow`](crate::params::ParamCow) to a term from this set.
            ///
            /// If no match is found, [`None`] is returned.
            ///
            /// # Note
            /// This method can be called in `const` contexts, requiring the type be [`ParamCow`](crate::params::ParamCow) with a `'static`
            /// lifetime parameter, but the regular [`From`] trait is implemented for all [`ParamLike`](crate::params::ParamLike) types.
            pub const fn from_param(p: &$crate::params::ParamCow<'static>) -> Option<Self> {
                if let Some(acc) = p.accession {
                    Self::from_accession(acc)
                } else {
                    None
                }
            }

            /// Retrieve a term set specific set of flags
            pub fn flags(&self) -> $flag_type {
                match self {
                    $(Self::$variant => $flags.into(),)*
                }
            }

            /// Retrieve the list of zero or more terms in the set which are
            /// parents of this term.
            pub fn parents(&self) -> Vec<Self> {
                match self {
                    $(Self::$variant => $parents.iter().flat_map(|s: &&str| {
                        let curie = s.parse::<$crate::params::CURIE>().unwrap();
                        Self::from_accession(curie.accession)
                    }).collect(),)*
                }
            }

        }

        impl<P> From<P> for $enum_name where P: $crate::params::ParamLike {
            fn from(value: P) -> Self {
                Self::from_accession(
                    value.accession().expect(
                        concat!("Cannot convert an uncontrolled parameter to ", stringify!($enum_name)))
                ).unwrap_or_else(
                    || panic!(
                        "Could not map {:?}:{} to {}",
                        value.controlled_vocabulary().unwrap(),
                        value.accession().unwrap(),
                        stringify!($enum_name)
                    )
                )
            }
        }

        impl From<$enum_name> for $crate::params::ParamCow<'static> {
            fn from(value: $enum_name) -> Self {
                value.to_param()
            }
        }

        impl From<$enum_name> for $crate::params::Param {
            fn from(value: $enum_name) -> Self {
                value.to_param().into()
            }
        }

        impl From<&$enum_name> for $crate::params::ParamCow<'static> {
            fn from(value: &$enum_name) -> Self {
                value.to_param()
            }
        }

        impl From<&$enum_name> for $crate::params::Param {
            fn from(value: &$enum_name) -> Self {
                value.to_param().into()
            }
        }
    };

    (
        #[value_type=$value_type:ty]
        #[flag_type=$flag_type:ty]
        $(#[$enum_attrs:meta])*
        $vis:vis enum $enum_name:ident {
            $(
                #[term(cv=$cv:ident, accession=$accession:literal, name=$term_name:literal, flags=$flags:tt, parents=$parents:tt)]
                $(#[$variant_attrs:meta])*
                $variant:ident($value_type2:ty)
            ),*
            $(,)?
        }
    ) => {
        $(#[$enum_attrs])*
        $vis enum $enum_name {
            $(
                $(#[$variant_attrs])*
                $variant($value_type),
            )*
        }

        /// These methods are part of the controlled vocabulary mapping
        impl $enum_name {

            /// Retrieve the accession number for this term, independent of its controlled vocabulary
            pub const fn accession(&self) -> $crate::params::AccessionIntCode {
                match self {
                    $(Self::$variant(_) => $accession,)*
                }
            }

            /// Retrieve the controlled vocabulary this term belongs to
            pub const fn controlled_vocabulary(&self) -> $crate::params::ControlledVocabulary {
                match self {
                    $(Self::$variant(_) => $crate::params::ControlledVocabulary::$cv,)*
                }
            }

            /// Retrieve the plain text human readable name for this term
            pub const fn name(&self) -> &'static str {
                match self {
                    $(Self::$variant(_) => $term_name,)*
                }
            }

            /// Attempt to map a string by name to retrieve one of the terms from this
            /// set.
            ///
            /// If no match is found, [`None`] is returned.
            pub fn from_name(name: &str) -> Option<Self> {
                match name {
                    $($term_name => Some(Self::$variant(Default::default())),)*
                    _ => None
                }
            }


            /// Attempt to map the numeric accession number to retrieve one of the terms from this
            /// set.
            ///
            /// If no match is found, [`None`] is returned.
            pub const fn from_accession(accession: $crate::params::AccessionIntCode, value: $value_type) -> Option<Self> {
                match accession {
                    $($accession => Some(Self::$variant(value)),)*
                    _ => None
                }
            }

            /// Convert this term into a [`ParamCow`](crate::params::ParamCow) without a value.
            pub const fn to_param(self) -> $crate::params::ParamCow<'static> {
                $crate::params::ParamCow::const_new(
                    self.name(),
                    $crate::params::ValueRef::Empty,
                    Some(self.accession()),
                    Some(self.controlled_vocabulary()),
                    $crate::params::Unit::Unknown
                )
            }

            /// Convert this term into a [`ParamCow`](crate::params::ParamCow) a specific boxed value.
            pub const fn to_param_value(self, value: $crate::params::ValueRef<'static>) -> $crate::params::ParamCow<'static> {
                $crate::params::ParamCow::const_new(
                    self.name(),
                    value,
                    Some(self.accession()),
                    Some(self.controlled_vocabulary()),
                    $crate::params::Unit::Unknown
                )
            }

            /// Convert a [`CURIE`]($crate::params::CURIE) by accession.
            ///
            /// If no match is found, [`None`] is returned.
            pub fn from_curie(curie: &$crate::params::CURIE, value: f32) -> Option<Self> {
                if matches!(curie.controlled_vocabulary, $crate::params::ControlledVocabulary::MS) {
                    Self::from_accession(curie.accession, value)
                } else {
                    None
                }
            }

            /// Attempt to convert a [`ParamCow`](crate::params::ParamCow) to a term from this set.
            ///
            /// If no match is found, [`None`] is returned.
            ///
            /// # Note
            /// This method can be called in `const` contexts, requiring the type be [`ParamCow`](crate::params::ParamCow) with a `'static`
            /// lifetime parameter, but the regular [`From`] trait is implemented for all [`ParamLike`](crate::params::ParamLike) types.
            pub fn from_param(p: &$crate::params::ParamCow<'static>) -> Option<Self> {
                if let Some(acc) = p.accession {
                    Self::from_accession(acc, p.value.clone().into())
                } else {
                    None
                }
            }

            /// Retrieve a term set specific set of flags
            pub fn flags(&self) -> $flag_type {
                match self {
                    $(Self::$variant(_) => $flags.into(),)*
                }
            }

            /// Retrieve the list of zero or more terms in the set which are
            /// parents of this term.
            pub fn parents(&self) -> Vec<Self> {
                match self {
                    $(Self::$variant(_) => $parents.iter().flat_map(|s: &&str| {
                        let curie = s.parse::<$crate::params::CURIE>().unwrap();
                        Self::from_accession(curie.accession, Default::default())
                    }).collect(),)*
                }
            }

        }

        impl<P> From<P> for $enum_name where P: $crate::params::ParamLike {
            fn from(value: P) -> Self {
                Self::from_accession(
                    value.accession().expect(
                        concat!("Cannot convert an uncontrolled parameter to ", stringify!($enum_name))),
                    value.value().into()
                ).unwrap_or_else(
                    || panic!(
                        "Could not map {:?}:{} to {}",
                        value.controlled_vocabulary().unwrap(),
                        value.accession().unwrap(),
                        stringify!($enum_name)
                    )
                )
            }
        }

        impl From<$enum_name> for $crate::params::ParamCow<'static> {
            fn from(value: $enum_name) -> Self {
                value.to_param()
            }
        }

        impl From<$enum_name> for $crate::params::Param {
            fn from(value: $enum_name) -> Self {
                value.to_param().into()
            }
        }

        impl From<&$enum_name> for $crate::params::ParamCow<'static> {
            fn from(value: &$enum_name) -> Self {
                value.to_param()
            }
        }

        impl From<&$enum_name> for $crate::params::Param {
            fn from(value: &$enum_name) -> Self {
                value.to_param().into()
            }
        }
    };
}

bitflags::bitflags! {
    #[doc="A bit mask encoding the different value types a [`Param`](crate::params::Param) can take on, associated with a term."]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct ValueType: u16 {
        const NoType = 0;
        const String = 0b00000001;
        const Integer = 0b00000010;
        const Float = 0b00000100;
        const Double = 0b00001000;
        const NonNegativeInteger = 0b00010000;
        const PositiveInteger = 0b00100000;
        const DateTime = 0b01000000;
        const Boolean = 0b10000000;

        const ListOf = 0b1000000000000000;
    }
}

impl From<u16> for ValueType {
    fn from(value: u16) -> Self {
        Self::from_bits_retain(value)
    }
}

impl ValueType {
    pub fn parse(&self, s: String) -> Result<Value, ParamValueParseError> {
        let v = match *self {
            Self::Integer | Self::NonNegativeInteger | Self::PositiveInteger => Value::Int(
                s.parse()
                    .map_err(|_| ParamValueParseError::FailedToExtractInt(Some(s)))?,
            ),
            Self::Float | Self::Double => Value::Float(
                s.parse()
                    .map_err(|_| ParamValueParseError::FailedToExtractFloat(Some(s)))?,
            ),
            _ => Value::String(s),
        };
        Ok(v)
    }

    pub fn parse_str<'a>(&self, s: &'a str) -> Result<ValueRef<'a>, ParamValueParseError> {
        let v = match *self {
            Self::Integer | Self::NonNegativeInteger | Self::PositiveInteger => ValueRef::Int(
                s.parse()
                    .map_err(|_| ParamValueParseError::FailedToExtractInt(Some(s.to_string())))?,
            ),
            Self::Float | Self::Double => ValueRef::Float(
                s.parse()
                    .map_err(|_| ParamValueParseError::FailedToExtractFloat(Some(s.to_string())))?,
            ),
            _ => ValueRef::String(Cow::Borrowed(s)),
        };
        Ok(v)
    }
}
