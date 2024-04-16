/*! Metadata describing mass spectrometry data files and their contents.
 */
#[macro_use]
mod file_description;
mod data_processing;
mod instrument;
mod software;
mod run;
#[macro_use]
mod traits;

pub use data_processing::{DataProcessing, ProcessingMethod, FormatConversion, DataTransformationAction, DataProcessingAction};
pub use software::{Software, custom_software_name};

pub use file_description::{FileDescription, SourceFile};

pub use instrument::{Component, ComponentType, InstrumentConfiguration};

pub use traits::MSDataFileMetadata;
pub use run::MassSpectrometryRun;

#[macro_export]
macro_rules! cvmap {
    (
        #[flag_type=$flag_type:ty]
        $(#[$enum_attrs:meta])*
        $vis:vis enum $enum_name:ident {
            $(
                #[term(cv=$cv:ident, accession=$accession:literal, name=$term_name:literal, flags=$flags:tt)]
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

        impl $enum_name {
            pub const fn accession(&self) -> u32 {
                match self {
                    $(Self::$variant => $accession,)*
                }
            }

            pub const fn controlled_vocabulary(&self) -> $crate::params::ControlledVocabulary {
                match self {
                    $(Self::$variant => $crate::params::ControlledVocabulary::$cv,)*
                }
            }

            pub const fn name(&self) -> &'static str {
                match self {
                    $(Self::$variant => $term_name,)*
                }
            }

            pub fn from_name(name: &str) -> Option<Self> {
                match name {
                    $($term_name => Some(Self::$variant),)*
                    _ => None
                }
            }

            pub const fn from_accession(accession: u32) -> Option<Self> {
                match accession {
                    $($accession => Some(Self::$variant),)*
                    _ => None
                }
            }

            pub const fn to_param(&self) -> $crate::params::ParamCow<'static> {
                $crate::params::ParamCow::const_new(
                    self.name(),
                    $crate::params::ValueRef::Empty,
                    Some(self.accession()),
                    Some(self.controlled_vocabulary()),
                    $crate::params::Unit::Unknown
                )
            }

            pub const fn from_param(p: &$crate::params::ParamCow<'static>) -> Option<Self> {
                if let Some(acc) = p.accession {
                    Self::from_accession(acc)
                } else {
                    None
                }
            }

            pub fn flags(&self) -> $flag_type {
                match self {
                    $(Self::$variant => $flags.into(),)*
                }
            }

        }

        impl<P> From<P> for $enum_name where P: $crate::params::ParamLike {
            fn from(value: P) -> Self {
                Self::from_accession(value.accession().unwrap()).unwrap()
            }
        }

        impl From<$enum_name> for $crate::params::ParamCow<'static> {
            fn from(value: $enum_name) -> Self {
                value.to_param()
            }
        }
    };
}