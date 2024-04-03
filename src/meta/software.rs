use std::collections::HashSet;

use crate::params::{ControlledVocabulary, ParamList};
use crate::{impl_param_described, Param};

/// A piece of software that was associated with the acquisition, transformation or otherwise
/// processing of mass spectrometry data. See <https://peptideatlas.org/tmp/mzML1.1.0.html#software>
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Software {
    /// A unique identifier for the software within processing metadata
    pub id: String,
    /// A string denoting a particular software version, but does no guarantee is given for its format
    pub version: String,
    /// Any associated vocabulary terms, including actual software name and type
    pub params: ParamList,
}

impl Software {
    pub fn new(id: String, version: String, params: ParamList) -> Self {
        Self {
            id,
            version,
            params,
        }
    }

    /// Find a unique identifier from an iterator over software IDs
    pub fn find_unique_id<'a>(id_stem: &str, softwares: impl IntoIterator<Item = &'a Self>) -> String {
        let software_ids: HashSet<_> = softwares.into_iter().map(|sw| &sw.id).collect();
        (0..)
            .into_iter()
            .map(|i| format!("{id_stem}_{i}"))
            .filter(|s| !software_ids.contains(s))
            .next()
            .unwrap()
    }
}

/// Create an instance of "custom unreleased software tool" with name `name`
pub fn custom_software_name(name: &str) -> Param {
    ControlledVocabulary::MS.param_val(1000799, "custom unreleased software tool", name)
}

impl_param_described!(Software);
