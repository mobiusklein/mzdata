use crate::impl_param_described;
use crate::params::ParamList;

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

impl_param_described!(Software);
