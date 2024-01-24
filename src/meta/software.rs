use crate::impl_param_described;
use crate::params::ParamList;

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
