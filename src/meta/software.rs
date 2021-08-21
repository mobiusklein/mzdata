use crate::impl_param_described;
use crate::params::ParamList;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Software {
    pub id: String,
    pub version: String,
    pub params: ParamList,
}

impl_param_described!(Software);
