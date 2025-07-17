use crate::impl_param_described;
use crate::params::{ParamDescribed, ParamList};


/// A sample that was used to generate a MS dataset
#[derive(Debug, Default, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Sample {
    /// A "unique" identifier
    pub id: String,
    /// The sample's human-friendly name
    pub name: Option<String>,
    /// The controlled vocabulary and user parameters of the sample
    pub params: ParamList
}

impl Sample {
    pub fn new(id: String, name: Option<String>, params: ParamList) -> Self {
        Self { id, name, params }
    }

    crate::find_param_method!(number, &crate::curie!(MS:1000001), "Find the sample number, if it is present");
    crate::find_param_method!(batch, &crate::curie!(MS:1000053), "Find the sample batch, if it is present");
}


impl_param_described!(Sample);