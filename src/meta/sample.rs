use crate::impl_param_described;
use crate::params::{ParamDescribed, ParamList, ParamValue, ValueRef};


#[derive(Debug, Default, Clone, PartialEq)]
pub struct Sample {
    pub id: String,
    pub name: Option<String>,
    pub params: ParamList
}

impl Sample {
    pub fn new(id: String, name: Option<String>, params: ParamList) -> Self {
        Self { id, name, params }
    }

    pub fn number(&self) -> Option<ValueRef> {
        self.get_param_by_curie(&crate::curie!(MS:1000001)).map(|p| p.value.as_ref())
    }
}


impl_param_described!(Sample);