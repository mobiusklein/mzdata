use crate::impl_param_described;
use crate::params::ParamList;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ProcessingMethod {
    pub order: i8,
    pub software_reference: String,
    pub params: ParamList,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DataProcessing {
    pub id: String,
    pub methods: Vec<ProcessingMethod>,
}

impl_param_described!(ProcessingMethod);

impl DataProcessing {
    pub fn push(&mut self, method: ProcessingMethod) {
        self.methods.push(method)
    }

    pub fn iter(&self) -> std::slice::Iter<ProcessingMethod> {
        self.methods.iter()
    }

    pub fn len(&self) -> usize {
        self.methods.len()
    }

    pub fn is_empty(&self) -> bool {
        self.methods.is_empty()
    }
}
