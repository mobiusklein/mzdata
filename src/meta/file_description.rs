use crate::params::{ParamList, Param, ParamDescribed};
use crate::impl_param_described;


#[derive(Debug, Clone, Default)]
pub struct SourceFile {
    pub name: String,
    pub location: String,
    pub id: String,
    pub file_format: Option<Param>,
    pub id_format: Option<Param>,
    pub params: ParamList,
}

#[derive(Debug, Clone, Default)]
pub struct FileDescription {
    pub contents: ParamList,
    pub source_files: Vec<SourceFile>,
}

impl_param_described!(SourceFile);


impl ParamDescribed for FileDescription {
    fn params(&self) -> &ParamList {
        &self.contents
    }

    fn params_mut(&mut self) -> &mut ParamList {
        &mut self.contents
    }
}