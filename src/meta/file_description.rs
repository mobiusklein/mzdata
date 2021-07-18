
use crate::params::{Param, ParamList};


pub struct SourceFile {
    pub name: String,
    pub location: String,
    pub id: String,
    pub file_format: Option<Param>,
    pub id_format: Option<Param>,
    pub params: ParamList
}


pub struct FileDescription {
    pub contents: ParamList,
    pub source_files: Vec<SourceFile>
}
