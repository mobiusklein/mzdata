use crate::impl_param_described;
use crate::params::{Param, ParamDescribed, ParamList, CURIE, ControlledVocabulary};

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SourceFile {
    pub name: String,
    pub location: String,
    pub id: String,
    pub file_format: Option<Param>,
    pub id_format: Option<Param>,
    pub params: ParamList,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FileDescription {
    pub contents: ParamList,
    pub source_files: Vec<SourceFile>,
}

impl FileDescription {
    pub fn new(contents: ParamList, source_files: Vec<SourceFile>) -> Self {
        Self {
            contents,
            source_files,
        }
    }

    /// Checks to see if the "MS1 spectrum" term is present in the file contents
    ///
    /// **Note**: This does not actually inspect the spectra in the file, only the metadata,
    /// which may be incorrect/missing.
    pub fn has_ms1_spectra(&self) -> bool {
        self.get_param_by_curie(&CURIE::new(ControlledVocabulary::MS, 1000579)).is_some()
    }

    /// Checks to see if the "MSn spectrum" term is present in the file contents.
    ///
    /// **Note**: This does not actually inspect the spectra in the file, only the metadata,
    /// which may be incorrect/missing.
    pub fn has_msn_spectra(&self) -> bool {
        self.get_param_by_curie(&CURIE::new(ControlledVocabulary::MS, 1000580)).is_some()
    }

    pub fn has_contents(&self) -> bool {
        !self.contents.is_empty()
    }
}

impl_param_described!(SourceFile);

impl ParamDescribed for FileDescription {
    fn params(&self) -> &[Param] {
        &self.contents
    }

    fn params_mut(&mut self) -> &mut ParamList {
        &mut self.contents
    }
}
