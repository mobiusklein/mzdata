use crate::params::{Param, ParamDescribed};

/// Description of the acquisition settings of the instrument prior to the start of the run
#[derive(Debug, Default, Clone)]
pub struct ScanSettings {
    /// A unique identifier
    pub id: String,
    /// List with the source files containing the acquisition settings
    pub source_file_refs: Vec<String>,
    /// Target list (or 'inclusion list') configured prior to the run
    pub targets: Vec<Vec<Param>>,
    /// The controlled vocabulary and user parameters of the settings
    pub params: Vec<Param>,
}

impl ScanSettings {
    pub fn new(
        id: String,
        params: Vec<Param>,
        source_file_refs: Vec<String>,
        targets: Vec<Vec<Param>>,
    ) -> Self {
        Self {
            id,
            params,
            source_file_refs,
            targets,
        }
    }
}

impl ParamDescribed for ScanSettings {
    fn params(&self) -> &[Param] {
        &self.params
    }

    fn params_mut(&mut self) -> &mut crate::ParamList {
        &mut self.params
    }
}
