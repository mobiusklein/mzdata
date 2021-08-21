use std::collections::HashMap;

use super::{DataProcessing, FileDescription, InstrumentConfiguration, Software};

pub trait MSDataFileMetadata {
    fn data_processings(&self) -> &Vec<DataProcessing>;
    fn instrument_configurations(&self) -> &HashMap<String, InstrumentConfiguration>;
    fn file_description(&self) -> &FileDescription;
    fn softwares(&self) -> &Vec<Software>;
}

#[macro_export]
macro_rules! impl_metadata_trait {
    ($t:ident, $R:ident, $C:ident, $D:ident) => {
        impl<R: $R, C: $C, D: $D> $crate::meta::MSDataFileMetadata for $t<R, C, D> {
            fn data_processings(&self) -> &Vec<$crate::meta::DataProcessing> {
                &self.data_processings
            }

            fn instrument_configurations(&self) -> &HashMap<String, InstrumentConfiguration> {
                &self.instrument_configurations
            }
            fn file_description(&self) -> &FileDescription {
                &self.file_description
            }
            fn softwares(&self) -> &Vec<Software> {
                &self.softwares
            }
        }
    };
}
