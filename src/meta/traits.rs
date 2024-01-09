use std::collections::HashMap;

use super::{
    DataProcessing, FileDescription, InstrumentConfiguration, MassSpectrometryRun, Software,
};

pub trait MSDataFileMetadata {
    fn data_processings(&self) -> &Vec<DataProcessing>;
    // TODO: Switch this to an IndexMap
    fn instrument_configurations(&self) -> &HashMap<u32, InstrumentConfiguration>;
    fn file_description(&self) -> &FileDescription;
    fn softwares(&self) -> &Vec<Software>;

    fn data_processings_mut(&mut self) -> &mut Vec<DataProcessing>;
    // TODO: Switch this to an IndexMap
    fn instrument_configurations_mut(&mut self) -> &mut HashMap<u32, InstrumentConfiguration>;
    fn file_description_mut(&mut self) -> &mut FileDescription;
    fn softwares_mut(&mut self) -> &mut Vec<Software>;

    fn copy_metadata_from<T: MSDataFileMetadata>(&mut self, source: &T) {
        *self.data_processings_mut() = source.data_processings().clone();
        *self.instrument_configurations_mut() = source.instrument_configurations().clone();
        *self.file_description_mut() = source.file_description().clone();
        *self.softwares_mut() = source.softwares().clone();

        match source.run_description() {
            Some(run) => {
                let desc = self.run_description_mut();
                desc.and_then(|r| {
                    *r = run.clone();
                    Some(())
                });
            }
            None => {
                let mut desc = self.run_description_mut();
                desc.take();
            }
        }
    }

    fn spectrum_count_hint(&self) -> Option<u64> {
        None
    }

    fn run_description(&self) -> Option<&MassSpectrometryRun> {
        None
    }

    fn run_description_mut(&mut self) -> Option<&mut MassSpectrometryRun> {
        None
    }
}

#[macro_export]
macro_rules! impl_metadata_trait {
    () => {
        fn data_processings(&self) -> &Vec<DataProcessing> {
            &self.data_processings
        }

        fn instrument_configurations(&self) -> &HashMap<u32, InstrumentConfiguration> {
            &self.instrument_configurations
        }
        fn file_description(&self) -> &FileDescription {
            &self.file_description
        }
        fn softwares(&self) -> &Vec<Software> {
            &self.softwares
        }

        fn data_processings_mut(&mut self) -> &mut Vec<DataProcessing> {
            &mut self.data_processings
        }

        fn instrument_configurations_mut(&mut self) -> &mut HashMap<u32, InstrumentConfiguration> {
            &mut self.instrument_configurations
        }

        fn file_description_mut(&mut self) -> &mut FileDescription {
            &mut self.file_description
        }

        fn softwares_mut(&mut self) -> &mut Vec<Software> {
            &mut self.softwares
        }
    };
}
