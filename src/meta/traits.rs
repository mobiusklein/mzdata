use std::collections::HashMap;

use super::{
    DataProcessing, FileDescription, InstrumentConfiguration, MassSpectrometryRun, Software,
};

/// Mass spectrometry data files have several facets of descriptive metadata
pub trait MSDataFileMetadata {
    /// The series of [`DataProcessing`] workflows applied to spectra in
    /// this data file.
    fn data_processings(&self) -> &Vec<DataProcessing>;
    /// A mapping over different [`InstrumentConfiguration`] modes that spectra
    /// were acquired under.
    fn instrument_configurations(&self) -> &HashMap<u32, InstrumentConfiguration>;
    /// A description of the contents and the sources for this mass spectrometry
    /// data file.
    fn file_description(&self) -> &FileDescription;
    /// The series of [`Software`] applied to the data file to apply different
    /// [`DataProcessing`] methods.
    fn softwares(&self) -> &Vec<Software>;

    /// Mutably access the [`DataProcessing`] list for this data file
    fn data_processings_mut(&mut self) -> &mut Vec<DataProcessing>;
    /// Mutably access the [`InstrumentConfiguration`] mapping for this data file
    fn instrument_configurations_mut(&mut self) -> &mut HashMap<u32, InstrumentConfiguration>;
    /// Mutably access the [`FileDescription`] description of the contents and the
    /// sources for this mass spectrometry data file.
    fn file_description_mut(&mut self) -> &mut FileDescription;
    /// Mutably access the list of [`Software`] of this data file.
    fn softwares_mut(&mut self) -> &mut Vec<Software>;

    /// Copy the metadata from another [`MSDataFileMetadata`] implementation into
    /// this one.
    fn copy_metadata_from<T: MSDataFileMetadata>(&mut self, source: &T) {
        *self.data_processings_mut() = source.data_processings().clone();
        *self.instrument_configurations_mut() = source.instrument_configurations().clone();
        *self.file_description_mut() = source.file_description().clone();
        *self.softwares_mut() = source.softwares().clone();

        match source.run_description() {
            Some(run) => {
                let desc = self.run_description_mut();
                if let Some(r) = desc {
                    *r = run.clone();
                };
            }
            None => {
                let mut desc = self.run_description_mut();
                desc.take();
            }
        }
    }

    /// A hint about how many spectra are in this data file
    fn spectrum_count_hint(&self) -> Option<u64> {
        None
    }

    /// Access the [`MassSpectrometryRun`] metadata record if it is available
    fn run_description(&self) -> Option<&MassSpectrometryRun> {
        None
    }

    /// Mutably access the [`MassSpectrometryRun`] metadata record if it is available
    fn run_description_mut(&mut self) -> Option<&mut MassSpectrometryRun> {
        None
    }

    /// Get the name of the primary source file, if available
    fn source_file_name(&self) -> Option<&str> {
        self.file_description().source_files.first().map(|s| s.name.as_str())
    }
}


#[macro_export]
/// Assumes a field for the non-`Option` facets of the [`MSDataFileMetadata`]
/// implementation are present.
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


#[macro_export]
/// Delegates the implementation of [`MSDataFileMetadata`] to an field
macro_rules! delegate_impl_metadata_trait {
    ($src:tt) => {

        fn data_processings(&self) -> &Vec<crate::meta::DataProcessing> {
            self.$src.data_processings()
        }

        fn instrument_configurations(&self) -> &std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
            self.$src.instrument_configurations()
        }

        fn file_description(&self) -> &crate::meta::FileDescription {
            self.$src.file_description()
        }

        fn softwares(&self) -> &Vec<crate::meta::Software> {
            self.$src.softwares()
        }

        fn data_processings_mut(&mut self) -> &mut Vec<crate::meta::DataProcessing> {
            self.$src.data_processings_mut()
        }

        fn instrument_configurations_mut(&mut self) -> &mut std::collections::HashMap<u32, crate::meta::InstrumentConfiguration> {
            self.$src.instrument_configurations_mut()
        }

        fn file_description_mut(&mut self) -> &mut crate::meta::FileDescription {
            self.$src.file_description_mut()
        }

        fn softwares_mut(&mut self) -> &mut Vec<crate::meta::Software> {
            self.$src.softwares_mut()
        }

        fn spectrum_count_hint(&self) -> Option<u64> {
            self.$src.spectrum_count_hint()
        }

        fn run_description(&self) -> Option<&crate::meta::MassSpectrometryRun> {
            self.$src.run_description()
        }

        fn run_description_mut(&mut self) -> Option<&mut crate::meta::MassSpectrometryRun> {
            self.$src.run_description_mut()
        }

    };
}