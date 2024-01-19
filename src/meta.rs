/*! Metadata describing mass spectrometry data files and their contents.
 */
#[macro_use]
mod file_description;
mod data_processing;
mod instrument;
mod software;
mod run;
#[macro_use]
mod traits;

pub use data_processing::{DataProcessing, ProcessingMethod};
pub use file_description::{FileDescription, SourceFile};
pub use instrument::{Component, ComponentType, InstrumentConfiguration};
pub use software::Software;
pub use traits::MSDataFileMetadata;
pub use run::MassSpectrometryRun;
