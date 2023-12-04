#[macro_use]
mod file_description;
mod data_processing;
mod instrument;
mod software;
#[macro_use]
mod traits;

pub use crate::meta::data_processing::{DataProcessing, ProcessingMethod};
pub use crate::meta::file_description::{FileDescription, SourceFile};
pub use crate::meta::instrument::{Component, ComponentType, InstrumentConfiguration};
pub use crate::meta::software::Software;
pub use crate::meta::traits::MSDataFileMetadata;
