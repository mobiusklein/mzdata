#[macro_use]
pub mod file_description;
pub mod data_processing;
pub mod instrument;
pub mod software;
#[macro_use]
pub mod traits;

pub use crate::meta::data_processing::{DataProcessing, ProcessingMethod};
pub use crate::meta::file_description::{FileDescription, SourceFile};
pub use crate::meta::instrument::{Component, ComponentType, InstrumentConfiguration};
pub use crate::meta::software::Software;
pub use crate::meta::traits::MSDataFileMetadata;
