#[macro_use]
pub mod file_description;
pub mod instrument;
pub mod software;
pub mod data_processing;

pub use crate::meta::file_description::{SourceFile, FileDescription};
pub use crate::meta::instrument::{ComponentType, Component, InstrumentConfiguration};
pub use crate::meta::software::Software;
pub use crate::meta::data_processing::{ProcessingMethod, DataProcessing};