use chrono::{DateTime, FixedOffset};

/// Metadata describing the experiment that does not belong in any other section
/// that covers some default options.
#[derive(Debug, Default, PartialEq, Hash, Eq, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MassSpectrometryRun {
    /// A unique identifier for the entire MS run.
    pub id: Option<String>,
    /// The default [`DataProcessing`](crate::DataProcessing) identifier to use, if any, for entities that do not
    /// list their own method.
    pub default_data_processing_id: Option<String>,
    /// The default [`InstrumentConfiguration`](crate::InstrumentConfiguration) identifier to use, if any, for entities that do
    /// not list their own configuration.
    pub default_instrument_id: Option<u32>,
    /// The default [`SourceFile`](crate::SourceFile) identifier that entries come from, for entities that do not list their own
    /// source.
    pub default_source_file_id: Option<String>,
    /// The date and time acquisition began.
    pub start_time: Option<DateTime<FixedOffset>>,
}

impl MassSpectrometryRun {
    pub fn new(
        id: Option<String>,
        default_data_processing_id: Option<String>,
        default_instrument_id: Option<u32>,
        default_source_file_id: Option<String>,
        start_time: Option<DateTime<FixedOffset>>,
    ) -> Self {
        Self {
            id,
            default_data_processing_id,
            default_instrument_id,
            default_source_file_id,
            start_time,
        }
    }
}
