/// Metadata describing the experiment that does not belong in any other section
/// that covers some default options.
#[derive(Debug, Default, PartialEq, Hash, Eq, Clone)]
pub struct MassSpectrometryRun {
    pub id: Option<String>,
    pub default_data_processing_id: Option<String>,
    pub default_instrument_id: Option<u32>,
    pub default_source_file_id: Option<String>,
    pub start_time: Option<String>,
}

impl MassSpectrometryRun {
    pub fn new(
        id: Option<String>,
        default_data_processing_id: Option<String>,
        default_instrument_id: Option<u32>,
        default_source_file_id: Option<String>,
        start_time: Option<String>,
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
