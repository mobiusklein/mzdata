use std::fmt::Display;

use crate::impl_param_described;
use crate::params::ParamList;

/// A distinguishing tag describing the part of an instrument a [`Component`] refers to
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum ComponentType {
    /// A mass analyzer
    Analyzer,
    /// A source for ions
    IonSource,
    /// An abundance measuring device
    Detector,
    #[default]
    Unknown,
}

impl Display for ComponentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// A description of a combination of parts that are described as part of an [`InstrumentConfiguration`]
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct Component {
    pub component_type: ComponentType,
    pub order: u8,
    pub params: ParamList,
}

/// A series of mass spectrometer components that together were engaged to acquire a mass spectrum
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct InstrumentConfiguration {
    /// The set of components involved
    pub components: Vec<Component>,
    /// A set of parameters that describe the instrument such as the model name or serial number
    pub params: ParamList,
    /// A reference to the data acquisition software involved in processing this configuration
    pub software_reference: String,
    /// A unique identifier translated to an ordinal identifying this configuration
    pub id: u32,
}

impl_param_described!(InstrumentConfiguration, Component);
