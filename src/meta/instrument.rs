use crate::impl_param_described;
use crate::params::ParamList;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[derive(Default)]
pub enum ComponentType {
    Analyzer,
    IonSource,
    Detector,
    #[default]
    Unknown,
}



#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct Component {
    pub component_type: ComponentType,
    pub order: u8,
    pub params: ParamList,
}

#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct InstrumentConfiguration {
    pub components: Vec<Component>,
    pub params: ParamList,
    pub software_reference: String,
    pub id: u32,
}

impl_param_described!(InstrumentConfiguration, Component);
