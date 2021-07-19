use crate::params::ParamList;

#[derive(Debug)]
pub enum ComponentType {
    Analyzer,
    IonSource,
    Detector,
    Unknown,
}

impl Default for ComponentType {
    fn default() -> ComponentType {
        ComponentType::Unknown
    }
}

#[derive(Default, Debug)]
pub struct Component {
    pub component_type: ComponentType,
    pub order: u8,
    pub params: ParamList,
}

#[derive(Default, Debug)]
pub struct InstrumentConfiguration {
    pub components: Vec<Component>,
    pub params: ParamList,
    pub software_reference: String,
    pub id: String,
}
