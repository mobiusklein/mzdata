use crate::{
    meta::{Component, ComponentType, InletTypeTerm, IonizationTypeTerm},
    params::ParamDescribed,
};

#[allow(non_camel_case_types)]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InstrumentSource {
    AlsoUnknown = 0,
    ESI = 1,
    APCI = 2,
    NanoESIOffline = 3,
    NanoESIOnline = 4,
    APPI = 5,
    AP_MALDI = 6,
    MALDI = 7,
    MultiMode = 8,
    NanoFlowESI = 9,
    Ultraspray = 10,
    CaptiveSpray = 11,
    EI = 16,
    GC_APCI = 17,
    VIP_HESI = 18,
    VIP_APCI = 19,
    #[default]
    Unknown = 255,
}

impl InstrumentSource {
    pub fn to_component(&self) -> Component {
        let mut comp = Component::default();
        comp.order = 1;
        comp.component_type = ComponentType::IonSource;

        match self {
            Self::ESI | Self::MultiMode | Self::Ultraspray | Self::VIP_HESI => {
                comp.add_param(IonizationTypeTerm::ElectrosprayIonization.into());
                comp.add_param(InletTypeTerm::ElectrosprayInlet.into());
            }
            Self::NanoESIOffline | Self::NanoESIOnline | Self::NanoFlowESI | Self::CaptiveSpray => {
                comp.add_param(IonizationTypeTerm::Nanoelectrospray.into());
                comp.add_param(InletTypeTerm::NanosprayInlet.into());
            }
            Self::APCI | Self::GC_APCI | Self::VIP_APCI => {
                comp.add_param(IonizationTypeTerm::AtmosphericPressureChemicalIonization.into());
            }
            Self::APPI => {
                comp.add_param(IonizationTypeTerm::AtmosphericPressurePhotoionization.into());
            }
            Self::EI => {
                comp.add_param(IonizationTypeTerm::ElectronIonization.into());
            }
            Self::MALDI => {
                comp.add_param(IonizationTypeTerm::MatrixAssistedLaserDesorptionIonization.into());
            }
            Self::AP_MALDI => {
                comp.add_param(
                    IonizationTypeTerm::AtmosphericPressureMatrixAssistedLaserDesorptionIonization
                        .into(),
                );
            }
            Self::AlsoUnknown | Self::Unknown => {}
        }

        comp
    }
}

impl From<u8> for InstrumentSource {
    fn from(value: u8) -> Self {
        match value {
            0 => Self::AlsoUnknown,
            1 => Self::ESI,
            2 => Self::APCI,
            3 => Self::NanoESIOffline,
            4 => Self::NanoESIOnline,
            5 => Self::APPI,
            6 => Self::AP_MALDI,
            7 => Self::MALDI,
            8 => Self::MultiMode,
            9 => Self::NanoFlowESI,
            10 => Self::Ultraspray,
            11 => Self::CaptiveSpray,
            16 => Self::EI,
            17 => Self::GC_APCI,
            18 => Self::VIP_HESI,
            19 => Self::VIP_APCI,
            255 => Self::Unknown,
            _ => Self::Unknown,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MsMsType {
    MS1 = 0,
    MRM = 2,
    DDAPASEF = 8,
    DIAPASEF = 9,
    PRMPASEF = 10,

    Unknown = -1,
}

impl MsMsType {
    pub const fn ms_level(&self) -> u8 {
        match self {
            MsMsType::MS1 => 1,
            MsMsType::MRM => 2,
            MsMsType::DDAPASEF => 2,
            MsMsType::DIAPASEF => 2,
            MsMsType::PRMPASEF => 2,
            MsMsType::Unknown => 0,
        }
    }
}

impl From<u8> for MsMsType {
    fn from(value: u8) -> Self {
        match value {
            0 => Self::MS1,
            2 => Self::MRM,
            8 => Self::DDAPASEF,
            9 => Self::DIAPASEF,
            10 => Self::PRMPASEF,
            _ => Self::Unknown,
        }
    }
}
