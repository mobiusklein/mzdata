#![allow(unused)]
//! Thank you ProteoWizard for building all of these look up tables.
use std::fmt::Display;

use crate::meta::{Component, ComponentType, DetectorTypeTerm, InstrumentConfiguration, IonizationTypeTerm, MassAnalyzerTerm};
use crate::params::ParamDescribed;
use crate::{params::ControlledVocabulary, Param};

macro_rules! param {
    ($name:expr, $acc:expr) => {
        ControlledVocabulary::MS.const_param_ident($name, $acc)
    };
}

#[allow(unused)]
#[derive(Debug)]
enum MatchType {
    Exact,
    Contains,
    StartsWith,
    EndsWith,
    ExactNoSpaces,
}

#[allow(unused, non_camel_case_types, clippy::upper_case_acronyms)]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InstrumentModelType {
    #[default]
    Unknown = -1,

    // Finnigan MAT
    MAT253,
    MAT900XP,
    MAT900XP_Trap,
    MAT95XP,
    MAT95XP_Trap,
    SSQ_7000,
    TSQ_7000,
    TSQ,

    // Thermo Electron
    Element_2,

    // Thermo Finnigan
    Delta_Plus_Advantage,
    Delta_Plus_XP,
    LCQ_Advantage,
    LCQ_Classic,
    LCQ_Deca,
    LCQ_Deca_XP_Plus,
    Neptune,
    DSQ,
    PolarisQ,
    Surveyor_MSQ,
    Tempus_TOF,
    Trace_DSQ,
    Triton,

    // Thermo Scientific
    LTQ,
    LTQ_Velos,
    LTQ_Velos_ETD,
    LTQ_Velos_Plus,
    LTQ_FT,
    LTQ_FT_Ultra,
    LTQ_Orbitrap,
    LTQ_Orbitrap_Classic,
    LTQ_Orbitrap_Discovery,
    LTQ_Orbitrap_XL,
    LTQ_Orbitrap_Velos,
    LTQ_Orbitrap_Velos_Pro,
    LTQ_Orbitrap_Elite,
    LXQ,
    LCQ_Fleet,
    ITQ_700,
    ITQ_900,
    ITQ_1100,
    GC_Quantum,
    LTQ_XL,
    LTQ_XL_ETD,
    LTQ_Orbitrap_XL_ETD,
    DFS,
    DSQ_II,
    ISQ,
    MALDI_LTQ_XL,
    MALDI_LTQ_Orbitrap,
    TSQ_Quantum,
    TSQ_Quantum_Access,
    TSQ_Quantum_Ultra,
    TSQ_Quantum_Ultra_AM,
    TSQ_Vantage,
    Element_XR,
    Element_GD,
    GC_IsoLink,
    Exactive,
    Exactive_Plus,
    Q_Exactive,
    Q_Exactive_Plus,
    Q_Exactive_HF,
    Q_Exactive_HF_X,
    Q_Exactive_UHMR,
    Surveyor_PDA,
    Accela_PDA,
    Orbitrap_Fusion,
    Orbitrap_Fusion_Lumos,
    Orbitrap_Fusion_ETD,
    Orbitrap_Ascend,
    Orbitrap_ID_X,
    TSQ_Quantiva,
    TSQ_Endura,
    TSQ_Altis,
    TSQ_Altis_Plus,
    TSQ_Quantis,
    TSQ_8000_Evo,
    TSQ_9000,
    Orbitrap_Exploris_120,
    Orbitrap_Exploris_240,
    Orbitrap_Exploris_480,
    Orbitrap_Exploris_GC_240,
    Orbitrap_Eclipse,
    Orbitrap_GC,
    Orbitrap_Astral,
}

impl Display for InstrumentModelType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let k = format!("{:?}", self).replace("_", " ");
        write!(f, "{k}")
    }
}

#[allow(unused)]
impl InstrumentModelType {
    pub fn to_param(self) -> Param {
        let val = match self {
            InstrumentModelType::Unknown => {
                param!("Thermo Fisher Scientific instrument model", 1000483)
            }
            InstrumentModelType::Accela_PDA => param!("Accela PDA", 1000623),
            InstrumentModelType::Delta_Plus_Advantage => param!("DELTA plusAdvantage", 1000153),
            InstrumentModelType::Delta_Plus_XP => param!("DELTAplusXP", 1000154),
            InstrumentModelType::DFS => param!("DFS", 1000640),
            InstrumentModelType::DSQ => param!("DSQ", 1000634),
            InstrumentModelType::DSQ_II => param!("DSQ II", 1000641),
            InstrumentModelType::Element_2 => param!("Element 2", 1000646),
            InstrumentModelType::Element_GD => param!("Element GD", 1000647),
            InstrumentModelType::Element_XR => param!("Element XR", 1000645),
            InstrumentModelType::Exactive => param!("Exactive", 1000649),
            InstrumentModelType::Exactive_Plus => param!("Exactive Plus", 1002526),
            InstrumentModelType::GC_IsoLink => param!("GC IsoLink", 1000648),
            InstrumentModelType::GC_Quantum => param!("GC Quantum", 1000558),
            InstrumentModelType::ISQ => param!("ISQ", 1001908),
            InstrumentModelType::ITQ_1100 => param!("ITQ 1100", 1000637),
            InstrumentModelType::ITQ_700 => param!("ITQ 700", 1000635),
            InstrumentModelType::ITQ_900 => param!("ITQ 900", 1000636),
            InstrumentModelType::LCQ_Advantage => param!("LCQ Advantage", 1000167),
            InstrumentModelType::LCQ_Classic => param!("LCQ Classic", 1000168),
            InstrumentModelType::LCQ_Deca => param!("LCQ Deca", 1000554),
            InstrumentModelType::LCQ_Deca_XP_Plus => param!("LCQ Deca XP Plus", 1000169),
            InstrumentModelType::LCQ_Fleet => param!("LCQ Fleet", 1000578),
            InstrumentModelType::LTQ => param!("LTQ", 1000447),
            InstrumentModelType::LTQ_FT => param!("LTQ FT", 1000448),
            InstrumentModelType::LTQ_FT_Ultra => param!("LTQ FT Ultra", 1000557),
            InstrumentModelType::LTQ_Orbitrap => param!("LTQ Orbitrap", 1000449),
            InstrumentModelType::LTQ_Orbitrap_Classic => param!("LTQ Orbitrap Classic", 1002835),
            InstrumentModelType::LTQ_Orbitrap_Discovery => {
                param!("LTQ Orbitrap Discovery", 1000555)
            }
            InstrumentModelType::LTQ_Orbitrap_Elite => param!("LTQ Orbitrap Elite", 1001910),
            InstrumentModelType::LTQ_Orbitrap_Velos => param!("LTQ Orbitrap Velos", 1001742),
            InstrumentModelType::LTQ_Orbitrap_Velos_Pro => {
                param!("LTQ Orbitrap Velos Pro", 1003096)
            }
            InstrumentModelType::LTQ_Orbitrap_XL => param!("LTQ Orbitrap XL", 1000556),
            InstrumentModelType::LTQ_Orbitrap_XL_ETD => param!("LTQ Orbitrap XL ETD", 1000639),
            InstrumentModelType::LTQ_Velos => param!("LTQ Velos", 1000855),
            InstrumentModelType::LTQ_Velos_ETD => param!("LTQ Velos ETD", 1000856),
            InstrumentModelType::LTQ_XL => param!("LTQ XL", 1000854),
            InstrumentModelType::LTQ_XL_ETD => param!("LTQ XL ETD", 1000638),
            InstrumentModelType::LXQ => param!("LXQ", 1000450),
            InstrumentModelType::MALDI_LTQ_Orbitrap => param!("MALDI LTQ Orbitrap", 1000643),
            InstrumentModelType::MALDI_LTQ_XL => param!("MALDI LTQ XL", 1000642),
            InstrumentModelType::MAT253 => param!("MAT253", 1000172),
            InstrumentModelType::MAT900XP => param!("MAT900XP", 1000173),
            InstrumentModelType::MAT900XP_Trap => param!("MAT900XP Trap", 1000174),
            InstrumentModelType::MAT95XP => param!("MAT95XP", 1000175),
            InstrumentModelType::MAT95XP_Trap => param!("MAT95XP Trap", 1000176),
            InstrumentModelType::Orbitrap_Ascend => param!("Orbitrap Ascend", 1003356),
            InstrumentModelType::Orbitrap_Astral => param!("Orbitrap Astral", 1003378),
            InstrumentModelType::Orbitrap_Eclipse => param!("Orbitrap Eclipse", 1003029),
            InstrumentModelType::Orbitrap_Exploris_120 => param!("Orbitrap Exploris 120", 1003095),
            InstrumentModelType::Orbitrap_Exploris_240 => param!("Orbitrap Exploris 240", 1003094),
            InstrumentModelType::Orbitrap_Exploris_480 => param!("Orbitrap Exploris 480", 1003028),
            InstrumentModelType::Orbitrap_Exploris_GC_240 => param!("Orbitrap Exploris GC 240", 1003423),
            InstrumentModelType::Orbitrap_Fusion => param!("Orbitrap Fusion", 1002416),
            InstrumentModelType::Orbitrap_Fusion_ETD => param!("Orbitrap Fusion ETD", 1002417),
            InstrumentModelType::Orbitrap_Fusion_Lumos => param!("Orbitrap Fusion Lumos", 1002732),
            InstrumentModelType::Orbitrap_ID_X => param!("Orbitrap ID-X", 1003112),
            InstrumentModelType::Orbitrap_GC => param!("Q Exactive GC Orbitrap", 1003395),
            InstrumentModelType::PolarisQ => param!("PolarisQ", 1000185),
            InstrumentModelType::Q_Exactive => param!("Q Exactive", 1001911),
            InstrumentModelType::Q_Exactive_HF => param!("Q Exactive HF", 1002523),
            InstrumentModelType::Q_Exactive_HF_X => param!("Q Exactive HF-X", 1002877),
            InstrumentModelType::Q_Exactive_Plus => param!("Q Exactive Plus", 1002634),
            InstrumentModelType::Q_Exactive_UHMR => param!("Q Exactive UHMR", 1003245),
            InstrumentModelType::SSQ_7000 => param!("SSQ 7000", 1000748),
            InstrumentModelType::Surveyor_MSQ => param!("Surveyor MSQ", 1000193),
            InstrumentModelType::Surveyor_PDA => param!("Surveyor PDA", 1000622),
            InstrumentModelType::Tempus_TOF => param!("TEMPUS TOF", 1000196),
            InstrumentModelType::Trace_DSQ => param!("TRACE DSQ", 1000197),
            InstrumentModelType::Triton => param!("TRITON", 1000198),
            InstrumentModelType::TSQ => param!("TSQ", 1000750),
            InstrumentModelType::TSQ_7000 => param!("TSQ 7000", 1000749),
            InstrumentModelType::TSQ_8000_Evo => param!("TSQ 8000 Evo", 1002525),
            InstrumentModelType::TSQ_9000 => param!("TSQ 9000", 1002876),
            InstrumentModelType::TSQ_Altis => param!("TSQ Altis", 1002874),
            InstrumentModelType::TSQ_Altis_Plus => param!("TSQ Altis Plus", 1003292),
            InstrumentModelType::TSQ_Endura => param!("TSQ Endura", 1002419),
            InstrumentModelType::TSQ_Quantis => param!("TSQ Quantis", 1002875),
            InstrumentModelType::TSQ_Quantiva => param!("TSQ Quantiva", 1002418),
            InstrumentModelType::TSQ_Quantum => param!("TSQ Quantum", 1000199),
            InstrumentModelType::TSQ_Quantum_Access => param!("TSQ Quantum Access", 1000644),
            InstrumentModelType::TSQ_Quantum_Ultra => param!("TSQ Quantum Ultra", 1000751),
            InstrumentModelType::TSQ_Quantum_Ultra_AM => param!("TSQ Quantum Ultra AM", 1000743),
            InstrumentModelType::TSQ_Vantage => param!("TSQ Vantage", 1001510),
            InstrumentModelType::LTQ_Velos_Plus => param!("Velos Plus", 1001909),
            InstrumentModelType::Neptune => param!("neptune", 1000179),
        };
        val.into()
    }
}

static INSTRUMENT_MODEL_TYPE_MATCH: [(&str, InstrumentModelType, MatchType); 89] = [
    (
        "MAT253",
        InstrumentModelType::MAT253,
        MatchType::ExactNoSpaces,
    ),
    (
        "MAT900XP",
        InstrumentModelType::MAT900XP,
        MatchType::ExactNoSpaces,
    ),
    (
        "MAT900XPTRAP",
        InstrumentModelType::MAT900XP_Trap,
        MatchType::ExactNoSpaces,
    ),
    (
        "MAT95XP",
        InstrumentModelType::MAT95XP,
        MatchType::ExactNoSpaces,
    ),
    (
        "MAT95XPTRAP",
        InstrumentModelType::MAT95XP_Trap,
        MatchType::ExactNoSpaces,
    ),
    (
        "SSQ7000",
        InstrumentModelType::SSQ_7000,
        MatchType::ExactNoSpaces,
    ),
    (
        "TSQ7000",
        InstrumentModelType::TSQ_7000,
        MatchType::ExactNoSpaces,
    ),
    (
        "TSQ8000EVO",
        InstrumentModelType::TSQ_8000_Evo,
        MatchType::ExactNoSpaces,
    ),
    (
        "TSQ9000",
        InstrumentModelType::TSQ_9000,
        MatchType::ExactNoSpaces,
    ),
    ("TSQ", InstrumentModelType::TSQ, MatchType::Exact),
    (
        "ELEMENT2",
        InstrumentModelType::Element_2,
        MatchType::ExactNoSpaces,
    ),
    (
        "DELTA PLUSADVANTAGE",
        InstrumentModelType::Delta_Plus_Advantage,
        MatchType::Exact,
    ),
    (
        "DELTAPLUSXP",
        InstrumentModelType::Delta_Plus_XP,
        MatchType::Exact,
    ),
    (
        "LCQ ADVANTAGE",
        InstrumentModelType::LCQ_Advantage,
        MatchType::Exact,
    ),
    (
        "LCQ CLASSIC",
        InstrumentModelType::LCQ_Classic,
        MatchType::Exact,
    ),
    ("LCQ DECA", InstrumentModelType::LCQ_Deca, MatchType::Exact),
    (
        "LCQ DECA XP",
        InstrumentModelType::LCQ_Deca_XP_Plus,
        MatchType::Exact,
    ),
    (
        "LCQ DECA XP PLUS",
        InstrumentModelType::LCQ_Deca_XP_Plus,
        MatchType::Exact,
    ),
    ("NEPTUNE", InstrumentModelType::Neptune, MatchType::Exact),
    ("DSQ", InstrumentModelType::DSQ, MatchType::Exact),
    ("POLARISQ", InstrumentModelType::PolarisQ, MatchType::Exact),
    (
        "SURVEYOR MSQ",
        InstrumentModelType::Surveyor_MSQ,
        MatchType::Exact,
    ),
    (
        "MSQ PLUS",
        InstrumentModelType::Surveyor_MSQ,
        MatchType::Exact,
    ),
    (
        "TEMPUS TOF",
        InstrumentModelType::Tempus_TOF,
        MatchType::Exact,
    ),
    (
        "TRACE DSQ",
        InstrumentModelType::Trace_DSQ,
        MatchType::Exact,
    ),
    ("TRITON", InstrumentModelType::Triton, MatchType::Exact),
    ("LTQ", InstrumentModelType::LTQ, MatchType::Exact),
    ("LTQ XL", InstrumentModelType::LTQ_XL, MatchType::Exact),
    ("LTQ FT", InstrumentModelType::LTQ_FT, MatchType::Exact),
    ("LTQ-FT", InstrumentModelType::LTQ_FT, MatchType::Exact),
    (
        "LTQ FT ULTRA",
        InstrumentModelType::LTQ_FT_Ultra,
        MatchType::Exact,
    ),
    (
        "LTQ ORBITRAP",
        InstrumentModelType::LTQ_Orbitrap,
        MatchType::Exact,
    ),
    (
        "LTQ ORBITRAP CLASSIC",
        InstrumentModelType::LTQ_Orbitrap_Classic,
        MatchType::Exact,
    ),
    (
        "LTQ ORBITRAP DISCOVERY",
        InstrumentModelType::LTQ_Orbitrap_Discovery,
        MatchType::Exact,
    ),
    (
        "LTQ ORBITRAP XL",
        InstrumentModelType::LTQ_Orbitrap_XL,
        MatchType::Exact,
    ),
    (
        "ORBITRAP VELOS PRO",
        InstrumentModelType::LTQ_Orbitrap_Velos_Pro,
        MatchType::Contains,
    ),
    (
        "ORBITRAP VELOS",
        InstrumentModelType::LTQ_Orbitrap_Velos,
        MatchType::Contains,
    ),
    (
        "ORBITRAP ELITE",
        InstrumentModelType::LTQ_Orbitrap_Elite,
        MatchType::Contains,
    ),
    (
        "VELOS PLUS",
        InstrumentModelType::LTQ_Velos_Plus,
        MatchType::Contains,
    ),
    (
        "VELOS PRO",
        InstrumentModelType::LTQ_Velos_Plus,
        MatchType::Contains,
    ),
    (
        "LTQ VELOS",
        InstrumentModelType::LTQ_Velos,
        MatchType::Exact,
    ),
    (
        "LTQ VELOS ETD",
        InstrumentModelType::LTQ_Velos_ETD,
        MatchType::Exact,
    ),
    ("LXQ", InstrumentModelType::LXQ, MatchType::Exact),
    (
        "LCQ FLEET",
        InstrumentModelType::LCQ_Fleet,
        MatchType::Exact,
    ),
    ("ITQ 700", InstrumentModelType::ITQ_700, MatchType::Exact),
    ("ITQ 900", InstrumentModelType::ITQ_900, MatchType::Exact),
    ("ITQ 1100", InstrumentModelType::ITQ_1100, MatchType::Exact),
    (
        "GC QUANTUM",
        InstrumentModelType::GC_Quantum,
        MatchType::Exact,
    ),
    (
        "LTQ XL ETD",
        InstrumentModelType::LTQ_XL_ETD,
        MatchType::Exact,
    ),
    (
        "LTQ ORBITRAP XL ETD",
        InstrumentModelType::LTQ_Orbitrap_XL_ETD,
        MatchType::Exact,
    ),
    ("DFS", InstrumentModelType::DFS, MatchType::Exact),
    ("DSQ II", InstrumentModelType::DSQ_II, MatchType::Exact),
    ("ISQ SERIES", InstrumentModelType::ISQ, MatchType::Exact),
    (
        "MALDI LTQ XL",
        InstrumentModelType::MALDI_LTQ_XL,
        MatchType::Exact,
    ),
    (
        "MALDI LTQ ORBITRAP",
        InstrumentModelType::MALDI_LTQ_Orbitrap,
        MatchType::Exact,
    ),
    (
        "TSQ QUANTUM",
        InstrumentModelType::TSQ_Quantum,
        MatchType::Exact,
    ),
    (
        "TSQ QUANTUM ACCESS",
        InstrumentModelType::TSQ_Quantum_Access,
        MatchType::Contains,
    ),
    (
        "TSQ QUANTUM ULTRA",
        InstrumentModelType::TSQ_Quantum_Ultra,
        MatchType::Exact,
    ),
    (
        "TSQ QUANTUM ULTRA AM",
        InstrumentModelType::TSQ_Quantum_Ultra_AM,
        MatchType::Exact,
    ),
    (
        "TSQ VANTAGE",
        InstrumentModelType::TSQ_Vantage,
        MatchType::StartsWith,
    ),
    (
        "TSQ QUANTIVA",
        InstrumentModelType::TSQ_Quantiva,
        MatchType::Exact,
    ),
    (
        "TSQ ENDURA",
        InstrumentModelType::TSQ_Endura,
        MatchType::Exact,
    ),
    (
        "TSQ ALTIS",
        InstrumentModelType::TSQ_Altis,
        MatchType::Exact,
    ),
    (
        "TSQ ALTIS PLUS",
        InstrumentModelType::TSQ_Altis_Plus,
        MatchType::Exact,
    ),
    (
        "TSQ QUANTIS",
        InstrumentModelType::TSQ_Quantis,
        MatchType::Exact,
    ),
    (
        "ELEMENT XR",
        InstrumentModelType::Element_XR,
        MatchType::Exact,
    ),
    (
        "ELEMENT GD",
        InstrumentModelType::Element_GD,
        MatchType::Exact,
    ),
    (
        "GC ISOLINK",
        InstrumentModelType::GC_IsoLink,
        MatchType::Exact,
    ),
    (
        "ORBITRAP ID-X",
        InstrumentModelType::Orbitrap_ID_X,
        MatchType::Exact,
    ),
    (
        "Q EXACTIVE PLUS",
        InstrumentModelType::Q_Exactive_Plus,
        MatchType::Contains,
    ),
    (
        "Q EXACTIVE HF-X",
        InstrumentModelType::Q_Exactive_HF_X,
        MatchType::Contains,
    ),
    (
        "Q EXACTIVE HF",
        InstrumentModelType::Q_Exactive_HF,
        MatchType::Contains,
    ),
    (
        "Q EXACTIVE UHMR",
        InstrumentModelType::Q_Exactive_UHMR,
        MatchType::Contains,
    ),
    (
        "Q EXACTIVE",
        InstrumentModelType::Q_Exactive,
        MatchType::Contains,
    ),
    (
        "EXACTIVE PLUS",
        InstrumentModelType::Exactive_Plus,
        MatchType::Contains,
    ),
    (
        "EXACTIVE",
        InstrumentModelType::Exactive,
        MatchType::Contains,
    ),
    (
        "ORBITRAP EXPLORIS 120",
        InstrumentModelType::Orbitrap_Exploris_120,
        MatchType::Exact,
    ),
    (
        "ORBITRAP EXPLORIS 240",
        InstrumentModelType::Orbitrap_Exploris_240,
        MatchType::Exact,
    ),
    (
        "ORBITRAP EXPLORIS 480",
        InstrumentModelType::Orbitrap_Exploris_480,
        MatchType::Exact,
    ),
    (
        "ORBITRAP EXPLORIS GC 240",
        InstrumentModelType::Orbitrap_Exploris_GC_240,
        MatchType::Exact,
    ),
    (
        "ORBITRAP GC",
        InstrumentModelType::Orbitrap_GC,
        MatchType::Contains,
    ),
    (
        "ECLIPSE",
        InstrumentModelType::Orbitrap_Eclipse,
        MatchType::Contains,
    ),
    (
        "ASTRAL",
        InstrumentModelType::Orbitrap_Astral,
        MatchType::Contains,
    ),
    (
        "FUSION ETD",
        InstrumentModelType::Orbitrap_Fusion_ETD,
        MatchType::Contains,
    ),
    (
        "FUSION LUMOS",
        InstrumentModelType::Orbitrap_Fusion_Lumos,
        MatchType::Contains,
    ),
    (
        "FUSION",
        InstrumentModelType::Orbitrap_Fusion,
        MatchType::Contains,
    ),
    (
        "ASCEND",
        InstrumentModelType::Orbitrap_Ascend,
        MatchType::Contains,
    ),
    (
        "SURVEYOR PDA",
        InstrumentModelType::Surveyor_PDA,
        MatchType::Exact,
    ),
    (
        "ACCELA PDA",
        InstrumentModelType::Accela_PDA,
        MatchType::Exact,
    ),
];

pub fn parse_instrument_model(instrument_model: &str) -> InstrumentModelType {
    let model_type = instrument_model.to_uppercase();
    let model_type_no_spaces = model_type.replace(" ", "");
    log::debug!("Parsing instrument model: '{}' -> '{}' (no spaces: '{}')",
        instrument_model, model_type, model_type_no_spaces);
    
    for (key, model_enum, match_type) in INSTRUMENT_MODEL_TYPE_MATCH.iter() {
        let hit = match match_type {
            MatchType::Exact => **key == model_type,
            MatchType::Contains => model_type.contains(key),
            MatchType::StartsWith => model_type.starts_with(key),
            MatchType::EndsWith => model_type.ends_with(key),
            MatchType::ExactNoSpaces => **key == model_type_no_spaces,
        };
        if hit {
            log::debug!("Matched instrument model '{}' with pattern '{}' ({:?}) -> {:?}",
                instrument_model, key, match_type, model_enum);
            return *model_enum;
        }
    }
    log::warn!("Failed to infer instrument model from name string '{}' - will use Unknown", instrument_model);
    InstrumentModelType::Unknown
}


pub fn instrument_model_to_mass_analyzers(model: InstrumentModelType) -> Vec<MassAnalyzerTerm> {
    match model {
        InstrumentModelType::Exactive
        | InstrumentModelType::Exactive_Plus
        | InstrumentModelType::Q_Exactive
        | InstrumentModelType::Q_Exactive_Plus
        | InstrumentModelType::Q_Exactive_HF_X
        | InstrumentModelType::Q_Exactive_HF
        | InstrumentModelType::Q_Exactive_UHMR
        | InstrumentModelType::Orbitrap_Exploris_120
        | InstrumentModelType::Orbitrap_Exploris_240
        | InstrumentModelType::Orbitrap_Exploris_480
        | InstrumentModelType::Orbitrap_Exploris_GC_240
        | InstrumentModelType::Orbitrap_GC => {
            vec![MassAnalyzerTerm::Orbitrap]
        }
        InstrumentModelType::LTQ_Orbitrap
        | InstrumentModelType::LTQ_Orbitrap_Classic
        | InstrumentModelType::LTQ_Orbitrap_Discovery
        | InstrumentModelType::LTQ_Orbitrap_XL
        | InstrumentModelType::MALDI_LTQ_Orbitrap
        | InstrumentModelType::LTQ_Orbitrap_Velos
        | InstrumentModelType::LTQ_Orbitrap_Velos_Pro
        | InstrumentModelType::LTQ_Orbitrap_Elite
        | InstrumentModelType::Orbitrap_Fusion
        | InstrumentModelType::Orbitrap_Fusion_Lumos
        | InstrumentModelType::Orbitrap_Fusion_ETD
        | InstrumentModelType::Orbitrap_Ascend
        | InstrumentModelType::Orbitrap_ID_X
        | InstrumentModelType::Orbitrap_Eclipse => {
            vec![MassAnalyzerTerm::Orbitrap, MassAnalyzerTerm::LinearIonTrap]
        }
        InstrumentModelType::Orbitrap_Astral => {
            vec![
                MassAnalyzerTerm::Orbitrap,
                MassAnalyzerTerm::AsymmetricTrackLosslessTimeOfFlightAnalyzer,
            ]
        }
        InstrumentModelType::LTQ_FT | InstrumentModelType::LTQ_FT_Ultra => {
            vec![
                MassAnalyzerTerm::FourierTransformIonCyclotronResonanceMassSpectrometer,
                MassAnalyzerTerm::LinearIonTrap,
            ]
        }
        InstrumentModelType::SSQ_7000
        | InstrumentModelType::Surveyor_MSQ
        | InstrumentModelType::DSQ
        | InstrumentModelType::DSQ_II
        | InstrumentModelType::ISQ
        | InstrumentModelType::Trace_DSQ
        | InstrumentModelType::GC_IsoLink => {
            vec![MassAnalyzerTerm::Quadrupole]
        }
        InstrumentModelType::TSQ_7000
        | InstrumentModelType::TSQ_8000_Evo
        | InstrumentModelType::TSQ_9000
        | InstrumentModelType::TSQ
        | InstrumentModelType::TSQ_Quantum
        | InstrumentModelType::TSQ_Quantum_Access
        | InstrumentModelType::TSQ_Quantum_Ultra
        | InstrumentModelType::TSQ_Quantum_Ultra_AM
        | InstrumentModelType::GC_Quantum
        | InstrumentModelType::TSQ_Quantiva
        | InstrumentModelType::TSQ_Endura
        | InstrumentModelType::TSQ_Altis
        | InstrumentModelType::TSQ_Altis_Plus
        | InstrumentModelType::TSQ_Quantis => {
            vec![MassAnalyzerTerm::Quadrupole]
        }

        InstrumentModelType::LCQ_Advantage
        | InstrumentModelType::LCQ_Classic
        | InstrumentModelType::LCQ_Deca
        | InstrumentModelType::LCQ_Deca_XP_Plus
        | InstrumentModelType::LCQ_Fleet
        | InstrumentModelType::PolarisQ
        | InstrumentModelType::ITQ_700
        | InstrumentModelType::ITQ_900 => {
            vec![MassAnalyzerTerm::QuadrupoleIonTrap]
        }

        InstrumentModelType::LTQ
        | InstrumentModelType::LXQ
        | InstrumentModelType::LTQ_XL
        | InstrumentModelType::LTQ_XL_ETD
        | InstrumentModelType::LTQ_Orbitrap_XL_ETD
        | InstrumentModelType::ITQ_1100
        | InstrumentModelType::MALDI_LTQ_XL
        | InstrumentModelType::LTQ_Velos
        | InstrumentModelType::LTQ_Velos_ETD
        | InstrumentModelType::LTQ_Velos_Plus => {
            vec![MassAnalyzerTerm::LinearIonTrap]
        }

        InstrumentModelType::DFS
        | InstrumentModelType::MAT253
        | InstrumentModelType::MAT900XP
        | InstrumentModelType::MAT900XP_Trap
        | InstrumentModelType::MAT95XP
        | InstrumentModelType::MAT95XP_Trap => {
            vec![MassAnalyzerTerm::MagneticSector]
        }

        InstrumentModelType::Tempus_TOF => vec![MassAnalyzerTerm::TimeOfFlight],
        _ => Vec::default(),
    }
}

pub fn instrument_model_to_detector(model: InstrumentModelType) -> Vec<DetectorTypeTerm> {
    match model {
        InstrumentModelType::Q_Exactive |
        InstrumentModelType::Q_Exactive_Plus |
        InstrumentModelType::Q_Exactive_HF |
        InstrumentModelType::Q_Exactive_HF_X |
        InstrumentModelType::Q_Exactive_UHMR |
        InstrumentModelType::Orbitrap_Exploris_120 |
        InstrumentModelType::Orbitrap_Exploris_240 |
        InstrumentModelType::Orbitrap_Exploris_480 |
        InstrumentModelType::Orbitrap_Exploris_GC_240 |
        InstrumentModelType::Orbitrap_GC => {
            vec![DetectorTypeTerm::InductiveDetector]
        },

        InstrumentModelType::Exactive |
        InstrumentModelType::Exactive_Plus => {
            vec![DetectorTypeTerm::InductiveDetector]
        },

        InstrumentModelType::Orbitrap_Astral => {
            vec![DetectorTypeTerm::InductiveDetector, DetectorTypeTerm::ElectronMultiplier]
        },

        InstrumentModelType::LTQ_FT |
        InstrumentModelType::LTQ_FT_Ultra => {
            vec![DetectorTypeTerm::InductiveDetector, DetectorTypeTerm::ElectronMultiplier]
        }

        InstrumentModelType::Orbitrap_Fusion |
        InstrumentModelType::Orbitrap_Fusion_Lumos |
        InstrumentModelType::Orbitrap_Fusion_ETD |
        InstrumentModelType::Orbitrap_Ascend |
        InstrumentModelType::Orbitrap_ID_X |
        InstrumentModelType::Orbitrap_Eclipse => {
            vec![DetectorTypeTerm::InductiveDetector, DetectorTypeTerm::ElectronMultiplier]
        }


        InstrumentModelType::LTQ_Orbitrap |
        InstrumentModelType::LTQ_Orbitrap_Classic |
        InstrumentModelType::LTQ_Orbitrap_Discovery |
        InstrumentModelType::LTQ_Orbitrap_XL |
        InstrumentModelType::LTQ_Orbitrap_XL_ETD |
        InstrumentModelType::MALDI_LTQ_Orbitrap |
        InstrumentModelType::LTQ_Orbitrap_Velos |
        InstrumentModelType::LTQ_Orbitrap_Velos_Pro |
        InstrumentModelType::LTQ_Orbitrap_Elite => {
            vec![DetectorTypeTerm::InductiveDetector, DetectorTypeTerm::ElectronMultiplier]
        }

        InstrumentModelType::LCQ_Advantage |
        InstrumentModelType::LCQ_Classic |
        InstrumentModelType::LCQ_Deca |
        InstrumentModelType::LCQ_Deca_XP_Plus |
        InstrumentModelType::LCQ_Fleet |
        InstrumentModelType::PolarisQ |
        InstrumentModelType::ITQ_700 |
        InstrumentModelType::ITQ_900 => {
            vec![DetectorTypeTerm::ElectronMultiplier]
        }

        InstrumentModelType::LTQ |
        InstrumentModelType::LXQ |
        InstrumentModelType::LTQ_XL |
        InstrumentModelType::LTQ_XL_ETD |
        InstrumentModelType::ITQ_1100 |
        InstrumentModelType::MALDI_LTQ_XL |
        InstrumentModelType::LTQ_Velos |
        InstrumentModelType::LTQ_Velos_ETD |
        InstrumentModelType::LTQ_Velos_Plus => {
            vec![DetectorTypeTerm::ElectronMultiplier]
        }

        InstrumentModelType::SSQ_7000 |
        InstrumentModelType::Surveyor_MSQ |
        InstrumentModelType::DSQ |
        InstrumentModelType::DSQ_II |
        InstrumentModelType::ISQ |
        InstrumentModelType::Trace_DSQ |
        InstrumentModelType::GC_IsoLink => {
            vec![DetectorTypeTerm::ElectronMultiplier]
        },

        InstrumentModelType::TSQ_7000 |
        InstrumentModelType::TSQ_8000_Evo |
        InstrumentModelType::TSQ_9000 |
        InstrumentModelType::TSQ |
        InstrumentModelType::TSQ_Quantum |
        InstrumentModelType::TSQ_Quantum_Access |
        InstrumentModelType::TSQ_Quantum_Ultra |
        InstrumentModelType::TSQ_Quantum_Ultra_AM |
        InstrumentModelType::TSQ_Vantage |
        // InstrumentModelType::TSQ_Vantage_EMR |
        // InstrumentModelType::TSQ_Vantage_AM |
        InstrumentModelType::GC_Quantum |
        InstrumentModelType::TSQ_Quantiva |
        InstrumentModelType::TSQ_Endura |
        InstrumentModelType::TSQ_Altis |
        InstrumentModelType::TSQ_Altis_Plus |
        InstrumentModelType::TSQ_Quantis => {
            vec![DetectorTypeTerm::ElectronMultiplier]
        }

        InstrumentModelType::DFS |
        InstrumentModelType::MAT253 |
        InstrumentModelType::MAT900XP |
        InstrumentModelType::MAT900XP_Trap |
        InstrumentModelType::MAT95XP |
        InstrumentModelType::MAT95XP_Trap => {
            vec![DetectorTypeTerm::ElectronMultiplier]
        }

        InstrumentModelType::Tempus_TOF |
        InstrumentModelType::Element_2 |
        InstrumentModelType::Element_XR |
        InstrumentModelType::Element_GD |
        InstrumentModelType::Delta_Plus_Advantage |
        InstrumentModelType::Delta_Plus_XP |
        InstrumentModelType::Neptune |
        InstrumentModelType::Triton => {
            vec![]
        }

        InstrumentModelType::Surveyor_PDA |
        InstrumentModelType::Accela_PDA => {
            vec![]
        }

        InstrumentModelType::Unknown => {
            vec![]
        },
    }
}

pub fn instrument_model_to_ion_sources(model: InstrumentModelType) -> Vec<IonizationTypeTerm> {
    match model {
        InstrumentModelType::SSQ_7000
        | InstrumentModelType::TSQ_7000
        | InstrumentModelType::TSQ_8000_Evo
        | InstrumentModelType::TSQ_9000
        | InstrumentModelType::Surveyor_MSQ
        | InstrumentModelType::LCQ_Advantage
        | InstrumentModelType::LCQ_Classic
        | InstrumentModelType::LCQ_Deca
        | InstrumentModelType::LCQ_Deca_XP_Plus
        | InstrumentModelType::LCQ_Fleet
        | InstrumentModelType::LXQ
        | InstrumentModelType::LTQ
        | InstrumentModelType::LTQ_XL
        | InstrumentModelType::LTQ_XL_ETD
        | InstrumentModelType::LTQ_Velos
        | InstrumentModelType::LTQ_Velos_ETD
        | InstrumentModelType::LTQ_Velos_Plus
        | InstrumentModelType::LTQ_FT
        | InstrumentModelType::LTQ_FT_Ultra
        | InstrumentModelType::LTQ_Orbitrap
        | InstrumentModelType::LTQ_Orbitrap_Classic
        | InstrumentModelType::LTQ_Orbitrap_Discovery
        | InstrumentModelType::LTQ_Orbitrap_XL
        | InstrumentModelType::LTQ_Orbitrap_XL_ETD
        | InstrumentModelType::LTQ_Orbitrap_Velos
        | InstrumentModelType::LTQ_Orbitrap_Velos_Pro
        | InstrumentModelType::LTQ_Orbitrap_Elite
        | InstrumentModelType::Exactive
        | InstrumentModelType::Exactive_Plus
        | InstrumentModelType::Q_Exactive
        | InstrumentModelType::Q_Exactive_Plus
        | InstrumentModelType::Q_Exactive_HF
        | InstrumentModelType::Q_Exactive_HF_X
        | InstrumentModelType::Q_Exactive_UHMR
        | InstrumentModelType::Orbitrap_Exploris_120
        | InstrumentModelType::Orbitrap_Exploris_240
        | InstrumentModelType::Orbitrap_Exploris_480
        | InstrumentModelType::Orbitrap_Eclipse
        | InstrumentModelType::Orbitrap_Fusion
        | InstrumentModelType::Orbitrap_Fusion_Lumos
        | InstrumentModelType::Orbitrap_Fusion_ETD
        | InstrumentModelType::Orbitrap_Ascend
        | InstrumentModelType::Orbitrap_ID_X
        | InstrumentModelType::Orbitrap_Astral
        | InstrumentModelType::TSQ
        | InstrumentModelType::TSQ_Quantum
        | InstrumentModelType::TSQ_Quantum_Access
        | InstrumentModelType::TSQ_Quantum_Ultra
        | InstrumentModelType::TSQ_Quantum_Ultra_AM
        | InstrumentModelType::TSQ_Quantiva
        | InstrumentModelType::TSQ_Endura
        | InstrumentModelType::TSQ_Altis
        | InstrumentModelType::TSQ_Altis_Plus
        | InstrumentModelType::TSQ_Quantis => {
            vec![IonizationTypeTerm::ElectrosprayIonization]
        }

        InstrumentModelType::DSQ
        | InstrumentModelType::PolarisQ
        | InstrumentModelType::ITQ_700
        | InstrumentModelType::ITQ_900
        | InstrumentModelType::ITQ_1100
        | InstrumentModelType::Trace_DSQ
        | InstrumentModelType::GC_Quantum
        | InstrumentModelType::DFS
        | InstrumentModelType::DSQ_II
        | InstrumentModelType::ISQ
        | InstrumentModelType::GC_IsoLink
        | InstrumentModelType::Orbitrap_Exploris_GC_240
        | InstrumentModelType::Orbitrap_GC => {
            vec![IonizationTypeTerm::ElectronIonization]
        }
        InstrumentModelType::MALDI_LTQ_XL | InstrumentModelType::MALDI_LTQ_Orbitrap => {
            vec![IonizationTypeTerm::MatrixAssistedLaserDesorptionIonization]
        }
        _ => Vec::default(),
    }
}


macro_rules! comp {
    ($config:ident, $comptype:expr, $term:expr) => {
        $config.new_component($comptype).add_param($term.into())
    };
}

macro_rules! analyzer {
    ($config:ident, $term:expr) => {
        comp!($config, ComponentType::Analyzer, $term)
    };
    ($config:ident quadrupole) => {
        analyzer!($config, MassAnalyzerTerm::Quadrupole)
    };
    ($config:ident orbitrap) => {
        analyzer!($config, MassAnalyzerTerm::Orbitrap)
    };
    ($config:ident radial) => {
        analyzer!($config, MassAnalyzerTerm::RadialEjectionLinearIonTrap)
    }
}

macro_rules! detector {
    ($config:ident, $term:expr) => {
        comp!($config, ComponentType::Detector, $term)
    };
    ($config:ident inductive) => {
        detector!($config, DetectorTypeTerm::InductiveDetector)
    };
    ($config:ident electron) => {
        detector!($config, DetectorTypeTerm::ElectronMultiplier)
    }
}


pub fn create_instrument_configurations(model: InstrumentModelType, source: Component) -> Vec<InstrumentConfiguration> {
    let mut configs = Vec::new();

    match model {
        InstrumentModelType::Q_Exactive |
        InstrumentModelType::Q_Exactive_Plus |
        InstrumentModelType::Q_Exactive_HF |
        InstrumentModelType::Q_Exactive_HF_X |
        InstrumentModelType::Q_Exactive_UHMR |
        InstrumentModelType::Orbitrap_Exploris_120 |
        InstrumentModelType::Orbitrap_Exploris_240 |
        InstrumentModelType::Orbitrap_Exploris_480 |
        InstrumentModelType::Orbitrap_Exploris_GC_240|
		InstrumentModelType::Orbitrap_GC  => {
            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());
            comp!(config, ComponentType::Analyzer, MassAnalyzerTerm::Quadrupole);
            comp!(config, ComponentType::Analyzer, MassAnalyzerTerm::Orbitrap);
            comp!(config, ComponentType::Detector, DetectorTypeTerm::InductiveDetector);
        },

        InstrumentModelType::Exactive |
        InstrumentModelType::Exactive_Plus => {
            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());
            comp!(config, ComponentType::Analyzer, MassAnalyzerTerm::Orbitrap);
            comp!(config, ComponentType::Detector, DetectorTypeTerm::InductiveDetector);

        }

        InstrumentModelType::Orbitrap_Astral => {
            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());
            comp!(config, ComponentType::Analyzer, MassAnalyzerTerm::Quadrupole);
            comp!(config, ComponentType::Analyzer, MassAnalyzerTerm::Orbitrap);
            comp!(config, ComponentType::Detector, DetectorTypeTerm::InductiveDetector);

            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());
            comp!(config, ComponentType::Analyzer, MassAnalyzerTerm::Quadrupole);
            comp!(config, ComponentType::Analyzer, MassAnalyzerTerm::AsymmetricTrackLosslessTimeOfFlightAnalyzer);
            comp!(config, ComponentType::Detector, DetectorTypeTerm::InductiveDetector);
        }
        InstrumentModelType::LTQ_FT |
        InstrumentModelType::LTQ_FT_Ultra => {
            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());
            analyzer!(config, MassAnalyzerTerm::FourierTransformIonCyclotronResonanceMassSpectrometer);
            detector!(config, DetectorTypeTerm::InductiveDetector);

            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());
            analyzer!(config, MassAnalyzerTerm::RadialEjectionLinearIonTrap);
            detector!(config, DetectorTypeTerm::ElectronMultiplier);
        },
        InstrumentModelType::Orbitrap_Fusion |
        InstrumentModelType::Orbitrap_Fusion_Lumos |
        InstrumentModelType::Orbitrap_Fusion_ETD |
        InstrumentModelType::Orbitrap_Ascend |
        InstrumentModelType::Orbitrap_ID_X |
        InstrumentModelType::Orbitrap_Eclipse => {
            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());

            analyzer!(config quadrupole);
            analyzer!(config orbitrap);
            detector!(config inductive);

            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());

            analyzer!(config quadrupole);
            analyzer!(config radial);
            detector!(config electron);

        }

        InstrumentModelType::LTQ_Orbitrap |
        InstrumentModelType::LTQ_Orbitrap_Classic |
        InstrumentModelType::LTQ_Orbitrap_Discovery |
        InstrumentModelType::LTQ_Orbitrap_XL |
        InstrumentModelType::LTQ_Orbitrap_XL_ETD |
        InstrumentModelType::MALDI_LTQ_Orbitrap |
        InstrumentModelType::LTQ_Orbitrap_Velos |
        InstrumentModelType::LTQ_Orbitrap_Velos_Pro |
        InstrumentModelType::LTQ_Orbitrap_Elite => {
            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());

            analyzer!(config orbitrap);
            detector!(config inductive);

            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());

            analyzer!(config radial);
            detector!(config electron);
        },
        InstrumentModelType::LCQ_Advantage |
        InstrumentModelType::LCQ_Classic |
        InstrumentModelType::LCQ_Deca |
        InstrumentModelType::LCQ_Deca_XP_Plus |
        InstrumentModelType::LCQ_Fleet |
        InstrumentModelType::PolarisQ |
        InstrumentModelType::ITQ_700 |
        InstrumentModelType::ITQ_900 => {
            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());

            analyzer!(config, MassAnalyzerTerm::QuadrupoleIonTrap);
            detector!(config electron);
        },


        InstrumentModelType::LTQ |
        InstrumentModelType::LXQ |
        InstrumentModelType::LTQ_XL |
        InstrumentModelType::LTQ_XL_ETD |
        InstrumentModelType::ITQ_1100 |
        InstrumentModelType::MALDI_LTQ_XL |
        InstrumentModelType::LTQ_Velos |
        InstrumentModelType::LTQ_Velos_ETD |
        InstrumentModelType::LTQ_Velos_Plus => {
            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());

            analyzer!(config radial);
            detector!(config electron);
        },
        InstrumentModelType::SSQ_7000 |
        InstrumentModelType::Surveyor_MSQ |
        InstrumentModelType::DSQ |
        InstrumentModelType::DSQ_II |
        InstrumentModelType::ISQ |
        InstrumentModelType::Trace_DSQ |
        InstrumentModelType::GC_IsoLink => {
            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());

            analyzer!(config quadrupole);
            detector!(config electron);

        },
        InstrumentModelType::TSQ_7000 |
        InstrumentModelType::TSQ_8000_Evo |
        InstrumentModelType::TSQ_9000 |
        InstrumentModelType::TSQ |
        InstrumentModelType::TSQ_Quantum |
        InstrumentModelType::TSQ_Quantum_Access |
        InstrumentModelType::TSQ_Quantum_Ultra |
        InstrumentModelType::TSQ_Quantum_Ultra_AM |
        InstrumentModelType::GC_Quantum |
        InstrumentModelType::TSQ_Quantiva |
        InstrumentModelType::TSQ_Endura |
        InstrumentModelType::TSQ_Altis |
        InstrumentModelType::TSQ_Altis_Plus |
        InstrumentModelType::TSQ_Quantis => {
            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());
            analyzer!(config quadrupole);
            analyzer!(config quadrupole);
            analyzer!(config quadrupole);
            detector!(config electron);
        },
        InstrumentModelType::DFS |
        InstrumentModelType::MAT253 |
        InstrumentModelType::MAT900XP |
        InstrumentModelType::MAT900XP_Trap |
        InstrumentModelType::MAT95XP |
        InstrumentModelType::MAT95XP_Trap => {
            configs.push(InstrumentConfiguration::default());
            let config = configs.last_mut().unwrap();
            config.push(source.clone());

            analyzer!(config, MassAnalyzerTerm::MagneticSector);
            detector!(config, DetectorTypeTerm::ElectronMultiplier);
        },
        _ => {}
    }

    configs
}