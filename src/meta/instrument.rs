use std::fmt::Display;

use crate::impl_param_described;
use crate::params::{ParamCow, ParamList};

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

impl Component {
    pub fn name(&self) -> Option<&str> {
        match self.component_type {
            ComponentType::Analyzer => self
                .params
                .iter()
                .flat_map(|p| {
                    if let Some(u) = p.accession {
                        Some(MassAnalyzerTerm::from_accession(u).unwrap().name())
                    } else {
                        None
                    }
                })
                .next(),
            ComponentType::IonSource => self
                .params
                .iter()
                .flat_map(|p| {
                    if let Some(u) = p.accession {
                        if let Some(t) = IonizationTypeTerm::from_accession(u) {
                            Some(t.name())
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                })
                .next(),
            ComponentType::Detector => self
                .params
                .iter()
                .flat_map(|p| {
                    if let Some(u) = p.accession {
                        Some(DetectorTypeTerm::from_accession(u).unwrap().name())
                    } else {
                        None
                    }
                })
                .next(),
            ComponentType::Unknown => None,
        }
    }

    pub fn parent_types(&self) -> Vec<ParamCow<'static>> {
        match self.component_type {
            ComponentType::Analyzer => self
                .params
                .iter()
                .flat_map(|p| {
                    if let Some(u) = p.accession {
                        Some(
                            MassAnalyzerTerm::from_accession(u)
                                .unwrap()
                                .parents()
                                .into_iter()
                                .map(|t| t.to_param())
                                .collect(),
                        )
                    } else {
                        None
                    }
                })
                .next()
                .unwrap_or_default(),
            ComponentType::IonSource => self
                .params
                .iter()
                .flat_map(|p| {
                    if let Some(u) = p.accession {
                        if let Some(t) = IonizationTypeTerm::from_accession(u) {
                            Some(t.parents().into_iter().map(|t| t.to_param()).collect())
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                })
                .next()
                .unwrap_or_default(),
            ComponentType::Detector => self
                .params
                .iter()
                .flat_map(|p| {
                    if let Some(u) = p.accession {
                        Some(
                            DetectorTypeTerm::from_accession(u)
                                .unwrap()
                                .parents()
                                .into_iter()
                                .map(|t| t.to_param())
                                .collect(),
                        )
                    } else {
                        None
                    }
                })
                .next()
                .unwrap_or_default(),
            ComponentType::Unknown => vec![],
        }
    }
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

crate::cvmap! {
    #[flag_type=i32]
    #[allow(unused)]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_component.py', "mass-analyzer"]).decode('utf8')
    for line in buf.splitlines():
        cog.outl(line)
    ]]]*/
    pub enum MassAnalyzerTerm {
        #[term(cv=MS, accession=1000079, name="fourier transform ion cyclotron resonance mass spectrometer", flags={0}, parents={["MS:1000443"]})]
        FourierTransformIonCyclotronResonanceMassSpectrometer,
        #[term(cv=MS, accession=1000080, name="magnetic sector", flags={0}, parents={["MS:1000443"]})]
        MagneticSector,
        #[term(cv=MS, accession=1000081, name="quadrupole", flags={0}, parents={["MS:1000443"]})]
        Quadrupole,
        #[term(cv=MS, accession=1000084, name="time-of-flight", flags={0}, parents={["MS:1000443"]})]
        TimeOfFlight,
        #[term(cv=MS, accession=1000254, name="electrostatic energy analyzer", flags={0}, parents={["MS:1000443"]})]
        ElectrostaticEnergyAnalyzer,
        #[term(cv=MS, accession=1000264, name="ion trap", flags={0}, parents={["MS:1000443"]})]
        IonTrap,
        #[term(cv=MS, accession=1000284, name="stored waveform inverse fourier transform", flags={0}, parents={["MS:1000443"]})]
        StoredWaveformInverseFourierTransform,
        #[term(cv=MS, accession=1000288, name="cyclotron", flags={0}, parents={["MS:1000443"]})]
        Cyclotron,
        #[term(cv=MS, accession=1000291, name="linear ion trap", flags={0}, parents={["MS:1000264"]})]
        LinearIonTrap,
        #[term(cv=MS, accession=1000443, name="mass analyzer type", flags={0}, parents={[]})]
        MassAnalyzerType,
        #[term(cv=MS, accession=1000484, name="orbitrap", flags={0}, parents={["MS:1000443"]})]
        Orbitrap,
        #[term(cv=MS, accession=1003379, name="asymmetric track lossless time-of-flight analyzer", flags={0}, parents={["MS:1000084"]})]
        AsymmetricTrackLosslessTimeOfFlightAnalyzer,
    }
    //[[[end]]] (checksum: 6db5fba0f613bb1ce44fa5312b4383b3)
}

crate::cvmap! {
    #[flag_type=i32]
    #[allow(unused)]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_component.py', "ionization-type"]).decode('utf8')
    for line in buf.splitlines():
        cog.outl(line)
    ]]]*/
    pub enum IonizationTypeTerm {
        #[term(cv=MS, accession=1000008, name="ionization type", flags={0}, parents={[]})]
        IonizationType,
        #[term(cv=MS, accession=1000071, name="chemical ionization", flags={0}, parents={["MS:1000008"]})]
        ChemicalIonization,
        #[term(cv=MS, accession=1000073, name="electrospray ionization", flags={0}, parents={["MS:1000008"]})]
        ElectrosprayIonization,
        #[term(cv=MS, accession=1000074, name="fast atom bombardment ionization", flags={0}, parents={["MS:1000008"]})]
        FastAtomBombardmentIonization,
        #[term(cv=MS, accession=1000227, name="multiphoton ionization", flags={0}, parents={["MS:1000008"]})]
        MultiphotonIonization,
        #[term(cv=MS, accession=1000240, name="atmospheric pressure ionization", flags={0}, parents={["MS:1000008"]})]
        AtmosphericPressureIonization,
        #[term(cv=MS, accession=1000247, name="desorption ionization", flags={0}, parents={["MS:1000008"]})]
        DesorptionIonization,
        #[term(cv=MS, accession=1000255, name="flowing afterglow", flags={0}, parents={["MS:1000008"]})]
        FlowingAfterglow,
        #[term(cv=MS, accession=1000257, name="field desorption", flags={0}, parents={["MS:1000247"]})]
        FieldDesorption,
        #[term(cv=MS, accession=1000258, name="field ionization", flags={0}, parents={["MS:1000008"]})]
        FieldIonization,
        #[term(cv=MS, accession=1000259, name="glow discharge ionization", flags={0}, parents={["MS:1000008"]})]
        GlowDischargeIonization,
        #[term(cv=MS, accession=1000271, name="Negative Ion chemical ionization", flags={0}, parents={["MS:1000008"]})]
        NegativeIonChemicalIonization,
        #[term(cv=MS, accession=1000272, name="neutralization reionization mass spectrometry", flags={0}, parents={["MS:1000008"]})]
        NeutralizationReionizationMassSpectrometry,
        #[term(cv=MS, accession=1000273, name="photoionization", flags={0}, parents={["MS:1000008"]})]
        Photoionization,
        #[term(cv=MS, accession=1000274, name="pyrolysis mass spectrometry", flags={0}, parents={["MS:1000008"]})]
        PyrolysisMassSpectrometry,
        #[term(cv=MS, accession=1000276, name="resonance enhanced multiphoton ionization", flags={0}, parents={["MS:1000008"]})]
        ResonanceEnhancedMultiphotonIonization,
        #[term(cv=MS, accession=1000380, name="adiabatic ionization", flags={0}, parents={["MS:1000008"]})]
        AdiabaticIonization,
        #[term(cv=MS, accession=1000381, name="associative ionization", flags={0}, parents={["MS:1000008"]})]
        AssociativeIonization,
        #[term(cv=MS, accession=1000382, name="atmospheric pressure photoionization", flags={0}, parents={["MS:1000240"]})]
        AtmosphericPressurePhotoionization,
        #[term(cv=MS, accession=1000383, name="autodetachment", flags={0}, parents={["MS:1000008"]})]
        Autodetachment,
        #[term(cv=MS, accession=1000384, name="autoionization", flags={0}, parents={["MS:1000008"]})]
        Autoionization,
        #[term(cv=MS, accession=1000385, name="charge exchange ionization", flags={0}, parents={["MS:1000008"]})]
        ChargeExchangeIonization,
        #[term(cv=MS, accession=1000386, name="chemi-ionization", flags={0}, parents={["MS:1000008"]})]
        ChemiIonization,
        #[term(cv=MS, accession=1000387, name="desorption/ionization on silicon", flags={0}, parents={["MS:1000247"]})]
        DesorptionIonizationOnSilicon,
        #[term(cv=MS, accession=1000388, name="dissociative ionization", flags={0}, parents={["MS:1000008"]})]
        DissociativeIonization,
        #[term(cv=MS, accession=1000389, name="electron ionization", flags={0}, parents={["MS:1000008"]})]
        ElectronIonization,
        #[term(cv=MS, accession=1000393, name="laser desorption ionization", flags={0}, parents={["MS:1000247"]})]
        LaserDesorptionIonization,
        #[term(cv=MS, accession=1000395, name="liquid secondary ionization", flags={0}, parents={["MS:1000008"]})]
        LiquidSecondaryIonization,
        #[term(cv=MS, accession=1000397, name="microelectrospray", flags={0}, parents={["MS:1000073"]})]
        Microelectrospray,
        #[term(cv=MS, accession=1000398, name="nanoelectrospray", flags={0}, parents={["MS:1000073"]})]
        Nanoelectrospray,
        #[term(cv=MS, accession=1000399, name="penning ionization", flags={0}, parents={["MS:1000008"]})]
        PenningIonization,
        #[term(cv=MS, accession=1000400, name="plasma desorption ionization", flags={0}, parents={["MS:1000008"]})]
        PlasmaDesorptionIonization,
        #[term(cv=MS, accession=1000402, name="secondary ionization", flags={0}, parents={["MS:1000008"]})]
        SecondaryIonization,
        #[term(cv=MS, accession=1000403, name="soft ionization", flags={0}, parents={["MS:1000008"]})]
        SoftIonization,
        #[term(cv=MS, accession=1000404, name="spark ionization", flags={0}, parents={["MS:1000008"]})]
        SparkIonization,
        #[term(cv=MS, accession=1000405, name="surface-assisted laser desorption ionization", flags={0}, parents={["MS:1000247"]})]
        SurfaceAssistedLaserDesorptionIonization,
        #[term(cv=MS, accession=1000406, name="surface ionization", flags={0}, parents={["MS:1000008"]})]
        SurfaceIonization,
        #[term(cv=MS, accession=1000407, name="thermal ionization", flags={0}, parents={["MS:1000008"]})]
        ThermalIonization,
        #[term(cv=MS, accession=1000408, name="vertical ionization", flags={0}, parents={["MS:1000008"]})]
        VerticalIonization,
        #[term(cv=MS, accession=1000446, name="fast ion bombardment", flags={0}, parents={["MS:1000008"]})]
        FastIonBombardment,
        #[term(cv=MS, accession=1002011, name="desorption electrospray ionization", flags={0}, parents={["MS:1000240"]})]
        DesorptionElectrosprayIonization,
        #[term(cv=MS, accession=1003235, name="paper spray ionization", flags={0}, parents={["MS:1000008"]})]
        PaperSprayIonization,
        #[term(cv=MS, accession=1003248, name="proton transfer reaction", flags={0}, parents={["MS:1000008"]})]
        ProtonTransferReaction,
        #[term(cv=MS, accession=1003249, name="proton transfer charge reduction", flags={0}, parents={["MS:1000008"]})]
        ProtonTransferChargeReduction,
    }
    // [[[end]]] (checksum: 75410a9cb2735cb0ac1665d65cbf8eb7)
}

crate::cvmap! {
    #[flag_type=i32]
    #[allow(unused)]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_component.py', "inlet-type"]).decode('utf8')
    for line in buf.splitlines():
        cog.outl(line)
    ]]]*/
    pub enum InletTypeTerm {
        #[term(cv=MS, accession=1000007, name="inlet type", flags={0}, parents={[]})]
        InletType,
        #[term(cv=MS, accession=1000055, name="continuous flow fast atom bombardment", flags={0}, parents={["MS:1000007"]})]
        ContinuousFlowFastAtomBombardment,
        #[term(cv=MS, accession=1000056, name="direct inlet", flags={0}, parents={["MS:1000007"]})]
        DirectInlet,
        #[term(cv=MS, accession=1000057, name="electrospray inlet", flags={0}, parents={["MS:1000007"]})]
        ElectrosprayInlet,
        #[term(cv=MS, accession=1000058, name="flow injection analysis", flags={0}, parents={["MS:1000007"]})]
        FlowInjectionAnalysis,
        #[term(cv=MS, accession=1000059, name="inductively coupled plasma", flags={0}, parents={["MS:1000007"]})]
        InductivelyCoupledPlasma,
        #[term(cv=MS, accession=1000060, name="infusion", flags={0}, parents={["MS:1000007"]})]
        Infusion,
        #[term(cv=MS, accession=1000061, name="jet separator", flags={0}, parents={["MS:1000007"]})]
        JetSeparator,
        #[term(cv=MS, accession=1000062, name="membrane separator", flags={0}, parents={["MS:1000007"]})]
        MembraneSeparator,
        #[term(cv=MS, accession=1000063, name="moving belt", flags={0}, parents={["MS:1000007"]})]
        MovingBelt,
        #[term(cv=MS, accession=1000064, name="moving wire", flags={0}, parents={["MS:1000007"]})]
        MovingWire,
        #[term(cv=MS, accession=1000065, name="open split", flags={0}, parents={["MS:1000007"]})]
        OpenSplit,
        #[term(cv=MS, accession=1000066, name="particle beam", flags={0}, parents={["MS:1000007"]})]
        ParticleBeam,
        #[term(cv=MS, accession=1000067, name="reservoir", flags={0}, parents={["MS:1000007"]})]
        Reservoir,
        #[term(cv=MS, accession=1000068, name="septum", flags={0}, parents={["MS:1000007"]})]
        Septum,
        #[term(cv=MS, accession=1000069, name="thermospray inlet", flags={0}, parents={["MS:1000007"]})]
        ThermosprayInlet,
        #[term(cv=MS, accession=1000248, name="direct insertion probe", flags={0}, parents={["MS:1000007"]})]
        DirectInsertionProbe,
        #[term(cv=MS, accession=1000249, name="direct liquid introduction", flags={0}, parents={["MS:1000007"]})]
        DirectLiquidIntroduction,
        #[term(cv=MS, accession=1000396, name="membrane inlet", flags={0}, parents={["MS:1000007"]})]
        MembraneInlet,
        #[term(cv=MS, accession=1000485, name="nanospray inlet", flags={0}, parents={["MS:1000057"]})]
        NanosprayInlet,
    }
    // [[[end]]] (checksum: e7a44857303f45b80298f18523df0088)
}

crate::cvmap! {
    #[flag_type=i32]
    #[allow(unused)]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_component.py', "detector-type"]).decode('utf8')
    for line in buf.splitlines():
        cog.outl(line)
    ]]]*/
    pub enum DetectorTypeTerm {
        #[term(cv=MS, accession=1000026, name="detector type", flags={0}, parents={[]})]
        DetectorType,
        #[term(cv=MS, accession=1000107, name="channeltron", flags={0}, parents={["MS:1000026"]})]
        Channeltron,
        #[term(cv=MS, accession=1000110, name="daly detector", flags={0}, parents={["MS:1000026"]})]
        DalyDetector,
        #[term(cv=MS, accession=1000112, name="faraday cup", flags={0}, parents={["MS:1000026"]})]
        FaradayCup,
        #[term(cv=MS, accession=1000115, name="multi-collector", flags={0}, parents={["MS:1000026"]})]
        MultiCollector,
        #[term(cv=MS, accession=1000116, name="photomultiplier", flags={0}, parents={["MS:1000026"]})]
        Photomultiplier,
        #[term(cv=MS, accession=1000253, name="electron multiplier", flags={0}, parents={["MS:1000026"]})]
        ElectronMultiplier,
        #[term(cv=MS, accession=1000345, name="array detector", flags={0}, parents={["MS:1000026"]})]
        ArrayDetector,
        #[term(cv=MS, accession=1000346, name="conversion dynode", flags={0}, parents={["MS:1000026"]})]
        ConversionDynode,
        #[term(cv=MS, accession=1000347, name="dynode", flags={0}, parents={["MS:1000026"]})]
        Dynode,
        #[term(cv=MS, accession=1000348, name="focal plane collector", flags={0}, parents={["MS:1000026"]})]
        FocalPlaneCollector,
        #[term(cv=MS, accession=1000349, name="ion-to-photon detector", flags={0}, parents={["MS:1000026"]})]
        IonToPhotonDetector,
        #[term(cv=MS, accession=1000350, name="point collector", flags={0}, parents={["MS:1000026"]})]
        PointCollector,
        #[term(cv=MS, accession=1000351, name="postacceleration detector", flags={0}, parents={["MS:1000026"]})]
        PostaccelerationDetector,
        #[term(cv=MS, accession=1000621, name="photodiode array detector", flags={0}, parents={["MS:1000345"]})]
        PhotodiodeArrayDetector,
        #[term(cv=MS, accession=1000624, name="inductive detector", flags={0}, parents={["MS:1000026"]})]
        InductiveDetector,
        #[term(cv=MS, accession=1000818, name="Acquity UPLC PDA", flags={0}, parents={["MS:1000126", "MS:1000621"]})]
        AcquityUPLCPDA,
        #[term(cv=MS, accession=1002308, name="fluorescence detector", flags={0}, parents={["MS:1000026"]})]
        FluorescenceDetector,
    }
    //[[[end]]] (checksum: 5b1d0baa7e3919e011a06b4061265fab)
}
