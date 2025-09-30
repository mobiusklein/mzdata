use std::fmt::Display;

use crate::impl_param_described;
use crate::params::{ParamCow, ParamLike, ParamList};

/// A distinguishing tag describing the part of an instrument a [`Component`] refers to
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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

/// A description of a combination of parts that are described as part of an [`InstrumentConfiguration`].
/// There may be more than one component of the same type in a singel configuration, e.g. a triple-quad instrument
/// can have three separate [`ComponentType::Analyzer`] components.
///
/// A component may also be described by more than one [`Param`](crate::params::Param), such as the
#[derive(Default, Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Component {
    /// The kind of component this describes
    pub component_type: ComponentType,
    /// The order in the sequence of components that the analytes interact with
    pub order: u8,
    pub params: ParamList,
}

impl Component {
    pub fn mass_analyzer(&self) -> Option<MassAnalyzerTerm> {
        self.params
            .iter()
            .filter(|p| p.is_ms())
            .flat_map(|p| {
                if let Some(u) = p.accession {
                    MassAnalyzerTerm::from_accession(u)
                } else {
                    None
                }
            })
            .next()
    }

    pub fn detector(&self) -> Option<DetectorTypeTerm> {
        self.params
            .iter()
            .filter(|p| p.is_ms())
            .flat_map(|p| {
                if let Some(u) = p.accession {
                    DetectorTypeTerm::from_accession(u)
                } else {
                    None
                }
            })
            .next()
    }

    pub fn ionization_type(&self) -> Option<IonizationTypeTerm> {
        self.params
            .iter()
            .filter(|p| p.is_ms())
            .flat_map(|p| {
                if let Some(u) = p.accession {
                    IonizationTypeTerm::from_accession(u)
                } else {
                    None
                }
            })
            .next()
    }

    pub fn name(&self) -> Option<&str> {
        let it = self.params.iter().filter(|p| p.is_ms());
        match self.component_type {
            ComponentType::Analyzer => it
                .flat_map(|p| {
                    p.accession
                        .and_then(|u| MassAnalyzerTerm::from_accession(u))
                        .map(|u| u.name())
                })
                .next(),
            ComponentType::IonSource => it
                .flat_map(|p| {
                    p.accession
                        .and_then(|u| IonizationTypeTerm::from_accession(u))
                        .map(|u| u.name())
                })
                .next(),
            ComponentType::Detector => it
                .flat_map(|p| {
                    p.accession
                        .and_then(|u| DetectorTypeTerm::from_accession(u))
                        .map(|u| u.name())
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
                    p.accession.and_then(|u| {
                        MassAnalyzerTerm::from_accession(u)
                            .map(|t| t.parents().into_iter().map(|t| t.to_param()).collect())
                    })
                })
                .next()
                .unwrap_or_default(),
            ComponentType::IonSource => self
                .params
                .iter()
                .flat_map(|p| {
                    p.accession.and_then(|u| {
                        IonizationTypeTerm::from_accession(u)
                            .map(|t| t.parents().into_iter().map(|t| t.to_param()).collect())
                    })
                })
                .next()
                .unwrap_or_default(),
            ComponentType::Detector => self
                .params
                .iter()
                .flat_map(|p| {
                    p.accession.and_then(|u| {
                        DetectorTypeTerm::from_accession(u)
                            .map(|t| t.parents().into_iter().map(|t| t.to_param()).collect())
                    })
                })
                .next()
                .unwrap_or_default(),
            ComponentType::Unknown => vec![],
        }
    }
}

/// A series of mass spectrometer components that together were engaged to acquire a mass spectrum
#[derive(Default, Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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

impl InstrumentConfiguration {
    /// Add a new [`Component`] to the configuration, added at the end of the list
    pub fn new_component(&mut self, component_type: ComponentType) -> &mut Component {
        let component = Component {
            component_type,
            ..Default::default()
        };
        self.push(component);
        self.components.last_mut().unwrap()
    }

    pub fn len(&self) -> usize {
        self.components.len()
    }

    pub fn is_empty(&self) -> bool {
        self.components.is_empty()
    }

    /// Add a new [`Component`] to the end of the list, setting the [`Component::order`] field
    /// accordingly.
    pub fn push(&mut self, mut value: Component) {
        let n = self.len();
        value.order = n as u8;
        self.components.push(value)
    }

    pub fn iter(&self) -> std::slice::Iter<'_, Component> {
        self.components.iter()
    }

    pub fn last(&self) -> Option<&Component> {
        self.components.last()
    }

    pub fn last_mut(&mut self) -> Option<&mut Component> {
        self.components.last_mut()
    }
}

impl_param_described!(InstrumentConfiguration, Component);

crate::cvmap! {
    #[flag_type=i32]
    #[allow(unused)]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_component.py', "mass-analyzer"]).decode('utf8')
    for line in buf.splitlines():
        cog.outl(line)
    ]]]*/
    pub enum MassAnalyzerTerm {
        #[term(cv=MS, accession=1000078, name="axial ejection linear ion trap", flags={0}, parents={["MS:1000291"]})]
        #[doc="axial ejection linear ion trap - A linear ion trap mass spectrometer where ions are ejected along the axis of the analyzer."]
        AxialEjectionLinearIonTrap,
        #[term(cv=MS, accession=1000079, name="fourier transform ion cyclotron resonance mass spectrometer", flags={0}, parents={["MS:1000443"]})]
        #[doc="fourier transform ion cyclotron resonance mass spectrometer - A mass spectrometer based on the principle of ion cyclotron resonance in which an ion in a magnetic field moves in a circular orbit at a frequency characteristic of its m/z value. Ions are coherently excited to a larger radius orbit using a pulse of radio frequency energy and their image charge is detected on receiver plates as a time domain signal. Fourier transformation of the time domain signal results in a frequency domain signal which is converted to a mass spectrum based in the inverse relationship between frequency and m/z."]
        FourierTransformIonCyclotronResonanceMassSpectrometer,
        #[term(cv=MS, accession=1000080, name="magnetic sector", flags={0}, parents={["MS:1000443"]})]
        #[doc="magnetic sector - A device that produces a magnetic field perpendicular to a charged particle beam that deflects the beam to an extent that is proportional to the particle momentum per unit charge. For a monoenergetic beam, the deflection is proportional to m/z."]
        MagneticSector,
        #[term(cv=MS, accession=1000081, name="quadrupole", flags={0}, parents={["MS:1000443"]})]
        #[doc="quadrupole - A mass spectrometer that consists of four parallel rods whose centers form the corners of a square and whose opposing poles are connected. The voltage applied to the rods is a superposition of a static potential and a sinusoidal radio frequency potential. The motion of an ion in the x and y dimensions is described by the Matthieu equation whose solutions show that ions in a particular m/z range can be transmitted along the z axis."]
        Quadrupole,
        #[term(cv=MS, accession=1000082, name="quadrupole ion trap", flags={0}, parents={["MS:1000264"]})]
        #[doc="quadrupole ion trap - Quadrupole Ion Trap mass analyzer captures the ions in a three dimensional ion trap and then selectively ejects them by varying the RF and DC potentials."]
        QuadrupoleIonTrap,
        #[term(cv=MS, accession=1000083, name="radial ejection linear ion trap", flags={0}, parents={["MS:1000291"]})]
        #[doc="radial ejection linear ion trap - A linear ion trap mass spectrometer where ions are ejected along the radius of the analyzer."]
        RadialEjectionLinearIonTrap,
        #[term(cv=MS, accession=1000084, name="time-of-flight", flags={0}, parents={["MS:1000443"]})]
        #[doc="time-of-flight - Instrument that separates ions by m/z in a field-free region after acceleration to a fixed acceleration energy."]
        TimeOfFlight,
        #[term(cv=MS, accession=1000254, name="electrostatic energy analyzer", flags={0}, parents={["MS:1000443"]})]
        #[doc="electrostatic energy analyzer - A device consisting of conducting parallel plates, concentric cylinders or concentric spheres that separates charged particles according to their kinetic energy by means of an electric field that is constant in time."]
        ElectrostaticEnergyAnalyzer,
        #[term(cv=MS, accession=1000264, name="ion trap", flags={0}, parents={["MS:1000443"]})]
        #[doc="ion trap - A device for spatially confining ions using electric and magnetic fields alone or in combination."]
        IonTrap,
        #[term(cv=MS, accession=1000284, name="stored waveform inverse fourier transform", flags={0}, parents={["MS:1000443"]})]
        #[doc="stored waveform inverse fourier transform - A technique to create excitation waveforms for ions in FT-ICR mass spectrometer or Paul ion trap. An excitation waveform in the time-domain is generated by taking the inverse Fourier transform of an appropriate frequency-domain programmed excitation spectrum, in which the resonance frequencies of ions to be excited are included. This technique may be used for selection of precursor ions in MS2 experiments."]
        StoredWaveformInverseFourierTransform,
        #[term(cv=MS, accession=1000288, name="cyclotron", flags={0}, parents={["MS:1000443"]})]
        #[doc="cyclotron - A device that uses an oscillating electric field and magnetic field to accelerate charged particles."]
        Cyclotron,
        #[term(cv=MS, accession=1000291, name="linear ion trap", flags={0}, parents={["MS:1000264"]})]
        #[doc="linear ion trap - A two dimensional Paul ion trap in which ions are confined in the axial dimension by means of an electric field at the ends of the trap."]
        LinearIonTrap,
        #[term(cv=MS, accession=1000443, name="mass analyzer type", flags={0}, parents={[]})]
        #[doc="mass analyzer type - Mass analyzer separates the ions according to their mass-to-charge ratio."]
        MassAnalyzerType,
        #[term(cv=MS, accession=1000484, name="orbitrap", flags={0}, parents={["MS:1000443"]})]
        #[doc="orbitrap - An ion trapping device that consists of an outer barrel-like electrode and a coaxial inner spindle-like electrode that form an electrostatic field with quadro-logarithmic potential distribution. The frequency of harmonic oscillations of the orbitally trapped ions along the axis of the electrostatic field is independent of the ion velocity and is inversely proportional to the square root of m/z so that the trap can be used as a mass analyzer."]
        Orbitrap,
        #[term(cv=MS, accession=1003379, name="asymmetric track lossless time-of-flight analyzer", flags={0}, parents={["MS:1000084"]})]
        #[doc="asymmetric track lossless time-of-flight analyzer - A TOF-like mass analyzer with asymmetric ion mirrors to direct ions into transversal asymmetric oscillations and ion foil shapes and maintains ion packet for transmission and resolution."]
        AsymmetricTrackLosslessTimeOfFlightAnalyzer,
    }
    //[[[end]]] (sum: xfYRYQf4J4)
}

crate::cvmap! {
    #[flag_type=i32]
    #[allow(unused)]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_component.py', "ionization-type"]).decode('utf8')
    for line in buf.splitlines():
        cog.outl(line)
    ]]]*/
    pub enum IonizationTypeTerm {
        #[term(cv=MS, accession=1000008, name="ionization type", flags={0}, parents={[]})]
        #[doc="ionization type - The method by which gas phase ions are generated from the sample."]
        IonizationType,
        #[term(cv=MS, accession=1000070, name="atmospheric pressure chemical ionization", flags={0}, parents={["MS:1000240"]})]
        #[doc="atmospheric pressure chemical ionization - Chemical ionization that takes place at atmospheric pressure as opposed to the reduced pressure is normally used for chemical ionization."]
        AtmosphericPressureChemicalIonization,
        #[term(cv=MS, accession=1000071, name="chemical ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="chemical ionization - The formation of a new ion by the reaction of a neutral species with an ion. The process may involve transfer of an electron, a proton or other charged species between the reactants. When a positive ion results from chemical ionization the term may be used without qualification. When a negative ion results the term negative ion chemical ionization should be used. Note that this term is not synonymous with chemi-ionization."]
        ChemicalIonization,
        #[term(cv=MS, accession=1000073, name="electrospray ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="electrospray ionization - A process in which ionized species in the gas phase are produced from an analyte-containing solution via highly charged fine droplets, by means of spraying the solution from a narrow-bore needle tip at atmospheric pressure in the presence of a high electric field. When a pressurized gas is used to aid in the formation of a stable spray, the term pneumatically assisted electrospray ionization is used. The term ion spray is not recommended."]
        ElectrosprayIonization,
        #[term(cv=MS, accession=1000074, name="fast atom bombardment ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="fast atom bombardment ionization - The ionization of any species by the interaction of a focused beam of neutral atoms having a translational energy of several thousand eV with a sample that is typically dissolved in a solvent matrix. See also secondary ionization."]
        FastAtomBombardmentIonization,
        #[term(cv=MS, accession=1000075, name="matrix-assisted laser desorption ionization", flags={0}, parents={["MS:1000247"]})]
        #[doc="matrix-assisted laser desorption ionization - The formation of gas-phase ions from molecules that are present in a solid or solvent matrix that is irradiated with a pulsed laser. See also laser desorption/ionization."]
        MatrixAssistedLaserDesorptionIonization,
        #[term(cv=MS, accession=1000227, name="multiphoton ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="multiphoton ionization - Photoionization of an atom or molecule in which in two or more photons are absorbed."]
        MultiphotonIonization,
        #[term(cv=MS, accession=1000239, name="atmospheric pressure matrix-assisted laser desorption ionization", flags={0}, parents={["MS:1000240"]})]
        #[doc="atmospheric pressure matrix-assisted laser desorption ionization - Matrix-assisted laser desorption ionization in which the sample target is at atmospheric pressure and the ions formed by the pulsed laser are sampled through a small aperture into the mass spectrometer."]
        AtmosphericPressureMatrixAssistedLaserDesorptionIonization,
        #[term(cv=MS, accession=1000240, name="atmospheric pressure ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="atmospheric pressure ionization - Any ionization process in which ions are formed in the gas phase at atmospheric pressure."]
        AtmosphericPressureIonization,
        #[term(cv=MS, accession=1000247, name="desorption ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="desorption ionization - The formation of ions from a solid or liquid material after the rapid vaporization of that sample."]
        DesorptionIonization,
        #[term(cv=MS, accession=1000255, name="flowing afterglow", flags={0}, parents={["MS:1000008"]})]
        #[doc="flowing afterglow - An ion source immersed in a flow of helium or other inert buffer gas that carries the ions through a meter-long reactor at pressures around 100 Pa."]
        FlowingAfterglow,
        #[term(cv=MS, accession=1000257, name="field desorption", flags={0}, parents={["MS:1000247"]})]
        #[doc="field desorption - The formation of gas-phase ions from a material deposited on a solid surface in the presence of a high electric field. Because this process may encompass ionization by field ionization or other mechanisms, it is not recommended as a synonym for field desorption ionization."]
        FieldDesorption,
        #[term(cv=MS, accession=1000258, name="field ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="field ionization - The removal of electrons from any species by interaction with a high electric field."]
        FieldIonization,
        #[term(cv=MS, accession=1000259, name="glow discharge ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="glow discharge ionization - The formation of ions in the gas phase and from solid samples at the cathode by application of a voltage to a low pressure gas."]
        GlowDischargeIonization,
        #[term(cv=MS, accession=1000271, name="Negative Ion chemical ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="Negative Ion chemical ionization - Chemical ionization that results in the formation of negative ions."]
        NegativeIonChemicalIonization,
        #[term(cv=MS, accession=1000272, name="neutralization reionization mass spectrometry", flags={0}, parents={["MS:1000008"]})]
        #[doc="neutralization reionization mass spectrometry - With this technique, m/z selected ions form neutrals by charge transfer to a collision gas or by dissociation. The neutrals are separated from the remaining ions and ionized in collisions with a second gas. This method is used to investigate reaction intermediates and other unstable species."]
        NeutralizationReionizationMassSpectrometry,
        #[term(cv=MS, accession=1000273, name="photoionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="photoionization - The ionization of an atom or molecule by a photon, written M + h? ? M^+ + e. The term photon impact is not recommended."]
        Photoionization,
        #[term(cv=MS, accession=1000274, name="pyrolysis mass spectrometry", flags={0}, parents={["MS:1000008"]})]
        #[doc="pyrolysis mass spectrometry - A mass spectrometry technique in which the sample is heated to the point of decomposition and the gaseous decomposition products are introduced into the ion source."]
        PyrolysisMassSpectrometry,
        #[term(cv=MS, accession=1000276, name="resonance enhanced multiphoton ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="resonance enhanced multiphoton ionization - Multiphoton ionization in which the ionization cross section is significantly enhanced because the energy of the incident photons is resonant with an intermediate excited state of the neutral species."]
        ResonanceEnhancedMultiphotonIonization,
        #[term(cv=MS, accession=1000278, name="surface enhanced laser desorption ionization", flags={0}, parents={["MS:1000406"]})]
        #[doc="surface enhanced laser desorption ionization - The formation of ionized species in the gas phase from analytes deposited on a particular surface substrate which is irradiated with a laser beam of which wavelength is absorbed by the surface. See also desorption/ionization on silicon and laser desorption/ionization."]
        SurfaceEnhancedLaserDesorptionIonization,
        #[term(cv=MS, accession=1000279, name="surface enhanced neat desorption", flags={0}, parents={["MS:1000406"]})]
        #[doc="surface enhanced neat desorption - Matrix-assisted laser desorption ionization in which the matrix is covalently linked to the target surface."]
        SurfaceEnhancedNeatDesorption,
        #[term(cv=MS, accession=1000380, name="adiabatic ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="adiabatic ionization - A process whereby an electron is removed from an atom, ion, or molecule to produce an ion in its lowest energy state."]
        AdiabaticIonization,
        #[term(cv=MS, accession=1000381, name="associative ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="associative ionization - An ionization process in which two excited atoms or molecules react to form a single positive ion and an electron."]
        AssociativeIonization,
        #[term(cv=MS, accession=1000382, name="atmospheric pressure photoionization", flags={0}, parents={["MS:1000240"]})]
        #[doc="atmospheric pressure photoionization - Atmospheric pressure chemical ionization in which the reactant ions are generated by photo-ionization."]
        AtmosphericPressurePhotoionization,
        #[term(cv=MS, accession=1000383, name="autodetachment", flags={0}, parents={["MS:1000008"]})]
        #[doc="autodetachment - The formation of a neutral when a negative ion in a discrete state with an energy greater than the detachment threshold loses an electron spontaneously without further interaction with an energy source."]
        Autodetachment,
        #[term(cv=MS, accession=1000384, name="autoionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="autoionization - The formation of an ion when an atom or molecule in a discrete state with an energy greater than the ionization threshold loses an electron spontaneously without further interaction with an energy source."]
        Autoionization,
        #[term(cv=MS, accession=1000385, name="charge exchange ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="charge exchange ionization - The interaction of an ion with an atom or molecule in which the charge on the ion is transferred to the neutral without the dissociation of either. Synonymous with charge transfer ionization."]
        ChargeExchangeIonization,
        #[term(cv=MS, accession=1000386, name="chemi-ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="chemi-ionization - The reaction of a neutral molecule with an internally excited molecule to form an ion. Note that this term is not synonymous with chemical ionization."]
        ChemiIonization,
        #[term(cv=MS, accession=1000387, name="desorption/ionization on silicon", flags={0}, parents={["MS:1000247"]})]
        #[doc="desorption/ionization on silicon - The formation of ions by laser desorption ionization of a sample deposited on a porous silicon surface."]
        DesorptionIonizationOnSilicon,
        #[term(cv=MS, accession=1000388, name="dissociative ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="dissociative ionization - The reaction of a gas-phase molecule that results in its decomposition to form products, one of which is an ion."]
        DissociativeIonization,
        #[term(cv=MS, accession=1000389, name="electron ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="electron ionization - The ionization of an atom or molecule by electrons that are typically accelerated to energies between 50 and 150 eV. Usually 70 eV electrons are used to produce positive ions. The term 'electron impact' is not recommended."]
        ElectronIonization,
        #[term(cv=MS, accession=1000393, name="laser desorption ionization", flags={0}, parents={["MS:1000247"]})]
        #[doc="laser desorption ionization - The formation of gas-phase ions by the interaction of a pulsed laser with a solid or liquid material."]
        LaserDesorptionIonization,
        #[term(cv=MS, accession=1000395, name="liquid secondary ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="liquid secondary ionization - The ionization of any species by the interaction of a focused beam of ions with a sample that is dissolved in a solvent matrix. See also fast atom bombardment and secondary ionization."]
        LiquidSecondaryIonization,
        #[term(cv=MS, accession=1000397, name="microelectrospray", flags={0}, parents={["MS:1000073"]})]
        #[doc="microelectrospray - Electrospray ionization at a solvent flow rate of 300-800 nL/min where the flow is a result of a mechanical pump. See nanoelectrospray."]
        Microelectrospray,
        #[term(cv=MS, accession=1000398, name="nanoelectrospray", flags={0}, parents={["MS:1000073"]})]
        #[doc="nanoelectrospray - Electrospray ionization at a flow rate less than ~25 nL/min. Nanoelectrospray is synonymous with nanospray. The flow is dependent on the potenial on the tip of the electrospray needle and/or a gas presure to push the sample through the needle. See also electrospray ionization and microelectrospray."]
        Nanoelectrospray,
        #[term(cv=MS, accession=1000399, name="penning ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="penning ionization - Ionization that occurs through the interaction of two or more neutral gaseous species, at least one of which is internally excited."]
        PenningIonization,
        #[term(cv=MS, accession=1000400, name="plasma desorption ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="plasma desorption ionization - The ionization of material in a solid sample by bombarding it with ionic or neutral atoms formed as a result of the fission of a suitable nuclide, typically 252Cf. Synonymous with fission fragment ionization."]
        PlasmaDesorptionIonization,
        #[term(cv=MS, accession=1000402, name="secondary ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="secondary ionization - The process in which ions are ejected from a sample surface as a result of bombardment by a primary beam of atoms or ions."]
        SecondaryIonization,
        #[term(cv=MS, accession=1000403, name="soft ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="soft ionization - The formation of gas-phase ions without extensive fragmentation."]
        SoftIonization,
        #[term(cv=MS, accession=1000404, name="spark ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="spark ionization - The formation of ions from a solid material by an intermittent electrical discharge."]
        SparkIonization,
        #[term(cv=MS, accession=1000405, name="surface-assisted laser desorption ionization", flags={0}, parents={["MS:1000247"]})]
        #[doc="surface-assisted laser desorption ionization - The formation of gas-phase ions from molecules that are deposited on a particular surface substrate that is irradiated with a pulsed laser. See also matrix-assisted laser desorption ionization."]
        SurfaceAssistedLaserDesorptionIonization,
        #[term(cv=MS, accession=1000406, name="surface ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="surface ionization - The ionization of a neutral species when it interacts with a solid surface with an appropriate work function and temperature."]
        SurfaceIonization,
        #[term(cv=MS, accession=1000407, name="thermal ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="thermal ionization - The ionization of a neutral species through contact with a high temperature surface."]
        ThermalIonization,
        #[term(cv=MS, accession=1000408, name="vertical ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="vertical ionization - A process in which an electron is removed from or added to a molecule without a change in the positions of the atoms. The resulting ion is typically in an excited vibrational state."]
        VerticalIonization,
        #[term(cv=MS, accession=1000446, name="fast ion bombardment", flags={0}, parents={["MS:1000008"]})]
        #[doc="fast ion bombardment - The ionization of any species by the interaction of a focused beam of ions having a translational energy of several thousand eV with a solid sample."]
        FastIonBombardment,
        #[term(cv=MS, accession=1002011, name="desorption electrospray ionization", flags={0}, parents={["MS:1000240"]})]
        #[doc="desorption electrospray ionization - Combination of electrospray and desorption ionization method that ionizes gases, liquids and solids in open air under atmospheric pressure."]
        DesorptionElectrosprayIonization,
        #[term(cv=MS, accession=1003235, name="paper spray ionization", flags={0}, parents={["MS:1000008"]})]
        #[doc="paper spray ionization - The ionization of analytes from a piece of paper by applying a solvent and voltage."]
        PaperSprayIonization,
        #[term(cv=MS, accession=1003248, name="proton transfer reaction", flags={0}, parents={["MS:1000008"]})]
        #[doc="proton transfer reaction - Process to transfer a proton from a hydronium ion (H3O+) to neutral analyte, leading to a protonated analyte, which typically does not lead to fragmentation."]
        ProtonTransferReaction,
        #[term(cv=MS, accession=1003249, name="proton transfer charge reduction", flags={0}, parents={["MS:1000008"]})]
        #[doc="proton transfer charge reduction - Process to transfer one or more protons from a multiply charged cation (peptide or protein ion) to a proton acceptor anion or neutral basic compound, thereby reducing the charge of the original analyte."]
        ProtonTransferChargeReduction,
    }
    // [[[end]]] (sum: sgB0zq2JzK)
}

crate::cvmap! {
    #[flag_type=i32]
    #[allow(unused)]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_component.py', "inlet-type"]).decode('utf8')
    for line in buf.splitlines():
        cog.outl(line)
    ]]]*/
    pub enum InletTypeTerm {
        #[term(cv=MS, accession=1000007, name="inlet type", flags={0}, parents={[]})]
        #[doc="inlet type - The nature of the sample inlet."]
        InletType,
        #[term(cv=MS, accession=1000055, name="continuous flow fast atom bombardment", flags={0}, parents={["MS:1000007"]})]
        #[doc="continuous flow fast atom bombardment - Fast atom bombardment ionization in which the analyte in solution is entrained in a flowing liquid matrix."]
        ContinuousFlowFastAtomBombardment,
        #[term(cv=MS, accession=1000056, name="direct inlet", flags={0}, parents={["MS:1000007"]})]
        #[doc="direct inlet - The sample is directly inserted into the ion source, usually on the end of a heatable probe."]
        DirectInlet,
        #[term(cv=MS, accession=1000057, name="electrospray inlet", flags={0}, parents={["MS:1000007"]})]
        #[doc="electrospray inlet - Inlet used for introducing the liquid sample into an electrospray ionization source."]
        ElectrosprayInlet,
        #[term(cv=MS, accession=1000058, name="flow injection analysis", flags={0}, parents={["MS:1000007"]})]
        #[doc="flow injection analysis - Sample is directly injected or infused into the ionization source."]
        FlowInjectionAnalysis,
        #[term(cv=MS, accession=1000059, name="inductively coupled plasma", flags={0}, parents={["MS:1000007"]})]
        #[doc="inductively coupled plasma - A gas discharge ion source in which the energy to the plasma is supplied by electromagnetic induction."]
        InductivelyCoupledPlasma,
        #[term(cv=MS, accession=1000060, name="infusion", flags={0}, parents={["MS:1000007"]})]
        #[doc="infusion - The continuous flow of solution of a sample into the ionization source."]
        Infusion,
        #[term(cv=MS, accession=1000061, name="jet separator", flags={0}, parents={["MS:1000007"]})]
        #[doc="jet separator - A device that separates carrier gas from gaseous analyte molecules on the basis of diffusivity."]
        JetSeparator,
        #[term(cv=MS, accession=1000062, name="membrane separator", flags={0}, parents={["MS:1000007"]})]
        #[doc="membrane separator - A device to separate carrier molecules from analyte molecules on the basis of ease of diffusion across a semipermeable membrane."]
        MembraneSeparator,
        #[term(cv=MS, accession=1000063, name="moving belt", flags={0}, parents={["MS:1000007"]})]
        #[doc="moving belt - Continuous moving surface in the form of a belt which passes through an ion source carrying analyte molecules."]
        MovingBelt,
        #[term(cv=MS, accession=1000064, name="moving wire", flags={0}, parents={["MS:1000007"]})]
        #[doc="moving wire - Continuous moving surface in the form of a wire which passes through an ion source carrying analyte molecules."]
        MovingWire,
        #[term(cv=MS, accession=1000065, name="open split", flags={0}, parents={["MS:1000007"]})]
        #[doc="open split - A division of flowing stream of liquid into two streams."]
        OpenSplit,
        #[term(cv=MS, accession=1000066, name="particle beam", flags={0}, parents={["MS:1000007"]})]
        #[doc="particle beam - Method for generating ions from a solution of an analyte."]
        ParticleBeam,
        #[term(cv=MS, accession=1000067, name="reservoir", flags={0}, parents={["MS:1000007"]})]
        #[doc="reservoir - A sample inlet method involving a reservoir."]
        Reservoir,
        #[term(cv=MS, accession=1000068, name="septum", flags={0}, parents={["MS:1000007"]})]
        #[doc="septum - A disc composed of a flexible material that seals the entrance to the reservoir. Can also be entrance to the vacuum chamber."]
        Septum,
        #[term(cv=MS, accession=1000069, name="thermospray inlet", flags={0}, parents={["MS:1000007"]})]
        #[doc="thermospray inlet - A method for generating gas phase ions from a solution of an analyte by rapid heating of the sample."]
        ThermosprayInlet,
        #[term(cv=MS, accession=1000248, name="direct insertion probe", flags={0}, parents={["MS:1000007"]})]
        #[doc="direct insertion probe - A device for introducing a solid or liquid sample into a mass spectrometer ion source for desorption ionization."]
        DirectInsertionProbe,
        #[term(cv=MS, accession=1000249, name="direct liquid introduction", flags={0}, parents={["MS:1000007"]})]
        #[doc="direct liquid introduction - The delivery of a liquid sample into a mass spectrometer for spray or desorption ionization."]
        DirectLiquidIntroduction,
        #[term(cv=MS, accession=1000396, name="membrane inlet", flags={0}, parents={["MS:1000007"]})]
        #[doc="membrane inlet - A semi-permeable membrane separator that permits the passage of gas sample directly to the mass spectrometer ion source."]
        MembraneInlet,
        #[term(cv=MS, accession=1000485, name="nanospray inlet", flags={0}, parents={["MS:1000057"]})]
        #[doc="nanospray inlet - Nanospray Inlet."]
        NanosprayInlet,
    }
    // [[[end]]] (sum: TlXp3x9TQi)
}

crate::cvmap! {
    #[flag_type=i32]
    #[allow(unused)]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_component.py', "detector-type"]).decode('utf8')
    for line in buf.splitlines():
        cog.outl(line)
    ]]]*/
    pub enum DetectorTypeTerm {
        #[term(cv=MS, accession=1000026, name="detector type", flags={0}, parents={[]})]
        #[doc="detector type - Type of detector used in the mass spectrometer."]
        DetectorType,
        #[term(cv=MS, accession=1000107, name="channeltron", flags={0}, parents={["MS:1000026"]})]
        #[doc="channeltron - A horn-shaped (or cone-shaped) continuous dynode particle multiplier. The ion strikes the inner surface of the device and induces the production of secondary electrons that in turn impinge on the inner surfaces to produce more secondary electrons. This avalanche effect produces an increase in signal in the final measured current pulse."]
        Channeltron,
        #[term(cv=MS, accession=1000108, name="conversion dynode electron multiplier", flags={0}, parents={["MS:1000346"]})]
        #[doc="conversion dynode electron multiplier - A surface that is held at high potential so that ions striking the surface produce electrons that are subsequently detected."]
        ConversionDynodeElectronMultiplier,
        #[term(cv=MS, accession=1000109, name="conversion dynode photomultiplier", flags={0}, parents={["MS:1000346"]})]
        #[doc="conversion dynode photomultiplier - A detector in which ions strike a conversion dynode to produce electrons that in turn generate photons through a phosphorescent screen that are detected by a photomultiplier."]
        ConversionDynodePhotomultiplier,
        #[term(cv=MS, accession=1000110, name="daly detector", flags={0}, parents={["MS:1000026"]})]
        #[doc="daly detector - Detector consisting of a conversion dynode, scintillator and photomultiplier. The metal knob at high potential emits secondary electrons when ions impinge on the surface. The secondary electrons are accelerated onto the scintillator that produces light that is then detected by the photomultiplier detector."]
        DalyDetector,
        #[term(cv=MS, accession=1000111, name="electron multiplier tube", flags={0}, parents={["MS:1000253"]})]
        #[doc="electron multiplier tube - A device to amplify the current of a beam or packet of charged particles or photons by incidence upon the surface of an electrode to produce secondary electrons."]
        ElectronMultiplierTube,
        #[term(cv=MS, accession=1000112, name="faraday cup", flags={0}, parents={["MS:1000026"]})]
        #[doc="faraday cup - A conducting cup or chamber that intercepts a charged particle beam and is electrically connected to a current measuring device."]
        FaradayCup,
        #[term(cv=MS, accession=1000113, name="focal plane array", flags={0}, parents={["MS:1000348"]})]
        #[doc="focal plane array - An array of detectors for spatially disperse ion beams in which all ions simultaneously impinge on the detector plane."]
        FocalPlaneArray,
        #[term(cv=MS, accession=1000114, name="microchannel plate detector", flags={0}, parents={["MS:1000345"]})]
        #[doc="microchannel plate detector - A thin plate that contains a closely spaced array of channels that each act as a continuous dynode particle multiplier. A charged particle, fast neutral particle, or photon striking the plate causes a cascade of secondary electrons that ultimately exits the opposite side of the plate."]
        MicrochannelPlateDetector,
        #[term(cv=MS, accession=1000115, name="multi-collector", flags={0}, parents={["MS:1000026"]})]
        #[doc="multi-collector - A detector system commonly used in inductively coupled plasma mass spectrometers."]
        MultiCollector,
        #[term(cv=MS, accession=1000116, name="photomultiplier", flags={0}, parents={["MS:1000026"]})]
        #[doc="photomultiplier - A detector for conversion of the ion/electron signal into photon(s) which are then amplified and detected."]
        Photomultiplier,
        #[term(cv=MS, accession=1000253, name="electron multiplier", flags={0}, parents={["MS:1000026"]})]
        #[doc="electron multiplier - A device to amplify the current of a beam or packet of charged particles or photons by incidence upon the surface of an electrode to produce secondary electrons. The secondary electrons are then accelerated to other electrodes or parts of a continuous electrode to produce further secondary electrons."]
        ElectronMultiplier,
        #[term(cv=MS, accession=1000345, name="array detector", flags={0}, parents={["MS:1000026"]})]
        #[doc="array detector - Detector comprising several ion collection elements, arranged in a line or grid where each element is an individual detector."]
        ArrayDetector,
        #[term(cv=MS, accession=1000346, name="conversion dynode", flags={0}, parents={["MS:1000026"]})]
        #[doc="conversion dynode - A surface that is held at high potential such that ions striking the surface produce electrons that are subsequently detected."]
        ConversionDynode,
        #[term(cv=MS, accession=1000347, name="dynode", flags={0}, parents={["MS:1000026"]})]
        #[doc="dynode - One of a series of electrodes in a photomultiplier tube. Such an arrangement is able to amplify the current emitted by the photocathode."]
        Dynode,
        #[term(cv=MS, accession=1000348, name="focal plane collector", flags={0}, parents={["MS:1000026"]})]
        #[doc="focal plane collector - A detector for spatially disperse ion beams in which all ions simultaneously impinge on the detector plane."]
        FocalPlaneCollector,
        #[term(cv=MS, accession=1000349, name="ion-to-photon detector", flags={0}, parents={["MS:1000026"]})]
        #[doc="ion-to-photon detector - A detector in which ions strike a conversion dynode to produce electrons that in turn strike a phosphor and the resulting photons are detected by a photomultiplier."]
        IonToPhotonDetector,
        #[term(cv=MS, accession=1000350, name="point collector", flags={0}, parents={["MS:1000026"]})]
        #[doc="point collector - A detector in which the ion beam is focused onto a point and the individual ions arrive sequentially."]
        PointCollector,
        #[term(cv=MS, accession=1000351, name="postacceleration detector", flags={0}, parents={["MS:1000026"]})]
        #[doc="postacceleration detector - A detector in which the charged particles are accelerated to a high velocity and impinge on a conversion dynode, emitting secondary electrons. The electrons are accelerated onto a phosphor screen, which emits photons that are in turn detected using a photomultiplier or other photon detector."]
        PostaccelerationDetector,
        #[term(cv=MS, accession=1000621, name="photodiode array detector", flags={0}, parents={["MS:1000345"]})]
        #[doc="photodiode array detector - An array detector used to record spectra in the ultraviolet and visible region of light."]
        PhotodiodeArrayDetector,
        #[term(cv=MS, accession=1000624, name="inductive detector", flags={0}, parents={["MS:1000026"]})]
        #[doc="inductive detector - Inductive detector."]
        InductiveDetector,
        #[term(cv=MS, accession=1000818, name="Acquity UPLC PDA", flags={0}, parents={["MS:1000126", "MS:1000621"]})]
        #[doc="Acquity UPLC PDA - Acquity UPLC Photodiode Array Detector."]
        AcquityUPLCPDA,
        #[term(cv=MS, accession=1000819, name="Acquity UPLC FLR", flags={0}, parents={["MS:1000126", "MS:1002308"]})]
        #[doc="Acquity UPLC FLR - Acquity UPLC Fluorescence Detector."]
        AcquityUPLCFLR,
        #[term(cv=MS, accession=1002308, name="fluorescence detector", flags={0}, parents={["MS:1000026"]})]
        #[doc="fluorescence detector - A detector using a fluorescent signal after excitation with light."]
        FluorescenceDetector,
    }
    //[[[end]]] (sum: 6NYSHJaWvt)
}
