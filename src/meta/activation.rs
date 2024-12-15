use std::fmt::Display;



crate::cvmap! {
    #[flag_type=i32]
    #[allow(unused)]
    #[doc = "A method used for dissociation or fragmentation."]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_activation.py']).decode('utf8')
    for line in buf.splitlines():
        cog.outl(line)
    ]]]*/
    pub enum DissociationMethodTerm {
        #[term(cv=MS, accession=1000044, name="dissociation method", flags={0}, parents={[]})]
        #[doc = "dissociation method - Fragmentation method used for dissociation or fragmentation."]
        DissociationMethod,
        #[term(cv=MS, accession=1000133, name="collision-induced dissociation", flags={0}, parents={["MS:1000044"]})]
        #[doc = "collision-induced dissociation - The dissociation of an ion after collisional excitation. The term collisional-activated dissociation is not recommended."]
        CollisionInducedDissociation,
        #[term(cv=MS, accession=1000134, name="plasma desorption", flags={0}, parents={["MS:1000044"]})]
        #[doc = "plasma desorption - The ionization of material in a solid sample by bombarding it with ionic or neutral atoms formed as a result of the fission of a suitable nuclide, typically 252Cf. Synonymous with fission fragment ionization."]
        PlasmaDesorption,
        #[term(cv=MS, accession=1000135, name="post-source decay", flags={0}, parents={["MS:1000044"]})]
        #[doc = "post-source decay - A technique specific to reflectron time-of-flight mass spectrometers where product ions of metastable transitions or collision-induced dissociations generated in the drift tube prior to entering the reflectron are m/z separated to yield product ion spectra."]
        PostSourceDecay,
        #[term(cv=MS, accession=1000136, name="surface-induced dissociation", flags={0}, parents={["MS:1000044"]})]
        #[doc = "surface-induced dissociation - Fragmentation that results from the collision of an ion with a surface."]
        SurfaceInducedDissociation,
        #[term(cv=MS, accession=1000242, name="blackbody infrared radiative dissociation", flags={0}, parents={["MS:1000044"]})]
        #[doc = "blackbody infrared radiative dissociation - A special case of infrared multiphoton dissociation wherein excitation of the reactant ion is caused by absorption of infrared photons radiating from heated blackbody surroundings, which are usually the walls of a vacuum chamber. See also infrared multiphoton dissociation."]
        BlackbodyInfraredRadiativeDissociation,
        #[term(cv=MS, accession=1000250, name="electron capture dissociation", flags={0}, parents={["MS:1000044"]})]
        #[doc = "electron capture dissociation - A process in which a multiply protonated molecules interacts with a low energy electrons. Capture of the electron leads the liberation of energy and a reduction in charge state of the ion with the production of the (M + nH) (n-1)+ odd electron ion, which readily fragments."]
        ElectronCaptureDissociation,
        #[term(cv=MS, accession=1000262, name="infrared multiphoton dissociation", flags={0}, parents={["MS:1000435"]})]
        #[doc = "infrared multiphoton dissociation - Multiphoton ionization where the reactant ion dissociates as a result of the absorption of multiple infrared photons."]
        InfraredMultiphotonDissociation,
        #[term(cv=MS, accession=1000282, name="sustained off-resonance irradiation", flags={0}, parents={["MS:1000044"]})]
        #[doc = "sustained off-resonance irradiation - A technique associated with Fourier transform ion cyclotron resonance (FT-ICR) mass spectrometry to carry out ion/neutral reactions such as low-energy collision-induced dissociation. A radio-frequency electric field of slightly off-resonance to the cyclotron frequency of the reactant ion cyclically accelerates and decelerates the reactant ion that is confined in the Penning ion trap. The ion's orbit does not exceed the dimensions of ion trap while the ion undergoes an ion/neutral species process that produces a high average translational energy for an extended time."]
        SustainedOffResonanceIrradiation,
        #[term(cv=MS, accession=1000422, name="beam-type collision-induced dissociation", flags={0}, parents={["MS:1000133"]})]
        #[doc = "beam-type collision-induced dissociation - A collision-induced dissociation process that occurs in a beam-type collision cell."]
        BeamTypeCollisionInducedDissociation,
        #[term(cv=MS, accession=1000433, name="low-energy collision-induced dissociation", flags={0}, parents={["MS:1000044"]})]
        #[doc = "low-energy collision-induced dissociation - A collision-induced dissociation process wherein the precursor ion has the translational energy lower than approximately 1000 eV. This process typically requires multiple collisions and the collisional excitation is cumulative."]
        LowEnergyCollisionInducedDissociation,
        #[term(cv=MS, accession=1000435, name="photodissociation", flags={0}, parents={["MS:1000044"]})]
        #[doc = "photodissociation - A process wherein the reactant ion is dissociated as a result of absorption of one or more photons."]
        Photodissociation,
        #[term(cv=MS, accession=1000598, name="electron transfer dissociation", flags={0}, parents={["MS:1000044"]})]
        #[doc = "electron transfer dissociation - A process to fragment ions in a mass spectrometer by inducing fragmentation of cations (e.g. peptides or proteins) by transferring electrons from radical-anions."]
        ElectronTransferDissociation,
        #[term(cv=MS, accession=1000599, name="pulsed q dissociation", flags={0}, parents={["MS:1000044"]})]
        #[doc = "pulsed q dissociation - A process that involves precursor ion activation at high Q, a time delay to allow the precursor to fragment, then a rapid pulse to low Q where all fragment ions are trapped. The product ions can then be scanned out of the ion trap and detected."]
        PulsedQDissociation,
        #[term(cv=MS, accession=1001880, name="in-source collision-induced dissociation", flags={0}, parents={["MS:1000044"]})]
        #[doc = "in-source collision-induced dissociation - The dissociation of an ion as a result of collisional excitation during ion transfer from an atmospheric pressure ion source and the mass spectrometer vacuum."]
        InSourceCollisionInducedDissociation,
        #[term(cv=MS, accession=1002000, name="LIFT", flags={0}, parents={["MS:1000044"]})]
        #[doc = "LIFT - A Bruker's proprietary technique where molecular ions are initially accelerated at lower energy, then collide with inert gas in a collision cell that is then 'lifted' to high potential. The use of inert gas is optional, as it could lift also fragments provided by LID."]
        LIFT,
        #[term(cv=MS, accession=1002472, name="trap-type collision-induced dissociation", flags={0}, parents={["MS:1000133"]})]
        #[doc = "trap-type collision-induced dissociation - A collision-induced dissociation process that occurs in a trap-type collision cell."]
        TrapTypeCollisionInducedDissociation,
        #[term(cv=MS, accession=1002481, name="higher energy beam-type collision-induced dissociation", flags={0}, parents={["MS:1000422"]})]
        #[doc = "higher energy beam-type collision-induced dissociation - A collision-induced dissociation process wherein the projectile ion has the translational energy higher than approximately 1000 eV."]
        HigherEnergyBeamTypeCollisionInducedDissociation,
        #[term(cv=MS, accession=1002678, name="supplemental beam-type collision-induced dissociation", flags={0}, parents={["MS:1000422"]})]
        #[doc = "supplemental beam-type collision-induced dissociation - A supplemental collision-induced dissociation process that occurs in a beam-type collision cell in addition to another primary type of dissociation."]
        SupplementalBeamTypeCollisionInducedDissociation,
        #[term(cv=MS, accession=1002679, name="supplemental collision-induced dissociation", flags={0}, parents={["MS:1000133"]})]
        #[doc = "supplemental collision-induced dissociation - The dissociation of an ion after supplemental collisional excitation."]
        SupplementalCollisionInducedDissociation,
        #[term(cv=MS, accession=1003246, name="ultraviolet photodissociation", flags={0}, parents={["MS:1000435"]})]
        #[doc = "ultraviolet photodissociation - Multiphoton ionization where the reactant ion dissociates as a result of the absorption of multiple UV photons."]
        UltravioletPhotodissociation,
        #[term(cv=MS, accession=1003247, name="negative electron transfer dissociation", flags={0}, parents={["MS:1000044"]})]
        #[doc = "negative electron transfer dissociation - A process to fragment ions in a mass spectrometer by inducing fragmentation of anions (e.g. peptides or proteins) by transferring electrons to a radical-cation."]
        NegativeElectronTransferDissociation,
        #[term(cv=MS, accession=1003294, name="electron activated dissociation", flags={0}, parents={["MS:1000250"]})]
        #[doc = "electron activated dissociation - A process to fragment ions in a high intensity electron beam which results in a dissociation of various analytes ranging from singly charged small molecules to multiply protonated proteins."]
        ElectronActivatedDissociation,
    }
    //[[[end]]] (checksum: c0f6c7cfcb5d43054288f14702549985)
}

impl DissociationMethodTerm {

    pub fn is_electronic(&self) -> bool {
        match self {
            Self::ElectronActivatedDissociation
            | Self::ElectronCaptureDissociation
            | Self::ElectronTransferDissociation
            | Self::NegativeElectronTransferDissociation => true,
            _ => false
        }
    }

    pub fn is_collisional(&self) -> bool {
        match self {
            Self::CollisionInducedDissociation
            | Self::LowEnergyCollisionInducedDissociation
            | Self::BeamTypeCollisionInducedDissociation
            | Self::TrapTypeCollisionInducedDissociation
            | Self::InSourceCollisionInducedDissociation
            | Self::SupplementalBeamTypeCollisionInducedDissociation
            | Self::SupplementalCollisionInducedDissociation
            | Self::HigherEnergyBeamTypeCollisionInducedDissociation => true,
            _ => false,
        }
    }
}

crate::cvmap! {
    #[value_type=f32]
    #[flag_type=i32]
    #[doc = "Energy for an ion experiencing some form of dissociation"]
    #[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_energy.py']).decode('utf8')
    for line in buf.splitlines():
        cog.outl(line)
    ]]]*/
    pub enum DissociationEnergyTerm {
        #[term(cv=MS, accession=1000045, name="collision energy", flags={0}, parents={["MS:1000510"]})]
        #[doc = "collision energy - Energy for an ion experiencing collision with a stationary gas particle resulting in dissociation of the ion."]
        CollisionEnergy(f32),
        #[term(cv=MS, accession=1000138, name="normalized collision energy", flags={0}, parents={["MS:1000510"]})]
        #[doc = "normalized collision energy - Instrument setting, expressed in percent, for adjusting collisional energies of ions in an effort to provide equivalent excitation of all ions."]
        NormalizedCollisionEnergy(f32),
        #[term(cv=MS, accession=1002013, name="collision energy ramp start", flags={0}, parents={["MS:1000045"]})]
        #[doc = "collision energy ramp start - Collision energy at the start of the collision energy ramp."]
        CollisionEnergyRampStart(f32),
        #[term(cv=MS, accession=1002014, name="collision energy ramp end", flags={0}, parents={["MS:1000045"]})]
        #[doc = "collision energy ramp end - Collision energy at the end of the collision energy ramp."]
        CollisionEnergyRampEnd(f32),
        #[term(cv=MS, accession=1002218, name="percent collision energy ramp start", flags={0}, parents={["MS:1000138"]})]
        #[doc = "percent collision energy ramp start - Collision energy at the start of the collision energy ramp in percent, normalized to the mass of the ion."]
        PercentCollisionEnergyRampStart(f32),
        #[term(cv=MS, accession=1002219, name="percent collision energy ramp end", flags={0}, parents={["MS:1000138"]})]
        #[doc = "percent collision energy ramp end - Collision energy at the end of the collision energy ramp in percent, normalized to the mass of the ion."]
        PercentCollisionEnergyRampEnd(f32),
        #[term(cv=MS, accession=1002680, name="supplemental collision energy", flags={0}, parents={["MS:1000510"]})]
        #[doc = "supplemental collision energy - Energy for an ion experiencing supplemental collision with a stationary gas particle resulting in dissociation of the ion."]
        SupplementalCollisionEnergy(f32),
        #[term(cv=MS, accession=1003410, name="electron beam energy", flags={0}, parents={["MS:1000510"]})]
        #[doc = "electron beam energy - The kinetic energy of the electron beam used in dissociation methods induced by a free electron beam, such as electron-capture dissociation (ECD), electron-detachment dissociation (EDD), and electron-activated dissociation (EAD)."]
        ElectronBeamEnergy(f32),
    }
    //[[[end]]] (checksum: 252597ff6067b90ae0bc37f63da81021)
}

impl Display for DissociationEnergyTerm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl DissociationEnergyTerm {
    pub const fn is_supplemental(&self) -> bool {
        matches!(self, Self::SupplementalCollisionEnergy(_))
    }

    pub const fn is_ramp_start(&self) -> bool {
        matches!(self, Self::CollisionEnergyRampStart(_) | Self::PercentCollisionEnergyRampStart(_))
    }

    pub const fn is_ramp_end(&self) -> bool {
        matches!(self, Self::CollisionEnergyRampEnd(_) | Self::PercentCollisionEnergyRampEnd(_))
    }

    pub const fn energy(&self) -> f32 {
        *match self {
            DissociationEnergyTerm::CollisionEnergy(x) => x,
            DissociationEnergyTerm::NormalizedCollisionEnergy(x) => x,
            DissociationEnergyTerm::CollisionEnergyRampStart(x) => x,
            DissociationEnergyTerm::CollisionEnergyRampEnd(x) => x,
            DissociationEnergyTerm::PercentCollisionEnergyRampStart(x) => x,
            DissociationEnergyTerm::PercentCollisionEnergyRampEnd(x) => x,
            DissociationEnergyTerm::SupplementalCollisionEnergy(x) => x,
            DissociationEnergyTerm::ElectronBeamEnergy(x) => x,
        }
    }

    pub fn energy_mut(&mut self) -> &mut f32 {
        match self {
            DissociationEnergyTerm::CollisionEnergy(x) => x,
            DissociationEnergyTerm::NormalizedCollisionEnergy(x) => x,
            DissociationEnergyTerm::CollisionEnergyRampStart(x) => x,
            DissociationEnergyTerm::CollisionEnergyRampEnd(x) => x,
            DissociationEnergyTerm::PercentCollisionEnergyRampStart(x) => x,
            DissociationEnergyTerm::PercentCollisionEnergyRampEnd(x) => x,
            DissociationEnergyTerm::SupplementalCollisionEnergy(x) => x,
            DissociationEnergyTerm::ElectronBeamEnergy(x) => x,
        }
    }
}

impl AsRef<f32> for DissociationEnergyTerm {
    fn as_ref(&self) -> &f32 {
        match self {
            DissociationEnergyTerm::CollisionEnergy(x) => x,
            DissociationEnergyTerm::NormalizedCollisionEnergy(x) => x,
            DissociationEnergyTerm::CollisionEnergyRampStart(x) => x,
            DissociationEnergyTerm::CollisionEnergyRampEnd(x) => x,
            DissociationEnergyTerm::PercentCollisionEnergyRampStart(x) => x,
            DissociationEnergyTerm::PercentCollisionEnergyRampEnd(x) => x,
            DissociationEnergyTerm::SupplementalCollisionEnergy(x) => x,
            DissociationEnergyTerm::ElectronBeamEnergy(x) => x,
        }
    }
}

#[allow(unused)]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum DissociationEnergy {
    Energy(DissociationEnergyTerm),
    Ramp { start: DissociationEnergyTerm, end: DissociationEnergyTerm },
    Combined { primary: DissociationEnergyTerm, supplementary: DissociationEnergyTerm}
}

#[cfg(test)]
mod test {
    use crate::{params::{ControlledVocabulary, ParamCow, ValueRef}, Param};

    use super::*;

    #[test]
    fn test_categories() {
        assert!(DissociationMethodTerm::CollisionInducedDissociation.is_collisional());
        assert!(DissociationMethodTerm::BeamTypeCollisionInducedDissociation.is_collisional());
        assert!(!DissociationMethodTerm::BeamTypeCollisionInducedDissociation.is_electronic());
        assert!(!DissociationMethodTerm::ElectronTransferDissociation.is_collisional());
        assert!(DissociationMethodTerm::ElectronTransferDissociation.is_electronic());

        assert_eq!(DissociationEnergyTerm::CollisionEnergy(30.0).to_string(), "CollisionEnergy(30.0)");

        assert!(!DissociationEnergyTerm::CollisionEnergy(30.0).is_ramp_end());
        assert!(!DissociationEnergyTerm::NormalizedCollisionEnergy(30.0).is_ramp_start());
        assert!(*DissociationEnergyTerm::NormalizedCollisionEnergy(30.0).as_ref() == 30.0);

        let mut energies = [
            DissociationEnergyTerm::CollisionEnergy(10.0),
            DissociationEnergyTerm::NormalizedCollisionEnergy(10.0),
            DissociationEnergyTerm::CollisionEnergyRampStart(10.0),
            DissociationEnergyTerm::CollisionEnergyRampEnd(10.0),
            DissociationEnergyTerm::PercentCollisionEnergyRampStart(10.0),
            DissociationEnergyTerm::PercentCollisionEnergyRampEnd(10.0),
            DissociationEnergyTerm::SupplementalCollisionEnergy(10.0),
            DissociationEnergyTerm::ElectronBeamEnergy(10.0),
        ];
        for e in energies.iter_mut() {
            assert_eq!(e.energy(), 10.0);
            assert_eq!(*e.energy_mut(), 10.0);
            assert_eq!(*e.as_ref(), 10.0);
        }
    }

    #[test]
    fn test_meta() {
        // #[term(cv=MS, accession=1000138, name="normalized collision energy", flags={0}, parents={["MS:1000510"]})]
        // #[doc = "normalized collision energy - Instrument setting, expressed in percent, for adjusting collisional energies of ions in an effort to provide equivalent excitation of all ions."]
        // NormalizedCollisionEnergy(f32),

        assert_eq!(
            DissociationEnergyTerm::from_accession(1000138, 0.0),
            Some(DissociationEnergyTerm::NormalizedCollisionEnergy(0.0))
        );

        let mut param = DissociationEnergyTerm::NormalizedCollisionEnergy(30.0).to_param();
        param.value = 30.0.into();
        assert_eq!(
            DissociationEnergyTerm::from_param(&param),
            Some(DissociationEnergyTerm::NormalizedCollisionEnergy(30.0))
        );

        let term: DissociationEnergyTerm = param.clone().into();
        assert_eq!(
            term,
            DissociationEnergyTerm::NormalizedCollisionEnergy(30.0)
        );

        let p: Param = DissociationEnergyTerm::NormalizedCollisionEnergy(30.0).into();
        param.value = ValueRef::Empty;
        assert_eq!(
            p,
            param
        );

        let p: Param = (&DissociationEnergyTerm::NormalizedCollisionEnergy(30.0)).into();
        assert_eq!(
            p,
            param
        );

        let p: ParamCow = DissociationEnergyTerm::NormalizedCollisionEnergy(30.0).into();
        assert_eq!(
            p,
            param
        );

        let p: ParamCow = (&DissociationEnergyTerm::NormalizedCollisionEnergy(30.0)).into();
        assert_eq!(
            p,
            param
        );

        param.value = 30.0.into();
        let term: DissociationEnergyTerm = Param::from(param.clone()).into();
        assert_eq!(
            term,
            DissociationEnergyTerm::NormalizedCollisionEnergy(30.0)
        );

        assert_eq!(
            DissociationEnergyTerm::NormalizedCollisionEnergy(30.0).accession(),
            1000138
        );

        assert_eq!(
            DissociationEnergyTerm::NormalizedCollisionEnergy(30.0).parents(),
            []
        );

        assert_eq!(
            DissociationEnergyTerm::NormalizedCollisionEnergy(30.0).flags(),
            0,
        );

        assert_eq!(
            DissociationEnergyTerm::NormalizedCollisionEnergy(30.0).controlled_vocabulary(),
            ControlledVocabulary::MS,
        );

        assert_eq!(
            DissociationEnergyTerm::NormalizedCollisionEnergy(30.0).name(),
            "normalized collision energy"
        );

        assert_eq!(
            DissociationEnergyTerm::from_name("normalized collision energy"),
            Some(DissociationEnergyTerm::NormalizedCollisionEnergy(Default::default()))
        );
    }
}