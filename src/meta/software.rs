use std::collections::HashSet;

use crate::params::{ControlledVocabulary, ParamList};
use crate::{impl_param_described, Param};

/// A piece of software that was associated with the acquisition, transformation or otherwise
/// processing of mass spectrometry data.
///
/// See <https://peptideatlas.org/tmp/mzML1.1.0.html#software>
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Software {
    /// A unique identifier for the software within processing metadata
    pub id: String,
    /// A string denoting a particular software version, but does no guarantee is given for its format
    pub version: String,
    /// Any associated vocabulary terms, including actual software name and type
    pub params: ParamList,
}

bitflags::bitflags! {
    #[doc="A bit mask encoding the different kinds of software."]
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct SoftwareType: u8 {
        const None = 0;
        const Analysis = 0b00000001;
        const DataProcessing = 0b00000010;
        const Acquisition = 0b00000100;
    }
}

#[allow(unused)]
impl SoftwareType {
    /// Is this software for analysis?
    pub const fn is_analysis(&self) -> bool {
        self.contains(Self::Analysis)
    }

    /// Is this software for data processing?
    pub const fn is_data_processing(&self) -> bool {
        self.contains(Self::DataProcessing)
    }

    /// Is this software for data acquisition?
    pub const fn is_acquisition(&self) -> bool {
        self.contains(Self::Acquisition)
    }
}

impl From<u8> for SoftwareType {
    fn from(value: u8) -> Self {
        Self::from_bits_retain(value)
    }
}

impl Software {
    pub fn new(id: String, version: String, params: ParamList) -> Self {
        Self {
            id,
            version,
            params,
        }
    }

    /// Find the term defining the software
    pub fn find_software_term(&self) -> Option<SoftwareTerm> {
        self.params
            .iter()
            .flat_map(|p| {
                if let Some(i) = p.accession {
                    SoftwareTerm::from_accession(i)
                } else {
                    None
                }
            })
            .next()
    }

    /// Is this software for analysis?
    pub fn is_analysis(&self) -> bool {
        self.find_software_term().map(|s| s.flags().is_analysis()).unwrap_or(false)
    }

    /// Is this software for data processing?
    pub fn is_data_processing(&self) -> bool {
        self.find_software_term().map(|s| s.flags().is_data_processing()).unwrap_or(false)
    }

    /// Is this software for data acquisition?
    pub fn is_acquisition(&self) -> bool {
        self.find_software_term().map(|s| s.flags().is_acquisition()).unwrap_or(false)
    }

    /// Find a unique identifier from an iterator over software IDs
    pub fn find_unique_id<'a>(
        id_stem: &str,
        softwares: impl IntoIterator<Item = &'a Self>,
    ) -> String {
        let software_ids: HashSet<_> = softwares.into_iter().map(|sw| &sw.id).collect();
        (0..)
            .map(|i| format!("{id_stem}_{i}"))
            .find(|s| !software_ids.contains(s))
            .unwrap()
    }
}

/// Create an instance of "custom unreleased software tool" with name `name`
pub fn custom_software_name(name: &str) -> Param {
    ControlledVocabulary::MS.param_val(1000799, "custom unreleased software tool", name)
}

impl_param_described!(Software);

crate::cvmap! {
    #[flag_type=SoftwareType]
    #[allow(unused, clippy::upper_case_acronyms)]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_software.py']).decode('utf8')
    for line in buf.splitlines():
        cog.outl(line)
    ]]]*/
    pub enum SoftwareTerm {
        #[term(cv=MS, accession=1000531, name="software", flags={0}, parents={[]})]
        Software,
        #[term(cv=MS, accession=1000532, name="Xcalibur", flags={7}, parents={["MS:1000693", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        Xcalibur,
        #[term(cv=MS, accession=1000533, name="Bioworks", flags={3}, parents={["MS:1000693", "MS:1001456", "MS:1001457"]})]
        Bioworks,
        #[term(cv=MS, accession=1000534, name="MassLynx", flags={7}, parents={["MS:1000694", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        MassLynx,
        #[term(cv=MS, accession=1000535, name="FlexAnalysis", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        FlexAnalysis,
        #[term(cv=MS, accession=1000536, name="Data Explorer", flags={7}, parents={["MS:1000691", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        DataExplorer,
        #[term(cv=MS, accession=1000537, name="4700 Explorer", flags={7}, parents={["MS:1000691", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        _4700Explorer,
        #[term(cv=MS, accession=1000538, name="massWolf", flags={2}, parents={["MS:1001457"]})]
        MassWolf,
        #[term(cv=MS, accession=1000539, name="Voyager Biospectrometry Workstation System", flags={7}, parents={["MS:1000691", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        VoyagerBiospectrometryWorkstationSystem,
        #[term(cv=MS, accession=1000540, name="FlexControl", flags={4}, parents={["MS:1000692", "MS:1001455"]})]
        FlexControl,
        #[term(cv=MS, accession=1000541, name="ReAdW", flags={2}, parents={["MS:1001457"]})]
        ReAdW,
        #[term(cv=MS, accession=1000542, name="MzStar", flags={2}, parents={["MS:1001457"]})]
        MzStar,
        #[term(cv=MS, accession=1000551, name="Analyst", flags={7}, parents={["MS:1000690", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        Analyst,
        #[term(cv=MS, accession=1000553, name="Trapper", flags={2}, parents={["MS:1001457"]})]
        Trapper,
        #[term(cv=MS, accession=1000591, name="MzWiff", flags={2}, parents={["MS:1001457"]})]
        MzWiff,
        #[term(cv=MS, accession=1000600, name="Proteios", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        Proteios,
        #[term(cv=MS, accession=1000601, name="ProteinLynx Global Server", flags={3}, parents={["MS:1000694", "MS:1001456", "MS:1001457"]})]
        ProteinLynxGlobalServer,
        #[term(cv=MS, accession=1000615, name="ProteoWizard software", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        ProteoWizardSoftware,
        #[term(cv=MS, accession=1000650, name="Proteome Discoverer", flags={3}, parents={["MS:1000693", "MS:1001456", "MS:1001457"]})]
        ProteomeDiscoverer,
        #[term(cv=MS, accession=1000659, name="4000 Series Explorer Software", flags={7}, parents={["MS:1000690", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        _4000SeriesExplorerSoftware,
        #[term(cv=MS, accession=1000661, name="GPS Explorer", flags={7}, parents={["MS:1000690", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        GPSExplorer,
        #[term(cv=MS, accession=1000662, name="LightSight Software", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        LightSightSoftware,
        #[term(cv=MS, accession=1000663, name="ProteinPilot Software", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        ProteinPilotSoftware,
        #[term(cv=MS, accession=1000664, name="TissueView Software", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        TissueViewSoftware,
        #[term(cv=MS, accession=1000665, name="MarkerView Software", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        MarkerViewSoftware,
        #[term(cv=MS, accession=1000666, name="MRMPilot Software", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        MRMPilotSoftware,
        #[term(cv=MS, accession=1000667, name="BioAnalyst", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        BioAnalyst,
        #[term(cv=MS, accession=1000668, name="Pro ID", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        ProID,
        #[term(cv=MS, accession=1000669, name="Pro ICAT", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        ProICAT,
        #[term(cv=MS, accession=1000670, name="Pro Quant", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        ProQuant,
        #[term(cv=MS, accession=1000671, name="Pro BLAST", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        ProBLAST,
        #[term(cv=MS, accession=1000672, name="Cliquid", flags={0}, parents={["MS:1000690"]})]
        Cliquid,
        #[term(cv=MS, accession=1000673, name="MIDAS Workflow Designer", flags={0}, parents={["MS:1000690"]})]
        MIDASWorkflowDesigner,
        #[term(cv=MS, accession=1000674, name="MultiQuant", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        MultiQuant,
        #[term(cv=MS, accession=1000678, name="MassHunter Data Acquisition", flags={4}, parents={["MS:1000689", "MS:1001455"]})]
        MassHunterDataAcquisition,
        #[term(cv=MS, accession=1000679, name="MassHunter Easy Access", flags={4}, parents={["MS:1000689", "MS:1001455"]})]
        MassHunterEasyAccess,
        #[term(cv=MS, accession=1000680, name="MassHunter Qualitative Analysis", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        MassHunterQualitativeAnalysis,
        #[term(cv=MS, accession=1000681, name="MassHunter Quantitative Analysis", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        MassHunterQuantitativeAnalysis,
        #[term(cv=MS, accession=1000682, name="MassHunter Metabolite ID", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        MassHunterMetaboliteID,
        #[term(cv=MS, accession=1000683, name="MassHunter BioConfirm", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        MassHunterBioConfirm,
        #[term(cv=MS, accession=1000684, name="Genespring MS", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        GenespringMS,
        #[term(cv=MS, accession=1000685, name="MassHunter Mass Profiler", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        MassHunterMassProfiler,
        #[term(cv=MS, accession=1000686, name="METLIN", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        METLIN,
        #[term(cv=MS, accession=1000687, name="Spectrum Mill for MassHunter Workstation", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        SpectrumMillForMassHunterWorkstation,
        #[term(cv=MS, accession=1000688, name="6300 Series Ion Trap Data Analysis Software", flags={7}, parents={["MS:1000689", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        _6300SeriesIonTrapDataAnalysisSoftware,
        #[term(cv=MS, accession=1000689, name="Agilent software", flags={0}, parents={["MS:1000531"]})]
        AgilentSoftware,
        #[term(cv=MS, accession=1000690, name="SCIEX software", flags={0}, parents={["MS:1000531"]})]
        SCIEXSoftware,
        #[term(cv=MS, accession=1000691, name="Applied Biosystems software", flags={0}, parents={["MS:1000531"]})]
        AppliedBiosystemsSoftware,
        #[term(cv=MS, accession=1000692, name="Bruker software", flags={0}, parents={["MS:1000531"]})]
        BrukerSoftware,
        #[term(cv=MS, accession=1000693, name="Thermo Finnigan software", flags={0}, parents={["MS:1000531"]})]
        ThermoFinniganSoftware,
        #[term(cv=MS, accession=1000694, name="Waters software", flags={0}, parents={["MS:1000531"]})]
        WatersSoftware,
        #[term(cv=MS, accession=1000706, name="apexControl", flags={4}, parents={["MS:1000692", "MS:1001455"]})]
        ApexControl,
        #[term(cv=MS, accession=1000707, name="BioTools", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        BioTools,
        #[term(cv=MS, accession=1000708, name="CLINPROT", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        CLINPROT,
        #[term(cv=MS, accession=1000709, name="CLINPROT micro", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        CLINPROTMicro,
        #[term(cv=MS, accession=1000710, name="CLINPROT robot", flags={0}, parents={["MS:1000692"]})]
        CLINPROTRobot,
        #[term(cv=MS, accession=1000711, name="ClinProTools", flags={0}, parents={["MS:1000692"]})]
        ClinProTools,
        #[term(cv=MS, accession=1000712, name="Compass", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        Compass,
        #[term(cv=MS, accession=1000713, name="Compass for HCT/esquire", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        CompassForHCTEsquire,
        #[term(cv=MS, accession=1000714, name="Compass for micrOTOF", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        CompassForMicrOTOF,
        #[term(cv=MS, accession=1000715, name="Compass OpenAccess", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        CompassOpenAccess,
        #[term(cv=MS, accession=1000716, name="Compass Security Pack", flags={0}, parents={["MS:1000692"]})]
        CompassSecurityPack,
        #[term(cv=MS, accession=1000717, name="CompassXport", flags={2}, parents={["MS:1000692", "MS:1001457"]})]
        CompassXport,
        #[term(cv=MS, accession=1000718, name="CompassXtract", flags={2}, parents={["MS:1000692", "MS:1001457"]})]
        CompassXtract,
        #[term(cv=MS, accession=1000719, name="DataAnalysis", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        DataAnalysis,
        #[term(cv=MS, accession=1000720, name="dpControl", flags={4}, parents={["MS:1000692", "MS:1001455"]})]
        DpControl,
        #[term(cv=MS, accession=1000721, name="esquireControl", flags={4}, parents={["MS:1000692", "MS:1001455"]})]
        EsquireControl,
        #[term(cv=MS, accession=1000722, name="flexImaging", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        FlexImaging,
        #[term(cv=MS, accession=1000723, name="GENOLINK", flags={0}, parents={["MS:1000692"]})]
        GENOLINK,
        #[term(cv=MS, accession=1000724, name="GenoTools", flags={0}, parents={["MS:1000692"]})]
        GenoTools,
        #[term(cv=MS, accession=1000725, name="HCTcontrol", flags={4}, parents={["MS:1000692", "MS:1001455"]})]
        HCTcontrol,
        #[term(cv=MS, accession=1000726, name="micrOTOFcontrol", flags={4}, parents={["MS:1000692", "MS:1001455"]})]
        MicrOTOFcontrol,
        #[term(cv=MS, accession=1000727, name="PolyTools", flags={0}, parents={["MS:1000692"]})]
        PolyTools,
        #[term(cv=MS, accession=1000728, name="ProfileAnalysis", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        ProfileAnalysis,
        #[term(cv=MS, accession=1000729, name="PROTEINEER", flags={0}, parents={["MS:1000692"]})]
        PROTEINEER,
        #[term(cv=MS, accession=1000730, name="PROTEINEER dp", flags={0}, parents={["MS:1000692"]})]
        PROTEINEERDp,
        #[term(cv=MS, accession=1000731, name="PROTEINEER fc", flags={0}, parents={["MS:1000692"]})]
        PROTEINEERFc,
        #[term(cv=MS, accession=1000732, name="PROTEINEER spII", flags={0}, parents={["MS:1000692"]})]
        PROTEINEERSpII,
        #[term(cv=MS, accession=1000733, name="PROTEINEER-LC", flags={0}, parents={["MS:1000692"]})]
        PROTEINEERLC,
        #[term(cv=MS, accession=1000734, name="ProteinScape", flags={1}, parents={["MS:1000692", "MS:1001456"]})]
        ProteinScape,
        #[term(cv=MS, accession=1000735, name="PureDisk", flags={0}, parents={["MS:1000692"]})]
        PureDisk,
        #[term(cv=MS, accession=1000736, name="QuantAnalysis", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        QuantAnalysis,
        #[term(cv=MS, accession=1000737, name="spControl", flags={4}, parents={["MS:1000692", "MS:1001455"]})]
        SpControl,
        #[term(cv=MS, accession=1000738, name="TargetAnalysis", flags={0}, parents={["MS:1000692"]})]
        TargetAnalysis,
        #[term(cv=MS, accession=1000739, name="WARP-LC", flags={0}, parents={["MS:1000692", "MS:1001139"]})]
        WARPLC,
        #[term(cv=MS, accession=1000752, name="TOPP software", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        TOPPSoftware,
        #[term(cv=MS, accession=1000753, name="BaselineFilter", flags={0}, parents={["MS:1000752"]})]
        BaselineFilter,
        #[term(cv=MS, accession=1000754, name="DBExporter", flags={0}, parents={["MS:1000752"]})]
        DBExporter,
        #[term(cv=MS, accession=1000755, name="DBImporter", flags={0}, parents={["MS:1000752"]})]
        DBImporter,
        #[term(cv=MS, accession=1000756, name="FileConverter", flags={0}, parents={["MS:1000752"]})]
        FileConverter,
        #[term(cv=MS, accession=1000757, name="FileFilter", flags={0}, parents={["MS:1000752"]})]
        FileFilter,
        #[term(cv=MS, accession=1000758, name="FileMerger", flags={0}, parents={["MS:1000752"]})]
        FileMerger,
        #[term(cv=MS, accession=1000759, name="InternalCalibration", flags={0}, parents={["MS:1000752"]})]
        InternalCalibration,
        #[term(cv=MS, accession=1000760, name="MapAligner", flags={0}, parents={["MS:1000752"]})]
        MapAligner,
        #[term(cv=MS, accession=1000761, name="MapNormalizer", flags={0}, parents={["MS:1000752"]})]
        MapNormalizer,
        #[term(cv=MS, accession=1000762, name="NoiseFilter", flags={0}, parents={["MS:1000752"]})]
        NoiseFilter,
        #[term(cv=MS, accession=1000763, name="PeakPicker", flags={0}, parents={["MS:1000752"]})]
        PeakPicker,
        #[term(cv=MS, accession=1000764, name="Resampler", flags={0}, parents={["MS:1000752"]})]
        Resampler,
        #[term(cv=MS, accession=1000765, name="SpectraFilter", flags={0}, parents={["MS:1000752"]})]
        SpectraFilter,
        #[term(cv=MS, accession=1000766, name="TOFCalibration", flags={0}, parents={["MS:1000752"]})]
        TOFCalibration,
        #[term(cv=MS, accession=1000799, name="custom unreleased software tool", flags={0}, parents={["MS:1000531"]})]
        CustomUnreleasedSoftwareTool,
        #[term(cv=MS, accession=1000817, name="HyStar", flags={0}, parents={["MS:1000692"]})]
        HyStar,
        #[term(cv=MS, accession=1000871, name="SRM software", flags={0}, parents={["MS:1000531"]})]
        SRMSoftware,
        #[term(cv=MS, accession=1000872, name="MaRiMba", flags={0}, parents={["MS:1000871"]})]
        MaRiMba,
        #[term(cv=MS, accession=1000873, name="peptide attribute calculation software", flags={0}, parents={["MS:1000531"]})]
        PeptideAttributeCalculationSoftware,
        #[term(cv=MS, accession=1000874, name="SSRCalc", flags={0}, parents={["MS:1000873"]})]
        SSRCalc,
        #[term(cv=MS, accession=1000922, name="Skyline", flags={0}, parents={["MS:1000871", "MS:1001139"]})]
        Skyline,
        #[term(cv=MS, accession=1000923, name="TIQAM", flags={0}, parents={["MS:1000871"]})]
        TIQAM,
        #[term(cv=MS, accession=1000925, name="ATAQS", flags={0}, parents={["MS:1000871"]})]
        ATAQS,
        #[term(cv=MS, accession=1001139, name="quantitation software name", flags={0}, parents={["MS:1000531", "MS:1001129"]})]
        QuantitationSoftwareName,
        #[term(cv=MS, accession=1001207, name="Mascot", flags={1}, parents={["MS:1001456"]})]
        Mascot,
        #[term(cv=MS, accession=1001208, name="SEQUEST", flags={1}, parents={["MS:1001456"]})]
        SEQUEST,
        #[term(cv=MS, accession=1001209, name="Phenyx", flags={1}, parents={["MS:1001456"]})]
        Phenyx,
        #[term(cv=MS, accession=1001327, name="Spectronaut", flags={1}, parents={["MS:1001456", "MS:1003207"]})]
        Spectronaut,
        #[term(cv=MS, accession=1001455, name="acquisition software", flags={0}, parents={["MS:1000531"]})]
        AcquisitionSoftware,
        #[term(cv=MS, accession=1001456, name="analysis software", flags={0}, parents={["MS:1000531"]})]
        AnalysisSoftware,
        #[term(cv=MS, accession=1001457, name="data processing software", flags={0}, parents={["MS:1000531"]})]
        DataProcessingSoftware,
        #[term(cv=MS, accession=1001461, name="greylag", flags={1}, parents={["MS:1001456"]})]
        Greylag,
        #[term(cv=MS, accession=1001475, name="OMSSA", flags={1}, parents={["MS:1001456"]})]
        OMSSA,
        #[term(cv=MS, accession=1001476, name="X!Tandem", flags={1}, parents={["MS:1001456"]})]
        XTandem,
        #[term(cv=MS, accession=1001477, name="SpectraST", flags={1}, parents={["MS:1001456", "MS:1003207", "MS:1003406"]})]
        SpectraST,
        #[term(cv=MS, accession=1001478, name="Mascot Parser", flags={1}, parents={["MS:1001456"]})]
        MascotParser,
        #[term(cv=MS, accession=1001483, name="SCIEX TOF/TOF Series Explorer Software", flags={7}, parents={["MS:1000690", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        SCIEXTOFTOFSeriesExplorerSoftware,
        #[term(cv=MS, accession=1001487, name="ProteinExtractor", flags={1}, parents={["MS:1000692", "MS:1001456"]})]
        ProteinExtractor,
        #[term(cv=MS, accession=1001488, name="Mascot Distiller", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        MascotDistiller,
        #[term(cv=MS, accession=1001489, name="Mascot Integra", flags={1}, parents={["MS:1001456"]})]
        MascotIntegra,
        #[term(cv=MS, accession=1001490, name="Percolator", flags={1}, parents={["MS:1001456"]})]
        Percolator,
        #[term(cv=MS, accession=1001557, name="Shimadzu Corporation software", flags={0}, parents={["MS:1000531"]})]
        ShimadzuCorporationSoftware,
        #[term(cv=MS, accession=1001558, name="MALDI Solutions", flags={7}, parents={["MS:1001455", "MS:1001456", "MS:1001457", "MS:1001557"]})]
        MALDISolutions,
        #[term(cv=MS, accession=1001561, name="Scaffold", flags={1}, parents={["MS:1001456"]})]
        Scaffold,
        #[term(cv=MS, accession=1001582, name="XCMS", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        XCMS,
        #[term(cv=MS, accession=1001583, name="MaxQuant", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        MaxQuant,
        #[term(cv=MS, accession=1001585, name="MyriMatch", flags={1}, parents={["MS:1001456"]})]
        MyriMatch,
        #[term(cv=MS, accession=1001586, name="DirecTag", flags={1}, parents={["MS:1001456"]})]
        DirecTag,
        #[term(cv=MS, accession=1001587, name="TagRecon", flags={1}, parents={["MS:1001456"]})]
        TagRecon,
        #[term(cv=MS, accession=1001588, name="Pepitome", flags={1}, parents={["MS:1001456"]})]
        Pepitome,
        #[term(cv=MS, accession=1001795, name="Empower", flags={3}, parents={["MS:1000694", "MS:1001456", "MS:1001457"]})]
        Empower,
        #[term(cv=MS, accession=1001796, name="UNIFY", flags={3}, parents={["MS:1000694", "MS:1001456", "MS:1001457"]})]
        UNIFY,
        #[term(cv=MS, accession=1001798, name="LECO software", flags={0}, parents={["MS:1000531"]})]
        LECOSoftware,
        #[term(cv=MS, accession=1001799, name="ChromaTOF software", flags={7}, parents={["MS:1001455", "MS:1001456", "MS:1001457", "MS:1001798"]})]
        ChromaTOFSoftware,
        #[term(cv=MS, accession=1001830, name="Progenesis LC-MS", flags={0}, parents={["MS:1001139"]})]
        ProgenesisLCMS,
        #[term(cv=MS, accession=1001831, name="SILACAnalyzer", flags={0}, parents={["MS:1001139", "MS:1000752"]})]
        SILACAnalyzer,
        #[term(cv=MS, accession=1001877, name="ChromaTOF HRT software", flags={7}, parents={["MS:1001455", "MS:1001456", "MS:1001457", "MS:1001798"]})]
        ChromaTOFHRTSoftware,
        #[term(cv=MS, accession=1001878, name="MALDI Solutions Microbial Identification", flags={0}, parents={["MS:1001558"]})]
        MALDISolutionsMicrobialIdentification,
        #[term(cv=MS, accession=1001886, name="SQID", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        SQID,
        #[term(cv=MS, accession=1001912, name="PinPoint", flags={3}, parents={["MS:1000693", "MS:1001456", "MS:1001457"]})]
        PinPoint,
        #[term(cv=MS, accession=1001914, name="pymzML", flags={2}, parents={["MS:1001457"]})]
        PymzML,
        #[term(cv=MS, accession=1001946, name="PEAKS Studio", flags={3}, parents={["MS:1001139", "MS:1001456", "MS:1001457"]})]
        PEAKSStudio,
        #[term(cv=MS, accession=1001947, name="PEAKS Online", flags={3}, parents={["MS:1001139", "MS:1001456", "MS:1001457"]})]
        PEAKSOnline,
        #[term(cv=MS, accession=1001948, name="PEAKS Node", flags={3}, parents={["MS:1001139", "MS:1001456", "MS:1001457"]})]
        PEAKSNode,
        #[term(cv=MS, accession=1001949, name="BSI software", flags={0}, parents={["MS:1000531"]})]
        BSISoftware,
        #[term(cv=MS, accession=1001973, name="DeBunker", flags={1}, parents={["MS:1001456"]})]
        DeBunker,
        #[term(cv=MS, accession=1001977, name="MSQuant", flags={1}, parents={["MS:1001456"]})]
        MSQuant,
        #[term(cv=MS, accession=1001984, name="Ascore software", flags={1}, parents={["MS:1001456"]})]
        AscoreSoftware,
        #[term(cv=MS, accession=1002043, name="ProteinProspector", flags={1}, parents={["MS:1001456"]})]
        ProteinProspector,
        #[term(cv=MS, accession=1002047, name="MS-GF", flags={1}, parents={["MS:1001456"]})]
        MSGF,
        #[term(cv=MS, accession=1002048, name="MS-GF+", flags={1}, parents={["MS:1001456"]})]
        MSGFplus,
        #[term(cv=MS, accession=1002059, name="Microsoft Excel", flags={0}, parents={["MS:1001139"]})]
        MicrosoftExcel,
        #[term(cv=MS, accession=1002063, name="FindPairs", flags={0}, parents={["MS:1001139"]})]
        FindPairs,
        #[term(cv=MS, accession=1002076, name="PAnalyzer", flags={1}, parents={["MS:1001456"]})]
        PAnalyzer,
        #[term(cv=MS, accession=1002123, name="x-Tracker", flags={0}, parents={["MS:1001139"]})]
        XTracker,
        #[term(cv=MS, accession=1002124, name="ProteoSuite", flags={0}, parents={["MS:1001139"]})]
        ProteoSuite,
        #[term(cv=MS, accession=1002129, name="ITRAQAnalyzer", flags={0}, parents={["MS:1001139", "MS:1000752"]})]
        ITRAQAnalyzer,
        #[term(cv=MS, accession=1002131, name="TOPP noise filter", flags={0}, parents={["MS:1000752"]})]
        TOPPNoiseFilter,
        #[term(cv=MS, accession=1002132, name="TOPP NoiseFilterGaussian", flags={0}, parents={["MS:1002131"]})]
        TOPPNoiseFilterGaussian,
        #[term(cv=MS, accession=1002133, name="TOPP NoiseFilterSGolay", flags={0}, parents={["MS:1002131"]})]
        TOPPNoiseFilterSGolay,
        #[term(cv=MS, accession=1002134, name="TOPP peak picker", flags={0}, parents={["MS:1000752"]})]
        TOPPPeakPicker,
        #[term(cv=MS, accession=1002135, name="TOPP PeakPickerHiRes", flags={0}, parents={["MS:1002134"]})]
        TOPPPeakPickerHiRes,
        #[term(cv=MS, accession=1002136, name="TOPP PeakPickerWavelet", flags={0}, parents={["MS:1002134"]})]
        TOPPPeakPickerWavelet,
        #[term(cv=MS, accession=1002137, name="TOPP spectra filter", flags={0}, parents={["MS:1000752"]})]
        TOPPSpectraFilter,
        #[term(cv=MS, accession=1002138, name="TOPP SpectraFilterBernNorm", flags={0}, parents={["MS:1002137"]})]
        TOPPSpectraFilterBernNorm,
        #[term(cv=MS, accession=1002139, name="TOPP SpectraFilterMarkerMower", flags={0}, parents={["MS:1002137"]})]
        TOPPSpectraFilterMarkerMower,
        #[term(cv=MS, accession=1002140, name="TOPP SpectraFilterNLargest", flags={0}, parents={["MS:1002137"]})]
        TOPPSpectraFilterNLargest,
        #[term(cv=MS, accession=1002141, name="TOPP SpectraFilterNormalizer", flags={0}, parents={["MS:1002137"]})]
        TOPPSpectraFilterNormalizer,
        #[term(cv=MS, accession=1002142, name="TOPP SpectraFilterParentPeakMower", flags={0}, parents={["MS:1002137"]})]
        TOPPSpectraFilterParentPeakMower,
        #[term(cv=MS, accession=1002143, name="TOPP SpectraFilterScaler", flags={0}, parents={["MS:1002137"]})]
        TOPPSpectraFilterScaler,
        #[term(cv=MS, accession=1002144, name="TOPP SpectraFilterSqrtMower", flags={0}, parents={["MS:1002137"]})]
        TOPPSpectraFilterSqrtMower,
        #[term(cv=MS, accession=1002145, name="TOPP SpectraFilterThresholdMower", flags={0}, parents={["MS:1002137"]})]
        TOPPSpectraFilterThresholdMower,
        #[term(cv=MS, accession=1002146, name="TOPP SpectraFilterWindowMower", flags={0}, parents={["MS:1002137"]})]
        TOPPSpectraFilterWindowMower,
        #[term(cv=MS, accession=1002147, name="TOPP map aligner", flags={0}, parents={["MS:1000752"]})]
        TOPPMapAligner,
        #[term(cv=MS, accession=1002148, name="TOPP MapAlignerIdentification", flags={0}, parents={["MS:1002147"]})]
        TOPPMapAlignerIdentification,
        #[term(cv=MS, accession=1002149, name="TOPP MapAlignerPoseClustering", flags={0}, parents={["MS:1002147"]})]
        TOPPMapAlignerPoseClustering,
        #[term(cv=MS, accession=1002150, name="TOPP MapAlignerSpectrum", flags={0}, parents={["MS:1002147"]})]
        TOPPMapAlignerSpectrum,
        #[term(cv=MS, accession=1002154, name="TOPP DTAExtractor", flags={0}, parents={["MS:1000752"]})]
        TOPPDTAExtractor,
        #[term(cv=MS, accession=1002155, name="TOPP IDMerger", flags={0}, parents={["MS:1000752"]})]
        TOPPIDMerger,
        #[term(cv=MS, accession=1002156, name="TOPP IDFileConverter", flags={0}, parents={["MS:1000752"]})]
        TOPPIDFileConverter,
        #[term(cv=MS, accession=1002157, name="TOPP SpectraMerger", flags={0}, parents={["MS:1000752"]})]
        TOPPSpectraMerger,
        #[term(cv=MS, accession=1002158, name="TOPP MzTabExporter", flags={0}, parents={["MS:1000752"]})]
        TOPPMzTabExporter,
        #[term(cv=MS, accession=1002159, name="TOPP MassTraceExtractor", flags={0}, parents={["MS:1000752"]})]
        TOPPMassTraceExtractor,
        #[term(cv=MS, accession=1002160, name="TOPP PrecursorMassCorrector", flags={0}, parents={["MS:1000752"]})]
        TOPPPrecursorMassCorrector,
        #[term(cv=MS, accession=1002161, name="TOPP HighResPrecursorMassCorrector", flags={0}, parents={["MS:1000752"]})]
        TOPPHighResPrecursorMassCorrector,
        #[term(cv=MS, accession=1002162, name="TOPP AdditiveSeries", flags={0}, parents={["MS:1000752"]})]
        TOPPAdditiveSeries,
        #[term(cv=MS, accession=1002163, name="TOPP Decharger", flags={0}, parents={["MS:1000752"]})]
        TOPPDecharger,
        #[term(cv=MS, accession=1002164, name="TOPP EICExtractor", flags={0}, parents={["MS:1000752"]})]
        TOPPEICExtractor,
        #[term(cv=MS, accession=1002165, name="TOPP feature finder", flags={0}, parents={["MS:1000752"]})]
        TOPPFeatureFinder,
        #[term(cv=MS, accession=1002166, name="TOPP FeatureFinderCentroided", flags={0}, parents={["MS:1002165"]})]
        TOPPFeatureFinderCentroided,
        #[term(cv=MS, accession=1002167, name="TOPP FeatureFinderRaw", flags={0}, parents={["MS:1002165"]})]
        TOPPFeatureFinderRaw,
        #[term(cv=MS, accession=1002168, name="TOPP FeatureFinderIsotopeWavelet", flags={0}, parents={["MS:1002165"]})]
        TOPPFeatureFinderIsotopeWavelet,
        #[term(cv=MS, accession=1002169, name="TOPP FeatureFinderMetabo", flags={0}, parents={["MS:1002165"]})]
        TOPPFeatureFinderMetabo,
        #[term(cv=MS, accession=1002170, name="TOPP FeatureFinderMRM", flags={0}, parents={["MS:1002165"]})]
        TOPPFeatureFinderMRM,
        #[term(cv=MS, accession=1002171, name="TOPP ProteinQuantifier", flags={0}, parents={["MS:1000752"]})]
        TOPPProteinQuantifier,
        #[term(cv=MS, accession=1002172, name="TOPP ConsensusMapNormalizer", flags={0}, parents={["MS:1000752"]})]
        TOPPConsensusMapNormalizer,
        #[term(cv=MS, accession=1002173, name="TOPP MapRTTransformer", flags={0}, parents={["MS:1000752"]})]
        TOPPMapRTTransformer,
        #[term(cv=MS, accession=1002174, name="TOPP feature linker", flags={0}, parents={["MS:1000752"]})]
        TOPPFeatureLinker,
        #[term(cv=MS, accession=1002175, name="TOPP FeatureLinkerLabeled", flags={0}, parents={["MS:1002174"]})]
        TOPPFeatureLinkerLabeled,
        #[term(cv=MS, accession=1002176, name="TOPP FeatureLinkerUnlabeled", flags={0}, parents={["MS:1002174"]})]
        TOPPFeatureLinkerUnlabeled,
        #[term(cv=MS, accession=1002177, name="TOPP FeatureLinkerUnlabeledQT", flags={0}, parents={["MS:1002174"]})]
        TOPPFeatureLinkerUnlabeledQT,
        #[term(cv=MS, accession=1002178, name="TOPP CompNovo", flags={0}, parents={["MS:1000752"]})]
        TOPPCompNovo,
        #[term(cv=MS, accession=1002179, name="TOPP CompNovoCID", flags={0}, parents={["MS:1000752"]})]
        TOPPCompNovoCID,
        #[term(cv=MS, accession=1002180, name="TOPP software adaptor", flags={0}, parents={["MS:1000752"]})]
        TOPPSoftwareAdaptor,
        #[term(cv=MS, accession=1002181, name="TOPP InspectAdapter", flags={0}, parents={["MS:1002180"]})]
        TOPPInspectAdapter,
        #[term(cv=MS, accession=1002182, name="TOPP MascotAdapter", flags={0}, parents={["MS:1002180"]})]
        TOPPMascotAdapter,
        #[term(cv=MS, accession=1002183, name="TOPP MascotAdapterOnline", flags={0}, parents={["MS:1002180"]})]
        TOPPMascotAdapterOnline,
        #[term(cv=MS, accession=1002184, name="TOPP OMSSAAdapter", flags={0}, parents={["MS:1002180"]})]
        TOPPOMSSAAdapter,
        #[term(cv=MS, accession=1002185, name="TOPP PepNovoAdapter", flags={0}, parents={["MS:1002180"]})]
        TOPPPepNovoAdapter,
        #[term(cv=MS, accession=1002186, name="TOPP XTandemAdapter", flags={0}, parents={["MS:1002180"]})]
        TOPPXTandemAdapter,
        #[term(cv=MS, accession=1002187, name="TOPP SpecLibSearcher", flags={0}, parents={["MS:1000752"]})]
        TOPPSpecLibSearcher,
        #[term(cv=MS, accession=1002188, name="TOPP ConsensusID", flags={0}, parents={["MS:1000752"]})]
        TOPPConsensusID,
        #[term(cv=MS, accession=1002189, name="TOPP IDConflictResolver", flags={0}, parents={["MS:1000752"]})]
        TOPPIDConflictResolver,
        #[term(cv=MS, accession=1002190, name="TOPP IDFilter", flags={0}, parents={["MS:1000752"]})]
        TOPPIDFilter,
        #[term(cv=MS, accession=1002191, name="TOPP IDMapper", flags={0}, parents={["MS:1000752"]})]
        TOPPIDMapper,
        #[term(cv=MS, accession=1002192, name="TOPP IDPosteriorErrorProbability", flags={0}, parents={["MS:1000752"]})]
        TOPPIDPosteriorErrorProbability,
        #[term(cv=MS, accession=1002193, name="TOPP IDRTCalibration", flags={0}, parents={["MS:1000752"]})]
        TOPPIDRTCalibration,
        #[term(cv=MS, accession=1002194, name="TOPP PeptideIndexer", flags={0}, parents={["MS:1000752"]})]
        TOPPPeptideIndexer,
        #[term(cv=MS, accession=1002195, name="TOPP PrecursorIonSelector", flags={0}, parents={["MS:1000752"]})]
        TOPPPrecursorIonSelector,
        #[term(cv=MS, accession=1002196, name="TOPP MRMMapper", flags={0}, parents={["MS:1000752"]})]
        TOPPMRMMapper,
        #[term(cv=MS, accession=1002197, name="TOPP OpenSwath component", flags={0}, parents={["MS:1000752"]})]
        TOPPOpenSwathComponent,
        #[term(cv=MS, accession=1002198, name="TOPP OpenSwathAnalyzer", flags={0}, parents={["MS:1002197"]})]
        TOPPOpenSwathAnalyzer,
        #[term(cv=MS, accession=1002199, name="TOPP OpenSwathChromatogramExtractor", flags={0}, parents={["MS:1002197"]})]
        TOPPOpenSwathChromatogramExtractor,
        #[term(cv=MS, accession=1002200, name="TOPP OpenSwathDecoyGenerator", flags={0}, parents={["MS:1002197"]})]
        TOPPOpenSwathDecoyGenerator,
        #[term(cv=MS, accession=1002201, name="TOPP OpenSwathFeatureXMLToTSV", flags={0}, parents={["MS:1002197"]})]
        TOPPOpenSwathFeatureXMLToTSV,
        #[term(cv=MS, accession=1002202, name="TOPP OpenSwathRTNormalizer", flags={0}, parents={["MS:1002197"]})]
        TOPPOpenSwathRTNormalizer,
        #[term(cv=MS, accession=1002203, name="TOPP ProteinInference", flags={0}, parents={["MS:1000752"]})]
        TOPPProteinInference,
        #[term(cv=MS, accession=1002204, name="TOPP FalseDiscoveryRate", flags={0}, parents={["MS:1000752"]})]
        TOPPFalseDiscoveryRate,
        #[term(cv=MS, accession=1002205, name="ProteoWizard msconvert", flags={0}, parents={["MS:1000615"]})]
        ProteoWizardMsconvert,
        #[term(cv=MS, accession=1002206, name="ProteoWizard idconvert", flags={0}, parents={["MS:1000615"]})]
        ProteoWizardIdconvert,
        #[term(cv=MS, accession=1002207, name="ProteoWizard chainsaw", flags={0}, parents={["MS:1000615"]})]
        ProteoWizardChainsaw,
        #[term(cv=MS, accession=1002208, name="ProteoWizard msaccess", flags={0}, parents={["MS:1000615"]})]
        ProteoWizardMsaccess,
        #[term(cv=MS, accession=1002209, name="ProteoWizard SeeMS", flags={0}, parents={["MS:1000615"]})]
        ProteoWizardSeeMS,
        #[term(cv=MS, accession=1002210, name="IsobariQ", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        IsobariQ,
        #[term(cv=MS, accession=1002220, name="MRMaid", flags={0}, parents={["MS:1000871"]})]
        MRMaid,
        #[term(cv=MS, accession=1002237, name="mzidLib", flags={1}, parents={["MS:1001456"]})]
        MzidLib,
        #[term(cv=MS, accession=1002238, name="mzidLib:Omssa2Mzid", flags={0}, parents={["MS:1002237"]})]
        MzidLibOmssa2Mzid,
        #[term(cv=MS, accession=1002239, name="mzidLib:Tandem2Mzid", flags={0}, parents={["MS:1002237"]})]
        MzidLibTandem2Mzid,
        #[term(cv=MS, accession=1002240, name="mzidLib:Csv2Mzid", flags={0}, parents={["MS:1002237"]})]
        MzidLibCsv2Mzid,
        #[term(cv=MS, accession=1002241, name="mzidLib:ProteoGrouper", flags={0}, parents={["MS:1002237"]})]
        MzidLibProteoGrouper,
        #[term(cv=MS, accession=1002242, name="mzidLib:Thresholder", flags={0}, parents={["MS:1002237"]})]
        MzidLibThresholder,
        #[term(cv=MS, accession=1002243, name="mzidLib:Perform emPAI on mzid", flags={0}, parents={["MS:1002237"]})]
        MzidLibPerformEmPAIOnMzid,
        #[term(cv=MS, accession=1002244, name="mzidLib:FalseDiscoveryRate", flags={0}, parents={["MS:1002237"]})]
        MzidLibFalseDiscoveryRate,
        #[term(cv=MS, accession=1002245, name="mzidLib:Mzidentml2Csv", flags={0}, parents={["MS:1002237"]})]
        MzidLibMzidentml2Csv,
        #[term(cv=MS, accession=1002246, name="mzidLib:CombineSearchEngines", flags={0}, parents={["MS:1002237"]})]
        MzidLibCombineSearchEngines,
        #[term(cv=MS, accession=1002247, name="mzidLib:InsertMetaDataFromFasta", flags={0}, parents={["MS:1002237"]})]
        MzidLibInsertMetaDataFromFasta,
        #[term(cv=MS, accession=1002251, name="Comet", flags={1}, parents={["MS:1001456"]})]
        Comet,
        #[term(cv=MS, accession=1002261, name="Byonic", flags={1}, parents={["MS:1001456"]})]
        Byonic,
        #[term(cv=MS, accession=1002285, name="Trans-Proteomic Pipeline", flags={1}, parents={["MS:1001456"]})]
        TransProteomicPipeline,
        #[term(cv=MS, accession=1002286, name="Trans-Proteomic Pipeline software", flags={1}, parents={["MS:1001456"]})]
        TransProteomicPipelineSoftware,
        #[term(cv=MS, accession=1002287, name="PeptideProphet", flags={0}, parents={["MS:1002286"]})]
        PeptideProphet,
        #[term(cv=MS, accession=1002288, name="iProphet", flags={0}, parents={["MS:1002286"]})]
        IProphet,
        #[term(cv=MS, accession=1002289, name="ProteinProphet", flags={0}, parents={["MS:1002286"]})]
        ProteinProphet,
        #[term(cv=MS, accession=1002290, name="XPRESS", flags={0}, parents={["MS:1002286"]})]
        XPRESS,
        #[term(cv=MS, accession=1002291, name="Libra", flags={0}, parents={["MS:1002286"]})]
        Libra,
        #[term(cv=MS, accession=1002292, name="PTMProphet", flags={0}, parents={["MS:1002286"]})]
        PTMProphet,
        #[term(cv=MS, accession=1002333, name="conversion software", flags={2}, parents={["MS:1001457"]})]
        ConversionSoftware,
        #[term(cv=MS, accession=1002334, name="ProCon", flags={0}, parents={["MS:1002333"]})]
        ProCon,
        #[term(cv=MS, accession=1002335, name="PRIDE Converter2", flags={0}, parents={["MS:1002333"]})]
        PRIDEConverter2,
        #[term(cv=MS, accession=1002336, name="Amanda", flags={1}, parents={["MS:1001456"]})]
        Amanda,
        #[term(cv=MS, accession=1002337, name="Andromeda", flags={1}, parents={["MS:1001456"]})]
        Andromeda,
        #[term(cv=MS, accession=1002342, name="mzmine", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        Mzmine,
        #[term(cv=MS, accession=1002344, name="Maltcms", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        Maltcms,
        #[term(cv=MS, accession=1002381, name="MALDI Solutions LC-MALDI", flags={7}, parents={["MS:1001455", "MS:1001456", "MS:1001457", "MS:1001557"]})]
        MALDISolutionsLCMALDI,
        #[term(cv=MS, accession=1002383, name="SCiLS software", flags={0}, parents={["MS:1000531"]})]
        SCiLSSoftware,
        #[term(cv=MS, accession=1002384, name="SCiLS Lab", flags={3}, parents={["MS:1002383", "MS:1001456", "MS:1001457"]})]
        SCiLSLab,
        #[term(cv=MS, accession=1002386, name="preprocessing software", flags={2}, parents={["MS:1001457"]})]
        PreprocessingSoftware,
        #[term(cv=MS, accession=1002387, name="PIA", flags={1}, parents={["MS:1002414", "MS:1001456"]})]
        PIA,
        #[term(cv=MS, accession=1002410, name="Anubis", flags={0}, parents={["MS:1000871", "MS:1001139"]})]
        Anubis,
        #[term(cv=MS, accession=1002414, name="postprocessing software", flags={2}, parents={["MS:1001457"]})]
        PostprocessingSoftware,
        #[term(cv=MS, accession=1002452, name="Maui", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        Maui,
        #[term(cv=MS, accession=1002458, name="PeptideShaker", flags={1}, parents={["MS:1001456"]})]
        PeptideShaker,
        #[term(cv=MS, accession=1002524, name="PepFinder", flags={2}, parents={["MS:1001457"]})]
        PepFinder,
        #[term(cv=MS, accession=1002543, name="xiFDR", flags={1}, parents={["MS:1001456"]})]
        XiFDR,
        #[term(cv=MS, accession=1002544, name="xi", flags={1}, parents={["MS:1001456"]})]
        Xi,
        #[term(cv=MS, accession=1002546, name="Skyline mzQuantML converter", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        SkylineMzQuantMLConverter,
        #[term(cv=MS, accession=1002574, name="ASAPRatio", flags={0}, parents={["MS:1002286"]})]
        ASAPRatio,
        #[term(cv=MS, accession=1002575, name="Tide", flags={1}, parents={["MS:1001456"]})]
        Tide,
        #[term(cv=MS, accession=1002596, name="ProLuCID", flags={1}, parents={["MS:1001456"]})]
        ProLuCID,
        #[term(cv=MS, accession=1002598, name="DTASelect", flags={1}, parents={["MS:1001456"]})]
        DTASelect,
        #[term(cv=MS, accession=1002645, name="MSDK", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        MSDK,
        #[term(cv=MS, accession=1002661, name="Morpheus", flags={1}, parents={["MS:1001456"]})]
        Morpheus,
        #[term(cv=MS, accession=1002673, name="OpenXQuest", flags={0}, parents={["MS:1000752"]})]
        OpenXQuest,
        #[term(cv=MS, accession=1002714, name="FLASHDeconv", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        FLASHDeconv,
        #[term(cv=MS, accession=1002720, name="MSPathFinder", flags={1}, parents={["MS:1001456"]})]
        MSPathFinder,
        #[term(cv=MS, accession=1002750, name="NIST MSPepSearch", flags={1}, parents={["MS:1001456"]})]
        NISTMSPepSearch,
        #[term(cv=MS, accession=1002826, name="MetaMorpheus", flags={1}, parents={["MS:1001456"]})]
        MetaMorpheus,
        #[term(cv=MS, accession=1002869, name="mzR", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        MzR,
        #[term(cv=MS, accession=1002870, name="MSnbase", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        MSnbase,
        #[term(cv=MS, accession=1002871, name="CAMERA", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        CAMERA,
        #[term(cv=MS, accession=1002878, name="small molecule analysis software", flags={1}, parents={["MS:1001456"]})]
        SmallMoleculeAnalysisSoftware,
        #[term(cv=MS, accession=1002879, name="Progenesis QI", flags={0}, parents={["MS:1002878"]})]
        ProgenesisQI,
        #[term(cv=MS, accession=1002880, name="Compound Discoverer", flags={0}, parents={["MS:1002878"]})]
        CompoundDiscoverer,
        #[term(cv=MS, accession=1002881, name="MyCompoundID", flags={0}, parents={["MS:1002878"]})]
        MyCompoundID,
        #[term(cv=MS, accession=1002901, name="TopPIC", flags={1}, parents={["MS:1001456"]})]
        TopPIC,
        #[term(cv=MS, accession=1002902, name="TopFD", flags={1}, parents={["MS:1001456"]})]
        TopFD,
        #[term(cv=MS, accession=1002903, name="TopMG", flags={1}, parents={["MS:1001456"]})]
        TopMG,
        #[term(cv=MS, accession=1002964, name="lipidomics analysis software", flags={0}, parents={["MS:1002878"]})]
        LipidomicsAnalysisSoftware,
        #[term(cv=MS, accession=1002965, name="Lipid Data Analyzer", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        LipidDataAnalyzer,
        #[term(cv=MS, accession=1002967, name="LipidHunter", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        LipidHunter,
        #[term(cv=MS, accession=1002968, name="LipidXplorer", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        LipidXplorer,
        #[term(cv=MS, accession=1002969, name="LipidMatch", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        LipidMatch,
        #[term(cv=MS, accession=1002970, name="Greazy", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        Greazy,
        #[term(cv=MS, accession=1002971, name="LipidBlast", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        LipidBlast,
        #[term(cv=MS, accession=1002972, name="Lipid-Pro", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        LipidPro,
        #[term(cv=MS, accession=1002973, name="LipidFinder", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        LipidFinder,
        #[term(cv=MS, accession=1002974, name="LipiDex", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        LipiDex,
        #[term(cv=MS, accession=1002975, name="LIQUID", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        LIQUID,
        #[term(cv=MS, accession=1002976, name="ALEX", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        ALEX,
        #[term(cv=MS, accession=1002977, name="ALEX123", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        ALEX123,
        #[term(cv=MS, accession=1002978, name="LIMSA", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        LIMSA,
        #[term(cv=MS, accession=1002979, name="LOBSTAHS", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        LOBSTAHS,
        #[term(cv=MS, accession=1002980, name="LipidQA", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        LipidQA,
        #[term(cv=MS, accession=1002981, name="Proline", flags={1}, parents={["MS:1001456"]})]
        Proline,
        #[term(cv=MS, accession=1002982, name="PepNovo", flags={1}, parents={["MS:1001456"]})]
        PepNovo,
        #[term(cv=MS, accession=1002983, name="pNovo", flags={1}, parents={["MS:1001456"]})]
        PNovo,
        #[term(cv=MS, accession=1002984, name="Novor", flags={1}, parents={["MS:1001456"]})]
        Novor,
        #[term(cv=MS, accession=1002987, name="IdentiPy", flags={1}, parents={["MS:1001456"]})]
        IdentiPy,
        #[term(cv=MS, accession=1002990, name="ms_deisotope", flags={2}, parents={["MS:1001457"]})]
        MsDeisotope,
        #[term(cv=MS, accession=1002991, name="python-psims", flags={0}, parents={["MS:1002333"]})]
        PythonPsims,
        #[term(cv=MS, accession=1003010, name="LPPtiger", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        LPPtiger,
        #[term(cv=MS, accession=1003011, name="pFind", flags={1}, parents={["MS:1001456"]})]
        PFind,
        #[term(cv=MS, accession=1003013, name="i3tms", flags={1}, parents={["MS:1001456"]})]
        I3tms,
        #[term(cv=MS, accession=1003014, name="MSFragger", flags={1}, parents={["MS:1001456"]})]
        MSFragger,
        #[term(cv=MS, accession=1003018, name="Philosopher", flags={1}, parents={["MS:1001456"]})]
        Philosopher,
        #[term(cv=MS, accession=1003023, name="OpenPepXL", flags={0}, parents={["MS:1000752"]})]
        OpenPepXL,
        #[term(cv=MS, accession=1003082, name="MS-DIAL", flags={2}, parents={["MS:1002878", "MS:1001457"]})]
        MSDIAL,
        #[term(cv=MS, accession=1003108, name="PatternLab", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        PatternLab,
        #[term(cv=MS, accession=1003109, name="SIM-XL", flags={1}, parents={["MS:1001456"]})]
        SIMXL,
        #[term(cv=MS, accession=1003111, name="QUIN-XL", flags={0}, parents={["MS:1001139"]})]
        QUINXL,
        #[term(cv=MS, accession=1003118, name="EPIFANY", flags={1}, parents={["MS:1001456", "MS:1000752"]})]
        EPIFANY,
        #[term(cv=MS, accession=1003141, name="ProSight", flags={1}, parents={["MS:1001456"]})]
        ProSight,
        #[term(cv=MS, accession=1003142, name="TDPortal", flags={1}, parents={["MS:1001456"]})]
        TDPortal,
        #[term(cv=MS, accession=1003145, name="ThermoRawFileParser", flags={2}, parents={["MS:1001457"]})]
        ThermoRawFileParser,
        #[term(cv=MS, accession=1003146, name="pyteomics", flags={1}, parents={["MS:1001456"]})]
        Pyteomics,
        #[term(cv=MS, accession=1003162, name="PTX-QC", flags={1}, parents={["MS:1001456"]})]
        PTXQC,
        #[term(cv=MS, accession=1003164, name="QuaMeter IDFree", flags={1}, parents={["MS:1001456"]})]
        QuaMeterIDFree,
        #[term(cv=MS, accession=1003165, name="iMonDB", flags={1}, parents={["MS:1001456"]})]
        IMonDB,
        #[term(cv=MS, accession=1003202, name="BiblioSpec", flags={1}, parents={["MS:1001456", "MS:1003207"]})]
        BiblioSpec,
        #[term(cv=MS, accession=1003207, name="library creation software", flags={0}, parents={["MS:1000531", "MS:1003171"]})]
        LibraryCreationSoftware,
        #[term(cv=MS, accession=1003232, name="PeakForest", flags={1}, parents={["MS:1001456", "MS:1003207", "MS:1002878"]})]
        PeakForest,
        #[term(cv=MS, accession=1003253, name="DIA-NN", flags={1}, parents={["MS:1001139", "MS:1001456", "MS:1003207"]})]
        DIANN,
        #[term(cv=MS, accession=1003281, name="Casanovo", flags={1}, parents={["MS:1001456"]})]
        Casanovo,
        #[term(cv=MS, accession=1003309, name="Goslin", flags={2}, parents={["MS:1001457", "MS:1002414", "MS:1002964"]})]
        Goslin,
        #[term(cv=MS, accession=1003357, name="ANN-SoLo", flags={1}, parents={["MS:1001456"]})]
        ANNSoLo,
        #[term(cv=MS, accession=1003376, name="ChemClipse", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        ChemClipse,
        #[term(cv=MS, accession=1003377, name="OpenChrom", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        OpenChrom,
        #[term(cv=MS, accession=1003382, name="waters_connect", flags={7}, parents={["MS:1000694", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        WatersConnect,
        #[term(cv=MS, accession=1003386, name="Spectra", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        Spectra,
        #[term(cv=MS, accession=1003387, name="MetaboAnnotation", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        MetaboAnnotation,
        #[term(cv=MS, accession=1003388, name="CompoundDb", flags={3}, parents={["MS:1001456", "MS:1001457", "MS:1003207"]})]
        CompoundDb,
        #[term(cv=MS, accession=1003399, name="quality control software", flags={2}, parents={["MS:1001457"]})]
        QualityControlSoftware,
        #[term(cv=MS, accession=1003400, name="rmzqc", flags={0}, parents={["MS:1003399"]})]
        Rmzqc,
        #[term(cv=MS, accession=1003401, name="jmzqc", flags={0}, parents={["MS:1003399"]})]
        Jmzqc,
        #[term(cv=MS, accession=1003402, name="pymzqc", flags={0}, parents={["MS:1003399"]})]
        Pymzqc,
        #[term(cv=MS, accession=1003405, name="mzRecal", flags={2}, parents={["MS:1001457"]})]
        MzRecal,
        #[term(cv=MS, accession=1003406, name="spectrum clustering software", flags={0}, parents={["MS:1000531"]})]
        SpectrumClusteringSoftware,
        #[term(cv=MS, accession=1003407, name="Scout", flags={1}, parents={["MS:1001456"]})]
        Scout,
        #[term(cv=MS, accession=1003413, name="Kojak", flags={1}, parents={["MS:1001456"]})]
        Kojak,
        #[term(cv=MS, accession=1003425, name="quantms", flags={1}, parents={["MS:1001456", "MS:1001139"]})]
        Quantms,
        #[term(cv=MS, accession=1003426, name="xQuest/xProphet", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        XQuestXProphet,
        #[term(cv=MS, accession=1003427, name="PeakView", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        PeakView,
        #[term(cv=MS, accession=1003428, name="Perseus", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        Perseus,
        #[term(cv=MS, accession=1003429, name="FragPipe", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        FragPipe,
        #[term(cv=MS, accession=1003430, name="OpenMS", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        OpenMS,
        #[term(cv=MS, accession=1003431, name="pLink", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        PLink,
        #[term(cv=MS, accession=1003432, name="pLink2", flags={1}, parents={["MS:1003431", "MS:1001139", "MS:1001456"]})]
        PLink2,
        #[term(cv=MS, accession=4000151, name="MsQuality", flags={1}, parents={["MS:1001456"]})]
        MsQuality,
    }
    //[[[end]]] (checksum: ef73f9f08114185c6c9e5395d215df54)
}

#[cfg(test)]
mod test {
    use crate::params::ParamDescribed;

    use super::*;

    #[test]
    fn cvmap_test() {
        assert_eq!(
            SoftwareTerm::SCIEXTOFTOFSeriesExplorerSoftware.accession(),
            1001483
        );
        assert_eq!(
            SoftwareTerm::SCIEXTOFTOFSeriesExplorerSoftware.name(),
            "SCIEX TOF/TOF Series Explorer Software"
        );
        assert_eq!(
            SoftwareTerm::SCIEXTOFTOFSeriesExplorerSoftware.flags(),
            SoftwareType::Analysis | SoftwareType::Acquisition | SoftwareType::DataProcessing
        );
        assert!(
            SoftwareTerm::SCIEXTOFTOFSeriesExplorerSoftware.flags().is_analysis(),
        );
    }

    #[test]
    fn sw_test() {
        let mut sw = Software::new("foo".into(), "v0.1.0".into(), vec![custom_software_name("foo")]);
        assert_eq!(sw.id, "foo");
        assert_eq!(sw.version, "v0.1.0");
        assert!(sw.find_software_term().is_some());
        sw.add_param(SoftwareTerm::ANNSoLo.into());
        assert!(!sw.is_analysis());
        sw.params_mut().reverse();
        assert!(sw.is_analysis());
        assert!(sw.find_software_term().is_some());


        let id = Software::find_unique_id("foo", [sw].iter());
        assert_eq!(id, "foo_0");
    }
}
