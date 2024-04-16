use std::collections::HashSet;

use crate::params::{ControlledVocabulary, ParamList};
use crate::{impl_param_described, Param};

/// A piece of software that was associated with the acquisition, transformation or otherwise
/// processing of mass spectrometry data. See <https://peptideatlas.org/tmp/mzML1.1.0.html#software>
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Software {
    /// A unique identifier for the software within processing metadata
    pub id: String,
    /// A string denoting a particular software version, but does no guarantee is given for its format
    pub version: String,
    /// Any associated vocabulary terms, including actual software name and type
    pub params: ParamList,
}

bitflags::bitflags! {
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

    pub const fn is_analysis(&self) -> bool {
        self.contains(Self::Analysis)
    }

    pub const fn is_data_processing(&self) -> bool {
        self.contains(Self::DataProcessing)
    }

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

    /// Find a unique identifier from an iterator over software IDs
    pub fn find_unique_id<'a>(id_stem: &str, softwares: impl IntoIterator<Item = &'a Self>) -> String {
        let software_ids: HashSet<_> = softwares.into_iter().map(|sw| &sw.id).collect();
        (0..)
            .into_iter()
            .map(|i| format!("{id_stem}_{i}"))
            .filter(|s| !software_ids.contains(s))
            .next()
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
    #[allow(unused)]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_software.py']).decode('utf8')
    for line in buf.splitlines():
        cog.outl(line)
    ]]]*/
    pub enum SoftwareTerm {
        #[term(cv=MS, accession=1000531, name="software", flags={0})]
        Software,
        #[term(cv=MS, accession=1000689, name="Agilent software", flags={0})]
        AgilentSoftware,
        #[term(cv=MS, accession=1000690, name="SCIEX software", flags={0})]
        SCIEXSoftware,
        #[term(cv=MS, accession=1000691, name="Applied Biosystems software", flags={0})]
        AppliedBiosystemsSoftware,
        #[term(cv=MS, accession=1000692, name="Bruker software", flags={0})]
        BrukerSoftware,
        #[term(cv=MS, accession=1000693, name="Thermo Finnigan software", flags={0})]
        ThermoFinniganSoftware,
        #[term(cv=MS, accession=1000694, name="Waters software", flags={0})]
        WatersSoftware,
        #[term(cv=MS, accession=1000706, name="apexControl", flags={4})]
        ApexControl,
        #[term(cv=MS, accession=1000707, name="BioTools", flags={3})]
        BioTools,
        #[term(cv=MS, accession=1000708, name="CLINPROT", flags={3})]
        CLINPROT,
        #[term(cv=MS, accession=1000709, name="CLINPROT micro", flags={3})]
        CLINPROTMicro,
        #[term(cv=MS, accession=1000710, name="CLINPROT robot", flags={0})]
        CLINPROTRobot,
        #[term(cv=MS, accession=1000711, name="ClinProTools", flags={0})]
        ClinProTools,
        #[term(cv=MS, accession=1000712, name="Compass", flags={3})]
        Compass,
        #[term(cv=MS, accession=1000713, name="Compass for HCT/esquire", flags={3})]
        CompassForHCTEsquire,
        #[term(cv=MS, accession=1000714, name="Compass for micrOTOF", flags={3})]
        CompassForMicrOTOF,
        #[term(cv=MS, accession=1000715, name="Compass OpenAccess", flags={3})]
        CompassOpenAccess,
        #[term(cv=MS, accession=1000716, name="Compass Security Pack", flags={0})]
        CompassSecurityPack,
        #[term(cv=MS, accession=1000717, name="CompassXport", flags={2})]
        CompassXport,
        #[term(cv=MS, accession=1000718, name="CompassXtract", flags={2})]
        CompassXtract,
        #[term(cv=MS, accession=1000719, name="DataAnalysis", flags={3})]
        DataAnalysis,
        #[term(cv=MS, accession=1000720, name="dpControl", flags={4})]
        DpControl,
        #[term(cv=MS, accession=1000721, name="esquireControl", flags={4})]
        EsquireControl,
        #[term(cv=MS, accession=1000722, name="flexImaging", flags={3})]
        FlexImaging,
        #[term(cv=MS, accession=1000723, name="GENOLINK", flags={0})]
        GENOLINK,
        #[term(cv=MS, accession=1000724, name="GenoTools", flags={0})]
        GenoTools,
        #[term(cv=MS, accession=1000725, name="HCTcontrol", flags={4})]
        HCTcontrol,
        #[term(cv=MS, accession=1000726, name="micrOTOFcontrol", flags={4})]
        MicrOTOFcontrol,
        #[term(cv=MS, accession=1000727, name="PolyTools", flags={0})]
        PolyTools,
        #[term(cv=MS, accession=1000728, name="ProfileAnalysis", flags={3})]
        ProfileAnalysis,
        #[term(cv=MS, accession=1000729, name="PROTEINEER", flags={0})]
        PROTEINEER,
        #[term(cv=MS, accession=1000730, name="PROTEINEER dp", flags={0})]
        PROTEINEERDp,
        #[term(cv=MS, accession=1000731, name="PROTEINEER fc", flags={0})]
        PROTEINEERFc,
        #[term(cv=MS, accession=1000732, name="PROTEINEER spII", flags={0})]
        PROTEINEERSpII,
        #[term(cv=MS, accession=1000733, name="PROTEINEER-LC", flags={0})]
        PROTEINEERLC,
        #[term(cv=MS, accession=1000734, name="ProteinScape", flags={1})]
        ProteinScape,
        #[term(cv=MS, accession=1000735, name="PureDisk", flags={0})]
        PureDisk,
        #[term(cv=MS, accession=1000736, name="QuantAnalysis", flags={3})]
        QuantAnalysis,
        #[term(cv=MS, accession=1000737, name="spControl", flags={4})]
        SpControl,
        #[term(cv=MS, accession=1000738, name="TargetAnalysis", flags={0})]
        TargetAnalysis,
        #[term(cv=MS, accession=1000739, name="WARP-LC", flags={0})]
        WARPLC,
        #[term(cv=MS, accession=1000799, name="custom unreleased software tool", flags={0})]
        CustomUnreleasedSoftwareTool,
        #[term(cv=MS, accession=1000817, name="HyStar", flags={0})]
        HyStar,
        #[term(cv=MS, accession=1000871, name="SRM software", flags={0})]
        SRMSoftware,
        #[term(cv=MS, accession=1000872, name="MaRiMba", flags={0})]
        MaRiMba,
        #[term(cv=MS, accession=1000873, name="peptide attribute calculation software", flags={0})]
        PeptideAttributeCalculationSoftware,
        #[term(cv=MS, accession=1000874, name="SSRCalc", flags={0})]
        SSRCalc,
        #[term(cv=MS, accession=1000922, name="Skyline", flags={0})]
        Skyline,
        #[term(cv=MS, accession=1000923, name="TIQAM", flags={0})]
        TIQAM,
        #[term(cv=MS, accession=1000925, name="ATAQS", flags={0})]
        ATAQS,
        #[term(cv=MS, accession=1001139, name="quantitation software name", flags={0})]
        QuantitationSoftwareName,
        #[term(cv=MS, accession=1001455, name="acquisition software", flags={0})]
        AcquisitionSoftware,
        #[term(cv=MS, accession=1001456, name="analysis software", flags={0})]
        AnalysisSoftware,
        #[term(cv=MS, accession=1001457, name="data processing software", flags={0})]
        DataProcessingSoftware,
        #[term(cv=MS, accession=1001461, name="greylag", flags={1})]
        Greylag,
        #[term(cv=MS, accession=1001475, name="OMSSA", flags={1})]
        OMSSA,
        #[term(cv=MS, accession=1001476, name="X!Tandem", flags={1})]
        XTandem,
        #[term(cv=MS, accession=1001477, name="SpectraST", flags={1})]
        SpectraST,
        #[term(cv=MS, accession=1001478, name="Mascot Parser", flags={1})]
        MascotParser,
        #[term(cv=MS, accession=1001483, name="SCIEX TOF/TOF Series Explorer Software", flags={7})]
        SCIEXTOFTOFSeriesExplorerSoftware,
        #[term(cv=MS, accession=1001487, name="ProteinExtractor", flags={1})]
        ProteinExtractor,
        #[term(cv=MS, accession=1001488, name="Mascot Distiller", flags={1})]
        MascotDistiller,
        #[term(cv=MS, accession=1001489, name="Mascot Integra", flags={1})]
        MascotIntegra,
        #[term(cv=MS, accession=1001490, name="Percolator", flags={1})]
        Percolator,
        #[term(cv=MS, accession=1001557, name="Shimadzu Corporation software", flags={0})]
        ShimadzuCorporationSoftware,
        #[term(cv=MS, accession=1001558, name="MALDI Solutions", flags={7})]
        MALDISolutions,
        #[term(cv=MS, accession=1001561, name="Scaffold", flags={1})]
        Scaffold,
        #[term(cv=MS, accession=1001582, name="XCMS", flags={3})]
        XCMS,
        #[term(cv=MS, accession=1001583, name="MaxQuant", flags={1})]
        MaxQuant,
        #[term(cv=MS, accession=1001585, name="MyriMatch", flags={1})]
        MyriMatch,
        #[term(cv=MS, accession=1001586, name="DirecTag", flags={1})]
        DirecTag,
        #[term(cv=MS, accession=1001587, name="TagRecon", flags={1})]
        TagRecon,
        #[term(cv=MS, accession=1001588, name="Pepitome", flags={1})]
        Pepitome,
        #[term(cv=MS, accession=1001795, name="Empower", flags={3})]
        Empower,
        #[term(cv=MS, accession=1001796, name="UNIFY", flags={3})]
        UNIFY,
        #[term(cv=MS, accession=1001798, name="LECO software", flags={0})]
        LECOSoftware,
        #[term(cv=MS, accession=1001799, name="ChromaTOF software", flags={7})]
        ChromaTOFSoftware,
        #[term(cv=MS, accession=1001830, name="Progenesis LC-MS", flags={0})]
        ProgenesisLCMS,
        #[term(cv=MS, accession=1001831, name="SILACAnalyzer", flags={0})]
        SILACAnalyzer,
        #[term(cv=MS, accession=1001877, name="ChromaTOF HRT software", flags={7})]
        ChromaTOFHRTSoftware,
        #[term(cv=MS, accession=1001878, name="MALDI Solutions Microbial Identification", flags={0})]
        MALDISolutionsMicrobialIdentification,
        #[term(cv=MS, accession=1001886, name="SQID", flags={3})]
        SQID,
        #[term(cv=MS, accession=1001912, name="PinPoint", flags={3})]
        PinPoint,
        #[term(cv=MS, accession=1001914, name="pymzML", flags={2})]
        PymzML,
        #[term(cv=MS, accession=1001946, name="PEAKS Studio", flags={3})]
        PEAKSStudio,
        #[term(cv=MS, accession=1001947, name="PEAKS Online", flags={3})]
        PEAKSOnline,
        #[term(cv=MS, accession=1001948, name="PEAKS Node", flags={3})]
        PEAKSNode,
        #[term(cv=MS, accession=1001949, name="BSI software", flags={0})]
        BSISoftware,
        #[term(cv=MS, accession=1001973, name="DeBunker", flags={1})]
        DeBunker,
        #[term(cv=MS, accession=1001977, name="MSQuant", flags={1})]
        MSQuant,
        #[term(cv=MS, accession=1001984, name="Ascore software", flags={1})]
        AscoreSoftware,
        #[term(cv=MS, accession=1002043, name="ProteinProspector", flags={1})]
        ProteinProspector,
        #[term(cv=MS, accession=1002047, name="MS-GF", flags={1})]
        MSGF,
        #[term(cv=MS, accession=1002048, name="MS-GF+", flags={1})]
        MSGFplus,
        #[term(cv=MS, accession=1002059, name="Microsoft Excel", flags={0})]
        MicrosoftExcel,
        #[term(cv=MS, accession=1002063, name="FindPairs", flags={0})]
        FindPairs,
        #[term(cv=MS, accession=1002076, name="PAnalyzer", flags={1})]
        PAnalyzer,
        #[term(cv=MS, accession=1002123, name="x-Tracker", flags={0})]
        XTracker,
        #[term(cv=MS, accession=1002124, name="ProteoSuite", flags={0})]
        ProteoSuite,
        #[term(cv=MS, accession=1002129, name="ITRAQAnalyzer", flags={0})]
        ITRAQAnalyzer,
        #[term(cv=MS, accession=1002210, name="IsobariQ", flags={1})]
        IsobariQ,
        #[term(cv=MS, accession=1002220, name="MRMaid", flags={0})]
        MRMaid,
        #[term(cv=MS, accession=1002237, name="mzidLib", flags={1})]
        MzidLib,
        #[term(cv=MS, accession=1002238, name="mzidLib:Omssa2Mzid", flags={0})]
        MzidLibOmssa2Mzid,
        #[term(cv=MS, accession=1002239, name="mzidLib:Tandem2Mzid", flags={0})]
        MzidLibTandem2Mzid,
        #[term(cv=MS, accession=1002240, name="mzidLib:Csv2Mzid", flags={0})]
        MzidLibCsv2Mzid,
        #[term(cv=MS, accession=1002241, name="mzidLib:ProteoGrouper", flags={0})]
        MzidLibProteoGrouper,
        #[term(cv=MS, accession=1002242, name="mzidLib:Thresholder", flags={0})]
        MzidLibThresholder,
        #[term(cv=MS, accession=1002243, name="mzidLib:Perform emPAI on mzid", flags={0})]
        MzidLibPerformEmPAIOnMzid,
        #[term(cv=MS, accession=1002244, name="mzidLib:FalseDiscoveryRate", flags={0})]
        MzidLibFalseDiscoveryRate,
        #[term(cv=MS, accession=1002245, name="mzidLib:Mzidentml2Csv", flags={0})]
        MzidLibMzidentml2Csv,
        #[term(cv=MS, accession=1002246, name="mzidLib:CombineSearchEngines", flags={0})]
        MzidLibCombineSearchEngines,
        #[term(cv=MS, accession=1002247, name="mzidLib:InsertMetaDataFromFasta", flags={0})]
        MzidLibInsertMetaDataFromFasta,
        #[term(cv=MS, accession=1002251, name="Comet", flags={1})]
        Comet,
        #[term(cv=MS, accession=1002261, name="Byonic", flags={1})]
        Byonic,
        #[term(cv=MS, accession=1002285, name="Trans-Proteomic Pipeline", flags={1})]
        TransProteomicPipeline,
        #[term(cv=MS, accession=1002286, name="Trans-Proteomic Pipeline software", flags={1})]
        TransProteomicPipelineSoftware,
        #[term(cv=MS, accession=1002287, name="PeptideProphet", flags={0})]
        PeptideProphet,
        #[term(cv=MS, accession=1002288, name="iProphet", flags={0})]
        IProphet,
        #[term(cv=MS, accession=1002289, name="ProteinProphet", flags={0})]
        ProteinProphet,
        #[term(cv=MS, accession=1002290, name="XPRESS", flags={0})]
        XPRESS,
        #[term(cv=MS, accession=1002291, name="Libra", flags={0})]
        Libra,
        #[term(cv=MS, accession=1002292, name="PTMProphet", flags={0})]
        PTMProphet,
        #[term(cv=MS, accession=1002333, name="conversion software", flags={2})]
        ConversionSoftware,
        #[term(cv=MS, accession=1002334, name="ProCon", flags={0})]
        ProCon,
        #[term(cv=MS, accession=1002335, name="PRIDE Converter2", flags={0})]
        PRIDEConverter2,
        #[term(cv=MS, accession=1002336, name="Amanda", flags={1})]
        Amanda,
        #[term(cv=MS, accession=1002337, name="Andromeda", flags={1})]
        Andromeda,
        #[term(cv=MS, accession=1002342, name="MZmine", flags={3})]
        MZmine,
        #[term(cv=MS, accession=1002344, name="Maltcms", flags={3})]
        Maltcms,
        #[term(cv=MS, accession=1002381, name="MALDI Solutions LC-MALDI", flags={7})]
        MALDISolutionsLCMALDI,
        #[term(cv=MS, accession=1002383, name="SCiLS software", flags={0})]
        SCiLSSoftware,
        #[term(cv=MS, accession=1002384, name="SCiLS Lab", flags={3})]
        SCiLSLab,
        #[term(cv=MS, accession=1002386, name="preprocessing software", flags={2})]
        PreprocessingSoftware,
        #[term(cv=MS, accession=1002387, name="PIA", flags={1})]
        PIA,
        #[term(cv=MS, accession=1002410, name="Anubis", flags={0})]
        Anubis,
        #[term(cv=MS, accession=1002414, name="postprocessing software", flags={2})]
        PostprocessingSoftware,
        #[term(cv=MS, accession=1002452, name="Maui", flags={3})]
        Maui,
        #[term(cv=MS, accession=1002458, name="PeptideShaker", flags={1})]
        PeptideShaker,
        #[term(cv=MS, accession=1002524, name="PepFinder", flags={2})]
        PepFinder,
        #[term(cv=MS, accession=1002543, name="xiFDR", flags={1})]
        XiFDR,
        #[term(cv=MS, accession=1002544, name="xi", flags={1})]
        Xi,
        #[term(cv=MS, accession=1002546, name="Skyline mzQuantML converter", flags={1})]
        SkylineMzQuantMLConverter,
        #[term(cv=MS, accession=1002574, name="ASAPRatio", flags={0})]
        ASAPRatio,
        #[term(cv=MS, accession=1002575, name="Tide", flags={1})]
        Tide,
        #[term(cv=MS, accession=1002596, name="ProLuCID", flags={1})]
        ProLuCID,
        #[term(cv=MS, accession=1002598, name="DTASelect", flags={1})]
        DTASelect,
        #[term(cv=MS, accession=1002645, name="MSDK", flags={3})]
        MSDK,
        #[term(cv=MS, accession=1002661, name="Morpheus", flags={1})]
        Morpheus,
        #[term(cv=MS, accession=1002714, name="FLASHDeconv", flags={3})]
        FLASHDeconv,
        #[term(cv=MS, accession=1002720, name="MSPathFinder", flags={1})]
        MSPathFinder,
        #[term(cv=MS, accession=1002750, name="NIST MSPepSearch", flags={1})]
        NISTMSPepSearch,
        #[term(cv=MS, accession=1002826, name="MetaMorpheus", flags={1})]
        MetaMorpheus,
        #[term(cv=MS, accession=1002869, name="mzR", flags={3})]
        MzR,
        #[term(cv=MS, accession=1002870, name="MSnbase", flags={3})]
        MSnbase,
        #[term(cv=MS, accession=1002871, name="CAMERA", flags={3})]
        CAMERA,
        #[term(cv=MS, accession=1002878, name="small molecule analysis software", flags={1})]
        SmallMoleculeAnalysisSoftware,
        #[term(cv=MS, accession=1002879, name="Progenesis QI", flags={0})]
        ProgenesisQI,
        #[term(cv=MS, accession=1002880, name="Compound Discoverer", flags={0})]
        CompoundDiscoverer,
        #[term(cv=MS, accession=1002881, name="MyCompoundID", flags={0})]
        MyCompoundID,
        #[term(cv=MS, accession=1002901, name="TopPIC", flags={1})]
        TopPIC,
        #[term(cv=MS, accession=1002902, name="TopFD", flags={1})]
        TopFD,
        #[term(cv=MS, accession=1002903, name="TopMG", flags={1})]
        TopMG,
        #[term(cv=MS, accession=1002964, name="lipidomics analysis software", flags={0})]
        LipidomicsAnalysisSoftware,
        #[term(cv=MS, accession=1002965, name="Lipid Data Analyzer", flags={2})]
        LipidDataAnalyzer,
        #[term(cv=MS, accession=1002967, name="LipidHunter", flags={2})]
        LipidHunter,
        #[term(cv=MS, accession=1002968, name="LipidXplorer", flags={2})]
        LipidXplorer,
        #[term(cv=MS, accession=1002969, name="LipidMatch", flags={2})]
        LipidMatch,
        #[term(cv=MS, accession=1002970, name="Greazy", flags={2})]
        Greazy,
        #[term(cv=MS, accession=1002971, name="LipidBlast", flags={2})]
        LipidBlast,
        #[term(cv=MS, accession=1002972, name="Lipid-Pro", flags={2})]
        LipidPro,
        #[term(cv=MS, accession=1002973, name="LipidFinder", flags={2})]
        LipidFinder,
        #[term(cv=MS, accession=1002974, name="LipiDex", flags={2})]
        LipiDex,
        #[term(cv=MS, accession=1002975, name="LIQUID", flags={2})]
        LIQUID,
        #[term(cv=MS, accession=1002976, name="ALEX", flags={2})]
        ALEX,
        #[term(cv=MS, accession=1002977, name="ALEX123", flags={2})]
        ALEX123,
        #[term(cv=MS, accession=1002978, name="LIMSA", flags={2})]
        LIMSA,
        #[term(cv=MS, accession=1002979, name="LOBSTAHS", flags={2})]
        LOBSTAHS,
        #[term(cv=MS, accession=1002980, name="LipidQA", flags={2})]
        LipidQA,
        #[term(cv=MS, accession=1002981, name="Proline", flags={1})]
        Proline,
        #[term(cv=MS, accession=1002982, name="PepNovo", flags={1})]
        PepNovo,
        #[term(cv=MS, accession=1002983, name="pNovo", flags={1})]
        PNovo,
        #[term(cv=MS, accession=1002984, name="Novor", flags={1})]
        Novor,
        #[term(cv=MS, accession=1002987, name="IdentiPy", flags={1})]
        IdentiPy,
        #[term(cv=MS, accession=1002990, name="ms_deisotope", flags={2})]
        MsDeisotope,
        #[term(cv=MS, accession=1002991, name="python-psims", flags={0})]
        PythonPsims,
        #[term(cv=MS, accession=1003010, name="LPPtiger", flags={2})]
        LPPtiger,
        #[term(cv=MS, accession=1003011, name="pFind", flags={1})]
        PFind,
        #[term(cv=MS, accession=1003013, name="i3tms", flags={1})]
        I3tms,
        #[term(cv=MS, accession=1003014, name="MSFragger", flags={1})]
        MSFragger,
        #[term(cv=MS, accession=1003018, name="Philosopher", flags={1})]
        Philosopher,
        #[term(cv=MS, accession=1003082, name="MS-DIAL", flags={2})]
        MSDIAL,
        #[term(cv=MS, accession=1003108, name="PatternLab", flags={1})]
        PatternLab,
        #[term(cv=MS, accession=1003109, name="SIM-XL", flags={1})]
        SIMXL,
        #[term(cv=MS, accession=1003111, name="QUIN-XL", flags={0})]
        QUINXL,
        #[term(cv=MS, accession=1003118, name="EPIFANY", flags={1})]
        EPIFANY,
        #[term(cv=MS, accession=1003141, name="ProSight", flags={1})]
        ProSight,
        #[term(cv=MS, accession=1003142, name="TDPortal", flags={1})]
        TDPortal,
        #[term(cv=MS, accession=1003145, name="ThermoRawFileParser", flags={2})]
        ThermoRawFileParser,
        #[term(cv=MS, accession=1003146, name="pyteomics", flags={1})]
        Pyteomics,
        #[term(cv=MS, accession=1003162, name="PTX-QC", flags={1})]
        PTXQC,
        #[term(cv=MS, accession=1003164, name="QuaMeter IDFree", flags={1})]
        QuaMeterIDFree,
        #[term(cv=MS, accession=1003165, name="iMonDB", flags={1})]
        IMonDB,
        #[term(cv=MS, accession=1003202, name="BiblioSpec", flags={1})]
        BiblioSpec,
        #[term(cv=MS, accession=1003207, name="library creation software", flags={0})]
        LibraryCreationSoftware,
        #[term(cv=MS, accession=1003232, name="PeakForest", flags={1})]
        PeakForest,
        #[term(cv=MS, accession=1003253, name="DIA-NN", flags={1})]
        DIANN,
        #[term(cv=MS, accession=1003281, name="Casanovo", flags={1})]
        Casanovo,
        #[term(cv=MS, accession=1003309, name="Goslin", flags={2})]
        Goslin,
        #[term(cv=MS, accession=1003357, name="ANN-SoLo", flags={1})]
        ANNSoLo,
        #[term(cv=MS, accession=1003376, name="ChemClipse", flags={3})]
        ChemClipse,
        #[term(cv=MS, accession=1003377, name="OpenChrom", flags={3})]
        OpenChrom,
        #[term(cv=MS, accession=1003382, name="waters_connect", flags={7})]
        WatersConnect,
        #[term(cv=MS, accession=1003386, name="Spectra", flags={3})]
        Spectra,
        #[term(cv=MS, accession=1003387, name="MetaboAnnotation", flags={3})]
        MetaboAnnotation,
        #[term(cv=MS, accession=1003388, name="CompoundDb", flags={3})]
        CompoundDb,
        #[term(cv=MS, accession=1003399, name="quality control software", flags={2})]
        QualityControlSoftware,
        #[term(cv=MS, accession=1003400, name="rmzqc", flags={0})]
        Rmzqc,
        #[term(cv=MS, accession=1003401, name="jmzqc", flags={0})]
        Jmzqc,
        #[term(cv=MS, accession=1003402, name="pymzqc", flags={0})]
        Pymzqc,
        #[term(cv=MS, accession=1003405, name="mzRecal", flags={2})]
        MzRecal,
        #[term(cv=MS, accession=4000151, name="MsQuality", flags={1})]
        MsQuality,
    }
    //[[[end]]] (checksum: cd5ca38fc06b4604402d65497de7b0a9)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn cvmap_test() {
        assert_eq!(SoftwareTerm::SCIEXTOFTOFSeriesExplorerSoftware.accession(), 1001483);
        assert_eq!(SoftwareTerm::SCIEXTOFTOFSeriesExplorerSoftware.name(), "SCIEX TOF/TOF Series Explorer Software");
        assert_eq!(SoftwareTerm::SCIEXTOFTOFSeriesExplorerSoftware.flags(), SoftwareType::Analysis | SoftwareType::Acquisition | SoftwareType::DataProcessing);
    }
}