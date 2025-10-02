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
        #[doc="software - Software related to the recording or transformation of spectra."]
        Software,
        #[term(cv=MS, accession=1000532, name="Xcalibur", flags={7}, parents={["MS:1000693", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        #[doc="Xcalibur - Thermo Finnigan software for data acquisition and analysis."]
        Xcalibur,
        #[term(cv=MS, accession=1000533, name="Bioworks", flags={3}, parents={["MS:1000693", "MS:1001456", "MS:1001457"]})]
        #[doc="Bioworks - Thermo Finnigan software for data analysis of peptides and proteins."]
        Bioworks,
        #[term(cv=MS, accession=1000534, name="MassLynx", flags={7}, parents={["MS:1000694", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        #[doc="MassLynx - Micromass software for data acquisition and analysis."]
        MassLynx,
        #[term(cv=MS, accession=1000535, name="FlexAnalysis", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        #[doc="FlexAnalysis - Bruker software for data analysis."]
        FlexAnalysis,
        #[term(cv=MS, accession=1000536, name="Data Explorer", flags={7}, parents={["MS:1000691", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        #[doc="Data Explorer - Applied Biosystems software for data acquisition and analysis."]
        DataExplorer,
        #[term(cv=MS, accession=1000537, name="4700 Explorer", flags={7}, parents={["MS:1000691", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        #[doc="4700 Explorer - Applied Biosystems software for data acquisition and analysis."]
        _4700Explorer,
        #[term(cv=MS, accession=1000538, name="massWolf", flags={2}, parents={["MS:1001457"]})]
        #[doc="massWolf - A software for converting Waters raw directory format to mzXML or mzML. MassWolf was originally developed at the Institute for Systems Biology."]
        MassWolf,
        #[term(cv=MS, accession=1000539, name="Voyager Biospectrometry Workstation System", flags={7}, parents={["MS:1000691", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        #[doc="Voyager Biospectrometry Workstation System - Applied Biosystems MALDI-TOF data acquisition and analysis system."]
        VoyagerBiospectrometryWorkstationSystem,
        #[term(cv=MS, accession=1000540, name="FlexControl", flags={4}, parents={["MS:1000692", "MS:1001455"]})]
        #[doc="FlexControl - Bruker software for data acquisition."]
        FlexControl,
        #[term(cv=MS, accession=1000541, name="ReAdW", flags={2}, parents={["MS:1001457"]})]
        #[doc="ReAdW - A software program for converting Thermo Finnigan RAW file format to mzXML or mzML. ReAdW was originally developed at the Institute for Systems Biology. Its whimsical interleaved spelling and capitalization is pronounced 'readraw'."]
        ReAdW,
        #[term(cv=MS, accession=1000542, name="MzStar", flags={2}, parents={["MS:1001457"]})]
        #[doc="MzStar - A software program for converting Applied Biosystems wiff file format to mzXML format. MzStar was originally developed at the Institute for Systems Biology. It is now obsoleted by the MzWiff program."]
        MzStar,
        #[term(cv=MS, accession=1000551, name="Analyst", flags={7}, parents={["MS:1000690", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        #[doc="Analyst - SCIEX or Applied Biosystems|MDS SCIEX software for data acquisition."]
        Analyst,
        #[term(cv=MS, accession=1000553, name="Trapper", flags={2}, parents={["MS:1001457"]})]
        #[doc="Trapper - A software program for converting Agilent MassHunter format to mzXML or mzML. Trapper was originally developed at the Institute for Systems Biology."]
        Trapper,
        #[term(cv=MS, accession=1000591, name="MzWiff", flags={2}, parents={["MS:1001457"]})]
        #[doc="MzWiff - A software program for converting Applied Biosystems wiff file format to the mzXML or mzML format. MzWiff is currently maintained at the Institute for Systems Biology. It replaces the slower mzStar program."]
        MzWiff,
        #[term(cv=MS, accession=1000600, name="Proteios", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="Proteios - Database application and analysis platform for proteomics."]
        Proteios,
        #[term(cv=MS, accession=1000601, name="ProteinLynx Global Server", flags={3}, parents={["MS:1000694", "MS:1001456", "MS:1001457"]})]
        #[doc="ProteinLynx Global Server - Waters software for data analysis."]
        ProteinLynxGlobalServer,
        #[term(cv=MS, accession=1000615, name="ProteoWizard software", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="ProteoWizard software - ProteoWizard software for data processing and analysis. Primarily developed by the labs of P. Malick and D. Tabb."]
        ProteoWizardSoftware,
        #[term(cv=MS, accession=1000650, name="Proteome Discoverer", flags={3}, parents={["MS:1000693", "MS:1001456", "MS:1001457"]})]
        #[doc="Proteome Discoverer - Thermo Scientific software for data analysis of peptides and proteins."]
        ProteomeDiscoverer,
        #[term(cv=MS, accession=1000659, name="4000 Series Explorer Software", flags={7}, parents={["MS:1000690", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        #[doc="4000 Series Explorer Software - SCIEX or Applied Biosystems software for data acquisition and analysis."]
        _4000SeriesExplorerSoftware,
        #[term(cv=MS, accession=1000661, name="GPS Explorer", flags={7}, parents={["MS:1000690", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        #[doc="GPS Explorer - SCIEX or Applied Biosystems software for data acquisition and analysis."]
        GPSExplorer,
        #[term(cv=MS, accession=1000662, name="LightSight Software", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        #[doc="LightSight Software - SCIEX or Applied Biosystems|MDS SCIEX software metabolite identification."]
        LightSightSoftware,
        #[term(cv=MS, accession=1000663, name="ProteinPilot Software", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        #[doc="ProteinPilot Software - SCIEX or Applied Biosystems|MDS SCIEX software for protein ID and quant."]
        ProteinPilotSoftware,
        #[term(cv=MS, accession=1000664, name="TissueView Software", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        #[doc="TissueView Software - Applied Biosystems|MDS SCIEX software for tissue imaging."]
        TissueViewSoftware,
        #[term(cv=MS, accession=1000665, name="MarkerView Software", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        #[doc="MarkerView Software - Applied Biosystems|MDS SCIEX software for metabolomics and biomarker profiling."]
        MarkerViewSoftware,
        #[term(cv=MS, accession=1000666, name="MRMPilot Software", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        #[doc="MRMPilot Software - Applied Biosystems|MDS SCIEX software for MRM assay development."]
        MRMPilotSoftware,
        #[term(cv=MS, accession=1000667, name="BioAnalyst", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        #[doc="BioAnalyst - Applied Biosystems|MDS SCIEX software for bio-related data exploration."]
        BioAnalyst,
        #[term(cv=MS, accession=1000668, name="Pro ID", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        #[doc="Pro ID - Applied Biosystems|MDS SCIEX software for protein identification."]
        ProID,
        #[term(cv=MS, accession=1000669, name="Pro ICAT", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        #[doc="Pro ICAT - Applied Biosystems|MDS SCIEX software for protein ID and quant by ICAT."]
        ProICAT,
        #[term(cv=MS, accession=1000670, name="Pro Quant", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        #[doc="Pro Quant - Applied Biosystems|MDS SCIEX software for protein ID and quant by iTRAQ."]
        ProQuant,
        #[term(cv=MS, accession=1000671, name="Pro BLAST", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        #[doc="Pro BLAST - Applied Biosystems|MDS SCIEX software for MS-BLAST identification."]
        ProBLAST,
        #[term(cv=MS, accession=1000672, name="Cliquid", flags={0}, parents={["MS:1000690"]})]
        #[doc="Cliquid - SCIEX Cliquid software for data analysis and quantitation."]
        Cliquid,
        #[term(cv=MS, accession=1000673, name="MIDAS Workflow Designer", flags={0}, parents={["MS:1000690"]})]
        #[doc="MIDAS Workflow Designer - Applied Biosystems|MDS SCIEX software for MRM assay development."]
        MIDASWorkflowDesigner,
        #[term(cv=MS, accession=1000674, name="MultiQuant", flags={3}, parents={["MS:1000690", "MS:1001456", "MS:1001457"]})]
        #[doc="MultiQuant - Applied Biosystems|MDS SCIEX software for MRM-based quantitation."]
        MultiQuant,
        #[term(cv=MS, accession=1000678, name="MassHunter Data Acquisition", flags={4}, parents={["MS:1000689", "MS:1001455"]})]
        #[doc="MassHunter Data Acquisition - Software for data acquisition of 6000 series instruments."]
        MassHunterDataAcquisition,
        #[term(cv=MS, accession=1000679, name="MassHunter Easy Access", flags={4}, parents={["MS:1000689", "MS:1001455"]})]
        #[doc="MassHunter Easy Access - Software for open access data acquisition."]
        MassHunterEasyAccess,
        #[term(cv=MS, accession=1000680, name="MassHunter Qualitative Analysis", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        #[doc="MassHunter Qualitative Analysis - Software for data analysis of data from 6000 series instruments."]
        MassHunterQualitativeAnalysis,
        #[term(cv=MS, accession=1000681, name="MassHunter Quantitative Analysis", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        #[doc="MassHunter Quantitative Analysis - Software for quantitation of Triple Quadrupole and Quadrupole Time-of-Flight data."]
        MassHunterQuantitativeAnalysis,
        #[term(cv=MS, accession=1000682, name="MassHunter Metabolite ID", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        #[doc="MassHunter Metabolite ID - Software for identification of metabolites."]
        MassHunterMetaboliteID,
        #[term(cv=MS, accession=1000683, name="MassHunter BioConfirm", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        #[doc="MassHunter BioConfirm - Software for protein characterization."]
        MassHunterBioConfirm,
        #[term(cv=MS, accession=1000684, name="Genespring MS", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        #[doc="Genespring MS - Software for quantitation and statistical analysis of TOF and Q-TOF LC/MS data."]
        GenespringMS,
        #[term(cv=MS, accession=1000685, name="MassHunter Mass Profiler", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        #[doc="MassHunter Mass Profiler - Software for quantitation and statistical analysis of TOF and Q-TOF LC/MS data."]
        MassHunterMassProfiler,
        #[term(cv=MS, accession=1000686, name="METLIN", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        #[doc="METLIN - Personal Metabolite Database for MassHunter Workstation. Software for identification of human metabolites."]
        METLIN,
        #[term(cv=MS, accession=1000687, name="Spectrum Mill for MassHunter Workstation", flags={3}, parents={["MS:1000689", "MS:1001456", "MS:1001457"]})]
        #[doc="Spectrum Mill for MassHunter Workstation - Software for protein identification and characterization of complex protein digest mixtures."]
        SpectrumMillForMassHunterWorkstation,
        #[term(cv=MS, accession=1000688, name="6300 Series Ion Trap Data Analysis Software", flags={7}, parents={["MS:1000689", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        #[doc="6300 Series Ion Trap Data Analysis Software - Software for data analysis of 6300 series ion trap mass spectrometers."]
        _6300SeriesIonTrapDataAnalysisSoftware,
        #[term(cv=MS, accession=1000689, name="Agilent software", flags={0}, parents={["MS:1000531"]})]
        #[doc="Agilent software - Agilent software for data acquisition and analysis."]
        AgilentSoftware,
        #[term(cv=MS, accession=1000690, name="SCIEX software", flags={0}, parents={["MS:1000531"]})]
        #[doc="SCIEX software - SCIEX or Applied Biosystems software for data acquisition and analysis."]
        SCIEXSoftware,
        #[term(cv=MS, accession=1000691, name="Applied Biosystems software", flags={0}, parents={["MS:1000531"]})]
        #[doc="Applied Biosystems software - Applied Biosystems|MDS SCIEX software for data acquisition and analysis."]
        AppliedBiosystemsSoftware,
        #[term(cv=MS, accession=1000692, name="Bruker software", flags={0}, parents={["MS:1000531"]})]
        #[doc="Bruker software - Bruker software for data acquisition and analysis."]
        BrukerSoftware,
        #[term(cv=MS, accession=1000693, name="Thermo Finnigan software", flags={0}, parents={["MS:1000531"]})]
        #[doc="Thermo Finnigan software - Thermo Finnigan software for data acquisition and analysis."]
        ThermoFinniganSoftware,
        #[term(cv=MS, accession=1000694, name="Waters software", flags={0}, parents={["MS:1000531"]})]
        #[doc="Waters software - Waters software for data acquisition and analysis."]
        WatersSoftware,
        #[term(cv=MS, accession=1000706, name="apexControl", flags={4}, parents={["MS:1000692", "MS:1001455"]})]
        #[doc="apexControl - Bruker software for data acquisition."]
        ApexControl,
        #[term(cv=MS, accession=1000707, name="BioTools", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        #[doc="BioTools - Bruker software for data analysis."]
        BioTools,
        #[term(cv=MS, accession=1000708, name="CLINPROT", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        #[doc="CLINPROT - Bruker CLINPROT software."]
        CLINPROT,
        #[term(cv=MS, accession=1000709, name="CLINPROT micro", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        #[doc="CLINPROT micro - Bruker CLINPROT micro software."]
        CLINPROTMicro,
        #[term(cv=MS, accession=1000710, name="CLINPROT robot", flags={0}, parents={["MS:1000692"]})]
        #[doc="CLINPROT robot - Bruker CLINPROT robot software."]
        CLINPROTRobot,
        #[term(cv=MS, accession=1000711, name="ClinProTools", flags={0}, parents={["MS:1000692"]})]
        #[doc="ClinProTools - Bruker ClinProTools software."]
        ClinProTools,
        #[term(cv=MS, accession=1000712, name="Compass", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        #[doc="Compass - Bruker Compass software."]
        Compass,
        #[term(cv=MS, accession=1000713, name="Compass for HCT/esquire", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        #[doc="Compass for HCT/esquire - Bruker Compass for HCT/esquire software."]
        CompassForHCTEsquire,
        #[term(cv=MS, accession=1000714, name="Compass for micrOTOF", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        #[doc="Compass for micrOTOF - Bruker Compass for micrOTOF software."]
        CompassForMicrOTOF,
        #[term(cv=MS, accession=1000715, name="Compass OpenAccess", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        #[doc="Compass OpenAccess - Bruker compass OpenAccess software."]
        CompassOpenAccess,
        #[term(cv=MS, accession=1000716, name="Compass Security Pack", flags={0}, parents={["MS:1000692"]})]
        #[doc="Compass Security Pack - Bruker compass Security Pack software."]
        CompassSecurityPack,
        #[term(cv=MS, accession=1000717, name="CompassXport", flags={2}, parents={["MS:1000692", "MS:1001457"]})]
        #[doc="CompassXport - Bruker stand-alone software for data conversion."]
        CompassXport,
        #[term(cv=MS, accession=1000718, name="CompassXtract", flags={2}, parents={["MS:1000692", "MS:1001457"]})]
        #[doc="CompassXtract - Bruker software library for data access."]
        CompassXtract,
        #[term(cv=MS, accession=1000719, name="DataAnalysis", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        #[doc="DataAnalysis - Bruker software for data analysis."]
        DataAnalysis,
        #[term(cv=MS, accession=1000720, name="dpControl", flags={4}, parents={["MS:1000692", "MS:1001455"]})]
        #[doc="dpControl - Bruker software for data acquisition."]
        DpControl,
        #[term(cv=MS, accession=1000721, name="esquireControl", flags={4}, parents={["MS:1000692", "MS:1001455"]})]
        #[doc="esquireControl - Bruker software for data acquisition."]
        EsquireControl,
        #[term(cv=MS, accession=1000722, name="flexImaging", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        #[doc="flexImaging - Bruker software for data analysis."]
        FlexImaging,
        #[term(cv=MS, accession=1000723, name="GENOLINK", flags={0}, parents={["MS:1000692"]})]
        #[doc="GENOLINK - Bruker GENOLINK software."]
        GENOLINK,
        #[term(cv=MS, accession=1000724, name="GenoTools", flags={0}, parents={["MS:1000692"]})]
        #[doc="GenoTools - Bruker GenoTools software."]
        GenoTools,
        #[term(cv=MS, accession=1000725, name="HCTcontrol", flags={4}, parents={["MS:1000692", "MS:1001455"]})]
        #[doc="HCTcontrol - Bruker software for data acquisition."]
        HCTcontrol,
        #[term(cv=MS, accession=1000726, name="micrOTOFcontrol", flags={4}, parents={["MS:1000692", "MS:1001455"]})]
        #[doc="micrOTOFcontrol - Bruker software for data acquisition."]
        MicrOTOFcontrol,
        #[term(cv=MS, accession=1000727, name="PolyTools", flags={0}, parents={["MS:1000692"]})]
        #[doc="PolyTools - Bruker PolyTools software."]
        PolyTools,
        #[term(cv=MS, accession=1000728, name="ProfileAnalysis", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        #[doc="ProfileAnalysis - Bruker software for data analysis."]
        ProfileAnalysis,
        #[term(cv=MS, accession=1000729, name="PROTEINEER", flags={0}, parents={["MS:1000692"]})]
        #[doc="PROTEINEER - Bruker PROTEINEER software."]
        PROTEINEER,
        #[term(cv=MS, accession=1000730, name="PROTEINEER dp", flags={0}, parents={["MS:1000692"]})]
        #[doc="PROTEINEER dp - Bruker PROTEINEER dp software."]
        PROTEINEERDp,
        #[term(cv=MS, accession=1000731, name="PROTEINEER fc", flags={0}, parents={["MS:1000692"]})]
        #[doc="PROTEINEER fc - Bruker PROTEINEER fc software."]
        PROTEINEERFc,
        #[term(cv=MS, accession=1000732, name="PROTEINEER spII", flags={0}, parents={["MS:1000692"]})]
        #[doc="PROTEINEER spII - Bruker PROTEINEER spII software."]
        PROTEINEERSpII,
        #[term(cv=MS, accession=1000733, name="PROTEINEER-LC", flags={0}, parents={["MS:1000692"]})]
        #[doc="PROTEINEER-LC - Bruker PROTEINEER-LC software."]
        PROTEINEERLC,
        #[term(cv=MS, accession=1000734, name="ProteinScape", flags={1}, parents={["MS:1000692", "MS:1001456"]})]
        #[doc="ProteinScape - Bruker ProteinScape software."]
        ProteinScape,
        #[term(cv=MS, accession=1000735, name="PureDisk", flags={0}, parents={["MS:1000692"]})]
        #[doc="PureDisk - BrukerPureDisk software."]
        PureDisk,
        #[term(cv=MS, accession=1000736, name="QuantAnalysis", flags={3}, parents={["MS:1000692", "MS:1001456", "MS:1001457"]})]
        #[doc="QuantAnalysis - Bruker software for data analysis."]
        QuantAnalysis,
        #[term(cv=MS, accession=1000737, name="spControl", flags={4}, parents={["MS:1000692", "MS:1001455"]})]
        #[doc="spControl - Bruker software for data acquisition."]
        SpControl,
        #[term(cv=MS, accession=1000738, name="TargetAnalysis", flags={0}, parents={["MS:1000692"]})]
        #[doc="TargetAnalysis - Bruker TargetAnalysis software."]
        TargetAnalysis,
        #[term(cv=MS, accession=1000739, name="WARP-LC", flags={0}, parents={["MS:1000692", "MS:1001139"]})]
        #[doc="WARP-LC - Bruker WARP-LC software."]
        WARPLC,
        #[term(cv=MS, accession=1000752, name="TOPP software", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="TOPP software - TOPP (The OpenMS proteomics pipeline) software."]
        TOPPSoftware,
        #[term(cv=MS, accession=1000753, name="BaselineFilter", flags={0}, parents={["MS:1000752"]})]
        #[doc="BaselineFilter - Removes the baseline from profile spectra using a top-hat filter."]
        BaselineFilter,
        #[term(cv=MS, accession=1000754, name="DBExporter", flags={0}, parents={["MS:1000752"]})]
        #[doc="DBExporter - Exports data from an OpenMS database to a file."]
        DBExporter,
        #[term(cv=MS, accession=1000755, name="DBImporter", flags={0}, parents={["MS:1000752"]})]
        #[doc="DBImporter - Imports data to an OpenMS database."]
        DBImporter,
        #[term(cv=MS, accession=1000756, name="FileConverter", flags={0}, parents={["MS:1000752"]})]
        #[doc="FileConverter - Converts between different MS file formats."]
        FileConverter,
        #[term(cv=MS, accession=1000757, name="FileFilter", flags={0}, parents={["MS:1000752"]})]
        #[doc="FileFilter - Extracts or manipulates portions of data from peak, feature or consensus feature files."]
        FileFilter,
        #[term(cv=MS, accession=1000758, name="FileMerger", flags={0}, parents={["MS:1000752"]})]
        #[doc="FileMerger - Merges several MS files into one file."]
        FileMerger,
        #[term(cv=MS, accession=1000759, name="InternalCalibration", flags={0}, parents={["MS:1000752"]})]
        #[doc="InternalCalibration - Applies an internal calibration."]
        InternalCalibration,
        #[term(cv=MS, accession=1000760, name="MapAligner", flags={0}, parents={["MS:1000752"]})]
        #[doc="MapAligner - OBSOLETE Corrects retention time distortions between maps."]
        MapAligner,
        #[term(cv=MS, accession=1000761, name="MapNormalizer", flags={0}, parents={["MS:1000752"]})]
        #[doc="MapNormalizer - Normalizes peak intensities in an MS run."]
        MapNormalizer,
        #[term(cv=MS, accession=1000762, name="NoiseFilter", flags={0}, parents={["MS:1000752"]})]
        #[doc="NoiseFilter - OBSOLETE Removes noise from profile spectra by using different smoothing techniques."]
        NoiseFilter,
        #[term(cv=MS, accession=1000763, name="PeakPicker", flags={0}, parents={["MS:1000752"]})]
        #[doc="PeakPicker - OBSOLETE Finds mass spectrometric peaks in profile mass spectra."]
        PeakPicker,
        #[term(cv=MS, accession=1000764, name="Resampler", flags={0}, parents={["MS:1000752"]})]
        #[doc="Resampler - Transforms an LC/MS map into a resampled map or a png image."]
        Resampler,
        #[term(cv=MS, accession=1000765, name="SpectraFilter", flags={0}, parents={["MS:1000752"]})]
        #[doc="SpectraFilter - OBSOLETE Applies a filter to peak spectra."]
        SpectraFilter,
        #[term(cv=MS, accession=1000766, name="TOFCalibration", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOFCalibration - Applies time of flight calibration."]
        TOFCalibration,
        #[term(cv=MS, accession=1000799, name="custom unreleased software tool", flags={0}, parents={["MS:1000531"]})]
        #[doc="custom unreleased software tool - A software tool that has not yet been released. The value should describe the software. Please do not use this term for publicly available software - contact the PSI-MS working group in order to have another CV term added."]
        CustomUnreleasedSoftwareTool,
        #[term(cv=MS, accession=1000817, name="HyStar", flags={0}, parents={["MS:1000692"]})]
        #[doc="HyStar - Bruker software for hyphenated experiments."]
        HyStar,
        #[term(cv=MS, accession=1000871, name="SRM software", flags={0}, parents={["MS:1000531"]})]
        #[doc="SRM software - Software used to predict, select, or optimize transitions or analyze the results of selected reaction monitoring runs."]
        SRMSoftware,
        #[term(cv=MS, accession=1000872, name="MaRiMba", flags={0}, parents={["MS:1000871"]})]
        #[doc="MaRiMba - Software used to predict transitions for selected reaction monitoring experiments based on observed spectrum libraries developed and distributed by the Institute for Systems Biology."]
        MaRiMba,
        #[term(cv=MS, accession=1000873, name="peptide attribute calculation software", flags={0}, parents={["MS:1000531"]})]
        #[doc="peptide attribute calculation software - Software used to predict or calculate numerical attributes of peptides."]
        PeptideAttributeCalculationSoftware,
        #[term(cv=MS, accession=1000874, name="SSRCalc", flags={0}, parents={["MS:1000873"]})]
        #[doc="SSRCalc - Sequence Specific Retention Calculator estimates the retention time of peptides based on their sequence."]
        SSRCalc,
        #[term(cv=MS, accession=1000922, name="Skyline", flags={0}, parents={["MS:1000871", "MS:1001139"]})]
        #[doc="Skyline - Software used to predict, select, and optimize transitions as well as analyze the results of selected reaction monitoring runs developed and distributed by the MacCoss lab at the University of Washington."]
        Skyline,
        #[term(cv=MS, accession=1000923, name="TIQAM", flags={0}, parents={["MS:1000871"]})]
        #[doc="TIQAM - Software used to predict, select, and optimize transitions for selected reaction monitoring experiments developed and distributed by the Institute for Systems Biology."]
        TIQAM,
        #[term(cv=MS, accession=1000925, name="ATAQS", flags={0}, parents={["MS:1000871"]})]
        #[doc="ATAQS - Software suite used to predict, select, and optimize transitions as well as analyze the results of selected reaction monitoring runs developed and distributed by the Institute for Systems Biology."]
        ATAQS,
        #[term(cv=MS, accession=1001139, name="quantitation software name", flags={0}, parents={["MS:1000531", "MS:1001129"]})]
        #[doc="quantitation software name - Quantitation software name."]
        QuantitationSoftwareName,
        #[term(cv=MS, accession=1001207, name="Mascot", flags={1}, parents={["MS:1001456"]})]
        #[doc="Mascot - The name of the Mascot search engine."]
        Mascot,
        #[term(cv=MS, accession=1001208, name="SEQUEST", flags={1}, parents={["MS:1001456"]})]
        #[doc="SEQUEST - The name of the SEQUEST search engine."]
        SEQUEST,
        #[term(cv=MS, accession=1001209, name="Phenyx", flags={1}, parents={["MS:1001456"]})]
        #[doc="Phenyx - The name of the Phenyx search engine."]
        Phenyx,
        #[term(cv=MS, accession=1001327, name="Spectronaut", flags={1}, parents={["MS:1001456", "MS:1003207"]})]
        #[doc="Spectronaut - Commercial cross-vendor software for library (peptide centric), and library-free (spectrum centric) analysis and quantification of DIA data."]
        Spectronaut,
        #[term(cv=MS, accession=1001455, name="acquisition software", flags={0}, parents={["MS:1000531"]})]
        #[doc="acquisition software - Acquisition software."]
        AcquisitionSoftware,
        #[term(cv=MS, accession=1001456, name="analysis software", flags={0}, parents={["MS:1000531"]})]
        #[doc="analysis software - Analysis software."]
        AnalysisSoftware,
        #[term(cv=MS, accession=1001457, name="data processing software", flags={0}, parents={["MS:1000531"]})]
        #[doc="data processing software - Data processing software."]
        DataProcessingSoftware,
        #[term(cv=MS, accession=1001461, name="greylag", flags={1}, parents={["MS:1001456"]})]
        #[doc="greylag - Greylag identification software."]
        Greylag,
        #[term(cv=MS, accession=1001475, name="OMSSA", flags={1}, parents={["MS:1001456"]})]
        #[doc="OMSSA - Open Mass Spectrometry Search Algorithm was used to analyze the spectra."]
        OMSSA,
        #[term(cv=MS, accession=1001476, name="X!Tandem", flags={1}, parents={["MS:1001456"]})]
        #[doc="X!Tandem - X!Tandem was used to analyze the spectra."]
        XTandem,
        #[term(cv=MS, accession=1001477, name="SpectraST", flags={1}, parents={["MS:1001456", "MS:1003207", "MS:1003406"]})]
        #[doc="SpectraST - Open-source software for mass spectral library creation and searching, developed at the Institute for Systems Biology and the Hong Kong University of Science and Technology. Part of the Trans-Proteomic Pipeline."]
        SpectraST,
        #[term(cv=MS, accession=1001478, name="Mascot Parser", flags={1}, parents={["MS:1001456"]})]
        #[doc="Mascot Parser - Mascot Parser was used to analyze the spectra."]
        MascotParser,
        #[term(cv=MS, accession=1001483, name="SCIEX TOF/TOF Series Explorer Software", flags={7}, parents={["MS:1000690", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        #[doc="SCIEX TOF/TOF Series Explorer Software - SCIEX or Applied Biosystems software for TOF/TOF data acquisition and analysis."]
        SCIEXTOFTOFSeriesExplorerSoftware,
        #[term(cv=MS, accession=1001487, name="ProteinExtractor", flags={1}, parents={["MS:1000692", "MS:1001456"]})]
        #[doc="ProteinExtractor - An algorithm for protein determination/assembly integrated into Bruker's ProteinScape."]
        ProteinExtractor,
        #[term(cv=MS, accession=1001488, name="Mascot Distiller", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        #[doc="Mascot Distiller - Mascot Distiller."]
        MascotDistiller,
        #[term(cv=MS, accession=1001489, name="Mascot Integra", flags={1}, parents={["MS:1001456"]})]
        #[doc="Mascot Integra - Mascot Integra."]
        MascotIntegra,
        #[term(cv=MS, accession=1001490, name="Percolator", flags={1}, parents={["MS:1001456"]})]
        #[doc="Percolator - Percolator."]
        Percolator,
        #[term(cv=MS, accession=1001557, name="Shimadzu Corporation software", flags={0}, parents={["MS:1000531"]})]
        #[doc="Shimadzu Corporation software - Shimadzu Corporation software."]
        ShimadzuCorporationSoftware,
        #[term(cv=MS, accession=1001558, name="MALDI Solutions", flags={7}, parents={["MS:1001455", "MS:1001456", "MS:1001457", "MS:1001557"]})]
        #[doc="MALDI Solutions - Shimadzu Biotech software for data acquisition, processing, and analysis."]
        MALDISolutions,
        #[term(cv=MS, accession=1001561, name="Scaffold", flags={1}, parents={["MS:1001456"]})]
        #[doc="Scaffold - Scaffold analysis software."]
        Scaffold,
        #[term(cv=MS, accession=1001582, name="XCMS", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="XCMS - Bioconductor package XCMS for preprocessing high-throughput, untargeted analyte profiling data."]
        XCMS,
        #[term(cv=MS, accession=1001583, name="MaxQuant", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        #[doc="MaxQuant - MaxQuant is a quantitative proteomics software package designed for analyzing large mass spectrometric data sets. It is specifically aimed at high resolution MS data."]
        MaxQuant,
        #[term(cv=MS, accession=1001585, name="MyriMatch", flags={1}, parents={["MS:1001456"]})]
        #[doc="MyriMatch - Tabb Lab software for directly comparing peptides in a database to tandem mass spectra."]
        MyriMatch,
        #[term(cv=MS, accession=1001586, name="DirecTag", flags={1}, parents={["MS:1001456"]})]
        #[doc="DirecTag - Tabb Lab software for generating sequence tags from tandem mass spectra."]
        DirecTag,
        #[term(cv=MS, accession=1001587, name="TagRecon", flags={1}, parents={["MS:1001456"]})]
        #[doc="TagRecon - Tabb Lab software for reconciling sequence tags to a protein database."]
        TagRecon,
        #[term(cv=MS, accession=1001588, name="Pepitome", flags={1}, parents={["MS:1001456"]})]
        #[doc="Pepitome - Tabb Lab software for spectral library searches on tandem mass spectra."]
        Pepitome,
        #[term(cv=MS, accession=1001795, name="Empower", flags={3}, parents={["MS:1000694", "MS:1001456", "MS:1001457"]})]
        #[doc="Empower - Waters Empower software for liquid chromatography and mass spectrometry acquisition."]
        Empower,
        #[term(cv=MS, accession=1001796, name="UNIFY", flags={3}, parents={["MS:1000694", "MS:1001456", "MS:1001457"]})]
        #[doc="UNIFY - Waters UNIFY software for liquid chromatography and mass spectrometry acquisition."]
        UNIFY,
        #[term(cv=MS, accession=1001798, name="LECO software", flags={0}, parents={["MS:1000531"]})]
        #[doc="LECO software - LECO software for data acquisition and analysis."]
        LECOSoftware,
        #[term(cv=MS, accession=1001799, name="ChromaTOF software", flags={7}, parents={["MS:1001455", "MS:1001456", "MS:1001457", "MS:1001798"]})]
        #[doc="ChromaTOF software - Software for acquisition, processing and analysis of data for LECO instruments."]
        ChromaTOFSoftware,
        #[term(cv=MS, accession=1001830, name="Progenesis LC-MS", flags={0}, parents={["MS:1001139"]})]
        #[doc="Progenesis LC-MS - Software from Nonlinear Dynamics for LC-MS label-free workflow."]
        ProgenesisLCMS,
        #[term(cv=MS, accession=1001831, name="SILACAnalyzer", flags={0}, parents={["MS:1001139", "MS:1000752"]})]
        #[doc="SILACAnalyzer - Software for SILAC workflow."]
        SILACAnalyzer,
        #[term(cv=MS, accession=1001877, name="ChromaTOF HRT software", flags={7}, parents={["MS:1001455", "MS:1001456", "MS:1001457", "MS:1001798"]})]
        #[doc="ChromaTOF HRT software - Software for acquisition, processing and analysis of data for LECO instruments."]
        ChromaTOFHRTSoftware,
        #[term(cv=MS, accession=1001878, name="MALDI Solutions Microbial Identification", flags={0}, parents={["MS:1001558"]})]
        #[doc="MALDI Solutions Microbial Identification - Shimadzu Biotech software for data acquisition, processing, and analysis."]
        MALDISolutionsMicrobialIdentification,
        #[term(cv=MS, accession=1001886, name="SQID", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="SQID - Software for data analysis of peptides and proteins."]
        SQID,
        #[term(cv=MS, accession=1001912, name="PinPoint", flags={3}, parents={["MS:1000693", "MS:1001456", "MS:1001457"]})]
        #[doc="PinPoint - Thermo Scientific PinPoint SRM analysis software."]
        PinPoint,
        #[term(cv=MS, accession=1001914, name="pymzML", flags={2}, parents={["MS:1001457"]})]
        #[doc="pymzML - Python module to interface mzML Data."]
        PymzML,
        #[term(cv=MS, accession=1001946, name="PEAKS Studio", flags={3}, parents={["MS:1001139", "MS:1001456", "MS:1001457"]})]
        #[doc="PEAKS Studio - PEAKS Studio software for data analysis."]
        PEAKSStudio,
        #[term(cv=MS, accession=1001947, name="PEAKS Online", flags={3}, parents={["MS:1001139", "MS:1001456", "MS:1001457"]})]
        #[doc="PEAKS Online - PEAKS Online software for high throughput data analysis."]
        PEAKSOnline,
        #[term(cv=MS, accession=1001948, name="PEAKS Node", flags={3}, parents={["MS:1001139", "MS:1001456", "MS:1001457"]})]
        #[doc="PEAKS Node - PEAKS Node software for high throughput data analysis."]
        PEAKSNode,
        #[term(cv=MS, accession=1001949, name="BSI software", flags={0}, parents={["MS:1000531"]})]
        #[doc="BSI software - Bioinformatics Solutions Inc. Software for data processing and analysis."]
        BSISoftware,
        #[term(cv=MS, accession=1001973, name="DeBunker", flags={1}, parents={["MS:1001456"]})]
        #[doc="DeBunker - DeBunker software."]
        DeBunker,
        #[term(cv=MS, accession=1001977, name="MSQuant", flags={1}, parents={["MS:1001456"]})]
        #[doc="MSQuant - MSQuant software."]
        MSQuant,
        #[term(cv=MS, accession=1001984, name="Ascore software", flags={1}, parents={["MS:1001456"]})]
        #[doc="Ascore software - Ascore software."]
        AscoreSoftware,
        #[term(cv=MS, accession=1002043, name="ProteinProspector", flags={1}, parents={["MS:1001456"]})]
        #[doc="ProteinProspector - ProteinProspector software for data acquisition and analysis."]
        ProteinProspector,
        #[term(cv=MS, accession=1002047, name="MS-GF", flags={1}, parents={["MS:1001456"]})]
        #[doc="MS-GF - MS-GF software used to re-score the peptide-spectrum matches."]
        MSGF,
        #[term(cv=MS, accession=1002048, name="MS-GF+", flags={1}, parents={["MS:1001456"]})]
        #[doc="MS-GF+ - MS-GF+ software used to analyze the spectra."]
        MSGFplus,
        #[term(cv=MS, accession=1002059, name="Microsoft Excel", flags={0}, parents={["MS:1001139"]})]
        #[doc="Microsoft Excel - Microsoft Excel (can be used for spectral counting)."]
        MicrosoftExcel,
        #[term(cv=MS, accession=1002063, name="FindPairs", flags={0}, parents={["MS:1001139"]})]
        #[doc="FindPairs - Software e.g. for SILAC and 14N/15N workflow, part of the PeakQuant suite."]
        FindPairs,
        #[term(cv=MS, accession=1002076, name="PAnalyzer", flags={1}, parents={["MS:1001456"]})]
        #[doc="PAnalyzer - PAnalyzer software for getting protein evidence categories."]
        PAnalyzer,
        #[term(cv=MS, accession=1002123, name="x-Tracker", flags={0}, parents={["MS:1001139"]})]
        #[doc="x-Tracker - X-Tracker generic tool for quantitative proteomics."]
        XTracker,
        #[term(cv=MS, accession=1002124, name="ProteoSuite", flags={0}, parents={["MS:1001139"]})]
        #[doc="ProteoSuite - ProteoSuite software for the analysis of quantitative proteomics data."]
        ProteoSuite,
        #[term(cv=MS, accession=1002129, name="ITRAQAnalyzer", flags={0}, parents={["MS:1001139", "MS:1000752"]})]
        #[doc="ITRAQAnalyzer - Software for iTRAQ workflow. Extracts and normalizes iTRAQ information from an MS experiment."]
        ITRAQAnalyzer,
        #[term(cv=MS, accession=1002131, name="TOPP noise filter", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP noise filter - Noise filter component of the TOPP software."]
        TOPPNoiseFilter,
        #[term(cv=MS, accession=1002132, name="TOPP NoiseFilterGaussian", flags={0}, parents={["MS:1002131"]})]
        #[doc="TOPP NoiseFilterGaussian - Removes noise from profile spectra by using a gaussian smoothing."]
        TOPPNoiseFilterGaussian,
        #[term(cv=MS, accession=1002133, name="TOPP NoiseFilterSGolay", flags={0}, parents={["MS:1002131"]})]
        #[doc="TOPP NoiseFilterSGolay - Removes noise from profile spectra by using a Savitzky-Golay smoothing."]
        TOPPNoiseFilterSGolay,
        #[term(cv=MS, accession=1002134, name="TOPP peak picker", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP peak picker - Peak picker component of the TOPP software."]
        TOPPPeakPicker,
        #[term(cv=MS, accession=1002135, name="TOPP PeakPickerHiRes", flags={0}, parents={["MS:1002134"]})]
        #[doc="TOPP PeakPickerHiRes - Finds mass spectrometric peaks in high-resoluted profile mass spectra."]
        TOPPPeakPickerHiRes,
        #[term(cv=MS, accession=1002136, name="TOPP PeakPickerWavelet", flags={0}, parents={["MS:1002134"]})]
        #[doc="TOPP PeakPickerWavelet - Finds mass spectrometric peaks with a wavelet algorithm in low-resoluted profile mass spectra."]
        TOPPPeakPickerWavelet,
        #[term(cv=MS, accession=1002137, name="TOPP spectra filter", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP spectra filter - Spectra filter component of the TOPP software."]
        TOPPSpectraFilter,
        #[term(cv=MS, accession=1002138, name="TOPP SpectraFilterBernNorm", flags={0}, parents={["MS:1002137"]})]
        #[doc="TOPP SpectraFilterBernNorm - Applies a Bern et al normalization to peak spectra."]
        TOPPSpectraFilterBernNorm,
        #[term(cv=MS, accession=1002139, name="TOPP SpectraFilterMarkerMower", flags={0}, parents={["MS:1002137"]})]
        #[doc="TOPP SpectraFilterMarkerMower - Applies a filter to peak spectra for marked peaks."]
        TOPPSpectraFilterMarkerMower,
        #[term(cv=MS, accession=1002140, name="TOPP SpectraFilterNLargest", flags={0}, parents={["MS:1002137"]})]
        #[doc="TOPP SpectraFilterNLargest - Retains the n largest peaks of a peak spectra."]
        TOPPSpectraFilterNLargest,
        #[term(cv=MS, accession=1002141, name="TOPP SpectraFilterNormalizer", flags={0}, parents={["MS:1002137"]})]
        #[doc="TOPP SpectraFilterNormalizer - Applies a TIC/maximal intensity normalization to peak spectra."]
        TOPPSpectraFilterNormalizer,
        #[term(cv=MS, accession=1002142, name="TOPP SpectraFilterParentPeakMower", flags={0}, parents={["MS:1002137"]})]
        #[doc="TOPP SpectraFilterParentPeakMower - Filters putative unfragmented precursor ions from tandem spectra."]
        TOPPSpectraFilterParentPeakMower,
        #[term(cv=MS, accession=1002143, name="TOPP SpectraFilterScaler", flags={0}, parents={["MS:1002137"]})]
        #[doc="TOPP SpectraFilterScaler - Applies a filter to peak spectra after intensity scaling according to rank."]
        TOPPSpectraFilterScaler,
        #[term(cv=MS, accession=1002144, name="TOPP SpectraFilterSqrtMower", flags={0}, parents={["MS:1002137"]})]
        #[doc="TOPP SpectraFilterSqrtMower - Applies a filter to peak spectra after intensity scaling to the square root."]
        TOPPSpectraFilterSqrtMower,
        #[term(cv=MS, accession=1002145, name="TOPP SpectraFilterThresholdMower", flags={0}, parents={["MS:1002137"]})]
        #[doc="TOPP SpectraFilterThresholdMower - Applies a filter of peaks below a given threshold to peak spectra."]
        TOPPSpectraFilterThresholdMower,
        #[term(cv=MS, accession=1002146, name="TOPP SpectraFilterWindowMower", flags={0}, parents={["MS:1002137"]})]
        #[doc="TOPP SpectraFilterWindowMower - Applies a filter of the largest peaks in a sliding window over a peak spectrum."]
        TOPPSpectraFilterWindowMower,
        #[term(cv=MS, accession=1002147, name="TOPP map aligner", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP map aligner - Map aligner component of the TOPP software."]
        TOPPMapAligner,
        #[term(cv=MS, accession=1002148, name="TOPP MapAlignerIdentification", flags={0}, parents={["MS:1002147"]})]
        #[doc="TOPP MapAlignerIdentification - Corrects retention time distortions between maps based on common peptide identifications."]
        TOPPMapAlignerIdentification,
        #[term(cv=MS, accession=1002149, name="TOPP MapAlignerPoseClustering", flags={0}, parents={["MS:1002147"]})]
        #[doc="TOPP MapAlignerPoseClustering - Corrects retention time distortions between maps using a pose clustering approach."]
        TOPPMapAlignerPoseClustering,
        #[term(cv=MS, accession=1002150, name="TOPP MapAlignerSpectrum", flags={0}, parents={["MS:1002147"]})]
        #[doc="TOPP MapAlignerSpectrum - Corrects retention time distortions between maps by spectrum alignment."]
        TOPPMapAlignerSpectrum,
        #[term(cv=MS, accession=1002154, name="TOPP DTAExtractor", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP DTAExtractor - Extracts spectra of an MS run file to several files in DTA format."]
        TOPPDTAExtractor,
        #[term(cv=MS, accession=1002155, name="TOPP IDMerger", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP IDMerger - Merges several protein/peptide identification files into one file."]
        TOPPIDMerger,
        #[term(cv=MS, accession=1002156, name="TOPP IDFileConverter", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP IDFileConverter - Converts identification engine file formats."]
        TOPPIDFileConverter,
        #[term(cv=MS, accession=1002157, name="TOPP SpectraMerger", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP SpectraMerger - Merges spectra from an LC/MS map, either by precursor or by RT blocks."]
        TOPPSpectraMerger,
        #[term(cv=MS, accession=1002158, name="TOPP MzTabExporter", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP MzTabExporter - Exports various XML formats to an mzTab file."]
        TOPPMzTabExporter,
        #[term(cv=MS, accession=1002159, name="TOPP MassTraceExtractor", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP MassTraceExtractor - Annotates mass traces in centroided LC/MS maps."]
        TOPPMassTraceExtractor,
        #[term(cv=MS, accession=1002160, name="TOPP PrecursorMassCorrector", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP PrecursorMassCorrector - Correct the precursor entries of tandem MS scans."]
        TOPPPrecursorMassCorrector,
        #[term(cv=MS, accession=1002161, name="TOPP HighResPrecursorMassCorrector", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP HighResPrecursorMassCorrector - Performs precursor mz correction on centroided high resolution data."]
        TOPPHighResPrecursorMassCorrector,
        #[term(cv=MS, accession=1002162, name="TOPP AdditiveSeries", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP AdditiveSeries - Computes an additive series to quantify a peptide in a set of samples."]
        TOPPAdditiveSeries,
        #[term(cv=MS, accession=1002163, name="TOPP Decharger", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP Decharger - Decharges and merges different feature charge variants of the same chemical entity."]
        TOPPDecharger,
        #[term(cv=MS, accession=1002164, name="TOPP EICExtractor", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP EICExtractor - Quantifies signals at given positions in (raw or picked) LC/MS maps."]
        TOPPEICExtractor,
        #[term(cv=MS, accession=1002165, name="TOPP feature finder", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP feature finder - Feature finder component of the TOPP software."]
        TOPPFeatureFinder,
        #[term(cv=MS, accession=1002166, name="TOPP FeatureFinderCentroided", flags={0}, parents={["MS:1002165"]})]
        #[doc="TOPP FeatureFinderCentroided - Detects two-dimensional features in centroided LC-MS data."]
        TOPPFeatureFinderCentroided,
        #[term(cv=MS, accession=1002167, name="TOPP FeatureFinderRaw", flags={0}, parents={["MS:1002165"]})]
        #[doc="TOPP FeatureFinderRaw - Detects two-dimensional features in uncentroided LC-MS data."]
        TOPPFeatureFinderRaw,
        #[term(cv=MS, accession=1002168, name="TOPP FeatureFinderIsotopeWavelet", flags={0}, parents={["MS:1002165"]})]
        #[doc="TOPP FeatureFinderIsotopeWavelet - Detects two-dimensional features in uncentroided LC-MS data with a wavelet algorithm."]
        TOPPFeatureFinderIsotopeWavelet,
        #[term(cv=MS, accession=1002169, name="TOPP FeatureFinderMetabo", flags={0}, parents={["MS:1002165"]})]
        #[doc="TOPP FeatureFinderMetabo - Detects two-dimensional features in centroided LC-MS data of metabolites."]
        TOPPFeatureFinderMetabo,
        #[term(cv=MS, accession=1002170, name="TOPP FeatureFinderMRM", flags={0}, parents={["MS:1002165"]})]
        #[doc="TOPP FeatureFinderMRM - Quantifies features LC-MS/MS MRM data."]
        TOPPFeatureFinderMRM,
        #[term(cv=MS, accession=1002171, name="TOPP ProteinQuantifier", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP ProteinQuantifier - Computes protein abundances from annotated feature/consensus maps."]
        TOPPProteinQuantifier,
        #[term(cv=MS, accession=1002172, name="TOPP ConsensusMapNormalizer", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP ConsensusMapNormalizer - Normalizes maps of one consensus XML file (after linking)."]
        TOPPConsensusMapNormalizer,
        #[term(cv=MS, accession=1002173, name="TOPP MapRTTransformer", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP MapRTTransformer - Applies retention time transformations to maps."]
        TOPPMapRTTransformer,
        #[term(cv=MS, accession=1002174, name="TOPP feature linker", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP feature linker - Feature linker component of the TOPP software."]
        TOPPFeatureLinker,
        #[term(cv=MS, accession=1002175, name="TOPP FeatureLinkerLabeled", flags={0}, parents={["MS:1002174"]})]
        #[doc="TOPP FeatureLinkerLabeled - Groups corresponding isotope-labeled features in a feature map."]
        TOPPFeatureLinkerLabeled,
        #[term(cv=MS, accession=1002176, name="TOPP FeatureLinkerUnlabeled", flags={0}, parents={["MS:1002174"]})]
        #[doc="TOPP FeatureLinkerUnlabeled - Groups corresponding features from multiple maps."]
        TOPPFeatureLinkerUnlabeled,
        #[term(cv=MS, accession=1002177, name="TOPP FeatureLinkerUnlabeledQT", flags={0}, parents={["MS:1002174"]})]
        #[doc="TOPP FeatureLinkerUnlabeledQT - Groups corresponding features from multiple maps using a quality threshold clustering approach."]
        TOPPFeatureLinkerUnlabeledQT,
        #[term(cv=MS, accession=1002178, name="TOPP CompNovo", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP CompNovo - Performs a peptide/protein identification with the CompNovo engine."]
        TOPPCompNovo,
        #[term(cv=MS, accession=1002179, name="TOPP CompNovoCID", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP CompNovoCID - Performs a peptide/protein identification with the CompNovo engine in collision-induced dissociation (CID) mode."]
        TOPPCompNovoCID,
        #[term(cv=MS, accession=1002180, name="TOPP software adaptor", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP software adaptor - Software adaptor to an external program in the TOPP software."]
        TOPPSoftwareAdaptor,
        #[term(cv=MS, accession=1002181, name="TOPP InspectAdapter", flags={0}, parents={["MS:1002180"]})]
        #[doc="TOPP InspectAdapter - Identifies MS2 spectra using the external program Inspect."]
        TOPPInspectAdapter,
        #[term(cv=MS, accession=1002182, name="TOPP MascotAdapter", flags={0}, parents={["MS:1002180"]})]
        #[doc="TOPP MascotAdapter - Identifies MS2 spectra using the external program Mascot."]
        TOPPMascotAdapter,
        #[term(cv=MS, accession=1002183, name="TOPP MascotAdapterOnline", flags={0}, parents={["MS:1002180"]})]
        #[doc="TOPP MascotAdapterOnline - Identifies MS2 spectra using the online version of the external program Mascot."]
        TOPPMascotAdapterOnline,
        #[term(cv=MS, accession=1002184, name="TOPP OMSSAAdapter", flags={0}, parents={["MS:1002180"]})]
        #[doc="TOPP OMSSAAdapter - Identifies MS2 spectra using the external program OMSSA."]
        TOPPOMSSAAdapter,
        #[term(cv=MS, accession=1002185, name="TOPP PepNovoAdapter", flags={0}, parents={["MS:1002180"]})]
        #[doc="TOPP PepNovoAdapter - Identifies MS2 spectra using the external program PepNovo."]
        TOPPPepNovoAdapter,
        #[term(cv=MS, accession=1002186, name="TOPP XTandemAdapter", flags={0}, parents={["MS:1002180"]})]
        #[doc="TOPP XTandemAdapter - Identifies MS2 spectra using the external program XTandem."]
        TOPPXTandemAdapter,
        #[term(cv=MS, accession=1002187, name="TOPP SpecLibSearcher", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP SpecLibSearcher - Identifies peptide MS2 spectra by spectral matching with a searchable spectral library."]
        TOPPSpecLibSearcher,
        #[term(cv=MS, accession=1002188, name="TOPP ConsensusID", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP ConsensusID - Computes a consensus identification from peptide identifications of several identification engines."]
        TOPPConsensusID,
        #[term(cv=MS, accession=1002189, name="TOPP IDConflictResolver", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP IDConflictResolver - Resolves ambiguous annotations of features with peptide identifications."]
        TOPPIDConflictResolver,
        #[term(cv=MS, accession=1002190, name="TOPP IDFilter", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP IDFilter - Filters results from protein or peptide identification engines based on different criteria."]
        TOPPIDFilter,
        #[term(cv=MS, accession=1002191, name="TOPP IDMapper", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP IDMapper - Assigns protein/peptide identifications to feature or consensus features."]
        TOPPIDMapper,
        #[term(cv=MS, accession=1002192, name="TOPP IDPosteriorErrorProbability", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP IDPosteriorErrorProbability - Estimates posterior error probabilities using a mixture model."]
        TOPPIDPosteriorErrorProbability,
        #[term(cv=MS, accession=1002193, name="TOPP IDRTCalibration", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP IDRTCalibration - Calibrate Retention times of peptide hits to standards."]
        TOPPIDRTCalibration,
        #[term(cv=MS, accession=1002194, name="TOPP PeptideIndexer", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP PeptideIndexer - Refreshes the protein references for all peptide hits."]
        TOPPPeptideIndexer,
        #[term(cv=MS, accession=1002195, name="TOPP PrecursorIonSelector", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP PrecursorIonSelector - A tool for precursor ion selection based on identification results."]
        TOPPPrecursorIonSelector,
        #[term(cv=MS, accession=1002196, name="TOPP MRMMapper", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP MRMMapper - MRMMapper maps measured chromatograms (mzML) and the transitions used (TraML)."]
        TOPPMRMMapper,
        #[term(cv=MS, accession=1002197, name="TOPP OpenSwath component", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP OpenSwath component - OpenSwath component of the TOPP software."]
        TOPPOpenSwathComponent,
        #[term(cv=MS, accession=1002198, name="TOPP OpenSwathAnalyzer", flags={0}, parents={["MS:1002197"]})]
        #[doc="TOPP OpenSwathAnalyzer - Picks peaks and finds features in an SRM experiment."]
        TOPPOpenSwathAnalyzer,
        #[term(cv=MS, accession=1002199, name="TOPP OpenSwathChromatogramExtractor", flags={0}, parents={["MS:1002197"]})]
        #[doc="TOPP OpenSwathChromatogramExtractor - Extract chromatograms (XIC) from a MS2 map file."]
        TOPPOpenSwathChromatogramExtractor,
        #[term(cv=MS, accession=1002200, name="TOPP OpenSwathDecoyGenerator", flags={0}, parents={["MS:1002197"]})]
        #[doc="TOPP OpenSwathDecoyGenerator - Generates decoys according to different models for a specific TraML."]
        TOPPOpenSwathDecoyGenerator,
        #[term(cv=MS, accession=1002201, name="TOPP OpenSwathFeatureXMLToTSV", flags={0}, parents={["MS:1002197"]})]
        #[doc="TOPP OpenSwathFeatureXMLToTSV - Converts a featureXML to a mProphet tsv (tab separated values)."]
        TOPPOpenSwathFeatureXMLToTSV,
        #[term(cv=MS, accession=1002202, name="TOPP OpenSwathRTNormalizer", flags={0}, parents={["MS:1002197"]})]
        #[doc="TOPP OpenSwathRTNormalizer - Generates a transformation file for retention time space into normalized space."]
        TOPPOpenSwathRTNormalizer,
        #[term(cv=MS, accession=1002203, name="TOPP ProteinInference", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP ProteinInference - Infer proteins from a list of (high-confidence) peptides."]
        TOPPProteinInference,
        #[term(cv=MS, accession=1002204, name="TOPP FalseDiscoveryRate", flags={0}, parents={["MS:1000752"]})]
        #[doc="TOPP FalseDiscoveryRate - Estimates the false discovery rate on peptide and protein level using decoy searches."]
        TOPPFalseDiscoveryRate,
        #[term(cv=MS, accession=1002205, name="ProteoWizard msconvert", flags={0}, parents={["MS:1000615"]})]
        #[doc="ProteoWizard msconvert - Converts, filters, and processes mass spectrometry data in variety of formats."]
        ProteoWizardMsconvert,
        #[term(cv=MS, accession=1002206, name="ProteoWizard idconvert", flags={0}, parents={["MS:1000615"]})]
        #[doc="ProteoWizard idconvert - Converts, filters, and processes identifications from shotgun proteomics experiments."]
        ProteoWizardIdconvert,
        #[term(cv=MS, accession=1002207, name="ProteoWizard chainsaw", flags={0}, parents={["MS:1000615"]})]
        #[doc="ProteoWizard chainsaw - Filters and processes protein sequence databases."]
        ProteoWizardChainsaw,
        #[term(cv=MS, accession=1002208, name="ProteoWizard msaccess", flags={0}, parents={["MS:1000615"]})]
        #[doc="ProteoWizard msaccess - Filters, processes, and displays mass spectrometry data in a variety of ways."]
        ProteoWizardMsaccess,
        #[term(cv=MS, accession=1002209, name="ProteoWizard SeeMS", flags={0}, parents={["MS:1000615"]})]
        #[doc="ProteoWizard SeeMS - An interactive GUI application to view and filter mass spectrometry data in a variety of formats."]
        ProteoWizardSeeMS,
        #[term(cv=MS, accession=1002210, name="IsobariQ", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        #[doc="IsobariQ - A quantitative software package designed for analysis of IPTL, TMT and iTRAQ data."]
        IsobariQ,
        #[term(cv=MS, accession=1002220, name="MRMaid", flags={0}, parents={["MS:1000871"]})]
        #[doc="MRMaid - A web-based SRM assay design tool whose transitions are generated by mining the millions of identified peptide spectra held in the EBI's PRIDE database."]
        MRMaid,
        #[term(cv=MS, accession=1002237, name="mzidLib", flags={1}, parents={["MS:1001456"]})]
        #[doc="mzidLib - A library of Java routines for manipulating mzIdentML files."]
        MzidLib,
        #[term(cv=MS, accession=1002238, name="mzidLib:Omssa2Mzid", flags={0}, parents={["MS:1002237"]})]
        #[doc="mzidLib:Omssa2Mzid - A converter for OMSSA OMX to mzIdentML."]
        MzidLibOmssa2Mzid,
        #[term(cv=MS, accession=1002239, name="mzidLib:Tandem2Mzid", flags={0}, parents={["MS:1002237"]})]
        #[doc="mzidLib:Tandem2Mzid - A converter for Tandem XML to mzIdentML."]
        MzidLibTandem2Mzid,
        #[term(cv=MS, accession=1002240, name="mzidLib:Csv2Mzid", flags={0}, parents={["MS:1002237"]})]
        #[doc="mzidLib:Csv2Mzid - A converter for CSV files (following OMSSA CSV style) to mzIdentML."]
        MzidLibCsv2Mzid,
        #[term(cv=MS, accession=1002241, name="mzidLib:ProteoGrouper", flags={0}, parents={["MS:1002237"]})]
        #[doc="mzidLib:ProteoGrouper - A generic and parameterizable protein inference algorithm for mzIdentML files."]
        MzidLibProteoGrouper,
        #[term(cv=MS, accession=1002242, name="mzidLib:Thresholder", flags={0}, parents={["MS:1002237"]})]
        #[doc="mzidLib:Thresholder - A routine for keeping only identifications passing a given threshold or setting passThreshold to true or false for SpectrumIdentificationItem or ProteinDetectionHypothesis in mzIdentML files."]
        MzidLibThresholder,
        #[term(cv=MS, accession=1002243, name="mzidLib:Perform emPAI on mzid", flags={0}, parents={["MS:1002237"]})]
        #[doc="mzidLib:Perform emPAI on mzid - A routine for adding emPAI quantitative values to an mzIdentML file."]
        MzidLibPerformEmPAIOnMzid,
        #[term(cv=MS, accession=1002244, name="mzidLib:FalseDiscoveryRate", flags={0}, parents={["MS:1002237"]})]
        #[doc="mzidLib:FalseDiscoveryRate - A routine for calculating local FDR, q-value and FDRScore for mzIdentML files, based on a decoy search."]
        MzidLibFalseDiscoveryRate,
        #[term(cv=MS, accession=1002245, name="mzidLib:Mzidentml2Csv", flags={0}, parents={["MS:1002237"]})]
        #[doc="mzidLib:Mzidentml2Csv - A tool for converting mzIdentML files to CSV format."]
        MzidLibMzidentml2Csv,
        #[term(cv=MS, accession=1002246, name="mzidLib:CombineSearchEngines", flags={0}, parents={["MS:1002237"]})]
        #[doc="mzidLib:CombineSearchEngines - A tool for combining results analysed in parallel in two or three search engines into a single mzIdentML file."]
        MzidLibCombineSearchEngines,
        #[term(cv=MS, accession=1002247, name="mzidLib:InsertMetaDataFromFasta", flags={0}, parents={["MS:1002237"]})]
        #[doc="mzidLib:InsertMetaDataFromFasta - A tool for adding additional meta data from a FASTA file to DBSequence entries (sequence and description) in mzIdentML files."]
        MzidLibInsertMetaDataFromFasta,
        #[term(cv=MS, accession=1002251, name="Comet", flags={1}, parents={["MS:1001456"]})]
        #[doc="Comet - Comet open-source sequence search engine developed at the University of Washington."]
        Comet,
        #[term(cv=MS, accession=1002261, name="Byonic", flags={1}, parents={["MS:1001456"]})]
        #[doc="Byonic - Byonic search engine from Protein Metrics."]
        Byonic,
        #[term(cv=MS, accession=1002285, name="Trans-Proteomic Pipeline", flags={1}, parents={["MS:1001456"]})]
        #[doc="Trans-Proteomic Pipeline - A suite of open source tools for the processing of MS2 proteomics data developed by the Seattle Proteome Center at the Institute for Systems Biology."]
        TransProteomicPipeline,
        #[term(cv=MS, accession=1002286, name="Trans-Proteomic Pipeline software", flags={1}, parents={["MS:1001456"]})]
        #[doc="Trans-Proteomic Pipeline software - A software program that is a component of the Trans-Proteomic Pipeline."]
        TransProteomicPipelineSoftware,
        #[term(cv=MS, accession=1002287, name="PeptideProphet", flags={0}, parents={["MS:1002286"]})]
        #[doc="PeptideProphet - A program in the TPP that calculates PSM probabilities for MS2 proteomics data searched with any of the supported sequence or spectral library search engines via the pepXML format."]
        PeptideProphet,
        #[term(cv=MS, accession=1002288, name="iProphet", flags={0}, parents={["MS:1002286"]})]
        #[doc="iProphet - A program in the TPP that calculates distinct peptide probabilities based on several lines of corroborating evidence including search results from multiple search engines via the pepXML format."]
        IProphet,
        #[term(cv=MS, accession=1002289, name="ProteinProphet", flags={0}, parents={["MS:1002286"]})]
        #[doc="ProteinProphet - A program in the TPP that calculates protein-level probabilities based on input PSM or peptide-level probabilities from PeptideProphet or iProphet. The output is written in the protXML format."]
        ProteinProphet,
        #[term(cv=MS, accession=1002290, name="XPRESS", flags={0}, parents={["MS:1002286"]})]
        #[doc="XPRESS - A program in the TPP that calculates PSM-level abundances based on 2-channel isotope-labelled data such as ICAT, SILAC, etc."]
        XPRESS,
        #[term(cv=MS, accession=1002291, name="Libra", flags={0}, parents={["MS:1002286"]})]
        #[doc="Libra - A program in the TPP that calculates PSM, peptide, and protein-level abundances based on N-channel isobaric label peptide data such as iTRAQ, TMT, etc."]
        Libra,
        #[term(cv=MS, accession=1002292, name="PTMProphet", flags={0}, parents={["MS:1002286"]})]
        #[doc="PTMProphet - A program in the TPP that calculates PTM localization probabilities by re-analyzing the peaks that are available to distinguish between possible modification sites."]
        PTMProphet,
        #[term(cv=MS, accession=1002333, name="conversion software", flags={2}, parents={["MS:1001457"]})]
        #[doc="conversion software - Computer software primarily designed to convert data represented in one format to another format, sometimes with minor data alterations in the process."]
        ConversionSoftware,
        #[term(cv=MS, accession=1002334, name="ProCon", flags={0}, parents={["MS:1002333"]})]
        #[doc="ProCon - Java software designed to convert one of several proteomics identification results formats into mzIdentML or PRIDE XML."]
        ProCon,
        #[term(cv=MS, accession=1002335, name="PRIDE Converter2", flags={0}, parents={["MS:1002333"]})]
        #[doc="PRIDE Converter2 - Java software designed to convert one of several proteomics identification results formats into PRIDE XML."]
        PRIDEConverter2,
        #[term(cv=MS, accession=1002336, name="Amanda", flags={1}, parents={["MS:1001456"]})]
        #[doc="Amanda - Amanda scoring system for PSM identification."]
        Amanda,
        #[term(cv=MS, accession=1002337, name="Andromeda", flags={1}, parents={["MS:1001456"]})]
        #[doc="Andromeda - Andromeda is a peptide search engine."]
        Andromeda,
        #[term(cv=MS, accession=1002342, name="mzmine", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="mzmine - A framework for differential analysis of mass spectrometry data."]
        Mzmine,
        #[term(cv=MS, accession=1002344, name="Maltcms", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="Maltcms - Modular Application Toolkit for Chromatography Mass-Spectrometry is an application framework mainly for developers."]
        Maltcms,
        #[term(cv=MS, accession=1002381, name="MALDI Solutions LC-MALDI", flags={7}, parents={["MS:1001455", "MS:1001456", "MS:1001457", "MS:1001557"]})]
        #[doc="MALDI Solutions LC-MALDI - Software for automated LC-MALDI analysis and reporting."]
        MALDISolutionsLCMALDI,
        #[term(cv=MS, accession=1002383, name="SCiLS software", flags={0}, parents={["MS:1000531"]})]
        #[doc="SCiLS software - SCiLS software for data acquisition and analysis."]
        SCiLSSoftware,
        #[term(cv=MS, accession=1002384, name="SCiLS Lab", flags={3}, parents={["MS:1002383", "MS:1001456", "MS:1001457"]})]
        #[doc="SCiLS Lab - SCiLS Lab software."]
        SCiLSLab,
        #[term(cv=MS, accession=1002386, name="preprocessing software", flags={2}, parents={["MS:1001457"]})]
        #[doc="preprocessing software - Preprocessing software."]
        PreprocessingSoftware,
        #[term(cv=MS, accession=1002387, name="PIA", flags={1}, parents={["MS:1002414", "MS:1001456"]})]
        #[doc="PIA - PIA - Protein Inference Algorithms, a toolbox for protein inference and identification analysis."]
        PIA,
        #[term(cv=MS, accession=1002410, name="Anubis", flags={0}, parents={["MS:1000871", "MS:1001139"]})]
        #[doc="Anubis - Anubis software for selected reaction monitoring data."]
        Anubis,
        #[term(cv=MS, accession=1002414, name="postprocessing software", flags={2}, parents={["MS:1001457"]})]
        #[doc="postprocessing software - Postprocessing software."]
        PostprocessingSoftware,
        #[term(cv=MS, accession=1002452, name="Maui", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="Maui - The Maltcms Graphical User Interface."]
        Maui,
        #[term(cv=MS, accession=1002458, name="PeptideShaker", flags={1}, parents={["MS:1001456"]})]
        #[doc="PeptideShaker - PeptideShaker is a software for the interpretation of proteomics identification results."]
        PeptideShaker,
        #[term(cv=MS, accession=1002524, name="PepFinder", flags={2}, parents={["MS:1001457"]})]
        #[doc="PepFinder - Thermo Scientific PepFinder BioPharma analysis software."]
        PepFinder,
        #[term(cv=MS, accession=1002543, name="xiFDR", flags={1}, parents={["MS:1001456"]})]
        #[doc="xiFDR - Target/Decoy based FDR estimation for crosslinking peptide-identifications."]
        XiFDR,
        #[term(cv=MS, accession=1002544, name="xi", flags={1}, parents={["MS:1001456"]})]
        #[doc="xi - Search engine for crosslinked peptides."]
        Xi,
        #[term(cv=MS, accession=1002546, name="Skyline mzQuantML converter", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        #[doc="Skyline mzQuantML converter - A software package to convert Skyline report to mzQuantML."]
        SkylineMzQuantMLConverter,
        #[term(cv=MS, accession=1002574, name="ASAPRatio", flags={0}, parents={["MS:1002286"]})]
        #[doc="ASAPRatio - A program in the TPP that calculates PSM, peptide, and protein-level abundances based on 2-channel isotope-labelled data such as ICAT, SILAC, etc."]
        ASAPRatio,
        #[term(cv=MS, accession=1002575, name="Tide", flags={1}, parents={["MS:1001456"]})]
        #[doc="Tide - Tide open-source sequence search program developed at the University of Washington."]
        Tide,
        #[term(cv=MS, accession=1002596, name="ProLuCID", flags={1}, parents={["MS:1001456"]})]
        #[doc="ProLuCID - The SEQUEST-like sequence search engine ProLuCID, developed in the Yates Lab at the Scripps Research Institute."]
        ProLuCID,
        #[term(cv=MS, accession=1002598, name="DTASelect", flags={1}, parents={["MS:1001456"]})]
        #[doc="DTASelect - Analysis software designed to reassemble the SEQUEST peptide identifications and to highlight the most significant matches."]
        DTASelect,
        #[term(cv=MS, accession=1002645, name="MSDK", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="MSDK - Mass Spectrometry Development Kit (MSDK) is a Java library of algorithms for processing of mass spectrometry data."]
        MSDK,
        #[term(cv=MS, accession=1002661, name="Morpheus", flags={1}, parents={["MS:1001456"]})]
        #[doc="Morpheus - Morpheus search engine."]
        Morpheus,
        #[term(cv=MS, accession=1002673, name="OpenXQuest", flags={0}, parents={["MS:1000752"]})]
        #[doc="OpenXQuest - Cross-Linking MS search engine."]
        OpenXQuest,
        #[term(cv=MS, accession=1002714, name="FLASHDeconv", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="FLASHDeconv - Ultrafast, High-Quality Feature Deconvolution for Top-Down Proteomics."]
        FLASHDeconv,
        #[term(cv=MS, accession=1002717, name="Waters DATA Convert", flags={3}, parents={["MS:1000694", "MS:1001456", "MS:1001457"]})]
        #[doc="Waters DATA Convert - Software for processing and exporting Waters MassLynx and waters_connect data."]
        WatersDATAConvert,
        #[term(cv=MS, accession=1002720, name="MSPathFinder", flags={1}, parents={["MS:1001456"]})]
        #[doc="MSPathFinder - PNNL top-down/bottom-up analysis software for identifying peptides and proteoforms in fragmentation mass spectra."]
        MSPathFinder,
        #[term(cv=MS, accession=1002750, name="NIST MSPepSearch", flags={1}, parents={["MS:1001456"]})]
        #[doc="NIST MSPepSearch - Search tool of the NIST (National Institute of Standards and Technology) for spectral library searches."]
        NISTMSPepSearch,
        #[term(cv=MS, accession=1002826, name="MetaMorpheus", flags={1}, parents={["MS:1001456"]})]
        #[doc="MetaMorpheus - MetaMorpheus search engine."]
        MetaMorpheus,
        #[term(cv=MS, accession=1002869, name="mzR", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="mzR - Bioconductor package mzR for reading and writing mass spectrometry data files."]
        MzR,
        #[term(cv=MS, accession=1002870, name="MSnbase", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="MSnbase - Bioconductor package MSnbase provides infrastructure for manipulation, processing and visualization of mass spectrometry and proteomics data, ranging from raw to quantitative and annotated data."]
        MSnbase,
        #[term(cv=MS, accession=1002871, name="CAMERA", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="CAMERA - Bioconductor package CAMERA for annotation of peak lists generated by xcms, rule based annotation of isotopes and adducts, isotope validation, EIC correlation based tagging of unknown adducts and fragments."]
        CAMERA,
        #[term(cv=MS, accession=1002878, name="small molecule analysis software", flags={1}, parents={["MS:1001456"]})]
        #[doc="small molecule analysis software - Software for the analysis of small molecules."]
        SmallMoleculeAnalysisSoftware,
        #[term(cv=MS, accession=1002879, name="Progenesis QI", flags={0}, parents={["MS:1002878"]})]
        #[doc="Progenesis QI - Metabolomics analysis software for LC-MS data from Nonlinear Dynamics."]
        ProgenesisQI,
        #[term(cv=MS, accession=1002880, name="Compound Discoverer", flags={0}, parents={["MS:1002878"]})]
        #[doc="Compound Discoverer - Metabolomics analysis software from Thermo Fisher Scientific."]
        CompoundDiscoverer,
        #[term(cv=MS, accession=1002881, name="MyCompoundID", flags={0}, parents={["MS:1002878"]})]
        #[doc="MyCompoundID - Metabolite identification tool MyCompoundID."]
        MyCompoundID,
        #[term(cv=MS, accession=1002901, name="TopPIC", flags={1}, parents={["MS:1001456"]})]
        #[doc="TopPIC - TopPIC: a software tool for top-down mass spectrometry-based proteoform identification and characterization."]
        TopPIC,
        #[term(cv=MS, accession=1002902, name="TopFD", flags={1}, parents={["MS:1001456"]})]
        #[doc="TopFD - Top-down mass spectral feature detection."]
        TopFD,
        #[term(cv=MS, accession=1002903, name="TopMG", flags={1}, parents={["MS:1001456"]})]
        #[doc="TopMG - A mass graph-based approach for the identification of modified proteoforms using top-down tandem mass spectra."]
        TopMG,
        #[term(cv=MS, accession=1002964, name="lipidomics analysis software", flags={0}, parents={["MS:1002878"]})]
        #[doc="lipidomics analysis software - Lipidomics analysis software."]
        LipidomicsAnalysisSoftware,
        #[term(cv=MS, accession=1002965, name="Lipid Data Analyzer", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="Lipid Data Analyzer - Lipid Data Analyzer software for lipid quantification."]
        LipidDataAnalyzer,
        #[term(cv=MS, accession=1002967, name="LipidHunter", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="LipidHunter - Software for identification of phospholipids by high-throughput processing of LC-MS and shotgun lipidomics datasets."]
        LipidHunter,
        #[term(cv=MS, accession=1002968, name="LipidXplorer", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="LipidXplorer - Software for consensual cross-platform lipidomics."]
        LipidXplorer,
        #[term(cv=MS, accession=1002969, name="LipidMatch", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="LipidMatch - An automated workflow for rule-based lipid identification using untargeted high-resolution tandem mass spectrometry data."]
        LipidMatch,
        #[term(cv=MS, accession=1002970, name="Greazy", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="Greazy - Open-source software for automated phospholipid tandem mass spectrometry identification."]
        Greazy,
        #[term(cv=MS, accession=1002971, name="LipidBlast", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="LipidBlast - LC-MS-based lipidomics and automated identification of lipids using the LipidBlast in-silico MS/MS library."]
        LipidBlast,
        #[term(cv=MS, accession=1002972, name="Lipid-Pro", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="Lipid-Pro - A computational lipid identification solution for untargeted lipidomics on data-independent acquisition tandem mass spectrometry platforms."]
        LipidPro,
        #[term(cv=MS, accession=1002973, name="LipidFinder", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="LipidFinder - A computational workflow for the discovery of lipids for the identification of eicosanoid-phosphoinositides in platelets."]
        LipidFinder,
        #[term(cv=MS, accession=1002974, name="LipiDex", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="LipiDex - An integrated software package for high-confidence lipid identification."]
        LipiDex,
        #[term(cv=MS, accession=1002975, name="LIQUID", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="LIQUID - An-open source software for identifying lipids in LC-MS/MS-based lipidomics data."]
        LIQUID,
        #[term(cv=MS, accession=1002976, name="ALEX", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="ALEX - Analysis of lipid experiments, a calculator for m/z values of intact lipid molecules (MS1)."]
        ALEX,
        #[term(cv=MS, accession=1002977, name="ALEX123", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="ALEX123 - Analysis of lipid experiments 123, a calculator with m/z values of intact lipid molecules (MS1) and their fragment ions at the MS2 and MS3 level."]
        ALEX123,
        #[term(cv=MS, accession=1002978, name="LIMSA", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="LIMSA - Software tool for the quantitative analysis of mass spectrometric lipidome data."]
        LIMSA,
        #[term(cv=MS, accession=1002979, name="LOBSTAHS", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="LOBSTAHS - Adduct-Based lipidomics software for the discovery and identification of oxidative stress biomarkers."]
        LOBSTAHS,
        #[term(cv=MS, accession=1002980, name="LipidQA", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="LipidQA - Lipid qualitative/quantitative analysis software for identification and quantitation of complex lipid molecular species."]
        LipidQA,
        #[term(cv=MS, accession=1002981, name="Proline", flags={1}, parents={["MS:1001456"]})]
        #[doc="Proline - The Proline software suite for mass spectrometry based proteomics."]
        Proline,
        #[term(cv=MS, accession=1002982, name="PepNovo", flags={1}, parents={["MS:1001456"]})]
        #[doc="PepNovo - PepNovo tool for de novo peptide sequencing."]
        PepNovo,
        #[term(cv=MS, accession=1002983, name="pNovo", flags={1}, parents={["MS:1001456"]})]
        #[doc="pNovo - pNovo tool for de novo peptide sequencing and identification using HCD spectra."]
        PNovo,
        #[term(cv=MS, accession=1002984, name="Novor", flags={1}, parents={["MS:1001456"]})]
        #[doc="Novor - Novor real-time peptide de novo sequencing software tool."]
        Novor,
        #[term(cv=MS, accession=1002987, name="IdentiPy", flags={1}, parents={["MS:1001456"]})]
        #[doc="IdentiPy - IdentiPy."]
        IdentiPy,
        #[term(cv=MS, accession=1002990, name="ms_deisotope", flags={2}, parents={["MS:1001457"]})]
        #[doc="ms_deisotope - ms_deisotope, a library for deisotoping and charge state deconvolution of mass spectra."]
        MsDeisotope,
        #[term(cv=MS, accession=1002991, name="python-psims", flags={0}, parents={["MS:1002333"]})]
        #[doc="python-psims - python-psims, a library for generating mzML and mzIdentML."]
        PythonPsims,
        #[term(cv=MS, accession=1003010, name="LPPtiger", flags={2}, parents={["MS:1002964", "MS:1001457"]})]
        #[doc="LPPtiger - Software for lipidome-specific prediction and identification of oxidized phospholipids from LC-MS datasets."]
        LPPtiger,
        #[term(cv=MS, accession=1003011, name="pFind", flags={1}, parents={["MS:1001456"]})]
        #[doc="pFind - Sequence-tag-based search engine pFind."]
        PFind,
        #[term(cv=MS, accession=1003013, name="i3tms", flags={1}, parents={["MS:1001456"]})]
        #[doc="i3tms - i3-tms search engine and data-analysis software."]
        I3tms,
        #[term(cv=MS, accession=1003014, name="MSFragger", flags={1}, parents={["MS:1001456"]})]
        #[doc="MSFragger - A database search-based peptide identification tool."]
        MSFragger,
        #[term(cv=MS, accession=1003018, name="Philosopher", flags={1}, parents={["MS:1001456"]})]
        #[doc="Philosopher - General proteomics processing toolkit for shotgun proteomics."]
        Philosopher,
        #[term(cv=MS, accession=1003023, name="OpenPepXL", flags={0}, parents={["MS:1000752"]})]
        #[doc="OpenPepXL - Cross-Linking MS search engine."]
        OpenPepXL,
        #[term(cv=MS, accession=1003082, name="MS-DIAL", flags={2}, parents={["MS:1002878", "MS:1001457"]})]
        #[doc="MS-DIAL - Data processing software for untargeted metabolomics and lipidomics that supports multiple instruments and MS vendors."]
        MSDIAL,
        #[term(cv=MS, accession=1003108, name="PatternLab", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        #[doc="PatternLab - PatternLab for Proteomics is an integrated computational environment for analyzing shotgun proteomic data."]
        PatternLab,
        #[term(cv=MS, accession=1003109, name="SIM-XL", flags={1}, parents={["MS:1001456"]})]
        #[doc="SIM-XL - Identifying crosslinked peptides in complex protein mixtures"]
        SIMXL,
        #[term(cv=MS, accession=1003111, name="QUIN-XL", flags={0}, parents={["MS:1001139"]})]
        #[doc="QUIN-XL - Quantification of crosslinked peptides in complex protein mixtures"]
        QUINXL,
        #[term(cv=MS, accession=1003118, name="EPIFANY", flags={1}, parents={["MS:1001456", "MS:1000752"]})]
        #[doc="EPIFANY - A Method for Efficient High-Confidence Protein Inference. The tool is part of the OpenMS framework"]
        EPIFANY,
        #[term(cv=MS, accession=1003141, name="ProSight", flags={1}, parents={["MS:1001456"]})]
        #[doc="ProSight - ProSight: Database search engine for top-down proteomics."]
        ProSight,
        #[term(cv=MS, accession=1003142, name="TDPortal", flags={1}, parents={["MS:1001456"]})]
        #[doc="TDPortal - TDPortal: Database search engine for top-down proteomics."]
        TDPortal,
        #[term(cv=MS, accession=1003145, name="ThermoRawFileParser", flags={2}, parents={["MS:1001457"]})]
        #[doc="ThermoRawFileParser - Cross-platform software to convert Thermo RAW files to a number of open formats."]
        ThermoRawFileParser,
        #[term(cv=MS, accession=1003146, name="pyteomics", flags={1}, parents={["MS:1001456"]})]
        #[doc="pyteomics - Python module that helps handling various proteomics data analysis tasks."]
        Pyteomics,
        #[term(cv=MS, accession=1003162, name="PTX-QC", flags={1}, parents={["MS:1001456"]})]
        #[doc="PTX-QC - Proteomics (PTX) - QualityControl (QC) software for QC report generation and visualization."]
        PTXQC,
        #[term(cv=MS, accession=1003164, name="QuaMeter IDFree", flags={1}, parents={["MS:1001456"]})]
        #[doc="QuaMeter IDFree - QuaMeter IDFree software for QC metric calculation."]
        QuaMeterIDFree,
        #[term(cv=MS, accession=1003165, name="iMonDB", flags={1}, parents={["MS:1001456"]})]
        #[doc="iMonDB - iMonDB software to extract, store, and manage mass spectrometry instrument parameters from raw data files."]
        IMonDB,
        #[term(cv=MS, accession=1003202, name="BiblioSpec", flags={1}, parents={["MS:1001456", "MS:1003207"]})]
        #[doc="BiblioSpec - A suite of software tools for creating and searching MS/MS peptide spectrum libraries, developed at the University of Washington"]
        BiblioSpec,
        #[term(cv=MS, accession=1003207, name="library creation software", flags={0}, parents={["MS:1000531", "MS:1003171"]})]
        #[doc="library creation software - Library creation software"]
        LibraryCreationSoftware,
        #[term(cv=MS, accession=1003232, name="PeakForest", flags={1}, parents={["MS:1001456", "MS:1003207", "MS:1002878"]})]
        #[doc="PeakForest - comprehensive infrastructure to organize, curate and share a multi- instrument spectral library for metabolomics data annotation developed and distributed by the French National infrastructure in metabolomics and fluxomics (MetaboHUB)."]
        PeakForest,
        #[term(cv=MS, accession=1003253, name="DIA-NN", flags={1}, parents={["MS:1001139", "MS:1001456", "MS:1003207"]})]
        #[doc="DIA-NN - A universal software for data-independent acquisition (DIA) proteomics data processing"]
        DIANN,
        #[term(cv=MS, accession=1003281, name="Casanovo", flags={1}, parents={["MS:1001456"]})]
        #[doc="Casanovo - Casanovo is a deep learning-based de novo spectrum identification tool. Official website https://github.com/Noble-Lab/casanovo/."]
        Casanovo,
        #[term(cv=MS, accession=1003309, name="Goslin", flags={2}, parents={["MS:1001457", "MS:1002414", "MS:1002964"]})]
        #[doc="Goslin - The Goslin implementations parse common lipid name dialects and normalize them to the recent lipid shorthand nomenclature based on grammars on succinct lipid nomenclature."]
        Goslin,
        #[term(cv=MS, accession=1003357, name="ANN-SoLo", flags={1}, parents={["MS:1001456"]})]
        #[doc="ANN-SoLo - ANN-SoLo (Approximate Nearest Neighbor Spectral Library) is a spectral library search engine for fast and accurate open modification searching. ANN-SoLo uses approximate nearest neighbor indexing to speed up open modification searching by selecting only a limited number of the most relevant library spectra to compare to an unknown query spectrum. This is combined with a cascade search strategy to maximize the number of identified unmodified and modified spectra while strictly controlling the false discovery rate and the shifted dot product score to sensitively match modified spectra to their unmodified counterpart."]
        ANNSoLo,
        #[term(cv=MS, accession=1003376, name="ChemClipse", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="ChemClipse - ChemClipse is part of the Eclipse Science project. Primarily developed by Lablicate GmbH."]
        ChemClipse,
        #[term(cv=MS, accession=1003377, name="OpenChrom", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="OpenChrom - OpenChrom is an Open Source software for data processing and analysis. Based upon Eclipse ChemClipse."]
        OpenChrom,
        #[term(cv=MS, accession=1003382, name="waters_connect", flags={7}, parents={["MS:1000694", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        #[doc="waters_connect - Waters Corporation waters_connect software for liquid chromatography and mass spectrometry acquisition and processing."]
        WatersConnect,
        #[term(cv=MS, accession=1003386, name="Spectra", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="Spectra - Bioconductor package Spectra for mass spectrometry data representation and processing."]
        Spectra,
        #[term(cv=MS, accession=1003387, name="MetaboAnnotation", flags={3}, parents={["MS:1001456", "MS:1001457"]})]
        #[doc="MetaboAnnotation - Bioconductor package MetaboAnnotation for annotation of untargeted metabolomics data."]
        MetaboAnnotation,
        #[term(cv=MS, accession=1003388, name="CompoundDb", flags={3}, parents={["MS:1001456", "MS:1001457", "MS:1003207"]})]
        #[doc="CompoundDb - Bioconductor package CompoundDb for creation, usage and maintenance of public or library-specific annotation databases and spectra libraries."]
        CompoundDb,
        #[term(cv=MS, accession=1003399, name="quality control software", flags={2}, parents={["MS:1001457"]})]
        #[doc="quality control software - Software that creates or manipulates QC-related data."]
        QualityControlSoftware,
        #[term(cv=MS, accession=1003400, name="rmzqc", flags={0}, parents={["MS:1003399"]})]
        #[doc="rmzqc - An R package for reading, validating, and writing mzQC files."]
        Rmzqc,
        #[term(cv=MS, accession=1003401, name="jmzqc", flags={0}, parents={["MS:1003399"]})]
        #[doc="jmzqc - A Java package for reading, validating, and writing mzQC files."]
        Jmzqc,
        #[term(cv=MS, accession=1003402, name="pymzqc", flags={0}, parents={["MS:1003399"]})]
        #[doc="pymzqc - A Python package for reading, validating, and writing mzQC files."]
        Pymzqc,
        #[term(cv=MS, accession=1003405, name="mzRecal", flags={2}, parents={["MS:1001457"]})]
        #[doc="mzRecal - MS1 recalibration using identified peptides as internal calibrants."]
        MzRecal,
        #[term(cv=MS, accession=1003406, name="spectrum clustering software", flags={0}, parents={["MS:1000531"]})]
        #[doc="spectrum clustering software - Software designed to group multiple mass spectra by high similarity, generally with the goal of grouping replicate spectra derived from the same analyte."]
        SpectrumClusteringSoftware,
        #[term(cv=MS, accession=1003407, name="Scout", flags={1}, parents={["MS:1001456"]})]
        #[doc="Scout - Identifying crosslinked peptides in complex protein mixtures"]
        Scout,
        #[term(cv=MS, accession=1003413, name="Kojak", flags={1}, parents={["MS:1001456"]})]
        #[doc="Kojak - Kojak open-source crosslinked peptide sequence search engine developed at the Institute for Systems Biology."]
        Kojak,
        #[term(cv=MS, accession=1003425, name="quantms", flags={1}, parents={["MS:1001456", "MS:1001139"]})]
        #[doc="quantms - Cloud-based pipeline for quantitative proteomics that enables the reanalysis of public proteomics data."]
        Quantms,
        #[term(cv=MS, accession=1003426, name="xQuest/xProphet", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        #[doc="xQuest/xProphet - A software to identify cross-linked peptides from LC-MS/MS spectra."]
        XQuestXProphet,
        #[term(cv=MS, accession=1003427, name="PeakView", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        #[doc="PeakView - A software for spectral analysis and data interrogation in qualitative review of LC-MS and MS/MS data."]
        PeakView,
        #[term(cv=MS, accession=1003428, name="Perseus", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        #[doc="Perseus - A software for interpreting protein quantification, interaction and post-translational modification data."]
        Perseus,
        #[term(cv=MS, accession=1003429, name="FragPipe", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        #[doc="FragPipe - A computational platform for analyzing mass spectrometry-based proteomics data."]
        FragPipe,
        #[term(cv=MS, accession=1003430, name="OpenMS", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        #[doc="OpenMS - A software for LC-MS data management and analysis."]
        OpenMS,
        #[term(cv=MS, accession=1003431, name="pLink", flags={1}, parents={["MS:1001139", "MS:1001456"]})]
        #[doc="pLink - A tool for the analysis of chemically cross-linked proteins using mass spectrometry."]
        PLink,
        #[term(cv=MS, accession=1003432, name="pLink2", flags={1}, parents={["MS:1003431", "MS:1001139", "MS:1001456"]})]
        #[doc="pLink2 - Upgraded version of pLink tool, provides a graphical user interface and faster with newly designed index structure."]
        PLink2,
        #[term(cv=MS, accession=1003446, name="SCIEX OS", flags={7}, parents={["MS:1000690", "MS:1001455", "MS:1001456", "MS:1001457"]})]
        #[doc="SCIEX OS - SCIEX OS software."]
        SCIEXOS,
        #[term(cv=MS, accession=1003447, name="SCIEX MS Data Converter", flags={2}, parents={["MS:1000690", "MS:1001457"]})]
        #[doc="SCIEX MS Data Converter - A software for converting SCIEX wiff or wiff2 format to mzML."]
        SCIEXMSDataConverter,
        #[term(cv=MS, accession=4000151, name="MsQuality", flags={1}, parents={["MS:1001456"]})]
        #[doc="MsQuality - MsQuality – an interoperable open-source package for the calculation of standardized quality metrics of mass spectrometry data."]
        MsQuality,
        #[term(cv=MS, accession=4000189, name="DIAMetric", flags={1}, parents={["MS:1001456"]})]
        #[doc="DIAMetric - DIAMetric is a Data-Independent Acquisition Quality Metric Generator."]
        DIAMetric,
    }
    //[[[end]]] (sum: hwFwvphdPv)
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
