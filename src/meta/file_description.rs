use regex::{self, Regex};
use std::io;
use std::path::Path;

use crate::impl_param_described;
use crate::io::infer_format;
use crate::params::{
    ControlledVocabulary, Param, ParamDescribed, ParamList, ParamValue, ValueRef, CURIE,
};

/// Description of a source file, including location and type.
///
/// This is usually another mass spectrometry data file. For some
/// vendor raw formats that are directories, there can be many files.
/// See <https://peptideatlas.org/tmp/mzML1.1.0.html#sourceFile>.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SourceFile {
    /// The name of the source file without any path or location information
    pub name: String,
    /// The URI-formatted location where the file was retrieved.
    pub location: String,
    /// A unique identifier for this source file
    pub id: String,
    /// The [`MassSpectrometerFileFormatTerm`]-defined parameter for this file
    pub file_format: Option<Param>,
    /// The [`NativeSpectrumIdentifierFormatTerm`]-defined parameter for this file
    pub id_format: Option<Param>,
    /// The rest of the parameters for this file.
    pub params: ParamList,
}

impl SourceFile {
    /// Create a new [`SourceFile`] from a path.
    ///
    /// This function makes a minimal effort to infer information about the file,
    /// using [`infer_format`] to populate [`SourceFile::file_format`]
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let path = path.as_ref();
        let format = infer_format(path)
            .ok()
            .and_then(|(format, _)| format.as_param());
        let inst = Self {
            name: path
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_default(),
            location: path
                .canonicalize()?
                .parent()
                .map(|s| format!("file://{}", s.to_string_lossy()))
                .unwrap_or_else(|| "file://".to_string()),
            file_format: format,
            ..Default::default()
        };
        Ok(inst)
    }

    /// Convert [`SourceFile::id_format`] into a [`NativeSpectrumIDFormat`] carrying its own
    /// parser machinery, if such a term mapping exists
    pub fn native_id_format(&self) -> Option<NativeSpectrumIDFormat> {
        self.id_format
            .as_ref()
            .and_then(|p| p.curie())
            .and_then(|p| NativeSpectrumIdentifierFormatTerm::from_curie(&p))
            .map(|t| t.build())
    }
}

/// A description of the file data file and its contents
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FileDescription {
    /// Descriptors of the content spectra
    pub contents: ParamList,
    /// Any source files involved in producing the current file, such as vendor raw files.
    pub source_files: Vec<SourceFile>,
}

impl FileDescription {
    pub fn new(contents: ParamList, source_files: Vec<SourceFile>) -> Self {
        Self {
            contents,
            source_files,
        }
    }

    /// Checks to see if the "MS1 spectrum" term is present in the file contents
    ///
    /// **Note**: This does not actually inspect the spectra in the file, only the metadata,
    /// which may be incorrect/missing.
    pub fn has_ms1_spectra(&self) -> bool {
        self.get_param_by_curie(&CURIE::new(ControlledVocabulary::MS, 1000579))
            .is_some()
    }

    /// Checks to see if the "MSn spectrum" term is present in the file contents.
    ///
    /// **Note**: This does not actually inspect the spectra in the file, only the metadata,
    /// which may be incorrect/missing.
    pub fn has_msn_spectra(&self) -> bool {
        self.get_param_by_curie(&CURIE::new(ControlledVocabulary::MS, 1000580))
            .is_some()
    }

    pub fn has_contents(&self) -> bool {
        !self.contents.is_empty()
    }
}

impl_param_described!(SourceFile);

impl ParamDescribed for FileDescription {
    fn params(&self) -> &[Param] {
        &self.contents
    }

    fn params_mut(&mut self) -> &mut ParamList {
        &mut self.contents
    }
}

crate::cvmap! {
    #[flag_type=&str]
    #[allow(unused)]
    #[doc = "A text-based schema that defines how native spectrum identifiers are formatted.
    These patterns are often found in mzML-compatible formats."]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_native_ids.py']).decode('utf8')
    for line in buf.splitlines():
        cog.outl(line)
    ]]]*/
    pub enum NativeSpectrumIdentifierFormatTerm {
        #[term(cv=MS, accession=1000767, name="native spectrum identifier format", flags={r"(.+)"}, parents={[]})]
        #[doc = r"native spectrum identifier format - `(.+)`"]
        NativeSpectrumIdentifierFormat,
        #[term(cv=MS, accession=1000768, name="Thermo nativeID format", flags={r"controllerType=(?<controllerType>\d+) controllerNumber=(?<controllerNumber>\d+) scan=(?<scan>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"Thermo nativeID format - `controllerType=(?<controllerType>\d+) controllerNumber=(?<controllerNumber>\d+) scan=(?<scan>\d+)`"]
        ThermoNativeIDFormat,
        #[term(cv=MS, accession=1000769, name="Waters nativeID format", flags={r"function=(?<function>\d+) process=(?<process>\d+) scan=(?<scan>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"Waters nativeID format - `function=(?<function>\d+) process=(?<process>\d+) scan=(?<scan>\d+)`"]
        WatersNativeIDFormat,
        #[term(cv=MS, accession=1000770, name="WIFF nativeID format", flags={r"sample=(?<sample>\d+) period=(?<period>\d+) cycle=(?<cycle>\d+) experiment=(?<experiment>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"WIFF nativeID format - `sample=(?<sample>\d+) period=(?<period>\d+) cycle=(?<cycle>\d+) experiment=(?<experiment>\d+)`"]
        WIFFNativeIDFormat,
        #[term(cv=MS, accession=1000771, name="Bruker/Agilent YEP nativeID format", flags={r"scan=(?<scan>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"Bruker/Agilent YEP nativeID format - `scan=(?<scan>\d+)`"]
        BrukerAgilentYEPNativeIDFormat,
        #[term(cv=MS, accession=1000772, name="Bruker BAF nativeID format", flags={r"scan=(?<scan>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"Bruker BAF nativeID format - `scan=(?<scan>\d+)`"]
        BrukerBAFNativeIDFormat,
        #[term(cv=MS, accession=1000773, name="Bruker FID nativeID format", flags={r"file=(?<file>\S+)"}, parents={["MS:1000767"]})]
        #[doc = r"Bruker FID nativeID format - `file=(?<file>\S+)`"]
        BrukerFIDNativeIDFormat,
        #[term(cv=MS, accession=1000774, name="multiple peak list nativeID format", flags={r"index=(?<index>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"multiple peak list nativeID format - `index=(?<index>\d+)`"]
        MultiplePeakListNativeIDFormat,
        #[term(cv=MS, accession=1000775, name="single peak list nativeID format", flags={r"file=(?<file>\S+)"}, parents={["MS:1000767"]})]
        #[doc = r"single peak list nativeID format - `file=(?<file>\S+)`"]
        SinglePeakListNativeIDFormat,
        #[term(cv=MS, accession=1000776, name="scan number only nativeID format", flags={r"scan=(?<scan>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"scan number only nativeID format - `scan=(?<scan>\d+)`"]
        ScanNumberOnlyNativeIDFormat,
        #[term(cv=MS, accession=1000777, name="spectrum identifier nativeID format", flags={r"spectrum=(?<spectrum>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"spectrum identifier nativeID format - `spectrum=(?<spectrum>\d+)`"]
        SpectrumIdentifierNativeIDFormat,
        #[term(cv=MS, accession=1000823, name="Bruker U2 nativeID format", flags={r"declaration=(?<declaration>\d+) collection=(?<collection>\d+) scan=(?<scan>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"Bruker U2 nativeID format - `declaration=(?<declaration>\d+) collection=(?<collection>\d+) scan=(?<scan>\d+)`"]
        BrukerU2NativeIDFormat,
        #[term(cv=MS, accession=1000824, name="no nativeID format", flags={r"(.+)"}, parents={["MS:1000767"]})]
        #[doc = r"no nativeID format - `(.+)`"]
        NoNativeIDFormat,
        #[term(cv=MS, accession=1000929, name="Shimadzu Biotech nativeID format", flags={r"source=(?<source>\S+) start=(?<start>\d+) end=(?<end>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"Shimadzu Biotech nativeID format - `source=(?<source>\S+) start=(?<start>\d+) end=(?<end>\d+)`"]
        ShimadzuBiotechNativeIDFormat,
        #[term(cv=MS, accession=1001186, name="Mobilion MBI nativeID format", flags={r"frame=(?<frame>\d+) scan=(?<scan>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"Mobilion MBI nativeID format - `frame=(?<frame>\d+) scan=(?<scan>\d+)`"]
        MobilionMBINativeIDFormat,
        #[term(cv=MS, accession=1001480, name="SCIEX TOF/TOF nativeID format", flags={r"jobRun=(?<jobRun>\d+) spotLabel=(?<spotLabel>\S+) spectrum=(?<spectrum>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"SCIEX TOF/TOF nativeID format - `jobRun=(?<jobRun>\d+) spotLabel=(?<spotLabel>\S+) spectrum=(?<spectrum>\d+)`"]
        SCIEXTOFTOFNativeIDFormat,
        #[term(cv=MS, accession=1001508, name="Agilent MassHunter nativeID format", flags={r"scanId=(?<scanId>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"Agilent MassHunter nativeID format - `scanId=(?<scanId>\d+)`"]
        AgilentMassHunterNativeIDFormat,
        #[term(cv=MS, accession=1001526, name="spectrum from database integer nativeID format", flags={r"databasekey=(?<databasekey>-?\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"spectrum from database integer nativeID format - `databasekey=(?<databasekey>-?\d+)`"]
        SpectrumFromDatabaseIntegerNativeIDFormat,
        #[term(cv=MS, accession=1001528, name="Mascot query number", flags={r"query=(?<query>\d+)"}, parents={["MS:1000767", "MS:1001405"]})]
        #[doc = r"Mascot query number - `query=(?<query>\d+)`"]
        MascotQueryNumber,
        #[term(cv=MS, accession=1001531, name="spectrum from ProteinScape database nativeID format", flags={r"databasekey=(?<databasekey>-?\d+)"}, parents={["MS:1000767", "MS:1001529"]})]
        #[doc = r"spectrum from ProteinScape database nativeID format - `databasekey=(?<databasekey>-?\d+)`"]
        SpectrumFromProteinScapeDatabaseNativeIDFormat,
        #[term(cv=MS, accession=1001532, name="spectrum from database string nativeID format", flags={r"databasekey=(?<databasekey>\S+)"}, parents={["MS:1000767", "MS:1001529"]})]
        #[doc = r"spectrum from database string nativeID format - `databasekey=(?<databasekey>\S+)`"]
        SpectrumFromDatabaseStringNativeIDFormat,
        #[term(cv=MS, accession=1001559, name="SCIEX TOF/TOF T2D nativeID format", flags={r"file=(?<file>\S+)"}, parents={["MS:1000767"]})]
        #[doc = r"SCIEX TOF/TOF T2D nativeID format - `file=(?<file>\S+)`"]
        SCIEXTOFTOFT2DNativeIDFormat,
        #[term(cv=MS, accession=1001562, name="Scaffold nativeID format", flags={r"(.+)"}, parents={["MS:1000767"]})]
        #[doc = r"Scaffold nativeID format - `(.+)`"]
        ScaffoldNativeIDFormat,
        #[term(cv=MS, accession=1002303, name="Bruker Container nativeID format", flags={r"(.+)"}, parents={["MS:1000767"]})]
        #[doc = r"Bruker Container nativeID format - `(.+)`"]
        BrukerContainerNativeIDFormat,
        #[term(cv=MS, accession=1002532, name="UIMF nativeID format", flags={r"frame=(?<frame>\d+) scan=(?<scan>\d+) frameType=(?<frameType>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"UIMF nativeID format - `frame=(?<frame>\d+) scan=(?<scan>\d+) frameType=(?<frameType>\d+)`"]
        UIMFNativeIDFormat,
        #[term(cv=MS, accession=1002818, name="Bruker TDF nativeID format", flags={r"frame=(?<frame>\d+) scan=(?<scan>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"Bruker TDF nativeID format - `frame=(?<frame>\d+) scan=(?<scan>\d+)`"]
        BrukerTDFNativeIDFormat,
        #[term(cv=MS, accession=1002898, name="Shimadzu Biotech QTOF nativeID format", flags={r"scan=(?<scan>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"Shimadzu Biotech QTOF nativeID format - `scan=(?<scan>\d+)`"]
        ShimadzuBiotechQTOFNativeIDFormat,
        #[term(cv=MS, accession=1003283, name="Bruker TSF nativeID format", flags={r"frame=(?<frame>\d+)"}, parents={["MS:1000767"]})]
        #[doc = r"Bruker TSF nativeID format - `frame=(?<frame>\d+)`"]
        BrukerTSFNativeIDFormat,
    }
    //[[[end]]] (sum: q+gQ/J50Sc)
}

/// A text-based schema that defines how native spectrum identifiers are formatted.
///
/// These patterns are often found in mzML-compatible formats.
#[derive(Debug, Clone)]
pub struct NativeSpectrumIDFormat {
    pub term: NativeSpectrumIdentifierFormatTerm,
    parser: Regex,
    field_names: Vec<Option<String>>,
}

impl PartialEq for NativeSpectrumIDFormat {
    fn eq(&self, other: &Self) -> bool {
        self.term == other.term
    }
}

impl From<NativeSpectrumIdentifierFormatTerm> for NativeSpectrumIDFormat {
    fn from(value: NativeSpectrumIdentifierFormatTerm) -> Self {
        Self::new(value)
    }
}

#[derive(Debug, Clone, thiserror::Error, PartialEq)]
pub enum NativeIDFormatError {
    /// The ID format required a different number of groups than what were found
    #[error("{term:?} required {expected} arguments, but received {received} arguments")]
    IncorrectArgumentNumber {
        term: NativeSpectrumIdentifierFormatTerm,
        expected: usize,
        received: usize,
    },
    /// The ID format's pattern did not match the provided string
    #[error("{term:?} did not match {text}")]
    PatternMismatch {
        term: NativeSpectrumIdentifierFormatTerm,
        text: String,
    },
}

impl NativeSpectrumIDFormat {
    pub fn new(term: NativeSpectrumIdentifierFormatTerm) -> Self {
        let parser = term.parser();
        let field_names = parser
            .capture_names()
            .skip(1)
            .map(|s| s.map(|i| i.to_string()))
            .collect();
        Self {
            term,
            parser,
            field_names,
        }
    }

    pub const fn name(&self) -> &str {
        self.term.name()
    }

    pub const fn curie(&self) -> CURIE {
        CURIE::new(self.term.controlled_vocabulary(), self.term.accession())
    }

    /// This parses the provided string, returning the captured groups of the ID pattern if they are present
    /// in a [`regex::Captures`] structure that can be indexed by group number.
    pub fn parse<'h>(&self, ident: &'h str) -> Option<regex::Captures<'h>> {
        self.parser.captures(ident)
    }

    /// This parses the provided string, returning the capture groups as (name, value) pairs they are present
    pub fn parse_named<'h>(
        &self,
        ident: &'h str,
    ) -> Result<Vec<(Option<String>, &'h str)>, NativeIDFormatError> {
        if let Some(hits) = self.parser.captures(ident) {
            Ok(self
                .parser
                .capture_names()
                .enumerate()
                .map(|(i, name)| {
                    let m = if let Some(name_) = name {
                        hits.name(name_).unwrap()
                    } else {
                        hits.get(i).unwrap()
                    };
                    (name.map(|s| s.to_string()), m.as_str())
                })
                .collect())
        } else {
            Err(NativeIDFormatError::PatternMismatch {
                term: self.term,
                text: ident.to_string(),
            })
        }
    }

    /// Given the field values of a nativeID format, create string in that format
    pub fn format<'h>(
        &self,
        values: impl IntoIterator<Item = ValueRef<'h>>,
    ) -> Result<String, NativeIDFormatError> {
        let mut buffer = String::with_capacity(64);
        let names = &self.field_names;
        let n_names = names.len().saturating_sub(1);
        let values: Vec<_> = values.into_iter().collect();
        if values.len() != names.len() {
            return Err(NativeIDFormatError::IncorrectArgumentNumber {
                term: self.term,
                expected: names.len(),
                received: values.len(),
            });
        }
        for (i, (k, v)) in names.iter().zip(values).enumerate() {
            match k {
                Some(k) => {
                    buffer.push_str(k);
                    buffer.push('=');
                    buffer.push_str(&v.as_str());
                }
                None => {
                    buffer.push_str(&v.as_str());
                }
            };
            if i < n_names {
                buffer.push(' ');
            }
        }
        Ok(buffer)
    }
}

impl NativeSpectrumIdentifierFormatTerm {
    /// Create a new [`regex::Regex`] for this identifier format.
    pub fn parser(&self) -> regex::Regex {
        regex::Regex::new(self.flags()).unwrap()
    }

    /// This parses the provided string, returning the captured groups of the ID pattern if they are present
    /// in a [`regex::Captures`] structure that can be indexed by group number.
    ///
    /// # Note
    /// This method creates a new regular expression on each invocation, making it expensive to invoke.
    /// If you must call this repeatedly, instead use [`NativeSpectrumIdentifierFormatTerm::build`] to create a
    /// the regular expression once and re-use it directly.
    pub fn parse<'h>(&self, ident: &'h str) -> Option<regex::Captures<'h>> {
        let parser = self.parser();
        parser.captures(ident)
    }

    /// Create a [`NativeSpectrumIDFormat`] that owns the [`Regex`] produced
    /// by [`NativeSpectrumIdentifierFormatTerm::parser`]
    pub fn build(&self) -> NativeSpectrumIDFormat {
        (*self).into()
    }

    /// Given the field values of a nativeID format, create string in that format
    ///
    /// # Note
    /// This method creates a new expression formatter every time.
    /// Use [`NativeSpectrumIdentifierFormatTerm::build`] to create a re-useable
    /// parser/formatter.
    pub fn format<'h>(&self, values: impl IntoIterator<Item = ValueRef<'h>>) -> String {
        self.build().format(values).unwrap()
    }

    /// This parses the provided string, returning the capture groups as (name, value) pairs they are present.
    ///
    /// # Note
    /// This method creates a new regular expression on each invocation, making it expensive to invoke.
    /// If you must call this repeatedly, instead use [`NativeSpectrumIdentifierFormatTerm::build`] to create a
    /// the regular expression once and re-use it directly.
    pub fn parse_named<'h>(&self, ident: &'h str) -> Vec<(Option<String>, &'h str)> {
        self.build().parse_named(ident).unwrap()
    }
}

crate::cvmap! {
    #[flag_type=i32]
    #[allow(unused)]
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    /*[[[cog
    import cog
    import subprocess
    buf = subprocess.check_output(['python', 'cv/extract_file_formats.py']).decode('utf8')
    for line in buf.splitlines():
        cog.outl(line)
    ]]]*/
    pub enum MassSpectrometerFileFormatTerm {
        #[term(cv=MS, accession=1000526, name="Waters raw format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Waters raw format - Waters data file format found in a Waters RAW directory, generated from an MS acquisition."]
        WatersRaw,
        #[term(cv=MS, accession=1000560, name="mass spectrometer file format", flags={0}, parents={["MS:1001459"]})]
        #[doc = "mass spectrometer file format - The format of the file being used. This could be a instrument or vendor specific proprietary file format or a converted open file format."]
        MassSpectrometerFile,
        #[term(cv=MS, accession=1000562, name="ABI WIFF format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "ABI WIFF format - Applied Biosystems WIFF file format."]
        ABIWIFF,
        #[term(cv=MS, accession=1000563, name="Thermo RAW format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Thermo RAW format - Thermo Scientific RAW file format."]
        ThermoRAW,
        #[term(cv=MS, accession=1000564, name="PSI mzData format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "PSI mzData format - Proteomics Standards Inititative mzData file format."]
        PSIMzData,
        #[term(cv=MS, accession=1000565, name="Micromass PKL format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Micromass PKL format - Micromass PKL file format."]
        MicromassPKL,
        #[term(cv=MS, accession=1000566, name="ISB mzXML format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "ISB mzXML format - Institute of Systems Biology mzXML file format."]
        ISBMzXML,
        #[term(cv=MS, accession=1000567, name="Bruker/Agilent YEP format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Bruker/Agilent YEP format - Bruker/Agilent YEP file format."]
        BrukerAgilentYEP,
        #[term(cv=MS, accession=1000584, name="mzML format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "mzML format - Proteomics Standards Inititative mzML file format."]
        MzML,
        #[term(cv=MS, accession=1000613, name="DTA format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "DTA format - SEQUEST DTA file format."]
        DTA,
        #[term(cv=MS, accession=1000614, name="ProteinLynx Global Server mass spectrum XML format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "ProteinLynx Global Server mass spectrum XML format - Peak list file format used by ProteinLynx Global Server."]
        ProteinLynxGlobalServerMassSpectrumXML,
        #[term(cv=MS, accession=1000740, name="parameter file", flags={0}, parents={["MS:1000560"]})]
        #[doc = "parameter file - Parameter file used to configure the acquisition of raw data on the instrument."]
        ParameterFile,
        #[term(cv=MS, accession=1000742, name="Bioworks SRF format", flags={0}, parents={["MS:1000560", "MS:1001040"]})]
        #[doc = "Bioworks SRF format - Thermo Finnigan SRF file format."]
        BioworksSRF,
        #[term(cv=MS, accession=1000815, name="Bruker BAF format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Bruker BAF format - Bruker BAF raw file format."]
        BrukerBAF,
        #[term(cv=MS, accession=1000816, name="Bruker U2 format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Bruker U2 format - Bruker HyStar U2 file format."]
        BrukerU2,
        #[term(cv=MS, accession=1000825, name="Bruker FID format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Bruker FID format - Bruker FID file format."]
        BrukerFID,
        #[term(cv=MS, accession=1000930, name="Shimadzu Biotech database entity", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Shimadzu Biotech database entity - Shimadzu Biotech format."]
        ShimadzuBiotechDatabaseEntity,
        #[term(cv=MS, accession=1001062, name="Mascot MGF format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Mascot MGF format - Mascot MGF file format."]
        MascotMGF,
        #[term(cv=MS, accession=1001185, name="Mobilion MBI format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Mobilion MBI format - Mobilion MBI file format."]
        MobilionMBI,
        #[term(cv=MS, accession=1001245, name="PerSeptive PKS format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "PerSeptive PKS format - PerSeptive peak list file format."]
        PerSeptivePKS,
        #[term(cv=MS, accession=1001246, name="SCIEX API III format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "SCIEX API III format - PE SCIEX peak list file format."]
        SCIEXAPIIII,
        #[term(cv=MS, accession=1001247, name="Bruker XML format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Bruker XML format - Bruker data exchange XML format."]
        BrukerXML,
        #[term(cv=MS, accession=1001369, name="text format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "text format - Simple text file format of 'm/z \\[intensity\\]' values for a PMF (or single MS2) search."]
        Text,
        #[term(cv=MS, accession=1001463, name="Phenyx XML format", flags={0}, parents={["MS:1000560", "MS:1001040"]})]
        #[doc = "Phenyx XML format - Phenyx open XML file format."]
        PhenyxXML,
        #[term(cv=MS, accession=1001466, name="MS2 format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "MS2 format - MS2 file format for MS2 spectral data."]
        MS2,
        #[term(cv=MS, accession=1001481, name="SCIEX TOF/TOF database", flags={0}, parents={["MS:1000560"]})]
        #[doc = "SCIEX TOF/TOF database - Applied Biosystems/MDS Analytical Technologies TOF/TOF instrument database."]
        SCIEXTOFTOFDatabase,
        #[term(cv=MS, accession=1001509, name="Agilent MassHunter format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Agilent MassHunter format - A data file format found in an Agilent MassHunter directory which contains raw data acquired by an Agilent mass spectrometer."]
        AgilentMassHunter,
        #[term(cv=MS, accession=1001527, name="Proteinscape spectra", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Proteinscape spectra - Spectra from Bruker/Protagen Proteinscape database."]
        ProteinscapeSpectra,
        #[term(cv=MS, accession=1001560, name="SCIEX TOF/TOF T2D format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "SCIEX TOF/TOF T2D format - Applied Biosystems/MDS Analytical Technologies TOF/TOF instrument export format."]
        SCIEXTOFTOFT2D,
        #[term(cv=MS, accession=1001881, name="mz5 format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "mz5 format - mz5 file format, modelled after mzML."]
        Mz5,
        #[term(cv=MS, accession=1002302, name="Bruker Container format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Bruker Container format - Bruker Container raw file format."]
        BrukerContainer,
        #[term(cv=MS, accession=1002385, name="SCiLS Lab format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "SCiLS Lab format - SCiLS Lab file format."]
        SCiLSLab,
        #[term(cv=MS, accession=1002441, name="Andi-MS format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Andi-MS format - AIA Analytical Data Interchange file format for mass spectrometry data."]
        AndiMS,
        #[term(cv=MS, accession=1002531, name="UIMF format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "UIMF format - SQLite-based file format created at Pacific Northwest National Lab. It stores an intermediate analysis of ion-mobility mass spectrometry data."]
        UIMF,
        #[term(cv=MS, accession=1002597, name="MS1 format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "MS1 format - MS1 file format for MS1 spectral data."]
        MS1,
        #[term(cv=MS, accession=1002817, name="Bruker TDF format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Bruker TDF format - Bruker TDF raw file format."]
        BrukerTDF,
        #[term(cv=MS, accession=1002838, name="mzMLb format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "mzMLb format - mzMLb file format, mzML encapsulated within HDF5."]
        MzMLb,
        #[term(cv=MS, accession=1002899, name="msalign format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "msalign format - msalign file format."]
        Msalign,
        #[term(cv=MS, accession=1002900, name="feature format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "feature format - TopFD feature file format."]
        Feature,
        #[term(cv=MS, accession=1002966, name="chrom format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "chrom format - The Lipid Data Analyzer native chrom format."]
        Chrom,
        #[term(cv=MS, accession=1002996, name="Andromeda:apl file format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Andromeda:apl file format - Peak list file format of the Andromeda search engine."]
        AndromedaAplFile,
        #[term(cv=MS, accession=1003009, name="Shimadzu Biotech LCD format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Shimadzu Biotech LCD format - Shimadzu Biotech LCD file format."]
        ShimadzuBiotechLCD,
        #[term(cv=MS, accession=1003282, name="Bruker TSF format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Bruker TSF format - Bruker TSF raw file format."]
        BrukerTSF,
        #[term(cv=MS, accession=1003374, name="Open Chromatography Binary OCB format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "Open Chromatography Binary OCB format - ChemClipse/OpenChrom file format."]
        OpenChromatographyBinaryOCB,
        #[term(cv=MS, accession=1003448, name="SCIEX WIFF2 format", flags={0}, parents={["MS:1000560"]})]
        #[doc = "SCIEX WIFF2 format - SCIEX WIFF2 file format."]
        SCIEXWIFF2,
    }
    //[[[end]]] (sum: tevATyl5zS)
}

#[allow(unused, clippy::upper_case_acronyms)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[doc = "The kinds of spectrum that might be represented"]
pub enum SpectrumType {
    #[doc = "mass spectrum - A plot of the relative abundance of a beam or other collection of ions as a function of the mass-to-charge ratio (m/z)."]
    MassSpectrum,
    #[doc = "charge inversion mass spectrum - The measurement of the relative abundance of ions that result from a charge inversion reaction as a function of m/z."]
    ChargeInversionMassSpectrum,
    #[doc = "constant neutral gain spectrum - A spectrum formed of all product ions that have been produced by gain of a pre-selected neutral mass following the reaction with and addition of the gas in a collision cell."]
    ConstantNeutralGainSpectrum,
    #[doc = "constant neutral loss spectrum - A spectrum formed of all product ions that have been produced with a selected m/z decrement from any precursor ions. The spectrum shown correlates to the precursor ion spectrum. See also neutral loss spectrum."]
    ConstantNeutralLossSpectrum,
    #[doc = "e/2 mass spectrum - A mass spectrum obtained using a sector mass spectrometer in which the electric sector field E is set to half the value required to transmit the main ion-beam. This spectrum records the signal from doubly charged product ions of charge-stripping reactions."]
    E2MassSpectrum,
    #[doc = "precursor ion spectrum - Spectrum generated by scanning precursor m/z while monitoring a fixed product m/z."]
    PrecursorIonSpectrum,
    #[doc = "product ion spectrum - OBSOLETE A mass spectrum recorded from any spectrometer in which the appropriate m/z separation scan function is set to record the product ion or ions of selected precursor ions."]
    ProductIonSpectrum,
    #[doc = "MS1 spectrum - Mass spectrum created by a single-stage MS experiment or the first stage of a multi-stage experiment."]
    MS1Spectrum,
    #[doc = "MSn spectrum - MSn refers to multi-stage MS2 experiments designed to record product ion spectra where n is the number of product ion stages (progeny ions). For ion traps, sequential MS/MS experiments can be undertaken where n > 2 whereas for a simple triple quadrupole system n=2. Use the term ms level (MS:1000511) for specifying n."]
    MSnSpectrum,
    #[doc = "CRM spectrum - Spectrum generated from MSn experiment with three or more stages of m/z separation and in which a particular multi-step reaction path is monitored."]
    CRMSpectrum,
    #[doc = "SIM spectrum - Spectrum obtained with the operation of a mass spectrometer in which the abundances of one ion or several ions of specific m/z values are recorded rather than the entire mass spectrum (Selected Ion Monitoring)."]
    SIMSpectrum,
    #[doc = "SRM spectrum - Spectrum obtained when data are acquired from specific product ions corresponding to m/z values of selected precursor ions a recorded via two or more stages of mass spectrometry. The precursor/product ion pair is called a transition pair. Data can be obtained for a single transition pair or multiple transition pairs. Multiple time segments of different transition pairs can exist in a single file. Single precursor ions can have multiple product ions consitituting multiple transition pairs. Selected reaction monitoring can be performed as tandem mass spectrometry in time or tandem mass spectrometry in space."]
    SRMSpectrum,
    #[doc = "PDA spectrum - OBSOLETE Spectrum generated from a photodiode array detector (ultraviolet/visible spectrum)."]
    PDASpectrum,
    #[doc = "enhanced multiply charged spectrum - MS1 spectrum that is enriched in multiply-charged ions compared to singly-charged ions."]
    EnhancedMultiplyChargedSpectrum,
    #[doc = "time-delayed fragmentation spectrum - MSn spectrum in which the product ions are collected after a time delay, which allows the observation of lower energy fragmentation processes after precursor ion activation."]
    TimeDelayedFragmentationSpectrum,
    #[doc = "electromagnetic radiation spectrum - A plot of the relative intensity of electromagnetic radiation as a function of the wavelength."]
    ElectromagneticRadiationSpectrum,
    #[doc = "emission spectrum - A plot of the relative intensity of electromagnetic radiation emitted by atoms or molecules when excited."]
    EmissionSpectrum,
    #[doc = "absorption spectrum - A plot of the relative intensity of electromagnetic radiation absorbed by atoms or molecules when excited."]
    AbsorptionSpectrum,
}
#[doc = r" These methods are part of the controlled vocabulary mapping"]
impl SpectrumType {
    #[doc = r" Retrieve the accession number for this term, independent of its controlled vocabulary"]
    pub const fn accession(&self) -> crate::params::AccessionIntCode {
        match self {
            Self::MassSpectrum => 1000294,
            Self::ChargeInversionMassSpectrum => 1000322,
            Self::ConstantNeutralGainSpectrum => 1000325,
            Self::ConstantNeutralLossSpectrum => 1000326,
            Self::E2MassSpectrum => 1000328,
            Self::PrecursorIonSpectrum => 1000341,
            Self::ProductIonSpectrum => 1000343,
            Self::MS1Spectrum => 1000579,
            Self::MSnSpectrum => 1000580,
            Self::CRMSpectrum => 1000581,
            Self::SIMSpectrum => 1000582,
            Self::SRMSpectrum => 1000583,
            Self::PDASpectrum => 1000620,
            Self::EnhancedMultiplyChargedSpectrum => 1000789,
            Self::TimeDelayedFragmentationSpectrum => 1000790,
            Self::ElectromagneticRadiationSpectrum => 1000804,
            Self::EmissionSpectrum => 1000805,
            Self::AbsorptionSpectrum => 1000806,
        }
    }
    #[doc = r" Retrieve the controlled vocabulary this term belongs to"]
    pub const fn controlled_vocabulary(&self) -> crate::params::ControlledVocabulary {
        match self {
            Self::MassSpectrum => crate::params::ControlledVocabulary::MS,
            Self::ChargeInversionMassSpectrum => crate::params::ControlledVocabulary::MS,
            Self::ConstantNeutralGainSpectrum => crate::params::ControlledVocabulary::MS,
            Self::ConstantNeutralLossSpectrum => crate::params::ControlledVocabulary::MS,
            Self::E2MassSpectrum => crate::params::ControlledVocabulary::MS,
            Self::PrecursorIonSpectrum => crate::params::ControlledVocabulary::MS,
            Self::ProductIonSpectrum => crate::params::ControlledVocabulary::MS,
            Self::MS1Spectrum => crate::params::ControlledVocabulary::MS,
            Self::MSnSpectrum => crate::params::ControlledVocabulary::MS,
            Self::CRMSpectrum => crate::params::ControlledVocabulary::MS,
            Self::SIMSpectrum => crate::params::ControlledVocabulary::MS,
            Self::SRMSpectrum => crate::params::ControlledVocabulary::MS,
            Self::PDASpectrum => crate::params::ControlledVocabulary::MS,
            Self::EnhancedMultiplyChargedSpectrum => crate::params::ControlledVocabulary::MS,
            Self::TimeDelayedFragmentationSpectrum => crate::params::ControlledVocabulary::MS,
            Self::ElectromagneticRadiationSpectrum => crate::params::ControlledVocabulary::MS,
            Self::EmissionSpectrum => crate::params::ControlledVocabulary::MS,
            Self::AbsorptionSpectrum => crate::params::ControlledVocabulary::MS,
        }
    }
    #[doc = r" Retrieve the plain text human readable name for this term"]
    pub const fn name(&self) -> &'static str {
        match self {
            Self::MassSpectrum => "mass spectrum",
            Self::ChargeInversionMassSpectrum => "charge inversion mass spectrum",
            Self::ConstantNeutralGainSpectrum => "constant neutral gain spectrum",
            Self::ConstantNeutralLossSpectrum => "constant neutral loss spectrum",
            Self::E2MassSpectrum => "e/2 mass spectrum",
            Self::PrecursorIonSpectrum => "precursor ion spectrum",
            Self::ProductIonSpectrum => "product ion spectrum",
            Self::MS1Spectrum => "MS1 spectrum",
            Self::MSnSpectrum => "MSn spectrum",
            Self::CRMSpectrum => "CRM spectrum",
            Self::SIMSpectrum => "SIM spectrum",
            Self::SRMSpectrum => "SRM spectrum",
            Self::PDASpectrum => "PDA spectrum",
            Self::EnhancedMultiplyChargedSpectrum => "enhanced multiply charged spectrum",
            Self::TimeDelayedFragmentationSpectrum => "time-delayed fragmentation spectrum",
            Self::ElectromagneticRadiationSpectrum => "electromagnetic radiation spectrum",
            Self::EmissionSpectrum => "emission spectrum",
            Self::AbsorptionSpectrum => "absorption spectrum",
        }
    }
    #[doc = r" Attempt to map a string by name to retrieve one of the terms from this"]
    #[doc = r" set."]
    #[doc = r""]
    #[doc = r" If no match is found, [`None`] is returned."]
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "mass spectrum" => Some(Self::MassSpectrum),
            "charge inversion mass spectrum" => Some(Self::ChargeInversionMassSpectrum),
            "constant neutral gain spectrum" => Some(Self::ConstantNeutralGainSpectrum),
            "constant neutral loss spectrum" => Some(Self::ConstantNeutralLossSpectrum),
            "e/2 mass spectrum" => Some(Self::E2MassSpectrum),
            "precursor ion spectrum" => Some(Self::PrecursorIonSpectrum),
            "product ion spectrum" => Some(Self::ProductIonSpectrum),
            "MS1 spectrum" => Some(Self::MS1Spectrum),
            "MSn spectrum" => Some(Self::MSnSpectrum),
            "CRM spectrum" => Some(Self::CRMSpectrum),
            "SIM spectrum" => Some(Self::SIMSpectrum),
            "SRM spectrum" => Some(Self::SRMSpectrum),
            "PDA spectrum" => Some(Self::PDASpectrum),
            "enhanced multiply charged spectrum" => Some(Self::EnhancedMultiplyChargedSpectrum),
            "time-delayed fragmentation spectrum" => Some(Self::TimeDelayedFragmentationSpectrum),
            "electromagnetic radiation spectrum" => Some(Self::ElectromagneticRadiationSpectrum),
            "emission spectrum" => Some(Self::EmissionSpectrum),
            "absorption spectrum" => Some(Self::AbsorptionSpectrum),
            _ => None,
        }
    }
    #[doc = r" Attempt to map the numeric accession number to retrieve one of the terms from this"]
    #[doc = r" set."]
    #[doc = r""]
    #[doc = r" If no match is found, [`None`] is returned."]
    pub const fn from_accession(accession: crate::params::AccessionIntCode) -> Option<Self> {
        match accession {
            1000294 => Some(Self::MassSpectrum),
            1000322 => Some(Self::ChargeInversionMassSpectrum),
            1000325 => Some(Self::ConstantNeutralGainSpectrum),
            1000326 => Some(Self::ConstantNeutralLossSpectrum),
            1000328 => Some(Self::E2MassSpectrum),
            1000341 => Some(Self::PrecursorIonSpectrum),
            1000343 => Some(Self::ProductIonSpectrum),
            1000579 => Some(Self::MS1Spectrum),
            1000580 => Some(Self::MSnSpectrum),
            1000581 => Some(Self::CRMSpectrum),
            1000582 => Some(Self::SIMSpectrum),
            1000583 => Some(Self::SRMSpectrum),
            1000620 => Some(Self::PDASpectrum),
            1000789 => Some(Self::EnhancedMultiplyChargedSpectrum),
            1000790 => Some(Self::TimeDelayedFragmentationSpectrum),
            1000804 => Some(Self::ElectromagneticRadiationSpectrum),
            1000805 => Some(Self::EmissionSpectrum),
            1000806 => Some(Self::AbsorptionSpectrum),
            _ => None,
        }
    }
    #[doc = r" Convert this term into a [`ParamCow`](crate::params::ParamCow) without a value."]
    pub const fn to_param(self) -> crate::params::ParamCow<'static> {
        crate::params::ParamCow::const_new(
            self.name(),
            crate::params::ValueRef::Empty,
            Some(self.accession()),
            Some(self.controlled_vocabulary()),
            crate::params::Unit::Unknown,
        )
    }
    #[doc = r" Convert a [`CURIE`]($crate::params::CURIE) by accession."]
    #[doc = r""]
    #[doc = r" If no match is found, [`None`] is returned."]
    pub const fn from_curie(curie: &crate::params::CURIE) -> Option<Self> {
        if matches!(
            curie.controlled_vocabulary,
            crate::params::ControlledVocabulary::MS
        ) {
            Self::from_accession(curie.accession)
        } else {
            None
        }
    }
    #[doc = r" Attempt to convert a [`ParamCow`](crate::params::ParamCow) to a term from this set."]
    #[doc = r""]
    #[doc = r" If no match is found, [`None`] is returned."]
    #[doc = r""]
    #[doc = r" # Note"]
    #[doc = r" This method can be called in `const` contexts, requiring the type be [`ParamCow`](crate::params::ParamCow) with a `'static`"]
    #[doc = r" lifetime parameter, but the regular [`From`] trait is implemented for all [`ParamLike`](crate::params::ParamLike) types."]
    pub const fn from_param(p: &crate::params::ParamCow<'static>) -> Option<Self> {
        if let Some(acc) = p.accession {
            Self::from_accession(acc)
        } else {
            None
        }
    }
    #[doc = r" Retrieve a term set specific set of flags"]
    pub fn flags(&self) -> i32 {
        match self {
            Self::MassSpectrum => 0,
            Self::ChargeInversionMassSpectrum => 0,
            Self::ConstantNeutralGainSpectrum => 0,
            Self::ConstantNeutralLossSpectrum => 0,
            Self::E2MassSpectrum => 0,
            Self::PrecursorIonSpectrum => 0,
            Self::ProductIonSpectrum => 0,
            Self::MS1Spectrum => 0,
            Self::MSnSpectrum => 0,
            Self::CRMSpectrum => 0,
            Self::SIMSpectrum => 0,
            Self::SRMSpectrum => 0,
            Self::PDASpectrum => 0,
            Self::EnhancedMultiplyChargedSpectrum => 0,
            Self::TimeDelayedFragmentationSpectrum => 0,
            Self::ElectromagneticRadiationSpectrum => 0,
            Self::EmissionSpectrum => 0,
            Self::AbsorptionSpectrum => 0,
        }
    }
    #[doc = r" Retrieve the list of zero or more terms in the set which are"]
    #[doc = r" parents of this term."]
    pub fn parents(&self) -> Vec<Self> {
        match self {
            Self::MassSpectrum => { ["MS:1000524", "MS:1000559"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::ChargeInversionMassSpectrum => { ["MS:1000294"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::ConstantNeutralGainSpectrum => { ["MS:1000294"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::ConstantNeutralLossSpectrum => { ["MS:1000294"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::E2MassSpectrum => { ["MS:1000294"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::PrecursorIonSpectrum => { ["MS:1000294"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::ProductIonSpectrum => { ["MS:1000294"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::MS1Spectrum => { ["MS:1000294"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::MSnSpectrum => { ["MS:1000294"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::CRMSpectrum => { ["MS:1000294"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::SIMSpectrum => { ["MS:1000294"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::SRMSpectrum => { ["MS:1000294"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::PDASpectrum => { ["MS:1000524", "MS:1000559"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::EnhancedMultiplyChargedSpectrum => { ["MS:1000579"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::TimeDelayedFragmentationSpectrum => { ["MS:1000580"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::ElectromagneticRadiationSpectrum => { ["MS:1000524", "MS:1000559"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::EmissionSpectrum => { ["MS:1000524", "MS:1000559"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
            Self::AbsorptionSpectrum => { ["MS:1000524", "MS:1000559"] }
                .iter()
                .flat_map(|s: &&str| {
                    let curie = s.parse::<crate::params::CURIE>().unwrap();
                    Self::from_accession(curie.accession)
                })
                .collect(),
        }
    }
}
impl<P> From<P> for SpectrumType
where
    P: crate::params::ParamLike,
{
    fn from(value: P) -> Self {
        Self::from_accession(value.accession().expect(concat!(
            "Cannot convert an uncontrolled parameter to ",
            stringify!(SpectrumType)
        )))
        .unwrap_or_else(|| {
            panic!(
                "Could not map {:?}:{} to {}",
                value.controlled_vocabulary().unwrap(),
                value.accession().unwrap(),
                stringify!(SpectrumType)
            )
        })
    }
}
impl From<SpectrumType> for crate::params::ParamCow<'static> {
    fn from(value: SpectrumType) -> Self {
        value.to_param()
    }
}
impl From<SpectrumType> for crate::params::Param {
    fn from(value: SpectrumType) -> Self {
        value.to_param().into()
    }
}
impl From<&SpectrumType> for crate::params::ParamCow<'static> {
    fn from(value: &SpectrumType) -> Self {
        value.to_param()
    }
}
impl From<&SpectrumType> for crate::params::Param {
    fn from(value: &SpectrumType) -> Self {
        value.to_param().into()
    }
}

macro_rules! t {
    ($t:expr) => {
        ($t, $t.to_param())
    };
}
const SPECTRUM_TYPES: &[(crate::meta::SpectrumType, crate::params::ParamCow<'static>)] = &[
    t!(crate::meta::SpectrumType::MS1Spectrum),
    t!(crate::meta::SpectrumType::MSnSpectrum),
    t!(crate::meta::SpectrumType::MassSpectrum),
    t!(crate::meta::SpectrumType::ChargeInversionMassSpectrum),
    t!(crate::meta::SpectrumType::ConstantNeutralGainSpectrum),
    t!(crate::meta::SpectrumType::ConstantNeutralLossSpectrum),
    t!(crate::meta::SpectrumType::E2MassSpectrum),
    t!(crate::meta::SpectrumType::PrecursorIonSpectrum),
    t!(crate::meta::SpectrumType::ProductIonSpectrum),
    t!(crate::meta::SpectrumType::MS1Spectrum),
    t!(crate::meta::SpectrumType::MSnSpectrum),
    t!(crate::meta::SpectrumType::CRMSpectrum),
    t!(crate::meta::SpectrumType::SIMSpectrum),
    t!(crate::meta::SpectrumType::SRMSpectrum),
    t!(crate::meta::SpectrumType::PDASpectrum),
    t!(crate::meta::SpectrumType::EnhancedMultiplyChargedSpectrum),
    t!(crate::meta::SpectrumType::TimeDelayedFragmentationSpectrum),
    t!(crate::meta::SpectrumType::ElectromagneticRadiationSpectrum),
    t!(crate::meta::SpectrumType::EmissionSpectrum),
    t!(crate::meta::SpectrumType::AbsorptionSpectrum),
];

impl SpectrumType {
    /// Check if this a mass spectrum or some other kind of spectrum
    pub fn is_mass_spectrum(&self) -> bool {
        self.parents().contains(&Self::MassSpectrum) || *self == Self::MassSpectrum
    }

    /// Get the default measurement dimension for this kind of spectrum
    pub const fn default_main_axis(&self) -> crate::spectrum::ArrayType {
        match self {
            SpectrumType::PDASpectrum => crate::spectrum::ArrayType::WavelengthArray,
            SpectrumType::ElectromagneticRadiationSpectrum => {
                crate::spectrum::ArrayType::WavelengthArray
            }
            SpectrumType::EmissionSpectrum => crate::spectrum::ArrayType::WavelengthArray,
            SpectrumType::AbsorptionSpectrum => crate::spectrum::ArrayType::WavelengthArray,
            _ => crate::spectrum::ArrayType::MZArray,
        }
    }

    /// Get an array of all available types
    pub const fn all_types() -> &'static [(SpectrumType, crate::params::ParamCow<'static>)] {
        SPECTRUM_TYPES
    }
}

#[cfg(test)]
mod test {
    use crate::params::ParamCow;

    use super::*;

    #[test]
    fn test_source_file() -> io::Result<()> {
        let sf = SourceFile::from_path("test/data/small.mzML")?;
        assert_eq!(
            MassSpectrometerFileFormatTerm::MzML.to_param(),
            sf.file_format.clone().unwrap()
        );
        Ok(())
    }

    #[test]
    fn test_parser() {
        let ident = NativeSpectrumIdentifierFormatTerm::ThermoNativeIDFormat
            .parse("controllerType=0 controllerNumber=1 scan=25788")
            .unwrap();
        let scan_number = ident.name("scan").unwrap().as_str();
        assert_eq!(scan_number, "25788");
    }

    #[test]
    fn test_meta() {
        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat.flags(),
            r"frame=(?<frame>\d+) scan=(?<scan>\d+)"
        );
        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat.parents(),
            [NativeSpectrumIdentifierFormatTerm::NativeSpectrumIdentifierFormat]
        );
        let param: ParamCow<'static> =
            (NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat).into();
        assert_eq!(
            param.curie().unwrap(),
            CURIE::new(ControlledVocabulary::MS, 1002818)
        );
        let param: ParamCow<'static> =
            (&NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat).into();
        assert_eq!(
            param.curie().unwrap(),
            CURIE::new(ControlledVocabulary::MS, 1002818)
        );

        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::from_param(&param),
            Some(NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat)
        );

        let param: Param = (NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat).into();
        assert_eq!(
            param.curie().unwrap(),
            CURIE::new(ControlledVocabulary::MS, 1002818)
        );
        let param: Param = (&NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat).into();
        assert_eq!(
            param.curie().unwrap(),
            CURIE::new(ControlledVocabulary::MS, 1002818)
        );

        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat.accession(),
            1002818
        );
        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat.controlled_vocabulary(),
            ControlledVocabulary::MS
        );
        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat.name(),
            "Bruker TDF nativeID format"
        );

        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::from_name("Bruker TDF nativeID format"),
            Some(NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat)
        );

        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::from_curie(&CURIE::new(
                ControlledVocabulary::MS,
                1002818
            )),
            Some(NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat)
        )
    }

    #[test]
    fn test_format() {
        let fmt = NativeSpectrumIdentifierFormatTerm::ThermoNativeIDFormat.format([
            ValueRef::Int(0),
            ValueRef::Int(1),
            ValueRef::Int(25788),
        ]);
        assert_eq!(fmt, "controllerType=0 controllerNumber=1 scan=25788");

        let fmt = NativeSpectrumIdentifierFormatTerm::ThermoNativeIDFormat
            .build()
            .format([ValueRef::Int(0), ValueRef::Int(1)])
            .unwrap_err();
        assert_eq!(
            fmt,
            NativeIDFormatError::IncorrectArgumentNumber {
                term: NativeSpectrumIdentifierFormatTerm::ThermoNativeIDFormat,
                expected: 3,
                received: 2
            }
        );
    }
}
