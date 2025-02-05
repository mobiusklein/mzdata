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
    //[[[end]]] (checksum: abe810fc9e7449c83803cd458f5cdc6c)
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
    }
    //[[[end]]] (checksum: 7e00698426115a6a91ea38badd77a05d)
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
            NativeSpectrumIdentifierFormatTerm::from_curie(&CURIE::new(ControlledVocabulary::MS, 1002818)),
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
