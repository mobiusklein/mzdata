use super::*;

/// Controlled vocabularies used in mass spectrometry data files
#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(u8)]
pub enum ControlledVocabulary {
    /// The PSI-MS Controlled Vocabulary [https://www.ebi.ac.uk/ols4/ontologies/ms](https://www.ebi.ac.uk/ols4/ontologies/ms)
    MS = 1,
    /// The Unit Ontology [https://www.ebi.ac.uk/ols4/ontologies/uo](https://www.ebi.ac.uk/ols4/ontologies/uo)
    UO,
    /// The Experimental Factor Ontology <https://www.ebi.ac.uk/ols4/ontologies/efo>
    EFO,
    /// The Ontology for Biomedical Investigations <https://www.ebi.ac.uk/ols4/ontologies/obi>
    OBI,
    /// The Human Ancestry Ontology <https://www.ebi.ac.uk/ols4/ontologies/hancestro>
    HANCESTRO,
    /// The Basic Formal Ontology <https://www.ebi.ac.uk/ols4/ontologies/bfo>
    BFO,
    /// The NCI Thesaurus OBO Edition <https://www.ebi.ac.uk/ols4/ontologies/ncit>
    NCIT,
    /// The BRENDA Tissue Ontology <https://www.ebi.ac.uk/ols4/ontologies/bto>
    BTO,
    /// The PRIDE Controlled Vocabulary <https://www.ebi.ac.uk/ols4/ontologies/pride>
    PRIDE,
    /// The Imaging MS Controlled Vocabulary <https://www.ms-imaging.org/imzml/>
    #[cfg(feature = "imzml")]
    IMS,
    Unknown,
}

const MS_CV: &str = "MS";
const UO_CV: &str = "UO";
const EFO_CV: &str = "EFO";
const OBI_CV: &str = "OBI";
const HANCESTRO_CV: &str = "HANCESTRO";
const BFO_CV: &str = "BFO";
const BTO_CV: &str = "BTO";
const NCIT_CV: &str = "NCIT";
const PRIDE_CV: &str = "PRIDE";
#[cfg(feature = "imzml")]
const IMS_CV: &str = "IMS";

const MS_CV_BYTES: &[u8] = MS_CV.as_bytes();
const UO_CV_BYTES: &[u8] = UO_CV.as_bytes();
const EFO_CV_BYTES: &[u8] = EFO_CV.as_bytes();
const OBI_CV_BYTES: &[u8] = OBI_CV.as_bytes();
const HANCESTRO_CV_BYTES: &[u8] = HANCESTRO_CV.as_bytes();
const BFO_CV_BYTES: &[u8] = BFO_CV.as_bytes();
const BTO_CV_BYTES: &[u8] = BTO_CV.as_bytes();
const NCIT_CV_BYTES: &[u8] = NCIT_CV.as_bytes();
const PRIDE_CV_BYTES: &[u8] = PRIDE_CV.as_bytes();
#[cfg(feature = "imzml")]
const IMS_CV_BYTES: &[u8] = IMS_CV.as_bytes();

impl TryFrom<u8> for ControlledVocabulary {
    type Error = ControlledVocabularyResolutionError;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            1 => Ok(Self::MS),
            2 => Ok(Self::UO),
            3 => Ok(Self::EFO),
            4 => Ok(Self::OBI),
            5 => Ok(Self::HANCESTRO),
            6 => Ok(Self::BFO),
            7 => Ok(Self::NCIT),
            8 => Ok(Self::BTO),
            9 => Ok(Self::PRIDE),
            #[cfg(feature = "imzml")]
            10 => Ok(Self::IMS),
            _ => {
                Err(ControlledVocabularyResolutionError::UnknownControlledVocabularyCode(value))
            }
        }
    }
}

impl<'a> ControlledVocabulary {
    /// Get the CURIE namespace prefix for this controlled vocabulary
    pub const fn prefix(&self) -> Cow<'static, str> {
        match &self {
            Self::MS => Cow::Borrowed(MS_CV),
            Self::UO => Cow::Borrowed(UO_CV),
            Self::EFO => Cow::Borrowed(EFO_CV),
            Self::OBI => Cow::Borrowed(OBI_CV),
            Self::HANCESTRO => Cow::Borrowed(HANCESTRO_CV),
            Self::BFO => Cow::Borrowed(BFO_CV),
            Self::NCIT => Cow::Borrowed(NCIT_CV),
            Self::BTO => Cow::Borrowed(BTO_CV),
            Self::PRIDE => Cow::Borrowed(PRIDE_CV),
            #[cfg(feature = "imzml")]
            Self::IMS => Cow::Borrowed(IMS_CV),
            Self::Unknown => panic!("Cannot encode unknown CV"),
        }
    }

    /// Like [`ControlledVocabulary::prefix`], but obtain a byte string instead
    pub const fn as_bytes(&self) -> &'static [u8] {
        match &self {
            Self::MS => MS_CV_BYTES,
            Self::UO => UO_CV_BYTES,
            Self::EFO => EFO_CV_BYTES,
            Self::OBI => OBI_CV_BYTES,
            Self::HANCESTRO => HANCESTRO_CV_BYTES,
            Self::BFO => BFO_CV_BYTES,
            Self::NCIT => NCIT_CV_BYTES,
            Self::BTO => BTO_CV_BYTES,
            Self::PRIDE => PRIDE_CV_BYTES,
            #[cfg(feature = "imzml")]
            Self::IMS => IMS_CV_BYTES,
            Self::Unknown => panic!("Cannot encode unknown CV"),
        }
    }

    pub const fn as_option(&self) -> Option<Self> {
        match self {
            Self::Unknown => None,
            _ => Some(*self),
        }
    }

    /// Create a [`Param`] whose accession comes from this controlled vocabulary namespace with
    /// an empty value.
    ///
    /// # Arguments
    /// - `accession`: The accession code for the [`Param`]. If specified as a [`CURIE`] or a string-like type,
    ///   any namespace is ignored.
    /// - `name`: The name of the parameter
    /// # See Also
    /// - [`ControlledVocabulary::param_val`]
    pub fn param<A: Into<AccessionLike<'a>>, S: Into<String>>(
        &self,
        accession: A,
        name: S,
    ) -> Param {
        let mut param = Param::new();
        param.controlled_vocabulary = Some(*self);
        param.name = name.into();

        let accession: AccessionLike = accession.into();

        match accession {
            AccessionLike::Text(s) => {
                if let Some(nb) = s.split(':').next_back() {
                    param.accession = Some(nb.parse().unwrap_or_else(|_| {
                        panic!("Expected accession to be numeric, got {}", s)
                    }))
                }
            }
            AccessionLike::Number(n) => param.accession = Some(n),
            AccessionLike::CURIE(c) => param.accession = Some(c.accession),
        }
        param
    }

    pub const fn curie(&self, accession: AccessionIntCode) -> CURIE {
        CURIE::new(*self, accession)
    }

    /// Create a [`ParamCow`] from this namespace in a `const` context, useful for preparing
    /// global constants or inlined variables.
    ///
    /// All parameters must have a `'static` lifetime.
    ///
    /// # Arguments
    /// - `name`: The name of the controlled vocabulary term.
    /// - `value`: The wrapped value as a constant.
    /// - `accession`: The a priori determined accession code for the term
    /// - `unit`: The unit associated with the value
    pub const fn const_param(
        &self,
        name: &'static str,
        value: ValueRef<'static>,
        accession: AccessionIntCode,
        unit: Unit,
    ) -> ParamCow<'static> {
        ParamCow {
            name: Cow::Borrowed(name),
            value,
            accession: Some(accession),
            controlled_vocabulary: Some(*self),
            unit,
        }
    }

    /// Create a [`ParamCow`] from this namespace in a `const` context with an empty
    /// value and no unit.
    ///
    /// See [`ControlledVocabulary::const_param`] for more details.
    pub const fn const_param_ident(
        &self,
        name: &'static str,
        accession: AccessionIntCode,
    ) -> ParamCow<'static> {
        self.const_param(name, ValueRef::Empty, accession, Unit::Unknown)
    }

    /// Create a [`ParamCow`] from this namespace in a `const` context with an empty
    /// value but a specified unit.
    ///
    /// This is intended to create a "template" that will be copied and have a value specified.
    ///
    /// See [`ControlledVocabulary::const_param`] for more details.
    pub const fn const_param_ident_unit(
        &self,
        name: &'static str,
        accession: AccessionIntCode,
        unit: Unit,
    ) -> ParamCow<'static> {
        self.const_param(name, ValueRef::Empty, accession, unit)
    }

    /// Create a [`Param`] whose accession comes from this controlled vocabulary namespace with
    /// the given value.
    ///
    /// # Arguments
    /// - `accession`: The accession code for the [`Param`]. If specified as a [`CURIE`] or a string-like type,
    ///   any namespace is ignored.
    /// - `name`: The name of the parameter
    /// - `value`: The value of the parameter
    ///
    /// # See Also
    /// - [`ControlledVocabulary::param`]
    pub fn param_val<S: Into<String>, A: Into<AccessionLike<'a>>, V: Into<Value>>(
        &self,
        accession: A,
        name: S,
        value: V,
    ) -> Param {
        let mut param = self.param(accession, name);
        param.value = value.into();
        param
    }
}

/// An error describing a failure to map a controlled vocabulary identifier
/// to a known namespace
#[derive(Debug, Clone, Error)]
pub enum ControlledVocabularyResolutionError {
    #[error("Unrecognized controlled vocabulary {0}")]
    UnknownControlledVocabulary(String),
    #[error("Unrecognized controlled vocabulary code {0}")]
    UnknownControlledVocabularyCode(u8),
}

impl FromStr for ControlledVocabulary {
    type Err = ControlledVocabularyResolutionError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "MS" | "PSI-MS" => Ok(Self::MS),
            "UO" => Ok(Self::UO),
            EFO_CV => Ok(Self::EFO),
            OBI_CV => Ok(Self::OBI),
            BFO_CV => Ok(Self::BFO),
            HANCESTRO_CV => Ok(Self::HANCESTRO),
            #[cfg(feature = "imzml")]
            IMS_CV => Ok(Self::IMS),
            _ => Ok(Self::Unknown),
        }
    }
}

pub type AccessionIntCode = u32;
pub type AccessionByteCode7 = [u8; 7];

#[allow(unused)]
#[derive(Debug)]
#[repr(u8)]
enum AccessionCode {
    Int(AccessionIntCode),
    Byte7(AccessionByteCode7),
}

impl Display for AccessionCode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AccessionCode::Int(v) => write!(f, "{v:07}"),
            AccessionCode::Byte7(v) => write!(
                f,
                "{}",
                core::str::from_utf8(v)
                    .map_err(|e| format!("ERROR:{e}"))
                    .unwrap()
            ),
        }
    }
}

#[derive(Debug, thiserror::Error, Clone, PartialEq)]
pub enum AccessionCodeParseError {
    #[error("The acccession code was too long: {0}")]
    AccessionCodeTooLong(String),
    #[error("The acccession code was not in range: {0}")]
    AccessionCodeNotInRange(String),
}

impl FromStr for AccessionCode {
    type Err = AccessionCodeParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() > 7 {
            return Err(AccessionCodeParseError::AccessionCodeTooLong(s.to_string()));
        }
        if !s.is_ascii() {
            return Err(AccessionCodeParseError::AccessionCodeNotInRange(
                s.to_string(),
            ));
        }
        if let Ok(u) = s.parse::<AccessionIntCode>() {
            Ok(Self::Int(u))
        } else {
            let mut bytes = AccessionByteCode7::default();
            for (byte_from, byte_to) in s.as_bytes().iter().rev().zip(bytes.iter_mut().rev()) {
                *byte_to = *byte_from;
            }
            Ok(Self::Byte7(bytes))
        }
    }
}

/// A CURIE is a namespace + accession identifier
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CURIE {
    pub controlled_vocabulary: ControlledVocabulary,
    pub accession: AccessionIntCode,
}

impl CURIE {
    pub const fn new(cv_id: ControlledVocabulary, accession: AccessionIntCode) -> Self {
        Self {
            controlled_vocabulary: cv_id,
            accession,
        }
    }

    pub fn as_param(&self) -> Param {
        let mut param = Param::new();
        param.controlled_vocabulary = Some(self.controlled_vocabulary);
        param.accession = Some(self.accession);
        param
    }

    #[inline(always)]
    pub fn accession_int(&self) -> u32 {
        self.accession
    }

    #[inline(always)]
    pub fn controlled_vocabulary(&self) -> ControlledVocabulary {
        self.controlled_vocabulary
    }
}

impl Display for CURIE {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{:07}",
            self.controlled_vocabulary.prefix(),
            self.accession
        )
    }
}

impl<T: ParamLike> PartialEq<T> for CURIE {
    fn eq(&self, other: &T) -> bool {
        if !other.is_controlled()
            || other
                .controlled_vocabulary()
                .map(|c| c != self.controlled_vocabulary)
                .unwrap_or_default()
        {
            false
        } else {
            other
                .accession()
                .map(|a| a == self.accession)
                .unwrap_or_default()
        }
    }
}

#[derive(Debug, Error)]
pub enum CURIEParsingError {
    #[error("{0} is not a recognized controlled vocabulary")]
    UnknownControlledVocabulary(
        #[from]
        #[source]
        ControlledVocabularyResolutionError,
    ),
    #[error("Failed to parse accession number {0}")]
    AccessionParsingError(
        #[from]
        #[source]
        num::ParseIntError,
    ),
    #[error("Did not detect a namespace separator ':' token")]
    MissingNamespaceSeparator,
}

impl FromStr for CURIE {
    type Err = CURIEParsingError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut tokens = s.split(':');
        let cv = tokens
            .next()
            .ok_or(CURIEParsingError::MissingNamespaceSeparator)?;
        let accession = tokens.next();
        if accession.is_none() {
            Err(CURIEParsingError::MissingNamespaceSeparator)
        } else {
            let cv: ControlledVocabulary = cv.parse::<ControlledVocabulary>()?;

            let accession = accession.unwrap().parse()?;
            Ok(CURIE::new(cv, accession))
        }
    }
}

impl TryFrom<&Param> for CURIE {
    type Error = String;

    fn try_from(value: &Param) -> Result<Self, Self::Error> {
        match (value.controlled_vocabulary, value.accession) {
            (Some(cv), Some(acc)) => Ok(CURIE::new(cv, acc)),
            _ => Err(format!(
                "{} is missing controlled vocabulary or accession",
                value.name()
            )),
        }
    }
}

impl<'a> TryFrom<&ParamCow<'a>> for CURIE {
    type Error = String;

    fn try_from(value: &ParamCow<'a>) -> Result<Self, Self::Error> {
        match (value.controlled_vocabulary, value.accession) {
            (Some(cv), Some(acc)) => Ok(CURIE::new(cv, acc)),
            _ => Err(format!(
                "{} is missing controlled vocabulary or accession",
                value.name()
            )),
        }
    }
}

pub fn curie_to_num(curie: &str) -> (Option<ControlledVocabulary>, Option<AccessionIntCode>) {
    let mut parts = curie.split(':');
    let prefix = match parts.next() {
        Some(v) => v.parse::<ControlledVocabulary>().unwrap().as_option(),
        None => None,
    };
    if let Some(k) = curie.split(':').nth(1) {
        match k.parse() {
            Ok(v) => (prefix, Some(v)),
            Err(_) => (prefix, None),
        }
    } else {
        (prefix, None)
    }
}
