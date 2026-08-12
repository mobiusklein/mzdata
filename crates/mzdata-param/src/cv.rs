use std::{
    collections::{HashMap, HashSet, VecDeque},
    fs, path,
};

use super::*;

use bincode::serde::Compat;
use context_error::{BoxedError, CreateError, FullErrorContent, StaticErrorContent};
use mzcv::{CVData, CVIndex, CVSource, CVStructure, HashBufReader};
use serde::{Deserialize, Serialize};
use sha1::Digest;
use std::sync::{Arc, OnceLock};

impl TryFrom<mzcv::ControlledVocabulary> for ControlledVocabulary {
    type Error = ControlledVocabularyResolutionError;

    fn try_from(value: mzcv::ControlledVocabulary) -> Result<Self, Self::Error> {
        Ok(match value {
            mzcv::ControlledVocabulary::MS => Self::MS,
            mzcv::ControlledVocabulary::UO => Self::UO,
            mzcv::ControlledVocabulary::EFO => Self::EFO,
            mzcv::ControlledVocabulary::OBI => Self::OBI,
            mzcv::ControlledVocabulary::HANCESTRO => Self::HANCESTRO,
            mzcv::ControlledVocabulary::BFO => Self::BFO,
            mzcv::ControlledVocabulary::NCIT => Self::NCIT,
            mzcv::ControlledVocabulary::BTO => Self::BTO,
            mzcv::ControlledVocabulary::PRIDE => Self::PRIDE,
            _ => {
                return Err(
                    ControlledVocabularyResolutionError::UnknownControlledVocabulary(
                        value.to_string(),
                    ),
                );
            }
        })
    }
}

impl From<ControlledVocabulary> for mzcv::ControlledVocabulary {
    fn from(value: ControlledVocabulary) -> Self {
        match value {
            ControlledVocabulary::MS => mzcv::ControlledVocabulary::MS,
            ControlledVocabulary::UO => mzcv::ControlledVocabulary::UO,
            ControlledVocabulary::EFO => mzcv::ControlledVocabulary::EFO,
            ControlledVocabulary::OBI => mzcv::ControlledVocabulary::OBI,
            ControlledVocabulary::HANCESTRO => mzcv::ControlledVocabulary::HANCESTRO,
            ControlledVocabulary::BFO => mzcv::ControlledVocabulary::BFO,
            ControlledVocabulary::NCIT => mzcv::ControlledVocabulary::NCIT,
            ControlledVocabulary::BTO => mzcv::ControlledVocabulary::BTO,
            ControlledVocabulary::PRIDE => mzcv::ControlledVocabulary::PRIDE,
            #[cfg(feature = "imzml")]
            ControlledVocabulary::IMS => mzcv::ControlledVocabulary::Unknown,
            ControlledVocabulary::Unknown => mzcv::ControlledVocabulary::Unknown,
        }
    }
}

impl TryFrom<mzcv::Curie> for CURIE {
    type Error = CURIEParsingError;
    fn try_from(value: mzcv::Curie) -> Result<Self, Self::Error> {
        Ok(CURIE::new(
            value.cv.try_into()?,
            match value.accession {
                mzcv::AccessionCode::Numeric(v) => v,
                mzcv::AccessionCode::Alphanumeric(_, _) => "".parse()?,
            },
        ))
    }
}

impl TryFrom<mzcv::OboIdentifier> for CURIE {
    type Error = CURIEParsingError;

    fn try_from(value: mzcv::OboIdentifier) -> Result<Self, Self::Error> {
        match mzcv::Curie::try_from(value) {
            Ok(v) => v.try_into(),
            Err(e) => match e {
                mzcv::CURIEParsingError::UnknownControlledVocabulary => {
                    Err(CURIEParsingError::UnknownControlledVocabulary(
                        ControlledVocabularyResolutionError::UnknownControlledVocabulary(
                            "".to_string(),
                        ),
                    ))
                }
                mzcv::CURIEParsingError::AccessionParsingError(_) => {
                    let v = "".parse::<u32>()?;
                    Ok(CURIE::new(ControlledVocabulary::Unknown, v))
                }
                mzcv::CURIEParsingError::MissingNamespaceSeparator => {
                    Err(CURIEParsingError::MissingNamespaceSeparator)
                }
            },
        }
    }
}

#[derive(Clone, Debug, bincode::Encode, bincode::Decode)]
pub struct MSTerm {
    pub name: Box<str>,
    accession: bincode::serde::Compat<CURIE>,
    pub(crate) parents: HashSet<bincode::serde::Compat<CURIE>>,
    pub synonyms: Vec<Box<str>>,
    pub relationships: Vec<(Box<str>, Box<str>)>,
}

impl MSTerm {
    #[inline(always)]
    pub fn accession(&self) -> CURIE {
        self.accession.0
    }

    #[inline(always)]
    pub fn curie(&self) -> CURIE {
        self.accession()
    }

    #[inline(always)]
    pub fn parents(&self) -> impl ExactSizeIterator<Item = &CURIE> + '_ {
        self.parents.iter().map(|v| &v.0)
    }

    #[inline(always)]
    pub fn is_a(&self, accession: &CURIE) -> bool {
        self.parents.iter().any(|v| v.0 == *accession)
    }

    #[inline(always)]
    pub fn has_value_types(&self) -> impl Iterator<Item = &str> + '_ {
        self.relationships.iter().filter_map(|(r, v)| {
            if r.as_ref() == "has_value_type" {
                Some(v.as_ref())
            } else {
                None
            }
        })
    }
}

impl CVData for MSTerm {
    type Index = CURIE;

    fn index(&self) -> Option<Self::Index> {
        Some(self.accession.0)
    }

    fn name(&self) -> Option<std::borrow::Cow<'_, str>> {
        Some(Cow::Borrowed(self.name.as_ref()))
    }

    fn synonyms(&self) -> impl Iterator<Item = &str> {
        self.synonyms.iter().map(|s| s.as_ref())
    }

    fn parents(&self) -> impl Iterator<Item = &Self::Index> {
        self.parents.iter().map(|c| &c.0)
    }

    fn curie(&self) -> Option<mzcv::Curie> {
        Some(mzcv::Curie {
            cv: self.accession.0.controlled_vocabulary().into(),
            accession: mzcv::AccessionCode::Numeric(self.accession.0.accession_int()),
        })
    }
}

/// A global facade for interacting with the PSI-MS controlled vocabulary
#[derive(Debug)]
pub struct MSVocabulary();

static MS: OnceLock<mzcv::CVIndex<MSVocabulary>> = OnceLock::new();

impl MSVocabulary {
    /// Get the version metadata of the vocabulary
    pub fn version() -> &'static mzcv::CVVersion {
        Self::init().version()
    }

    /// Update the PSI-MS vocabulary to the version at the provided path. This will update
    /// the parsed cache on disk and attempt to set the global in-memory representation to
    /// use this version.
    ///
    /// # Warning
    /// If the controlled vocabulary is already loaded, this will fail to update the current
    /// version in memory, but should still update the version on disk so that it will be read
    /// from when the a new process starts.
    ///
    /// This *does not* update the compile time static version embedded in the binary used by [`Self::init_static`].
    /// That can only be upgraded at build time.
    pub fn update_from_obo<P: AsRef<path::Path>>(
        obo_path: P,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let infile = obo_path.as_ref();
        let cv_reader: HashBufReader<Box<dyn std::io::Read>, sha1::Sha1> =
            if infile.extension().is_some_and(|s| s == "gz") {
                HashBufReader::boxed(flate2::read::GzDecoder::new(fs::File::open(&infile)?))
            } else {
                HashBufReader::boxed(fs::File::open(&infile)?)
            };

        let (ver, data, errs) = match MSVocabulary::parse([cv_reader].into_iter()) {
            Ok(r) => r,
            Err(e) => {
                let root = BoxedError::default().add_underlying_errors(e);
                return Err(Box::new(root));
            }
        };
        for e in errs {
            log::warn!("While updating PSI-MS CV: {e}");
        }
        let mut this = CVIndex::<Self>::empty();
        this.update_from_structure(ver, data)?;
        if let Err(_) = MS.set(this) {
            log::error!("Failed to update global in-memory singleton, value was already initialized. On-disk cache updated.")
        }
        Ok(())
    }

    /// Get or create the global instance of the [`mzcv::CVIndex`] for [`MSVocabulary`]
    pub fn init() -> &'static mzcv::CVIndex<MSVocabulary> {
        MS.get_or_init(|| {
            let (cv, errs) = mzcv::CVIndex::<Self>::init();
            for e in errs {
                log::error!("Error while initializing MS vocabulary database: {e}");
            }
            cv
        })
    }

    /// Get or create the global instance of the [`mzcv::CVIndex`] for [`MSVocabulary`]
    /// using only the statically embedded data.
    pub fn init_static() -> &'static mzcv::CVIndex<MSVocabulary> {
        MS.get_or_init(|| {
            let cv = mzcv::CVIndex::<Self>::init_static();
            cv
        })
    }

    /// Get a term by its [`CURIE`]
    pub fn get(curie: CURIE) -> Option<Arc<MSTerm>> {
        Self::init().get_by_index(&curie)
    }

    /// Get a term by its human-readable name
    pub fn get_by_name(name: &str) -> Option<Arc<MSTerm>> {
        Self::init().get_by_name(name)
    }

    /// Get the immediate parents of this term
    ///
    /// This follows the `is_a` OBO relationship.
    pub fn parents_of(curie: CURIE) -> Option<HashSet<Compat<CURIE>>> {
        Self::get(curie).map(|v| v.parents.clone())
    }

    /// Get all parents of this term, recursively querying the graph
    ///
    /// This follows the `is_a` OBO relationship.
    pub fn parents_of_recursive(curie: CURIE) -> HashSet<Compat<CURIE>> {
        let cv = Self::init();
        let mut acc = HashSet::new();
        let mut queue = VecDeque::from([Compat(curie)]);
        while !queue.is_empty() {
            let idx = queue.pop_front().unwrap();
            if let Some(term) = cv.get_by_index(&idx.0) {
                for p in term.parents.iter() {
                    if acc.contains(p) {
                        continue;
                    }
                    queue.push_back(p.clone());
                }
            }
            acc.insert(idx);
        }
        acc.remove(&Compat(curie));
        acc
    }

    /// Get the immediate children of this term
    ///
    /// This follows the `is_a` OBO relationship.
    pub fn children_of(curie: CURIE) -> Option<&'static [Compat<CURIE>]> {
        let cv = Self::init();
        cv.data()
            .child_map()
            .get(&Compat(curie))
            .map(|v| v.as_slice())
    }

    /// Test if `child` is derived from `parent` along a transitive `is_a` relationship
    pub fn is_child_of(child: CURIE, parent: CURIE) -> bool {
        let cv = Self::init();

        let child_term = match cv.get_by_index(&child) {
            Some(t) => t,
            None => {
                return false;
            }
        };

        let parent = Compat(parent);
        if child_term.parents.contains(&parent) {
            return true;
        }

        let mut seen = HashSet::new();
        let mut queue = VecDeque::from_iter(child_term.parents().copied());

        while !queue.is_empty() {
            let index = queue.pop_front().unwrap();
            seen.insert(index);
            if let Some(term) = cv.get_by_index(&index) {
                for of in term.parents() {
                    if *of == parent.0 {
                        return true;
                    }
                    if seen.contains(of) {
                        continue;
                    }
                    queue.push_back(*of);
                }
            }
        }
        return false;
    }

    /// Get all children of this term, recursively querying the graph
    ///
    /// This follows the `is_a` OBO relationship.
    pub fn children_of_recursive(curie: CURIE) -> HashSet<Compat<CURIE>> {
        let cv = Self::init();
        let mut acc = HashSet::new();
        let mut queue = VecDeque::from([Compat(curie)]);
        while !queue.is_empty() {
            let idx = queue.pop_front().unwrap();
            if let Some(children) = cv.data().child_map().get(&idx).map(|v| v.as_slice()) {
                for p in children.iter() {
                    if acc.contains(p) {
                        continue;
                    }
                    queue.push_back(p.clone());
                }
            }
            acc.insert(idx);
        }
        acc.remove(&Compat(curie));
        acc
    }
}

fn convert_stanza_to_msterm(obj: mzcv::OboStanza) -> Option<Arc<MSTerm>> {
    let mut data = MSTerm {
        accession: bincode::serde::Compat(obj.id.try_into().ok()?),
        name: obj.lines["name"][0].0.clone(),
        synonyms: obj.synonyms.iter().map(|s| s.synonym.clone()).collect(),
        parents: HashSet::new(),
        relationships: Vec::new(),
    };

    for (rel, tp, _ns, _comment) in obj.relationship {
        match rel {
            mzcv::RelationType::IsA => {
                match mzcv::Curie::try_from(tp) {
                    Ok(val) => {
                        let parent = CURIE::try_from(val).unwrap();
                        data.parents.insert(bincode::serde::Compat(parent))
                    }
                    Err(e) => {
                        log::error!("Failed to parse CV identifier: {e}");
                        continue;
                    }
                };
            }
            mzcv::RelationType::Other(rel) => {
                data.relationships
                    .push((rel, tp.to_string().into_boxed_str()));
            }
        }
    }
    Some(Arc::new(data))
}

/// A data structure for walking up and down trees.
#[derive(Debug, bincode::Encode, bincode::Decode)]
pub struct VocabularyData<T: CVData>
where
    T::Index: Serialize + for<'a> Deserialize<'a>,
{
    terms: Vec<Arc<T>>,
    child_map: HashMap<Compat<T::Index>, Vec<Compat<T::Index>>>,
}

impl<T: CVData> IntoIterator for VocabularyData<T>
where
    T::Index: Serialize + for<'a> Deserialize<'a>,
{
    type Item = Arc<T>;

    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.terms.into_iter()
    }
}

impl<T: CVData> Default for VocabularyData<T>
where
    T::Index: Serialize + for<'a> Deserialize<'a>,
{
    fn default() -> Self {
        Self {
            terms: Vec::default(),
            child_map: HashMap::default(),
        }
    }
}

impl<T: CVData> VocabularyData<T>
where
    T::Index: Serialize + for<'a> Deserialize<'a>,
{
    pub fn new(
        terms: Vec<Arc<T>>,
        child_map: HashMap<Compat<T::Index>, Vec<Compat<T::Index>>>,
    ) -> Self {
        Self { terms, child_map }
    }

    pub fn terms(&self) -> &[Arc<T>] {
        &self.terms
    }

    pub fn child_map(&self) -> &HashMap<Compat<T::Index>, Vec<Compat<T::Index>>> {
        &self.child_map
    }
}

impl<T: CVData> CVStructure<T> for VocabularyData<T>
where
    T::Index: Serialize + for<'a> Deserialize<'a>,
{
    fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    fn len(&self) -> usize {
        self.terms.len()
    }

    fn clear(&mut self) {
        self.terms.clear();
        self.child_map.clear();
    }

    type Index = usize;

    fn iter_indexed(&self) -> Self::IterIndexed<'_> {
        self.terms.iter_indexed()
    }

    fn iter_data(&self) -> Self::IterData<'_> {
        self.terms.iter().cloned()
    }

    fn add(&mut self, data: std::sync::Arc<T>) {
        self.terms.push(data);
    }

    fn index(&self, index: Self::Index) -> Option<std::sync::Arc<T>> {
        self.terms.get(index).cloned()
    }

    fn remove(&mut self, index: Self::Index) {
        self.terms.remove(index);
    }

    type IterIndexed<'a>
        = <Vec<Arc<T>> as CVStructure<T>>::IterIndexed<'a>
    where
        Self: 'a;

    type IterData<'a>
        = <Vec<Arc<T>> as CVStructure<T>>::IterData<'a>
    where
        Self: 'a;
}

impl mzcv::CVSource for MSVocabulary {
    type Data = MSTerm;

    type Structure = VocabularyData<MSTerm>;

    fn cv_name() -> &'static str {
        "MS"
    }

    fn files() -> &'static [mzcv::CVFile] {
        &[mzcv::CVFile {
            name: "MS",
            extension: "obo",
            url: Some("http://purl.obolibrary.org/obo/ms.obo"),
            compression: mzcv::CVCompression::None,
        }]
    }

    fn static_data() -> Option<(mzcv::CVVersion, Self::Structure)> {
        #[cfg(feature = "static_data")]
        {
            use bincode::config::Configuration;
            use std::io;
            let buf = io::Cursor::new(include_bytes!("ms.dat"));
            let mut reader = flate2::bufread::GzDecoder::new(buf);
            let cache = bincode::decode_from_std_read::<
                (mzcv::CVVersion, Self::Structure),
                Configuration,
                _,
            >(&mut reader, Configuration::default())
            .unwrap();
            Some(cache)
        }
        #[cfg(not(feature = "static_data"))]
        {
            None
        }
    }

    fn parse(
        mut reader: impl Iterator<Item = mzcv::HashBufReader<Box<dyn std::io::Read>, impl Digest>>,
    ) -> Result<
        (
            mzcv::CVVersion,
            Self::Structure,
            Vec<BoxedError<'static, mzcv::CVError>>,
        ),
        Vec<BoxedError<'static, mzcv::CVError>>,
    > {
        use context_error::CreateError;
        let reader = reader.next().unwrap();
        let obo = mzcv::OboOntology::from_raw(reader).map_err(|e| {
            vec![
                context_error::BoxedError::small(
                    mzcv::CVError::FileCouldNotBeParsed,
                    e.get_short_description(),
                    e.get_long_description(),
                )
                .add_contexts(e.get_contexts().iter().cloned()),
            ]
        })?;
        let version = obo.version();
        let terms: Vec<Arc<MSTerm>> = obo
            .objects
            .into_iter()
            .filter(|o| o.stanza_type == mzcv::OboStanzaType::Term)
            .filter_map(convert_stanza_to_msterm)
            .collect();

        let mut child_map: HashMap<Compat<CURIE>, Vec<Compat<CURIE>>> =
            HashMap::with_capacity(terms.len());
        for term in terms.iter() {
            for parent in term.parents.iter() {
                child_map
                    .entry(parent.clone())
                    .or_default()
                    .push(term.accession.clone());
            }
        }

        let state = VocabularyData::new(terms, child_map);
        Ok((version, state, Vec::new()))
    }
}

pub trait CVTraversal<T: CVSource> {
    /// Helper method to mirror [`mzcv::CVIndex::get_by_index`]
    fn get_by_index(&self, key: <T::Data as CVData>::Index) -> Option<Arc<T::Data>>;

    fn child_map(
        &self,
    ) -> &HashMap<Compat<<T::Data as CVData>::Index>, Vec<Compat<<T::Data as CVData>::Index>>>;

    /// Get the immediate children of this term
    ///
    /// This follows the `is_a` OBO relationship.
    fn children_of(
        &self,
        curie: <T::Data as CVData>::Index,
    ) -> Option<&'_ [Compat<<T::Data as CVData>::Index>]>;

    /// Get the immediate parents of this term
    ///
    /// This follows the `is_a` OBO relationship.
    fn parents_of(
        &self,
        curie: <T::Data as CVData>::Index,
    ) -> Option<HashSet<Compat<<T::Data as CVData>::Index>>>;

    /// Get all parents of this term, recursively querying the graph
    ///
    /// This follows the `is_a` OBO relationship.
    fn parents_of_recursive(
        &self,
        curie: <T::Data as CVData>::Index,
    ) -> HashSet<Compat<<T::Data as CVData>::Index>> {
        let mut acc: HashSet<Compat<<<T as CVSource>::Data as CVData>::Index>> = HashSet::new();
        let mut queue = VecDeque::from([Compat(curie.clone())]);
        while !queue.is_empty() {
            let idx = queue.pop_front().unwrap();
            if let Some(term) = self.get_by_index(idx.0.clone()) {
                for p in term.parents() {
                    let p = Compat(p.clone());
                    if acc.contains(&p) {
                        continue;
                    }
                    queue.push_back(p);
                }
            }
            acc.insert(idx);
        }
        acc.remove(&Compat(curie));
        acc
    }

    /// Get all children of this term, recursively querying the graph
    ///
    /// This follows the `is_a` OBO relationship.
    fn children_of_recursive(
        &self,
        curie: <T::Data as CVData>::Index,
    ) -> HashSet<Compat<<T::Data as CVData>::Index>> {
        let cv = self;
        let mut acc = HashSet::new();
        let mut queue = VecDeque::from([Compat(curie.clone())]);
        while !queue.is_empty() {
            let idx = queue.pop_front().unwrap();
            if let Some(children) = cv.child_map().get(&idx).map(|v| v.as_slice()) {
                for p in children.iter() {
                    if acc.contains(p) {
                        continue;
                    }
                    queue.push_back(p.clone());
                }
            }
            acc.insert(idx);
        }
        acc.remove(&Compat(curie));
        acc
    }

    /// Test if `child` is derived from `parent` along a transitive `is_a` relationship
    fn is_child_of(
        &self,
        child: <T::Data as CVData>::Index,
        parent: <T::Data as CVData>::Index,
    ) -> bool {
        let cv = self;

        let child_term = match cv.get_by_index(child) {
            Some(t) => t,
            None => {
                return false;
            }
        };

        let parent = Compat(parent);
        if child_term.parents().any(|v| *v == parent.0) {
            return true;
        }

        let mut seen = HashSet::new();
        let mut queue = VecDeque::from_iter(child_term.parents().cloned());

        while !queue.is_empty() {
            let index = queue.pop_front().unwrap();
            seen.insert(index.clone());
            if let Some(term) = cv.get_by_index(index) {
                for of in term.parents() {
                    if *of == parent.0 {
                        return true;
                    }
                    if seen.contains(of) {
                        continue;
                    }
                    queue.push_back(of.clone());
                }
            }
        }
        return false;
    }
}

impl CVTraversal<MSVocabulary> for CVIndex<MSVocabulary> {
    fn child_map(
        &self,
    ) -> &HashMap<Compat<<MSTerm as CVData>::Index>, Vec<Compat<<MSTerm as CVData>::Index>>> {
        self.data().child_map()
    }

    fn children_of(
        &self,
        curie: CURIE,
    ) -> Option<&'_ [Compat<<<MSVocabulary as CVSource>::Data as CVData>::Index>]> {
        let state = self.data();
        state.child_map().get(&Compat(curie)).map(|v| &**v)
    }

    fn parents_of(
        &self,
        curie: CURIE,
    ) -> Option<HashSet<Compat<<<MSVocabulary as CVSource>::Data as CVData>::Index>>> {
        self.get_by_index(&curie)
            .as_ref()
            .map(|v| v.parents.clone())
    }

    fn get_by_index(
        &self,
        key: <<MSVocabulary as CVSource>::Data as CVData>::Index,
    ) -> Option<Arc<<MSVocabulary as CVSource>::Data>> {
        self.get_by_index(&key)
    }
}
