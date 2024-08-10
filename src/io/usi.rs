use std::fmt::Display;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum USIParseError {
    #[error("Protocol {0} is not recognized. Expected `mzspec`, got {1}")]
    UnknownProtocol(String, String),
    #[error("Did not find a dataset section in {0}")]
    MissingDataset(String),
    #[error("Did not find a run name in {0}")]
    MissingRun(String)
}


#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Protocol {
    MZSpec
}

impl Display for Protocol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MZSpec => write!(f, "mzspec")
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Identifier {
    Scan(u64),
    Index(u64),
}


#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct USI {
    pub protocol: Protocol,
    pub dataset: String,
    pub run_name: String,
    pub identifier: Option<Identifier>,
}

impl Display for USI {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(ident) = self.identifier {
            let (ident_class, ident_val) = match ident {
                Identifier::Scan(i) => ("scan", i),
                Identifier::Index(i) => ("index", i),
            };
            write!(f, "{}:{}:{}:{ident_class}:{ident_val}", self.protocol, self.dataset, self.run_name)
        } else {
            write!(f, "{}:{}:{}", self.protocol, self.dataset, self.run_name)
        }
    }
}