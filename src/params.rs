use std::collections::HashMap;
use std::str;

#[derive(Debug, Clone, Default)]
pub struct Param {
    pub name: String,
    pub value: String,
    pub accession: String,
    pub controlled_vocabulary: Option<String>,
    pub unit_info: Option<String>,
}

impl Param {
    pub fn new() -> Param {
        Param {
            ..Default::default()
        }
    }

    pub fn coerce<T: str::FromStr>(&self) -> Result<T, T::Err> {
        self.value.parse::<T>()
    }
}

pub type ParamList = Vec<Param>;

pub type ParamMap = HashMap<String, Param>;
