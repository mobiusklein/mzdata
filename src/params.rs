use std::collections::HashMap;
use std::str;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Param {
    pub name: String,
    pub value: String,
    pub accession: String,
    pub controlled_vocabulary: Option<String>,
    pub unit_name: Option<String>,
    pub unit_accession: Option<String>,
}

impl Param {
    pub fn new() -> Param {
        Param {
            ..Default::default()
        }
    }

    pub fn new_key_value(name: String, value: String) -> Param {
        let mut inst = Self::new();
        inst.name = name;
        inst.value = value;
        inst
    }

    pub fn coerce<T: str::FromStr>(&self) -> Result<T, T::Err> {
        self.value.parse::<T>()
    }

    pub fn is_controlled(&self) -> bool {
        self.accession.is_empty()
    }

    pub fn with_unit<S: Into<String>, A: Into<String>>(mut self, accession: S, name: A) -> Param {
        self.unit_accession = Some(accession.into());
        self.unit_name = Some(name.into());
        self
    }
}

#[derive(Debug, Clone)]
pub struct ControlledVocabulary {
    pub prefix: String,
}

impl ControlledVocabulary {
    pub fn new(prefix: String) -> ControlledVocabulary {
        ControlledVocabulary { prefix }
    }

    pub fn param<S: Into<String>, A: Into<String>>(&self, accession: A, name: S) -> Param {
        let mut param = Param::new();
        param.controlled_vocabulary = Some(self.prefix.clone());
        param.name = name.into();
        param.accession = accession.into();
        param
    }

    pub fn param_val<S: Into<String>, A: Into<String>, V: ToString>(
        &self,
        accession: A,
        name: S,
        value: V,
    ) -> Param {
        let mut param = self.param(accession, name);
        param.value = value.to_string();
        param
    }
}

pub type ParamList = Vec<Param>;

pub type ParamMap = HashMap<String, Param>;

pub trait ParamDescribed {
    fn params(&self) -> &ParamList;
    fn params_mut(&mut self) -> &mut ParamList;

    fn add_param(&mut self, param: Param) {
        self.params_mut().push(param);
    }

    fn get_param_by_name(&self, name: &str) -> Option<&Param> {
        for param in self.params().iter() {
            if param.name == name {
                return Some(param);
            }
        }
        None
    }

    fn get_param_by_accession(&self, accession: &str) -> Option<&Param> {
        for param in self.params().iter() {
            if param.accession == accession {
                return Some(param);
            }
        }
        None
    }
}

#[macro_export]
macro_rules! impl_param_described {
    ($($t:ty), +) => {$(

        impl $crate::params::ParamDescribed for $t {
            fn params(&self) -> &$crate::params::ParamList {
                return &self.params
            }

            fn params_mut(&mut self) -> &mut $crate::params::ParamList {
                return &mut self.params
            }
        }
    )+};
}


/// Units that a term's value might have
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Unit {
    // Mass
    MZ,
    Mass,
    PartsPerMillion,
    Nanometer,

    // Time
    Minute,
    Second,
    Hour,
    Millisecond,
    VoltSecondPerSquareCentimeter,

    // Intensity
    DetectorCounts,
    PercentBasePeak,
    PercentBasePeakTimes100,
    AbsorbanceUnit,
    CountsPerSecond,

    // Collision Energy
    Electronvolt,
    PercentElectronVolt,
    Volt
}