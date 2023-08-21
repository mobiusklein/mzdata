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

    fn remove_param(&mut self, index: usize) -> Param {
        self.params_mut().remove(index)
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Unit {
    // Mass
    MZ,
    Mass,
    PartsPerMillion,
    Nanometer,

    // Time
    Minute,
    Second,
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
    Volt,

    Unknown,
}

impl Unit {
    pub fn for_param(&self) -> (&'static str, &'static str) {
        match self {
            Self::Millisecond => ("UO:0000028", "millisecond"),
            Self::Second => ("UO:0000010", "second"),
            Self::Minute => ("UO:0000031", "minute"),

            Self::MZ => ("MS:1000040", "m/z"),
            Self::Mass => ("UO:000221", "dalton"),

            Self::DetectorCounts => ("MS:1000131", "number of detector counts"),
            Self::PercentBasePeak => ("MS:1000132", "percent of base peak"),
            Self::PercentBasePeakTimes100 => ("MS:1000905", "percent of base peak times 100"),
            Self::AbsorbanceUnit => ("UO:0000269", "absorbance unit"),
            Self::CountsPerSecond => ("MS:1000814", "counts per second"),

            Self::Electronvolt => ("UO:0000266", "electronvolt"),
            Self::PercentElectronVolt => ("UO:0000187", "percent"),

            _ => ("", ""),
        }
    }
}
