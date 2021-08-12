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
