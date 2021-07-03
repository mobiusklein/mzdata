use std::fmt;
use std::cmp;
use std::hash;
use std::option::Option;

use crate::coordinate::{CoordinateLike, MZ};

#[derive(Default, Copy, Clone, Debug)]
pub struct Peak {
    pub mz: f64,
    pub intensity: f32,
    pub index: u32
}

impl fmt::Display for Peak {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Peak({}, {}, {})", self.mz, self.intensity, self.index)
    }
}

impl hash::Hash for Peak {
    fn hash<H: hash::Hasher>(&self, state: &mut H) {
        let mz_val: i64 = self.mz.round() as i64;
        mz_val.hash(state);
    }
}

impl cmp::PartialOrd<Peak> for Peak {
    fn partial_cmp(&self, other: &Peak) -> Option<cmp::Ordering> {
        return self.mz.partial_cmp(&other.mz);
    }
}

impl cmp::PartialEq<Peak> for Peak {
    fn eq(&self, other: &Peak) -> bool {
        if (self.mz - other.mz).abs() > 1e-3 {
            return false;
        }
        else if (self.intensity - other.intensity).abs() > 1e-3 {
            return false;
        }
        return true;
    }

    fn ne(&self, other: &Peak) -> bool {
        return !(self == other)
    }
}

impl CoordinateLike<MZ> for Peak {
    fn get_coordinate(&self) -> f64 {
        return self.mz
    }
}
