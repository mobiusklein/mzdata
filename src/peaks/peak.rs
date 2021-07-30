use std::cmp;
use std::fmt;
use std::hash;

use crate::peaks::coordinate::{CoordinateLike, IndexedCoordinate, Mass, MZ};

#[derive(Default, Clone, Debug)]
pub struct CentroidPeak {
    pub mz: f64,
    pub intensity: f32,
    pub index: u32,
}

impl fmt::Display for CentroidPeak {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "CentroidPeak({}, {}, {})",
            self.mz, self.intensity, self.index
        )
    }
}

impl hash::Hash for CentroidPeak {
    fn hash<H: hash::Hasher>(&self, state: &mut H) {
        let mz_val: i64 = self.mz.round() as i64;
        mz_val.hash(state);
    }
}

impl cmp::PartialOrd<CentroidPeak> for CentroidPeak {
    fn partial_cmp(&self, other: &CentroidPeak) -> Option<cmp::Ordering> {
        self.mz.partial_cmp(&other.mz)
    }
}

impl cmp::PartialEq<CentroidPeak> for CentroidPeak {
    fn eq(&self, other: &CentroidPeak) -> bool {
        if (self.mz - other.mz).abs() > 1e-3 || (self.intensity - other.intensity).abs() > 1e-3 {
            return false;
        }
        true
    }
}

impl CoordinateLike<MZ> for CentroidPeak {
    #[inline]
    fn coordinate(&self) -> f64 {
        self.mz
    }
}

impl IndexedCoordinate<MZ> for CentroidPeak {
    fn get_index(&self) -> u32 {
        self.index
    }

    fn set_index(&mut self, index: u32) {
        self.index = index;
    }
}

pub trait IntensityMeasurement {
    fn intensity(&self) -> f32;
}

impl IntensityMeasurement for CentroidPeak {
    fn intensity(&self) -> f32 {self.intensity}
}

pub trait CentroidLike : IndexedCoordinate<MZ> + IntensityMeasurement {
    fn as_centroid(&self) -> CentroidPeak {
        CentroidPeak { mz: self.coordinate(), intensity: self.intensity(), index: self.get_index() }
    }
}

impl<T: IndexedCoordinate<MZ> + IntensityMeasurement> CentroidLike for T {}



#[derive(Default, Clone, Debug)]
pub struct DeconvolutedPeak {
    pub neutral_mass: f64,
    pub intensity: f32,
    pub charge: i32,
    pub index: u32,
}

impl fmt::Display for DeconvolutedPeak {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "DeconvolutedPeak({}, {}, {}, {})",
            self.neutral_mass, self.intensity, self.charge, self.index
        )
    }
}

impl hash::Hash for DeconvolutedPeak {
    fn hash<H: hash::Hasher>(&self, state: &mut H) {
        let neutral_mass: i64 = self.neutral_mass.round() as i64;
        neutral_mass.hash(state);
    }
}

impl cmp::PartialOrd<DeconvolutedPeak> for DeconvolutedPeak {
    fn partial_cmp(&self, other: &DeconvolutedPeak) -> Option<cmp::Ordering> {
        self.neutral_mass.partial_cmp(&other.neutral_mass)
    }
}

impl cmp::PartialEq<DeconvolutedPeak> for DeconvolutedPeak {
    fn eq(&self, other: &DeconvolutedPeak) -> bool {
        if (self.neutral_mass - other.neutral_mass).abs() > 1e-3
            || self.charge != other.charge
            || (self.intensity - other.intensity).abs() > 1e-3
        {
            return false;
        }
        true
    }
}

impl CoordinateLike<Mass> for DeconvolutedPeak {
    #[inline]
    fn coordinate(&self) -> f64 {
        self.neutral_mass
    }
}

impl IndexedCoordinate<Mass> for DeconvolutedPeak {
    fn get_index(&self) -> u32 {
        self.index
    }

    fn set_index(&mut self, index: u32) {
        self.index = index;
    }
}

impl CoordinateLike<MZ> for DeconvolutedPeak {
    #[inline]
    fn coordinate(&self) -> f64 {
        let charge_carrier: f64 = 1.007276;
        let charge = self.charge as f64;
        (self.neutral_mass - charge_carrier * charge) / charge
    }
}

impl IntensityMeasurement for DeconvolutedPeak {
    fn intensity(&self) -> f32 {self.intensity}
}

pub trait KnownCharge {
    fn charge(&self) -> i32;
}

impl KnownCharge for DeconvolutedPeak {
    fn charge(&self) -> i32 {self.charge}
}

pub trait DeconvolutedCentroid : IndexedCoordinate<Mass> + IntensityMeasurement + KnownCharge {
    fn as_centroid(&self) -> DeconvolutedPeak {
        DeconvolutedPeak {
            neutral_mass: self.coordinate(),
            intensity: self.intensity(),
            charge: self.charge(),
            index: self.get_index()
        }
    }
}

impl<T: IndexedCoordinate<Mass> + IntensityMeasurement + KnownCharge> DeconvolutedCentroid for T {}