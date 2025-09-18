use std::borrow::{Borrow, Cow};

use super::bindata::{ArrayRetrievalError, ArrayType, BinaryArrayMap, ByteArrayView};
use crate::params::{Param, ParamDescribed};
use crate::spectrum::scan_properties::{
    ChromatogramDescription, ChromatogramType, Precursor, ScanPolarity,
};
use mzpeaks::coordinate::{Time, MZ};
use mzpeaks::feature::{FeatureView, SimpleFeature, TimeInterval};

#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Chromatogram {
    description: ChromatogramDescription,
    pub arrays: BinaryArrayMap,
}

macro_rules! as_feature_view {
    ($chromatogram:ident, $view:ident => $then:tt) => {
        if let Ok(t) = $chromatogram.time() {
            if let Ok(i) = $chromatogram.intensity() {
                let $view = FeatureView::<Time, Time>::new(t.borrow(), t.borrow(), i.borrow());
                Some($then)
            } else {
                None
            }
        } else {
            None
        }
    };
}

#[allow(unused)]
pub(crate) fn as_simple_feature(chromatogram: &Chromatogram) -> Option<SimpleFeature<MZ, Time>> {
    if let Ok(t) = chromatogram.time() {
        if let Ok(i) = chromatogram.intensity() {
            let mut f = SimpleFeature::<MZ, Time>::empty(0.0);
            f.extend(t.iter().zip(i.iter()).map(|(y, z)| (0.0f64, *y, *z)));
            return Some(f);
        }
    }
    None
}

impl TimeInterval<Time> for Chromatogram {
    fn start_time(&self) -> Option<f64> {
        if let Ok(t) = self.time() {
            t.first().copied()
        } else {
            None
        }
    }

    fn end_time(&self) -> Option<f64> {
        if let Ok(t) = self.time() {
            t.last().copied()
        } else {
            None
        }
    }

    fn apex_time(&self) -> Option<f64> {
        as_feature_view!(self, view => {
            view.apex_time()
        })?
    }

    fn area(&self) -> f32 {
        as_feature_view!(self, view => {
            view.area()
        })
        .unwrap_or_default()
    }

    fn iter_time(&self) -> impl Iterator<Item = f64> {
        self.arrays
            .get(&ArrayType::TimeArray)
            .map(|a| a.iter_f64().unwrap())
            .into_iter()
            .flatten()
    }
}


/// Analog of [`SpectrumLike`](crate::spectrum::SpectrumLike) for chromatograms or
/// other measures over time.
pub trait ChromatogramLike: ParamDescribed {
    /// The method to access the spectrum description itself, which supplies
    /// the data for most other methods on this trait.
    fn description(&self) -> &ChromatogramDescription;

    fn description_mut(&mut self) -> &mut ChromatogramDescription;

    /// Access the (first) precursor information, if it exists.
    #[inline]
    fn precursor(&self) -> Option<&Precursor> {
        let desc = self.description();
        desc.precursor.first()
    }

    /// Iterate over all precursors of the spectrum
    fn precursor_iter(&self) -> impl Iterator<Item = &Precursor> {
        let desc = self.description();
        desc.precursor.iter()
    }

    /// Mutably access the (first) precursor information, if it exists
    fn precursor_mut(&mut self) -> Option<&mut Precursor> {
        let desc = self.description_mut();
        desc.precursor.first_mut()
    }

    /// Iterate over all precursors of the spectrum mutably
    fn precursor_iter_mut(&mut self) -> impl Iterator<Item = &mut Precursor> {
        let desc = self.description_mut();
        desc.precursor.iter_mut()
    }

    #[inline]
    fn start_time(&self) -> Option<f64> {
        if let Ok(t) = self.time() {
            t.first().copied()
        } else {
            None
        }
    }

    #[inline]
    fn end_time(&self) -> Option<f64> {
        if let Ok(t) = self.time() {
            t.last().copied()
        } else {
            None
        }
    }

    /// Access the MS exponentiation level
    #[inline]
    fn ms_level(&self) -> Option<u8> {
        self.description().ms_level
    }

    /// Access the native ID string for the spectrum
    #[inline]
    fn id(&self) -> &str {
        &self.description().id
    }

    #[inline]
    fn chromatogram_type(&self) -> ChromatogramType {
        self.description().chromatogram_type
    }

    /// Access the index of the spectrum in the source file
    #[inline]
    fn index(&self) -> usize {
        self.description().index
    }

    #[inline]
    fn polarity(&self) -> ScanPolarity {
        self.description().polarity
    }

    fn is_aggregate(&self) -> bool {
        self.description().is_aggregate()
    }

    fn is_electromagnetic_radiation(&self) -> bool {
        self.description().is_electromagnetic_radiation()
    }

    fn is_ion_current(&self) -> bool {
        self.description().is_ion_current()
    }

    fn time(&self) -> Result<Cow<'_, [f64]>, ArrayRetrievalError>;
    fn intensity(&self) -> Result<Cow<'_, [f32]>, ArrayRetrievalError>;
}

impl Chromatogram {
    pub fn new(description: ChromatogramDescription, arrays: BinaryArrayMap) -> Self {
        Self {
            description,
            arrays,
        }
    }

    pub fn time(&self) -> Result<Cow<'_, [f64]>, ArrayRetrievalError> {
        if let Some(a) = self.arrays.get(&ArrayType::TimeArray) {
            a.to_f64()
        } else {
            Err(ArrayRetrievalError::NotFound(ArrayType::TimeArray))
        }
    }

    pub fn intensity(&self) -> Result<Cow<'_, [f32]>, ArrayRetrievalError> {
        if let Some(a) = self.arrays.get(&ArrayType::IntensityArray) {
            a.to_f32()
        } else {
            Err(ArrayRetrievalError::NotFound(ArrayType::IntensityArray))
        }
    }

    pub fn apex_time(&self) -> Option<f64> {
        TimeInterval::apex_time(&self)
    }

    pub fn area(&self) -> f32 {
        TimeInterval::area(&self)
    }
}

impl ChromatogramLike for Chromatogram {
    fn description(&self) -> &ChromatogramDescription {
        &self.description
    }

    fn time(&self) -> Result<Cow<'_, [f64]>, ArrayRetrievalError> {
        self.time()
    }

    fn intensity(&self) -> Result<Cow<'_, [f32]>, ArrayRetrievalError> {
        self.intensity()
    }

    fn description_mut(&mut self) -> &mut ChromatogramDescription {
        &mut self.description
    }
}

impl ParamDescribed for Chromatogram {
    fn params(&self) -> &[Param] {
        self.description.params()
    }

    fn params_mut(&mut self) -> &mut crate::ParamList {
        self.description.params_mut()
    }
}
