use std::borrow::Cow;

use crate::params::Param;
use crate::spectrum::scan_properties::{
    Precursor, ScanPolarity, ChromatogramType, ChromatogramDescription,
};
use super::signal::{ArrayType, BinaryArrayMap, ArrayRetrievalError, ByteArrayView};



#[derive(Debug, Default, Clone)]
pub struct Chromatogram {
    description: ChromatogramDescription,
    pub arrays: BinaryArrayMap
}

pub trait ChromatogramBehavior {
    /// The method to access the spectrum description itself, which supplies
    /// the data for most other methods on this trait.
    fn description(&self) -> &ChromatogramDescription;

    /// Access the precursor information, if it exists.
    #[inline]
    fn precursor(&self) -> Option<&Precursor> {
        let desc = self.description();
        if let Some(precursor) = &desc.precursor {
            Some(precursor)
        } else {
            None
        }
    }

    /// A shortcut method to retrieve the scan start time
    /// of a spectrum.
    #[inline]
    fn start_time(&self) -> f64 {
        todo!()
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
    fn chromatogram_typ(&self) -> ChromatogramType {
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

    #[inline]
    fn params(&self) -> &[Param] {
        &self.description().params
    }

    fn time(&self) -> Result<Cow<'_, [f64]>, ArrayRetrievalError>;
    fn intensity(&self) -> Result<Cow<'_, [f32]>, ArrayRetrievalError>;
}

impl Chromatogram {
    pub fn new(description: ChromatogramDescription, arrays: BinaryArrayMap) -> Self { Self { description, arrays } }

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
}


impl ChromatogramBehavior for Chromatogram {
    fn description(&self) -> &ChromatogramDescription {
        &self.description
    }

    fn time(&self) -> Result<Cow<'_, [f64]>, ArrayRetrievalError> {
        self.time()
    }

    fn intensity(&self) -> Result<Cow<'_, [f32]>, ArrayRetrievalError> {
        self.intensity()
    }
}


