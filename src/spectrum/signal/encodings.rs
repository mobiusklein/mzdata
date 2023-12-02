use std::{ops::{Mul, AddAssign}, io, fmt::Display};
use bytemuck::{self, Pod};
use thiserror::{self, Error};

use num_traits::Float;
#[cfg(feature = "numpress")]
use numpress;

use crate::params::{ParamCow, ControlledVocabulary, Unit};

pub type Bytes = Vec<u8>;

pub fn to_bytes<T: Pod>(data: &[T]) -> Bytes {
    bytemuck::cast_slice(data).to_vec()
}

pub fn as_bytes<T: Pod>(data: &[T]) -> &[u8] {
    bytemuck::cast_slice(data)
}

pub fn vec_as_bytes<T: Pod>(data: Vec<T>) -> Bytes {
    bytemuck::cast_vec(data)
}

/// The kinds of data arrays found in mass spectrometry data files governed
/// by the PSI-MS controlled vocabulary.
#[derive(Debug, Clone, PartialEq, Hash, Eq, PartialOrd, Ord, Default)]
pub enum ArrayType {
    #[default]
    Unknown,
    MZArray,
    IntensityArray,
    ChargeArray,
    SignalToNoiseArray,
    TimeArray,
    WavelengthArray,
    IonMobilityArray,
    MeanIonMobilityArray,
    RawIonMobilityArray,
    DeconvolutedIonMobilityArray,
    NonStandardDataArray { name: Box<String> },
}

impl Display for ArrayType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl ArrayType {
    pub const fn preferred_dtype(&self) -> BinaryDataArrayType {
        match self {
            ArrayType::MZArray => BinaryDataArrayType::Float64,
            ArrayType::IntensityArray => BinaryDataArrayType::Float32,
            ArrayType::ChargeArray => BinaryDataArrayType::Int32,
            _ => BinaryDataArrayType::Float32,
        }
    }

    pub const fn as_param_const(&self) -> ParamCow<'static> {
        const CV: ControlledVocabulary = ControlledVocabulary::MS;
        match self {
            ArrayType::MZArray => CV.const_param_ident_unit("m/z array", 1000514, Unit::MZ),
            ArrayType::IntensityArray => CV.const_param_ident_unit("intensity array", 1000515, Unit::DetectorCounts),
            ArrayType::ChargeArray => CV.const_param_ident("charge array", 1000516),
            ArrayType::TimeArray => CV.const_param_ident_unit("time array", 1000595, Unit::Minute),
            ArrayType::RawIonMobilityArray => CV.const_param_ident("raw ion mobility array", 1003007),
            ArrayType::MeanIonMobilityArray => CV.const_param_ident("mean ion mobility array", 1002816),
            ArrayType::DeconvolutedIonMobilityArray => {
                CV.const_param_ident("deconvoluted ion mobility array", 1003154)
            }
            ArrayType::NonStandardDataArray { name: _name } => {
                panic!("Cannot format NonStandardDataArray in a const context");
            }
            _ => {
                panic!("Could not determine how to name for array");
            }
        }
    }

    pub const fn as_param_with_unit_const(&self, unit: Unit) -> ParamCow<'static> {
        const CV: ControlledVocabulary = ControlledVocabulary::MS;
        match self {
            ArrayType::MZArray => CV.const_param_ident_unit("m/z array", 1000514, unit),
            ArrayType::IntensityArray => CV.const_param_ident_unit("intensity array", 1000515, unit),
            ArrayType::ChargeArray => CV.const_param_ident_unit("charge array", 1000516, unit),
            ArrayType::TimeArray => CV.const_param_ident_unit("time array", 1000595, unit),
            ArrayType::RawIonMobilityArray => CV.const_param_ident_unit("raw ion mobility array", 1003007, unit),
            ArrayType::MeanIonMobilityArray => CV.const_param_ident_unit("mean ion mobility array", 1002816, unit),
            ArrayType::DeconvolutedIonMobilityArray => {
                CV.const_param_ident_unit("deconvoluted ion mobility array", 1003154, unit)
            }
            ArrayType::NonStandardDataArray { name: _name } => {
                panic!("Cannot format NonStandardDataArray in a const context");
            }
            _ => {
                panic!("Could not determine how to name for array");
            }
        }
    }
}


/// The canonical primitive data types found in MS data file formats
/// supported by the PSI-MS controlled vocabulary
#[derive(Debug, Clone, Copy, PartialEq, Hash, Eq, Default)]
pub enum BinaryDataArrayType {
    #[default]
    Unknown,
    Float64,
    Float32,
    Int64,
    Int32,
    ASCII,
}

impl Display for BinaryDataArrayType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl BinaryDataArrayType {
    pub const fn size_of(&self) -> usize {
        match self {
            BinaryDataArrayType::Unknown | BinaryDataArrayType::ASCII => 1,
            BinaryDataArrayType::Float32 | BinaryDataArrayType::Int32 => 4,
            BinaryDataArrayType::Float64 | BinaryDataArrayType::Int64 => 8,
        }
    }
}


/// The range of compression and encoding states that a raw byte buffer
/// might be in during different stages of decoding. Other than `[BinaryCompressionType::Decoded]`,
/// these states may or may not include intermediate base64 encoding.
#[derive(Debug, Clone, Copy, PartialEq, Hash, Default)]
pub enum BinaryCompressionType {
    #[default]
    NoCompression,
    Zlib,
    NumpressLinear,
    NumpressSLOF,
    NumpressPIC,
    NumpressLinearZlib,
    NumpressSLOFZlib,
    NumpressPICZlib,
    LinearPrediction,
    DeltaPrediction,
    Decoded,
}

impl BinaryCompressionType {
    /// Generate a user-understandable message about why a compression conversion operation failed
    pub fn unsupported_msg(&self, context: Option<&str>) -> String {
        match context {
            Some(ctx) => format!("Cannot decode array compressed with {:?} ({})", self, ctx),
            None => format!("Cannot decode array compressed with {:?}", self)
        }
    }
}

impl Display for BinaryCompressionType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}


/// A high level set of failure modes that an operation to retrieve a typed memory buffer
/// from a `[BinaryArrayMap]` might encounter. May also be used to represented conversion
/// during reading or writing.
#[derive(Debug, Clone, Error, PartialEq)]
pub enum ArrayRetrievalError {
    #[error("Array type {0:?} not found")]
    NotFound(ArrayType),
    #[error("An error occurred while decompressing: {0}")]
    DecompressionError(String),
    #[error("The requested data type does not match the number of bytes available in the buffer")]
    DataTypeSizeMismatch,
}

impl From<ArrayRetrievalError> for io::Error {
    fn from(value: ArrayRetrievalError) -> Self {
        match value {
            ArrayRetrievalError::NotFound(_) => io::Error::new(io::ErrorKind::NotFound, value),
            ArrayRetrievalError::DecompressionError(e) => io::Error::new(io::ErrorKind::InvalidData, e),
            ArrayRetrievalError::DataTypeSizeMismatch => io::Error::new(io::ErrorKind::InvalidData, value),
        }
    }
}

#[cfg(feature = "numpress")]
impl From<numpress::Error> for ArrayRetrievalError {
    fn from(value: numpress::Error) -> Self {
        ArrayRetrievalError::DecompressionError(value.to_string())
    }
}


pub fn linear_prediction_decoding<F: Float + Mul + AddAssign>(values: &mut [F]) -> &mut [F] {
    if values.len() < 2 {
        return values
    }
    let two = F::from(2.0).unwrap();

    let prev2 = values[1];
    let prev1 = values[2];
    let offset = values[1];

    values.iter_mut().skip(2).fold((prev1, prev2), |(prev1, prev2), current| {
        let tmp = *current + two * prev1 - prev2 - offset;
        let prev1 = *current;
        let prev2 = prev1;
        *current = tmp;
        (prev1, prev2)
    });

    for i in 0..values.len() {
        if i < 2 {
            continue;
        }
        let v = values[i] + two * values[i - 1] - values[i - 2] - values[1];
        values[i] = v;
    }
    values
}


pub fn linear_prediction_encoding<F: Float + Mul<F> + AddAssign>(values: &mut [F]) -> &mut [F] {
    let n = values.len();
    if n < 3 {
        return values
    }
    let offset = values[1];
    let prev2 = values[0];
    let prev1 = values[1];
    let two = F::from(2.0).unwrap();

    values.iter_mut().fold((prev1, prev2), |(prev1, prev2), val| {
        *val += offset - two * prev1 + prev2;
        let tmp = prev1;
        let prev1 = *val + two * prev1 - prev2 - offset;
        let prev2 = tmp;
        (prev1, prev2)
    });
    values
}


pub fn delta_decoding<F: Float + Mul + AddAssign>(values: &mut [F]) -> &mut [F] {
    if values.len() < 2 {
        return values
    }

    let offset = values[0];
    let prev = values[1];

    values.iter_mut().skip(2).fold(prev, |prev, current| {
        *current += prev - offset;
        *current
    });
    values
}


pub fn delta_encoding<F: Float + Mul + AddAssign>(values: &mut [F]) -> &mut [F] {
    let n = values.len();
    if n < 2 {
        return values
    }
    let prev = values[0];
    let offset = values[0];

    let it = values.iter_mut();
    it.skip(1).fold(prev, |prev, current| {
        let tmp = *current;
        *current += offset - prev;
        tmp
    });
    values
}
