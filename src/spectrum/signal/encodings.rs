use std::{fmt::{self, Formatter}, ops::Mul, io};
use bytemuck::{self, Pod};

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
#[derive(Debug, Clone, PartialEq, Hash, Eq, PartialOrd, Ord)]
pub enum ArrayType {
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

impl Default for ArrayType {
    fn default() -> ArrayType {
        ArrayType::Unknown
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
#[derive(Debug, Clone, Copy, PartialEq, Hash)]
pub enum BinaryDataArrayType {
    Unknown,
    Float64,
    Float32,
    Int64,
    Int32,
    ASCII,
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

impl Default for BinaryDataArrayType {
    fn default() -> BinaryDataArrayType {
        BinaryDataArrayType::Unknown
    }
}


/// The range of compression and encoding states that a raw byte buffer
/// might be in during different stages of decoding. Other than `[BinaryCompressionType::Decoded]`,
/// these states may or may not include intermediate base64 encoding.
#[derive(Debug, Clone, Copy, PartialEq, Hash)]
pub enum BinaryCompressionType {
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

impl Default for BinaryCompressionType {
    fn default() -> BinaryCompressionType {
        BinaryCompressionType::NoCompression
    }
}

/// A high level set of failure modes that an operation to retrieve a typed memory buffer
/// from a `[BinaryArrayMap]` might encounter. May also be used to represented conversion
/// during reading or writing.
#[derive(Debug, Clone)]
pub enum ArrayRetrievalError {
    NotFound,
    DecompressionError(String),
    DecodeError,
    DataTypeSizeMismatch,
}

impl std::fmt::Display for ArrayRetrievalError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::error::Error for ArrayRetrievalError {}


impl From<ArrayRetrievalError> for io::Error {
    fn from(value: ArrayRetrievalError) -> Self {
        match value {
            ArrayRetrievalError::NotFound => io::Error::new(io::ErrorKind::NotFound, value),
            ArrayRetrievalError::DecompressionError(e) => io::Error::new(io::ErrorKind::InvalidData, e),
            ArrayRetrievalError::DecodeError => io::Error::new(io::ErrorKind::InvalidData, value),
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


pub fn linear_prediction_decoding<F: Float + Mul<F>>(values: &mut [F]) -> &mut [F] {
    let two = F::from(2.0).unwrap();

    for i in 0..values.len() {
        if i < 2 {
            continue;
        }
        let v = values[i] + two * values[i - 1] - values[i - 2] - values[1];
        values[i] = v;
    }
    values
}


pub fn linear_prediction_encoding<F: Float + Mul<F>>(values: &mut [F]) -> &mut [F] {
    let n = values.len();
    if n < 3 {
        return values
    }
    let offset = values[1];
    let mut prev2 = values[0];
    let mut prev1 = values[1];
    let two = F::from(2.0).unwrap();
    for i in 0..values.len() {
        values[i] = offset + values[i] - two * prev1 + prev2;
        let tmp = prev1;
        prev1 = values[i] + two * prev1 - prev2 - offset;
        prev2 = tmp;
    }
    values
}


pub fn delta_decoding<F: Float + Mul<F>>(values: &mut [F]) -> &mut [F] {
    for i in 0..values.len() {
        if i < 2 {
            continue;
        }
        values[i] = values[i] + values[i - 1] - values[0];
    }
    values
}


pub fn delta_encoding<F: Float + Mul<F>>(values: &mut [F]) -> &mut [F] {
    let n = values.len();
    if n < 2 {
        return values
    }
    let mut prev = values[0];
    let offset = values[0];
    for i in 1..n {
        let tmp = values[i];
        values[i] = offset + values[i] - prev;
        prev = tmp;
    }
    values
}
