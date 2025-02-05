use bytemuck::{self, Pod};
use std::{
    fmt::Display,
    io,
    ops::{AddAssign, Mul},
};
use thiserror::{self, Error};

use num_traits::Float;
#[cfg(feature = "numpress")]
use numpress;

use crate::{
    params::{ControlledVocabulary, ParamCow, Unit},
    Param,
};

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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
    BaselineArray,
    ResolutionArray,
    PressureArray,
    TemperatureArray,
    FlowRateArray,
    NonStandardDataArray {
        name: Box<String>,
    },
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

    /// Create a [`ArrayType::NonStandardDataArray`] with the provided name.
    pub fn nonstandard<S: ToString>(name: S) -> ArrayType {
        ArrayType::NonStandardDataArray {
            name: name.to_string().into(),
        }
    }

    /// Test if the the array describes an ion mobility quantity.
    pub const fn is_ion_mobility(&self) -> bool {
        matches!(
            self,
            Self::IonMobilityArray
                | Self::MeanIonMobilityArray
                | Self::RawIonMobilityArray
                | Self::DeconvolutedIonMobilityArray
        )
    }

    pub fn as_param(&self, unit: Option<Unit>) -> Param {
        const CV: ControlledVocabulary = ControlledVocabulary::MS;
        match self {
            ArrayType::MZArray => CV
                .const_param_ident_unit("m/z array", 1000514, unit.unwrap_or(Unit::MZ))
                .into(),
            ArrayType::IntensityArray => CV
                .const_param_ident_unit(
                    "intensity array",
                    1000515,
                    unit.unwrap_or(Unit::DetectorCounts),
                )
                .into(),
            ArrayType::ChargeArray => CV.const_param_ident("charge array", 1000516).into(),
            ArrayType::TimeArray => CV
                .const_param_ident_unit("time array", 1000595, unit.unwrap_or(Unit::Minute))
                .into(),
            ArrayType::WavelengthArray => CV.const_param_ident_unit("wavelength array", 1000617, Unit::Nanometer).into(),
            ArrayType::SignalToNoiseArray => CV.const_param_ident("signal to noise array", 1000517).into(),
            ArrayType::IonMobilityArray => CV.const_param_ident_unit("ion mobility array", 1002893, unit.unwrap_or_default()).into(),
            ArrayType::RawIonMobilityArray => CV
                .const_param_ident_unit("raw ion mobility array", 1003007, unit.unwrap_or_default())
                .into(),
            ArrayType::MeanIonMobilityArray => CV
                .const_param_ident_unit(
                    "mean ion mobility array",
                    1002816,
                    unit.unwrap_or_default(),
                )
                .into(),
            ArrayType::DeconvolutedIonMobilityArray => CV
                .const_param_ident_unit(
                    "deconvoluted ion mobility array",
                    1003154,
                    unit.unwrap_or_default(),
                )
                .into(),
            ArrayType::NonStandardDataArray { name } => {
                let mut p = CV.param_val(1000786, "non-standard data array", name.to_string());
                p.unit = unit.unwrap_or_default();
                p
            }
            ArrayType::BaselineArray => {
                CV.const_param_ident("baseline array", 1002530).into()
            },
            ArrayType::ResolutionArray => {
                CV.const_param_ident("resolution array", 1002529).into()
            },
            ArrayType::PressureArray => {
                let mut p = CV.const_param_ident("pressure array", 1000821);
                p.unit = unit.unwrap_or_default();
                p.into()
            },
            ArrayType::TemperatureArray => {
                let mut p = CV.const_param_ident("temperature array", 1000822);
                p.unit = unit.unwrap_or_default();
                p.into()
            },
            ArrayType::FlowRateArray => {
                let mut p = CV.const_param_ident("flow rate array", 1000820);
                p.unit = unit.unwrap_or_default();
                p.into()
            },
            _ => {
                panic!("Could not determine how to name for array {}", self);
            }
        }
    }

    pub const fn as_param_const(&self) -> ParamCow<'static> {
        const CV: ControlledVocabulary = ControlledVocabulary::MS;
        match self {
            ArrayType::MZArray => CV.const_param_ident_unit("m/z array", 1000514, Unit::MZ),
            ArrayType::IntensityArray => {
                CV.const_param_ident_unit("intensity array", 1000515, Unit::DetectorCounts)
            }
            ArrayType::ChargeArray => CV.const_param_ident("charge array", 1000516),
            ArrayType::TimeArray => CV.const_param_ident_unit("time array", 1000595, Unit::Minute),
            ArrayType::WavelengthArray => CV.const_param_ident_unit("wavelength array", 1000617, Unit::Nanometer),
            ArrayType::SignalToNoiseArray => CV.const_param_ident("signal to noise array", 1000517),
            ArrayType::IonMobilityArray => CV.const_param_ident("ion mobility array", 1002893),
            ArrayType::RawIonMobilityArray => {
                CV.const_param_ident("raw ion mobility array", 1003007)
            }
            ArrayType::MeanIonMobilityArray => {
                CV.const_param_ident("mean ion mobility array", 1002816)
            }
            ArrayType::DeconvolutedIonMobilityArray => {
                CV.const_param_ident("deconvoluted ion mobility array", 1003154)
            }
            ArrayType::NonStandardDataArray { name: _name } => {
                panic!(
                    "Cannot format NonStandardDataArray in a const context, please use `as_param`"
                );
            }
            ArrayType::BaselineArray => {
                CV.const_param_ident("baseline array", 1002530)
            },
            ArrayType::ResolutionArray => {
                CV.const_param_ident("resolution array", 1002529)
            },
            ArrayType::PressureArray => {
                CV.const_param_ident("pressure array", 1000821)
            },
            ArrayType::TemperatureArray => {
                CV.const_param_ident("temperature array", 1000822)
            },
            ArrayType::FlowRateArray => {
                CV.const_param_ident("flow rate array", 1000820)
            },
            _ => {
                panic!("Could not determine how to name for array");
            }
        }
    }

    pub const fn as_param_with_unit_const(&self, unit: Unit) -> ParamCow<'static> {
        const CV: ControlledVocabulary = ControlledVocabulary::MS;
        match self {
            ArrayType::MZArray => CV.const_param_ident_unit("m/z array", 1000514, unit),
            ArrayType::IntensityArray => {
                CV.const_param_ident_unit("intensity array", 1000515, unit)
            }
            ArrayType::ChargeArray => CV.const_param_ident_unit("charge array", 1000516, unit),
            ArrayType::TimeArray => CV.const_param_ident_unit("time array", 1000595, unit),
            ArrayType::RawIonMobilityArray => {
                CV.const_param_ident_unit("raw ion mobility array", 1003007, unit)
            }
            ArrayType::MeanIonMobilityArray => {
                CV.const_param_ident_unit("mean ion mobility array", 1002816, unit)
            }
            ArrayType::DeconvolutedIonMobilityArray => {
                CV.const_param_ident_unit("deconvoluted ion mobility array", 1003154, unit)
            }
            ArrayType::NonStandardDataArray { name: _name } => {
                panic!(
                    "Cannot format NonStandardDataArray in a const context, please use `as_param`"
                );
            }
            ArrayType::BaselineArray => {
                CV.const_param_ident_unit("baseline array", 1002530, unit)
            },
            ArrayType::ResolutionArray => {
                CV.const_param_ident_unit("resolution array", 1002529, unit)
            },
            ArrayType::PressureArray => {
                CV.const_param_ident_unit("pressure array", 1000821, unit)
            },
            ArrayType::TemperatureArray => {
                CV.const_param_ident_unit("temperature array", 1000822, unit)
            },
            ArrayType::FlowRateArray => {
                CV.const_param_ident_unit("flow rate array", 1000820, unit)
            },
            _ => {
                panic!("Could not determine how to name for array");
            }
        }
    }
}

/// The canonical primitive data types found in MS data file formats
/// supported by the PSI-MS controlled vocabulary
#[derive(Debug, Clone, Copy, PartialEq, Hash, Eq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
/// might be in during different stages of decoding. Other than `Decoded`,
/// these states may or may not include intermediate base64 encoding.
#[derive(Debug, Clone, Copy, PartialEq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
            None => format!("Cannot decode array compressed with {:?}", self),
        }
    }

    pub const fn as_param(&self) -> Option<ParamCow> {
        let (name, accession) = match self {
            BinaryCompressionType::NoCompression => ("no compression", 1000576),
            BinaryCompressionType::Zlib => ("zlib compression", 1000574),
            BinaryCompressionType::NumpressLinear => {
                ("MS-Numpress linear prediction compression", 1002312)
            }
            BinaryCompressionType::NumpressSLOF => {
                ("MS-Numpress positive integer compression", 1002313)
            }
            BinaryCompressionType::NumpressPIC => {
                ("MS-Numpress short logged float compression", 1002314)
            }
            BinaryCompressionType::NumpressLinearZlib => (
                "MS-Numpress linear prediction compression followed by zlib compression",
                1002746,
            ),
            BinaryCompressionType::NumpressSLOFZlib => (
                "MS-Numpress positive integer compression followed by zlib compression",
                1002477,
            ),
            BinaryCompressionType::NumpressPICZlib => (
                "MS-Numpress short logged float compression followed by zlib compression",
                1002478,
            ),
            BinaryCompressionType::LinearPrediction => todo!(),
            BinaryCompressionType::DeltaPrediction => todo!(),
            BinaryCompressionType::Decoded => return None,
        };
        Some(ControlledVocabulary::MS.const_param_ident(name, accession))
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

impl From<bytemuck::PodCastError> for ArrayRetrievalError {
    fn from(value: bytemuck::PodCastError) -> Self {
        match value {
            bytemuck::PodCastError::TargetAlignmentGreaterAndInputNotAligned => {
                Self::DataTypeSizeMismatch
            }
            bytemuck::PodCastError::OutputSliceWouldHaveSlop => Self::DataTypeSizeMismatch,
            bytemuck::PodCastError::SizeMismatch => Self::DataTypeSizeMismatch,
            bytemuck::PodCastError::AlignmentMismatch => Self::DataTypeSizeMismatch,
        }
    }
}

impl From<ArrayRetrievalError> for io::Error {
    fn from(value: ArrayRetrievalError) -> Self {
        match value {
            ArrayRetrievalError::NotFound(_) => io::Error::new(io::ErrorKind::NotFound, value),
            ArrayRetrievalError::DecompressionError(e) => {
                io::Error::new(io::ErrorKind::InvalidData, e)
            }
            ArrayRetrievalError::DataTypeSizeMismatch => {
                io::Error::new(io::ErrorKind::InvalidData, value)
            }
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
        return values;
    }
    let two = F::from(2.0).unwrap();

    let prev2 = values[1];
    let prev1 = values[2];
    let offset = values[1];

    values
        .iter_mut()
        .skip(2)
        .fold((prev1, prev2), |(prev1, prev2), current| {
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
        return values;
    }
    let offset = values[1];
    let prev2 = values[0];
    let prev1 = values[1];
    let two = F::from(2.0).unwrap();

    values
        .iter_mut()
        .fold((prev1, prev2), |(prev1, prev2), val| {
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
        return values;
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
        return values;
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_dtype_size() {
        assert_eq!(BinaryDataArrayType::ASCII.size_of(), 1);
        assert_eq!(BinaryDataArrayType::Float32.size_of(), 4);
        assert_eq!(BinaryDataArrayType::Int32.size_of(), 4);
        assert_eq!(BinaryDataArrayType::Float64.size_of(), 8);
        assert_eq!(BinaryDataArrayType::Int64.size_of(), 8);
    }

    #[test]
    fn test_array_type_param() {
        let array_types = [
            ArrayType::MZArray,
            ArrayType::IntensityArray,
            ArrayType::ChargeArray,
            ArrayType::SignalToNoiseArray,
            ArrayType::TimeArray,
            ArrayType::WavelengthArray,
            ArrayType::IonMobilityArray,
            ArrayType::MeanIonMobilityArray,
            ArrayType::RawIonMobilityArray,
            ArrayType::DeconvolutedIonMobilityArray,
        ];

        for at in array_types {
            assert_eq!(at.as_param_const().name, at.as_param(None).name)
        }
    }

    #[test]
    fn test_binary_encoding_conv() {
        let encodings = [
            BinaryCompressionType::Decoded,
            BinaryCompressionType::NoCompression,
            BinaryCompressionType::NumpressLinear,
            BinaryCompressionType::NumpressLinearZlib,
            BinaryCompressionType::NumpressPIC,
            BinaryCompressionType::NumpressPICZlib,
            BinaryCompressionType::NumpressSLOF,
            BinaryCompressionType::NumpressSLOFZlib,
            BinaryCompressionType::Zlib,
        ];

        for enc in encodings {
            let reps = match enc {
                BinaryCompressionType::NoCompression => ("no compression", 1000576),
                BinaryCompressionType::Zlib => ("zlib compression", 1000574),
                BinaryCompressionType::NumpressLinear => {
                    ("MS-Numpress linear prediction compression", 1002312)
                }
                BinaryCompressionType::NumpressSLOF => {
                    ("MS-Numpress positive integer compression", 1002313)
                }
                BinaryCompressionType::NumpressPIC => {
                    ("MS-Numpress short logged float compression", 1002314)
                }
                BinaryCompressionType::NumpressLinearZlib => (
                    "MS-Numpress linear prediction compression followed by zlib compression",
                    1002746,
                ),
                BinaryCompressionType::NumpressSLOFZlib => (
                    "MS-Numpress positive integer compression followed by zlib compression",
                    1002477,
                ),
                BinaryCompressionType::NumpressPICZlib => (
                    "MS-Numpress short logged float compression followed by zlib compression",
                    1002478,
                ),
                _ =>  ("", 0),
            };
            if let Some(p) = enc.as_param() {
                assert_eq!(p.name, reps.0);
                assert_eq!(p.accession.unwrap(), reps.1);
            }
        }
    }
}
