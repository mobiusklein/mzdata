#![allow(unused)]
use std::fmt::Display;

use crate::impl_param_described;
use crate::params::{ControlledVocabulary, Param, ParamCow, ParamList};

use super::Software;

/// Describe a data processing method stage tied to a specific piece of [`Software`]
///
/// See <https://peptideatlas.org/tmp/mzML1.1.0.html#processingMethod>
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ProcessingMethod {
    /// The order of application of this method in the processing pipeline
    pub order: i8,
    /// The software used
    pub software_reference: String,
    /// Controlled vocabulary and user parameters
    pub params: ParamList,
}

/// Describe a complete data processing method, a series of [`ProcessingMethod`] transformations
/// through a pipeline of [`Software`].
///
/// See <https://peptideatlas.org/tmp/mzML1.1.0.html#dataProcessing>
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DataProcessing {
    /// The identifier for this data processing pipeline
    pub id: String,
    /// The set of processing steps applied
    pub methods: Vec<ProcessingMethod>,
}

impl_param_described!(ProcessingMethod);

impl DataProcessing {
    pub fn push(&mut self, method: ProcessingMethod) {
        self.methods.push(method)
    }

    pub fn iter(&self) -> std::slice::Iter<ProcessingMethod> {
        self.methods.iter()
    }

    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, ProcessingMethod> {
        self.methods.iter_mut()
    }

    pub fn len(&self) -> usize {
        self.methods.len()
    }

    pub fn is_empty(&self) -> bool {
        self.methods.is_empty()
    }

    pub fn highest_order(&self) -> i8 {
        self.iter().map(|p| p.order).max().unwrap_or_default()
    }

    pub fn remove(&mut self, index: usize) -> ProcessingMethod {
        self.methods.remove(index)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum DataTransformationAction {
    FormatConversion(FormatConversion),
    DataProcessingAction(DataProcessingAction),
    Other(Param),
}

impl Display for DataTransformationAction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Hash)]
pub enum DataProcessingAction {
    Deisotoping,
    ChargeDeconvolution,
    PeakPicking,
    Smoothing,
    BaselineReduction,
    ChargeStateCalculation,
    PrecursorRecalculation,
    IntensityNormalization,
    MZCalibration,
    DataFiltering,
    AdductDeconvolution,
    IonMobilityDeconvolution,
}

impl From<DataProcessingAction> for Param {
    fn from(value: DataProcessingAction) -> Self {
        value.as_param_const().into()
    }
}

impl From<DataProcessingAction> for ParamCow<'static> {
    fn from(value: DataProcessingAction) -> Self {
        value.as_param_const()
    }
}

impl Display for DataProcessingAction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl DataProcessingAction {
    pub const fn as_param_const(&self) -> ParamCow<'static> {
        const CV: ControlledVocabulary = ControlledVocabulary::MS;

        match self {
            DataProcessingAction::Deisotoping => CV.const_param_ident("deisotoping", 1000033),
            DataProcessingAction::ChargeDeconvolution => {
                CV.const_param_ident("charge deconvolution", 1000034)
            }
            DataProcessingAction::PeakPicking => CV.const_param_ident("peak picking", 1000035),
            DataProcessingAction::Smoothing => CV.const_param_ident("smoothing", 1000592),
            DataProcessingAction::BaselineReduction => {
                CV.const_param_ident("baseline reduction", 1000593)
            }
            DataProcessingAction::ChargeStateCalculation => {
                CV.const_param_ident("charge state calculation", 1000778)
            }
            DataProcessingAction::PrecursorRecalculation => {
                CV.const_param_ident("precursor recalculation", 1000780)
            }
            DataProcessingAction::IntensityNormalization => {
                CV.const_param_ident("intensity normalization", 1001484)
            }
            DataProcessingAction::MZCalibration => CV.const_param_ident("m/z calibration", 1001485),
            DataProcessingAction::DataFiltering => CV.const_param_ident("data filtering", 1001486),
            DataProcessingAction::AdductDeconvolution => {
                CV.const_param_ident("adduct deconvolution", 1003220)
            }
            DataProcessingAction::IonMobilityDeconvolution => {
                CV.const_param_ident("ion mobility deconvolution", 1003222)
            }
        }
    }

    pub const fn from_accession(accession: u32) -> Option<DataProcessingAction> {
        match accession {
            1000033 => Some(DataProcessingAction::Deisotoping),
            1000034 => Some(DataProcessingAction::ChargeDeconvolution),
            1000035 => Some(DataProcessingAction::PeakPicking),
            1000592 => Some(DataProcessingAction::Smoothing),
            1000593 => Some(DataProcessingAction::BaselineReduction),
            1000778 => Some(DataProcessingAction::ChargeStateCalculation),
            1000780 => Some(DataProcessingAction::PrecursorRecalculation),
            1001484 => Some(DataProcessingAction::IntensityNormalization),
            1001485 => Some(DataProcessingAction::MZCalibration),
            1001486 => Some(DataProcessingAction::DataFiltering),
            1003220 => Some(DataProcessingAction::AdductDeconvolution),
            1003222 => Some(DataProcessingAction::IonMobilityDeconvolution),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Hash)]
pub enum FormatConversion {
    ConversionToMzML,
    ConversionToMzMLb,
}

impl Display for FormatConversion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl From<FormatConversion> for Param {
    fn from(value: FormatConversion) -> Self {
        value.to_param_const().into()
    }
}

impl From<FormatConversion> for ParamCow<'static> {
    fn from(value: FormatConversion) -> Self {
        value.to_param_const()
    }
}

impl FormatConversion {
    pub const fn to_param_const(&self) -> ParamCow<'static> {
        const CV: ControlledVocabulary = ControlledVocabulary::MS;

        match self {
            FormatConversion::ConversionToMzML => {
                CV.const_param_ident("Conversion to mzML", 1000544)
            }
            FormatConversion::ConversionToMzMLb => {
                CV.const_param_ident("Conversion to mzMLb", 1002839)
            }
        }
    }

    pub const fn from_accession(accession: u32) -> Option<FormatConversion> {
        match accession {
            1000544 => Some(Self::ConversionToMzML),
            1002839 => Some(Self::ConversionToMzMLb),
            _ => None,
        }
    }
}
