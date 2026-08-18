//! The majority of this code is adapted from https://github.com/jspaezp/timsrust-calibration.
//!
//! It is replicated to be compatible with `timsrust` v0.4.1 instead of v0.5+ which introduces
//! greater complexity, and to avoid adding a second SQLite3 implementation.
use rusqlite::Connection;
use thiserror::Error;
use timsrust::converters::{ConvertableDomain, Tof2MzConverter, Scan2ImConverter};

use super::sql::{FromSQL, SQLFrame};

fn require_at<T: rusqlite::types::FromSql>(
    row: &rusqlite::Row<'_>,
    index: usize,
    table: &str,
    column: &str,
) -> Result<T, rusqlite::Error> {
    match row.get::<usize, T>(index) {
        Ok(value) => Ok(value),
        Err(_) => Err(rusqlite::Error::InvalidColumnName(format!(
            "{table} did not contain {column} at index {index}"
        ))),
    }
}

/// One row of the `MzCalibration` table (physical TOF->m/z model).
///
/// Only the columns this crate's converters actually consume are kept; see
/// the module docs for the columns dropped as unused (`T2`, `dC2`, `C2`,
/// `C3`, `C4`).
#[derive(Clone, Debug, PartialEq)]
pub struct MzCalibration {
    /// `MzCalibration.Id` (primary key), referenced by `Frames.MzCalibration`.
    pub id: u32,
    /// `MzCalibration.ModelType`; only `1` is supported.
    pub model_type: u8,
    /// `MzCalibration.DigitizerTimebase`.
    pub digitizer_timebase: f64,
    /// `MzCalibration.DigitizerDelay`.
    pub digitizer_delay: f64,
    /// `MzCalibration.T1`, the reference digitizer temperature the
    /// calibration was fit at.
    pub t1: f64,
    /// `MzCalibration.dC1`, the per-degree drift coefficient for `C1`.
    pub dc1: f64,
    /// `MzCalibration.C0`
    pub c0: Option<f64>,
    /// `MzCalibration.C1`
    pub c1: Option<f64>,
}

impl MzCalibration {
    pub fn new(
        id: u32,
        model_type: u8,
        digitizer_timebase: f64,
        digitizer_delay: f64,
        t1: f64,
        dc1: f64,
        c0: Option<f64>,
        c1: Option<f64>,
    ) -> Self {
        Self {
            id,
            model_type,
            digitizer_timebase,
            digitizer_delay,
            t1,
            dc1,
            c0,
            c1,
        }
    }
}

impl FromSQL for MzCalibration {
    fn from_row(row: &rusqlite::Row<'_>) -> Result<Self, rusqlite::Error> {
        const TABLE_NAME: &str = "MzCalibration";
        Ok(MzCalibration::new(
            row.get(0).unwrap_or_default(),
            row.get(1).unwrap_or_default(),
            require_at(row, 2, TABLE_NAME, "DigitizerTimebase")?,
            require_at(row, 3, TABLE_NAME, "DigitizerDelay")?,
            require_at(row, 4, TABLE_NAME, "T1")?,
            require_at(row, 5, TABLE_NAME, "dC1")?,
            require_at(row, 6, TABLE_NAME, "C0")?,
            require_at(row, 7, TABLE_NAME, "C1")?,
        ))
    }

    fn get_sql() -> String {
        "SELECT Id, ModelType, DigitizerTimebase, DigitizerDelay, T1, dC1, C0, C1 FROM MzCalibration".into()
    }
}

/// One row of the `TimsCalibration` table (physical scan->1/K0 mobility model).
///
/// Only the columns this crate's converters actually consume are kept; see
/// the module docs for the columns dropped as unused (`C5`, `C8`, `C9`).
#[derive(Clone, Debug, PartialEq)]
pub struct TimsCalibration {
    /// `TimsCalibration.Id` (primary key), referenced by
    /// `Frames.TimsCalibration`.
    pub id: u32,
    /// `TimsCalibration.ModelType`; only `2` is supported
    pub model_type: u8,
    /// `TimsCalibration.C0`
    pub c0: Option<f64>,
    /// `TimsCalibration.C1`
    pub c1: Option<f64>,
    /// `TimsCalibration.C2`
    pub c2: Option<f64>,
    /// `TimsCalibration.C3`
    pub c3: Option<f64>,
    /// `TimsCalibration.C4`
    pub c4: Option<f64>,
    /// `TimsCalibration.C6`
    pub c6: Option<f64>,
    /// `TimsCalibration.C7`
    pub c7: Option<f64>,
}

impl TimsCalibration {
    pub fn new(
        id: u32,
        model_type: u8,
        c0: Option<f64>,
        c1: Option<f64>,
        c2: Option<f64>,
        c3: Option<f64>,
        c4: Option<f64>,
        c6: Option<f64>,
        c7: Option<f64>,
    ) -> Self {
        Self {
            id,
            model_type,
            c0,
            c1,
            c2,
            c3,
            c4,
            c6,
            c7,
        }
    }
}

impl FromSQL for TimsCalibration {
    fn from_row(row: &rusqlite::Row<'_>) -> Result<Self, rusqlite::Error> {
        const TABLE_NAME: &str = "TimsCalibration";
        Ok(Self::new(
            row.get(0).unwrap_or_default(),
            require_at(row, 1, TABLE_NAME, "ModelType")?,
            row.get(2).ok(),
            row.get(3).ok(),
            row.get(4).ok(),
            row.get(5).ok(),
            row.get(6).ok(),
            row.get(7).ok(),
            row.get(8).ok(),
        ))
    }

    fn get_sql() -> String {
        "SELECT Id, ModelType, C0, C1, C2, C3, C4, C6, C7 FROM TimsCalibration".into()
    }
}

/// The parameters models for converting indices to m/z or ion mobility.
///
/// There may be multiple models for each dimension, but any given frame
/// uses a single model for each dimension. Not all models use all parameters
/// and not all model types are supported.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CalibrationParameters {
    pub mz: Vec<MzCalibration>,
    pub tims: Vec<TimsCalibration>,
}

impl CalibrationParameters {
    pub fn new(mz: Vec<MzCalibration>, tims: Vec<TimsCalibration>) -> Self {
        Self { mz, tims }
    }

    pub fn from_sql(connection: &Connection) -> Result<Self, rusqlite::Error> {
        let mz = MzCalibration::read_from(connection, [])?;
        let tims = TimsCalibration::read_from(connection, [])?;
        Ok(Self::new(mz, tims))
    }

    pub fn find_mz_model_for_frame(&self, frame: &SQLFrame) -> Result<MzCalibrationModel, MzCalibrationError> {
        match self.mz.iter().find(|m| m.id == frame.mz_calibration).map(|v| {
            MzCalibrationModel::try_from((v, frame.t1))
        }) {
            Some(value) => value,
            None => Err(MzCalibrationError::ModelNotFound(frame.mz_calibration)),
        }
    }

    pub fn find_tims_model_for_frame(&self, frame: &SQLFrame) -> Result<TimsCalibrationModel, IonMobilityCalibrationError> {
        match self.tims.iter().find(|m| m.id == frame.tims_calibration).map(|v| {
            TimsCalibrationModel::try_from(v)
        }) {
            Some(value) => value,
            None => Err(IonMobilityCalibrationError::ModelNotFound(frame.mz_calibration)),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct TimsCalibrationModel1 {
    pub c6: f64,
    pub c7: f64,
    pub offset: f64,
    pub slope: f64,
}

impl TimsCalibrationModel1 {
    pub fn new(c6: f64, c7: f64, offset: f64, slope: f64) -> Self {
        Self { c6, c7, offset, slope }
    }
}

#[derive(Debug, Error)]
pub enum IonMobilityCalibrationError {
    #[error("Ion mobility calibration model type {0} is not supported")]
    UnsupportedModel(u8),
    #[error("Missing model parameters: {0}")]
    MissingParameters(&'static str),
    #[error("Ion mobility model ID {0} not found")]
    ModelNotFound(u32),
}

impl TryFrom<&'_ TimsCalibration> for TimsCalibrationModel1 {
    type Error = IonMobilityCalibrationError;

    fn try_from(value: &'_ TimsCalibration) -> Result<Self, Self::Error> {
        if value.model_type != 2 {
            return Err(IonMobilityCalibrationError::UnsupportedModel(value.model_type));
        }
        let c0 = value.c0.ok_or( IonMobilityCalibrationError::MissingParameters("c0"))?;
        let c1 = value.c1.ok_or( IonMobilityCalibrationError::MissingParameters("c1"))?;
        let c2 = value.c2.ok_or( IonMobilityCalibrationError::MissingParameters("c2"))?;
        let c3 = value.c3.ok_or( IonMobilityCalibrationError::MissingParameters("c3"))?;
        let c4 = value.c4.ok_or( IonMobilityCalibrationError::MissingParameters("c4"))?;
        let c6 = value.c6.ok_or( IonMobilityCalibrationError::MissingParameters("c6"))?;
        let c7 = value.c7.ok_or( IonMobilityCalibrationError::MissingParameters("c7"))?;

        let slope = if c1 == 0.0 { 0.0 } else { (c3 - c2) / c1 };
        let offset = c2 - slope * (c4 + c0);
        Ok(Self::new(c6, c7, offset, slope))
    }
}

impl ConvertableDomain for TimsCalibrationModel1 {
    fn convert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        1.0 / (self.c6 + self.c7 / (self.offset + self.slope * value.into()))
    }

    fn invert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        // im = 1/(c6 + c7/(offset + slope*scan))
        // => scan = (c7/(1/im - c6) - offset) / slope
        let denom = (1.0 / value.into()) - self.c6;
        ((self.c7 / denom) - self.offset) / self.slope
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MzCalibrationModel1 {
    pub c0: f64,
    pub c1: f64,
    pub digitizer_timebase: f64,
    pub digitize_delay: f64,
}

impl MzCalibrationModel1 {
    pub fn new(c0: f64, c1: f64, digitizer_timebase: f64, digitize_delay: f64) -> Self {
        Self { c0, c1, digitizer_timebase, digitize_delay }
    }
}

#[derive(Debug, Error)]
pub enum MzCalibrationError {
    #[error("Mz calibration model type {0} is not supported")]
    UnsupportedModel(u8),
    #[error("Missing model parameters: {0}")]
    MissingParameters(&'static str),
    #[error("Mz model ID {0} not found")]
    ModelNotFound(u32),
    #[error("Mz calibration models are disabled")]
    Disabled,
}

impl TryFrom<(&'_ MzCalibration, f64)> for MzCalibrationModel1 {
    type Error = MzCalibrationError;

    fn try_from(value: (&'_ MzCalibration, f64)) -> Result<Self, Self::Error> {
        let (value, t1 ) = value;
        if value.model_type != 1 {
            return Err(MzCalibrationError::UnsupportedModel(value.model_type))
        }
        let c0 = value.c0.ok_or(MzCalibrationError::MissingParameters("c0"))?;
        let c1 = value.c1.ok_or(MzCalibrationError::MissingParameters("c1"))?;
        let cf = value.dc1 * (value.t1 - t1);
        let cf = 1.0 + (cf / 1.0e6);
        Ok(Self::new(c0, c1 * cf, value.digitizer_timebase, value.digitizer_delay))
    }
}

impl ConvertableDomain for MzCalibrationModel1 {
    fn convert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        let tof = (value.into() * self.digitizer_timebase) + self.digitize_delay;
        let inner = tof - self.c0;
        (self.c1 * inner.powi(2)) / 1e12
    }

    fn invert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        let mz = value.into();
        let tof = ((mz * 1e12) / self.c1).sqrt() + self.c0;
        (tof - self.digitize_delay) / self.digitizer_timebase
    }
}

pub fn clamp_u32(value: f64) -> u32 {
    const MAX_INDEX: f64 = (u32::MAX - 1) as f64;
    if value.is_nan() || value < 0.0 {
        0
    } else if value >= MAX_INDEX {
        u32::MAX - 1
    } else {
        value.round() as u32
    }
}


#[derive(Debug, Clone)]
pub enum TimsCalibrationModel {
    Basic(Scan2ImConverter),
    Model1(TimsCalibrationModel1),
}

impl ConvertableDomain for TimsCalibrationModel {
    fn convert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        match self {
            TimsCalibrationModel::Basic(scan2_im_converter) => scan2_im_converter.convert(value),
            TimsCalibrationModel::Model1(tims_calibration_model1) => tims_calibration_model1.convert(value),
        }
    }

    fn invert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        match self {
            TimsCalibrationModel::Basic(scan2_im_converter) => scan2_im_converter.invert(value),
            TimsCalibrationModel::Model1(tims_calibration_model1) => tims_calibration_model1.invert(value),
        }
    }
}

impl From<TimsCalibrationModel1> for TimsCalibrationModel {
    fn from(v: TimsCalibrationModel1) -> Self {
        Self::Model1(v)
    }
}

impl From<Scan2ImConverter> for TimsCalibrationModel {
    fn from(v: Scan2ImConverter) -> Self {
        Self::Basic(v)
    }
}

impl TryFrom<&'_ TimsCalibration> for TimsCalibrationModel {
    type Error = IonMobilityCalibrationError;

    fn try_from(value: &'_ TimsCalibration) -> Result<Self, Self::Error> {
        TimsCalibrationModel1::try_from(value).map(|v| TimsCalibrationModel::Model1(v))
    }
}


#[derive(Debug, Clone)]
pub enum MzCalibrationModel {
    Basic(Tof2MzConverter),
    Model1(MzCalibrationModel1)
}

impl TryFrom<(&'_ MzCalibration, f64)> for MzCalibrationModel {
    type Error = MzCalibrationError;

    fn try_from(value: (&'_ MzCalibration, f64)) -> Result<Self, Self::Error> {
        MzCalibrationModel1::try_from(value).map(|v| v.into())
    }
}

impl ConvertableDomain for MzCalibrationModel {
    fn convert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        match self {
            MzCalibrationModel::Basic(tof2_mz_converter) => tof2_mz_converter.convert(value),
            MzCalibrationModel::Model1(mz_calibration_model1) => mz_calibration_model1.convert(value),
        }
    }

    fn invert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        match self {
            MzCalibrationModel::Basic(tof2_mz_converter) => tof2_mz_converter.invert(value),
            MzCalibrationModel::Model1(mz_calibration_model1) => mz_calibration_model1.invert(value),
        }
    }
}

impl From<MzCalibrationModel1> for MzCalibrationModel {
    fn from(v: MzCalibrationModel1) -> Self {
        Self::Model1(v)
    }
}

impl From<Tof2MzConverter> for MzCalibrationModel {
    fn from(v: Tof2MzConverter) -> Self {
        Self::Basic(v)
    }
}


#[cfg(test)]
mod test_im {
    use super::*;

    #[test]
    fn scan2im_matches_fork_reference() {
        let cal = TimsCalibration {
            id: 1,
            model_type: 2,
            c0: Some(1.0),
            c1: Some(708.0),
            c2: Some(241.751905250524),
            c3: Some(99.2437539638487),
            c4: Some(33.9622641509434),
            c6: Some(0.0071422641733084),
            c7: Some(164.998795925213),
        };
        let conv = TimsCalibrationModel1::try_from(&cal).unwrap();
        const TOL: f64 = 5e-2;
        let im1 = conv.convert(1u32);
        assert!((im1 - 1.45).abs() < TOL, "im1={im1}");
        let im708 = conv.convert(708u32);
        assert!((im708 - 0.64).abs() < TOL, "im708={im708}");

        // round trip
        let back = conv.invert(im708) as u32;
        assert!((back as i64 - 708).abs() <= 1);
    }
}

#[cfg(test)]
mod test_mz {
    use super::*;

        #[test]
    fn tof2mz_matches_fork_reference() {
        let cal = MzCalibration {
            id: 1,
            model_type: 1,
            digitizer_timebase: 0.125,
            digitizer_delay: 25741.0,
            t1: 20.9410989491122,
            dc1: 20.0,
            c0: Some(286.065160463331),
            c1: Some(154317.348188993),
        };
        let real_t1 = 20.9455139021767;
        let conv = MzCalibrationModel1::try_from((&cal, real_t1)).unwrap();

        let mz0 = f64::from(conv.convert(0u32));
        let mz_max = f64::from(conv.convert(636029u32));
        const TOL: f64 = 1e-3;
        assert!((mz0 - 99.990834).abs() < TOL, "mz0={mz0}");
        assert!((mz_max - 1700.005).abs() < TOL, "mz_max={mz_max}");

        // round trip
        let back = conv.convert(mz_max);
        assert!(((back as u32) as i64 - 636029).abs() <= 1);
    }
}
