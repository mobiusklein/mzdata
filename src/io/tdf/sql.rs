use std::{collections::HashMap, path::Path, sync::Arc};

use parking_lot::{ReentrantMutex, ReentrantMutexGuard};
use rusqlite::{Connection, Error, Params, Row};

use crate::{params::Value, spectrum::ScanPolarity};

pub trait FromSQL: Sized {
    fn from_row(row: &Row<'_>) -> Result<Self, Error>;

    fn get_sql() -> String;

    fn read_from<I: Params>(connection: &Connection, params: I) -> Result<Vec<Self>, Error> {
        let sql = Self::get_sql();
        let mut stmt = connection.prepare(&sql)?;
        let out: Result<Vec<Self>, Error> = stmt
            .query_map(params, |row: &Row<'_>| Self::from_row(row))?
            .collect();
        out
    }

    fn read_from_where<I: Params>(
        connection: &Connection,
        params: I,
        condition_sql_fragment: &str,
    ) -> Result<Vec<Self>, Error> {
        let sql = Self::get_sql();
        let sql = format!("{sql} WHERE {condition_sql_fragment}");
        let mut stmt = connection.prepare(&sql)?;
        let out: Result<Vec<Self>, Error> = stmt
            .query_map(params, |row: &Row<'_>| Self::from_row(row))?
            .collect();
        out
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SQLPrecursor {
    pub id: usize,
    pub mz: f64,
    pub charge: i32,
    pub scan_average: f64,
    pub intensity: f64,
    pub precursor_frame: usize,
}

impl SQLPrecursor {
    #[allow(unused)]
    pub fn new(
        id: usize,
        mz: f64,
        charge: i32,
        scan_average: f64,
        intensity: f64,
        precursor_frame: usize,
    ) -> Self {
        Self {
            id,
            mz,
            charge,
            scan_average,
            intensity,
            precursor_frame,
        }
    }
}

impl FromSQL for SQLPrecursor {
    fn from_row(row: &Row<'_>) -> Result<Self, Error> {
        let this = Self {
            id: row.get(0).unwrap_or_default(),
            mz: row.get(1).unwrap_or_default(),
            charge: row.get(2).unwrap_or_default(),
            scan_average: row.get(3).unwrap_or_default(),
            intensity: row.get(4).unwrap_or_default(),
            precursor_frame: row.get(5).unwrap_or_default(),
        };

        Ok(this)
    }

    fn get_sql() -> String {
        "SELECT Id, MonoisotopicMz, Charge, ScanNumber, Intensity, Parent FROM Precursors".into()
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct SQLPasefFrameMsMs {
    pub frame: usize,
    pub scan_start: usize,
    pub scan_end: usize,
    pub isolation_mz: f64,
    pub isolation_width: f64,
    pub collision_energy: f64,
    pub precursor: usize,
}

impl FromSQL for SQLPasefFrameMsMs {
    fn get_sql() -> String {
        "SELECT Frame, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth, CollisionEnergy, Precursor FROM PasefFrameMsMsInfo".to_string()
    }

    fn from_row(row: &Row) -> Result<Self, Error> {
        let this = Self {
            frame: row.get(0).unwrap_or_default(),
            scan_start: row.get(1).unwrap_or_default(),
            scan_end: row.get(2).unwrap_or_default(),
            isolation_mz: row.get(3).unwrap_or_default(),
            isolation_width: row.get(4).unwrap_or_default(),
            collision_energy: row.get(5).unwrap_or_default(),
            precursor: row.get(6).unwrap_or_default(),
        };

        Ok(this)
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct SQLDIAFrameMsMsWindow {
    pub window_group: u32,
    pub scan_start: usize,
    pub scan_end: usize,
    pub isolation_mz: f64,
    pub isolation_width: f64,
    pub collision_energy: f32,
}

impl SQLDIAFrameMsMsWindow {
    pub fn new(
        window_group: u32,
        scan_number_begin: usize,
        scan_number_end: usize,
        isolation_mz: f64,
        isolation_width: f64,
        collision_energy: f32,
    ) -> Self {
        Self {
            window_group,
            scan_start: scan_number_begin,
            scan_end: scan_number_end,
            isolation_mz,
            isolation_width,
            collision_energy,
        }
    }
}

impl FromSQL for SQLDIAFrameMsMsWindow {
    fn from_row(row: &Row<'_>) -> Result<Self, Error> {
        let this = SQLDIAFrameMsMsWindow::new(
            row.get(0)?,
            row.get(1)?,
            row.get(2)?,
            row.get(3)?,
            row.get(4)?,
            row.get(5)?,
        );
        Ok(this)
    }

    fn get_sql() -> String {
        "SELECT WindowGroup, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth, CollisionEnergy FROM DiaFrameMsMsWindows".into()
    }
}

pub struct KeyValue(String, Value);

impl FromSQL for KeyValue {
    fn from_row(row: &Row<'_>) -> Result<Self, Error> {
        let key = row.get(0)?;
        let val = row.get(1)?;
        Ok(KeyValue(key, Value::new(val)))
    }

    fn get_sql() -> String {
        "SELECT Key, Value FROM GlobalMetadata".into()
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct SQLFrame {
    pub id: usize,
    pub time: f64,
    pub polarity: ScanPolarity,
    pub scan_mode: u8,
    pub msms_type: u8,
    pub tims_id: usize,
    pub max_intensity: f32,
    pub summed_intensities: f32,
    pub accumulation_time: f32,
    pub ramp_time: f32,
    pub property_group: usize,
    pub num_scans: usize,
    pub num_peaks: usize,
}

impl SQLFrame {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        id: usize,
        time: f64,
        polarity: ScanPolarity,
        scan_mode: u8,
        msms_type: u8,
        tims_id: usize,
        max_intensity: f32,
        summed_intensities: f32,
        accumulation_time: f32,
        ramp_time: f32,
        property_group: usize,
        num_scans: usize,
        num_peaks: usize,
    ) -> Self {
        Self {
            id,
            time,
            polarity,
            scan_mode,
            msms_type,
            tims_id,
            max_intensity,
            summed_intensities,
            accumulation_time,
            ramp_time,
            property_group,
            num_scans,
            num_peaks,
        }
    }
}

impl FromSQL for SQLFrame {
    fn from_row(row: &Row<'_>) -> Result<Self, Error> {
        let this = Self::new(
            row.get(0)?,
            row.get(1)?,
            row.get::<usize, String>(2).map(|c| match c.as_str() {
                "+" => ScanPolarity::Positive,
                "-" => ScanPolarity::Negative,
                _ => ScanPolarity::Unknown,
            })?,
            row.get(3)?,
            row.get(4)?,
            row.get(5)?,
            row.get(6)?,
            row.get(7)?,
            row.get(8)?,
            row.get(9)?,
            row.get(10)?,
            row.get(11)?,
            row.get(12)?,
        );
        Ok(this)
    }

    fn get_sql() -> String {
        "SELECT Id, Time, Polarity, ScanMode, MsMsType, TimsId, MaxIntensity, SummedIntensities, AccumulationTime, RampTime, PropertyGroup, NumScans, NumPeaks FROM Frames".into()
    }
}

#[derive(Debug)]
pub struct RawTDFSQLReader {
    pub connection: ReentrantMutex<Connection>,
}

#[allow(unused)]
impl RawTDFSQLReader {
    pub fn new(tdf_path: &Path) -> Result<Self, Error> {
        let connection = ReentrantMutex::new(Connection::open(tdf_path)?);
        Ok(Self { connection })
    }

    pub fn connection(&self) -> ReentrantMutexGuard<'_, Connection> {
        self.connection.try_lock().unwrap()
    }

    pub fn metadata(&self) -> Result<HashMap<String, Value>, Error> {
        Ok(KeyValue::read_from(&self.connection(), [])?
            .into_iter()
            .map(|KeyValue(k, v)| (k, v))
            .collect())
    }

    pub fn query<T: FromSQL>(&self, sql: &str, params: impl Params) -> Result<Vec<T>, Error> {
        let conn = self.connection();
        let mut stmt = conn.prepare(sql)?;
        let out: Result<Vec<T>, Error> = stmt
            .query_map(params, |row: &Row<'_>| T::from_row(row))?
            .collect();
        out
    }

    fn has_table(&self, table: &str) -> Result<bool, Error> {
        self.connection().query_row(
            "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
            [table],
            |_| Ok(true),
        )
    }

    pub fn is_dda_pasef(&self) -> Result<bool, Error> {
        if self.has_table("PasefFrameMsMsInfo")? {
            self.connection()
                .query_row("SELECT COUNT(*) FROM PasefFrameMsMsInfo", [], |row| {
                    let val: u64 = row.get(0)?;
                    Ok(val > 0)
                })
        } else {
            Ok(false)
        }
    }

    pub fn is_dia_pasef(&self) -> Result<bool, Error> {
        if self.has_table("DiaFrameMsMsInfo")? {
            self.connection()
                .query_row("SELECT COUNT(*) FROM DiaFrameMsMsInfo", [], |row| {
                    let val: u64 = row.get(0)?;
                    Ok(val > 0)
                })
        } else {
            Ok(false)
        }
    }

    pub fn is_prm_pasef(&self) -> Result<bool, Error> {
        if self.has_table("PrmFrameMsMsInfo")? {
            self.connection()
                .query_row("SELECT COUNT(*) FROM PrmFrameMsMsInfo", [], |row| {
                    let val: u64 = row.get(0)?;
                    Ok(val > 0)
                })
        } else {
            Ok(false)
        }
    }

    pub fn precursors_for(&self, frame_id: usize) -> Result<Vec<SQLPrecursor>, Error> {
        SQLPrecursor::read_from_where(&self.connection(), [frame_id], "Parent = ?")
    }

    pub fn pasef_frames_for(
        &self,
        frame_id: Option<usize>,
        precursor_id: Option<usize>,
    ) -> Result<Vec<SQLPasefFrameMsMs>, Error> {
        match (frame_id, precursor_id) {
            (Some(frame_id), Some(precursor_id)) => SQLPasefFrameMsMs::read_from_where(
                &self.connection(),
                [frame_id, precursor_id],
                "Frame = ? AND Precursor = ?",
            ),
            (Some(frame_id), None) => {
                SQLPasefFrameMsMs::read_from_where(&self.connection(), [frame_id], "Frame = ?")
            }
            (None, Some(precursor_id)) => SQLPasefFrameMsMs::read_from_where(
                &self.connection(),
                [precursor_id],
                "Precursor = ?",
            ),
            (None, None) => SQLPasefFrameMsMs::read_from(&self.connection(), []),
        }
    }

    pub fn dia_window_group_for(&self, frame_id: usize) -> Result<u32, Error> {
        let conn = self.connection();
        let mut stmt = conn.prepare("SELECT WindowGroup FROM DiaFrameMsMsInfo WHERE Frame = ?")?;
        let window_group: u32 = stmt.query_row([frame_id], |row| row.get(0))?;
        Ok(window_group)
    }
}

#[derive(Debug, Clone)]
pub struct PasefPrecursor {
    pub precursor: Arc<SQLPrecursor>,
    pub pasef_msms: Arc<SQLPasefFrameMsMs>,
}

impl PasefPrecursor {
    pub fn new(precursor: Arc<SQLPrecursor>, pasef_msms: Arc<SQLPasefFrameMsMs>) -> Self {
        Self {
            precursor,
            pasef_msms,
        }
    }
}

#[derive(Debug, Default, Clone)]
pub enum TDFMSnFacet {
    PasefDDA(PasefPrecursor),
    PasefDIA(Arc<SQLDIAFrameMsMsWindow>),
    #[allow(unused)]
    CompoundPasefDDA(Box<[PasefPrecursor]>),
    #[allow(unused)]
    CompoundPasefDIA(Arc<[Arc<SQLDIAFrameMsMsWindow>]>),
    #[default]
    None,
}

impl From<Arc<SQLDIAFrameMsMsWindow>> for TDFMSnFacet {
    fn from(value: Arc<SQLDIAFrameMsMsWindow>) -> Self {
        Self::PasefDIA(value)
    }
}

impl From<PasefPrecursor> for TDFMSnFacet {
    fn from(value: PasefPrecursor) -> Self {
        Self::PasefDDA(value)
    }
}

#[allow(unused)]
pub type PasefPrecursorIter<'a> = std::iter::Chain<
    std::slice::Iter<'a, PasefPrecursor>,
    std::option::IntoIter<&'a PasefPrecursor>,
>;

#[allow(unused)]
pub type DIAWindowIter<'a> = std::iter::Chain<
    std::slice::Iter<'a, Arc<SQLDIAFrameMsMsWindow>>,
    std::option::IntoIter<&'a Arc<SQLDIAFrameMsMsWindow>>,
>;

impl TDFMSnFacet {
    #[allow(unused)]
    pub fn scan_range(&self) -> Option<(usize, usize)> {
        match self {
            Self::PasefDDA(pasef_precursor) => {
                let pasef = &pasef_precursor.pasef_msms;
                Some((pasef.scan_start, pasef.scan_end))
            }
            Self::PasefDIA(dia) => Some((dia.scan_start, dia.scan_end)),
            Self::None => None,
            Self::CompoundPasefDDA(_) => None,
            Self::CompoundPasefDIA(_) => None,
        }
    }

    #[allow(unused)]
    pub fn is_compound(&self) -> bool {
        matches!(self, Self::CompoundPasefDDA(_) | Self::CompoundPasefDIA(_))
    }

    pub fn precursor(&self) -> Option<&SQLPrecursor> {
        if let Self::PasefDDA(prec) = self {
            Some(&prec.precursor)
        } else {
            None
        }
    }

    pub fn pasef_msms(&self) -> Option<&SQLPasefFrameMsMs> {
        if let Self::PasefDDA(prec) = self {
            Some(&prec.pasef_msms)
        } else {
            None
        }
    }

    pub fn dia_window(&self) -> Option<&SQLDIAFrameMsMsWindow> {
        if let Self::PasefDIA(dia) = self {
            Some(dia)
        } else {
            None
        }
    }

    #[allow(unused)]
    pub fn iter_pasef_msms(&self) -> PasefPrecursorIter<'_> {
        let c = if let Self::CompoundPasefDDA(dda) = self {
            dda.iter()
        } else {
            [].iter()
        };

        if let Self::PasefDDA(p) = self {
            c.chain(Some(p))
        } else {
            c.chain(None)
        }
    }

    #[allow(unused)]
    pub fn iter_dia_windows(&self) -> DIAWindowIter<'_> {
        let c = if let Self::CompoundPasefDIA(dda) = self {
            dda.iter()
        } else {
            [].iter()
        };

        if let Self::PasefDIA(p) = self {
            c.chain(Some(p))
        } else {
            c.chain(None)
        }
    }
}
