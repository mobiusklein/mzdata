use std::fmt::Display;
use std::str::FromStr;

use serde::{Deserialize, Serialize};

use super::usi::USI;
use crate::params::{ControlledVocabulary, Param, ParamCow, Value, CURIE};
use crate::spectrum::{ArrayType, BinaryArrayMap, BinaryDataArrayType, DataArray, IsolationWindowState, MultiLayerSpectrum, Precursor, ScanPolarity, SignalContinuity, SpectrumDescription};
use crate::{curie, prelude::*};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum Status {
    #[serde(rename = "READABLE")]
    Readable,
    #[serde(rename = "PEAK UNAVAILABLE")]
    PeakUnavailable,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PROXIValue(Value);

impl Default for PROXIValue {
    fn default() -> Self {
        Self(Value::Empty)
    }
}

impl PROXIValue {
    pub fn is_empty(&self) -> bool {
        matches!(self.0, Value::Empty)
    }
}

fn proxi_value_serialize<S>(proxi_value: &PROXIValue, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    match &proxi_value.0 {
        Value::String(v) => serializer.serialize_str(v),
        Value::Float(v) => serializer.serialize_f64(*v),
        Value::Int(v) => serializer.serialize_i64(*v),
        Value::Buffer(v) => serializer.serialize_bytes(v),
        Value::Boolean(v) => serializer.serialize_bool(*v),
        Value::Empty => serializer.serialize_unit(),
    }
}

fn proxi_value_deserialize<'de, D>(deserializer: D) -> Result<PROXIValue, D::Error>
where
    D: serde::Deserializer<'de>,
{
    struct PROXIValueVisit {}
    impl<'de> serde::de::Visitor<'de> for PROXIValueVisit {
        type Value = PROXIValue;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("PROXIValue string")
        }

        fn visit_str<E>(self, v: &str) -> Result<Self::Value, E>
        where
            E: serde::de::Error,
        {
            match v.parse::<PROXIValue>() {
                Ok(v) => Ok(v),
                Err(e) => Err(E::custom(e)),
            }
        }

        fn visit_bool<E>(self, v: bool) -> Result<Self::Value, E>
            where
                E: serde::de::Error, {
            Ok(Value::Boolean(v).into())
        }

        fn visit_i64<E>(self, v: i64) -> Result<Self::Value, E>
            where
                E: serde::de::Error, {
            Ok(Value::Int(v).into())
        }

        fn visit_u64<E>(self, v: u64) -> Result<Self::Value, E>
            where
                E: serde::de::Error, {
            Ok(Value::Int(v as i64).into())
        }

        fn visit_f64<E>(self, v: f64) -> Result<Self::Value, E>
            where
                E: serde::de::Error, {
            Ok(Value::Float(v).into())
        }

        fn visit_none<E>(self) -> Result<Self::Value, E>
            where
                E: serde::de::Error, {
            Ok(Value::Empty.into())
        }

        fn visit_unit<E>(self) -> Result<Self::Value, E>
            where
                E: serde::de::Error, {
            Ok(Value::Empty.into())
        }
    }

    deserializer.deserialize_any(PROXIValueVisit {})
}

impl ParamValue for PROXIValue {
    fn is_empty(&self) -> bool {
        <Value as ParamValue>::is_empty(&self.0)
    }

    fn is_i64(&self) -> bool {
        <Value as ParamValue>::is_i64(&self.0)
    }

    fn is_f64(&self) -> bool {
        <Value as ParamValue>::is_f64(&self.0)
    }

    fn is_buffer(&self) -> bool {
        <Value as ParamValue>::is_buffer(&self.0)
    }

    fn is_str(&self) -> bool {
        <Value as ParamValue>::is_str(&self.0)
    }

    fn to_f64(&self) -> Result<f64, crate::params::ParamValueParseError> {
        <Value as ParamValue>::to_f64(&self.0)
    }

    fn to_i64(&self) -> Result<i64, crate::params::ParamValueParseError> {
        <Value as ParamValue>::to_i64(&self.0)
    }

    fn to_str(&self) -> std::borrow::Cow<'_, str> {
        <Value as ParamValue>::to_str(&self.0)
    }

    fn to_buffer(&self) -> Result<std::borrow::Cow<'_, [u8]>, crate::params::ParamValueParseError> {
        <Value as ParamValue>::to_buffer(&self.0)
    }

    fn parse<T: FromStr>(&self) -> Result<T, T::Err> {
        <Value as ParamValue>::parse(&self.0)
    }

    fn as_bytes(&self) -> std::borrow::Cow<'_, [u8]> {
        <Value as ParamValue>::as_bytes(&self.0)
    }

    fn as_ref(&self) -> crate::params::ValueRef<'_> {
        <Value as ParamValue>::as_ref(&self.0)
    }

    fn data_len(&self) -> usize {
        <Value as ParamValue>::data_len(&self.0)
    }

    fn is_boolean(&self) -> bool {
        <Value as ParamValue>::is_boolean(&self.0)
    }

    fn to_bool(&self) -> Result<bool, crate::params::ParamValueParseError> {
        <Value as ParamValue>::to_bool(&self.0)
    }
}

impl Display for PROXIValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.0 {
            Value::String(v) => write!(f, "{v}"),
            Value::Float(v) => write!(f, "{v}"),
            Value::Int(v) => write!(f, "{v}"),
            Value::Buffer(v) => write!(f, "{v:?}"),
            Value::Boolean(v) => write!(f, "{v}"),
            Value::Empty => Ok(()),
        }
    }
}

impl FromStr for PROXIValue {
    type Err = <Value as FromStr>::Err;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let inner = s.parse()?;
        Ok(Self(inner))
    }
}

impl From<Value> for PROXIValue {
    fn from(value: Value) -> Self {
        Self(value)
    }
}

impl From<PROXIValue> for Value {
    fn from(value: PROXIValue) -> Self {
        value.0
    }
}

fn curie_serialize<S>(curie: &CURIE, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    serializer.serialize_str(&curie.to_string())
}

fn curie_deserialize<'de, D>(deserializer: D) -> Result<CURIE, D::Error>
where
    D: serde::Deserializer<'de>,
{
    struct CURIEVisit {}
    impl<'de> serde::de::Visitor<'de> for CURIEVisit {
        type Value = CURIE;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("expected CURIE string")
        }

        fn visit_str<E>(self, v: &str) -> Result<Self::Value, E>
        where
            E: serde::de::Error,
        {
            match v.parse::<CURIE>() {
                Ok(v) => Ok(v),
                Err(e) => Err(E::custom(e)),
            }
        }
    }

    deserializer.deserialize_str(CURIEVisit {})
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PROXIParam {
    #[serde(
        serialize_with = "curie_serialize",
        deserialize_with = "curie_deserialize"
    )]
    pub accession: CURIE,
    pub name: String,
    #[serde(
        serialize_with = "proxi_value_serialize",
        deserialize_with = "proxi_value_deserialize",
        skip_serializing_if = "PROXIValue::is_empty",
        default
    )]
    pub value: PROXIValue,
}

impl ParamValue for PROXIParam {
    fn is_empty(&self) -> bool {
        <PROXIValue as ParamValue>::is_empty(&self.value)
    }

    fn is_i64(&self) -> bool {
        <PROXIValue as ParamValue>::is_i64(&self.value)
    }

    fn is_f64(&self) -> bool {
        <PROXIValue as ParamValue>::is_f64(&self.value)
    }

    fn is_buffer(&self) -> bool {
        <PROXIValue as ParamValue>::is_buffer(&self.value)
    }

    fn is_str(&self) -> bool {
        <PROXIValue as ParamValue>::is_str(&self.value)
    }

    fn to_f64(&self) -> Result<f64, crate::params::ParamValueParseError> {
        <PROXIValue as ParamValue>::to_f64(&self.value)
    }

    fn to_i64(&self) -> Result<i64, crate::params::ParamValueParseError> {
        <PROXIValue as ParamValue>::to_i64(&self.value)
    }

    fn to_str(&self) -> std::borrow::Cow<'_, str> {
        <PROXIValue as ParamValue>::to_str(&self.value)
    }

    fn to_buffer(&self) -> Result<std::borrow::Cow<'_, [u8]>, crate::params::ParamValueParseError> {
        <PROXIValue as ParamValue>::to_buffer(&self.value)
    }

    fn parse<T: FromStr>(&self) -> Result<T, T::Err> {
        <PROXIValue as ParamValue>::parse(&self.value)
    }

    fn as_bytes(&self) -> std::borrow::Cow<'_, [u8]> {
        <PROXIValue as ParamValue>::as_bytes(&self.value)
    }

    fn as_ref(&self) -> crate::params::ValueRef<'_> {
        <PROXIValue as ParamValue>::as_ref(&self.value)
    }

    fn data_len(&self) -> usize {
        <PROXIValue as ParamValue>::data_len(&self.value)
    }

    fn is_boolean(&self) -> bool {
        <PROXIValue as ParamValue>::is_boolean(&self.value)
    }

    fn to_bool(&self) -> Result<bool, crate::params::ParamValueParseError> {
        <PROXIValue as ParamValue>::to_bool(&self.value)
    }
}

impl PROXIParam {
    pub fn new<S: ToString, V: Into<PROXIValue>>(accession: CURIE, name: S, value: V) -> Self {
        Self {
            accession,
            name: name.to_string(),
            value: value.into(),
        }
    }
}

impl From<Param> for PROXIParam {
    fn from(value: Param) -> Self {
        Self {
            accession: value.curie().unwrap(),
            name: value.name,
            value: value.value.into(),
        }
    }
}

impl<'a> From<ParamCow<'a>> for PROXIParam {
    fn from(value: ParamCow<'a>) -> Self {
        Self {
            accession: value.curie().unwrap(),
            name: value.name.to_string(),
            value: Value::from(value.value).into(),
        }
    }
}

const MS1_SPECTRUM: ParamCow = ControlledVocabulary::MS.const_param_ident("MS1 spectrum", 1000579);
const MSN_SPECTRUM: ParamCow = ControlledVocabulary::MS.const_param_ident("MSn spectrum", 1000580);
const NEGATIVE_SCAN: ParamCow =
    ControlledVocabulary::MS.const_param_ident("negative scan", 1000129);
const POSITIVE_SCAN: ParamCow =
    ControlledVocabulary::MS.const_param_ident("positive scan", 1000130);
const PROFILE_SPECTRUM: ParamCow =
    ControlledVocabulary::MS.const_param_ident("profile spectrum", 1000128);
const CENTROID_SPECTRUM: ParamCow =
    ControlledVocabulary::MS.const_param_ident("centroid spectrum", 1000127);

fn usi_serialize<S>(usi: &Option<USI>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    match usi {
        Some(usi) => serializer.serialize_str(&usi.to_string()),
        None => serializer.serialize_none()
    }

}

fn usi_deserialize<'de, D>(deserializer: D) -> Result<Option<USI>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    struct USIVisit {}
    impl<'de> serde::de::Visitor<'de> for USIVisit {
        type Value = Option<USI>;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("USI string or null")
        }

        fn visit_none<E>(self) -> Result<Self::Value, E> {
            Ok(None)
        }

        fn visit_unit<E>(self) -> Result<Self::Value, E>
            where
                E: serde::de::Error, {
            Ok(None)
        }

        fn visit_str<E>(self, v: &str) -> Result<Self::Value, E>
        where
            E: serde::de::Error,
        {
            if v == "null" {
                Ok(None)
            } else {
                match v.parse::<USI>() {
                    Ok(v) => Ok(Some(v)),
                    Err(e) => Err(E::custom(e)),
                }
            }
        }
    }

    deserializer.deserialize_any(USIVisit {})
}

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct PROXISpectrum {
    #[serde(serialize_with = "usi_serialize", deserialize_with = "usi_deserialize")]
    pub usi: Option<USI>,
    pub status: Option<Status>,
    pub attributes: Vec<PROXIParam>,
    #[serde(default)]
    pub mzs: Vec<f64>,
    #[serde(default)]
    pub intensities: Vec<f32>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub charges: Option<Vec<i32>>,
}

impl PROXISpectrum {
    pub fn add_attribute<P: Into<PROXIParam>>(&mut self, param: P) {
        self.attributes.push(param.into())
    }
}

impl From<&PROXISpectrum> for SpectrumDescription {
    fn from(value: &PROXISpectrum) -> Self {
        let mut this = SpectrumDescription::default();
        if let Some(usi) = value.usi.as_ref() {
            this.id = usi.to_string();
            if let Some(ident) = usi.identifier.as_ref() {
                match ident {
                    super::usi::Identifier::Scan(i) => this.index = (i.saturating_sub(1)) as usize,
                    super::usi::Identifier::Index(i) => this.index = *i as usize,
                    super::usi::Identifier::NativeID(_values) => todo!(),
                }
            } else {
                this.index = 0;
            }
        }

        let mut has_precursor = false;
        let mut precursor = Precursor::default();

        for param in value.attributes.iter() {
            if matches!(param.accession.controlled_vocabulary, ControlledVocabulary::UO) {
                continue;
            }
            match param.name.as_str() {
                "ms level" => {
                    this.ms_level = param.value.to_i32().expect("Failed to parse ms level") as u8;
                }
                "positive scan" => {
                    this.polarity = ScanPolarity::Positive;
                }
                "negative scan" => {
                    this.polarity = ScanPolarity::Negative;
                }
                "profile spectrum" => {
                    this.signal_continuity = SignalContinuity::Profile;
                }
                "centroid spectrum" => {
                    this.signal_continuity = SignalContinuity::Centroid;
                },

                "scan start time" => {
                    if let Some(s) = this.acquisition.first_scan_mut() {
                        s.start_time = param.value.to_f64().expect("Failed to extract scan start time") / 60.0;
                    }
                },
                "ion injection time" => {
                    if let Some(s) = this.acquisition.first_scan_mut() {
                        s.injection_time = param.to_f32().expect("Failed to extract ion injection time");
                    }
                }
                "filter string" => {
                    if let Some(s) = this.acquisition.first_scan_mut() {
                        let mut param = Param::new_key_value("filter string", param.value.to_string());
                        param.controlled_vocabulary = Some(ControlledVocabulary::MS);
                        param.accession = Some(1000512);
                        s.add_param(param);
                    }
                }

                "selected ion m/z" => {
                    has_precursor = true;
                    precursor.ion_mut().mz = param.to_f64().expect("Failed to parse ion m/z");
                }
                "peak intensity" => {
                    has_precursor = true;
                    precursor.ion_mut().intensity =
                        param.to_f32().expect("Failed to parse peak intensity");
                }
                "charge state" => {
                    has_precursor = true;
                    precursor.ion_mut().charge =
                        Some(param.to_i32().expect("Failed to parse ion charge"));
                }

                "isolation window target m/z" => {
                    has_precursor = true;
                    precursor.isolation_window.target = param
                        .to_f32()
                        .expect("Failed to parse isolation window target");
                    precursor.isolation_window.flags = match precursor.isolation_window.flags {
                        IsolationWindowState::Unknown => IsolationWindowState::Complete,
                        IsolationWindowState::Explicit => IsolationWindowState::Complete,
                        IsolationWindowState::Offset => {
                            precursor.isolation_window.lower_bound = precursor.isolation_window.target - precursor.isolation_window.lower_bound;
                            precursor.isolation_window.upper_bound += precursor.isolation_window.target;
                            IsolationWindowState::Complete
                        }
                        IsolationWindowState::Complete => IsolationWindowState::Complete,
                    };
                }
                "isolation window lower offset" => {
                    has_precursor = true;
                    let lower_bound = param
                        .to_f32()
                        .expect("Failed to parse isolation window limit");
                    match precursor.isolation_window.flags {
                        IsolationWindowState::Unknown => {
                            precursor.isolation_window.flags = IsolationWindowState::Offset;
                            precursor.isolation_window.lower_bound = lower_bound;
                        }
                        IsolationWindowState::Complete => {
                            precursor.isolation_window.lower_bound = precursor.isolation_window.target - lower_bound;
                        }
                        _ => {}
                    }
                }
                "isolation window upper offset" => {
                    has_precursor = true;
                    let upper_bound = param
                        .to_f32()
                        .expect("Failed to parse isolation window limit");
                    match precursor.isolation_window.flags {
                        IsolationWindowState::Unknown => {
                            precursor.isolation_window.flags = IsolationWindowState::Offset;
                            precursor.isolation_window.upper_bound = upper_bound;
                        }
                        IsolationWindowState::Complete => {
                            precursor.isolation_window.upper_bound = precursor.isolation_window.target + upper_bound;
                        }
                        _ => {}
                    }
                }
                "isolation window lower limit" => {
                    has_precursor = true;
                    let lower_bound = param
                        .to_f32()
                        .expect("Failed to parse isolation window limit");
                    if let IsolationWindowState::Unknown = precursor.isolation_window.flags {
                        precursor.isolation_window.flags = IsolationWindowState::Explicit;
                        precursor.isolation_window.lower_bound = lower_bound;
                    }
                }
                "isolation window upper limit" => {
                    has_precursor = true;
                    let upper_bound = param
                        .to_f32()
                        .expect("Failed to parse isolation window limit");
                    if let IsolationWindowState::Unknown = precursor.isolation_window.flags {
                        precursor.isolation_window.flags = IsolationWindowState::Explicit;
                        precursor.isolation_window.upper_bound = upper_bound;
                    }
                }
                _ => {
                    let mut p = Param::new_key_value(param.name.clone(), param.value.as_ref().to_owned());
                    p.accession = Some(param.accession.accession);
                    p.controlled_vocabulary = Some(param.accession.controlled_vocabulary);
                    this.add_param(p);
                }
            }
        }

        if has_precursor {
            this.precursor = Some(precursor);
        }
        this
    }
}

impl<C: CentroidLike + Default + BuildFromArrayMap + BuildArrayMapFrom, D: DeconvolutedCentroidLike + Default + BuildFromArrayMap + BuildArrayMapFrom> From<PROXISpectrum> for MultiLayerSpectrum<C, D> {
    fn from(value: PROXISpectrum) -> Self {
        let descr: SpectrumDescription = (&value).into();
        let mut arrays = BinaryArrayMap::default();

        let mut mz_array = DataArray::from_name_and_type(&ArrayType::MZArray, BinaryDataArrayType::Float64);
        mz_array.extend(&value.mzs).unwrap();
        arrays.add(mz_array);

        let mut intensity_array = DataArray::from_name_and_type(&ArrayType::IntensityArray, BinaryDataArrayType::Float32);
        intensity_array.extend(&value.intensities).unwrap();
        arrays.add(intensity_array);

        if let Some(charges) = value.charges.as_ref() {
            let mut charge_arrays = DataArray::from_name_and_type(&ArrayType::ChargeArray, BinaryDataArrayType::Int32);
            charge_arrays.extend(&charges).unwrap();
            arrays.add(charge_arrays);
        };

        MultiLayerSpectrum::from_arrays_and_description(arrays, descr)
    }
}

impl<T> From<&T> for PROXISpectrum
where
    T: SpectrumLike,
{
    fn from(value: &T) -> Self {
        let mut this = PROXISpectrum::default();
        this.status = Some(Status::Readable);

        let ms_level = value.ms_level();
        if ms_level == 1 {
            this.add_attribute(MS1_SPECTRUM.clone());
        } else if ms_level > 1 {
            this.add_attribute(MSN_SPECTRUM.clone());
        }

        for param in value.params().iter().filter(|p| p.is_controlled()) {
            if param.curie().unwrap() == curie!(MS:1003063) {
                this.usi = Some(param.value.as_str().parse().unwrap());
            } else {
                this.add_attribute(param.clone())
            }
        }

        this.add_attribute(PROXIParam::new(
            curie!(MS:10003057),
            "scan number",
            Value::Int((value.index() + 1) as i64),
        ));

        this.add_attribute(PROXIParam {
            name: "ms level".to_string(),
            value: PROXIValue(Value::Int(ms_level as i64)),
            accession: curie!(MS:1000511),
        });

        match value.polarity() {
            crate::spectrum::ScanPolarity::Unknown => {}
            crate::spectrum::ScanPolarity::Positive => {
                this.attributes.push(Param::from(POSITIVE_SCAN).into())
            }
            crate::spectrum::ScanPolarity::Negative => {
                this.attributes.push(Param::from(NEGATIVE_SCAN).into())
            }
        }

        match value.signal_continuity() {
            crate::spectrum::SignalContinuity::Unknown => {}
            crate::spectrum::SignalContinuity::Centroid => {
                this.attributes.push(Param::from(CENTROID_SPECTRUM).into())
            }
            crate::spectrum::SignalContinuity::Profile => {
                this.attributes.push(Param::from(PROFILE_SPECTRUM).into())
            }
        }

        if let Some(precursor) = value.precursor() {
            let iw = precursor.isolation_window();
            this.add_attribute(ControlledVocabulary::MS.param_val(
                "MS:1000827",
                "isolation window target m/z",
                iw.target,
            ));
            this.add_attribute(ControlledVocabulary::MS.param_val(
                "MS:1000828",
                "isolation window lower offset",
                iw.target - iw.lower_bound,
            ));
            this.add_attribute(ControlledVocabulary::MS.param_val(
                "MS:1000829",
                "isolation window upper offset",
                iw.upper_bound - iw.target,
            ));

            let ion = precursor.ion();
            this.add_attribute(PROXIParam::new(
                curie!(MS:1000744),
                "selected ion m/z",
                Value::Float(ion.mz),
            ));

            if let Some(z) = ion.charge {
                this.add_attribute(PROXIParam::new(
                    curie!(MS:1000041),
                    "charge state",
                    Value::Int(z as i64),
                ));
            }
        }

        for param in value.acquisition().params() {
            if param.is_controlled() {
                this.add_attribute(param.clone());
            }
        }

        for event in value.acquisition().iter() {
            let p = PROXIParam {
                name: "scan start time".into(),
                accession: curie!(MS:1000016),
                value: Value::Float(event.start_time * 60.0).into(),
            };
            this.add_attribute(p);
            for param in event.params() {
                if param.is_controlled() {
                    this.add_attribute(param.clone());
                }
            }
        }

        match value.peaks() {
            crate::spectrum::RefPeakDataLevel::Missing => {}
            crate::spectrum::RefPeakDataLevel::RawData(arrays) => {
                this.mzs = arrays.mzs().unwrap().to_vec();
                this.intensities = arrays.intensities().unwrap().to_vec();
                if let Ok(arr) = arrays.charges() {
                    this.charges = Some(arr.to_vec())
                }
            }
            crate::spectrum::RefPeakDataLevel::Centroid(peaks) => {
                (this.mzs, this.intensities) =
                    peaks.iter().map(|p| (p.mz(), p.intensity())).unzip();
            }
            crate::spectrum::RefPeakDataLevel::Deconvoluted(peaks) => {
                (this.mzs, this.intensities) =
                    peaks.iter().map(|p| (p.mz(), p.intensity())).unzip();
                this.charges = Some(peaks.iter().map(|p| p.charge()).collect::<Vec<_>>());
            }
        }

        if this.mzs.is_empty() {
            this.status = Some(Status::PeakUnavailable);
        }

        this
    }
}

#[cfg(test)]
mod test {
    use std::io;

    use super::*;
    use serde_json;
    use crate::MZReader;

    #[test]
    fn test_convert() -> io::Result<()> {
        let mut reader = MZReader::open_path("./test/data/batching_test.mzML")?;

        let scan = reader.get_spectrum_by_index(0).unwrap();
        let scan_message = PROXISpectrum::from(&scan);

        let message = serde_json::to_string(&scan_message)?;
         let dup: PROXISpectrum = serde_json::from_str(&message)?;
        assert_eq!(dup.usi, scan_message.usi);
        assert_eq!(dup.attributes, scan_message.attributes);
        Ok(())
    }
}