use std::{cmp::Ordering, fmt::Display, marker::PhantomData, str::FromStr};

use num_traits::AsPrimitive;
use serde::{Deserialize, Serialize};

use crate::{
    curie,
    io::usi::USI,
    params::{ControlledVocabulary, Param, ParamCow, Value, CURIE},
    prelude::*,
    spectrum::{
        ArrayType, BinaryArrayMap, BinaryDataArrayType, DataArray, IsolationWindowState,
        MultiLayerSpectrum, Precursor, ScanPolarity, SignalContinuity, SpectrumDescription,
    },
};

/// The possible PROXI backends
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub enum PROXIBackend {
    PeptideAtlas,
    MassIVE,
    Pride,
    Jpost,
    ProteomeXchange,
}

impl PROXIBackend {
    const ALL: &[Self] = &[
        Self::PeptideAtlas,
        Self::MassIVE,
        Self::Pride,
        Self::Jpost,
        Self::ProteomeXchange,
    ];

    /// The PROXI server base url which needs concatenating of the USI at the end
    const fn base_url(self) -> &'static str {
        match self {
            Self::PeptideAtlas => "http://www.peptideatlas.org/api/proxi/v0.1/spectra?resultType=full&usi=",
            Self::MassIVE => "http://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/spectra?resultType=full&usi=",
            Self::Pride => "http://www.ebi.ac.uk/pride/proxi/archive/v0.1/spectra?resultType=full&usi=",
            Self::Jpost => "https://repository.jpostdb.org/proxi/spectra?resultType=full&usi=",
            Self::ProteomeXchange => "http://proteomecentral.proteomexchange.org/api/proxi/v0.1/spectra?resultType=full&usi=",
        }
    }
}

impl USI {
    /// Retrieve this USI from the given PROXI backend. If no PROXI backend is indicated it will
    /// aggregate the results from all known backends and return the first successful spectrum.
    ///
    /// This function is only available with the feature `proxi`.
    pub fn get_spectrum_blocking(
        &self,
        backend: Option<PROXIBackend>,
    ) -> Result<(PROXIBackend, Vec<PROXISpectrum>), PROXIError> {
        backend.map_or_else(
            || {
                let client = reqwest::blocking::Client::new();
                let mut last_error = None;
                PROXIBackend::ALL
                    .iter()
                    .find_map(|backend| {
                        transform_response(
                            *backend,
                            client
                                .get(backend.base_url().to_string() + &self.to_string())
                                .send()
                                .and_then(reqwest::blocking::Response::json),
                        )
                        .map_err(|err| {
                            last_error = Some(err);
                        })
                        .ok()
                    })
                    .ok_or(last_error.unwrap_or(PROXIError::NotFound))
            },
            |backend| {
                transform_response(
                    backend,
                    reqwest::blocking::get(backend.base_url().to_string() + &self.to_string())
                        .and_then(reqwest::blocking::Response::json),
                )
            },
        )
    }

    /// Retrieve this USI from the given PROXI backend. If no PROXI backend is indicated it will
    /// aggregate the results from all known backends and return the first successful spectrum.
    ///
    /// This function is only available with the feature `proxi-async`.
    #[cfg(feature = "proxi-async")]
    pub async fn get_spectrum_async(
        &self,
        backend: Option<PROXIBackend>,
    ) -> Result<(PROXIBackend, Vec<PROXISpectrum>), PROXIError> {
        async fn get_response(
            client: &reqwest::Client,
            backend: PROXIBackend,
            usi: &str,
        ) -> Result<(PROXIBackend, Vec<PROXISpectrum>), PROXIError> {
            transform_response(
                backend,
                match client
                    .get(backend.base_url().to_string() + usi)
                    .send()
                    .await
                {
                    Ok(r) => r.json::<PROXIResponse>().await,
                    Err(e) => Err(e),
                },
            )
        }

        let client = reqwest::Client::new();
        if let Some(backend) = backend {
            get_response(&client, backend, &self.to_string()).await
        } else {
            use futures::StreamExt;

            let mut requests = futures::stream::FuturesUnordered::new();
            let mut last_error = None;
            for backend in PROXIBackend::ALL {
                requests.push(get_response(&client, *backend, &self.to_string()));
            }

            while let Some(res) = requests.next().await {
                match res {
                    Ok(s) => return Ok(s),
                    Err(e) => last_error = Some(e),
                }
            }

            Err(last_error.unwrap_or(PROXIError::NotFound))
        }
    }
}

fn transform_response(
    backend: PROXIBackend,
    response: Result<PROXIResponse, reqwest::Error>,
) -> Result<(PROXIBackend, Vec<PROXISpectrum>), PROXIError> {
    match response {
        Ok(PROXIResponse::Spectra(s)) => Ok((backend, s)),
        Ok(PROXIResponse::Error {
            detail,
            status,
            title,
            kind,
        }) => Err(PROXIError::Error {
            backend,
            detail,
            status,
            title,
            kind,
        }),
        Err(err) => Err(PROXIError::IO(backend, err)),
    }
}

#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum PROXIResponse {
    Spectra(Vec<PROXISpectrum>),
    Error {
        detail: String,
        status: usize,
        title: PROXIErrorType,
        #[serde(rename = "type")]
        kind: String,
    },
}

/// An error returned when accessing a PROXI server
#[derive(Debug)]
pub enum PROXIError {
    /// An error during the network request or decoding of the JSON response
    IO(PROXIBackend, reqwest::Error),
    /// A returned error by the server
    Error {
        /// Which backend failed
        backend: PROXIBackend,
        /// The type of error
        title: PROXIErrorType,
        /// Detailed explanation on the error
        detail: String,
        /// HTTP status code
        status: usize,
        /// The error kind, often "about:blank"
        kind: String,
    },
    /// An error when none of the aggregated backends returned a positive or negative result
    NotFound,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PROXIErrorType {
    /// The dataset is not present on this PROXI backend, it could be present on a different backend or the identifier does not exist
    #[serde(rename = "DatasetNotHere")]
    DataSetNotHere,
    /// The dataset is present, but the ms run is not, likely a mistake in the ms run name
    #[serde(rename = "MsRunNotFound")]
    MsRunNotFound,
    /// The dataset and the ms run are available, but the scan number does not exist in this file
    #[serde(rename = "ScanNotFound")]
    ScanNotFound,
    /// The dataset identifier is not of a recognisable format, commonly PXD identifiers are used eg 'PXD004939', but see the USI spec for more details
    #[serde(rename = "UnrecognizedIdentifierFormat")]
    UnrecognizedIdentifierFormat,
    /// The interpretation part of the USI is unable to be parsed, note that some PROXI backends require the addition of a charge to all peptides
    #[serde(rename = "MalformedInterpretation")]
    MalformedInterpretation,
    /// The index flag (scan/index/nativeid) is malformed
    #[serde(rename = "UnrecognizedIndexFlag")]
    UnrecognizedIndexFlag,
    /// Mandatory 'mzspec:' preamble is missing from the USI
    #[serde(rename = "MissingPreamble")]
    MissingPreamble,
    /// The USI is malformed and has too few fields
    #[serde(rename = "TooFewFields")]
    TooFewFields,
    #[serde(untagged)]
    Other(String),
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum Status {
    #[serde(rename = "READABLE")]
    Readable,
    #[serde(rename = "PEAK UNAVAILABLE")]
    PeakUnavailable,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PROXIValue(Value);

impl Default for PROXIValue {
    fn default() -> Self {
        Self(Value::Empty)
    }
}

impl PROXIValue {
    pub const fn is_empty(&self) -> bool {
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

        fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<Self::Value, E> {
            match v.parse::<PROXIValue>() {
                Ok(v) => Ok(v),
                Err(e) => Err(E::custom(e)),
            }
        }

        fn visit_bool<E: serde::de::Error>(self, v: bool) -> Result<Self::Value, E> {
            Ok(Value::Boolean(v).into())
        }

        fn visit_i64<E: serde::de::Error>(self, v: i64) -> Result<Self::Value, E> {
            Ok(Value::Int(v).into())
        }

        fn visit_u64<E: serde::de::Error>(self, v: u64) -> Result<Self::Value, E> {
            Ok(Value::Int(v as i64).into())
        }

        fn visit_f64<E: serde::de::Error>(self, v: f64) -> Result<Self::Value, E> {
            Ok(Value::Float(v).into())
        }

        fn visit_none<E: serde::de::Error>(self) -> Result<Self::Value, E> {
            Ok(Value::Empty.into())
        }

        fn visit_unit<E: serde::de::Error>(self) -> Result<Self::Value, E> {
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

        fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<Self::Value, E> {
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
        None => serializer.serialize_none(),
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

        fn visit_unit<E: serde::de::Error>(self) -> Result<Self::Value, E> {
            Ok(None)
        }

        fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<Self::Value, E> {
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

use serde::de::{self, Visitor};

/// MassIVE returns a list of strings instead of a list of numbers, this type can be deserialized if a number of string is given in the JSON
#[derive(Debug, Default, Clone)]
pub struct Wrapped<T>(T);

impl<T> std::ops::Deref for Wrapped<T> {
    type Target = T;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'de, T: 'static + Default + Copy + FromStr> serde::Deserialize<'de> for Wrapped<T>
where
    T::Err: Display,
    f64: AsPrimitive<T>,
    u64: AsPrimitive<T>,
    i64: AsPrimitive<T>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        deserializer
            .deserialize_any(PotentiallyWrappedNumberVisitor::<T>::default())
            .map(Wrapped)
    }
}

impl serde::Serialize for Wrapped<f64> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_f64(self.0)
    }
}

impl serde::Serialize for Wrapped<f32> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_f32(self.0)
    }
}

#[derive(Debug, Default)]
struct PotentiallyWrappedNumberVisitor<T> {
    marker: PhantomData<T>,
}

impl<'de, T: 'static + Copy + FromStr> Visitor<'de> for PotentiallyWrappedNumberVisitor<T>
where
    T::Err: Display,
    f64: AsPrimitive<T>,
    u64: AsPrimitive<T>,
    i64: AsPrimitive<T>,
{
    type Value = T;

    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        formatter.write_str("a wrapped number")
    }

    fn visit_f64<E: de::Error>(self, value: f64) -> Result<Self::Value, E> {
        Ok(value.as_())
    }

    fn visit_i64<E: de::Error>(self, value: i64) -> Result<Self::Value, E> {
        Ok(value.as_())
    }

    fn visit_u64<E: de::Error>(self, value: u64) -> Result<Self::Value, E> {
        Ok(value.as_())
    }

    fn visit_str<E: de::Error>(self, value: &str) -> Result<Self::Value, E> {
        value.parse().map_err(|e| serde::de::Error::custom(e))
    }
}

/// A spectrum returnd by a PROXI server
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct PROXISpectrum {
    /// The USI as returned by the PROXI server
    #[serde(
        default,
        serialize_with = "usi_serialize",
        deserialize_with = "usi_deserialize"
    )]
    pub usi: Option<USI>,
    /// The status of this request
    pub status: Option<Status>,
    /// Metadata for this spectrum
    #[serde(default)]
    pub attributes: Vec<PROXIParam>,
    #[serde(default)]
    pub mzs: Vec<Wrapped<f64>>,
    #[serde(default)]
    pub intensities: Vec<Wrapped<f32>>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub charges: Option<Vec<i32>>,
}

impl PROXISpectrum {
    pub fn add_attribute<P: Into<PROXIParam>>(&mut self, param: P) {
        self.attributes.push(param.into());
    }
}

impl From<&PROXISpectrum> for SpectrumDescription {
    fn from(value: &PROXISpectrum) -> Self {
        let mut this = Self::default();
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

        for param in &value.attributes {
            if matches!(
                param.accession.controlled_vocabulary,
                ControlledVocabulary::UO
            ) {
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
                }

                "scan start time" => {
                    if let Some(s) = this.acquisition.first_scan_mut() {
                        s.start_time = param
                            .value
                            .to_f64()
                            .expect("Failed to extract scan start time")
                            / 60.0;
                    }
                }
                "ion injection time" => {
                    if let Some(s) = this.acquisition.first_scan_mut() {
                        s.injection_time = param
                            .to_f32()
                            .expect("Failed to extract ion injection time");
                    }
                }
                "filter string" => {
                    if let Some(s) = this.acquisition.first_scan_mut() {
                        let mut param =
                            Param::new_key_value("filter string", param.value.to_string());
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
                        IsolationWindowState::Explicit | IsolationWindowState::Complete => {
                            IsolationWindowState::Complete
                        }
                        IsolationWindowState::Offset => {
                            precursor.isolation_window.lower_bound =
                                precursor.isolation_window.target
                                    - precursor.isolation_window.lower_bound;
                            precursor.isolation_window.upper_bound +=
                                precursor.isolation_window.target;
                            IsolationWindowState::Complete
                        }
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
                            precursor.isolation_window.lower_bound =
                                precursor.isolation_window.target - lower_bound;
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
                            precursor.isolation_window.upper_bound =
                                precursor.isolation_window.target + upper_bound;
                        }
                        _ => {}
                    }
                }
                "isolation window lower limit" => {
                    has_precursor = true;
                    let lower_bound = param
                        .to_f32()
                        .expect("Failed to parse isolation window limit");
                    if matches!(
                        precursor.isolation_window.flags,
                        IsolationWindowState::Unknown
                    ) {
                        precursor.isolation_window.flags = IsolationWindowState::Explicit;
                        precursor.isolation_window.lower_bound = lower_bound;
                    }
                }
                "isolation window upper limit" => {
                    has_precursor = true;
                    let upper_bound = param
                        .to_f32()
                        .expect("Failed to parse isolation window limit");
                    if matches!(
                        precursor.isolation_window.flags,
                        IsolationWindowState::Unknown
                    ) {
                        precursor.isolation_window.flags = IsolationWindowState::Explicit;
                        precursor.isolation_window.upper_bound = upper_bound;
                    }
                }
                _ => {
                    let mut p = Param::new_key_value(param.name.clone(), param.value.clone());
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

impl<
        C: CentroidLike + Default + BuildFromArrayMap + BuildArrayMapFrom,
        D: DeconvolutedCentroidLike + Default + BuildFromArrayMap + BuildArrayMapFrom,
    > From<PROXISpectrum> for MultiLayerSpectrum<C, D>
{
    fn from(value: PROXISpectrum) -> Self {
        let descr: SpectrumDescription = (&value).into();
        let mut arrays = BinaryArrayMap::default();

        let mut mz_array =
            DataArray::from_name_and_type(&ArrayType::MZArray, BinaryDataArrayType::Float64);
        mz_array
            .extend(&value.mzs.into_iter().map(|v| v.0).collect::<Vec<_>>())
            .unwrap();
        arrays.add(mz_array);

        let mut intensity_array =
            DataArray::from_name_and_type(&ArrayType::IntensityArray, BinaryDataArrayType::Float32);
        intensity_array
            .extend(
                &value
                    .intensities
                    .into_iter()
                    .map(|v| v.0)
                    .collect::<Vec<_>>(),
            )
            .unwrap();
        arrays.add(intensity_array);

        if let Some(charges) = value.charges.as_ref() {
            let mut charge_arrays =
                DataArray::from_name_and_type(&ArrayType::ChargeArray, BinaryDataArrayType::Int32);
            charge_arrays.extend(charges).unwrap();
            arrays.add(charge_arrays);
        };

        Self::from_arrays_and_description(arrays, descr)
    }
}

impl<T> From<&T> for PROXISpectrum
where
    T: SpectrumLike,
{
    fn from(value: &T) -> Self {
        let mut this = Self {
            status: Some(Status::Readable),
            ..Default::default()
        };

        match value.ms_level().cmp(&1) {
            Ordering::Equal => this.add_attribute(MS1_SPECTRUM.clone()),
            Ordering::Greater => this.add_attribute(MSN_SPECTRUM.clone()),
            Ordering::Less => (),
        }

        for param in value.params().iter().filter(|p| p.is_controlled()) {
            if param.curie().unwrap() == curie!(MS:1003063) {
                this.usi = Some(param.value.as_str().parse().unwrap());
            } else {
                this.add_attribute(param.clone());
            }
        }

        this.add_attribute(PROXIParam::new(
            curie!(MS:10003057),
            "scan number",
            Value::Int((value.index() + 1) as i64),
        ));

        this.add_attribute(PROXIParam {
            name: "ms level".to_string(),
            value: PROXIValue(Value::Int(value.ms_level() as i64)),
            accession: curie!(MS:1000511),
        });

        match value.polarity() {
            crate::spectrum::ScanPolarity::Unknown => {}
            crate::spectrum::ScanPolarity::Positive => {
                this.attributes.push(Param::from(POSITIVE_SCAN).into());
            }
            crate::spectrum::ScanPolarity::Negative => {
                this.attributes.push(Param::from(NEGATIVE_SCAN).into());
            }
        }

        match value.signal_continuity() {
            crate::spectrum::SignalContinuity::Unknown => {}
            crate::spectrum::SignalContinuity::Centroid => {
                this.attributes.push(Param::from(CENTROID_SPECTRUM).into());
            }
            crate::spectrum::SignalContinuity::Profile => {
                this.attributes.push(Param::from(PROFILE_SPECTRUM).into());
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
                this.mzs = arrays.mzs().unwrap().iter().copied().map(Wrapped).collect();
                this.intensities = arrays
                    .intensities()
                    .unwrap()
                    .iter()
                    .copied()
                    .map(Wrapped)
                    .collect();
                if let Ok(arr) = arrays.charges() {
                    this.charges = Some(arr.to_vec());
                }
            }
            crate::spectrum::RefPeakDataLevel::Centroid(peaks) => {
                (this.mzs, this.intensities) = peaks
                    .iter()
                    .map(|p| (Wrapped(p.mz()), Wrapped(p.intensity())))
                    .unzip();
            }
            crate::spectrum::RefPeakDataLevel::Deconvoluted(peaks) => {
                (this.mzs, this.intensities) = peaks
                    .iter()
                    .map(|p| (Wrapped(p.mz()), Wrapped(p.intensity())))
                    .unzip();
                this.charges = Some(
                    peaks
                        .iter()
                        .map(mzpeaks::KnownCharge::charge)
                        .collect::<Vec<_>>(),
                );
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
    use crate::MZReader;
    use serde_json;

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

    #[test]
    fn get_peptide_atlas() {
        let usi: USI =
            "mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555000000:VLHPLEGAVVIIFK/2"
                .parse()
                .unwrap();
        let (_, response) = usi
            .get_spectrum_blocking(Some(PROXIBackend::PeptideAtlas))
            .unwrap();
        dbg!(&response);
        assert!(!response.is_empty());
        todo!();
    }

    #[test]
    fn get_massive() {
        let usi: USI = "mzspec:MSV000078547:120228_nbut_3610_it_it_take2:scan:389"
            .parse()
            .unwrap();
        let (_, response) = usi
            .get_spectrum_blocking(Some(PROXIBackend::MassIVE))
            .unwrap();
        assert!(!response.is_empty());
    }

    #[test]
    fn get_pride() {
        let usi: USI =
            "mzspec:PXD043489:20201103_F1_UM5_Peng0013_SA_139H2_InS_Elastase.raw:scan:11809:VSLFPPSSEQLTSNASVV"
                .parse()
                .unwrap();
        let (_, response) = usi
            .get_spectrum_blocking(Some(PROXIBackend::Pride))
            .unwrap();
        assert!(!response.is_empty());
    }

    #[test]
    fn get_proteomexchange() {
        let usi: USI =
            "mzspec:PXD004939:Rice_phos_ABA_3h_20per_F1_R2:scan:2648:DAEKS[UNIMOD:21]PIN[UNIMOD:7]GR/2"
                .parse()
                .unwrap();
        let (_, response) = usi
            .get_spectrum_blocking(Some(PROXIBackend::ProteomeXchange))
            .unwrap();
        assert!(!response.is_empty());
    }

    #[test]
    fn get_aggregate() {
        for usi in [
            "mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2",
            "mzspec:MSV000078547:120228_nbut_3610_it_it_take2:scan:389",
            "mzspec:PXD043489:20201103_F1_UM5_Peng0013_SA_139H2_InS_Elastase.raw:scan:11809:VSLFPPSSEQLTSNASVV",
            "mzspec:PXD004939:Rice_phos_ABA_3h_20per_F1_R2:scan:2648:DAEKS[UNIMOD:21]PIN[UNIMOD:7]GR/2"] {
            println!("Trying: {usi}");
            let usi: USI = usi.parse().unwrap();
            let (_, response) = usi.get_spectrum_blocking(None).unwrap();
            assert!(!response.is_empty());
        }
    }
}

#[cfg(all(feature = "proxi-async", feature = "tokio"))]
mod tests {
    use super::*;

    #[tokio::test]
    async fn get_aggregate_async() {
        for usi in [
            "mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2",
            "mzspec:MSV000078547:120228_nbut_3610_it_it_take2:scan:389",
            "mzspec:PXD043489:20201103_F1_UM5_Peng0013_SA_139H2_InS_Elastase.raw:scan:11809:VSLFPPSSEQLTSNASVV",
            "mzspec:PXD004939:Rice_phos_ABA_3h_20per_F1_R2:scan:2648:DAEKS[UNIMOD:21]PIN[UNIMOD:7]GR/2"] {
            println!("Trying: {usi}");
            let usi: USI = usi.parse().unwrap();
            let (_, response) = usi.get_spectrum_async(None).await.unwrap();
            assert!(!response.is_empty());
        }
    }
}
