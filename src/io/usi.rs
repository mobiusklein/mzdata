use std::{fmt::Display, str::FromStr};

use serde::Deserialize;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum USIParseError {
    #[error("Protocol {0} is not recognized")]
    UnknownProtocol(String, String),
    #[error("Did not find a dataset section in {0}")]
    MissingDataset(String),
    #[error("Did not find a run name in {0}")]
    MissingRun(String),
    #[error("Index type {0} is not recognized")]
    UnknownIndexType(String, String),
    #[error("Malformed spectrum index {0}: {1}")]
    MalformedIndex(String, String, String),
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Protocol {
    #[default]
    MZSpec,
}

impl Display for Protocol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MZSpec => write!(f, "mzspec"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Identifier {
    Scan(u64),
    Index(u64),
    NativeID(Box<Vec<u64>>),
}

#[derive(Debug, Default, Clone, PartialEq, Eq, Hash)]
pub struct USI {
    pub protocol: Protocol,
    pub dataset: String,
    pub run_name: String,
    pub identifier: Option<Identifier>,
    pub interpretation: Option<String>,
    pub provenance: Option<String>,
}

pub enum PROXIBackend {
    PeptideAtlas,
    MassIVE,
    Pride,
    Jpost,
    ProteomeXchange,
}

impl PROXIBackend {
    fn base_url(&self) -> &'static str {
        match self {
            Self::PeptideAtlas => "http://www.peptideatlas.org/api/proxi/v0.1/spectra?resultType=full&usi=",
            Self::MassIVE => "http://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/spectra?resultType=full&usi=",
            Self::Pride => "http://www.ebi.ac.uk/pride/proxi/archive/v0.1/spectra?resultType=full&usi=",
            Self::Jpost => "https://repository.jpostdb.org/proxi/spectra?resultType=full&usi=",
            Self::ProteomeXchange => "http://proteomecentral.proteomexchange.org/api/proxi/v0.1/spectra?resultType=full&usi=",
        }
    }
}
use crate::{params::ControlledVocabulary, CentroidSpectrum, Param};

impl USI {
    #[cfg(feature = "reqwest")]
    pub fn get_spectrum(
        &self,
        backend: Option<PROXIBackend>,
    ) -> Result<CentroidSpectrum, USIRetrieveError> {
        use mzpeaks::CentroidPeak;

        use crate::spectrum::{Acquisition, SpectrumDescription};

        let response = if let Some(backend) = backend {
            reqwest::blocking::get(backend.base_url().to_string() + &self.to_string())?
                .json::<PROXIResponse>()
                .map_err(USIRetrieveError::IO)
        } else {
            let client = reqwest::blocking::Client::new();
            [
                PROXIBackend::PeptideAtlas,
                PROXIBackend::MassIVE,
                PROXIBackend::Pride,
                PROXIBackend::Jpost,
                PROXIBackend::ProteomeXchange,
            ]
            .iter()
            .find_map(|backend| {
                client
                    .get(backend.base_url().to_string() + &self.to_string())
                    .send()
                    .and_then(|r| r.json::<PROXIResponse>())
                    .ok()
            })
            .ok_or(USIRetrieveError::NotFound)
        }?
        .into_iter()
        .find(|r| r.status.as_ref().is_none_or(|s| s == "READABLE"))
        .ok_or(USIRetrieveError::NotFound)?;

        dbg!(&response);

        let description = SpectrumDescription::new(
            if let Some(Identifier::NativeID(id)) = &self.identifier {
                id.iter()
                    .map(|i| i.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            } else {
                String::new()
            },
            if let Some(Identifier::Index(i)) = &self.identifier {
                *i as usize
            } else {
                response
                    .attributes
                    .iter()
                    .find(|a| a.name == "scan number")
                    .and_then(|a| a.value.as_ref())
                    .and_then(|v| v.parse::<usize>().ok())
                    .unwrap_or(0)
            },
            0,
            crate::spectrum::ScanPolarity::Unknown,
            crate::spectrum::SignalContinuity::Unknown,
            response.attributes.iter().map(|a| a.into()).collect(),
            Acquisition::default(),
            None,
        ); // TODO: use more of the attributes

        let peaks = response
            .mzs
            .into_iter()
            .zip(response.intensities)
            .enumerate()
            .map(|(index, (mz, intensity))| {
                CentroidPeak::new(mz.0, intensity.0 as f32, index as u32)
            })
            .collect();

        Ok(CentroidSpectrum::new(description, peaks))
    }
}

#[derive(Debug)]
pub enum USIRetrieveError {
    IO(reqwest::Error),
    NotFound,
    // TODO: handle spectrum not found cases (will be in resulting json)
}

impl From<reqwest::Error> for USIRetrieveError {
    fn from(value: reqwest::Error) -> Self {
        Self::IO(value)
    }
}

type PROXIResponse = Vec<Spectrum>;

#[derive(Deserialize, Debug)]
struct Spectrum {
    #[serde(default)]
    attributes: Vec<Attribute>,
    #[serde(default)]
    mzs: Vec<PotentiallyWrappedNumber>,
    #[serde(default)]
    intensities: Vec<PotentiallyWrappedNumber>,
    usi: Option<String>,
    status: Option<String>,
}

use serde::de::{self, Visitor};

/// MassIVE returns a list of strings instead of a list of numbers, this type can be deserialized if a number of string is given in the JSON
#[derive(Debug)]
struct PotentiallyWrappedNumber(f64);

impl<'de> serde::Deserialize<'de> for PotentiallyWrappedNumber {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        deserializer
            .deserialize_any(PotentiallyWrappedNumberVisitor)
            .map(PotentiallyWrappedNumber)
    }
}

#[derive(Debug)]
struct PotentiallyWrappedNumberVisitor;

impl<'de> Visitor<'de> for PotentiallyWrappedNumberVisitor {
    type Value = f64;

    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        formatter.write_str("a f64")
    }

    fn visit_f64<E: de::Error>(self, value: f64) -> Result<Self::Value, E> {
        Ok(value)
    }

    fn visit_i64<E: de::Error>(self, value: i64) -> Result<Self::Value, E> {
        Ok(value as f64)
    }

    fn visit_u64<E: de::Error>(self, value: u64) -> Result<Self::Value, E> {
        Ok(value as f64)
    }

    fn visit_str<E: de::Error>(self, value: &str) -> Result<Self::Value, E> {
        value.parse().map_err(|e| serde::de::Error::custom(e))
    }
}

#[derive(Deserialize, Debug)]
struct Attribute {
    accession: String,
    name: String,
    value: Option<String>,
}

impl From<&Attribute> for Param {
    fn from(value: &Attribute) -> Self {
        let (controlled_vocabulary, accession) = if let Some((Ok(cv), Ok(accession))) = value
            .accession
            .split_once(':')
            .map(|(cv, accession)| (ControlledVocabulary::from_str(cv), accession.parse::<u32>()))
        {
            (Some(cv), Some(accession))
        } else {
            (None, None)
        };
        Param {
            name: value.name.to_string(),
            value: value
                .value
                .as_ref()
                .map(|v| crate::params::Value::new(v.to_string()))
                .unwrap_or(crate::params::Value::Empty),
            accession,
            controlled_vocabulary,
            unit: crate::params::Unit::Unknown,
        }
    }
}

impl FromStr for USI {
    type Err = USIParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut this = Self::default();
        let mut tokens = s.split(':');
        if let Some(protocol) = tokens.next() {
            match protocol {
                "mzspec" => {
                    this.protocol = Protocol::MZSpec;
                }
                _ => {
                    return Err(USIParseError::UnknownProtocol(
                        protocol.to_string(),
                        s.to_string(),
                    ))
                }
            };

            if let Some(dataset) = tokens.next() {
                this.dataset = dataset.to_string();

                if let Some(run) = tokens.next() {
                    this.run_name = run.to_string();

                    if let (Some(ident_type), Some(ident_value)) = (tokens.next(), tokens.next()) {
                        match ident_type {
                            "scan" => match ident_value.parse() {
                                Ok(v) => {
                                    this.identifier = Some(Identifier::Scan(v));
                                }
                                Err(e) => {
                                    return Err(USIParseError::MalformedIndex(
                                        ident_value.to_string(),
                                        e.to_string(),
                                        s.to_string(),
                                    ))
                                }
                            },
                            "index" => match ident_value.parse() {
                                Ok(v) => {
                                    this.identifier = Some(Identifier::Index(v));
                                }
                                Err(e) => {
                                    return Err(USIParseError::MalformedIndex(
                                        ident_value.to_string(),
                                        e.to_string(),
                                        s.to_string(),
                                    ))
                                }
                            },
                            "nativeId" => {
                                let res: Result<Vec<u64>, _> =
                                    ident_value.split(',').map(|t| t.parse()).collect();
                                match res {
                                    Ok(vals) => {
                                        this.identifier = Some(Identifier::NativeID(vals.into()))
                                    }
                                    Err(e) => {
                                        return Err(USIParseError::MalformedIndex(
                                            ident_value.to_string(),
                                            e.to_string(),
                                            s.to_string(),
                                        ))
                                    }
                                }
                            }
                            _ => {
                                return Err(USIParseError::UnknownIndexType(
                                    ident_type.to_string(),
                                    s.to_string(),
                                ))
                            }
                        };

                        this.interpretation = tokens.next().map(|s| s.to_string());
                        this.provenance = tokens.next_back().map(|s| s.to_string());
                    }
                    Ok(this)
                } else {
                    Err(USIParseError::MissingRun(s.to_string()))
                }
            } else {
                Err(USIParseError::MissingDataset(s.to_string()))
            }
        } else {
            Err(USIParseError::UnknownProtocol(
                "".to_string(),
                s.to_string(),
            ))
        }
    }
}

impl Display for USI {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(ident) = self.identifier.as_ref() {
            let (ident_class, ident_val) = match ident {
                Identifier::Scan(i) => ("scan", i.to_string()),
                Identifier::Index(i) => ("index", i.to_string()),
                Identifier::NativeID(parts) => (
                    "nativeId",
                    parts
                        .iter()
                        .map(|i| i.to_string())
                        .collect::<Vec<_>>()
                        .join(","),
                ),
            };
            write!(
                f,
                "{}:{}:{}:{ident_class}:{ident_val}",
                self.protocol, self.dataset, self.run_name
            )?;
            if let Some(interp) = self.interpretation.as_ref() {
                write!(f, ":{}", interp)?;
                if let Some(provenance) = self.provenance.as_ref() {
                    write!(f, ":{}", provenance)?;
                }
            }
            Ok(())
        } else {
            write!(f, "{}:{}:{}", self.protocol, self.dataset, self.run_name)
        }
    }
}

#[cfg(test)]
mod test {
    use mzpeaks::PeakCollection;

    use super::*;

    #[test]
    fn test_usi_parse() -> Result<(), USIParseError> {
        let usi: USI = "mzspec:PXD000001:foo".parse()?;

        assert_eq!(usi.protocol, Protocol::MZSpec);
        assert_eq!(usi.dataset, "PXD000001");
        assert_eq!(usi.run_name, "foo");

        let usi: USI = "mzspec:PXD000001:foo:index:100".parse()?;
        assert_eq!(usi.identifier, Some(Identifier::Index(100)));

        let usi: USI = "mzspec:PXD000001:foo:scan:100".parse()?;
        assert_eq!(usi.identifier, Some(Identifier::Scan(100)));
        Ok(())
    }

    #[test]
    fn test_example() -> Result<(), USIParseError> {
        let usi: USI = "mzspec:PXD019909:20180914_QE8_nLC0_BDA_SA_DIA_Skin_Dendritic_cells_DC_MT_600000:scan:62396:SAGQGEVLVYVEDPAGHQEEAK/3".parse()?;
        assert_eq!(usi.dataset, "PXD019909");
        assert_eq!(
            usi.run_name,
            "20180914_QE8_nLC0_BDA_SA_DIA_Skin_Dendritic_cells_DC_MT_600000"
        );
        assert_eq!(usi.identifier, Some(Identifier::Scan(62396)));
        assert_eq!(
            usi.interpretation,
            Some("SAGQGEVLVYVEDPAGHQEEAK/3".to_string())
        );
        Ok(())
    }

    #[test]
    fn get_peptide_atlas() {
        let usi: USI =
            "mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2"
                .parse()
                .unwrap();
        let response = usi.get_spectrum(Some(PROXIBackend::PeptideAtlas)).unwrap();
        assert!(!response.is_empty())
    }

    #[test]
    fn get_massive() {
        let usi: USI = "mzspec:MSV000078547:120228_nbut_3610_it_it_take2:scan:389"
            .parse()
            .unwrap();
        let response = usi.get_spectrum(Some(PROXIBackend::MassIVE)).unwrap();
        assert!(!response.is_empty())
    }

    #[test]
    fn get_pride() {
        let usi: USI =
            "mzspec:PXD043489:20201103_F1_UM5_Peng0013_SA_139H2_InS_Elastase.raw:scan:11809:VSLFPPSSEQLTSNASVV"
                .parse()
                .unwrap();
        let response = usi.get_spectrum(Some(PROXIBackend::Pride)).unwrap();
        assert!(!response.is_empty())
    }

    #[test]
    fn get_proteomexchange() {
        let usi: USI =
            "mzspec:PXD004939:Rice_phos_ABA_3h_20per_F1_R2:scan:2648:DAEKS[UNIMOD:21]PIN[UNIMOD:7]GR/2"
                .parse()
                .unwrap();
        let response = usi
            .get_spectrum(Some(PROXIBackend::ProteomeXchange))
            .unwrap();
        assert!(!response.is_empty())
    }

    #[test]
    fn get_aggregate() {
        for usi in ["mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2", "mzspec:MSV000078547:120228_nbut_3610_it_it_take2:scan:389","mzspec:PXD043489:20201103_F1_UM5_Peng0013_SA_139H2_InS_Elastase.raw:scan:11809:VSLFPPSSEQLTSNASVV","mzspec:PXD004939:Rice_phos_ABA_3h_20per_F1_R2:scan:2648:DAEKS[UNIMOD:21]PIN[UNIMOD:7]GR/2"] {
            let usi: USI = usi.parse().unwrap();
            let response = usi.get_spectrum(None).unwrap();
            assert!(!response.is_empty())
        }
    }
}
