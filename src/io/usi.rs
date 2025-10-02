use std::{fmt::Display, str::FromStr};

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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Identifier {
    Scan(u64),
    Index(u64),
    NativeID(Box<Vec<u64>>),
}

#[derive(Debug, Default, Clone, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct USI {
    pub protocol: Protocol,
    pub dataset: String,
    pub run_name: String,
    pub identifier: Option<Identifier>,
    pub interpretation: Option<String>,
    pub provenance: Option<(Repository, String)>,
}

impl FromStr for USI {
    type Err = USIParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut this = Self::default();
        let mut tokens = s.splitn(6, ':');
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

                        let tail = tokens.next().map_or((None, None), |tail| {
                            for repository in Repository::ALL {
                                if let Some((i, p)) =
                                    tail.rsplit_once(&format!(":{}-", repository.code()))
                                {
                                    return (
                                        Some(i.to_string()),
                                        Some((*repository, p.to_string())),
                                    );
                                }
                            }
                            (Some(tail.to_string()), None)
                        });
                        this.interpretation = tail.0;
                        this.provenance = tail.1;
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
                    write!(f, ":{}-{}", provenance.0, provenance.1)?;
                }
            }
            Ok(())
        } else {
            write!(f, "{}:{}:{}", self.protocol, self.dataset, self.run_name)
        }
    }
}

/// A repository that can be used for provenance IDs.
#[allow(non_camel_case_types)]
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, Ord, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Repository {
    #[default]
    Pride,
    PeptideAtlas,
    MassIVE,
    jPOST,
    iProx,
    PanoramaPublic,
}

impl Repository {
    /// All repositories listed for programmatic access
    const ALL: &[Self] = &[
        Self::Pride,
        Self::PeptideAtlas,
        Self::MassIVE,
        Self::jPOST,
        Self::iProx,
        Self::PanoramaPublic,
    ];

    /// Get the repository code as used in a USI provenance ID
    pub const fn code(&self) -> &'static str {
        match self {
            Self::Pride => "PR",
            Self::PeptideAtlas => "PA",
            Self::MassIVE => "MA",
            Self::jPOST => "JP",
            Self::iProx => "IP",
            Self::PanoramaPublic => "PP",
        }
    }
}

impl std::fmt::Display for Repository {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.code())
    }
}

#[cfg(test)]
mod test {
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
        assert_eq!(usi.provenance, None);
        Ok(())
    }

    #[test]
    fn test_example_modification_shielded() -> Result<(), USIParseError> {
        let usi: USI = "mzspec:PXD047679:20221122_EX3_UM7_Kadav001_SA_EXT00_16_DSS_11.raw:scan:40889:GGK[xlink:dss[138]#XLDSS]IEVQLK//KVESELIK[#XLDSS]PINPR/4".parse()?;
        assert_eq!(usi.dataset, "PXD047679");
        assert_eq!(
            usi.run_name,
            "20221122_EX3_UM7_Kadav001_SA_EXT00_16_DSS_11.raw"
        );
        assert_eq!(usi.identifier, Some(Identifier::Scan(40889)));
        assert_eq!(
            usi.interpretation,
            Some("GGK[xlink:dss[138]#XLDSS]IEVQLK//KVESELIK[#XLDSS]PINPR/4".to_string())
        );
        assert_eq!(usi.provenance, None);
        Ok(())
    }

    #[test]
    fn test_example_labile_shielded() -> Result<(), USIParseError> {
        let usi: USI = "mzspec:PXD047679:20221122_EX3_UM7_Kadav001_SA_EXT00_16_DSS_11.raw:scan:40889:{Cation:K}GGKIEVQLK/4".parse()?;
        assert_eq!(usi.dataset, "PXD047679");
        assert_eq!(
            usi.run_name,
            "20221122_EX3_UM7_Kadav001_SA_EXT00_16_DSS_11.raw"
        );
        assert_eq!(usi.identifier, Some(Identifier::Scan(40889)));
        assert_eq!(
            usi.interpretation,
            Some("{Cation:K}GGKIEVQLK/4".to_string())
        );
        assert_eq!(usi.provenance, None);
        Ok(())
    }

    #[test]
    fn test_example_global_shielded() -> Result<(), USIParseError> {
        let usi: USI = "mzspec:PXD047679:20221122_EX3_UM7_Kadav001_SA_EXT00_16_DSS_11.raw:scan:40889:<[Glu->pyro-Glu]@N-term:E>GGKIEVQLK/4".parse()?;
        assert_eq!(usi.dataset, "PXD047679");
        assert_eq!(
            usi.run_name,
            "20221122_EX3_UM7_Kadav001_SA_EXT00_16_DSS_11.raw"
        );
        assert_eq!(usi.identifier, Some(Identifier::Scan(40889)));
        assert_eq!(
            usi.interpretation,
            Some("<[Glu->pyro-Glu]@N-term:E>GGKIEVQLK/4".to_string())
        );
        assert_eq!(usi.provenance, None);
        Ok(())
    }

    #[test]
    fn test_example_many_shielded() -> Result<(), USIParseError> {
        let usi: USI = "mzspec:PXD047679:20221122_EX3_UM7_Kadav001_SA_EXT00_16_DSS_11.raw:scan:40889:<[Cation:K]@N-term:E>GGKIE[Cation:K|Info:(((((([]}}}}}}}}<<<<<<{{{)(>>[((((){>>>><>>))}}{}]}]VQLK/4".parse()?;
        assert_eq!(usi.dataset, "PXD047679");
        assert_eq!(
            usi.run_name,
            "20221122_EX3_UM7_Kadav001_SA_EXT00_16_DSS_11.raw"
        );
        assert_eq!(usi.identifier, Some(Identifier::Scan(40889)));
        assert_eq!(
            usi.interpretation,
            Some("<[Cation:K]@N-term:E>GGKIE[Cation:K|Info:(((((([]}}}}}}}}<<<<<<{{{)(>>[((((){>>>><>>))}}{}]}]VQLK/4".to_string())
        );
        assert_eq!(usi.provenance, None);
        Ok(())
    }

    #[test]
    fn test_provenance() -> Result<(), USIParseError> {
        let usi: USI = "mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLH[UNIMOD:1]PLEGAVVIIFK/2:PR-G47".parse()?;
        assert_eq!(usi.dataset, "PXD000561");
        assert_eq!(usi.run_name, "Adult_Frontalcortex_bRP_Elite_85_f09");
        assert_eq!(usi.identifier, Some(Identifier::Scan(17555)));
        assert_eq!(
            usi.interpretation,
            Some("VLH[UNIMOD:1]PLEGAVVIIFK/2".to_string())
        );
        assert_eq!(usi.provenance, Some((Repository::Pride, "G47".to_string())));
        Ok(())
    }
}
