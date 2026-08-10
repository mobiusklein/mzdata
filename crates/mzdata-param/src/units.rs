use std::fmt::Display;

use crate::{CURIE, ControlledVocabulary, Param};

/// Easily define all units. See the usage for more details. The accession and name both have to be
/// included in `&'static str` and `byte slice` to allow for constant time matching.
macro_rules! units {
    [$($unit:ident, $accession:literal, $baccession:literal, $name:literal, $bname:literal, $cv:ident, $id:literal);*;] => {
        /// Units that a term's value might have
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
        #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
        pub enum Unit {
            Unknown,
            $($unit,)*
        }

        impl Unit {
            pub const fn for_param(&self) -> (&'static str, &'static str) {
                match self {
                    $(Self::$unit => ($accession, $name),)*
                    Self::Unknown => ("","")
                }
            }

            pub const fn from_name(name: &str) -> Unit {
                match name.as_bytes() {
                    $($bname => Self::$unit,)*
                    _ => Self::Unknown
                }
            }

            pub const fn from_accession(acc: &str) -> Unit {
                match acc.as_bytes() {
                    $($baccession => Self::$unit,)*
                    _ => Self::Unknown
                }
            }

            pub const fn from_curie(acc: &CURIE) -> Unit {
                match acc {
                    $(CURIE {
                        controlled_vocabulary: ControlledVocabulary::$cv,
                        accession: $id,
                    } => Self::$unit,)*
                    _ => Self::Unknown
                }
            }

            pub const fn to_curie(&self) -> Option<CURIE> {
                match self {
                    $(Self::$unit => Some(CURIE {
                        controlled_vocabulary: ControlledVocabulary::$cv,
                        accession: $id,
                    }),)*
                    Self::Unknown => None
                }
            }

            pub const fn from_param(param: &Param) -> Unit {
                param.unit
            }

            pub const fn is_unknown(&self) -> bool {
                matches!(self, Self::Unknown)
            }
        }
    };
}

units![
    AbsorbanceUnit, "UO:0000269", b"UO:0000269", "absorbance unit", b"absorbance unit", UO, 269;
    Celsius, "UO:0000027", b"UO:0000027", "degree Celsius", b"degree Celsius", UO, 27;
    CountsPerSecond, "MS:1000814", b"MS:1000814", "counts per second", b"counts per second", MS, 1000814;
    DetectorCounts, "MS:1000131", b"MS:1000131", "number of detector counts", b"number of detector counts", MS, 1000131;
    Dimensionless, "UO:0000186", b"UO:0000186", "dimensionless unit", b"dimensionless unit", UO, 186;
    Electronvolt, "UO:0000266", b"UO:0000266", "electronvolt", b"electronvolt", UO, 266;
    Kelvin, "UO:0000012", b"UO:0000012", "kelvin", b"kelvin", UO, 12;
    Mass, "UO:000221", b"UO:000221", "dalton", b"dalton", UO, 221;
    MicrolitersPerMinute, "UO:0000271", b"UO:0000271", "microliters per minute", b"microliters per minute", UO, 271;
    Millisecond, "UO:0000028", b"UO:0000028", "millisecond", b"millisecond", UO, 28;
    Minute, "UO:0000031", b"UO:0000031", "minute", b"minute", UO, 31;
    MZ, "MS:1000040", b"MS:1000040", "m/z", b"m/z", MS, 1000040;
    Nanometer, "UO:0000018", b"UO:0000018", "nanometer", b"nanometer", UO, 18;
    Micrometer, "UO:0000017", b"UO:0000017", "micrometer", b"micrometer", UO, 17;
    Millimeter, "UO:0000016", b"UO:0000016", "millimeter", b"millimeter", UO, 16;
    Centimeter, "UO:0000015", b"UO:0000015", "centimeter", b"centimeter", UO, 15;
    PartsPerMillion, "UO:0000169", b"UO:0000169", "parts per million ", b"parts per million ", UO, 169;
    Pascal, "UO:0000110", b"UO:0000110", "pascal", b"pascal", UO, 110;
    Percent, "UO:0000187", b"UO:0000187", "percent", b"percent", UO, 187;
    PercentBasePeak, "MS:1000132", b"MS:1000132", "percent of base peak", b"percent of base peak", MS, 1000132;
    PercentBasePeakTimes100, "MS:1000905", b"MS:1000905", "percent of base peak times 100", b"percent of base peak times 100", MS, 1000905;
    Psi, "UO:0010052", b"UO:0010052", "pounds per square inch", b"pounds per square inch", UO, 10052;
    Second, "UO:0000010", b"UO:0000010", "second", b"second", UO, 10;
    Volt, "UO:0000218", b"UO:0000218", "volt", b"volt", UO, 218;
    VoltSecondPerSquareCentimeter, "MS:1002814", b"MS:1002814", "volt-second per square centimeter", b"volt-second per square centimeter", MS, 1002814;
    Hertz, "UO:000106", b"UO:000106", "hertz", b"hertz", UO, 106;
    Liter, "UO:0000099", b"UO:0000099", "liter", b"liter", UO, 99;
    Milliliter, "UO:0000098", b"UO:0000098", "milliliter", b"milliliter", UO, 98;
    Microliter, "UO:0000101", b"UO:0000101", "microliter", b"microliter", UO, 101;
];

impl Default for Unit {
    fn default() -> Self {
        Self::Unknown
    }
}

impl Display for Unit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(format!("{:?}", self).as_str())
    }
}
