//! Reader implementation for Thermo RAW files, [`ThermoRawReaderType`].
//!
//! Depends upon the [`thermorawfilereader`] crate which manages the self-hosted `dotnet`
//! runtime.
//!
//! ```no_run
//! use std::io;
//!
//! use mzdata::prelude::*;
//! use mzdata::io::ThermoRawReader;
//!
//! # fn main() -> io::Result<()> {
//! let mut reader = ThermoRawReader::open_path("./test/data/small.RAW")?;
//! let scan = reader.get_spectrum_by_index(0).unwrap();
//! assert_eq!(scan.index(), 0);
//! assert_eq!(reader.len(), 48);
//! #    Ok(())
//! # }
//! ```
//! # Licensing
//! By using this library, you agree to the [RawFileReader License](https://github.com/thermofisherlsms/RawFileReader/blob/main/License.doc)
//!
mod reader;
mod instruments;

pub use reader::{ThermoRawReaderType, ThermoRawReader, is_thermo_raw_prefix};
