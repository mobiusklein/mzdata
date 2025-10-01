//! Reader implementation for Thermo RAW files, [`ThermoRawReaderType`].
//!
//! **Requires the `thermo` feature**
//!
//! Depends upon the [`thermorawfilereader`] crate which manages the self-hosted `.NET`
//! runtime. You must still have a working [`.NET 8`](https://dotnet.microsoft.com/en-us/download/dotnet/8.0) runtime installed on the machine you
//! wish to run this on until Thermo's library supports .NET ahead-of-time compilation. For scripted installation of the .NET runtime
//! see <https://dotnet.microsoft.com/en-us/download/dotnet/scripts>.
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
mod instruments;
mod reader;

pub use reader::{ThermoRawReader, ThermoRawReaderType, is_thermo_raw_prefix};

#[cfg(feature = "async")]
mod async_reader;
#[cfg(feature = "async")]
pub use async_reader::{
    ThermoRawReader as AsyncThermoRawReader, ThermoRawReaderType as AsyncThermoRawReaderType,
};
