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
//! # Instrument configuration database
//! Like other tooling that works with instrument vendor files, as new instruments are released, information about the model names and their
//! components takes time to percolate through the ecosystem. If `mzdata` encounters an instrument it doesn't recognize, it won't know how to
//! fill in the [`InstrumentConfiguration`](crate::meta::InstrumentConfiguration) for that instrument, and will `panic` by default. Set the
//! `MZDATA_IGNORE_UNKNOWN_INSTRUMENT` environment variable to `ignore` to tell the library to just assume this is the default (empty) instrument
//! configuration and carry on with a warning, or `silent` to suppress the warning too. Not all APIs may know how to cope with this!
mod instruments;
mod reader;

pub use reader::{is_thermo_raw_prefix, ThermoRawReader, ThermoRawReaderType};

#[cfg(feature = "async_partial")]
mod async_reader;
#[cfg(feature = "async_partial")]
pub use async_reader::{
    ThermoRawReader as AsyncThermoRawReader, ThermoRawReaderType as AsyncThermoRawReaderType,
};
