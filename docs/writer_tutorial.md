# Writing mass spectrometry data files with `mzdata`

## Table of contents
- [Writing mass spectrometry data files with `mzdata`](#writing-mass-spectrometry-data-files-with-mzdata)
  - [Table of contents](#table-of-contents)
  - [Creating a `SpectrumWriter`](#creating-a-spectrumwriter)
  - [Copying across metadata](#copying-across-metadata)
  - [Writing spectra](#writing-spectra)

## Creating a `SpectrumWriter`

`mzdata` uses the [`SpectrumWriter`] trait to define shared writing behavior across different
data file formats. The [`MzMLWriter`] type writes spectra to mzML files, and [`MGFWriter`] writes
spectra to MGF files, though other formats may be available as well. [`SpectrumWriter`] is agnostic
to the thing being written to, and [`MzMLWriter`] and [`MGFWriter`] can accept an arbitrary [`Write`]
implementation, though not all implementations are.

```rust
use std::{io, fs};

use mzdata::io::MzMLWriter;
use mzdata::prelude::*;

fn main() -> io::Result<()> {
    let fh = io::BufWriter::new(fs::File::create("tmp.mzML")?);
    let mut writer = MzMLWriter::new(fh);

    Ok(())
}
```

## Copying across metadata

If you have a source which implements [`MSDataFileMetadata`], you can easily copy metadata from your
source file to your output file, preserving the trace of information about where your data came from.

```rust
use std::{io, fs};

use mzdata::io::MzMLWriter;
use mzdata::prelude::*;

fn main() -> io::Result<()> {
    let reader = mzdata::MZReader::open_path("test/data/batching_test.mzML")?;

    let fh = io::BufWriter::new(fs::File::create("tmp.mzML")?);
    let mut writer = MzMLWriter::new(fh);

    writer.copy_metadata_from(&reader);

    // mzML files want to list how many spectra they contain
    if let Some(n_spectra) = reader.spectrum_count_hint() {
        writer.set_spectrum_count(n_spectra)
    } else {
        writer.set_spectrum_count(reader.len() as u64)
    }

    Ok(())
}
```

## Writing spectra

```rust
# use std::{io, fs};
# use mzdata::io::MzMLWriter;
# use mzdata::prelude::*;
#
#
# fn main() -> io::Result<()> {
#     let reader = mzdata::MZReader::open_path("test/data/batching_test.mzML")?;
#
#     let fh = io::BufWriter::new(fs::File::create("tmp.mzML")?);
#     let mut writer = MzMLWriter::new(fh);
#
#     writer.copy_metadata_from(&reader);
#
#     // mzML files want to list how many spectra they contain
#     if let Some(n_spectra) = reader.spectrum_count_hint() {
#         writer.set_spectrum_count(n_spectra)
#     } else {
#         writer.set_spectrum_count(reader.len() as u64)
#     }
#
    // Write spectra out one at a time, by reference
    for spec in reader {
        writer.write(&spec)?;
    }
#    Ok(())
# }
```

```rust
# use std::{io, fs};
# use mzdata::io::MzMLWriter;
# use mzdata::prelude::*;
#
#
# fn main() -> io::Result<()> {
#     let reader = mzdata::MZReader::open_path("test/data/batching_test.mzML")?;
#
#     let fh = io::BufWriter::new(fs::File::create("tmp.mzML")?);
#     let mut writer = MzMLWriter::new(fh);
#
#     writer.copy_metadata_from(&reader);
#
#     // mzML files want to list how many spectra they contain
#     if let Some(n_spectra) = reader.spectrum_count_hint() {
#         writer.set_spectrum_count(n_spectra)
#     } else {
#         writer.set_spectrum_count(reader.len() as u64)
#     }
#
    // Write out an iterator over spectra, in this case
    // using the owning variant for an iterator over owned
    // instances.
    writer.write_all_owned(reader)?;
#    Ok(())
# }
```


