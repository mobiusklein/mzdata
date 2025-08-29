# Accessing mass spectrometry data files with `mzdata`

## Table of contents

- [Accessing mass spectrometry data files with `mzdata`](#accessing-mass-spectrometry-data-files-with-mzdata)
  - [Table of contents](#table-of-contents)
  - [Open a `SpectrumSource`:](#open-a-spectrumsource)
    - [From a file path](#from-a-file-path)
    - [From an `io::Read`](#from-an-ioread)
      - [When you know the format](#when-you-know-the-format)
        - [No `io::Seek` support](#no-ioseek-support)
      - [When you don't know the format](#when-you-dont-know-the-format)
      - [Using `MZReader`](#using-mzreader)
      - [Manually handling format dispatch](#manually-handling-format-dispatch)
      - [Using `mz_read!`](#using-mz_read)
      - [`MassSpectrometryReadWriteProcess`](#massspectrometryreadwriteprocess)
  - [Getting a spectrum](#getting-a-spectrum)
    - [By native ID](#by-native-id)
    - [By index](#by-index)
    - [By scan start time](#by-scan-start-time)
  - [What can you do with a spectrum?](#what-can-you-do-with-a-spectrum)
  - [Fancy iterators](#fancy-iterators)
    - [Random access iterators](#random-access-iterators)
    - [Grouping iterators](#grouping-iterators)

## Open a `SpectrumSource`:

All the spectrum reading operations in `mzdata` rely on a core trait called [`SpectrumSource`] which provides the shared set of behaviors for a mass spectrometry data file. [`SpectrumSource`] requires that the type also be an [`Iterator`], and provides some basic features for random access when supported. See [Fancy iterators](#fancy-iterators) for further examples of specializing iterators.

There are many different ways to access data, on disk, in memory, or from some other kind of byte stream
transformation. They can change the ways we can interact with that data, be it the type of the data source,
whether it is random access, or what format it is stored in.

### From a file path

`mzdata` can read files on disk, accessed by path, easily. In this example we'll use the [`MZReader`] type
to figure out which reader to use for us automatically.

```rust
use std::io;

use mzdata;
use mzdata::prelude::*;

fn from_path() -> io::Result<()> {
    let file_path = "test/data/batching_test.mzML";
    let reader = mzdata::MZReader::open_path(file_path)?;
    println!("Found {} spectra", reader.len());
    Ok(())
}
```

When reading a file from disk, `mzdata` can make certain assumptions like that the file supports
the [`io::Seek`](std::io::Seek) trait and can read or build indices over the file quickly and guarantee
that the file supports full random access, like [`RandomAccessSpectrumIterator`].

Additionally, some binary formats like [`ThermoRawReader`], [`MzMLbReader`], or [`TDFFrameReader`] _require_ that there be a file on disk that exists outside of the Rust model of the file system in order to read it.

### From an `io::Read`

When you have a readable type that makes no guarantees about what backs it, `mzdata` can still read it, but the
code may be more complicated. If you already know what kind of file you're dealing with, you can directly create
the reader type you want:

#### When you know the format

```rust
use std::{io, fs};
use mzdata::prelude::*;

// If the readable source supports `io::Seek`
fn from_read_seek() -> io::Result<()> {
    let mut handle = fs::File::open("test/data/batching_test.mzML")?;
    let reader = mzdata::MzMLReader::new_indexed(&mut handle);

    println!("Found {} spectra", reader.len());
    Ok(())
}
```

##### No `io::Seek` support

When the source doesn't support [`io::Seek`](std::io::Seek), most reader types don't
even implement [`SpectrumSource`], although they still implement
[`Iterator`]. The [`StreamingSpectrumIterator`]
wrapper does implement parts of [`SpectrumSource`] using less efficient
mechanism, but in situations where it cannot satisfy the request, it will `panic` instead. It also,
naturally doesn't support reading spectra that have already been seen as the stream cannot be reversed.

```rust
use std::{io, fs};
use mzdata::prelude::*;
use mzdata::io::StreamingSpectrumIterator;

fn from_read() -> io::Result<()> {
    let mut handle = fs::File::open("test/data/batching_test.mzML")?;
    let reader = StreamingSpectrumIterator::new(mzdata::MzMLReader::new(&mut handle));

    println!("File reports {:?} spectra, but without an index this isn't guaranteed.", reader.spectrum_count_hint());
    Ok(())
}
```

#### When you don't know the format

If you don't know the format, it's still possible to do so if your source supports [`io::Seek`](std::io::Seek), and this lets you
be more flexible about things like file compression where the higher level features shown previously do not.

#### Using `MZReader`

The [`MZReader`] type, or more generally [`MZReaderType`] is an algebraic data type that wraps a data source, and proxies [`SpectrumSource`] and related traits to one of the underlying implementations.

```rust
use std::{io, fs};

use mzdata;
use mzdata::io::MZReader;
use mzdata::prelude::*;

fn dispatching() -> io::Result<()> {
    let reader = MZReader::open_path("test/data/batching_test.mzML")?;
    println!("Found {} spectra in {:?}", reader.len(), reader.as_format());
    Ok(())
}

```

If you need to customize the how the underlying reader works, [`MZReader`] is just an `enum` so you are free to use pattern matching to access it. For instance, if you have the `thermo` feature enabled, you can instruct the reader to use the Thermo peak centroiding algorithm on the backend before passing the data from .NET to Rust:

```rust

use std::{io, fs};

use mzdata;
#[cfg(feature = "thermo")]
use mzdata::io::{MZReader, thermo::ThermoRawReader};
use mzdata::prelude::*;

#[cfg(feature = "thermo")]
fn dispatching_features() -> io::Result<()> {
    let mut reader = MZReader::open_path("test/data/small.RAW")?;
    if let MZReader::ThermoRaw(thermo) = &mut reader {
        // Do the specializing operation here
        thermo.set_centroiding(true);
    } else {
        // Otherwise no peak picking specialization is done
    }
    Ok(())
}

```

#### Manually handling format dispatch

```rust
use std::{io, fs};

use mzdata;
use mzdata::io::{MassSpectrometryFormat, infer_from_stream, RestartableGzDecoder};
use mzdata::prelude::*;

fn infer_format_from_read() -> io::Result<()> {
    let mut handle = fs::File::open("test/data/batching_test.mzML")?;

    let (format_type, gzipped) = infer_from_stream(&mut handle)?;
    match format_type {
        MassSpectrometryFormat::MzML => {
            if !gzipped {
                let reader = mzdata::MzMLReader::new_indexed(&mut handle);
                println!("Found {} spectra", reader.len());
            } else {
                let mut handle = RestartableGzDecoder::new(io::BufReader::new(handle));
                let reader = mzdata::MzMLReader::new_indexed(&mut handle);
                println!("Found {} spectra", reader.len());
            }
        }
        MassSpectrometryFormat::MGF => {
            if !gzipped {
                let reader = mzdata::MGFReader::new_indexed(&mut handle);
                println!("Found {} spectra", reader.len());
            } else {
                let mut handle = RestartableGzDecoder::new(io::BufReader::new(handle));
                let reader = mzdata::MGFReader::new_indexed(&mut handle);
                println!("Found {} spectra", reader.len());
            }
        }

        // Handle more formats if they're appropriate here
        _ => {
            // ...
        }
    }
    Ok(())
}
```

Another drawback of this approach is that a reader type for a compressed file is fundamentally a different
type, so you either must use trait object over `dyn SpectrumSource`, replicate the logic using the reader
instance for each combination of formats and compression states, or move all logic operating on the reader
instance to another function that is generic over the `SpectrumSource` plus whatever other required traits
you want to use.

#### Using `mz_read!`

When you don't _need_ to keep the reader around beyond the current scope, the [`mz_read!`](crate::mz_read) macro can substantially simplify. It is like [`MZReader`], but without the trivial dispatch overhead, provided that the reader instance only lives as long as the enclosing scope:

```rust
use std::io;

use mzdata::Spectrum;
use mzdata::prelude::*;

fn mz_read_path() -> io::Result<()> {
    mzdata::mz_read!("test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz".as_ref(), reader => {
        println!("{} spectra found", reader.len());
    })?;
    Ok(())
}
```

Again, the drawback is that the reader itself cannot leave the scope created inside the macro, so you cannot hold onto the instance, nor do you know what kind of reader you have so you cannot customize what to do in each case. The indirection incurred by [`MZReader`] is trivial compared to the cost of all the other operations involved in doing I/O. This method is kept for the sake of brevity.

#### `MassSpectrometryReadWriteProcess`

[`MassSpectrometryReadWriteProcess`] can help work around this limitation, but it carries its own restrictions. For an example of that, see [mzconvert](https://github.com/mobiusklein/mzdata/blob/main/examples/mzconvert.rs) example program.

## Getting a spectrum

Once you've opened a [`SpectrumSource`], you can access spectra in several ways:

1. Incrementally using [`Iterator::next`]
2. By native ID using [`get_spectrum_by_id`](SpectrumSource::get_spectrum_by_id)
3. By index position using [`get_spectrum_by_index`](SpectrumSource::get_spectrum_by_index)
4. By scan start time using [`get_spectrum_by_time`](SpectrumSource::get_spectrum_by_time)

Using iterator methods is very flexible, and already well described in other places.

### By native ID

Most mass spectrometry instrument vendors have some form of textual notation for specifically identifying
a spectrum. Often this just a function of the the "scan number", which is synonymous with "index + 1", but
not always. The native ID is what is written in mzML files' `id` attribute for each spectrum.

```rust
# use std::{io, fs};
# use mzdata::prelude::*;
# fn from_read_seek() -> io::Result<()> {
#    let mut reader = mzdata::MZReader::open_path("test/data/batching_test.mzML")?;
if let Some(spec) = reader.get_spectrum_by_id("controllerType=0 controllerNumber=1 scan=25788") {
    println!("{} @ MS level {} has {} data points with mode {}",
        spec.id(),
        spec.ms_level(),
        spec.peaks().len(),
        spec.signal_continuity()
    );
}
#    Ok(())
# }
```

### By index

A spectrum's index is just that, its base-0 ordinal in a collection of spectra.

```rust
# use std::{io, fs};
# use mzdata::prelude::*;
# fn from_read_seek() -> io::Result<()> {
#    let mut reader = mzdata::MZReader::open_path("test/data/batching_test.mzML")?;
if let Some(spec) = reader.get_spectrum_by_index(48) {
    println!("{} @ MS level {} has {} data points with mode {}",
        spec.id(),
        spec.ms_level(),
        spec.peaks().len(),
        spec.signal_continuity()
    );
}
#    Ok(())
# }
```

### By scan start time

A spectrum's "scan start time" is the point in time when the spectrum was acquired. This assumes
that the source MS run is stored in chronological order, which should be true for most files from
raw MS runs, but not all, nor for programmatically generated files, especially if they are merged
across multiple MS runs.

```rust
# use std::{io, fs};
# use mzdata::prelude::*;
# fn from_read_seek() -> io::Result<()> {
#    let mut reader = mzdata::MZReader::open_path("test/data/batching_test.mzML")?;
if let Some(spec) = reader.get_spectrum_by_time(120.218212668) {
    println!("{} @ MS level {} has {} data points with mode {}",
        spec.id(),
        spec.ms_level(),
        spec.peaks().len(),
        spec.signal_continuity()
    );
}
#    Ok(())
# }
```

## What can you do with a spectrum?

A mass spectrum can come in a variety of shapes, but most of them are a combination of "spectrum acquisition details" and "signal data". `mzdata` provides several types for dealing with this variation, but uses the [`SpectrumLike`] trait to provide a common interface for most of their facets.

`mzdata` represents the "spectrum acquisition details" with the [`SpectrumDescription`] type. Most of the [`SpectrumLike`] interface can be treated as just an abstraction over [`SpectrumDescription`]-containing types.

The "signal data" component is more varied. There are raw spectra which contain profile or continuous spectra, similar to what the instrument sees ([`RawSpectrum`]). There are centroid spectra which contain discrete points in the m/z dimension with a measured intensity ([`CentroidSpectrum`]), often produced from profile spectra post-acquisition or centroided by the instrument itself during data acquisition. A spectrum may be centroided, but also deisotoped and charge state deconvolved, where the discrete points are in the neutral mass dimension with a known charge ([`DeconvolutedSpectrum`]). A spectrum might be in any of these three states, or multiple of them when transforming from one to the other ([`MultiLayerSpectrum`]).

Because a single mass spectrometry data file may contain spectra in different states, `mzdata` always reads [`MultiLayerSpectrum`] instances, but it is possible to convert between these four types. With the [`mzsignal`](https://crates.io/crates/mzsignal) library `mzdata` can also perform peak picking to centroid profile spectra when dealing with raw signal [`MultiLayerSpectrum::pick_peaks`] and [`RawSpectrum::pick_peaks_into`].

For more details, see the [spectrum tutorial](crate::tutorial::spectrum).

## Fancy iterators

### Random access iterators

The [`SpectrumSource`] trait requires that the type already be an [`Iterator`], which can be convenient for reading all the spectra in a data source sequentially. Sometimes we might want to iterate over a subset of spectra which appear part-way through an MS run. The natural solution would be to use [`Iterator::skip`] or [`Iterator::filter`], but these would *require* that the [`SpectrumSource`] consume all the intervening spectra until it reached some you were interested in. The [`RandomAccessSpectrumIterator`] trait extends [`SpectrumSource`] iterators with methods to immediately jump to the requested location without needing to do the intermediate parsing that mirror the `get_spectrum_by_` methods described earlier:

1. [`RandomAccessSpectrumIterator::start_from_id`]
2. [`RandomAccessSpectrumIterator::start_from_index`]
3. [`RandomAccessSpectrumIterator::start_from_time`]

These alter the [`Iterator`] state, potentially encountering an error if the requested coordinate couldn't be located.

### Grouping iterators

Spectra naturally form units of hierarchy with MS levels and precursor-product relationships connecting them. A complete unit of
that hierarchy is described as a [`SpectrumGroup`], a simple implementation of the [`SpectrumingGrouping`](crate::spectrum::SpectrumGrouping)
trait. The [`SpectrumSource::groups`] or [`SpectrumSource::into_groups`] methods produce [`SpectrumGroupingIterator`], iterators
over [`SpectrumGroup`]. These work by buffering to build trees of spectra before producing them and provide many of the same
metadata and fancy iterator methods that their underlying [`SpectrumSource`] do.