/*!
The data structures and components that represent a mass spectrum and how
to access their data.

A mass spectrum is made up of multiple components, the spectrum's signal data
itself, plus all the metadata that describes how that data was acquired by
the instrument.

# Components
- [`bindata`] includes structures for dealing with raw binary data arrays that may or may not
be byte-encoded but not strongly typed, though it does not include signal processing as that is outside the scope of
this crate.

# Spectra

Represent the collection of attributes and data that compose a single mass spectrum.

Because a mass spectrum may be obtained from sources with varying levels of detail,
several alternative structures are provided with a common set of trait-based methods
to unify access:

1. [`RawSpectrum`] for representing a spectrum that has not been decoded into distinct
       peaks yet, but whose data may be continuous or discrete.
2. [`CentroidSpectrum`] for representing spectra from sources which are guaranteed to
       be pre-centroided, like those from MGF files or other simple text representations.
3. [`Spectrum`] for representing a multi-layer representation of a spectrum where both
       raw data and a distinct peak list are available.

These structures all implement the [`SpectrumDescription`] trait

The [`SpectrumDescription`] trait is included in the crate prelude, and gives the caller
read-only access to components that describe a spectrum's metadata.

```rust
use std::fs::File;
use mzdata::prelude::*;
use mzpeaks::{Tolerance, prelude::*};
use mzdata::io::MzMLReader;
use mzdata::spectrum::{SignalContinuity};

let reader = MzMLReader::new(File::open("./test/data/small.mzML").unwrap());
for spectrum in reader {
    println!("Scan {} => BP {}", spectrum.id(), spectrum.peaks().base_peak().mz);
    if spectrum.signal_continuity() < SignalContinuity::Profile {
        let peak_picked = spectrum.into_centroid().unwrap();
        println!("Matches for 579.155: {:?}",
                 peak_picked.peaks.all_peaks_for(
                     579.155, Tolerance::Da(0.02)));
    }
}
```

*/

pub mod bindata;
pub(crate) mod chromatogram;
pub(crate) mod group;
pub(crate) mod scan_properties;
pub(crate) mod spectrum;
pub mod utils;

pub use crate::spectrum::bindata::{ArrayType, BinaryArrayMap, BinaryDataArrayType, DataArray};
pub use crate::spectrum::chromatogram::{Chromatogram, ChromatogramLike};
pub use crate::spectrum::scan_properties::*;
pub use crate::spectrum::spectrum::{
    CentroidSpectrum, CentroidSpectrumType, DeconvolutedSpectrum, DeconvolutedSpectrumType,
    MultiLayerSpectrum, PeakDataLevel, RawSpectrum, Spectrum, SpectrumConversionError,
    SpectrumLike, SpectrumProcessingError,
};

pub use group::{
    SpectrumGroup, SpectrumGroupIntoIter, SpectrumGroupIter, SpectrumGroupingIterator,
};

#[cfg(feature = "mzsignal")]
pub use group::{
    average_spectra, DeferredSpectrumAveragingIterator, SpectrumAveragingIterator,
    SpectrumGroupAveraging,
};
