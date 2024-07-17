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
3. [`MultiLayerSpectrum`] for representing a multi-layer representation of a spectrum where both
       raw data and a distinct peak list are available.
4. [`DeconvolutedSpectrum`] for representing spectra from sources which are guaranteed to be
    pre-centroided, deisotoped and charge state deconvoluted.

These structures all implement the [`SpectrumLike`] trait

The [`SpectrumLike`] trait is included in the crate prelude, and gives the caller
read-only access to components that describe a spectrum's metadata.

```rust
use mzpeaks::Tolerance;
use mzdata::MzMLReader;
use mzdata::prelude::*;
use mzdata::spectrum::SignalContinuity;

let reader = MzMLReader::open_path("./test/data/small.mzML").unwrap();
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
pub(crate) mod spectrum_types;
pub(crate) mod frame;
pub mod utils;

pub use crate::spectrum::bindata::{ArrayType, BinaryArrayMap, BinaryDataArrayType, DataArray};
pub use crate::spectrum::chromatogram::{Chromatogram, ChromatogramLike};
pub use crate::spectrum::scan_properties::*;
pub use crate::spectrum::spectrum_types::{
    CentroidSpectrum, CentroidSpectrumType, DeconvolutedSpectrum, DeconvolutedSpectrumType,
    MultiLayerSpectrum, RefPeakDataLevel, RawSpectrum, Spectrum, SpectrumConversionError,
    SpectrumLike, SpectrumProcessingError, PeakDataLevel, CentroidPeakAdapting,
    DeconvolutedPeakAdapting
};

pub use frame::{IonMobilityFrameDescription, MultiLayerIonMobilityFrame, IonMobilityFrameLike};

pub use group::{
    SpectrumGroup, SpectrumGroupIntoIter, SpectrumGroupIter, SpectrumGroupingIterator
};

#[cfg(feature = "mzsignal")]
pub use group::{
    average_spectra, DeferredSpectrumAveragingIterator, SpectrumAveragingIterator,
    SpectrumGroupAveraging,
};
