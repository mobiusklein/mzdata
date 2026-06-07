# Working with Ion Mobility Frames

## Table of contents

- [Working with Ion Mobility Frames](#working-with-ion-mobility-frames)
  - [Table of contents](#table-of-contents)
  - [What you will learn](#what-you-will-learn)
  - [From spectra to frames](#from-spectra-to-frames)
    - [What makes a frame different](#what-makes-a-frame-different)
    - [Opening an ion mobility source](#opening-an-ion-mobility-source)
  - [The IonMobilityFrameLike trait](#the-ionmobilityframelike-trait)
    - [Familiar metadata](#familiar-metadata)
    - [The ion mobility dimension](#the-ion-mobility-dimension)
    - [Raw 3D arrays](#raw-3d-arrays)
  - [Feature extraction with mzsignal](#feature-extraction-with-mzsignal)
    - [Why features, not peaks](#why-features-not-peaks)
    - [Extracting features](#extracting-features)
    - [Querying the feature map](#querying-the-feature-map)
  - [MS hierarchy and frame groups](#ms-hierarchy-and-frame-groups)
  - [Putting it together](#putting-it-together)

## What you will learn

By the end of this tutorial you will be able to:

- Open a file with ion mobility data and obtain an [`IonMobilityFrameSource`]
- Inspect frame metadata using the [`IonMobilityFrameLike`] trait
- Extract ion mobility features from an MS1 frame using `mzsignal`
- Query a `FeatureMap` by m/z to find matching features
- Iterate over MS1/MS2 frame groups

This tutorial assumes you have read the [spectrum tutorial](crate::tutorial::spectrum) and are
comfortable with the core `mzdata` spectrum API. Ion mobility frames map closely onto that
API — the table in [From spectra to frames](#from-spectra-to-frames) shows the correspondence.

## From spectra to frames

### What makes a frame different

A standard mass spectrum is two-dimensional: m/z versus intensity. An ion mobility frame adds
a third axis — **ion mobility** — measured as either drift time (milliseconds) or inverse
reduced mobility (1/K₀, V·s/cm²). Each point in a frame represents a (m/z, ion mobility,
intensity) triple, forming a surface rather than a curve.

Processing that surface yields **features** rather than peaks. A
[`Feature<MZ, IonMobility>`](mzpeaks::feature::Feature) is a contiguous region of the
(m/z, ion mobility) plane: it has a centroid m/z, a start and end on the mobility axis, and
a total integrated intensity.

The table below maps the spectrum concepts you already know to their frame equivalents:

| Spectrum concept | Frame equivalent |
|---|---|
| [`MultiLayerSpectrum`] | [`MultiLayerIonMobilityFrame`] |
| [`SpectrumLike`] | [`IonMobilityFrameLike`] |
| [`SpectrumSource`] | [`IonMobilityFrameSource`] |
| [`SpectrumGroup`] | [`IonMobilityFrameGroup`] |
| [`BinaryArrayMap`] (2D) | [`BinaryArrayMap3D`] (3D) |
| [`MZPeakSetType`](mzpeaks::peak_set::MZPeakSetType) | [`FeatureMap<MZ, IonMobility>`](mzpeaks::feature_map::FeatureMap) |
| [`MultiLayerSpectrum::pick_peaks`] | [`MultiLayerIonMobilityFrame::extract_features_simple`] |

### Opening an ion mobility source

Ion mobility data can come from many places, but the common :

- **Bruker TDF** files (`.d/` directories), accessed with the `bruker_tdf` feature
- **mzML** files that store 3D ion mobility data as separate scan events per mobility bin or that combine all mobility bins into a single array.

In both cases the workflow is: open an [`MZReader`], then convert it into an
[`IonMobilityFrameSource`] using [`IntoIonMobilityFrameSource::try_into_frame_source`].
The [`IntoIonMobilityFrameSource::has_ion_mobility`] method lets you check at runtime
whether the source contains ion mobility data before committing to the conversion.

The two generic parameters to `try_into_frame_source` are the centroid-feature type and the
deconvoluted-feature type. For most work, use `Feature<MZ, IonMobility>` and
`ChargedFeature<Mass, IonMobility>` from the prelude:

```rust
# use std::io;
# use mzdata::prelude::*;
# fn open_frame_source() -> io::Result<()> {
let mut reader = mzdata::MZReader::open_path(
    "test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz"
)?;

if reader.has_ion_mobility().is_some() {
    let mut frame_source = reader
        .try_into_frame_source::<
            Feature<MZ, IonMobility>,
            ChargedFeature<Mass, IonMobility>,
        >()
        .expect("file has ion mobility data");
    println!("Opened ion mobility source with {} frames", frame_source.len());
}
#   Ok(())
# }
```

## The IonMobilityFrameLike trait

### Familiar metadata

[`IonMobilityFrameLike`] mirrors [`SpectrumLike`] for frame metadata. The same access
pattern applies:

```rust
# use std::io;
# use mzdata::prelude::*;
# fn frame_metadata() -> io::Result<()> {
# let mut reader = mzdata::MZReader::open_path(
#     "test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz"
# )?;
# let mut frame_source = reader
#     .try_into_frame_source::<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>>()
#     .expect("file has ion mobility data");
# let frame = frame_source.get_frame_by_index(0).unwrap();
let id    = frame.id();
let level = frame.ms_level();
let time  = frame.start_time();
println!("{id} | MS{level} | {time:.2} min");
#   Ok(())
# }
```

Random-access retrieval mirrors the spectrum API:

```rust
# use std::io;
# use mzdata::prelude::*;
# fn frame_access() -> io::Result<()> {
# let mut reader = mzdata::MZReader::open_path(
#     "test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz"
# )?;
# let mut frame_source = reader
#     .try_into_frame_source::<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>>()
#     .expect("file has ion mobility data");
// By index
let _frame = frame_source.get_frame_by_index(0).unwrap();
// By native ID
let _frame = frame_source
    .get_frame_by_id("merged=42926 frame=9728 scanStart=1 scanEnd=705");
// By retention time (minutes)
let _frame = frame_source.get_frame_by_time(43.5);
#   Ok(())
# }
```

### The ion mobility dimension

The key new metadata field is the **ion mobility unit**, which tells you how to interpret the
mobility axis of the raw data:

```rust
# use std::io;
# use mzdata::prelude::*;
# fn im_unit() -> io::Result<()> {
# let mut reader = mzdata::MZReader::open_path(
#     "test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz"
# )?;
# let mut frame_source = reader
#     .try_into_frame_source::<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>>()
#     .expect("file has ion mobility data");
# let frame = frame_source.get_frame_by_index(0).unwrap();
// Unit::VoltSecondPerSquareCentimeter (1/K0) for PASEF, or Unit::Millisecond for drift tubes
let unit = frame.ion_mobility_unit();
println!("Ion mobility unit: {unit:?}");
#   Ok(())
# }
```

### Raw 3D arrays

For profile-mode frames, [`IonMobilityFrameLike::raw_arrays`] returns the
[`BinaryArrayMap3D`] — a three-axis array holding m/z, ion mobility, and intensity. This is
the frame analogue of the two-axis [`BinaryArrayMap`] from [`SpectrumLike::raw_arrays`].

## Feature extraction with mzsignal

The following requires the `mzsignal` feature:

```toml
[dependencies]
mzdata = { version = "*", features = ["mzsignal", "nalgebra"] }
```

### Why features, not peaks

A centroid peak is a single (m/z, intensity) point. In a frame, the equivalent is a
**feature**: a contiguous ridge in the (m/z, ion mobility) plane spanning multiple
mobility bins. A feature records the m/z centroid, the ion mobility range it covers, and
the area (integrated intensity) under that ridge.

### Extracting features

`MultiLayerIonMobilityFrame::extract_features_simple` finds features in a frame:

```ignore
# use std::io;
# use mzdata::prelude::*;
# use mzsignal::Tolerance;
# fn extract(path: &str) -> io::Result<()> {
let mut reader = mzdata::MZReader::open_path(path)?;
let mut frame_source = reader
    .try_into_frame_source::<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>>()
    .expect("file has ion mobility data");

let mut frame = frame_source.get_frame_by_index(0).unwrap();

// 15 ppm m/z tolerance, at least 2 consecutive bins, max 0.1-unit gap allowed
frame.extract_features_simple(Tolerance::PPM(15.0), 2, 0.1, None).unwrap();

let feature_map = frame.features.as_ref().unwrap();
println!("Extracted {} features", feature_map.len());
#   Ok(())
# }
```

The four parameters are:
- `error_tolerance` — m/z tolerance for grouping signal across mobility bins
- `min_length` — minimum number of consecutive bins a feature must span
- `maximum_gap_size` — maximum allowed gap (in ion mobility units) within a feature
- `peak_picker` — optional custom peak picker; `None` uses the default

### Querying the feature map

[`FeatureMap::all_features_for`] returns all features whose m/z centroid falls within the
requested tolerance of a query m/z. Each feature exposes its properties through the
[`FeatureLike`] trait:

- `coordinate()` — centroid m/z (via [`CoordinateLike<MZ>`](mzpeaks::prelude::CoordinateLike))
- `start_time()` / `end_time()` — the ion mobility range (return `Option<f64>`)
- `apex_time()` — the ion mobility at peak signal
- `area()` — integrated intensity across the feature

```ignore
# use std::io;
# use mzdata::prelude::*;
# use mzsignal::Tolerance;
# fn query(path: &str) -> io::Result<()> {
# let mut reader = mzdata::MZReader::open_path(path)?;
# let mut frame_source = reader
#     .try_into_frame_source::<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>>()
#     .expect("file has ion mobility data");
# let mut frame = frame_source.get_frame_by_index(0).unwrap();
# frame.extract_features_simple(Tolerance::PPM(15.0), 2, 0.1, None).unwrap();
let feature_map = frame.features.as_ref().unwrap();

let query_mz = 1456.95;
let hits = feature_map.all_features_for(query_mz, Tolerance::PPM(15.0));
for feature in &hits {
    println!(
        "m/z {:.4}  IM [{:.4}, {:.4}]  area {:.0}",
        feature.coordinate(),
        feature.start_time().unwrap_or(0.0),
        feature.end_time().unwrap_or(0.0),
        feature.area(),
    );
}
#   Ok(())
# }
```

## MS hierarchy and frame groups

Just as spectra form MS1/MS2 trees captured in [`SpectrumGroup`], frames form
[`IonMobilityFrameGroup`] instances. The grouping iterator pairs each MS1 precursor frame
with its MS2 product frames:

```rust
# use std::io;
# use mzdata::prelude::*;
# fn frame_groups() -> io::Result<()> {
# let mut reader = mzdata::MZReader::open_path(
#     "test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz"
# )?;
# let mut frame_source = reader
#     .try_into_frame_source::<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>>()
#     .expect("file has ion mobility data");
for group in frame_source.groups() {
    if let Some(precursor) = group.precursor() {
        println!(
            "MS1 frame {} @ {:.2} min — {} MS2 products",
            precursor.id(),
            precursor.start_time(),
            group.products().len()
        );
    }
}
#   Ok(())
# }
```

[`IonMobilityFrameGrouping`] provides the same accessors as [`SpectrumGrouping`]:
`precursor()`, `products()`, `earliest_frame()`, `highest_ms_level()`, and so on.

## Putting it together

Here is a complete example that opens a 3D mzML file, iterates its MS1 frame groups,
extracts features from each MS1 frame, and prints the ten most intense features.
This requires the `mzsignal` feature.

```ignore
use std::io;

use mzdata::prelude::*;
use mzsignal::Tolerance;

fn summarize_ms1_features(path: &str) -> io::Result<()> {
    let mut reader = mzdata::MZReader::open_path(path)?;
    let mut frame_source = reader
        .try_into_frame_source::<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>>()
        .expect("file has ion mobility data");

    for group in frame_source.groups() {
        let (precursor_frame, _products) = group.into_parts();
        let Some(mut frame) = precursor_frame else { continue };
        if frame.ms_level() != 1 { continue }

        frame.extract_features_simple(Tolerance::PPM(15.0), 2, 0.1, None).unwrap();

        let fm = frame.features.as_ref().unwrap();
        let mut features: Vec<_> = fm.iter().collect();
        features.sort_by(|a, b| b.area().partial_cmp(&a.area()).unwrap());

        println!(
            "Frame {} @ {:.2} min — {} features",
            frame.id(), frame.start_time(), fm.len()
        );
        for feat in features.iter().take(10) {
            println!(
                "  m/z {:.4}  IM [{:.4}, {:.4}]  area {:.0}",
                feat.coordinate(),
                feat.start_time().unwrap_or(0.0),
                feat.end_time().unwrap_or(0.0),
                feat.area(),
            );
        }
    }
    Ok(())
}
```
