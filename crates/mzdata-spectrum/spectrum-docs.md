# Documentation Review: mzdata-spectrum

Scope: `crates/mzdata-spectrum/src/` (lib.rs, spectrum_types.rs, peaks.rs, scan_properties.rs,
chromatogram.rs, frame.rs, utils.rs, group/*). No code was changed; this is a review only.

## How to read this report

- "Missing" means a public item has no `///` doc comment at all.
- "Thin" means a doc comment exists but omits information a caller would need to use the
  API correctly (units, algorithm, edge cases, panics, errors).
- The "Signal processing" section at the end collects every operation that actually
  transforms spectral signal (denoising, peak picking, reprofiling, feature extraction,
  summary statistics) because these deserve deeper treatment than ordinary accessors and
  are called out separately per your request.

---

## 1. Crate root (lib.rs)

- There is no crate-level (`//!`) documentation at all. `cargo doc` for this crate currently
  opens on an empty landing page. This is the highest-leverage gap in the whole crate: a
  short module-level doc would let readers on docs.rs immediately understand:
  - The three-tier spectrum model (`RawSpectrum` -> `CentroidSpectrumType` ->
    `DeconvolutedSpectrumType`) and how `MultiLayerSpectrum` unifies them.
  - How this crate relates to `mzdata-bindata` (arrays), `mzdata-param` (CV terms),
    `mzdata-meta` (run-level metadata), and the parent `mzdata` crate (I/O). None of that
    relationship is visible from inside this crate today; it currently only lives in the
    repository's top-level `CLAUDE.md`, which is not shipped with the published crate docs.
  - A minimal worked example (construct a spectrum, inspect `peaks()`, convert between
    layers) as a doctest. There are currently zero doctests anywhere in this crate.
- The re-export list itself has no grouping/commentary, so a docs.rs visitor sees a flat,
  undifferentiated list of ~30 names with no indication of which ones are the "entry
  points" (`MultiLayerSpectrum`, `SpectrumLike`, `SpectrumSource`-adjacent grouping types)
  versus supporting/internal types (`PeakDataIter`, `RawIter`, etc.).

## 2. spectrum_types.rs

- `SpectrumLike` (the central trait) has good per-method docs, but the trait itself has
  only a one-line summary. Given this is the primary abstraction of the whole crate, it
  deserves a paragraph explaining why `peaks()` returns an enum (`RefPeakDataLevel`) rather
  than a concrete type, and how that connects to the raw/centroid/deconvoluted model.
- `SpectrumLike::update_summaries`:
  - The doc does not say which four CV params it writes (total ion current, base peak m/z,
    base peak intensity, lowest/highest observed m/z), nor that it updates existing params
    in place if present and appends new ones otherwise. A caller cannot tell from the doc
    whether calling this twice is safe/idempotent (it is, except see below).
  - Possible accuracy problem worth a second look: in the "update an existing param"
    branch for `highest_mz_curie` (spectrum_types.rs around line 269), the value written is
    `mz_range.0` (the minimum), not `mz_range.1` (the maximum), while the "insert new"
    branch two lines later correctly uses `mz_range.1`. If that is really a logic bug and
    not intentional, then no doc comment can be "correct" until it's fixed; if it is
    intentional, the doc should explain why. Either way this is currently undocumented
    behavior that a reader would not expect from the method name. Flagging for
    verification rather than fixing, per instructions not to modify code.
- `RawSpectrum::mzs_mut` and `RawSpectrum::intensities_mut` have no doc comments at all
  (contrast with `mzs`/`intensities` just above them, which document their panic
  behavior).
- `MultiLayerSpectrum` struct fields `peaks` and `deconvoluted_peaks` use plain `//`
  comments instead of `///` (lines ~1004 and ~1007). Plain comments are invisible to
  `cargo doc`, so these two fields currently render with no description on docs.rs, unlike
  every other field in the same struct.
- `MultiLayerSpectrum::try_build_centroids` and `try_build_deconvoluted_centroids` have no
  doc comments at all. Their sibling `try_build_peaks` has a doc comment with `# Errors`.
  All three should say what "try build" means (attempt to materialize a peak list from
  `arrays` without discarding the arrays), and that they require
  `signal_continuity() == Centroid`, returning `NotCentroided` otherwise.
- `CentroidSpectrumType::reinterpret` / `DeconvolutedSpectrumType::reinterpret` say
  "potentially losing information" but do not say what information: fitted peak shape
  parameters, signal-to-noise, or any custom fields on `C`/`D` that are not representable
  as a `DataArray` are silently dropped when round-tripping through `BuildArrayMapFrom`.
- `SpectrumConversionError` / `SpectrumProcessingError`: each variant has a `#[error(...)]`
  message (used for `Display`), but none has a `///` doc comment. The `#[error]` message is
  not rustdoc, so docs.rs currently shows these variants with no description at all, only
  their bare names.
- Minor: `into_centroid` (on `RawSpectrum`) documents the `NotCentroided` error but does
  not clarify that this call trusts the recorded `signal_continuity` and does not itself
  inspect the array contents, i.e. it is not a substitute for peak picking. There is a "See
  also" pointing at `pick_peaks_into`, which helps, but a one-line explicit contrast
  ("does not pick peaks; only reinterprets already-centroided arrays") would remove
  ambiguity.

## 3. peaks.rs

- `PeakDataLevel` and `RefPeakDataLevel` enums have a doc comment on the type itself but
  none of their four variants (`Missing`, `RawData`, `Centroid`, `Deconvoluted`) are
  individually documented. Since these are the two most-matched-on enums in the crate,
  documenting what state each variant represents (and when a caller should expect to see
  it, e.g. `Missing` only occurs before any array/peak data has been attached) would help.
- `PeakDataLevel::try_build_peaks` (`pub(crate)`) says it "attempts to reconstruct one of
  the peak levels" but does not document the precedence rule: it checks
  `D::has_arrays_for` (deconvoluted) before `C::has_arrays_for` (centroid), so if arrays
  happen to satisfy both, the deconvoluted interpretation always wins silently. The same
  precedence exists in `MultiLayerIonMobilityFrame::try_build_features` (frame.rs) and
  should be documented in both places since it affects which peak type a caller gets back.
- `RefPeakDataLevel::search` has a `**NOTE**` about the `Deconvoluted` variant using "a
  different coordinate system" but does not say what that means concretely (the
  deconvoluted layer is indexed/searched by neutral mass produced via `mass_charge_ratio`,
  not by raw m/z), nor give guidance on how to convert a query m/z into that space before
  searching.
- `SummaryOps` (private trait) is reasonably documented per-method, but the "no guarantee
  of m/z order" warning only appears on `get`; it would be worth repeating (or centralizing
  in one place referenced by `#[doc = ...]` or a `# Note` cross-reference) since
  `PeakDataIter`/`RefPeakDataIter` both rely on it and a reader may only see one of the
  three near-duplicate types.

## 4. scan_properties.rs

- `IonMobilityMeasure::ion_mobility` has no doc comment, and this is the most significant
  correctness-adjacent gap in the crate. The method scans four different CURIEs and
  returns whichever value it finds first as a bare `f64`:
  - `MS:1002476` ion mobility drift time (milliseconds)
  - `MS:1002815` inverse reduced ion mobility (Vs/cm^2, i.e. 1/K0)
  - `MS:1001581` FAIMS compensation voltage (volts)
  - `MS:1003371` SELEXION compensation voltage (volts)

  These are physically different quantities with different units and different meanings
  (a drift time is not comparable to a compensation voltage). A caller who calls
  `ion_mobility()` without also calling `ion_mobility_type()` to check which CURIE matched
  has no way to know what the returned number means. This needs an explicit doc comment
  (and probably a `# Note`/warning admonition) explaining that the value's unit and
  interpretation depend on `ion_mobility_type()`, and that `SpectrumLike::ion_mobility`,
  `SelectedIon`, and `ScanEvent` all inherit this ambiguity because they all implement
  `IonMobilityMeasure`.
- `IsolationWindowBuilder` is a small state machine and the type-level doc explains its
  purpose, but the arithmetic performed when transitioning states is not documented and is
  easy to get wrong:
  - `.target(v)`: if the window is currently in `Offset` state, calling this rewrites
    `lower_bound`/`upper_bound` from "radius relative to target" into "absolute bound",
    via `lower_bound = target - lower_bound` and `upper_bound = target + upper_bound`. This
    is the crux of how the builder resolves offsets, and it is invisible from the method's
    one-line doc ("Set the target value").
  - `lower_offset`/`upper_offset` accept a delta from the target; `lower_limit`/
    `upper_limit` accept an absolute m/z value. The two families are easy to confuse by
    name alone and neither doc comment says "delta" or "absolute" explicitly.
  - A subtle, currently-silent hazard: `lower_limit`/`upper_limit` only act when
    `self.0.flags` is `Unknown`, `Explicit`, or `Complete`. If the window is in `Offset`
    state (because `lower_offset`/`upper_offset` was called first), calling `lower_limit`
    or `upper_limit` is a silent no-op, it neither errors nor panics, it just does nothing.
    Mixing offset-style and limit-style builder calls on the same window can therefore
    silently drop a call. This deserves an explicit warning in both method docs since
    nothing else would tip a caller off.
  - A worked example in the type-level doc (e.g. `.lower_offset(2.0).upper_offset(3.0)
    .target(500.0)` producing bounds `[498.0, 503.0]`) would make the state machine much
    easier to use correctly.
- `IsolationWindow::is_empty` doc says it is "distinct from `has_isolation`" but does not
  spell out the actual (slightly counter-intuitive) truth table: a window with
  `NoIsolation` flags and zero-valued bounds returns `false` from `is_empty`, because
  `is_empty` is defined as `lower_bound == 0.0 && upper_bound == 0.0 && has_isolation()`,
  i.e. "empty" only applies to a window that claims to isolate something but has no
  recorded width. Worth stating this outright with an example.
- `ScanWindow::contains` and `ScanWindow::is_empty` have no doc comments at all, while the
  parallel methods on `IsolationWindow` (`contains`, `is_empty`) do. These should at least
  get the one-line doc their `IsolationWindow` counterparts have.
- `Activation::energy` (`pub energy: f32`) has no doc comment and no stated unit. Whether
  this is eV, normalized collision energy, or some instrument-specific arbitrary scale is
  not documented anywhere in this file, and different activation methods in real-world
  files use different units for "energy" (a very common source of downstream bugs).
- `SelectedIon::intensity` (`pub intensity: f32`) has no doc comment and no stated unit
  (compare to `mz`, which is documented).
- `Acquisition`'s three public fields (`scans`, `combination`, `params`) have no per-field
  doc comments, unlike `ScanEvent`'s fields just below, which are all documented.
- `ChromatogramDescription`'s fields (`id`, `index`, `ms_level`, `polarity`,
  `chromatogram_type`, `params`, `precursor`, `products`) have zero per-field doc comments.
  This stands out because the near-identical `SpectrumDescription` struct in the same file
  documents every one of its fields. Bringing `ChromatogramDescription` up to the same
  standard would remove an inconsistency a reader is likely to notice.
- `ChromatogramType::is_aggregate`, `is_ion_current`, `is_electromagnetic_radiation` have
  no doc comments explaining the semantic grouping (e.g. why total-ion-current and
  base-peak chromatograms count as "aggregate" but selected-ion-monitoring does not).
  `ChromatogramType`'s individual variants are also undocumented; several (e.g.
  `FlowRateChromatogram`, `PressureChromatogram`) are not self-explanatory in an MS
  context to a non-specialist reader.
- `ScanCombination`'s variants carry `// MS:1000795`-style plain comments instead of `///`
  doc comments, so they render without any description on docs.rs.

## 5. chromatogram.rs

- `Chromatogram` (the struct) has no type-level doc comment at all. This is the most
  visible gap in the file, since every comparable struct in `spectrum_types.rs`
  (`RawSpectrum`, `CentroidSpectrumType`, etc.) has one. Its `arrays` field is also
  undocumented (its `description` field is private, so that omission is lower priority).
- `ChromatogramLike::time` / `ChromatogramLike::intensity` (the trait methods) and their
  concrete implementations on `Chromatogram` have no doc comments. It is not stated
  whether the time array is expected/required to be sorted ascending, what unit `time` is
  in (minutes, to match `ScanEvent::start_time`, presumably, but this is never stated in
  this file), or what happens when the arrays are absent versus empty.
- `TimeInterval<Time> for Chromatogram`'s `apex_time` and `area` implementations have no
  doc comments explaining what they compute or how. Concretely: `area` builds a
  `FeatureView` from the time/intensity arrays and delegates to
  `mzpeaks::feature::FeatureView::area`, i.e. this is (presumably) a trapezoidal
  integration over the recorded points, but nothing in this file says so, and nothing
  states the units of the result (intensity-seconds? intensity-minutes?) or that `area`
  silently returns `0.0` (not an error) when the arrays cannot be read.
- The `as_feature_view!` macro and `as_simple_feature` function are undocumented. They are
  private/`pub(crate)`, so this is lower priority, but a one-line comment on
  `as_feature_view!` would help future maintainers, since its early-return-via-`?`-like
  control flow inside a macro is not obvious at a glance.

## 6. frame.rs

- This file has the best signal-processing documentation in the crate
  (`extract_features_with`/`extract_features_simple` both have `# Arguments` sections),
  but still has real gaps, listed in the Signal Processing section below.
- `FeatureDataLevel` / `RefFeatureDataLevel`: same gap as `PeakDataLevel`/
  `RefPeakDataLevel` in peaks.rs, the enum has a summary doc but none of its four variants
  do.
- `MultiLayerIonMobilityFrame::description` (the struct field) has no doc comment, unlike
  its siblings `arrays`, `features`, `deconvoluted_features` in the same struct.
- `MultiLayerIonMobilityFrame::try_build_features` has no doc comment at all. It shares the
  same undocumented precedence rule as `PeakDataLevel::try_build_peaks` (checks
  deconvoluted-capable arrays before centroid-capable arrays) and the same implicit
  precondition (`signal_continuity() == Centroid`) that its sibling `try_build_peaks` in
  spectrum_types.rs also leaves undocumented.
- None of the `TryFrom`/`From` conversions between `RawSpectrum`, `MultiLayerSpectrum`, and
  `MultiLayerIonMobilityFrame` (roughly lines 458-542) have doc comments. These are not
  trivial field copies:
  - `From<MultiLayerIonMobilityFrame<C, D>> for RawSpectrum` "unstacks" the 3D array map
    into a flat 2D one and then writes the frame's ion mobility unit onto every
    ion-mobility-typed `DataArray`. None of that is explained, and the two internal
    `.unwrap()` calls on `unstack()` mean this conversion can panic for array
    configurations that fail to unstack; that possibility deserves a `# Panics` section.
  - `TryFrom<MultiLayerSpectrum<CPeak, DPeak>> for MultiLayerIonMobilityFrame<CFeat, DFeat>`
    only accepts frames built from raw arrays (returns `NotFound(IonMobilityArray)` if
    `arrays` is `None`, silently discarding any existing `peaks`/`deconvoluted_peaks` on
    the input) and this constraint is not documented on the impl.

## 7. group/ (spectrum.rs, frame.rs, utils.rs, mod.rs)

- `group/mod.rs` has no module-level doc. A short note that `IonMobilityFrameGrouping` is a
  parallel, near-duplicate implementation of `SpectrumGrouping` adapted for ion mobility
  frames (rather than an independent design) would help a reader who opens one file and
  wonders why the other looks identical.
- `SpectrumGroupIntoIter` / `SpectrumGroupIter` (and their `IonMobilityFrameGroup*`
  counterparts) have no doc comments describing their iteration order/contract: they always
  yield the precursor first (if present), then products in stored order. This ordering is
  a real behavioral guarantee, visible only by reading `next()`, and any caller relying on
  "precursor comes first" is relying on undocumented behavior.
- `SpectrumGrouping::lowest_ms_level` / `highest_ms_level`: the "no data" sentinel is
  `Option::None`, returned whenever the computed level is `0`. This conflates "the group is
  empty" with "every spectrum in the group legitimately has ms_level 0", and the crate
  elsewhere (`SpectrumLike::spectrum_type` doc in spectrum_types.rs) acknowledges that
  `mzdata` can hold non-MS spectra, which would plausibly report `ms_level() == 0`. The doc
  comment does not mention this edge case, and a caller cannot distinguish "empty group"
  from "group full of non-MS spectra" from the `None` return alone.
- `earliest_spectrum`/`latest_spectrum`: tie-breaking on equal `start_time` is not
  documented (uses `min_by`/`max_by`, whose tie-break behavior differs: `min_by` keeps the
  first equal element, `max_by` keeps the last). This is a minor point but worth a
  one-line note if determinism ever matters to a caller.
- `group/utils.rs`: `GroupIterState` has no doc comment on the enum or its variants. It is
  `pub(crate)`, so this is low priority, but one line per variant would make the iterator
  state machine easier to follow for future maintainers.

## 8. utils.rs

- `HasIonMobility` is the best-documented enum in the crate: it has a summary and every
  variant is documented. This is a good template to bring the other undocumented enums
  (`PeakDataLevel`, `RefPeakDataLevel`, `FeatureDataLevel`, `RefFeatureDataLevel`,
  `ChromatogramType`) up to.

---

## 9. Signal processing operations that need more detail

These are the operations that actually transform or summarize spectral signal. Because
they encode domain-specific algorithmic choices (not just data shuffling), under-describing
them is a bigger risk than the general accessor-documentation gaps above: a caller can
misuse these without any type error.

1. **`denoise` / `RawSpectrum::denoise`, `MultiLayerSpectrum::denoise`**
   (spectrum_types.rs). The doc is one line: "Apply a local denoising algorithm ... using
   `mzsignal::denoise`." Missing:
   - What the `scale` parameter means, its units, and what a reasonable range is (there is
     no `# Arguments` section at all).
   - A description (even a short one) of what "local denoising" does conceptually, e.g.
     whether it is a local noise-floor/background subtraction, a smoothing filter, or
     something else, so a reader is not forced to leave this crate and read `mzsignal`
     source to understand what will happen to their data.
   - That it mutates the intensity array in place and coerces its storage type to
     `Float32` (`store_as(BinaryDataArrayType::Float32)`), which is a side effect a caller
     working with a different original encoding (e.g. Float64) should know about.
   - Whether output intensities can go negative, or are clamped at zero.
   - No `# Errors` section despite returning `Result`.

2. **Peak picking: `pick_peaks`, `pick_peaks_with`, `pick_peaks_in_intervals`
   (`MultiLayerSpectrum`), and their `RawSpectrum` wrappers `pick_peaks_into`,
   `pick_peaks_with_into`, `pick_peaks_in_intervals_into`.** These have the best docs of
   the signal-processing group (the "if already centroided, no filtering is performed"
   behavior is called out), but still omit:
   - Any description of the underlying peak-detection algorithm (local maxima plus a
     curve fit controlled by `PeakFitType`, currently hardcoded to `PeakFitType::Quadratic`
     in every convenience method) or how the signal-to-noise threshold is estimated from
     the surrounding signal. A reader has to go to `mzsignal::peak_picker::PeakPicker` to
     find out, and even the one-line summary of what SNR estimation method is used would
     help a reader decide whether their data is a good fit.
   - The fact that `signal_to_noise_threshold` is unitless/relative rather than an absolute
     intensity cutoff is never stated.
   - No `# Panics` or `# Errors` sections on any of these methods even though they return
     `Result<_, SpectrumProcessingError>` and can fail via `PeakPickerError` or array
     access errors.
   - `pick_peaks_in_intervals` silently produces an empty peak set for any interval that
     matches no data; this normal-but-easy-to-miss behavior is not called out.

3. **`reprofile_with_shape` / `reprofile_with_shape_into`** (spectrum_types.rs). `dx` and
   `fwhm` are reasonably explained, but:
   - The doc never states that the result is a synthetic/theoretical profile (every peak
     is rendered as an idealized Gaussian; this is a simulation for visualization or
     downstream algorithms that require profile-shaped input, not a reconstruction of the
     original instrument signal). A reader could easily mistake the output for "the
     original profile, recovered."
   - The peak shape is hardcoded to `PeakShape::Gaussian` in this convenience method; the
     doc mentions the option to use `PeakSetReprofiler` directly for "non-uniform shapes"
     but does not mention that *this* method is Gaussian-only, i.e. it does not surface
     the constraint it is itself subject to.
   - Grid memory scales as `(max_mz - min_mz) / dx`; the doc gives a suggested `dx` range
     but does not connect that to the actual memory/array-size consequence in concrete
     terms (e.g. "a 1:2000 m/z range at dx=0.001 produces ~2,000,000 points per array").
   - No `# Errors` section despite returning `Result`.

4. **`extract_features_with` / `extract_features_simple`** (frame.rs). Best-documented
   signal processing code in the crate (has `# Arguments` and `# Type Arguments`), but:
   - Still no conceptual description of the linking/extraction algorithm itself, i.e. how
     points across ion mobility scans are associated into a single feature/trace beyond
     naming the three tunable parameters. A one-paragraph description (or a link to the
     relevant `mzsignal::feature_mapping` type-level docs) would let a reader reason about
     what `error_tolerance`/`maximum_gap_size` actually do to their data instead of
     treating them as opaque knobs.
   - `extract_features_with` calls `panic!("Cannot extract feature, unknown signal
     continuities")` when `signal_continuity()` is `Unknown`. This panic is completely
     undocumented; there is no `# Panics` section, and nothing in the caller-facing
     `IonMobilityFrameLike::signal_continuity` doc warns that this value must be
     Centroid or Profile before calling this method.
   - The `MapState`/`E: MapState<...>` type parameter is described as opening "a wide
     range of customization options" without saying what any of those options are or
     giving even one example beyond the default `PeakMapState` used by
     `extract_features_simple`.

5. **`SpectrumLike::update_summaries`** (spectrum_types.rs). Computes and writes back four
   CV params (TIC, base peak m/z, base peak intensity, m/z range) via a full linear scan
   over `self.peaks()`. The doc does not name which params are written, does not mention
   that it updates in place if the param CURIE already exists, and (see section 2 above)
   the branch that updates an existing "highest observed m/z" param appears to write the
   wrong tuple element. Whatever the resolution of that discrepancy, the doc should state
   precisely which four CURIEs are affected and the update-vs-insert behavior so a caller
   is not surprised by duplicate or stale params.

6. **`try_build_peaks` / `try_build_centroids` / `try_build_deconvoluted_centroids`**
   (spectrum_types.rs) **and `try_build_features`** (frame.rs). These convert raw
   `BinaryArrayMap`(3D) data into typed peak/feature collections. All four share an
   undocumented precedence rule (deconvoluted-capable arrays are preferred over
   centroid-capable arrays whenever both are technically satisfiable) and an undocumented
   precondition (`signal_continuity() == Centroid`; on `Profile` data these calls either
   return `NotCentroided` or, for `try_build_peaks`, silently no-op and return the current,
   possibly `Missing`, peak level). Two of the four (`try_build_centroids`,
   `try_build_deconvoluted_centroids`) have no doc comment at all.

7. **`IonMobilityMeasure::ion_mobility`** (scan_properties.rs). Not a transform, but it
   feeds directly into ion-mobility-aware processing elsewhere in the crate
   (`SpectrumLike::ion_mobility`/`has_ion_mobility`/`has_ion_mobility_class`), so getting
   its unit ambiguity documented (see section 4 above) matters for anyone building
   signal-processing logic on top of it, e.g. code that buckets or filters spectra by ion
   mobility value without checking `ion_mobility_type()` first would silently mix drift
   times and compensation voltages.

---

## Suggested priority order

If addressing incrementally, the following order maximizes correctness value per change:

1. `IonMobilityMeasure::ion_mobility` unit ambiguity (scan_properties.rs) - safety-relevant.
2. Verify and document (or fix) the `update_summaries` highest-m/z branch
   (spectrum_types.rs) - correctness-relevant.
3. `IsolationWindowBuilder` offset/limit mixing hazard (scan_properties.rs) -
   correctness-relevant, silent no-op.
4. Crate-level (`//!`) documentation in lib.rs - highest onboarding value.
5. Signal processing methods (`denoise`, `pick_peaks*`, `reprofile_with_shape*`,
   `extract_features_*`) - fill in `# Arguments`, `# Panics`, `# Errors`, and one paragraph
   of algorithm context each.
6. Fix the two `//` (non-doc) field comments on `MultiLayerSpectrum` in spectrum_types.rs.
7. Bring `ChromatogramDescription`, `ChromatogramType`, `Chromatogram` up to the same
   per-field documentation standard already used by `SpectrumDescription`.
8. Document the four peak/feature-data enums' variants (`PeakDataLevel`,
   `RefPeakDataLevel`, `FeatureDataLevel`, `RefFeatureDataLevel`) using `HasIonMobility` in
   utils.rs as the template.
9. Remaining missing doc comments listed above (mut accessors, `try_build_*` methods,
   conversions in frame.rs, group iterator ordering guarantees).
