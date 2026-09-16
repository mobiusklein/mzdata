# Documentation Review: mzdata-bindata

Scope: `crates/mzdata-bindata/src/` (lib.rs, encodings.rs, array.rs, conversion.rs, map.rs,
traits.rs, utils.rs). No code was changed; this is a review only.

## How to read this report

- "Missing" means a public item has no `///` doc comment at all.
- "Thin" means a doc comment exists but omits information a caller needs (units, byte vs.
  element offsets, panics, errors, algorithm behavior).
- Several findings below go beyond documentation and describe places where the code's
  observable behavior looks like it does not match what a reasonable reader would expect
  from its name or its (thin) doc comment. These are flagged explicitly as "worth
  verifying" rather than asserted as confirmed bugs, since no code was run or changed as
  part of this review, per your instructions. They are included because they were found by
  trying to document the code accurately, which is exactly the exercise you asked for.

---

## 1. Crate root (lib.rs)

- This is the one crate in the workspace with a real `//!` crate-level doc comment, and it
  is a reasonable summary. It could still be improved:
  - It does not mention that several compression schemes are lossy (Numpress family) while
    others are lossless (Zlib, Zstd, dictionary encoding), which is the single most
    important fact a new user of this crate needs before picking a `BinaryCompressionType`.
  - It does not mention that a meaningful chunk of the public surface
    (`BinaryCompressionType` variants `NumpressPIC`, `NumpressPICZlib`, `NumpressPICZstd`,
    `LinearPrediction`, `DeltaPrediction`) is declared (round-trips through CV accessions)
    but not actually implemented by `DataArray::encode_bytestring`/`decode` (see section 9
    below). A one-line pointer to "supported vs. reserved" compression variants at the
    crate level would save readers from discovering this by trial and error (a panic).
  - No doctest/worked example exists anywhere in the crate showing the intended encode ->
    store -> decode round trip, or the recommended way to add data to a fresh `DataArray`.

## 2. encodings.rs

- `as_bytes` has no doc comment, and its behavior is a real hazard next to its neighbor
  `to_bytes`, which is documented as converting "into an owned buffer of little-endian
  bytes." `as_bytes` does not convert anything: it is a zero-copy `bytemuck::cast_slice`
  reinterpretation of whatever the host's native byte order happens to be. On a
  little-endian host the two happen to agree; on a big-endian host they do not. Since
  mzML/most MS formats are little-endian on disk, using `as_bytes` where `to_bytes` was
  intended would silently produce wrong-endian output on big-endian systems. This
  distinction needs to be stated explicitly on `as_bytes` (and cross-referenced from
  `to_bytes`), not left to the reader to infer from the name.
- The `byte_rotation` module (`transpose_f32/f64/i32/i64`, `reverse_transpose_f32/f64/
  i32/i64`) is re-exported (`pub use byte_rotation::*;`) but has zero doc comments on any
  function. These implement a byte-plane transpose ("shuffle"), the same idea as
  blosc's/HDF5's shuffle filter: for N-byte values, all first bytes are grouped together,
  then all second bytes, and so on, which tends to make the output more compressible by a
  general-purpose compressor (zstd) when the underlying values vary smoothly (e.g. a
  sorted m/z array). None of this is explained anywhere; a reader has to reverse-engineer
  the purpose from the fact that `ShuffleZstd`/`DeltaShuffleZstd` call it before
  compressing. There is also a known-incomplete part flagged only in an inline (non-doc)
  comment inside `transpose_bytes_into`: `// This isn't endian-correct yet, see note in
  'transpose_bytes_into'`. That caveat needs to be promoted to a real `///` doc comment
  (with a `# Note` or similar) on the public functions that call it
  (`encode_dict_indices` in the dictionary encoder, and by extension anything using
  `ShuffleZstd`/`DeltaShuffleZstd` on a big-endian host), since right now the warning is
  invisible to anyone who does not read the private implementation module.
- The `dictionary_encoding` module (`DictionaryEncoder`, `DictionaryDecoder`,
  `dictionary_encoding`, `dictionary_decoding`) implements a full custom on-disk format
  (an 8-byte data offset, an 8-byte value count, a table of unique values, then an
  index buffer whose element width is chosen from the table's cardinality) with zero
  documentation anywhere: no description of the wire format, no explanation of when this
  encoding is a good choice (low-cardinality/repeated-value arrays, e.g. integer charge
  states or categorical flags) versus a poor one (high-cardinality continuous floating
  point data, where the value table itself approaches the size of the original data).
  Since this module happens to sit inside a private (`mod encodings;`) parent module, it
  is not reachable from outside the crate today, but it is exercised from `array.rs`
  (`compress_dict_zstd`/`decompress_dict_zstd`, which back the public
  `BinaryCompressionType::ZstdDict` variant), so this is still a real, user-facing feature
  that is entirely undocumented, just one level removed from the public surface.
- `ArrayType` (the enum naming what an array measures) has a good type-level doc, but none
  of its ~25 variants have per-variant docs. The ion mobility variants are the most in
  need of this: there are three parallel families (`Raw*`, `Mean*`, `Deconvoluted*`)
  crossed with three IM representations (`IonMobility`, `DriftTime`,
  `InverseReducedIonMobility`), and nothing in this file explains what "Raw" (per-scan,
  before any averaging), "Mean" (averaged, e.g. the value assigned to an extracted ion
  mobility feature, produced by `mzdata-spectrum`'s feature extraction), and
  "Deconvoluted" (assigned after charge/isotope deconvolution) mean in practice, nor when
  code should convert between them via `as_mean_ion_mobility`/`as_raw_ion_mobility`/
  `as_deconvoluted_ion_mobility` (which themselves have one-line docs that just restate
  the method name).
- `ArrayType::as_param`, `as_param_const`, and `as_param_with_unit_const` are three
  overlapping ways to get a CV `Param`/`ParamCow` for an array type. Each has a short doc,
  and two of them correctly warn that `NonStandardDataArray` panics in a `const` context,
  but there is no top-level guidance on which of the three a caller should reach for
  (allocating vs. `const`-context vs. explicit-unit-override). There is also an
  undocumented inconsistency between them: `as_param_const` bakes in
  `Unit::VoltSecondPerSquareCentimeter` as the default unit for the
  "InverseReducedIonMobility" array family but leaves the "DriftTime" array family
  defaulted to `Unit::Unknown`. A caller who does not read every match arm would not know
  that drift-time arrays need an explicit unit to be set while inverse-reduced-mobility
  arrays do not.
- `BinaryDataArrayType::swap_bytes` has a doc comment that is grammatically incomplete:
  "Byte order swap the data stored in" (the sentence has no object). This should read
  something like "Byte-order swap the data stored in `data`, interpreting it as elements
  of `self`'s size."
- `BinaryCompressionType` has a reasonable type-level doc but zero per-variant
  documentation across 17 variants, several of which are materially different algorithms
  with different tradeoffs (see section 9, which covers this in depth since it is
  squarely a signal-processing documentation gap, not just an API polish item).
- `ArrayRetrievalError::DataTypeSizeMismatch` is documented as "The requested data type
  does not match the number of bytes available in the buffer," but it is also returned by
  `ByteArrayView::to_f32`/`to_f64`/`to_i32`/`to_i64` (traits.rs) for the unrelated case of
  an unsupported/`Unknown` source dtype, where there is no size mismatch at all, just an
  unsupported conversion. The doc comment should be broadened, or a distinct variant
  should exist (not a code change requested here, just flagging that the current doc does
  not cover one of the variant's real uses).
- `linear_prediction_decoding`, `linear_prediction_encoding`, `delta_decoding`,
  `delta_encoding` are public (re-exported from lib.rs) and have zero doc comments. These
  implement predictive coding transforms (the wire formats behind
  `BinaryCompressionType::LinearPrediction`/`DeltaPrediction`, see section 9). Given they
  are public API, at minimum each needs: what the transform does conceptually, what the
  minimum input length is (both silently no-op below a length threshold: `delta_decoding`/
  `delta_encoding` below 2 elements, `linear_prediction_decoding`/`linear_prediction_
  encoding` below 3 elements for the encoder and below 2 for the decoder, but the decoder
  reads `values[2]` before checking that length, see section 9), and the round-trip
  contract (`decoding(encoding(x)) == x`). See section 9 for a more serious concern about
  whether `linear_prediction_decoding` is even self-consistent as currently written.

## 3. array.rs

- `DataArray`'s type-level doc comment is the best in the crate: it explains the
  performance/misuse tradeoff clearly and points at `decode_and_store`. This is a good
  template; most other types in this crate (and in `mzdata-spectrum`) fall well short of
  it.
- The four basic constructors `new()`, `from_name()`, `from_name_and_type()`, and
  `from_name_type_size()` all have zero doc comments. `from_name_type_size`'s `size`
  parameter is not documented as a byte count (it is passed straight to
  `Bytes::with_capacity(size)`), which matters because `size` is not "number of elements."
- `slice()` and `slice_buffer()` have no doc comments, and in both cases `start`/`end` are
  **byte offsets into the encoded buffer**, not element indices. This is easy to get
  wrong (`slice(0, 10)` does not return the first 10 `f64` values; for an 8-byte-wide
  dtype it would fail the `(end - start) % dtype.size_of() != 0` check unless `end` is a
  multiple of 8), and nothing in either method's signature or documentation says so. The
  same applies to the crate-private `decoded_slice`.
- `wrap()` documents that "the data are already in native byte order," but this is
  inconsistent terminology next to `to_bytes` in encodings.rs, which promises
  little-endian output specifically. "Native" (meaning "whatever this CPU uses") and
  "little-endian" are the same thing only on little-endian hosts. Since MS file formats
  are little-endian on disk, a big-endian-host caller who takes this doc literally could
  wrap already-big-endian bytes and end up with silently corrupted data on decode. This
  should say "little-endian" explicitly if that is really the intended contract (which
  the rest of the crate's byte-order handling, e.g. `BinaryDataArrayType::swap_bytes`
  being invoked under `#[cfg(target_endian = "big")]` elsewhere, strongly suggests it is).
- `encode_bytestring` has a one-line doc for what is, by line count, the largest and most
  branch-heavy method in the crate (about 150 lines, 15+ match arms). Missing:
  - A `# Panics` section. The method panics (not returns `Err`) in several concrete,
    reachable situations: requesting `BinaryCompressionType::Decoded` as the target
    ("Should never happen"); requesting any Numpress variant on non-`Float64` /
    non-`Float32`-or-`Float64` data depending on the sub-variant; and requesting any
    compression variant with no implemented codec (`NumpressPIC*`, `LinearPrediction`,
    `DeltaPrediction`), which falls into the generic `_ => panic!("Compression type {:?}
    is unsupported", compression)` arm. None of this is discoverable without reading the
    full match.
  - The "fast path" (if the array is already stored in the requested compression, the
    existing buffer is cloned and returned unchanged) is only visible via a `log::trace!`
    call, not documented for callers who care about that optimization.
- `decode()` has a good doc comment covering the `Cow`-borrow-when-possible behavior and
  documents that it may fail, but the internal `base64_decode!` macro
  (`unwrap_or_else(|e| panic!(...))`) means malformed base64 input causes a panic, not the
  `Err(ArrayRetrievalError)` the function's `Result` return type implies. This should
  either be called out as a documented panic condition, or (separately, as a code
  question worth asking, not something this review is fixing) reconsidered so decode
  errors are consistently returned as `Err`.
- `decode_mut()` has no doc comment at all (contrast with `decode()` and
  `decode_and_store()` just above it, both of which are documented).
- `clear()` has no doc comment, and its behavior is broader than the name suggests: it
  clears `self.data` **and** discards `self.params` (sets it to `None`) **and** resets the
  cached item count. A caller reading only the name "clear" would reasonably expect only
  the byte buffer to be affected (as with, say, `Vec::clear`), not the array's metadata
  parameters.
- `store_as()`'s doc ("Recode the stored data as the requested binary data type") does not
  mention that narrowing conversions are lossy (Float64 -> Float32 loses precision; any
  float -> integer conversion truncates rather than rounds, per `AsPrimitive`'s semantics
  in traits.rs) and gives no guidance on when this is safe to call.
- In the "compression codecs" `impl DataArray` block, only `compress_zstd` has a doc
  comment (and a good one, calling out the `MZDATA_ZSTD_LEVEL` environment variable). None
  of `compress_numpress_linear`, `compress_numpress_slof`, `decompress_numpress_linear`,
  `decompress_numpress_slof`, `compress_delta_zstd`, `compress_dict_zstd`,
  `decompress_zstd`, `decompress_delta_zstd`, or `decompress_dict_zstd` have any doc
  comment, despite several being `pub`. These are exactly the operations this review was
  asked to look at closely; see section 9.

## 4. conversion.rs

- `ArraysAvailable` (the enum returned by `has_arrays_for`/`has_arrays_3d_for`) has no doc
  comment at all, on the type or any of its three variants (`Unknown`, `Ok`,
  `MissingArrays`). This matters because `Unknown` is a distinct, real outcome (the target
  type never declared `arrays_required()`, so the check is inconclusive) that callers
  must not treat the same as either "definitely available" or "definitely missing" and
  is easy to mishandle in a `match` that only distinguishes `Ok` from "not Ok."
- `BuildFromArrayMap::try_from_arrays` and `from_arrays` have no doc comments.
  `from_arrays` is a thin `.unwrap()` wrapper around `try_from_arrays` with no `# Panics`
  section, even though it panics on any array-retrieval error. This "try_* is fallible,
  bare name panics" naming convention recurs across the crate (`try_from_arrays_3d`/
  `from_arrays_3d` here, `try_build_peaks`/etc. in `mzdata-spectrum`) but is never stated
  once anywhere as a crate-wide convention.
- `has_arrays_for`'s doc ("A pre-emptive check for the presence of the required arrays")
  does not mention that it returns `ArraysAvailable::Unknown` whenever the implementing
  type does not override `arrays_required()` (returns `None`), which is the single most
  important thing to know about this method's contract.
- `BuildArrayMapFrom::arrays_included`/`as_arrays` have no doc comments. `as_arrays` is the
  serialization half of the round trip and, for several impls in this file, encodes
  nontrivial, undocumented conventions:
  - `Feature<MZ, IonMobility>::as_arrays` (and the `ChargedFeature<Mass, IonMobility>`
    equivalent) flattens a `Vec<Feature>` into parallel m/z/intensity/ion-mobility arrays
    plus a synthetic, non-standard `"feature identifier array"` (an `Int32` array
    recording which source feature each point came from), and globally sorts every point
    across all features by `(m/z, ion_mobility, source_feature_index)` before writing.
    None of this, the marker-array convention or the global sort, is written down
    anywhere. Anyone hand-assembling a `BinaryArrayMap` to round-trip through
    `try_from_arrays` has to reverse-engineer the exact convention from source.
  - The corresponding `try_from_arrays` for `Feature<MZ, IonMobility>` assumes the marker
    array's values are small, zero-based, contiguous integers (it sizes the output `Vec`
    as `marker_array.max() `, then indexes `features[key as usize]` directly). A
    marker array that violates this (negative values are impossible since it is read via
    `to_i32` into `usize`, but non-contiguous or non-zero-based values are entirely
    possible if hand-constructed) will either allocate a larger-than-necessary `Vec` or
    panic on out-of-bounds indexing. This precondition is not documented as a `# Panics`
    note or otherwise.
- `BuildFromArrayMap3D::has_arrays_3d_for`'s algorithm is subtle and undocumented: it
  iterates the per-ion-mobility-value 2D slices and returns `ArraysAvailable::Ok` as soon
  as it finds **one** slice where all required arrays are present, rather than requiring
  all slices to satisfy the requirement. A 3D map where the required arrays are only
  present in some ion mobility slices (e.g. a partially-populated import) will report
  `Ok` based on the first satisfying slice. This "first match wins" semantics is not
  mentioned in the (nonexistent) doc comment on this method.
- `BuildArrayMap3DFrom::as_arrays_3d`'s default implementation calls
  `BuildArrayMapFrom::as_arrays(source).try_into().unwrap()`, which can panic (inherits
  every panic condition of `BinaryArrayMap3D::stack`, e.g. a NaN ion mobility value, see
  section 5) and has no doc comment mentioning this.
- None of the concrete `impl BuildFromArrayMap for X` / `impl BuildArrayMapFrom for X`
  blocks (`CentroidPeak`, `DeconvolutedPeak`, `Feature<MZ, IonMobility>`,
  `ChargedFeature<Mass, IonMobility>`, `IonMobilityAwareCentroidPeak`,
  `IonMobilityAwareDeconvolutedPeak`) have doc comments on the impl. Each one is the
  authoritative statement of "these are the required arrays, and this is how each field
  maps to which array," and none of that is written in prose anywhere; a reader has to
  read every `arrays_required()`/`as_arrays()` pair by hand to reconstruct the contract.
- Every deconvoluted-peak conversion in this file (`DeconvolutedPeak`,
  `ChargedFeature<Mass, IonMobility>`, `IonMobilityAwareDeconvolutedPeak`, and the
  `From<&BinaryArrayMap> for DeconvolutedPeakSet` impl) computes neutral mass via
  `crate::utils::neutral_mass(mz, charge)`, which hardcodes a single-proton-per-charge
  adduct model (see section 5). None of the call sites mention this assumption, so a
  reader working from `conversion.rs` alone would not know that non-proton adducts are
  not representable through this path.

## 5. utils.rs (mzdata-bindata)

- `PROTON`, `mass_charge_ratio`, and `neutral_mass` have zero documentation of any kind,
  not even a one-line comment. These three items are the single most consequential piece
  of domain logic in the crate: every deconvoluted-peak and charged-feature array
  conversion in `conversion.rs` goes through them. Missing:
  - What `PROTON` represents (the proton mass in Daltons, used as the mass of the charge
    carrier).
  - That the formula assumes every charge is carried by exactly one proton adduct (`z as
    f64 * PROTON`). This is a real, consequential modeling assumption: it cannot
    represent charge carriers other than protons (e.g. sodium or potassium adducts, or
    electron-based negative-mode ionization models used by some deconvolution tools). A
    reader relying on this crate for negative-mode or adduct-aware work needs to know this
    up front, not discover it by reading the arithmetic.
  - The sign convention: for negative `z`, `neutral_mass` computes `(mz * |z|) - z *
    PROTON`, which, because `z` is negative, is a subtraction of a negative quantity
    (equivalent to `+ |z| * PROTON`). Whether this is the intended physical model for
    negative-mode spectra or an artifact of applying the same formula symmetrically is not
    stated anywhere and would benefit from an explicit note either way.

## 6. map.rs

- `BinaryArrayMap`'s type doc is a single line; its one public field, `byte_buffer_map`,
  has no field-level doc (worth stating explicitly that callers can insert/rearrange keys
  directly, and that a `HashMap` gives no iteration-order guarantee, which matters for
  anything relying on insertion or m/z order from `iter()`/`iter_mut()`).
- Most accessor methods (`encode_array`, `decode_all_arrays`, `add`, `get`, `has_array`,
  `mzs`, `intensities`, `charges`, `ion_mobility`, `stack_ion_mobility`) are reasonably
  documented, better than the crate average. The following stand out as needing more:
  - `search()`'s doc is one line ("Search for a specific m/z"). See section 9: this
    method's logic looks like it may not actually perform the tolerance search it is
    named for in the common case, which is exactly the kind of thing a fuller doc comment
    (spelling out preconditions and the algorithm) would have caught before it shipped.
  - `sort_from_indices()` has no doc comment. Its algorithm (a cycle-following in-place
    permutation apply, using `usize::MAX` as a "visited" tombstone) needs an explanation
    of what `mask` is required to represent (an index permutation; index semantics -
    "source index for position i" vs. the reverse - are not stated), and it silently
    skips any array in the map whose length does not equal the mask length, without any
    documentation of that skip.
  - `sort_by_array()`'s one-line doc does not mention: (a) it is a no-op if the sort key
    is already sorted (early return), (b) arrays of a different length than the sort key
    are left untouched (same silent skip as above), and (c) calling it with an `ASCII` or
    `Unknown`-typed sort key panics via `todo!()` (twice, once in this method and once in
    the `sort_from_indices` swap loop). "This panics on unsupported key dtypes" is
    important enough to be a `# Panics` section, not left to be discovered as a runtime
    `not yet implemented` message.
  - `mzs_mut`, `intensities_mut`, `charge_mut`, and `ion_mobility_mut` each silently
    coerce the underlying `DataArray`'s stored dtype (via `store_as`) as a side effect of
    producing a mutable typed slice. None of the four docs mention that requesting a
    mutable view can change how the array is stored on disk/in memory as a byproduct.
    Separately, and more seriously: `ion_mobility_mut` calls `store_as(BinaryDataArrayType
    ::Float32)` and then returns the coerced buffer as `&mut [f64]` (the function's
    declared return type). `mzs_mut`, right above it, does the analogous thing correctly
    (`store_as(BinaryDataArrayType::Float64)` before producing an `&mut [f64]`).
    Reinterpreting a buffer that was just downcast to 4-byte `Float32` storage as 8-byte
    `f64` elements would read pairs of adjacent floats as single nonsense doubles. This
    looks like a real bug (an `f32` vs. `f64` mismatch), not just a documentation gap, and
    is worth testing/verifying directly; flagging it here because writing the doc comment
    for this method is what surfaces the inconsistency.
- `BinaryArrayMap3D::unstack` and `::stack` are the best-documented nontrivial methods in
  the crate (both have `# Errors` sections). Remaining gaps:
  - `unstack`'s fallback behavior when there is no m/z array to sort by is only logged
    (`log::debug!("Unsorted unstack")`), not documented: the flattened output is emitted
    in "ion-mobility-slice-then-original-order" with no other ordering guarantee. A
    caller flattening an array map whose primary axis is not m/z would want this stated in
    the doc, not discovered via a debug log line.
  - `stack`'s binning of points into per-unique-ion-mobility-value slices uses exact `f64`
    equality (through `NonNaNF64`'s `Eq`/`Ord`, which is `total_cmp`-based, i.e. exact,
    bit-pattern-aware equality). Two ion mobility values that differ by even one ULP due
    to floating point noise will be treated as different slices, potentially producing far
    more (nearly-duplicate) slices than a user expects from an instrument whose IM axis is
    not bit-exact across a frame. This is an important, non-obvious behavior for anyone
    stacking data from a real acquisition and is not documented.
  - `NonNaNF64`'s `Hash` implementation (`((self.0 * 10000.0) as i64).hash(state)`) hashes
    at a fixed 4-decimal-digit precision while its `Eq`/`Ord` implementations use full
    `f64` precision. This is not unsound for `HashMap` use (colliding hashes for unequal
    values are always allowed), but it is unusual enough, and specific enough to this
    crate's ion mobility value ranges, that it deserves at least one explanatory comment
    for future maintainers, since it looks at first glance like an attempt to bucket
    "close" values together that does not actually do that (equality is still exact).

## 7. traits.rs

- `ByteArrayView`'s conversion methods `to_f32`, `to_f64`, `to_i32`, `to_i64` have zero
  doc comments despite being, almost certainly, the most frequently called methods in the
  entire crate (every m/z/intensity/charge/ion-mobility accessor in `map.rs` goes through
  one of these). None of the four document that:
  - Narrowing conversions are lossy (`Float64 -> f32`, any float `-> ` integer) and that
    the integer conversions truncate toward zero rather than round (inherited from
    `num_traits::AsPrimitive`'s semantics, never restated here).
  - They convert (reinterpreting **and** numerically casting between representations),
    as opposed to `iter_f32`/`iter_f64`/`iter_i32`/`iter_i64`/`iter_u8`/`iter_type`, which
    only reinterpret the raw bytes as the requested type without checking or converting
    from the array's actual stored `dtype`. Mixing these two families up (e.g. calling
    `iter_f32()` on data actually stored as `Float64`) will silently misinterpret bytes
    rather than converting them, and nothing in either family's documentation (there is
    none on the `iter_*` methods) warns of this.
- `coerce_from`, `coerce`, and `convert` (the generic machinery behind all four `to_*`
  methods) have no doc comments at all, despite containing the crate's central
  endian-handling logic (`#[cfg(target_endian = "big")]` byte-swapping) and Cow-based
  zero-copy-when-possible behavior. A short explanation of the contract (`T` must be
  `Pod` and must match the source dtype's width; behavior on a big-endian host) would
  help future maintainers who touch this code without having to re-derive it from the
  `cfg` attributes.
- `ByteArrayViewMut::coerce_from_mut` has no doc comment and no `# Safety` note, despite
  using `unsafe { slice::from_raw_parts_mut(...) }` to reinterpret a raw mutable byte
  buffer as `&mut [T]`. Notably, its bound is `T: Clone + Sized`, not `T: Pod` as the
  immutable `coerce_from` requires. That asymmetry means the compiler does not enforce
  that `T` has a valid bit pattern for arbitrary bytes or has compatible alignment; the
  soundness of every call site currently rests entirely on undocumented caller discipline
  (always calling this with a `T` matching the array's real `dtype`). This is worth an
  explicit `# Safety` section and, separately, worth asking whether the `Pod` bound should
  have been required here too.
- `DataSliceIter` (the type behind all the `iter_*` methods) has no type-level doc
  comment, and its `next_value` method (which contains the actual big-endian handling)
  has none either.

## 8. Cross-cutting: error handling and panics

Several of the findings above share a pattern worth naming once: a number of functions
whose signature returns `Result<_, ArrayRetrievalError>` actually panic instead of
returning `Err` for at least one input:

- `DataArray::decode`'s internal base64 decode step panics on malformed input.
- `DataArray::encode_bytestring` panics for unsupported/unimplemented compression targets,
  for `Decoded` as a target, and for numpress requested on the wrong source dtype.
- `BinaryArrayMap::sort_by_array` (and the shared `sort_from_indices` helper) panic via
  `todo!()` for `ASCII`/`Unknown`-typed sort keys.
- `BuildArrayMap3DFrom::as_arrays_3d`'s default body unwraps a `Result`.

None of these are documented with `# Panics` sections. Since the crate otherwise commits
to `Result`-based error propagation (and `ArrayRetrievalError` exists specifically to
carry these failures), each of these panic paths is a point where the documented contract
(the function signature) and the actual behavior diverge, which is worth calling out as a
group, not just as isolated missing doc comments.

---

## 9. Signal processing operations that need more detail

This section collects the operations that actually transform the numeric content of an
array (as opposed to just moving bytes around), since these encode algorithmic choices
that a caller cannot infer from the type signature.

1. **`BinaryCompressionType` variants, in general.** The enum documents itself as "the
   range of compression and encoding states," but gives no guidance on what each of the
   17 variants actually does to the data or when to choose it. In particular:
   - **Zlib**: lossless, general-purpose, no numeric-specific tricks. Baseline choice.
   - **NumpressLinear**: lossy. Applies linear prediction plus a per-array
     optimally-chosen fixed-point scaling factor (`numpress::optimal_scaling`), then
     encodes as truncated integers. Suited to smoothly-varying `Float64` data (only
     `Float64` is accepted; `encode_bytestring` panics otherwise). No mention anywhere
     of the resulting error bound, or that "optimal scaling" is a per-call search over
     the data rather than a fixed constant.
   - **NumpressSLOF**: lossy, "short logged float," a logarithmic fixed-point
     quantization better suited to data with large dynamic range (typically intensity
     arrays) than NumpressLinear's linear-domain quantization. Works on either `Float32`
     or `Float64` (unlike NumpressLinear). None of this distinction between the two
     Numpress variants, or when to prefer one over the other, is documented anywhere.
   - **NumpressPIC**, **NumpressPICZlib**, **NumpressPICZstd**: declared (they have CV
     accessions, `Display`/`from_accession`/`as_param` support, and appear in error
     messages), but **not implemented**. Neither `encode_bytestring` nor `decode` has a
     match arm for them; encoding falls through to a generic panic, decoding falls
     through to a generic `Err`. Nothing documents that these three variants exist only
     for round-tripping CV terms and cannot actually be used to compress or decompress
     data through this crate today.
   - **LinearPrediction**, **DeltaPrediction**: also declared but not wired into
     `encode_bytestring`/`decode` at all (same fate as the PIC variants for encoding; for
     decoding they fall into the generic "cannot decode" `Err`). Confusingly, the
     underlying algorithms these two variants are presumably meant to represent **do**
     exist as free functions (`linear_prediction_encoding`/`_decoding`,
     `delta_encoding`/`_decoding` in encodings.rs) and are public, but a caller has to
     know to call them directly on a decoded buffer; there is no path from
     `DataArray::store_compressed(BinaryCompressionType::LinearPrediction)` to them. This
     split (enum variant exists, free functions exist, but they are not connected) is a
     significant, undocumented gap between what the type system advertises as possible
     and what actually works.
   - **Zstd / ShuffleZstd / DeltaShuffleZstd / ZstdDict**: lossless. Each applies a
     different lossless pre-processing step (none, byte-plane shuffle, delta-then-shuffle,
     dictionary-encode) before generic zstd compression. No documentation anywhere
     compares these or gives guidance on which suits which kind of array (e.g. dictionary
     encoding suits low-cardinality integer arrays; shuffle/delta-shuffle suit smoothly
     varying floating point arrays like m/z; plain Zstd is the safe default for anything
     else).
   - Recommendation: at minimum, the `BinaryCompressionType` type-level doc should state
     which variants are actually implemented today, which are lossy vs. lossless, and
     point to a short comparison; ideally each variant gets its own one- or two-line doc.

2. **`linear_prediction_decoding` / `linear_prediction_encoding`** (encodings.rs). Beyond
   having no doc comment (section 2 above), `linear_prediction_decoding`'s body looks
   internally inconsistent and is worth verifying directly rather than trusting a
   docstring written after the fact: it first runs a `fold` over `values.iter_mut().
   skip(2)` that overwrites each element in place via `*current = tmp` (using a
   `(prev1, prev2)` pair threaded through the fold), and then immediately runs a *second*,
   separate `for i in 0..values.len()` loop that recomputes `values[i]` again from
   `values[i - 1]` and `values[i - 2]` -- but by this point those neighbors have already
   been overwritten by the first pass, so the second pass is reading already-transformed
   data rather than the original encoded values. There is no test in this crate's test
   module exercising `linear_prediction_encoding`/`_decoding` (the closest counterpart,
   `test_decode_delta_zstd`, only exercises `delta_encoding`/`_decoding`, a different pair
   of functions). Given both the missing documentation and the missing test coverage,
   this function's correctness cannot be confirmed by reading alone and should be
   verified with a round-trip test before anyone relies on `BinaryCompressionType::
   LinearPrediction`-flavored data (or the raw functions) elsewhere.

3. **`BinaryArrayMap::search`** (map.rs). Documented only as "Search for a specific m/z."
   Reading the implementation: it computes the tolerance window's lower bound, then calls
   `binary_search_by` comparing array elements against that lower bound, and **only
   proceeds to the tolerance-checking linear scan if `binary_search_by` returns `Ok`**
   (i.e., some array element is exactly, bit-for-bit, equal to the lower bound). If it
   returns `Err` (the insertion-point case, which is what `binary_search_by` returns for
   essentially every real-world query against a continuous m/z array, since floats
   essentially never land exactly on the computed tolerance boundary), the function
   returns `None` immediately without ever inspecting the actual candidate region. By
   contrast, `BinaryArrayMap3D::search_ion_mobility` a little further down in the same
   file correctly handles both the `Ok(i)` and `Err(i)` arms of essentially the same
   binary-search pattern. This asymmetry strongly suggests `search` does not perform the
   tolerance-based nearest-neighbor search its (thin) documentation promises, and instead
   returns `None` for nearly all queries in practice. This is worth a direct test
   (`search` against a known array with a query that should match but does not fall
   exactly on a tolerance boundary) before relying on it; if confirmed, it is a
   correctness bug, not just a documentation gap, but it was found here specifically
   because the one-line doc comment invited checking what "search for a specific m/z"
   actually guarantees.

4. **`BinaryArrayMap::ion_mobility_mut`** (map.rs). As described in section 6: this method
   stores the target array as `BinaryDataArrayType::Float32` and then returns it coerced
   as `&mut [f64]`, whereas its sibling `mzs_mut` correctly stores as `Float64` before
   returning `&mut [f64]`. If this is not intentional, it means every mutable ion-mobility
   accessor silently corrupts the array's numeric content on first use (each pair of
   4-byte floats reinterpreted as one 8-byte double). Recommend a direct test of
   `ion_mobility_mut` (write a known value, read it back via `ion_mobility()`) before
   relying on this path in any signal-processing code.

5. **Numpress helper functions** (`DataArray::compress_numpress_linear`,
   `compress_numpress_slof`, `decompress_numpress_linear`, `decompress_numpress_slof`;
   array.rs). All four are undocumented (section 3). Beyond restating what each algorithm
   does (covered in item 1 above), the docs should state: `compress_numpress_linear` only
   accepts `Float64` input (enforced by its caller in `encode_bytestring`, not by the
   function's own signature, since it takes `&[f64]` directly, so this constraint is only
   visible at the call site); `compress_numpress_slof` is generic over
   `numpress::AsFloat64`, accepting either `f32` or `f64`, unlike its Linear counterpart,
   with no note anywhere explaining why the two Numpress modes differ in which source
   types they accept.

6. **Delta and shuffle pre-processing for Zstd** (`compress_delta_zstd`,
   `compress_dict_zstd`; array.rs). Both apply a lossless numeric transform
   (delta-encoding, or dictionary/categorical encoding) before generic compression, and
   neither documents what kind of array benefits from which transform. As a concrete,
   fillable gap: `compress_delta_zstd` is a good fit for arrays that are monotonic or
   smoothly varying (m/z arrays are a natural fit, since deltas between adjacent sorted
   m/z values are small and repetitive); `compress_dict_zstd` is a good fit for
   low-cardinality arrays (charge state arrays, or any small integer-coded categorical
   array) and a poor fit for high-cardinality continuous data (the value table itself
   would approach the size of the source data). None of this guidance exists in the crate
   today, and a caller has no way to know which compression scheme suits their array type
   without reading (and understanding) the codec implementations directly.

7. **`mass_charge_ratio` / `neutral_mass`** (utils.rs). As covered in section 5, these
   are the crate's core mass-to-m/z conversion functions and have zero documentation.
   Because every deconvoluted-peak and charged-feature conversion in `conversion.rs`
   depends on them, their single-proton-adduct assumption is effectively baked into the
   whole crate's notion of "neutral mass," and this needs to be stated explicitly
   somewhere a reader would find it (ideally on the functions themselves, and
   cross-referenced from every `impl BuildFromArrayMap` that calls them).

8. **`Feature<MZ, IonMobility>` / `ChargedFeature<Mass, IonMobility>` array round-trip**
   (conversion.rs, `as_arrays`/`try_from_arrays`/`try_from_arrays_3d`). This is the most
   algorithmically involved undocumented code in the crate: flattening multiple ion
   mobility features into parallel arrays, sorting all points globally by
   `(m/z, ion_mobility, feature_index)`, and recovering the grouping on the way back in
   via a synthetic, crate-invented `"feature identifier array"`. This deserves a
   dedicated paragraph (ideally on the `BuildArrayMapFrom`/`BuildFromArrayMap` trait impls
   themselves) describing the wire convention, since nothing else in the public API
   documents that this non-standard array even exists, let alone its meaning or its
   contiguous-zero-based-index precondition.

---

## Suggested priority order

1. Verify (write a focused test for) `BinaryArrayMap::search` (section 9, item 3) and
   `BinaryArrayMap::ion_mobility_mut` (section 9, item 4). Both look like they may not do
   what their names/doc comments claim, and both sit on a common signal-processing path
   (finding a data point near a query value; mutating ion mobility data in place).
2. Verify `linear_prediction_decoding`/`linear_prediction_encoding` (section 9, item 2)
   with a round-trip test; currently unverified and undocumented.
3. Document, at the `BinaryCompressionType` type level, which variants are actually
   implemented today versus merely declared for CV round-tripping (section 9, item 1);
   this prevents a caller from being surprised by a panic or a decode error for
   `NumpressPIC*`, `LinearPrediction`, or `DeltaPrediction`.
4. Add `# Panics` sections to `DataArray::encode_bytestring`, `DataArray::decode`, and
   `BinaryArrayMap::sort_by_array`/`sort_from_indices` (section 8).
5. Document `mass_charge_ratio`/`neutral_mass` and their single-proton-adduct assumption
   (section 5), since this quietly limits what the whole crate can represent.
6. Document the `Feature`/`ChargedFeature` array round-trip convention (section 9, item 8)
   and the general precedence/"required arrays" contract on `BuildFromArrayMap`/
   `BuildArrayMapFrom` (section 4).
7. Clarify the `as_bytes` vs. `to_bytes` and `wrap()`'s "native" vs. "little-endian" byte
   order language (section 2 and section 3), since both are correctness-adjacent for
   big-endian hosts.
8. Document `slice`/`slice_buffer`'s byte-offset (not element-index) semantics
   (section 3).
9. Fill in the remaining missing doc comments (constructors in array.rs, the `to_*`/
   `iter_*` families in traits.rs, `ArraysAvailable`, per-variant docs on `ArrayType`/
   `BinaryCompressionType`).
