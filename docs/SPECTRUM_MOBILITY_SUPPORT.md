# MSSpectrum ion mobility support

Source: `src/openms/include/OpenMS/KERNEL/MSSpectrum.h` (784 lines) and
`src/openms/source/KERNEL/MSSpectrum.cpp` (968 lines) at Core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`, together with
`src/openms/source/IONMOBILITY/IMDataArrayUtils.cpp`, which
`MSSpectrum::containsIMData` and `getIMData` delegate the array-name
interpretation to, and `src/openms/source/IONMOBILITY/IMTypes.cpp`, which owns
the `DriftTimeUnit` names and the `-1` unset sentinel.

Rust: `src/kernel/spectrum_mobility.rs`. Tests: `tests/spectrum_mobility.rs`.
Manifest: `tests/data/spectrum_mobility_provenance.json`.

This work package closes the ion-mobility half of the header. The value type
itself, its container and search operations, `getType`, the record metadata and
the range manager were ported by earlier waves and live in `src/kernel.rs`,
`src/kernel/spectrum_type.rs`, `src/kernel/acquisition_fields.rs` and
`src/kernel/ranges.rs`; the table below covers the whole header so that no
member is unaccounted for, and names the file that owns each row.

## API mapping

### Nested types

| C++ member | Rust counterpart |
| --- | --- |
| `struct RTLess` | not ported as a type: `a.rt < b.rt`, or `sort_by(\|a, b\| a.rt.total_cmp(&b.rt))`. The source comparator is a function object only because C++ needs one. |
| `struct IMLess` | not ported as a type: `a.drift_time < b.drift_time`. It has no class-test section upstream either. |
| `struct Chunk` (`start`, `end`, `is_sorted`, ctor) | `spectrum_mobility::Chunk` with the same three public fields and `Chunk::new` |
| `struct Chunks` | `spectrum_mobility::Chunks` |
| `Chunks::Chunks(const MSSpectrum&)` | `Chunks::new()`; the spectrum is passed to `add` instead — see *Native differences* |
| `Chunks::add(bool)` | `Chunks::add(&MSSpectrum, bool) -> Result<()>` |
| `Chunks::getChunks()` | `Chunks::chunks() -> &[Chunk]`, plus native `len` / `is_empty` |
| `enum class RasterAggregation { SUM, MAX }` | `spectrum_mobility::RasterAggregation { Sum, Max }`, `Sum` is the `Default` as in the source's default argument |
| `PeakType`, `CoordinateType`, `ContainerType` | `Peak1D`, `f64`, `Vec<Peak1D>` (`MSSpectrum::peaks`) — `src/kernel.rs` |
| `RangeManagerContainerType`, `RangeManagerType` | `kernel::ranges::RangeManager`, via `MSSpectrum::range_manager()` — `src/kernel/ranges.rs` |
| `FloatDataArray(s)`, `StringDataArray(s)`, `IntegerDataArray(s)` | `DataArray<f32>` / `<String>` / `<i32>` and `Vec<_>` of each — `src/kernel.rs` |
| `Iterator`, `ConstIterator`, `ReverseIterator`, `ConstReverseIterator` and the `using typename ContainerType::…` aliases | not ported: `&mut [Peak1D]` / `&[Peak1D]` and their standard iterators carry all of it |

### Exported `std::vector<Peak1D>` surface

| C++ member | Rust counterpart |
| --- | --- |
| `operator[]`, `front`, `back`, `data` | index and slice methods on `MSSpectrum::peaks` — `src/kernel.rs` |
| `begin`, `end`, `cbegin`, `cend`, `rbegin`, `rend` | `peaks.iter()`, `.iter_mut()`, `.rev()` |
| `size`, `empty` | `MSSpectrum::len`, `MSSpectrum::is_empty` |
| `resize`, `reserve`, `shrink_to_fit`, `push_back`, `emplace_back`, `pop_back`, `insert`, `erase`, `swap` | the same operations on `MSSpectrum::peaks`; the field is public, so no wrapper is added |

### Construction, assignment and equality

| C++ member | Rust counterpart |
| --- | --- |
| `MSSpectrum()` | `MSSpectrum::new()` / `Default` — `src/kernel.rs` |
| `MSSpectrum(const std::initializer_list<Peak1D>&)` | `MSSpectrum::from_peaks(Vec<Peak1D>)`, `From<Vec<Peak1D>>` |
| copy constructor, copy assignment | `Clone` |
| move constructor, move assignment | Rust move semantics; `std::mem::take` where the source needs a named moved-from value |
| `~MSSpectrum()` | compiler-generated `Drop` |
| `operator=(const SpectrumSettings&)` | not ported: the port has no `SpectrumSettings` base; the settings are named fields, assigned individually. `kernel::spectrum_helper::copy_spectrum_meta` is the closest operation. |
| `operator==`, `operator!=` | derived `PartialEq` — see *Native differences* for the `name` field |

### Scalar metadata

| C++ member | Rust counterpart |
| --- | --- |
| `getRT`, `setRT` | public field `MSSpectrum::rt`, seconds, `-1` unset — `src/kernel.rs` |
| `getDriftTime` | public field `MSSpectrum::drift_time` (raw, `-1` unset) plus `MSSpectrum::drift_time_if_set() -> Option<f64>` and `has_drift_time()` |
| `setDriftTime` | `MSSpectrum::set_drift_time(Option<f64>) -> Result<()>`, or the field directly |
| `getDriftTimeUnit` | public field `MSSpectrum::drift_time_unit`. No `drift_time_unit()` method is added: a method with a field's name makes every existing `[`MSSpectrum::drift_time_unit`]` doc link ambiguous, and `cargo doc -D warnings` then fails. |
| `setDriftTimeUnit` | the same public field |
| `getDriftTimeUnitAsString` | `MSSpectrum::drift_time_unit_as_string() -> &'static str` |
| `getMSLevel`, `setMSLevel` | public field `MSSpectrum::ms_level` — `src/kernel.rs` |
| `getName`, `setName` | public field `MSSpectrum::name` — `src/kernel.rs` |

### Data arrays

| C++ member | Rust counterpart |
| --- | --- |
| `getFloatDataArrays` (const and mutable), `setFloatDataArrays` | public field `MSSpectrum::float_data_arrays` — `src/kernel.rs` |
| `getStringDataArrays` (const and mutable), `setStringDataArrays` | public field `MSSpectrum::string_data_arrays` |
| `getIntegerDataArrays` (const and mutable), `setIntegerDataArrays` | public field `MSSpectrum::integer_data_arrays` |

### Sorting

| C++ member | Rust counterpart |
| --- | --- |
| `sortByIntensity(bool reverse)` | `MSSpectrum::sort_by_intensity(bool) -> Result<()>` — `src/kernel.rs` |
| `sortByPosition()` | `MSSpectrum::sort_by_position() -> Result<()>` — `src/kernel.rs` |
| `sortByIonMobility()` | `MSSpectrum::sort_by_ion_mobility() -> Result<()>` |
| `sortByPositionPresorted(const std::vector<Chunk>&)` | `MSSpectrum::sort_by_position_presorted(&[Chunk]) -> Result<()>` |
| `isSorted()` | `MSSpectrum::is_sorted() -> bool` — `src/kernel.rs` |
| `isSortedByIM()` | `MSSpectrum::is_sorted_by_im() -> Result<bool>` |
| `template<Predicate> isSorted(const Predicate&)` | not ported: the caller writes the comparison over its own indices. The source template exists to smuggle an index predicate through `std::is_sorted`, which Rust does not need. |
| `template<Predicate> sort(const Predicate&)` | not ported as a generic: the caller computes a permutation and applies it with `MSSpectrum::select`, which performs the same up-front data-array size check before any peak moves. `sort_by_ion_mobility` and `sort_by_position_presorted` do exactly that internally. |

### Searching

| C++ member | Rust counterpart |
| --- | --- |
| `findNearest(mz)` | `MSSpectrum::find_nearest` — `src/kernel.rs` |
| `findNearest(mz, tolerance)` | `MSSpectrum::find_nearest_with_tolerance` |
| `findNearest(mz, left, right)` | `MSSpectrum::find_nearest_in_window` |
| `findHighestInWindow(mz, left, right)` | `MSSpectrum::find_highest_in_window` |
| `MZBegin(mz)` / `PosBegin(mz)`, const and mutable | `MSSpectrum::mz_begin(f64) -> Result<usize>` |
| `MZEnd(mz)` / `PosEnd(mz)`, const and mutable | `MSSpectrum::mz_end(f64) -> Result<usize>` |
| the eight `(begin, mz, end)` subrange overloads | not ported: slice the peaks (`&peaks[a..b]`) and use `partition_point`; an index is a position in the whole list, so no iterator-pair overload is needed |
| `containsIMData()` | `MSSpectrum::contains_im_data() -> bool` |
| `getIMData()` | `MSSpectrum::im_data() -> Result<(usize, DriftTimeUnit)>` |
| `maybeGetIMData()` | `MSSpectrum::maybe_im_data() -> Option<(DriftTimeUnit, &[f32])>` |

### Remaining operations

| C++ member | Rust counterpart |
| --- | --- |
| `updateRanges()` | `MSSpectrum::range_manager() -> Result<RangeManager>`, recomputed on demand — `src/kernel/ranges.rs` |
| `clear(bool clear_meta_data)` | `MSSpectrum::clear(bool)` — `src/kernel.rs` |
| `select(const std::vector<Size>&)` | `MSSpectrum::select(&[usize]) -> Result<()>` — `src/kernel.rs` |
| `selectUnchecked(const std::vector<Size>&)` | not ported: `select` is the only entry point. Its checks are a bounds scan and a duplicate scan over a `vec![bool]`, which is the same order as the reorder itself, and skipping them is what makes the source overload undefined for duplicate indices. |
| `getType(const bool query_data)` | `MSSpectrum::get_type(bool) -> Result<SpectrumType>` — `src/kernel/spectrum_type.rs` |
| `using SpectrumSettings::getType` | public field `MSSpectrum::spectrum_type` |
| `getBasePeak()` (const and mutable) | `MSSpectrum::base_peak() -> Option<&Peak1D>` — `src/kernel.rs`; the mutable overload is an index into the public `peaks` field |
| `calculateTIC()` | `MSSpectrum::calculate_tic() -> f32` — `src/kernel.rs` |
| `rasterizeIMFrame(float*, …)` | `MSSpectrum::rasterize_im_frame(&ImFrameRaster) -> Result<Vec<f32>>` |
| `checkDataArraySizes_()` (protected) | `MSSpectrum::validate_data_arrays()` — `src/kernel.rs` |
| protected fields `retention_time_`, `drift_time_`, `drift_time_unit_`, `ms_level_`, `name_`, the three array vectors | the public fields of the same meaning — `src/kernel.rs` |
| `operator<<(std::ostream&, const MSSpectrum&)` | not ported for the whole spectrum. `Peak1D: Display` reproduces the per-peak `POS: … INT: …` lines the source prints; the `-- MSSPECTRUM BEGIN --` framing and the `SpectrumSettings` block have no counterpart. |
| `setIMFormat` (inherited from `SpectrumSettings`) | `metadata::SpectrumSettings::set_im_format` — `src/metadata/acquisition.rs`. Added by the IM-types package, which ported the quartet on the class that declares it; see [IM_TYPES_SUPPORT.md](IM_TYPES_SUPPORT.md) |
| `getIMFormat` (inherited from `SpectrumSettings`) | `metadata::SpectrumSettings::im_format` — `src/metadata/acquisition.rs`, over the public `ion_mobility_format` field |
| `setIMPeakType` (inherited from `SpectrumSettings`) | `metadata::SpectrumSettings::set_im_peak_type` — `src/metadata/acquisition.rs` |
| `getIMPeakType` (inherited from `SpectrumSettings`) | `metadata::SpectrumSettings::im_peak_type` — `src/metadata/acquisition.rs`, over the public `ion_mobility_peak_type` field |
| calling the four above **on an `MSSpectrum`** | still not available: the Rust `MSSpectrum` flattens `SpectrumSettings` into named fields rather than inheriting it, and it has no `ion_mobility_format` / `ion_mobility_peak_type` field of its own. `ImTypes::determine_im_format_with_stored` takes the stored value as an argument so a caller can supply it. This is the one part of this package's deferral that survives; see *Deferrals* |

## Preserved source conventions

- **The `-1` unset drift time.** `MSSpectrum::drift_time` keeps the source's
  `IMTypes::DRIFTTIME_NOT_SET` representation, so a spectrum round-trips through
  mzML unchanged. `drift_time_if_set()` is a view over it, not a second storage.
- **Name-only IM detection.** `contains_im_data` inspects float-array *names*
  and nothing else, exactly as `containsIMData` does. An array with the right
  name and no entries still marks the spectrum as an IM frame; the operations
  that read values check the length themselves.
- **First match wins.** A spectrum with two ion-mobility arrays uses the earlier
  one, because the source scans in storage order and returns on the first hit.
- **The `DriftTimeUnit` name table.** `"<NONE>"`, `"ms"`, `"1/K0"`, `"FAIMS_CV"`,
  `"CCS"` — the source's `NamesOfDriftTimeUnit`, in its order.
- **The unit rules of `IMDataArrayUtils::getIMUnit`,** including their
  precedence: exact PSI-MS term, then the two inverse-reduced `UserParam`
  prefixes as `1/K0`, then the generic `"Ion Mobility"` prefix whose embedded
  accession chooses `1/K0` (`MS:1002815`, `MS:1003006`), CCS (`MS:1002954`) or
  milliseconds. Prefix matching is case sensitive, as `StringUtils::hasPrefix` is.
- **A sorted IM array short-circuits.** `sort_by_ion_mobility` returns without
  touching anything when the array is already ordered, as the source's
  `std::is_sorted` guard does.
- **Stable ordering.** Both sorts are stable, so equal keys keep their input
  order, matching `std::stable_sort` and `std::inplace_merge`.
- **Chunk semantics.** `[start, end)` half-open, `Chunks::add` closing a run at
  the spectrum's current length, and an empty chunk list sorting nothing.
- **Raster geometry.** Row-major with m/z slowest, `mz_bin * im_bins + im_bin`;
  peaks outside `[min, max]` on either axis skipped; a value exactly at the
  maximum clamped into the last bin; `SUM` the default aggregation. The frame is
  not required to be sorted by m/z, and all peaks are visited linearly.
- **An empty IM frame rasterizes to zeros,** without a length check — the source
  returns early for an empty spectrum before comparing the IM array's length
  against the peak count.

## Native differences

- **`Chunks` does not hold the spectrum.** The source stores a
  `const MSSpectrum&` for the builder's lifetime and reads `size()` inside
  `add()`, i.e. it observes a spectrum that is being mutated through another
  path. Rust cannot express that, so the spectrum is an argument to `add`. The
  recorded chunks are identical; the class-test loop ports line for line.
- **`sort_by_position_presorted` requires the chunks to tile `[0, len)`** and
  rejects a run that claims to be sorted but is not. Neither is checked
  upstream, and both are observable there: a chunk list that stops short of the
  peak count leaves the tail unsorted when the spectrum has data arrays but
  sorts the whole list when it has none, because the two branches treat `chunks`
  differently; a dishonest `is_sorted` feeds `std::inplace_merge` an unsorted
  range, whose result is unspecified. Checking costs one linear scan, which the
  merge pays anyway. Because the port validates, its two branches would produce
  the same permutation, so it keeps only the chunk-aware one.
- **The merge is bottom-up rather than a balanced recursion.** The source merges
  adjacent runs through a recursive `std::inplace_merge`; the port merges
  adjacent runs pairwise until one remains. A stable merge of adjacent runs is
  associative in its result, so the permutation is the same, and the recursion
  depth is gone.
- **`sort_by_ion_mobility` and `is_sorted_by_im` require one ion-mobility value
  per peak.** The source relates the array to nothing: `sortByIonMobility` tests
  `std::is_sorted` on it and an empty range is sorted, so a spectrum whose
  ion-mobility array has no entries — which `checkDataArraySizes_` explicitly
  permits — is left in whatever order it was in, silently, and `isSortedByIM`
  reports `true` for it. The port refuses the mismatch.
- **Both refuse a non-finite ion-mobility array** rather than reporting it
  sorted. `std::is_sorted` uses `<`, under which a NaN neighbour compares
  unordered and the array is declared sorted. This is the one place the port
  differs from `MSSpectrum::is_sorted`, which folds an invalid coordinate into
  `false`.
- **`maybe_im_data` returns `Option<(DriftTimeUnit, &[f32])>`.** The source's
  `{DriftTimeUnit::NONE, {}}` conflates three states: no ion-mobility array, an
  ion-mobility array with no entries, and an ion-mobility array whose CV term
  declares no usable unit. The `Option` separates all three, and borrowing the
  values copies nothing where the source copies the whole array.
- **`rasterize_im_frame` returns an owned `Vec<f32>`.** The `float*` output
  parameter, the `Exception::NullPointer` that guards it and the caller's
  obligation to pre-allocate `im_bins * mz_bins` all disappear. Nothing is
  allocated before every check has passed, so a rejected call cannot clear an
  image the way the source's ordering does.
- **Non-finite values are rejected** by the rasterizer and both sorts. The
  source has no guard; a NaN m/z passes both range tests and
  `static_cast<Int64>(NaN)` is undefined behaviour in C++.
- **`set_drift_time` validates.** The source assigns whatever it is handed, so a
  NaN drift time reaches the range manager and every later comparison.
- **The pinned CV table.** `IMDataArrayUtils::getIMUnit` resolves array names
  against the PSI-MS controlled vocabulary at run time. The port pins the nine
  children of `MS:1002893` and their units, read out of `resources/cv/psi-ms.obo`,
  so this module needs neither the controlled-vocabulary feature nor a CV file,
  and `tests/spectrum_mobility.rs` runs under `--no-default-features`. At the
  pinned CV, three children declare `MS:1002814` (`1/K0`) and six declare
  `UO:0000028` (milliseconds); none declares `UO:0000324`, so the CV can never
  produce `CCS` and only the `"Ion Mobility (MS:1002954)"` fallback can.
- **`operator==` compares the name.** The source deliberately excludes `name_`
  (`MSSpectrum.cpp:516`) and its class-test section asserts that two spectra
  differing only in name are equal. `MSSpectrum` derives `PartialEq` over every
  field, so they are not equal here. `tests/spectrum_mobility.rs` asserts the
  Rust behaviour and names the divergence at the assertion. Changing it would
  mean hand-writing `PartialEq` on a struct this work package may not touch, and
  silently ignoring a field is worse than comparing it.
- **Serial.** Neither sort nor the rasterizer carries `#pragma omp` upstream, so
  this module has no parallelism gap to declare; the crate is serial throughout.

### Agreement with `src/kernel/ranges.rs`

`src/kernel/ranges.rs` contains a private `ion_mobility_array_index`, used to
decide whether a spectrum's mobility range comes from a float array or from the
scalar drift time. It applies the same rule as `contains_im_data`: the nine
pinned PSI-MS children of `MS:1002893`, plus the `"mean inverse reduced ion
mobility array"`, `"inverse reduced ion mobility"` and `"Ion Mobility"` prefixes.
**The two must stay in step.** If they diverge, `range_manager()` and
`contains_im_data()` disagree about whether a spectrum is an IM frame, and a
mobility range silently comes from the wrong place. The public entry point lives
here to avoid a method-name collision between the two work packages; the private
helper is the one to change if the rule ever moves, and this table is the record
of that dependency. `tests/spectrum_mobility.rs` asserts every recognized name
and unit, and `update_ranges_includes_mobility_from_the_im_array` pins the two
together on the class test's own fixture.

## Checked boundaries and evidence

| Boundary | Value | Where |
| --- | --- | --- |
| `MSSpectrum::MAX_MOBILITY_ITEMS` | 100 000 000 peaks, chunks or ion-mobility values | `sort_by_ion_mobility`, `sort_by_position_presorted`, `rasterize_im_frame`, `Chunks::add` |
| `MSSpectrum::MAX_RASTER_PIXELS` | 16 777 216 pixels (64 MiB of `f32`) | `ImFrameRaster::validate` |
| pixel-count overflow | `checked_mul` | `ImFrameRaster::pixels` |

`MAX_MOBILITY_ITEMS` matches `RangeManager::MAX_ITEMS`, so a spectrum whose
ranges can be computed can also be sorted by ion mobility.

Every preflight runs before anything is allocated or mutated, and every mutating
operation commits through `MSSpectrum::select`, which validates the index
permutation and the data-array lengths before it touches storage. A rejected
call therefore leaves the spectrum byte-identical; `tests/spectrum_mobility.rs`
asserts that for every error path.

**Evidence tiers.** Every expectation carrying a section number is tier 3
(source review): the literals of the 71 `START_SECTION`s of
`MSSpectrum_test.cpp` transcribed into `tests/spectrum_mobility.rs`, including
`getPrefilledSpec`'s ion-mobility permutation and the five-chunk presorted
fixture. `rasterizeIMFrame` has no upstream section, so its bin indices,
aggregation and max-bound clamp are tier 4 (independently derived from
`MSSpectrum.cpp:857-967`), as are the native boundaries. No C++ was built or
executed. `tests/data/spectrum_mobility_provenance.json` hashes all eight source
files and pins 21 source anchors.

## Section accounting

`MSSpectrum_test.cpp` has 71 `START_SECTION`s. 68 are ported into
`tests/spectrum_mobility.rs`, none of them merely mapped to another file's test.
The remaining three had no `MSSpectrum` surface to assert when this package was
written, because the port had no field for either property; this file maps them
onto the ported enums, each with one concrete value from the section, in
`im_format_and_peak_type_enums_exist_but_no_spectrum_field_does`:

| Section | Asserts | Mapped to here | Now ported in |
| --- | --- | --- | --- |
| `void setIMFormat(IMFormat)` | 2 | `IonMobilityFormat::PerPeak.name() == "im_peak"` and `IonMobilityFormat::None.name() == "none"` | `tests/im_types.rs::set_im_format_round_trips`, on `SpectrumSettings` |
| `IMPeakType getIMPeakType() const` | 1 | `IonMobilityPeakType::default() == Unknown` | `tests/im_types.rs::get_im_peak_type_defaults_to_unknown` |
| `void setIMPeakType(IMPeakType)` | 2 | `IonMobilityPeakType::Centroid.name() == "im_centroided"` and `Profile.name() == "im_profile"` | `tests/im_types.rs::set_im_peak_type_round_trips` |

None is unaccounted for. The three enum assertions in this file are kept as they
are: they still hold, and they pin the enum names that the round-trip tests in
`tests/im_types.rs` rely on. The accounting for those three sections now belongs
to [IM_TYPES_SUPPORT.md](IM_TYPES_SUPPORT.md), which asserts them on
`SpectrumSettings`, the class that declares the members.

Two further sections assert C++-only mechanics and are ported to their nearest
Rust meaning, with the substitution stated at the assertion: `~MSSpectrum()`
(compiler-generated `Drop`) and the `noexcept` move-constructor check of
`MSSpectrum(const MSSpectrum&&)` (Rust moves are unconditional).

## Deferrals

- **Partly closed.** `MSSpectrum` needs an `im_format` and an `im_peak_type`
  field to complete `SpectrumSettings`' ion-mobility surface. Adding them means
  editing the struct in `src/kernel.rs`, which this work package may not do.
  The four accessors themselves are ported — on `SpectrumSettings`, the class
  that declares them — by the IM-types package
  ([IM_TYPES_SUPPORT.md](IM_TYPES_SUPPORT.md)); what survives is only that an
  `MSSpectrum` cannot answer them itself.
- **Closed.** `tests/data/spectrum_mobility_provenance.json` was not yet listed
  in `SOURCE_PROVENANCE.json`, so `tools/check_core_sdk.py` did not verify its
  hashes. It is now registered under both `current_sdk_reference_manifests` and
  `kernel_wave2_reference_manifests`.
- `MSSpectrum::operator==` excluding `name_` is left as a stated divergence
  rather than reproduced; reproducing it needs a hand-written `PartialEq` on a
  struct this work package does not own.
