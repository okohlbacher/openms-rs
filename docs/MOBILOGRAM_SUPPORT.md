# Mobility peaks and mobilograms

This native group is based on OpenMS Core SDK commit
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. It requires no new dependency.
The implementation is [mobilogram.rs](../src/kernel/mobilogram.rs); source
fixtures and independent boundaries are in [mobilogram.rs](../tests/mobilogram.rs).
The exact source hashes are recorded in
[mobilogram_provenance.json](../tests/data/mobilogram_provenance.json).

## Value and container mapping

`MobilityPeak1D` has public `mobility: f64` and `intensity: f32`, with zero
defaults, `DIMENSION = 1`, a constructor, `Copy`, `Clone`, floating `PartialEq`,
`Display` and native `Hash`. The coordinate is in the mobilogram's unit; no unit
conversion is inferred. Source position/mobility aliases refer to this one
scalar. Source comparator overloads become ordinary scalar comparisons of
these public fields. There is no `Eq` or invented total floating ordering.

Display retains `POS: … INT: …`. Rust's numeric formatting is used, including an
optional formatter precision; C++ stream locale/precision state is not emulated.
Hash covers both values and normalizes signed zero. Hash values are not promised
stable across processes, Rust versions or languages.

`Mobilogram` owns `peaks`, `rt` (default -1 seconds), `drift_time_unit` (default
`None`), and float, integer and string `DataArray` vectors. The existing
`DriftTimeUnit` covers every source unit and its literal name. Public fields and
the standard `Vec` API replace mechanical accessors, iteration, indexing,
reserve/resize, insertion, erasure and move/copy operations. RT comparisons use
the public `rt` field. There is no source mobilogram-level MetaInfo inheritance.

| Source operation | Native operation |
| --- | --- |
| `updateRanges`, mobility/intensity extrema | `ranges()` returns current optional `NumericRange` values |
| `isSorted`, custom predicate overload | `is_sorted`, `is_sorted_by` |
| `sortByPosition`, `sortByIntensity` | `sort_by_position`, `sort_by_intensity(reverse)` |
| stable predicate `sort` | `sort_by`, callback receives the original container and two indices |
| `MBBegin` / `PosBegin` and iterator subranges | `mobility_begin`, `mobility_begin_in` |
| `MBEnd` / `PosEnd` and iterator subranges | `mobility_end`, `mobility_end_in` |
| three `findNearest` overloads | `find_nearest`, `find_nearest_with_tolerance`, `find_nearest_in_window` |
| `findHighestInWindow` | `find_highest_in_window` |
| mutable / const `getBasePeak` | `base_peak_mut`, `base_peak`, `base_peak_index` |
| `calculateTIC` | `calculate_tic` |
| `select`, `selectUnchecked` | checked `select`; no unchecked duplicate-index API |
| `clear` | `clear` removes peaks and all arrays, retaining RT/unit |
| `operator==` | `source_equal` ignores arrays, exactly as source |
| source partial `swap` | `swap_peak_data` exchanges peaks/RT/unit only |
| full native value equality / swap | `PartialEq` includes arrays; `std::mem::swap` exchanges everything |

Source `swap` can leave arrays misaligned; `swap_peak_data` deliberately exposes
that behavior under an explicit name. Sorting/selection subsequently reject
misaligned arrays when they are consumed. Source position/intensity sort has an
already-sorted early return, so that no-op still ignores unused malformed arrays.
Custom `sort_by` checks lengths before invoking the first callback. Empty arrays
remain placeholders; duplicates and out-of-bounds selection indices are errors.

Ranges are recomputed on demand, consistent with existing native spectra.
They cannot be stale after public mutation or subset selection. Source cached
range state, direct cache manipulation, generic inherited RangeManager algebra
and its debug warning are not represented by this group. This is not a claim of
a complete general RangeManager implementation, mobility-bearing MSExperiment,
mobility peak picking, or mobilogram file transport.

## Arithmetic and checked boundaries

Sorts are stable, with all nonempty array values moving in the same permutation.
Array names, metadata, processing handles and string buffers retain ownership;
sorting does not deep-clone their payload. Predicate callbacks see original
indices throughout sorting and must supply a strict weak ordering. A panic in
caller code has ordinary Rust unwinding semantics; mutations start only after
all comparisons and fallible preflight have completed.

Nearest search preserves source lower-bound and tie behavior: exact queries
select the first duplicate, midpoint ties select the predecessor, and queries
above the final coordinate select the final duplicate. Tolerance borders are
inclusive. Symmetric and asymmetric searches retain their separate source
branches, including defined finite negative-tolerance behavior. Asymmetric
search checks only the bound on the chosen side before optionally trying the
immediate opposite neighbor; it is not rewritten as a generic interval filter.
Highest-in-window rejects a reversed interval rather than constructing an
invalid source iterator range. Empty nearest returns `None` (the source
no-tolerance overload throws); other empty searches also return `None`.

Searches validate the requested coordinate subrange before binary search.
Other operations inspect only the fields they use: TIC/base-peak need intensity,
position sorting needs mobility, and custom sorting checks array shapes without
inspecting metadata. Finite signed values are permitted. Nonfinite consumed
scalars, invalid subranges, invalid selection indices and nonfinite TIC results
are checked errors. Bounds computed from finite extreme tolerances may overflow
to infinity; their source comparison behavior is retained. TIC accumulates in
storage order as f32 after every addition, including source rounding of small
increments following a large intensity. Maximum-intensity ties select the first
peak.

`MobilogramLimits` defaults are 10,000,000 peaks, 100,000 consumed arrays,
50,000,000 visits/comparisons and 256 MiB cumulative temporary vector storage.
Each checked operation has a limits variant (binary bounds share
`bound_with_limits`). Iterative merge sorting and checked permutation building
charge actual traversal/comparison work; final swaps and discarded string
destructors are precharged before mutation. Allocation uses fallible reservation.
Caller predicate work is external to these counters. Standard `Vec`/`Clone`,
`PartialEq`, `Display`, `clear`, and direct field operations retain normal Rust
ownership/allocation behavior; these counters apply to the documented checked
scientific operations, not every standard-library operation.

## Shared array descriptions

The existing `DataArray<T>` now also carries `metadata: MetaInfo` and
`data_processing: Vec<Arc<DataProcessing>>`, matching the source
MetaInfoDescription payload without a duplicate mobility-only array model.
`DataArray::new` and `Default` initialize them empty. Source processing records
are shared handles; native clones preserve that sharing. Native equality compares
the complete represented values rather than emulating pointer addresses.
Existing public `DataArray { name, data }` literals must add
`..Default::default()` or use `DataArray::new(name, data)`.

Existing spectrum/chromatogram sort/select and theoretical generators preserve
the fields. HiRes picking deliberately preserves them on its retained mobility
array as a native retention correction: source `PeakPickerHiRes.cpp:134–136`
constructs a fresh array and copies only its name, losing its description.
Peak values and weighted mobility arithmetic are unchanged. Arrays reported as
omitted remain omitted. Spectrum and chromatogram validation meters
description traversal against a separate fixed 50M-work/256-MiB description
budget before inspecting nested records. Shared processing payload is counted
conservatively even when cloning only copies Arc handles. The existing processing
algorithms' own scientific work limits are otherwise unchanged; this does not
claim a new single budget across their historical whole-container clones.
The summary array-clear operation and annotation replacements charge description
destruction using their existing shared operation budgets. RNA append stages only
array values, so it neither clones nor destroys descriptions.

Name-and-values-only XML projections reject nonempty descriptions before output;
see [DATA_ARRAY_XML_SUPPORT.md](DATA_ARRAY_XML_SUPPORT.md). Native storage does
not silently discard metadata merely because an existing transport cannot yet
represent it. This additive change does not claim completion of source
DataArrays comparison operators or generic MetaInfoDescription ordering.

## Verification

Tests include the literal 21-peak source fixture, source nearest/highest/TIC
assertions, source tie permutations, full array metadata and shared ownership,
empty cases, subrange searches, negative tolerance branches, signed values,
f32 rounding, checked limits, unchanged state after errors and no-copy string
movement. [data_array_descriptions.rs](../tests/data_array_descriptions.rs) checks
existing containers, peptide/RNA append, HiRes mobility picking and clear budgets.
No C++ build or execution is used to obtain expected values.
