# Mobility peaks and mobilograms

This native group is based on OpenMS Core SDK commit
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. It requires no new dependency.
The implementation is [mobilogram.rs](../src/kernel/mobilogram.rs); source
fixtures and independent boundaries are in [mobilogram.rs](../tests/mobilogram.rs).
The exact source hashes are recorded in
[mobilogram_provenance.json](../tests/data/mobilogram_provenance.json).
The inherited range surface lives in [ranges.rs](../src/kernel/ranges.rs) as
`Mobilogram::range_manager()`, and the class test's 48 sections are ported in
[chromatogram_merge.rs](../tests/chromatogram_merge.rs), audited at the end of
this document with
[chromatogram_merge_provenance.json](../tests/data/chromatogram_merge_provenance.json).

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

The table below now lists every public member of `Mobilogram.h`, so that a member
without a native counterpart cannot go unnoticed. Members marked "not ported"
state why.

| Source operation | Native operation |
| --- | --- |
| `struct RTLess` / `RTLess::operator()` | compare the public `rt` fields: `a.rt < b.rt` |
| `PeakType` (`MobilityPeak1D`) | `kernel::MobilityPeak1D` |
| `CoordinateType` | `f64` |
| `ContainerType` (`std::vector<MobilityPeak1D>`) | `Vec<MobilityPeak1D>`, the public `peaks` field |
| `RangeManagerContainerType`, `RangeManagerType` | `kernel::ranges::RangeManager::mobilogram()` |
| `FloatDataArray(s)`, `StringDataArray(s)`, `IntegerDataArray(s)` | `DataArray<f32>`, `DataArray<String>`, `DataArray<i32>` and their `Vec`s |
| `Iterator`/`iterator`, `ConstIterator`/`const_iterator`, `ReverseIterator`/`reverse_iterator`, `ConstReverseIterator`/`const_reverse_iterator` | slice iterators over `peaks`; not ported as type aliases |
| `Mobilogram()`, copy and move constructors, `~Mobilogram()` | `Mobilogram::new()` / `Default`, `Clone`, a Rust move, `Drop` |
| `operator=(const Mobilogram&)`, `operator=(Mobilogram&&)` | assignment of a `Clone` or of a move |
| `operator==` | `source_equal` ignores arrays, exactly as source |
| `operator!=` | `!source_equal(..)` / `!=` |
| `operator[]`, `front`, `back`, `begin`, `cbegin`, `end`, `cend`, `rbegin`, `crbegin`, `rend`, `crend`, `empty`, `erase`, `push_back`, `emplace_back`, `pop_back`, `insert`, `resize`, `reserve`, `size` | the standard `Vec`/slice API on the public `peaks` field |
| `swap` | `swap_peak_data` (see below) |
| `updateRanges`, mobility/intensity extrema | `ranges()` returns current optional `NumericRange` values; `range_manager()` returns the equivalent `RangeManager` value |
| `getRT`, `setRT` | public `rt` field |
| `getDriftTimeUnit`, `setDriftTimeUnit` | public `drift_time_unit` field |
| `getDriftTimeUnitAsString` | `drift_time_unit_as_str()` |
| `getFloatDataArrays` const/mutable, `setFloatDataArrays` | public `float_data_arrays` field |
| `getStringDataArrays` const/mutable, `setStringDataArrays` | public `string_data_arrays` field |
| `getIntegerDataArrays` const/mutable, `setIntegerDataArrays` | public `integer_data_arrays` field |
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
| protected `checkDataArraySizes_` | `validate()`, and the array check inside `sort_by` / `select` |
| protected `data_`, `retention_time_`, `drift_time_unit_`, `float_data_arrays_`, `string_data_arrays_`, `integer_data_arrays_` | the public fields above |
| `operator<<(std::ostream&, const Mobilogram&)` | `Display for Mobilogram`, with an optional formatter precision |
| source partial `swap` | `swap_peak_data` exchanges peaks/RT/unit only |
| full native value equality / swap | `PartialEq` includes arrays; `std::mem::swap` exchanges everything |

### Inherited `RangeManagerContainer<RangeMobility, RangeIntensity>`

`Mobilogram` publicly inherits the range container, so every member of
`RangeManager.h` for the mobility and intensity dimensions is part of its public
surface. All of them map onto one accessor.

| Inherited member | Native counterpart |
| --- | --- |
| `updateRanges()` (pure virtual, overridden by `Mobilogram`) | `Mobilogram::range_manager()` in [ranges.rs](../src/kernel/ranges.rs) |
| `getRange()` const and mutable | `Mobilogram::range_manager()`, which returns an owned value rather than a reference into the mobilogram |
| `getMinMobility`, `getMaxMobility`, `getMinIntensity`, `getMaxIntensity` | `min_mobility()`, `max_mobility()`, `min_intensity()`, `max_intensity()` on the returned `RangeManager` |
| `setMinMobility`, `setMaxMobility`, `setMinIntensity`, `setMaxIntensity`, `extendMobility`, `extendIntensity`, `containsMobility`, `containsIntensity`, `clearRanges`, `clearRange`, and the whole `RangeManager` algebra (`assign`, `extend`, `scaleBy`, `minSpanIfSingular`, `pushInto`, `clampTo`, `hasRange`, `containsAll`) | the same operations on the `RangeManager` value; see [RANGES_SUPPORT.md](RANGES_SUPPORT.md) |

The source caches these ranges inside the mobilogram and requires an explicit
`updateRanges()` call to refresh them, which is why `select()` and `sort()` carry
`@note`s about what the cache does and does not survive. This port computes the
ranges on demand: `range_manager()` walks the current peaks each time, so a stale
or too-wide cache cannot exist, no mutation can invalidate anything, and the
cache-mutating members have no counterpart because there is no cache to mutate.
That divergence is deliberate and was ratified when `ranges.rs` was written; it
is what closes this header's recorded residual.

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

## Section audit: `Mobilogram_test.cpp`

All 48 `START_SECTION`s of `Mobilogram_test.cpp`, in file order, are ported into
[chromatogram_merge.rs](../tests/chromatogram_merge.rs); none is merely mapped
onto an existing test. The class-test literals are unchanged between the pinned
revision `54a232fe2cae9c590d5c997fa49d20e7769860fb` that this group was written
against and the current target `bc9cc12514c768385ce121d6ca4bb710fe1983c4`, which
[mobilogram_provenance.json](../tests/data/mobilogram_provenance.json) records as
an empty `source_changes` list. [mobilogram.rs](../tests/mobilogram.rs) remains
the home of the native boundary coverage — ceilings, atomicity, string-buffer
movement and shared array descriptions — that has no class-test section.

| # | Line | Section | Rust test |
| --- | --- | --- | --- |
| 1 | 48 | `Mobilogram()` | `mobilogram_source_construction_and_scalar_accessors` |
| 2 | 55 | `~Mobilogram()` | `mobilogram_source_construction_and_scalar_accessors` (explicit `drop`; the C++ section has no assertion) |
| 3 | 61 | `[EXTRA] Mobilogram()` | `mobilogram_source_construction_and_scalar_accessors` |
| 4 | 75 | `double getRT() const` | `mobilogram_source_construction_and_scalar_accessors` |
| 5 | 82 | `void setRT(double)` | `mobilogram_source_construction_and_scalar_accessors` |
| 6 | 91 | `double getDriftTimeUnit() const` | `mobilogram_source_construction_and_scalar_accessors` |
| 7 | 98 | `double getDriftTimeUnitAsString() const` | `mobilogram_source_construction_and_scalar_accessors` |
| 8 | 105 | `void setDriftTimeUnit(double)` | `mobilogram_source_construction_and_scalar_accessors` |
| 9 | 118 | `virtual void updateRanges()` | `mobilogram_source_update_ranges` |
| 10 | 148 | `Mobilogram(const Mobilogram&)` | `mobilogram_source_copy_move_and_assignment` |
| 11 | 167 | `Mobilogram(const Mobilogram&&)` | `mobilogram_source_copy_move_and_assignment` |
| 12 | 200 | `Mobilogram& operator=(const Mobilogram&)` | `mobilogram_source_copy_move_and_assignment` |
| 13 | 227 | `Mobilogram& operator=(const Mobilogram&&)` | `mobilogram_source_copy_move_and_assignment` |
| 14 | 274 | `bool operator==(const Mobilogram&) const` | `mobilogram_source_equality` |
| 15 | 301 | `bool operator!=(const Mobilogram&) const` | `mobilogram_source_equality` |
| 16 | 332 | `void sortByIntensity(bool reverse = false)` | `mobilogram_source_sorts` |
| 17 | 402 | `void sortByPosition()` | `mobilogram_source_sorts` |
| 18 | 442 | `[EXTRA] sorting reorders all data arrays alongside the peaks` | `mobilogram_source_sorting_regression` |
| 19 | 567 | `bool isSorted() const` | `mobilogram_source_is_sorted_and_predicate_is_sorted` |
| 20 | 591 | `template<class Predicate> bool isSorted(const Predicate&) const` | `mobilogram_source_is_sorted_and_predicate_is_sorted` |
| 21 | 614 | `template<class Predicate> void sort(const Predicate&)` | `mobilogram_source_predicate_sort` |
| 22 | 648 | `Mobilogram& select(const std::vector<Size>&)` | `mobilogram_source_select` |
| 23 | 731 | `Iterator MBEnd(CoordinateType)` | `mobilogram_source_mb_bounds` |
| 24 | 744 | `Iterator MBBegin(CoordinateType)` | `mobilogram_source_mb_bounds` |
| 25 | 758 | `Iterator MBBegin(Iterator, CoordinateType, Iterator)` | `mobilogram_source_mb_bounds` |
| 26 | 772 | `ConstIterator MBBegin(ConstIterator, CoordinateType, ConstIterator) const` | `mobilogram_source_mb_bounds` |
| 27 | 786 | `Iterator MBEnd(Iterator, CoordinateType, Iterator)` | `mobilogram_source_mb_bounds` |
| 28 | 800 | `ConstIterator MBEnd(ConstIterator, CoordinateType, ConstIterator) const` | `mobilogram_source_mb_bounds` |
| 29 | 813 | `ConstIterator MBEnd(CoordinateType) const` | `mobilogram_source_mb_bounds` |
| 30 | 826 | `ConstIterator MBBegin(CoordinateType) const` | `mobilogram_source_mb_bounds` |
| 31 | 841 | `Iterator PosBegin(CoordinateType)` | `mobilogram_source_pos_bounds_are_mb_bounds` |
| 32 | 853 | `Iterator PosBegin(Iterator, CoordinateType, Iterator)` | `mobilogram_source_pos_bounds_are_mb_bounds` |
| 33 | 867 | `ConstIterator PosBegin(CoordinateType) const` | `mobilogram_source_pos_bounds_are_mb_bounds` |
| 34 | 879 | `ConstIterator PosBegin(ConstIterator, CoordinateType, ConstIterator) const` | `mobilogram_source_pos_bounds_are_mb_bounds` |
| 35 | 893 | `Iterator PosEnd(CoordinateType)` | `mobilogram_source_pos_bounds_are_mb_bounds` |
| 36 | 905 | `Iterator PosEnd(Iterator, CoordinateType, Iterator)` | `mobilogram_source_pos_bounds_are_mb_bounds` |
| 37 | 919 | `ConstIterator PosEnd(CoordinateType) const` | `mobilogram_source_pos_bounds_are_mb_bounds` |
| 38 | 931 | `ConstIterator PosEnd(ConstIterator, CoordinateType, ConstIterator) const` | `mobilogram_source_pos_bounds_are_mb_bounds` |
| 39 | 971 | `Size findNearest(CoordinateType) const` | `mobilogram_source_searches` |
| 40 | 993 | `Size findNearest(CoordinateType, CoordinateType) const` | `mobilogram_source_searches` |
| 41 | 1017 | `Size findNearest(CoordinateType, CoordinateType, CoordinateType) const` | `mobilogram_source_searches` |
| 42 | 1045 | `Size findHighestInWindow(CoordinateType, CoordinateType, CoordinateType) const` | `mobilogram_source_searches` |
| 43 | 1077 | `ConstIterator getBasePeak() const` | `mobilogram_source_base_peak_and_tic` |
| 44 | 1088 | `Iterator getBasePeak()` | `mobilogram_source_base_peak_and_tic` |
| 45 | 1099 | `PeakType::IntensityType calculateTIC() const` | `mobilogram_source_base_peak_and_tic` |
| 46 | 1108 | `void clear()` | `mobilogram_source_clear` |
| 47 | 1130 | `[Mobilogram::RTLess] bool operator()(const Mobilogram&, const Mobilogram&) const` | `mobilogram_source_rt_less` |
| 48 | 1165 | `[EXTRA] std::ostream& operator<<(std::ostream&, const Mobilogram&)` | `mobilogram_source_stream_layout` |

Three assertions have no runtime counterpart and are recorded rather than
asserted: the `noexcept` move-constructor check and the two moved-from emptiness
checks in sections 11 and 13. Rust moves never panic and the compiler rejects any
use of a moved-from value, so both properties hold statically. The file's four
`static_assert`s on rule-of-5/6 and fast-vector semantics are likewise compile-time
C++ properties with no Rust counterpart: `Clone`, `Default`, `Drop` and moves are
derived or built in.
