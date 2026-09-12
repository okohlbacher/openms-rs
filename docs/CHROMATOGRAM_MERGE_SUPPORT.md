# MSChromatogram: merging, subrange search and the remaining members

This group closes `KERNEL/MSChromatogram.h` at OpenMS Core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. It requires no new dependency.

The new implementation is
[chromatogram_merge.rs](../src/kernel/chromatogram_merge.rs); the rest of the
header was already covered by [kernel.rs](../src/kernel.rs),
[acquisition_fields.rs](../src/kernel/acquisition_fields.rs) and
[ranges.rs](../src/kernel/ranges.rs). The tests are
[chromatogram_merge.rs](../tests/chromatogram_merge.rs) and the exact source
hashes are in
[chromatogram_merge_provenance.json](../tests/data/chromatogram_merge_provenance.json).

Four members of the header had no native counterpart before this group:
`mergePeaks`, the eight subrange `RTBegin`/`RTEnd`/`PosBegin`/`PosEnd`
overloads, `template<class Predicate> void sort`, and the value predicates
`MZLess` / `operator==` / `operator<<`. Everything else is a re-statement of an
existing mapping, recorded here so the header's member table is complete in one
place.

## API mapping

Every public member of `MSChromatogram.h`. "Public field" means the native type
exposes the value directly instead of a getter/setter pair, which is the crate's
established choice for the kernel containers.

### Types and base-type definitions

| Source member | Native counterpart |
| --- | --- |
| `struct MZLess` / `MZLess::operator()` | `chromatogram_merge::chromatogram_mz_less(a, b)` |
| `PeakType` (`ChromatogramPeak`) | `kernel::ChromatogramPeak` |
| `CoordinateType` | `f64` |
| `ContainerType` (`std::vector<ChromatogramPeak>`) | `Vec<ChromatogramPeak>`, the public `peaks` field |
| `RangeManagerType` (`RangeManager<RangeRT, RangeIntensity>`) | `kernel::ranges::RangeManager::chromatogram()` |
| `FloatDataArray` / `FloatDataArrays` | `DataArray<f32>` / `Vec<DataArray<f32>>` |
| `StringDataArray` / `StringDataArrays` | `DataArray<String>` / `Vec<DataArray<String>>` |
| `IntegerDataArray` / `IntegerDataArrays` | `DataArray<i32>` / `Vec<DataArray<i32>>` |
| `Iterator`, `ConstIterator`, `ReverseIterator`, `ConstReverseIterator` | `std::slice::IterMut`, `Iter`, `Rev<IterMut>`, `Rev<Iter>` via `peaks` |
| `iterator`, `const_iterator`, `size_type`, `value_type`, `reference`, `const_reference`, `pointer`, `difference_type` | not ported as names: the `Vec`/slice API supplies each one, and a Rust alias for `usize` or `&T` carries no information |

### Container surface exported from `std::vector`

All of these are the standard `Vec`/slice API on the public `peaks` field. The
source re-exports them with `using` so that `MSChromatogram` looks like a vector
while privately inheriting one; the native type simply owns a `Vec`.

| Source member | Native counterpart |
| --- | --- |
| `operator[]` | `peaks[i]` |
| `begin`, `cbegin`, `end`, `cend` | `peaks.iter()`, `peaks.iter_mut()` |
| `rbegin`, `rend` | `peaks.iter().rev()` |
| `size` | `len()` (also `peaks.len()`) |
| `empty` | `is_empty()` |
| `front`, `back` | `peaks.first()`, `peaks.last()` |
| `resize`, `reserve` | `peaks.resize(..)`, `peaks.reserve(..)` |
| `push_back`, `emplace_back`, `pop_back` | `peaks.push(..)`, `peaks.pop()` |
| `insert`, `erase` | `peaks.insert(..)`, `peaks.remove(..)` / `peaks.drain(..)` |
| `swap` | `std::mem::swap(&mut a.peaks, &mut b.peaks)` |

### Construction, assignment and equality

| Source member | Native counterpart |
| --- | --- |
| `MSChromatogram()` | `MSChromatogram::new()` / `Default` |
| copy constructor | `Clone` |
| move constructor | a Rust move; no separate operation |
| `~MSChromatogram()` | `Drop` (automatic) |
| `operator=(const MSChromatogram&)` | assignment of a `Clone` |
| `operator=(MSChromatogram&&)` | assignment of a move |
| `operator==` | `MSChromatogram::source_equal` (new); the derived `PartialEq` differs only in also comparing the name |
| `operator!=` | `!source_equal(..)` / `!=` |
| `updateRanges()` | `MSChromatogram::range_manager()` and `MSChromatogram::ranges()`, computed on demand |
| `operator<<(std::ostream&, const MSChromatogram&)` | `Display for MSChromatogram` (new) |

### Meta information

| Source member | Native counterpart |
| --- | --- |
| `getName()` / `setName(name)` | public `name` field |
| `getMZ()` | `product.mz` (source `getProduct().getMZ()`) |
| `getFloatDataArrays()` const / mutable, `setFloatDataArrays(fda)` | public `float_data_arrays` field |
| `getStringDataArrays()` const / mutable, `setStringDataArrays(sda)` | public `string_data_arrays` field |
| `getIntegerDataArrays()` const / mutable, `setIntegerDataArrays(ida)` | public `integer_data_arrays` field |

### Sorting

| Source member | Native counterpart |
| --- | --- |
| `sortByIntensity(bool reverse)` | `MSChromatogram::sort_by_intensity(reverse)` |
| `sortByPosition()` | `MSChromatogram::sort_by_position()` |
| `isSorted()` | `MSChromatogram::is_sorted()` |
| `template<class Predicate> void sort(lambda)` | `MSChromatogram::sort_by(less)` and `sort_by_with_limits(less, limits)` (new) |

### Searching a peak or peak range

| Source member | Native counterpart |
| --- | --- |
| `findNearest(rt)` | `MSChromatogram::find_nearest(rt)` |
| `RTBegin(rt)` (mutable and const) | `MSChromatogram::rt_begin(rt)` |
| `RTEnd(rt)` (mutable and const) | `MSChromatogram::rt_end(rt)` |
| `RTBegin(begin, rt, end)` (mutable and const) | `MSChromatogram::rt_begin_in(rt, range)` (new) |
| `RTEnd(begin, rt, end)` (mutable and const) | `MSChromatogram::rt_end_in(rt, range)` (new) |
| `PosBegin(rt)` / `PosBegin(begin, rt, end)`, both const-qualifications | the same `rt_begin` / `rt_begin_in`; the source declares `PosBegin` an alias of `RTBegin` and its body is a forwarding call |
| `PosEnd(rt)` / `PosEnd(begin, rt, end)`, both const-qualifications | the same `rt_end` / `rt_end_in`, for the same reason |

The source needs four overloads per bound because it returns iterators, which
must come in mutable and const flavours, and because a subrange is expressed as
an iterator pair. The native functions return an index, which is neither
mutable nor const and which addresses the whole chromatogram, so one function
covers the four. The index a subrange search returns is a chromatogram index,
not an offset into the subrange.

### Subsetting and merging

| Source member | Native counterpart |
| --- | --- |
| `clear(bool clear_meta_data)` | `MSChromatogram::clear(clear_metadata)` |
| `select(indices)` | `MSChromatogram::select(&indices)` |
| `selectUnchecked(indices)` | not ported: the checked `select` is the only entry point. The source's unchecked variant exists to skip a bounds re-scan on a permutation it built itself and documents a repeated index as undefined behaviour; the native `select` always checks range and uniqueness, so no API can reach that state |
| `mergePeaks(other, add_meta)` | `MSChromatogram::merge_peaks(&other, add_meta)` and `merge_peaks_with_options(&other, options)` (new) |

### Protected members

| Source member | Native counterpart |
| --- | --- |
| `checkDataArraySizes_()` | `MSChromatogram::validate_data_arrays()`, which is public because a caller who edits the public array fields needs it |
| `name_`, `float_data_arrays_`, `string_data_arrays_`, `integer_data_arrays_` | the public fields above |

### Inherited surfaces

`MSChromatogram` privately inherits `std::vector<ChromatogramPeak>` and publicly
inherits `RangeManagerContainer<RangeRT, RangeIntensity>` and
`ChromatogramSettings`.

| Inherited surface | Native counterpart |
| --- | --- |
| `RangeManagerContainer::getRange()`, `updateRanges()`, and the `RangeRT`/`RangeIntensity` minimum, maximum, `extend*`, `contains*` and `clearRanges` algebra | `MSChromatogram::range_manager()` returns a fresh `RangeManager` value carrying the same algebra; see [RANGES_SUPPORT.md](RANGES_SUPPORT.md). The source's cached ranges are replaced by on-demand computation |
| `ChromatogramSettings` (native id, instrument settings, acquisition info, source file, precursor, product, data processing, chromatogram type, `MetaInfoInterface`) | the public `native_id`, `instrument_settings`, `acquisition_info`, `source_file`, `precursor`, `product`, `data_processing`, `chromatogram_type` and `metadata` fields; see [CHROMATOGRAM_TOOLS_SUPPORT.md](CHROMATOGRAM_TOOLS_SUPPORT.md) and [RECORD_METADATA_MIGRATION.md](RECORD_METADATA_MIGRATION.md) |
| `ChromatogramSettings::getComment()` / `setComment(comment)` | not ported: no native field carries it. It participates in the source `operator==`, so `source_equal` cannot observe a difference that does not exist here. `ChromatogramSettings.h` is a separate header and its own ledger entry |

## Preserved source conventions

**The merge key is a millisecond bucket.** `setSumSimilarUnion` compares
`round(rt * 1000.0)`; the comment above it calls this "within 1/1000 seconds",
which is looser than the code. Two points 0.00002 s apart (`0.00149` and
`0.00151`) fall in buckets 1 and 2 and stay separate, while two points 0.00099 s
apart (`0.00050` and `0.00149`) share bucket 1 and are summed. The bucketing is
reproduced rather than replaced by a distance test, because it decides which
points are summed. `f64::round` and C's `round` both break halves away from
zero, so the bucket boundaries agree exactly.

**Merge order and tie handling.** The union drains the destination first when
the other side is exhausted and vice versa, in the source's order of checks. A
tie keeps the destination's retention time, not the other's, and sets the
intensity to the `f32` sum of the two.

**The destination's product m/z is not touched**, as the header states, so a
merged chromatogram keeps advertising its own transition.

**`add_meta` accumulates.** The source reads any existing
`merged_chromatogram_mzs` value, appends `other.getMZ()` and writes the list
back, so repeated merges build a list. `merged_chromatogram_mzs()` reads it.

**Search semantics.** `RTBegin`/`PosBegin` are `std::lower_bound`, so an exact
query selects the first duplicate; `RTEnd`/`PosEnd` are `std::upper_bound`. A
subrange whose begin equals its end returns that position, matching
`lower_bound(begin, .., begin)`. `findNearest` resolves a midpoint tie towards
the lower retention time, because the source compares with `<` and falls through
to the predecessor.

**`sort(Predicate)` validates before it compares.** The source `@exception`
requires that a mis-sized data array is rejected before the lambda is invoked,
so the chromatogram is left unchanged and the predicate never indexes out of
bounds. `sort_by` checks the arrays first for the same reason.

**Sorting is stable** and moves every non-empty parallel array in the same
permutation. An array with no entries is a declared placeholder and is left
alone, here as everywhere else in the kernel.

**`clear(false)` still drops the data arrays**, because they are parallel to the
points; only the descriptive metadata is governed by the flag.

**The stream layout, including its empty settings block.**
`operator<<(std::ostream&, const ChromatogramSettings&)` takes its argument
unnamed and writes only `-- CHROMATOGRAMSETTINGS BEGIN --` and
`-- CHROMATOGRAMSETTINGS END --`, so no setting has ever appeared in a
chromatogram dump. `Display for MSChromatogram` reproduces that, delimiters
included, because callers and the class test match on the surrounding text.

## Native differences

**Sorted input is checked, not assumed.** The source `@note` documents unsorted
input to `mergePeaks` and to every range search as undefined and checks nothing.
Both check and return `Error::UnsortedData`: the cost is one pass over data the
operation reads anyway, and an unnoticed unsorted input silently produces a
wrongly summed, wrongly ordered result.

**Annotation arrays are a decision, not an accident.** The source updates the
points and leaves the float, string and integer arrays at their pre-merge
length. The header's `@note` says they are "not guaranteed to be correct"; the
consequence is stronger, because the chromatogram then fails its own
`checkDataArraySizes_` and a later `sortByPosition()` or `select()` throws.
`MergedDataArrays` makes the three possibilities explicit: `Reject` (the
default) refuses a merge when either side carries a non-empty array, `Drop`
removes them so the result is consistent, and `Source` reproduces the source's
untouched arrays. `ChromatogramMergeOptions::source()` selects the last.

**A non-finite summed intensity is an error.** Two large `f32` intensities can
sum to infinity; the source stores it. The port reports
`Error::InvalidValue` and leaves both chromatograms unchanged.

**No range cache to go stale.** The source's `mergePeaks` changes both the
retention-time span and the maximum intensity without calling `updateRanges()`
and without documenting that it does not, so a destination whose ranges were
current becomes one whose cached range is too *narrow* — the opposite of the
documented `select()` case, and a range that no longer contains its own points.
This port has no cache: `range_manager()` recomputes, so the state cannot exist.
Recorded as a C++ finding.

**`other` is borrowed, not mutably referenced.** The source signature takes
`MSChromatogram&` although the body only reads it, which forces a caller holding
a `const` chromatogram to copy. `merge_peaks` takes `&MSChromatogram`. A side
effect is that `a.merge_peaks(&a, ..)` is a compile error, where the source
permits the aliased call.

**Equality is spelled out twice.** `source_equal` is the source predicate, which
ignores `name_`. The derived `PartialEq` also compares the name, because a value
type that silently ignores one of its own fields is a trap for a Rust reader.
The source also excludes the cached ranges from equality; there is no cache
here, so nothing is excluded for that reason.

**No unchecked selection.** `selectUnchecked` is not ported; see the table.

**Predicate sorting can invoke the predicate twice per comparison**, because the
standard-library sort is driven by a three-way ordering while the source
predicate answers only "less". A predicate that is not a strict weak ordering
panics here instead of being undefined behaviour.

**Serial only.** Neither the source nor this module parallelises, so there is no
performance gap to record for the merge itself. The crate is serial by design;
no `#pragma omp` appears in `MSChromatogram.cpp` or `Mobilogram.cpp`.

## Checked boundaries and evidence

`ChromatogramMergeLimits` defaults to 10,000,000 combined points, 1,000,000
entries in the `merged_chromatogram_mzs` list and 256 MiB of temporary storage.
`ChromatogramSortLimits` defaults to 10,000,000 points and 256 MiB. Both are
checked before anything is allocated or mutated. The merge builds its result in
a temporary with fallible reservation and prepares the metadata value before
either is committed, so a rejection — an exceeded ceiling, an unsorted input, an
invalid record, a non-finite sum, a non-finite product m/z, a
`merged_chromatogram_mzs` value that is not a float list, or a refused
annotation array — leaves this chromatogram byte-for-byte as it was and never
writes to `other` at all. `sort_by` checks its ceilings and the array sizes
before the first predicate call. Evaluating a caller predicate is the caller's
cost and is not metered, as in [MOBILOGRAM_SUPPORT.md](MOBILOGRAM_SUPPORT.md).

Evidence is **tier 3 (source review)**: every expected value in
[chromatogram_merge.rs](../tests/chromatogram_merge.rs) is a literal transcribed
from `MSChromatogram_test.cpp` or `Mobilogram_test.cpp` at the pinned revision,
or is read off the implementation lines recorded as `source_anchors` in the
manifest. No C++ was built, linked or executed and no retained C++ output is
claimed. The native-only expectations — the millisecond-bucket demonstration,
the three annotation-array policies, the ceilings and the transactional
rejections — are tier 4.

## Section audit: `MSChromatogram_test.cpp`

All 43 `START_SECTION`s, in file order. Every one is ported into
[chromatogram_merge.rs](../tests/chromatogram_merge.rs); none is mapped onto a
pre-existing test.

| # | Line | Section | Rust test |
| --- | --- | --- | --- |
| 1 | 34 | `MSChromatogram()` | `source_construction_destruction_and_name` |
| 2 | 41 | `virtual ~MSChromatogram()` | `source_construction_destruction_and_name` (explicit `drop`; the C++ section has no assertion) |
| 3 | 69 | `const std::string& getName() const` | `source_construction_destruction_and_name` |
| 4 | 78 | `void setName(const std::string&)` | `source_construction_destruction_and_name` (the C++ section is `NOT_TESTABLE`) |
| 5 | 84 | `const FloatDataArrays& getFloatDataArrays() const` | `source_data_array_accessors` |
| 6 | 91 | `FloatDataArrays& getFloatDataArrays()` | `source_data_array_accessors` |
| 7 | 99 | `const StringDataArrays& getStringDataArrays() const` | `source_data_array_accessors` |
| 8 | 106 | `StringDataArrays& getStringDataArrays()` | `source_data_array_accessors` |
| 9 | 114 | `const IntegerDataArrays& getIntegerDataArrays() const` | `source_data_array_accessors` |
| 10 | 121 | `IntegerDataArrays& getIntegerDataArrays()` | `source_data_array_accessors` |
| 11 | 129 | `void sortByIntensity(bool reverse=false)` | `source_sort_by_intensity_permutes_every_array` |
| 12 | 229 | `void sortByPosition()` | `source_sort_by_position_permutes_every_array` |
| 13 | 320 | `bool isSorted() const` | `source_is_sorted` |
| 14 | 345 | `[EXTRA] sorting reorders string and integer data arrays even when no float data array is present` | `source_sorting_and_selection_regression` |
| 15 | 504 | `Size findNearest(CoordinateType rt) const` | `source_find_nearest` |
| 16 | 549 | `Iterator RTBegin(CoordinateType)` | `source_rt_bounds_and_subranges` |
| 17 | 579 | `Iterator RTBegin(Iterator, CoordinateType, Iterator)` | `source_rt_bounds_and_subranges` |
| 18 | 610 | `Iterator RTEnd(CoordinateType)` | `source_rt_bounds_and_subranges` |
| 19 | 640 | `Iterator RTEnd(Iterator, CoordinateType, Iterator)` | `source_rt_bounds_and_subranges` |
| 20 | 670 | `ConstIterator RTBegin(CoordinateType) const` | `source_rt_bounds_and_subranges` |
| 21 | 701 | `ConstIterator RTBegin(ConstIterator, CoordinateType, ConstIterator) const` | `source_rt_bounds_and_subranges` |
| 22 | 731 | `ConstIterator RTEnd(CoordinateType) const` | `source_rt_bounds_and_subranges` |
| 23 | 762 | `ConstIterator RTEnd(ConstIterator, CoordinateType, ConstIterator) const` | `source_rt_bounds_and_subranges` |
| 24 | 801 | `Iterator PosBegin(CoordinateType)` | `source_pos_bounds_are_rt_bounds` |
| 25 | 813 | `Iterator PosBegin(Iterator, CoordinateType, Iterator)` | `source_pos_bounds_are_rt_bounds` |
| 26 | 827 | `ConstIterator PosBegin(CoordinateType) const` | `source_pos_bounds_are_rt_bounds` |
| 27 | 839 | `ConstIterator PosBegin(ConstIterator, CoordinateType, ConstIterator) const` | `source_pos_bounds_are_rt_bounds` |
| 28 | 853 | `Iterator PosEnd(CoordinateType)` | `source_pos_bounds_are_rt_bounds` |
| 29 | 865 | `Iterator PosEnd(Iterator, CoordinateType, Iterator)` | `source_pos_bounds_are_rt_bounds` |
| 30 | 879 | `ConstIterator PosEnd(CoordinateType) const` | `source_pos_bounds_are_rt_bounds` |
| 31 | 891 | `ConstIterator PosEnd(ConstIterator, CoordinateType, ConstIterator) const` | `source_pos_bounds_are_rt_bounds` |
| 32 | 908 | `MSChromatogram(const MSChromatogram&)` | `source_copy_move_and_assignment` |
| 33 | 934 | `MSChromatogram(const MSChromatogram&&)` | `source_copy_move_and_assignment` |
| 34 | 972 | `MSChromatogram& operator=(const MSChromatogram&)` | `source_copy_move_and_assignment` |
| 35 | 1005 | `MSChromatogram& operator=(const MSChromatogram&&)` | `source_copy_move_and_assignment` |
| 36 | 1057 | `bool operator==(const MSChromatogram&) const` | `source_equality_ignores_only_the_name` |
| 37 | 1103 | `bool operator!=(const MSChromatogram&) const` | `source_equality_ignores_only_the_name` |
| 38 | 1152 | `virtual void updateRanges()` | `source_update_ranges` |
| 39 | 1179 | `void clear(bool clear_meta_data)` | `source_clear` |
| 40 | 1207 | `double getMZ() const` | `source_product_mz_and_mz_less` |
| 41 | 1218 | `[MSChromatogram::MZLess] bool operator()(...)` | `source_product_mz_and_mz_less` |
| 42 | 1237 | `void mergePeaks(MSChromatogram& other)` | `source_merge_peaks` |
| 43 | 1257 | `std::ostream& operator<<(std::ostream&, const MSChromatogram&)` | `source_stream_layout` |

Three assertions have no runtime counterpart and are recorded rather than
asserted: the `noexcept` move-constructor check (section 33) and the two
moved-from emptiness checks (sections 33 and 35). Rust moves never panic, and
the compiler rejects any use of a moved-from value, so both properties hold
statically.

Beyond the sections, `predicate_sort_follows_a_data_array_and_is_bounded`
exercises the ordering half of `sort(Predicate)`, which the class test covers
only through its rejection case, using the equivalent fixture from
`Mobilogram_test.cpp:614-646`. Four further tests —
`merge_key_is_a_millisecond_bucket_not_a_distance`,
`merge_drains_both_sides_and_keeps_the_destination_product`,
`merge_metadata_list_accumulates_and_rejects_a_wrong_type`,
`merge_annotation_array_policies` and `merge_rejections_are_transactional` —
cover merge behaviour that the single-assertion class-test section does not
reach.

The 48 sections of `Mobilogram_test.cpp` are audited in
[MOBILOGRAM_SUPPORT.md](MOBILOGRAM_SUPPORT.md) and are ported in the same test
binary.
