# Ranges: RangeManager, SpectrumRangeManager, ChromatogramRangeManager

This native group ports `OpenMS/KERNEL/RangeManager.h`,
`OpenMS/KERNEL/SpectrumRangeManager.h` and
`OpenMS/KERNEL/ChromatogramRangeManager.h` at core SDK commit
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The implementation is
[ranges.rs](../src/kernel/ranges.rs); the ported class tests and native
boundaries are in [ranges.rs](../tests/ranges.rs); the source hashes and
anchors are in [ranges_provenance.json](../tests/data/ranges_provenance.json).
It adds no dependency and does not build or execute C++.

## Design: pure value types, ranges computed on demand

The source mixes the range **algebra** (`RangeBase`, `RangeManager`) into every
peak container through `RangeManagerContainer`, which stores a cached range and
asks the caller to refresh it with `updateRanges()`. Every one of the 21 TOPP
`updateRanges()` call sites is an update-then-read; the mutable `getRange()` has
no mutating caller in core. This port therefore keeps the algebra as pure value
types and replaces the cache with on-demand accessors:

| Source (cached) | Native (on demand) |
| --- | --- |
| `MSSpectrum::updateRanges(); getRange()` / `getMinMZ()` … | [`MSSpectrum::range_manager`] |
| `MSChromatogram::updateRanges(); getRange()` | [`MSChromatogram::range_manager`] |
| `Mobilogram::updateRanges(); getRange()` | [`Mobilogram::range_manager`] |
| `MSExperiment::updateRanges(); spectrumRanges()` | [`MSExperiment::spectrum_range_manager`] (global + per MS level) |
| `MSExperiment::updateRanges(); chromatogramRanges()` | [`MSExperiment::chromatogram_range_manager`] |
| `MSExperiment::updateRanges(); combinedRanges()` / `getRange()` | [`MSExperiment::combined_range_manager`] |
| `clearRanges()` on a container | nothing: the next call recomputes |

The experiment accessors are named `*_range_manager` because
`MSExperiment::chromatogram_ranges()` and `combined_ranges()` already exist in
`experiment_summary.rs` with the three-dimensional `ExperimentRanges` result;
the new accessors are the four-dimensional, MS-level-aware counterparts and
agree with the legacy summaries where they overlap (tested).

The variadic template `RangeManager<RangeBases...>` becomes a run-time
dimension set: [`RangeManager`] holds an ordered set of [`MSDim`] values with
one [`RangeBase`] each. `RangeManager::new(&[MSDim])` builds any set; presets
name the source instantiations (`experiment()` = `RangeAllType` /
`MSExperiment::RangeManagerType`, `spectrum()`, `chromatogram()`,
`mobilogram()`, `chromatogram_manager()` = `ChromatogramRangeManager`,
`spectrum_manager()` = `SpectrumRangeManager::BaseType`). Cross-manager
operations act on the dimensions in common, exactly as the source's
`for_each_base_` folds with `is_base_of_v`; "no dimension in common" is the
source's `Exception::InvalidRange` → `Error::InvalidRange`. What the source
rejects at compile time — calling `getMinRT()` on a manager without RT — is
`Error::InvalidValue` here. Equality compares the ranges of each dimension and
ignores declaration order; order affects only `Display`.

## API mapping

### RangeManager.h

| Source member | Native counterpart |
| --- | --- |
| `enum MSDim { RT, MZ, INT, IM }` | `MSDim::{Rt, Mz, Intensity, Mobility}`, `MSDim::ALL`, `MSDim::label` |
| `RangeBase()` | `RangeBase::new` / `Default` |
| `RangeBase(const double single)` | `RangeBase::singular` (`Result`: rejects non-finite) |
| `RangeBase(const double min, const double max)` (`@throws InvalidRange`) | `RangeBase::from_min_max` → `Error::InvalidRange` when `min > max` |
| copy / move ctor, copy / move assignment, d'tor | `Clone`, `Copy`, ordinary assignment and drop |
| `operator RangeRT() / RangeMZ() / RangeIntensity() / RangeMobility()` | not ported: the typed subclasses do not exist; a `RangeBase` is placed into a dimension with `RangeManager::extend_range(dim, &range)` or `range_for_dim_mut` |
| `clear()` | `RangeBase::clear` |
| `isEmpty()` | `RangeBase::is_empty` |
| `contains(const double)` | `RangeBase::contains` |
| `contains(const RangeBase&)` | `RangeBase::contains_range` |
| `setMin` / `setMax` | `RangeBase::set_min` / `set_max` (`Result`: rejects non-finite) |
| `getMin` / `getMax` (`@throws InvalidRange` when empty) | `RangeBase::min` / `max` → `Error::InvalidRange` when empty |
| `extend(const RangeBase&)` | `RangeBase::extend` |
| `extend(const double)` | `RangeBase::extend_value` (`Result`: rejects non-finite) |
| `extendLeftRight` | `RangeBase::extend_left_right` (`Result`) |
| `minSpanIfSingular` | `RangeBase::min_span_if_singular` (`Result`) |
| `clampTo` (`@throw InvalidRange` if other empty) | `RangeBase::clamp_to` → `Error::InvalidRange` |
| `pushInto` (`@throw InvalidRange` if sandbox empty) | `RangeBase::push_into` → `Error::InvalidRange` |
| `scaleBy` | `RangeBase::scale_by` (`Result`) |
| `shift` | `RangeBase::shift` (`Result`) |
| `center()` (NaN when empty) | `RangeBase::center` → `Option<f64>`, `None` when empty |
| `getSpan()` (NaN when empty) | `RangeBase::span` → `Option<f64>`, `None` when empty |
| `operator==` | `PartialEq` on the stored endpoints |
| `getNonEmptyRange()` | `RangeBase::non_empty_range` → `(f64::MIN, f64::MAX)` when empty |
| protected `min_`, `max_` | private fields |
| `operator<<(ostream&, const RangeBase&)` | `Display`: `[min, max]` or `[, ]` |
| `RangeRT` / `RangeMZ` / `RangeIntensity` / `RangeMobility` structs, `DIM` | not ported as types: each is a `RangeBase` stored under `MSDim::Rt` / `Mz` / `Intensity` / `Mobility` in a `RangeManager` |
| `setMinRT` / `setMaxRT` / `getMinRT` / `getMaxRT` / `extendRT` / `containsRT(double)` / `containsRT(RangeBase)` | `RangeManager::set_min_rt` / `set_max_rt` / `min_rt` / `max_rt` / `extend_rt` / `contains_rt` / `contains_rt_range` |
| the same seven for MZ | `RangeManager::set_min_mz` / `set_max_mz` / `min_mz` / `max_mz` / `extend_mz` / `contains_mz` / `contains_mz_range` |
| the same seven for Intensity | `RangeManager::set_min_intensity` / `set_max_intensity` / `min_intensity` / `max_intensity` / `extend_intensity` / `contains_intensity` / `contains_intensity_range` |
| the same seven for Mobility | `RangeManager::set_min_mobility` / `set_max_mobility` / `min_mobility` / `max_mobility` / `extend_mobility` / `contains_mobility` / `contains_mobility_range` |
| `RangeRT::isEmpty()` and siblings (base-qualified calls) | `RangeManager::is_dim_empty(dim)` |
| `RangeRT::clear()` / `RangeRT::scaleBy()` and siblings (base-qualified calls) | `RangeManager::clear_dim(dim)`; `range_for_dim_mut(dim)?.scale_by(...)` |
| `operator<<` for each typed range | `Display` of `RangeManager` prints `label: [min, max]` per dimension |
| `enum HasRangeType { ALL, SOME, NONE }` | `HasRangeType::{All, Some, None}` |
| `RangeManager<RangeBases...>` | `RangeManager` with a run-time dimension set: `RangeManager::new(&[MSDim])`, presets `experiment`, `spectrum`, `chromatogram`, `mobilogram`, `chromatogram_manager`, `spectrum_manager`; `dims`, `has_dim` |
| `ThisRangeType` | not needed: one type |
| `operator==` / `operator!=` | `PartialEq` (per-dimension ranges, order-insensitive) |
| `assignUnsafe` (returns bool) | `RangeManager::assign_unsafe` → `bool` |
| `assign` (`@throw InvalidRange`) | `RangeManager::assign` → `Error::InvalidRange` |
| `extendUnsafe` / `extend` | `RangeManager::extend_unsafe` / `extend` |
| `scaleBy` | `RangeManager::scale_by` (`Result`, atomic) |
| `minSpanIfSingular` | `RangeManager::min_span_if_singular` (`Result`, atomic) |
| `pushIntoUnsafe` / `pushInto` | `RangeManager::push_into_unsafe` → `Result<bool>` / `push_into` |
| `clampToUnsafe` / `clampTo` | `RangeManager::clamp_to_unsafe` → `bool` / `clamp_to` |
| `getRangeForDim(MSDim)` const / mutable (asserts presence) | `RangeManager::range_for_dim` / `range_for_dim_mut` → `Error::InvalidValue` when absent |
| `hasRange()` | `RangeManager::has_range` |
| `containsAll` (`@throws InvalidRange` if no overlap) | `RangeManager::contains_all` → `Result<bool>` |
| `clearRanges()` | `RangeManager::clear_ranges` |
| `clear(DIM_UNIT)` | `RangeManager::clear_dim(MSDim)`; the three IM units collapse to `MSDim::Mobility`, as the source's switch does |
| `printRange(ostream&)`, `operator<<(RangeManager)` | `Display` / `to_string()` |
| protected `for_each_base_` (2), `static_for_each_base_` | private iteration over the dimension array |
| `RangeManagerContainer::updateRanges()` (pure virtual) | not needed: containers compute ranges on demand |
| `RangeManagerContainer::getRange()` const | `MSSpectrum::range_manager`, `MSChromatogram::range_manager`, `Mobilogram::range_manager`, `MSExperiment::combined_range_manager` |
| `RangeManagerContainer::getRange()` mutable ("avoid updateRanges() for minor changes") | not ported: there is no cache to edit; the returned `RangeManager` is an owned value |
| `RangeManagerContainer` virtual d'tor | not needed |
| `using RangeAllType` | `RangeManager::experiment()` |

### SpectrumRangeManager.h

| Source member | Native counterpart |
| --- | --- |
| `SpectrumRangeManager : RangeManager<RangeMZ, RangeIntensity, RangeMobility, RangeRT>` | `SpectrumRangeManager` holding a global `RangeManager::spectrum_manager()` and a `BTreeMap<u32, RangeManager>` |
| `BaseType` | `RangeManager` (`RangeManager::spectrum_manager()` shape) |
| default / copy / move ctor, copy / move assignment, d'tor | `new` / `Default`, `Clone`, assignment, drop |
| inherited base members (`getMinRT()`, `clear(DIM_UNIT)`, …) on the global ranges | `SpectrumRangeManager::global` / `global_mut` |
| `clearRanges()` | `SpectrumRangeManager::clear_ranges` |
| `extend(const BaseType& other, UInt ms_level = 0)` | `SpectrumRangeManager::extend(&RangeManager, ms_level)` |
| `byMSLevel(UInt ms_level = 0)` (`@throw InvalidValue` if absent) | `SpectrumRangeManager::by_ms_level` → `Option<&RangeManager>`; `None` for an unknown level **and for level 0**, which the source never stores in the map (it throws) — read `global()` instead |
| `getMSLevels()` | `SpectrumRangeManager::ms_levels` → `BTreeSet<u32>` |
| `extendRT(double rt, UInt ms_level = 0)` | `SpectrumRangeManager::extend_rt(rt, ms_level)` |
| `extendMZ(double mz, UInt ms_level = 0)` | `SpectrumRangeManager::extend_mz(mz, ms_level)` |
| `extendUnsafe(const MSSpectrum&, UInt ms_level = 0)` | `SpectrumRangeManager::extend_spectrum(&MSSpectrum, ms_level)` (computes `MSSpectrum::range_manager` instead of reading a cache; `Result`) |
| protected `ms_level_ranges_` | private map |

The pinned header declares no `clear(dim)`, `insert` or `getRange` of its own;
`clear(DIM_UNIT)` is the inherited base member and, as in the source, acts on
the global ranges only (`global_mut().clear_dim(dim)`).

### ChromatogramRangeManager.h

| Source member | Native counterpart |
| --- | --- |
| `ChromatogramRangeManager : RangeManager<RangeRT, RangeIntensity, RangeMZ>` | `RangeManager::chromatogram_manager()` |
| `BaseType` | `RangeManager` |

The class adds nothing beyond the alias; no separate Rust type is warranted.

## Preserved source conventions

- **Empty representation.** `min > max` marks an empty range; the default is
  `(f64::MAX, f64::MIN)`. `extend_left_right(-19)` on `[2, 8]` leaves the
  inverted pair `(23, -13)`, which is empty but not equal to a default-empty
  range — the class test's `TEST_EQUAL(empty1, empty2)` compares two defaults.
- **Setter symmetry.** `set_min(7)` on `[4, 6]` yields `[7, 7]`; `set_max(2)`
  yields `[2, 2]` (class test `setMin` / `setMax`).
- **Tie rule.** `extend` keeps the existing endpoint on equality, as
  `std::min`/`std::max` return their first argument; this preserves the sign
  of a zero endpoint, and the experiment merges spectra before chromatograms
  (`MSExperiment.cpp:718-719`), so a `0.0` spectrum endpoint beats a `-0.0`
  chromatogram endpoint. Tested.
- **`pushInto`.** Contained → unchanged; wider than the sandbox → cut to the
  sandbox span keeping the minimum; then shifted right or left. All eleven
  class-test assertions reproduce.
- **`clampTo`** may empty a dimension (`[1, 47110]` clamped to `[-10, -9]` is
  empty), and the manager-level variants skip dimensions that are empty in
  `rhs`, so `rm.clampTo(rmi)` leaves RT untouched when `rmi.RT` is empty.
- **`containsAll`.** An empty `rhs` dimension is contained; a non-empty `rhs`
  dimension is not contained in an empty `self` dimension; only overlapping
  dimensions count; no overlap is an error.
- **`hasRange`.** `NONE` / `SOME` / `ALL` over the dimensions the manager
  carries, so an experiment manager whose mobility is empty reports `Some`.
- **Scaling** a singular or empty dimension is a no-op
  (`scaleBy` class test, `rtmz == copy`).
- **Mobility source** (`MSSpectrum.cpp:586-604`): when the spectrum has an
  ion-mobility float data array (source `containsIMData()`), every value of the
  first such array extends mobility and the scalar drift time is ignored;
  otherwise the scalar `drift_time` extends mobility whenever it is not the
  sentinel `-1` (`IMTypes::DRIFTTIME_NOT_SET`) — the comparison is exact, so
  `-2` counts, as in the source. The array is recognised as the source
  `IMDataArrayUtils::getIMUnit` does: an exact PSI-MS name that is a child of
  `MS:1002893 ! ion mobility array` at the pinned CV (nine names), or a name
  starting with one of the `Constants::UserParam` fallbacks `"Ion Mobility"`,
  `"inverse reduced ion mobility"`, `"mean inverse reduced ion mobility array"`.
  The unit is not needed for a range and is not derived.
- **Experiment roles** (`MSExperiment.cpp:670-720`):
  - *spectrum-only*: every spectrum's m/z, intensity and mobility plus its RT
    extend the global ranges and its MS level's ranges
    (`KERNEL/MSExperiment.cpp:693-699`); a spectrum
    without peaks still contributes its RT (`:696`); a level-0 spectrum extends
    the global ranges twice and registers no level, because level 0 addresses
    the global ranges (`SpectrumRangeManager.h:84,126,137,148`).
  - *chromatogram-only*: every chromatogram's RT and intensity from its points,
    and its `getMZ()` — the **product** m/z (`MSChromatogram.cpp:81-84`) — even
    when it has no points (`KERNEL/MSExperiment.cpp:713`), so a chromatogram
    without a product
    contributes m/z `0`.
  - *combined*: the global spectrum ranges merged first, then the chromatogram
    ranges (`:718-719`), so combined RT and intensity include chromatogram
    points and combined m/z includes product m/z. The pre-existing
    `MSExperiment::ranges(ms_level)` iterates spectra only and is left as it is;
    the combined role lives in `combined_range_manager`.
  - The source's `updateRanges()` also refreshes each spectrum's and
    chromatogram's own cache as a side effect (`:691`, `:709`); on-demand
    computation covers this implicitly.
- **Display** prints one `label: [min, max]` line per dimension in declaration
  order, with `[, ]` for an empty dimension, reproducing the `printRange`
  class-test string exactly.

## Native differences

- **No cache.** Nothing on any container is mutated by a range query; every
  accessor returns an owned `RangeManager` computed from the current data. The
  mutable `getRange()` and the `updateRanges()` obligation have no counterpart.
- **Run-time dimension set.** Calling a typed accessor on a manager that lacks
  the dimension is `Error::InvalidValue` rather than a compile error;
  `RangeManager::new` rejects an empty or duplicated set the same way.
- **Finite domain.** Every value entering a `RangeBase` (constructors, setters,
  `extend_value`, arithmetic arguments) must be finite, and an arithmetic result
  that leaves the finite domain is `Error::InvalidRange` with the value
  unchanged. The source silently ignores a NaN in `extend` and lets NaN or
  infinities poison a range through `setMin`, `extendLeftRight`, `scaleBy` or
  `shift`. Manager-level `scale_by`, `min_span_if_singular` and
  `push_into_unsafe` are atomic: on error no dimension changes.
- **`Option` for absence.** `center()` and `span()` return `None` where the
  source returns NaN; `by_ms_level` returns `None` where the source throws
  `Exception::InvalidValue`; `min()`/`max()` keep the source's error (an empty
  range is a caller mistake, not an absence).
- **`SpectrumRangeManager::extend` registers a level only on success.** The
  source's `ms_level_ranges_[ms_level].extend(other)` default-constructs the
  entry before `extend` could throw; with the source's fixed `BaseType` the
  throw cannot occur, so this is unobservable there.
- **Validation.** `range_manager()` runs the container's `validate()` first
  (finite coordinates, positive MS level, aligned arrays) and additionally
  rejects a non-finite drift time or ion-mobility array value.
- **Formatting.** `Display` uses Rust's shortest round-trip number formatting;
  the source stream default prints six significant digits and may use
  exponent notation. The class-test integers format identically.
- **Bounds.** `RangeManager::MAX_ITEMS` (peaks, mobility values, spectra and
  chromatograms visited by one call) and `RangeManager::MAX_BYTES` (per-level
  managers an experiment query may allocate) are checked in a preflight before
  any allocation; exceeding either is `Error::InvalidValue`. The per-level
  byte preflight counts distinct levels incrementally and stops at the budget,
  so its own temporary set cannot exceed the budget it enforces. These bounds
  are far above any class-test input and are not exercised by a test.
- **Serial.** The source range code carries no `#pragma omp`; nothing is
  parallel here either.
- **Not added.** `SpectrumRanges` / `ExperimentRanges` in `src/kernel.rs`
  still carry no mobility dimension; this work package could not edit that
  file beyond registering the module. The four-dimensional managers returned
  here carry mobility, so the gap is in the legacy summary structs only.

## Checked boundaries and evidence

Evidence is **tier 3 (source review)**: every expected literal in
`tests/ranges.rs` is transcribed from `RangeManager_test.cpp` and
`SpectrumRangeManager_test.cpp` at `bc9cc12`, and the container behaviour is
read from `MSSpectrum.cpp`, `MSChromatogram.cpp`, `Mobilogram.cpp`,
`MSExperiment.cpp` and `IMDataArrayUtils.cpp` at the line anchors recorded in
the manifest. No C++ was compiled or executed; no retained C++ output exists
for these headers, and no oracle driver was written.

Section accounting (`RangeManager_test.cpp`, 38 sections → 38 Rust tests;
`SpectrumRangeManager_test.cpp`, 10 sections → 10 Rust tests). The three
`NOT_TESTABLE` sections (`isEmpty`, `getMin`, `getMax`) and the two
pointer-only sections (`RangeMType()`, `~RangeMType()`) are ported as small
tests of the same member so that no section is unaccounted for. Class-test
`RM::updateRanges()` / `updateRanges2()` are the test helpers
`rm_update_ranges` / `rm_update_ranges2`. `TEST_REAL_SIMILAR` values are all
exactly representable (`-47`, `1700`, `-23553.5`, …) and are asserted exactly.

Native tests additionally cover: every preset dimension set; every typed
accessor family; `extend_range`, `clear_dim`, `min_span_if_singular`; the
`*_unsafe` overlap reports; non-finite rejection with unchanged state for all
`RangeBase` mutators and for manager-level atomic operations; the signed-zero
tie rule; `Display` of empty and fractional ranges; spectrum ranges from peaks,
from the scalar drift time (including `-2`, NaN, and the `-1` sentinel), from
each recognised ion-mobility array name, from an empty placeholder array and a
NaN array value; chromatogram and mobilogram ranges; the three experiment roles
on one fixture with a peakless spectrum, a pointless chromatogram, three MS
levels, a scalar and an array mobility source; level-0 spectra; empty
experiments; a chromatogram without a product; agreement with the legacy
`combined_ranges()`; and failure propagation from a malformed spectrum or
chromatogram.

Self-audit: 48 ported, 0 mapped-with-evidence, 0 mapped-without-evidence,
0 unaccounted.
