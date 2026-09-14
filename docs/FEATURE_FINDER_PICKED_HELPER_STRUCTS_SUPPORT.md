# FeatureFinderAlgorithmPicked helper structures

[`src/analysis/feature_finder_picked/helper_structs.rs`](../src/analysis/feature_finder_picked/helper_structs.rs)
ports `FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h` and its
implementation at core SDK `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. These
are the seeds, mass traces and isotope patterns that `FeatureFinderAlgorithmPicked`
and the trace fitters (`TraceFitter`, `GaussTraceFitter`, `EGHTraceFitter`)
share. This is the first step of the FeatureFinderCentroided dependency chain
in [`EARLY_TOPP_BUILD_PLAN.md`](EARLY_TOPP_BUILD_PLAN.md).

Tests: [`tests/feature_finder_picked_helper_structs.rs`](../tests/feature_finder_picked_helper_structs.rs).
Manifest: [`tests/data/feature_finder_picked_helper_structs_provenance.json`](../tests/data/feature_finder_picked_helper_structs_provenance.json).

The header is one Rust file. The source wrapper struct
`FeatureFinderAlgorithmPickedHelperStructs` has no state or members of its own,
so it becomes the module. Nothing in the module is feature-gated, and it names
only `crate::kernel`, an edge `analysis` already had, so the module graph is
unchanged.

## API mapping

Every public member of the header, and every member of the two headers that
share its typedefs, is listed.

### `Seed`

| Source | Rust |
| --- | --- |
| `Size spectrum` | `Seed::spectrum: usize` |
| `Size peak` | `Seed::peak: usize` |
| `float intensity` | `Seed::intensity: f32` |
| `bool operator<(const Seed&) const` | `Seed::is_less_intense_than(&self, &Seed) -> bool`; no `PartialOrd` (see native differences) |
| aggregate initialisation | `Seed::new(spectrum, peak, intensity)`, `Default` |

### `MassTrace`

| Source | Rust |
| --- | --- |
| `const Peak1D* max_peak = nullptr` | `MassTrace::max_peak: Option<TracePeak>` |
| `double max_rt` | `MassTrace::max_rt: f64` |
| `double theoretical_int` | `MassTrace::theoretical_int: f64` |
| `std::vector<std::pair<double, const Peak1D*>> peaks` | `MassTrace::peaks: Vec<TracePeak>`, where `TracePeak { spectrum, peak, rt, mz, intensity }` |
| `ConvexHull2D getConvexhull() const` | `MassTrace::convex_hull(&self) -> Result<ConvexHull2D>` |
| `void updateMaximum()` | `MassTrace::update_maximum(&mut self)` |
| `double getAvgMZ() const` | `MassTrace::avg_mz(&self) -> f64` |
| `bool isValid() const` | `MassTrace::is_valid(&self) -> bool` |

### `MassTraces`

| Source | Rust |
| --- | --- |
| private base `std::vector<MassTrace>` (`privvec`) | private field `traces: Vec<MassTrace>` |
| `using privvec::size` | `MassTraces::len`, plus `is_empty` |
| `using privvec::at` | `MassTraces::get` / `get_mut` -> `Option` |
| `using privvec::reserve` | `MassTraces::reserve(capacity) -> Result<()>` |
| `using privvec::push_back` | `MassTraces::push` |
| `using privvec::operator[]` | `Index<usize>` / `IndexMut<usize>` |
| `using privvec::back` | `MassTraces::last` / `last_mut` -> `Option` |
| `using privvec::clear` | `MassTraces::clear` |
| `using privvec::begin`, `end`, `iterator`, `const_iterator` | `MassTraces::iter` / `iter_mut`, `IntoIterator` for `&MassTraces` and `&mut MassTraces`, `as_slice` |
| `MassTraces()` | `MassTraces::new`, `Default` |
| `Size getPeakCount() const` | `MassTraces::peak_count(&self) -> usize` |
| `bool isValid(double seed_mz, double trace_tolerance)` | `MassTraces::is_valid(&self, seed_mz, trace_tolerance) -> bool` |
| `Size getTheoreticalmaxPosition() const` | `MassTraces::theoretical_max_position(&self) -> Result<usize>` |
| `void updateBaseline()` | `MassTraces::update_baseline(&mut self)` |
| `std::pair<double,double> getRTBounds() const` | `MassTraces::rt_bounds(&self) -> Result<(f64, f64)>` |
| `void computeIntensityProfile(std::list<std::pair<double,double>>&) const` | `MassTraces::intensity_profile(&self) -> Result<Vec<(f64, f64)>>` |
| `Size max_trace` | `MassTraces::max_trace: usize` |
| `double baseline` | `MassTraces::baseline: f64` |
| (native) | ceilings `MassTraces::MAX_TRACES`, `MAX_PEAKS`, `MAX_PROFILE_STEPS` |

### `TheoreticalIsotopePattern`

| Source | Rust |
| --- | --- |
| `std::vector<double> intensity` | `TheoreticalIsotopePattern::intensity: Vec<f64>` |
| `Size optional_begin` | `optional_begin: usize` |
| `Size optional_end` | `optional_end: usize` |
| `double max` | `max: f64` |
| `Size trimmed_left` | `trimmed_left: usize` |
| `Size size() const` | `TheoreticalIsotopePattern::len`, plus `is_empty` |

### `IsotopePattern`

| Source | Rust |
| --- | --- |
| `std::vector<SignedSize> peak` (-1 not found, -2 removed) | `IsotopePattern::peak: Vec<PatternPeak>` with `PatternPeak::{NotFound, Removed, Found(usize)}` and `PatternPeak::index` |
| `std::vector<Size> spectrum` | `spectrum: Vec<usize>` |
| `std::vector<double> intensity` | `intensity: Vec<f64>` |
| `std::vector<double> mz_score` | `mz_score: Vec<f64>` |
| `std::vector<double> theoretical_mz` | `theoretical_mz: Vec<f64>` |
| `TheoreticalIsotopePattern theoretical_pattern` | `theoretical_pattern: TheoreticalIsotopePattern` |
| `explicit IsotopePattern(Size size)` | `IsotopePattern::new(size) -> Result<IsotopePattern>`; ceiling `IsotopePattern::MAX_SIZE` |

### Not ported: `FeatureFinderDefs` and `IsotopeCluster`

None of these members is on the FeatureFinderCentroided path.

| Source | Status |
| --- | --- |
| `FEATUREFINDER/FeatureFinderDefs.h` `struct FeatureFinderDefs` | Not ported. The header has no includer. Its struct duplicates the `FeatureFinderDefs` struct defined in `FEATUREFINDER/FeatureFinderAlgorithmPicked.h`, so the two headers cannot be included together. |
| `FeatureFinderDefs::IndexPair` (`IsotopeCluster::IndexPair`) | Not ported. Neither copy is used anywhere in `source/FEATUREFINDER/`. |
| `FeatureFinderDefs::ChargedIndexSet` (`IsotopeCluster::ChargedIndexSet`) | Not ported, unused. |
| `FeatureFinderDefs::IndexSet` (`IsotopeCluster::IndexSet`) | Not ported, unused. |
| `FeatureFinderDefs::Flag { UNUSED, USED }` | Not ported, unused. |
| `FeatureFinderDefs::NoSuccessor` exception | Not ported; nothing throws it. |
| `DATASTRUCTURES/IsotopeCluster.h` `IndexPair`, `IndexSet`, `ChargedIndexSet` (with `charge`), `peaks`, `scans` | Not ported. It is reached through the two `FeatureFinderDefs` copies above and the legacy `Fitter1D.h` typedefs `IndexSet`/`ChargedIndexSet`, which the `*Fitter1D` family (unmapped) declares. The picked algorithm uses none of them. The ledger's candidate `src/processing/deisotoping.rs` holds an unrelated `Deisotoper` result struct that only shares the name. |

`FeatureFinderAlgorithm.h` (the abstract base with `Summary`) is not included by
this header and is outside this module.

## Preserved source conventions

- **Peak representation.** `MassTrace` peaks are `(spectrum index, peak index)`
  plus copied `rt: f64`, `mz: f64` and `intensity: f32`, where the source holds
  `(RT, const Peak1D*)`. `FeatureFinderAlgorithmPicked` sorts its map before
  building traces. Afterwards it writes only float data arrays and never a peak's
  m/z or intensity (source review of `FeatureFinderAlgorithmPicked.cpp`), so the
  copies read what the pointers would. The pointer identity that the class test
  checks (`max_peak == &p1_10`) becomes equality of the indices.
- **Precision.** Intensities are `f32`. `avg_mz` sums `mz * intensity` and
  `intensity` in `f64` after promoting each `f32`. `update_baseline` compares in
  `f64`. `intensity_profile` adds the promoted intensities into `f64` entries.
  This is what the source does. The oracle confirms that ten `f32` 0.1 peaks at
  one retention time give `0x3ff0000004000000`, which an `f32` accumulation would
  not.
- **First maximum wins.** `update_maximum` and `theoretical_max_position` use a
  strict `>`, so ties keep the first and a leading NaN is never replaced.
  `update_baseline` uses a strict `<` after seeding from the first peak in trace
  order.
- **Defaults.** `max_trace` is 0 (the source constructor). `IsotopePattern::new`
  sets every `peak` to `NotFound` (source `-1`) and every other vector entry to
  zero.
- **Validity rules.** `MassTrace::is_valid` needs at least 3 peaks.
  `MassTraces::is_valid` needs at least 2 traces and some trace whose average m/z
  is within the tolerance, inclusive, of the seed m/z. A NaN average never
  matches.
- **RT bounds** start from `(f64::MAX, -f64::MAX)` and skip NaN, so traces
  without peaks return those start values.
- **Profile merge.** `intensity_profile` reproduces the source's single forward
  `std::list` walk exactly. That includes unsorted traces and repeated retention
  times within a trace, where the result is not sorted (oracle case
  `boundary_profile_unsorted_duplicates`).
- **Errors where the source throws.** `theoretical_max_position` and
  `rt_bounds` return `Error::InvalidValue` on an empty collection, where the
  source throws `Exception::Precondition` in every build mode. The messages are
  the source's.

## Native differences

| Source behaviour | This port | Reason |
| --- | --- | --- |
| `Seed::operator<` has two uses in `FeatureFinderAlgorithmPicked`: `std::sort` of the seeds (`.cpp:548`), and in debug mode the order of `std::map<Seed, std::string> abort_reasons_` (`.h:153`), filled at `.cpp:1138` and written to `debug/abort_reasons.featureXML` at `.cpp:1028-1045` | `is_less_intense_than`, no `PartialOrd` | An intensity-only ordering would contradict the structural `PartialEq`. The sort order, including ties, is the caller's (B6) decision. A caller that reproduces the debug abort-reason output (B7, B10 debug 5) must key its map by `is_less_intense_than`. Seeds of equal intensity then collapse into one entry, which keeps the first seed's position and the last reason, and entries come out in ascending `f32` intensity. |
| `MassTraces::isValid` is non-`const` | `&self` | It modifies nothing. |
| `at` throws `std::out_of_range`; `back` on empty and out-of-range `operator[]` are undefined | `get`/`last` return `Option`; `Index` panics like slice indexing | Checked accessors exist for every access. |
| `reserve` fails only at the allocator's limit | `reserve` refuses more than `MAX_TRACES` (100,000) and maps allocation failure to `Error::InvalidValue` | Bounded pre-allocation; the collection is unchanged on error. |
| `IsotopePattern(Size)` always allocates | `IsotopePattern::new` refuses more than `MAX_SIZE` (100,000) isotopes | Bounded allocation. |
| `getConvexhull` accepts any coordinates into a `std::map` | `convex_hull` refuses non-finite coordinates (from `ConvexHull2D::from_points`) and more than `MAX_PEAKS` peaks | A NaN key breaks the map's ordering; infinities are not hull coordinates. |
| `ConvexHull2D::addPoint` keeps the first m/z of a retention time when a later one is equal: it skips a point its `DBoundingBox<1>` encloses, and `enlarge` replaces only on strict `<` or `>` | `convex_hull` merges the m/z range of one retention time with `f64::min` and `f64::max` in kernel `ConvexHull2D::add_points` | Known difference, not a design choice. For `-0.0` and `+0.0` m/z at one retention time, Rust leaves the sign that `min` and `max` return unspecified. An independent review fuzz against the product-SDK libOpenMS (705 cases, not retained) matched every hull bit for bit except a hand-built signed-zero case: 2 differing rows on Linux x86-64, 1 on macOS arm64. m/z is positive on the FeatureFinderCentroided path, and the retained oracle rows hold no signed zero. The fix belongs in the kernel: strict comparisons that keep the first value. |
| `computeIntensityProfile` on an empty collection dereferences `begin()` | returns an empty profile | Undefined behaviour in the source. |
| `computeIntensityProfile` never terminates when a NaN retention time meets a profile entry | returns `Error::InvalidValue` | A hang is not a result. A NaN that is only copied or appended passes through as in the source. |
| `computeIntensityProfile` fills a caller's list, documented as empty | returns a new `Vec` | That contract as a return value; an index-linked list keeps the source's insertion cost. |
| `computeIntensityProfile` has no bound | refuses more than `MAX_PEAKS` (1,000,000) total peaks or more than `MAX_PROFILE_STEPS` (100,000,000) worst-case merge steps, before allocating | Bounded work. The step bound is the sum, over every trace after the first, of all earlier peaks plus the trace's own peaks. |
| `MassTrace::max_rt`, `theoretical_int`, `MassTraces::baseline` and `TheoreticalIsotopePattern`'s scalars start uninitialised | start at zero | Rust has no indeterminate values. `update_baseline` on traces without peaks keeps the previous value, as the source does; that value is 0 unless set. |
| `peak` indices are signed `-1`/`-2`/`>= 0` | `PatternPeak` enum | The source assigns no other negative value, so no information is lost. |

## Checked boundaries and evidence

**Tier 1, executed differential.** `../oracle/feature-finder-picked-helper-structs/`
has three parts:
- `driver.cpp` replays all 14 class-test sections in source order with the
  source literals, plus 13 boundary cases;
- the driver links the product SDK's `libOpenMS.dylib` (Debug, core `4fdec46`);
  that revision's helper-struct, `ConvexHull2D` and `Peak1D` sources are
  identical to `bc9cc12`;
- `run.sh` writes 164 rows, with floats as IEEE-754 bit patterns. It was run
  twice with byte-identical output.

The retained output is `tests/data/feature_finder_picked_helper_structs_oracle.tsv`.
The test `oracle_replay_matches_the_product_sdk_output_row_for_row` rebuilds
every input and requires every row, with floats compared bit for bit. The
comparison is exact, so no tolerance is declared.

That bitwise contract applies on every platform. It is deliberately stricter
than the early-bundle plan, which asserts bit equality only on macOS arm64
against the same-platform oracle and 1e-9 relative elsewhere. The reasons:
- the members use only basic arithmetic, comparisons and `fabs`, so no libm
  result enters;
- the one expression a compiler may contract into a fused multiply-add,
  `getAvgMZ`'s `sum += mz * intensity`, gives the same bits either way for the
  recorded `mt_avg` and `inexact` cases (checked with Python `math.fma`; the
  products of `mt1` are exact), so the retained rows do not depend on whether
  the Debug oracle contracts;
- all 164 rows matched on Linux x86-64 (IBMI kim; stable and 1.85.0, with and
  without default features) and, in an independent review run, on macOS arm64.

Regenerating the oracle from a Release or FMA-contracting build, or adding
`getAvgMZ` cases, requires re-measuring before the bitwise assertion is kept.

Tier 1 covers every computational member: `getConvexhull`, `updateMaximum`,
`getAvgMZ`, both `isValid`, `getPeakCount`, `getTheoreticalmaxPosition`,
`updateBaseline`, `getRTBounds`, `computeIntensityProfile`, `Seed::operator<`,
`IsotopePattern(Size)` and `TheoreticalIsotopePattern::size`. The accessors
have tier 4 evidence only (below). No hull row holds a signed-zero m/z, where
the port has a known difference (see native differences).

The rows cover:
- `avg_mz` exactly 1000 for `mt1`, the `mt_avg` value, and an inexact
  three-peak case;
- the 18 hull points of `mt1` and the hull of an unsorted trace with repeated
  retention times;
- the 12-entry profile;
- merges with unsorted traces, duplicate retention times, prepended and appended
  peaks, and an empty first trace;
- `f64` accumulation;
- tie and NaN rules for `update_maximum`, `theoretical_max_position` and
  `update_baseline`;
- RT bounds for traces without peaks or with NaN retention times;
- NaN averages;
- the inclusive tolerance in `MassTraces::is_valid`;
- the `Seed` comparison for equal intensities, NaN and signed zeros;
- the two `Precondition` errors.

No recorded value depends on a Debug-only precondition. Source undefined
behaviour is not probed.

**Tier 3, source review.** Each `START_SECTION` of
`FeatureFinderAlgorithmPickedHelperStructs_test.cpp` (14 sections) is its own
test, with the literals and comparison semantics unchanged.
- `TEST_EQUAL` is exact equality after converting to the actual value's type:
  `getAvgMZ() == 1000`, `baseline == p2_4`, RT bounds 677.1 and 679.8.
- `TEST_REAL_SIMILAR` is ClassTest `isRealSimilar` with the default absolute
  tolerance 1e-5 and ratio 1 + 1e-5. It applies to `mt_avg` (10.4459) and to the
  12 profile entries, which are compared with the source's `f32` literal sums.
- `TEST_EXCEPTION(Exception::Precondition, ...)` becomes `Err(Error::InvalidValue(_))`.

**Tier 4, native boundaries.**
- the ceilings of `intensity_profile`, `convex_hull`, `IsotopePattern::new` and
  `reserve`, checked before allocation, with the collection unchanged on error;
- the NaN merge error, plus the NaN copy and append pass-through;
- the empty profile;
- non-finite hull coordinates;
- `clear` keeping `max_trace` and `baseline`;
- the accessors, which the oracle driver does not compare: `MassTraces::reserve`,
  `clear`, `get`/`get_mut`, `last`/`last_mut`, `iter`/`iter_mut`, `as_slice`,
  `is_empty`, the `IntoIterator` impls and `PatternPeak::index`;
- `Seed` ordering ignoring position;
- `update_baseline` keeping its value without peaks.

## C++ issue candidates (source review, not executed)

Recorded for the integrator's `OpenMS_CPP_ISSUES.md`; unconfirmed until
reproduced.

1. `MassTraces::computeIntensityProfile` dereferences `this->begin()` without
   checking for an empty collection (FeatureFinderAlgorithmPickedHelperStructs.cpp:196-198).
   This is undefined behaviour for an empty `MassTraces`. The shipped callers pass
   non-empty traces.
2. `MassTraces::computeIntensityProfile` never terminates when a NaN retention
   time meets a profile entry (FeatureFinderAlgorithmPickedHelperStructs.cpp:210-236):
   none of `>`, `<` and `==` holds, and no branch advances. Fix: treat unordered
   comparisons as an error, or advance the profile iterator.
3. `MassTraces::updateBaseline` leaves `baseline` uninitialised when there are
   traces but none holds a peak (FeatureFinderAlgorithmPickedHelperStructs.cpp:135-158;
   the constructor does not initialise it).

The empty `MassTraces::getTheoreticalmaxPosition`/`getRTBounds` preconditions
are documented and thrown in every build mode, so they are not defects.

## Ledger notes

- `FeatureFinderAlgorithmPickedHelperStructs.h`: this file covers every public
  member. It is at 100% rustdoc coverage, with tier 1 evidence for every
  computational member and tier 4 for the accessors. Suggested status is
  `complete`, with `rust` `src/analysis/feature_finder_picked/helper_structs.rs`,
  `tests` `tests/feature_finder_picked_helper_structs.rs` and this document. The
  integrator decides.
- `FeatureFinderDefs.h`: not ported; no includer, and a duplicate of the struct
  in `FeatureFinderAlgorithmPicked.h`.
- `IsotopeCluster.h`: the ledger's candidate mapping to
  `src/processing/deisotoping.rs` is wrong. No Rust file covers it; see the
  not-ported table above. Its `evidence_requires_review` status comes from the
  `pub struct IsotopeCluster` name match. A review entry can only set
  `complete`, `partial` or `native_equivalent`, so the remap needs a generator
  or naming change, not a review.
- The provenance manifest lists the context-only sources (`FeatureFinderDefs.h`,
  `IsotopeCluster.h`, `Fitter1D.h` and `FeatureFinderAlgorithmPicked.h`/`.cpp`)
  under `context_sources`. Their paths are relative to the include or source
  directory and do not start with the core repository prefix.
  `tools/core_sdk_coverage.py` counts every prefixed path string in a manifest
  as reference evidence, and this spelling keeps the unported headers
  `unmapped`. Regenerating the ledger with this manifest moves only this header
  from `unmapped` to `evidence_requires_review`. It also adds the manifest to the
  reference manifests of `ConvexHull2D.h` and `KERNEL/Peak1D.h`. Separately,
  `helper_structs.rs` joins the candidate files of `KERNEL/MassTrace.h`, because
  both declare `pub struct MassTrace`; that is a name collision, not coverage.
