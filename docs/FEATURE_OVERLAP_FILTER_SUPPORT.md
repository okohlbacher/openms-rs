# Feature overlap filter support

Native coverage of `PROCESSING/FEATURE/FeatureOverlapFilter.h` and
`PROCESSING/FEATURE/FeatureOverlapFilter.cpp`, and of the quadtree the source
bundles as `src/openms/extern/Quadtree` (`Quadtree.h`, `Box.h`, `Vector2.h`), at
Core SDK `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Work package B9-OVERLAP of
the early TOPP bundle, source mode only (decision D5).

| Artifact | Path |
| --- | --- |
| Implementation | `src/processing/feature_overlap_filter.rs`, `src/processing/feature_overlap_filter/quadtree.rs` |
| Tests | `tests/feature_overlap_filter.rs` |
| Fixtures | `tests/data/feature_overlap_filter/feature_overlap_filter_oracle.tsv`, `quadtree_grid.tsv`, `quadtree_nan.tsv` |
| Manifest | `tests/data/feature_overlap_filter_provenance.json` |
| Oracle | `../oracle/feature-overlap-filter/` (outside this repository); the C2 facts in `../oracle/featurefinder-picked/` |

Consumers in the pinned core: `FeatureFinderCentroided.cpp` (topp `174b576`,
lines 304-316) calls `mergeFAIMSFeatures(features, 5.0, 0.05)` when
`faims_merge_features` is set on FAIMS input, and so do
`Biosaur2Algorithm.cpp:361`, `FeatureFinderIdentificationAlgorithm.cpp:621` and
`FeatureFinderAlgorithmMetaboIdent.cpp:256`;
`FeatureFinderAlgorithmMetaboIdent.cpp:490` calls `filter` with a trace-level
check. The FeatureFinderCentroided FAIMS closure itself is deferred (D5, package
B11); this package ports the library only.

The module uses the edges `processing -> kernel` and `processing -> metadata`,
which exist, and adds `processing -> concept` for the `FAIMS_CV` key. That edge
closes no cycle: `concept` names no other top-level module
(`tools/check_module_cycles.py` reports it as a new acyclic edge).

## Quadtree: crates first

The source takes its quadtree from a vendored third-party library, so a crate
was looked for before porting it. A crate qualifies only if its query order
equals `extern/Quadtree`: 16 values per leaf before a split, no split at depth 8,
`f32` boxes, strict intersection, a node's own values reported before its
children, children north-west, north-east, south-west, south-east, and the
source's split, removal and merge rules. The order is observable, because the
filter's callbacks change state in it.

| Crate (version, crates.io 2026-09-14) | Why it does not qualify |
| --- | --- |
| `aabb-quadtree` 0.2.0 (2019, no MSRV declared, `euclid`) | Inclusive intersection plus `contains` of either corner and an epsilon closeness test; an item spanning a node's midpoint goes to `in_all`, any other item into every intersecting child (duplicates); `query` sorts results by item id and removes duplicates; near-duplicate boxes are dropped unless `allow_duplicates`; panics when an insert fails |
| `quadtree-f32` 0.5.0 | Four items per node, inclusive `overlaps_rect`, items kept in a `BTreeMap`, query order not the insertion order; no MSRV declared |
| `quadtree` 0.5.0 (benbaarber) | A point quadtree over `glam::Vec2` with a Barnes-Hut variant; stores points, not boxes |
| `quadtree_rs` 0.1.3 | Integer (`num::PrimInt`) anchors and region sizes; no floating-point boxes |
| `rectutils` 0.7.0 (MSRV 1.72) | Built once from all entries (no incremental `add`/`remove`), depth up to 64, an entry copied into every intersecting leaf |
| `spart` 0.6.1 (MSRV 1.85) | Point-based trees with a capacity split; no box values and a different split rule |
| `rstar` 0.13 (MSRV 1.85) | An R*-tree: different partitioning and a different query order |

None qualifies, so `quadtree.rs` ports the vendored header (MIT, Pierre Vigier
2019). No dependency is requested. The port is 1:1 in every rule that decides
placement and order, and it is checked against the pinned header executed
(below).

## API mapping

### `FeatureOverlapFilter.h`

| C++ | Rust |
| --- | --- |
| `enum class FeatureOverlapMode { CONVEX_HULL, TRACE_LEVEL, CENTROID_BASED }` | `FeatureOverlapMode::{ConvexHull, TraceLevel, CentroidBased}`, `ConvexHull` the `Default` as documented |
| `enum class MergeIntensityMode { SUM, MAX }` | `MergeIntensityMode::{Sum, Max}`, `Sum` the `Default` |
| `struct CentroidTolerances` and its defaults 5.0, 0.05, `true`, `false` | `CentroidTolerances` with public fields and the same `Default` |
| `class FeatureOverlapFilter` | `FeatureOverlapFilter`, a unit struct holding the associated functions |
| `using OverlapMode = FeatureOverlapMode` | `feature_overlap_filter::OverlapMode` (a module-level alias; Rust has no inherent associated type aliases) |
| `static void filter(FeatureMap&, comparator = quality >, callback = true, bool check_overlap_at_trace_level = true)` | `FeatureOverlapFilter::filter(&mut FeatureMap, comparator, on_overlap, bool) -> Result<()>`; the default arguments are `FeatureOverlapFilter::higher_overall_quality`, `FeatureOverlapFilter::always_overlapping` and `true`, passed explicitly |
| `static void filter(FeatureMap&, comparator, callback, FeatureOverlapMode, const CentroidTolerances& = CentroidTolerances())` | `FeatureOverlapFilter::filter_with_mode(.., FeatureOverlapMode, &CentroidTolerances) -> Result<()>` |
| `std::function<bool(const Feature&, const Feature&)>` comparator | `FnMut(&Feature, &Feature) -> bool` |
| `std::function<bool(Feature&, Feature&)>` callback | `FnMut(&mut Feature, &mut Feature) -> bool`; a callback that fails goes to `FeatureOverlapFilter::filter_with_fallible_callback` (`FnMut(..) -> Result<bool>`), the counterpart of a throwing callback |
| `static std::function<bool(Feature&, Feature&)> createFAIMSMergeCallback(MergeIntensityMode = SUM, bool write_meta_values = true)` | `FeatureOverlapFilter::create_faims_merge_callback(MergeIntensityMode, bool) -> impl FnMut(&mut Feature, &mut Feature) -> Result<bool>`, built on `FaimsMergeCallback::merge` |
| the lambda it returns | `FaimsMergeCallback { intensity_mode, write_meta_values }` with `merge(&self, &mut Feature, &Feature) -> Result<bool>` |
| `static void mergeOverlappingFeatures(FeatureMap&, double = 5.0, double = 0.05, bool = true, bool = false, MergeIntensityMode = SUM, bool = true)` | `FeatureOverlapFilter::merge_overlapping_features(&mut FeatureMap, f64, f64, bool, bool, MergeIntensityMode, bool) -> Result<()>` |
| `static void mergeFAIMSFeatures(FeatureMap&, double = 5.0, double = 0.05)` | `FeatureOverlapFilter::merge_faims_features(&mut FeatureMap, f64, f64) -> Result<()>` |
| meta value names `"merged_centroid_rts"`, `"merged_centroid_mzs"`, `"merged_centroid_IMs"`, `"FAIMS_merge_count"` | `MERGED_CENTROID_RTS`, `MERGED_CENTROID_MZS`, `MERGED_CENTROID_IMS`, `FAIMS_MERGE_COUNT` |
| used: `Constants::UserParam::FAIMS_CV` | `concept::constants::user_param::FAIMS_CV` (ported earlier) |
| file-local `struct MassTraceBounds`, `FeatureBoundsMap`, `getFeatureBounds`, `hasOverlappingBounds`, `tracesOverlap` (`.cpp:20-120`) | private `TraceBounds`, `trace_bounds` (keyed by unique ID in a `BTreeMap`) and `Plan::overlaps`/`Plan::traces_of`; `MassTraceBounds::sub_index` is written but never read in the source and is not kept |
| the `getBox` lambdas (`.cpp:143-169`) | private `Plan::feature_box` |
| `Feature::getConvexHull().getBoundingBox()` as used there | private `SourceBox::of_feature`, which keeps `DBoundingBox::enlarge` and `isEmpty` exactly |
| `FeatureMap::updateRanges()` / `getMinMZ()` ... as used there | private `Range` in `Plan::new` |
| new, no source counterpart | `FeatureOverlapFilter::MAX_FEATURES`, `FeatureOverlapFilter::MAX_CANDIDATE_VISITS` |

### `extern/Quadtree`

| C++ | Rust |
| --- | --- |
| `template<typename T> class Vector2` with `x`, `y`, `operator+=`, `operator/=`, `operator+`, `operator/` | `quadtree::Vector2` (`f32`), `Add` and `Div<f32>`; the compound forms are covered by the operators |
| `template<typename T> class Box` with `left`, `top`, `width`, `height` | `quadtree::QuadBox` (`f32`) with the same public fields; renamed so it does not shadow `std::boxed::Box` |
| `Box(Left, Top, Width, Height)`, `Box(position, size)` | `QuadBox::new`, `QuadBox::from_position_size` |
| `getRight`, `getBottom`, `getTopLeft`, `getCenter`, `getSize` | `right`, `bottom`, `top_left`, `center`, `size` |
| `contains`, `intersects` | `contains`, `intersects` |
| `template<typename T, typename GetBox, typename Equal = std::equal_to<T>, typename Float = float> class Quadtree` | `quadtree::Quadtree<T>`; `Float` is `f32`, `Equal` is `PartialEq`, and `GetBox` is an argument of each operation |
| `Quadtree(const Box& box, const GetBox&, const Equal&)` | `Quadtree::new(QuadBox)` |
| `void add(const T&)` | `add(T, &G) -> Result<()>` |
| `void remove(const T&)` | `remove(&T, &G) -> Result<()>` |
| `std::vector<T> query(const Box&) const` | `query(QuadBox, &G) -> Vec<T>`; `query_into` reuses a buffer |
| `std::vector<std::pair<T, T>> findAllIntersections() const` | `find_all_intersections(&G) -> Result<Vec<(T, T)>>` |
| `Box getBox() const` | `bounds()` |
| `Threshold = 16`, `MaxDepth = 8` | `Quadtree::THRESHOLD`, `Quadtree::MAX_DEPTH` |
| private `Node`, `isLeaf`, `computeBox`, `getQuadrant`, `split`, `removeValue`, `tryMerge`, the recursive helpers | private arena nodes, `child_box`, `quadrant`, `split`, `try_merge`; `removeValue` is the `swap_remove` in `remove` |
| new, no source counterpart | `len`, `is_empty`, `MAX_VALUES`, `MAX_INTERSECTIONS` |

## Preserved source conventions

Each item is covered by the executed oracle unless marked otherwise.

- **Order.** Features are sorted stably by the comparator; candidates are handed
  to the callback in quadtree query order; survivors keep the sorted order. The
  merge functions sort by intensity, descending, ties in input order
  (`bound_intensity_ties_stable`).
- **`f32` boxes.** Centroid boxes are `float(float(mz) - tol)` wide
  `float(2 * tol)`; hull boxes convert the `double` bounding box corners and
  extents separately; the extent is `float(float(min) - tol) - 1` wide
  `max - min + 2` in `float`. Every operation is in the source's order and
  precision (`bound_float_rt_rounded_before_tolerance`,
  `bound_float_mz_rounded_before_tolerance`, `bound_extent_rounding_found`).
- **Strict intersection before the exact test.** A pair whose boxes do not
  intersect strictly is never tested: a zero tolerance merges nothing, not even
  identical features (`bound_zero_tolerance`); a tolerance below the float
  spacing merges nothing (`bound_tiny_tolerance`); at 2^24 s the boxes collapse
  (`bound_float_collapse_rt`, `bound_merge_float_collapse`). The exact test is
  inclusive in `double` (`bound_exact_tolerance`, `bound_just_outside`,
  `bound_float_same_box_exact_rejects`).
- **Boxes at query time.** Boxes are recomputed from the features at every
  query, as the source's pointer-based `getBox` does, so a callback that moves a
  feature changes later queries (`order_centroid_loose_rule3`, whose callback
  moves the survivor by 0.5 s).
- **Extent from `updateRanges`.** Centroids, then each feature's hull bounding
  box unless `DBoundingBox::isEmpty`, which also calls a box of zero width or
  height empty (`edge_zero_extent_hull_many`, `edge_zero_extent_hull_within_margin`).
- **Removal by unique ID.** Every feature whose unique ID was marked is erased,
  and a feature is skipped as a querier when its ID is marked. Features sharing
  ID 0 go together (`c2_uid0_wipe`, `edge_trace_uid0_pooling`); in the trace mode
  they also pool their trace bounds.
- **Double merge.** A removed feature can still be merged into a later survivor
  (`c2_three_cvs`: 1000/900/800 become 1900/1700).
- **`mergeFAIMSFeatures`.** Nothing happens without a `FAIMS_CV`; with exactly
  one FAIMS feature the map is only reordered (`edge_single_faims`); the survivor
  of a first merge loses `FAIMS_CV` and absorbs nothing more; the map holds the
  FAIMS features in merge order, then the others in input order
  (`edge_faims_reorder`); map-level data is kept (`edge_map_metadata_preserved`);
  an integer `FAIMS_CV` is converted to `double` (`edge_faims_cv_int`).
- **Merge callback.** Existing merged lists are extended; `FAIMS_CV` is removed
  from the survivor only when it has no `merged_centroid_IMs` yet
  (`edge_existing_merged_lists`); `FAIMS_merge_count` is an integer; `MAX` keeps
  the first intensity on a tie, as `std::max`.
- **Trace bounds.** `getFeatureBounds` reads feature hull `i` for subordinate
  `i`, the m/z bounds from its first and last outline points, and the retention
  times from the subordinate's first hull: the start is the first outline point
  with m/z above zero; the end walks the outline backwards and, for any hull
  whose first scan has an m/z extent, stops at that first scan again. Traces
  therefore compare by their start (`edge_trace_bounds_collapse_trace` keeps two
  features whose traces overlap for 1.5 s; `edge_trace_bounds_collapse_hull`
  merges them). A trace whose computed start lies after its end is skipped
  (`.cpp:84-87`): an outline without m/z above zero whose first scan is a single
  point, or one whose first scan has m/z 0 at its lower edge, gives a start at
  the second scan and an end at the first; such a trace takes no part in the
  overlap test (`edge_trace_inverted_bounds_skipped` keeps a feature whose only
  overlapping traces are of this kind, `edge_trace_inverted_bounds_control`
  removes it when one of them is regular).
- **Multi-hull features.** The hull box of a feature with several hulls is the
  box `DBoundingBox::enlarge` builds from every hull's corners
  (`edge_multi_hull_bounding_box`).
- **Quadtree rules.** Threshold, depth limit, strict quadrant tests (a box
  touching a centre line from the west stays west only when its right edge is
  strictly left of it), query order, `swap_remove` on removal, merging of
  children after a leaf removal, and the pair order of `findAllIntersections`
  (`quadtree_grid.tsv`); NaN and infinite boxes, including the source's
  `(+inf, +inf, -inf, -inf)` sentinel, follow IEEE comparisons
  (`quadtree_nan.tsv`).

## Native differences

1. **Atomic.** An error leaves the map exactly as it was. The source sorts the map
   before `getFeatureBounds` can throw (`edge_trace_missing_sub_hull` leaves the
   C++ map sorted), and `mergeFAIMSFeatures` has moved the features into its two
   temporary maps when its callback throws, so the C++ map is left holding
   features stripped of their metadata (`edge_faims_cv_empty`). The merge
   functions keep a journal of the values their callbacks overwrote and undo
   it. The journal records each field of each feature once, with its value
   before the call, so it holds at most six entries per feature however many
   features a survivor absorbs (the source itself keeps only the current
   merged lists); the
   generic filters copy the features first when a failure after a callback is
   possible (trace mode, `require_same_im`, more than
   `MAX_CANDIDATE_VISITS` possible candidates), and
   `filter_with_fallible_callback` always does.
2. **Undefined source behaviour is refused**, before anything changes:
   - hull modes: a feature without a convex hull, or with an empty one, is
     `Error::MissingInformation`. The source converts the `±DBL_MAX` sentinel box
     to `float` (undefined before C++23); the Debug library aborts on the
     quadtree's containment assertion, and the Release replica silently ignores
     the feature (`edge_hull_less_convex_hull`, `debug_only`). In the trace mode
     this refusal also covers the empty feature hull a subordinate's m/z bounds
     would be read from (the source reads `front()` of an empty vector);
   - trace mode: a subordinate without a matching feature hull (the source reads
     past the hull vector) is `Error::InvalidValue`; a candidate pair of
     which one feature has no trace bounds is `Error::InvalidValue` when it is
     reached (the source dereferences `std::map::end()`);
   - a box or extent that does not fit `f32` is `Error::InvalidValue`;
   - a `FAIMS_CV` that is read and is not numeric, or an existing merged list that
     is not a float list, is `Error::InvalidValue`. The source throws
     `ConversionError` for an empty value or a non-list, but reads the `double`
     member of the value's union for a string or a list
     (`edge_faims_cv_empty`, `edge_merged_list_wrong_type`).
3. **Validated input.** Every feature must pass `Feature::validate` (finite
   coordinates, intensity, quality and meta values), nonzero unique IDs must be
   distinct (a `FeatureMap` invariant; the source would remove every feature with
   a marked ID), and centroid-mode tolerances must be finite and nonnegative
   (`Box.h` requires a positive size and the source checks nothing).
   `merge_faims_features` validates only the FAIMS features, and only when there
   are some, as the source reads nothing otherwise.
4. **Errors instead of silent values.** A summed intensity that is not a finite
   `f32` is `Error::InvalidValue`; the source stores infinity.
5. **Bounds.** `MAX_FEATURES` (`FeatureMap::MAX_ITEMS`, ten million) and
   `MAX_CANDIDATE_VISITS` (2^32 candidates per call, the querier included);
   `Quadtree::MAX_VALUES` and `Quadtree::MAX_INTERSECTIONS`. The source is
   unbounded and its worst case is quadratic.
6. **Debug assertions are not emulated.** The port follows a Release build: a
   value box outside the root, or a query box that misses it, is handled by the
   quadrant and intersection tests. The four oracle cases the Debug library
   aborts on are compared with the Release replica
   (`bound_float_collapse_rt`, `bound_merge_float_collapse`,
   `edge_zero_extent_hull_outside_margin`; `edge_hull_less_convex_hull` is
   refused by difference 2).
7. **Callbacks.** The comparator may be called twice per comparison (once in the
   source's merge sort); a comparator that is not a strict weak ordering is
   undefined in the source and may panic in Rust's sort. A Rust callback cannot
   throw: a failing callback returns `Err` through
   `filter_with_fallible_callback`.
8. **The quadtree takes the box function per call** and removal of an absent
   value is an error (the source asserts, then writes through `end()`).
9. **Serial**, as the source; no `parallel` path, so no thread-count contract
   applies.

## Checked boundaries and evidence

In `tests/feature_overlap_filter.rs` (33 tests) and the unit tests of
`src/processing/feature_overlap_filter.rs` (3 tests).

**Tier 3, class test.** The 14 START_SECTIONs of `FeatureOverlapFilter_test.cpp`,
transcribed with their literals, comparisons and `TEST_REAL_SIMILAR` default
tolerance unchanged; see the accounting below.

**Tier 1, executed differential.** `../oracle/feature-overlap-filter/` links the
product-sdk libOpenMS (Debug, core `4fdec46`; the installed
`FeatureOverlapFilter.h` equals the pin, and the traced `PROCESSING/FEATURE`,
`KERNEL` and `DATASTRUCTURES` sources are unchanged between `4fdec46` and the pin
per the C2 manifest). It runs 78 cases, each in its own process and twice, and
records inputs, every callback invocation, exceptions and outputs with doubles
and floats as their bits. `oracle_cases_replay_bit_for_bit` rebuilds every input
through the port's model (the round trip of every input is checked first), runs
the same call, and compares:

- the callback sequence `(best, other, return)`, in order (16,917 invocations);
- every output feature bit for bit: position, intensity, quality, charge,
  unique ID, every meta value, every hull and subordinate hull;
- the map-level fields for the `mergeFAIMSFeatures` cases;
- for the six cases where C++ throws, the error class
  (`InvalidRange`, `MissingInformation`, `ConversionError` → `InvalidValue`) and
  an unchanged map.

71 cases compare equal, six are errors matched to the C++ exception, and one is
refused natively (`edge_hull_less_convex_hull`). The families: the 14 class-test
sections as 16 cases, two sections calling twice (`class_*`); the three C2 facts
(`c2_*`); 19 random-map cases of 240 or 700 features in all three modes with
callbacks that keep everything, remove everything, remove by a rule, and move the
survivor (`order_*`, which split the tree over several levels); 16 float-boundary
cases (`bound_*`); and 24 edge cases (`edge_*`).

`c2_*` tests also transcribe the C2 literals directly
(`../oracle/featurefinder-picked/results/omp1/faims_facts.jsonl`, sha256
`5451c697...`): all FAIMS features with unique ID 0 are removed; the valid-ID
control merges to 1500 with `merged_centroid_IMs` `[-45, -60]`; three voltages
give 1900 and 1700.

**Tier 2, executed probes.** The pinned `FeatureOverlapFilter.cpp` and
`extern/Quadtree` headers compiled into a replica with `-DNDEBUG` agree byte for
byte with the library on all 74 cases the library completes, which also
confirms that the library behaves as the pinned source; three of the four
`debug_only` cases are compared with the replica, and the fourth is refused
(native difference 2). The pinned quadtree header, driven directly with
359 boxes through 359 adds, 145 queries, 186 removals and two intersection lists,
gives the same records with and without assertions (`quadtree_grid.tsv`), and a
NaN/infinity scenario runs with `-DNDEBUG` (`quadtree_nan.tsv`). Both replay
bit for bit.

**Sensitivity.** Four deliberate mutations of the port were each caught by the
replay before being reverted: reversing the child order of a query; computing a
centroid box from `rt - tol` in `f64` without the source's intermediate `float`;
computing the extent without that intermediate `float`; and counting
zero-extent hull boxes in the extent. The last two needed dedicated cases
(`bound_extent_rounding_found`, found by a float32 simulation kept in the oracle's
`tools/extent_search2.py`, and `edge_zero_extent_hull_many`). Two later
mutations were caught the same way: deleting the skip of a trace whose start
lies after its end fails the replay of `edge_trace_inverted_bounds_skipped`
(it had survived before that case existed), and journaling every change instead
of each field's first fails the three journal unit tests.

**Tier 4, native contracts.** Atomicity after a merge has happened
(`faims_merge_error_after_a_merge_restores_the_map`), after a callback has run in
trace mode (`trace_mode_candidate_without_bounds_restores_the_map_after_a_callback`),
after a fallible callback fails, on `f32` overflow, and after a survivor has
changed every journaled field several times
(`merge_error_after_repeated_changes_of_the_same_fields_restores_the_map`); a
journal of at most one entry per field per feature, with 2,000 co-located
features merged into one and, on crafted input, 500 FAIMS features that all
carry `merged_centroid_IMs` (unit tests
`merge_journal_records_each_field_once_however_many_merges`,
`faims_journal_records_each_field_once_on_crafted_input`,
`journal_rollback_restores_the_values_before_the_run`); each refusal of
difference 2 and 3 with the map unchanged; `merge_faims_features` leaving a map
without FAIMS features untouched even when its features are invalid;
`create_faims_merge_callback` through `filter_with_fallible_callback` equal to
`merge_overlapping_features` on a 240-feature map; `FaimsMergeCallback::merge`
leaving the survivor unchanged on a conversion error; quadtree constants and box
predicates; removal of an absent value.

**Not tested.** `MAX_CANDIDATE_VISITS` and the value and intersection ceilings
are not reached by a test: reaching 2^32 candidates takes tens of seconds even in
an optimised build.

**Platform.** The oracle ran on macOS arm64 (AppleClang 21, `-ffp-contract=off`
as libOpenMS); the replay passes bit for bit on the Linux gate host, on Rust
1.85.0 and on the current stable toolchain. The port uses only IEEE 754 basic
operations and `f64` to `f32` conversions, which Rust defines as
round-to-nearest-even on every target.

## Class-test section accounting

| `FeatureOverlapFilter_test.cpp` section | Rust test |
| --- | --- |
| `(Filter FeatureMap)` | `section_filter_feature_map` |
| `mergeOverlappingFeatures - basic merging with SUM intensity` | `section_merge_overlapping_features_basic_sum` |
| `mergeOverlappingFeatures - MAX intensity mode` | `section_merge_overlapping_features_max_intensity` |
| `mergeOverlappingFeatures - require_same_charge` | `section_merge_overlapping_features_require_same_charge` |
| `mergeOverlappingFeatures - require_same_im with FAIMS_CV` | `section_merge_overlapping_features_require_same_im_with_faims_cv` |
| `mergeOverlappingFeatures - require_same_im with same FAIMS_CV` | `section_merge_overlapping_features_require_same_im_with_same_faims_cv` |
| `mergeOverlappingFeatures - features without FAIMS_CV` | `section_merge_overlapping_features_without_faims_cv` |
| `mergeOverlappingFeatures - mixed FAIMS_CV presence with require_same_im` | `section_merge_overlapping_features_mixed_faims_cv_presence` |
| `mergeOverlappingFeatures - write_meta_values=false` | `section_merge_overlapping_features_write_meta_values_false` |
| `mergeOverlappingFeatures - no merge when outside tolerance` | `section_merge_overlapping_features_no_merge_outside_tolerance` |
| `mergeOverlappingFeatures - multiple features merging` | `section_merge_overlapping_features_multiple_features` |
| `mergeFAIMSFeatures - only merges features with different FAIMS_CV` | `section_merge_faims_features_only_merges_different_faims_cv` |
| `mergeFAIMSFeatures - does NOT merge features with same FAIMS_CV` | `section_merge_faims_features_does_not_merge_same_faims_cv` |
| `mergeFAIMSFeatures - no-op on non-FAIMS data` | `section_merge_faims_features_no_op_on_non_faims_data` |

14 sections, 14 ported, none unaccounted. The same 14 run in the oracle as
`class_*` and replay bit for bit. The Quadtree has no class test in the core.

## Source defects observed

Candidates for `OpenMS_CPP_ISSUES.md`, with executed evidence in the oracle:

1. `mergeFAIMSFeatures` on features without unique IDs (as
   `FeatureFinderAlgorithmPicked` returns them before the tool assigns IDs)
   removes every FAIMS feature, because removal is keyed by unique ID
   (`FeatureOverlapFilter.cpp:269-281`; `c2_uid0_wipe`, also C2).
2. `filter` does not skip candidates that were already removed, so a removed
   feature is merged again into a later survivor and its intensity counted twice;
   in `mergeFAIMSFeatures` a survivor absorbs at most one feature because its
   `FAIMS_CV` is removed after the first merge (`c2_three_cvs`).
3. The quadtree boxes are `float` and intersect strictly, while the tolerance
   test is inclusive `double`: a zero tolerance, a tolerance below the float
   spacing, or coordinates near 2^24 never merge features that are within the
   tolerances (`bound_zero_tolerance`, `bound_tiny_tolerance`,
   `bound_float_collapse_rt`).
4. `getFeatureBounds` takes a trace's retention-time end from the end of the hull
   outline, which returns to the first scan, so trace-level overlap compares only
   start times (`.cpp:69-83`; `edge_trace_bounds_collapse_trace`).
5. `FeatureMap::updateRanges` skips zero-width hull boxes
   (`DBoundingBox::isEmpty`), so the quadtree extent can exclude a feature's box
   and a Debug build aborts on a valid map (`edge_zero_extent_hull_outside_margin`);
   a feature without a hull aborts the same way in the hull modes
   (`edge_hull_less_convex_hull`).
6. Exception safety: `filter` throws `MissingInformation` after sorting the map,
   and an exception inside the `mergeFAIMSFeatures` callback leaves the caller's
   map with moved-from features (`edge_trace_missing_sub_hull`,
   `edge_faims_cv_empty`).
7. Undefined behaviour on reachable inputs: `tracesOverlap` dereferences
   `std::map::end()` for a feature without trace bounds, `getFeatureBounds`
   indexes past the hull vector when a feature has more subordinates than hulls,
   and `(double)` on a string `FAIMS_CV` reads the wrong union member (source
   review; not executed).

## Deferrals

- The FeatureFinderCentroided FAIMS closure (`faims_merge_features`) is B11's,
  after decision D5. No corrected mode exists here: the source mode is the only
  mode.
- Ledger status: `FeatureOverlapFilter.h` can be recorded as ported with the
  evidence above; the promotion is the integrator's.
