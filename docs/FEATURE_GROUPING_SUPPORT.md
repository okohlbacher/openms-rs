# Feature grouping support (`FeatureDistance`, `BaseGroupFinder`, `StablePairFinder`, `FeatureGroupingAlgorithm*`)

Source revision: OpenMS4-core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

| Source header | Rust | Status |
|---|---|---|
| `ANALYSIS/MAPMATCHING/FeatureDistance.h` | `src/analysis/mapmatching/feature_distance.rs` | complete |
| `ANALYSIS/MAPMATCHING/BaseGroupFinder.h` | `src/analysis/mapmatching/base_group_finder.rs` | complete |
| `ANALYSIS/MAPMATCHING/StablePairFinder.h` | `src/analysis/mapmatching/stable_pair_finder.rs` | complete |
| `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h` | `src/analysis/mapmatching/feature_grouping_algorithm.rs` | complete |
| `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h` | `src/analysis/mapmatching/feature_grouping_algorithm_unlabeled.rs` | complete |
| `ANALYSIS/MAPMATCHING/QTClusterFinder.h` | `src/analysis/mapmatching/qt_cluster_finder.rs` | partial (tie order of equal-quality clusters) |
| `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h` | `src/analysis/mapmatching/feature_grouping_algorithm_qt.rs` | complete (inherits the finder's tie order) |
| `FeatureGroupingAlgorithmKD.h` | none | blocked on `DATASTRUCTURES/KDTree.h` |
| `LabeledPairFinder.h`, `FeatureGroupingAlgorithmLabeled.h` | none | not ported |
| `WNetMatcher.h`, `FeatureGroupingAlgorithmWNet.h` | none | not ported |
| `FeatureMapping.h`, `PipEchoAlgorithm.h` | none | not ported |

Parameter handling follows `src/processing/peak_picking.rs`: `defaults()`,
`from_param`, `from_param_with_warnings`, `to_param`, validated through
`DefaultParamHandler`. Parameter-tree construction shares
`src/analysis/mapmatching/params.rs` (crate-private).

## API mapping

### FeatureDistance

| Source | Rust |
|---|---|
| `static const double infinity` | `FeatureDistance::INFINITY` |
| `FeatureDistance(double max_intensity = 1.0, bool force_constraints = false)` | `new(max_intensity, force_constraints)`, `Default` (1.0, false) |
| `~FeatureDistance()` | drop |
| `operator=` | `Clone` (all state is copied, including `max_intensity` and `force_constraints`; the source recomputes derived members, which here are computed per call) |
| `std::pair<bool,double> operator()(left, right)` | `distance(&left, &right) -> Result<(bool, f64)>` |
| `DefaultParamHandler` `getDefaults` / `setParameters` / `getParameters` | `defaults()`, `from_param`, `from_param_with_warnings`, `with_param(max_intensity, force, &Param)`, `to_param()` |
| protected `DistanceParams_` (`max_difference`, `exponent`, `weight`, `norm_factor`, `max_diff_ppm`, `relevant`) | `DistanceComponent` (`max_difference`, `exponent`, `weight`); `norm_factor` and `relevant` derived per call; `max_diff_ppm` is `FeatureDistance::mz_ppm` |
| protected `updateMembers_`, `distance_` | `from_complete_param`, `DistanceComponent::distance` (private) |
| protected `params_rt_`, `params_mz_`, `params_intensity_`, `total_weight_reciprocal_`, `max_intensity_`, `ignore_charge_`, `ignore_adduct_`, `force_constraints_`, `log_transform_` | public fields `rt`, `mz`, `intensity_exponent` / `intensity_weight` (with `intensity_component()`), `max_intensity`, `ignore_charge`, `ignore_adduct`, `force_constraints`, `log_transform`; the reciprocal is computed per call |

### BaseGroupFinder

| Source | Rust |
|---|---|
| abstract class (`DefaultParamHandler`, `ProgressLogger`) | trait `BaseGroupFinder` |
| `virtual void run(const vector<ConsensusMap>&, ConsensusMap&) = 0` | `fn run(&self, &[ConsensusMap], &mut ConsensusMap) -> Result<()>` |
| protected `checkIds_(maps)` | `check_ids(&[ConsensusMap]) -> Result<()>` |
| private copy and assignment | not applicable |

### StablePairFinder

| Source | Rust |
|---|---|
| `StablePairFinder()`, destructor | `Default`, `from_param`, `from_param_with_warnings`, `defaults`, `to_param` |
| `run(input_maps, result_map)` | `BaseGroupFinder::run`; `run_with_max_intensity(input, result, max_intensity)` for callers that model a stale range; `max_intensity(&[&ConsensusMap])` |
| protected `updateMembers_` | `from_complete_param` |
| protected `compatibleIDs_(feat1, feat2)` | `compatible_ids(&feat1, &feat2) -> Result<bool>` |
| protected `getBestHitSequence_(pep_id)` | `best_hit_sequence(&pep_id) -> Result<Option<String>>` |
| protected `second_nearest_gap_`, `use_IDs_` | public fields `second_nearest_gap`, `use_identifications`; the embedded distance parameters are the field `distance` |
| protected `enum { RT, MZ }` | not needed |

### FeatureGroupingAlgorithm

| Source | Rust |
|---|---|
| abstract class, constructor, destructor, private copy and assignment | trait `FeatureGroupingAlgorithm` |
| `virtual void group(const vector<FeatureMap>&, ConsensusMap&) = 0` | `fn group(&mut self, &[FeatureMap], &mut ConsensusMap) -> Result<()>` |
| `virtual void group(const vector<ConsensusMap>&, ConsensusMap&)` | provided method `group_consensus` (converts keeping unique IDs, then forwards) |
| `void transferSubelements(maps, out) const` | `transfer_subelements(&[ConsensusMap], &mut ConsensusMap) -> Result<()>` |
| protected template `postprocess_(maps, out)` | `postprocess(&[InputIdentifications], &mut ConsensusMap)`, with `feature_map_identifications` and `consensus_map_identifications` |

### FeatureGroupingAlgorithmUnlabeled

| Source | Rust |
|---|---|
| constructor, destructor | `Default`, `from_param`, `from_param_with_warnings`, `defaults`, `to_param` |
| `template setReference(int map_id, const MapType&)` | `set_reference(u64, &FeatureMap)`, `set_reference_peak_map(u64, &MSExperiment, &mut UniqueIdGenerator)` |
| `ConsensusMap& getResultMap()` | `result_map()`, `result_map_mut()` |
| `group(maps, out)` | `FeatureGroupingAlgorithm::group` |
| `addToGroup(int map_id, const FeatureMap&)` | `add_to_group(u64, &FeatureMap) -> Result<()>` |
| private `pairfinder_input_` | private field, plus `reference_range_current` (see below) |

### QTClusterFinder

| Source | Rust |
|---|---|
| typedefs `PairDistances`, `ElementMapping`, `Heap`, `Grid` | internal: element mapping `BTreeMap<(map, feature), BTreeSet<id>>`, a versioned `BinaryHeap`, `HashGrid<usize>`; `PairDistances` is unused by the source too |
| `QTClusterFinder()`, destructor | `Default`, `from_param`, `from_param_with_warnings`, `defaults`, `to_param` |
| `run(const vector<ConsensusMap>&, ConsensusMap&)` | `BaseGroupFinder::run` |
| `run(const vector<FeatureMap>&, ConsensusMap&)` | `run_feature_maps`; both through the generic `run_maps<M: QtInputMap>` (source template `run_`) |
| private `run_internal_`, `setParameters_`, `computeClustering_`, `addClusterElements_`, `makeConsensusFeature_`, `removeFromElementMapping_`, `createConsensusFeature_`, `updateClustering_`, `getDistance_`, `distIsOutlier_` | private functions of the same names and bodies |
| private `num_maps_`, `use_IDs_`, `min_nr_diffs_per_bin_`, `min_score_`, `noID_penalty_`, `max_diff_rt_`, `max_diff_mz_`, `nr_partitions_`, `feature_distance_`, `already_used_` | public parameter fields `use_identifications`, `min_nr_diffs_per_bin`, `min_id_score_for_tol_calc`, `no_id_penalty`, `nr_partitions`, `distance`; per-run state is local |
| private `bin_tolerances_` | private field, read through `bin_tolerances()`; never cleared, as the source member |
| protected `enum { RT, MZ }` | not needed |

### FeatureGroupingAlgorithmQT

| Source | Rust |
|---|---|
| constructor, destructor, private copy and assignment | `Default`, `from_param`, `from_param_with_warnings`, `defaults`, `to_param` |
| `group(const vector<FeatureMap>&, ConsensusMap&)` | `FeatureGroupingAlgorithm::group` |
| `group(const vector<ConsensusMap>&, ConsensusMap&)` | `FeatureGroupingAlgorithm::group_consensus` (overridden: no conversion) |
| private template `group_` | private `group_impl` |

## Preserved source conventions

- `FeatureDistance`: charges must agree unless one is 0 or `ignore_charge` is
  set; adduct formulas (`dc_charge_adducts`) must agree when both features are
  annotated and `ignore_adduct` is off; a ppm m/z tolerance is relative to the
  **left** feature; an exceeded maximum invalidates the pair, or returns
  infinity under `force_constraints`; a component with weight or exponent 0 is
  dropped from the distance and from the total weight; exponents 1 and 2 use
  the source's multiplication shortcuts; the result is
  `(rt + mz + intensity) * (1 / total_weight)`; the non-log intensity
  difference is taken in `f32` (`fabs(float - float)`).
- `StablePairFinder`: nearest and second-nearest bookkeeping exactly as the
  source (an invalid best match only updates the second-nearest distance);
  the quality arithmetic, including `getQuality() * (size - 1)` evaluated in
  `f32`; merged features take the union of handles and append both features'
  identifications; singletons with fewer than two handles get quality 0; the
  result is stably sorted by m/z; `result` keeps its metadata, headers and
  identifications (`clear(false)`).
- The normalising intensity comes from *stored* ranges in the source. They are
  current right after `MapConversion`, and a pair-finder result never has them
  updated. `FeatureGroupingAlgorithmUnlabeled` therefore uses both maps for its
  first pairing and only the added map afterwards, in `group` and in the
  incremental path. `TOPP_FeatureLinkerUnlabeled_2` (intensity weight 0.5)
  runs through this and matches.
- `group`: the reference is the largest map (the first of equal ones); the
  others are added in input order; `out` keeps its column headers and metadata
  (the pin's `ConsensusMap::swap` does not swap `MetaInfo`) and takes
  everything else from a fresh map; `postprocess` appends identifications,
  stamps `map_index` on unassigned peptide identifications, and applies the
  three stable sorts by quality, maps and size.
- `transferSubelements`: columns renumbered consecutively in input order; the
  first of several consensus features with one unique ID wins the lookup;
  `old_map_index` renumbering and `map_index` removal as in the source.

- `QTClusterFinder`: an empty input returns before the partition and the
  maximum-m/z check; m/z partitions at gaps wider than the tolerance with the
  source's `j >= boundaries * (n / nr_partitions)` rule and a final boundary
  of `max + 1`; per partition `FeatureDistance(max_intensity, true)`, the
  `[1e-16, 1e16]` m/z and `1e16` intensity checks, and the ppm grid width
  `max_difference * max_mz * 1e-6`; clusters for every feature, one heap
  round, invalid clusters dropped at the top, updates re-filling the affected
  clusters; `createConsensusFeature_` with its adduct meta values, including
  the source's never-advanced `best_quality_index`, so the first element's
  adduct string (`Adduct::toAdductString`) overrides the best ion.
- With IDs: RT tolerances per bin by the Double-MAD rule
  `max(20, q2 + 2 * 1.4826 * (q3 - q2))` with the quantile of
  `MATH/MathFunctions.h` (`id = max(0, n*q - 1)`), bin starts at
  `minRT` and midway between medians, the diffs of the median that closes a
  bin dropped as in the source, the last bin merged with the previous one,
  `std::map::insert` semantics (an existing bin start keeps its tolerance) and
  the member never cleared between runs; unidentified neighbours outside
  their bin's tolerance skipped, others penalised by `noID_penalty`.
- `FeatureGroupingAlgorithmQT`: a fresh finder per call, `out`'s features
  replaced (`clear(false)`), then `postprocess`.

## Native differences

- Errors instead of undefined behaviour or silent insertion in
  `transferSubelements` (an unknown input map or unique ID, a missing
  `(map, old_map_index)` column that `operator[]` would add), and `None`
  instead of dereferencing `end()` in `getBestHitSequence_` for an
  identification without hits.
- Peptide identifications are validated before best-hit comparison
  (`PeptideIdentification::best_hit`); an invalid record is an error.
- `add_to_group` and `run` are transactional; the source replaces
  `pairfinder_input_[1]` before running.
- `FeatureDistance::to_param` leaves out the derived
  `distance_intensity:max_difference` that the source stores in its current
  parameters; `intensity_component()` reports it.
- Map ids are unsigned; the source's `int map_id` wraps when negative.
- `MAX_COMPARISONS` (2^36 pairs) bounds a `StablePairFinder` run.
- Progress logging and the `V_` debug macro are not ported.
- `postprocess` sorts the features directly instead of through
  `ConsensusMap::sort_by_*`, which validate the whole map and would reject two
  singletons from different inputs that share a unique ID.

- QT order: the source's Fibonacci heap resolves equal qualities by its
  internal layout; here ties go to the lower cluster id, cells are visited in
  ascending index and neighbours in map order (see `docs/CLUSTERING_SUPPORT.md`
  for `HashGrid` and `QTCluster`). Mirror-image pair clusters, the usual tie,
  hold the same elements either way. Consequently the "first element" whose
  adduct string becomes the best ion can differ from the source's when three
  or more maps carry adducts.
- QT does not sort the input features' peptide hits in place (the source
  `const_cast`s its `const` input to sort them); the grouping never reads
  them afterwards.
- `Adduct::toAdductString` (`DATASTRUCTURES/Adduct.h`, not ported) is
  reproduced privately for plain elements; an isotope-labelled adduct formula
  is `Error::Unsupported`. Non-string `Group` meta values are an error (the
  source stringifies them).
- `MAX_FEATURES` (2^24) bounds a QT run; the partition count must fit the
  source's `int`.

## Checked boundaries and evidence

- Tier 1: `tests/mapmatching_feature_grouping_topp.rs` reproduces the
  algorithm path of `TOPP_FeatureLinkerUnlabeled_1` to `_3` (test-data
  0cb15f2, `topp/CMakeLists.txt:2879-2887`; the tool path from
  `FeatureLinkerUnlabeled.cpp` at topp 174b576: `setReference`, `addToGroup`,
  `getResultMap`, `sortByMZ`) and matches the retained C++ outputs: count and
  order of consensus features, every handle's map index and unique ID exactly,
  and positions, intensities, charges and qualities within the six
  significant digits stored. The upstream whitelist (`id=`) skips handle
  lines, so this is stricter than the upstream test.
- Tier 3: every section of `FeatureDistance_test.cpp`,
  `StablePairFinder_test.cpp`, `FeatureGroupingAlgorithm_test.cpp` and
  `BaseGroupFinder_test.cpp` (`tests/mapmatching_feature_grouping.rs`).
- Tier 4: defaults and restrictions, the log-transformed intensity distance,
  input checks, `use_identifications`, `old_map_index` renumbering, and `group`
  plus the incremental path on a hand-built three-map case (upstream
  `NOT_TESTABLE`).
- Tier 1 for QT: `tests/mapmatching_qt_cluster_finder_topp.rs` reproduces
  `TOPP_FeatureLinkerUnlabeledQT_1`, `_3`, `_4` (IDs, bin tolerances) and `_6`
  (IDs with default parameters, RT 200) (test-data 0cb15f2,
  `topp/CMakeLists.txt:2907-2925`; tool path `FeatureLinkerBase.cpp`: load
  without hulls and subordinates, clear meta values except
  `dc_charge_adducts`/`Group`, `group`) and matches the retained outputs
  with the same strictness. Tier 3: every section of `QTClusterFinder_test.cpp`
  (`tests/mapmatching_qt_cluster_finder.rs`).
- Not reproduced: `TOPP_FeatureLinkerUnlabeledQT_5`, whose three clusters have
  equal quality and equal handle keys, so its order is the source heap's tie
  order (the named gap). `TOPP_FeatureLinkerUnlabeledQT_2` links consensusXML
  input through the tool and is not attempted.
- Reader defect found on the way: the featureXML and consensusXML readers
  refuse `IdentificationRun date="0000-00-00T00:00:00"`, which the source's
  `DateTime::fromString` turns into a null date. The QT differential replaces
  that attribute only when a reader refuses it, and says so at the helper.
- Not covered by tier 1: `group` of the unlabeled algorithm (no TOPP test
  calls it for featureXML input), `transferSubelements` on real data, and
  `StablePairFinder`'s `use_identifications`.

## Backlog

- `QTClusterFinder` tie order: reproducing the source for clusters of exactly
  equal quality needs `boost::heap::fibonacci_heap`'s root-list and
  consolidation order and Boost's `unordered_map` iteration order of the grid
  cells (the `HashGrid` port iterates cells in ascending order). The retained
  `FeatureLinkerUnlabeledQT_5` output is the test for it.
- `FeatureGroupingAlgorithmKD`: needs `DATASTRUCTURES/KDTree.h` (libkdtree++),
  not ported.
- `LabeledPairFinder` and `FeatureGroupingAlgorithmLabeled`, `WNetMatcher`
  and `FeatureGroupingAlgorithmWNet`, `FeatureMapping`, `PipEchoAlgorithm`:
  not started.
