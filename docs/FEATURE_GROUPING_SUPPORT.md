# Feature grouping support (`FeatureDistance`, `BaseGroupFinder`, `StablePairFinder`, `FeatureGroupingAlgorithm*`)

Source revision: OpenMS4-core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

| Source header | Rust | Status |
|---|---|---|
| `ANALYSIS/MAPMATCHING/FeatureDistance.h` | `src/analysis/mapmatching/feature_distance.rs` | complete |
| `ANALYSIS/MAPMATCHING/BaseGroupFinder.h` | `src/analysis/mapmatching/base_group_finder.rs` | complete |
| `ANALYSIS/MAPMATCHING/StablePairFinder.h` | `src/analysis/mapmatching/stable_pair_finder.rs` | complete |
| `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h` | `src/analysis/mapmatching/feature_grouping_algorithm.rs` | complete |
| `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h` | `src/analysis/mapmatching/feature_grouping_algorithm_unlabeled.rs` | complete |
| `QTClusterFinder.h`, `FeatureGroupingAlgorithmQT.h` | none | not ported (see Backlog) |
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
- Not covered by tier 1: `group` itself (no TOPP test calls it for featureXML
  input), `transferSubelements` on real data, and `use_identifications`.

## Backlog

- `QTClusterFinder` and `FeatureGroupingAlgorithmQT`: depend on the `HashGrid`
  iteration order (ported `partial`), the tie order of
  `boost::heap::fibonacci_heap` with `update_lazy`, and `unordered_set`
  iteration of cluster ids; the retained `FeatureLinkerUnlabeledQT_*` outputs
  are the tier-1 target.
- `FeatureGroupingAlgorithmKD`: needs `DATASTRUCTURES/KDTree.h` (libkdtree++),
  not ported.
- `LabeledPairFinder` and `FeatureGroupingAlgorithmLabeled`, `WNetMatcher`
  and `FeatureGroupingAlgorithmWNet`, `FeatureMapping`, `PipEchoAlgorithm`:
  not started.
