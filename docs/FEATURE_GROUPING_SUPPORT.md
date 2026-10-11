# Feature grouping support (`FeatureDistance`, `BaseGroupFinder`, `StablePairFinder`, `LabeledPairFinder`, `QTClusterFinder`, `KDTreeFeatureMaps`, `FeatureMapping`, `MapAlignmentAlgorithmKD`, `WNetMatcher`, `FeatureGroupingAlgorithm*`)

Source revision: OpenMS4-core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

| Source header | Rust | Status |
|---|---|---|
| `ANALYSIS/MAPMATCHING/FeatureDistance.h` | `src/analysis/mapmatching/feature_distance.rs` | complete |
| `ANALYSIS/MAPMATCHING/BaseGroupFinder.h` | `src/analysis/mapmatching/base_group_finder.rs` | complete |
| `ANALYSIS/MAPMATCHING/StablePairFinder.h` | `src/analysis/mapmatching/stable_pair_finder.rs` | complete |
| `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h` | `src/analysis/mapmatching/feature_grouping_algorithm.rs` | complete |
| `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h` | `src/analysis/mapmatching/feature_grouping_algorithm_unlabeled.rs` | complete |
| `ANALYSIS/MAPMATCHING/QTClusterFinder.h` | `src/analysis/mapmatching/qt_cluster_finder.rs` | complete (extraction order of equal-quality clusters reproduced) |
| `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h` | `src/analysis/mapmatching/feature_grouping_algorithm_qt.rs` | complete |
| `ANALYSIS/MAPMATCHING/LabeledPairFinder.h` | `src/analysis/mapmatching/labeled_pair_finder.rs` | complete |
| `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h` | `src/analysis/mapmatching/feature_grouping_algorithm_labeled.rs` | complete |
| `ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h` | `src/analysis/mapmatching/kd_tree_feature_maps.rs` | complete |
| `ANALYSIS/QUANTITATION/KDTreeFeatureNode.h` | `src/analysis/mapmatching/kd_tree_feature_maps.rs` | complete |
| `DATASTRUCTURES/KDTree.h` (libkdtree++) | private subset in `kd_tree_feature_maps.rs` (insert, `optimise`, `find_within_range`) | partial: only what `KDTreeFeatureMaps` uses |
| `ANALYSIS/MAPMATCHING/FeatureMapping.h` | `src/analysis/mapmatching/feature_mapping.rs` | complete (reproduces the source's map-1 exclusion) |
| `ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmKD.h` | `src/analysis/mapmatching/map_alignment_algorithm_kd.rs` | partial (LOWESS without `auto_span`, `akima`) |
| `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h` (with `ClusterProxyKD`) | `src/analysis/mapmatching/feature_grouping_algorithm_kd.rs` | partial (inherits the LOWESS gap when warping) |
| `ANALYSIS/MAPMATCHING/WNetMatcher.h` | `src/analysis/mapmatching/wnet_matcher.rs` | complete |
| `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmWNet.h` | `src/analysis/mapmatching/feature_grouping_algorithm_wnet.rs` | complete |
| `wnetalign` `WNetAligner<2>`, `wnet` `WassersteinNetwork`, `pylmcf`/LEMON `NetworkSimplex` (third-party, not OpenMS headers) | `wnet_aligner.rs`, `lemon_network_simplex.rs` (crate-private solver) | the subset both OpenMS classes reach: `DIM = 2`, dense factory, no chain optimisation, no derivatives |
| `PipEchoAlgorithm.h` | none | not started |

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
| `virtual void run(const vector<ConsensusMap>&, ConsensusMap&) = 0` | `fn run(&mut self, &[ConsensusMap], &mut ConsensusMap) -> Result<()>` (mutable: `QTClusterFinder` keeps `bin_tolerances_` across runs) |
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

### LabeledPairFinder

| Source | Rust |
|---|---|
| constructor, destructor, private copy and assignment | `Default`, `from_param`, `from_param_with_warnings`, `defaults`, `to_param` |
| `run(input_maps, result_map)` | `BaseGroupFinder::run` (clock-seeded IDs); `run_with_generator(input, result, &mut UniqueIdGenerator) -> Result<RtPairEstimate>` |
| protected `PValue_(x, m, sig1, sig2)` | private `p_value` (with the source's literal `0.707106781`) |
| the `cout` estimation messages | `RtPairEstimate` (`rt_pair_dist`, `rt_dev_low`, `rt_dev_high`, `estimation_pairs`, `estimated`) and `RtPairEstimate::warnings` |

### FeatureGroupingAlgorithmLabeled

| Source | Rust |
|---|---|
| constructor, destructor, private copy and assignment | `Default`, `from_param`, `from_param_with_warnings`, `defaults`, `to_param`; the parameters are those of `LabeledPairFinder` (field `pair_finder`) |
| `group(maps, out)` | `FeatureGroupingAlgorithm::group`; `group_with_generator(maps, out, &mut UniqueIdGenerator) -> Result<RtPairEstimate>` |

### KDTreeFeatureMaps and KDTreeFeatureNode

| Source | Rust |
|---|---|
| `KDTreeFeatureMaps()`, `KDTreeFeatureMaps(maps, param)`, destructor, copy | `new()`, `from_maps(&maps)` (the parameters have no effect in the source), `Clone` |
| `typedef FeatureKDTree` | private `Tree` (libkdtree++ subset) |
| `addMaps`, `addFeature`, `optimizeTree`, `clear` | `add_maps`, `add_feature`, `optimize_tree`, `clear` (fallible for the ceilings) |
| `feature`, `rt`, `mz`, `intensity`, `charge`, `mapIndex` | same names, `Result` for an index out of range (unchecked in the source) |
| `size`, `treeSize`, `numMaps` | `size`, `tree_size`, `num_maps` |
| `getNeighborhood(index, result, rt_tol, mz_tol, mz_ppm, include_same_map, max_pairwise_log_fc)` | `get_neighborhood(...)`, appends |
| `queryRegion(rt_low, rt_high, mz_low, mz_high, result, ignored_map_index)` | `query_region(..., Option<usize>)`; `None` is `numeric_limits<Size>::max()` |
| `applyTransformations(vector<TransformationModelLowess*>)` | `apply_transformations(Fn(map, rt) -> Result<f64>)` |
| protected `updateMembers_`, `features_`, `map_index_`, `rt_`, `num_maps_`, `kd_tree_` | private fields; `updateMembers_` is empty in the source |
| `KDTreeFeatureNode(data, idx)`, copy, `getIndex`, `operator[]` | `KDTreeFeatureNode::new(&data, idx)`, `Copy`, `index`, `get(i) -> Result<f64>` (`ElementNotFound` for `i > 1`) |

### FeatureMapping

| Source | Rust |
|---|---|
| `FeatureMappingInfo` (`feature_maps`, `kd_tree`) | `FeatureMappingInfo<'a>` (`feature_maps: &'a [FeatureMap]`, `kd_tree`), built by `FeatureMappingInfo::new` |
| `FeatureToMs2Indices` (`assignedMS2` keyed by feature pointer, `unassignedMS2`) | `FeatureToMs2Indices` (`assigned_ms2` keyed by global kd-tree index, `unassigned_ms2`) |
| `static assignMS2IndexToFeature(spectra, info, mz_tol, rt_tol, ppm)` | `assign_ms2_index_to_feature(&spectra, &info, mz_tol, rt_tol, ppm) -> Result<FeatureToMs2Indices>` |

### MapAlignmentAlgorithmKD

| Source | Rust |
|---|---|
| `MapAlignmentAlgorithmKD(num_maps, param)`, destructor | `new(num_maps, &Param) -> Result<Self>`, `with_settings(num_maps, Option<KdWarpSettings>)` |
| `addRTFitData(kd_data)`, `fitLOWESS()`, `transform(kd_data)` | `add_rt_fit_data`, `fit_lowess() -> Result<Vec<String>>` (the warnings), `transform` |
| protected `updateMembers_`, `computeCCs_`, `getCCs_`, `filterCCs_` | `KdWarpSettings::from_param`, private `compute_ccs`, `get_ccs`, `filter_ccs` |
| private `fit_data_`, `transformations_`, `param_`, tolerances | private fields; `fit_data()` reads the fit data |
| `TransformationModelLowess` parameters (`LOWESS:`) | `LowessParameters` (`to_param`, `fit`) over `LowessModel` |

### FeatureGroupingAlgorithmKD and ClusterProxyKD

| Source | Rust |
|---|---|
| constructor, destructor, private copy and assignment | `Default`, `from_param`, `from_param_with_warnings`, `defaults`, `to_param` (built with the source's call sequence) |
| `group(vector<FeatureMap>)`, `group(vector<ConsensusMap>)` | `FeatureGroupingAlgorithm::group`, `group_consensus` |
| private `group_`, `runClustering_`, `updateClusterProxies_`, `computeBestClusterForCenter_`, `addConsensusFeature_` | private functions of the same names |
| private `rt_tol_secs_`, `mz_tol_`, `mz_ppm_`, `feature_distance_`, `progress_` | public parameter fields; the distance functor is built per run |
| `ClusterProxyKD` (constructors, `operator<`, `==`, `!=`, `getSize`, `isValid`, `getAvgDistance`, `getCenterIndex`) | `ClusterProxyKD` (`new`, `Default`, `Ord`, `Eq`, `size`, `is_valid`, `avg_distance`, `center_index`) |

### WNetMatcher

| Source | Rust |
|---|---|
| `enum class DistanceMetric { L1, L2, LINF }` | `DistanceMetric { L1, L2, Linf }` (defined with the aligner, re-exported) |
| `struct MatchResult { matched_pairs, cost }` | `MatchResult { matched_pairs: Vec<(usize, usize)>, cost: f64 }` |
| `static MatchResult match(positions_a, intensities_a, positions_b, intensities_b, metric = LINF, max_distance = 100, trash_cost = 100)` | `WNetMatcher::match_points(..) -> Result<MatchResult>` (`match` is a Rust keyword), with `match_with_defaults` and the constants `DEFAULT_METRIC`, `DEFAULT_MAX_DISTANCE`, `DEFAULT_TRASH_COST` |
| `static DistanceMetric metricFromString(const std::string&)` (warns and returns LINF for other strings) | `metric_from_string`; `metric_from_string_with_warning` also returns the warning text |

### FeatureGroupingAlgorithmWNet and the Wasserstein aligner

| Source | Rust |
|---|---|
| constructor (`defaults_`: `distance_metric`, `max_rt_shift`, `mz_unit`, `max_mz_shift_da`, `max_mz_shift_ppm`, `trash_cost`, `normalize_intensities`, `decharge_mz`), destructor, deleted copy and move | `Default`, `from_param`, `from_param_with_warnings`, `defaults`, `to_param` (built with the source's call sequence) |
| `group(vector<FeatureMap>, ConsensusMap&)`, `group(vector<ConsensusMap>, ConsensusMap&)` | `FeatureGroupingAlgorithm::group`, `group_consensus` |
| private `group_<MapType>`, anonymous `transformMz`, `dechargeMz`, `mapToSpectrum`, `parseDistanceMetric`, `buildAlignParams`, `convertAndNormalizeSpectra`, `collectCandidates`, `UnionFind` | private functions and types of the same names |
| `Spectrum<2>` (`sum_intensities`, `normalized`, `size`) | `Spectrum` (`new`, `sum_intensities`, `normalized`, `len`) |
| `WNetAligner<2>` (constructor, `set_point`, `total_cost`, `scale_factor`, `no_subgraphs`, `flows_for_target`, `consensus_for_target`) | `WNetAligner` (`new`, `set_point`, `total_cost`, `scale_factor`, `no_subgraphs`, `flows_for_target`, `consensus_for_target`) |
| `lemon::NetworkSimplex<StaticDigraph, int64_t, int64_t>::run(BLOCK_SEARCH)` | crate-private `lemon_network_simplex::solve` |

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

- `LabeledPairFinder`: exactly one input map and two column headers with one
  file name, labeled `light` and `heavy`; RT estimation from all
  charge- and m/z-compatible pairs: median, a 100-bin histogram of the
  `|map|/2` distances around it, the median bin count as background, the
  highest bin as centre, the sigma walk down to the background, and
  `GaussFitter` (bit-identical to the Release build on both oracle inputs);
  every candidate pair gets a fresh unique ID, the light features' meta
  values are cleared, both features' identifications are appended with their
  `map_index`; greedy selection by quality (stable) on the handles' unique
  IDs; protein and unassigned identifications appended; stable sort by m/z.
- `KDTreeFeatureMaps`: the libkdtree++ tree is reproduced node for node:
  insertion left when strictly less on the level's axis, `optimise()` from the
  in-order sequence with libstdc++'s introselect `nth_element`, and range
  queries in pre-order with the region-bounds pruning; coordinates are read
  live, so queries after `apply_transformations` use the new RTs on the old
  tree. The query order is checked against the Release build.
- `FeatureMapping`: the source passes `true` as `queryRegion`'s
  `ignored_map_index`, so features of map 1 are never assigned; reproduced.
  Ties in m/z distance go to the first feature in query order; the starting
  distance is 1e11 and the starting index 0.
- `MapAlignmentAlgorithmKD`: BFS connected components over
  `getNeighborhood(..., false, max_pairwise_log_fc)`; the size filter
  `max(2, min_rel_cc_size * num_maps)` truncated to an integer; the conflict
  filter; the source's charge filter is a no-op (its `continue` only continues
  the inner loop), reproduced; fewer than 50 fit points give the identity
  through (0,0), (1,1), (1e6,1e6).
- `FeatureGroupingAlgorithmKD`: partition boundaries at gaps above the larger
  m/z tolerance with the `j >= boundaries * (n / nr_partitions)` rule and a
  final `max + 1`; the aligner collects fit data from every partition before
  one LOWESS fit per map; linking per partition on a fresh kd-tree (warped
  and re-optimised when enabled); the distance functor from the link
  tolerances, `ignore_charge` only for `Any`; clusters chosen by size, then
  average distance, then the higher centre index; average quality summed in
  `f32`; the `best ion` index starting at global index 0 (a source quirk,
  see the module documentation); `postprocess` as for the other algorithms.

- Wasserstein aligner (`wnet_aligner.rs`): intensities and positions are scaled
  by `s = sqrt(2^60 / (max(sum A, sum B) * trash_cost))` and truncated to
  `int64`; distances are computed in `double` from the scaled positions and
  truncated; an arc exists when the truncated distance is at most
  `trunc(max_distance * s)` and the first coordinates lie within the window
  `x +- (double) trunc(max_distance * s)`. The flow network of a connected
  component is source -> empirical peak (capacity intensity), empirical ->
  theoretical (capacity `min(theoretical * point, empirical)`, cost distance),
  theoretical -> sink, and the trash arc source -> sink (capacity the total,
  cost `trunc(trash_cost * s)`), with the total the larger of the two
  intensity sums *of that component*. Components are numbered by their lowest
  peak, peaks without an arc belong to none, and arcs sit in the order of
  `(start, end)`, which is the `StaticDigraph` arc numbering.
- LEMON's `StaticDigraph` visits nodes and arcs from the highest id down, so
  the simplex's internal node `i` is graph node `N - 1 - i` and its arcs are
  stored in the `arc_mixing` order (`skip = max(M / N, 3)`). Block search with
  block size `max(floor(sqrt(M)), 10)`, the heuristic initial pivots (reverse
  search from the sink over in-arcs by descending arc id), and the
  strongly-feasible-tree leaving-arc rule (`d < delta` on the first path,
  `d <= delta` on the second) pick *which* optimum the solver stops at when
  costs tie; the port transcribes all of it.
- The consensus is `std::sort` (the C++ Release build's introsort, through
  `math::source_sort`) of the matching arcs by descending flow, then a greedy
  pass taking an arc while neither peak is used. Equal flows therefore end in
  the libstdc++ order, which the lattice instances of the oracle exercise.
- `FeatureGroupingAlgorithmWNet`: `trash_cost <= 0` means `max_rt_shift`; the
  m/z axis is `ln(m/z)` for ppm with `mz_scale = max_rt_shift / log1p(ppm *
  1e-6)` and `m/z` for Da with `max_rt_shift / max_mz_shift_da`; the maximum
  distance is `max_rt_shift`, `2 * max_rt_shift` (L1) or `sqrt(2) *
  max_rt_shift` (L2); normalisation divides by the intensity sum only when a
  map is non-empty; map pairs with an empty map are skipped; a candidate's
  cost is the pair's total cost divided by the number of consensus pairs;
  candidates are sorted by `(cost, lower, higher flat feature number)` and
  merged by a union-find that refuses merges putting two features of one map in
  a group; groups become consensus features in ascending order of their root,
  members by map and feature, then `postprocess`. Consensus features are
  appended to `out`'s (the source does not clear it).

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

- QT order is reproduced, not approximated: the heap is
  `math::fibonacci_heap::FibonacciHeap` (push, lazy update, pop with
  consolidation, as `boost/heap/fibonacci_heap.hpp` 1.89), the grid
  (`HashGrid`) iterates cells and pairs in Boost.Unordered 1.89's order, and a
  cluster's per-map neighbours and a feature's cluster set are
  `math::hash_order::StdHashMap`s in libstdc++ 14's `std::unordered_map`
  order (identity hash, prime rehash policy, bucket-begin insertion). Which of
  several equal-quality clusters is extracted first, which of several equally
  near features becomes a neighbour, the order a cluster lists its elements
  (and so the "first element" whose adduct string becomes the best ion), and
  the order the neighbour distances are summed in (floating point addition
  is not associative) therefore follow the executed C++. The orders are
  validated by `tests/hash_order_oracle.rs` (14 operation scripts against
  `std::unordered_map<size_t, int>` and `OpenMS::HashGrid<int>`) and end to
  end by `tests/mapmatching_qt_tie_order_oracle.rs`. Not ported from
  `boost::heap::fibonacci_heap`, because nothing reaches them: `increase`,
  `decrease`, eager `update`, `erase`, `merge` and the cascading-cut marks.
  The pointer-keyed `std::unordered_map<const GridFeature*, ...>` of the
  source (`ElementMapping`, `tmp_element_mapping`) is only looked up here:
  its iteration order depends on object addresses, the source's loops over
  it touch each feature's own set independently, so the result does not
  depend on it.
- QT does not sort the input features' peptide hits in place (the source
  `const_cast`s its `const` input to sort them); the grouping never reads
  them afterwards.
- `Adduct::toAdductString` (`DATASTRUCTURES/Adduct.h`, not ported) is
  reproduced privately for plain elements; an isotope-labelled adduct formula
  is `Error::Unsupported`. Non-string `Group` meta values are an error (the
  source stringifies them).
- `MAX_FEATURES` (2^24) bounds a QT run; the partition count must fit the
  source's `int`.

- `LabeledPairFinder`: the input is sorted by position stably (the source's
  `std::sort` leaves equal positions unordered); a candidate score that is
  not a number (a zero deviation meeting an exact match) is an error; the
  estimation's histogram walk is bounded; `MAX_COMPARISONS` and
  `MAX_CANDIDATES` bound a run; the RT estimate is returned instead of
  printed.
- `KDTreeFeatureMaps`: `MAX_FEATURES` and `MAX_DEPTH` (a degenerate
  insertion chain) are checked; accessors report an out-of-range index; the
  partition scan of `nth_element` is bounded, which only matters for
  coordinates that are not numbers (undefined behaviour in the source).
- `FeatureMapping`: `assigned_ms2` is keyed by the feature's global kd-tree
  index; the source keys by address, whose order follows memory layout.
- `FeatureGroupingAlgorithmKD`: a failed LOWESS fit is returned as an error
  (the source logs it and returns an empty, unprocessed result); a candidate
  without a finite distance is an error (the source records `SIZE_MAX`);
  input features must have finite coordinates.
- `MapAlignmentAlgorithmKD`: `TransformationModelLowess` is not ported; the
  fit runs through `LowessModel` without `auto_span` or `akima`
  (`Error::Unsupported`); an aligner built from an empty tree is
  unconfigured and its operations fail (the source leaves the tolerances
  uninitialised).
- QT, labeled and KD consensus features take their elements' peptide
  identifications (with `map_index`), as `ConsensusFeature::insert(map,
  BaseFeature)` does in the source (`insert_base_feature`). The QT port
  previously inserted only the handles; fixed, and the TOPP differentials now
  compare identifications.

- WNet: the source converts floating-point values to `int64_t` without a range
  check, which is undefined behaviour for NaN, infinite and out-of-range
  values. A non-finite scaled position, an intensity, capacity or distance that
  does not fit 64 bits, and a scale factor that is not finite (`trash_cost` 0
  with `max_rt_shift` 0) are `Error::InvalidValue` here. An m/z tolerance of
  exactly 0 (infinite `mz_scale`, every scaled distance `NaN`) is therefore an
  error rather than garbage costs.
- WNet: `std::sort` ties are reproduced (see above); `unordered_set` use in
  `consensus_for_target` has no order effect and is a `HashSet`.
- WNet ceilings the source lacks: `wnet_aligner::MAX_PEAKS` (2^25) and
  `MAX_MATCHING_ARCS` (2^24) per aligner, `feature_grouping_algorithm_wnet::
  MAX_FEATURES` (2^24) over all maps. The source logs per-pair counts at INFO;
  no logging is done here. A non-numeric `dc_charge_adduct_mass` meta value is
  an error (the source reads a union member of the wrong type).
- WNet: the caller-supplied `scale_factor` constructor argument, the
  one-dimensional chain factory (`create_1d`) and the derivative members
  (`signal_part_derivatives`, `spectrum_proportion_derivatives`) of the
  `wnet` library are not ported because no OpenMS code reaches them.

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
- `TOPP_FeatureLinkerUnlabeledQT_5` (three clusters of equal quality, RT 200)
  is reproduced exactly against the *executed* Release build
  (`tests/mapmatching_qt_tie_order_oracle.rs`, case `topp_qt5`, both views),
  but its retained output is not what that build produces: in the middle
  cluster the build takes the map-1 feature at RT 310.1 where the retained file
  has the one at RT 310.0, two equally near candidates. The retained file was
  not generated by the pinned Linux build (which of its differences causes the
  tie to go the other way is not established); the upstream test therefore
  depends on the platform. `topp_feature_linker_unlabeled_qt_5_agrees_with_
  retained_except_one_tie` states the difference exactly. (An earlier
  incarnation of this lane attributed the QT_5 mismatch to the tie order of
  the heap; the heap order is now exact, and the remaining difference is this
  candidate choice.)
- Tier 1 against the Release build for QT (oracle `M2-qt`, ibminode06,
  `openms4-release-bc9cc12-c19e494-174b576`, two identical runs):
  `tests/mapmatching_qt_tie_order_oracle.rs` compares 17 generated instances
  (noisy copies of a base map under Da and ppm tolerances with up to 20
  partitions, dense competition, integer lattices with equal intensities,
  lattices with sub-tolerance jitter, exact duplicates, mirror pairs,
  equidistant candidates) and the five `TOPP_FeatureLinkerUnlabeledQT_*`
  paths, each in two views: the consensus features in the order
  `QTClusterFinder::run` extracts them (handles, quality exactly) and after
  `FeatureGroupingAlgorithmQT::group`'s canonical sorts. With ties broken by
  cluster id, ascending cells and ascending map order, 7 of the 22 grouped
  comparisons (dense3, dense4, dups3, lat3, lat4, lat5 and QT_5) and 21 of
  the 22 extraction-order comparisons differed from the executed C++; all 44
  agree now.
- Reader defect found on the way: the featureXML and consensusXML readers
  refuse `IdentificationRun date="0000-00-00T00:00:00"`, which the source's
  `DateTime::fromString` turns into a null date. The QT differential replaces
  that attribute only when a reader refuses it, and says so at the helper.
- Not covered by tier 1: `group` of the unlabeled algorithm (no TOPP test
  calls it for featureXML input), `transferSubelements` on real data, and
  `StablePairFinder`'s `use_identifications`.

- Tier 1 against the Release build (oracle `M2-labeled`, ibminode06,
  `openms4-release-bc9cc12-c19e494-174b576`, two identical runs):
  `tests/mapmatching_labeled_grouping_oracle.rs` (FeatureLinkerLabeled_1/_2
  path and the class test's estimation input, fitted RT distance and
  deviation to 1e-12), `tests/mapmatching_kd_tree_feature_maps_oracle.rs`
  (1589 features in two maps: neighbourhood and region queries in the exact
  order returned; FeatureMapping's map-1 exclusion),
  `tests/mapmatching_kd_grouping_oracle.rs` (KD defaults tree;
  FeatureLinkerUnlabeledKD_1 to _8 paths including consensus input, LOWESS
  warp, charge and adduct merging, `best ion` and `LinkedGroups`).
- Tier 1 against retained outputs: `tests/mapmatching_labeled_grouping_topp.rs`
  (FeatureLinkerLabeled_1, _2; case 2's retained qualities predate the pin's
  `GaussFitter` and differ from the Release build by up to 3.7e-4 relative,
  so they are compared with the oracle instead) and
  `tests/mapmatching_kd_grouping_topp.rs` (FeatureLinkerUnlabeledKD_1 to _8,
  case 2 through `keep_subelements` and `transferSubelements`).
- Tier 3: every section of `LabeledPairFinder_test.cpp`,
  `FeatureGroupingAlgorithmLabeled_test.cpp`, `KDTreeFeatureMaps_test.cpp`,
  `KDTreeFeatureNode_test.cpp`, `FeatureMapping_test.cpp`,
  `FeatureGroupingAlgorithmKD_test.cpp` and `MapAlignmentAlgorithmKD_test.cpp`.

- Tier 1 against executed C++ for WNet (oracle `M2-wnet`, ibminode06, two
  identical runs): the Release install `openms4-release-bc9cc12-c19e494-174b576`
  was configured without WNet (no `WNetMatcher` or
  `FeatureGroupingAlgorithmWNet` symbol in `libOpenMS.so`), so the oracle
  compiles the pinned `WNetMatcher.cpp` and `FeatureGroupingAlgorithmWNet.cpp`
  with that install's g++ 14.4.0, headers and `libOpenMS` flags against the
  header-only libraries at the commits `cmake_findExternalLibs.cmake` pins
  (wnetalign v0.9.8 `cff6a19`, wnet v0.9.11 `18a1525`, pylmcf v0.9.8
  `d2c9c52`). `tests/mapmatching_wnet_oracle.rs`: 31 generated instances
  (random under L1/L2/LINF, equal-weight integer lattices and a full grid with
  massive cost and flow ties, duplicates, zero intensities, 10-component
  cluster sets, trash from 0.5 to 1000, zero and huge maximum distance,
  realistic scaled `ln(m/z)` positions with normalised intensities up to 400
  points per side); scale factor and total cost compared bit for bit, every
  flow triple in order, the consensus in selection order, and
  `WNetMatcher::match`. `tests/mapmatching_wnet_grouping_oracle.rs`: 13 tool-path
  cases (default, Da, L2, L1 without normalisation, trash 500, tight
  tolerances, the unlabeled 2 and 3 inputs, charge and adduct inputs with and
  without `decharge_mz`, consensus-map input) and the four class-test
  sections; the defaults tree equals the executed C++'s stored `ParamXML`.
- Tier 1 against the retained output: `tests/mapmatching_wnet_topp.rs`
  reproduces `TOPP_FeatureLinkerWNet_1` (test-data 0cb15f2,
  `topp/CMakeLists.txt:2955-2958`: the FeatureLinkerUnlabeled_1 inputs and
  the default parameters) with every handle exact; the retained file was
  generated by tool version `3.6.0-pre-wnetalign-2026-03-25` and agrees with
  the pin here.
- Tier 3: every section of `WNetMatcher_test.cpp` and
  `FeatureGroupingAlgorithmWNet_test.cpp` (`tests/mapmatching_wnet.rs`).
  Tier 4: the single-pair cost, the competition for one theoretical peak,
  unconnected pairs, and every refusal listed under native differences.

## Backlog

- `TransformationModelLowess` (unmapped): automatic span selection and
  `akima` interpolation, which `MapAlignmentAlgorithmKD` and
  `FeatureGroupingAlgorithmKD` need for those parameter values.
- `PipEchoAlgorithm` (and the private `PIPECHO/` sources, about 3100 lines):
  untouched. Its `TransferFDRModel` trains a linear-kernel SVM through
  `ML/SVM/SimpleSVM.h` (libsvm: cross-validated grid search over `log2_C`,
  probability estimates), which is not ported (`docs/core-sdk-coverage.json`
  lists `SimpleSVM.h` as unmapped; only fixed-model RBF prediction exists, in
  `feature_finding_metabo/predictor.rs`). The Release install does contain
  `PipEchoAlgorithm`, so an oracle is available once the SVM trainer exists.
