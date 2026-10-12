# Map alignment

`analysis::mapmatching` holds the map-alignment headers of `ANALYSIS/MAPMATCHING`: the superimposers, the alignment algorithms and the alignment evaluation. The transformation models they produce are in `analysis::transformations` ([TRANSFORMATIONS_SUPPORT.md](TRANSFORMATIONS_SUPPORT.md)), and applying a transformation to a map is `analysis::alignment_transformer`. The reference is OpenMS4-core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`; hashes are in [map_alignment_provenance.json](../tests/data/map_alignment_provenance.json). No dependency is added.

| Header | Rust | Status |
|---|---|---|
| `BaseSuperimposer.h` | `mapmatching/superimposer.rs` | complete |
| `PoseClusteringShiftSuperimposer.h` | `mapmatching/pose_clustering_shift.rs` | complete |
| `PoseClusteringAffineSuperimposer.h` | `mapmatching/pose_clustering_affine.rs` | complete |
| `MapAlignmentEvaluationAlgorithm.h` | `mapmatching/alignment_evaluation.rs` | complete |
| `MapAlignmentEvaluationAlgorithmPrecision.h` | `mapmatching/alignment_evaluation.rs` | complete |
| `MapAlignmentEvaluationAlgorithmRecall.h` | `mapmatching/alignment_evaluation.rs` | complete |
| `MapAlignmentAlgorithmIdentification.h` | `mapmatching/alignment_identification.rs` | complete |
| `MapAlignmentAlgorithmPoseClustering.h` | `mapmatching/alignment_pose_clustering.rs` | complete (composes `mapmatching/stable_pair_finder.rs`) |
| `MapAlignmentAlgorithmKD.h` | `mapmatching/alignment_kd.rs` | complete; runs on the trait `KdFeatureMaps`, which `KDTreeFeatureMaps` (`mapmatching/kd_tree_feature_maps.rs`) implements |
| `MapAlignmentAlgorithmTreeGuided.h` | `mapmatching/alignment_tree_guided.rs` | complete (`buildTree` on `ml/clustering`) |

## API mapping

### `BaseSuperimposer`

| Source | Rust |
|---|---|
| abstract class (`DefaultParamHandler` named `BaseSuperimposer`, `ProgressLogger`) | trait `BaseSuperimposer` |
| pure virtual `run(map_model, map_scene, TransformationDescription&)` | `superimpose(model, scene, &mut ProgressReporter) -> Result<Superimposition>`; provided `run(model, scene) -> Result<TransformationDescription>` |
| private copy constructor and assignment | none: implementors are plain values |

`Superimposition` carries the transformation and the lines the source writes to `OPENMS_LOG_WARN`.

### `PoseClusteringShiftSuperimposer`

| Source | Rust |
|---|---|
| constructor with its six defaults | `Default`; `defaults()`, `to_param()`, `from_param()`, `from_param_with_warnings()` |
| `run(const ConsensusMap&, const ConsensusMap&, TransformationDescription&)` | `BaseSuperimposer::superimpose` / `run` |
| parameters `mz_pair_max_distance`, `num_used_points`, `shift_bucket_size`, `max_shift`, `dump_buckets`, `dump_pairs` | fields of the same names; `max_buckets` is a native limit |

### `PoseClusteringAffineSuperimposer`

| Source | Rust |
|---|---|
| constructor with its nine defaults | `Default`; `defaults()`, `to_param()`, `from_param()`, `from_param_with_warnings()` |
| `run(const ConsensusMap&, const ConsensusMap&, TransformationDescription&)` | `BaseSuperimposer::superimpose` / `run` |
| `run(const std::vector<Peak2D>&, const std::vector<Peak2D>&, TransformationDescription&)` | `superimpose_points(model, scene, &mut ProgressReporter)` |
| parameters `mz_pair_max_distance`, `rt_pair_distance_fraction`, `num_used_points`, `scaling_bucket_size`, `shift_bucket_size`, `max_shift`, `max_scaling`, `dump_buckets`, `dump_pairs` | fields of the same names; `max_buckets`, `max_quadruples` are native limits |
| file-local `initializeHashTables`, `affineTransformationHashing`, `scalingEstimate`, `shiftEstimate`, `computeIntensityRatio` | private: `Hashes`, `Hashing::hash`, and the shared `Histogram` methods |

### `MapAlignmentEvaluationAlgorithm` and subclasses

| Source | Rust |
|---|---|
| abstract `MapAlignmentEvaluationAlgorithm` | trait `MapAlignmentEvaluationAlgorithm` |
| pure virtual `evaluate(in, gt, rt_dev, mz_dev, int_dev, use_charge, double& out)` | `evaluate(...) -> Result<f64>` |
| `isSameHandle` | provided trait method `is_same_handle`; free function `is_same_handle` |
| protected `HandleIterator` typedef | none: handles are a slice |
| `MapAlignmentEvaluationAlgorithmPrecision` | `MapAlignmentEvaluationPrecision` |
| `MapAlignmentEvaluationAlgorithmRecall` | `MapAlignmentEvaluationRecall` |

### `MapAlignmentAlgorithmIdentification`

| Source | Rust |
|---|---|
| constructor with its eight defaults | `Default`; `defaults()`, `to_param()`, `set_parameters()`, `from_param()` |
| `setReference<DataType>(data)` | `set_reference<D: AlignmentData>(&data)`, returning its info lines |
| `align<DataType>(data, transformations, reference_index = -1)` | `align(&data, Option<usize>)` and `align_with_progress`, returning `IdentificationAlignment` (transformations, info lines, warnings) |
| the `IsFCMap` concept and the four `getRetentionTimes_` overloads (`PeptideIdentificationList`, `IdentificationData`, feature and consensus maps) | trait `AlignmentData`, implemented for `Vec<PeptideIdentification>`, `IdentificationData`, `FeatureMap`, `ConsensusMap` |
| typedefs `SeqToList`, `SeqToValue` | `SeqToList`, `SeqToValue` (`BTreeMap`, ordered as `std::map<std::string, ...>`) |
| protected `reference_index_`, `reference_`, `min_run_occur_`, `use_feature_rt_`, `use_adducts_`, `min_score_`, `score_cutoff_`, `score_type_` | private state; `reference()` reads the reference |
| protected `better_` | private state of type `ScoreFilter`; `score_filter()` reads it |
| `computeMedians_`, `computeTransformations_`, `checkParameters_` | private `medians`, `compute_transformations`, `check_parameters` |
| `getReference_` (declared, never defined in the pinned `.cpp`) | none |
| `handleIdDataScoreType_` | private `identification_score_type` |
| `getBestScoringHit` | private `best_hit` |

### `MapAlignmentAlgorithmPoseClustering`

| Source | Rust |
|---|---|
| constructor: `superimposer:` and `pairfinder:` sections, `max_num_peaks_considered` | `Default`; `defaults()`, `to_param()`, `set_parameters()`; the `pairfinder:` section configures the composed `StablePairFinder` |
| `setReference<MapType>(map)` | `set_reference_features`, `set_reference_peaks`, `set_reference_consensus` |
| `align(const FeatureMap&, ...)`, `align(const PeakMap&, ...)`, `align(const ConsensusMap&, ...)` | `align_features`, `align_peaks`, `align_consensus` (with the composed `StablePairFinder`), and `*_with_pair_finder` variants that take the pairing step as a closure |
| `updateMembers_` | inside `set_parameters` |
| protected `superimposer_`, `pairfinder_`, `reference_`, `max_num_peaks_considered_` | `superimposer`, `pair_finder` (a `StablePairFinder`), `reference()`, `max_num_peaks_considered` |

### `MapAlignmentAlgorithmKD`

| Source | Rust |
|---|---|
| constructor `(Size num_maps, const Param& param)` | `new(num_maps, &Param) -> Result<Self>` |
| `addRTFitData(const KDTreeFeatureMaps&)` | `add_rt_fit_data(&K)` for `K: KdFeatureMaps` |
| `fitLOWESS()` | `fit_lowess() -> Result<Vec<String>>`, returning the source's warnings |
| `transform(KDTreeFeatureMaps&) const` | `transform(&mut K)` |
| protected `computeCCs_`, `getCCs_`, `filterCCs_` | public `compute_ccs` (labels and count), `get_ccs`, `filter_ccs` |
| protected `updateMembers_` | inside `new` |
| private `fit_data_`, `transformations_` | read-only `fit_data()`, `transformations()` (`LowessModel`s) |
| private `param_`, `max_pairwise_log_fc_`, `rt_tol_secs_`, `mz_tol_`, `mz_ppm_` | private state |
| private default constructor | none |
| the `KDTreeFeatureMaps` members it calls: `size`, `rt`, `charge`, `mapIndex`, `getNeighborhood`, `applyTransformations` with `optimizeTree` | trait `KdFeatureMaps`: `len`, `rt`, `charge`, `map_index`, `neighborhood`, `apply_transformations` |

### `MapAlignmentAlgorithmTreeGuided`

| Source | Rust |
|---|---|
| constructor: `model:` (`TransformationModelDefaults::getDefaults("b_spline")`), `model_type`, `align_algorithm:` with `use_feature_rt` `true` | `new()`; `defaults()`, `to_param()`, `set_parameters()`, `from_param()`, `from_param_with_warnings()` |
| static `buildTree(feature_maps, tree, maps_ranges)` | `build_tree(&[FeatureMap]) -> Result<(Vec<BinaryTreeNode>, Vec<Vec<f64>>)>`: `ClusterHierarchical` with `AverageLinkage` on `peptide_identifications_pearson_similarity` |
| `treeGuidedAlignment(tree, feature_maps_transformed, maps_ranges, map_transformed, trafo_order)` | `tree_guided_alignment(&[TreeMerge], &mut [FeatureMap], &[Vec<f64>], &mut UniqueIdGenerator) -> Result<TreeGuidedMerge>` |
| `align(data, transformations)` | `align(&mut [FeatureMap], &mut UniqueIdGenerator) -> Result<TreeGuidedAlignment>` (builds the tree, logs `  Alignment follows Newick tree: ...` with `ClusterAnalyzer::newick_tree` and the distances, then `align_with_tree`); `align_with_tree(&mut [FeatureMap], &[TreeMerge], &mut UniqueIdGenerator)` takes the tree from the caller |
| `computeTrafosByOriginalRT(feature_maps, map_transformed, transformations, trafo_order)` | `compute_trafos_by_original_rt(&[FeatureMap], &FeatureMap)`; the unused `trafo_order` is dropped |
| static `computeTransformedFeatureMaps(feature_maps, transformations)` | `compute_transformed_feature_maps(&mut [FeatureMap], &[TransformationDescription])` |
| protected typedef `SeqAndRTList` | `SeqAndRtList` (`BTreeMap<String, Vec<f64>>`) |
| protected `updateMembers_`, `model_type_`, `model_param_`, `align_algorithm_` | inside `set_parameters`; `model_type`, `model_param()`, `align_algorithm` |
| protected nested functor `PeptideIdentificationsPearsonDistance_` | `peptide_identifications_pearson_similarity` |
| protected static `addPeptideSequences_`, `extractSeqAndRt_` | `extract_seq_and_rt` |
| `BinaryTreeNode` (`left_child`, `right_child`, `distance`) | `TreeMerge = (usize, usize)`; the alignment does not use the distance |
| private copy constructor and assignment | none: a plain value |

## Preserved source conventions

**Histograms.** Both superimposers vote into `Math::LinearInterpolation` histograms: a vote at a fractional bucket is split linearly between the two neighbouring buckets, a vote within one bucket outside either end lands on that end, and a vote further out is dropped. The histogram is filtered with a 21-sample top-hat (`MorphologicalFilter::filterRange`, instantiated on `double`), its noise level taken from the sorted values against the line through the first bucket with slope `(last - first) / n / 3`, and the buckets below it zeroed; three rounds of a weighted mean and `1.5 * sd` cut-off then give the centroid. The affine superimposer's high end uses a range one bucket wider on each side, as in the source.

**Selection of points.** The shift superimposer keeps the `num_used_points` most intense elements with `ConstRefVector::sortByIntensity(true)` and sorts them by m/z, the affine one with `std::nth_element` on the reversed vector; both orders are the Release build's libstdc++ permutations (`math::source_sort` for `std::sort`, a generic port of `std::nth_element` for the other). `-1` keeps all points. The affine superimposer takes its RT ends from all points, not only the kept ones.

**Hashing windows.** Pairs are formed only within `mz_pair_max_distance`; each m/z window gets unit weight shared by its elements, less 0.1, so windows of ten or more elements are skipped. The last model element is never a first element of a pair (`i < size - 1`), and the affine superimposer's second-point window is centred on the first point's m/z, both as in the source.

**Results.** Each superimposer returns a `linear` transformation fitted from explicit `slope` and `intercept` without data, so its `model_parameters()` holds exactly those two, as the class tests assert.

**Evaluation.** Handles are the same feature when RT, m/z and intensity are within the deviations (the intensity difference in `float`) and, with `use_charge`, the charges agree. Only ground-truth consensus features of at least two elements count. Precision averages `shared / tool size` over them, recall `shared / (matched tool features * ground-truth size)`.

**Identification alignment.** Only the best hit of each identification is used. Retention times are grouped per sequence (for identification data: per molecule, with `+[adduct]` under `use_adducts`), the median taken per input, and the reference is the given one or the median of medians over the sequences occurring in at least `min_run_occur` runs. Data points are `(input median, reference)`, annotated with the sequence, within `max_rt_shift`. For identification data without a named score type, the score type of the first scored match's most recent score is used (`IdentificationData::pickScoreType`, by the processing software's score priority) and kept. The info lines (`Using score type: <name>`, "Alignment based on:", "- N data points for sample I", "(reference)", "(K outliers removed)") and the warnings are returned as the source logs them; `OPENMS_LOG_DEBUG` lines are not. The object's state is the source's:

- `set_reference` reads `use_feature_rt`, `score_cutoff` and `score_type` from the parameters but not `min_score`, which only an alignment copies, so a reference collected before the first alignment compares scores with 0;
- a reference given by `set_reference` is filtered by each alignment and stays filtered; an internal reference is cleared at the start of the next internal-reference alignment;
- the score filter is set only when a feature or consensus map is collected, and identification lists use whatever filter is current, so `score_cutoff` has no effect on identification lists aligned before any map;
- with `use_feature_rt`, a feature whose identifications all fail the score filter is still recorded, under the empty sequence.

**Pose-clustering alignment.** The superimposer's transformation is applied to every consensus feature and to its first handle only; the pairs' scene retention times are taken back through the algebraic inverse; the final model is `linear` with default parameters, which fails without pairs as the source's linear constructor does.

**KD alignment.** Components are found by breadth-first search seeded at the lowest unvisited feature, so the labels, the component order and the member order do not depend on the order in which the kd-tree returns a neighbourhood; the fit data are `(feature RT, component average RT)` in component order, the average summed in member order. A component is kept when it has at least `max(2, warp:min_rel_cc_size * number of maps)` features (`std::max`, then truncated to `Size`) and, unless `warp:max_nr_conflicts` is -1, at most that many repeated maps. The class documentation also demands one charge per component, but the source's `continue` on a second charge continues only the loop over the members, so charges never remove a component; this is kept. A map with fewer than 50 data points gets a LOWESS model of `(0, 0), (1, 1), (1e6, 1e6)` and the source's warning. `warp:min_rel_cc_size` and `warp:max_nr_conflicts` are read at each filtering, the tolerances at construction.

**Tree-guided alignment.** Per map, every identification of a feature with at least one hit adds the feature RT under the sequence of its *first* hit (not the best) and to the map's RT list, which is sorted. The similarity of two maps walks both sequence lists together; the "median" RT of a shared sequence is `Math::median(..., sorted = true)` on the *unsorted* per-sequence list, the middle element by position; the Pearson coefficient of the two median lists is cast to `float` and multiplied by the shared count over the number of walk steps, which the source calls the union size but which stops at the end of the shorter list. `ClusterHierarchical` turns it into the distance `1 - similarity`. Along the tree, the child with the larger 10/90 percentile range (`v[size * 0.9] - v[size * 0.1]`, indices truncated; the right child on a tie) is the reference; the ranges stay those of the original maps, as the source never merges them. The merged map lives at the smaller index. Merging is `FeatureMap::operator+=`: document identifier, file path and type and the map's unique ID reset; proteins, unassigned identifications, processing records and features appended; colliding unique IDs replaced. The final transformations align each original map to the merged map (reference index 1) and fit `model_type` with its `model:` subsection.

## Native differences

- **Debug dumps.** A non-empty `dump_buckets` or `dump_pairs` writes the source's text files with the source's names and lines (numbers at six significant digits): `<base>_<serial>` and `<base><serial>` for the shift superimposer; `<base>_scale_<serial>`, `_low_`, `_high_` and `<base>_phase_two_<serial>` for the affine one, whose first hashing round truncates the pairs file the second then fills, as in the source. The serial is the source's function-local `static` invocation counter, one per class, kept as a process-wide atomic; it is the one piece of global state in these modules. A file that cannot be created or written is an `Error::Io`, where the source's `std::ofstream` fails silently and the run continues.
- **Undefined behaviour is refused.** An empty map (the source loops past the end of the model map in the shift superimposer), a zero bucket size (an infinite bucket count converted to `Int`), a histogram whose noise level the source can only find by reading past its buffer, a pair similarity or intensity ratio that is not finite (two zero intensities give `0/0`), a first feature without identification under `score_cutoff` (the source reads `features[0].getPeptideIdentifications()[0]` unchecked), and a consensus feature without handles in pose-clustering alignment (the source writes through an end iterator) are errors.
- **Identifications without retention time are refused** in the identification alignment, where the source records NaN, which then decides medians and sorts arbitrarily. With `use_feature_rt` a missing RT never becomes the closest identification, as NaN never is in the source.
- **Ranges are computed, not trusted.** The source states `updateRanges()` as a precondition of the superimposers "for performance reasons"; the ranges are computed from the maps here (with handles, as `ConsensusMap::updateRanges` does).
- **A failed call ends its progress section** (`ProgressReporter::section`), where the source's exception skips `endProgress`.
- **Precision of an unusable ground truth** is the source's NaN (`1.0 / 0 * 0`); recall refuses it, as the source does since its guard was added.
- **Resource limits.** Histogram sizes and the number of quadruples hashed by the affine superimposer are bounded (`max_buckets`, `max_quadruples`); the source bounds neither.
- **The pair finder is composed, as in the source.** `MapAlignmentAlgorithmPoseClustering` holds a `StablePairFinder` (`pairfinder:` section, defaults included); the `*_with_pair_finder` variants of the three `align_*` methods take the pairing step as a closure for tests that pin the surrounding algorithm without the finder.

- **The kd-tree is a trait.** `MapAlignmentAlgorithmKD` takes any `KdFeatureMaps`; `KDTreeFeatureMaps` implements it (`apply_transformations` followed by `optimize_tree`), and `FeatureGroupingAlgorithmKD` drives the aligner on it. The class and oracle tests of this header use a brute-force stand-in (`tests/support/kd_pool.rs`), whose answers are those of the tree because the algorithm does not depend on the order of a neighbourhood; `tests/mapmatching_kd_grouping*.rs` run the same algorithm on the real tree. A second, parallel port of the header (`map_alignment_algorithm_kd.rs` of lane M2, which fitted LOWESS without `auto_span` and `akima`) was collapsed into this one at integration. The fit-point count per map is bounded by `MAX_FIT_POINTS`.
- **KD refusals.** An algorithm built from an empty `Param` (whose tolerances the source leaves uninitialised), a feature of a map index past `num_maps` (the source writes past `fit_data_`), `transform` before `fit_lowess` (null model pointers) and a `warp:max_nr_conflicts` outside `int` are errors.
- **The guide tree.** `align` builds it with `build_tree`; `align_with_tree` takes the merges a caller supplies (`(left_child, right_child)` of each `BinaryTreeNode`, dummy nodes included) and does not log the Newick line.
- **Tree-guided refusals.** No maps, a range list of another length than the maps, a merge naming a map out of range or the same map twice (undefined in the source), and merging features that refer to an identification graph (the source merges `id_data_`, which the native feature map does not own) are errors. The similarity of two maps without a shared sequence is the source's `Exception::InvalidRange`; a single shared sequence gives NaN, positive here where the source's `0/0` is negative.
- **Unique-ID conflicts** are resolved with IDs from the caller's `UniqueIdGenerator`, where the source draws from its process-wide random generator.

## Checked boundaries and evidence

The superimposers, the identification alignment, the evaluation and the transformation models are tier 1: `tests/map_alignment_oracle.rs` (6 tests, features `idxml` and `consensusxml`) compares 3780 values with records of the executed Linux x86_64 Release build (`openms4-release-bc9cc12-c19e494-174b576`, conda-forge gcc 14.4.0, run on ibminode06 in two invocations of two fresh-process runs, all byte-identical), and on dax all 3780 are bit-identical. The cases are the shift superimposer on five synthetic map pairs (default, dense, narrow buckets, swapped, 100 used points) and the affine one on seven, through both overloads; the identification alignment of the two upstream idXML files without reference, with either as reference, with a score cutoff (`min_score` 30) and with `max_rt_shift` 0 and 2; precision and recall at three tolerance settings on the upstream consensusXML pair, including one where the source's precision is 1.5; and the models listed in `docs/TRANSFORMATIONS_SUPPORT.md`. The records are tests/data/map_alignment_oracle_records.tsv and _summaries.txt; the driver, the input generator and the hashes of the build are recorded in [map_alignment_provenance.json](../tests/data/map_alignment_provenance.json). Off Linux x86_64 the comparison allows `1e-12` relative, for the platform libm's `ln` and `exp` in the affine superimposer's log-scaled histogram.

The remaining checks are tier 3 (the pinned class tests' literals and data files) or tier 4 (derived).

- `tests/pose_clustering.rs` (7 tests): the shift superimposer's slope 1 and intercept -20.4; the affine superimposer's slope 1 and intercept -0.4 through both overloads; its `num_used_points` 2 and 3 cases on the reversed five-point maps (the 3-point case, slope 0.977273 and intercept -0.368182, depends on the `nth_element` selection); the wrong-parameter case and its four warnings; the parameter layouts and restrictions; the refusals; the shift superimposer's two dump files, line by line, derived from the source's dump statements; and a derived check of the pose-clustering alignment with a positional pair finder.
- `tests/map_alignment_evaluation.rs` (5 tests, feature `consensusxml`): precision 0.757143 and recall 0.5 on the upstream consensusXML pair, recall 0 for an empty tool map and the refusal of an empty ground truth, `isSameHandle` with and without charge, and the implementable interface. The precision NaN of an unusable ground truth is derived.
- `tests/map_alignment_identification.rs` (5 tests, feature `idxml`): alignment of the two upstream idXML files without reference (10 points each, common scale), with internal reference (identity, 10 points on the reference scale), with external reference (1 transformation, 10 points), the repeated-align regression with `min_run_occur` 3; and derived cases for feature maps (`use_feature_rt`, the score filter and the empty-sequence record) and for identification data (the automatically picked score type and its `Using score type` line, the best match per observation, a missing named score type).

`MapAlignmentAlgorithmKD` and `MapAlignmentAlgorithmTreeGuided` are tier 1 too: `tests/map_alignment_oracle_kd_tree.rs` (3 tests, features `paramxml` and `featurexml`) compares with `../oracle/M1-kdtree`, run on ibminode06 against the same Release install in two fresh-process runs, byte-identical. KD: eight cases on four synthetic pools (4 maps of 140 peptides with noise; with 30% duplicates under `max_nr_conflicts` 0, 1 and -1; a Da tolerance without fold-change filter; 30 peptides, below the 50-point LOWESS minimum; LOWESS span 0.3; two m/z partitions), run in C++ on the real `KDTreeFeatureMaps` with `FeatureGroupingAlgorithmKD`'s defaults: every neighbourhood size, component label, filtered component size and transformed RT, and the LOWESS warnings in order. Tree-guided: the class-test maps (`b_spline` and `linear`) and synthetic sets of 4 and 5 maps (`b_spline`, `linear`, `lowess`, `interpolated`): every pair's two-map clustering distance against `1 - similarity`, the first merge, the RT lists, and, following the C++ tree, the merge order, every merged feature's RT and `original_RT`, every transformation's data points and notes, and every final RT; and the defaults, written as paramXML. This oracle found that the least-squares `linear` model rounded its means differently from GTE's `ApprHeightLine2` (`sum * (1 / n)`), one ulp in the fitted line; `LinearModel::fit` now divides as GTE does.

- `tests/map_alignment_kd.rs` (10 tests): the class test's construction with 42 maps and an empty `Param` (tier 3; its other sections are `NOT_TESTABLE`); derived components, filters, fit data, LOWESS fallback and refusals.
- `tests/map_alignment_tree_guided.rs` (7 tests, feature `featurexml`): the class test on its three featureXML inputs (tier 3): identification counts 6, 5 and 5; the tree distances 1.84834e-04 and 0.505752 from the similarity (the second through the average-linkage update); following the tree `(0, 2), (0, 1)`, 15 merged features in the order 1, 2, 0 with unique IDs 20 first and 14 last and `original_RT` on all but the first five; data points at each map's median feature RT per sequence, not all identities; `original_RT` after `computeTransformedFeatureMaps`; and `TransformationModelDefaults_test.cpp`'s second section, the `model:` subtree equal to `getDefaults("b_spline")`. Derived: the similarity's arithmetic, refusals and unique-ID conflicts.

- `tests/map_alignment_pose_clustering_class.rs` (1 test, feature `mzml`): `MapAlignmentAlgorithmPoseClustering_test.cpp`'s `align(PeakMap)` section end to end with the composed `StablePairFinder`: 307 data points, slope 1.01164, intercept -32.0912 (tier 3; its other two sections are `NOT_TESTABLE` upstream).
- `tests/map_alignment_oracle_kd_tree.rs` also replays `build_tree` (every merge and distance of the Release build's tree) and `align` end to end (tier 1).
