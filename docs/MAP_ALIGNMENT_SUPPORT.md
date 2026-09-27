# Map alignment

`analysis::mapmatching` holds the map-alignment headers of `ANALYSIS/MAPMATCHING`: the superimposers, the alignment algorithms and the alignment evaluation. The transformation models they produce are in `analysis::transformations` ([TRANSFORMATIONS_SUPPORT.md](TRANSFORMATIONS_SUPPORT.md)), and applying a transformation to a map is `analysis::alignment_transformer`. The reference is OpenMS4-core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`; hashes are in [map_alignment_provenance.json](../tests/data/map_alignment_provenance.json). No dependency is added.

| Header | Rust | Status |
|---|---|---|
| `BaseSuperimposer.h` | `mapmatching/superimposer.rs` | complete |
| `PoseClusteringShiftSuperimposer.h` | `mapmatching/pose_clustering_shift.rs` | partial: the `[DEBUG]` dump files |
| `PoseClusteringAffineSuperimposer.h` | `mapmatching/pose_clustering_affine.rs` | partial: the `[DEBUG]` dump files |
| `MapAlignmentEvaluationAlgorithm.h` | `mapmatching/alignment_evaluation.rs` | complete |
| `MapAlignmentEvaluationAlgorithmPrecision.h` | `mapmatching/alignment_evaluation.rs` | complete |
| `MapAlignmentEvaluationAlgorithmRecall.h` | `mapmatching/alignment_evaluation.rs` | complete |
| `MapAlignmentAlgorithmIdentification.h` | `mapmatching/alignment_identification.rs` | complete |
| `MapAlignmentAlgorithmPoseClustering.h` | `mapmatching/alignment_pose_clustering.rs` | partial: `StablePairFinder` |
| `MapAlignmentAlgorithmKD.h` | none | not started: `KDTreeFeatureMaps` |
| `MapAlignmentAlgorithmTreeGuided.h` | none | not started: `ClusterHierarchical`, `AverageLinkage`, `ClusterAnalyzer::newickTree` |

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
| constructor: `superimposer:` and `pairfinder:` sections, `max_num_peaks_considered` | `Default`; `defaults()`, `to_param()`, `set_parameters()`; the `pairfinder:` section is empty here and carried through unread |
| `setReference<MapType>(map)` | `set_reference_features`, `set_reference_peaks`, `set_reference_consensus` |
| `align(const FeatureMap&, ...)`, `align(const PeakMap&, ...)`, `align(const ConsensusMap&, ...)` | `align_features`, `align_peaks`, `align_consensus`, each taking the pair finder as a closure |
| `updateMembers_` | inside `set_parameters` |
| protected `superimposer_`, `pairfinder_`, `reference_`, `max_num_peaks_considered_` | `superimposer`, `pair_finder` (its parameters), `reference()`, `max_num_peaks_considered` |

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

## Native differences

- **The debug dumps are refused.** A non-empty `dump_buckets` or `dump_pairs` returns `Error::Unsupported`. The source writes several text files per call, named with a process-global invocation counter; they are diagnostic output only and do not affect the result. This is the remaining gap of both superimposers.
- **Undefined behaviour is refused.** An empty map (the source loops past the end of the model map in the shift superimposer), a zero bucket size (an infinite bucket count converted to `Int`), a histogram whose noise level the source can only find by reading past its buffer, a pair similarity or intensity ratio that is not finite (two zero intensities give `0/0`), a first feature without identification under `score_cutoff` (the source reads `features[0].getPeptideIdentifications()[0]` unchecked), and a consensus feature without handles in pose-clustering alignment (the source writes through an end iterator) are errors.
- **Identifications without retention time are refused** in the identification alignment, where the source records NaN, which then decides medians and sorts arbitrarily. With `use_feature_rt` a missing RT never becomes the closest identification, as NaN never is in the source.
- **Ranges are computed, not trusted.** The source states `updateRanges()` as a precondition of the superimposers "for performance reasons"; the ranges are computed from the maps here (with handles, as `ConsensusMap::updateRanges` does).
- **A failed call ends its progress section** (`ProgressReporter::section`), where the source's exception skips `endProgress`.
- **Precision of an unusable ground truth** is the source's NaN (`1.0 / 0 * 0`); recall refuses it, as the source does since its guard was added.
- **Resource limits.** Histogram sizes and the number of quadruples hashed by the affine superimposer are bounded (`max_buckets`, `max_quadruples`); the source bounds neither.
- **The pair finder is a parameter.** `StablePairFinder` belongs to the feature-grouping headers, which this lane does not port; `MapAlignmentAlgorithmPoseClustering` takes the pairing step as a closure `&[ConsensusMap; 2] -> Result<ConsensusMap>` so that it plugs in unchanged when the finder exists, and its `pairfinder:` defaults are empty until then.

## Checked boundaries and evidence

All evidence is tier 3 (the pinned class tests' literals and data files) or tier 4 (derived); no C++ was executed.

- `tests/pose_clustering.rs` (6 tests): the shift superimposer's slope 1 and intercept -20.4; the affine superimposer's slope 1 and intercept -0.4 through both overloads; its `num_used_points` 2 and 3 cases on the reversed five-point maps (the 3-point case, slope 0.977273 and intercept -0.368182, depends on the `nth_element` selection); the wrong-parameter case and its four warnings; the parameter layouts and restrictions; the refusals; and a derived check of the pose-clustering alignment with a positional pair finder.
- `tests/map_alignment_evaluation.rs` (5 tests, feature `consensusxml`): precision 0.757143 and recall 0.5 on the upstream consensusXML pair, recall 0 for an empty tool map and the refusal of an empty ground truth, `isSameHandle` with and without charge, and the implementable interface. The precision NaN of an unusable ground truth is derived.
- `tests/map_alignment_identification.rs` (5 tests, feature `idxml`): alignment of the two upstream idXML files without reference (10 points each, common scale), with internal reference (identity, 10 points on the reference scale), with external reference (1 transformation, 10 points), the repeated-align regression with `min_run_occur` 3; and derived cases for feature maps (`use_feature_rt`, the score filter and the empty-sequence record) and for identification data (the automatically picked score type and its `Using score type` line, the best match per observation, a missing named score type).

Not run: `MapAlignmentAlgorithmPoseClustering_test.cpp` (307 pairs, slope 1.01164, intercept -32.0912), which needs `StablePairFinder`; the `TransformationModelDefaults` test's second section, which needs `MapAlignmentAlgorithmTreeGuided`.
