# Consensus-map normalisation support (`ConsensusMapNormalizerAlgorithm{Median,Quantile,Threshold}`)

Source revision: OpenMS4-core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

| Source header | Rust | Status |
|---|---|---|
| `ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h` | `src/analysis/mapmatching/consensus_map_normalizer_median.rs` | complete |
| `ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h` | `src/analysis/mapmatching/consensus_map_normalizer_quantile.rs` | complete |
| `ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmThreshold.h` | `src/analysis/mapmatching/consensus_map_normalizer_threshold.rs` | complete |

The three classes are stateless bags of static functions with suppressed copy
and assignment. Each becomes a unit struct with associated functions. The one
TOPP consumer, `ConsensusMapNormalizer`, is not ported yet.

## API mapping

### ConsensusMapNormalizerAlgorithmMedian

| Source | Rust |
|---|---|
| default constructor, destructor; private copy constructor and assignment | unit struct `ConsensusMapNormalizerAlgorithmMedian` (`Clone`/`Copy` are harmless: no state) |
| `enum NormalizationMethod { NM_SCALE, NM_SHIFT }` | `NormalizationMethod::{Scale, Shift}` |
| `static void normalizeMaps(map, method, acc_filter, desc_filter)` | `normalize_maps(&mut map, method, &IdentificationFilter) -> Result<Medians>` |
| `static Size computeMedians(map, medians&, acc_filter, desc_filter)` | `compute_medians(&map, &IdentificationFilter) -> Result<Medians>` (`reference_map` is the return value, `medians` the out-parameter) |
| `static bool passesFilters_(cf_it, map, acc_filter, desc_filter)` | `IdentificationFilter::new(acc, desc)?.passes(&feature, &map) -> Result<bool>` |
| `friend class ConsensusMapNormalizerAlgorithmThreshold` | `IdentificationFilter` is public and shared |

### ConsensusMapNormalizerAlgorithmQuantile

| Source | Rust |
|---|---|
| default constructor, destructor; private copy and assignment | unit struct |
| `static void normalizeMaps(map)` | `normalize_maps(&mut map) -> Result<()>` |
| `static void resample(data_in, data_out&, n)` | `resample(&[f64], usize) -> Result<Vec<f64>>` |
| `static void extractIntensityVectors(map, out&)` | `extract_intensity_vectors(&map) -> Result<Vec<Vec<f64>>>` |
| `static void setNormalizedIntensityValues(feature_ints, map&)` | `set_normalized_intensity_values(&[Vec<f64>], &mut map) -> Result<()>` |

### ConsensusMapNormalizerAlgorithmThreshold

| Source | Rust |
|---|---|
| default constructor, destructor; private copy and assignment | unit struct |
| `static vector<double> computeCorrelation(map, ratio_threshold, acc_filter, desc_filter)` | `compute_correlation(&map, f64, &IdentificationFilter) -> Result<ThresholdFactors>` (`ratios` is the source result) |
| `static void normalizeMaps(map, ratios)` | `normalize_maps(&mut map, &[f64]) -> Result<()>` |

## Preserved source conventions

- The maps are the column headers `0..n`, `n` their count; the reference is
  the argmax of the stored header `size` (strict `>`, first on ties, map 0 when
  all sizes are 0 — OpenMS issue #9864 finding A3, pinned by the class test).
- The filters: an empty expression, or one matching the empty string, passes
  every feature; otherwise accessions (empty ones skipped) are searched, and
  the description of the first protein hit with that accession in each
  protein identification. Boost perl syntax through `BoostRegex`. Text after an
  embedded NUL is invisible, as with the source's `c_str()`.
- Median bail-out: if any map has no passing handle every median is 1.0 and
  the reference is 0. Threshold bail-out: if any map keeps no ratio every factor
  is 1.0.
- Shift targets the largest median, searched from
  `numeric_limits<double>::min()` (the smallest *positive* double), so with all
  medians non-positive the target is map 0.
- Arithmetic: `(intensity * medians[ref]) / medians[j]`,
  `intensity + medians[max] - medians[j]`, `intensity * ratio`, all in `f64`
  and stored as `f32`; quantile reference built as `sum_i value / n_maps` in map
  order; resampling weights `1 - (p - floor)` and `1 - (ceil - p)`; threshold
  factor `Math::mean` (sum in feature order divided by count); a zero
  intensity counts as absent from a ratio; ratio window open at both ends.
- Only handle intensities change; consensus-feature intensities are left
  alone (pinned by the median class test).
- Handle iteration order is (map index, unique id), the source's `IndexLess`
  handle set order, which fixes the quantile write-back order.

## Native differences

- Errors instead of undefined behaviour: a handle whose map index has no
  column header (the source indexes a vector with it); a missing column header
  in the threshold path (dereferenced `end()`); a map without handles in
  quantile normalisation (`front()` of an empty vector); a `ratios` vector too
  short for `normalizeMaps`; a malformed regular expression (uncaught
  `regex_error`).
- A non-finite normalised intensity (a zero median or factor) is refused; the
  source stores the infinity or NaN. Every mutation is transactional.
- The regular expressions are compiled once per call rather than once per
  feature.
- The informational and warning log lines become fields of `Medians` and
  `ThresholdFactors` (`passing_features`, `enough_features`); nothing is
  logged. The `ProgressLogger` of `normalizeMaps` is not ported.
- Ceilings: `MAX_HANDLES` (2^27 handles) and, for the threshold's dense
  matrix, `MAX_MATRIX_CELLS` (2^27 maps x features), checked before
  allocation. The source has none.
- Quantile ranks ties by handle position, as the source's sort of
  `(value, index)` pairs does; `-0.0` and `0.0` compare equal as under `<`.

## Checked boundaries and evidence

- Tier 1: `tests/mapmatching_consensus_map_normalizer_topp.rs` reproduces the
  algorithm path of `TOPP_ConsensusMapNormalizer_1` (test-data 0cb15f2,
  `topp/CMakeLists.txt:337-339`; default `robust_regression`: `sortBySize`,
  `computeCorrelation(0.67, "", "")`, `normalizeMaps`) and matches every
  handle of the retained C++ output within the six significant digits the file
  stores (relative 5e-6 plus one `f32` unit).
- Tier 3: every section of `ConsensusMapNormalizerAlgorithmMedian_test.cpp`
  (`tests/mapmatching_consensus_map_normalizer.rs`).
- Tier 4: the quantile (`TODO` upstream) and threshold (`NOT_TESTABLE`
  upstream) sections, with hand-derived values shown at each assertion; error
  paths and transactional failure.
- Not executed against the C++ Release build: the quantile algorithm has no
  retained output and no oracle yet (a backlog item: an
  `../oracle/M2-normalizer` driver on ibminode06).
