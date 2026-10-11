# OPENSWATHALGO algorithms: Scoring, StatsHelpers, DataFrameWriter

Source pin: OpenMS4-core `bc9cc12`. Rust: `src/analysis/openswath/{scoring,stats_helpers,data_frame_writer}.rs`.
Tests: `tests/openswath_algo.rs` (plus a unit test in `data_frame_writer.rs`). Provenance:
`tests/data/openswath_foundation_provenance.json`. The OPENSWATHALGO data structures and
interfaces are documented in [OPENSWATH_DATAACCESS_SUPPORT.md](OPENSWATH_DATAACCESS_SUPPORT.md).

## API mapping

### `OPENSWATHALGO/ALGO/Scoring.h` (namespace `OpenSwath::Scoring`)

| Source | Rust (`analysis::openswath::scoring`) | Notes |
|---|---|---|
| `typedef pos2D` | `Pos2D = (u32, u32)` | |
| `struct pair_hash` | not ported | Rust tuples implement `Hash` |
| `typedef XCorrEntry` | `XCorrEntry = (i32, f64)` | |
| `struct XCorrArrayType` (`data`, `begin`, `end`) | `XCorrArrayType { data }`, `iter`, `iter_mut`, `IntoIterator for &` | |
| `NormalizedManhattanDist` | `normalized_manhattan_dist(&mut, &mut) -> Result<f64>` | normalises both inputs in place, as the source |
| `RootMeanSquareDeviation` | `root_mean_square_deviation` | |
| `SpectralAngle` | `spectral_angle` | zero-length vector gives pi/2 |
| `normalizedCrossCorrelation` | `normalized_cross_correlation(&mut, &mut, maxdelay, lag)` | standardises in place |
| `normalizedCrossCorrelationPost` | `normalized_cross_correlation_post` | |
| `calculateCrossCorrelation` | `calculate_cross_correlation` | |
| `xcorrArrayGetMaxPeak` | `xcorr_array_get_max_peak -> Result<usize>` | returns an index instead of an iterator |
| `standardize_data` | `standardize_data` | population standard deviation |
| `normalize_sum` | `normalize_sum` | |
| `computeAndAppendRank` | `compute_and_append_rank -> Result<u32>` | replaces `ranks_out`, as the source |
| `computeRankVector` | `compute_rank_vector -> Result<Vec<u32>>` | appends one vector per row |
| `rankedMutualInformation` | `ranked_mutual_information -> Result<f64>` | |

### `OPENSWATHALGO/ALGO/StatsHelpers.h`

| Source | Rust (`analysis::openswath::stats_helpers`) | Notes |
|---|---|---|
| `normalize(in, factor, out)` | `normalize_into` (exact out-parameter form), `normalize` (fresh vector) | |
| `norm<T>` (2 instantiations) | `norm(&[f64])` | |
| `dotProd<Texp,Ttheo>` (double, float, int) | `dot_prod`, `dot_prod_f32`, `dot_prod_i32` | float accumulates in `f32`, int in checked `i32`, as Eigen does in the element type |
| `dotprodScoring` | `dotprod_scoring` | |
| `manhattanDist<Texp,Ttheo>` | `manhattan_dist` | |
| `manhattanScoring` | `manhattan_scoring` | |
| `cor_pearson` (`[[deprecated]]`) | `cor_pearson` (`#[deprecated]`, points at `math::statistic_functions::pearson_correlation_coefficient`) | `double` form only |
| class `mean_and_stddev` | `MeanAndStddev` (`push`, `Extend`, `FromIterator`, `sample_variance`, `standard_variance`, `sample_stddev`, `standard_stddev`, `mean`, `count`, `variance`, `stddev`) | nullary `operator()` is `stddev`; typedefs `argument_type`/`result_type` dropped |

### `OPENSWATHALGO/DATAACCESS/DataFrameWriter.h`

| Source | Rust (`analysis::openswath::data_frame_writer`) | Notes |
|---|---|---|
| `IDataFrameWriter` (`colnames`, `store`) | trait `DataFrameWriter` | both return `Result` |
| `DataMatrix` | `DataMatrix` + native read accessors `column_names`, `row_names`, `rows` | the source has no getters |
| `CSVWriter(filename)`, `~CSVWriter` | `CsvWriter::create`, `CsvWriter::finish` | see differences |

## Preserved source conventions

- Every precondition-free edge case keeps its source result: zero vectors in `spectral_angle` (pi/2),
  all-zero and constant vectors in `standardize_data`, a zero sum in `normalize_sum` (unchanged),
  delays without overlap in the cross-correlation (0), ties in `xcorr_array_get_max_peak` (first wins),
  empty input to `normalized_cross_correlation_post` (0/0 = NaN entries).
- `NormalizedManhattanDist` follows the code, not its Doxygen formula (which shows a square root the
  code does not take); the class test pins the code.
- `compute_and_append_rank` gives equal values the position of their first occurrence in sorted order
  (ranks can skip), and the name notwithstanding replaces its output.
- `ranked_mutual_information` visits joint states in ascending order, so the floating-point sum is the
  source's sum.
- `normalize`: a factor that is not positive leaves the (resized) output untouched, so the in-place use
  in `dotprodScoring`/`manhattanScoring` leaves an all-zero vector as its square roots.
- The TSV writer formats values as `std::setprecision(5)` in the default float format (`%.5g`).

## Native differences

- **Summation order.** Eigen reduces in vector blocks. The `f64` reductions of `StatsHelpers.cpp`
  (`norm`, `dot_prod`, `manhattan_dist` and the sums of `manhattanScoring`) and of `Scoring.cpp`
  (`standardize_data`, `calculate_cross_correlation`) add in Eigen's own order, with the Eigen 5.0.1
  two-lane kernels of the Linux x86-64 Release build that `math::fitters::levenberg_marquardt` carries
  and measured; they agree with that build bit for bit (`tests/openswath_scoring_oracle.rs`,
  `tests/openswath_peak_picking.rs`) and are identical on every machine. The `f32` and `i32` dot
  products accumulate in plain order. No crate was taken for these reductions: `nalgebra` (considered)
  does not reproduce Eigen's block order.
- **Checked preconditions.** `OPENSWATH_PRECONDITION` is compiled out of release builds; the port checks
  every one (equal non-zero lengths, non-empty cross-correlation array, ranks within their maximum) and
  returns `Error::InvalidValue`. The iterator-pair functions check that the second sequence is at least
  as long as the first.
- **Resource ceiling.** A cross-correlation with more than `MAX_XCORR_ENTRIES` (2^24) entries is refused
  before allocation; `lag <= 0` (an endless loop in the source) is refused; a negative `maxdelay` yields
  an empty array.
- **NaN ranks.** The source sorts with `<`, undefined for NaN; the port sorts by IEEE total order (same
  result for every other input), so each NaN gets its own rank.
- **Integer dot product.** Overflow is undefined in the source and an error here.
- **`CSVWriter`.** Creation fails if the file cannot be created (the source silently writes nothing);
  write errors are reported; `finish` replaces the destructor's flush and its "have flushed and closed
  the file stream" line on standard output, which the port does not print.

## Checked boundaries and evidence

| Behaviour | Evidence |
|---|---|
| Scoring: all 13 sections of `openswathalgo/Scoring_test.cpp` (distances, spectral angle incl. 4 NumPy cases, `normalize_sum`, `standardize_data`, cross-correlations incl. full range and the 3-size equivalence against the reference loop, properties, edge cases, ranks, mutual information) | tier 3, `tests/openswath_algo.rs` |
| StatsHelpers: `openms/source/StatsHelpers_test.cpp` (norm, normalise, dot product, Manhattan) | tier 3 |
| Normalisation rules, scoring on hand-derived vectors, `MeanAndStddev` against the textbook data set, `cor_pearson` against the crate's Pearson | tier 4 |
| Guards: lag 0, unequal lengths, entry ceiling, empty arrays, rank above maximum, integer overflow | tier 4 |
| `DataMatrix`, TSV layout and `%.5g` formatting | tier 4 (derived from the C stream rule) |

Measured: `cargo test --all-features --test openswath_algo` on dax, 30 passed.
