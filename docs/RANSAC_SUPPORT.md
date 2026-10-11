# Regression, RANSAC and RT normalization: native header equivalents

This document covers eight headers of core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4` that form one dependency chain:

| Header | Rust | Status |
| --- | --- | --- |
| `ML/REGRESSION/LinearRegression.h` (+ `.cpp`) | [`src/math/ml/linear_regression.rs`](../src/math/ml/linear_regression.rs) | every member; `getTValue`/`getLower`/`getUpper` not bit-identical (see below) |
| `ML/REGRESSION/LinearRegressionWithoutIntercept.h` (+ `.cpp`) | [`src/math/ml/linear_regression_without_intercept.rs`](../src/math/ml/linear_regression_without_intercept.rs) | every member |
| `ML/REGRESSION/QuadraticRegression.h` (+ `.cpp`) | [`src/math/ml/quadratic_regression.rs`](../src/math/ml/quadratic_regression.rs) | every member |
| `ML/RANSAC/RANSACModel.h` | [`src/math/ml/ransac.rs`](../src/math/ml/ransac.rs) | every member |
| `ML/RANSAC/RANSACModelLinear.h` (+ `.cpp`) | `src/math/ml/ransac.rs` | every member |
| `ML/RANSAC/RANSACModelQuadratic.h` (+ `.cpp`) | `src/math/ml/ransac.rs` | every member |
| `ML/RANSAC/RANSAC.h` (header-only) | `src/math/ml/ransac.rs` | every member |
| `ANALYSIS/OPENSWATH/MRMRTNormalizer.h` (+ `.cpp`) | [`src/analysis/openswath/mrm_rt_normalizer.rs`](../src/analysis/openswath/mrm_rt_normalizer.rs) | every member |

`RANSAC` draws its samples through `Math::RandomShuffler` (`MATH/MathFunctions.h`),
ported as [`src/math/random_shuffler.rs`](../src/math/random_shuffler.rs).

Tests: [`tests/ransac.rs`](../tests/ransac.rs) (class-test literals),
[`tests/mrm_rt_normalizer.rs`](../tests/mrm_rt_normalizer.rs) (class-test
literals plus independently derived cases) and
[`tests/ransac_oracle.rs`](../tests/ransac_oracle.rs) (executed differential).
Manifest: [`tests/data/ransac_provenance.json`](../tests/data/ransac_provenance.json).

The regression and RANSAC types live in `math::ml` because they need nothing
above `math`, which keeps `math` a sink of the module graph and lets
`processing` (`InternalCalibration`, `MZTrafoModel`) use them later without a
cycle. `MRMRTNormalizer` needs only `math` too, but belongs to OpenSWATH.

## API mapping

### `LinearRegression`

| Source member | Native representation |
| --- | --- |
| `LinearRegression()` | `LinearRegression::new`, `Default` (every value `0.0`) |
| `~LinearRegression()` | drop |
| copy constructor / `operator=` (declared private, not implemented) | `Clone` and `PartialEq` are derived; the source forbids copies only because it never needed them |
| `computeRegression(P, x_begin, x_end, y_begin, compute_goodness = true)` | `compute_regression(&mut self, p, &[f64], &[f64], bool) -> Result<()>` |
| `computeRegressionWeighted(P, x_begin, x_end, y_begin, w_begin, compute_goodness = true)` | `compute_regression_weighted(&mut self, p, x, y, w, bool) -> Result<()>` |
| `getIntercept` / `getSlope` / `getXIntercept` | `intercept` / `slope` / `x_intercept` |
| `getLower` / `getUpper` / `getTValue` | `lower` / `upper` / `t_value` |
| `getRSquared` / `getStandDevRes` / `getMeanRes` | `r_squared` / `stand_dev_res` / `mean_res` |
| `getStandErrSlope` / `getChiSquared` / `getRSD` | `stand_err_slope` / `chi_squared` / `rsd` |
| `static computePointY(x, slope, intercept)` | `LinearRegression::compute_point_y` |
| protected `computeGoodness_` | private `compute_goodness` |
| protected `computeChiSquare` / `computeWeightedChiSquare` (templates) | private `chi_square` / `weighted_chi_square` |
| protected data members `intercept_` ... `rsd_` | private fields behind the getters |

### `LinearRegressionWithoutIntercept`

| Source member | Native representation |
| --- | --- |
| `LinearRegressionWithoutIntercept()` | `new`, `Default` |
| `addData(double x, double y)` | `add_data(x, y)` |
| `addData(std::vector<double>& x, std::vector<double>& y)` | `add_data_slices(&[f64], &[f64]) -> Result<()>` (distinct name for the overload) |
| `getSlope()` | `slope() -> Option<f64>` |
| private `sum_xx_`, `sum_xy_`, `n_` | private fields; native `len` / `is_empty` read the count |

### `QuadraticRegression`

| Source member | Native representation |
| --- | --- |
| `QuadraticRegression()` | `new`, `Default` |
| `computeRegression(x_begin, x_end, y_begin)` | `compute_regression(&mut self, x, y) -> Result<()>` |
| `computeRegressionWeighted(x_begin, x_end, y_begin, w_begin)` | `compute_regression_weighted(&mut self, x, y, w) -> Result<()>` |
| `eval(double x) const` | `eval(&self, x)` |
| `static eval(A, B, C, x)` | `QuadraticRegression::eval_with(a, b, c, x)` (distinct name for the overload) |
| `getA` / `getB` / `getC` / `getChiSquared` | `a` / `b` / `c` / `chi_squared` |
| free `computeChiSquareWeighted` (in the `.cpp`, not the header) | inlined in `compute_regression_weighted` |

### `RansacModel`, `RansacModelLinear`, `RansacModelQuadratic`

| Source member | Native representation |
| --- | --- |
| `template<class ModelT = int> class RansacModel` (CRTP base) | trait `RansacModel` |
| typedef `DPair` | `type DPair = (f64, f64)` |
| typedefs `DVec`, `DVecIt` | `&[DPair]` in, `Vec<DPair>` out |
| typedef `ModelParameters` | `type ModelParameters = Vec<f64>` |
| `rm_fit(begin, end)` / `rm_fit_impl` | `RansacModel::rm_fit(&self, &[DPair]) -> Result<ModelParameters>` |
| `rm_rsq(begin, end)` / `rm_rsq_impl` | `RansacModel::rm_rsq(&self, &[DPair]) -> Result<f64>` |
| `rm_rss(begin, end, coefficients)` / `rm_rss_impl` | `RansacModel::rm_rss(&self, &[DPair], &[f64]) -> Result<f64>` |
| `rm_inliers(begin, end, coefficients, max_threshold)` / `rm_inliers_impl` | `RansacModel::rm_inliers(&self, &[DPair], &[f64], f64) -> Result<Vec<DPair>>` |
| `class RansacModelLinear : public RansacModel<RansacModelLinear>` | unit struct `RansacModelLinear` implementing the trait |
| `class RansacModelQuadratic : public RansacModel<RansacModelQuadratic>` | unit struct `RansacModelQuadratic` implementing the trait |

The source's `*_impl` members are `static`; the unit structs are stateless, so
`RansacModelLinear.rm_fit(..)` is the same call.

### `RANSACParam`, `RANSAC`

| Source member | Native representation |
| --- | --- |
| `struct RANSACParam` | `RansacParam` (public fields `n`, `k`, `t`, `d`, `relative_d`) |
| `RANSACParam()` | `Default` |
| `RANSACParam(n, k, t, d, relative_d = false)` | `RansacParam::new(..) -> Result<Self>` |
| `RANSACParam::toString()` | `Display` |
| `template<typename TModelType = RansacModelLinear> class RANSAC` | `Ransac<M = RansacModelLinear>` |
| `explicit RANSAC(uint64_t seed = time(nullptr))` | `Ransac::new(seed: u64)`; no time-seeded default (see below) |
| `~RANSAC()` | drop |
| `setSeed(uint64_t)` | `set_seed(u64)` |
| `ransac(pairs, const RANSACParam&)` | `ransac_with_param(&mut self, &[DPair], &RansacParam)` |
| `ransac(pairs, n, k, t, d, relative_d = false)` | `ransac(&mut self, &[DPair], n, k, t, d, relative_d) -> Result<Vec<DPair>>` |
| private `shuffler_` | private field; native `Ransac::with_model(model, shuffler)` sets it explicitly |

### `Math::RandomShuffler` (the class of `MATH/MathFunctions.h`)

| Source member | Native representation |
| --- | --- |
| `explicit RandomShuffler(int seed)` | `RandomShuffler::new(i32)` |
| `explicit RandomShuffler(const boost::mt19937_64&)` | `RandomShuffler::from_engine(Mt64)` |
| `RandomShuffler()` | `Default` (seed 5489, `RandomShuffler::DEFAULT_SEED`) |
| public member `rng_` | `engine` / `engine_mut` / `into_engine`, and `next_u64` for `rng_()` |
| `portable_random_shuffle(first, last)` | `portable_random_shuffle(&mut [T])` |
| `seed(uint64_t)` | `seed(u64)` |
| — | native `with_seed(u64)`: any 64-bit seed |

### `MRMRTNormalizer`

| Source member | Native representation |
| --- | --- |
| `static removeOutliersRANSAC(pairs, rsq_limit, coverage_limit, max_iterations, max_rt_threshold, sampling_size)` | `remove_outliers_ransac(.., seed: u64)` |
| `static removeOutliersIterative(pairs, rsq_limit, coverage_limit, use_chauvenet, method)` | `remove_outliers_iterative(.., method: OutlierMethod)`; `OutlierMethod: FromStr` parses the source's names |
| `static chauvenet_probability(residuals, pos)` | `chauvenet_probability(&[f64], usize) -> Result<f64>` |
| `static chauvenet(residuals, pos)` | `chauvenet(&[f64], usize) -> Result<bool>` |
| `static computeBinnedCoverage(rtRange, pairs, nrBins, minPeptidesPerBin, minBinsFilled)` | `compute_binned_coverage((f64, f64), .., i32, i32, i32) -> Result<bool>` |
| protected `static jackknifeOutlierCandidate_(x, y)` | public `jackknife_outlier_candidate` (the class test reaches it through a subclass) |
| protected `static residualOutlierCandidate_(x, y)` | public `residual_outlier_candidate` |

## Preserved source conventions

- **Fit arithmetic.** The unweighted line is `gte::ApprHeightLine2`: centroid by
  summing and multiplying by the reciprocal of the count, slope `Sxy / Sxx`
  about that centroid, intercept `-slope * mean_x + mean_y`; failure when the
  centroid is not finite or `Sxx` is not positive. The weighted line and the
  quadratic solve their normal equations through `gte::LinearSystem`'s
  explicit 2x2 and cofactor 3x3 inverses, with the same statement order,
  a zero determinant as the only singularity test, and each matrix-vector row
  sum started at `0.0`. Chi-squared, the goodness statistics and every model
  residual follow the source's operation order. All of this is bit-identical
  to the Release build on every oracle case.
- **Goodness statistics** are computed only when requested *and* for more than
  two points; otherwise they keep the previous fit's values. A weighted fit's
  statistics are those of the unweighted points (the source passes bare `x`
  and `y` to `computeGoodness_`). `RansacModelLinear::rm_rsq` of exactly two
  points is therefore `0.0`.
- **`RansacModelQuadratic::rm_rsq` returns the chi-squared sum**, not an
  R-squared, as the source does.
- **Inliers** are points with squared residual strictly *below* the threshold
  (the base class comment says `<=`; both models use `<`).
- **RANSAC.** The whole input is reshuffled cumulatively every iteration; the
  first `n` points are the sample; a set is accepted when more than `d` other
  points, or all of them, are inliers; a larger set always wins and equal sizes
  are decided by the residual sum of squares; iterations stop once a set holds
  every point. Only the initial sample's fit is guarded (the source's
  `catch (...)`); refitting an accepted set propagates its error.
- **Seeds.** `RANSAC(uint64_t seed)` hands the seed to `RandomShuffler(int)`, so
  only the low 32 bits count, read as signed and sign-extended; `setSeed`
  keeps all 64 bits. Both are reproduced, and the oracle's seeds
  `4294967297` and `18446744073709551615` exercise the narrowing.
- **`removeOutliersIterative`** evaluates Chauvenet's criterion on the residuals
  of **all input pairs**, at the candidate's index **among the remaining
  points**. After the first removal these refer to different peptides. The
  class test's expectations depend on it; it is reproduced and listed below as
  a defect candidate.
- **`computeBinnedCoverage`** clamps positions at or beyond the range end into
  the last bin and truncates toward zero, so positions just below the range
  fall into the first bin.

## Native differences

- **Student's t quantile.** Boost computes `quantile(students_t(df), p)` in
  `long double` (Hill's approximation plus a Halley step on the incomplete beta
  function). The port starts from `statrs`'s inversion (AS 64/109, `statrs`
  is already a dependency) and polishes it with Newton steps on the tail
  probability. Unpolished, `statrs` was measured up to `2.4e-13` relative from
  Boost on `t_value` and `3.2e-12` on `lower`, which failed the `1e-12`
  relative bound the differential fixed before its first run; polished, the
  maxima are `5.6e-14` (`t_value`) and `5.7e-13` (`lower`/`upper`). The
  remaining difference is `statrs`'s incomplete beta function, which limits the
  tail probability to about `1e-14` relative. `getLower` and `getUpper` inherit the
  difference, amplified by the interval's `1 / (1 - g)`. Nothing in the crate
  consumes these three values; slope, intercept and R-squared are exact.
- **No time-seeded RANSAC.** The source's default seed is `time(nullptr)`, and
  `removeOutliersRANSAC` always uses it, so its result can change from one
  second to the next. `Ransac::new` requires a seed, and
  `remove_outliers_ransac` takes one; passing the value the source's clock
  returned reproduces the source (the oracle interposes `time` to prove it).
- **Checked instead of undefined.** Mismatched slice lengths
  (`InvalidRange`), a Chauvenet position outside the residuals, an empty
  residual vector, a non-positive bin count, a position below the first bin or
  a non-finite normalized position, a negative or non-finite coverage limit
  and model coefficient vectors that are too short are refused; the source
  reads out of bounds, divides by zero or converts with undefined behaviour.
- **Non-finite input** is refused by the weighted linear fit, the quadratic
  fit and RANSAC (`InvalidValue`); the source returns NaN coefficients there.
  The unweighted linear fit refuses it through the source's own rule (a
  non-finite centroid).
- **Atomic failures.** A failed fit leaves the regression unchanged; the source
  has already overwritten slope, intercept and chi-squared when it throws.
- **Method names** parse up front (`OutlierMethod`), so an unknown name is
  refused even for data that already meets the R-squared limit, where the
  source never looks at it.
- **Work bound.** `Ransac::ransac` refuses `k * pairs.len()` above
  `MAX_RANSAC_WORK` (`2^32`) before allocating.
- **Diagnostics.** The source's `OPENMS_LOG_DEBUG` lines, the
  `computeBinnedCoverage` bin warning on standard error and the negative-RSD
  dump on standard output are not reproduced.
- `LinearRegressionWithoutIntercept::slope` returns `None` where the source
  returns NaN for fewer than two observations; its observation count saturates
  instead of overflowing an `int`.
- `RansacParam`'s `Display` formats `t` as Rust does; the source's
  `std::stringstream` uses six significant digits.

## Checked boundaries and evidence

| Tier | Source | What it covers |
| --- | --- | --- |
| 1, executed differential | `../oracle/O4b-ransac` (driver sha256 in the manifest), Release install `openms4-release-bc9cc12-c19e494-174b576`, GCC 14, run twice on ibminode06 with byte-identical output | 128 cases, about 2,100 records: every regression getter, both RANSAC models, RANSAC with 16 seed/parameter combinations (full output order), relative `d`, `setSeed`, both preconditions, jackknife and residual candidates, 28 `removeOutliersIterative` runs, 16 Chauvenet probabilities, 5 coverage cases, 9 `removeOutliersRANSAC` runs with interposed `time()` |
| 3, source review | the class tests `LinearRegression_test.cpp`, `LinearRegressionWithoutIntercept_test.cpp`, `QuadraticRegression_test.cpp`, `RANSACModelLinear_test.cpp`, `RANSACModelQuadratic_test.cpp`, `MRMRTNormalizer_test.cpp` | every literal, with the class tests' `TEST_REAL_SIMILAR` tolerance |
| 4, independently derived | published two-sided 95 % t table; hand-derived coverage bins; a 30-point line with three gross outliers for `removeOutliersRANSAC` | the quantile's order of magnitude, bin arithmetic, RANSAC acceptance logic |

All doubles in the differential compare bit for bit except `t_value`, `lower`
and `upper` (see above).

## Candidate C++ issues

1. `MRMRTNormalizer::removeOutliersIterative` indexes the all-pairs residual
   vector with the candidate's position among the remaining points
   (`MRMRTNormalizer.cpp`, the `chauvenet(residuals, pos)` call): after one
   removal, Chauvenet's criterion is evaluated for a different peptide than the
   one removed.
2. `MRMRTNormalizer::removeOutliersRANSAC` seeds RANSAC with `time(nullptr)`, so
   the RT normalization of `OpenSwathWorkflow` is not reproducible between runs
   that straddle a second boundary.
3. `MRMRTNormalizer::computeBinnedCoverage` writes `binCounter[bin]` with a
   negative `bin` for positions more than one bin below the range, and with an
   arbitrary index for a zero-width range (division by zero); only the upper
   side is clamped.
4. `RANSAC(uint64_t seed)` silently narrows the seed to `int`, so seeds that
   differ only above bit 31 give the same permutations; `setSeed` does not
   narrow.
5. `RansacModel::rm_inliers`' comment says "error <= max_threshold"; both
   models test `<`.
6. `LinearRegression::computeRegression` and `QuadraticRegression`'s callers
   pass `&points.front()` of an empty vector to GTE for empty input (undefined
   behaviour before the `UnableToFit`).
