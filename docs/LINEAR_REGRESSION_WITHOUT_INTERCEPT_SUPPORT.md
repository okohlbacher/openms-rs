# LinearRegressionWithoutIntercept support (`ML/REGRESSION`)

Source revision: OpenMS4-core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`; hashes
in `tests/data/ml_linear_regression_without_intercept_provenance.json`.

| Source header | Rust | Status |
|---|---|---|
| `ML/REGRESSION/LinearRegressionWithoutIntercept.h` | `src/math/ml/linear_regression_without_intercept.rs` | complete |

The header is a least-squares fit of `y = c1 * x` through the origin. It needs
nothing above `math`, so it lives under `src/math/ml/` (`math` is a sink in the
module graph). `LinearRegression.h` and `QuadraticRegression.h`, the other two
headers of `ML/REGRESSION`, are not ported (backlog; `LinearRegression.h` has a
direct TOPP consumer in FeatureFinderMultiplex's source tree, but the
multiplex port reaches only this header's computation).

## API mapping

| Source | Rust |
|---|---|
| `LinearRegressionWithoutIntercept()` (`sum_xx_ = sum_xy_ = 0`, `n_ = 0`) | `new()`, `Default` |
| `addData(double x, double y)` | `add_data(x, y)` |
| `addData(std::vector<double>& x, std::vector<double>& y)` | `add_data_slices(&[f64], &[f64]) -> Result<()>` |
| `getSlope() const` (NaN for `n_ < 2`, else `sum_xy_ / sum_xx_`) | `slope() -> f64` |
| private `sum_xx_`, `sum_xy_`, `n_` | private fields; `n` is a `u64`, the source's `int` |
| destructor, implicit copy | ordinary drop, `Clone`, `Copy` |

## Preserved source conventions

- Fewer than two observations give NaN, whatever they are.
- `sum(x * x) = 0` with two or more observations divides by zero: NaN
  (`0 / 0`) or an infinity, as the source does.
- The observations are summed in insertion order, in `f64`, with the products
  formed before the additions (`sum_xx_ += x * x`).

## Native differences

- `add_data_slices` refuses a dependent vector shorter than the independent
  one; the source loops over `x` and reads `y[i]` unchecked (undefined
  behaviour). The accumulator is unchanged on that error. A longer `y` is
  ignored beyond `x.len()`, as in the source.
- The observation count cannot overflow (`u64`); the source's `int` would.

## Evidence

Tier 3: `tests/ml_linear_regression_without_intercept.rs` transcribes every
section of `LinearRegressionWithoutIntercept_test.cpp`. The sums and the
quotient are one operation each, so no executed differential was run.

## Consumers

`FeatureFinderMultiplexAlgorithm` estimates the intensity ratio of a peptide
pair as this slope (`src/analysis/featurefinder/multiplex/algorithm.rs`); a
private copy of the formula there was collapsed onto this type at integration.
