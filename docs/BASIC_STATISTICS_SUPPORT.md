# BasicStatistics: native header equivalent

[`src/math/basic_statistics.rs`](../src/math/basic_statistics.rs) provides a
native equivalent for the complete public surface of core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`
`MATH/STATISTICS/BasicStatistics.h`. The header is header-only — a template
class with no accompanying translation unit — so every body below is inline in
it.

Tests: [`tests/basic_statistics.rs`](../tests/basic_statistics.rs). Manifest:
[`tests/data/math_statistics_provenance.json`](../tests/data/math_statistics_provenance.json).

`BasicStatistics` accumulates **weighted** moments of a distribution: the
weights are the probability values, the abscissae are either the positions
`0, 1, ..., n-1` or an explicit coordinate vector, and the divisor is the
probability mass. It is not a sample-statistics class, so the degrees-of-freedom
question that dominates [`StatisticFunctions`](STATISTIC_FUNCTIONS_SUPPORT.md)
does not arise here.

## API mapping

Every public and protected member of the header appears here.

| Source member | Native representation |
| --- | --- |
| `template <typename RealT = double> class BasicStatistics` | non-generic `BasicStatistics` over `f64`; see differences |
| `typedef RealT RealType` | `f64` |
| `typedef std::vector<RealType> probability_container` | `Vec<f64>` / `&[f64]` |
| `typedef std::vector<RealType> coordinate_container` | `Vec<f64>` / `&[f64]` |
| `BasicStatistics()` | `BasicStatistics::new()`, `Default` |
| `BasicStatistics(BasicStatistics const&)` | `Copy` / `Clone` |
| `operator=(BasicStatistics const&)` | ordinary assignment (`Copy`) |
| `void clear()` | `clear(&mut self)` |
| `update(probability_begin, probability_end)` | `update(&mut self, probabilities: &[f64]) -> Result<()>` |
| `update(probability_begin, probability_end, coordinate_begin)` | `update_with_coordinates(&mut self, probabilities: &[f64], coordinates: &[f64]) -> Result<()>` |
| `RealType mean() const` | `mean(&self) -> f64` |
| `void setMean(RealType const&)` | `set_mean(&mut self, f64)` |
| `RealType variance() const` | `variance(&self) -> f64` |
| `void setVariance(RealType const&)` | `set_variance(&mut self, f64)` |
| `RealType sum() const` | `sum(&self) -> f64` |
| `void setSum(RealType const&)` | `set_sum(&mut self, f64)` |
| `static RealType sqrt2pi()` | `BasicStatistics::sqrt2pi()` and the associated constant `BasicStatistics::SQRT_2PI` |
| `RealType normalDensity_sqrt2pi(RealType) const` | `normal_density_sqrt2pi(&self, f64) -> f64` |
| `RealType normalDensity(RealType) const` | `normal_density(&self, f64) -> f64` |
| `void normalApproximation(probability_container&)` | `normal_approximation(&self, size: usize) -> Result<Vec<f64>>`, called with the existing length |
| `void normalApproximation(probability_container&, size_type size)` | the same `normal_approximation(size)` |
| `void normalApproximation(probability_container&, coordinate_container const&)` | `normal_approximation_at(&self, coordinates: &[f64]) -> Result<Vec<f64>>` |
| `friend operator<<(ostream&, BasicStatistics&)` | `impl Display` with the same field layout |
| protected `mean_`, `variance_`, `sum_` | private fields behind the accessors above |
| private `normalApproximationHelper_(probability, size)` | the body of `normal_approximation` |
| private `normalApproximationHelper_(probability, coordinate)` | the body of `normal_approximation_at` |
| Not in source but added | `MAX_ITEMS` resource ceiling |

The source's three `normalApproximation` overloads become two functions because
the difference between the first two is only whether the output vector is
resized first; a Rust function that returns the vector covers both, and the
caller that wanted "same size as before" passes that length.

## Preserved source conventions

- **Two passes, in order.** `update` accumulates `sum_` and `sum(p_i * i)`
  together in the first pass and divides to get the mean, then accumulates
  `sum(p_i * (i - mean)^2)` in a second pass and divides by the same sum. The
  passes are not fused and no shifted-moment identity is used, because either
  changes the last bits.
- **The variance divisor is the probability mass**, not `n` and not `n - 1`.
- **The zero-mass guard is asymmetric.** The position-based `update` ends with
  `if (sum_ == 0 && (isnan(mean_) || isinf(mean_))) { mean_ = 0; variance_ = 0; }`.
  The coordinate-based overload has no such guard, so a zero probability sum
  leaves both moments NaN there. Both behaviours are reproduced; see the C++
  issue below.
- **`normalDensity_sqrt2pi` is not scaled by `1 / sqrt(variance)`.** It is
  `exp` of the standardised square, so it is exactly `1` at the mean for every
  variance. The class test pins that.
- **`normalApproximation` scales by `sum()`.** Each entry is
  `density(i) / gauss_sum * sum()`, evaluated left to right, so the entries
  carry the distribution's own mass rather than integrating to one.
- **`sqrt2pi()` is the source's literal**, `2.50662827463100050240`. That
  literal parses to the `f64` `2.5066282746310007`, which is **one unit in the
  last place above** the correctly rounded `sqrt(2 * pi)` = `2.5066282746310002`.
  The port reproduces the source's value rather than correcting it, because the
  class test compares against it exactly and every density the library has ever
  produced carries it. `tests/basic_statistics.rs` asserts the one-ulp
  relationship directly, so a future "cleanup" to the true constant fails.

## Native differences

| Difference | Reason |
| --- | --- |
| Not generic over `RealT` | Only the `double` instantiation exists in the library and in the class test. An `f32` instantiation would only lose precision in an accumulation this port performs in `f64` deliberately. |
| `update` and `update_with_coordinates` return `Result` | They carry the `MAX_ITEMS` ceiling, and the coordinate overload checks the two lengths. |
| `update_with_coordinates` requires equal lengths | The source takes only a coordinate *begin* iterator and advances it once per probability, so a short coordinate range is read out of bounds. |
| Both refuse before clearing | A rejected call leaves the previous state intact; the source has nothing to reject. |
| `normal_density_sqrt2pi` branches on a zero variance | The source divides by it. The branch returns exactly what the division produces — `0.0` away from the mean, NaN at it — so this is a spelled-out equivalence, not a behaviour change, and it satisfies the crate's rule that every division is checked. |
| `normal_approximation*` refuse a zero or non-finite normalising sum | The source divides by it and fills the output with NaN or infinities. |
| `normal_approximation*` return a new `Vec` | The source writes into a caller-supplied container and resizes it. |
| `Display` prints full `f64` precision | The source's `operator<<` inherits the stream's default six significant digits. This is a debugging helper on both sides; no caller parses it. |
| `MAX_ITEMS` ceiling | The source resizes to whatever it is handed. |
| Serial only | The header carries no `#pragma omp`. |

## Checked boundaries and evidence

Resource boundaries: `MAX_ITEMS = 50_000_000` probabilities, coordinates or
approximation points, checked before the state is cleared or a vector is
allocated.

Numeric boundaries: an empty probability vector (all three parameters end at
zero, through the same guard as an all-zero vector); a zero probability mass
with and without coordinates; a zero variance in both density functions; a
normalising density sum that is zero or non-finite; mismatched probability and
coordinate lengths.

Evidence, per `docs/DIFFERENTIAL_VALIDATION.md`:

- **Tier 3 (source review, transcribed class-test literals)** for all 18
  `START_SECTION`s of `BasicStatistics_test.cpp`, in
  `tests/basic_statistics.rs`. A plain `grep -c START_SECTION` of that file
  reports 20; two of those occurrences are inside comments at lines 341 and 356,
  which refer to a section rather than opening one. The 18 include the
  195-value sample and its
  `sum = 15228.2`, `mean = 96.4639`, `variance = 3276.51` at
  `TOLERANCE_ABSOLUTE(0.1)`, the six `good_probs` of the normal approximation,
  and the four `normalDensity_sqrt2pi` values at `TOLERANCE_ABSOLUTE(0.0001)`.
- **Tier 4 (independently derived)**: the weighted mean `13/6` and variance
  `17/36` of `{0, 1, 3, 2, 0}`; the reflection identity
  `mean(coordinates 1000 - i) == 1000 - mean(positions)`; the symmetry of
  `normalDensity_sqrt2pi` about the mean, asserted as exact equality rather
  than a tolerance; the fact that the approximation's entries sum to `sum()`;
  and the one-ulp relationship between the source's `sqrt2pi` literal and
  `sqrt(2 * pi)`.
- **Tier 4 (Rust-only)** for the length check, the ceiling, and the two
  refusals above.

No tier 1 or tier 2 evidence exists for this group.

## Candidate C++ issue

Reported to the integrating agent rather than written into
`OpenMS_CPP_ISSUES.md`, which that agent owns:

- `BasicStatistics::update(probability_begin, probability_end, coordinate_begin)`
  (`BasicStatistics.h`, lines 113-141) takes no coordinate end iterator and
  advances it once per probability, reading past the end of a shorter
  coordinate range; and, unlike the position-based overload at lines 82-111, it
  has no zero-mass guard, so a zero probability sum silently leaves `mean()` and
  `variance()` NaN.
