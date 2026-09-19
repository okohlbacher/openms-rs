# StatisticFunctions: native header equivalent

[`src/math/statistic_functions.rs`](../src/math/statistic_functions.rs) provides
a native equivalent for the complete public surface of core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4` `MATH/StatisticFunctions.h`. The
accompanying `src/openms/source/MATH/StatisticFunctions.cpp` is an empty
translation unit: the header is header-only and every body below is inline in
it.

The module is registered from [`src/math/mod.rs`](../src/math/mod.rs), together
with [`BasicStatistics`](BASIC_STATISTICS_SUPPORT.md),
[`RankData`](RANK_DATA_SUPPORT.md) and [`Histogram`](HISTOGRAM_SUPPORT.md). Its
tests are [`tests/statistic_functions.rs`](../tests/statistic_functions.rs) and
its provenance manifest is
[`tests/data/math_statistics_provenance.json`](../tests/data/math_statistics_provenance.json).

The header is consumed by eight TOPP tools (`FileInfo`, `MRMPairFinder`,
`MapStatistics`, `MetaProSIP`, `NucleicAcidSearchEngine`, `OpenNuXL`,
`ProteomicsLFQ`, `QCMerger`) and by much of the library, so its conventions
propagate. The two that matter most are recorded below: the **degrees of
freedom**, which are not uniform, and the **median convention** for an even
count.

## API mapping

Every public member of the header appears here.

| Source member | Native representation |
| --- | --- |
| `struct AdaptiveQuantileResult` (`blended`, `half_raw`, `half_rob`, `upper_fence`, `tail_fraction`, `weight`) | `AdaptiveQuantileResult` with the same six public `f64` fields; `Default` reproduces the member initialisers, including `upper_fence = +inf` |
| `checkIteratorsNotNULL(begin, end)` | `check_not_empty(&[T]) -> Result<()>` |
| `checkIteratorsEqual(begin, end)` | `check_exhausted(&[T]) -> Result<()>` |
| `checkIteratorsAreValid(begin_b, end_b, begin_a, end_a)` | `check_ranges_end_together(&[T], &[U]) -> Result<()>`, argument order preserved |
| `sum(begin, end)` | `sum(&[f64]) -> f64` |
| `mean(begin, end)` | `mean(&[f64]) -> Result<f64>` |
| `median(begin, end, sorted = false)` | split: `median_sorted(&[f64])` for `sorted = true`, `median(&mut [f64])` for `sorted = false` (which sorts, as the source does) |
| `MAD(begin, end, median_of_numbers)` | `mad(&[f64], f64) -> Result<f64>` |
| `MeanAbsoluteDeviation(begin, end, mean_of_numbers)` | `mean_absolute_deviation(&[f64], f64) -> f64` |
| `absdev(begin, end, mean = DBL_MAX)` | split: `absdev(&[f64])` computes the mean, `absdev_with_mean(&[f64], f64)` takes it; the `DBL_MAX` sentinel is gone |
| `quantile1st(begin, end, sorted = false)` | split: `quantile1st_sorted(&[f64])`, `quantile1st(&mut [f64])` |
| `quantile3rd(begin, end, sorted = false)` | split: `quantile3rd_sorted(&[f64])`, `quantile3rd(&mut [f64])` |
| `quantile(begin, end, q)` | `quantile(&[f64], f64) -> Result<f64>` |
| `tukeyUpperFence(begin, end, k = 1.5)` | `tukey_upper_fence(&[f64], f64) -> Result<f64>`; the default is the constant `DEFAULT_TUKEY_FACTOR` |
| `tailFractionAbove(begin, end, threshold)` | `tail_fraction_above(&[f64], f64) -> f64` |
| `winsorizedQuantile(begin, end, q, upper_fence)` | `winsorized_quantile(&[f64], f64, f64) -> Result<f64>` |
| `adaptiveQuantile(begin, end, q, k = 1.5, r_sparse = 0.01, r_dense = 0.10)` | `adaptive_quantile(&[f64], f64, f64, f64, f64) -> Result<AdaptiveQuantileResult>`; defaults are `DEFAULT_TUKEY_FACTOR`, `DEFAULT_R_SPARSE`, `DEFAULT_R_DENSE` |
| `variance(begin, end, mean = DBL_MAX)` | split: `variance(&[f64])`, `variance_with_mean(&[f64], f64)` |
| `sd(begin, end, mean = DBL_MAX)` | split: `sd(&[f64])`, `sd_with_mean(&[f64], f64)` |
| `covariance(begin_a, end_a, begin_b, end_b)` | `covariance(&[f64], &[f64]) -> Result<f64>` |
| `meanSquareError(begin_a, end_a, begin_b, end_b)` | `mean_square_error(&[f64], &[f64]) -> Result<f64>` |
| `rootMeanSquareError(begin_a, end_a, begin_b, end_b)` | `root_mean_square_error(&[f64], &[f64]) -> Result<f64>` |
| `classificationRate(begin_a, end_a, begin_b, end_b)` | `classification_rate(&[f64], &[f64]) -> Result<f64>` |
| `matthewsCorrelationCoefficient(begin_a, end_a, begin_b, end_b)` | `matthews_correlation_coefficient(&[f64], &[f64]) -> Result<f64>` |
| `pearsonCorrelationCoefficient(begin_a, end_a, begin_b, end_b)` | `pearson_correlation_coefficient(&[f64], &[f64]) -> Result<f64>` |
| `computeRank(std::vector<Value>& w)` | `compute_rank(&mut [f64]) -> Result<()>`; its tie tolerance is the public `COMPUTE_RANK_TIE_TOLERANCE` |
| `rankCorrelationCoefficient(begin_a, end_a, begin_b, end_b)` | `rank_correlation_coefficient(&[f64], &[f64]) -> Result<f64>` |
| `SummaryStatistics<T>` default constructor | `SummaryStatistics::default()` (all zero, `count = 0`) |
| `SummaryStatistics<T>(T& data)` | `SummaryStatistics::new(&mut [f64]) -> Result<Self>`, which sorts as the source does |
| `SummaryStatistics` fields `mean`, `variance`, `lowerq`, `median`, `upperq`, `min`, `max`, `count` | public fields of the same names; `min`/`max` are `f64` rather than `T::value_type` |
| Not in source but added | `MAX_ITEMS`, `MAX_BYTES` resource ceilings; `DEFAULT_TUKEY_FACTOR`, `DEFAULT_R_SPARSE`, `DEFAULT_R_DENSE`, `COMPUTE_RANK_TIE_TOLERANCE` as named constants for the source's default arguments and literals |

The template parameters (`IteratorType`, `IteratorType1`, `IteratorType2`,
`Value`, `T`) are not carried: the port takes `&[f64]`. The three
`checkIterators*` helpers stay generic over the element type, because their
whole content is a length question and the class test calls one of them on
`std::vector<int>`.

## Preserved source conventions

- **Degrees of freedom, exactly as written.** `variance`, `sd` and `covariance`
  divide by `n - 1`. `meanSquareError` and `MeanAbsoluteDeviation` divide by
  `n`. Nothing was made uniform.
- **Median for an even count** averages the two middle values,
  `(x[n/2 - 1] + x[n/2]) / 2.0`; for an odd count it returns `x[(n-1)/2]`
  unchanged. That is why the return type is floating point.
- **Two different quantile conventions coexist.** `quantile1st`/`quantile3rd`
  are medians of halves and interpolate nothing; `quantile(q)` is
  Hyndman-Fan type 7 with linear interpolation. `tukeyUpperFence` uses the
  latter, so `Q1`/`Q3` inside a Tukey fence are **not** `quantile1st`/
  `quantile3rd`. Both are ported and both keep their own rule.
- **The even-count halves drop an element.** `quantile1st` takes the median of
  `[0, n/2 - 1)` and `quantile3rd` the median of `[n/2 + 1, n)`. The `-1` and
  `+1` are the source's, commented "to exclude median values" — an exclusion
  that is correct for an odd count, where the single middle element really is
  the median, and drops a real observation for an even count, where there is no
  middle element. For `n = 10` the lower half is four values, not five. The port
  reproduces this; see the C++ issue below.
- **Sizes below three** collapse `quantile1st` to the minimum and `quantile3rd`
  to the maximum, which is what the size-3 and size-4 cases produce anyway.
- **Accumulation order.** Pearson forms both means first and then accumulates
  the numerator and both denominator sums together in one pass. `variance`
  accumulates `diff * diff` in slice order. `BasicStatistics`-style fused or
  shifted-moment rewrites are not used anywhere.
- **`computeRank`'s relative tie tolerance.** Two neighbours tie when
  `|x[i+1] - x[i]| <= 1e-7 * |x[i+1]|`. The reference is the *later* value, so
  the tolerance collapses to zero next to a zero; and the test is pairwise
  while scanning, so a chain of near-neighbours becomes one block even when its
  ends are far apart. Both are reproduced.
- **Spearman correlates about the theoretical mean rank** `(n + 1) / 2`, not the
  observed mean of the ranks. The two differ when ties are present. The source's
  comment records the earlier integer division this replaced.
- **`rankCorrelationCoefficient` returns `0.0`**, not NaN, when either range's
  ranks are constant — the source's explicit `if (!sqsum_data || !sqsum_model)`.
  `pearsonCorrelationCoefficient` returns NaN in the same situation, and its
  `@brief` says so. The inconsistency is source behaviour and is preserved.
- **`MeanAbsoluteDeviation` of an empty range is NaN.** The source neither
  checks nor documents it, but the class test asserts the NaN, so the value is
  pinned rather than quietly changed. `absdev` checks and is the safe entry.
- **`winsorizedQuantile` also caps below at `0.0`.** The source calls that
  defensive and useful for absolute residuals; it makes the function unsuitable
  for signed data, and that is documented at the item rather than removed.
- **`adaptiveQuantile` with `r_dense <= r_sparse`** degenerates to the source's
  step function rather than erroring.

## Native differences

| Difference | Reason |
| --- | --- |
| `Result` instead of exceptions | `Exception::InvalidRange` maps to `Error::InvalidRange` (empty range, mismatched lengths); `Exception::InvalidValue` to `Error::InvalidValue` (`q` outside `[0, 1]`). |
| The `sorted` boolean becomes two functions | The `false` case mutates the caller's range. `&mut [f64]` states that in the signature; `&[f64]` proves the other does not. |
| The `mean = DBL_MAX` sentinel becomes two functions | A caller that genuinely wants the deviation about `DBL_MAX` cannot express it in the source. |
| Sortedness is always checked | `median_sorted`, `quantile1st_sorted`, `quantile3rd_sorted` and `quantile` return `Error::UnsortedData` for a non-ascending input. The source states the precondition as `@pre` and checks it only through `OPENMS_PRECONDITION`, compiled out of a release build. A NaN also fails this check, including in a one-element range, which has no adjacent pair to disagree and which `std::is_sorted` accepts. |
| A NaN is refused wherever ordering it would decide the answer | See "NaN policy" below. The four functions above report `Error::UnsortedData`, because a NaN makes the caller's sortedness claim false; the six that sort or stage a buffer themselves — `median`, `quantile1st`, `quantile3rd`, `mad`, `compute_rank`, `rank_correlation_coefficient` — report `Error::InvalidValue`. `SummaryStatistics::new` reports it too, except for the two sample shapes in which the permutation cannot be observed. |
| `variance`, `sd` and `covariance` refuse `n < 2` | The `n - 1` divisor is zero there and the source returns NaN. A NaN variance is exactly what this layer would propagate into everything built on it. `SummaryStatistics` keeps the source's `0.0` for `n <= 1`, because the source's own comment fixes that value. |
| A zero correlation denominator returns NaN explicitly instead of dividing | For both Matthews and Pearson a zero denominator forces a zero numerator (proved below), so `0 / 0 = NaN` is the source's value in every reachable case. The one divergence is a denominator that *underflows* to zero from non-zero deviations, where the source yields an infinity and the port yields NaN. |
| `adaptive_quantile` rejects non-finite `k`, `r_sparse`, `r_dense` | The source accepts a NaN threshold and lets it decide the blend weight through comparisons that are all false. |
| Sorting uses `f64::total_cmp` | Every sorting entry point has already refused a NaN, so the total order and `std::sort`'s `<` agree on everything that reaches the sort; `total_cmp` is kept because it also orders `-0.0` before `0.0` deterministically. |
| `compute_rank` on an empty slice is a no-op | The source computes `w.size() - 1` in unsigned arithmetic, which wraps to `SIZE_MAX`. |
| `MAX_ITEMS` / `MAX_BYTES` preflight | Every function that stages an owned buffer (`mad`, `tukey_upper_fence`, `winsorized_quantile`, `adaptive_quantile`, `compute_rank`, `rank_correlation_coefficient`) checks the ceiling before allocating, so a refusal leaves the input unchanged. The source has no ceiling. Each preflight is stated per buffer, not per call: `rank_correlation_coefficient` checks one `f64` buffer and then copies two, and `compute_rank` runs its own, wider preflight for the `(usize, f64)` pairs it stages, so the widest single buffer is what `MAX_BYTES` actually bounds. |
| Lengths are compared up front | `covariance` checks its second range only by re-testing the two *begin* iterators inside the loop — a test whose value never changes — and by comparing the second iterator to its end afterwards, so a short second range is dereferenced out of bounds before the mismatch is noticed. |
| `matthews_correlation_coefficient` compares each range with itself | The source's emptiness check is `checkIteratorsNotNULL(begin_a, end_b)`: iterators into two different containers, which is undefined behaviour and is not the check it intends. |
| Serial only | The header carries no `#pragma omp`, so there is no OpenMP gap to record for this file. |

### NaN policy

The source sorts with `std::sort`. Under `operator<` a NaN is incomparable with
every value, itself included, so `std::sort` is free to return any permutation
of the elements it cannot tell apart; and once the range holds two or more
distinct numbers as well, transitivity of incomparability fails (`1 ~ NaN` and
`NaN ~ 3` while `1 < 3`) and the strict-weak-ordering precondition is violated
outright. Either way a C++ call that sorts a NaN-bearing range has no single
answer to reproduce — the permutation, and therefore the statistic, is whatever
the library's introsort happens to do. The port draws the line at *ordering*:

| Group | Functions | NaN input |
| --- | --- | --- |
| Sorts or stages a buffer it then sorts | `median`, `quantile1st`, `quantile3rd`, `mad`, `compute_rank`, `rank_correlation_coefficient` | `Error::InvalidValue`, raised before the sort, so the caller's range is left in its original order |
| Sorts, but answers the shapes whose permutation is unobservable | `SummaryStatistics::new` | `Error::InvalidValue` as above, **except** for a one-value sample and an all-NaN sample; see below |
| Requires the caller to have sorted | `median_sorted`, `quantile1st_sorted`, `quantile3rd_sorted`, `quantile` | `Error::UnsortedData`, the same error any other order violation gets |
| Drops non-finite values first | `tukey_upper_fence`, `tail_fraction_above`, `winsorized_quantile`, `adaptive_quantile` | filtered out, exactly as the source's `std::isfinite` filter does; never reaches an ordering |
| Neither sorts nor buffers | `sum`, `mean`, `variance`, `sd`, `covariance`, `mean_square_error`, `root_mean_square_error`, `absdev`, `mean_absolute_deviation`, `pearson_correlation_coefficient` | propagates to a NaN result, as in the source |
| Classifies by comparison | `classification_rate`, `matthews_correlation_coefficient` | every comparison against a NaN is false, so the pair counts as agreeing (`classification_rate`) or increments none of the four confusion counts (`matthews`). Well defined in C++ and identical there, so it is reproduced rather than refused, and documented at both items |

`mad` additionally refuses a NaN `median_of_numbers`, and refuses the staged
differences when they contain a NaN neither input had — `inf - inf` is the only
way that happens.

`SummaryStatistics::new` is the one exception, and it is a narrow one. The test
it applies is whether the set of outputs the source may produce has exactly one
member. Two sample shapes pass it, and in both cases that is a proof rather than
an observation that the permutation happened not to matter:

- **one value.** A one-element range has exactly one permutation, so there is
  nothing for `std::sort` to choose. Every positional field is that value, and
  `variance` is the `0.0` the source substitutes for `n <= 1`;
- **every value a NaN.** `std::sort` may permute freely, but every permutation
  of an all-NaN range produces the same eight fields, because every field is
  read from, or computed out of, values that are all NaN.

`SummaryStatistics::of_nan_sample` computes both without sorting, since there is
nothing to order, and a NaN next to a number is still refused. Both shapes are
reached from real input and both are pinned against the Release C++ build:
`FileInfo`'s consensusXML `-s` blocks divide, and a pair of sub-features of
intensity `-0.0` and `0.0` under one centroid contributes `(-inf) + (+inf)` to
the per-consensus-feature sample (`../oracle/a7-fileinfo`, cases `c_nan_one_s`
and `c_nan_two_s`).

The refusal that remains is a **deferral, not a decision-D1 refusal**, and is
raised for the lead. Nothing is out of bounds and the Release build's values are
stable per input, so D1 would have the port reproduce them; it does not, because
libstdc++ compares every pair involving a NaN false and therefore moves nothing,
which makes the `minimum`, quartile and `maximum` lines positional reads of a
range whose elements `std::sort` was free to leave in any order. Measured: the oracle's `c_nan_then_finite_s`
and `c_finite_then_nan_s` hold the same two consensus features in opposite file
order and disagree on exactly those four lines. Reproducing them means porting
libstdc++'s `std::sort` permutation into `sort_ascending`, which every
`SummaryStatistics` caller consumes — its own wave. See CPP-347 and section 5.2
of `docs/FILE_INFO_A7_SUPPORT.md`.

An infinity is refused nowhere. Both `std::sort` and `f64::total_cmp` order it
consistently, so the source's answer is well defined and is the port's answer.

This policy is what makes the resource-ceiling and sortedness claims above
total rather than partial. Before it, five entry points accepted a NaN, sorted
it to one end by `total_cmp` and returned a plausible finite number:
`median(&mut [1.0, NaN, 3.0])` answered `3.0` where the median of the real
values is `2.0`, and `quantile1st(&mut [1.0, NaN, 3.0, 4.0, 5.0])` answered
`2.0`. `tests/statistic_functions.rs` pins every group in the table.

### Why a zero denominator implies a zero numerator

For Matthews, the denominator is
`sqrt((tp+fp)(tp+fn)(tn+fp)(tn+fn))` over non-negative counts. A zero factor
forces two counts to zero, and each of the four cases kills one term of
`tp*tn - fp*fn` and zeroes the other: `tp = fp = 0` gives `0*tn - 0*fn`;
`tp = fn = 0` gives `0*tn - fp*0`; `tn = fp = 0` gives `tp*0 - 0*fn`;
`tn = fn = 0` gives `tp*0 - fp*0`. For Pearson, a zero `denominator_a` means
every `temp_a` is zero, which makes every term of the numerator zero. So
returning NaN is the source's value, not a substitution — this is an
independently derived equivalence, not a transcribed literal.

## Checked boundaries and evidence

Resource boundaries: `MAX_ITEMS = 50_000_000` values per staged buffer and
`MAX_BYTES = 512 MiB` of staging, both checked in a preflight before any
allocation. `usize` arithmetic in the preflight is `checked_mul`.

Numeric boundaries: empty and single-element ranges; mismatched lengths;
`q` outside `[0, 1]` and NaN `q`; a NaN input at every entry point that orders
values, in both of the groups of the NaN-policy table, and the one-element
`[NaN]` range that has no adjacent pair; infinities, which are ordered rather
than refused, including the `inf - inf` NaN that `mad` can manufacture from two
non-NaN inputs; non-finite values, which `tukey_upper_fence`,
`tail_fraction_above`, `winsorized_quantile` and `adaptive_quantile` drop as the
source does; constant ranges in all three correlation coefficients.

Evidence, per `docs/DIFFERENTIAL_VALIDATION.md`:

- **Tier 3 (source review, transcribed class-test literals)** for all 20
  `START_SECTION`s of `StatisticFunctions_test.cpp`, in
  `tests/statistic_functions.rs`. The mapping from section to test function is
  the `// L<line>` comment above each test.
- **Tier 4 (independently derived)** where the expectation follows from a closed
  form or an invariant rather than the class test's printed digits. These are
  marked "Derived" in the test: the exact `6/5` and `6/7` of
  `MeanAbsoluteDeviation` against the class test's truncated `1.2` and
  `0.857142`; the sample variance `2.0` of `SummaryStatistics` for `{1, 3}`;
  the `±1` Spearman correlation of a strictly increasing sequence with itself
  and its reversal, which follows from the symmetry of the ranks about `mu`;
  the tie rank `0.5 * (i + z + 1) = 5.5` of `computeRank`; the Tukey fence
  `14.5` and tail fraction `1/201` of adaptive-quantile case A, and both its
  interpolated quantiles `0.6*9 + 0.4*1000` and `0.6*9 + 0.4*14.5`; the
  Matthews values `1`, `-1` and `2/sqrt(12)`; and the zero-denominator
  equivalence above.
- **Tier 4 (Rust-only)** for the resource ceilings, the sortedness and NaN
  refusals, and the `n < 2` variance refusal. The NaN refusals are three tests —
  `the_sorting_entry_points_refuse_a_nan`,
  `the_buffering_entry_points_refuse_a_nan` and
  `the_sorted_entry_points_reject_a_lone_nan` — which together cover every
  function in the first two rows of the NaN-policy table, assert that a refused
  call leaves the caller's range unmodified, and assert that an infinity is
  still ordered.

No tier 1 or tier 2 evidence exists for this group: no retained C++ output
corresponds to these free functions, and no oracle driver was built.

A separate reading of one class-test line: `Math::MAD(x2, x2 + 6, true)` at
`StatisticFunctions_test.cpp:88` passes the literal `true`, which converts to
the median `1.0`, while the line's own comment states the real median `1.5`.
Both give `1.5`, so the assertion passes either way; the port's test asserts
both and says why.

## Candidate C++ issues

Reported to the integrating agent rather than written here, since
`OpenMS_CPP_ISSUES.md` is owned by that agent:

1. `matthewsCorrelationCoefficient` (line 734) calls
   `checkIteratorsNotNULL(begin_a, end_b)`, mixing iterators from two different
   containers. Comparing them is undefined behaviour, and the intended
   emptiness check is not performed.
2. `covariance` (line 618) calls `checkIteratorsAreValid(begin_b, end_b,
   begin_a, end_a)` with the *begin* iterators inside the loop, so the test is
   loop-invariant and a shorter second range is dereferenced past its end before
   the trailing `checkIteratorsEqual` notices.
3. `computeRank` (line 821) computes `Size n = (w.size() - 1)` in unsigned
   arithmetic; an empty vector wraps `n` to `SIZE_MAX`.
4. `quantile1st`/`quantile3rd` (lines 254 and 289) exclude one element from each
   half for an even count, where no middle element exists to exclude.
5. `variance` (line 555) divides by `std::distance(begin, end) - 1` without
   checking, returning NaN for a single value; `SummaryStatistics` works around
   it locally instead of fixing it.
