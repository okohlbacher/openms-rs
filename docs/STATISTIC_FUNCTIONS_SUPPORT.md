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
| A NaN is **reproduced** wherever the source's own `std::sort` decides the answer | Lead decision D16; see "NaN policy" below. `median`, `quantile1st`, `quantile3rd`, `mad` and `SummaryStatistics::new` sort with the Release build's own permutation and read their order statistics positionally out of it. The four `_sorted` functions above still report `Error::UnsortedData`, because a NaN makes the *caller's* sortedness claim false; `compute_rank` and `rank_correlation_coefficient` still report `Error::InvalidValue`, for the two reasons under "Where a NaN lands". |
| `variance`, `sd` and `covariance` refuse `n < 2` | The `n - 1` divisor is zero there and the source returns NaN. A NaN variance is exactly what this layer would propagate into everything built on it. `SummaryStatistics` keeps the source's `0.0` for `n <= 1`, because the source's own comment fixes that value. |
| A zero correlation denominator returns NaN explicitly instead of dividing | For both Matthews and Pearson a zero denominator forces a zero numerator (proved below), so `0 / 0 = NaN` is the source's value in every reachable case. The one divergence is a denominator that *underflows* to zero from non-zero deviations, where the source yields an infinity and the port yields NaN. |
| `adaptive_quantile` rejects non-finite `k`, `r_sparse`, `r_dense` | The source accepts a NaN threshold and lets it decide the blend weight through comparisons that are all false. |
| Sorting is the Release build's own `std::sort` | `sort_ascending` is `source_sort_by(&mut values, \|a, b\| a < b)`, the libstdc++ introsort of `crate::math::source_sort` under the default `operator<` — the call at `:140`, `:244`, `:281` and `:948`, which `MAD` reaches through its own `median` call at `:189`. It replaced an `f64::total_cmp` sort, which differed from `operator<` exactly on a NaN and on a signed zero, and both differences are now closed. `compute_rank` is the one sort in this file still on `total_cmp`; "Where a NaN lands" says why. |
| A generated NaN carries the Release build's bits | IEEE 754 does not fix which NaN an operation produces from non-NaN operands. Every function whose arithmetic can generate one is built on `crate::math::x86_64`, so `inf - inf` is `0xfff8000000000000` on every host and not the arm64 default `0x7ff8000000000000`. No finite value changes; see "NaN bit patterns" below. |
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
outright.

**Unspecified is not unknowable.** The Release build runs one particular
algorithm — the conda-forge GCC 14.4.0 libstdc++ introsort — and it runs it
deterministically, so the permutation it leaves is a measurable fact about that
build. Lead decision **D16** (shared-math wave, 2026-09-19; `docs/VALIDATION.md`)
puts reproducing it in scope, and
`src/math/source_sort.rs` is the comparison-by-comparison port of it, validated
tier 1 against two oracle drivers over 2,272 inputs. This module's private
`sort_ascending` is therefore `std::sort(begin, end)` itself:
`source_sort_by(&mut values, |a, b| a < b)`.

| Group | Functions | NaN input |
| --- | --- | --- |
| Sorts, or stages a buffer it then sorts | `median`, `quantile1st`, `quantile3rd`, `mad`, `SummaryStatistics::new` | **reproduced**: the range is sorted into the permutation the Release build's `std::sort` leaves, and the order statistics are read positionally out of it, exactly as `StatisticFunctions.h:952-956` reads them |
| Sorts with a *different* `std::sort` | `compute_rank`, `rank_correlation_coefficient` | `Error::InvalidValue`, raised before the sort, so the caller's range is left in its original order; see 5.2 below for why this one was not moved |
| Requires the caller to have sorted | `median_sorted`, `quantile1st_sorted`, `quantile3rd_sorted`, `quantile` | `Error::UnsortedData`, the same error any other order violation gets. These do not sort; they verify a claim, and a NaN makes that claim false |
| Drops non-finite values first | `tukey_upper_fence`, `tail_fraction_above`, `winsorized_quantile`, `adaptive_quantile` | filtered out, exactly as the source's `std::isfinite` filter does; never reaches an ordering |
| Neither sorts nor buffers | `sum`, `mean`, `variance`, `sd`, `covariance`, `mean_square_error`, `root_mean_square_error`, `absdev`, `mean_absolute_deviation`, `pearson_correlation_coefficient` | propagates to a NaN result, as in the source |
| Classifies by comparison | `classification_rate`, `matthews_correlation_coefficient` | every comparison against a NaN is false, so the pair counts as agreeing (`classification_rate`) or increments none of the four confusion counts (`matthews`). Well defined in C++ and identical there, so it is reproduced rather than refused, and documented at both items |

**The one refusal `sort_ascending` still has** is decision D1's, and it is
unreachable in practice: `Error::InvalidValue` at the step where the introsort's
unbounded partition or final-insertion loop would read *outside* the vector,
which is undefined behaviour with no reproducible result. Only a comparison that
is not asymmetric can provoke it, and `<` on `f64` keys is asymmetric with NaN
keys too. `src/math/source_sort.rs`'s module documentation carries the proof,
and none of the 2,272 oracle inputs reached the guard.

`SummaryStatistics::new` no longer has a NaN special case at all. It used to
answer two shapes — a one-value sample and an all-NaN sample — whose set of
possible outputs has exactly one member, and refuse everything else; both are
now just ordinary sorts, and `SummaryStatistics::of_nan_sample` is gone.

An infinity is refused nowhere: `operator<` orders it consistently, so the
source's answer is well defined and is the port's answer.

#### Where a NaN lands, and the two refusals that remain

`std::sort` is `__introsort_loop` followed by `__final_insertion_sort`
(`stl_algo.h:1899-1910`), and `__introsort_loop` runs only
`while (__last - __first > int(_S_threshold))` with `_S_threshold` enumerated as
16 (`stl_algo.h:1806`, `stl_algo.h:1880`). `__final_insertion_sort`
(`:1812-1823`) is one `__insertion_sort` pass (`:1770-1788`) below the
threshold, and above it one over the first 16 elements plus an
`__unguarded_insertion_sort` over the rest. So at **16 elements or fewer** the
whole sort is a single `__insertion_sort` pass, and at **17 or more** the
partitioning is free to move a NaN, which for the `{NaN, 2..n}` family it does.
Neither bound says the NaN stands still below the threshold: a block move
carries it without any comparison involving it ever being true.

Measured with the reference compiler, three byte-stable runs: `{NaN, 2}` and
`{2, NaN}` come back unchanged; `{NaN, 2..16}` keeps the NaN at index 0, while
`{NaN, 2..17}` moves it to index 8 and `{NaN, 2..20}` to index 10 — which is why
that sample prints `minimum: 2` and `median: -nan`. `{3, NaN, 2}` sorts to
`{2, 3, NaN}` at three elements, because `__insertion_sort` relocates a whole
block when a later element belongs before `*__first`.

**The port reproduces all five**, and pins the whole permutation rather than the
index alone: `the_introsort_threshold_decides_where_a_nan_lands` and
`the_sorting_entry_points_reproduce_the_release_builds_permutation` in
`tests/statistic_functions.rs`. The 17-element case is worth spelling out
because it is the one the algorithm decides rather than leaves alone:
`__unguarded_partition_pivot` takes the median of `*(first+1)`, `*mid` and
`*(last-1)` — 2, 9 and 17 — and `__move_median_to_first` swaps the middle one to
the front, which is what carries the NaN to index 8; the partition then stops on
both sides at that NaN, and `__final_insertion_sort` bubbles the displaced 9
back into place without touching it again.

**The open question this section used to carry is answered.** It asked whether
reproducing an unspecified `std::sort` permutation is in scope at all. Decision
D16 says yes, on four grounds set out in `docs/VALIDATION.md`. The
five frozen oracle cases the question rested on — `c_nan_one_s`, `c_nan_two_s`,
`c_nan_then_finite_s`, `c_finite_then_nan_s` and `c_zero_swapped_s` of
`../oracle/a7-fileinfo` — are now reproduced and compared line for line by
`tests/file_info_a7.rs`; each differs from its retained Release report only in
native difference 5's `nan` / `-nan` spelling, which belongs to the FileInfo
text layer and not to this module. `CPP-347` is rewritten accordingly: the port
reproduces the defect, and the defect stands.

**Two refusals remain, deliberately.**

1. `compute_rank` and `rank_correlation_coefficient`. Their `std::sort`
   (`StatisticFunctions.h:829-830`) is a lambda comparing `std::pair::second`,
   not the default `operator<`, so it is a different call from the five
   `sort_ascending` covers; and a NaN additionally defeats their **relative tie
   test**, whose two comparisons are both false against a NaN and which would
   therefore merge every block the NaN touches. That is a second, independent
   behaviour that no oracle row measures, so the refusal stands rather than
   being guessed at. For a NaN-free range nothing is lost: `operator<` and
   `total_cmp` differ only on `±0.0`, which the tie test makes one block either
   way, and the ranks are written back by origin index, so the two sorts give
   the same result.
2. The `_sorted` entry points. They do not sort; they verify a *caller's* claim
   to have sorted, and there is no way to know which permutation a caller who
   asserts "already sorted" about a NaN-bearing range meant.
   `SummaryStatistics::new` does its own sorting and so reads through the
   private `_of_sorted` helpers, which do not re-check.

#### The cost of the faithful sort, and the fast path that removes it

`sort_ascending` is no longer a library sort. It builds a permutation of
`0..n` with the libstdc++ introsort reproduced in Rust, calling a closure for
every comparison, where it used to call `slice::sort_by`. That is **16.5x**
slower than a library sort of the same ten million values, and it cost the
public entry point **19.1x** wall clock and **2.9x** peak memory before this
repair. The samples it is handed are not small:
`src/format/file_info/peaks.rs:706` and `:714` give `summarize` every MS1 peak
intensity in the file, bounded only by
`FileInfo::MAX_STATISTICS_VALUES = 1 << 27`. On a routine LC-MS run that was a
real regression on `FileInfo -s`, not a theoretical one.

Lead decision **D17** (`docs/VALIDATION.md`) closes it without giving anything
up. Where the sample holds **no NaN** and **not both spellings of zero**,
`sort_ascending` sorts in place with `f64::total_cmp` and allocates nothing;
otherwise it runs the libstdc++ permutation, unchanged. That is not a
compromise, because in exactly that case the permutation cannot be observed:
without a NaN, `operator<` is a strict weak ordering whose equivalence relation
is numeric equality, and two numerically equal non-NaN doubles are
bit-identical — with `-0.0 == +0.0` the one exception in the whole format. Every
equivalence class is then a set of identical bytes, the sorted sequence is a
function of the multiset alone, and any correct sort writes what the Release
build writes. The guard is one O(n) pass testing those two things and nothing
else.

The argument is not what the port rests on.
`both_paths_agree_bit_for_bit_wherever_the_fast_one_is_taken` runs **both**
paths over the same adversarial samples — both zeros, both infinities,
subnormals, `DBL_MAX`, signalling and negative NaNs, heavy duplication,
ascending, descending, organ-pipe and sawtooth shapes, and raw random bit
patterns, at 21 lengths spanning libstdc++'s 16-element `_S_threshold` and its
heapsort fallback — and compares the results bit for bit, NaN payloads included.
Deleting either half of the guard makes it fail.
`the_public_entry_points_agree_with_the_release_builds_permutation` shows
`median`'s public surface reaching both paths, and
`the_guard_is_exactly_a_nan_or_both_zero_spellings` pins the boundary.

**The measurement is in [BENCHMARKS](BENCHMARKS.md) §8**, with the host, the
load, the command and the committed harness
(`tools/bench_sort_ascending.sh`, driving `sort_ascending_benchmark` in
`tests/statistic_functions.rs`); the numbers are not repeated here. In one
line: at ten million values the public entry point went from **2.23 s and
464 MiB** to **118 ms and 159 MiB**, and over every mzML fixture in
`tests/data` **809 of 809** statistics samples take the fast path.

Who calls this at all is still narrow. The consumers of the sorting entry points
are `FileInfo`'s `summarize` (`src/format/file_info/report.rs`) and `fasta.rs`'s
sequence-length summary. `mass_trace_detection.rs` has a `median` of its own and
does not reach this one, and the picked feature finder already called
`crate::math::source_sort` directly and is unchanged by all of this.

### Signed zeros

`-0.0` and `0.0` are the other pair `operator<` calls *equivalent*: both
`-0.0 < 0.0` and `0.0 < -0.0` are false. Unlike a NaN they do not break
`std::sort`'s precondition, so the call is well formed — but every permutation
is a conforming result, and libstdc++ leaves a short range as it found it.

`sort_ascending` used to order them by `f64::total_cmp`, which puts `-0.0` first
regardless of input order, so a caller reading order statistics positionally out
of such a sample disagreed with the source in the sign of a printed zero. One
did: `FileInfo`'s consensusXML `-s` `Intensity ratios` block, recorded as
**native difference 6** of `docs/FILE_INFO_A7_SUPPORT.md` and measured in both
orders (oracle cases `c_nan_one_s` and `c_zero_swapped_s`, the same consensus
feature with its two sub-feature intensities exchanged).

Running libstdc++'s own permutation closes that difference. The port keeps the
input order the Release build keeps, reproduces **both** members of the measured
pair rather than collapsing them onto one, and
`consensus_a_signed_zero_sample_keeps_the_release_builds_order` asserts equality
with both retained reports.
`a_signed_zero_keeps_the_order_the_release_build_keeps` pins the same thing at
the `SummaryStatistics` level. Native difference 6 is closed.

### NaN bit patterns

IEEE 754 fixes every finite result of `+`, `-`, `*`, `/` and `sqrt`, but not
which NaN bit pattern an operation *generates* out of non-NaN operands.
`inf - inf`, `0 * inf`, `0 / 0`, `inf / inf` and the square root of a negative
number all yield "some" NaN, and the answer is the instruction set's: the Linux
x86_64 Release build's SSE2 produces the "real indefinite" QNaN
`0xfff8000000000000`, which glibc spells `-nan`, while an arm64 host produces
the positive default NaN `0x7ff8000000000000`.

That difference is observable. `FileInfo`'s consensusXML `-s` block prints the
variance of `{1, +inf}`, which is `(1 - inf)^2 + (inf - inf)^2` over one — a
generated NaN — and native difference 5 of `docs/FILE_INFO_SUPPORT.md` is the
spelling of exactly that value. The *value* had to become host-independent
before the spelling could be fixed.

Every function here whose own arithmetic can generate a NaN is therefore built
on `crate::math::x86_64`'s `add`, `sub`, `mul`, `div`, `sqrt` and `abs`:

The table below is meant to be read as exhaustive over the module's public
surface, so what it leaves out is named here rather than left to be noticed.
`check_not_empty`, `check_exhausted` and `check_ranges_end_together` do no
arithmetic. `median`, `median_sorted`, `quantile1st`, `quantile1st_sorted`,
`quantile3rd`, `quantile3rd_sorted` and `SummaryStatistics::new` do none of
their own either: they order a range and read out of it, and the only
arithmetic under them is `median_of_sorted`'s even-size average and, for
`SummaryStatistics`, `mean` and `variance` — all three of which are listed.

| | Functions | Why |
| --- | --- | --- |
| **Rebuilt** | `sum`, `mean`, `variance`, `variance_with_mean`, `sd`, `sd_with_mean`, `covariance`, `mean_square_error`, `root_mean_square_error`, `mean_absolute_deviation`, `absdev`, `absdev_with_mean`, `mad`'s `fabs`, `median_of_sorted`'s even-size average, `quantile`'s linear blend | the source's own `double` arithmetic can produce a NaN from operands that are not NaN. `absdev` and `absdev_with_mean` inherit it: both are `mean_absolute_deviation` behind an emptiness check, and `absdev`'s centre is the rebuilt `mean` |
| **Not rebuilt** | `classification_rate`, `matthews_correlation_coefficient` (their counting), `compute_rank`, `rank_correlation_coefficient`, `tukey_upper_fence`, `tail_fraction_above`, `winsorized_quantile`, `adaptive_quantile` | no NaN can be generated: the first two count with comparisons, `compute_rank` averages small integers and `rank_correlation_coefficient` hands those ranks to `pearson_correlation_coefficient`, and the Tukey family drops every non-finite value before it computes anything, as the source's `std::isfinite` filter does |
| **Not rebuilt, a known gap** | `pearson_correlation_coefficient`, `matthews_correlation_coefficient` | both substitute an explicit `f64::NAN` for a division the source actually performs — a divergence that predates this work and is documented at each item — so making only their *other* operations bit-faithful would leave that substituted NaN as the single host-shaped value in the result. The crate's bit-faithful Pearson is `analysis::feature_finder_picked::scoring::source_pearson`; converging the two needs an oracle row of its own |

**Two operand orders in this module are not measured, and are named rather than
glossed.** SSE2's two-operand instructions return the **destination** operand
quieted when both operands are NaN, so where the source writes a *commutative*
operation on two *distinct* temporaries, which one GCC leaves in the destination
register decides the answer's payload — and that is a register-allocation fact,
not something the source fixes.

| item | the source expression | reachable when | what the port answers |
| --- | --- | --- | --- |
| `covariance` | `(*iter_a - mean_a) * (*iter_b - mean_b)` (`:619`) | `a[i]` a payload NaN (so `mean_a` is one too) while `b[i] == mean_b == +inf` | the `a`-derived NaN's payload |
| `median_of_sorted`, even size | `(*(it + n/2 - 1) + *(it + n/2)) / 2.0` | `median(&mut [nan_a, nan_b])` with two distinct payloads, which since D16 is summarised rather than refused | `nan_a`'s payload |

No oracle row measures either, and none is invented. What would settle them:
one call each on exactly those inputs, run on the Linux x86_64 Release build
with the result read back as bits; or, without running anything, disassembling
the two loops in the reference `libOpenMS.so` and reading which operand the
emitted `mulsd` and `addsd` write to. Everywhere else in the module the order is
forced: `subsd` and `divsd` are not commutative, `mul(diff, diff)` has one
operand twice, an accumulating `sum += x` makes the accumulator the destination,
and `quantile`'s blend cannot have both operands NaN because it refuses a NaN
range.

The helpers return the IEEE result whenever it is not a NaN, so **no finite
value changes**. That is asserted directly rather than assumed:
`the_x86_64_helpers_change_no_finite_result` compares bit for bit against the
plain-Rust arithmetic the module used before, over a battery reaching
subnormals, both zeros, `DBL_MAX` and ranges whose squared deviations overflow
to an infinity, across every function and every equally long pair. No frozen
expectation in `tests/statistic_functions.rs` moved.
`a_generated_nan_carries_the_release_builds_bits` then pins fifteen bit
patterns derived from the SSE2 rules of Intel SDM vol. 1 rather than from this
crate's output — including the *positive* NaN `andpd` leaves behind in
`mean_absolute_deviation`, because the absolute-value mask clears the sign bit
of a NaN like any other value, and the quieted payload a signalling NaN input
keeps.

One shape worth naming because it reads the other way round:
`variance_with_mean(&[+inf, -inf], 0.0)` generates no NaN at all. Both squared
deviations are `+inf` and they add, so the result is `+inf`. The indefinite NaN
appears where a *subtraction* cancels two infinities — `variance(&[1.0, +inf])`,
which is the `FileInfo` path, and `variance_with_mean(&[+inf, -inf], +inf)`.

The *spelling* is a separate layer and is not decided here:
`format::file_info::text_format` writes every NaN as `nan` where glibc writes
`-nan` for one whose sign bit is set. That is native difference 5, and the
shared-math wave deliberately did not touch it; it also needs A2's oracle row
re-captured against the Linux Release build.

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
values, in every group of the NaN-policy table, and the one-element `[NaN]`
range that has no adjacent pair; the `_S_threshold` boundary at 16, 17 and 20
elements; both signed zeros, in both input orders; infinities, which are ordered
rather than refused, including the `inf - inf` NaN that `mad` can manufacture
from two non-NaN inputs; the bit pattern of every NaN this module can generate;
non-finite values, which `tukey_upper_fence`,
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
- **Tier 4 (Rust-only)** for the resource ceilings, the sortedness refusals,
  the two NaN refusals 5.2 keeps, and the `n < 2` variance refusal:
  `the_sorted_entry_points_reject_a_lone_nan` and
  `mad_reproduces_a_nan_and_the_ranking_still_refuses_one`, which also assert
  that a refused call leaves the caller's range unmodified.
- **Tier 4 (independently derived, from the algorithm rather than a run)** for
  the `std::sort` permutation itself, which is the shared-math wave's own
  evidence. `the_sorting_entry_points_reproduce_the_release_builds_permutation`,
  `the_introsort_threshold_decides_where_a_nan_lands` and
  `a_signed_zero_keeps_the_order_the_release_build_keeps` derive every expected
  array from the libstdc++ algorithm — `__insertion_sort`'s two branches,
  `__move_median_to_first`'s swap and `_S_threshold` — and not from this crate's
  output. The three NaN positions they pin (index 0, 8 and 10 for
  `{NaN, 2..16}`, `{NaN, 2..17}` and `{NaN, 2..20}`) are the ones 5.2 measured
  against the reference compiler, so the derivation and the measurement agree.
- **Tier 4 (independently derived)** for the NaN bit patterns:
  `a_generated_nan_carries_the_release_builds_bits` applies the SSE2 rules of
  Intel SDM vol. 1 to the instruction the Release build emits, and
  `the_x86_64_helpers_change_no_finite_result` asserts the invariant those
  helpers rest on — that no non-NaN result moves — bit for bit against the
  plain-Rust arithmetic.

No tier 1 or tier 2 evidence exists for this group *as free functions*: no
retained C++ output corresponds to them and no oracle driver was built for this
header. The sorting behaviour is the exception and has tier 1 evidence at one
remove, through two consumers: `../oracle/ffap-complete-fix1` and
`../oracle/ffap-instr-completion` validate `crate::math::source_sort` itself
over 2,272 inputs, and `../oracle/a7-fileinfo`'s five frozen `-s` reports
(`c_nan_one_s`, `c_nan_two_s`, `c_nan_then_finite_s`, `c_finite_then_nan_s`,
`c_zero_swapped_s`) are what `tests/file_info_a7.rs` compares these functions'
output against, line for line.

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
