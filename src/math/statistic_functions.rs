// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Means, medians, quantiles, deviations and correlation coefficients.
//!
//! Port of `src/openms/include/OpenMS/MATH/StatisticFunctions.h` (the
//! accompanying `StatisticFunctions.cpp` is an empty translation unit; the
//! header is header-only). See `docs/STATISTIC_FUNCTIONS_SUPPORT.md`.
//!
//! The source takes iterator pairs; this port takes slices. An operation the
//! source expresses with a `bool sorted` default argument becomes two
//! functions, because the `false` case mutates the caller's range by sorting
//! it: [`crate::math::statistic_functions::median`] takes `&mut [f64]` and
//! sorts, [`crate::math::statistic_functions::median_sorted`] takes `&[f64]`
//! and does not. The same split applies to `quantile1st` and `quantile3rd`.
//!
//! Accumulation order is the source's throughout. Where the source's divisor
//! can be zero the port either refuses or returns the source's value
//! explicitly; each such place is documented at the item.
//!
//! # NaN bit patterns
//!
//! IEEE 754 fixes every finite result of `+`, `-`, `*`, `/` and `sqrt`, but not
//! which NaN bit pattern an operation *generates* out of non-NaN operands.
//! `inf - inf`, `0 * inf`, `0 / 0`, `inf / inf` and the square root of a
//! negative number all yield "some" NaN, and the answer is the instruction
//! set's: the Linux x86_64 Release build's SSE2 produces
//! `0xfff8000000000000`, which glibc spells `-nan`, while an arm64 host
//! produces the positive default NaN `0x7ff8000000000000`. That difference is
//! observable: `FileInfo`'s consensusXML `-s` `Intensity ratios` block prints
//! the variance of `{1, +inf}`, which is `(1 - inf)^2 + (inf - inf)^2` divided
//! by one — a generated NaN.
//!
//! Every function here whose own arithmetic can generate a NaN is therefore
//! built on `crate::math::x86_64`'s `add`, `sub`, `mul`, `div`, `sqrt` and
//! `abs` rather than on plain Rust operators. Those helpers return the IEEE
//! result whenever it is not a NaN, so **no finite value changes** — that is
//! asserted directly in `tests/statistic_functions.rs`
//! (`the_x86_64_helpers_change_no_finite_result`), and every frozen expectation
//! in that file is unchanged.
//!
//! Which functions, and why:
//!
//! - **Rebuilt**, because the source's own `double` arithmetic can produce a
//!   NaN from operands that are not NaN: [`sum`](crate::math::statistic_functions::sum) (`std::accumulate`'s `addsd`),
//!   [`mean`](crate::math::statistic_functions::mean) (its `divsd`), [`variance`](crate::math::statistic_functions::variance), [`variance_with_mean`](crate::math::statistic_functions::variance_with_mean), [`sd`](crate::math::statistic_functions::sd),
//!   [`sd_with_mean`](crate::math::statistic_functions::sd_with_mean), [`covariance`](crate::math::statistic_functions::covariance), [`mean_square_error`](crate::math::statistic_functions::mean_square_error),
//!   [`root_mean_square_error`](crate::math::statistic_functions::root_mean_square_error), [`mean_absolute_deviation`](crate::math::statistic_functions::mean_absolute_deviation) and the `fabs` of
//!   [`mad`](crate::math::statistic_functions::mad)'s staged differences, plus the two places an order statistic is
//!   *interpolated* rather than read: the even-size average in
//!   `median_of_sorted` (and so [`median`](crate::math::statistic_functions::median), [`median_sorted`](crate::math::statistic_functions::median_sorted),
//!   [`quantile1st`](crate::math::statistic_functions::quantile1st), [`quantile3rd`](crate::math::statistic_functions::quantile3rd) and their `_sorted` forms, which reduce
//!   to it) and the linear blend in [`quantile`](crate::math::statistic_functions::quantile).
//! - **Not rebuilt**, because no NaN can be generated:
//!   [`classification_rate`](crate::math::statistic_functions::classification_rate) and [`matthews_correlation_coefficient`](crate::math::statistic_functions::matthews_correlation_coefficient) count
//!   with comparisons; `compute_rank` averages small integers; and
//!   [`tukey_upper_fence`](crate::math::statistic_functions::tukey_upper_fence), [`tail_fraction_above`](crate::math::statistic_functions::tail_fraction_above), [`winsorized_quantile`](crate::math::statistic_functions::winsorized_quantile)
//!   and [`adaptive_quantile`](crate::math::statistic_functions::adaptive_quantile) drop every non-finite value before they compute
//!   anything, as the source's `std::isfinite` filter does, so no infinity
//!   reaches their arithmetic at all.
//! - **Not rebuilt, and a known gap**: [`pearson_correlation_coefficient`](crate::math::statistic_functions::pearson_correlation_coefficient) and
//!   [`matthews_correlation_coefficient`](crate::math::statistic_functions::matthews_correlation_coefficient) substitute an explicit `f64::NAN` for
//!   a division the source actually performs — a divergence that predates this
//!   work and is documented at each item. Making only their *other* operations
//!   bit-faithful would leave that substituted NaN as the single host-shaped
//!   value in the result and would be more misleading than the present state.
//!   The crate's bit-faithful Pearson is
//!   `analysis::feature_finder_picked::scoring::source_pearson`, which the
//!   picked feature finder measured against the Release build; converging the
//!   two needs an oracle row of its own.
//!
//! The *spelling* of a NaN is a separate layer and is not decided here:
//! `format::file_info::text_format` writes every NaN as `nan`, where glibc
//! writes `-nan` for one whose sign bit is set. That is native difference 5 of
//! `docs/FILE_INFO_SUPPORT.md`; making the value host-independent, which is
//! what this section is about, is its precondition and not its fix.
//!
//! # NaN
//!
//! A NaN has no place in a total order. Under `operator<` it is incomparable
//! with every value, itself included, so `std::sort`'s strict-weak-ordering
//! precondition fails as soon as the range holds two distinct numbers as well
//! (`1 ~ NaN` and `NaN ~ 3` while `1 < 3`), and the standard leaves the result
//! unspecified.
//!
//! Unspecified is not unknowable. The Release build runs one particular
//! algorithm — the conda-forge GCC 14.4.0 libstdc++ introsort — and it runs it
//! deterministically, so the permutation it leaves is a measurable fact about
//! that build rather than a coin toss. Decision **D16** of
//! `docs/EARLY_TOPP_WORK_PACKAGES.md` puts reproducing it in scope, and
//! [`crate::math::source_sort`] is the comparison-by-comparison port of it,
//! validated tier 1 against two oracle drivers over 2,272 inputs carrying ties,
//! signed zeros, infinities and four NaN bit patterns.
//!
//! So this module's private `sort_ascending` **is** `std::sort(begin, end)`:
//! `source_sort_by(&mut values, |a, b| a < b)`. Every entry point that orders
//! values therefore *reproduces* a NaN-bearing sample instead of refusing it:
//!
//! - [`median`](crate::math::statistic_functions::median),
//!   [`quantile1st`](crate::math::statistic_functions::quantile1st),
//!   [`quantile3rd`](crate::math::statistic_functions::quantile3rd),
//!   [`mad`](crate::math::statistic_functions::mad) and
//!   [`SummaryStatistics::new`](crate::math::statistic_functions::SummaryStatistics::new)
//!   sort the range and read their order statistics positionally out of it, as
//!   `StatisticFunctions.h:140`, `:244`, `:281`, `:189` and `:948-956` do. The
//!   positions the NaN lands in are libstdc++'s: `{NaN, 2..16}` leaves it at
//!   index 0, `{NaN, 2..17}` moves it to index 8 and `{NaN, 2..20}` to index
//!   10, because `__introsort_loop` only runs while the range is longer than
//!   `_S_threshold == 16` (`stl_algo.h:1806`, `:1880`, `:1899-1910`, with
//!   `__insertion_sort` at `:1770-1788` and `__final_insertion_sort` at
//!   `:1812-1823`). That is why `{NaN, 2..20}` prints `minimum: 2` and
//!   `median: -nan`.
//! - The only refusal left is
//!   [`Error::InvalidValue`](crate::Error::InvalidValue) at the step where the
//!   introsort's unbounded partition or final-insertion loop would read
//!   *outside* the vector — decision D1, undefined behaviour with no
//!   reproducible result. No asymmetric comparison can reach it, and `<` on
//!   `f64` keys is asymmetric with NaN keys too, so in practice no NaN input is
//!   turned away. [`crate::math::source_sort`]'s module documentation carries
//!   the proof.
//! - [`median_sorted`](crate::math::statistic_functions::median_sorted),
//!   [`quantile1st_sorted`](crate::math::statistic_functions::quantile1st_sorted),
//!   [`quantile3rd_sorted`](crate::math::statistic_functions::quantile3rd_sorted)
//!   and [`quantile`](crate::math::statistic_functions::quantile) are the
//!   exception, and deliberately so. They do not sort; they require the
//!   *caller* to have sorted, and a NaN makes that claim false. They keep
//!   returning [`Error::UnsortedData`](crate::Error::UnsortedData), as they do
//!   for any other order violation, because there is no way to know which
//!   permutation a caller who asserts "already sorted" about a NaN-bearing
//!   range meant. `SummaryStatistics::new` does its own sorting and so reads
//!   through the private `_of_sorted` helpers, which do not re-check.
//! - [`compute_rank`](crate::math::statistic_functions::compute_rank) and
//!   [`rank_correlation_coefficient`](crate::math::statistic_functions::rank_correlation_coefficient)
//!   still refuse. Their `std::sort` (`:829-830`) is a different call — a
//!   lambda comparing `std::pair::second`, not the default `operator<` — and a
//!   NaN additionally defeats their *tie test*, whose two comparisons are both
//!   false against a NaN, so every block it touches would be merged. That is a
//!   second, independent behaviour that no oracle row measures, so this wave
//!   left the refusal standing rather than guess at it; see section 5.2 of
//!   `docs/STATISTIC_FUNCTIONS_SUPPORT.md`. For a NaN-free range the two sorts
//!   agree on the ranks anyway: `operator<` and `total_cmp` differ only on
//!   `±0.0`, which the tie test makes one block either way.
//! - [`tukey_upper_fence`](crate::math::statistic_functions::tukey_upper_fence),
//!   [`tail_fraction_above`](crate::math::statistic_functions::tail_fraction_above),
//!   [`winsorized_quantile`](crate::math::statistic_functions::winsorized_quantile)
//!   and [`adaptive_quantile`](crate::math::statistic_functions::adaptive_quantile)
//!   **drop** non-finite values before ordering anything, exactly as the
//!   source's `std::isfinite` filter does, and so never see a NaN at all.
//!
//! What this does **not** make defensible is the source's own behaviour.
//! `Math::SummaryStatistics` reads `front()`, three quantiles and `back()`
//! positionally out of a range whose order the comparison did not determine, so
//! the same multiset in a different input order gives different printed lines —
//! `../oracle/a7-fileinfo`'s `c_nan_then_finite_s` and `c_finite_then_nan_s`
//! measure exactly that. The port now reproduces the defect rather than
//! refusing it; **CPP-347** of `docs/OpenMS_CPP_ISSUES.md` keeps the defect
//! itself on the record for upstream.
//!
//! # Signed zeros
//!
//! `-0.0` and `0.0` are the other pair `operator<` calls *equivalent*: both
//! `-0.0 < 0.0` and `0.0 < -0.0` are false. Unlike a NaN they do not break
//! `std::sort`'s precondition, so the call is well formed — but every
//! permutation is a conforming result, and libstdc++ leaves a short range as it
//! found it.
//!
//! `sort_ascending` used to order them by the IEEE-754 total order, which puts
//! `-0.0` first regardless of input order, and a caller reading order
//! statistics positionally out of such a sample then disagreed with the source
//! in the sign of a printed zero. One did: `FileInfo`'s consensusXML `-s`
//! `Intensity ratios` block, native difference 6 of
//! `docs/FILE_INFO_A7_SUPPORT.md`. Running libstdc++'s own permutation closes
//! that difference — the port now keeps the input order the Release build keeps
//! — and `tests/file_info_a7.rs` asserts the two retained Release reports byte
//! for byte rather than recording a divergence.
//!
//! An infinity is *not* refused anywhere: `operator<` orders it consistently,
//! so the source's answer is well defined and is reproduced. The functions that
//! neither sort nor buffer — `sum`, `mean`,
//! `variance`, `covariance`, `mean_square_error`, `mean_absolute_deviation`,
//! `pearson_correlation_coefficient` — let a NaN propagate into the result,
//! which is what the source does and is honest about the input. The two label
//! functions are the exception and are documented at the item:
//! [`classification_rate`](crate::math::statistic_functions::classification_rate)
//! and
//! [`matthews_correlation_coefficient`](crate::math::statistic_functions::matthews_correlation_coefficient)
//! classify by comparison, and every comparison against a NaN is false, so a
//! NaN pair silently falls through — source behaviour, reproduced.

use crate::math::source_sort::source_sort_by;
use crate::math::x86_64;
use crate::{Error, Result};

/// Maximum number of values a single call may stage into an owned buffer.
///
/// The source allocates without a ceiling; a hostile or accidental input is
/// therefore only bounded by the allocator. Native-only guard.
pub const MAX_ITEMS: usize = 50_000_000;

/// Maximum owned bytes a single call may stage before it does any work.
///
/// Checked in a preflight, so a refusal leaves the input untouched.
/// Native-only guard.
pub const MAX_BYTES: usize = 512 * 1024 * 1024;

/// Default Tukey factor `k` of the source's `tukeyUpperFence` and
/// `adaptiveQuantile`.
pub const DEFAULT_TUKEY_FACTOR: f64 = 1.5;

/// Default tail density below which `adaptiveQuantile` uses the robust value.
pub const DEFAULT_R_SPARSE: f64 = 0.01;

/// Default tail density above which `adaptiveQuantile` uses the raw value.
pub const DEFAULT_R_DENSE: f64 = 0.10;

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.to_string())
}

fn empty(message: &str) -> Error {
    Error::InvalidRange(message.to_string())
}

/// Refuse an owned staging buffer of `count` elements of `width` bytes.
fn preflight(count: usize, width: usize) -> Result<()> {
    if count > MAX_ITEMS {
        return Err(empty(
            "statistics input exceeds the supported element count",
        ));
    }
    let bytes = count
        .checked_mul(width)
        .ok_or_else(|| empty("statistics input exceeds the supported staging size"))?;
    if bytes > MAX_BYTES {
        return Err(empty("statistics input exceeds the supported staging size"));
    }
    Ok(())
}

/// True when no value is NaN.
fn is_free_of_nan(values: &[f64]) -> bool {
    !values.iter().any(|value| value.is_nan())
}

/// True when every adjacent pair is non-decreasing and no value is NaN.
///
/// The NaN test is written out rather than left to the pairwise comparison: a
/// NaN does make `pair[0] <= pair[1]` false for every pair it takes part in,
/// but a one-element range has no pair at all, and `[NaN]` would otherwise pass
/// for "sorted" the way `std::is_sorted` passes it.
fn is_ascending(values: &[f64]) -> bool {
    is_free_of_nan(values) && values.windows(2).all(|pair| pair[0] <= pair[1])
}

/// Fail when any value is NaN.
///
/// Since decision D16 the functions that *sort* no longer use this: they
/// reproduce the Release build's own permutation instead. Its two remaining
/// callers are [`compute_rank`](crate::math::statistic_functions::compute_rank)
/// and
/// [`rank_correlation_coefficient`](crate::math::statistic_functions::rank_correlation_coefficient),
/// whose `std::sort` is a different call — a lambda on `std::pair::second`
/// rather than the default `operator<` — and whose relative tie test a NaN
/// defeats in a second, independent way that no oracle row measures. See the
/// module's NaN section.
fn check_no_nan(values: &[f64]) -> Result<()> {
    if !is_free_of_nan(values) {
        return Err(bad("statistics input must not contain NaN"));
    }
    Ok(())
}

/// Sort in place into the order the C++ Release build's `std::sort` leaves.
///
/// Every sorting entry point of this header calls `std::sort(begin, end)` with
/// the default `operator<` and no comparator — `median` at
/// `MATH/StatisticFunctions.h:140`, `quantile1st` at `:244`, `quantile3rd` at
/// `:281`, `MAD` through its own `median` call at `:189`, and
/// `SummaryStatistics`'s unqualified `sort(data.begin(), data.end())` at `:948`
/// (core `bc9cc12`). This function is that call:
/// [`source_sort_by`] with
/// `|a, b| a < b`, which reproduces the libstdc++ introsort comparison by
/// comparison and move by move.
///
/// It replaces an `f64::total_cmp` sort, and the difference is visible in
/// exactly the two places `operator<` is not a total order:
///
/// - **a NaN**, which `operator<` makes incomparable with everything. The
///   permutation is then unspecified rather than undefined *as an operation* —
///   libstdc++ still runs, deterministically, and decision D16 of
///   `docs/EARLY_TOPP_WORK_PACKAGES.md` puts reproducing it in scope;
/// - **a signed zero**, which `operator<` makes *equivalent* to a plain zero.
///   `total_cmp` put `-0.0` first; libstdc++ leaves a short range as it found
///   it, and that is what the Release build prints.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the introsort's unbounded partition or
/// final-insertion loop would read outside the vector, which is the one thing
/// [`source_sort_by`] refuses, and
/// which no asymmetric comparison — `<` on `f64` keys, NaN keys included — can
/// provoke; when the owned copy the permutation is applied through cannot be
/// allocated; and — unreachably, but checked rather than left to panic — when
/// the sort returns something that is not a permutation of the input. `values`
/// is left in its original order in every case.
fn sort_ascending(values: &mut [f64]) -> Result<()> {
    let mut owned: Vec<f64> = Vec::new();
    owned
        .try_reserve_exact(values.len())
        .map_err(|_| bad("cannot allocate the sort buffer"))?;
    owned.extend_from_slice(values);
    source_sort_by(&mut owned, |a, b| a < b)?;
    // `source_sort_by` applies a permutation of `0..len`, so the length cannot
    // change; `copy_from_slice` would panic rather than refuse if it ever did,
    // and this module may not panic on an input.
    if owned.len() != values.len() {
        return Err(bad("the sort did not return a permutation"));
    }
    values.copy_from_slice(&owned);
    Ok(())
}

/// Fail when a range is empty.
///
/// Port of `checkIteratorsNotNULL`, which throws `Exception::InvalidRange` when
/// `begin == end`.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] when `values` is empty.
pub fn check_not_empty<T>(values: &[T]) -> Result<()> {
    if values.is_empty() {
        return Err(empty("statistics range must not be empty"));
    }
    Ok(())
}

/// Fail when a range still has elements left.
///
/// Port of `checkIteratorsEqual`, which throws `Exception::InvalidRange` when
/// the two iterators are *not* equal. The source uses it after a lock-step loop
/// to prove the second range was consumed exactly.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] when `values` is not empty.
pub fn check_exhausted<T>(values: &[T]) -> Result<()> {
    if !values.is_empty() {
        return Err(empty("statistics ranges must have the same length"));
    }
    Ok(())
}

/// Fail when exactly one of two ranges is exhausted.
///
/// Port of `checkIteratorsAreValid(begin_b, end_b, begin_a, end_a)`, whose test
/// is the exclusive-or `(begin_b == end_b) ^ (begin_a == end_a)`. Both empty is
/// accepted and both non-empty is accepted; only a disagreement is an error.
/// The argument order is the source's, `b` first.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] when one slice is empty and the other is not.
pub fn check_ranges_end_together<T, U>(b: &[T], a: &[U]) -> Result<()> {
    if b.is_empty() != a.is_empty() {
        return Err(empty("statistics ranges do not end simultaneously"));
    }
    Ok(())
}

/// Sum of a range of values.
///
/// `std::accumulate(begin, end, 0.0)`: a left fold from `0.0`, in slice order.
/// An empty range sums to `0.0`, as in the source, which does not check.
///
/// Each addition is `x86_64::add`, so `(+inf) + (-inf)` is the Release
/// build's `0xfff8000000000000` rather than the host's default NaN; see the
/// module's "NaN bit patterns" section. No finite sum changes.
pub fn sum(values: &[f64]) -> f64 {
    values
        .iter()
        .fold(0.0, |total, value| x86_64::add(total, *value))
}

/// Arithmetic mean of a range of values.
///
/// Computed as [`sum`] divided by the element count, in that order.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range, as the source's
/// `checkIteratorsNotNULL` throws `Exception::InvalidRange`.
pub fn mean(values: &[f64]) -> Result<f64> {
    check_not_empty(values)?;
    Ok(x86_64::div(sum(values), values.len() as f64))
}

/// Median of an already ascending range.
///
/// The source's `median(begin, end, /*sorted=*/true)`. For an even count the
/// two middle values are averaged as `(a + b) / 2.0`; for an odd count the
/// single middle value is returned unchanged. The return type is floating point
/// precisely because the even case has to average.
///
/// The source does not verify the claimed sortedness in a release build; this
/// port does, because the alternative is a silently wrong quantile.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range and
/// [`Error::UnsortedData`] when `values` is not ascending or contains a NaN,
/// including the single-value range `[NaN]`.
pub fn median_sorted(values: &[f64]) -> Result<f64> {
    check_not_empty(values)?;
    if !is_ascending(values) {
        return Err(Error::UnsortedData);
    }
    Ok(median_of_sorted(values))
}

/// Median of an ascending, non-empty slice. Callers have already checked both.
///
/// The even case averages the two middle values, `(*it1 + *it2) / 2.0` in the
/// source. That is `addsd` then `divsd`, so it is built on
/// [`x86_64::add`]/[`x86_64::div`]: `{-inf, +inf}` generates a NaN here, and
/// its bits are the Release build's.
fn median_of_sorted(values: &[f64]) -> f64 {
    let size = values.len();
    if size % 2 == 0 {
        x86_64::div(x86_64::add(values[size / 2 - 1], values[size / 2]), 2.0)
    } else {
        values[(size - 1) / 2]
    }
}

/// Median of a range, sorting it in place first.
///
/// The source's `median(begin, end, /*sorted=*/false)`, which likewise sorts
/// the caller's range; the `&mut` receiver makes that visible in the signature.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range. A NaN is **not**
/// refused: the range is sorted into the permutation the Release build's
/// `std::sort` leaves and the statistic is read positionally out of it, as the
/// source does; see the module's NaN section and decision D16. The only
/// [`Error::InvalidValue`] left is the one
/// [`source_sort_by`] raises when the
/// introsort would read outside the vector, which `<` on `f64` keys cannot
/// provoke, and it is raised before anything is written back, so a rejected
/// call leaves the caller's range in its original order.
pub fn median(values: &mut [f64]) -> Result<f64> {
    check_not_empty(values)?;
    sort_ascending(values)?;
    Ok(median_of_sorted(values))
}

/// Median absolute deviation, `median(|x_i - median_of_numbers|)`.
///
/// Sortedness of the input is neither required nor exploited. The median must
/// be supplied because the caller has usually computed it already; the source
/// says so explicitly and this port keeps the same contract rather than
/// recomputing it.
///
/// The absolute differences are staged into an owned buffer and sorted, exactly
/// as the source does.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range, and also when the range
/// exceeds [`MAX_ITEMS`] or [`MAX_BYTES`]; the ceiling is checked before the
/// buffer is allocated, and the source has no such ceiling.
///
/// Returns [`Error::InvalidValue`] only where
/// [`sort_ascending`](crate::math::statistic_functions) does: an introsort read
/// outside the vector, which `<` on `f64` keys cannot provoke. A NaN — in the
/// input, in `median_of_numbers`, or produced by an `inf - inf` difference that
/// neither input had — is **reproduced**, not refused: the source stages the
/// same differences and hands them to the same `std::sort`, and under decision
/// D16 that permutation is what this port reproduces.
pub fn mad(values: &[f64], median_of_numbers: f64) -> Result<f64> {
    check_not_empty(values)?;
    preflight(values.len(), size_of::<f64>())?;
    let mut diffs = Vec::with_capacity(values.len());
    for value in values {
        diffs.push(x86_64::abs(x86_64::sub(*value, median_of_numbers)));
    }
    sort_ascending(&mut diffs)?;
    Ok(median_of_sorted(&diffs))
}

/// Mean absolute deviation, `mean(|x_i - mean_of_numbers|)`.
///
/// The mean must be supplied, as in the source. The divisor is the element
/// count `n`, not `n - 1`.
///
/// An empty range returns NaN (`0.0 / 0.0`): the source neither checks nor
/// documents this, but its class test asserts the NaN, so the value is pinned
/// behaviour rather than an oversight this port may quietly change. Callers
/// that need an error should use [`absdev`], which checks.
pub fn mean_absolute_deviation(values: &[f64], mean_of_numbers: f64) -> f64 {
    let mut total = 0.0;
    for value in values {
        total = x86_64::add(total, x86_64::abs(x86_64::sub(*value, mean_of_numbers)));
    }
    x86_64::div(total, values.len() as f64)
}

/// Mean absolute deviation about the range's own mean.
///
/// The source's `absdev(begin, end)` with the default `mean` argument. That
/// default is the sentinel `std::numeric_limits<double>::max()`, so a caller
/// that genuinely wants the deviation about `DBL_MAX` cannot express it; the
/// port splits the overload into two functions instead, which removes the
/// sentinel entirely.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range.
pub fn absdev(values: &[f64]) -> Result<f64> {
    check_not_empty(values)?;
    let mean_value = mean(values)?;
    Ok(mean_absolute_deviation(values, mean_value))
}

/// Mean absolute deviation about an explicitly supplied mean.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range.
pub fn absdev_with_mean(values: &[f64], mean_of_numbers: f64) -> Result<f64> {
    check_not_empty(values)?;
    Ok(mean_absolute_deviation(values, mean_of_numbers))
}

/// First quartile of an already ascending range, by the median-of-halves rule.
///
/// The range is halved and the median of the lower half is returned. For fewer
/// than three values the lower half is empty and the minimum is returned, which
/// is what the size-3 and size-4 cases produce anyway.
///
/// For an **even** count the source takes the median of `[0, n/2 - 1)`, one
/// element short of the lower half, because its `-1` is written to "exclude the
/// median value" — an exclusion that only makes sense for an odd count. The
/// port reproduces the asymmetry: for `n = 10` the lower half is
/// `[0, 4)`, not `[0, 5)`. See `docs/STATISTIC_FUNCTIONS_SUPPORT.md`.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range and
/// [`Error::UnsortedData`] when `values` is not ascending or contains a NaN,
/// including the single-value range `[NaN]`.
pub fn quantile1st_sorted(values: &[f64]) -> Result<f64> {
    check_not_empty(values)?;
    if !is_ascending(values) {
        return Err(Error::UnsortedData);
    }
    Ok(quantile1st_of_sorted(values))
}

/// First quartile of an ascending, non-empty slice. Callers have already
/// checked both, exactly as for [`median_of_sorted`].
fn quantile1st_of_sorted(values: &[f64]) -> f64 {
    let size = values.len();
    if size < 3 {
        return values[0];
    }
    if size % 2 == 0 {
        return median_of_sorted(&values[..size / 2 - 1]);
    }
    median_of_sorted(&values[..size / 2])
}

/// First quartile of a range, sorting it in place first.
///
/// The source's `quantile1st(begin, end, /*sorted=*/false)`.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range. A NaN is **not**
/// refused: the range is sorted into the permutation the Release build's
/// `std::sort` leaves and the statistic is read positionally out of it, as the
/// source does; see the module's NaN section and decision D16. The only
/// [`Error::InvalidValue`] left is the one
/// [`source_sort_by`] raises when the
/// introsort would read outside the vector, which `<` on `f64` keys cannot
/// provoke, and it is raised before anything is written back, so a rejected
/// call leaves the caller's range in its original order.
pub fn quantile1st(values: &mut [f64]) -> Result<f64> {
    check_not_empty(values)?;
    sort_ascending(values)?;
    Ok(quantile1st_of_sorted(values))
}

/// Third quartile of an already ascending range, by the median-of-halves rule.
///
/// The range is halved and the median of the upper half is returned. For fewer
/// than three values the upper half is empty and the maximum is returned.
///
/// The upper half is `[n/2 + 1, n)`, whose `+1` again excludes one element for
/// an even count; see [`quantile1st_sorted`] for the same asymmetry.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range and
/// [`Error::UnsortedData`] when `values` is not ascending or contains a NaN,
/// including the single-value range `[NaN]`.
pub fn quantile3rd_sorted(values: &[f64]) -> Result<f64> {
    check_not_empty(values)?;
    if !is_ascending(values) {
        return Err(Error::UnsortedData);
    }
    Ok(quantile3rd_of_sorted(values))
}

/// Third quartile of an ascending, non-empty slice. Callers have already
/// checked both, exactly as for [`median_of_sorted`].
fn quantile3rd_of_sorted(values: &[f64]) -> f64 {
    let size = values.len();
    if size < 3 {
        return values[size - 1];
    }
    median_of_sorted(&values[size / 2 + 1..])
}

/// Third quartile of a range, sorting it in place first.
///
/// The source's `quantile3rd(begin, end, /*sorted=*/false)`.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range. A NaN is **not**
/// refused: the range is sorted into the permutation the Release build's
/// `std::sort` leaves and the statistic is read positionally out of it, as the
/// source does; see the module's NaN section and decision D16. The only
/// [`Error::InvalidValue`] left is the one
/// [`source_sort_by`] raises when the
/// introsort would read outside the vector, which `<` on `f64` keys cannot
/// provoke, and it is raised before anything is written back, so a rejected
/// call leaves the caller's range in its original order.
pub fn quantile3rd(values: &mut [f64]) -> Result<f64> {
    check_not_empty(values)?;
    sort_ascending(values)?;
    Ok(quantile3rd_of_sorted(values))
}

/// The `q`-quantile of an already ascending range, Hyndman-Fan type 7.
///
/// ```text
/// pos      = q * (n - 1)
/// idx      = floor(pos),  frac = pos - idx
/// quantile = (1 - frac) * x[idx] + frac * x[idx + 1]
/// ```
///
/// `q == 0` returns the minimum and `q == 1` the maximum exactly; a zero
/// fractional part short-circuits to `x[idx]` without touching `x[idx + 1]`,
/// which is what keeps `q == 1` from reading past the end.
///
/// A single value is returned unchanged for every `q`.
///
/// This is a different convention from [`quantile1st_sorted`] and
/// [`quantile3rd_sorted`], which interpolate nothing; the source carries both
/// and so does the port.
///
/// # Arguments
///
/// * `values` — ascending, non-empty.
/// * `q` — in `[0, 1]`.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range, [`Error::InvalidValue`]
/// when `q` is outside `[0, 1]` or is NaN, and [`Error::UnsortedData`] when
/// `values` is not ascending **or contains a NaN**, including the single-value
/// range `[NaN]`, which has no adjacent pair to disagree. The source states the
/// sortedness precondition as `@pre` and checks it only through
/// `OPENMS_PRECONDITION`, which is compiled out of a release build; this port
/// always checks.
pub fn quantile(values: &[f64], q: f64) -> Result<f64> {
    check_not_empty(values)?;
    if !(0.0..=1.0).contains(&q) {
        return Err(bad("quantile q must be in [0,1]"));
    }
    if !is_ascending(values) {
        return Err(Error::UnsortedData);
    }
    let n = values.len();
    if n == 1 {
        return Ok(values[0]);
    }
    let pos = q * (n - 1) as f64;
    let index = pos.floor();
    // `pos` is at most `n - 1`, so the truncation is exact and in range.
    let i = index as usize;
    let frac = pos - index;
    if frac == 0.0 {
        return Ok(values[i]);
    }
    Ok(x86_64::add(
        x86_64::mul(x86_64::sub(1.0, frac), values[i]),
        x86_64::mul(frac, values[i + 1]),
    ))
}

/// Tukey upper fence `Q3 + k * IQR` over the finite values of a range.
///
/// Non-finite values are dropped. With fewer than four finite values, or a
/// non-positive interquartile range, the fence is `+inf`, which every consumer
/// reads as "no fence". The quartiles are the interpolating [`quantile`] at
/// `0.25` and `0.75`, not [`quantile1st_sorted`]/[`quantile3rd_sorted`].
///
/// Reference: J. W. Tukey (1977). Exploratory Data Analysis.
///
/// # Arguments
///
/// * `k` — Tukey factor; the source's default is [`DEFAULT_TUKEY_FACTOR`].
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] when the range exceeds [`MAX_ITEMS`] or
/// [`MAX_BYTES`]; the source has no ceiling on the copy it makes. An empty
/// range is *not* an error here — the source does not check it, and fewer than
/// four finite values already yields `+inf`.
pub fn tukey_upper_fence(values: &[f64], k: f64) -> Result<f64> {
    let finite = finite_sorted_copy(values)?;
    if finite.len() < 4 {
        return Ok(f64::INFINITY);
    }
    let q1 = quantile(&finite, 0.25)?;
    let q3 = quantile(&finite, 0.75)?;
    let iqr = q3 - q1;
    // The source writes `if (!(iqr > 0.0))`, so a NaN interquartile range also
    // yields the "no fence" answer rather than falling through.
    if iqr.is_nan() || iqr <= 0.0 {
        return Ok(f64::INFINITY);
    }
    Ok(q3 + k * iqr)
}

/// Fraction of the finite values that are strictly greater than `threshold`.
///
/// Non-finite values are excluded from both the numerator and the denominator.
/// With no finite value at all the fraction is `0.0`, as in the source, which
/// checks that divisor explicitly.
pub fn tail_fraction_above(values: &[f64], threshold: f64) -> f64 {
    let mut n: usize = 0;
    let mut n_tail: usize = 0;
    for &value in values {
        if !value.is_finite() {
            continue;
        }
        n += 1;
        if value > threshold {
            n_tail += 1;
        }
    }
    if n == 0 {
        return 0.0;
    }
    n_tail as f64 / n as f64
}

/// The `q`-quantile after winsorizing the finite values at an upper fence.
///
/// Finite values are copied, capped above at `upper_fence` and below at `0.0`,
/// then sorted and passed to [`quantile`]. The lower cap is the source's, whose
/// comment calls it defensive and useful for absolute residuals; it means the
/// function is **not** a general winsorizer for signed data.
///
/// A non-finite `upper_fence` disables the capping and the raw quantile is
/// returned. An input with no finite value returns `0.0`.
///
/// Reference: J. W. Tukey (1962). The Future of Data Analysis.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `q` is outside `[0, 1]`, and
/// [`Error::InvalidRange`] when the range exceeds [`MAX_ITEMS`] or
/// [`MAX_BYTES`].
pub fn winsorized_quantile(values: &[f64], q: f64, upper_fence: f64) -> Result<f64> {
    preflight(values.len(), size_of::<f64>())?;
    let mut kept = Vec::with_capacity(values.len());
    for &value in values {
        if !value.is_finite() {
            continue;
        }
        kept.push(value);
    }
    if kept.is_empty() {
        return Ok(0.0);
    }
    if upper_fence.is_finite() {
        for value in &mut kept {
            if *value > upper_fence {
                *value = upper_fence;
            }
            if *value < 0.0 {
                *value = 0.0;
            }
        }
    }
    sort_ascending(&mut kept)?;
    quantile(&kept, q)
}

/// Everything [`adaptive_quantile`] computed on the way to its blended value.
///
/// Field names and defaults are the source's `AdaptiveQuantileResult`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveQuantileResult {
    /// The final blended (adaptive) quantile.
    pub blended: f64,
    /// The raw `q`-quantile of the finite values.
    pub half_raw: f64,
    /// The `q`-quantile after interquartile-range winsorization.
    pub half_rob: f64,
    /// Tukey upper fence `Q3 + k * IQR`, `+inf` when undefined.
    pub upper_fence: f64,
    /// Fraction of values above `upper_fence`.
    pub tail_fraction: f64,
    /// Blend weight in `[0, 1]`: `0` is fully robust, `1` fully raw.
    pub weight: f64,
}

impl Default for AdaptiveQuantileResult {
    /// The source's member initialisers: zero everywhere except `upper_fence`,
    /// which starts at `+inf` because "no fence" is not the same as a fence at
    /// zero. This is the value [`adaptive_quantile`] returns for an input with
    /// no finite value.
    fn default() -> Self {
        Self {
            blended: 0.0,
            half_raw: 0.0,
            half_rob: 0.0,
            upper_fence: f64::INFINITY,
            tail_fraction: 0.0,
            weight: 0.0,
        }
    }
}

/// A quantile that blends the raw and winsorized values by tail density.
///
/// With `UF` the Tukey upper fence of the finite inputs:
///
/// ```text
/// half_raw = quantile(values, q)
/// half_rob = winsorized_quantile(values, q, UF)
/// r        = tail_fraction_above(values, UF)      (0 when UF is not finite)
/// w        = 0                     for r <= r_sparse
///          = 1                     for r >= r_dense
///          = (r - r_sparse) / (r_dense - r_sparse) in between, clamped
/// blended  = (1 - w) * half_rob + w * half_raw
/// ```
///
/// Sparse outliers therefore leave the window where the robust estimate puts
/// it, while a genuinely broad tail pulls it back to the raw quantile.
///
/// A degenerate configuration with `r_dense <= r_sparse` is the source's own
/// step function: `w` is `1` when `r > r_sparse` and `0` otherwise.
///
/// References: J. W. Tukey (1962, 1977); R. J. Hyndman, Y. Fan (1996), Sample
/// Quantiles in Statistical Packages.
///
/// # Arguments
///
/// * `q` — target quantile in `[0, 1]`, e.g. `0.99` for a 99% half-width.
/// * `k` — Tukey factor; source default [`DEFAULT_TUKEY_FACTOR`].
/// * `r_sparse` — source default [`DEFAULT_R_SPARSE`].
/// * `r_dense` — source default [`DEFAULT_R_DENSE`].
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `q` is outside `[0, 1]` or when `k`,
/// `r_sparse` or `r_dense` is not finite, and [`Error::InvalidRange`] when the
/// range exceeds [`MAX_ITEMS`] or [`MAX_BYTES`]. The finiteness check is
/// native: the source accepts a NaN threshold and lets it decide the blend
/// weight through comparisons that are all false. An input with no finite value
/// yields `AdaptiveQuantileResult::default()`, not an error.
pub fn adaptive_quantile(
    values: &[f64],
    q: f64,
    k: f64,
    r_sparse: f64,
    r_dense: f64,
) -> Result<AdaptiveQuantileResult> {
    if !k.is_finite() || !r_sparse.is_finite() || !r_dense.is_finite() {
        return Err(bad("adaptive quantile parameters must be finite"));
    }
    let finite = finite_sorted_copy(values)?;
    if finite.is_empty() {
        return Ok(AdaptiveQuantileResult::default());
    }
    let half_raw = quantile(&finite, q)?;
    let upper_fence = tukey_upper_fence(&finite, k)?;
    let tail_fraction = if upper_fence.is_finite() {
        tail_fraction_above(&finite, upper_fence)
    } else {
        0.0
    };
    let half_rob = winsorized_quantile(&finite, q, upper_fence)?;
    let weight = if r_dense <= r_sparse {
        if tail_fraction > r_sparse { 1.0 } else { 0.0 }
    } else {
        let t = (tail_fraction - r_sparse) / (r_dense - r_sparse);
        // `std::max(0.0, std::min(1.0, t))`, spelled out so the clamp keeps the
        // source's branch order rather than `f64::clamp`'s NaN propagation.
        let capped = if t < 1.0 { t } else { 1.0 };
        if 0.0 < capped { capped } else { 0.0 }
    };
    Ok(AdaptiveQuantileResult {
        blended: (1.0 - weight) * half_rob + weight * half_raw,
        half_raw,
        half_rob,
        upper_fence,
        tail_fraction,
        weight,
    })
}

/// Copy the finite values of a range and sort them ascending.
fn finite_sorted_copy(values: &[f64]) -> Result<Vec<f64>> {
    preflight(values.len(), size_of::<f64>())?;
    let mut kept = Vec::with_capacity(values.len());
    for &value in values {
        if value.is_finite() {
            kept.push(value);
        }
    }
    sort_ascending(&mut kept)?;
    Ok(kept)
}

/// Sample variance about the range's own mean, with `n - 1` degrees of freedom.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range and for a single value.
/// The source divides by `n - 1` unchecked and so returns NaN for `n == 1`;
/// this port refuses instead, because a NaN variance is the exact failure this
/// layer propagates into everything built on it. The source's own
/// `SummaryStatistics` acknowledges the problem in a comment and substitutes
/// `0.0`, which [`SummaryStatistics`] reproduces.
pub fn variance(values: &[f64]) -> Result<f64> {
    let mean_value = mean(values)?;
    variance_with_mean(values, mean_value)
}

/// Sample variance about an explicitly supplied mean, `n - 1` degrees of
/// freedom.
///
/// The source expresses this as one function whose `mean` argument defaults to
/// the sentinel `std::numeric_limits<double>::max()`, meaning "compute it"; a
/// caller that passes `DBL_MAX` deliberately silently gets the wrong answer.
/// The port splits the overload in two and has no sentinel.
///
/// Squared deviations accumulate in slice order, `diff * diff` per element,
/// exactly as the source writes them.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range and for a single value;
/// see [`variance`].
pub fn variance_with_mean(values: &[f64], mean_of_numbers: f64) -> Result<f64> {
    check_not_empty(values)?;
    if values.len() < 2 {
        return Err(empty("variance needs at least two values"));
    }
    let mut sum_value = 0.0;
    for value in values {
        let diff = x86_64::sub(*value, mean_of_numbers);
        sum_value = x86_64::add(sum_value, x86_64::mul(diff, diff));
    }
    Ok(x86_64::div(sum_value, (values.len() - 1) as f64))
}

/// Sample standard deviation about the range's own mean.
///
/// The square root of [`variance`]; the source computes it the same way, so the
/// rounding of the intermediate variance is preserved.
///
/// # Errors
///
/// As [`variance`].
pub fn sd(values: &[f64]) -> Result<f64> {
    Ok(x86_64::sqrt(variance(values)?))
}

/// Sample standard deviation about an explicitly supplied mean.
///
/// # Errors
///
/// As [`variance_with_mean`].
pub fn sd_with_mean(values: &[f64], mean_of_numbers: f64) -> Result<f64> {
    Ok(x86_64::sqrt(variance_with_mean(values, mean_of_numbers)?))
}

/// Sample covariance of two equally long ranges, `n - 1` degrees of freedom.
///
/// Each range's mean is computed over that range, then the cross products
/// accumulate in slice order.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range, for ranges of different
/// length, and for a single pair. The source checks the lengths only by
/// re-testing the two *begin* iterators inside the loop — a test that never
/// changes — and by comparing the second iterator to its end afterwards, so a
/// short second range is read out of bounds before the mismatch is noticed;
/// this port compares the lengths up front. `n == 1` divides by zero in the
/// source and is refused here, as in [`variance`].
pub fn covariance(a: &[f64], b: &[f64]) -> Result<f64> {
    check_not_empty(a)?;
    if a.len() != b.len() {
        return Err(empty("covariance ranges must have the same length"));
    }
    if a.len() < 2 {
        return Err(empty("covariance needs at least two pairs"));
    }
    let mean_a = mean(a)?;
    let mean_b = mean(b)?;
    let mut sum_value = 0.0;
    for (value_a, value_b) in a.iter().zip(b.iter()) {
        let product = x86_64::mul(x86_64::sub(*value_a, mean_a), x86_64::sub(*value_b, mean_b));
        sum_value = x86_64::add(sum_value, product);
    }
    Ok(x86_64::div(sum_value, (a.len() - 1) as f64))
}

/// Mean square error between two equally long ranges.
///
/// Squared differences accumulate in slice order and are divided by the element
/// count `n` — not `n - 1`, unlike [`variance`] and [`covariance`].
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range or ranges of different
/// length, matching the source's `Exception::InvalidRange`.
pub fn mean_square_error(a: &[f64], b: &[f64]) -> Result<f64> {
    check_not_empty(a)?;
    if a.len() != b.len() {
        return Err(empty("mean square error ranges must have the same length"));
    }
    let mut error = 0.0;
    for (value_a, value_b) in a.iter().zip(b.iter()) {
        let tmp = x86_64::sub(*value_a, *value_b);
        error = x86_64::add(error, x86_64::mul(tmp, tmp));
    }
    Ok(x86_64::div(error, a.len() as f64))
}

/// Root mean square error, the square root of [`mean_square_error`].
///
/// # Errors
///
/// As [`mean_square_error`].
pub fn root_mean_square_error(a: &[f64], b: &[f64]) -> Result<f64> {
    Ok(x86_64::sqrt(mean_square_error(a, b)?))
}

/// Fraction of positions where two label ranges agree in sign.
///
/// A position counts as wrong when exactly one of the two values is negative;
/// zero counts as non-negative on both sides. The result is the number of
/// agreeing positions divided by `n`.
///
/// A NaN on either side makes both of the source's comparisons false, so the
/// position is counted as *agreeing*. That is the source's arithmetic, not a
/// substitution — the comparisons are well defined in C++ and produce the same
/// answer — so it is reproduced rather than refused; nothing here sorts, and
/// the module's NaN section says why that is the dividing line.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range or ranges of different
/// length.
pub fn classification_rate(a: &[f64], b: &[f64]) -> Result<f64> {
    check_not_empty(a)?;
    if a.len() != b.len() {
        return Err(empty("classification ranges must have the same length"));
    }
    let mut correct = a.len();
    for (&value_a, &value_b) in a.iter().zip(b.iter()) {
        if (value_a < 0.0 && value_b >= 0.0) || (value_a >= 0.0 && value_b < 0.0) {
            correct -= 1;
        }
    }
    Ok(correct as f64 / a.len() as f64)
}

/// Matthews correlation coefficient of predicted against real labels.
///
/// `a` holds the predicted labels and `b` the real ones; a value is positive
/// when it is `>= 0`. With the four counts of the confusion matrix,
/// `(tp*tn - fp*fn) / sqrt((tp+fp) * (tp+fn) * (tn+fp) * (tn+fn))`.
///
/// A zero denominator yields NaN. That is not a substitution: a zero factor
/// forces two of the four counts to zero, which makes the numerator exactly
/// zero as well, so the source's unchecked `0 / 0` is NaN in every case. The
/// port returns it explicitly rather than dividing.
///
/// A NaN on either side of a pair fails all four comparisons, so that pair
/// increments none of the counts and is silently excluded from the confusion
/// matrix. As in [`classification_rate`] the comparisons are well defined in
/// C++ and give the same answer, so the behaviour is reproduced rather than
/// refused.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range or ranges of different
/// length. The source's own emptiness check compares `begin_a` against `end_b`
/// — iterators into two different containers, which is undefined behaviour and
/// is not the check it intends; the port compares each range with itself.
pub fn matthews_correlation_coefficient(a: &[f64], b: &[f64]) -> Result<f64> {
    check_not_empty(a)?;
    if a.len() != b.len() {
        return Err(empty("Matthews ranges must have the same length"));
    }
    let mut tp: f64 = 0.0;
    let mut fp: f64 = 0.0;
    let mut tn: f64 = 0.0;
    let mut fn_: f64 = 0.0;
    for (&value_a, &value_b) in a.iter().zip(b.iter()) {
        if value_a < 0.0 && value_b >= 0.0 {
            fn_ += 1.0;
        } else if value_a < 0.0 && value_b < 0.0 {
            tn += 1.0;
        } else if value_a >= 0.0 && value_b >= 0.0 {
            tp += 1.0;
        } else if value_a >= 0.0 && value_b < 0.0 {
            fp += 1.0;
        }
    }
    let denominator = ((tp + fp) * (tp + fn_) * (tn + fp) * (tn + fn_)).sqrt();
    if denominator == 0.0 {
        return Ok(f64::NAN);
    }
    Ok((tp * tn - fp * fn_) / denominator)
}

/// Pearson (linear) correlation coefficient of two equally long ranges.
///
/// Both means are formed first, then the numerator and the two denominator sums
/// accumulate together in one pass, in slice order — the source's arrangement,
/// kept because a two-pass or fused rewrite rounds differently.
///
/// A range whose values are all equal makes its denominator zero and the result
/// NaN; the source documents exactly that and its class test pins it. As in
/// [`matthews_correlation_coefficient`] a zero denominator also forces a zero
/// numerator, so the NaN is returned explicitly instead of dividing. The one
/// case this treats differently from the source is a denominator that underflows
/// to zero from non-zero deviations, where the source yields an infinity; NaN
/// is the port's answer there.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range or ranges of different
/// length.
pub fn pearson_correlation_coefficient(a: &[f64], b: &[f64]) -> Result<f64> {
    check_not_empty(a)?;
    if a.len() != b.len() {
        return Err(empty("Pearson ranges must have the same length"));
    }
    let dist = a.len() as f64;
    let avg_a = sum(a) / dist;
    let avg_b = sum(b) / dist;
    let mut numerator = 0.0;
    let mut denominator_a = 0.0;
    let mut denominator_b = 0.0;
    for (&value_a, &value_b) in a.iter().zip(b.iter()) {
        let temp_a = value_a - avg_a;
        let temp_b = value_b - avg_b;
        numerator += temp_a * temp_b;
        denominator_a += temp_a * temp_a;
        denominator_b += temp_b * temp_b;
    }
    let denominator = (denominator_a * denominator_b).sqrt();
    if denominator == 0.0 {
        return Ok(f64::NAN);
    }
    Ok(numerator / denominator)
}

/// Relative tolerance of the source's `computeRank` tie test.
///
/// Two neighbouring values are a tie when they differ by no more than this
/// fraction of the *later* value's magnitude.
pub const COMPUTE_RANK_TIE_TOLERANCE: f64 = 0.000_000_1;

/// Replace every value by its rank, ties sharing the mean of their ranks.
///
/// Ranks are one-based. Values are ordered ascending together with their
/// original positions, and the ranks are written back in the original order.
///
/// Ties are **not** decided by exact equality: neighbours count as tied when
/// `|x[i+1] - x[i]| <= COMPUTE_RANK_TIE_TOLERANCE * |x[i+1]|`, a *relative*
/// test against the later value. Two consequences follow from the source's
/// formulation and are reproduced here. The tolerance collapses to zero when
/// the later value is zero, so a pair straddling zero is only a tie if both are
/// exactly zero; and the test is applied pairwise while scanning, so a run of
/// values each within tolerance of its neighbour is one tie block even when its
/// ends are far apart.
///
/// A tie block spanning sorted positions `i..z` receives the rank
/// `0.5 * (i + z + 1)`, which is the mean of the one-based ranks `i+1 ..= z`.
///
/// This is the ranking [`rank_correlation_coefficient`] is built on, and it is
/// a different function from [`crate::math::rank_data::rankdata`], which
/// reproduces SciPy: that one compares by exact equality and offers five tie
/// rules. Neither is a replacement for the other.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] when the slice exceeds [`MAX_ITEMS`] or
/// [`MAX_BYTES`]; the staging buffer is allocated only after that check. An
/// empty slice is a no-op: the source computes `w.size() - 1` in unsigned
/// arithmetic and wraps to `SIZE_MAX` on an empty vector, which this port
/// cannot and does not reproduce.
///
/// Returns [`Error::InvalidValue`] when any value is NaN. The ranking sorts,
/// and a NaN would additionally make the tie test — whose two comparisons are
/// both false against a NaN — split every block that touches it; see the
/// module's NaN section. The refusal happens before anything is written, so a
/// rejected call leaves `w` unchanged.
pub fn compute_rank(w: &mut [f64]) -> Result<()> {
    if w.is_empty() {
        return Ok(());
    }
    check_no_nan(w)?;
    preflight(w.len(), size_of::<(usize, f64)>())?;
    let mut w_idx: Vec<(usize, f64)> = Vec::with_capacity(w.len());
    for (index, &value) in w.iter().enumerate() {
        w_idx.push((index, value));
    }
    w_idx.sort_by(|left, right| left.1.total_cmp(&right.1));

    let n = w.len() - 1;
    let mut i = 0usize;
    while i < n {
        let tied = (w_idx[i + 1].1 - w_idx[i].1).abs()
            <= COMPUTE_RANK_TIE_TOLERANCE * w_idx[i + 1].1.abs();
        if !tied {
            w_idx[i].1 = (i + 1) as f64;
            i += 1;
        } else {
            let mut z = i + 1;
            while z <= n
                && (w_idx[z].1 - w_idx[i].1).abs() <= COMPUTE_RANK_TIE_TOLERANCE * w_idx[z].1.abs()
            {
                z += 1;
            }
            let rank = 0.5 * (i + z + 1) as f64;
            for entry in w_idx.iter_mut().take(z).skip(i) {
                entry.1 = rank;
            }
            i = z;
        }
    }
    if i == n {
        w_idx[n].1 = (n + 1) as f64;
    }
    for &(origin, rank) in &w_idx {
        w[origin] = rank;
    }
    Ok(())
}

/// Spearman rank correlation coefficient of two equally long ranges.
///
/// Both ranges are replaced by their [`compute_rank`] ranks and correlated
/// about the theoretical mean rank `(n + 1) / 2` rather than the observed mean
/// of the ranks. The two differ whenever ties are present, and the source keeps
/// the theoretical value deliberately — a comment records the earlier integer
/// division `(n + 1) / 2` it replaced.
///
/// A range whose ranks are all equal gives a zero sum of squares, and the
/// source returns `0.0` for it rather than NaN. That is reproduced: a constant
/// range correlates to zero here, unlike
/// [`pearson_correlation_coefficient`], which returns NaN.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range, for ranges of different
/// length, and when a range exceeds [`MAX_ITEMS`] or [`MAX_BYTES`]; and
/// [`Error::InvalidValue`] when either range contains a NaN, which
/// [`compute_rank`] cannot order. Both ranges are checked before either is
/// copied.
pub fn rank_correlation_coefficient(a: &[f64], b: &[f64]) -> Result<f64> {
    check_not_empty(a)?;
    if a.len() != b.len() {
        return Err(empty("rank correlation ranges must have the same length"));
    }
    check_no_nan(a)?;
    check_no_nan(b)?;
    preflight(a.len(), size_of::<f64>())?;
    let mut ranks_model = a.to_vec();
    let mut ranks_data = b.to_vec();
    compute_rank(&mut ranks_data)?;
    compute_rank(&mut ranks_model)?;

    let mu = (ranks_data.len() + 1) as f64 / 2.0;
    let mut sum_model_data = 0.0;
    let mut sqsum_data = 0.0;
    let mut sqsum_model = 0.0;
    for (&rank_data, &rank_model) in ranks_data.iter().zip(ranks_model.iter()) {
        sum_model_data += (rank_data - mu) * (rank_model - mu);
        sqsum_data += (rank_data - mu) * (rank_data - mu);
        sqsum_model += (rank_model - mu) * (rank_model - mu);
    }
    if sqsum_data == 0.0 || sqsum_model == 0.0 {
        return Ok(0.0);
    }
    Ok(sum_model_data / (sqsum_data.sqrt() * sqsum_model.sqrt()))
}

/// Count, mean, variance, extrema, quartiles and median of one sample.
///
/// Port of the source's `SummaryStatistics<T>` helper. Every field is the value
/// the corresponding free function would return, with two documented
/// substitutions the source makes itself.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SummaryStatistics {
    /// Number of values summarised.
    pub count: usize,
    /// Arithmetic mean, `0.0` when there are no values.
    pub mean: f64,
    /// Sample variance with `n - 1` degrees of freedom, `0.0` for `n <= 1`.
    pub variance: f64,
    /// Smallest value, `0.0` when there are no values.
    pub min: f64,
    /// First quartile by [`quantile1st_sorted`], `0.0` when there are no values.
    pub lowerq: f64,
    /// Median, `0.0` when there are no values.
    pub median: f64,
    /// Third quartile by [`quantile3rd_sorted`], `0.0` when there are no values.
    pub upperq: f64,
    /// Largest value, `0.0` when there are no values.
    pub max: f64,
}

impl SummaryStatistics {
    /// Summarise a sample, sorting it in place first.
    ///
    /// The source's constructor takes its container by mutable reference and
    /// sorts it, so the `&mut` receiver is the faithful signature.
    ///
    /// An empty sample yields all-zero fields and `count == 0` — the source
    /// calls this a sanity check against a core dump. A single value yields
    /// `variance == 0.0` for the same reason: the `n - 1` divisor is zero and
    /// the source substitutes the empty case's value rather than propagating
    /// NaN.
    ///
    /// A sample holding a NaN is summarised like any other, since decision
    /// D16: the sort is the Release build's own, so the permutation it leaves
    /// is reproduced rather than guessed at, and `min`, the three quantiles and
    /// `max` are read out of it positionally exactly as
    /// `StatisticFunctions.h:952-956` reads them. That the source reads order
    /// statistics out of a range whose order the comparison did not determine
    /// is a C++ defect and is recorded as `CPP-347`; reproducing it is a
    /// decision about this port, not a defence of the source. See the module's
    /// NaN section and section 5.2 of `docs/FILE_INFO_A7_SUPPORT.md`.
    ///
    /// The order statistics are taken through the private `_of_sorted` helpers
    /// rather than the public `_sorted` entry points, because those verify a
    /// caller's sortedness claim and `std::sort`'s own output is not ascending
    /// when a NaN is in it.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] only where
    /// [`source_sort_by`] does: an
    /// introsort read outside the vector, which `<` on `f64` keys — NaN keys
    /// included — cannot provoke, and which is raised before anything is
    /// written back, so a rejected call leaves the caller's sample in its
    /// original order.
    ///
    /// Returns [`Error::InvalidRange`] only when the internal mean refuses,
    /// which the emptiness check above makes unreachable; the signature keeps
    /// the error path rather than asserting the impossibility.
    pub fn new(data: &mut [f64]) -> Result<Self> {
        let count = data.len();
        if data.is_empty() {
            return Ok(Self::default());
        }
        sort_ascending(data)?;
        let mean_value = mean(data)?;
        let variance_value = if count > 1 {
            variance_with_mean(data, mean_value)?
        } else {
            0.0
        };
        // The order statistics are read out of the sorted range positionally,
        // as `StatisticFunctions.h:952-956` reads them, through the private
        // `_of_sorted` helpers rather than the public `_sorted` entry points:
        // those verify the caller's sortedness claim, and `std::sort`'s own
        // output is not ascending when a NaN is in it.
        Ok(Self {
            count,
            mean: mean_value,
            variance: variance_value,
            min: data[0],
            lowerq: quantile1st_of_sorted(data),
            median: median_of_sorted(data),
            upperq: quantile3rd_of_sorted(data),
            max: data[count - 1],
        })
    }
}
