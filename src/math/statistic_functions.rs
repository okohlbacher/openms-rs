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

/// True when every adjacent pair is non-decreasing and no value is NaN.
fn is_ascending(values: &[f64]) -> bool {
    values.windows(2).all(|pair| pair[0] <= pair[1])
}

/// Sort in place by the IEEE-754 total order.
///
/// The source calls `std::sort`, whose strict-weak-ordering precondition a NaN
/// violates; this port orders NaN deterministically instead, so a NaN cannot
/// make the result depend on the sort implementation.
fn sort_ascending(values: &mut [f64]) {
    values.sort_by(f64::total_cmp);
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
pub fn sum(values: &[f64]) -> f64 {
    values.iter().fold(0.0, |total, value| total + value)
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
    Ok(sum(values) / values.len() as f64)
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
/// [`Error::UnsortedData`] when `values` is not ascending or contains NaN.
pub fn median_sorted(values: &[f64]) -> Result<f64> {
    check_not_empty(values)?;
    if !is_ascending(values) {
        return Err(Error::UnsortedData);
    }
    Ok(median_of_sorted(values))
}

/// Median of an ascending, non-empty slice. Callers have already checked both.
fn median_of_sorted(values: &[f64]) -> f64 {
    let size = values.len();
    if size % 2 == 0 {
        (values[size / 2 - 1] + values[size / 2]) / 2.0
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
/// Returns [`Error::InvalidRange`] for an empty range.
pub fn median(values: &mut [f64]) -> Result<f64> {
    check_not_empty(values)?;
    sort_ascending(values);
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
pub fn mad(values: &[f64], median_of_numbers: f64) -> Result<f64> {
    check_not_empty(values)?;
    preflight(values.len(), size_of::<f64>())?;
    let mut diffs = Vec::with_capacity(values.len());
    for value in values {
        diffs.push((value - median_of_numbers).abs());
    }
    sort_ascending(&mut diffs);
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
        total += (value - mean_of_numbers).abs();
    }
    total / values.len() as f64
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
/// [`Error::UnsortedData`] when `values` is not ascending or contains NaN.
pub fn quantile1st_sorted(values: &[f64]) -> Result<f64> {
    check_not_empty(values)?;
    if !is_ascending(values) {
        return Err(Error::UnsortedData);
    }
    let size = values.len();
    if size < 3 {
        return Ok(values[0]);
    }
    if size % 2 == 0 {
        return Ok(median_of_sorted(&values[..size / 2 - 1]));
    }
    Ok(median_of_sorted(&values[..size / 2]))
}

/// First quartile of a range, sorting it in place first.
///
/// The source's `quantile1st(begin, end, /*sorted=*/false)`.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range.
pub fn quantile1st(values: &mut [f64]) -> Result<f64> {
    check_not_empty(values)?;
    sort_ascending(values);
    let size = values.len();
    if size < 3 {
        return Ok(values[0]);
    }
    if size % 2 == 0 {
        return Ok(median_of_sorted(&values[..size / 2 - 1]));
    }
    Ok(median_of_sorted(&values[..size / 2]))
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
/// [`Error::UnsortedData`] when `values` is not ascending or contains NaN.
pub fn quantile3rd_sorted(values: &[f64]) -> Result<f64> {
    check_not_empty(values)?;
    if !is_ascending(values) {
        return Err(Error::UnsortedData);
    }
    let size = values.len();
    if size < 3 {
        return Ok(values[size - 1]);
    }
    Ok(median_of_sorted(&values[size / 2 + 1..]))
}

/// Third quartile of a range, sorting it in place first.
///
/// The source's `quantile3rd(begin, end, /*sorted=*/false)`.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty range.
pub fn quantile3rd(values: &mut [f64]) -> Result<f64> {
    check_not_empty(values)?;
    sort_ascending(values);
    let size = values.len();
    if size < 3 {
        return Ok(values[size - 1]);
    }
    Ok(median_of_sorted(&values[size / 2 + 1..]))
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
/// `values` is not ascending. The source states the sortedness precondition as
/// `@pre` and checks it only through `OPENMS_PRECONDITION`, which is compiled
/// out of a release build; this port always checks.
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
    Ok((1.0 - frac) * values[i] + frac * values[i + 1])
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
    sort_ascending(&mut kept);
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
    sort_ascending(&mut kept);
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
        let diff = value - mean_of_numbers;
        sum_value += diff * diff;
    }
    Ok(sum_value / (values.len() - 1) as f64)
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
    Ok(variance(values)?.sqrt())
}

/// Sample standard deviation about an explicitly supplied mean.
///
/// # Errors
///
/// As [`variance_with_mean`].
pub fn sd_with_mean(values: &[f64], mean_of_numbers: f64) -> Result<f64> {
    Ok(variance_with_mean(values, mean_of_numbers)?.sqrt())
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
        sum_value += (value_a - mean_a) * (value_b - mean_b);
    }
    Ok(sum_value / (a.len() - 1) as f64)
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
        let tmp = value_a - value_b;
        error += tmp * tmp;
    }
    Ok(error / a.len() as f64)
}

/// Root mean square error, the square root of [`mean_square_error`].
///
/// # Errors
///
/// As [`mean_square_error`].
pub fn root_mean_square_error(a: &[f64], b: &[f64]) -> Result<f64> {
    Ok(mean_square_error(a, b)?.sqrt())
}

/// Fraction of positions where two label ranges agree in sign.
///
/// A position counts as wrong when exactly one of the two values is negative;
/// zero counts as non-negative on both sides. The result is the number of
/// agreeing positions divided by `n`.
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
pub fn compute_rank(w: &mut [f64]) -> Result<()> {
    if w.is_empty() {
        return Ok(());
    }
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
/// length, and when a range exceeds [`MAX_ITEMS`] or [`MAX_BYTES`].
pub fn rank_correlation_coefficient(a: &[f64], b: &[f64]) -> Result<f64> {
    check_not_empty(a)?;
    if a.len() != b.len() {
        return Err(empty("rank correlation ranges must have the same length"));
    }
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
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] only when an internal quantile refuses,
    /// which the sort makes unreachable for a non-empty sample; the signature
    /// keeps the error path rather than asserting the impossibility.
    pub fn new(data: &mut [f64]) -> Result<Self> {
        let count = data.len();
        if data.is_empty() {
            return Ok(Self::default());
        }
        sort_ascending(data);
        let mean_value = mean(data)?;
        let variance_value = if count > 1 {
            variance_with_mean(data, mean_value)?
        } else {
            0.0
        };
        Ok(Self {
            count,
            mean: mean_value,
            variance: variance_value,
            min: data[0],
            lowerq: quantile1st_sorted(data)?,
            median: median_sorted(data)?,
            upperq: quantile3rd_sorted(data)?,
            max: data[count - 1],
        })
    }
}
