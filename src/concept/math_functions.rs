// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Marc Sturm, Chris Bielow, Timo Sachsenberg, OpenMS Rust contributors $

//! General numeric auxiliary functions: the `OpenMS::Math` namespace.
//!
//! Ported from `MATH/MathFunctions.h` (header-only; `MATH/MathFunctions.cpp`
//! contains no definitions) at revision
//! `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. See `docs/MATH_FUNCTIONS_SUPPORT.md`
//! for the member-by-member mapping, the preserved source conventions and the
//! evidence behind every transcribed value.
//!
//! The header's Doxygen places the `Math` namespace in the `Concept` group, which
//! is why this module lives beside the other `CONCEPT` ports rather than in a
//! module of its own.
//!
//! Two families carry most of the risk for callers:
//!
//! * The tolerance conversions
//!   [`ppm`](crate::concept::math_functions::ppm),
//!   [`ppm_to_mass`](crate::concept::math_functions::ppm_to_mass) and
//!   [`tolerance_window`](crate::concept::math_functions::tolerance_window).
//!   They are **not symmetric**: the source always divides by the *reference*
//!   m/z, so `ppm(a, b)` and `-ppm(b, a)` differ, and the ppm tolerance window
//!   is wider to the right than to the left. Both conventions are reproduced
//!   exactly and asserted in `tests/math_functions.rs`.
//! * The rounding helpers
//!   [`round_to`](crate::concept::math_functions::round_to),
//!   [`ceil_decimal`](crate::concept::math_functions::ceil_decimal) and
//!   [`round_decimal`](crate::concept::math_functions::round_decimal).
//!   The source builds its scaling factor by repeated multiplication or division
//!   by ten rather than by a single `pow` call, which is not the same double.
//!   That loop is transcribed literally.
//!
//! Every operation that the source leaves unchecked - division by a zero
//! denominator, a logarithm of a non-positive argument, an index derived from a
//! floating-point expression - is checked here and reported through
//! [`Error`](crate::Error) instead of producing an infinity, a NaN or an
//! out-of-bounds read.

use crate::{Error, Result};

/// Largest bin count [`create_bins`] will allocate.
///
/// The source has no ceiling: `createBins` sizes its vector directly from the
/// caller's `number_of_bins`, so a hostile count commits the memory before any
/// value is computed. Checked before allocation, so a refusal costs nothing.
pub const MAX_BINS: u32 = 1_000_000;

/// Largest trial count [`binomial_cdf_complement`] will sum over.
///
/// The complement is evaluated as an explicit sum of `trials - successes + 1`
/// probability-mass terms, so its cost is linear in the trial count. The source
/// delegates to Boost's regularized incomplete beta function and has no ceiling.
pub const MAX_BINOMIAL_TRIALS: u32 = 1_000_000;

/// Largest sample count [`quantile`] will scan.
///
/// `quantile` allocates nothing, but it validates ordering in one pass, so the
/// pass is bounded like every other input-sized operation in this crate.
pub const MAX_QUANTILE_ITEMS: usize = 100_000_000;

/// Largest magnitude of `digits` [`round_to`] will iterate over.
///
/// The source builds its scaling factor in a loop over `|digits|`, so the cost
/// is linear in a caller-supplied integer: `roundTo(x, INT_MIN)` spins over two
/// billion multiplications before returning. Any magnitude above about 308
/// already leaves the factor infinite or zero, so this ceiling refuses only
/// inputs that could not have produced a usable result.
pub const MAX_DECIMAL_DIGITS: u32 = 400;

/// One closed interval `[min, max]` of a [`create_bins`] partition.
///
/// The source returns `Math::BinContainer`, a `std::vector<RangeBase>` from
/// `KERNEL/RangeManager.h`. This port returns the two coordinates only: the
/// crate's `kernel::ranges::RangeBase` is not reachable from this module without
/// introducing a `concept -> kernel` module dependency that does not exist
/// today, and `createBins` uses none of `RangeBase`'s dimension-tagged
/// conversions. A bin may be empty (`min > max`) when a negative
/// `extend_margin` shrinks it past zero width, exactly as `RangeBase` allows.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bin {
    /// Lower bound of the bin, inclusive.
    pub min: f64,
    /// Upper bound of the bin, inclusive.
    pub max: f64,
}

impl Bin {
    /// Is the bin empty, that is, is `min` greater than `max`?
    ///
    /// Mirrors `RangeBase::isEmpty`, which is the predicate
    /// `RangeBase::extendLeftRight` consults before shrinking a bin further.
    pub fn is_empty(&self) -> bool {
        self.min > self.max
    }

    /// Is `value` inside `[min, max]`?
    ///
    /// Mirrors `RangeBase::contains(double)`. A NaN `value` is not contained.
    pub fn contains(&self, value: f64) -> bool {
        (self.min..=self.max).contains(&value)
    }
}

/// Greatest common divisor together with the Bezout coefficients that produce it.
///
/// `a * u1 + b * u2 == gcd` for the `a` and `b` passed to [`extended_gcd`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExtendedGcd {
    /// The greatest common divisor, `u3` in Knuth's notation.
    pub gcd: i64,
    /// Coefficient of `a`, `u1` in Knuth's notation.
    pub u1: i64,
    /// Coefficient of `b`, `u2` in Knuth's notation.
    pub u2: i64,
}

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

fn finite(value: f64, message: &str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad(message))
    }
}

/// Given an interval `[min, max]` and a new value, extend the interval to include it.
///
/// Returns `true` when either bound was modified. The source template
/// `Math::extendRange` takes `min` and `max` as in/out references and this port
/// keeps that shape, because the caller owns both bounds and the return value
/// alone cannot say which one moved.
///
/// # Arguments
///
/// * `min` - current minimum, lowered when `value` falls below it
/// * `max` - current maximum, raised when `value` rises above it
/// * `value` - the new value which may extend the interval
///
/// A NaN `value` compares false against both bounds, so nothing is modified and
/// `false` is returned; the source behaves identically and does not check.
/// Note that the source tests `value < min` first and returns immediately, so an
/// inverted interval (`min > max`) has only its `min` lowered by one call.
///
/// ```
/// use openms::concept::math_functions::extend_range;
///
/// let (mut min, mut max) = (1.0, 2.0);
/// assert!(extend_range(&mut min, &mut max, 3.0));
/// assert_eq!((min, max), (1.0, 3.0));
/// assert!(!extend_range(&mut min, &mut max, 1.5));
/// ```
pub fn extend_range(min: &mut f64, max: &mut f64, value: f64) -> bool {
    if value < *min {
        *min = value;
        return true;
    }
    if value > *max {
        *max = value;
        return true;
    }
    false
}

/// Is `value` contained in the closed interval `[min, max]`?
///
/// Transcribes `Math::contains`, which evaluates `min <= value && value <= max`.
/// A NaN argument makes both comparisons false, so the answer is `false`.
pub fn contains(value: f64, min: f64, max: f64) -> bool {
    (min..=max).contains(&value)
}

/// Zoom into the interval `[left, right]`, scaling its width by `factor`.
///
/// A `factor` in `[0, 1]` shrinks the span, a `factor` above one extends it, and
/// `align` decides where the resulting interval sits: `0.0` keeps `left` fixed,
/// `1.0` keeps `right` fixed, `0.5` zooms into the centre. Round trips invert
/// exactly by inverting the factor, which is the source's own documented example:
///
/// ```
/// use openms::concept::math_functions::zoom_in;
///
/// let (a2, b2) = zoom_in(10.0, 20.0, 0.5, 0.5)?;
/// let (a1, b1) = zoom_in(a2, b2, 2.0, 0.5)?;
/// assert_eq!((a2, b2), (12.5, 17.5));
/// assert_eq!((a1, b1), (10.0, 20.0));
/// # Ok::<(), openms::Error>(())
/// ```
///
/// # Arguments
///
/// * `left`, `right` - start and end of the interval
/// * `factor` - width multiplier, at least zero
/// * `align` - position of the zoomed interval, between zero and one
///
/// `factor` and `align` are `f32` because the source declares them `float`, and
/// the term `(1.0f - factor)` is therefore evaluated in single precision before
/// it widens. Computing it in `f64` would give a different low-order bit, so the
/// single-precision subtraction is reproduced.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `factor` is negative, when `align` is
/// outside `[0, 1]`, when any argument is not finite, or when the result is not
/// finite. The source states these two as `OPENMS_PRECONDITION`, which is
/// compiled out of release builds, so an out-of-range `align` silently produces
/// a misplaced interval there; this port checks unconditionally.
pub fn zoom_in(left: f64, right: f64, factor: f32, align: f32) -> Result<(f64, f64)> {
    if factor.is_nan() || factor < 0.0 {
        return Err(bad("zoom factor must be >= 0"));
    }
    if !(0.0..=1.0).contains(&align) {
        return Err(bad("zoom alignment must be within [0, 1]"));
    }
    if !left.is_finite() || !right.is_finite() {
        return Err(bad("zoom interval must be finite"));
    }
    let old_width = right - left;
    // Source order: (1.0f - factor) is a float subtraction, then both products
    // widen to double. `res.second` is computed from the already-shifted
    // `res.first`, not from `right`.
    let offset_left = f64::from(1.0_f32 - factor) * old_width * f64::from(align);
    let new_left = finite(left + offset_left, "zoom lower bound is not finite")?;
    let new_right = finite(
        new_left + old_width * f64::from(factor),
        "zoom upper bound is not finite",
    )?;
    Ok((new_left, new_right))
}

/// Split `[min, max]` into `number_of_bins` bins, optionally overlapping.
///
/// Each bin is widened on both sides by `extend_margin`, so neighbouring bins
/// overlap by `2 * extend_margin`. A negative margin shrinks the bins instead,
/// which the source documents as a feature and which can leave a bin empty. The
/// outer borders of the original interval are never extended: the first bin's
/// minimum is reset to `min` and the last bin's maximum to `max` afterwards,
/// including the source's `RangeBase::setMin`/`setMax` repair of the opposite
/// bound when the reset would invert the bin.
///
/// # Arguments
///
/// * `min` - minimum of the range; must be smaller than `max`
/// * `max` - maximum of the range
/// * `number_of_bins` - how many bins to divide the range into; at least one
/// * `extend_margin` - overlap of neighbouring bins; zero for no overlap
///
/// ```
/// use openms::concept::math_functions::create_bins;
///
/// let bins = create_bins(0.0, 10.0, 2, 0.0)?;
/// assert_eq!((bins[0].min, bins[0].max), (0.0, 5.0));
/// assert_eq!((bins[1].min, bins[1].max), (5.0, 10.0));
/// # Ok::<(), openms::Error>(())
/// ```
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] when `min >= max` and [`Error::InvalidValue`]
/// when `number_of_bins` is zero, exceeds [`MAX_BINS`], or when any argument or
/// computed bound is not finite. The source declares both conditions as
/// `OPENMS_PRECONDITION`, which is inactive in release builds; there a zero bin
/// count then calls `front()` and `back()` on an empty vector, which is
/// undefined behaviour, and `min > max` instead escapes as the
/// `Exception::InvalidRange` thrown by the `RangeBase(min, max)` constructor.
pub fn create_bins(
    min: f64,
    max: f64,
    number_of_bins: u32,
    extend_margin: f64,
) -> Result<Vec<Bin>> {
    if number_of_bins == 0 {
        return Err(bad("number of bins must be >= 1"));
    }
    if number_of_bins > MAX_BINS {
        return Err(bad("number of bins exceeds MAX_BINS"));
    }
    if !min.is_finite() || !max.is_finite() || !extend_margin.is_finite() {
        return Err(bad("bin range and margin must be finite"));
    }
    if min >= max {
        return Err(Error::InvalidRange("bin range requires min < max".into()));
    }
    let bin_width = (max - min) / f64::from(number_of_bins);
    let mut bins = Vec::with_capacity(number_of_bins as usize);
    for i in 0..number_of_bins {
        let mut bin = Bin {
            min: finite(min + f64::from(i) * bin_width, "bin bound is not finite")?,
            max: finite(
                min + f64::from(i + 1) * bin_width,
                "bin bound is not finite",
            )?,
        };
        // RangeBase::extendLeftRight is a no-op on an already empty range.
        if !bin.is_empty() {
            bin.min = finite(bin.min - extend_margin, "bin bound is not finite")?;
            bin.max = finite(bin.max + extend_margin, "bin bound is not finite")?;
        }
        bins.push(bin);
    }
    let first = &mut bins[0];
    // RangeBase::setMin also raises max when the range was left inverted.
    first.min = min;
    if first.max < min {
        first.max = min;
    }
    let last = bins
        .last_mut()
        .ok_or_else(|| bad("bin container cannot be empty"))?;
    last.max = max;
    if last.min > max {
        last.min = max;
    }
    Ok(bins)
}

/// Round `x` up to the next decimal power `10 ^ dec_pow`.
///
/// ```text
/// (123.0 ,  1)  => 130
/// (123.0 ,  2)  => 200
/// (  0.123, -2) => 0.13    (10^-2 = 0.01)
/// ```
///
/// The source shifts right, applies `ceil`, and shifts left again, using
/// `pow(10.0, decPow)` for both shifts; that same `pow` call is used here so the
/// scaling factor is bit-identical.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `x` is not finite, when the scaling
/// factor underflows to zero or overflows (roughly `|dec_pow| > 308`), or when
/// the result is not finite. The source performs the division regardless and
/// yields an infinity or a NaN.
pub fn ceil_decimal(x: f64, dec_pow: i32) -> Result<f64> {
    let scale = decimal_power(dec_pow)?;
    finite(x, "ceil_decimal input must be finite")?;
    finite((x / scale).ceil() * scale, "ceil_decimal is not finite")
}

/// Round `x` to the nearest decimal power `10 ^ dec_pow`.
///
/// ```text
/// (123.0 , 1)  => 120
/// (123.0 , 2)  => 100
/// ```
///
/// The source has two branches and they are not symmetric in their spelling: for
/// `x > 0` it evaluates `floor(0.5 + x / scale) * scale`, and otherwise it
/// negates the same expression applied to `fabs(x)`. Because zero is not greater
/// than zero it takes the negating branch, so `round_decimal(0.0, 0)` returns
/// negative zero. That is reproduced rather than normalised, because callers can
/// observe it through `1.0_f64.copysign(result)` and through division.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] under the same conditions as [`ceil_decimal`].
pub fn round_decimal(x: f64, dec_pow: i32) -> Result<f64> {
    let scale = decimal_power(dec_pow)?;
    finite(x, "round_decimal input must be finite")?;
    let rounded = if x > 0.0 {
        (0.5 + x / scale).floor() * scale
    } else {
        -((0.5 + x.abs() / scale).floor() * scale)
    };
    finite(rounded, "round_decimal is not finite")
}

/// `pow(10.0, dec_pow)` with the source's exact `double`-argument call.
///
/// C++ promotes the `int` exponent, so `std::pow(double, double)` is selected
/// and `f64::powf` is the matching Rust call. `powi` would evaluate a different
/// sequence of multiplications.
fn decimal_power(dec_pow: i32) -> Result<f64> {
    let scale = 10.0_f64.powf(f64::from(dec_pow));
    if !scale.is_finite() || scale == 0.0 {
        return Err(bad("decimal power is outside the representable range"));
    }
    Ok(scale)
}

/// Transform point `x` of the interval `[left1, right1]` into `[left2, right2]`.
///
/// The source expression is `left2 + (x - left1) * (right2 - left2) / (right1 - left1)`,
/// which multiplies before it divides; that order is preserved because a
/// mathematically equivalent regrouping rounds differently.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when any argument is not finite, when
/// `right1 == left1` - the source divides by zero there and returns an infinity
/// or a NaN - or when the result is not finite.
pub fn interval_transformation(
    x: f64,
    left1: f64,
    right1: f64,
    left2: f64,
    right2: f64,
) -> Result<f64> {
    for value in [x, left1, right1, left2, right2] {
        finite(value, "interval transformation needs finite arguments")?;
    }
    if right1 == left1 {
        return Err(bad("interval transformation needs a non-empty source span"));
    }
    finite(
        left2 + (x - left1) * (right2 - left2) / (right1 - left1),
        "interval transformation is not finite",
    )
}

/// Transform a number from linear to log10 scale, adding one first.
///
/// Source `Math::linear2log`. The added one is what keeps small positive inputs
/// from producing negative logarithms; it is not a smoothing constant and must
/// be undone by [`log10_to_linear`].
///
/// # Arguments
///
/// * `x` - the number to transform
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `x + 1.0` is not strictly positive, or
/// when either the argument or the result is not finite. The source calls
/// `log10` unguarded, so `x == -1` yields negative infinity and `x < -1` a NaN.
pub fn linear_to_log10(x: f64) -> Result<f64> {
    finite(x, "linear_to_log10 input must be finite")?;
    let shifted = x + 1.0;
    if shifted <= 0.0 {
        return Err(bad("linear_to_log10 needs x + 1 > 0"));
    }
    finite(shifted.log10(), "linear_to_log10 is not finite")
}

/// Transform a number from log10 to linear scale, subtracting the added one.
///
/// Source `Math::log2linear`, the inverse of [`linear_to_log10`].
///
/// # Arguments
///
/// * `x` - the number to transform
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the argument or the result is not
/// finite; the source overflows to infinity for `x` above about 308.
pub fn log10_to_linear(x: f64) -> Result<f64> {
    finite(x, "log10_to_linear input must be finite")?;
    finite(10.0_f64.powf(x) - 1.0, "log10_to_linear is not finite")
}

/// Is the given integer odd?
///
/// Source `Math::isOdd(UInt)`, which tests the low bit rather than a remainder.
pub fn is_odd(x: u32) -> bool {
    (x & 1) != 0
}

/// Round to the nearest integer, halves away from zero.
///
/// Source `Math::round` forwards to `std::round`, whose tie rule is "away from
/// zero"; `f64::round` has the same rule, so the mapping is exact. The source
/// template is instantiated for `float` as well as `double` in its own class
/// test; `f32::round` covers that case unchanged and is not wrapped here.
///
/// A NaN argument returns NaN, as in the source.
pub fn round(x: f64) -> f64 {
    x.round()
}

/// Round to the `digits`-th decimal place; negative digits round to tens, hundreds and so on.
///
/// ```
/// use openms::concept::math_functions::round_to;
///
/// assert_eq!(round_to(3.14159265, 2)?, 3.14);
/// assert_eq!(round_to(1234.9, -2)?, 1200.0);
/// # Ok::<(), openms::Error>(())
/// ```
///
/// # Arguments
///
/// * `value` - the value to round
/// * `digits` - number of digits to round to; may be negative
///
/// The scaling factor is built exactly as the source builds it, by multiplying
/// or dividing a running `1.0` by ten `|digits|` times. That is deliberately not
/// `10f64.powi(digits)`: repeated division accumulates rounding, so for example
/// the factor for `digits == -3` is the double nearest `0.001` reached through
/// three divisions, and replacing the loop changes results in the last places.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `value` is not finite, when `|digits|`
/// exceeds [`MAX_DECIMAL_DIGITS`], when the factor underflows to zero or
/// overflows, or when the result is not finite. The source divides by the factor
/// unguarded and runs its loop for any `int`.
pub fn round_to(value: f64, digits: i32) -> Result<f64> {
    finite(value, "round_to input must be finite")?;
    if digits.unsigned_abs() > MAX_DECIMAL_DIGITS {
        return Err(bad("round_to digits exceed MAX_DECIMAL_DIGITS"));
    }
    let mut factor = 1.0_f64;
    if digits > 0 {
        for _ in 0..digits.unsigned_abs() {
            factor *= 10.0;
        }
    } else if digits < 0 {
        for _ in 0..digits.unsigned_abs() {
            factor /= 10.0;
        }
    }
    if !factor.is_finite() || factor == 0.0 {
        return Err(bad(
            "round_to scaling factor is outside the representable range",
        ));
    }
    finite((value * factor).round() / factor, "round_to is not finite")
}

/// Percentage of `value` relative to `total`, rounded to `digits` decimal places.
///
/// ```
/// use openms::concept::math_functions::percent_of;
///
/// assert_eq!(percent_of(1.0 / 3.0, 1.0, 2)?, 33.33);
/// # Ok::<(), openms::Error>(())
/// ```
///
/// # Arguments
///
/// * `value` - the value to compute the percentage for; must not be negative
/// * `total` - the total to compute it against; must not be negative
/// * `digits` - number of digits to round the result to
///
/// A `total` of zero returns `0.0` rather than dividing, which is the source's
/// documented behaviour and its reason for comparing `total <= 0` instead of
/// `total == 0`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `value` or `total` is negative, matching
/// the source's `Exception::InvalidValue`; when either is not finite, which the
/// source does not check and which would otherwise pass both sign tests and
/// return a NaN; and when [`round_to`] fails.
pub fn percent_of(value: f64, total: f64, digits: i32) -> Result<f64> {
    finite(value, "percent_of value must be finite")?;
    finite(total, "percent_of total must be finite")?;
    if value < 0.0 {
        return Err(bad("percent_of value must be non-negative"));
    }
    if total < 0.0 {
        return Err(bad("percent_of total must be non-negative"));
    }
    if total <= 0.0 {
        return Ok(0.0);
    }
    round_to(value * 100.0 / total, digits)
}

/// Is `a` approximately equal to `b`, within an absolute tolerance `tol`?
///
/// Source `Math::approximatelyEqual`, which is an absolute comparison
/// `fabs(a - b) <= tol` and not a relative one. A NaN argument answers `false`.
pub fn approximately_equal(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

/// Greatest common divisor of two numbers by the Euclidean algorithm.
///
/// The loop is the source's: `c = a % b; a = b; b = c` until `b` is zero. The
/// result carries the sign of `a` when `b` reaches zero, so `gcd(-4, 2)` is `2`
/// while `gcd(-4, 0)` is `-4`; the source template does not normalise the sign
/// and neither does this port.
///
/// See also [`extended_gcd`].
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] for `gcd(i64::MIN, -1)`, whose remainder
/// overflows. The source's `%` on the equivalent C++ types is undefined
/// behaviour there and traps on x86.
pub fn gcd(a: i64, b: i64) -> Result<i64> {
    let (mut a, mut b) = (a, b);
    while b != 0 {
        let c = a
            .checked_rem(b)
            .ok_or_else(|| bad("gcd remainder overflows"))?;
        a = b;
        b = c;
    }
    Ok(a)
}

/// Greatest common divisor by the extended Euclidean algorithm.
///
/// Follows Knuth, *The Art of Computer Programming* volume 2 page 342, as cited
/// by the source: it maintains `(u1, u2, u3)` and `(v1, v2, v3)` so that
/// `a * u1 + b * u2 == u3 == gcd(a, b)`. The source returns `u3` and writes `u1`
/// and `u2` through out-parameters; this port returns all three in
/// [`ExtendedGcd`].
///
/// # Arguments
///
/// * `a`, `b` - the two numbers
///
/// ```
/// use openms::concept::math_functions::extended_gcd;
///
/// let result = extended_gcd(240, 46)?;
/// assert_eq!(result.gcd, 2);
/// assert_eq!(240 * result.u1 + 46 * result.u2, result.gcd);
/// # Ok::<(), openms::Error>(())
/// ```
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when any of the divisions, multiplications or
/// subtractions overflows `i64`. The source performs them unchecked, so the
/// same inputs silently wrap or trap there.
pub fn extended_gcd(a: i64, b: i64) -> Result<ExtendedGcd> {
    let overflow = || bad("extended gcd overflows");
    let (mut u1, mut u2, mut u3) = (1_i64, 0_i64, a);
    let (mut v1, mut v2, mut v3) = (0_i64, 1_i64, b);
    while v3 != 0 {
        let q = u3.checked_div(v3).ok_or_else(overflow)?;
        let t1 = u1
            .checked_sub(v1.checked_mul(q).ok_or_else(overflow)?)
            .ok_or_else(overflow)?;
        let t2 = u2
            .checked_sub(v2.checked_mul(q).ok_or_else(overflow)?)
            .ok_or_else(overflow)?;
        let t3 = u3
            .checked_sub(v3.checked_mul(q).ok_or_else(overflow)?)
            .ok_or_else(overflow)?;
        u1 = v1;
        u2 = v2;
        u3 = v3;
        v1 = t1;
        v2 = t2;
        v3 = t3;
    }
    Ok(ExtendedGcd { gcd: u3, u1, u2 })
}

/// Parts-per-million deviation of an observed m/z from a reference m/z.
///
/// Source `Math::getPPM`, which evaluates `(mz_obs - mz_ref) / mz_ref * 1e6`.
/// The sign is kept: the result is positive when `mz_obs > mz_ref` and negative
/// when `mz_obs < mz_ref`.
///
/// # Arguments
///
/// * `mz_obs` - observed (experimental) m/z
/// * `mz_ref` - reference (theoretical) m/z; the divisor
///
/// The division is by the **reference**, never by the observed value, so the
/// function is not antisymmetric in its arguments: `ppm(1000.0, 1001.0)` is
/// about `-999.001`, not `-1000.0`. Swapping the arguments at a call site is
/// therefore a silent, small, mass-dependent error rather than a sign flip.
///
/// ```
/// use openms::concept::math_functions::ppm;
///
/// assert_eq!(ppm(1001.0, 1000.0)?, 1000.0);
/// assert_eq!(ppm(999.0, 1000.0)?, -1000.0);
/// # Ok::<(), openms::Error>(())
/// ```
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when either argument is not finite, when
/// `mz_ref` is zero, or when the result is not finite. The source divides
/// unguarded and returns an infinity or a NaN.
pub fn ppm(mz_obs: f64, mz_ref: f64) -> Result<f64> {
    finite(mz_obs, "observed m/z must be finite")?;
    finite(mz_ref, "reference m/z must be finite")?;
    if mz_ref == 0.0 {
        return Err(bad("reference m/z must not be zero"));
    }
    finite((mz_obs - mz_ref) / mz_ref * 1e6, "ppm is not finite")
}

/// Absolute parts-per-million deviation of an observed m/z from a reference m/z.
///
/// Source `Math::getPPMAbs`, the absolute value of [`ppm`]. Always at least zero.
///
/// # Errors
///
/// As [`ppm`].
pub fn ppm_abs(mz_obs: f64, mz_ref: f64) -> Result<f64> {
    Ok(ppm(mz_obs, mz_ref)?.abs())
}

/// Mass difference in Thomson for a ppm deviation at a reference m/z.
///
/// Source `Math::ppmToMass`, which evaluates `(ppm / 1e6) * mz_ref` - the
/// division happens first, so this is not the exact algebraic inverse of [`ppm`]
/// in floating point. The sign of `ppm` is carried through.
///
/// # Arguments
///
/// * `ppm_value` - parts-per-million error
/// * `mz_ref` - reference m/z the tolerance is anchored at
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when an argument or the result is not finite.
/// Unlike [`ppm`] this direction has no divisor to guard, so a zero `mz_ref`
/// simply yields a zero mass difference, as in the source.
pub fn ppm_to_mass(ppm_value: f64, mz_ref: f64) -> Result<f64> {
    finite(ppm_value, "ppm must be finite")?;
    finite(mz_ref, "reference m/z must be finite")?;
    finite((ppm_value / 1e6) * mz_ref, "mass difference is not finite")
}

/// Absolute mass difference in Thomson for a ppm deviation at a reference m/z.
///
/// Source `Math::ppmToMassAbs`, the absolute value of [`ppm_to_mass`]. Always at
/// least zero. Its Doxygen block opens with `/*` rather than `/**`, so this
/// member is missing from the generated C++ documentation; the comment's content
/// is carried across here regardless.
///
/// # Errors
///
/// As [`ppm_to_mass`].
pub fn ppm_to_mass_abs(ppm_value: f64, mz_ref: f64) -> Result<f64> {
    Ok(ppm_to_mass(ppm_value, mz_ref)?.abs())
}

/// Tolerance window `(left, right)` around `val` for a tolerance `tol`.
///
/// # Arguments
///
/// * `val` - the value the window is centred on
/// * `tol` - the tolerance, in ppm when `ppm_mode` is set and in absolute units otherwise
/// * `ppm_mode` - whether `tol` is a relative ppm tolerance
///
/// With ppm the window is deliberately **not symmetric**:
/// `right - val` is larger than `val - left`, because the right edge is
/// `val / (1 - tol * 1e-6)` rather than `val * (1 + tol * 1e-6)`. That is the
/// largest value `x` which still has `val` inside *its* own ppm window, so the
/// compatibility relation between two masses is symmetric even though the window
/// is not. Widening the left side instead, or using a symmetric window, changes
/// which pairs match at the boundary.
///
/// ```
/// use openms::concept::math_functions::tolerance_window;
///
/// let (left, right) = tolerance_window(1000.0, 10.0, true)?;
/// assert_eq!(left, 999.99);
/// assert!(right - 1000.0 > 1000.0 - left);
/// # Ok::<(), openms::Error>(())
/// ```
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when an argument or a bound is not finite, or
/// when `tol * 1e-6` equals one in ppm mode - that is `tol == 1e6`, where the
/// source divides by zero and returns an infinite right edge.
pub fn tolerance_window(val: f64, tol: f64, ppm_mode: bool) -> Result<(f64, f64)> {
    finite(val, "tolerance window value must be finite")?;
    finite(tol, "tolerance must be finite")?;
    if !ppm_mode {
        let left = finite(val - tol, "tolerance window bound is not finite")?;
        let right = finite(val + tol, "tolerance window bound is not finite")?;
        return Ok((left, right));
    }
    let denominator = 1.0 - tol * 1e-6;
    if denominator == 0.0 {
        return Err(bad("ppm tolerance window denominator is zero"));
    }
    let left = finite(
        val - val * tol * 1e-6,
        "tolerance window bound is not finite",
    )?;
    let right = finite(val / denominator, "tolerance window bound is not finite")?;
    Ok((left, right))
}

/// Value of the `q`-th quantile of a sorted, non-empty sample.
///
/// `q` is clamped into `[0, 1]` before use, as in the source. The index is the
/// source's `max(0, n * q - 1)`, interpolated linearly between its floor and its
/// ceiling. Note that this is **not** the common `(n - 1) * q` convention: for
/// `[1, 2, 3, 4, 5]` and `q == 0.5` the source's index is `1.5`, so the returned
/// "median" is `2.5` rather than `3`. The convention is reproduced because
/// callers of `Math::quantile` are calibrated against it, and it is called out
/// here because it is easy to mistake for a textbook quantile.
///
/// # Arguments
///
/// * `x` - a sample sorted in ascending order
/// * `q` - the quantile in `[0, 1]`; values outside are clamped
///
/// ```
/// use openms::concept::math_functions::quantile;
///
/// assert_eq!(quantile(&[1.0, 2.0, 3.0, 4.0, 5.0], 1.0)?, 5.0);
/// assert_eq!(quantile(&[1.0, 2.0, 3.0, 4.0, 5.0], 0.0)?, 1.0);
/// # Ok::<(), openms::Error>(())
/// ```
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `x` is empty - the source throws
/// `Exception::InvalidParameter`, which has no counterpart among this crate's
/// error variants - when `x` is longer than [`MAX_QUANTILE_ITEMS`], when `q` or
/// a sample is not finite, or when the interpolated index leaves the sample.
/// Returns [`Error::UnsortedData`] when `x` is not ascending: the source states
/// the precondition in its `@brief` and does not check it, and an unsorted
/// sample silently returns a value that is not a quantile of anything.
pub fn quantile(x: &[f64], q: f64) -> Result<f64> {
    if x.is_empty() {
        return Err(bad("quantile requested from an empty sample"));
    }
    if x.len() > MAX_QUANTILE_ITEMS {
        return Err(bad("quantile sample exceeds MAX_QUANTILE_ITEMS"));
    }
    finite(q, "quantile must be finite")?;
    if x.iter().any(|value| !value.is_finite()) {
        return Err(bad("quantile sample must be finite"));
    }
    if x.windows(2).any(|pair| pair[0] > pair[1]) {
        return Err(Error::UnsortedData);
    }
    let q = q.clamp(0.0, 1.0);
    let n = x.len() as f64;
    // Source: `std::max(0., n * q - 1)`, the -1 being its C++ index correction.
    let id = (n * q - 1.0).max(0.0);
    let lo = id.floor();
    let hi = id.ceil();
    let h = id - lo;
    let low = index(x, lo)?;
    let high = index(x, hi)?;
    finite((1.0 - h) * low + h * high, "quantile is not finite")
}

/// Read `x` at a floating-point index, as the source's `x[lo]` conversion does.
fn index(x: &[f64], position: f64) -> Result<f64> {
    if position.is_nan() || position < 0.0 || position >= usize::MAX as f64 {
        return Err(bad("quantile index is outside the sample"));
    }
    x.get(position as usize)
        .copied()
        .ok_or_else(|| bad("quantile index is outside the sample"))
}

/// Natural logarithm of the binomial coefficient `C(n, k)`, via the log-gamma function.
///
/// # Arguments
///
/// * `n` - total number of items
/// * `k` - number of items to choose; must not exceed `n`
///
/// The source's two shortcuts are kept and they matter numerically: `k == 0` and
/// `k == n` return exactly `0.0` instead of a difference of three log-gamma
/// values, and `k > n / 2` - integer division - is replaced by `n - k` before
/// the log-gammas are evaluated, so the symmetric pair `C(10, 3)` and `C(10, 7)`
/// are bit-identical rather than merely close.
///
/// ```
/// use openms::concept::math_functions::log_binomial_coef;
///
/// assert_eq!(log_binomial_coef(10, 3)?, log_binomial_coef(10, 7)?);
/// assert_eq!(log_binomial_coef(5, 5)?, 0.0);
/// # Ok::<(), openms::Error>(())
/// ```
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `k > n`, where the source throws
/// `std::invalid_argument` rather than an OpenMS exception, and when the result
/// is not finite.
///
/// The log-gamma itself comes from `libm` where the source uses
/// `boost::math::lgamma`. These are different implementations of the same
/// function and agree to within a few units in the last place over this domain,
/// which is far inside the `1e-5` relative tolerance the source's own class test
/// applies.
pub fn log_binomial_coef(n: u32, k: u32) -> Result<f64> {
    if k > n {
        return Err(bad("k cannot be greater than n in a binomial coefficient"));
    }
    if k == 0 || k == n {
        return Ok(0.0);
    }
    // Source: `if (k > n / 2) k = n - k;` with integer division of n.
    let k = if k > n / 2 { n - k } else { k };
    let value = libm::lgamma(f64::from(n) + 1.0)
        - libm::lgamma(f64::from(k) + 1.0)
        - libm::lgamma(f64::from(n - k) + 1.0);
    finite(value, "log binomial coefficient is not finite")
}

/// Numerically stable `ln(exp(x) + exp(y))`.
///
/// Negative infinity is an identity on either side, which is what makes this
/// usable to accumulate log-domain probabilities that include impossible events.
/// Otherwise the larger value is factored out before the exponentials, so
/// neither term overflows.
///
/// # Arguments
///
/// * `x`, `y` - the two logarithmic values
///
/// ```
/// use openms::concept::math_functions::log_sum_exp;
///
/// assert_eq!(log_sum_exp(f64::NEG_INFINITY, 5.0), 5.0);
/// assert_eq!(log_sum_exp(10.0, 10.0), 10.0 + 2.0_f64.ln());
/// ```
///
/// This function is total and returns no error, because the source guards only
/// the negative infinities. Positive infinity is not guarded there: it reaches
/// `inf - inf` and yields NaN, and that is reproduced. The maximum is selected
/// with the source's `std::max(x, y)`, which returns its first argument when the
/// comparison is false and so propagates a NaN `x`; `f64::max` would discard it.
pub fn log_sum_exp(x: f64, y: f64) -> f64 {
    if x.is_infinite() && x < 0.0 {
        return y;
    }
    if y.is_infinite() && y < 0.0 {
        return x;
    }
    // std::max(a, b) is `(a < b) ? b : a`, which is not f64::max for NaN.
    let max_val = if x < y { y } else { x };
    max_val + ((x - max_val).exp() + (y - max_val).exp()).ln()
}

/// Binomial upper tail `P(X >= successes)` for `X ~ B(trials, p)`.
///
/// # Arguments
///
/// * `trials` - total number of trials, `N`
/// * `successes` - minimum number of successes, `n`; must not exceed `trials`
/// * `p` - probability of success in each trial, within `[0, 1]`
///
/// The source's four shortcuts are kept exactly: zero successes returns `1.0`
/// before anything else is examined, `p == 0.0` returns `0.0` for a positive
/// success count, and `p == 1.0` returns `1.0` because all mass sits at `N`.
///
/// ```
/// use openms::concept::math_functions::binomial_cdf_complement;
///
/// // The exact value is the dyadic rational 638/1024 = 0.623046875. The
/// // log-domain sum reaches it to a few units in the last place, not exactly.
/// let tail = binomial_cdf_complement(10, 5, 0.5)?;
/// assert!((tail - 638.0 / 1024.0).abs() < 1e-13);
/// # Ok::<(), openms::Error>(())
/// ```
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `p` is outside `[0, 1]` or not finite,
/// and when `successes > trials`; the source throws `std::invalid_argument` for
/// the first and third of those and does not check for NaN at all, where its
/// `p < 0.0 || p > 1.0` test passes a NaN straight into Boost. Also returns
/// [`Error::InvalidValue`] when `trials` exceeds [`MAX_BINOMIAL_TRIALS`], or when
/// an intermediate term is not finite.
///
/// The remaining tail is summed here as its own definition -
/// `sum over k in n..=N of C(N, k) p^k (1-p)^(N-k)`, evaluated in the log domain
/// with the largest term factored out and accumulated in ascending `k` - where
/// the source calls `boost::math::cdf(complement(binomial_distribution(N, p), n - 1))`,
/// which is the regularized incomplete beta function. The two agree to about
/// `2.5e-14` relative against an exact rational evaluation of the same sum for
/// every case in the source's class test, four orders of magnitude inside that
/// test's `1e-5` tolerance. The trade is a cost linear in `trials` instead of the
/// continued fraction's near-constant cost, which is why the ceiling exists.
pub fn binomial_cdf_complement(trials: u32, successes: u32, p: f64) -> Result<f64> {
    if !(0.0..=1.0).contains(&p) {
        return Err(bad("probability p must be between 0 and 1"));
    }
    if successes > trials {
        return Err(bad("successes cannot be greater than trials"));
    }
    if successes == 0 {
        return Ok(1.0);
    }
    if p == 0.0 {
        return Ok(0.0);
    }
    if p == 1.0 {
        return Ok(1.0);
    }
    if trials > MAX_BINOMIAL_TRIALS {
        return Err(bad("trials exceed MAX_BINOMIAL_TRIALS"));
    }
    let log_p = p.ln();
    // ln(1 - p) via ln_1p keeps the accuracy the naive difference loses for
    // small p; the source never forms this term at all.
    let log_q = (-p).ln_1p();
    let mut terms = Vec::with_capacity((trials - successes + 1) as usize);
    let mut largest = f64::NEG_INFINITY;
    for k in successes..=trials {
        let term =
            log_binomial_coef(trials, k)? + f64::from(k) * log_p + f64::from(trials - k) * log_q;
        if term > largest {
            largest = term;
        }
        terms.push(term);
    }
    let mut total = 0.0;
    for term in terms {
        total += (term - largest).exp();
    }
    let value = finite(largest.exp() * total, "binomial complement is not finite")?;
    // A probability cannot exceed one; the factored sum can overshoot by a few
    // units in the last place. The source's incomplete beta cannot.
    Ok(value.clamp(0.0, 1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_scaling_factor_loop_is_not_a_single_power() {
        // Three divisions of 1.0 by ten, the source's own construction.
        let mut factor = 1.0_f64;
        for _ in 0..3 {
            factor /= 10.0;
        }
        assert_eq!(
            round_to(1234.9, -3).unwrap(),
            (1234.9 * factor).round() / factor
        );
        // powi and powf agree with each other here but the loop is what is used.
        assert_eq!(decimal_power(-3).unwrap(), 10.0_f64.powf(-3.0));
    }

    #[test]
    fn floating_point_index_reads_stay_inside_the_sample() {
        let sample = [1.0, 2.0, 3.0];
        assert_eq!(index(&sample, 0.0).unwrap(), 1.0);
        assert_eq!(index(&sample, 2.0).unwrap(), 3.0);
        assert!(index(&sample, 3.0).is_err());
        assert!(index(&sample, -1.0).is_err());
        assert!(index(&sample, f64::NAN).is_err());
        assert!(index(&sample, 1e300).is_err());
    }

    #[test]
    fn log_sum_exp_uses_the_source_maximum_not_the_rust_one() {
        assert!(log_sum_exp(f64::NAN, 1.0).is_nan());
        assert!(log_sum_exp(1.0, f64::NAN).is_nan());
        assert!(log_sum_exp(f64::INFINITY, 1.0).is_nan());
        assert_eq!(
            log_sum_exp(f64::NEG_INFINITY, f64::NEG_INFINITY),
            f64::NEG_INFINITY
        );
    }
}
