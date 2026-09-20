// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Emulation of the instructions the Linux x86_64 Release build of
//! `libOpenMS.so` (`openms4-release-bc9cc12-c19e494-174b576`, GCC 14.4, `-O3
//! -mssse3 -ffp-contract=off`) emits, where IEEE 754 alone does not fix the
//! result.
//!
//! The arithmetic helpers return the IEEE result whenever it is not NaN, so
//! they never change a number. For NaN they follow the SSE2 rule (Intel SDM
//! vol. 1, "Rules for handling NaNs"): the first NaN operand of the
//! instruction, quieted, or the default NaN `0xfff8000000000000` when the
//! operation itself is invalid. Rust's own arithmetic leaves those bits to the
//! host, and an arm64 host produces the positive default NaN.
//!
//! The rules are the instruction set's, not any one header's: the picked
//! feature finder's `intensityScore_` was the first caller to need them, and
//! [`crate::math::statistic_functions`] is the second — a `-nan` the Release
//! build prints out of `Math::SummaryStatistics` has the same provenance as a
//! `-nan` it stores in a feature's `FWHM`. The module was promoted out of
//! `analysis::feature_finder_picked::scoring` for that second caller; not one
//! line of the emulation changed in the move.

/// SSE2's default ("real indefinite") NaN.
pub(crate) const DEFAULT_NAN: f64 = f64::from_bits(0xfff8_0000_0000_0000);

fn quiet(x: f64) -> f64 {
    f64::from_bits(x.to_bits() | 0x0008_0000_0000_0000)
}

/// The NaN rule of a two-operand SSE2 instruction whose destination
/// register holds `first`.
fn nan_rule(first: f64, second: f64, result: f64) -> f64 {
    if !result.is_nan() {
        result
    } else if first.is_nan() {
        quiet(first)
    } else if second.is_nan() {
        quiet(second)
    } else {
        DEFAULT_NAN
    }
}

/// `addsd`: `first + second`.
pub(crate) fn add(first: f64, second: f64) -> f64 {
    nan_rule(first, second, first + second)
}

/// `subsd`: `first - second`.
pub(crate) fn sub(first: f64, second: f64) -> f64 {
    nan_rule(first, second, first - second)
}

/// `mulsd`: `first * second`.
pub(crate) fn mul(first: f64, second: f64) -> f64 {
    nan_rule(first, second, first * second)
}

/// `divsd`: `first / second`.
pub(crate) fn div(first: f64, second: f64) -> f64 {
    nan_rule(first, second, first / second)
}

/// `sqrtsd`.
pub(crate) fn sqrt(x: f64) -> f64 {
    nan_rule(x, x, x.sqrt())
}

/// `andpd` with the absolute-value mask: clears the sign bit, of a NaN too.
pub(crate) fn abs(x: f64) -> f64 {
    f64::from_bits(x.to_bits() & !(1 << 63))
}

/// `cvttsd2si %xmm, %r64`: truncation towards zero, and the "integer
/// indefinite" value `0x8000000000000000` for NaN and for every value
/// outside the signed 64-bit range.
pub(crate) fn cvttsd2si(x: f64) -> i64 {
    const TWO_TO_63: f64 = 9_223_372_036_854_775_808.0;
    if x.is_nan() || !(-TWO_TO_63..TWO_TO_63).contains(&x) {
        i64::MIN
    } else {
        x as i64
    }
}

/// `cvttsd2si %xmm, %r64` followed by the low 32 bits of the register: how
/// the Release build converts a `double` to `UInt` in `intensityScore_`.
///
/// The low 32 bits of the indefinite value are 0. An in-range value keeps
/// its low 32 bits, so `-1.0` becomes `0xffffffff` and `2^32 + 2` becomes 2.
pub(crate) fn truncate_to_u32(x: f64) -> u32 {
    cvttsd2si(x) as u32
}

/// How the Release build converts a `double` to `Size` (`size_t`):
/// `comisd` against `2^63`; below it (or NaN, which compares unordered)
/// `cvttsd2si`; at or above it `subsd 2^63`, `cvttsd2si` and `btc $63`.
///
/// Values in `[0, 2^64)` convert exactly. NaN gives `2^63`, `+inf` and
/// every value of `2^64` and above give 0 (the indefinite value with its
/// top bit flipped), and a negative value in `(-2^63, 0)` wraps modulo
/// `2^64`. This is the instruction sequence of step 2.5 of `run_`
/// (`libOpenMS.so` `0x18e46f4`-`0x18e46fe` and `0x18e6c3b`-`0x18e6c44`)
/// and of `getIsotopeDistribution_`.
pub(crate) fn truncate_to_u64(x: f64) -> u64 {
    const TWO_TO_63: f64 = 9_223_372_036_854_775_808.0;
    if x >= TWO_TO_63 {
        (cvttsd2si(x - TWO_TO_63) as u64) ^ (1 << 63)
    } else {
        cvttsd2si(x) as u64
    }
}

/// `cvtsd2ss`: `x` narrowed to `f32`. A NaN keeps its sign and the top 22
/// bits of its payload and is quieted, as the instruction does; Rust's `as`
/// leaves NaN bits to the host.
pub(crate) fn narrow(x: f64) -> f32 {
    if x.is_nan() {
        let bits = x.to_bits();
        let sign = ((bits >> 32) as u32) & 0x8000_0000;
        let payload = ((bits >> 29) as u32) & 0x003f_ffff;
        f32::from_bits(sign | 0x7fc0_0000 | payload)
    } else {
        x as f32
    }
}

/// `cvtss2sd`: `x` widened to `f64`. A NaN keeps its sign and its payload,
/// shifted to the top of the wider payload, and is quieted, as the
/// instruction does (`DataValue(float)` of `setMetaValue("FWHM", fwhm)`).
pub(crate) fn widen(x: f32) -> f64 {
    if x.is_nan() {
        let bits = x.to_bits();
        let sign = u64::from(bits & 0x8000_0000) << 32;
        let payload = u64::from(bits & 0x003f_ffff) << 29;
        f64::from_bits(sign | 0x7ff8_0000_0000_0000 | payload)
    } else {
        f64::from(x)
    }
}
