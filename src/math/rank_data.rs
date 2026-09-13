// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! SciPy-compatible ranking with selectable tie and NaN handling.
//!
//! Port of `src/openms/include/OpenMS/MATH/STATISTICS/RankData.h`, a
//! header-only struct of static templates with no translation unit. See
//! `docs/RANK_DATA_SUPPORT.md`.
//!
//! The source replicates `scipy.stats.rankdata(a, method=..., nan_policy=...)`
//! with one-dimensional semantics, and its class test pins whole vectors
//! against SciPy 1.17.1. Ranks are one-based, as in SciPy.
//!
//! This is **not** the ranking that
//! [`crate::math::statistic_functions::rank_correlation_coefficient`] uses.
//! [`crate::math::statistic_functions::compute_rank`] detects ties with a
//! relative tolerance and only ever averages them; the function here compares
//! by exact equality and offers five tie rules. The two disagree on inputs with
//! near-ties, and both are kept because both have callers.

use crate::{Error, Result};

/// Maximum number of values one call may rank.
///
/// The call stages two owned buffers of this length. Native-only guard; the
/// source allocates whatever it is handed.
pub const MAX_ITEMS: usize = 50_000_000;

/// How tied values share their ranks.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum RankMethod {
    /// The mean of the ranks the tied values would have received. SciPy's
    /// default, and the source's.
    #[default]
    Average,
    /// The lowest of those ranks, sometimes called competition ranking.
    Min,
    /// The highest of those ranks.
    Max,
    /// Like [`RankMethod::Min`], but the next distinct value takes the next
    /// integer rank rather than skipping the tied block.
    Dense,
    /// No tie at all: ties are broken by original position, so every rank is
    /// distinct.
    Ordinal,
}

/// What ranking does when the input contains a NaN.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum NanPolicy {
    /// A single NaN makes every output NaN. SciPy's documented behaviour and
    /// the source's default.
    #[default]
    Propagate,
    /// Rank the non-NaN values among themselves; NaN positions stay NaN.
    Omit,
    /// Refuse the input.
    Raise,
}

/// A value that can be ranked.
///
/// Implemented for `f64`, `f32` and `i32`, matching the source's
/// `rankdata_double`, `rankdata_float` and `rankdata_int` forwarders. The
/// source's template promotes every type to `double` before comparing, and so
/// does this trait, which is what makes an `f32` tie and an `f64` tie agree.
pub trait Rankable: Copy {
    /// The value widened to `f64` for comparison.
    fn as_f64(self) -> f64;

    /// Whether the value is a NaN. Always false for an integer type, as the
    /// source's `if constexpr` branch makes it.
    fn is_nan_value(self) -> bool;
}

impl Rankable for f64 {
    fn as_f64(self) -> f64 {
        self
    }
    fn is_nan_value(self) -> bool {
        self.is_nan()
    }
}

impl Rankable for f32 {
    fn as_f64(self) -> f64 {
        f64::from(self)
    }
    fn is_nan_value(self) -> bool {
        self.is_nan()
    }
}

impl Rankable for i32 {
    fn as_f64(self) -> f64 {
        f64::from(self)
    }
    fn is_nan_value(self) -> bool {
        false
    }
}

/// Rank the values of `a`, one-based, with the given tie and NaN handling.
///
/// The returned vector has one entry per input position, in input order.
/// Positions are ordered by value with a **stable** sort, so
/// [`RankMethod::Ordinal`] breaks ties by original position and every other
/// method sees tied values as one contiguous block.
///
/// Ties are decided by exact equality of the values widened to `f64`. That is
/// SciPy's rule and differs from
/// [`crate::math::statistic_functions::compute_rank`], which uses a relative
/// tolerance.
///
/// An empty input returns an empty vector rather than an error.
///
/// # Arguments
///
/// * `method` — see [`RankMethod`]; the source's default is
///   [`RankMethod::Average`].
/// * `nan_policy` — see [`NanPolicy`]; the source's default is
///   [`NanPolicy::Propagate`].
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `nan_policy` is [`NanPolicy::Raise`]
/// and the input contains a NaN — the source throws `std::invalid_argument`
/// there — and [`Error::InvalidRange`] when the input exceeds [`MAX_ITEMS`].
///
/// ```
/// use openms::math::rank_data::{NanPolicy, RankMethod, rankdata};
///
/// let ranks = rankdata(&[0.0, 2.0, 3.0, 2.0], RankMethod::Average, NanPolicy::Propagate)?;
/// assert_eq!(ranks, vec![1.0, 2.5, 4.0, 2.5]);
/// # Ok::<(), openms::Error>(())
/// ```
pub fn rankdata<T: Rankable>(
    a: &[T],
    method: RankMethod,
    nan_policy: NanPolicy,
) -> Result<Vec<f64>> {
    let n = a.len();
    if n > MAX_ITEMS {
        return Err(Error::InvalidRange(
            "rank input exceeds the supported element count".to_string(),
        ));
    }
    let mut ranks = vec![f64::NAN; n];
    if n == 0 {
        return Ok(ranks);
    }

    let any_nan = a.iter().any(|value| value.is_nan_value());
    if any_nan {
        match nan_policy {
            NanPolicy::Propagate => return Ok(ranks),
            NanPolicy::Raise => {
                return Err(Error::InvalidValue(
                    "NaN present but the NaN policy is Raise".to_string(),
                ));
            }
            NanPolicy::Omit => {}
        }
    }

    let mut idx: Vec<usize> = Vec::with_capacity(n);
    for (index, value) in a.iter().enumerate() {
        if nan_policy == NanPolicy::Omit && value.is_nan_value() {
            continue;
        }
        idx.push(index);
    }
    if idx.is_empty() {
        return Ok(ranks);
    }
    // Stable, so equal values keep their input order; no NaN can reach the
    // comparison, because Propagate and Raise have already returned and Omit
    // has filtered them out. The comparison is the source's `<` rather than the
    // IEEE-754 total order, so that `-0.0` and `0.0` stay interchangeable and
    // `Ordinal` still ranks them by input position.
    idx.sort_by(|&left, &right| {
        a[left]
            .as_f64()
            .partial_cmp(&a[right].as_f64())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    match method {
        RankMethod::Ordinal => {
            for (position, &origin) in idx.iter().enumerate() {
                ranks[origin] = (position + 1) as f64;
            }
        }
        RankMethod::Dense => {
            let mut lo = 0usize;
            let mut dense_rank = 1usize;
            while lo < idx.len() {
                let hi = tie_block_end(a, &idx, lo);
                for &origin in &idx[lo..hi] {
                    ranks[origin] = dense_rank as f64;
                }
                dense_rank += 1;
                lo = hi;
            }
        }
        RankMethod::Min | RankMethod::Max | RankMethod::Average => {
            let mut lo = 0usize;
            while lo < idx.len() {
                let hi = tie_block_end(a, &idx, lo);
                // One-based ranks of this block span [lo + 1, hi].
                let r_min = (lo + 1) as f64;
                let r_max = hi as f64;
                let value = match method {
                    RankMethod::Min => r_min,
                    RankMethod::Max => r_max,
                    _ => 0.5 * (r_min + r_max),
                };
                for &origin in &idx[lo..hi] {
                    ranks[origin] = value;
                }
                lo = hi;
            }
        }
    }
    Ok(ranks)
}

/// End of the block of values equal to `a[idx[lo]]`, starting at `lo`.
fn tie_block_end<T: Rankable>(a: &[T], idx: &[usize], lo: usize) -> usize {
    let value = a[idx[lo]].as_f64();
    let mut hi = lo + 1;
    while hi < idx.len() && a[idx[hi]].as_f64() == value {
        hi += 1;
    }
    hi
}

/// [`rankdata`] over `f64`, the source's `rankdata_double`.
///
/// # Errors
///
/// As [`rankdata`].
pub fn rankdata_f64(a: &[f64], method: RankMethod, nan_policy: NanPolicy) -> Result<Vec<f64>> {
    rankdata(a, method, nan_policy)
}

/// [`rankdata`] over `f32`, the source's `rankdata_float`.
///
/// # Errors
///
/// As [`rankdata`].
pub fn rankdata_f32(a: &[f32], method: RankMethod, nan_policy: NanPolicy) -> Result<Vec<f64>> {
    rankdata(a, method, nan_policy)
}

/// [`rankdata`] over `i32`, the source's `rankdata_int`.
///
/// An integer input never contains a NaN, so `nan_policy` cannot change the
/// result.
///
/// # Errors
///
/// As [`rankdata`].
pub fn rankdata_i32(a: &[i32], method: RankMethod, nan_policy: NanPolicy) -> Result<Vec<f64>> {
    rankdata(a, method, nan_policy)
}
