// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Scoring functions used by MRM scoring: distances, spectral angle,
//! cross-correlation, ranks and ranked mutual information.
//!
//! Ports `OPENSWATHALGO/ALGO/Scoring.h` with `Scoring.cpp` (the source
//! namespace `OpenSwath::Scoring`). See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_ALGO_SUPPORT.md).
//!
//! # Preconditions
//!
//! The source states its preconditions with `OPENSWATH_PRECONDITION`, which is
//! compiled out of release builds, where a violation reads out of bounds. The
//! port checks every one and returns [`Error::InvalidValue`](crate::Error::InvalidValue).
//!
//! # Summation order
//!
//! As in [`stats_helpers`](crate::analysis::openswath::stats_helpers), the source's Eigen
//! reductions add in blocks and most of the port adds left to right; results
//! agree to rounding and are identical on every machine. The upstream class
//! test pins the source's own Eigen cross-correlation against a left-to-right
//! reference loop with a relative tolerance.
//!
//! [`standardize_data`] and [`calculate_cross_correlation`] (and so
//! [`normalized_cross_correlation`]) add in Eigen's own order instead: the
//! vectorised `redux` of Eigen 3.4 with two-lane SSE2 `double` packets, as the
//! Linux x86-64 Release build compiles it. `MRMTransitionGroupPicker`'s peak
//! quality is built from them and is compared with that build bit for bit
//! (`tests/openswath_peak_picking.rs`); a left-to-right sum moved it by one
//! unit in the last place.

use crate::{Error, Result};

/// A 2-D position `(row, column)`. Source typedef `Scoring::pos2D`.
///
/// The source's `pair_hash` functor (hashing a pair as `h1 ^ h2`) has no port:
/// Rust tuples implement `Hash` already.
pub type Pos2D = (u32, u32);

/// One cross-correlation entry, `(lag, correlation)`. Source typedef `XCorrEntry`.
pub type XCorrEntry = (i32, f64);

/// Upper bound on the number of entries a cross-correlation may produce,
/// checked before anything is allocated (16 Mi entries, 256 MiB).
pub const MAX_XCORR_ENTRIES: usize = 1 << 24;

/// A cross-correlation array of `(lag, correlation)` pairs in ascending lag
/// order. Source struct `XCorrArrayType`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct XCorrArrayType {
    /// The entries.
    pub data: Vec<XCorrEntry>,
}

impl XCorrArrayType {
    /// Iterate over the entries (source `begin()`/`end()`).
    pub fn iter(&self) -> std::slice::Iter<'_, XCorrEntry> {
        self.data.iter()
    }

    /// Iterate mutably over the entries.
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, XCorrEntry> {
        self.data.iter_mut()
    }
}

impl<'a> IntoIterator for &'a XCorrArrayType {
    type Item = &'a XCorrEntry;
    type IntoIter = std::slice::Iter<'a, XCorrEntry>;
    fn into_iter(self) -> Self::IntoIter {
        self.data.iter()
    }
}

fn check_same_nonempty(x: usize, y: usize) -> Result<()> {
    if x != y || x == 0 {
        return Err(Error::InvalidValue(format!(
            "both vectors need to have the same non-zero length, got {x} and {y}"
        )));
    }
    Ok(())
}

/// Divide each element of `x` by the sum of the vector, in place; `x` is
/// unchanged when the sum is zero.
///
/// A NaN sum is not zero, so it spreads NaN, as in the source.
pub fn normalize_sum(x: &mut [f64]) {
    let s: f64 = x.iter().sum();
    if s != 0.0 {
        for v in x {
            *v /= s;
        }
    }
}

/// The normalised Manhattan distance between two vectors.
///
/// Equivalent to the function `delta_ratio_sum` from mQuest, used to compare
/// library and experimental intensities. Both vectors are first normalised
/// **in place** by [`normalize_sum`]; the result is
///
/// ```text
/// d = 1/N * sum_i |x_i / sum(x) - y_i / sum(y)|
/// ```
///
/// The source documents the formula with a square root that the code does not
/// take; the port follows the code, which the class test pins.
///
/// # Errors
///
/// [`Error::InvalidValue`] unless both vectors have the same non-zero length;
/// both are then unchanged.
pub fn normalized_manhattan_dist(x: &mut [f64], y: &mut [f64]) -> Result<f64> {
    check_same_nonempty(x.len(), y.len())?;
    normalize_sum(x);
    normalize_sum(y);
    let sum: f64 = x.iter().zip(y.iter()).map(|(a, b)| (a - b).abs()).sum();
    Ok(sum / x.len() as f64)
}

/// The root mean square deviation between two vectors,
/// `sqrt(1/N * sum_i (x_i - y_i)^2)`.
///
/// # Errors
///
/// [`Error::InvalidValue`] unless both vectors have the same non-zero length.
pub fn root_mean_square_deviation(x: &[f64], y: &[f64]) -> Result<f64> {
    check_same_nonempty(x.len(), y.len())?;
    let squared: f64 = x.iter().zip(y).map(|(a, b)| (a - b) * (a - b)).sum();
    Ok((squared / x.len() as f64).sqrt())
}

/// The spectral angle (arc cosine of the normalised dot product) between two
/// vectors:
///
/// ```text
/// theta = acos( sum_i x_i*y_i / sqrt(sum_i x_i^2 * sum_i y_i^2) )
/// ```
///
/// When either vector has zero length the cosine is taken as 0, so the angle is
/// pi/2; the cosine is clamped to `[-1, 1]` before `acos`.
///
/// # Errors
///
/// [`Error::InvalidValue`] unless both vectors have the same non-zero length.
pub fn spectral_angle(x: &[f64], y: &[f64]) -> Result<f64> {
    check_same_nonempty(x.len(), y.len())?;
    let x_len = x.iter().map(|v| v * v).sum::<f64>().sqrt();
    let y_len = y.iter().map(|v| v * v).sum::<f64>().sqrt();
    let denominator = x_len * y_len;
    let dot: f64 = x.iter().zip(y).map(|(a, b)| a * b).sum();
    let theta = if denominator == 0.0 {
        0.0
    } else {
        dot / denominator
    };
    Ok(theta.clamp(-1.0, 1.0).acos())
}

/// Sum `value(0) .. value(n - 1)` in the order of Eigen 3.4's vectorised
/// `redux` (`LinearVectorizedTraversal`, no unrolling) with two-lane SSE2
/// `double` packets and a first aligned index of zero.
///
/// The first aligned index is zero for every expression the source sums
/// here: a `Map` over `std::vector` storage, which glibc's allocator aligns to
/// 16 bytes, and element-wise expressions, for which Eigen always starts at
/// zero. Two packet accumulators take alternating pairs of values, are added,
/// take one more pair when one is left, and are reduced lane 0 plus lane 1; an
/// odd last value is added after that. A single value is returned as it is,
/// and an empty sum is zero.
fn eigen_sum(n: usize, value: impl Fn(usize) -> f64) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let aligned_size = n / 2 * 2;
    if aligned_size == 0 {
        return value(0);
    }
    let aligned_size2 = n / 4 * 4;
    let mut first = [value(0), value(1)];
    if aligned_size > 2 {
        let mut second = [value(2), value(3)];
        let mut index = 4;
        while index < aligned_size2 {
            first[0] += value(index);
            first[1] += value(index + 1);
            second[0] += value(index + 2);
            second[1] += value(index + 3);
            index += 4;
        }
        first = [first[0] + second[0], first[1] + second[1]];
        if aligned_size > aligned_size2 {
            first[0] += value(aligned_size2);
            first[1] += value(aligned_size2 + 1);
        }
    }
    let mut sum = first[0] + first[1];
    for index in aligned_size..n {
        sum += value(index);
    }
    sum
}

/// Standardise a vector in place: subtract the mean and divide by the
/// population standard deviation (divisor `N`).
///
/// An empty vector and an all-zero vector are left unchanged; a constant
/// non-zero vector (standard deviation 0) is divided by 1 instead, so it
/// becomes all zeros. Both sums are added in Eigen's order (see the module
/// documentation).
pub fn standardize_data(data: &mut [f64]) {
    if data.is_empty() {
        return;
    }
    let n = data.len() as f64;
    let mean = eigen_sum(data.len(), |i| data[i]) / n;
    let mut stdev = (eigen_sum(data.len(), |i| {
        let centred = data[i] - mean;
        centred * centred
    }) / n)
        .sqrt();
    if mean == 0.0 && stdev == 0.0 {
        return;
    }
    if stdev == 0.0 {
        stdev = 1.0;
    }
    for v in data {
        *v = (*v - mean) / stdev;
    }
}

fn xcorr_entries(maxdelay: i32, lag: i32) -> Result<usize> {
    if lag <= 0 {
        return Err(Error::InvalidValue(format!(
            "cross-correlation lag must be positive, got {lag}"
        )));
    }
    if maxdelay < 0 {
        return Ok(0);
    }
    let entries = (2 * i64::from(maxdelay)) / i64::from(lag) + 1;
    let entries = usize::try_from(entries).unwrap_or(usize::MAX);
    if entries > MAX_XCORR_ENTRIES {
        return Err(Error::InvalidValue(format!(
            "cross-correlation would produce {entries} entries, more than {MAX_XCORR_ENTRIES}"
        )));
    }
    Ok(entries)
}

/// Cross-correlation of two equally long vectors without normalisation.
///
/// For every delay `d` from `-maxdelay` to `maxdelay` in steps of `lag`, the
/// entry is `(d, sum_i data1[i] * data2[i + d])` over the indices where both
/// exist; a delay with no overlap yields 0. A negative `maxdelay` yields an
/// empty array. Each dot product is added in Eigen's order (see the module
/// documentation).
///
/// # Errors
///
/// [`Error::InvalidValue`] when the vectors differ in length, when `lag` is not
/// positive (the source loops forever for 0), or when the array would exceed
/// [`MAX_XCORR_ENTRIES`].
pub fn calculate_cross_correlation(
    data1: &[f64],
    data2: &[f64],
    maxdelay: i32,
    lag: i32,
) -> Result<XCorrArrayType> {
    if data1.len() != data2.len() {
        return Err(Error::InvalidValue(format!(
            "both data vectors need to have the same length, got {} and {}",
            data1.len(),
            data2.len()
        )));
    }
    let entries = xcorr_entries(maxdelay, lag)?;
    let size = i64::try_from(data1.len()).unwrap_or(i64::MAX);
    let mut result = XCorrArrayType {
        data: Vec::with_capacity(entries),
    };
    let mut delay = -i64::from(maxdelay);
    while delay <= i64::from(maxdelay) {
        let start = 0.max(-delay);
        let end = size.min(size - delay);
        let mut sxy = 0.0;
        if end > start {
            // Both bounds lie in [0, size], so the conversions cannot fail.
            let (s, e) = (start as usize, end as usize);
            let shifted = (start + delay) as usize;
            sxy = eigen_sum(e - s, |i| data1[s + i] * data2[shifted + i]);
        }
        // |delay| <= maxdelay, which fits i32.
        result.data.push((delay as i32, sxy));
        delay += i64::from(lag);
    }
    Ok(result)
}

/// Cross-correlation of two already normalised vectors: the entries of
/// [`calculate_cross_correlation`] divided by the vector length.
///
/// Empty vectors divide 0 by 0 and yield NaN entries, as in the source.
///
/// # Errors
///
/// Those of [`calculate_cross_correlation`].
pub fn normalized_cross_correlation_post(
    normalized_data1: &[f64],
    normalized_data2: &[f64],
    maxdelay: i32,
    lag: i32,
) -> Result<XCorrArrayType> {
    let mut result =
        calculate_cross_correlation(normalized_data1, normalized_data2, maxdelay, lag)?;
    let n = normalized_data1.len() as f64;
    for entry in result.iter_mut() {
        entry.1 /= n;
    }
    Ok(result)
}

/// Normalised cross-correlation: standardise both vectors **in place** with
/// [`standardize_data`], then [`normalized_cross_correlation_post`].
///
/// This replaces the source's older `calcxcorr`. The source declares a default
/// `lag = 1` on the definition only, where it is unusable; pass 1.
///
/// # Errors
///
/// [`Error::InvalidValue`] unless both vectors have the same non-zero length
/// (checked before either is modified), and those of
/// [`calculate_cross_correlation`].
pub fn normalized_cross_correlation(
    data1: &mut [f64],
    data2: &mut [f64],
    maxdelay: i32,
    lag: i32,
) -> Result<XCorrArrayType> {
    check_same_nonempty(data1.len(), data2.len())?;
    xcorr_entries(maxdelay, lag)?;
    standardize_data(data1);
    standardize_data(data2);
    normalized_cross_correlation_post(data1, data2, maxdelay, lag)
}

/// Index of the best peak of a cross-correlation: the entry with the highest
/// correlation, the first one on ties. Source `xcorrArrayGetMaxPeak`, which
/// returns an iterator.
///
/// # Errors
///
/// [`Error::InvalidValue`] for an empty array.
pub fn xcorr_array_get_max_peak(array: &XCorrArrayType) -> Result<usize> {
    let first = array.data.first().ok_or_else(|| {
        Error::InvalidValue("cannot get the highest apex from an empty array".into())
    })?;
    let mut max = first.1;
    let mut max_index = 0;
    for (i, entry) in array.data.iter().enumerate() {
        if entry.1 > max {
            max = entry.1;
            max_index = i;
        }
    }
    Ok(max_index)
}

/// Compute the ranks of the elements of `v` into `ranks_out` and return the
/// highest rank.
///
/// The rank of a value is the position, in ascending order, of the first
/// element equal to it: equal values share a rank and the next distinct value
/// takes its own position, so ranks can skip: `[5, 3, 3, 8]` ranks
/// `[2, 0, 0, 3]` and the highest rank is 3. Despite the name,
/// `ranks_out` is **replaced**, not appended to, as in the source; the
/// appending happens in [`compute_rank_vector`].
///
/// The source sorts with `<`, which is undefined for NaN input. The port sorts
/// by IEEE total order: identical for every other input, and each NaN then gets
/// a rank of its own.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `v` has more than `u32::MAX` elements (the
/// source truncates the index).
pub fn compute_and_append_rank(v: &[f64], ranks_out: &mut Vec<u32>) -> Result<u32> {
    if u32::try_from(v.len()).is_err() {
        return Err(Error::InvalidValue(
            "cannot rank more than u32::MAX values".into(),
        ));
    }
    let mut order: Vec<usize> = (0..v.len()).collect();
    order.sort_by(|&i, &j| v[i].total_cmp(&v[j]));
    ranks_out.clear();
    ranks_out.resize(v.len(), 0);
    let mut x = 0.0;
    let mut y: u32 = 0;
    for (i, &index) in order.iter().enumerate() {
        if v[index] != x {
            x = v[index];
            // i < v.len() <= u32::MAX, checked above.
            y = i as u32;
        }
        ranks_out[index] = y;
    }
    Ok(y)
}

/// Rank every row of `intensity`, append the rank vectors to `ranks`, and
/// return the highest rank of each row.
///
/// # Errors
///
/// Those of [`compute_and_append_rank`]; `ranks` is unchanged on error.
pub fn compute_rank_vector(intensity: &[Vec<f64>], ranks: &mut Vec<Vec<u32>>) -> Result<Vec<u32>> {
    let mut new_ranks = Vec::with_capacity(intensity.len());
    let mut max_rank_vec = Vec::with_capacity(intensity.len());
    for row in intensity {
        let mut row_ranks = Vec::new();
        max_rank_vec.push(compute_and_append_rank(row, &mut row_ranks)?);
        new_ranks.push(row_ranks);
    }
    ranks.extend(new_ranks);
    Ok(max_rank_vec)
}

/// Estimate the mutual information, in bits, between two vectors of ranks.
///
/// With `n` pairs, joint-state counts `c_ab` and marginal counts `c_a`, `c_b`:
///
/// ```text
/// MI = ( sum_ab c_ab * ln(c_ab / c_a / c_b) / n + ln(n) ) / ln(2)
/// ```
///
/// Joint states are visited in ascending `(a, b)` order, as in the source, so
/// the floating-point sum is the source's sum.
///
/// # Errors
///
/// [`Error::InvalidValue`] unless both vectors have the same non-zero length,
/// and when a rank exceeds its `max_rank` (the source then writes out of
/// bounds) or the vectors are longer than `u32::MAX`.
pub fn ranked_mutual_information(
    ranked_data1: &[u32],
    ranked_data2: &[u32],
    max_rank1: u32,
    max_rank2: u32,
) -> Result<f64> {
    check_same_nonempty(ranked_data1.len(), ranked_data2.len())?;
    if u32::try_from(ranked_data1.len()).is_err() {
        return Err(Error::InvalidValue(
            "cannot compute mutual information over more than u32::MAX values".into(),
        ));
    }
    if let Some(&bad) = ranked_data1.iter().find(|&&r| r > max_rank1) {
        return Err(Error::InvalidValue(format!(
            "rank {bad} exceeds the first maximum rank {max_rank1}"
        )));
    }
    if let Some(&bad) = ranked_data2.iter().find(|&&r| r > max_rank2) {
        return Err(Error::InvalidValue(format!(
            "rank {bad} exceeds the second maximum rank {max_rank2}"
        )));
    }
    let n = ranked_data1.len();
    let first_states = max_rank1 as usize + 1;
    let second_states = u64::from(max_rank2) + 1;
    let mut first_counts = vec![0u32; first_states];
    let mut second_counts = vec![0u32; max_rank2 as usize + 1];
    let mut joint: Vec<u64> = Vec::with_capacity(n);
    for (&a, &b) in ranked_data1.iter().zip(ranked_data2) {
        first_counts[a as usize] += 1;
        second_counts[b as usize] += 1;
        joint.push(u64::from(a) * second_states + u64::from(b));
    }
    joint.sort_unstable();

    let mut mutual_information = 0.0;
    let mut i = 0;
    while i < joint.len() {
        let state = joint[i];
        let mut j = i + 1;
        while j < joint.len() && joint[j] == state {
            j += 1;
        }
        let count = (j - i) as f64;
        let a = (state / second_states) as usize;
        let b = (state % second_states) as usize;
        mutual_information +=
            count * (count / f64::from(first_counts[a]) / f64::from(second_counts[b])).ln();
        i = j;
    }
    let n = n as f64;
    mutual_information /= n;
    mutual_information += n.ln();
    mutual_information /= 2.0f64.ln();
    Ok(mutual_information)
}

#[cfg(test)]
mod tests {
    use super::eigen_sum;

    #[test]
    fn eigen_sum_follows_the_two_packet_order() {
        // Left to right, 1e16 + 1 rounds back to 1e16 and the sum is 2; Eigen
        // pairs lane 0 (1e16, -1e16) and lane 1 (1, 1) first and gets 3.
        let values = [1e16, 1.0, -1e16, 1.0, 1.0];
        assert_eq!(eigen_sum(values.len(), |i| values[i]), 3.0);
        assert_eq!(values.iter().sum::<f64>(), 2.0);
        assert_eq!(eigen_sum(0, |_| unreachable!()), 0.0);
        assert_eq!(eigen_sum(1, |_| 7.5), 7.5);
        // Six values: the two accumulators, then one more pair, then lanes.
        let six = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        assert_eq!(
            eigen_sum(6, |i| six[i]),
            ((1.0 + 3.0) + 5.0) + ((2.0 + 4.0) + 6.0)
        );
    }
}
