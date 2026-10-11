// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Small vector statistics used by the OpenSWATH scores: normalisation, norms,
//! dot products, Manhattan distances and a running mean/standard deviation.
//!
//! Ports `OPENSWATHALGO/ALGO/StatsHelpers.h` with `StatsHelpers.cpp`. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_ALGO_SUPPORT.md).
//!
//! # Summation order
//!
//! The source reduces with Eigen (`norm`, `dot`, `sum`), whose vectorised
//! reductions add in blocks, not left to right. The `f64` reductions
//! ([`norm`], [`dot_prod`], [`manhattan_dist`] and the sums of
//! [`manhattan_scoring`]) add in Eigen's own order, with the Eigen 5.0.1
//! two-lane kernels of the Linux x86-64 Release build that
//! `crate::math::fitters::levenberg_marquardt` carries and measured, so they
//! agree with that build to the last bit (`tests/openswath_scoring_oracle.rs`)
//! and are identical on every machine. The `f32` and `i32` dot products keep
//! their element-type accumulation in plain order.
//!
//! # Lengths
//!
//! The source's iterator-pair functions read the second sequence for as many
//! elements as the first has, unchecked. The port takes slices and returns
//! [`Error::InvalidValue`](crate::Error::InvalidValue) when the second is
//! shorter; a longer second slice is allowed, as in the source.

use crate::math::fitters::levenberg_marquardt::{eigen_dot, eigen_squared_norm, eigen_sum};
use crate::{Error, Result};

/// `Eigen::Map<VectorXd>::sum()`: zero for an empty vector.
fn eigen_total(values: &[f64]) -> f64 {
    if values.is_empty() {
        0.0
    } else {
        eigen_sum(values.len(), |i| values[i])
    }
}

fn check_second_len(first: usize, second: usize) -> Result<()> {
    if second < first {
        return Err(Error::InvalidValue(format!(
            "the second vector has {second} elements, fewer than the first's {first}"
        )));
    }
    Ok(())
}

/// Normalise `intensities` by `normalization_factor` into `normalized`, exactly
/// as the source.
///
/// `normalized` is resized to the length of `intensities`. Only when
/// `normalization_factor > 0` is every element overwritten with
/// `intensity / normalization_factor`; otherwise (zero, negative or NaN) the
/// resized vector is left as it is, so a fresh vector holds zeros and a reused
/// one keeps its old leading values. [`normalize`] is the fresh-vector form.
pub fn normalize_into(intensities: &[f64], normalization_factor: f64, normalized: &mut Vec<f64>) {
    normalized.resize(intensities.len(), 0.0);
    if normalization_factor > 0.0 {
        for (out, &value) in normalized.iter_mut().zip(intensities) {
            *out = value / normalization_factor;
        }
    }
}

/// Normalise intensities by `normalization_factor`.
///
/// Each value is divided by the factor when it is positive; otherwise the
/// result is all zeros, which is what the source produces into a fresh output
/// vector. See [`normalize_into`] for the exact source form.
pub fn normalize(intensities: &[f64], normalization_factor: f64) -> Vec<f64> {
    let mut out = Vec::new();
    normalize_into(intensities, normalization_factor, &mut out);
    out
}

pub(super) fn normalize_in_place(values: &mut [f64], normalization_factor: f64) {
    if normalization_factor > 0.0 {
        for value in values {
            *value /= normalization_factor;
        }
    }
}

/// The Euclidean norm of `values`; 0 for an empty slice.
pub fn norm(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    eigen_squared_norm(values.len(), |i| values[i]).sqrt()
}

/// The dot product of `exp` and the first `exp.len()` elements of `theo`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `theo` is shorter than `exp`.
pub fn dot_prod(exp: &[f64], theo: &[f64]) -> Result<f64> {
    check_second_len(exp.len(), theo.len())?;
    if exp.is_empty() {
        return Ok(0.0);
    }
    Ok(eigen_dot(exp.len(), |i| exp[i], |i| theo[i]))
}

/// The dot product of two `f32` sequences, accumulated in `f32` and widened
/// to `f64`, as the source's `float` instantiation (Eigen computes the dot
/// product in the element type and casts the result).
///
/// # Errors
///
/// [`Error::InvalidValue`] when `theo` is shorter than `exp`.
pub fn dot_prod_f32(exp: &[f32], theo: &[f32]) -> Result<f64> {
    check_second_len(exp.len(), theo.len())?;
    Ok(f64::from(
        exp.iter().zip(theo).map(|(a, b)| a * b).sum::<f32>(),
    ))
}

/// The dot product of two `i32` sequences, accumulated in integers and
/// converted to `f64`, as the source's `int` instantiation.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `theo` is shorter than `exp`, or when a
/// product or the running sum overflows `i32` (undefined behaviour in the
/// source).
pub fn dot_prod_i32(exp: &[i32], theo: &[i32]) -> Result<f64> {
    check_second_len(exp.len(), theo.len())?;
    let mut sum: i32 = 0;
    for (a, b) in exp.iter().zip(theo) {
        sum = a
            .checked_mul(*b)
            .and_then(|product| sum.checked_add(product))
            .ok_or_else(|| Error::InvalidValue("integer dot product overflows i32".into()))?;
    }
    Ok(f64::from(sum))
}

/// The dot-product score of an experimental against a theoretical intensity
/// vector.
///
/// Take the square root of every intensity, normalise each vector by its
/// Euclidean norm, and return their dot product. A vector whose norm is zero is
/// left unnormalised (its square roots are used as they are), as in the source.
/// A negative intensity has no real square root and makes the score NaN.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `theorint` is shorter than `int_exp`.
pub fn dotprod_scoring(mut int_exp: Vec<f64>, mut theorint: Vec<f64>) -> Result<f64> {
    check_second_len(int_exp.len(), theorint.len())?;
    int_exp.iter_mut().for_each(|v| *v = v.sqrt());
    theorint.iter_mut().for_each(|v| *v = v.sqrt());
    let exp_total = norm(&int_exp);
    let theo_total = norm(&theorint);
    normalize_in_place(&mut int_exp, exp_total);
    normalize_in_place(&mut theorint, theo_total);
    dot_prod(&int_exp, &theorint)
}

/// The Manhattan distance, `sum |exp[i] - theo[i]|`, over the elements of
/// `exp`; 0 for an empty `exp`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `theo` is shorter than `exp`.
pub fn manhattan_dist(exp: &[f64], theo: &[f64]) -> Result<f64> {
    check_second_len(exp.len(), theo.len())?;
    if exp.is_empty() {
        return Ok(0.0);
    }
    Ok(eigen_sum(exp.len(), |i| (exp[i] - theo[i]).abs()))
}

/// The Manhattan score of an experimental against a theoretical intensity
/// vector.
///
/// Take the square root of every intensity, normalise each vector by its sum
/// (its TIC), and return their Manhattan distance. A vector whose sum is not
/// positive is left unnormalised, as in the source.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `theorint` is shorter than `int_exp`.
pub fn manhattan_scoring(mut int_exp: Vec<f64>, mut theorint: Vec<f64>) -> Result<f64> {
    check_second_len(int_exp.len(), theorint.len())?;
    int_exp.iter_mut().for_each(|v| *v = v.sqrt());
    theorint.iter_mut().for_each(|v| *v = v.sqrt());
    let exp_total = eigen_total(&int_exp);
    let theo_total = eigen_total(&theorint);
    normalize_in_place(&mut int_exp, exp_total);
    normalize_in_place(&mut theorint, theo_total);
    manhattan_dist(&int_exp, &theorint)
}

/// The Pearson correlation of `x` and the first `x.len()` elements of `y`.
///
/// When either centred sum of squares is below 1e-12 the result is 0. An empty
/// `x` divides by zero and yields NaN, as in the source. The source template
/// computes in the element type of `x`; this port is the `double` form.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `y` is shorter than `x`.
#[deprecated(note = "use math::statistic_functions::pearson_correlation_coefficient instead")]
pub fn cor_pearson(x: &[f64], y: &[f64]) -> Result<f64> {
    check_second_len(x.len(), y.len())?;
    let (mut m1, mut m2, mut s1, mut s2, mut corr) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let n = x.len();
    let nd = n as f64;
    for (a, b) in x.iter().zip(y) {
        corr += a * b;
        m1 += a;
        m2 += b;
        s1 += a * a;
        s2 += b * b;
    }
    m1 /= nd;
    m2 /= nd;
    s1 -= m1 * m1 * nd;
    s2 -= m2 * m2 * nd;
    if s1 < 1.0e-12 || s2 < 1.0e-12 {
        Ok(0.0)
    } else {
        corr -= m1 * m2 * nd;
        corr /= (s1 * s2).sqrt();
        Ok(corr)
    }
}

/// Running mean and standard deviation (Welford's algorithm). Source functor
/// `mean_and_stddev`, meant for `std::for_each`; here feed it with
/// [`Self::push`], [`Extend`] or [`FromIterator`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MeanAndStddev {
    m: f64,
    q: f64,
    c: u64,
}

impl MeanAndStddev {
    /// An empty accumulator.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add one sample (source `operator()(double)`).
    pub fn push(&mut self, sample: f64) {
        let delta = sample - self.m;
        self.c += 1;
        self.m += delta / self.c as f64;
        self.q += delta * (sample - self.m);
    }

    /// The sample variance (divided by `n - 1`); 0 for fewer than two samples.
    pub fn sample_variance(&self) -> f64 {
        if self.c > 1 {
            self.q / (self.c - 1) as f64
        } else {
            0.0
        }
    }

    /// The population variance (divided by `n`); 0 for fewer than two samples.
    pub fn standard_variance(&self) -> f64 {
        if self.c > 1 {
            self.q / self.c as f64
        } else {
            0.0
        }
    }

    /// The sample standard deviation.
    pub fn sample_stddev(&self) -> f64 {
        self.sample_variance().sqrt()
    }

    /// The population standard deviation.
    pub fn standard_stddev(&self) -> f64 {
        self.standard_variance().sqrt()
    }

    /// The mean; 0 before any sample.
    pub fn mean(&self) -> f64 {
        self.m
    }

    /// Number of samples.
    pub fn count(&self) -> u64 {
        self.c
    }

    /// The sample variance, as [`Self::sample_variance`].
    pub fn variance(&self) -> f64 {
        self.sample_variance()
    }

    /// The sample standard deviation, as [`Self::sample_stddev`]; also the
    /// source's nullary `operator()()`.
    pub fn stddev(&self) -> f64 {
        self.sample_stddev()
    }
}

impl Extend<f64> for MeanAndStddev {
    fn extend<I: IntoIterator<Item = f64>>(&mut self, iter: I) {
        for sample in iter {
            self.push(sample);
        }
    }
}

impl FromIterator<f64> for MeanAndStddev {
    fn from_iter<I: IntoIterator<Item = f64>>(iter: I) -> Self {
        let mut out = Self::new();
        out.extend(iter);
        out
    }
}
