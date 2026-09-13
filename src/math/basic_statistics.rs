// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! An accumulating weighted distribution with a normal approximation.
//!
//! Port of `src/openms/include/OpenMS/MATH/STATISTICS/BasicStatistics.h`, a
//! header-only template with no accompanying translation unit. See
//! `docs/BASIC_STATISTICS_SUPPORT.md`.
//!
//! The source is `template <typename RealT = double> class BasicStatistics`.
//! Only the `double` instantiation exists in the library and in its class test,
//! so [`crate::math::basic_statistics::BasicStatistics`] is not generic; an
//! `f32` instantiation would only lose precision in an accumulation this port
//! deliberately performs in `f64`.
//!
//! These are *weighted* moments over a probability vector, not the sample
//! moments of [`crate::math::statistic_functions`]. The weights are the values
//! themselves and the divisor is their sum, so the degrees of freedom question
//! does not arise here.

use crate::{Error, Result};

/// Maximum number of probabilities or coordinates one call may allocate.
///
/// Native-only guard: the source resizes to whatever size it is handed.
pub const MAX_ITEMS: usize = 50_000_000;

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.to_string())
}

/// Weighted sum, mean and variance of a distribution, plus its normal
/// approximation.
///
/// Intended usage follows the source: create an instance, set the three
/// parameters by [`BasicStatistics::update`],
/// [`BasicStatistics::update_with_coordinates`] or the individual setters, then
/// read them back or draw a normal approximation from them.
///
/// The three parameters are independent once set: the setters exist so that a
/// caller can describe a distribution it never measured, and
/// [`BasicStatistics::normal_density`] uses only `mean` and `variance`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BasicStatistics {
    mean: f64,
    variance: f64,
    sum: f64,
}

impl BasicStatistics {
    /// `sqrt(2 * pi)`, which normalises
    /// [`BasicStatistics::normal_density_sqrt2pi`].
    ///
    /// This is the source's literal `2.50662827463100050240`, spelled with the
    /// digits that round-trip to the same `f64`. That value is one unit in the
    /// last place *above* the correctly rounded `sqrt(2 * pi)`, which is
    /// `2.5066282746310002`; the source's constant is reproduced rather than
    /// corrected, because its class test compares against it exactly and every
    /// density the library has ever produced carries it.
    pub const SQRT_2PI: f64 = 2.506_628_274_631_000_7;

    /// A distribution with zero sum, mean and variance.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set sum, mean and variance to zero.
    pub fn clear(&mut self) {
        self.mean = 0.0;
        self.variance = 0.0;
        self.sum = 0.0;
    }

    /// Recompute the parameters from a probability vector at positions
    /// `0, 1, ..., n-1`.
    ///
    /// May be called as often as needed with different inputs; each call clears
    /// first. The accumulation is the source's two passes, both in slice order:
    /// the first forms `sum` and `sum(p_i * i)` together and divides to get the
    /// mean, the second forms `sum(p_i * (i - mean)^2)` and divides by the same
    /// sum. Fusing the passes or using a shifted-moment identity would change
    /// the last bits, so neither is done.
    ///
    /// A zero probability mass makes both quotients NaN; the source detects
    /// that case afterwards and substitutes zeros, and so does this. An empty
    /// input therefore leaves all three parameters at zero.
    ///
    /// The source counts positions in `unsigned`, which wraps after 2^32
    /// entries; `usize` here does not, and the count is capped by
    /// [`MAX_ITEMS`] anyway.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the input exceeds [`MAX_ITEMS`].
    /// The parameters are cleared only after that check, so a refused call
    /// leaves the previous state intact.
    pub fn update(&mut self, probabilities: &[f64]) -> Result<()> {
        if probabilities.len() > MAX_ITEMS {
            return Err(bad("probability vector exceeds the supported length"));
        }
        self.clear();
        for (pos, probability) in probabilities.iter().enumerate() {
            self.sum += probability;
            self.mean += probability * pos as f64;
        }
        self.mean /= self.sum;

        for (pos, probability) in probabilities.iter().enumerate() {
            let mut diff = pos as f64 - self.mean;
            diff *= diff;
            self.variance += probability * diff;
        }
        self.variance /= self.sum;

        if self.sum == 0.0 && (self.mean.is_nan() || self.mean.is_infinite()) {
            self.mean = 0.0;
            self.variance = 0.0;
        }
        Ok(())
    }

    /// Recompute the parameters from a probability vector evaluated at explicit
    /// coordinates.
    ///
    /// The accumulation mirrors [`BasicStatistics::update`] with `coordinates`
    /// in place of the positions.
    ///
    /// This overload has **no** zero-mass guard in the source, so a zero
    /// probability sum leaves `mean` and `variance` NaN where the
    /// position-based overload would report zeros. That asymmetry is source
    /// behaviour and is preserved; callers that accept untrusted weights should
    /// test [`BasicStatistics::sum`] before reading the moments.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the two slices have different
    /// lengths, or when the input exceeds [`MAX_ITEMS`]. The source takes only
    /// a coordinate *begin* iterator and advances it once per probability, so a
    /// short coordinate range is read out of bounds; the port compares the
    /// lengths instead. The state is cleared only after both checks.
    pub fn update_with_coordinates(
        &mut self,
        probabilities: &[f64],
        coordinates: &[f64],
    ) -> Result<()> {
        if probabilities.len() != coordinates.len() {
            return Err(bad(
                "probability and coordinate vectors must have the same length",
            ));
        }
        if probabilities.len() > MAX_ITEMS {
            return Err(bad("probability vector exceeds the supported length"));
        }
        self.clear();
        for (probability, coordinate) in probabilities.iter().zip(coordinates.iter()) {
            self.sum += probability;
            self.mean += probability * coordinate;
        }
        self.mean /= self.sum;

        for (probability, coordinate) in probabilities.iter().zip(coordinates.iter()) {
            let mut diff = coordinate - self.mean;
            diff *= diff;
            self.variance += probability * diff;
        }
        self.variance /= self.sum;
        Ok(())
    }

    /// The mean.
    pub fn mean(&self) -> f64 {
        self.mean
    }

    /// Set the mean, without touching the sum or the variance.
    pub fn set_mean(&mut self, mean: f64) {
        self.mean = mean;
    }

    /// The variance.
    pub fn variance(&self) -> f64 {
        self.variance
    }

    /// Set the variance, without touching the sum or the mean.
    pub fn set_variance(&mut self, variance: f64) {
        self.variance = variance;
    }

    /// The probability mass.
    pub fn sum(&self) -> f64 {
        self.sum
    }

    /// Set the probability mass, without touching the mean or the variance.
    pub fn set_sum(&mut self, sum: f64) {
        self.sum = sum;
    }

    /// `sqrt(2 * pi)`; see [`BasicStatistics::SQRT_2PI`].
    ///
    /// Present because the source exposes it as a static member function that
    /// callers use to normalise
    /// [`BasicStatistics::normal_density_sqrt2pi`] themselves.
    pub fn sqrt2pi() -> f64 {
        Self::SQRT_2PI
    }

    /// Density of the normal approximation at `coordinate`, times
    /// `sqrt(2 * pi)`.
    ///
    /// `exp(-(coordinate - mean)^2 / 2 / variance)`, saving the division that
    /// [`BasicStatistics::normal_density`] performs. Note that this is not
    /// scaled by `1 / sqrt(variance)` either: it is `exp` of the standardised
    /// square, so it equals `1` at the mean for every variance.
    ///
    /// A zero variance is the source's unchecked `x / 0`: the density is `0`
    /// away from the mean and NaN exactly at it, and this port returns those
    /// two values directly instead of dividing. A negative variance — only
    /// reachable through [`BasicStatistics::set_variance`] — is left to the
    /// source's arithmetic and grows without bound.
    pub fn normal_density_sqrt2pi(&self, coordinate: f64) -> f64 {
        let mut shifted = coordinate - self.mean;
        shifted *= shifted;
        if self.variance == 0.0 {
            return if shifted == 0.0 { f64::NAN } else { 0.0 };
        }
        (-shifted / 2.0 / self.variance).exp()
    }

    /// Density of the normal approximation at `coordinate`.
    ///
    /// [`BasicStatistics::normal_density_sqrt2pi`] divided by
    /// [`BasicStatistics::SQRT_2PI`].
    pub fn normal_density(&self, coordinate: f64) -> f64 {
        self.normal_density_sqrt2pi(coordinate) / Self::SQRT_2PI
    }

    /// Normal approximation sampled at coordinate positions
    /// `0, 1, ..., size-1`.
    ///
    /// Each entry is `density(i) / sum(density) * sum()`, in that order, so the
    /// returned values carry the distribution's own mass rather than
    /// integrating to one.
    ///
    /// This covers both of the source's position-based overloads: the one that
    /// takes an explicit `size` and resizes, and the one that reuses
    /// `probability.size()`. Pass the existing length for the latter.
    ///
    /// A `size` of zero succeeds and returns an empty vector. The source's
    /// `normalApproximationHelper_` runs two `for (i = 0; i < size; ++i)` loops,
    /// neither of which executes: `gaussSum` stays zero but nothing ever divides
    /// by it, and the container is left empty without an exception. Both
    /// position-based overloads therefore succeed on an empty request, and so
    /// does this. The zero-sum refusal below applies only where the source would
    /// actually have divided.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `size` exceeds [`MAX_ITEMS`], and
    /// when `size` is non-zero and the densities sum to zero or to a non-finite
    /// value, which would otherwise make every entry NaN or infinite. The source
    /// divides unchecked.
    pub fn normal_approximation(&self, size: usize) -> Result<Vec<f64>> {
        if size > MAX_ITEMS {
            return Err(bad(
                "normal approximation size exceeds the supported length",
            ));
        }
        if size == 0 {
            return Ok(Vec::new());
        }
        let mut gauss_sum = 0.0;
        for index in 0..size {
            gauss_sum += self.normal_density_sqrt2pi(index as f64);
        }
        if gauss_sum == 0.0 || !gauss_sum.is_finite() {
            return Err(bad("normal approximation has no usable density mass"));
        }
        let mut probability = Vec::with_capacity(size);
        for index in 0..size {
            probability.push(self.normal_density_sqrt2pi(index as f64) / gauss_sum * self.sum);
        }
        Ok(probability)
    }

    /// Normal approximation sampled at explicit coordinates.
    ///
    /// As [`BasicStatistics::normal_approximation`], with `coordinates` in
    /// place of `0, 1, ..., size-1`; the result has one entry per coordinate.
    ///
    /// An empty coordinate vector succeeds and returns an empty vector, for the
    /// same reason as a zero `size` there: the source's helper never reaches a
    /// division.
    ///
    /// # Errors
    ///
    /// As [`BasicStatistics::normal_approximation`].
    pub fn normal_approximation_at(&self, coordinates: &[f64]) -> Result<Vec<f64>> {
        if coordinates.len() > MAX_ITEMS {
            return Err(bad(
                "normal approximation size exceeds the supported length",
            ));
        }
        if coordinates.is_empty() {
            return Ok(Vec::new());
        }
        let mut gauss_sum = 0.0;
        for &coordinate in coordinates {
            gauss_sum += self.normal_density_sqrt2pi(coordinate);
        }
        if gauss_sum == 0.0 || !gauss_sum.is_finite() {
            return Err(bad("normal approximation has no usable density mass"));
        }
        let mut probability = Vec::with_capacity(coordinates.len());
        for &coordinate in coordinates {
            probability.push(self.normal_density_sqrt2pi(coordinate) / gauss_sum * self.sum);
        }
        Ok(probability)
    }
}

impl std::fmt::Display for BasicStatistics {
    /// The source's debugging `operator<<`, field for field.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "BasicStatistics:  mean={}  variance={}  sum={}",
            self.mean, self.variance, self.sum
        )
    }
}
