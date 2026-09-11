// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Owned masses and scaled integer weights from source OpenMS::ims::Weights.
//! This value utility permits zero weights; decomposition solvers impose their
//! own stricter alphabet rules. See `docs/IMS_WEIGHTS_SUPPORT.md`.

use crate::{Error, Result};
use std::fmt::Write;

pub const MAX_IMS_WEIGHTS: usize = 1_000_000;
pub const MAX_IMS_WEIGHTS_WORK: usize = 50_000_000;
pub const MAX_IMS_WEIGHTS_OUTPUT_BYTES: usize = 8 * 1024 * 1024;
// Exclusive: unlike u64::MAX as f64, the power of two is an exact boundary.
const U64_END: f64 = 18_446_744_073_709_551_616.0;

/// An immutable-length, owned alphabet and its current integer scaling.
/// `new`/`Default` use `None` for the source default's uninitialized precision.
/// Standard cloning is bounded by the private one-million-element invariant.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IMSWeights {
    masses: Vec<f64>,
    weights: Vec<u64>,
    precision: Option<f64>,
}

impl IMSWeights {
    pub fn new() -> Self {
        Self::default()
    }

    /// Preserve finite signed original masses and finite precision whenever the
    /// literal source floor(mass/precision+0.5) is representable as u64. Zero
    /// precision is therefore accepted only for an empty alphabet.
    pub fn from_masses(masses: &[f64], precision: f64) -> Result<Self> {
        if masses.len() > MAX_IMS_WEIGHTS {
            return Err(invalid("alphabet element limit exceeded"));
        }
        let mut work = Work::default();
        let weights = scaled(masses, precision, &mut work)?;
        work.consume(masses.len())?;
        let mut owned = vector(masses.len())?;
        owned.extend_from_slice(masses);
        Ok(Self {
            masses: owned,
            weights,
            precision: Some(precision),
        })
    }
    pub fn len(&self) -> usize {
        self.weights.len()
    }
    pub fn is_empty(&self) -> bool {
        self.weights.is_empty()
    }
    pub fn weights(&self) -> &[u64] {
        &self.weights
    }
    pub fn masses(&self) -> &[f64] {
        &self.masses
    }
    pub fn precision(&self) -> Option<f64> {
        self.precision
    }
    /// Checked equivalent of both source getWeight and operator[].
    pub fn weight(&self, index: usize) -> Result<u64> {
        self.weights
            .get(index)
            .copied()
            .ok_or_else(|| invalid("weight index out of range"))
    }
    pub fn alphabet_mass(&self, index: usize) -> Result<f64> {
        self.masses
            .get(index)
            .copied()
            .ok_or_else(|| invalid("mass index out of range"))
    }
    pub fn back(&self) -> Result<u64> {
        self.weights
            .last()
            .copied()
            .ok_or_else(|| invalid("empty alphabet has no last weight"))
    }

    /// Recompute from the original masses, including after GCD scaling. Invalid
    /// precision or a later overflowing weight leaves all old state unchanged.
    pub fn set_precision(&mut self, precision: f64) -> Result<()> {
        let weights = scaled(&self.masses, precision, &mut Work::default())?;
        self.weights = weights;
        self.precision = Some(precision);
        Ok(())
    }
    /// Swap both members of the mass/weight pair after checking both indices.
    pub fn swap(&mut self, first: usize, second: usize) -> Result<()> {
        if first >= self.len() || second >= self.len() {
            return Err(invalid("swap index out of range"));
        }
        self.weights.swap(first, second);
        self.masses.swap(first, second);
        Ok(())
    }
    /// Ordered f64 sum of original mass * u32 count. Finite signed results are
    /// valid; the current precision and integer weights do not affect this sum.
    pub fn parent_mass(&self, decomposition: &[u32]) -> Result<f64> {
        if decomposition.len() != self.len() {
            return Err(Error::InvalidValue(format!(
                "The passed decomposition has the wrong size. Expected {} but got {}.",
                self.len(),
                decomposition.len()
            )));
        }
        let mut work = Work::default();
        work.consume(self.len())?;
        let mut parent_mass = 0.0;
        for (&mass, &count) in self.masses.iter().zip(decomposition) {
            parent_mass += mass * f64::from(count);
        }
        finite(parent_mass, "parent mass is not finite")
    }

    /// GCD scaling without rerounding. Source returns true for exactly two
    /// coprime entries (despite no numerical change), but false for a coprime
    /// alphabet of three or more. All-zero alphabets of size>=2 are errors.
    pub fn divide_by_gcd(&mut self) -> Result<bool> {
        self.divide_with_work(&mut Work::default())
    }
    fn divide_with_work(&mut self, work: &mut Work) -> Result<bool> {
        if self.len() < 2 {
            return Ok(false);
        }
        let mut divisor = gcd(self.weights[0], self.weights[1], work)?;
        for &weight in &self.weights[2..] {
            divisor = gcd(divisor, weight, work)?;
            if divisor == 1 {
                return Ok(false);
            }
        }
        if divisor == 0 {
            return Err(invalid("all-zero weights have no positive GCD"));
        }
        let precision = finite(
            self.precision
                .ok_or_else(|| invalid("precision is unset"))?
                * divisor as f64,
            "GCD-scaled precision is not finite",
        )?;
        work.consume(self.len())?;
        self.precision = Some(precision);
        for weight in &mut self.weights {
            *weight /= divisor;
        }
        Ok(true)
    }
    pub fn min_rounding_error(&self) -> Result<f64> {
        self.rounding_error_bound(true)
    }
    pub fn max_rounding_error(&self) -> Result<f64> {
        self.rounding_error_bound(false)
    }
    fn rounding_error_bound(&self, minimum: bool) -> Result<f64> {
        if self.is_empty() {
            return Ok(0.0);
        }
        let precision = self
            .precision
            .ok_or_else(|| invalid("precision is unset"))?;
        Work::default().consume(self.len())?;
        let mut bound = 0.0;
        for (&weight, &mass) in self.weights.iter().zip(&self.masses) {
            let error = (precision * weight as f64 - mass) / mass;
            // Literal comparisons retain source handling of 0/0 (ignored), and
            // ignore an overflowing error in the opposite direction. Only a
            // nonfinite returned bound is a native checked error.
            if (minimum && error < 0.0 && error < bound)
                || (!minimum && error > 0.0 && error > bound)
            {
                bound = error;
            }
        }
        finite(bound, "rounding-error bound is not finite")
    }
    /// Source stream representation: every integer weight followed by '\n'.
    /// Complete length/work preflight occurs before allocating the string.
    pub fn to_text(&self) -> Result<String> {
        let mut work = Work::default();
        work.consume(self.len())?;
        let mut length = 0usize;
        for &weight in &self.weights {
            let digits = if weight == 0 {
                1
            } else {
                weight.ilog10() as usize + 1
            };
            length = length
                .checked_add(digits + 1)
                .filter(|&n| n <= MAX_IMS_WEIGHTS_OUTPUT_BYTES)
                .ok_or_else(|| invalid("formatted output byte limit exceeded"))?;
        }
        work.consume(length)?;
        let mut text = String::new();
        text.try_reserve_exact(length)
            .map_err(|_| invalid("formatted output allocation failed"))?;
        for weight in &self.weights {
            writeln!(text, "{weight}").expect("String writing cannot fail");
        }
        Ok(text)
    }
}

fn scaled(masses: &[f64], precision: f64, work: &mut Work) -> Result<Vec<u64>> {
    finite(precision, "precision must be finite")?;
    work.consume(masses.len())?;
    let mut weights = vector(masses.len())?;
    for &mass in masses {
        finite(mass, "original mass must be finite")?;
        let value = (mass / precision + 0.5).floor();
        if !value.is_finite() || !(0.0..U64_END).contains(&value) {
            return Err(invalid("rounded weight is outside u64 range"));
        }
        weights.push(value as u64);
    }
    Ok(weights)
}
fn vector<T>(length: usize) -> Result<Vec<T>> {
    let mut vector = Vec::new();
    vector
        .try_reserve_exact(length)
        .map_err(|_| invalid("weight vector allocation failed"))?;
    Ok(vector)
}
fn gcd(mut a: u64, mut b: u64, work: &mut Work) -> Result<u64> {
    work.consume(1)?;
    while b != 0 {
        work.consume(1)?;
        (a, b) = (b, a % b);
    }
    Ok(a)
}
fn finite(value: f64, message: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(message))
    }
}
fn invalid(message: &'static str) -> Error {
    Error::InvalidValue(format!("IMS weights: {message}"))
}
struct Work {
    remaining: usize,
}
impl Default for Work {
    fn default() -> Self {
        Self {
            remaining: MAX_IMS_WEIGHTS_WORK,
        }
    }
}
impl Work {
    fn consume(&mut self, count: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(count)
            .ok_or_else(|| invalid("work limit exceeded"))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn late_shared_work_failure_precedes_gcd_commit() {
        let mut weights = IMSWeights::from_masses(&[30.0, 50.0, 80.0], 1.0).unwrap();
        let old = weights.clone();
        let mut work = Work { remaining: 10 };
        assert!(weights.divide_with_work(&mut work).is_err());
        assert_eq!(weights, old);
    }
}
