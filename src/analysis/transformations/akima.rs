// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// The algorithm follows gte::IntpAkima1 and gte::IntpAkimaNonuniform1 of the
// Geometric Tools Engine (David Eberly, Boost Software License 1.0), which the
// pinned OpenMS vendors at src/openms/extern/GTE/Mathematics/.

//! The Akima spline behind the `akima` option of
//! `TransformationModelInterpolated`.
//!
//! The source's `AkimaInterpolator` (a private class in
//! `TransformationModelInterpolated.cpp`) wraps
//! `gte::IntpAkimaNonuniform1<double>` from the vendored Geometric Tools
//! Engine. This module ports that template and the part of its base
//! `gte::IntpAkima1` it reaches: the slope table, the derivative rule and the
//! cubic pieces, operation for operation. No crate was used; see
//! `docs/THIRD_PARTY_CRATE_DECISIONS.md` (the Akima entry) for the crates that
//! were considered and why none reproduces the GTE end conditions.

use super::{bad, finite};
use crate::Result;

/// Piecewise-cubic Akima interpolant over strictly increasing abscissae.
#[derive(Clone, Debug)]
pub(super) struct Akima {
    /// The abscissae, strictly increasing.
    x: Vec<f64>,
    /// `c0 + c1*t + c2*t^2 + c3*t^3` for the interval starting at `x[i]`,
    /// with `t = x - x[i]`.
    polynomials: Vec<[f64; 4]>,
}

impl Akima {
    /// Build the interpolant, as the `IntpAkimaNonuniform1` constructor does.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for fewer than three points,
    /// arrays of different length, abscissae that are not strictly increasing
    /// (both `LogAssert` failures in GTE, which throw `std::runtime_error`),
    /// or a non-finite slope or coefficient.
    pub(super) fn new(x: &[f64], f: &[f64]) -> Result<Self> {
        let n = x.len();
        if n < 3 || f.len() != n {
            return Err(bad("Akima interpolation needs at least three points"));
        }
        if x.windows(2).any(|pair| pair[1] <= pair[0]) {
            return Err(bad(
                "Akima interpolation needs strictly increasing x values",
            ));
        }
        // slope[i + 2] is the chord slope of interval i; two extrapolated
        // slopes pad each end.
        let mut slope = vec![0.0; n + 3];
        for i in 0..n - 1 {
            slope[i + 2] = finite((f[i + 1] - f[i]) / (x[i + 1] - x[i]))?;
        }
        slope[1] = 2.0 * slope[2] - slope[3];
        slope[0] = 2.0 * slope[1] - slope[2];
        slope[n + 1] = 2.0 * slope[n] - slope[n - 1];
        slope[n + 2] = 2.0 * slope[n + 1] - slope[n];
        let derivatives = (0..n)
            .map(|i| finite(derivative(&slope[i..i + 4])))
            .collect::<Result<Vec<_>>>()?;
        let mut polynomials = Vec::with_capacity(n - 1);
        for i in 0..n - 1 {
            let (f0, f1) = (f[i], f[i + 1]);
            let (d0, d1) = (derivatives[i], derivatives[i + 1]);
            let df = f1 - f0;
            let dx = x[i + 1] - x[i];
            let dx2 = dx * dx;
            let dx3 = dx2 * dx;
            polynomials.push([
                f0,
                d0,
                finite((3.0 * df - dx * (d1 + 2.0 * d0)) / dx2)?,
                finite((dx * (d0 + d1) - 2.0 * df) / dx3)?,
            ]);
        }
        Ok(Self {
            x: x.to_vec(),
            polynomials,
        })
    }

    /// Evaluate at `value`, clamped to the abscissa range first as
    /// `IntpAkima1::operator()` does.
    pub(super) fn eval(&self, value: f64) -> f64 {
        let first = self.x[0];
        let last = self.x[self.x.len() - 1];
        // std::min(std::max(x, xmin), xmax)
        let value = if value < first { first } else { value };
        let value = if last < value { last } else { value };
        // IntpAkimaNonuniform1::Lookup: a linear scan for the first interval
        // whose right end lies beyond the value; the last interval otherwise.
        let mut index = self.x.len() - 2;
        for i in 0..self.x.len() - 1 {
            if value < self.x[i + 1] {
                index = i;
                break;
            }
        }
        let t = value - self.x[index];
        let c = &self.polynomials[index];
        c[0] + t * (c[1] + t * (c[2] + t * c[3]))
    }
}

/// `IntpAkima1::ComputeDerivative` over four consecutive slopes.
fn derivative(slope: &[f64]) -> f64 {
    if slope[1] != slope[2] {
        if slope[0] != slope[1] {
            if slope[2] != slope[3] {
                let ad0 = (slope[3] - slope[2]).abs();
                let ad1 = (slope[0] - slope[1]).abs();
                (ad0 * slope[1] + ad1 * slope[2]) / (ad0 + ad1)
            } else {
                slope[2]
            }
        } else if slope[2] != slope[3] {
            slope[1]
        } else {
            0.5 * (slope[1] + slope[2])
        }
    } else {
        slope[1]
    }
}
