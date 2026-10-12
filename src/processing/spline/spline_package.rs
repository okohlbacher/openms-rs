// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The spline fit of one package of profile data points.
//!
//! Native equivalent of `PROCESSING/MISC/SplinePackage.h`. See
//! `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md`, which covers this header as a
//! dependency of the multiplex feature finder.

use super::cubic::CubicSpline2d;
use crate::{Error, Result};

/// The spline fit of a single package of data points (source
/// `SplinePackage`), the building block of
/// [`SplineInterpolatedPeaks`](super::spline_interpolated_peaks::SplineInterpolatedPeaks).
///
/// Data points in MS spectra (or chromatograms) are often not equidistant in
/// m/z (or RT) but come in packages separated by wide zero-intensity ranges.
/// A package holds the natural cubic spline through one such set of points.
#[derive(Clone, Debug)]
pub struct SplinePackage {
    /// Smallest position for which the fit is valid (source `pos_min_`).
    pos_min: f64,
    /// Largest position for which the fit is valid (source `pos_max_`).
    pos_max: f64,
    /// Sensible step width to scan through the package: the mean spacing of
    /// the input positions (source `pos_step_width_`). Each navigator rescales
    /// it.
    pos_step_width: f64,
    /// The interpolating spline (source `spline_`).
    spline: CubicSpline2d,
}

impl SplinePackage {
    /// Fit the spline through `pos` and `intensity` (source constructor).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the vectors differ in size or hold
    /// fewer than two points (source `Exception::IllegalArgument`, "m/z (or
    /// RT) and intensity vectors either not of the same size or too short."),
    /// and the errors of [`CubicSpline2d::new`], which also refuses positions
    /// that are not strictly increasing or not finite, where the source builds
    /// a spline of NaN coefficients.
    pub fn new(pos: &[f64], intensity: &[f64]) -> Result<Self> {
        if !(pos.len() == intensity.len() && pos.len() > 1) {
            return Err(Error::InvalidValue(
                "m/z (or RT) and intensity vectors either not of the same size or too short."
                    .into(),
            ));
        }
        let spline = CubicSpline2d::new(pos, intensity)?;
        let pos_min = pos[0];
        let pos_max = pos[pos.len() - 1];
        Ok(Self {
            pos_min,
            pos_max,
            pos_step_width: (pos_max - pos_min) / (pos.len() - 1) as f64,
            spline,
        })
    }

    /// The smallest position for which the fit is valid (source `getPosMin`).
    pub fn pos_min(&self) -> f64 {
        self.pos_min
    }

    /// The largest position for which the fit is valid (source `getPosMax`).
    pub fn pos_max(&self) -> f64 {
        self.pos_max
    }

    /// A sensible step width for the package (source `getPosStepWidth`).
    pub fn pos_step_width(&self) -> f64 {
        self.pos_step_width
    }

    /// Whether `pos` lies in `[pos_min, pos_max]` (source `isInPackage`).
    pub fn is_in_package(&self, pos: f64) -> bool {
        pos >= self.pos_min && pos <= self.pos_max
    }

    /// The interpolated intensity at `pos`, clamped below at zero, or zero
    /// outside the package (source `eval`). The clamp is `std::max(0.0, v)`,
    /// which yields `0.0` for a NaN `v`.
    ///
    /// # Errors
    ///
    /// The overflow error of [`CubicSpline2d::eval`], which the source does not
    /// check.
    pub fn eval(&self, pos: f64) -> Result<f64> {
        if self.is_in_package(pos) {
            let value = self.spline.eval(pos)?;
            Ok(if 0.0 < value { value } else { 0.0 })
        } else {
            Ok(0.0)
        }
    }
}
