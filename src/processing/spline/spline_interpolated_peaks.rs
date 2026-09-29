// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Spline interpolation of MS1 spectra and chromatograms, package by package.
//!
//! Native equivalent of `PROCESSING/MISC/SplineInterpolatedPeaks.h`. See
//! `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md`, which covers this header as a
//! dependency of the multiplex feature finder.

use super::spline_package::SplinePackage;
use crate::kernel::{MSChromatogram, MSSpectrum};
use crate::{Error, Result};

/// The step width scaling source `getNavigator` defaults to.
pub const DEFAULT_NAVIGATOR_SCALING: f64 = 0.7;

/// Ceiling on the input data points, checked before any work vector is
/// allocated. The source has no bound.
pub const MAX_POINTS: usize = 100_000_000;

/// Spline interpolation of an MS1 spectrum or a chromatogram (source
/// `SplineInterpolatedPeaks`).
///
/// A set of splines, each interpolating the data in one m/z (or RT) range.
/// Between these ranges there are no raw data points and the intensity is
/// zero: a spline through non-equidistant data swings wildly where there are
/// no points and cannot be used there, so such regions are taken to be zero.
///
/// Construction first drops zero-intensity points whose neighbours on both
/// sides are zero too, then the "Thermo bug" zeros — a zero following its
/// predecessor within 1/50 of the previous spacing — and finally splits the
/// points wherever a spacing exceeds twice the previous one. Every part with at
/// least two points becomes a [`SplinePackage`].
#[derive(Clone, Debug)]
pub struct SplineInterpolatedPeaks {
    /// Smallest input position (source `pos_min_`).
    pos_min: f64,
    /// Largest input position (source `pos_max_`).
    pos_max: f64,
    /// The packages, in position order (source `packages_`).
    packages: Vec<SplinePackage>,
}

impl SplineInterpolatedPeaks {
    /// Interpolate the points `(pos[i], intensity[i])` (source constructor
    /// taking two vectors). The positions are assumed sorted, as the source
    /// notes.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the vectors differ in size or hold
    /// fewer than three points (source `Exception::IllegalArgument`, "m/z and
    /// intensity vectors either not of the same size or too short."), when a
    /// package's positions are not strictly increasing and finite (see
    /// [`SplinePackage::new`]), and [`Error::InvalidRange`] for more than
    /// [`MAX_POINTS`] points.
    ///
    /// When fewer than two points survive the removal of surplus zeros (all
    /// intensities zero, for instance), the source reads past its work vector;
    /// this port builds no package, so [`Self::navigator`] then fails as it
    /// does for any spectrum without packages.
    pub fn new(pos: &[f64], intensity: &[f64]) -> Result<Self> {
        if !(pos.len() == intensity.len() && pos.len() > 2) {
            return Err(Error::InvalidValue(
                "m/z and intensity vectors either not of the same size or too short.".into(),
            ));
        }
        if pos.len() > MAX_POINTS {
            return Err(Error::InvalidRange(format!(
                "SplineInterpolatedPeaks: more than {MAX_POINTS} data points"
            )));
        }
        let n = pos.len();
        // start a new package if delta m/z is greater than new_package times previous one
        let new_package = 2.0;

        // remove unnecessary zeros, i.e. zero intensity data points with zeros to the left and right
        let mut pos_slim1 = Vec::with_capacity(n);
        let mut intensity_slim1 = Vec::with_capacity(n);
        if intensity[0] != 0.0 || intensity[1] != 0.0 {
            pos_slim1.push(pos[0]);
            intensity_slim1.push(intensity[0]);
        }
        let mut current_zero = intensity[0] == 0.0;
        let mut next_zero = intensity[1] == 0.0;
        for i in 1..n - 1 {
            let last_zero = current_zero;
            current_zero = next_zero;
            next_zero = intensity[i + 1] == 0.0;
            if !last_zero || !current_zero || !next_zero {
                pos_slim1.push(pos[i]);
                intensity_slim1.push(intensity[i]);
            }
        }
        if intensity[n - 1] != 0.0 || intensity[n - 2] != 0.0 {
            pos_slim1.push(pos[n - 1]);
            intensity_slim1.push(intensity[n - 1]);
        }

        let mut packages = Vec::new();
        if pos_slim1.len() >= 2 {
            // remove Thermo bug zeros: a zero quickly following a data point
            let scaling_thermo_bug = 1.0 / 50.0;
            let mut pos_slim2 = Vec::with_capacity(pos_slim1.len());
            let mut intensity_slim2 = Vec::with_capacity(pos_slim1.len());
            pos_slim2.extend_from_slice(&pos_slim1[..2]);
            intensity_slim2.extend_from_slice(&intensity_slim1[..2]);
            for i in 2..pos_slim1.len() {
                if intensity_slim1[i] == 0.0
                    && (pos_slim1[i] - pos_slim1[i - 1]) < (pos_slim1[i - 1] - pos_slim1[i - 2]) * scaling_thermo_bug
                {
                    continue;
                }
                pos_slim2.push(pos_slim1[i]);
                intensity_slim2.push(intensity_slim1[i]);
            }

            // subdivide spectrum into packages
            let mut start_package = vec![true, false];
            for i in 2..pos_slim2.len() {
                start_package.push(
                    (pos_slim2[i] - pos_slim2[i - 1]) / (pos_slim2[i - 1] - pos_slim2[i - 2]) > new_package,
                );
            }

            // fill the packages
            let mut pos_package = Vec::new();
            let mut intensity_package = Vec::new();
            for i in 0..pos_slim2.len() {
                if start_package[i] && i > 0 {
                    if intensity_package.len() > 1 {
                        packages.push(SplinePackage::new(&pos_package, &intensity_package)?);
                    }
                    pos_package.clear();
                    intensity_package.clear();
                }
                pos_package.push(pos_slim2[i]);
                intensity_package.push(intensity_slim2[i]);
            }
            if intensity_package.len() > 1 {
                packages.push(SplinePackage::new(&pos_package, &intensity_package)?);
            }
        }
        Ok(Self {
            pos_min: pos[0],
            pos_max: pos[n - 1],
            packages,
        })
    }

    /// Interpolate a spectrum's peaks, m/z against intensity widened to
    /// `f64` (source constructor taking an `MSSpectrum`).
    ///
    /// # Errors
    ///
    /// As [`Self::new`].
    pub fn from_spectrum(raw_spectrum: &MSSpectrum) -> Result<Self> {
        let pos: Vec<f64> = raw_spectrum.peaks.iter().map(|p| p.mz).collect();
        let intensity: Vec<f64> = raw_spectrum.peaks.iter().map(|p| f64::from(p.intensity)).collect();
        Self::new(&pos, &intensity)
    }

    /// Interpolate a chromatogram's peaks, RT against intensity widened to
    /// `f64` (source constructor taking an `MSChromatogram`).
    ///
    /// # Errors
    ///
    /// As [`Self::new`].
    pub fn from_chromatogram(raw_chromatogram: &MSChromatogram) -> Result<Self> {
        let pos: Vec<f64> = raw_chromatogram.peaks.iter().map(|p| p.rt).collect();
        let intensity: Vec<f64> = raw_chromatogram.peaks.iter().map(|p| f64::from(p.intensity)).collect();
        Self::new(&pos, &intensity)
    }

    /// The smallest input position (source `getPosMin`).
    pub fn pos_min(&self) -> f64 {
        self.pos_min
    }

    /// The largest input position (source `getPosMax`).
    pub fn pos_max(&self) -> f64 {
        self.pos_max
    }

    /// The number of spline packages (source `size`). Check it before asking
    /// for a [`Navigator`]: without packages, [`Self::navigator`] fails.
    pub fn size(&self) -> usize {
        self.packages.len()
    }

    /// The packages, in position order.
    pub fn packages(&self) -> &[SplinePackage] {
        &self.packages
    }

    /// A navigator over the packages whose step widths are scaled by `scaling`
    /// (source `getNavigator`, whose default scaling is
    /// [`DEFAULT_NAVIGATOR_SCALING`]). The navigator borrows this value, which
    /// the source can only ask its callers to keep alive.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] "packages_ must not be empty" when there
    /// is no package (source `Exception::InvalidSize`).
    pub fn navigator(&self, scaling: f64) -> Result<Navigator<'_>> {
        if self.packages.is_empty() {
            return Err(Error::InvalidValue("packages_ must not be empty".into()));
        }
        Ok(Navigator {
            packages: &self.packages,
            last_package: 0,
            pos_max: self.pos_max,
            pos_step_width_scaling: scaling,
        })
    }
}

/// Iterator-like access to the spline packages (source
/// `SplineInterpolatedPeaks::Navigator`).
///
/// It remembers the last package it used, so that successive nearby queries
/// start their search there. The source's default constructor (for pyOpenMS)
/// leaves the package pointer null, so every call on it is undefined; it is
/// not ported.
#[derive(Clone, Debug)]
pub struct Navigator<'a> {
    /// The packages (source `packages_`); never empty.
    packages: &'a [SplinePackage],
    /// Index of the package last accessed (source `last_package_`).
    last_package: usize,
    /// Largest position of the spectrum (source `pos_max_`).
    pos_max: f64,
    /// Scaling of the packages' step widths (source
    /// `pos_step_width_scaling_`). A step somewhat smaller than the mean raw
    /// spacing is often advantageous.
    pos_step_width_scaling: f64,
}

impl Navigator<'_> {
    /// The spline-interpolated intensity at `pos`, zero between and outside
    /// the packages (source `Navigator::eval`).
    ///
    /// The search starts at the last package used and walks left or right; a
    /// position before every package leaves the remembered package unchanged,
    /// as in the source.
    ///
    /// # Errors
    ///
    /// The errors of [`SplinePackage::eval`].
    pub fn eval(&mut self, pos: f64) -> Result<f64> {
        let packages = self.packages;
        if pos < packages[self.last_package].pos_min() {
            // look left
            for i in (0..=self.last_package).rev() {
                if pos > packages[i].pos_max() {
                    self.last_package = i;
                    return Ok(0.0);
                }
                if pos >= packages[i].pos_min() {
                    self.last_package = i;
                    return packages[i].eval(pos);
                }
            }
        } else {
            // look right
            for (i, package) in packages.iter().enumerate().skip(self.last_package) {
                if pos < package.pos_min() {
                    self.last_package = i;
                    return Ok(0.0);
                }
                if pos <= package.pos_max() {
                    self.last_package = i;
                    return package.eval(pos);
                }
            }
        }
        Ok(0.0)
    }

    /// The next sensible position for scanning through the data (source
    /// `Navigator::getNextPos`).
    ///
    /// Inside a package, the position advances by the package's step width
    /// times the scaling; a step that would leave the package jumps to the
    /// start of the next package, and past the last package the result is the
    /// spectrum's largest position. A position in a gap between packages moves
    /// to the start of the package to its right, and one before the first
    /// package to the first package's start.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a NaN position, on which the
    /// source's package search never terminates.
    pub fn next_pos(&mut self, pos: f64) -> Result<f64> {
        if pos.is_nan() {
            return Err(Error::InvalidValue(
                "SplineInterpolatedPeaks::Navigator: position is NaN".into(),
            ));
        }
        let packages = self.packages;
        let max_index = packages.len() - 1;
        let mut i = self.last_package;
        // find correct package
        while !packages[i].is_in_package(pos) {
            if pos < packages[i].pos_min() {
                if i == 0 {
                    self.last_package = 0;
                    return Ok(packages[0].pos_min());
                }
                i -= 1;
                // m/z in the gap?
                if pos > packages[i].pos_max() {
                    self.last_package = i + 1;
                    return Ok(packages[i + 1].pos_min());
                }
            } else {
                if i == max_index {
                    self.last_package = max_index;
                    return Ok(self.pos_max);
                }
                i += 1;
                // m/z in the gap?
                if pos < packages[i].pos_min() {
                    self.last_package = i;
                    return Ok(packages[i].pos_min());
                }
            }
        }
        // find m/z in the package
        let step = self.pos_step_width_scaling * packages[i].pos_step_width();
        if pos + step > packages[i].pos_max() {
            // The next step leaves the package: move to the package to the right.
            if i == max_index {
                self.last_package = max_index;
                return Ok(self.pos_max);
            }
            self.last_package = i + 1;
            Ok(packages[i + 1].pos_min())
        } else {
            self.last_package = i;
            Ok(pos + step)
        }
    }
}
