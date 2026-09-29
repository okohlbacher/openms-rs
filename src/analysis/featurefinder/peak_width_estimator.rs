// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Rough estimation of the peak width at any m/z.
//!
//! Native equivalent of `FEATUREFINDER/PeakWidthEstimator.h`, ported as a
//! dependency of the multiplex clustering. See
//! `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md`.

use crate::kernel::MSExperiment;
use crate::processing::peak_picking::PeakBoundary;
use crate::processing::spline::{BSpline2d, BoundaryCondition};
use crate::{Error, Result};

/// Rough estimation of the peak width at m/z (source `PeakWidthEstimator`).
///
/// From the picked peaks of a data set and the peak boundaries
/// `PeakPickerHiRes` reports for them, a B-spline through (m/z, boundary width)
/// estimates the typical peak width at any m/z.
#[derive(Clone, Debug)]
pub struct PeakWidthEstimator {
    /// B-spline for peak width interpolation (source `bspline_`).
    bspline: BSpline2d,
    /// First m/z of the interpolation range (source `mz_min_`).
    mz_min: f64,
    /// Last m/z of the interpolation range (source `mz_max_`).
    mz_max: f64,
}

impl PeakWidthEstimator {
    /// Fit the width model to `exp_picked` and its `boundaries` (source
    /// constructor).
    ///
    /// Spectra and boundary lists are paired in order, and within a spectrum
    /// peaks and boundaries, each pairing stopping at the shorter side, as the
    /// source's loops do. A peak's width is `max - min` of its boundary. The
    /// interpolation range runs from the **first** collected m/z to the
    /// **last**, not from the smallest to the largest: the two agree for one
    /// sorted spectrum but not across spectra, and the source uses the first
    /// and last (`peaks_mz.front()`, `peaks_mz.back()`). The spline is a
    /// `BSpline2d` with wavelength `min(500, (mz_max - mz_min) / 2)`, zero second
    /// derivative at the ends and the source's `num_nodes` of 1 (automatic).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when no peak is paired with a boundary
    /// (the source reads `front()` of an empty vector), and "Unable to fit
    /// B-spline to data." when the B-spline cannot be set up (source
    /// `Exception::UnableToFit`), including every error [`BSpline2d::with_options`]
    /// reports.
    pub fn new(exp_picked: &MSExperiment, boundaries: &[Vec<PeakBoundary>]) -> Result<Self> {
        let mut peaks_mz = Vec::new();
        let mut peaks_width = Vec::new();
        for (spectrum, spectrum_boundaries) in exp_picked.spectra.iter().zip(boundaries) {
            for (peak, boundary) in spectrum.peaks.iter().zip(spectrum_boundaries) {
                peaks_mz.push(peak.mz);
                peaks_width.push(boundary.max - boundary.min);
            }
        }
        let (Some(&mz_min), Some(&mz_max)) = (peaks_mz.first(), peaks_mz.last()) else {
            return Err(Error::InvalidValue(
                "PeakWidthEstimator: no picked peak with a boundary".into(),
            ));
        };
        let unable = |detail: String| {
            Error::InvalidValue(format!("Unable to fit B-spline to data. {detail}"))
        };
        let bspline = BSpline2d::with_options(
            &peaks_mz,
            &peaks_width,
            500.0_f64.min((mz_max - mz_min) / 2.0),
            BoundaryCondition::ZeroSecond,
            1,
        )
        .map_err(|e| unable(e.to_string()))?;
        if !bspline.ok() {
            return Err(unable(String::new()));
        }
        Ok(Self {
            bspline,
            mz_min,
            mz_max,
        })
    }

    /// The estimated peak width at `mz` (source `getPeakWidth`). Outside the
    /// interpolation range the width at the nearer end is used, as the source
    /// notes.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] "Estimated peak width is negative." for a
    /// negative estimate (source `Exception::InvalidValue`), and the errors of
    /// [`BSpline2d::eval`] (a non-finite `mz`).
    pub fn peak_width(&self, mz: f64) -> Result<f64> {
        let at = if mz < self.mz_min {
            self.mz_min
        } else if mz > self.mz_max {
            self.mz_max
        } else {
            mz
        };
        let width = self.bspline.eval(at)?;
        if width < 0.0 {
            return Err(Error::InvalidValue(
                "Estimated peak width is negative.".into(),
            ));
        }
        Ok(width)
    }
}
