// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Rough estimation of the peak width at an m/z.
//!
//! Ports `FEATUREFINDER/PeakWidthEstimator.h` and its `.cpp`. See
//! `docs/FEATUREFINDER_MODELS_SUPPORT.md`.

use crate::Result;
use crate::error::Error;
use crate::kernel::MSExperiment;
use crate::processing::peak_picking::PeakBoundary;
use crate::processing::spline::b_spline::{BSpline2d, BoundaryCondition};

/// The typical peak width at arbitrary m/z, interpolated by a smoothing
/// B-spline through the widths of picked peaks (source `PeakWidthEstimator`).
///
/// The widths come from the peak boundaries [`crate::processing::peak_picking::PeakPickerHiRes`]
/// reports: `max - min` of each peak's boundary.
#[derive(Clone, Debug)]
pub struct PeakWidthEstimator {
    /// Source `bspline_`.
    bspline: BSpline2d,
    /// m/z of the first peak (source `mz_min_`).
    mz_min: f64,
    /// m/z of the last peak (source `mz_max_`).
    mz_max: f64,
}

impl PeakWidthEstimator {
    /// Fit the width spline to the picked peaks of `exp_picked` and their
    /// `boundaries` (source constructor).
    ///
    /// Spectra and boundary lists are paired by position, and within them
    /// peaks and boundaries, each pair stopping at the shorter side, as the
    /// source's two parallel iterator loops do. (The picker reports boundaries
    /// for picked spectra only, so with spectra it skipped the pairing
    /// shifts; the source has the same behaviour.) The spline is a
    /// [`BSpline2d`] with cutoff wavelength `min(500, (mz_max - mz_min) / 2)`,
    /// second-derivative-zero boundaries and one node, where `mz_min` and
    /// `mz_max` are the m/z of the *first and last* peak visited.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] with the source's message "Unable to
    /// fit B-spline to data." when the spline cannot be set up, including for
    /// no peaks at all (where the source reads `front()` of an empty vector),
    /// non-finite values, or a wavelength the spline refuses.
    pub fn new(exp_picked: &MSExperiment, boundaries: &[Vec<PeakBoundary>]) -> Result<Self> {
        let mut peaks_mz = Vec::new();
        let mut peaks_width = Vec::new();
        for (spectrum, spectrum_boundaries) in exp_picked.spectra.iter().zip(boundaries) {
            for (peak, boundary) in spectrum.peaks.iter().zip(spectrum_boundaries) {
                peaks_mz.push(peak.mz);
                peaks_width.push(boundary.max - boundary.min);
            }
        }
        let unable = |detail: String| {
            Error::InvalidValue(format!(
                "UnableToFit: Unable to fit B-spline to data. {detail}"
            ))
        };
        let (Some(&mz_min), Some(&mz_max)) = (peaks_mz.first(), peaks_mz.last()) else {
            return Err(unable("There are no picked peaks.".into()));
        };
        let wavelength = 500.0f64.min((mz_max - mz_min) / 2.0);
        let bspline = BSpline2d::with_options(
            &peaks_mz,
            &peaks_width,
            wavelength,
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

    /// The estimated peak width at `mz` (source `getPeakWidth`).
    ///
    /// Outside the interpolation range the width at the nearer end (the first
    /// or last peak's m/z) is used.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] "Estimated peak width is negative."
    /// when the spline gives a negative width (the source's
    /// `Exception::InvalidValue`), and when the spline cannot be evaluated.
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
