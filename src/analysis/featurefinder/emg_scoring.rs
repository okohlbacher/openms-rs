// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Scoring of chromatographic elution peaks by an exponentially modified
//! Gaussian fit.
//!
//! Ports `FEATUREFINDER/EmgScoring.h` (all members are inline; `EmgScoring.cpp`
//! is empty). The class carries "the original ideas from
//! FeatureFinderAlgorithmMRM" as an interface for scoring chromatographic
//! peaks. See `docs/FEATUREFINDER_MODELS_SUPPORT.md`.

use crate::Result;
use crate::error::Error;
use crate::kernel::geometry::Point2D;
use crate::kernel::mrm::{MRMFeature, MRMTransitionGroup};
use crate::kernel::{MSChromatogram, MSSpectrum, Peak1D};
use crate::param::Param;
use crate::processing::SpectrumFilter;
use crate::processing::smoothing::{GaussFilter, GaussianWidth};

use super::base_model::Model1D;
use super::fitter1d::Fitter1D;
use super::fitters::EmgFitter1D;

/// Scores an elution peak by the quality of an [`EmgFitter1D`] fit (source
/// `EmgScoring`).
#[derive(Clone, Debug, Default)]
pub struct EmgScoring {
    /// The parameters handed to every [`EmgFitter1D`] (source
    /// `fitter_emg1D_params_`); unspecified ones stay at their defaults.
    fitter_emg1d_params: Param,
}

impl EmgScoring {
    /// A scorer with default fitter parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Overwrite the parameters for the EMG fitter; unspecified parameters
    /// stay default (source `setFitterParam`). See [`EmgScoring::defaults`].
    pub fn set_fitter_param(&mut self, param: Param) {
        self.fitter_emg1d_params = param;
    }

    /// The default parameters of the EMG fit (source `getDefaults`, which
    /// returns `EmgFitter1D().getDefaults()`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults(&self) -> Result<Param> {
        EmgFitter1D::defaults()
    }

    /// The elution profile fit score of a peak group: the mean of
    /// [`EmgScoring::elution_model_fit`] (without smoothing) over the first
    /// convex hull of the feature that `mrmfeature` holds for each
    /// chromatogram's native ID, in chromatogram order (source template
    /// `calcElutionFitScore`, here for chromatogram groups).
    ///
    /// An aborted fit contributes -1, "just the lowest pearson correlation to
    /// a fit that you can have" (a source TODO considers penalizing it more).
    /// A group without chromatograms yields NaN (`0 / 0`), as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when a feature has no convex hull
    /// (the source's debug-only precondition asks for exactly one and release
    /// builds read `getConvexHulls()[0]`), the error of
    /// [`MRMFeature::feature`] for a native ID without a feature (the source's
    /// `feature_map_[key]` then silently returns the first feature), and what
    /// the fit refuses.
    pub fn calc_elution_fit_score<T>(
        &self,
        mrmfeature: &MRMFeature,
        transition_group: &MRMTransitionGroup<MSChromatogram, T>,
    ) -> Result<f64> {
        let mut avg_score = 0.0f64;
        let smooth_data = false;
        let chromatograms = transition_group.chromatograms();
        for chromatogram in chromatograms {
            let feature = mrmfeature.feature(&chromatogram.native_id)?;
            let hull = feature.convex_hulls.first().ok_or_else(|| {
                Error::MissingInformation(format!(
                    "the feature of '{}' has no convex hull",
                    chromatogram.native_id
                ))
            })?;
            avg_score += self.elution_model_fit(&hull.hull_points(), smooth_data)?;
        }
        avg_score /= chromatograms.len() as f64;
        Ok(avg_score)
    }

    /// The elution model fit score of one section of hull points `(rt,
    /// intensity)`: the quality of an [`EmgFitter1D`] fit to the points padded
    /// by three zero-intensity points on either side at the mean RT spacing,
    /// optionally Gaussian-smoothed first (source `elutionModelFit`, "Fxn from
    /// FeatureFinderAlgorithmMRM").
    ///
    /// Fewer than two points score -1 ("we need at least 2 datapoints in order
    /// to create a fit"). The quality is the fitter's `float`, widened.
    ///
    /// # Errors
    ///
    /// Returns what smoothing (a non-positive mean spacing makes no Gaussian
    /// width) and the fit refuse.
    pub fn elution_model_fit(&self, current_section: &[Point2D], smooth_data: bool) -> Result<f64> {
        if current_section.len() < 2 {
            return Ok(-1.0);
        }
        let data_to_fit = Self::prepare_fit(current_section, smooth_data)?;
        let (quality, _model) = self.fit_rt(&data_to_fit)?;
        Ok(f64::from(quality))
    }

    /// Fit the RT profile with an [`EmgFitter1D`] configured by the stored
    /// parameters (source `fitRT_`). The fitter maps a NaN quality to -1.
    fn fit_rt(&self, rt_input_data: &[Peak1D]) -> Result<(f32, Model1D)> {
        let fitter = EmgFitter1D::from_param(&self.fitter_emg1d_params)?;
        fitter.fit1d(rt_input_data)
    }

    /// The points as peaks (position = RT, intensity narrowed to `float`),
    /// with three zero peaks appended and three prepended at the average RT
    /// distance, then Gaussian-smoothed with width `4 * average distance` when
    /// `smooth_data` is set (source `prepareFit_`).
    fn prepare_fit(current_section: &[Point2D], smooth_data: bool) -> Result<Vec<Peak1D>> {
        let mut filter_spec = MSSpectrum::default();
        for pa in current_section {
            filter_spec.peaks.push(Peak1D::new(pa.rt, pa.mz as f32));
        }
        let mut distances = Vec::with_capacity(filter_spec.peaks.len());
        for j in 1..filter_spec.peaks.len() {
            distances.push(filter_spec.peaks[j].mz - filter_spec.peaks[j - 1].mz);
        }
        let dist_average = distances.iter().fold(0.0, |acc, &d| acc + d) / distances.len() as f64;
        for _ in 0..3 {
            let last = filter_spec.peaks[filter_spec.peaks.len() - 1].mz;
            filter_spec
                .peaks
                .push(Peak1D::new(last + dist_average, 0.0));
        }
        for _ in 0..3 {
            let first = filter_spec.peaks[0].mz;
            filter_spec
                .peaks
                .insert(0, Peak1D::new(first - dist_average, 0.0));
        }
        // "To get an estimate of the peak quality, we probably should not
        // smooth and/or transform the data."
        if smooth_data {
            let filter = GaussFilter::new(GaussianWidth::Absolute(4.0 * dist_average))?;
            filter.filter_spectrum(&mut filter_spec)?;
        }
        Ok(filter_spec.peaks)
    }
}
