// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Clustering the results of the multiplex filtering.
//!
//! Native equivalent of `FEATUREFINDER/MultiplexClustering.h`. See
//! `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md`.

use std::collections::BTreeMap;

use super::filtered_ms_experiment::MultiplexFilteredMSExperiment;
use crate::analysis::featurefinder::peak_width_estimator::PeakWidthEstimator;
use crate::concept::progress_logger::ProgressReporter;
use crate::kernel::{MSExperiment, NumericRange};
use crate::math::ml::clustering::grid_based_cluster::{GridBasedCluster, Point};
use crate::math::ml::clustering::grid_based_clustering::GridBasedClustering;
use crate::math::source_sort::source_sort_by;
use crate::processing::peak_picking::PeakBoundary;
use crate::{Error, Result};

/// Ceiling on the grid lines in either direction, checked while the grid is
/// built. The source has none, and loops forever on a zero step.
pub const MAX_GRID_LINES: usize = 50_000_000;

/// Scaled Euclidean distance for the clustering (source
/// `MultiplexClustering::MultiplexDistance`).
///
/// Points are (m/z, RT); the RT difference is scaled by `rt_scaling` before the
/// Euclidean distance is taken.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultiplexDistance {
    /// Scaling of RT coordinates (source `rt_scaling_`).
    rt_scaling: f64,
}

impl Default for MultiplexDistance {
    /// Scaling 1 (source default constructor).
    fn default() -> Self {
        Self { rt_scaling: 1.0 }
    }
}

impl MultiplexDistance {
    /// A distance with the RT scaling `rt_scaling` (source constructor).
    pub const fn new(rt_scaling: f64) -> Self {
        Self { rt_scaling }
    }

    /// The RT scaling.
    pub const fn rt_scaling(&self) -> f64 {
        self.rt_scaling
    }

    /// The distance between `p1` and `p2` in the (m/z, RT) plane (source
    /// `operator()`): `sqrt(dx*dx + s*s*dy*dy)`, evaluated left to right.
    pub fn distance(&self, p1: &Point, p2: &Point) -> f64 {
        let dx = p1[0] - p2[0];
        let dy = p1[1] - p2[1];
        (dx * dx + self.rt_scaling * self.rt_scaling * dy * dy).sqrt()
    }
}

/// Clusters the results of the multiplex filtering (source
/// `MultiplexClustering`).
///
/// The filtering identified regions of the picked (and profile) data that
/// correspond to peptide features. This groups the data points that belong to
/// the same peptide feature, with the general-purpose
/// [`GridBasedClustering`] on a grid whose m/z spacing follows the peak width
/// (or the m/z tolerance) and whose RT spacing is the typical elution time.
///
/// The source derives from `ProgressLogger`; [`Self::cluster_with_progress`]
/// takes a reporter instead.
#[derive(Clone, Debug, PartialEq)]
pub struct MultiplexClustering {
    /// Grid lines in m/z (source `grid_spacing_mz_`).
    grid_spacing_mz: Vec<f64>,
    /// Grid lines in RT (source `grid_spacing_rt_`).
    grid_spacing_rt: Vec<f64>,
    /// Scaling of RT in the distance (source `rt_scaling_`).
    rt_scaling: f64,
    /// Typical retention time (source `rt_typical_`).
    rt_typical: f64,
}

/// Grid lines from `min` while below `max`, each the previous advanced by
/// `step`, then `max` itself.
fn grid(min: f64, max: f64, mut step: impl FnMut(f64) -> Result<f64>) -> Result<Vec<f64>> {
    let mut lines = Vec::new();
    let mut value = min;
    while value < max {
        if lines.len() >= MAX_GRID_LINES {
            return Err(Error::InvalidRange(format!(
                "MultiplexClustering: more than {MAX_GRID_LINES} grid lines"
            )));
        }
        lines.push(value);
        let next = step(value)?;
        if !(next > value) {
            return Err(Error::InvalidValue(
                "MultiplexClustering: the grid spacing is not positive (the source loops forever)"
                    .into(),
            ));
        }
        value = next;
    }
    lines.push(max);
    Ok(lines)
}

/// The median m/z of all peaks of `exp` as the source reads it: the element at
/// index `n / 2` of the ascending `std::sort` of every peak's m/z.
fn median_mz(exp: &MSExperiment) -> Result<f64> {
    let mut mz: Vec<f64> = exp
        .spectra
        .iter()
        .flat_map(|s| s.peaks.iter().map(|p| p.mz))
        .collect();
    if mz.is_empty() {
        return Err(Error::InvalidValue(
            "MultiplexClustering: no peak to determine the RT scaling from (the source reads an empty vector)".into(),
        ));
    }
    source_sort_by(&mut mz, |a, b| a < b)?;
    Ok(mz[mz.len() / 2])
}

fn required(range: Option<NumericRange>, what: &str) -> Result<NumericRange> {
    range.ok_or_else(|| {
        Error::InvalidValue(format!(
            "MinMZ,MaxMZ,MinRT,MaxRT values outside of sensible value ranges. Are they uninitialized? (no {what} range)"
        ))
    })
}

impl MultiplexClustering {
    /// Prepare clustering for profile input (source constructor taking the
    /// profile and centroided data, the boundaries and `rt_typical`).
    ///
    /// The grid spans the m/z and RT ranges of `exp_profile` (all MS levels)
    /// widened by 0.01 on each side. Its m/z lines advance by 0.4 times the
    /// peak width a [`PeakWidthEstimator`] of `exp_picked` and `boundaries`
    /// estimates there ("jitter of the peak centres is less than 0.4 peak
    /// widths", so neighbouring peaks at one RT never share a cluster), its RT
    /// lines by `rt_typical`, the elution time of a characteristic peptide.
    /// The RT scaling is the estimated peak width at the median m/z of all
    /// picked peaks divided by `rt_typical`.
    ///
    /// The source reads the ranges that the caller last stored with
    /// `updateRanges`; the port computes them from the data.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] "Centroided data and the corresponding list of
    ///   peak boundaries do not contain same number of spectra." (source
    ///   `Exception::IllegalArgument`).
    /// - [`Error::InvalidValue`] when the profile data holds no peak or
    ///   spectrum (the source uses its uninitialised ranges), when there is no
    ///   picked peak, and when a grid step is not positive; and
    ///   [`Error::InvalidRange`] beyond [`MAX_GRID_LINES`].
    /// - The errors of [`PeakWidthEstimator`].
    pub fn from_profile(
        exp_profile: &MSExperiment,
        exp_picked: &MSExperiment,
        boundaries: &[Vec<PeakBoundary>],
        rt_typical: f64,
    ) -> Result<Self> {
        if exp_picked.len() != boundaries.len() {
            return Err(Error::InvalidValue(
                "Centroided data and the corresponding list of peak boundaries do not contain same number of spectra."
                    .into(),
            ));
        }
        let ranges = exp_profile.ranges(0)?;
        let mz_range = required(ranges.mz, "m/z")?;
        let rt_range = required(ranges.rt, "RT")?;
        let margin = 1e-2;
        let (mz_min, mz_max) = (mz_range.min - margin, mz_range.max + margin);
        let (rt_min, rt_max) = (rt_range.min - margin, rt_range.max + margin);

        let estimator = PeakWidthEstimator::new(exp_picked, boundaries)?;
        let scaling = 0.4;
        let grid_spacing_mz = grid(mz_min, mz_max, |mz| {
            Ok(mz + scaling * estimator.peak_width(mz)?)
        })?;
        let grid_spacing_rt = grid(rt_min, rt_max, |rt| Ok(rt + rt_typical))?;
        let rt_scaling = estimator.peak_width(median_mz(exp_picked)?)? / rt_typical;
        Ok(Self {
            grid_spacing_mz,
            grid_spacing_rt,
            rt_scaling,
            rt_typical,
        })
    }

    /// Prepare clustering for centroided input (source constructor taking the
    /// data, the m/z tolerance, its unit and `rt_typical`).
    ///
    /// The grid spans the m/z and RT ranges of the MS1 spectra of `exp`,
    /// widened by 0.01 on each side. Its m/z lines advance by the tolerance
    /// (the jitter of peak centres is assumed below one tolerance): multiplied
    /// by `1 + tol / 1e6` in ppm mode (`mz_tolerance_unit_in_ppm`), increased
    /// by `tol` Th otherwise. Its RT lines advance by `rt_typical`. The RT
    /// scaling is the tolerance in Th at the median m/z of **all** peaks (every
    /// MS level, as the source collects them) divided by `rt_typical`.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] "MinMZ,MaxMZ,MinRT,MaxRT values outside of
    ///   sensible value ranges. Are they uninitialized? ..." when there are no
    ///   MS1 data or the m/z range is not inside `[0, 1e12]` or the RT range
    ///   not inside `[-1e12, 1e12]` (source `Exception::IllegalArgument`).
    /// - [`Error::InvalidValue`] when there is no peak, and when a grid step is
    ///   not positive; [`Error::InvalidRange`] beyond [`MAX_GRID_LINES`].
    pub fn from_centroided(
        exp: &MSExperiment,
        mz_tolerance: f64,
        mz_tolerance_unit_in_ppm: bool,
        rt_typical: f64,
    ) -> Result<Self> {
        let ranges = exp.ranges(1)?;
        let mz_range = required(ranges.mz, "m/z")?;
        let rt_range = required(ranges.rt, "RT")?;
        if !(mz_range.min >= 0.0
            && mz_range.max <= 1.0e12
            && rt_range.min >= -1.0e12
            && rt_range.max <= 1.0e12)
        {
            return Err(Error::InvalidValue(format!(
                "MinMZ,MaxMZ,MinRT,MaxRT values outside of sensible value ranges. Are they uninitialized? ({}/{}/{}/{}",
                mz_range.min, mz_range.max, rt_range.min, rt_range.max
            )));
        }
        let margin = 1e-2;
        let (mz_min, mz_max) = (mz_range.min - margin, mz_range.max + margin);
        let (rt_min, rt_max) = (rt_range.min - margin, rt_range.max + margin);
        let scaling = 1.0;
        let grid_spacing_mz = if mz_tolerance_unit_in_ppm {
            grid(mz_min, mz_max, |mz| {
                Ok(mz * (1.0 + scaling * mz_tolerance / 1_000_000.0))
            })?
        } else {
            grid(mz_min, mz_max, |mz| Ok(mz + scaling * mz_tolerance))?
        };
        let grid_spacing_rt = grid(rt_min, rt_max, |rt| Ok(rt + rt_typical))?;
        let median = median_mz(exp)?;
        let rt_scaling = if mz_tolerance_unit_in_ppm {
            (median * mz_tolerance / 1_000_000.0) / rt_typical
        } else {
            mz_tolerance / rt_typical
        };
        Ok(Self {
            grid_spacing_mz,
            grid_spacing_rt,
            rt_scaling,
            rt_typical,
        })
    }

    /// The m/z grid lines (source `grid_spacing_mz_`).
    pub fn grid_spacing_mz(&self) -> &[f64] {
        &self.grid_spacing_mz
    }

    /// The RT grid lines (source `grid_spacing_rt_`).
    pub fn grid_spacing_rt(&self) -> &[f64] {
        &self.grid_spacing_rt
    }

    /// The RT scaling of the distance (source `rt_scaling_`).
    pub fn rt_scaling(&self) -> f64 {
        self.rt_scaling
    }

    /// The typical retention time (source `rt_typical_`).
    pub fn rt_typical(&self) -> f64 {
        self.rt_typical
    }

    /// Cluster the filter results (source `cluster`).
    ///
    /// # Errors
    ///
    /// As [`Self::cluster_with_progress`].
    pub fn cluster(
        &self,
        filter_results: &[MultiplexFilteredMSExperiment],
    ) -> Result<Vec<BTreeMap<i32, GridBasedCluster>>> {
        self.cluster_with_progress(filter_results, &mut ProgressReporter::silent())
    }

    /// [`Self::cluster`], reporting progress with the label "clustering
    /// filtered LC-MS data", one step per filter result.
    ///
    /// Each filter result's peaks, at (m/z, RT), are clustered on the grid with
    /// the [`MultiplexDistance`] of this RT scaling; each map holds one
    /// pattern's clusters by cluster index. The source's commented-out
    /// `extendClustersY` step is not run.
    ///
    /// # Errors
    ///
    /// The errors of [`GridBasedClustering`] (a point outside the grid, among
    /// others) and of `progress`.
    pub fn cluster_with_progress(
        &self,
        filter_results: &[MultiplexFilteredMSExperiment],
        progress: &mut ProgressReporter<'_>,
    ) -> Result<Vec<BTreeMap<i32, GridBasedCluster>>> {
        progress.start_count(filter_results.len(), "clustering filtered LC-MS data")?;
        let distance = MultiplexDistance::new(self.rt_scaling);
        let mut cluster_results = Vec::with_capacity(filter_results.len());
        for (i, result) in filter_results.iter().enumerate() {
            progress.set_count(i + 1)?;
            let mut clustering = GridBasedClustering::without_properties(
                |p1: &Point, p2: &Point| distance.distance(p1, p2),
                &result.mz(),
                &result.rt(),
                self.grid_spacing_mz.clone(),
                self.grid_spacing_rt.clone(),
            )?;
            clustering.cluster()?;
            cluster_results.push(clustering.results().clone());
        }
        progress.end()?;
        Ok(cluster_results)
    }
}
