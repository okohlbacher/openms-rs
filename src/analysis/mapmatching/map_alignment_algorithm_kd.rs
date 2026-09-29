// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! RT alignment from conflict-free connected components of a kd-tree
//! compatibility graph, as KD feature grouping uses it to warp retention
//! times before linking.
//!
//! Native equivalent of `ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmKD.h`. See
//! `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! The algorithm computes connected components (CCs) of a compatibility graph
//! on the features of a [`KDTreeFeatureMaps`]: two features of different maps
//! are connected when each lies in the other's RT and m/z tolerance window
//! (and, optionally, their intensities differ by at most a log10 fold
//! change). CCs of sufficient size without too many features from the same
//! map are taken as sets of corresponding features; every feature of such a
//! CC contributes the point (its RT, the CC's mean RT) to its map's fit data,
//! and a LOWESS transformation per map is fitted to those points.
//!
//! # Differences from the source
//!
//! - The source documents that the CCs must also have one charge state, and
//!   collects the non-zero charges of each CC, but the `continue` that should
//!   reject a CC with two charges only continues the inner loop, so the check
//!   has no effect. The port reproduces that: charges are not checked.
//! - The LOWESS fit is `TransformationModelLowess`, which this port does not
//!   provide yet; [`LowessParameters`] reads its parameter section and fits
//!   with [`LowessModel`] (the same `FastLowessSmoothing` and interpolated
//!   model). Its automatic span selection by cross-validation (`auto_span`),
//!   `akima` interpolation, and more than 64 robustifying iterations are
//!   refused with [`Error::Unsupported`](crate::error::Error::Unsupported).
//! - The source's parameter-less constructor leaves the tolerances
//!   uninitialised; here an empty parameter tree leaves the algorithm
//!   unconfigured and every operation that needs them fails.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use super::kd_tree_feature_maps::KDTreeFeatureMaps;
use super::params::{add_float, add_int, add_string, bad, float, int, string};
use crate::analysis::transformations::{
    DataPoint, Extrapolation, Interpolation, LowessModel, LowessOptions,
};
use crate::error::{Error, Result};
use crate::param::Param;

/// Largest number of fit points one map collects, and the point ceiling of
/// the LOWESS fit.
pub const MAX_FIT_POINTS: usize = 1 << 24;

/// Work ceiling (local-regression sample operations) of one LOWESS fit.
pub const MAX_LOWESS_WORK: usize = 1 << 40;

/// The parameters of the source's `TransformationModelLowess` (its
/// `getDefaultParameters`), as KD feature grouping stores them under
/// `LOWESS:`.
#[derive(Clone, Debug, PartialEq)]
pub struct LowessParameters {
    /// Fraction of data points per local regression (`span`, default 2/3,
    /// 0 to 1).
    pub span: f64,
    /// Select the span by cross-validation (`auto_span`, default false; not
    /// supported by this port).
    pub auto_span: bool,
    /// Lower bound for the auto-selected span (`auto_span_min`, 0.15).
    pub auto_span_min: f64,
    /// Upper bound for the auto-selected span (`auto_span_max`, 0.80).
    pub auto_span_max: f64,
    /// Minimum neighbours in auto mode (`auto_min_neighbors`, 5).
    pub auto_min_neighbors: i64,
    /// Folds for cross-validation (`auto_k_folds`, 5).
    pub auto_k_folds: i64,
    /// Cross-validation metric (`auto_metric`, `mae`).
    pub auto_metric: String,
    /// Explicit span grid (`auto_span_grid`, empty).
    pub auto_span_grid: String,
    /// Robustifying iterations (`num_iterations`, default 3, minimum 0).
    pub num_iterations: i64,
    /// Skip distance; negative means 1% of the x range (`delta`, -1).
    pub delta: f64,
    /// Interpolation between the smoothed points (`interpolation_type`:
    /// `linear`, `cspline` (default) or `akima`).
    pub interpolation_type: String,
    /// Extrapolation outside the data (`extrapolation_type`:
    /// `two-point-linear`, `four-point-linear` (default) or
    /// `global-linear`).
    pub extrapolation_type: String,
}

impl Default for LowessParameters {
    fn default() -> Self {
        Self {
            span: 2.0 / 3.0,
            auto_span: false,
            auto_span_min: 0.15,
            auto_span_max: 0.80,
            auto_min_neighbors: 5,
            auto_k_folds: 5,
            auto_metric: "mae".into(),
            auto_span_grid: String::new(),
            num_iterations: 3,
            delta: -1.0,
            interpolation_type: "cspline".into(),
            extrapolation_type: "four-point-linear".into(),
        }
    }
}

impl LowessParameters {
    /// The source `TransformationModelLowess::getDefaultParameters` tree,
    /// with `advanced` on every entry when `advanced` is set (as KD feature
    /// grouping tags them).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self, advanced: bool) -> Result<Param> {
        let mut p = Param::new();
        add_float(
            &mut p,
            "span",
            self.span,
            "Fraction of datapoints (f) to use for each local regression (determines the amount of smoothing). Choosing this parameter in the range .2 to .8 usually results in a good fit.",
            advanced,
            Some(0.0),
            Some(1.0),
        )?;
        add_string(
            &mut p,
            "auto_span",
            if self.auto_span { "true" } else { "false" },
            "If true, or if 'span' is 0, automatically select LOWESS span by cross-validation.",
            advanced,
            &["true", "false"],
        )?;
        add_float(
            &mut p,
            "auto_span_min",
            self.auto_span_min,
            "Lower bound for auto-selected span.",
            advanced,
            Some(0.001),
            None,
        )?;
        add_float(
            &mut p,
            "auto_span_max",
            self.auto_span_max,
            "Upper bound for auto-selected span.",
            advanced,
            None,
            Some(0.99),
        )?;
        add_int(
            &mut p,
            "auto_min_neighbors",
            self.auto_min_neighbors,
            "Minimum number of neighbors (span*n) enforced in auto mode.",
            advanced,
            Some(3),
            None,
        )?;
        add_int(
            &mut p,
            "auto_k_folds",
            self.auto_k_folds,
            "K-folds for CV when n>50 (else LOO is used).",
            advanced,
            Some(2),
            None,
        )?;
        add_string(
            &mut p,
            "auto_metric",
            &self.auto_metric,
            "Metric for CV selection: one of {'p90','p95','p99','rmse','mae'}.",
            advanced,
            &["p90", "p95", "p99", "rmse", "mae"],
        )?;
        add_string(
            &mut p,
            "auto_span_grid",
            &self.auto_span_grid,
            "Optional explicit grid of span candidates in (0,1]. Comma-separated list, e.g. '0.2,0.3,0.5'.  If empty, a default grid is used.",
            advanced,
            &[],
        )?;
        add_int(
            &mut p,
            "num_iterations",
            self.num_iterations,
            "Number of robustifying iterations for lowess fitting.",
            advanced,
            Some(0),
            None,
        )?;
        add_float(
            &mut p,
            "delta",
            self.delta,
            "Nonnegative parameter which may be used to save computations (recommended value is 0.01 of the range of the input, e.g. for data ranging from 1000 seconds to 2000 seconds, it could be set to 10). Setting a negative value will automatically do this.",
            advanced,
            None,
            None,
        )?;
        add_string(
            &mut p,
            "interpolation_type",
            &self.interpolation_type,
            "Method to use for interpolation between datapoints computed by lowess. 'linear': Linear interpolation. 'cspline': Use the cubic spline for interpolation. 'akima': Use an akima spline for interpolation",
            advanced,
            &["linear", "cspline", "akima"],
        )?;
        add_string(
            &mut p,
            "extrapolation_type",
            &self.extrapolation_type,
            "Method to use for extrapolation outside the data range. 'two-point-linear': Uses a line through the first and last point to extrapolate. 'four-point-linear': Uses a line through the first and second point to extrapolate in front and and a line through the last and second-to-last point in the end. 'global-linear': Uses a linear regression to fit a line through all data points and use it for interpolation.",
            advanced,
            &["two-point-linear", "four-point-linear", "global-linear"],
        )?;
        Ok(p)
    }

    /// Read a validated tree whose every default is present.
    pub(crate) fn from_complete_param(p: &Param) -> Result<Self> {
        Ok(Self {
            span: float(p, "span")?,
            auto_span: string(p, "auto_span")? == "true",
            auto_span_min: float(p, "auto_span_min")?,
            auto_span_max: float(p, "auto_span_max")?,
            auto_min_neighbors: int(p, "auto_min_neighbors")?,
            auto_k_folds: int(p, "auto_k_folds")?,
            auto_metric: string(p, "auto_metric")?,
            auto_span_grid: string(p, "auto_span_grid")?,
            num_iterations: int(p, "num_iterations")?,
            delta: float(p, "delta")?,
            interpolation_type: string(p, "interpolation_type")?,
            extrapolation_type: string(p, "extrapolation_type")?,
        })
    }

    /// Fit the LOWESS transformation to `data`, as the source's
    /// `TransformationModelLowess(data, params)` does without automatic span
    /// selection: sort by x, smooth with `FastLowessSmoothing` (a negative
    /// `delta` becomes 1% of the x range), and interpolate the smoothed
    /// points.
    ///
    /// # Errors
    ///
    /// [`Error::Unsupported`] for `auto_span`, `akima` interpolation or more
    /// than 64 iterations; [`Error::InvalidValue`] for fewer than two points
    /// (source `Exception::IllegalArgument`, "'lowess' model requires more
    /// data"), an invalid option, and the fit's own errors.
    pub fn fit(&self, data: &[DataPoint]) -> Result<LowessModel> {
        if self.auto_span {
            return Err(Error::Unsupported(
                "LOWESS span selection by cross-validation (auto_span)".into(),
            ));
        }
        if data.len() < 2 {
            return Err(bad("'lowess' model requires more data"));
        }
        let interpolation = match self.interpolation_type.as_str() {
            "linear" => Interpolation::Linear,
            "cspline" => Interpolation::CubicSpline,
            other => {
                return Err(Error::Unsupported(format!(
                    "LOWESS interpolation type '{other}'"
                )));
            }
        };
        let extrapolation = match self.extrapolation_type.as_str() {
            "two-point-linear" => Extrapolation::TwoPointLinear,
            "four-point-linear" => Extrapolation::FourPointLinear,
            "global-linear" => Extrapolation::GlobalLinear,
            other => {
                return Err(bad(format!("invalid LOWESS extrapolation type '{other}'")));
            }
        };
        let iterations = usize::try_from(self.num_iterations)
            .map_err(|_| bad("LOWESS num_iterations must not be negative"))?;
        if iterations > 64 {
            return Err(Error::Unsupported(
                "more than 64 LOWESS robustifying iterations".into(),
            ));
        }
        let options = LowessOptions {
            span: self.span,
            iterations,
            delta: if self.delta < 0.0 {
                None
            } else {
                Some(self.delta)
            },
            interpolation,
            extrapolation,
            max_points: MAX_FIT_POINTS,
            max_work: MAX_LOWESS_WORK,
        };
        LowessModel::fit(data, options)
    }
}

/// The settings the source reads in `updateMembers_` (and `filterCCs_`)
/// from the KD grouping parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct KdWarpSettings {
    /// RT tolerance of the compatibility graph (`warp:rt_tol`).
    pub rt_tol: f64,
    /// m/z tolerance of the compatibility graph (`warp:mz_tol`).
    pub mz_tol: f64,
    /// Whether `mz_tol` is in ppm (`mz_unit` = `ppm`).
    pub mz_ppm: bool,
    /// Maximum absolute log10 fold change between connected features; a
    /// negative value disables the check (`warp:max_pairwise_log_fc`).
    pub max_pairwise_log_fc: f64,
    /// Minimum CC size relative to the number of maps
    /// (`warp:min_rel_cc_size`).
    pub min_rel_cc_size: f64,
    /// Features from an already represented map a CC may hold; -1 allows any
    /// number (`warp:max_nr_conflicts`).
    pub max_nr_conflicts: i64,
    /// The LOWESS parameters (`LOWESS:`).
    pub lowess: LowessParameters,
}

impl KdWarpSettings {
    /// Read the settings from the KD grouping parameters: `warp:rt_tol`,
    /// `warp:mz_tol`, `mz_unit`, `warp:max_pairwise_log_fc`,
    /// `warp:min_rel_cc_size`, `warp:max_nr_conflicts` and the `LOWESS:`
    /// section (missing LOWESS entries take their defaults, as the source's
    /// `setDefaults` fills them).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a missing or mistyped entry.
    pub fn from_param(p: &Param) -> Result<Self> {
        // Source: `TransformationModelLowess` receives `param_.copy("LOWESS:",
        // true)` and calls `setDefaults` on it.
        let mut lowess = p.copy("LOWESS:", true)?;
        lowess.set_defaults(&LowessParameters::default().to_param(false)?, "", false)?;
        Ok(Self {
            rt_tol: float(p, "warp:rt_tol")?,
            mz_tol: float(p, "warp:mz_tol")?,
            mz_ppm: string(p, "mz_unit")? == "ppm",
            max_pairwise_log_fc: float(p, "warp:max_pairwise_log_fc")?,
            min_rel_cc_size: float(p, "warp:min_rel_cc_size")?,
            max_nr_conflicts: int(p, "warp:max_nr_conflicts")?,
            lowess: LowessParameters::from_complete_param(&lowess)?,
        })
    }
}

/// KD-tree based RT alignment (source `MapAlignmentAlgorithmKD`).
#[derive(Clone, Debug)]
pub struct MapAlignmentAlgorithmKD {
    fit_data: Vec<Vec<DataPoint>>,
    transformations: Vec<Option<LowessModel>>,
    settings: Option<KdWarpSettings>,
}

impl MapAlignmentAlgorithmKD {
    /// An aligner for `num_maps` maps with the KD grouping parameters
    /// `param` (source constructor). An empty tree leaves it unconfigured.
    ///
    /// # Errors
    ///
    /// The [`KdWarpSettings::from_param`] errors for a non-empty tree.
    pub fn new(num_maps: usize, param: &Param) -> Result<Self> {
        let settings = if param.is_empty() {
            None
        } else {
            Some(KdWarpSettings::from_param(param)?)
        };
        Ok(Self::with_settings(num_maps, settings))
    }

    /// An aligner with already read settings.
    pub fn with_settings(num_maps: usize, settings: Option<KdWarpSettings>) -> Self {
        Self {
            fit_data: vec![Vec::new(); num_maps],
            transformations: vec![None; num_maps],
            settings,
        }
    }

    fn settings(&self) -> Result<&KdWarpSettings> {
        self.settings
            .as_ref()
            .ok_or_else(|| bad("MapAlignmentAlgorithmKD has no parameters"))
    }

    /// The fit data collected so far, one list of (RT, CC mean RT) per map.
    pub fn fit_data(&self) -> &[Vec<DataPoint>] {
        &self.fit_data
    }

    /// Collect the fit points of the conflict-free CCs of `kd_data` (source
    /// `addRTFitData`).
    ///
    /// # Errors
    ///
    /// When unconfigured, when a map index exceeds the number of maps, or
    /// when a map exceeds [`MAX_FIT_POINTS`].
    pub fn add_rt_fit_data(&mut self, kd_data: &KDTreeFeatureMaps<'_>) -> Result<()> {
        let ccs = self.get_ccs(kd_data)?;
        let filtered = self.filter_ccs(kd_data, &ccs)?;
        let mut avg_rts = BTreeMap::new();
        for (cc_index, cc) in &filtered {
            let mut avg_rt = 0.0;
            for &i in cc {
                avg_rt += kd_data.rt(i)?;
            }
            avg_rt /= cc.len() as f64;
            avg_rts.insert(*cc_index, avg_rt);
        }
        for (cc_index, cc) in &filtered {
            for &i in cc {
                let map = kd_data.map_index(i)?;
                let data = self
                    .fit_data
                    .get_mut(map)
                    .ok_or_else(|| bad(format!("map index {map} exceeds the aligner's maps")))?;
                if data.len() == MAX_FIT_POINTS {
                    return Err(bad(format!(
                        "RT fit data exceeds {MAX_FIT_POINTS} points for map {map}"
                    )));
                }
                data.push(DataPoint::new(kd_data.rt(i)?, avg_rts[cc_index]));
            }
        }
        Ok(())
    }

    /// Fit one LOWESS transformation per map (source `fitLOWESS`). A map
    /// with fewer than 50 fit points gets the identity, fitted through (0,
    /// 0), (1, 1) and (1e6, 1e6), and a warning (returned, one per map, in
    /// the source's wording).
    ///
    /// # Errors
    ///
    /// When unconfigured, and the [`LowessParameters::fit`] errors (the
    /// source's caller logs them and abandons the grouping).
    pub fn fit_lowess(&mut self) -> Result<Vec<String>> {
        let lowess = self.settings()?.lowess.clone();
        let mut warnings = Vec::new();
        for i in 0..self.fit_data.len() {
            let n = self.fit_data[i].len();
            let model = if n < 50 {
                warnings.push(format!(
                    "Warning: Only {n} data points for LOWESS fit of map {i}. Consider adjusting RT or m/z tolerance or max_pairwise_log_fc, decreasing min_rel_cc_size, or increasing max_nr_conflicts."
                ));
                let identity = [
                    DataPoint::new(0.0, 0.0),
                    DataPoint::new(1.0, 1.0),
                    DataPoint::new(1e6, 1e6),
                ];
                lowess.fit(&identity)?
            } else {
                lowess.fit(&self.fit_data[i])?
            };
            self.transformations[i] = Some(model);
        }
        Ok(warnings)
    }

    /// Apply the fitted transformations to `kd_data`'s retention times and
    /// rebalance its tree (source `transform`).
    ///
    /// # Errors
    ///
    /// When a map has no fitted transformation (the source dereferences a
    /// null pointer), and the transformation's evaluation errors.
    pub fn transform(&self, kd_data: &mut KDTreeFeatureMaps<'_>) -> Result<()> {
        let transformations = &self.transformations;
        kd_data.apply_transformations(|map, rt| {
            transformations
                .get(map)
                .and_then(Option::as_ref)
                .ok_or_else(|| bad(format!("no LOWESS transformation for map {map}")))?
                .apply(rt)
        })?;
        kd_data.optimize_tree()
    }

    /// The CC index of every feature, by repeated breadth-first search over
    /// the compatibility graph without storing it (source `computeCCs_`),
    /// and the number of CCs.
    fn compute_ccs(&self, kd_data: &KDTreeFeatureMaps<'_>) -> Result<(Vec<usize>, usize)> {
        let s = self.settings()?;
        let n = kd_data.size();
        let mut result = vec![usize::MAX; n];
        let mut visited = vec![false; n];
        let mut queue = VecDeque::new();
        let mut search_pos = 0;
        let mut cc_index = 0;
        loop {
            let Some(start) = (search_pos..n).find(|&i| !visited[i]) else {
                break;
            };
            queue.push_back(start);
            visited[start] = true;
            search_pos = start + 1;
            while let Some(i) = queue.pop_front() {
                result[i] = cc_index;
                let mut compatible = Vec::new();
                kd_data.get_neighborhood(
                    i,
                    &mut compatible,
                    s.rt_tol,
                    s.mz_tol,
                    s.mz_ppm,
                    false,
                    s.max_pairwise_log_fc,
                )?;
                for j in compatible {
                    if !visited[j] {
                        queue.push_back(j);
                        visited[j] = true;
                    }
                }
            }
            cc_index += 1;
        }
        Ok((result, cc_index))
    }

    /// The CCs, each the ascending list of its features (source `getCCs_`).
    fn get_ccs(&self, kd_data: &KDTreeFeatureMaps<'_>) -> Result<BTreeMap<usize, Vec<usize>>> {
        let (cc_index, _) = self.compute_ccs(kd_data)?;
        let mut result: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (i, cc) in cc_index.into_iter().enumerate() {
            result.entry(cc).or_default().push(i);
        }
        Ok(result)
    }

    /// The CCs of at least `max(2, min_rel_cc_size * num_maps)` features
    /// with at most `max_nr_conflicts` features from an already represented
    /// map (source `filterCCs_`; its charge check is a no-op, see the module
    /// documentation).
    fn filter_ccs(
        &self,
        kd_data: &KDTreeFeatureMaps<'_>,
        ccs: &BTreeMap<usize, Vec<usize>>,
    ) -> Result<BTreeMap<usize, Vec<usize>>> {
        let s = self.settings()?;
        let num_maps = self.fit_data.len();
        // `Size min_size = max(2.0, rel * num_maps)`: the double truncates.
        let min_size = 2.0_f64.max(s.min_rel_cc_size * num_maps as f64) as usize;
        let mut filtered = BTreeMap::new();
        for (index, cc) in ccs {
            if cc.len() < min_size {
                continue;
            }
            let mut passes = true;
            if s.max_nr_conflicts != -1 {
                let mut map_indices = BTreeSet::new();
                let mut nr_conflicts = 0i64;
                for &i in cc {
                    let map = kd_data.map_index(i)?;
                    if map_indices.contains(&map) {
                        nr_conflicts += 1;
                        if nr_conflicts > s.max_nr_conflicts {
                            passes = false;
                            break;
                        }
                    } else {
                        map_indices.insert(map);
                    }
                }
            }
            if passes {
                filtered.insert(*index, cc.clone());
            }
        }
        Ok(filtered)
    }
}
