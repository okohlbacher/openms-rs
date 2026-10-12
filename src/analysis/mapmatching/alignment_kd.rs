// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reference-free map alignment on a kd-tree of pooled features.
//!
//! Ports `ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmKD.h` and its `.cpp`.
//! The source works on a `KDTreeFeatureMaps` (`ANALYSIS/QUANTITATION`); the
//! part of it the algorithm uses is the trait
//! [`KdFeatureMaps`](crate::analysis::mapmatching::alignment_kd::KdFeatureMaps),
//! which [`KDTreeFeatureMaps`](crate::analysis::mapmatching::kd_tree_feature_maps::KDTreeFeatureMaps)
//! implements and which the brute-force pool of the tests implements too. This is the
//! one implementation of the header: KD feature grouping
//! ([`FeatureGroupingAlgorithmKD`](crate::analysis::mapmatching::feature_grouping_algorithm_kd::FeatureGroupingAlgorithmKD))
//! drives it on the real kd-tree. See `docs/MAP_ALIGNMENT_SUPPORT.md` and
//! `docs/FEATURE_GROUPING_SUPPORT.md`.

use super::kd_tree_feature_maps::KDTreeFeatureMaps;
use super::superimposer::float;
use crate::analysis::transformations::{DataPoint, LowessModel, ModelConfig};
use crate::param::{Param, ParamValue};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// Largest number of fit points one map collects. Checked on every insertion;
/// the source has no ceiling.
pub const MAX_FIT_POINTS: usize = 1 << 24;

/// The pooled features of several maps with the range queries of a kd-tree
/// on (RT, m/z): the members of the source's `KDTreeFeatureMaps` that
/// [`MapAlignmentAlgorithmKD`] calls.
///
/// Features are addressed by a global index `0..len()`; the accessors are
/// called with indices below `len()` only (an implementation may panic or
/// answer NaN beyond it). The retention time is the cached, possibly
/// transformed one; m/z, intensity and charge are the feature's own.
pub trait KdFeatureMaps {
    /// Number of pooled features (source `size`).
    fn len(&self) -> usize;

    /// Whether no feature is pooled.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The cached, possibly transformed, retention time of feature `index`
    /// (source `rt`).
    fn rt(&self, index: usize) -> f64;

    /// The charge of feature `index` (source `charge`).
    fn charge(&self, index: usize) -> i32;

    /// The input map feature `index` was added from (source `mapIndex`).
    fn map_index(&self, index: usize) -> usize;

    /// The features compatible with feature `index` (source
    /// `getNeighborhood`): those inside the tolerance window
    /// `Math::getTolWindow` gives around its RT (absolute `rt_tol`) and m/z
    /// (`mz_tol`, in ppm with `mz_ppm`), bounds included; without
    /// `include_features_from_same_map`, not from its own map; and with a
    /// non-negative `max_pairwise_log_fc`, only those whose
    /// `|log10(intensity / intensity of index)|` is at most that.
    ///
    /// [`MapAlignmentAlgorithmKD`] does not depend on the order of the
    /// result.
    ///
    /// # Errors
    ///
    /// As the implementation documents.
    fn neighborhood(
        &self,
        index: usize,
        rt_tol: f64,
        mz_tol: f64,
        mz_ppm: bool,
        include_features_from_same_map: bool,
        max_pairwise_log_fc: f64,
    ) -> Result<Vec<usize>>;

    /// Replace the cached retention time of every feature by
    /// `transformations[map_index].apply(original RT)` (source
    /// `applyTransformations`), then rebuild the tree on the new retention
    /// times (source `optimizeTree`).
    ///
    /// # Errors
    ///
    /// As the implementation documents.
    fn apply_transformations(&mut self, transformations: &[LowessModel]) -> Result<()>;
}

/// Tolerances read from the parameters (source `updateMembers_`).
#[derive(Clone, Copy, Debug, PartialEq)]
struct Tolerances {
    rt_tol_secs: f64,
    mz_tol: f64,
    mz_ppm: bool,
    max_pairwise_log_fc: f64,
}

/// An efficient reference-free feature map alignment algorithm for unlabeled
/// data (source `MapAlignmentAlgorithmKD`).
///
/// The algorithm computes the connected components of a compatibility graph
/// on the features of all maps, whose edges connect two features of different
/// maps when each lies in the other's RT and m/z tolerance window (and,
/// optionally, their intensities are within a fold change). Components that
/// are large enough and do not hold too many features of one map are taken as
/// sets of corresponding features. For each map, a LOWESS transformation is
/// then fitted from each of its features in such a component to the
/// component's average retention time, so that the deviation from the mean
/// retention time within the components is minimised.
///
/// The source's class documentation also requires the features of a
/// component to share one charge state. Its `filterCCs_` collects the charges
/// but its `continue` on a second charge only continues the inner loop, so
/// the charge never removes a component; this is reproduced.
///
/// The parameters are those `FeatureGroupingAlgorithmKD` passes: `warp:rt_tol`,
/// `warp:mz_tol`, `mz_unit`, `warp:max_pairwise_log_fc`,
/// `warp:min_rel_cc_size`, `warp:max_nr_conflicts` and the `LOWESS:` section.
/// The algorithm is not a `DefaultParamHandler` in the source and has no
/// defaults of its own.
#[derive(Clone, Debug)]
pub struct MapAlignmentAlgorithmKD {
    fit_data: Vec<Vec<DataPoint>>,
    transformations: Vec<LowessModel>,
    param: Param,
    tolerances: Option<Tolerances>,
}

impl MapAlignmentAlgorithmKD {
    /// A new algorithm for `num_maps` maps with the given parameters (source
    /// constructor).
    ///
    /// As in the source, an empty `param` is accepted and reads nothing; the
    /// source's tolerance members then stay uninitialised, and
    /// [`MapAlignmentAlgorithmKD::add_rt_fit_data`] refuses to run.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `param` is not empty and one of
    /// `warp:rt_tol`, `warp:mz_tol`, `mz_unit` and `warp:max_pairwise_log_fc`
    /// is missing (the source's `Exception::ElementNotFound`) or not a
    /// number where one is read (`Exception::ConversionError`).
    pub fn new(num_maps: usize, param: &Param) -> Result<Self> {
        let tolerances = if param.is_empty() {
            None
        } else {
            Some(Tolerances {
                rt_tol_secs: float(param, "warp:rt_tol")?,
                mz_tol: float(param, "warp:mz_tol")?,
                // `toString() == "ppm"`: only the string value `ppm` matches.
                mz_ppm: matches!(param.value("mz_unit")?, ParamValue::String(unit) if unit == "ppm"),
                max_pairwise_log_fc: float(param, "warp:max_pairwise_log_fc")?,
            })
        };
        Ok(Self {
            fit_data: vec![Vec::new(); num_maps],
            transformations: Vec::new(),
            param: param.clone(),
            tolerances,
        })
    }

    /// The number of maps the algorithm was built for.
    pub fn num_maps(&self) -> usize {
        self.fit_data.len()
    }

    /// The collected fit data per map: pairs of a feature's retention time
    /// and the average retention time of its component (source member
    /// `fit_data_`).
    pub fn fit_data(&self) -> &[Vec<DataPoint>] {
        &self.fit_data
    }

    /// The fitted LOWESS transformations, one per map, after
    /// [`MapAlignmentAlgorithmKD::fit_lowess`] (source member
    /// `transformations_`).
    pub fn transformations(&self) -> &[LowessModel] {
        &self.transformations
    }

    /// Compute the data points of the RT transformation in the current
    /// `kd_data` and add them to the fit data (source `addRTFitData`).
    ///
    /// The connected components of the compatibility graph are filtered
    /// (see [`MapAlignmentAlgorithmKD::filter_ccs`]); for each remaining
    /// component, the average of its members' retention times is computed,
    /// and each member adds `(its RT, the average)` to its map's fit data, in
    /// component order and, within a component, in feature order.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the algorithm was built without
    /// parameters, when `kd_data` holds a feature of a map past
    /// [`MapAlignmentAlgorithmKD::num_maps`] (the source writes past the end
    /// of `fit_data_`), and the errors of
    /// [`MapAlignmentAlgorithmKD::filter_ccs`] and of the queries.
    pub fn add_rt_fit_data<K: KdFeatureMaps + ?Sized>(&mut self, kd_data: &K) -> Result<()> {
        let ccs = self.get_ccs(kd_data)?;
        let filtered = self.filter_ccs(kd_data, &ccs)?;
        drop(ccs);
        let mut averages = BTreeMap::new();
        for (&index, cc) in &filtered {
            let mut average = 0.0;
            for &i in cc {
                average += kd_data.rt(i);
            }
            average /= cc.len() as f64;
            averages.insert(index, average);
        }
        let num_maps = self.fit_data.len();
        for (index, cc) in &filtered {
            let average = averages[index];
            for &i in cc {
                let map = kd_data.map_index(i);
                let data = self.fit_data.get_mut(map).ok_or_else(|| {
                    Error::InvalidValue(format!(
                        "feature {i} belongs to map {map}, but the alignment has {num_maps} maps"
                    ))
                })?;
                if data.len() >= MAX_FIT_POINTS {
                    return Err(Error::InvalidValue(format!(
                        "RT fit data exceeds {MAX_FIT_POINTS} points for map {map}"
                    )));
                }
                data.push(DataPoint::new(kd_data.rt(i), average));
            }
        }
        Ok(())
    }

    /// Fit a LOWESS model to each map's fit data (source `fitLOWESS`), with
    /// the parameters of the `LOWESS:` section, replacing any earlier fit.
    ///
    /// A map with fewer than 50 data points gets a model fitted to the
    /// identity points `(0, 0)`, `(1, 1)`, `(1e6, 1e6)` instead, and the
    /// source's warning, which is returned.
    ///
    /// # Errors
    ///
    /// Returns the errors of reading the `LOWESS:` section and of the LOWESS
    /// fit (the source's constructor exceptions); the previous models are
    /// kept then.
    pub fn fit_lowess(&mut self) -> Result<Vec<String>> {
        let mut warnings = Vec::new();
        let mut transformations = Vec::with_capacity(self.fit_data.len());
        for (i, data) in self.fit_data.iter().enumerate() {
            // The source reads the section inside the loop, so nothing is
            // read for zero maps.
            let lowess_param = self.param.copy("LOWESS:", true)?;
            let config = ModelConfig::from_param("lowess", &lowess_param)?;
            let n = data.len();
            let identity;
            let points = if n < 50 {
                warnings.push(format!(
                    "Warning: Only {n} data points for LOWESS fit of map {i}. Consider adjusting RT or m/z tolerance or max_pairwise_log_fc, decreasing min_rel_cc_size, or increasing max_nr_conflicts."
                ));
                identity = [
                    DataPoint::new(0.0, 0.0),
                    DataPoint::new(1.0, 1.0),
                    DataPoint::new(1e6, 1e6),
                ];
                &identity[..]
            } else {
                &data[..]
            };
            transformations.push(fit(points, &config)?);
        }
        self.transformations = transformations;
        Ok(warnings)
    }

    /// Transform the retention times of `kd_data` with the fitted models and
    /// rebuild its tree (source `transform`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] before
    /// [`MapAlignmentAlgorithmKD::fit_lowess`] has run (the source
    /// dereferences null model pointers), and the errors of
    /// [`KdFeatureMaps::apply_transformations`].
    pub fn transform<K: KdFeatureMaps + ?Sized>(&self, kd_data: &mut K) -> Result<()> {
        if self.transformations.len() != self.fit_data.len() {
            return Err(Error::InvalidValue(
                "the LOWESS transformations have not been fitted".into(),
            ));
        }
        kd_data.apply_transformations(&self.transformations)
    }

    /// Label every feature with its connected component (source
    /// `computeCCs_`); returns the labels and the number of components.
    ///
    /// Components are found by breadth-first search, seeded at the
    /// lowest-indexed feature not yet visited, so the labels do not depend on
    /// the order in which [`KdFeatureMaps::neighborhood`] returns its result.
    /// Neighbours are queried with the tolerances and `max_pairwise_log_fc`
    /// of the parameters, excluding features of the same map.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the algorithm was built without
    /// parameters or a neighbour index is out of range, and the errors of
    /// the queries.
    pub fn compute_ccs<K: KdFeatureMaps + ?Sized>(
        &self,
        kd_data: &K,
    ) -> Result<(Vec<usize>, usize)> {
        let t = self.tolerances.ok_or_else(|| {
            Error::InvalidValue(
                "MapAlignmentAlgorithmKD was built without parameters; its tolerances are unset"
                    .into(),
            )
        })?;
        let n = kd_data.len();
        let mut labels = vec![usize::MAX; n];
        let mut visited = vec![false; n];
        let mut queue = VecDeque::new();
        let mut cc_index = 0;
        for seed in 0..n {
            if visited[seed] {
                continue;
            }
            visited[seed] = true;
            queue.push_back(seed);
            while let Some(i) = queue.pop_front() {
                labels[i] = cc_index;
                let neighbours = kd_data.neighborhood(
                    i,
                    t.rt_tol_secs,
                    t.mz_tol,
                    t.mz_ppm,
                    false,
                    t.max_pairwise_log_fc,
                )?;
                for j in neighbours {
                    let seen = visited.get_mut(j).ok_or_else(|| {
                        Error::InvalidValue(format!("neighbour index {j} is out of range"))
                    })?;
                    if !*seen {
                        *seen = true;
                        queue.push_back(j);
                    }
                }
            }
            cc_index += 1;
        }
        Ok((labels, cc_index))
    }

    /// The connected components, each a list of feature indices in ascending
    /// order, keyed by component label (source `getCCs_`).
    ///
    /// # Errors
    ///
    /// As [`MapAlignmentAlgorithmKD::compute_ccs`].
    pub fn get_ccs<K: KdFeatureMaps + ?Sized>(
        &self,
        kd_data: &K,
    ) -> Result<BTreeMap<usize, Vec<usize>>> {
        let (labels, _) = self.compute_ccs(kd_data)?;
        let mut result: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (i, label) in labels.into_iter().enumerate() {
            result.entry(label).or_default().push(i);
        }
        Ok(result)
    }

    /// Keep the components with at least `max(2, warp:min_rel_cc_size *
    /// num_maps)` features (truncated to an integer) and, unless
    /// `warp:max_nr_conflicts` is -1, at most that many features from a map
    /// already seen in the component (source `filterCCs_`).
    ///
    /// The source's charge check does not remove anything (see the type's
    /// documentation), so the charges are not read.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `warp:min_rel_cc_size` or
    /// `warp:max_nr_conflicts` is missing or of the wrong type (the source
    /// reads both here, and throws), or `warp:max_nr_conflicts` does not fit
    /// an `int`.
    pub fn filter_ccs<K: KdFeatureMaps + ?Sized>(
        &self,
        kd_data: &K,
        ccs: &BTreeMap<usize, Vec<usize>>,
    ) -> Result<BTreeMap<usize, Vec<usize>>> {
        let relative = float(&self.param, "warp:min_rel_cc_size")? * self.fit_data.len() as f64;
        // `std::max(2.0, relative)` keeps 2.0 unless `2.0 < relative`, also
        // for NaN; the conversion to `Size` truncates.
        let min_size = if 2.0 < relative { relative } else { 2.0 } as usize;
        let max_nr_conflicts = i32::try_from(self.param.value("warp:max_nr_conflicts")?.to_i64()?)
            .map_err(|_| {
                Error::InvalidValue("parameter 'warp:max_nr_conflicts' does not fit an int".into())
            })?;
        let mut filtered = BTreeMap::new();
        for (&index, cc) in ccs {
            if cc.len() < min_size {
                continue;
            }
            let mut passes = true;
            if max_nr_conflicts != -1 {
                let mut maps = BTreeSet::new();
                let mut conflicts = 0;
                for &i in cc {
                    if !maps.insert(kd_data.map_index(i)) {
                        conflicts += 1;
                        if conflicts > max_nr_conflicts {
                            passes = false;
                            break;
                        }
                    }
                }
            }
            if passes {
                filtered.insert(index, cc.clone());
            }
        }
        Ok(filtered)
    }
}

/// Fit the `LOWESS:` configuration, as `new TransformationModelLowess(data,
/// lowess_param)` does.
fn fit(data: &[DataPoint], config: &ModelConfig) -> Result<LowessModel> {
    match config {
        ModelConfig::Lowess(options) => LowessModel::fit(data, *options),
        ModelConfig::LowessAuto(options, auto) => {
            Ok(LowessModel::fit_auto(data, *options, auto)?.0)
        }
        _ => Err(Error::InvalidValue(
            "the LOWESS parameters did not give a LOWESS model".into(),
        )),
    }
}

/// The kd-tree of the source: its retention times are the cached ones
/// ([`KDTreeFeatureMaps::rt`]), the neighbourhood is
/// [`KDTreeFeatureMaps::get_neighborhood`], and applying the transformations
/// replaces the cached times and rebalances the tree (`applyTransformations`
/// followed by `optimizeTree`, as `MapAlignmentAlgorithmKD::transform` calls
/// them).
impl KdFeatureMaps for KDTreeFeatureMaps<'_> {
    fn len(&self) -> usize {
        self.size()
    }

    fn rt(&self, index: usize) -> f64 {
        KDTreeFeatureMaps::rt(self, index).unwrap_or(f64::NAN)
    }

    fn charge(&self, index: usize) -> i32 {
        KDTreeFeatureMaps::charge(self, index).unwrap_or(0)
    }

    fn map_index(&self, index: usize) -> usize {
        KDTreeFeatureMaps::map_index(self, index).unwrap_or(usize::MAX)
    }

    fn neighborhood(
        &self,
        index: usize,
        rt_tol: f64,
        mz_tol: f64,
        mz_ppm: bool,
        include_features_from_same_map: bool,
        max_pairwise_log_fc: f64,
    ) -> Result<Vec<usize>> {
        let mut result = Vec::new();
        self.get_neighborhood(
            index,
            &mut result,
            rt_tol,
            mz_tol,
            mz_ppm,
            include_features_from_same_map,
            max_pairwise_log_fc,
        )?;
        Ok(result)
    }

    fn apply_transformations(&mut self, transformations: &[LowessModel]) -> Result<()> {
        KDTreeFeatureMaps::apply_transformations(self, |map, rt| {
            transformations
                .get(map)
                .ok_or_else(|| {
                    Error::InvalidValue(format!("no LOWESS transformation for map {map}"))
                })?
                .apply(rt)
        })?;
        self.optimize_tree()
    }
}
