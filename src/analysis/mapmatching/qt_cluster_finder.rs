// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A variant of QT clustering for the detection of feature groups.
//!
//! Native equivalent of `ANALYSIS/MAPMATCHING/QTClusterFinder.h`. See
//! `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! All features of all input maps are accumulated, and every feature is a
//! potential cluster centre. For every centre, its nearest neighbours from the
//! other maps are added to the potential cluster. Then the cluster of highest
//! quality is extracted repeatedly, and the remaining clusters are updated.
//!
//! # Properties affecting the grouping
//!
//! To be included in a cluster, a feature must
//!
//! - differ from the centre in RT and m/z by less than
//!   `distance_RT:max_difference` and `distance_MZ:max_difference`;
//! - match the centre's charge, unless `ignore_charge` is set;
//! - with `use_identifications`, when both it and the centre are annotated,
//!   share an annotation with the centre.
//!
//! Every cluster holds at most one feature per input map: the one closest to
//! the centre that meets the criteria and does not belong to a better cluster.
//! Closeness is the
//! [`FeatureDistance`](crate::analysis::mapmatching::feature_distance::FeatureDistance),
//! and the cluster quality comes from the number of elements and their
//! distances to the centre (see
//! [`QTCluster`](crate::analysis::mapmatching::qt_cluster::QTCluster)).
//!
//! # Optimisation
//!
//! Features are hashed into a two-dimensional grid, and clustering runs in one
//! round, updating the affected clusters after each extraction. The source's
//! look-up table of feature distances is not used by its implementation either.
//!
//! # Order
//!
//! The source keeps its clusters in a `boost::heap::fibonacci_heap` ordered by
//! quality, visits grid cells in Boost's hash order, and keeps neighbours in
//! hash maps. Here the heap is a binary heap of (quality, cluster id) entries
//! whose ties go to the **lower** cluster id, the grid is visited in ascending
//! cell order ([`HashGrid`](crate::math::ml::clustering::hash_grid::HashGrid)), and neighbours are ordered by map index. Clusters
//! of exactly equal quality that compete for the same features are therefore
//! resolved deterministically, where the source resolves them by its heap's
//! internal layout; the mirror-image clusters of a feature pair, the common case
//! of a tie, contain the same elements either way.

use std::collections::{BTreeMap, BTreeSet, BinaryHeap};

use super::base_group_finder::BaseGroupFinder;
use super::feature_distance::FeatureDistance;
use super::grid_feature::GridFeature;
use super::params::{add_float, add_int, add_string, bad, flag, flag_text, float, int};
use super::qt_cluster::{BulkData, Element, QTCluster};
use crate::chemistry::{EmpiricalFormula, element_table};
use crate::concept::constants::user_param::{
    ADDUCT_GROUP, DC_CHARGE_ADDUCTS, IIMN_BEST_ION, IIMN_LINKED_GROUPS,
};
use crate::error::{Error, Result};
use crate::kernel::features::{
    BaseFeature, ConsensusFeature, ConsensusMap, FeatureHandle, FeatureMap,
};
use crate::math::ml::clustering::hash_grid::HashGrid;
use crate::metadata::MetaValue;
use crate::param::{DefaultParamHandler, Param};

/// Parameter-handler name of the source class.
pub const QT_CLUSTER_FINDER_NAME: &str = "QTClusterFinder";

/// Largest number of features, summed over all input maps, one run accepts.
/// Checked before anything is allocated; the source has no ceiling.
pub const MAX_FEATURES: usize = 1 << 24;

/// An input map of the cluster finder: a feature map or a consensus map (the
/// source's `MapType` template parameter).
pub trait QtInputMap: Sized {
    /// The features, as base features, in map order.
    fn base_features(&self) -> Vec<&BaseFeature>;
    /// The maximum intensity and m/z of the map's stored ranges after
    /// `updateRanges()`, or `None` for an empty dimension.
    ///
    /// # Errors
    ///
    /// The map's range validation error.
    fn max_intensity_and_mz(&self) -> Result<(Option<f64>, Option<f64>)>;
    /// A new map of the features whose m/z lies in `[start, end)` (the
    /// source's per-partition `push_back` copy).
    fn partition(&self, start: f64, end: f64) -> Self;
}

impl QtInputMap for FeatureMap {
    fn base_features(&self) -> Vec<&BaseFeature> {
        self.features.iter().map(|feature| &feature.base).collect()
    }
    fn max_intensity_and_mz(&self) -> Result<(Option<f64>, Option<f64>)> {
        let ranges = self.ranges()?;
        Ok((ranges.intensity.map(|r| r.max), ranges.mz.map(|r| r.max)))
    }
    fn partition(&self, start: f64, end: f64) -> Self {
        FeatureMap::from_features(
            self.features
                .iter()
                .filter(|f| f.mz >= start && f.mz < end)
                .cloned()
                .collect(),
        )
    }
}

impl QtInputMap for ConsensusMap {
    fn base_features(&self) -> Vec<&BaseFeature> {
        self.features.iter().map(|feature| &feature.base).collect()
    }
    fn max_intensity_and_mz(&self) -> Result<(Option<f64>, Option<f64>)> {
        let ranges = self.ranges()?;
        Ok((ranges.intensity.map(|r| r.max), ranges.mz.map(|r| r.max)))
    }
    fn partition(&self, start: f64, end: f64) -> Self {
        ConsensusMap::from_features(
            self.features
                .iter()
                .filter(|f| f.mz >= start && f.mz < end)
                .cloned()
                .collect(),
        )
    }
}

/// QT cluster finder (source `QTClusterFinder`).
#[derive(Clone, Debug, PartialEq)]
pub struct QTClusterFinder {
    /// Never link features annotated with different peptides; only the best
    /// hit per identification counts (`use_identifications`, default false).
    pub use_identifications: bool,
    /// Number of m/z partitions (`nr_partitions`, default 100, minimum 1).
    pub nr_partitions: i64,
    /// With IDs: how many differences from matching IDs an RT bin collects
    /// before a linking tolerance is computed for it (`min_nr_diffs_per_bin`,
    /// default 50, minimum 5).
    pub min_nr_diffs_per_bin: i64,
    /// With IDs: the minimum score of an ID trusted for tolerance calculation
    /// (`min_IDscore_forTolCalc`, default 1).
    pub min_id_score_for_tol_calc: f64,
    /// With IDs: the distance penalty for unidentified features, in `[0, 1]`
    /// (`noID_penalty`, default 0).
    pub no_id_penalty: f64,
    /// The distance parameters. Its `max_intensity` and `force_constraints`
    /// are replaced for every partition (the source builds
    /// `FeatureDistance(max_intensity, true)`).
    pub distance: FeatureDistance,
    /// Start RT of each bin and its linking tolerance for unidentified
    /// features (source `bin_tolerances_`). Like the source member it is never
    /// cleared, so repeated runs of one finder keep earlier bins; a bin start
    /// already present keeps its tolerance.
    bin_tolerances: Vec<(f64, f64)>,
}

impl Default for QTClusterFinder {
    fn default() -> Self {
        Self {
            use_identifications: false,
            nr_partitions: 100,
            min_nr_diffs_per_bin: 50,
            min_id_score_for_tol_calc: 1.0,
            no_id_penalty: 0.0,
            distance: FeatureDistance::default(),
            bin_tolerances: Vec::new(),
        }
    }
}

/// The per-partition settings of source `setParameters_`.
struct RunSettings {
    max_diff_rt: f64,
    max_diff_mz: f64,
    distance: FeatureDistance,
}

/// One heap entry: a snapshot of a cluster's quality. Newer snapshots replace
/// older ones through `version`.
#[derive(Clone, Copy, Debug)]
struct HeapEntry {
    quality: f64,
    id: usize,
    version: u64,
}

impl PartialEq for HeapEntry {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}
impl Eq for HeapEntry {}
impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HeapEntry {
    /// Higher quality first, then the lower cluster id.
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.quality
            .total_cmp(&other.quality)
            .then_with(|| other.id.cmp(&self.id))
            .then_with(|| self.version.cmp(&other.version))
    }
}

/// Identity of a grid feature: `(map index, feature index)`.
type FeatureKey = (usize, usize);

fn key(element: &Element<'_>) -> FeatureKey {
    (element.feature.map_index(), element.feature.feature_index())
}

impl QTClusterFinder {
    /// The source defaults tree: `use_identifications`, `nr_partitions`,
    /// `min_nr_diffs_per_bin`, `min_IDscore_forTolCalc`, `noID_penalty` and
    /// the [`FeatureDistance::defaults`] at the root.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters in the layout of [`Self::defaults`], carrying the
    /// current values.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let mut p = Param::new();
        add_string(
            &mut p,
            "use_identifications",
            flag_text(self.use_identifications),
            "Never link features that are annotated with different peptides (only the best hit per peptide identification is taken into account).",
            false,
            &["true", "false"],
        )?;
        add_int(
            &mut p,
            "nr_partitions",
            self.nr_partitions,
            "How many partitions in m/z space should be used for the algorithm (more partitions means faster runtime and more memory efficient execution).",
            false,
            Some(1),
            None,
        )?;
        add_int(
            &mut p,
            "min_nr_diffs_per_bin",
            self.min_nr_diffs_per_bin,
            "If IDs are used: How many differences from matching IDs should be used to calculate a linking tolerance for unIDed features in an RT region. RT regions will be extended until that number is reached.",
            false,
            Some(5),
            None,
        )?;
        add_float(
            &mut p,
            "min_IDscore_forTolCalc",
            self.min_id_score_for_tol_calc,
            "If IDs are used: What is the minimum score of an ID to assume a reliable match for tolerance calculation. Check your current score type!",
            false,
            None,
            None,
        )?;
        add_float(
            &mut p,
            "noID_penalty",
            self.no_id_penalty,
            "If IDs are used: For the normalized distances, how high should the penalty for missing IDs be? 0 = no bias, 1 = IDs inside the max tolerances always preferred (even if much further away).",
            false,
            Some(0.0),
            Some(1.0),
        )?;
        p.insert("", &self.distance.to_param()?)?;
        Ok(p)
    }

    /// Build from parameters, as the constructor followed by `setParameters`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its restriction (source `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`Self::from_param`], also returning the unknown-parameter warnings.
    ///
    /// # Errors
    ///
    /// As [`Self::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(QT_CLUSTER_FINDER_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    /// Read a validated tree whose every default is present.
    pub(crate) fn from_complete_param(p: &Param) -> Result<Self> {
        Ok(Self {
            use_identifications: flag(p, "use_identifications")?,
            nr_partitions: int(p, "nr_partitions")?,
            min_nr_diffs_per_bin: int(p, "min_nr_diffs_per_bin")?,
            min_id_score_for_tol_calc: float(p, "min_IDscore_forTolCalc")?,
            no_id_penalty: float(p, "noID_penalty")?,
            distance: FeatureDistance::from_complete_param(1.0, true, p)?,
            bin_tolerances: Vec::new(),
        })
    }

    /// A finder with the same parameters and no RT bins, as a freshly
    /// constructed source finder has.
    pub(crate) fn fresh(&self) -> Self {
        Self {
            bin_tolerances: Vec::new(),
            ..self.clone()
        }
    }

    /// The RT bins and their linking tolerances computed so far (source
    /// `bin_tolerances_`), ascending by bin start.
    pub fn bin_tolerances(&self) -> &[(f64, f64)] {
        &self.bin_tolerances
    }

    /// Run on feature maps (source `run(const std::vector<FeatureMap>&, ...)`).
    ///
    /// The source requires the input maps' ranges to be current; here they
    /// are computed.
    ///
    /// # Errors
    ///
    /// As [`Self::run_maps`].
    pub fn run_feature_maps(
        &mut self,
        input: &[FeatureMap],
        result: &mut ConsensusMap,
    ) -> Result<()> {
        self.run_maps(input, result)
    }

    /// Run on feature or consensus maps (source template `run_`).
    ///
    /// With `use_identifications`, first derives RT linking tolerances for
    /// unidentified features from the RT spread of confidently identified
    /// peptides (see [`Self::bin_tolerances`]). Then, unless every map is
    /// empty, the features are split into m/z partitions at gaps wider than
    /// the m/z tolerance (at most about `nr_partitions`, never inside a
    /// possible cluster) and each partition is clustered separately; the
    /// consensus features are appended to `result` in extraction order, after
    /// its features are cleared (source `clear(false)`: metadata, headers and
    /// identifications stay).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for fewer than two input maps (source
    /// `Exception::IllegalArgument`, raised per partition), when a partition's
    /// maximum m/z is outside `[1e-16, 1e16]` or its maximum intensity above
    /// `1e16` (the source's "Has 'updateRanges' been called" check), when the
    /// feature count exceeds [`MAX_FEATURES`], when a grid cell index does not
    /// fit, and the [`FeatureDistance::distance`] and consensus errors. On
    /// error `result` is unchanged; the bin tolerances may have been extended.
    pub fn run_maps<M: QtInputMap>(
        &mut self,
        input: &[M],
        result: &mut ConsensusMap,
    ) -> Result<()> {
        let total: usize = input.iter().map(|m| m.base_features().len()).sum();
        if total > MAX_FEATURES {
            return Err(bad(format!(
                "{total} features exceed the {MAX_FEATURES}-feature ceiling"
            )));
        }
        if self.use_identifications {
            self.compute_bin_tolerances(input)?;
        }

        let mut massrange: Vec<f64> = input
            .iter()
            .flat_map(|m| m.base_features().into_iter().map(|f| f.mz))
            .collect();
        massrange.sort_by(f64::total_cmp);
        let mut features = Vec::new();
        if massrange.is_empty() {
            result.features = features;
            return Ok(());
        }

        if self.nr_partitions == 1 {
            self.run_internal(input, &mut features)?;
        } else {
            let unit_ppm = self.distance.mz_ppm;
            let mz_tol = self.distance.mz.max_difference;
            // Source `setParameters_(1, 1)`: in ppm the tolerance is taken at
            // m/z 1 until the loop below replaces it.
            let mut massrange_diff = if unit_ppm { mz_tol * 1e-6 } else { mz_tol };
            let count = i32::try_from(massrange.len())
                .map_err(|_| bad("too many features for the source's int partition count"))?;
            let partitions = i32::try_from(self.nr_partitions)
                .map_err(|_| bad("nr_partitions does not fit the source's int"))?;
            let pts_per_partition = (count / partitions) as usize;
            let mut boundaries = vec![massrange[0]];
            for j in 0..massrange.len() - 1 {
                if unit_ppm {
                    massrange_diff = mz_tol * 1e-6 * massrange[j + 1];
                }
                if (massrange[j] - massrange[j + 1]).abs() > massrange_diff
                    && j >= boundaries.len() * pts_per_partition
                {
                    boundaries.push((massrange[j] + massrange[j + 1]) / 2.0);
                }
            }
            boundaries.push(massrange[massrange.len() - 1] + 1.0);
            for window in boundaries.windows(2) {
                let parts: Vec<M> = input
                    .iter()
                    .map(|m| m.partition(window[0], window[1]))
                    .collect();
                self.run_internal(&parts, &mut features)?;
            }
        }
        result.features = features;
        Ok(())
    }

    /// Source `setParameters_`: the per-partition distance and grid sizes.
    fn settings(&self, max_intensity: f64, max_mz: f64) -> Result<RunSettings> {
        if !(1e-16..=1e16).contains(&max_mz) || max_intensity > 1e16 {
            return Err(bad(format!(
                "Maximum m/z or intensity out of range (m/z: {max_mz}, intensity: {max_intensity}). Has 'updateRanges' been called on the input maps?"
            )));
        }
        let mut max_diff_mz = self.distance.mz.max_difference;
        if self.distance.mz_ppm {
            max_diff_mz *= max_mz * 1e-6;
        }
        let mut distance = self.distance.clone();
        distance.max_intensity = max_intensity;
        distance.force_constraints = true;
        Ok(RunSettings {
            max_diff_rt: self.distance.rt.max_difference,
            max_diff_mz,
            distance,
        })
    }

    /// Source `run_internal_`: cluster one partition.
    fn run_internal<M: QtInputMap>(
        &self,
        input: &[M],
        out: &mut Vec<ConsensusFeature>,
    ) -> Result<()> {
        let num_maps = input.len();
        if num_maps < 2 {
            return Err(bad("At least two input maps required"));
        }
        let mut max_intensity = f64::MIN;
        let mut max_mz = f64::MIN;
        for map in input {
            let (intensity, mz) = map.max_intensity_and_mz()?;
            if let Some(value) = intensity {
                max_intensity = max_intensity.max(value);
            }
            if let Some(value) = mz {
                max_mz = max_mz.max(value);
            }
        }
        let settings = self.settings(max_intensity, max_mz)?;

        let features: Vec<Vec<&BaseFeature>> =
            input.iter().map(QtInputMap::base_features).collect();
        let mut grid_features = Vec::new();
        for (map_index, map) in features.iter().enumerate() {
            for (feature_index, feature) in map.iter().enumerate() {
                grid_features.push(GridFeature::new(feature, map_index, feature_index));
            }
        }
        let mut grid: HashGrid<usize> =
            HashGrid::new([settings.max_diff_rt, settings.max_diff_mz].into());
        for (index, feature) in grid_features.iter().enumerate() {
            grid.insert([feature.rt(), feature.mz()].into(), index)?;
        }

        let mut state = Clustering {
            finder: self,
            settings: &settings,
            grid: &grid,
            grid_features: &grid_features,
            already_used: BTreeSet::new(),
            clusters: Vec::new(),
            versions: Vec::new(),
            heap: BinaryHeap::new(),
            element_mapping: BTreeMap::new(),
        };
        state.compute_clustering(num_maps)?;
        while let Some(feature) = state.make_consensus_feature()? {
            out.push(feature);
        }
        Ok(())
    }

    /// Source `distIsOutlier_`: whether `dist` reaches the tolerance of the
    /// bin containing `rt` (the first bin for an RT before it). No bins means
    /// no outlier.
    fn dist_is_outlier(&self, dist: f64, rt: f64) -> bool {
        if self.bin_tolerances.is_empty() {
            return false;
        }
        // `upper_bound(rt)`: the first bin start greater than `rt`.
        let upper = self
            .bin_tolerances
            .partition_point(|&(start, _)| start <= rt);
        let index = if upper == 0 { 0 } else { upper - 1 };
        dist >= self.bin_tolerances[index].1
    }

    /// The RT-tolerance estimation of source `run_` with IDs.
    fn compute_bin_tolerances<M: QtInputMap>(&mut self, input: &[M]) -> Result<()> {
        let mut ided_feat_rts: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        let mut min_rt = f64::MAX;
        for map in input {
            for feature in map.base_features() {
                if feature.rt < min_rt {
                    min_rt = feature.rt;
                }
                if feature.peptide_identifications.is_empty() {
                    continue;
                }
                let mut copy = feature.clone();
                copy.sort_peptide_identifications()?;
                let first = &copy.peptide_identifications[0];
                let Some(hit) = first.hits.first() else {
                    continue;
                };
                let higher = first.higher_score_better;
                let score = hit.score;
                if (score > self.min_id_score_for_tol_calc && higher)
                    || (score < self.min_id_score_for_tol_calc && !higher)
                {
                    let key = format!("{}/{}", hit.sequence, copy.charge);
                    ided_feat_rts.entry(key).or_default().push(copy.rt);
                }
            }
        }

        let mut medians_diffs: Vec<(f64, Vec<f64>)> = Vec::with_capacity(ided_feat_rts.len());
        for (_, mut rts) in ided_feat_rts {
            rts.sort_by(f64::total_cmp);
            let middle = rts.len() / 2;
            let median = rts[middle];
            let diffs = rts
                .iter()
                .enumerate()
                .filter(|&(i, _)| i != middle)
                .map(|(_, rt)| (rt - median).abs())
                .collect();
            medians_diffs.push((median, diffs));
        }
        medians_diffs.sort_by(|a, b| a.0.total_cmp(&b.0).then_with(|| lexicographic(&a.1, &b.1)));

        let min_nr = usize::try_from(self.min_nr_diffs_per_bin).unwrap_or(0);
        let mut tmp_diffs: Vec<f64> = Vec::new();
        let mut last_tmp_diffs: Vec<f64> = Vec::new();
        let mut start_rt = min_rt;
        for (cnt, (median, diffs)) in medians_diffs.iter().enumerate() {
            if tmp_diffs.len() > min_nr && !diffs.is_empty() {
                tmp_diffs.sort_by(f64::total_cmp);
                self.insert_bin(start_rt, tolerance(&tmp_diffs)?);
                std::mem::swap(&mut tmp_diffs, &mut last_tmp_diffs);
                tmp_diffs.clear();
                if cnt > 0 {
                    start_rt = (median + medians_diffs[cnt - 1].0) / 2.0;
                }
            } else {
                tmp_diffs.extend_from_slice(diffs);
            }
        }
        tmp_diffs.sort_by(f64::total_cmp);
        // `std::merge` of two ascending ranges, the first range winning ties.
        let mut merged = Vec::with_capacity(tmp_diffs.len() + last_tmp_diffs.len());
        let (mut i, mut j) = (0, 0);
        while i < tmp_diffs.len() && j < last_tmp_diffs.len() {
            if last_tmp_diffs[j] < tmp_diffs[i] {
                merged.push(last_tmp_diffs[j]);
                j += 1;
            } else {
                merged.push(tmp_diffs[i]);
                i += 1;
            }
        }
        merged.extend_from_slice(&tmp_diffs[i..]);
        merged.extend_from_slice(&last_tmp_diffs[j..]);
        if !merged.is_empty() {
            self.insert_bin(start_rt, tolerance(&merged)?);
        }
        Ok(())
    }

    /// `std::map::insert`: a bin start already present keeps its tolerance.
    fn insert_bin(&mut self, start: f64, tol: f64) {
        let position = self.bin_tolerances.partition_point(|&(s, _)| s < start);
        if self
            .bin_tolerances
            .get(position)
            .is_some_and(|&(s, _)| s == start)
        {
            return;
        }
        self.bin_tolerances.insert(position, (start, tol));
    }
}

impl BaseGroupFinder for QTClusterFinder {
    /// Run on consensus maps (source `run(const std::vector<ConsensusMap>&, ...)`).
    ///
    /// # Errors
    ///
    /// As [`QTClusterFinder::run_maps`].
    fn run(&mut self, input: &[ConsensusMap], result: &mut ConsensusMap) -> Result<()> {
        self.run_maps(input, result)
    }
}

/// Lexicographic order of two `f64` vectors, as `std::vector`'s `operator<`.
fn lexicographic(a: &[f64], b: &[f64]) -> std::cmp::Ordering {
    for (x, y) in a.iter().zip(b) {
        let order = x.total_cmp(y);
        if order != std::cmp::Ordering::Equal {
            return order;
        }
    }
    a.len().cmp(&b.len())
}

/// Source Double-MAD tolerance: `max(20, q2 + 2 * 1.4826 * (q3 - q2))` with
/// the quantiles of `MATH/MathFunctions.h`.
fn tolerance(sorted: &[f64]) -> Result<f64> {
    let q2 = math_functions_quantile(sorted, 0.5)?;
    let q3 = math_functions_quantile(sorted, 0.75)?;
    let value = q2 + 2. * 1.4826 * (q3 - q2);
    // `std::max(min_tolerance, value)` returns the first unless it is less.
    Ok(if 20.0 < value { value } else { 20.0 })
}

/// `Math::quantile(const T1& x, double q)` of `MATH/MathFunctions.h` (not the
/// iterator overload of `StatisticFunctions.h`): `id = max(0, n*q - 1)`,
/// interpolating between `x[floor(id)]` and `x[ceil(id)]`.
fn math_functions_quantile(x: &[f64], q: f64) -> Result<f64> {
    if x.is_empty() {
        return Err(bad("Quantile requested from empty container."));
    }
    let q = q.clamp(0.0, 1.0);
    let n = x.len() as f64;
    let id = (n * q - 1.0).max(0.0);
    let lo = id.floor();
    let hi = id.ceil();
    let at = |i: f64| {
        x.get(i as usize)
            .copied()
            .ok_or_else(|| bad("quantile index outside the range"))
    };
    let h = id - lo;
    Ok((1.0 - h) * at(lo)? + h * at(hi)?)
}

/// `Adduct::toAdductString(ion_string, charge)`: `[M` then each element,
/// sorted by symbol, as sign, count when above one, and symbol, then `]`, the
/// absolute charge when above one, and its sign.
///
/// # Errors
///
/// The [`EmpiricalFormula::parse`] error, and [`Error::Unsupported`] for a
/// formula with isotope-labelled atoms, whose source symbols this port does
/// not enumerate.
fn to_adduct_string(ion_string: &str, charge: i32) -> Result<String> {
    let formula = EmpiricalFormula::parse(ion_string)?;
    let mut entries: BTreeMap<&str, i32> = BTreeMap::new();
    let mut covered = 0i64;
    for element in element_table() {
        let count = formula.count(element.symbol())?;
        if count != 0 {
            entries.insert(element.symbol(), count);
            covered += i64::from(count);
        }
    }
    if covered != formula.atom_count() {
        return Err(Error::Unsupported(format!(
            "adduct string of isotope-labelled formula '{ion_string}'"
        )));
    }
    let mut s = String::from("[M");
    for (symbol, count) in entries {
        s.push(if count > 0 { '+' } else { '-' });
        if count.unsigned_abs() > 1 {
            s.push_str(&count.unsigned_abs().to_string());
        }
        s.push_str(symbol);
    }
    s.push(']');
    if charge.unsigned_abs() > 1 {
        s.push_str(&charge.unsigned_abs().to_string());
    }
    s.push(if charge >= 0 { '+' } else { '-' });
    Ok(s)
}

/// The working state of one partition's clustering.
struct Clustering<'s, 'a> {
    finder: &'s QTClusterFinder,
    settings: &'s RunSettings,
    grid: &'s HashGrid<usize>,
    grid_features: &'a [GridFeature<'a>],
    /// Features already in a consensus feature (source `already_used_`).
    already_used: BTreeSet<FeatureKey>,
    /// All clusters by id (the source's `BulkData` vector plus heap handles).
    clusters: Vec<QTCluster<'a>>,
    /// The version of each cluster's current heap entry.
    versions: Vec<u64>,
    heap: BinaryHeap<HeapEntry>,
    /// Which clusters contain which feature (source `ElementMapping`).
    element_mapping: BTreeMap<FeatureKey, BTreeSet<usize>>,
}

impl<'a> Clustering<'_, 'a> {
    /// Source `computeClustering_`: one cluster per grid feature, in grid
    /// order, each filled from the neighbouring cells.
    fn compute_clustering(&mut self, num_maps: usize) -> Result<()> {
        let max_distance = 1.0 + self.finder.no_id_penalty;
        let grid_features = self.grid_features;
        let mut centres = Vec::new();
        for (cell, (_, index)) in self.grid.iter() {
            let x = i32::try_from(cell[0]).map_err(|_| bad("grid cell index does not fit"))?;
            let y = i32::try_from(cell[1]).map_err(|_| bad("grid cell index does not fit"))?;
            centres.push((x, y, *index));
        }
        for (id, (x, y, index)) in centres.into_iter().enumerate() {
            let data = BulkData::new(&grid_features[index], num_maps, max_distance, x, y, id);
            let mut cluster = QTCluster::new(data, self.finder.use_identifications);
            self.add_cluster_elements(&mut cluster)?;
            for element in cluster.elements() {
                self.element_mapping
                    .entry(key(&element))
                    .or_default()
                    .insert(id);
            }
            self.heap.push(HeapEntry {
                quality: cluster.current_quality(),
                id,
                version: 0,
            });
            self.clusters.push(cluster);
            self.versions.push(0);
        }
        Ok(())
    }

    /// Source `addClusterElements_`.
    fn add_cluster_elements(&self, cluster: &mut QTCluster<'a>) -> Result<()> {
        cluster.initialize_cluster();
        let x = i64::from(cluster.x_coord());
        let y = i64::from(cluster.y_coord());
        let center = cluster.center_point();
        let grid_features = self.grid_features;
        for i in x - 1..=x + 1 {
            for j in y - 1..=y + 1 {
                let Some(cell) = self.grid.grid_find(&[i, j]) else {
                    continue;
                };
                for (_, index) in cell {
                    let neighbor = &grid_features[*index];
                    if self
                        .already_used
                        .contains(&(neighbor.map_index(), neighbor.feature_index()))
                    {
                        continue;
                    }
                    if std::ptr::eq(center, neighbor) {
                        continue;
                    }
                    let mut dist = self
                        .settings
                        .distance
                        .distance(center.feature(), neighbor.feature())?
                        .1;
                    if dist == FeatureDistance::INFINITY {
                        continue;
                    }
                    if self.finder.use_identifications && neighbor.annotations().is_empty() {
                        let rt_dist = (neighbor.rt() - center.rt()).abs();
                        if self.finder.dist_is_outlier(rt_dist, center.rt()) {
                            continue;
                        }
                        dist += self.finder.no_id_penalty;
                    }
                    cluster.add(neighbor, dist);
                }
            }
        }
        cluster.finalize_cluster();
        Ok(())
    }

    /// The id of the best cluster still in the heap, dropping stale entries.
    fn pop_top(&mut self) -> Option<usize> {
        while let Some(entry) = self.heap.pop() {
            if self.versions[entry.id] == entry.version {
                return Some(entry.id);
            }
        }
        None
    }

    /// Source `makeConsensusFeature_` with `updateClustering_`: extract the
    /// best valid cluster as a consensus feature, or `None` when none is left.
    fn make_consensus_feature(&mut self) -> Result<Option<ConsensusFeature>> {
        let best_id = loop {
            let Some(id) = self.pop_top() else {
                return Ok(None);
            };
            if self.clusters[id].is_invalid() {
                self.remove_from_element_mapping(id);
                continue;
            }
            break id;
        };
        let elements = self.clusters[best_id].elements();
        let quality = self.clusters[best_id].current_quality();
        let feature = self.create_consensus_feature(quality, &elements)?;
        self.update_clustering(&elements, best_id)?;
        Ok(Some(feature))
    }

    /// Source `removeFromElementMapping_`.
    fn remove_from_element_mapping(&mut self, id: usize) {
        for element in self.clusters[id].elements() {
            if let Some(ids) = self.element_mapping.get_mut(&key(&element)) {
                ids.remove(&id);
            }
        }
    }

    /// Source `createConsensusFeature_`.
    fn create_consensus_feature(
        &mut self,
        quality: f64,
        elements: &[Element<'a>],
    ) -> Result<ConsensusFeature> {
        let mut feature = ConsensusFeature::new();
        feature.quality = quality as f32;
        let mut best_quality = 0.0f32;
        let mut linked_groups = Vec::new();
        for element in elements {
            self.already_used.insert(key(element));
            let elem = element.feature.feature();
            feature.insert(FeatureHandle::new(element.map_index as u64, elem))?;
            if let Some(adducts) = elem.metadata.get(DC_CHARGE_ADDUCTS) {
                feature
                    .metadata
                    .insert(elem.unique_id.to_string(), adducts.clone());
                if elem.quality > best_quality {
                    feature
                        .metadata
                        .insert(IIMN_BEST_ION.into(), adducts.clone());
                    best_quality = elem.quality;
                }
            }
            if let Some(group) = elem.metadata.get(ADDUCT_GROUP) {
                linked_groups.push(meta_text(group)?);
            }
        }
        // The source never advances `best_quality_index`, so the first
        // element's adduct, if any, overrides the best one above.
        if let Some(first) = elements.first() {
            let elem = first.feature.feature();
            if let Some(adducts) = elem.metadata.get(DC_CHARGE_ADDUCTS) {
                let text = to_adduct_string(adducts.as_str()?, elem.charge)?;
                feature
                    .metadata
                    .insert(IIMN_BEST_ION.into(), MetaValue::from(text));
            }
        }
        if !linked_groups.is_empty() {
            feature
                .metadata
                .insert(IIMN_LINKED_GROUPS.into(), MetaValue::from(linked_groups));
        }
        feature.compute_consensus()?;
        Ok(feature)
    }

    /// Source `updateClustering_`: remove the best cluster's elements from
    /// every other cluster holding them, refill those clusters, and invalidate
    /// clusters whose centre was taken.
    fn update_clustering(&mut self, elements: &[Element<'a>], best_id: usize) -> Result<()> {
        for element in elements {
            let feature_key = key(element);
            let ids: Vec<usize> = match self.element_mapping.get_mut(&feature_key) {
                Some(ids) => {
                    ids.remove(&best_id);
                    ids.iter().copied().collect()
                }
                None => Vec::new(),
            };
            let mut tmp_mapping: BTreeMap<FeatureKey, BTreeSet<usize>> = BTreeMap::new();
            for id in ids {
                if self.clusters[id].is_invalid() || !self.clusters[id].update(elements) {
                    continue;
                }
                self.remove_from_element_mapping(id);
                let mut cluster = self.clusters[id].clone();
                self.add_cluster_elements(&mut cluster)?;
                self.clusters[id] = cluster;
                self.versions[id] += 1;
                self.heap.push(HeapEntry {
                    quality: self.clusters[id].current_quality(),
                    id,
                    version: self.versions[id],
                });
                for neighbor in self.clusters[id].elements() {
                    tmp_mapping.entry(key(&neighbor)).or_default().insert(id);
                }
            }
            for (feature, ids) in tmp_mapping {
                self.element_mapping.entry(feature).or_default().extend(ids);
            }
        }
        Ok(())
    }
}

/// A meta value as the source's `DataValue::toString()` renders a string.
fn meta_text(value: &MetaValue) -> Result<String> {
    Ok(value.as_str()?.to_string())
}
