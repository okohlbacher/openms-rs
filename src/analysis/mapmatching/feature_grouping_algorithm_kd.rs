// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Feature grouping for unlabeled data with a kd-tree, optionally after an
//! internal LOWESS RT warp.
//!
//! Native equivalent of `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h`
//! (with its `ClusterProxyKD`). See `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! The algorithm takes feature or consensus maps and searches for
//! corresponding (consensus) features across them. The m/z range is cut into
//! partitions at gaps wider than the larger of the warp and link m/z
//! tolerances (about `nr_partitions` of them), so that no cluster can reach
//! across a boundary. With `warp:enabled`, a
//! [`MapAlignmentAlgorithmKD`](crate::analysis::mapmatching::alignment_kd::MapAlignmentAlgorithmKD) collects RT fit data from every partition and
//! fits one LOWESS transformation per map; the linking then runs on the
//! warped retention times (reported RTs stay the original ones).
//!
//! Linking is greedy: every feature is the centre of a potential cluster that
//! takes, from each map, the compatible unassigned feature closest to the
//! centre by [`FeatureDistance`](crate::analysis::mapmatching::feature_distance::FeatureDistance) (RT and m/z within `link:rt_tol` and
//! `link:mz_tol` of the centre, charges and adducts compatible per
//! `link:charge_merging` and `link:adduct_merging`). The largest cluster with
//! the smallest average distance (then the highest centre index) becomes a
//! consensus feature, its features are marked assigned, and the potential
//! clusters of their neighbours are recomputed, until every feature is
//! assigned.
//!
//! # Differences from the source
//!
//! - When the LOWESS fit fails, the source logs the error and returns with
//!   the output cleared and not post-processed; this port returns the error
//!   and leaves the output unchanged.
//! - A candidate whose distance to the centre is not a number cannot be
//!   chosen (the source then records the index `numeric_limits<Size>::max()`
//!   and dereferences it); this port returns an error.
//! - The warp is [`MapAlignmentAlgorithmKD`](crate::analysis::mapmatching::alignment_kd::MapAlignmentAlgorithmKD)
//!   on the real [`KDTreeFeatureMaps`](crate::analysis::mapmatching::kd_tree_feature_maps::KDTreeFeatureMaps), fitting `TransformationModelLowess` with
//!   its full parameter set (automatic span selection, `akima` interpolation).
//!
//! The consensus feature's `best ion` annotation reproduces a source quirk:
//! the best-quality index starts at 0, a global index into the partition's
//! features, so when no member qualifies the annotation comes from the
//! partition's first feature if that one carries an adduct.

use std::collections::{BTreeMap, BTreeSet};

use super::alignment_kd::MapAlignmentAlgorithmKD;
use super::feature_distance::FeatureDistance;
use super::feature_grouping_algorithm::{
    FeatureGroupingAlgorithm, InputIdentifications, consensus_map_identifications,
    feature_map_identifications, insert_base_feature, postprocess,
};
use super::kd_tree_feature_maps::KDTreeFeatureMaps;
use super::params::{add_float, add_int, add_string, bad, float, int, string};
use super::qt_cluster_finder::{QtInputMap, meta_text, to_adduct_string};
use crate::analysis::transformations::{LowessAutoSpan, LowessOptions};
use crate::chemistry::EmpiricalFormula;
use crate::concept::constants::user_param::{
    ADDUCT_GROUP, DC_CHARGE_ADDUCTS, IIMN_BEST_ION, IIMN_LINKED_GROUPS,
};
use crate::error::Result;
use crate::kernel::features::{ConsensusFeature, ConsensusMap, FeatureMap};
use crate::metadata::MetaValue;
use crate::param::{DefaultParamHandler, Param, ParamValue};

/// Parameter-handler name of the source class.
pub const FEATURE_GROUPING_ALGORITHM_KD_NAME: &str = "FeatureGroupingAlgorithmKD";

/// A proxy for a potential cluster: its size, the average distance of its
/// features to the centre, and the centre's index (source `ClusterProxyKD`).
///
/// The order is the source's `operator<`, which makes the first element of an
/// ordered set the cluster to take next: larger size first, then smaller
/// average distance, then **higher** centre index. The default proxy has
/// size 0 and is not valid.
#[derive(Clone, Copy, Debug, Default)]
pub struct ClusterProxyKD {
    size: usize,
    avg_distance: f64,
    center_index: usize,
}

impl ClusterProxyKD {
    /// A proxy (source constructor).
    pub fn new(size: usize, avg_distance: f64, center_index: usize) -> Self {
        Self {
            size,
            avg_distance,
            center_index,
        }
    }

    /// The cluster size (source `getSize`).
    pub fn size(&self) -> usize {
        self.size
    }

    /// Whether the proxy describes a cluster: a non-zero size (source
    /// `isValid`).
    pub fn is_valid(&self) -> bool {
        self.size != 0
    }

    /// The average distance to the centre (source `getAvgDistance`).
    pub fn avg_distance(&self) -> f64 {
        self.avg_distance
    }

    /// The centre's index (source `getCenterIndex`).
    pub fn center_index(&self) -> usize {
        self.center_index
    }
}

impl Ord for ClusterProxyKD {
    fn cmp(&self, rhs: &Self) -> std::cmp::Ordering {
        rhs.size
            .cmp(&self.size)
            .then_with(|| {
                self.avg_distance
                    .partial_cmp(&rhs.avg_distance)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| rhs.center_index.cmp(&self.center_index))
    }
}

impl PartialOrd for ClusterProxyKD {
    fn partial_cmp(&self, rhs: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(rhs))
    }
}

impl PartialEq for ClusterProxyKD {
    /// Neither proxy precedes the other (source `operator==`).
    fn eq(&self, rhs: &Self) -> bool {
        self.cmp(rhs) == std::cmp::Ordering::Equal
    }
}

impl Eq for ClusterProxyKD {}

/// Which charges may be linked (`link:charge_merging`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ChargeMerging {
    /// Only equal charges (`Identical`).
    Identical,
    /// Equal charges, or charge zero (unknown) with any (`With_charge_zero`,
    /// the default).
    #[default]
    WithChargeZero,
    /// Any charges (`Any`).
    Any,
}

impl ChargeMerging {
    fn as_str(self) -> &'static str {
        match self {
            Self::Identical => "Identical",
            Self::WithChargeZero => "With_charge_zero",
            Self::Any => "Any",
        }
    }
}

/// Which adducts may be linked (`link:adduct_merging`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AdductMerging {
    /// Only the same adduct, or none on both (`Identical`).
    Identical,
    /// Also features without adduct annotation (`With_unknown_adducts`).
    WithUnknownAdducts,
    /// Any adducts (`Any`, the default).
    #[default]
    Any,
}

impl AdductMerging {
    fn as_str(self) -> &'static str {
        match self {
            Self::Identical => "Identical",
            Self::WithUnknownAdducts => "With_unknown_adducts",
            Self::Any => "Any",
        }
    }
}

/// KD feature grouping (source `FeatureGroupingAlgorithmKD`).
#[derive(Clone, Debug, PartialEq)]
pub struct FeatureGroupingAlgorithmKD {
    /// Warp feature RTs with LOWESS before linking (`warp:enabled`, true).
    pub warp_enabled: bool,
    /// Width of the warp's RT tolerance window in seconds (`warp:rt_tol`,
    /// 100, minimum 0).
    pub warp_rt_tol: f64,
    /// The warp's m/z tolerance in ppm or Da (`warp:mz_tol`, 5, minimum 0).
    pub warp_mz_tol: f64,
    /// Maximum absolute log10 fold change between compatible signals in the
    /// warp's compatibility graph; negative disables the check
    /// (`warp:max_pairwise_log_fc`, advanced, 0.5).
    pub warp_max_pairwise_log_fc: f64,
    /// Minimum connected-component size relative to the number of maps
    /// (`warp:min_rel_cc_size`, advanced, 0.5, 0 to 1).
    pub warp_min_rel_cc_size: f64,
    /// Conflicts (features of one map) a component may hold, -1 for any
    /// (`warp:max_nr_conflicts`, advanced, 0, minimum -1).
    pub warp_max_nr_conflicts: i64,
    /// Width of the linking RT tolerance window in seconds (`link:rt_tol`,
    /// 30, minimum 0).
    pub link_rt_tol: f64,
    /// The linking m/z tolerance in ppm or Da (`link:mz_tol`, 10, minimum 0).
    pub link_mz_tol: f64,
    /// `link:charge_merging`.
    pub link_charge_merging: ChargeMerging,
    /// `link:adduct_merging`.
    pub link_adduct_merging: AdductMerging,
    /// Whether the m/z tolerances are in ppm (`mz_unit` = `ppm`, the
    /// default) rather than Da.
    pub mz_ppm: bool,
    /// Number of partitions in m/z space (`nr_partitions`, 100, minimum 1).
    pub nr_partitions: i64,
    /// The distance parameters: the exponents, weights and intensity log
    /// transformation of [`FeatureDistance`]. Its maximum differences, unit,
    /// `ignore_charge` and `ignore_adduct` are not parameters here; every run
    /// sets them from the link tolerances, `mz_unit` and
    /// `link:charge_merging`.
    pub distance: FeatureDistance,
    /// The LOWESS options of the warp (`LOWESS:span`, `num_iterations`,
    /// `delta`, `interpolation_type`, `extrapolation_type`; all advanced).
    pub lowess: LowessOptions,
    /// The automatic span selection of the warp's LOWESS (`LOWESS:auto_span`
    /// and the `auto_*` parameters), `None` when `auto_span` is `false`.
    pub lowess_auto_span: Option<LowessAutoSpan>,
}

impl Default for FeatureGroupingAlgorithmKD {
    fn default() -> Self {
        let distance = FeatureDistance {
            intensity_weight: 1.0,
            log_transform: true,
            ..FeatureDistance::default()
        };
        Self {
            warp_enabled: true,
            warp_rt_tol: 100.0,
            warp_mz_tol: 5.0,
            warp_max_pairwise_log_fc: 0.5,
            warp_min_rel_cc_size: 0.5,
            warp_max_nr_conflicts: 0,
            link_rt_tol: 30.0,
            link_mz_tol: 10.0,
            link_charge_merging: ChargeMerging::default(),
            link_adduct_merging: AdductMerging::default(),
            mz_ppm: true,
            nr_partitions: 100,
            distance,
            lowess: LowessOptions::default(),
            lowess_auto_span: None,
        }
    }
}

impl FeatureGroupingAlgorithmKD {
    /// The source defaults tree.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters in the layout of [`Self::defaults`], built with the
    /// source constructor's sequence of calls.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let mut p = Param::new();
        add_string(
            &mut p,
            "warp:enabled",
            if self.warp_enabled { "true" } else { "false" },
            "Whether or not to internally warp feature RTs using LOWESS transformation before linking (reported RTs in results will always be the original RTs)",
            false,
            &["true", "false"],
        )?;
        add_float(
            &mut p,
            "warp:rt_tol",
            self.warp_rt_tol,
            "Width of RT tolerance window (sec)",
            false,
            Some(0.0),
            None,
        )?;
        add_float(
            &mut p,
            "warp:mz_tol",
            self.warp_mz_tol,
            "m/z tolerance (in ppm or Da)",
            false,
            Some(0.0),
            None,
        )?;
        add_float(
            &mut p,
            "warp:max_pairwise_log_fc",
            self.warp_max_pairwise_log_fc,
            "Maximum absolute log10 fold change between two compatible signals during compatibility graph construction. Two signals from different maps will not be connected by an edge in the compatibility graph if absolute log fold change exceeds this limit (they might still end up in the same connected component, however). Note: this does not limit fold changes in the linking stage, only during RT alignment, where we try to find high-quality alignment anchor points. Setting this to a value < 0 disables the FC check.",
            true,
            None,
            None,
        )?;
        add_float(
            &mut p,
            "warp:min_rel_cc_size",
            self.warp_min_rel_cc_size,
            "Only connected components containing compatible features from at least max(2, (warp_min_occur * number_of_input_maps)) input maps are considered for computing the warping function",
            true,
            Some(0.0),
            Some(1.0),
        )?;
        add_int(
            &mut p,
            "warp:max_nr_conflicts",
            self.warp_max_nr_conflicts,
            "Allow up to this many conflicts (features from the same map) per connected component to be used for alignment (-1 means allow any number of conflicts)",
            true,
            Some(-1),
            None,
        )?;
        add_float(
            &mut p,
            "link:rt_tol",
            self.link_rt_tol,
            "Width of RT tolerance window (sec)",
            false,
            Some(0.0),
            None,
        )?;
        add_float(
            &mut p,
            "link:mz_tol",
            self.link_mz_tol,
            "m/z tolerance (in ppm or Da)",
            false,
            Some(0.0),
            None,
        )?;
        add_string(
            &mut p,
            "link:charge_merging",
            self.link_charge_merging.as_str(),
            "whether to disallow charge mismatches (Identical), allow to link charge zero (i.e., unknown charge state) with every charge state, or disregard charges (Any).",
            false,
            &["Identical", "With_charge_zero", "Any"],
        )?;
        add_string(
            &mut p,
            "link:adduct_merging",
            self.link_adduct_merging.as_str(),
            "whether to only allow the same adduct for linking (Identical), also allow linking features with adduct-free ones, or disregard adducts (Any).",
            false,
            &["Identical", "With_unknown_adducts", "Any"],
        )?;
        add_string(
            &mut p,
            "mz_unit",
            if self.mz_ppm { "ppm" } else { "Da" },
            "Unit of m/z tolerance",
            false,
            &["ppm", "Da"],
        )?;
        add_int(
            &mut p,
            "nr_partitions",
            self.nr_partitions,
            "Number of partitions in m/z space",
            false,
            Some(1),
            None,
        )?;

        // The FeatureDistance parameters, then the source's overrides.
        p.insert("", &self.distance.to_param()?)?;
        p.set_value(
            "distance_intensity:weight",
            ParamValue::Float(self.distance.intensity_weight),
            "",
            &[],
        )?;
        p.set_value(
            "distance_intensity:log_transform",
            ParamValue::String(
                if self.distance.log_transform {
                    "enabled"
                } else {
                    "disabled"
                }
                .into(),
            ),
            "",
            &[],
        )?;
        p.add_tag("distance_intensity:weight", "advanced")?;
        p.add_tag("distance_intensity:log_transform", "advanced")?;
        for key in [
            "distance_RT:max_difference",
            "distance_MZ:max_difference",
            "distance_MZ:unit",
            "ignore_charge",
            "ignore_adduct",
        ] {
            p.remove(key)?;
        }
        // The source inserts `TransformationModelLowess::getDefaultParameters()`
        // and tags every entry `advanced`.
        let mut lowess = self.lowess.to_param_with(self.lowess_auto_span.as_ref())?;
        let keys: Vec<String> = lowess.iter()?.map(|item| item.key).collect();
        for key in keys {
            lowess.add_tag(&key, "advanced")?;
        }
        p.insert("LOWESS:", &lowess)?;
        p.set_section_description(
            "LOWESS",
            "LOWESS parameters for internal RT transformations (only relevant if 'warp:enabled' is set to 'true')",
        )?;
        Ok(p)
    }

    /// Build from parameters, as the constructor followed by `setParameters`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`](crate::error::Error::InvalidValue) when
    /// a parameter has the wrong type or violates its restriction.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`Self::from_param`], also returning the unknown-parameter warnings.
    ///
    /// # Errors
    ///
    /// As [`Self::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(FEATURE_GROUPING_ALGORITHM_KD_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    /// Read a validated tree whose every default is present.
    pub(crate) fn from_complete_param(p: &Param) -> Result<Self> {
        let charge = match string(p, "link:charge_merging")?.as_str() {
            "Identical" => ChargeMerging::Identical,
            "With_charge_zero" => ChargeMerging::WithChargeZero,
            "Any" => ChargeMerging::Any,
            other => return Err(bad(format!("invalid link:charge_merging '{other}'"))),
        };
        let adduct = match string(p, "link:adduct_merging")?.as_str() {
            "Identical" => AdductMerging::Identical,
            "With_unknown_adducts" => AdductMerging::WithUnknownAdducts,
            "Any" => AdductMerging::Any,
            other => return Err(bad(format!("invalid link:adduct_merging '{other}'"))),
        };
        let mut distance = FeatureDistance::default();
        distance.rt.exponent = float(p, "distance_RT:exponent")?;
        distance.rt.weight = float(p, "distance_RT:weight")?;
        distance.mz.exponent = float(p, "distance_MZ:exponent")?;
        distance.mz.weight = float(p, "distance_MZ:weight")?;
        distance.intensity_exponent = float(p, "distance_intensity:exponent")?;
        distance.intensity_weight = float(p, "distance_intensity:weight")?;
        distance.log_transform = match string(p, "distance_intensity:log_transform")?.as_str() {
            "enabled" => true,
            "disabled" => false,
            other => {
                return Err(bad(format!(
                    "invalid distance_intensity:log_transform '{other}'"
                )));
            }
        };
        let (lowess, lowess_auto_span) = LowessOptions::from_param(&p.copy("LOWESS:", true)?)?;
        Ok(Self {
            warp_enabled: string(p, "warp:enabled")? == "true",
            warp_rt_tol: float(p, "warp:rt_tol")?,
            warp_mz_tol: float(p, "warp:mz_tol")?,
            warp_max_pairwise_log_fc: float(p, "warp:max_pairwise_log_fc")?,
            warp_min_rel_cc_size: float(p, "warp:min_rel_cc_size")?,
            warp_max_nr_conflicts: int(p, "warp:max_nr_conflicts")?,
            link_rt_tol: float(p, "link:rt_tol")?,
            link_mz_tol: float(p, "link:mz_tol")?,
            link_charge_merging: charge,
            link_adduct_merging: adduct,
            mz_ppm: string(p, "mz_unit")? == "ppm",
            nr_partitions: int(p, "nr_partitions")?,
            distance,
            lowess,
            lowess_auto_span,
        })
    }

    /// The distance functor of one run: the source copies the `distance_*`
    /// sections, sets the maximum differences to the link tolerances, the
    /// m/z unit, `ignore_charge` for `link:charge_merging` = `Any`, and
    /// normalises intensities by the largest input intensity.
    fn run_distance(&self, max_intensity: f64) -> Result<FeatureDistance> {
        let mut params = self.distance.to_param()?;
        for key in ["ignore_charge", "ignore_adduct"] {
            params.remove(key)?;
        }
        params.set_value(
            "distance_RT:max_difference",
            ParamValue::Float(self.link_rt_tol),
            "",
            &[],
        )?;
        params.set_value(
            "distance_MZ:max_difference",
            ParamValue::Float(self.link_mz_tol),
            "",
            &[],
        )?;
        params.set_value(
            "distance_MZ:unit",
            ParamValue::String(if self.mz_ppm { "ppm" } else { "Da" }.into()),
            "",
            &[],
        )?;
        if self.link_charge_merging == ChargeMerging::Any {
            params.set_value("ignore_charge", ParamValue::String("true".into()), "", &[])?;
        }
        Ok(FeatureDistance::with_param(max_intensity, false, &params)?.0)
    }

    /// Source template `group_`.
    fn group_impl<M: QtInputMap>(
        &self,
        input_maps: &[M],
        identifications: Vec<InputIdentifications<'_>>,
        out: &mut ConsensusMap,
    ) -> Result<()> {
        if input_maps.len() < 2 {
            return Err(bad("At least two maps must be given!"));
        }
        let mut result = out.clone();
        result.features.clear();

        // All m/z values for partitioning, and the intensity maximum.
        let mut massrange = Vec::new();
        let mut max_intensity = 0.0f64;
        for map in input_maps {
            for feature in map.base_features() {
                // Finite coordinates, so that the sorts and the tree are
                // well defined (the source does not check).
                feature.validate()?;
                massrange.push(feature.mz);
                let intensity = f64::from(feature.intensity);
                if intensity > max_intensity {
                    max_intensity = intensity;
                }
            }
        }
        let distance = self.run_distance(max_intensity)?;
        massrange.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        if massrange.is_empty() {
            postprocess(&identifications, &mut result)?;
            *out = result;
            return Ok(());
        }

        // Partition boundaries at gaps wider than the larger m/z tolerance.
        let nr_partitions = usize::try_from(self.nr_partitions)
            .ok()
            .filter(|&n| n >= 1)
            .ok_or_else(|| bad("nr_partitions must be at least 1"))?;
        let pts_per_partition = massrange.len() / nr_partitions;
        let max_mz_tol = self.link_mz_tol.max(self.warp_mz_tol);
        let mut boundaries = vec![massrange[0]];
        for j in 0..massrange.len() - 1 {
            let diff = if self.mz_ppm {
                max_mz_tol * 1e-6 * massrange[j + 1]
            } else {
                max_mz_tol
            };
            if (massrange[j] - massrange[j + 1]).abs() > diff
                && j >= boundaries.len() * pts_per_partition
            {
                boundaries.push((massrange[j] + massrange[j + 1]) / 2.0);
            }
        }
        boundaries.push(massrange[massrange.len() - 1] + 1.0);

        let partition = |j: usize| -> Vec<M> {
            input_maps
                .iter()
                .map(|map| map.partition(boundaries[j], boundaries[j + 1]))
                .collect()
        };

        // RT transformation models.
        let mut aligner = MapAlignmentAlgorithmKD::new(input_maps.len(), &self.to_param()?)?;
        if self.warp_enabled {
            for j in 0..boundaries.len() - 1 {
                let maps = partition(j);
                let kd_data = KDTreeFeatureMaps::from_maps(&maps)?;
                aligner.add_rt_fit_data(&kd_data)?;
            }
            aligner.fit_lowess()?;
        }

        // Alignment and linking per partition.
        for j in 0..boundaries.len() - 1 {
            let maps = partition(j);
            let mut kd_data = KDTreeFeatureMaps::from_maps(&maps)?;
            if self.warp_enabled {
                aligner.transform(&mut kd_data)?;
            }
            self.run_clustering(&kd_data, &distance, &mut result)?;
        }

        postprocess(&identifications, &mut result)?;
        *out = result;
        Ok(())
    }

    /// Source `runClustering_`.
    fn run_clustering(
        &self,
        kd_data: &KDTreeFeatureMaps<'_>,
        distance: &FeatureDistance,
        out: &mut ConsensusMap,
    ) -> Result<()> {
        let n = kd_data.size();
        let mut potential_clusters = BTreeSet::new();
        let mut cluster_for_idx = vec![ClusterProxyKD::default(); n];
        let mut assigned = vec![false; n];
        let all: BTreeSet<usize> = (0..n).collect();
        self.update_cluster_proxies(
            &mut potential_clusters,
            &mut cluster_for_idx,
            &all,
            &assigned,
            kd_data,
            distance,
        )?;
        while let Some(first) = potential_clusters.first() {
            let i = first.center_index();
            let mut cf_indices = Vec::new();
            self.compute_best_cluster_for_center(i, &mut cf_indices, &assigned, kd_data, distance)?;
            self.add_consensus_feature(&cf_indices, kd_data, out)?;
            for &f in &cf_indices {
                assigned[f] = true;
                potential_clusters.remove(&cluster_for_idx[f]);
            }
            let mut update_these = BTreeSet::new();
            for &f in &cf_indices {
                let mut neighbors = Vec::new();
                kd_data.get_neighborhood(
                    f,
                    &mut neighbors,
                    self.link_rt_tol,
                    self.link_mz_tol,
                    self.mz_ppm,
                    true,
                    -1.0,
                )?;
                for k in neighbors {
                    if !assigned[k] {
                        update_these.insert(k);
                    }
                }
            }
            self.update_cluster_proxies(
                &mut potential_clusters,
                &mut cluster_for_idx,
                &update_these,
                &assigned,
                kd_data,
                distance,
            )?;
        }
        Ok(())
    }

    /// Source `updateClusterProxies_`.
    fn update_cluster_proxies(
        &self,
        potential_clusters: &mut BTreeSet<ClusterProxyKD>,
        cluster_for_idx: &mut [ClusterProxyKD],
        update_these: &BTreeSet<usize>,
        assigned: &[bool],
        kd_data: &KDTreeFeatureMaps<'_>,
        distance: &FeatureDistance,
    ) -> Result<()> {
        for &i in update_these {
            let old_proxy = cluster_for_idx[i];
            let mut unused = Vec::new();
            let new_proxy =
                self.compute_best_cluster_for_center(i, &mut unused, assigned, kd_data, distance)?;
            if new_proxy != old_proxy {
                potential_clusters.remove(&old_proxy);
                cluster_for_idx[i] = new_proxy;
                potential_clusters.insert(new_proxy);
            }
        }
        Ok(())
    }

    /// Source `computeBestClusterForCenter_`.
    fn compute_best_cluster_for_center(
        &self,
        i: usize,
        cf_indices: &mut Vec<usize>,
        assigned: &[bool],
        kd_data: &KDTreeFeatureMaps<'_>,
        distance: &FeatureDistance,
    ) -> Result<ClusterProxyKD> {
        let mut points_for_map_index: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        let mut neighbors = Vec::new();
        kd_data.get_neighborhood(
            i,
            &mut neighbors,
            self.link_rt_tol,
            self.link_mz_tol,
            self.mz_ppm,
            true,
            -1.0,
        )?;
        let charge_i = kd_data.charge(i)?;
        let f_i = kd_data.feature(i)?;
        let adduct_i = f_i.metadata.get(DC_CHARGE_ADDUCTS);
        for &k in &neighbors {
            if assigned[k] {
                continue;
            }
            let charge_k = kd_data.charge(k)?;
            match self.link_charge_merging {
                ChargeMerging::Identical if charge_k != charge_i => continue,
                ChargeMerging::WithChargeZero if charge_k != charge_i && charge_k != 0 => continue,
                _ => {}
            }
            let adduct_k = kd_data.feature(k)?.metadata.get(DC_CHARGE_ADDUCTS);
            match self.link_adduct_merging {
                AdductMerging::Identical => {
                    if adduct_k.is_some() != adduct_i.is_some() {
                        continue;
                    }
                    if let (Some(a), Some(b)) = (adduct_k, adduct_i) {
                        if formula(a)? != formula(b)? {
                            continue;
                        }
                    }
                }
                AdductMerging::WithUnknownAdducts => {
                    if adduct_k.is_some() && adduct_i.is_none() {
                        continue;
                    }
                    if let (Some(a), Some(b)) = (adduct_k, adduct_i) {
                        if a != b && formula(a)? != formula(b)? {
                            continue;
                        }
                    }
                }
                AdductMerging::Any => {}
            }
            points_for_map_index
                .entry(kd_data.map_index(k)?)
                .or_default()
                .push(k);
        }
        // The centre is always part of the cluster, alone for its map.
        points_for_map_index.insert(kd_data.map_index(i)?, vec![i]);

        let mut avg_distance = 0.0;
        for candidates in points_for_map_index.values() {
            let mut min_dist = f64::MAX;
            let mut best_index = None;
            for &c in candidates {
                let dist = distance.distance(kd_data.feature(c)?, f_i)?.1;
                if dist < min_dist {
                    min_dist = dist;
                    best_index = Some(c);
                }
            }
            let best = best_index
                .ok_or_else(|| bad("no candidate has a finite distance to the cluster centre"))?;
            cf_indices.push(best);
            avg_distance += min_dist;
        }
        avg_distance /= cf_indices.len() as f64;
        Ok(ClusterProxyKD::new(cf_indices.len(), avg_distance, i))
    }

    /// Source `addConsensusFeature_`.
    fn add_consensus_feature(
        &self,
        indices: &[usize],
        kd_data: &KDTreeFeatureMaps<'_>,
        out: &mut ConsensusMap,
    ) -> Result<()> {
        let mut cf = ConsensusFeature::new();
        let mut avg_quality = 0.0f32;
        let mut best_quality = 0.0f32;
        let mut best_quality_index = 0usize;
        let mut linked_groups = Vec::new();
        for &i in indices {
            let feature = kd_data.feature(i)?;
            insert_base_feature(&mut cf, kd_data.map_index(i)? as u64, feature)?;
            avg_quality += feature.quality;
            if feature.metadata.contains_key(DC_CHARGE_ADDUCTS)
                && feature.quality > best_quality
                && feature.charge != 0
            {
                best_quality = feature.quality;
                best_quality_index = i;
            }
            if let Some(group) = feature.metadata.get(ADDUCT_GROUP) {
                linked_groups.push(meta_text(group)?);
            }
        }
        let best = kd_data.feature(best_quality_index)?;
        if let Some(adducts) = best.metadata.get(DC_CHARGE_ADDUCTS) {
            let text = to_adduct_string(adducts.as_str()?, best.charge)?;
            cf.metadata
                .insert(IIMN_BEST_ION.into(), MetaValue::from(text));
        }
        if !linked_groups.is_empty() {
            cf.metadata
                .insert(IIMN_LINKED_GROUPS.into(), MetaValue::from(linked_groups));
        }
        avg_quality /= indices.len() as f32;
        cf.quality = avg_quality;
        cf.compute_consensus()?;
        out.features.push(cf);
        Ok(())
    }
}

/// The empirical formula of a `dc_charge_adducts` value (the source's
/// `EmpiricalFormula(DataValue)`, which reads the value as a string).
fn formula(value: &MetaValue) -> Result<EmpiricalFormula> {
    EmpiricalFormula::parse(value.as_str()?)
}

impl FeatureGroupingAlgorithm for FeatureGroupingAlgorithmKD {
    /// Group feature maps (source `group(const std::vector<FeatureMap>&,
    /// ...)`).
    ///
    /// `out`'s features are replaced; its metadata, column headers and
    /// identifications stay, and the input maps' identifications are
    /// appended by [`postprocess`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`](crate::error::Error::InvalidValue) for fewer than two maps (source
    /// `Exception::IllegalArgument`), and the errors of the warp, the
    /// kd-tree and the distance. On error `out` is unchanged.
    fn group(&mut self, maps: &[FeatureMap], out: &mut ConsensusMap) -> Result<()> {
        let identifications = maps.iter().map(feature_map_identifications).collect();
        self.group_impl(maps, identifications, out)
    }

    /// Group consensus maps directly, as for feature maps.
    ///
    /// # Errors
    ///
    /// As [`FeatureGroupingAlgorithm::group`].
    fn group_consensus(&mut self, maps: &[ConsensusMap], out: &mut ConsensusMap) -> Result<()> {
        let identifications = maps.iter().map(consensus_map_identifications).collect();
        self.group_impl(maps, identifications, out)
    }
}
