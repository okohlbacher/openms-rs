// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Feature grouping by Wasserstein optimal transport.
//!
//! Native equivalent of `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmWNet.h`.
//! See `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! Corresponding features across several maps are found by solving a
//! minimum-cost network flow over their `(m/z, RT)` positions. Every pair of
//! maps is aligned independently by
//! [`WNetAligner`](crate::analysis::mapmatching::wnet_aligner::WNetAligner)
//! (2D: scaled m/z and RT), and the pairwise matches are merged into
//! consensus groups by a union-find that never puts two features of one map
//! into a group:
//!
//! 1. each map becomes a spectrum of `(transformed m/z * mz_scale, RT)` with
//!    the feature intensities as weights; with `mz_unit = ppm` the m/z axis is
//!    `ln(m/z)`, so a multiplicative tolerance is an additive one, and with
//!    `decharge_mz` a feature of charge `z > 1` is first projected to its
//!    `[M+H]+` m/z (through the feature's `dc_charge_adduct_mass` meta value
//!    when it has one, pure protons otherwise); `mz_scale` maps the m/z
//!    tolerance onto `max_rt_shift`, so both axes are in seconds;
//! 2. with `normalize_intensities` each spectrum's intensities are divided by
//!    their sum;
//! 3. for every map pair `i < j` with both maps non-empty the aligner's
//!    consensus gives the candidate pairs, each with the average cost per pair
//!    of that map pair;
//! 4. the candidates are sorted by `(cost, lower feature, higher feature)`
//!    over the flat feature numbering and merged greedily;
//! 5. a consensus feature is built for each union-find group, groups in
//!    ascending order of their root, members by map and feature.
//!
//! # Differences from the source
//!
//! - `max_mz_shift_da` or `max_mz_shift_ppm` of exactly 0 makes the source's
//!   m/z scale infinite; every scaled distance is then `NaN`, and the source
//!   converts that to a 64-bit integer, which is undefined behaviour. This
//!   port refuses such an input (the aligner's
//!   [`Error::InvalidValue`](crate::error::Error::InvalidValue) for a
//!   non-finite scaled position).
//! - `max_rt_shift` of 0 (with `trash_cost` 0, which then also resolves to 0)
//!   divides by zero inside the scale factor and is refused for the same
//!   reason.
//! - A feature with a non-numeric `dc_charge_adduct_mass` meta value makes the
//!   source read a union member of the wrong type; here it is an error.
//! - The source logs per-pair candidate counts and costs at INFO level; no
//!   logging is done here.

use std::collections::BTreeMap;

use super::feature_grouping_algorithm::{
    FeatureGroupingAlgorithm, InputIdentifications, consensus_map_identifications,
    feature_map_identifications, insert_base_feature, postprocess,
};
use super::params::{add_float, add_string, bad, flag, flag_text, float, string};
use super::qt_cluster_finder::QtInputMap;
use super::wnet_aligner::{DistanceMetric, Spectrum, WNetAligner};
use crate::concept::constants::PROTON_MASS_U;
use crate::error::Result;
use crate::kernel::features::{BaseFeature, ConsensusFeature, ConsensusMap, FeatureMap};
use crate::math::source_sort::source_sort_by;
use crate::param::{DefaultParamHandler, Param};

/// Parameter-handler name of the source class.
pub const FEATURE_GROUPING_ALGORITHM_WNET_NAME: &str = "FeatureGroupingAlgorithmWNet";

/// Meta value holding the adduct mass `MetaboliteFeatureDeconvolution` stored
/// on a feature.
const DC_CHARGE_ADDUCT_MASS: &str = "dc_charge_adduct_mass";

/// Largest number of features, summed over all input maps, one run accepts.
/// Checked before anything is allocated; the source has no ceiling.
pub const MAX_FEATURES: usize = 1 << 24;

/// Wasserstein feature grouping (source `FeatureGroupingAlgorithmWNet`).
#[derive(Clone, Debug, PartialEq)]
pub struct FeatureGroupingAlgorithmWNet {
    /// Distance metric for comparing feature positions (`distance_metric`,
    /// advanced, `LINF`).
    pub distance_metric: DistanceMetric,
    /// Maximum allowed RT difference in seconds (`max_rt_shift`, 100,
    /// minimum 0).
    pub max_rt_shift: f64,
    /// Whether the m/z tolerance is in ppm (`mz_unit` = `ppm`, the default)
    /// rather than Da.
    pub mz_ppm: bool,
    /// Maximum m/z difference in Da, used when `mz_unit` is `Da`
    /// (`max_mz_shift_da`, 0.3, minimum 0).
    pub max_mz_shift_da: f64,
    /// Maximum m/z difference in ppm, used when `mz_unit` is `ppm`
    /// (`max_mz_shift_ppm`, 10, minimum 0).
    pub max_mz_shift_ppm: f64,
    /// Cost of leaving a feature unmatched, in seconds; 0 derives it from
    /// `max_rt_shift` (`trash_cost`, advanced, 0, minimum 0).
    pub trash_cost: f64,
    /// Normalise each map's intensities to a sum of one before aligning
    /// (`normalize_intensities`, true).
    pub normalize_intensities: bool,
    /// Project observed m/z values to their singly-charged equivalent before
    /// aligning (`decharge_mz`, false).
    pub decharge_mz: bool,
}

impl Default for FeatureGroupingAlgorithmWNet {
    fn default() -> Self {
        Self {
            distance_metric: DistanceMetric::Linf,
            max_rt_shift: 100.0,
            mz_ppm: true,
            max_mz_shift_da: 0.3,
            max_mz_shift_ppm: 10.0,
            trash_cost: 0.0,
            normalize_intensities: true,
            decharge_mz: false,
        }
    }
}

/// The derived quantities of one run (source `AlignParams`).
struct AlignParams {
    metric: DistanceMetric,
    trash_cost: f64,
    mz_scale: f64,
    max_distance: f64,
    use_ppm: bool,
    normalize: bool,
    decharge: bool,
}

/// A candidate pair of flat feature numbers (source `Candidate`).
struct Candidate {
    flat_i: usize,
    flat_j: usize,
    cost: f64,
}

/// Union-find that refuses merges that put two features of one map into one
/// group (source `UnionFind`).
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
    group_maps: Vec<Vec<usize>>,
}

impl UnionFind {
    fn new(feat_map: &[usize]) -> Self {
        Self {
            parent: (0..feat_map.len()).collect(),
            rank: vec![0; feat_map.len()],
            group_maps: feat_map.iter().map(|map| vec![*map]).collect(),
        }
    }

    /// Find with path halving.
    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    /// True if merged (or already in one group), false on a map conflict.
    fn unite(&mut self, a: usize, b: usize) -> bool {
        let mut a = self.find(a);
        let mut b = self.find(b);
        if a == b {
            return true;
        }
        for map in &self.group_maps[b] {
            if self.group_maps[a].contains(map) {
                return false;
            }
        }
        if self.rank[a] < self.rank[b] {
            std::mem::swap(&mut a, &mut b);
        }
        self.parent[b] = a;
        if self.rank[a] == self.rank[b] {
            self.rank[a] += 1;
        }
        let moved = std::mem::take(&mut self.group_maps[b]);
        self.group_maps[a].extend(moved);
        true
    }
}

impl FeatureGroupingAlgorithmWNet {
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
            "distance_metric",
            match self.distance_metric {
                DistanceMetric::L1 => "L1",
                DistanceMetric::L2 => "L2",
                DistanceMetric::Linf => "LINF",
            },
            "Distance metric for comparing feature positions",
            true,
            &["L1", "L2", "LINF"],
        )?;
        add_float(
            &mut p,
            "max_rt_shift",
            self.max_rt_shift,
            "Maximum allowed RT difference (in seconds) for matching features.",
            false,
            Some(0.0),
            None,
        )?;
        add_string(
            &mut p,
            "mz_unit",
            if self.mz_ppm { "ppm" } else { "Da" },
            "Unit for the m/z tolerance. 'Da' uses max_mz_shift_da; 'ppm' uses max_mz_shift_ppm and log-transforms the m/z axis.",
            false,
            &["Da", "ppm"],
        )?;
        add_float(
            &mut p,
            "max_mz_shift_da",
            self.max_mz_shift_da,
            "Maximum allowed m/z difference in Daltons. Used when mz_unit is 'Da'.",
            false,
            Some(0.0),
            None,
        )?;
        add_float(
            &mut p,
            "max_mz_shift_ppm",
            self.max_mz_shift_ppm,
            "Maximum allowed m/z difference in ppm. Used when mz_unit is 'ppm'.",
            false,
            Some(0.0),
            None,
        )?;
        add_float(
            &mut p,
            "trash_cost",
            self.trash_cost,
            "Cost of leaving a feature unmatched (in seconds). Set to 0 to auto-derive from max_rt_shift (a match at the boundary costs the same as no match). Increase to prefer matching over leaving features unmatched.",
            true,
            Some(0.0),
            None,
        )?;
        add_string(
            &mut p,
            "normalize_intensities",
            flag_text(self.normalize_intensities),
            "Normalize feature intensities per map before alignment so that total intensity is equal across maps.",
            false,
            &["true", "false"],
        )?;
        add_string(
            &mut p,
            "decharge_mz",
            flag_text(self.decharge_mz),
            "Convert observed m/z values to their singly-charged (z=1) equivalent before alignment. Features with charge z>1 are projected to [M+H]+ m/z, making features of the same compound but different charge states directly comparable. If a feature carries a 'dc_charge_adduct_mass' MetaValue (set by MetaboliteFeatureDeconvolution), the stored adduct mass is used; otherwise pure H+ adducts are assumed. Features with unknown charge (z=0) are left unchanged.",
            false,
            &["true", "false"],
        )?;
        Ok(p)
    }

    /// Build from parameters, as the constructor followed by `setParameters`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`](crate::error::Error::InvalidValue) when a parameter has the wrong type or
    /// violates its restriction.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`Self::from_param`], also returning the unknown-parameter warnings.
    ///
    /// # Errors
    ///
    /// As [`Self::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(FEATURE_GROUPING_ALGORITHM_WNET_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    /// Read a validated tree whose every default is present.
    pub(crate) fn from_complete_param(p: &Param) -> Result<Self> {
        // The source's parseDistanceMetric: L1, L2, anything else LINF.
        let metric = match string(p, "distance_metric")?.as_str() {
            "L1" => DistanceMetric::L1,
            "L2" => DistanceMetric::L2,
            _ => DistanceMetric::Linf,
        };
        Ok(Self {
            distance_metric: metric,
            max_rt_shift: float(p, "max_rt_shift")?,
            mz_ppm: string(p, "mz_unit")? == "ppm",
            max_mz_shift_da: float(p, "max_mz_shift_da")?,
            max_mz_shift_ppm: float(p, "max_mz_shift_ppm")?,
            trash_cost: float(p, "trash_cost")?,
            normalize_intensities: flag(p, "normalize_intensities")?,
            decharge_mz: flag(p, "decharge_mz")?,
        })
    }

    /// Source `buildAlignParams`.
    fn align_params(&self) -> AlignParams {
        let mut trash_cost = self.trash_cost;
        if trash_cost <= 0.0 {
            trash_cost = self.max_rt_shift;
        }
        // The m/z tolerance in the units of the (possibly log) axis.
        let mz_shift_native = if self.mz_ppm {
            (self.max_mz_shift_ppm * 1e-6).ln_1p()
        } else {
            self.max_mz_shift_da
        };
        // An infinite scale maps any m/z difference to an infinite distance
        // when the tolerance is zero.
        let mz_scale = if mz_shift_native > 0.0 {
            self.max_rt_shift / mz_shift_native
        } else {
            f64::INFINITY
        };
        // After mz_scale both axes are in RT-equivalent seconds, and the
        // maximum distance reflects the norm so that a feature exactly at
        // the per-axis threshold lands on the boundary.
        let max_distance = match self.distance_metric {
            DistanceMetric::L1 => 2.0 * self.max_rt_shift,
            DistanceMetric::L2 => std::f64::consts::SQRT_2 * self.max_rt_shift,
            DistanceMetric::Linf => self.max_rt_shift,
        };
        AlignParams {
            metric: self.distance_metric,
            trash_cost,
            mz_scale,
            max_distance,
            use_ppm: self.mz_ppm,
            normalize: self.normalize_intensities,
            decharge: self.decharge_mz,
        }
    }

    /// Source `dechargeMz`: the `[M+H]+` m/z of a feature of charge above one.
    fn decharge_mz_of(feature: &BaseFeature) -> Result<f64> {
        let z = feature.charge;
        if z <= 1 {
            return Ok(feature.mz);
        }
        let mz = feature.mz;
        let neutral_mass = match feature.metadata.get(DC_CHARGE_ADDUCT_MASS) {
            Some(value) => mz * f64::from(z) - value.as_f64()?,
            None => mz * f64::from(z) - f64::from(z) * PROTON_MASS_U,
        };
        Ok(neutral_mass + PROTON_MASS_U)
    }

    /// Source `mapToSpectrum`.
    fn map_to_spectrum(features: &[&BaseFeature], p: &AlignParams) -> Result<Spectrum> {
        let mut positions = Vec::with_capacity(features.len());
        let mut intensities = Vec::with_capacity(features.len());
        for feature in features {
            let mz = if p.decharge {
                Self::decharge_mz_of(feature)?
            } else {
                feature.mz
            };
            let transformed = if p.use_ppm { mz.ln() } else { mz };
            positions.push([transformed * p.mz_scale, feature.rt]);
            intensities.push(f64::from(feature.intensity));
        }
        Spectrum::new(positions, intensities)
    }

    /// Source template `group_`.
    fn group_impl<M: QtInputMap>(
        &self,
        maps: &[M],
        identifications: Vec<InputIdentifications<'_>>,
        out: &mut ConsensusMap,
    ) -> Result<()> {
        if maps.len() < 2 {
            return Err(bad("At least two maps must be given!"));
        }
        let p = self.align_params();

        let features: Vec<Vec<&BaseFeature>> = maps.iter().map(QtInputMap::base_features).collect();
        let mut total_features = 0usize;
        for map in &features {
            total_features = total_features.saturating_add(map.len());
            if total_features > MAX_FEATURES {
                return Err(bad(format!("more than {MAX_FEATURES} features")));
            }
            for feature in map {
                feature.validate()?;
            }
        }

        // convertAndNormalizeSpectra
        let mut spectra: Vec<Spectrum> = Vec::with_capacity(features.len());
        for map in &features {
            let mut spectrum = Self::map_to_spectrum(map, &p)?;
            if p.normalize && !spectrum.is_empty() {
                spectrum = spectrum.normalized()?;
            }
            spectra.push(spectrum);
        }

        // Flat feature indexing: feature_offsets[map] + feature in map.
        let n = features.len();
        let mut feature_offsets = Vec::with_capacity(n);
        let mut offset = 0usize;
        for map in &features {
            feature_offsets.push(offset);
            offset += map.len();
        }
        let mut feat_map = vec![0usize; total_features];
        for (map_index, map) in features.iter().enumerate() {
            for feature_index in 0..map.len() {
                feat_map[feature_offsets[map_index] + feature_index] = map_index;
            }
        }

        // collectCandidates
        let mut candidates: Vec<Candidate> = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                if spectra[i].is_empty() || spectra[j].is_empty() {
                    continue;
                }
                let mut aligner = WNetAligner::new(
                    &spectra[i],
                    &[&spectra[j]],
                    p.metric,
                    p.max_distance,
                    p.trash_cost,
                )?;
                // Weight 1: solve a single unit-weight consensus.
                aligner.set_point(&[1.0])?;
                let (emp_ids, theo_ids) = aligner.consensus_for_target(0)?;
                // The average cost per pair: globally cheaper matches are
                // merged first whatever the map order.
                let pair_cost = if emp_ids.is_empty() {
                    0.0
                } else {
                    aligner.total_cost()? / emp_ids.len() as f64
                };
                for (e, t) in emp_ids.iter().zip(&theo_ids) {
                    candidates.push(Candidate {
                        flat_i: feature_offsets[i] + e,
                        flat_j: feature_offsets[j] + t,
                        cost: pair_cost,
                    });
                }
            }
        }

        // Ascending cost, ties broken by the ordered endpoint pair: the source
        // compares std::tuple<double, Size, Size> with `<`.
        source_sort_by(&mut candidates, |a, b| {
            let key_a = (a.cost, a.flat_i.min(a.flat_j), a.flat_i.max(a.flat_j));
            let key_b = (b.cost, b.flat_i.min(b.flat_j), b.flat_i.max(b.flat_j));
            if key_a.0 < key_b.0 {
                return true;
            }
            if key_b.0 < key_a.0 {
                return false;
            }
            if key_a.1 != key_b.1 {
                return key_a.1 < key_b.1;
            }
            key_a.2 < key_b.2
        })?;

        let mut union_find = UnionFind::new(&feat_map);
        for candidate in &candidates {
            union_find.unite(candidate.flat_i, candidate.flat_j);
        }

        // The groups by root, members by map and feature.
        let mut groups: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
        for (map_index, map) in features.iter().enumerate() {
            for feature_index in 0..map.len() {
                let root = union_find.find(feature_offsets[map_index] + feature_index);
                groups
                    .entry(root)
                    .or_default()
                    .push((map_index, feature_index));
            }
        }

        let mut result = out.clone();
        for members in groups.values() {
            let mut consensus = ConsensusFeature::new();
            for (map_index, feature_index) in members {
                insert_base_feature(
                    &mut consensus,
                    *map_index as u64,
                    features[*map_index][*feature_index],
                )?;
            }
            consensus.compute_consensus()?;
            result.features.push(consensus);
        }
        postprocess(&identifications, &mut result)?;
        *out = result;
        Ok(())
    }
}

impl FeatureGroupingAlgorithm for FeatureGroupingAlgorithmWNet {
    /// Group feature maps (source `group(const std::vector<FeatureMap>&, ...)`).
    ///
    /// The consensus features are appended to `out`'s, and the input maps'
    /// identifications are appended and everything sorted by
    /// [`postprocess`]. At least two maps are required.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`](crate::error::Error::InvalidValue) for fewer than two maps (source
    /// `Exception::IllegalArgument`), for the aligner's refusals, and for a
    /// feature with a non-finite position or intensity. On error `out` is
    /// unchanged.
    fn group(&mut self, maps: &[FeatureMap], out: &mut ConsensusMap) -> Result<()> {
        let identifications = maps.iter().map(feature_map_identifications).collect();
        self.group_impl(maps, identifications, out)
    }

    /// Group consensus maps directly (source
    /// `group(const std::vector<ConsensusMap>&, ...)`), as for feature maps.
    ///
    /// # Errors
    ///
    /// As [`FeatureGroupingAlgorithm::group`].
    fn group_consensus(&mut self, maps: &[ConsensusMap], out: &mut ConsensusMap) -> Result<()> {
        let identifications = maps.iter().map(consensus_map_identifications).collect();
        self.group_impl(maps, identifications, out)
    }
}
