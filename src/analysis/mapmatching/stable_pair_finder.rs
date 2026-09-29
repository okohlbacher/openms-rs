// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Pair finding across two consensus maps by mutual nearest neighbours with a
//! stability criterion.
//!
//! Native equivalent of `ANALYSIS/MAPMATCHING/StablePairFinder.h`. See
//! `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! The corresponding consensus features must be aligned, but may deviate
//! slightly in position. The distance is [`FeatureDistance`](crate::analysis::mapmatching::feature_distance::FeatureDistance).
//!
//! # Additional criteria for pairing
//!
//! With `use_identifications`, the peptide identifications annotated to two
//! features must be compatible (no annotation, or the same best-hit
//! sequences) for them to pair.
//!
//! Stability: the distance to the nearest neighbour must be smaller than the
//! distance to the second-nearest neighbour by the factor
//! `second_nearest_gap`. When that factor exceeds one, *lowering*
//! `max_difference` may yield *more* pairs, because it widens the gap between
//! nearest and second-nearest neighbour.
//!
//! # Quality
//!
//! The quality of pairing elements `i` and `j` at distance `d` with
//! second-nearest distances `d2_i`, `d2_j` and gap `g` is
//!
//! ```text
//! q_ij = (1 - d) * (1 - g*d/d2_i) * (1 - g*d/d2_j)
//! ```
//!
//! which lies in `[0, 1]` because pairing requires `d <= 1`, `g*d <= d2_i` and
//! `g*d <= d2_j`. The merged consensus feature gets the average of `q_ij` and
//! the existing qualities weighted by the elements' sizes `s_i`, `s_j`:
//!
//! ```text
//! q = (q_ij + (s_i - 1)*q_i + (s_j - 1)*q_j) / (s_i + s_j - 1)
//! ```
//!
//! so stepwise linking (as `FeatureGroupingAlgorithmUnlabeled` does) ends with
//! qualities that incorporate every pairing. Missing elements are not
//! penalised.

use super::base_group_finder::{BaseGroupFinder, check_ids};
use super::feature_distance::FeatureDistance;
use super::params::{add_float, add_string, flag, flag_text, float};
use crate::error::{Error, Result};
use crate::identification::PeptideIdentification;
use crate::kernel::features::{ConsensusFeature, ConsensusMap};
use crate::param::{DefaultParamHandler, Param};
use std::collections::BTreeSet;

/// Parameter-handler name of the source class.
pub const STABLE_PAIR_FINDER_NAME: &str = "StablePairFinder";

/// Largest number of feature pairs, `|map 0| * |map 1|`, whose distance one
/// run evaluates. Checked before any work; the source compares every pair
/// without a ceiling.
pub const MAX_COMPARISONS: usize = 1 << 36;

/// Stable pair finder (source `StablePairFinder`).
#[derive(Clone, Debug, PartialEq)]
pub struct StablePairFinder {
    /// The factor by which the distance to the second-nearest neighbours must
    /// exceed the distance of the pair (`second_nearest_gap`, default 2,
    /// minimum 1).
    pub second_nearest_gap: f64,
    /// Never link features annotated with different peptides
    /// (`use_identifications`, default false).
    pub use_identifications: bool,
    /// The distance parameters (the `distance_*`, `ignore_charge` and
    /// `ignore_adduct` entries). Its `max_intensity` and `force_constraints`
    /// are replaced at every run, as the source constructs a fresh functor.
    pub distance: FeatureDistance,
}

impl Default for StablePairFinder {
    fn default() -> Self {
        Self {
            second_nearest_gap: 2.0,
            use_identifications: false,
            distance: FeatureDistance::default(),
        }
    }
}

impl StablePairFinder {
    /// The source defaults tree: `second_nearest_gap`, `use_identifications`
    /// and the [`FeatureDistance::defaults`] at the root.
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
        add_float(
            &mut p,
            "second_nearest_gap",
            self.second_nearest_gap,
            "Only link features whose distance to the second nearest neighbors (for both sides) is larger by 'second_nearest_gap' than the distance between the matched pair itself.",
            false,
            Some(1.0),
            None,
        )?;
        add_string(
            &mut p,
            "use_identifications",
            flag_text(self.use_identifications),
            "Never link features that are annotated with different peptides (features without ID's always match; only the best hit per peptide identification is considered).",
            false,
            &["true", "false"],
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

    /// As [`Self::from_param`], also returning the unknown-parameter warnings
    /// the source logs.
    ///
    /// # Errors
    ///
    /// As [`Self::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(STABLE_PAIR_FINDER_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    /// Read a validated tree whose every default is present.
    pub(crate) fn from_complete_param(p: &Param) -> Result<Self> {
        Ok(Self {
            second_nearest_gap: float(p, "second_nearest_gap")?,
            use_identifications: flag(p, "use_identifications")?,
            distance: FeatureDistance::from_complete_param(1.0, false, p)?,
        })
    }

    /// The maximum intensity the source normalises by: the larger of the two
    /// maps' intensity maxima over consensus features **and** handles (the
    /// `ConsensusMap::updateRanges` range), skipping an empty map, and
    /// `f64::MIN` (source `numeric_limits<double>::lowest()`) when both are
    /// empty.
    ///
    /// The source reads the maps' *stored* ranges, which are only current
    /// when the caller ran `updateRanges()`; this computes them. A caller
    /// reproducing a stale range passes its own value to
    /// [`Self::run_with_max_intensity`], as
    /// [`FeatureGroupingAlgorithmUnlabeled`](super::feature_grouping_algorithm_unlabeled::FeatureGroupingAlgorithmUnlabeled)
    /// does.
    ///
    /// # Errors
    ///
    /// Returns the [`ConsensusMap::ranges`] validation error.
    pub fn max_intensity(maps: &[&ConsensusMap]) -> Result<f64> {
        let mut max_intensity = f64::MIN;
        for map in maps {
            if let Some(range) = map.ranges()?.intensity {
                max_intensity = max_intensity.max(range.max);
            }
        }
        Ok(max_intensity)
    }

    /// As [`BaseGroupFinder::run`], with the normalising maximum intensity
    /// supplied by the caller rather than computed by
    /// [`Self::max_intensity`].
    ///
    /// # Errors
    ///
    /// As [`BaseGroupFinder::run`].
    pub fn run_with_max_intensity(
        &self,
        input: &[ConsensusMap],
        result: &mut ConsensusMap,
        max_intensity: f64,
    ) -> Result<()> {
        if input.len() != 2 {
            return Err(Error::InvalidValue(
                "exactly two input maps required".into(),
            ));
        }
        check_ids(input)?;
        let (map0, map1) = (&input[0].features, &input[1].features);
        if map0.len().saturating_mul(map1.len()) > MAX_COMPARISONS {
            return Err(Error::InvalidValue(format!(
                "{} x {} feature pairs exceed the {MAX_COMPARISONS}-comparison ceiling",
                map0.len(),
                map1.len()
            )));
        }

        let mut distance = self.distance.clone();
        distance.max_intensity = max_intensity;
        distance.force_constraints = false;

        let init = (FeatureDistance::INFINITY, FeatureDistance::INFINITY);
        let mut nn_index_0: Vec<Option<usize>> = vec![None; map0.len()];
        let mut nn_distance_0 = vec![init; map0.len()];
        let mut nn_index_1: Vec<Option<usize>> = vec![None; map1.len()];
        let mut nn_distance_1 = vec![init; map1.len()];

        // Best hits are compared as sets of sequence strings; computing them
        // once per feature leaves the comparison order unchanged.
        let (best0, best1) = if self.use_identifications {
            (best_hit_sets(map0)?, best_hit_sets(map1)?)
        } else {
            (Vec::new(), Vec::new())
        };

        for (fi0, feat0) in map0.iter().enumerate() {
            for (fi1, feat1) in map1.iter().enumerate() {
                if self.use_identifications && !compatible(&best0[fi0], &best1[fi1]) {
                    continue;
                }
                let (valid, d) = distance.distance(&feat0.base, &feat1.base)?;
                // Only best matches must satisfy the constraints; a
                // second-best distance may become smaller than the best one.
                update(&mut nn_distance_0[fi0], &mut nn_index_0[fi0], d, valid, fi1);
                update(&mut nn_distance_1[fi1], &mut nn_index_1[fi1], d, valid, fi0);
            }
        }

        let gap = self.second_nearest_gap;
        let mut features = Vec::new();
        let mut singleton0 = vec![true; map0.len()];
        let mut singleton1 = vec![true; map1.len()];
        for fi0 in 0..map0.len() {
            let (d, d2_0) = nn_distance_0[fi0];
            let Some(fi1) = nn_index_0[fi0] else {
                continue;
            };
            if !(d < FeatureDistance::INFINITY && d * gap <= d2_0) {
                continue;
            }
            let (d1, d2_1) = nn_distance_1[fi1];
            // The source's `<=`, false for a NaN.
            #[allow(clippy::neg_cmp_op_on_partial_ord)]
            let rejected = nn_index_1[fi1] != Some(fi0) || !(d1 * gap <= d2_1);
            if rejected {
                continue;
            }
            let mut f = ConsensusFeature::new();
            insert_consensus(&mut f, &map0[fi0])?;
            insert_consensus(&mut f, &map1[fi1])?;
            f.compute_consensus()?;
            let mut quality = 1.0 - d;
            let quality0 = 1.0 - d * gap / d2_0;
            let quality1 = 1.0 - d1 * gap / d2_1;
            quality = quality * quality0 * quality1;
            let size0 = map0[fi0].len().max(1);
            let size1 = map1[fi1].len().max(1);
            // `getQuality() * (size - 1)` is `float * size_t`, evaluated in
            // `float` before it is widened.
            let quality0 = f64::from(map0[fi0].quality * (size0 - 1) as f32);
            let quality1 = f64::from(map1[fi1].quality * (size1 - 1) as f32);
            f.quality = ((quality + quality0 + quality1) / (size0 + size1 - 1) as f64) as f32;
            features.push(f);
            singleton0[fi0] = false;
            singleton1[fi1] = false;
        }

        for (map, singleton) in [(map0, &singleton0), (map1, &singleton1)] {
            for (feature, &single) in map.iter().zip(singleton.iter()) {
                if single {
                    let mut copy = feature.clone();
                    if copy.len() < 2 {
                        copy.quality = 0.0;
                    }
                    features.push(copy);
                }
            }
        }

        // Canonical ordering (source `sortByMZ`, a stable sort).
        features.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal));
        result.features = features;
        Ok(())
    }

    /// Whether two consensus features have compatible peptide identifications
    /// (source `compatibleIDs_`): either has none, or the sets of best-hit
    /// sequences of their identifications (identifications without hits
    /// skipped) are equal.
    ///
    /// # Errors
    ///
    /// Returns the [`PeptideIdentification::best_hit`] validation error.
    pub fn compatible_ids(feat1: &ConsensusFeature, feat2: &ConsensusFeature) -> Result<bool> {
        Ok(compatible(&best_hit_set(feat1)?, &best_hit_set(feat2)?))
    }

    /// The sequence of the best-scoring hit (source `getBestHitSequence_`):
    /// highest score when higher is better, lowest otherwise, the first of
    /// several equal ones. `None` without hits, where the source's
    /// `min_element` would return `end()`.
    ///
    /// # Errors
    ///
    /// Returns the [`PeptideIdentification::best_hit`] validation error.
    pub fn best_hit_sequence(identification: &PeptideIdentification) -> Result<Option<String>> {
        Ok(identification
            .best_hit()?
            .map(|hit| hit.sequence.to_string()))
    }
}

impl BaseGroupFinder for StablePairFinder {
    /// Pair the consensus features of exactly two maps (source `run`).
    ///
    /// For every element, the nearest valid and the second-nearest neighbour
    /// in the other map are found over all pairs. Two elements that are each
    /// other's nearest neighbour, at finite distance, with both second-nearest
    /// distances at least `second_nearest_gap` times larger, become one
    /// consensus feature holding both elements' handles and peptide
    /// identifications, with [`ConsensusFeature::compute_consensus`] applied and
    /// the quality of the module documentation. All other elements are copied
    /// unchanged, map 0 first, a copy with fewer than two handles getting
    /// quality 0. The result is sorted by m/z, stably.
    ///
    /// Only the features of `result` are replaced; its metadata, column
    /// headers and identifications are kept, as the source's `clear(false)`
    /// keeps them. Protein and unassigned peptide identifications are added by
    /// the feature grouping algorithm, not here.
    ///
    /// The normalising intensity is [`StablePairFinder::max_intensity`] of the
    /// two maps.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] unless there are exactly two input
    /// maps, when their column headers share a map identifier (source
    /// `Exception::IllegalArgument`), when `|map 0| * |map 1|` exceeds
    /// [`MAX_COMPARISONS`], and the [`FeatureDistance::distance`] and
    /// consensus-computation errors. On error `result` is unchanged.
    fn run(&self, input: &[ConsensusMap], result: &mut ConsensusMap) -> Result<()> {
        let refs: Vec<&ConsensusMap> = input.iter().take(2).collect();
        let max_intensity = Self::max_intensity(&refs)?;
        self.run_with_max_intensity(input, result, max_intensity)
    }
}

/// Source nearest/second-nearest update for one element.
fn update(
    nn: &mut (f64, f64),
    index: &mut Option<usize>,
    distance: f64,
    valid: bool,
    other: usize,
) {
    if distance < nn.1 {
        if valid && distance < nn.0 {
            nn.1 = nn.0;
            nn.0 = distance;
            *index = Some(other);
        } else {
            nn.1 = distance;
        }
    }
}

/// Source `ConsensusFeature::insert(const ConsensusFeature&)`: the union of
/// the handles, keeping existing ones, and the peptide identifications
/// appended.
fn insert_consensus(target: &mut ConsensusFeature, other: &ConsensusFeature) -> Result<()> {
    target.merge(other)?;
    target
        .peptide_identifications
        .extend(other.peptide_identifications.iter().cloned());
    Ok(())
}

/// Best-hit sequence sets, or `None` for a feature without identifications,
/// which is compatible with everything.
type BestHits = Option<BTreeSet<String>>;

fn best_hit_set(feature: &ConsensusFeature) -> Result<BestHits> {
    if feature.peptide_identifications.is_empty() {
        return Ok(None);
    }
    let mut set = BTreeSet::new();
    for identification in &feature.peptide_identifications {
        if let Some(sequence) = StablePairFinder::best_hit_sequence(identification)? {
            set.insert(sequence);
        }
    }
    Ok(Some(set))
}

fn best_hit_sets(features: &[ConsensusFeature]) -> Result<Vec<BestHits>> {
    features.iter().map(best_hit_set).collect()
}

fn compatible(a: &BestHits, b: &BestHits) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => a == b,
        _ => true,
    }
}
