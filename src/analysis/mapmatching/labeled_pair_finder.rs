// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Matching of isotope-labeled features: features at a fixed distance in RT
//! and m/z within the same map.
//!
//! Native equivalent of `ANALYSIS/MAPMATCHING/LabeledPairFinder.h`. See
//! `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! Exactly one input map is searched for pairs of features of equal charge
//! whose m/z distance is one of `mz_pair_dists` divided by the charge (within
//! `mz_dev`) and whose RT distance from light to heavy lies in
//! `[rt_pair_dist - rt_dev_low, rt_pair_dist + rt_dev_high]`. Every candidate
//! pair is scored by the geometric mean of two bi-Gaussian p-values (m/z and
//! RT), and pairs are then taken greedily, best score first, so that no
//! feature is used twice.
//!
//! With `rt_estimate`, the optimal RT distance and its allowed deviation are
//! first estimated by fitting a Gaussian to the histogram of the RT distances
//! of all charge- and m/z-compatible feature pairs; this works only when the
//! map holds a significant number of pairs.
//!
//! The source's `@todo` notes remain open: labeled MRM experiments (Q1 m/z
//! and charges) are only partly supported through the advanced `mrm` flag,
//! and more than one mass delta per pair (for example from missed cleavages)
//! is not.

use super::base_group_finder::{BaseGroupFinder, check_ids};
use super::feature_grouping_algorithm::insert_base_feature;
use super::params::{add_float, add_string, bad, flag, flag_text, float, string};
use crate::concept::UniqueIdGenerator;
use crate::error::Result;
use crate::kernel::features::{ConsensusFeature, ConsensusMap};
use crate::math::fitters::gauss::{GaussFitResult, GaussFitter};
use crate::math::histogram::Histogram;
use crate::param::{DefaultParamHandler, Param, ParamValue};
use std::cmp::Ordering;
use std::collections::BTreeSet;

/// Parameter-handler name of the source class.
pub const LABELED_PAIR_FINDER_NAME: &str = "LabeledPairFinder";

/// Largest number of candidate comparisons, `|map|^2 * |mz_pair_dists|`, one
/// RT estimation or pair search may evaluate, and the largest number of
/// candidate pairs it may keep. Checked before the work that could exceed it;
/// the source has no ceiling.
pub const MAX_COMPARISONS: usize = 1 << 36;

/// Largest number of candidate pairs (distances for the RT estimation, scored
/// matches for the pair search) one run keeps in memory.
pub const MAX_CANDIDATES: usize = 1 << 26;

/// Iterations the sigma estimate's histogram walk may take in one direction.
/// The walk moves one bin width per step over a histogram of one hundred
/// bins, so it ends after about a hundred steps; the ceiling only stops a walk
/// whose step is lost to rounding, which the source would repeat forever.
const MAX_WALK_STEPS: usize = 1 << 20;

/// The labeled pair finder (source `LabeledPairFinder`).
#[derive(Clone, Debug, PartialEq)]
pub struct LabeledPairFinder {
    /// Estimate the optimal RT distance and deviation by fitting a Gaussian to
    /// the histogram of pair distances (`rt_estimate`, default true). When
    /// false, or when no pair is found, `rt_pair_dist`, `rt_dev_low` and
    /// `rt_dev_high` are used.
    pub rt_estimate: bool,
    /// Optimal pair distance in RT (seconds) from light to heavy feature
    /// (`rt_pair_dist`, default -20).
    pub rt_pair_dist: f64,
    /// Maximum allowed deviation below the optimal RT distance
    /// (`rt_dev_low`, default 15, minimum 0).
    pub rt_dev_low: f64,
    /// Maximum allowed deviation above the optimal RT distance
    /// (`rt_dev_high`, default 15, minimum 0).
    pub rt_dev_high: f64,
    /// Optimal pair distances in m/z (Th) for charge +1, divided by the charge
    /// for higher charges (`mz_pair_dists`, default `[4.0]`).
    pub mz_pair_dists: Vec<f64>,
    /// Maximum allowed deviation from the optimal m/z distance (`mz_dev`,
    /// default 0.05, minimum 0).
    pub mz_dev: f64,
    /// Treat the features as MRM chromatograms, also comparing the precursor
    /// m/z stored in the `MZ` meta value (`mrm`, advanced, default false).
    pub mrm: bool,
}

impl Default for LabeledPairFinder {
    fn default() -> Self {
        Self {
            rt_estimate: true,
            rt_pair_dist: -20.0,
            rt_dev_low: 15.0,
            rt_dev_high: 15.0,
            mz_pair_dists: vec![4.0],
            mz_dev: 0.05,
            mrm: false,
        }
    }
}

/// The RT pair distance and deviations a run used, and how it got them.
#[derive(Clone, Debug, PartialEq)]
pub struct RtPairEstimate {
    /// The RT distance from light to heavy the pair search used.
    pub rt_pair_dist: f64,
    /// The deviation below `rt_pair_dist` the pair search allowed.
    pub rt_dev_low: f64,
    /// The deviation above `rt_pair_dist` the pair search allowed.
    pub rt_dev_high: f64,
    /// Number of candidate RT distances the estimation collected; zero when
    /// `rt_estimate` was off.
    pub estimation_pairs: usize,
    /// Whether the three values above come from the Gaussian fit (true) or
    /// from the parameters (false).
    pub estimated: bool,
}

impl RtPairEstimate {
    /// The warnings the source prints on standard output during estimation:
    /// no pairs found (the manual settings are used), or fewer than fifty
    /// (the estimate is probably not reliable).
    pub fn warnings(&self, rt_estimate: bool) -> Vec<String> {
        if !rt_estimate {
            return Vec::new();
        }
        if self.estimation_pairs == 0 {
            vec![
                "Warning: Could not find pairs for RT distance estimation. The manual settings are used!"
                    .to_string(),
            ]
        } else if self.estimation_pairs < 50 {
            vec![format!(
                "Warning: Found only {} pairs. The estimated shift and std deviation are probably not reliable!",
                self.estimation_pairs
            )]
        } else {
            Vec::new()
        }
    }
}

/// Source `PValue_`: the two-sided p-value at `x` of a bi-Gaussian with mean
/// `m`, left standard deviation `sig1` and right standard deviation `sig2`,
/// `1 - erf(|x - m| / sig / 0.707106781)`, with the source's truncated
/// constant.
fn p_value(x: f64, m: f64, sig1: f64, sig2: f64) -> f64 {
    if m < x {
        1.0 - libm::erf((x - m) / sig2 / 0.707_106_781)
    } else {
        1.0 - libm::erf((m - x) / sig1 / 0.707_106_781)
    }
}

impl LabeledPairFinder {
    /// The source defaults tree.
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
            "rt_estimate",
            flag_text(self.rt_estimate),
            "If 'true' the optimal RT pair distance and deviation are estimated by fitting a gaussian distribution to the histogram of pair distance. Note that this works only datasets with a significant amount of pairs! If 'false' the parameters 'rt_pair_dist', 'rt_dev_low' and 'rt_dev_high' define the optimal distance.",
            false,
            &["true", "false"],
        )?;
        add_float(
            &mut p,
            "rt_pair_dist",
            self.rt_pair_dist,
            "optimal pair distance in RT [sec] from light to heavy feature",
            false,
            None,
            None,
        )?;
        add_float(
            &mut p,
            "rt_dev_low",
            self.rt_dev_low,
            "maximum allowed deviation below optimal retention time distance",
            false,
            Some(0.0),
            None,
        )?;
        add_float(
            &mut p,
            "rt_dev_high",
            self.rt_dev_high,
            "maximum allowed deviation above optimal retention time distance",
            false,
            Some(0.0),
            None,
        )?;
        p.set_value(
            "mz_pair_dists",
            ParamValue::FloatList(self.mz_pair_dists.clone()),
            "optimal pair distances in m/z [Th] for features with charge +1 (adapted to +2, +3, .. by division through charge)",
            &[],
        )?;
        add_float(
            &mut p,
            "mz_dev",
            self.mz_dev,
            "maximum allowed deviation from optimal m/z distance\n",
            false,
            Some(0.0),
            None,
        )?;
        add_string(
            &mut p,
            "mrm",
            flag_text(self.mrm),
            "this option should be used if the features correspond mrm chromatograms (additionally the precursor is taken into account)",
            true,
            &["true", "false"],
        )?;
        Ok(p)
    }

    /// Build from parameters, as the constructor followed by `setParameters`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`](crate::error::Error::InvalidValue) when
    /// a parameter has the wrong type or violates its restriction (source
    /// `Exception::InvalidParameter`).
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
        let mut handler = DefaultParamHandler::new(LABELED_PAIR_FINDER_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    /// Read a validated tree whose every default is present.
    pub(crate) fn from_complete_param(p: &Param) -> Result<Self> {
        Ok(Self {
            // The source compares the string with "true".
            rt_estimate: string(p, "rt_estimate")? == "true",
            rt_pair_dist: float(p, "rt_pair_dist")?,
            rt_dev_low: float(p, "rt_dev_low")?,
            rt_dev_high: float(p, "rt_dev_high")?,
            mz_pair_dists: p.value("mz_pair_dists")?.to_double_vector()?,
            mz_dev: float(p, "mz_dev")?,
            mrm: flag(p, "mrm")?,
        })
    }

    /// Run the pair search, drawing the unique IDs of the candidate pairs
    /// from `generator` (the source's `setUniqueId()` on its singleton, one
    /// ID per candidate pair in search order), and report the RT parameters
    /// used.
    ///
    /// `result`'s features are replaced by the chosen pairs, sorted by m/z;
    /// its column headers, metadata and identifications stay, and the input
    /// map's protein and unassigned peptide identifications are appended.
    /// Each pair holds a handle to the light feature under the header labeled
    /// `light` and one to the heavy feature under the header labeled `heavy`;
    /// the handles' unique IDs are those of the input consensus features. The
    /// pair's peptide identifications are the light and then the heavy
    /// feature's, each stamped with its handle's map index in `map_index`.
    ///
    /// Differences from the source, all on input the source mishandles:
    ///
    /// - the input is sorted by position with a stable sort, where the source
    ///   uses `std::sort`, whose order of features at equal positions is
    ///   unspecified;
    /// - a candidate score that is not a number (from `mz_dev`, `rt_dev_low`
    ///   or `rt_dev_high` of zero with an exact match, where the source
    ///   divides zero by zero) is an error; the source sorts it into an
    ///   unspecified position;
    /// - the estimation's histogram walk is bounded by an iteration ceiling.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`](crate::error::Error::InvalidValue)
    /// (source `Exception::IllegalArgument`) unless there is exactly one input
    /// map, `result` has exactly two column headers with the same file name,
    /// one labeled `light` and one `heavy`, and the input's column headers
    /// have unique IDs; also for a work ceiling exceeded, in `mrm` mode for a
    /// feature without a numeric `MZ` meta value (the source's
    /// `Exception::ConversionError`), and for the histogram and fit errors of
    /// the RT estimation (the source's `Exception::OutOfRange` and
    /// `Exception::UnableToFit`). On error `result` is unchanged.
    pub fn run_with_generator(
        &self,
        input: &[ConsensusMap],
        result: &mut ConsensusMap,
        generator: &mut UniqueIdGenerator,
    ) -> Result<RtPairEstimate> {
        if input.len() != 1 {
            return Err(bad("exactly one input map required"));
        }
        if result.column_headers.len() != 2 {
            return Err(bad("two file descriptions required"));
        }
        let first = result.column_headers.values().next();
        let last = result.column_headers.values().next_back();
        if let (Some(first), Some(last)) = (first, last) {
            if first.filename != last.filename {
                return Err(bad(
                    "the two file descriptions have to contain the same file name",
                ));
            }
        }
        check_ids(input)?;

        // Look up the light and heavy index; a later header with the same
        // label replaces an earlier one, as in the source loop.
        let mut light_index = None;
        let mut heavy_index = None;
        for (index, header) in &result.column_headers {
            if header.label == "heavy" {
                heavy_index = Some(*index);
            } else if header.label == "light" {
                light_index = Some(*index);
            }
        }
        let (Some(light_index), Some(heavy_index)) = (light_index, heavy_index) else {
            return Err(bad("the input maps have to be labeled 'light' and 'heavy'"));
        };

        let map = &input[0];
        // Finite coordinates, so that the orderings below are total.
        for feature in &map.features {
            feature.validate()?;
        }
        // Source `ConstRefVector::sortByPosition`: by RT, then m/z.
        let mut model: Vec<&ConsensusFeature> = map.features.iter().collect();
        model.sort_by(|a, b| less(a.rt, b.rt).then_with(|| less(a.mz, b.mz)));

        let comparisons = model
            .len()
            .checked_mul(model.len())
            .and_then(|n| n.checked_mul(self.mz_pair_dists.len().max(1)));
        if comparisons.is_none_or(|n| n > MAX_COMPARISONS) {
            return Err(bad(format!(
                "labeled pair finding compares more than {MAX_COMPARISONS} feature pairs"
            )));
        }

        let mut estimate = RtPairEstimate {
            rt_pair_dist: self.rt_pair_dist,
            rt_dev_low: self.rt_dev_low,
            rt_dev_high: self.rt_dev_high,
            estimation_pairs: 0,
            estimated: false,
        };
        if self.rt_estimate {
            self.estimate_rt(&model, &mut estimate)?;
        }

        let matches =
            self.candidate_pairs(&model, &estimate, light_index, heavy_index, generator)?;

        // Greedy selection, best quality first (source `sortByQuality(true)`,
        // a stable sort), marking the handles' unique IDs as used.
        let mut order: Vec<usize> = (0..matches.len()).collect();
        order.sort_by(|&a, &b| less(f64::from(matches[b].quality), f64::from(matches[a].quality)));
        let mut used = BTreeSet::new();
        let mut chosen = Vec::new();
        for index in order {
            let candidate = &matches[index];
            let handles = candidate.handles();
            let (Some(first), Some(last)) = (handles.first(), handles.last()) else {
                continue;
            };
            if !used.contains(&first.unique_id) && !used.contains(&last.unique_id) {
                used.insert(first.unique_id);
                used.insert(last.unique_id);
                chosen.push(candidate.clone());
            }
        }

        let mut out = result.clone();
        out.features = chosen;
        out.protein_identifications
            .extend(map.protein_identifications.iter().cloned());
        out.unassigned_peptide_identifications
            .extend(map.unassigned_peptide_identifications.iter().cloned());
        // Source `sortByMZ`, a stable sort.
        out.features.sort_by(|a, b| less(a.mz, b.mz));
        *result = out;
        Ok(estimate)
    }

    /// The source's RT estimation block, updating `estimate` in place.
    fn estimate_rt(
        &self,
        model: &[&ConsensusFeature],
        estimate: &mut RtPairEstimate,
    ) -> Result<()> {
        let mut dists = Vec::new();
        for first in model {
            for second in model {
                for &mz_pair_dist in &self.mz_pair_dists {
                    let expected = first.mz + mz_pair_dist / f64::from(first.charge);
                    if second.charge == first.charge
                        && second.mz >= expected - self.mz_dev
                        && second.mz <= expected + self.mz_dev
                    {
                        if dists.len() == MAX_CANDIDATES {
                            return Err(bad(format!(
                                "RT distance estimation collects more than {MAX_CANDIDATES} pairs"
                            )));
                        }
                        dists.push(second.rt - first.rt);
                    }
                }
            }
        }
        estimate.estimation_pairs = dists.len();
        if dists.is_empty() {
            return Ok(());
        }

        // First estimate of the optimal shift: median of the distances.
        dists.sort_by(f64::total_cmp);
        let median_index = dists.len() / 2;
        // Histogram of the distances, only the maximum number of pairs,
        // centered around the median.
        let max_pairs = model.len() / 2;
        let start_index = median_index.saturating_sub(max_pairs / 2);
        let end_index = (dists.len() - 1).min(median_index + max_pairs / 2);
        let start_value = dists[start_index];
        let end_value = dists[end_index];
        // Ensure that we have 100 bins.
        let bin_step = (end_value - start_value).abs() / 99.999;
        let mut hist = Histogram::with_bounds(start_value, end_value, bin_step)?;
        for &dist in &dists[start_index..=end_index] {
            hist.inc(dist)?;
        }
        // Median of the bins (uniform background distribution).
        let mut bins = hist.bins().to_vec();
        bins.sort_by(f64::total_cmp);
        let bin_median = bins[bins.len() / 2];
        // Scale: the maximum of the histogram above the median.
        let max_value = hist
            .max_value()
            .ok_or_else(|| bad("RT distance histogram has no bins"))?;
        let mut initial = GaussFitResult::new(max_value - bin_median, dists[median_index], -1.0);
        // Overwrite the center with the position of the highest bin.
        for i in 0..hist.len() {
            if hist.bin(i)? == max_value {
                initial.x0 = hist.center_of_bin(i)?;
                break;
            }
        }
        // Sigma: the first position where the count drops to the median.
        let mut pos = initial.x0;
        let mut steps = 0;
        while pos > start_value && hist.bin_value(pos)? > bin_median {
            pos -= bin_step;
            steps += 1;
            if steps > MAX_WALK_STEPS {
                return Err(bad("RT distance histogram walk does not terminate"));
            }
        }
        let sigma_low = initial.x0 - pos;
        pos = initial.x0;
        steps = 0;
        while pos < end_value && hist.bin_value(pos)? > bin_median {
            pos += bin_step;
            steps += 1;
            if steps > MAX_WALK_STEPS {
                return Err(bad("RT distance histogram walk does not terminate"));
            }
        }
        let sigma_high = pos - initial.x0;
        initial.sigma = (sigma_high + sigma_low) / 6.0;

        let mut points = Vec::with_capacity(hist.len());
        for i in 0..hist.len() {
            points.push((hist.center_of_bin(i)?, hist.bin(i)?));
        }
        let mut fitter = GaussFitter::new();
        fitter.set_initial_parameters(initial);
        let fitted = fitter.fit(&points)?;
        estimate.rt_pair_dist = fitted.x0;
        estimate.rt_dev_low = fitted.sigma.abs() * 3.0;
        estimate.rt_dev_high = fitted.sigma.abs() * 3.0;
        estimate.estimated = true;
        Ok(())
    }

    /// Every candidate pair, scored, in the source's search order.
    fn candidate_pairs(
        &self,
        model: &[&ConsensusFeature],
        estimate: &RtPairEstimate,
        light_index: u64,
        heavy_index: u64,
        generator: &mut UniqueIdGenerator,
    ) -> Result<Vec<ConsensusFeature>> {
        let mut matches = Vec::new();
        for first in model {
            for &mz_pair_dist in &self.mz_pair_dists {
                let low = first.rt + estimate.rt_pair_dist - estimate.rt_dev_low;
                let high = first.rt + estimate.rt_pair_dist + estimate.rt_dev_high;
                // Source `lower_bound` with `ConsensusFeature::RTLess`.
                let begin = model.partition_point(|f| f.rt < low);
                for second in model[begin..].iter().take_while(|f| f.rt <= high) {
                    let charge = f64::from(first.charge);
                    let accepted = if self.mrm {
                        // Precursor and fragment m/z differences; charge stays.
                        let mut prec_mz_diff = (meta_mz(second)? - meta_mz(first)?).abs();
                        prec_mz_diff = if first.charge != 0 {
                            (prec_mz_diff - mz_pair_dist / charge).abs()
                        } else {
                            (prec_mz_diff - mz_pair_dist).abs()
                        };
                        let frag_mz_diff = (first.mz - second.mz).abs();
                        second.charge == first.charge
                            && prec_mz_diff < self.mz_dev
                            && (frag_mz_diff < self.mz_dev
                                || (frag_mz_diff - mz_pair_dist).abs() < self.mz_dev)
                    } else {
                        let expected = first.mz + mz_pair_dist / charge;
                        second.charge == first.charge
                            && second.mz >= expected - self.mz_dev
                            && second.mz <= expected + self.mz_dev
                    };
                    if !accepted {
                        continue;
                    }
                    let score = (p_value(
                        second.mz - first.mz,
                        mz_pair_dist / charge,
                        self.mz_dev,
                        self.mz_dev,
                    ) * p_value(
                        second.rt - first.rt,
                        estimate.rt_pair_dist,
                        estimate.rt_dev_low,
                        estimate.rt_dev_high,
                    ))
                    .sqrt();
                    if score.is_nan() {
                        return Err(bad(
                            "labeled pair score is not a number (a zero deviation with an exact match)",
                        ));
                    }
                    if matches.len() == MAX_CANDIDATES {
                        return Err(bad(format!(
                            "labeled pair finding keeps more than {MAX_CANDIDATES} candidate pairs"
                        )));
                    }
                    matches.push(pair(
                        first,
                        second,
                        light_index,
                        heavy_index,
                        score,
                        generator,
                    )?);
                }
            }
        }
        Ok(matches)
    }
}

/// The ordering `operator<` induces on finite values: `-0.0` and `0.0` are
/// equal, as they are to the source's comparators.
fn less(a: f64, b: f64) -> Ordering {
    a.partial_cmp(&b).unwrap_or(Ordering::Equal)
}

/// The `MZ` meta value of an MRM feature, as the source's `(double)
/// getMetaValue("MZ")`.
fn meta_mz(feature: &ConsensusFeature) -> Result<f64> {
    feature
        .metadata
        .get("MZ")
        .ok_or_else(|| bad("MRM feature has no 'MZ' meta value"))?
        .as_f64()
}

/// Source: a new consensus feature with a fresh unique ID, the light feature
/// inserted under `light_index` (meta values then cleared) and the heavy
/// feature under `heavy_index`, quality `score`, the light feature's charge,
/// and the monoisotopic consensus position.
fn pair(
    light: &ConsensusFeature,
    heavy: &ConsensusFeature,
    light_index: u64,
    heavy_index: u64,
    score: f64,
    generator: &mut UniqueIdGenerator,
) -> Result<ConsensusFeature> {
    let mut consensus = ConsensusFeature::new();
    consensus.unique_id = generator.get_unique_id();
    insert_base_feature(&mut consensus, light_index, &light.base)?;
    consensus.metadata.clear();
    insert_base_feature(&mut consensus, heavy_index, &heavy.base)?;
    consensus.quality = score as f32;
    consensus.charge = light.charge;
    consensus.compute_monoisotopic_consensus()?;
    Ok(consensus)
}

impl BaseGroupFinder for LabeledPairFinder {
    /// [`LabeledPairFinder::run_with_generator`] with a clock-seeded
    /// generator, discarding the RT report.
    fn run(&mut self, input: &[ConsensusMap], result: &mut ConsensusMap) -> Result<()> {
        let mut generator = UniqueIdGenerator::new();
        self.run_with_generator(input, result, &mut generator)
            .map(|_| ())
    }
}
