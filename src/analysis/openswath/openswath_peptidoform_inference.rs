// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OpenSWATH inference of peptidoforms (IPF).
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathPeptidoformInference.h` and its
//! translation unit. The class is file-format agnostic and works only on
//! compact typed rows. It follows PyProphet's IPF: a precursor layer combines
//! peak-group evidence with optional MS1/MS2 precursor evidence into a
//! precursor peak-group PEP, which is then the prior of a transition layer
//! whose hypotheses are the candidate peptidoforms plus H0 (no peptidoform
//! explains the peak group). See Rosenberger, G., Liu, Y., Rost, H. et al.
//! Inference and quantification of peptidoforms in large sample cohorts by
//! SWATH-MS. Nat Biotechnol 35, 781-788 (2017),
//! <https://doi.org/10.1038/nbt.3908>.
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_INFERENCE_SUPPORT.md).
//!
//! # Evaluation order
//!
//! [`OpenSwathPeptidoformInference::apply_bm`] multiplies the evidence of
//! each feature/hypothesis group in row order after sorting by feature and
//! hypothesis. The source sorts with the unstable `std::sort`, so the order of
//! rows *within* a group (and therefore the last bit of a product of three or
//! more factors) is unspecified there; the port sorts stably and multiplies in
//! input order. Across-run propagation likewise visits alignment groups in
//! ascending id where the source walks an `std::unordered_map`. A product of
//! at most two factors is commutative in IEEE arithmetic and unaffected.
//!
//! The source contains no `#pragma omp`; the port is serial.

use crate::analysis::openswath::openswath_inference_config::PeptidoformInferenceConfig;
use crate::analysis::openswath::openswath_inference_data::{
    IPFAlignmentRow, IPFPrecursorProbabilityRow, IPFPrecursorRow, IPFResultRow, IPFTransitionRow,
};
use crate::math::multiple_testing::{MAX_ITEMS, compute_model_fdr};
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};

/// Largest number of transition rows across-run propagation may produce
/// (native; the source has no ceiling). Propagation copies every confident row
/// of an alignment group to every feature of that group, so its output grows
/// with the product of the two; the size is computed before anything is
/// allocated. Equal to the multiple-testing ceiling
/// [`MAX_ITEMS`], which the propagated rows would meet next.
pub const MAX_PROPAGATED_ROWS: usize = MAX_ITEMS;

/// OpenSWATH inference of peptidoforms (IPF). Stateless.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OpenSwathPeptidoformInference;

/// Compact Bayesian model row used by the public helper methods.
///
/// The meaning of `hypothesis` depends on the inference layer: precursor
/// inference uses `1` for the true precursor state and `0` for the false
/// precursor state; transition / peptidoform inference uses peptide ids as
/// hypotheses and reserves `-1` for the H0/null hypothesis.
#[derive(Clone, Debug, PartialEq)]
pub struct BayesianModelRow {
    /// Feature (peak group) id; `-1` by default.
    pub feature_id: i64,
    /// Hypothesis, see the type notes; `-1` by default.
    pub hypothesis: i64,
    /// Prior probability of the hypothesis; `0.0` by default.
    pub prior: f64,
    /// Likelihood of one piece of evidence under the hypothesis; `0.0` by
    /// default.
    pub evidence: f64,
}

impl Default for BayesianModelRow {
    fn default() -> Self {
        Self {
            feature_id: -1,
            hypothesis: -1,
            prior: 0.0,
            evidence: 0.0,
        }
    }
}

/// Posterior probability row returned by
/// [`OpenSwathPeptidoformInference::apply_bm`]. `hypothesis` follows the
/// conventions of [`BayesianModelRow`].
#[derive(Clone, Debug, PartialEq)]
pub struct PosteriorRow {
    /// Feature (peak group) id; `-1` by default.
    pub feature_id: i64,
    /// Hypothesis; `-1` by default.
    pub hypothesis: i64,
    /// Posterior probability of the hypothesis; `0.0` by default.
    pub posterior: f64,
}

impl Default for PosteriorRow {
    fn default() -> Self {
        Self {
            feature_id: -1,
            hypothesis: -1,
            posterior: 0.0,
        }
    }
}

/// The source's anonymous `TransitionInferenceRow`: a transition row joined to
/// its precursor PEP and alignment group (`-1` for none).
#[derive(Clone, Debug)]
struct TransitionInferenceRow {
    feature_id: i64,
    transition_id: i64,
    hypothesis: i64,
    pep: f64,
    precursor_peakgroup_pep: f64,
    bmask: i32,
    num_peptidoforms: i32,
    alignment_group_id: i64,
}

impl TransitionInferenceRow {
    /// The tuple the source sorts and collapses propagated rows by.
    fn key(&self) -> (i64, i64, i64, i32, i32, i64) {
        (
            self.feature_id,
            self.transition_id,
            self.hypothesis,
            self.bmask,
            self.num_peptidoforms,
            self.alignment_group_id,
        )
    }
}

impl OpenSwathPeptidoformInference {
    /// Compute model FDR / q-values from posterior error probabilities.
    ///
    /// This reuses [`compute_model_fdr`] and follows the PyProphet IPF
    /// convention of computing the model FDR from sorted PEP values with
    /// maximum-rank tie handling. The result is in the original input order;
    /// a single `NaN` input makes every output `NaN`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when `pep_values` exceeds
    /// [`MAX_ITEMS`] (native).
    pub fn compute_model_fdr(pep_values: &[f64]) -> Result<Vec<f64>> {
        compute_model_fdr(pep_values)
    }

    /// Build precursor-layer Bayesian model rows.
    ///
    /// For each input row with MS1 precursor evidence, and again for each with
    /// MS2 precursor evidence, this emits a true (`1`) row with prior
    /// `1 - ms2_peakgroup_pep` and evidence `1 - precursor_pep`, then a false
    /// (`0`) row with prior `ms2_peakgroup_pep` and evidence `precursor_pep`.
    /// Afterwards, for each feature (once, in first-occurrence order) that
    /// carried no precursor evidence on any row, it emits the synthetic
    /// true/false pair PyProphet uses: evidence `0` for the true state and `1`
    /// for the false state.
    pub fn prepare_precursor_bm(rows: &[IPFPrecursorRow]) -> Vec<BayesianModelRow> {
        let mut bm_rows = Vec::with_capacity(rows.len().saturating_mul(4));
        let mut with_evidence: HashSet<i64> = HashSet::with_capacity(rows.len());
        for row in rows {
            for pep in [row.ms1_precursor_pep, row.ms2_precursor_pep]
                .into_iter()
                .flatten()
            {
                with_evidence.insert(row.feature_id);
                bm_rows.push(BayesianModelRow {
                    feature_id: row.feature_id,
                    hypothesis: 1,
                    prior: 1.0 - row.ms2_peakgroup_pep,
                    evidence: 1.0 - pep,
                });
                bm_rows.push(BayesianModelRow {
                    feature_id: row.feature_id,
                    hypothesis: 0,
                    prior: row.ms2_peakgroup_pep,
                    evidence: pep,
                });
            }
        }
        let mut missing_added: HashSet<i64> = HashSet::with_capacity(rows.len());
        for row in rows {
            if !with_evidence.contains(&row.feature_id) && missing_added.insert(row.feature_id) {
                bm_rows.push(BayesianModelRow {
                    feature_id: row.feature_id,
                    hypothesis: 1,
                    prior: 1.0 - row.ms2_peakgroup_pep,
                    evidence: 0.0,
                });
                bm_rows.push(BayesianModelRow {
                    feature_id: row.feature_id,
                    hypothesis: 0,
                    prior: row.ms2_peakgroup_pep,
                    evidence: 1.0,
                });
            }
        }
        bm_rows
    }

    /// Apply a compact Bayesian model and return posterior probabilities.
    ///
    /// Rows are grouped by feature and hypothesis. Within each group the
    /// evidence values are multiplied and the smallest prior in the group is
    /// retained, matching the PyProphet implementation being ported; the
    /// product times that prior is the group's likelihood-prior. Each
    /// feature's likelihood-priors are then normalised by their sum.
    ///
    /// # Returns
    ///
    /// One posterior per feature/hypothesis pair, ordered by feature id and
    /// then hypothesis. A feature whose likelihood-prior sum is zero or not
    /// finite gets a posterior of `0` for every hypothesis instead of a
    /// division by zero, as in the source. An empty input gives an empty
    /// result. See the module notes on the multiplication order.
    pub fn apply_bm(rows: &[BayesianModelRow]) -> Vec<PosteriorRow> {
        if rows.is_empty() {
            return Vec::new();
        }
        let mut sorted: Vec<&BayesianModelRow> = rows.iter().collect();
        sorted.sort_by_key(|r| (r.feature_id, r.hypothesis));

        let mut likelihoods: Vec<(i64, i64, f64)> = Vec::with_capacity(sorted.len());
        let mut i = 0;
        while i < sorted.len() {
            let feature_id = sorted[i].feature_id;
            let hypothesis = sorted[i].hypothesis;
            let mut evidence_product = 1.0;
            let mut prior = sorted[i].prior;
            let mut j = i;
            while j < sorted.len()
                && sorted[j].feature_id == feature_id
                && sorted[j].hypothesis == hypothesis
            {
                evidence_product *= sorted[j].evidence;
                prior = f64_min(prior, sorted[j].prior);
                j += 1;
            }
            likelihoods.push((feature_id, hypothesis, evidence_product * prior));
            i = j;
        }

        let mut posteriors = Vec::with_capacity(likelihoods.len());
        let mut i = 0;
        while i < likelihoods.len() {
            let feature_id = likelihoods[i].0;
            let mut sum = 0.0;
            let mut j = i;
            while j < likelihoods.len() && likelihoods[j].0 == feature_id {
                sum += likelihoods[j].2;
                j += 1;
            }
            let valid = sum.is_finite() && sum > 0.0;
            for entry in &likelihoods[i..j] {
                posteriors.push(PosteriorRow {
                    feature_id,
                    hypothesis: entry.1,
                    posterior: if valid { entry.2 / sum } else { 0.0 },
                });
            }
            i = j;
        }
        posteriors
    }

    /// Perform precursor-layer inference and return precursor peak-group PEPs.
    ///
    /// This step combines peak-group evidence with optional MS1/MS2 precursor
    /// evidence. The returned precursor peak-group PEP is then used as the
    /// prior for the transition/peptidoform layer.
    ///
    /// With neither [`PeptidoformInferenceConfig::ipf_ms1_scoring`] nor
    /// [`PeptidoformInferenceConfig::ipf_ms2_scoring`], every input row whose
    /// `ms2_peakgroup_pep` is strictly below
    /// [`PeptidoformInferenceConfig::ipf_max_precursor_peakgroup_pep`] is
    /// returned unchanged, in input order. Otherwise each evidence source is
    /// thresholded independently — a disabled source, or a PEP at or above
    /// [`PeptidoformInferenceConfig::ipf_max_precursor_pep`], is treated as
    /// absent — the Bayesian model is applied, and the PEP of the true state
    /// (`1 - posterior`) is returned for each feature where it is strictly
    /// below `ipf_max_precursor_peakgroup_pep`, ordered by feature id.
    pub fn precursor_inference(
        &self,
        precursor_rows: &[IPFPrecursorRow],
        config: &PeptidoformInferenceConfig,
    ) -> Vec<IPFPrecursorProbabilityRow> {
        if !(config.ipf_ms1_scoring || config.ipf_ms2_scoring) {
            return precursor_rows
                .iter()
                .filter(|row| row.ms2_peakgroup_pep < config.ipf_max_precursor_peakgroup_pep)
                .map(|row| IPFPrecursorProbabilityRow {
                    feature_id: row.feature_id,
                    precursor_peakgroup_pep: row.ms2_peakgroup_pep,
                })
                .collect();
        }
        // The source resets the evidence when `value >= ipf_max_precursor_pep`;
        // a `NaN` on either side compares false there and keeps the evidence.
        let keep = |enabled: bool, pep: Option<f64>| -> Option<f64> {
            match pep {
                Some(value)
                    if enabled
                        && !matches!(
                            value.partial_cmp(&config.ipf_max_precursor_pep),
                            Some(Ordering::Greater | Ordering::Equal)
                        ) =>
                {
                    Some(value)
                }
                _ => None,
            }
        };
        let filtered: Vec<IPFPrecursorRow> = precursor_rows
            .iter()
            .map(|row| IPFPrecursorRow {
                feature_id: row.feature_id,
                ms2_peakgroup_pep: row.ms2_peakgroup_pep,
                ms1_precursor_pep: keep(config.ipf_ms1_scoring, row.ms1_precursor_pep),
                ms2_precursor_pep: keep(config.ipf_ms2_scoring, row.ms2_precursor_pep),
            })
            .collect();
        let posteriors = Self::apply_bm(&Self::prepare_precursor_bm(&filtered));
        posteriors
            .into_iter()
            .filter(|row| row.hypothesis == 1)
            .map(|row| IPFPrecursorProbabilityRow {
                feature_id: row.feature_id,
                precursor_peakgroup_pep: 1.0 - row.posterior,
            })
            .filter(|row| row.precursor_peakgroup_pep < config.ipf_max_precursor_peakgroup_pep)
            .collect()
    }

    /// Perform complete peptidoform inference.
    ///
    /// Runs [`Self::precursor_inference`], joins each transition row of a
    /// retained feature to its precursor peak-group PEP (the last one when a
    /// feature id repeats) and to its alignment group (the row's own, else the
    /// last alignment row for the feature, else none), optionally propagates
    /// confident evidence across aligned runs, applies the transition-layer
    /// Bayesian model and computes model q-values from the PEPs
    /// (`1 - posterior`) of every hypothesis including H0 — pooled, or per
    /// number of peptidoforms with
    /// [`PeptidoformInferenceConfig::ipf_grouped_fdr`] (a feature's number is
    /// taken from its first transition row with a positive count, `0` if it
    /// has none).
    ///
    /// The transition-layer prior of H0 is the precursor peak-group PEP; each
    /// peptidoform hypothesis gets `(1 - precursor_peakgroup_pep) /
    /// num_peptidoforms`. The evidence is `1 - pep` when `bmask == 1`, `pep`
    /// otherwise.
    ///
    /// # Arguments
    ///
    /// * `precursor_rows` — peak-group/precursor evidence.
    /// * `transition_rows` — transition-level evidence.
    /// * `alignment_rows` — optional alignment-group memberships for
    ///   across-run propagation.
    /// * `config` — IPF configuration.
    ///
    /// # Returns
    ///
    /// The final result rows without the H0 hypothesis, ordered by feature id
    /// and then peptide id; empty when no precursor survives or there are no
    /// transitions.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] (the source's
    /// `Exception::Precondition`) when a non-null transition row reaching the
    /// Bayesian model has `num_peptidoforms <= 0`, and
    /// [`Error::InvalidRange`] when propagation would exceed
    /// [`MAX_PROPAGATED_ROWS`] or the model-FDR input exceeds [`MAX_ITEMS`]
    /// (both native).
    pub fn infer(
        &self,
        precursor_rows: &[IPFPrecursorRow],
        transition_rows: &[IPFTransitionRow],
        alignment_rows: &[IPFAlignmentRow],
        config: &PeptidoformInferenceConfig,
    ) -> Result<Vec<IPFResultRow>> {
        let precursor_probabilities = self.precursor_inference(precursor_rows, config);
        if precursor_probabilities.is_empty() || transition_rows.is_empty() {
            return Ok(Vec::new());
        }
        let mut precursor_pep: HashMap<i64, f64> =
            HashMap::with_capacity(precursor_probabilities.len());
        for row in &precursor_probabilities {
            precursor_pep.insert(row.feature_id, row.precursor_peakgroup_pep);
        }
        let mut alignment_group: HashMap<i64, i64> = HashMap::with_capacity(alignment_rows.len());
        for row in alignment_rows {
            alignment_group.insert(row.feature_id, row.alignment_group_id);
        }

        let mut merged: Vec<TransitionInferenceRow> = Vec::with_capacity(transition_rows.len());
        let mut num_peptidoforms_by_feature: HashMap<i64, i32> = HashMap::new();
        for row in transition_rows {
            let Some(&pep) = precursor_pep.get(&row.feature_id) else {
                continue;
            };
            merged.push(TransitionInferenceRow {
                feature_id: row.feature_id,
                transition_id: row.transition_id,
                hypothesis: row.peptide_id,
                pep: row.pep,
                precursor_peakgroup_pep: pep,
                bmask: row.bmask,
                num_peptidoforms: row.num_peptidoforms,
                alignment_group_id: row
                    .alignment_group_id
                    .unwrap_or_else(|| alignment_group.get(&row.feature_id).copied().unwrap_or(-1)),
            });
            if row.num_peptidoforms > 0 {
                num_peptidoforms_by_feature
                    .entry(row.feature_id)
                    .or_insert(row.num_peptidoforms);
            }
        }
        if merged.is_empty() {
            return Ok(Vec::new());
        }
        if config.propagate_signal_across_runs {
            merged = propagate_signal_across_runs(&merged, config.across_run_confidence_threshold)?;
        }

        let bm_rows = prepare_transition_bm(&merged)?;
        let posteriors = Self::apply_bm(&bm_rows);
        if posteriors.is_empty() {
            return Ok(Vec::new());
        }
        let pep_values: Vec<f64> = posteriors.iter().map(|row| 1.0 - row.posterior).collect();
        let qvalues = if config.ipf_grouped_fdr {
            let mut groups: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
            for (index, row) in posteriors.iter().enumerate() {
                let count = num_peptidoforms_by_feature
                    .get(&row.feature_id)
                    .copied()
                    .unwrap_or(0);
                groups.entry(count).or_default().push(index);
            }
            let mut qvalues = vec![1.0; posteriors.len()];
            for indices in groups.values() {
                let group_peps: Vec<f64> = indices.iter().map(|&i| pep_values[i]).collect();
                let group_q = Self::compute_model_fdr(&group_peps)?;
                for (slot, &index) in indices.iter().enumerate() {
                    qvalues[index] = group_q[slot];
                }
            }
            qvalues
        } else {
            Self::compute_model_fdr(&pep_values)?
        };

        let mut results = Vec::with_capacity(posteriors.len());
        for (index, row) in posteriors.iter().enumerate() {
            if row.hypothesis == -1 {
                continue;
            }
            let Some(&pep) = precursor_pep.get(&row.feature_id) else {
                continue;
            };
            results.push(IPFResultRow {
                feature_id: row.feature_id,
                peptide_id: row.hypothesis,
                precursor_peakgroup_pep: pep,
                qvalue: qvalues[index],
                pep: pep_values[index],
            });
        }
        results.sort_by_key(|row| (row.feature_id, row.peptide_id));
        Ok(results)
    }
}

/// `std::min(a, b)`: `b < a ? b : a`, so a `NaN` second argument keeps `a`.
fn f64_min(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}

/// The source's file-static `propagateSignalAcrossRuns_`.
///
/// Rows without an alignment group (negative id) pass through in input order.
/// For each alignment group (ascending id here, see the module notes), every
/// distinct feature of the group receives its own rows plus every other row of
/// the group whose PEP is at or below `confidence_threshold`, relabelled with
/// its feature id; duplicates by the source's sort key are collapsed into one
/// row carrying the minimum PEP and the minimum precursor PEP.
fn propagate_signal_across_runs(
    rows: &[TransitionInferenceRow],
    confidence_threshold: f64,
) -> Result<Vec<TransitionInferenceRow>> {
    let mut result: Vec<TransitionInferenceRow> = Vec::with_capacity(rows.len());
    let mut groups: BTreeMap<i64, Vec<&TransitionInferenceRow>> = BTreeMap::new();
    for row in rows {
        if row.alignment_group_id < 0 {
            result.push(row.clone());
        } else {
            groups.entry(row.alignment_group_id).or_default().push(row);
        }
    }

    let mut group_features: Vec<Vec<i64>> = Vec::with_capacity(groups.len());
    let mut total = result.len();
    for members in groups.values() {
        let mut features: Vec<i64> = members.iter().map(|r| r.feature_id).collect();
        features.sort_unstable();
        features.dedup();
        let confident = members
            .iter()
            .filter(|r| r.pep <= confidence_threshold)
            .count();
        // Upper bound of the rows this group emits: each feature keeps its own
        // rows and imports every confident row of the others.
        let bound = features
            .len()
            .checked_mul(confident)
            .and_then(|n| n.checked_add(members.len()))
            .ok_or_else(|| too_many_propagated(usize::MAX))?;
        total = total
            .checked_add(bound)
            .ok_or_else(|| too_many_propagated(usize::MAX))?;
        if total > MAX_PROPAGATED_ROWS {
            return Err(too_many_propagated(total));
        }
        group_features.push(features);
    }

    for (members, features) in groups.values().zip(&group_features) {
        let mut propagated: Vec<TransitionInferenceRow> = Vec::new();
        for &target in features {
            for row in members {
                if row.feature_id == target || row.pep <= confidence_threshold {
                    let mut copy = (*row).clone();
                    copy.feature_id = target;
                    propagated.push(copy);
                }
            }
        }
        propagated.sort_by_key(TransitionInferenceRow::key);
        for row in propagated {
            if let Some(back) = result.last_mut() {
                if back.key() == row.key() {
                    back.pep = f64_min(back.pep, row.pep);
                    back.precursor_peakgroup_pep =
                        f64_min(back.precursor_peakgroup_pep, row.precursor_peakgroup_pep);
                    continue;
                }
            }
            result.push(row);
        }
    }
    Ok(result)
}

fn too_many_propagated(total: usize) -> Error {
    Error::InvalidRange(format!(
        "OpenSwathPeptidoformInference: across-run propagation would produce {total} rows, \
         above the ceiling of {MAX_PROPAGATED_ROWS}"
    ))
}

/// The source's file-static `prepareTransitionBM_`.
fn prepare_transition_bm(rows: &[TransitionInferenceRow]) -> Result<Vec<BayesianModelRow>> {
    let mut bm_rows = Vec::with_capacity(rows.len());
    for row in rows {
        if row.hypothesis != -1 && row.num_peptidoforms <= 0 {
            return Err(Error::InvalidValue(
                "Transition inference requires num_peptidoforms > 0 for non-null hypotheses."
                    .to_string(),
            ));
        }
        let prior = if row.hypothesis == -1 {
            row.precursor_peakgroup_pep
        } else {
            (1.0 - row.precursor_peakgroup_pep) / f64::from(row.num_peptidoforms)
        };
        let evidence = if row.bmask == 1 {
            1.0 - row.pep
        } else {
            row.pep
        };
        bm_rows.push(BayesianModelRow {
            feature_id: row.feature_id,
            hypothesis: row.hypothesis,
            prior,
            evidence,
        });
    }
    Ok(bm_rows)
}
