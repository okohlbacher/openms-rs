// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `ANALYSIS/ID/ConsensusIDAlgorithm*.h`: consensus peptide
//! identification from several search runs of one spectrum or feature.
//!
//! | C++ class | Rust |
//! |---|---|
//! | `ConsensusIDAlgorithm` (abstract base) | [`ConsensusIDAlgorithm`] trait, [`ConsensusIDFilter`] |
//! | `ConsensusIDAlgorithmIdentity` (abstract) | the shared grouping behind the four identity types |
//! | `ConsensusIDAlgorithmSimilarity` (abstract) | the shared grouping behind the two similarity types |
//! | `ConsensusIDAlgorithmAverage` / `Best` / `Worst` / `Ranks` | [`ConsensusIDAlgorithmAverage`], [`ConsensusIDAlgorithmBest`], [`ConsensusIDAlgorithmWorst`], [`ConsensusIDAlgorithmRanks`] |
//! | `ConsensusIDAlgorithmPEPMatrix` / `PEPIons` | [`ConsensusIDAlgorithmPEPMatrix`], [`ConsensusIDAlgorithmPEPIons`] |
//!
//! The C++ abstract classes exist to share code between subclasses through
//! protected virtual hooks (`apply_`, `preprocess_`, `getAggregateScore_`,
//! `getSimilarity_`). Rust has no inheritance, so each concrete type
//! implements [`ConsensusIDAlgorithm`] and the hooks are private functions of
//! this module; the protected `HitInfo`/`SequenceGrouping` state is internal.
//!
//! Each type carries its parameters as fields and converts to and from the
//! source parameter tree with `defaults`, `from_param`,
//! `from_param_with_warnings` and `to_param`, under the source handler names
//! (`ConsensusIDAlgorithmBest`, ...).
//!
//! See the `ConsensusID` TOPP tool documentation for the `filter:` parameters.

use super::needleman_wunsch::{NeedlemanWunsch, ScoringMatrix};
use crate::analysis::id_filter::{DuplicatePolicy, remove_duplicate_peptide_hits};
use crate::chemistry::{AASequence, EmpiricalFormula, PeptideFragmentType, SequenceModification};
use crate::identification::{PeptideEvidence, PeptideHit, PeptideIdentification};
use crate::metadata::MetaValue;
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

/// Most peptide hits, summed over all input identifications, one consensus
/// computation accepts.
pub const MAX_CONSENSUS_HITS: usize = 1_000_000;

/// Most sequence-similarity evaluations one similarity-based consensus
/// (PEPMatrix, PEPIons) performs: for each distinct sequence, one per hit of
/// every other identification. Checked before any work is done.
pub const MAX_SIMILARITY_EVALUATIONS: usize = 1_000_000;

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

// ---------------------------------------------------------------------------
// Shared filter parameters (the base class)
// ---------------------------------------------------------------------------

/// The `filter:` parameters every consensus algorithm shares (source
/// `ConsensusIDAlgorithm` members `considered_hits_`, `min_support_`,
/// `count_empty_`, `keep_old_scores_`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConsensusIDFilter {
    /// `filter:considered_hits`: the number of top hits in each ID run that are
    /// considered for consensus scoring; 0 considers all hits.
    pub considered_hits: usize,
    /// `filter:min_support`: for each peptide hit from an ID run, the fraction
    /// (0 to 1) of other ID runs that must support that hit, otherwise it is
    /// removed.
    pub min_support: f64,
    /// `filter:count_empty`: count empty ID runs (those with no peptide hit
    /// for the current spectrum) when calculating `min_support`.
    pub count_empty: bool,
    /// `filter:keep_old_scores`: keep the original scores as meta values named
    /// `<score type>_score` on each consensus hit.
    pub keep_old_scores: bool,
}

impl Default for ConsensusIDFilter {
    fn default() -> Self {
        Self {
            considered_hits: 0,
            min_support: 0.0,
            count_empty: false,
            keep_old_scores: false,
        }
    }
}

impl ConsensusIDFilter {
    fn write(&self, param: &mut Param) -> Result<()> {
        let considered = i64::try_from(self.considered_hits)
            .ok()
            .filter(|v| *v <= i64::from(i32::MAX))
            .ok_or_else(|| {
                bad("'filter:considered_hits' does not fit the source 32-bit parameter")
            })?;
        param.set_value(
            "filter:considered_hits",
            ParamValue::Integer(considered),
            "The number of top hits in each ID run that are considered for consensus scoring ('0' for all hits).",
            &[],
        )?;
        param.set_min_int("filter:considered_hits", 0)?;
        param.set_value(
            "filter:min_support",
            ParamValue::Float(self.min_support),
            "For each peptide hit from an ID run, the fraction of other ID runs that must support that hit (otherwise it is removed).",
            &[],
        )?;
        param.set_min_float("filter:min_support", 0.0)?;
        param.set_max_float("filter:min_support", 1.0)?;
        flag_value(
            param,
            "filter:count_empty",
            self.count_empty,
            "Count empty ID runs (i.e. those containing no peptide hit for the current spectrum) when calculating 'min_support'?",
        )?;
        flag_value(
            param,
            "filter:keep_old_scores",
            self.keep_old_scores,
            "if set, keeps the original scores as user params",
        )
    }

    fn read(param: &Param) -> Result<Self> {
        let considered = param.value("filter:considered_hits")?.to_i64()?;
        Ok(Self {
            considered_hits: usize::try_from(considered)
                .map_err(|_| bad("'filter:considered_hits' must not be negative"))?,
            min_support: param.value("filter:min_support")?.to_f64()?,
            count_empty: param.value("filter:count_empty")?.as_str()? == "true",
            keep_old_scores: param.value("filter:keep_old_scores")?.as_str()? == "true",
        })
    }
}

fn flag_value(param: &mut Param, key: &str, value: bool, description: &str) -> Result<()> {
    param.set_value(
        key,
        ParamValue::String(if value { "true" } else { "false" }.into()),
        description,
        &[],
    )?;
    param.set_valid_strings(key, &["true".to_owned(), "false".to_owned()])
}

// ---------------------------------------------------------------------------
// The algorithm interface
// ---------------------------------------------------------------------------

/// Calculates a consensus peptide identification from several identification
/// runs of one spectrum or (consensus) feature (source abstract base
/// `ConsensusIDAlgorithm`).
///
/// Make sure the score type and the score orientation
/// ([`PeptideIdentification::score_type`],
/// [`PeptideIdentification::higher_score_better`]) are set properly on the
/// input.
pub trait ConsensusIDAlgorithm {
    /// The shared `filter:` parameters.
    fn filter(&self) -> &ConsensusIDFilter;

    /// Replace `ids` (input: normally more than one identification) by their
    /// single consensus identification.
    ///
    /// `se_info` maps a run identifier ([`PeptideIdentification::identifier`])
    /// to its search engine name; a matching engine prefixes the score type
    /// recorded for `filter:keep_old_scores` (`<engine>_<score type>`).
    /// `number_of_runs` is the number of ID runs, 0 meaning `ids.len()`.
    ///
    /// Each input identification is sorted by score, cut to
    /// `filter:considered_hits` and freed of duplicate sequences before the
    /// algorithm runs. The output identification carries only the score type
    /// and orientation of the first input (after any algorithm preprocessing);
    /// identifier, RT, m/z and meta values are not carried over. Its hits are
    /// sorted by score, and each carries the meta value `consensus_support`,
    /// `target_decoy` when the first hit of that sequence had one, and the
    /// peptide evidences of every hit of that sequence. Hits whose support is
    /// below `filter:min_support` are dropped, so the output identification
    /// may have no hits.
    ///
    /// Empty input is left unchanged. Returns the warnings the source logs.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when an identification or hit fails
    /// validation, the same sequence occurs with two different non-zero
    /// charges, the input exceeds [`MAX_CONSENSUS_HITS`], or the
    /// algorithm-specific conditions its type documents. On error `ids` is
    /// unchanged.
    fn apply_with_search_engines(
        &self,
        ids: &mut Vec<PeptideIdentification>,
        se_info: &BTreeMap<String, String>,
        number_of_runs: usize,
    ) -> Result<Vec<String>>;

    /// [`ConsensusIDAlgorithm::apply_with_search_engines`] with no search
    /// engine information.
    ///
    /// # Errors
    ///
    /// As [`ConsensusIDAlgorithm::apply_with_search_engines`].
    fn apply(
        &self,
        ids: &mut Vec<PeptideIdentification>,
        number_of_runs: usize,
    ) -> Result<Vec<String>> {
        self.apply_with_search_engines(ids, &BTreeMap::new(), number_of_runs)
    }
}

/// Source `ConsensusIDAlgorithm::HitInfo`.
#[derive(Clone, Debug)]
struct HitInfo {
    charge: i32,
    scores: Vec<f64>,
    types: Vec<String>,
    target_decoy: String,
    evidence: BTreeSet<PeptideEvidence>,
    final_score: f64,
    support: f64,
}

/// An [`AASequence`] ordered as the source `AASequence::operator<`, which
/// orders the source `SequenceGrouping` map and therefore the output of equal
/// scores.
#[derive(Clone, Debug, PartialEq, Eq)]
struct SourceOrdered(AASequence);

impl Ord for SourceOrdered {
    fn cmp(&self, other: &Self) -> Ordering {
        source_sequence_cmp(&self.0, &other.0)
    }
}

impl PartialOrd for SourceOrdered {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Source `SequenceGrouping`: sequence -> hit information.
type Grouping = BTreeMap<SourceOrdered, HitInfo>;

/// Source `AASequence::operator<`: length first, then the N-terminal
/// modification (none before any), then residue by residue the one-letter
/// code and the modification, then the C-terminal modification.
///
/// The source compares two residue modifications by *pointer*, an order that
/// depends on allocation; this compares them by name and full ID instead, and
/// finally by the native [`AASequence`] order, so that the order is total and
/// deterministic.
fn source_sequence_cmp(a: &AASequence, b: &AASequence) -> Ordering {
    fn modification(
        m: Option<&SequenceModification>,
        n: Option<&SequenceModification>,
    ) -> Ordering {
        match (m, n) {
            (None, None) => Ordering::Equal,
            (None, Some(_)) => Ordering::Less,
            (Some(_), None) => Ordering::Greater,
            (Some(m), Some(n)) => m
                .name()
                .cmp(n.name())
                .then_with(|| m.full_id().cmp(n.full_id())),
        }
    }
    a.len()
        .cmp(&b.len())
        .then_with(|| modification(a.n_terminal_modification(), b.n_terminal_modification()))
        .then_with(|| {
            let (x, y) = (a.as_str().as_bytes(), b.as_str().as_bytes());
            for index in 0..x.len() {
                let order = x[index].cmp(&y[index]).then_with(|| {
                    modification(
                        a.residue_modification(index).ok().flatten(),
                        b.residue_modification(index).ok().flatten(),
                    )
                });
                if order != Ordering::Equal {
                    return order;
                }
            }
            Ordering::Equal
        })
        .then_with(|| modification(a.c_terminal_modification(), b.c_terminal_modification()))
        .then_with(|| a.cmp(b))
}

/// Source `compareChargeStates_`: a zero charge adopts the new one; two
/// different non-zero charges are an error.
fn compare_charge_states(recorded: &mut i32, new: i32, peptide: &AASequence) -> Result<()> {
    if *recorded == 0 {
        *recorded = new;
    } else if new != 0 && *recorded != new {
        return Err(bad(format!(
            "Conflicting charge states found for peptide '{peptide}': {recorded}, {new}"
        )));
    }
    Ok(())
}

fn target_decoy_text(hit: &PeptideHit) -> String {
    hit.metadata
        .get("target_decoy")
        .map(ToString::to_string)
        .unwrap_or_default()
}

fn prefixed_score_type(id: &PeptideIdentification, se_info: &BTreeMap<String, String>) -> String {
    match se_info.get(&id.identifier) {
        Some(engine) => format!("{engine}_{}", id.score_type),
        None => id.score_type.clone(),
    }
}

/// The algorithm-specific step, source `apply_`.
enum Method<'a> {
    Identity(&'a dyn Fn(&[f64], bool) -> f64),
    Ranks,
    PepMatrix(NeedlemanWunsch),
    PepIons {
        mass_tolerance: f64,
        min_shared: usize,
    },
}

/// Source `ConsensusIDAlgorithm::apply`, run on a copy so that an error
/// leaves `ids` unchanged.
fn run(
    filter: &ConsensusIDFilter,
    method: Method<'_>,
    ids: &mut Vec<PeptideIdentification>,
    se_info: &BTreeMap<String, String>,
    number_of_runs: usize,
) -> Result<Vec<String>> {
    if ids.is_empty() {
        return Ok(Vec::new());
    }
    let mut total = 0_usize;
    for id in ids.iter() {
        id.validate()?;
        total = total.saturating_add(id.hits.len());
    }
    if total > MAX_CONSENSUS_HITS {
        return Err(bad(format!(
            "consensus input of {total} peptide hits exceeds {MAX_CONSENSUS_HITS}"
        )));
    }
    let number_of_runs = if number_of_runs != 0 {
        number_of_runs
    } else {
        ids.len()
    };
    let mut work = ids.clone();
    for id in &mut work {
        id.sort()?;
        if filter.considered_hits > 0 && id.hits.len() > filter.considered_hits {
            id.hits.truncate(filter.considered_hits);
        }
    }
    remove_duplicate_peptide_hits(&mut work, DuplicatePolicy::Sequence)?;

    let mut warnings = Vec::new();
    let results = match method {
        Method::Identity(aggregate) => {
            identity_preprocess(&work, &mut warnings)?;
            identity_group(filter, &work, se_info, number_of_runs, aggregate)?
        }
        Method::Ranks => {
            let (runs, considered) = ranks_preprocess(filter, &mut work, number_of_runs)?;
            let aggregate = |scores: &[f64], _: bool| -> f64 {
                match runs.checked_sub(scores.len()) {
                    Some(missing) => {
                        let sum: f64 = scores.iter().sum::<f64>() + (missing * considered) as f64;
                        1.0 - sum / (considered * runs) as f64
                    }
                    None => f64::NAN,
                }
            };
            let grouping = identity_group(filter, &work, se_info, number_of_runs, &aggregate)?;
            if let Some(found) = grouping
                .values()
                .find(|info| info.final_score.is_nan())
                .map(|info| info.scores.len())
            {
                return Err(bad(format!(
                    "a sequence was found in {found} identifications but only {runs} runs were given; the source's unsigned subtraction wraps here"
                )));
            }
            grouping
        }
        Method::PepMatrix(aligner) => {
            let mut similarity =
                |a: &AASequence, b: &AASequence| pep_matrix_similarity(&aligner, a, b);
            similarity_group(filter, &work, se_info, number_of_runs, &mut similarity)?
        }
        Method::PepIons {
            mass_tolerance,
            min_shared,
        } => {
            let water = EmpiricalFormula::parse("H2O")?.mono_mass();
            let mut cache = BTreeMap::new();
            let mut similarity = |a: &AASequence, b: &AASequence| {
                pep_ions_similarity(&mut cache, water, mass_tolerance, min_shared, a, b)
            };
            similarity_group(filter, &work, se_info, number_of_runs, &mut similarity)?
        }
    };

    let mut output = PeptideIdentification {
        score_type: work[0].score_type.clone(),
        higher_score_better: work[0].higher_score_better,
        ..PeptideIdentification::default()
    };
    for (sequence, info) in results {
        if info.support < filter.min_support {
            continue;
        }
        let mut hit = PeptideHit {
            sequence: sequence.0,
            charge: info.charge,
            score: info.final_score,
            evidences: info.evidence.into_iter().collect(),
            ..PeptideHit::default()
        };
        hit.metadata.insert(
            "consensus_support".into(),
            MetaValue::try_from(info.support)?,
        );
        if !info.target_decoy.is_empty() {
            hit.metadata
                .insert("target_decoy".into(), MetaValue::from(info.target_decoy));
        }
        if filter.keep_old_scores {
            for (kind, score) in info.types.iter().zip(&info.scores) {
                hit.metadata
                    .insert(format!("{kind}_score"), MetaValue::try_from(*score)?);
            }
        }
        output.hits.push(hit);
    }
    output.sort()?;
    *ids = vec![output];
    Ok(warnings)
}

// ---------------------------------------------------------------------------
// Identity-based grouping (source ConsensusIDAlgorithmIdentity)
// ---------------------------------------------------------------------------

/// Source `ConsensusIDAlgorithmIdentity::preprocess_`: score orientations
/// must agree (error), score types should (warning).
fn identity_preprocess(ids: &[PeptideIdentification], warnings: &mut Vec<String>) -> Result<()> {
    let higher_better = ids[0].higher_score_better;
    let mut types = BTreeSet::new();
    for id in ids {
        if id.higher_score_better != higher_better {
            let hi_lo = if higher_better {
                "higher/lower"
            } else {
                "lower/higher"
            };
            return Err(bad(format!(
                "Score types '{}' and '{}' have different orientations ({hi_lo} is better) and cannot be compared meaningfully.",
                ids[0].score_type, id.score_type
            )));
        }
        types.insert(id.score_type.as_str());
    }
    if types.len() > 1 {
        let joined = types.into_iter().collect::<Vec<_>>().join("'/'");
        warnings.push(format!(
            "Warning: Different score types for peptide hits found ('{joined}'). If the scores are not comparable, results will be meaningless."
        ));
    }
    Ok(())
}

/// Source `ConsensusIDAlgorithmIdentity::apply_` after `preprocess_`.
fn identity_group(
    filter: &ConsensusIDFilter,
    ids: &[PeptideIdentification],
    se_info: &BTreeMap<String, String>,
    number_of_runs: usize,
    aggregate: &dyn Fn(&[f64], bool) -> f64,
) -> Result<Grouping> {
    let mut results = Grouping::new();
    for id in ids {
        let score_type = prefixed_score_type(id, se_info);
        for hit in &id.hits {
            let key = SourceOrdered(hit.sequence.clone());
            if let Some(info) = results.get_mut(&key) {
                compare_charge_states(&mut info.charge, hit.charge, &hit.sequence)?;
                info.scores.push(hit.score);
                info.types.push(score_type.clone());
                info.evidence.extend(hit.evidences.iter().cloned());
            } else {
                results.insert(
                    key,
                    HitInfo {
                        charge: hit.charge,
                        scores: vec![hit.score],
                        types: vec![score_type.clone()],
                        target_decoy: target_decoy_text(hit),
                        evidence: hit.evidences.iter().cloned().collect(),
                        final_score: 0.0,
                        support: 0.0,
                    },
                );
            }
        }
    }
    let higher_better = ids[0].higher_score_better;
    // `number_of_runs` and `ids.len()` are both at least one here.
    let n_other_ids = if filter.count_empty {
        number_of_runs
    } else {
        ids.len()
    } - 1;
    for info in results.values_mut() {
        info.final_score = aggregate(&info.scores, higher_better);
        // With 'count_empty' false, 'n_other_ids' may be zero; the support is
        // then defined as one to avoid a NaN.
        info.support = if n_other_ids > 0 {
            (info.scores.len() as f64 - 1.0) / n_other_ids as f64
        } else {
            1.0
        };
    }
    Ok(results)
}

/// Source `ConsensusIDAlgorithmRanks::preprocess_`, returning the current
/// number of runs and of considered hits.
fn ranks_preprocess(
    filter: &ConsensusIDFilter,
    ids: &mut [PeptideIdentification],
    number_of_runs: usize,
) -> Result<(usize, usize)> {
    let runs = if number_of_runs > 0 {
        number_of_runs
    } else {
        ids.len()
    };
    let mut considered = filter.considered_hits;
    for id in ids.iter_mut() {
        id.sort()?;
        for (rank, hit) in id.hits.iter_mut().enumerate() {
            hit.score = rank as f64;
        }
        id.score_type = "ConsensusID_ranks".into();
        id.higher_score_better = true; // not true now, but after normalizing
        if filter.considered_hits == 0 && id.hits.len() > considered {
            considered = id.hits.len();
        }
    }
    Ok((runs, considered))
}

// ---------------------------------------------------------------------------
// Similarity-based grouping (source ConsensusIDAlgorithmSimilarity)
// ---------------------------------------------------------------------------

/// Source `ConsensusIDAlgorithmSimilarity::apply_`.
fn similarity_group(
    filter: &ConsensusIDFilter,
    ids: &[PeptideIdentification],
    se_info: &BTreeMap<String, String>,
    number_of_runs: usize,
    similarity: &mut dyn FnMut(&AASequence, &AASequence) -> Result<f64>,
) -> Result<Grouping> {
    for id in ids {
        if !matches!(
            id.score_type.as_str(),
            "Posterior Error Probability" | "pep" | "MS:1001493"
        ) {
            return Err(bad(format!(
                "Score type must be 'Posterior Error Probability' (was '{}')",
                id.score_type
            )));
        }
    }
    let total: usize = ids.iter().map(|id| id.hits.len()).sum();
    let evaluations = ids
        .iter()
        .map(|id| id.hits.len().saturating_mul(total - id.hits.len()))
        .fold(0_usize, usize::saturating_add);
    if evaluations > MAX_SIMILARITY_EVALUATIONS {
        return Err(bad(format!(
            "similarity consensus needs {evaluations} evaluations, more than {MAX_SIMILARITY_EVALUATIONS}"
        )));
    }

    let mut results = Grouping::new();
    for (index1, id1) in ids.iter().enumerate() {
        let score_type = prefixed_score_type(id1, se_info);
        for hit1 in &id1.hits {
            let key = SourceOrdered(hit1.sequence.clone());
            if let Some(info) = results.get_mut(&key) {
                // Already scored: the source records the unprefixed score type
                // on this path.
                compare_charge_states(&mut info.charge, hit1.charge, &hit1.sequence)?;
                info.scores.push(hit1.score);
                info.types.push(id1.score_type.clone());
                info.evidence.extend(hit1.evidences.iter().cloned());
                continue;
            }
            // Similarity and "1 - PEP" of the best match in every other run:
            // highest similarity, ties broken by the better PEP.
            let mut best_matches: Vec<(f64, f64)> = Vec::with_capacity(ids.len() - 1);
            for (index2, id2) in ids.iter().enumerate() {
                if index1 == index2 || id2.hits.is_empty() {
                    continue;
                }
                let mut best: Option<(f64, f64)> = None;
                for hit2 in &id2.hits {
                    let current = (
                        similarity(&hit1.sequence, &hit2.sequence)?,
                        1.0 - hit2.score,
                    );
                    // std::max_element: the first of equal maxima wins.
                    best = match best {
                        Some(b) if !(b.0 < current.0 || (b.0 == current.0 && b.1 < current.1)) => {
                            Some(b)
                        }
                        _ => Some(current),
                    };
                }
                best_matches.extend(best);
            }
            let mut score = hit1.score;
            let mut sum_sim = 1.0;
            for (sim, one_minus_pep) in &best_matches {
                score += sim * (1.0 - one_minus_pep); // undo "1 - PEP"
                sum_sim += sim;
            }
            score /= sum_sim * sum_sim;
            let n_other_ids = if filter.count_empty {
                number_of_runs - 1
            } else {
                best_matches.len()
            };
            let support = if n_other_ids == 0 {
                // Only one ID run: similarity is ill-defined.
                if filter.count_empty { 0.0 } else { 1.0 }
            } else {
                (sum_sim - 1.0) / n_other_ids as f64
            };
            if !score.is_finite() || !support.is_finite() {
                return Err(bad(format!(
                    "consensus score or support for '{}' is not finite",
                    hit1.sequence
                )));
            }
            results.insert(
                key,
                HitInfo {
                    charge: hit1.charge,
                    scores: vec![hit1.score],
                    types: vec![score_type.clone()],
                    target_decoy: target_decoy_text(hit1),
                    evidence: hit1.evidences.iter().cloned().collect(),
                    final_score: score,
                    support,
                },
            );
        }
    }
    Ok(results)
}

/// Source `ConsensusIDAlgorithmPEPMatrix::getSimilarity_`: alignment score of
/// the unmodified sequences, normalized by the smaller self-alignment score.
fn pep_matrix_similarity(aligner: &NeedlemanWunsch, a: &AASequence, b: &AASequence) -> Result<f64> {
    // Modifications cannot be taken into account here.
    let (x, y) = (a.as_str(), b.as_str());
    if x == y {
        return Ok(1.0);
    }
    let score = aligner.align(x, y)?;
    if score < 0 {
        return Ok(0.0);
    }
    let normalizer = aligner.align(x, x)?.min(aligner.align(y, y)?);
    let similarity = score as f64 / normalizer as f64;
    if !similarity.is_finite() {
        return Err(bad(format!(
            "PEPMatrix similarity of '{x}' and '{y}' is not finite: a self-alignment score is zero"
        )));
    }
    Ok(similarity)
}

/// b and y ion masses of one sequence, as the source builds them: the first
/// b ion is the full mass of the one-residue prefix (with N-terminal
/// modification), the first y ion that of the one-residue suffix (with
/// C-terminal modification), and each further ion adds the *full* residue
/// mass (internal mass plus water) of the next residue.
fn pep_ions_ladder(sequence: &AASequence, water: f64) -> Result<Vec<f64>> {
    let n = sequence.len();
    if n == 0 {
        return Err(bad(
            "PEPIons similarity needs non-empty sequences (the source writes out of bounds)",
        ));
    }
    let residue = |index: usize| -> Result<f64> {
        sequence
            .internal_mass(index)?
            .map(|mass| mass + water)
            .ok_or_else(|| bad(format!("'{sequence}' has a residue without a known mass")))
    };
    let mut ions = vec![0.0; 2 * n];
    ions[0] = sequence
        .prefix(1)?
        .mono_mass_for(PeptideFragmentType::Full, 0)?;
    ions[n] = sequence
        .suffix(1)?
        .mono_mass_for(PeptideFragmentType::Full, 0)?;
    for i in 1..n {
        ions[i] = ions[i - 1] + residue(i)?;
        ions[n + i] = ions[n + i - 1] + residue(n - i - 1)?;
    }
    ions.sort_by(f64::total_cmp);
    Ok(ions)
}

/// Source `ConsensusIDAlgorithmPEPIons::getSimilarity_`: shared peak count of
/// the b/y ion ladders.
fn pep_ions_similarity(
    cache: &mut BTreeMap<(SourceOrdered, SourceOrdered), f64>,
    water: f64,
    mass_tolerance: f64,
    min_shared: usize,
    a: &AASequence,
    b: &AASequence,
) -> Result<f64> {
    if a == b {
        return Ok(1.0);
    }
    // The measure is not symmetric; the source orders the pair first.
    let (a, b) = if source_sequence_cmp(b, a) == Ordering::Less {
        (b, a)
    } else {
        (a, b)
    };
    let key = (SourceOrdered(a.clone()), SourceOrdered(b.clone()));
    if let Some(&cached) = cache.get(&key) {
        return Ok(cached);
    }
    let ions1 = pep_ions_ladder(a, water)?;
    let ions2 = pep_ions_ladder(b, water)?;
    // Each best-matching ion of the second sequence counts once.
    let mut matches: Vec<f64> = Vec::new();
    let mut start = 0_usize;
    for &ion in &ions1 {
        let lower = start + ions2[start..].partition_point(|&x| x < ion - mass_tolerance);
        if lower == ions2.len() {
            break; // all remaining values are too low
        }
        let upper = lower + ions2[lower..].partition_point(|&x| x <= ion + mass_tolerance);
        let mut best_match = 0.0;
        let mut best_diff = mass_tolerance + 1.0;
        for &candidate in &ions2[lower..upper] {
            let diff = (ion - candidate).abs();
            if diff < best_diff {
                best_diff = diff;
                best_match = candidate;
            }
        }
        if best_diff <= mass_tolerance {
            matches.push(best_match);
        }
        start = lower;
    }
    matches.sort_by(f64::total_cmp);
    matches.dedup();
    let similarity = if matches.len() >= min_shared {
        // Source `matches.size() / float(min(...))`: single precision.
        f64::from(matches.len() as f32 / ions1.len().min(ions2.len()) as f32)
    } else {
        0.0
    };
    cache.insert(key, similarity);
    Ok(similarity)
}

// ---------------------------------------------------------------------------
// Concrete algorithms
// ---------------------------------------------------------------------------

macro_rules! handler_param_api {
    ($name:literal) => {
        /// The default parameters, in source declaration order: the
        /// `filter:` parameters of [`ConsensusIDFilter`], then this
        /// algorithm's own.
        ///
        /// # Errors
        ///
        /// Only a parameter-tree resource failure.
        pub fn defaults() -> Result<Param> {
            Self::default().to_param()
        }

        /// Build from parameters, as source `setParameters` followed by
        /// `updateMembers_`. Missing parameters take their defaults; unknown
        /// ones are ignored (see `from_param_with_warnings`).
        ///
        /// # Errors
        ///
        /// Returns [`Error::InvalidValue`] when a parameter has the wrong type
        /// or violates its restriction (source `Exception::InvalidParameter`).
        pub fn from_param(parameters: &Param) -> Result<Self> {
            Ok(Self::from_param_with_warnings(parameters)?.0)
        }

        /// As `from_param`, also returning the unknown-parameter warnings the
        /// source logs.
        ///
        /// # Errors
        ///
        /// As `from_param`.
        pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
            let mut handler = DefaultParamHandler::new($name)?;
            handler.set_defaults(Self::defaults()?)?;
            handler.set_parameters_with(parameters, Self::from_complete_param)
        }
    };
}

macro_rules! identity_algorithm {
    ($(#[$doc:meta])* $type:ident, $name:literal, $aggregate:expr) => {
        $(#[$doc])*
        ///
        /// Built on the identity grouping: see [`ConsensusIDAlgorithm`].
        #[derive(Clone, Copy, Debug, Default, PartialEq)]
        pub struct $type {
            /// The shared `filter:` parameters.
            pub filter: ConsensusIDFilter,
        }

        impl $type {
            handler_param_api!($name);

            /// The parameters in the layout of `defaults`, carrying the
            /// current values (source `getParameters`).
            ///
            /// # Errors
            ///
            /// Returns [`Error::InvalidValue`] when `considered_hits` does not
            /// fit the source 32-bit parameter.
            pub fn to_param(&self) -> Result<Param> {
                let mut param = Param::new();
                self.filter.write(&mut param)?;
                Ok(param)
            }

            fn from_complete_param(param: &Param) -> Result<Self> {
                Ok(Self {
                    filter: ConsensusIDFilter::read(param)?,
                })
            }
        }

        impl ConsensusIDAlgorithm for $type {
            fn filter(&self) -> &ConsensusIDFilter {
                &self.filter
            }

            fn apply_with_search_engines(
                &self,
                ids: &mut Vec<PeptideIdentification>,
                se_info: &BTreeMap<String, String>,
                number_of_runs: usize,
            ) -> Result<Vec<String>> {
                let aggregate: &dyn Fn(&[f64], bool) -> f64 = &$aggregate;
                run(
                    &self.filter,
                    Method::Identity(aggregate),
                    ids,
                    se_info,
                    number_of_runs,
                )
            }
        }
    };
}

identity_algorithm!(
    /// Consensus score is the average of the search-engine scores of each
    /// sequence (source `ConsensusIDAlgorithmAverage`).
    ///
    /// Only identical sequences are scored together. Score orientations of all
    /// inputs must agree, and the score types should.
    ///
    /// # Errors
    ///
    /// Besides those of [`ConsensusIDAlgorithm::apply_with_search_engines`],
    /// [`Error::InvalidValue`] when the score orientations disagree.
    ConsensusIDAlgorithmAverage,
    "ConsensusIDAlgorithmAverage",
    |scores: &[f64], _higher_better: bool| scores.iter().sum::<f64>() / scores.len() as f64
);

identity_algorithm!(
    /// Consensus score is the best search-engine score of each sequence
    /// (source `ConsensusIDAlgorithmBest`).
    ///
    /// Only identical sequences are scored together. Score orientations of all
    /// inputs must agree, and the score types should.
    ///
    /// # Errors
    ///
    /// Besides those of [`ConsensusIDAlgorithm::apply_with_search_engines`],
    /// [`Error::InvalidValue`] when the score orientations disagree.
    ConsensusIDAlgorithmBest,
    "ConsensusIDAlgorithmBest",
    |scores: &[f64], higher_better: bool| best_of(scores, higher_better)
);

identity_algorithm!(
    /// Consensus score is the worst search-engine score of each sequence
    /// (source `ConsensusIDAlgorithmWorst`).
    ///
    /// Only identical sequences are scored together. Score orientations of all
    /// inputs must agree, and the score types should.
    ///
    /// # Errors
    ///
    /// Besides those of [`ConsensusIDAlgorithm::apply_with_search_engines`],
    /// [`Error::InvalidValue`] when the score orientations disagree.
    ConsensusIDAlgorithmWorst,
    "ConsensusIDAlgorithmWorst",
    |scores: &[f64], higher_better: bool| best_of(scores, !higher_better)
);

/// `*max_element` when `highest`, otherwise `*min_element`; `scores` is never
/// empty and holds finite values.
fn best_of(scores: &[f64], highest: bool) -> f64 {
    scores[1..].iter().fold(scores[0], |best, &score| {
        if (highest && score > best) || (!highest && score < best) {
            score
        } else {
            best
        }
    })
}

/// Consensus from search ranks rather than scores (source
/// `ConsensusIDAlgorithmRanks`).
///
/// Each sequence gets a score from each ID run based on its rank there: 0 for
/// the best hit, 1 for the second, up to `N = considered_hits - 1`; a sequence
/// a run did not report scores `considered_hits` from that run. The scores of
/// a sequence are averaged and normalized to the range 0 (exclusive, worst) to
/// 1 (inclusive, best). With `filter:considered_hits` 0, the largest number of
/// hits of any run takes its place.
///
/// The output score type is `ConsensusID_ranks`, higher is better. Unlike the
/// other identity algorithms, score orientations are not checked, as in the
/// source, since only ranks are used.
///
/// # Errors
///
/// Besides those of [`ConsensusIDAlgorithm::apply_with_search_engines`],
/// [`Error::InvalidValue`] when a sequence occurs in more identifications than
/// `number_of_runs`. The source subtracts unsigned sizes there, which wraps
/// and yields a meaningless score.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ConsensusIDAlgorithmRanks {
    /// The shared `filter:` parameters.
    pub filter: ConsensusIDFilter,
}

impl ConsensusIDAlgorithmRanks {
    handler_param_api!("ConsensusIDAlgorithmRanks");

    /// The parameters in the layout of `defaults`, carrying the current values
    /// (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `considered_hits` does not fit the
    /// source 32-bit parameter.
    pub fn to_param(&self) -> Result<Param> {
        let mut param = Param::new();
        self.filter.write(&mut param)?;
        Ok(param)
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        Ok(Self {
            filter: ConsensusIDFilter::read(param)?,
        })
    }
}

impl ConsensusIDAlgorithm for ConsensusIDAlgorithmRanks {
    fn filter(&self) -> &ConsensusIDFilter {
        &self.filter
    }

    fn apply_with_search_engines(
        &self,
        ids: &mut Vec<PeptideIdentification>,
        se_info: &BTreeMap<String, String>,
        number_of_runs: usize,
    ) -> Result<Vec<String>> {
        run(&self.filter, Method::Ranks, ids, se_info, number_of_runs)
    }
}

/// Similarity-based consensus from posterior error probabilities and
/// substitution-matrix alignment of the sequences (source
/// `ConsensusIDAlgorithmPEPMatrix`).
///
/// Similarity-based algorithms need posterior error probabilities as scores
/// (score type `Posterior Error Probability`, `pep` or `MS:1001493`) to
/// combine scores and similarities into a consensus score; see Nahnsen et
/// al., "Probabilistic consensus scoring improves tandem mass spectrometry
/// peptide identification", J. Proteome Res. 2011, PMID 21644507.
///
/// For each sequence, first seen with PEP `p` in some run, the best match in
/// every other run is the hit of highest similarity `s_i` (ties: lower PEP
/// `p_i`), and the consensus score is `(p + sum s_i * p_i) / (1 + sum s_i)^2`.
/// The support is the average best-match similarity, over the other runs with
/// hits, or over `number_of_runs - 1` with `filter:count_empty`.
///
/// The similarity of two sequences is the [`NeedlemanWunsch`] score of their
/// unmodified forms, 0 when negative, else divided by the smaller of the two
/// self-alignment scores; identical unmodified sequences have similarity 1.
///
/// # Errors
///
/// Besides those of [`ConsensusIDAlgorithm::apply_with_search_engines`],
/// [`Error::InvalidValue`] when a score type is not a PEP, the input needs
/// more than [`MAX_SIMILARITY_EVALUATIONS`] similarity evaluations, an
/// alignment fails (see [`NeedlemanWunsch::align`]), or a self-alignment
/// score of zero makes a similarity non-finite (the source divides by zero).
///
/// # Notes
///
/// The source dereferences `max_element` of an empty range, which is
/// undefined behaviour, when another identification has no hits. This skips
/// such identifications: they contribute no best match, as an empty run
/// supports nothing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConsensusIDAlgorithmPEPMatrix {
    /// The shared `filter:` parameters.
    pub filter: ConsensusIDFilter,
    /// `matrix`: the substitution matrix for alignment-based similarity.
    pub matrix: ScoringMatrix,
    /// `penalty`: the alignment gap penalty, the same for gap opening and
    /// extension; at least 1.
    pub penalty: i32,
}

impl Default for ConsensusIDAlgorithmPEPMatrix {
    fn default() -> Self {
        Self {
            filter: ConsensusIDFilter::default(),
            matrix: ScoringMatrix::PAM30MS,
            penalty: 5,
        }
    }
}

impl ConsensusIDAlgorithmPEPMatrix {
    handler_param_api!("ConsensusIDAlgorithmPEPMatrix");

    /// The parameters in the layout of `defaults`, carrying the current values
    /// (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `considered_hits` does not fit the
    /// source 32-bit parameter.
    pub fn to_param(&self) -> Result<Param> {
        let mut param = Param::new();
        self.filter.write(&mut param)?;
        param.set_value(
            "matrix",
            ParamValue::String(self.matrix.name().into()),
            "Substitution matrix to use for alignment-based similarity scoring",
            &[],
        )?;
        param.set_valid_strings("matrix", &ScoringMatrix::NAMES.map(ToOwned::to_owned))?;
        param.set_value(
            "penalty",
            ParamValue::Integer(i64::from(self.penalty)),
            "Alignment gap penalty (the same value is used for gap opening and extension)",
            &[],
        )?;
        param.set_min_int("penalty", 1)?;
        Ok(param)
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        let penalty = param.value("penalty")?.to_i64()?;
        if penalty <= 0 {
            return Err(bad("Gap penalty should be positive"));
        }
        Ok(Self {
            filter: ConsensusIDFilter::read(param)?,
            matrix: ScoringMatrix::from_name(param.value("matrix")?.as_str()?)?,
            penalty: i32::try_from(penalty)
                .map_err(|_| bad("'penalty' does not fit the source 32-bit parameter"))?,
        })
    }
}

impl ConsensusIDAlgorithm for ConsensusIDAlgorithmPEPMatrix {
    fn filter(&self) -> &ConsensusIDFilter {
        &self.filter
    }

    fn apply_with_search_engines(
        &self,
        ids: &mut Vec<PeptideIdentification>,
        se_info: &BTreeMap<String, String>,
        number_of_runs: usize,
    ) -> Result<Vec<String>> {
        if self.penalty <= 0 {
            return Err(bad("Gap penalty should be positive"));
        }
        let aligner = NeedlemanWunsch::new(self.matrix, self.penalty);
        run(
            &self.filter,
            Method::PepMatrix(aligner),
            ids,
            se_info,
            number_of_runs,
        )
    }
}

/// Similarity-based consensus from posterior error probabilities and shared
/// b/y fragment ions (source `ConsensusIDAlgorithmPEPIons`).
///
/// The consensus score and support follow [`ConsensusIDAlgorithmPEPMatrix`];
/// only the similarity differs. It is the shared peak count of the two b/y
/// ion ladders: every ion of the first sequence is matched to the closest ion
/// of the second within `mass_tolerance` Da, each ion of the second counts
/// once, and the count, when at least `min_shared`, is divided by the smaller
/// ladder length. Identical sequences have similarity 1.
///
/// The ladders reproduce the source arithmetic, which is not a textbook ion
/// ladder: the first b and y ions are full (water-including) one-residue
/// masses, and each further ion adds a full residue mass, so every step
/// carries an extra water. The measure is not symmetric; the pair is put in
/// source `AASequence` order first. The ratio is computed in single precision,
/// as the source's `float` division, which is visible in the support values.
///
/// # Errors
///
/// Besides those of [`ConsensusIDAlgorithmPEPMatrix`] (other than alignment),
/// [`Error::InvalidValue`] for an empty sequence (source out-of-bounds write)
/// or a residue without a known mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConsensusIDAlgorithmPEPIons {
    /// The shared `filter:` parameters.
    pub filter: ConsensusIDFilter,
    /// `mass_tolerance`: maximum difference between fragment masses, in Da,
    /// for fragments to be considered shared; at least 0.
    pub mass_tolerance: f64,
    /// `min_shared`: the minimal number of shared fragments needed to evaluate
    /// the similarity from the shared peak count; at least 1.
    pub min_shared: usize,
}

impl Default for ConsensusIDAlgorithmPEPIons {
    fn default() -> Self {
        Self {
            filter: ConsensusIDFilter::default(),
            mass_tolerance: 0.5,
            min_shared: 2,
        }
    }
}

impl ConsensusIDAlgorithmPEPIons {
    handler_param_api!("ConsensusIDAlgorithmPEPIons");

    /// The parameters in the layout of `defaults`, carrying the current values
    /// (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `considered_hits` or `min_shared`
    /// does not fit the source 32-bit parameter.
    pub fn to_param(&self) -> Result<Param> {
        let mut param = Param::new();
        self.filter.write(&mut param)?;
        param.set_value(
            "mass_tolerance",
            ParamValue::Float(self.mass_tolerance),
            "Maximum difference between fragment masses (in Da) for fragments to be considered 'shared' between peptides .",
            &[],
        )?;
        param.set_min_float("mass_tolerance", 0.0)?;
        let min_shared = i32::try_from(self.min_shared)
            .map_err(|_| bad("'min_shared' does not fit the source 32-bit parameter"))?;
        param.set_value(
            "min_shared",
            ParamValue::Integer(i64::from(min_shared)),
            "The minimal number of 'shared' fragments (between two suggested peptides) that is necessary to evaluate the similarity based on shared peak count (SPC).",
            &[],
        )?;
        param.set_min_int("min_shared", 1)?;
        Ok(param)
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        Ok(Self {
            filter: ConsensusIDFilter::read(param)?,
            mass_tolerance: param.value("mass_tolerance")?.to_f64()?,
            min_shared: usize::try_from(param.value("min_shared")?.to_i64()?)
                .map_err(|_| bad("'min_shared' must not be negative"))?,
        })
    }
}

impl ConsensusIDAlgorithm for ConsensusIDAlgorithmPEPIons {
    fn filter(&self) -> &ConsensusIDFilter {
        &self.filter
    }

    fn apply_with_search_engines(
        &self,
        ids: &mut Vec<PeptideIdentification>,
        se_info: &BTreeMap<String, String>,
        number_of_runs: usize,
    ) -> Result<Vec<String>> {
        if !self.mass_tolerance.is_finite() || self.mass_tolerance < 0.0 {
            return Err(bad("'mass_tolerance' must be finite and non-negative"));
        }
        run(
            &self.filter,
            Method::PepIons {
                mass_tolerance: self.mass_tolerance,
                min_shared: self.min_shared,
            },
            ids,
            se_info,
            number_of_runs,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(text: &str) -> AASequence {
        AASequence::parse(text).unwrap()
    }

    #[test]
    fn source_order_is_length_first() {
        assert_eq!(source_sequence_cmp(&seq("K"), &seq("AA")), Ordering::Less);
        assert_eq!(
            source_sequence_cmp(&seq("AB"), &seq("AA")),
            Ordering::Greater
        );
        assert_eq!(
            source_sequence_cmp(&seq("PEPM"), &seq("PEPM(Oxidation)")),
            Ordering::Less
        );
    }

    #[test]
    fn charge_conflict_is_an_error() {
        let mut charge = 0;
        compare_charge_states(&mut charge, 2, &seq("PEPTIDE")).unwrap();
        assert_eq!(charge, 2);
        compare_charge_states(&mut charge, 0, &seq("PEPTIDE")).unwrap();
        assert!(compare_charge_states(&mut charge, 3, &seq("PEPTIDE")).is_err());
    }

    #[test]
    fn best_and_worst_pick_by_orientation() {
        assert_eq!(best_of(&[0.2, 0.1, 0.3], true), 0.3);
        assert_eq!(best_of(&[0.2, 0.1, 0.3], false), 0.1);
    }
}
