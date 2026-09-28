// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Target/decoy FDR, q-value and posterior-probability score calculations.
//! Legacy and Basic reproduce different source formulas; see FDR_SUPPORT.md.

use crate::identification::{PeptideIdentification, ProteinIdentification, TargetDecoyType};
use crate::kernel::ConsensusMap;
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::BTreeMap;

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("FDR score/arithmetic must be finite"))
    }
}
fn target(label: TargetDecoyType) -> Result<bool> {
    match label {
        TargetDecoyType::Target | TargetDecoyType::TargetAndDecoy => Ok(true),
        TargetDecoyType::Decoy => Ok(false),
        TargetDecoyType::Unknown => Err(bad("FDR requires target_decoy annotations")),
    }
}
fn better(a: f64, b: f64, higher: bool) -> bool {
    if higher { a > b } else { a < b }
}
// Finite scores only; partial_cmp deliberately equates positive and negative zero.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Score(f64);
impl Eq for Score {}
impl PartialOrd for Score {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Score {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0
            .partial_cmp(&other.0)
            .expect("validated finite score")
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum FdrOutput {
    #[default]
    QValue,
    Fdr,
}
impl FdrOutput {
    fn name(self) -> &'static str {
        match self {
            Self::QValue => "q-value",
            Self::Fdr => "FDR",
        }
    }
    fn peptide_name(self) -> &'static str {
        match self {
            Self::QValue => "peptide q-value",
            Self::Fdr => "peptide FDR",
        }
    }
}
/// Fractional target labels are supported by the source Basic calculation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScoreLabel {
    pub score: f64,
    pub target_fraction: f64,
}
impl ScoreLabel {
    pub fn new(score: f64, is_target: bool) -> Self {
        Self {
            score,
            target_fraction: if is_target { 1. } else { 0. },
        }
    }
}
#[derive(Clone, Debug)]
pub struct ScoreCurve {
    values: BTreeMap<Score, f64>,
    higher_score_better: bool,
    exact: bool,
}
impl ScoreCurve {
    /// Entries are in ascending original-score order.
    pub fn entries(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        self.values.iter().map(|(s, v)| (s.0, *v))
    }
    pub fn len(&self) -> usize {
        self.values.len()
    }
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
    /// Legacy requires a score present in its inputs. Basic uses source cutoff
    /// lookup: next greater/equal score for higher-better, next smaller/equal
    /// score for lower-better (clamped at the smallest stored score).
    pub fn value(&self, score: f64) -> Result<f64> {
        finite(score)?;
        let key = Score(score);
        let result = if self.exact {
            self.values.get(&key)
        } else if self.higher_score_better {
            self.values.range(key..).next().map(|(_, v)| v)
        } else {
            self.values
                .range(..=key)
                .next_back()
                .or_else(|| self.values.first_key_value())
                .map(|(_, v)| v)
        };
        result
            .copied()
            .ok_or_else(|| bad("score lies outside the FDR lookup curve"))
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FalseDiscoveryRate {
    pub output: FdrOutput,
    pub use_all_hits: bool,
    pub split_charge_variants: bool,
    pub treat_runs_separately: bool,
    pub add_decoy_peptides: bool,
    pub add_decoy_proteins: bool,
    /// Actual Basic source formula: true (D+1)/(T+1), false (D+1)/(T+D+1).
    /// Legacy calculations ignore this option.
    pub conservative: bool,
    pub max_hits: usize,
    pub max_records: usize,
    pub max_groups: usize,
}
impl Default for FalseDiscoveryRate {
    fn default() -> Self {
        Self {
            output: FdrOutput::QValue,
            use_all_hits: false,
            split_charge_variants: false,
            treat_runs_separately: false,
            add_decoy_peptides: false,
            add_decoy_proteins: false,
            conservative: true,
            max_hits: 1_000_000,
            max_records: 1_000_000,
            max_groups: 100_000,
        }
    }
}
impl FalseDiscoveryRate {
    fn limits(&self, records: usize, hits: usize, groups: usize) -> Result<()> {
        if self.max_records == 0
            || self.max_hits == 0
            || self.max_groups == 0
            || records > self.max_records
            || hits > self.max_hits
            || groups > self.max_groups
        {
            return Err(bad("FDR resource limit exceeded or zero"));
        }
        Ok(())
    }
    /// Source legacy calculation: inclusive threshold D/T, cumulative minimum
    /// for q-values, then its original decoy-to-target score assignment.
    pub fn calculate_legacy(
        &self,
        targets: &[f64],
        decoys: &[f64],
        higher: bool,
    ) -> Result<ScoreCurve> {
        self.limits(
            0,
            targets
                .len()
                .checked_add(decoys.len())
                .ok_or_else(|| bad("FDR count overflow"))?,
            0,
        )?;
        for &score in targets.iter().chain(decoys) {
            finite(score)?;
        }
        let mut targets = targets.to_vec();
        let mut decoys = decoys.to_vec();
        let order = |a: &f64, b: &f64| {
            if higher {
                b.partial_cmp(a).unwrap()
            } else {
                a.partial_cmp(b).unwrap()
            }
        };
        targets.sort_by(order);
        decoys.sort_by(order);
        let mut values = BTreeMap::new();
        let (mut i, mut d) = (0, 0);
        while i < targets.len() {
            let score = targets[i];
            let mut end = i + 1;
            while end < targets.len() && targets[end] == score {
                end += 1;
            }
            while d < decoys.len() && !better(score, decoys[d], higher) {
                d += 1;
            }
            values.insert(Score(score), d as f64 / end as f64);
            i = end;
        }
        if self.output == FdrOutput::QValue {
            cumulative_minimum(&mut values, higher);
            targets.reverse();
        }
        // Faithful source behavior: the raw-FDR target order remains best first.
        // Consequently its decoy lookup is not a general nearest-neighbor lookup.
        for score in decoys {
            let k = if self.output == FdrOutput::QValue {
                targets.partition_point(|&t| !better(t, score, higher))
            } else if targets.first().is_some_and(|&t| !better(t, score, higher)) {
                targets.len()
            } else {
                0
            };
            let chosen = if targets.is_empty() {
                None
            } else if k == 0 {
                Some(targets[0])
            } else if k == targets.len() {
                targets.last().copied()
            } else if (targets[k] - score).abs() < (targets[k - 1] - score).abs() {
                Some(targets[k])
            } else {
                Some(targets[k - 1])
            };
            let value = chosen.map_or(1., |s| values[&Score(s)]);
            values.insert(Score(score), value);
        }
        Ok(ScoreCurve {
            values,
            higher_score_better: higher,
            exact: true,
        })
    }
    /// Source Basic uses absolute 1e-12 tie groups and adds one to both its
    /// decoy count and denominator; raw FDR may exceed one.
    pub fn calculate_basic(&self, observations: &[ScoreLabel], higher: bool) -> Result<ScoreCurve> {
        self.limits(0, observations.len(), 0)?;
        for observation in observations {
            finite(observation.score)?;
            if !observation.target_fraction.is_finite()
                || !(0. ..=1.).contains(&observation.target_fraction)
            {
                return Err(bad("target fraction must be in 0..=1"));
            }
        }
        let mut observations = observations.to_vec();
        observations.sort_by(|a, b| {
            let order = a
                .score
                .partial_cmp(&b.score)
                .unwrap()
                .then(a.target_fraction.partial_cmp(&b.target_fraction).unwrap());
            if higher { order.reverse() } else { order }
        });
        let mut values = BTreeMap::new();
        let mut decoys = 0.;
        let mut i = 0;
        while i < observations.len() {
            let score = observations[i].score;
            let mut end = i;
            while end < observations.len() && (observations[end].score - score).abs() <= 1e-12 {
                decoys += 1. - observations[end].target_fraction;
                end += 1;
            }
            let denominator = end as f64 + 1. - if self.conservative { decoys } else { 0. };
            values.insert(Score(score), finite((decoys + 1.) / denominator)?);
            i = end;
        }
        if self.output == FdrOutput::QValue {
            cumulative_minimum(&mut values, higher);
        }
        Ok(ScoreCurve {
            values,
            higher_score_better: higher,
            exact: false,
        })
    }
    /// Intended estimated-q calculation: running mean PEP, or one minus running
    /// mean posterior probability, grouping exact ties at their final endpoint.
    /// Corrects the source's write into a reserved but unsized vector.
    pub fn calculate_estimated(&self, probabilities: &[f64], higher: bool) -> Result<ScoreCurve> {
        self.limits(0, probabilities.len(), 0)?;
        if probabilities
            .iter()
            .any(|v| !v.is_finite() || !(0. ..=1.).contains(v))
        {
            return Err(bad("posterior probabilities must be in 0..=1"));
        }
        let mut scores = probabilities.to_vec();
        scores.sort_by(|a, b| {
            if higher {
                b.partial_cmp(a).unwrap()
            } else {
                a.partial_cmp(b).unwrap()
            }
        });
        let mut values = BTreeMap::new();
        let mut sum = 0.;
        for (i, &score) in scores.iter().enumerate() {
            sum += if higher { 1. - score } else { score };
            values.insert(Score(score), sum / (i + 1) as f64);
        }
        Ok(ScoreCurve {
            values,
            higher_score_better: higher,
            exact: true,
        })
    }
}
fn cumulative_minimum(values: &mut BTreeMap<Score, f64>, higher: bool) {
    let mut minimum: f64 = 1.;
    if higher {
        for value in values.values_mut() {
            minimum = minimum.min(*value);
            *value = minimum;
        }
    } else {
        for value in values.values_mut().rev() {
            minimum = minimum.min(*value);
            *value = minimum;
        }
    }
}
fn old_score(metadata: &mut crate::metadata::MetaInfo, key: &str, value: f64) -> Result<()> {
    if key == "target_decoy" {
        return Err(bad(
            "FDR score backup would overwrite the target_decoy label",
        ));
    }
    metadata.insert(key.into(), MetaValue::try_from(value)?);
    Ok(())
}
fn count_hits(mut lengths: impl Iterator<Item = usize>) -> Result<usize> {
    lengths.try_fold(0usize, |n, v| {
        n.checked_add(v).ok_or_else(|| bad("FDR count overflow"))
    })
}
fn validate_peptides(options: &FalseDiscoveryRate, ids: &[PeptideIdentification]) -> Result<()> {
    options.limits(
        ids.len(),
        count_hits(ids.iter().map(|id| id.hits.len()))?,
        0,
    )?;
    for id in ids {
        id.validate()?;
    }
    Ok(())
}
fn validate_proteins(options: &FalseDiscoveryRate, ids: &[ProteinIdentification]) -> Result<()> {
    options.limits(
        ids.len(),
        count_hits(ids.iter().map(|id| id.hits.len()))?,
        count_hits(ids.iter().map(|id| id.indistinguishable_groups.len()))?,
    )?;
    for id in ids {
        id.validate()?;
    }
    Ok(())
}
fn consistent(current: &mut Option<(String, bool)>, score_type: &str, higher: bool) -> Result<()> {
    if let Some((kind, direction)) = current {
        if kind != score_type || *direction != higher {
            return Err(bad("FDR pool mixes score types or score directions"));
        }
    } else {
        *current = Some((score_type.into(), higher));
    }
    Ok(())
}
#[derive(Default)]
struct PeptidePool {
    kind: Option<(String, bool)>,
    positions: Vec<(usize, usize)>,
    observations: Vec<ScoreLabel>,
}
impl FalseDiscoveryRate {
    /// Concatenated-search source apply: sorts each record, keeps only its best
    /// hit unless use_all_hits, optionally pools by run and charge, removes decoys.
    pub fn apply_peptides(
        &self,
        ids: &mut [PeptideIdentification],
        annotate_peptide_fdr: bool,
    ) -> Result<()> {
        self.apply_peptides_inner(ids, false, annotate_peptide_fdr)
    }
    /// Source Basic counts the first hit unless use_all_hits, and annotates all
    /// hits. Records must already be sorted when the first hit is selected.
    /// Pool restrictions apply to both counting and writing in the native API.
    pub fn apply_basic_peptides(&self, ids: &mut [PeptideIdentification]) -> Result<()> {
        self.apply_peptides_inner(ids, true, false)
    }
    fn apply_peptides_inner(
        &self,
        ids: &mut [PeptideIdentification],
        basic: bool,
        annotate: bool,
    ) -> Result<()> {
        validate_peptides(self, ids)?;
        if ids.is_empty() {
            return Ok(());
        }
        let mut result = ids.to_vec();
        if !basic {
            for id in &mut result {
                id.sort()?;
                if !self.use_all_hits {
                    id.hits.truncate(1);
                }
            }
        }
        let mut pools: BTreeMap<(String, i32), PeptidePool> = BTreeMap::new();
        for (i, id) in result.iter().enumerate() {
            if basic
                && !self.use_all_hits
                && id
                    .hits
                    .windows(2)
                    .any(|w| better(w[1].score, w[0].score, id.higher_score_better))
            {
                return Err(bad(
                    "Basic FDR requires sorted peptide hits when selecting first hit",
                ));
            }
            for (j, hit) in id.hits.iter().enumerate() {
                let key = (
                    if self.treat_runs_separately {
                        id.identifier.clone()
                    } else {
                        String::new()
                    },
                    if self.split_charge_variants {
                        hit.charge
                    } else {
                        0
                    },
                );
                if !pools.contains_key(&key) && pools.len() >= self.max_groups {
                    return Err(bad("FDR pool limit exceeded"));
                }
                let pool = pools.entry(key).or_default();
                consistent(&mut pool.kind, &id.score_type, id.higher_score_better)?;
                let is_target = target(hit.target_decoy_type()?)?;
                pool.positions.push((i, j));
                if !basic || self.use_all_hits || j == 0 {
                    pool.observations
                        .push(ScoreLabel::new(hit.score, is_target));
                }
            }
        }
        self.limits(result.len(), 0, pools.len())?;
        let mut keep: Vec<Vec<bool>> = result.iter().map(|id| vec![true; id.hits.len()]).collect();
        for pool in pools.values() {
            let higher = pool.kind.as_ref().unwrap().1;
            if basic && pool.observations.is_empty() {
                return Err(bad("Basic FDR pool contains no selected scores"));
            }
            let targets: Vec<_> = pool
                .observations
                .iter()
                .filter(|p| p.target_fraction > 0.)
                .map(|p| p.score)
                .collect();
            let decoys: Vec<_> = pool
                .observations
                .iter()
                .filter(|p| p.target_fraction == 0.)
                .map(|p| p.score)
                .collect();
            let degenerate = !basic && (targets.is_empty() || decoys.is_empty());
            let curve = if basic {
                self.calculate_basic(&pool.observations, higher)?
            } else {
                self.calculate_legacy(&targets, &decoys, higher)?
            };
            let mut peptide_values: BTreeMap<(String, bool), f64> = BTreeMap::new();
            if annotate && !degenerate {
                for &(i, j) in &pool.positions {
                    let hit = &result[i].hits[j];
                    let label = target(hit.target_decoy_type()?)?;
                    let entry = peptide_values
                        .entry((hit.sequence.as_str().into(), label))
                        .or_insert(hit.score);
                    if better(hit.score, *entry, higher) {
                        *entry = hit.score;
                    }
                }
                let targets: Vec<_> = peptide_values
                    .iter()
                    .filter(|((_, t), _)| *t)
                    .map(|(_, s)| *s)
                    .collect();
                let decoys: Vec<_> = peptide_values
                    .iter()
                    .filter(|((_, t), _)| !*t)
                    .map(|(_, s)| *s)
                    .collect();
                let peptide_curve = self.calculate_legacy(&targets, &decoys, higher)?;
                for score in peptide_values.values_mut() {
                    *score = peptide_curve.value(*score)?;
                }
            }
            for &(i, j) in &pool.positions {
                let key = format!("{}_score", result[i].score_type);
                let hit = &mut result[i].hits[j];
                let label = target(hit.target_decoy_type()?)?;
                if !label && (!self.add_decoy_peptides || degenerate) {
                    keep[i][j] = false;
                    continue;
                }
                old_score(&mut hit.metadata, &key, hit.score)?;
                hit.score = if degenerate {
                    0.
                } else {
                    curve.value(hit.score)?
                };
                if annotate && !degenerate {
                    old_score(
                        &mut hit.metadata,
                        self.output.peptide_name(),
                        peptide_values[&(hit.sequence.as_str().into(), label)],
                    )?;
                }
            }
        }
        for (id, keep) in result.iter_mut().zip(keep) {
            let mut i = 0;
            id.hits.retain(|_| {
                let value = keep[i];
                i += 1;
                value
            });
            id.score_type = self.output.name().into();
            id.higher_score_better = false;
            if !basic {
                id.sort()?;
            }
        }
        ids.clone_from_slice(&result);
        Ok(())
    }
    /// Best first hit per unmodified sequence, with targets preferred on equal
    /// scores. Updates only each record's first hit; requires one hit per record
    /// to avoid the source producing records with mixed score units.
    pub fn apply_basic_peptide_level(&self, ids: &mut [PeptideIdentification]) -> Result<()> {
        validate_peptides(self, ids)?;
        if ids.is_empty() {
            return Ok(());
        }
        if self.split_charge_variants || self.treat_runs_separately {
            return Err(bad("peptide-level FDR does not split runs or charges"));
        }
        let mut kind = None;
        let mut representatives: BTreeMap<String, ScoreLabel> = BTreeMap::new();
        for id in ids.iter() {
            if id.hits.len() > 1 {
                return Err(bad("peptide-level FDR requires at most one hit per record"));
            }
            if let Some(hit) = id.hits.first() {
                consistent(&mut kind, &id.score_type, id.higher_score_better)?;
                let candidate = ScoreLabel::new(hit.score, target(hit.target_decoy_type()?)?);
                let entry = representatives
                    .entry(hit.sequence.as_str().into())
                    .or_insert(candidate);
                if better(candidate.score, entry.score, id.higher_score_better)
                    || candidate.score == entry.score
                        && candidate.target_fraction > entry.target_fraction
                {
                    *entry = candidate;
                }
            }
        }
        let Some((_, higher)) = kind else {
            return Ok(());
        };
        let curve = self.calculate_basic(
            &representatives.values().copied().collect::<Vec<_>>(),
            higher,
        )?;
        let mut result = ids.to_vec();
        for id in &mut result {
            if let Some(hit) = id.hits.first_mut() {
                if !self.add_decoy_peptides && !target(hit.target_decoy_type()?)? {
                    id.hits.clear();
                    continue;
                }
                let representative = representatives[hit.sequence.as_str()];
                old_score(&mut hit.metadata, &id.score_type, hit.score)?;
                hit.score = curve.value(representative.score)?;
                id.score_type = self.output.peptide_name().into();
                id.higher_score_better = false;
            }
        }
        ids.clone_from_slice(&result);
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DecoyAffix<'a> {
    Prefix(&'a str),
    Suffix(&'a str),
}
impl DecoyAffix<'_> {
    fn validate(self) -> Result<()> {
        let text = match self {
            Self::Prefix(s) | Self::Suffix(s) => s,
        };
        if text.is_empty() {
            Err(bad("picked FDR requires an explicit nonempty decoy affix"))
        } else {
            Ok(())
        }
    }
    fn strip(self, accession: &str) -> Option<&str> {
        match self {
            Self::Prefix(s) => accession.strip_prefix(s),
            Self::Suffix(s) => accession.strip_suffix(s),
        }
    }
}
impl FalseDiscoveryRate {
    /// Separate searches: counts all hits regardless of use_all_hits. Decoy IDs
    /// remain unchanged unless add_decoy_peptides is true. No labels are required.
    pub fn apply_separate_peptides(
        &self,
        targets: &mut [PeptideIdentification],
        decoys: &mut [PeptideIdentification],
    ) -> Result<()> {
        self.limits(count_hits([targets.len(), decoys.len()].into_iter())?, 0, 0)?;
        validate_peptides(self, targets)?;
        validate_peptides(self, decoys)?;
        if targets.is_empty() || decoys.is_empty() {
            return Ok(());
        }
        let mut kind = None;
        for id in targets.iter().chain(decoys.iter()) {
            consistent(&mut kind, &id.score_type, id.higher_score_better)?;
        }
        let higher = targets[0].higher_score_better;
        let target_scores: Vec<_> = targets
            .iter()
            .flat_map(|id| id.hits.iter().map(|h| h.score))
            .collect();
        let decoy_scores: Vec<_> = decoys
            .iter()
            .flat_map(|id| id.hits.iter().map(|h| h.score))
            .collect();
        let curve = self.calculate_legacy(&target_scores, &decoy_scores, higher)?;
        let mut forward = targets.to_vec();
        let mut reverse = if self.add_decoy_peptides {
            decoys.to_vec()
        } else {
            Vec::new()
        };
        for id in forward.iter_mut().chain(reverse.iter_mut()) {
            let key = format!("{}_score", id.score_type);
            for hit in &mut id.hits {
                old_score(&mut hit.metadata, &key, hit.score)?;
                hit.score = curve.value(hit.score)?;
            }
            id.score_type = self.output.name().into();
            id.higher_score_better = false;
        }
        targets.clone_from_slice(&forward);
        if self.add_decoy_peptides {
            decoys.clone_from_slice(&reverse);
        }
        Ok(())
    }
    /// Source concatenated protein apply pools all runs and updates protein hits;
    /// protein-group scores are unchanged. No hit sorting/rank reassignment occurs.
    pub fn apply_proteins(&self, ids: &mut [ProteinIdentification]) -> Result<()> {
        validate_proteins(self, ids)?;
        if ids.is_empty() {
            return Ok(());
        }
        let mut kind = None;
        let mut targets = Vec::new();
        let mut decoys = Vec::new();
        for id in ids.iter() {
            consistent(&mut kind, &id.score_type, id.higher_score_better)?;
            for hit in &id.hits {
                if target(hit.target_decoy_type()?)? {
                    targets.push(hit.score);
                } else {
                    decoys.push(hit.score);
                }
            }
        }
        let curve = self.calculate_legacy(&targets, &decoys, ids[0].higher_score_better)?;
        let mut result = ids.to_vec();
        for id in &mut result {
            replace_protein_scores(id, &curve, self.output.name(), self.add_decoy_proteins)?;
        }
        ids.clone_from_slice(&result);
        Ok(())
    }
    /// Source separate protein searches update only forward records. Decoy records
    /// and all protein groups remain unchanged, independently of retention options.
    pub fn apply_separate_proteins(
        &self,
        targets: &mut [ProteinIdentification],
        decoys: &[ProteinIdentification],
    ) -> Result<()> {
        self.limits(count_hits([targets.len(), decoys.len()].into_iter())?, 0, 0)?;
        validate_proteins(self, targets)?;
        validate_proteins(self, decoys)?;
        if targets.is_empty() || decoys.is_empty() {
            return Ok(());
        }
        let mut kind = None;
        for id in targets.iter().chain(decoys) {
            consistent(&mut kind, &id.score_type, id.higher_score_better)?;
        }
        let target_scores: Vec<_> = targets
            .iter()
            .flat_map(|id| id.hits.iter().map(|h| h.score))
            .collect();
        let decoy_scores: Vec<_> = decoys
            .iter()
            .flat_map(|id| id.hits.iter().map(|h| h.score))
            .collect();
        let curve = self.calculate_legacy(
            &target_scores,
            &decoy_scores,
            targets[0].higher_score_better,
        )?;
        let mut result = targets.to_vec();
        for id in &mut result {
            replace_protein_scores(id, &curve, self.output.name(), true)?;
        }
        targets.clone_from_slice(&result);
        Ok(())
    }
    /// Basic protein calculation, optionally including indistinguishable groups.
    /// A group is target if any listed accession has a target protein hit.
    /// All groups are retained, even when decoy hits are removed.
    pub fn apply_basic_protein(
        &self,
        id: &mut ProteinIdentification,
        groups_too: bool,
    ) -> Result<()> {
        validate_proteins(self, std::slice::from_ref(id))?;
        if id.hits.is_empty() {
            return Err(bad("no protein scores for Basic FDR"));
        }
        let mut labels = BTreeMap::new();
        let mut observations = Vec::new();
        for hit in &id.hits {
            let is_target = target(hit.target_decoy_type()?)?;
            if labels.insert(hit.accession.as_str(), is_target).is_some() {
                return Err(bad("duplicate protein accession in FDR pool"));
            }
            observations.push(ScoreLabel::new(hit.score, is_target));
        }
        let curve = self.calculate_basic(&observations, id.higher_score_better)?;
        let mut result = id.clone();
        if groups_too && !id.indistinguishable_groups.is_empty() {
            let groups = id
                .indistinguishable_groups
                .iter()
                .map(|group| {
                    let mut is_target = false;
                    for accession in &group.accessions {
                        is_target |= *labels
                            .get(accession.as_str())
                            .ok_or_else(|| bad("protein group references missing accession"))?;
                    }
                    Ok(ScoreLabel::new(group.probability, is_target))
                })
                .collect::<Result<Vec<_>>>()?;
            let group_curve = self.calculate_basic(&groups, id.higher_score_better)?;
            for group in &mut result.indistinguishable_groups {
                group.probability = group_curve
                    .values
                    .range(Score(group.probability)..)
                    .next()
                    .map(|(_, v)| *v)
                    .ok_or_else(|| bad("protein group score outside Basic lookup curve"))?;
            }
        }
        replace_protein_scores(
            &mut result,
            &curve,
            self.output.name(),
            self.add_decoy_proteins,
        )?;
        *id = result;
        Ok(())
    }
    /// Pairs accessions by a supplied affix, using target_decoy labels to identify
    /// decoy hits. Target wins score ties. Losing hits are still scored from the
    /// picked curve; decoy retention follows add_decoy_proteins.
    pub fn apply_picked_protein(
        &self,
        id: &mut ProteinIdentification,
        affix: DecoyAffix<'_>,
        groups_too: bool,
    ) -> Result<()> {
        affix.validate()?;
        validate_proteins(self, std::slice::from_ref(id))?;
        let mut picked: BTreeMap<String, ScoreLabel> = BTreeMap::new();
        let mut seen = BTreeMap::new();
        for hit in &id.hits {
            let label = target(hit.target_decoy_type()?)?;
            if seen.insert(hit.accession.as_str(), ()).is_some() {
                return Err(bad("duplicate protein accession in picked FDR"));
            }
            let key = if label {
                hit.accession.as_str()
            } else {
                affix
                    .strip(&hit.accession)
                    .ok_or_else(|| bad("decoy accession does not match configured affix"))?
            };
            if key.is_empty() {
                return Err(bad("empty accession after removing decoy affix"));
            }
            let candidate = ScoreLabel::new(hit.score, label);
            let entry = picked.entry(key.into()).or_insert(candidate);
            if better(candidate.score, entry.score, id.higher_score_better)
                || candidate.score == entry.score
                    && candidate.target_fraction > entry.target_fraction
            {
                *entry = candidate;
            }
        }
        if picked.is_empty() {
            return Err(bad("no protein scores for picked FDR"));
        }
        let curve = self.calculate_basic(
            &picked.values().copied().collect::<Vec<_>>(),
            id.higher_score_better,
        )?;
        let mut result = id.clone();
        if groups_too && !id.indistinguishable_groups.is_empty() {
            let mut groups = Vec::new();
            for group in &id.indistinguishable_groups {
                let mut decoy_picked = false;
                for accession in &group.accessions {
                    if !seen.contains_key(accession.as_str()) {
                        return Err(bad("picked protein group references missing accession"));
                    }
                    let stripped = affix.strip(accession);
                    let label = picked
                        .get(stripped.unwrap_or(accession))
                        .ok_or_else(|| bad("picked protein group has no accession pair"))?
                        .target_fraction;
                    if stripped.is_none() && label > 0. {
                        groups.push(ScoreLabel::new(group.probability, true));
                        break;
                    } else if stripped.is_some() && label == 0. {
                        decoy_picked = true;
                    }
                }
                // Retain the source's order-sensitive mixed-group contribution:
                // a previously encountered picked decoy remains after a target.
                if decoy_picked {
                    groups.push(ScoreLabel::new(group.probability, false));
                }
            }
            let group_curve = self.calculate_basic(&groups, id.higher_score_better)?;
            for group in &mut result.indistinguishable_groups {
                group.probability = group_curve
                    .values
                    .range(Score(group.probability)..)
                    .next()
                    .map(|(_, v)| *v)
                    .ok_or_else(|| bad("picked protein group outside calculated score curve"))?;
            }
        }
        replace_protein_scores(
            &mut result,
            &curve,
            self.output.name(),
            self.add_decoy_proteins,
        )?;
        *id = result;
        Ok(())
    }
    /// One record at a time, avoiding the source vector overload's implicit
    /// first-record-only behavior. PEP/PP type and direction must agree.
    pub fn apply_estimated_protein(&self, id: &mut ProteinIdentification) -> Result<()> {
        validate_proteins(self, std::slice::from_ref(id))?;
        if !matches!(
            (id.score_type.as_str(), id.higher_score_better),
            ("Posterior Error Probability", false) | ("Posterior Probability", true)
        ) {
            return Err(bad(
                "estimated FDR requires correctly oriented PEP or posterior probability scores",
            ));
        }
        let curve = self.calculate_estimated(
            &id.hits.iter().map(|h| h.score).collect::<Vec<_>>(),
            id.higher_score_better,
        )?;
        let mut result = id.clone();
        replace_protein_scores(
            &mut result,
            &curve,
            "Estimated Q-Values",
            self.add_decoy_proteins,
        )?;
        *id = result;
        Ok(())
    }
    /// Normalized target-decoy ROC area, including a whole tied-score batch at
    /// the cutoff. None means full area. Requires binary labels; no targets is
    /// an error, no decoys returns one, and an empty input returns zero.
    pub fn roc_n(
        &self,
        observations: &[ScoreLabel],
        higher: bool,
        false_positive_cutoff: Option<usize>,
    ) -> Result<f64> {
        self.limits(0, observations.len(), 0)?;
        if false_positive_cutoff == Some(0) {
            return Err(bad("ROC false-positive cutoff must be positive or None"));
        }
        for p in observations {
            finite(p.score)?;
            if p.target_fraction != 0. && p.target_fraction != 1. {
                return Err(bad("ROC requires binary target/decoy labels"));
            }
        }
        if observations.is_empty() {
            return Ok(0.);
        }
        if !observations.iter().any(|p| p.target_fraction == 1.) {
            return Err(bad("ROC area undefined without targets"));
        }
        let mut sorted = observations.to_vec();
        sorted.sort_by(|a, b| {
            if higher {
                b.score.partial_cmp(&a.score).unwrap()
            } else {
                a.score.partial_cmp(&b.score).unwrap()
            }
        });
        let (mut targets, mut decoys, mut previous_targets, mut previous_decoys) =
            (0usize, 0usize, 0usize, 0usize);
        let mut area = 0.;
        let mut i = 0;
        while i < sorted.len() {
            let score = sorted[i].score;
            while i < sorted.len() && sorted[i].score == score {
                if sorted[i].target_fraction == 1. {
                    targets += 1;
                } else {
                    decoys += 1;
                }
                i += 1;
            }
            area += (decoys - previous_decoys) as f64 * (targets + previous_targets) as f64 / 2.;
            if false_positive_cutoff.is_some_and(|limit| decoys >= limit) {
                break;
            }
            previous_targets = targets;
            previous_decoys = decoys;
        }
        if decoys == 0 {
            Ok(1.)
        } else if targets == 0 {
            Err(bad("ROC cutoff reached before any target"))
        } else {
            finite(area / (decoys as f64 * targets as f64))
        }
    }
}
fn replace_protein_scores(
    id: &mut ProteinIdentification,
    curve: &ScoreCurve,
    new_type: &str,
    keep_decoys: bool,
) -> Result<()> {
    let key = format!("{}_score", id.score_type);
    let mut hits = Vec::new();
    for hit in &id.hits {
        if !keep_decoys && !target(hit.target_decoy_type()?)? {
            continue;
        }
        let mut hit = hit.clone();
        old_score(&mut hit.metadata, &key, hit.score)?;
        hit.score = curve.value(hit.score)?;
        hits.push(hit);
    }
    id.hits = hits;
    id.score_type = new_type.into();
    id.higher_score_better = false;
    Ok(())
}

/// The source's `getTDLabel_` after `checkTDAnnotation_`: 1 when the
/// `target_decoy` text starts with `t` (`target`, `target+decoy`), else 0.
///
/// # Errors
///
/// [`Error::MissingInformation`] without a `target_decoy` meta value (the
/// source's `Exception::MissingInformation`).
fn td_label(metadata: &crate::metadata::MetaInfo) -> Result<f64> {
    let value = metadata.get("target_decoy").ok_or_else(|| {
        Error::MissingInformation(
            "Meta value 'target_decoy' does not exist in all hits! Reindex with 'PeptideIndexer'"
                .into(),
        )
    })?;
    Ok(if value.to_string().starts_with('t') {
        1.
    } else {
        0.
    })
}

/// `trapezoidal_area_xEqy`: the area between a segment and the diagonal,
/// split into two triangles where the segment crosses it.
fn trapezoidal_area_x_eq_y(x1: f64, x2: f64, y1: f64, y2: f64) -> f64 {
    let height = x2 - x1;
    let b1 = y1 - x1;
    let b2 = y2 - x2;
    if b1.is_sign_negative() == b2.is_sign_negative() {
        (b1.abs() + b2.abs()) * height / 2.0
    } else {
        (b1 * b1 + b2 * b2) * height / (2.0 * (b1.abs() + b2.abs()))
    }
}

/// `trapezoidal_area`: `|x1 - x2| * ((y1 + y2) / 2)`.
fn trapezoidal_area(x1: f64, x2: f64, y1: f64, y2: f64) -> f64 {
    (x1 - x2).abs() * ((y1 + y2) / 2.0)
}

/// Score/label pairs sorted descending as `std::sort(rbegin, rend)` sorts
/// them: by score, then by label. Equal pairs are interchangeable, so the
/// sort's stability does not matter.
fn sort_pairs_descending(pairs: &mut [(f64, f64)]) {
    pairs.sort_by(|a, b| {
        b.0.partial_cmp(&a.0)
            .unwrap_or(Ordering::Equal)
            .then(b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal))
    });
}

impl FalseDiscoveryRate {
    /// `diffEstimatedEmpirical`: the area between the estimated FDR (the
    /// running mean of `1 - score`, scores being posterior probabilities)
    /// and the empirical target-decoy FDR, over pairs sorted best first,
    /// divided by `min(final estimate, pep_cutoff)` (source default
    /// `pep_cutoff` 1).
    ///
    /// At each score change the empirical FDR is `decoys / targets` (one
    /// while no target was seen) with [`FalseDiscoveryRate::conservative`],
    /// else `decoys / (targets + decoys)`; for the last entry it is always
    /// `decoys / (targets + decoys)`, as in the source. A label counts as a
    /// target when it is non-zero. An empty input returns 1 (the source only
    /// warns).
    pub fn diff_estimated_empirical(&self, scores_labels: &[(f64, f64)], pep_cutoff: f64) -> f64 {
        let Some((last, rest)) = scores_labels.split_last() else {
            return 1.0;
        };
        let mut diff_area = 0.0;
        let (mut est_prev, mut emp_prev) = (0.0, 0.0);
        let mut pep_sum = 0.0;
        let (mut true_pos, mut false_pos) = (0u32, 0u32);
        for (index, pair) in rest.iter().enumerate() {
            if pair.1 != 0.0 {
                true_pos = true_pos.wrapping_add(1);
            } else {
                false_pos = false_pos.wrapping_add(1);
            }
            pep_sum += 1.0 - pair.0;
            if scores_labels[index + 1].0 != pair.0 {
                let est = pep_sum / f64::from(true_pos.wrapping_add(false_pos));
                let emp = if self.conservative {
                    if true_pos == 0 {
                        1.0
                    } else {
                        f64::from(false_pos) / f64::from(true_pos)
                    }
                } else {
                    f64::from(false_pos) / f64::from(true_pos.wrapping_add(false_pos))
                };
                diff_area += trapezoidal_area_x_eq_y(est_prev, est, emp_prev, emp);
                est_prev = est;
                emp_prev = emp;
            }
        }
        if last.1 != 0.0 {
            true_pos = true_pos.wrapping_add(1);
        } else {
            false_pos = false_pos.wrapping_add(1);
        }
        pep_sum += 1.0 - last.0;
        let total = f64::from(true_pos.wrapping_add(false_pos));
        let est = pep_sum / total;
        let emp = f64::from(false_pos) / total;
        diff_area += trapezoidal_area_x_eq_y(est_prev, est, emp_prev, emp);
        // std::min(est, pepCutoff)
        diff_area / if pep_cutoff < est { pep_cutoff } else { est }
    }

    /// `rocN(const ScoreToTgtDecLabelPairs&, Size)`: the target-decoy ROC
    /// area over pairs already sorted best first, up to the score batch in
    /// which the `fp_cutoff`-th decoy appears (source default 50),
    /// normalised by `decoys * targets` (a 32-bit product, as the source's
    /// `UInt`s). No decoys at the end gives 1, an empty input 0.
    ///
    /// Unlike [`FalseDiscoveryRate::roc_n`] nothing is validated or sorted:
    /// a cutoff reached before any target returns the source's division by
    /// zero (infinity or NaN).
    pub fn roc_n_sorted(&self, scores_labels: &[(f64, f64)], fp_cutoff: usize) -> f64 {
        let Some((last, rest)) = scores_labels.split_last() else {
            return 0.0;
        };
        let mut roc = 0.0;
        let (mut tp, mut fp, mut tp_prev, mut fp_prev) = (0u32, 0u32, 0u32, 0u32);
        for (index, pair) in rest.iter().enumerate() {
            if pair.1 != 0.0 {
                tp = tp.wrapping_add(1);
            } else {
                fp = fp.wrapping_add(1);
            }
            if scores_labels[index + 1].0 != pair.0 {
                roc += trapezoidal_area(
                    f64::from(fp),
                    f64::from(fp_prev),
                    f64::from(tp),
                    f64::from(tp_prev),
                );
                if fp as usize >= fp_cutoff {
                    return roc / f64::from(fp.wrapping_mul(tp));
                }
                tp_prev = tp;
                fp_prev = fp;
            }
        }
        if last.1 != 0.0 {
            tp = tp.wrapping_add(1);
        } else {
            fp = fp.wrapping_add(1);
        }
        roc += trapezoidal_area(
            f64::from(fp),
            f64::from(fp_prev),
            f64::from(tp),
            f64::from(tp_prev),
        );
        if fp == 0 {
            return 1.0;
        }
        roc / f64::from(fp.wrapping_mul(tp))
    }

    fn roc_n_from_pairs(
        &self,
        mut pairs: Vec<(f64, f64)>,
        higher: bool,
        fp_cutoff: usize,
    ) -> Result<f64> {
        if pairs.is_empty() {
            return Err(Error::MissingInformation(
                "No scores could be extracted!".into(),
            ));
        }
        if pairs.iter().any(|p| p.0.is_nan()) {
            return Err(bad("ROC scores must not be NaN"));
        }
        if higher {
            sort_pairs_descending(&mut pairs);
        } else {
            pairs.sort_by(|a, b| {
                a.0.partial_cmp(&b.0)
                    .unwrap_or(Ordering::Equal)
                    .then(a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal))
            });
        }
        let cutoff = if fp_cutoff == 0 {
            pairs.len()
        } else {
            fp_cutoff
        };
        Ok(self.roc_n_sorted(&pairs, cutoff))
    }

    /// The score/label pairs of peptide identifications as the source's
    /// `getScores_` collects them: every hit with
    /// [`FalseDiscoveryRate::use_all_hits`], else each record's first hit;
    /// optionally only records of one run identifier.
    fn peptide_pairs<'a>(
        &self,
        ids: impl IntoIterator<Item = &'a PeptideIdentification>,
        identifier: Option<&str>,
        pairs: &mut Vec<(f64, f64)>,
    ) -> Result<()> {
        for id in ids {
            if identifier.is_some_and(|run| run != id.identifier) {
                continue;
            }
            let hits: &[crate::identification::PeptideHit] = if self.use_all_hits {
                &id.hits
            } else {
                // The source reads getHits()[0] unchecked.
                std::slice::from_ref(id.hits.first().ok_or_else(|| {
                    bad("peptide identification without hits (the source reads its first hit unchecked)")
                })?)
            };
            for hit in hits {
                pairs.push((hit.score, td_label(&hit.metadata)?));
            }
        }
        Ok(())
    }

    /// `rocN(const PeptideIdentificationList&, Size)` and the overload with
    /// a run identifier: the ROC area of the peptide hits (first hits unless
    /// [`FalseDiscoveryRate::use_all_hits`]) in the first record's score
    /// direction, up to `fp_cutoff` decoys (0 for the whole curve), computed
    /// by [`FalseDiscoveryRate::roc_n_sorted`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for no records (the source dereferences the
    /// first), a record without hits when first hits are used, or a NaN
    /// score; [`Error::MissingInformation`] for a hit without
    /// `target_decoy` or when no score is extracted (the source's
    /// `MissingInformation`).
    pub fn roc_n_peptides(
        &self,
        ids: &[PeptideIdentification],
        fp_cutoff: usize,
        identifier: Option<&str>,
    ) -> Result<f64> {
        let first = ids
            .first()
            .ok_or_else(|| bad("rocN needs at least one peptide identification"))?;
        let mut pairs = Vec::new();
        self.peptide_pairs(ids, identifier, &mut pairs)?;
        self.roc_n_from_pairs(pairs, first.higher_score_better, fp_cutoff)
    }

    /// `rocN(const ConsensusMap&, Size, bool)` and the overload with a run
    /// identifier: as [`FalseDiscoveryRate::roc_n_peptides`] over the
    /// features' peptide identifications, then (optionally) the unassigned
    /// ones.
    ///
    /// The whole-map overload takes the score direction from the first
    /// feature identification (lower-better when none has one); the source's
    /// identifier overload reads the first feature's first identification
    /// unchecked, which is an error here when absent.
    ///
    /// # Errors
    ///
    /// As [`FalseDiscoveryRate::roc_n_peptides`].
    pub fn roc_n_consensus(
        &self,
        map: &ConsensusMap,
        fp_cutoff: usize,
        identifier: Option<&str>,
        include_unassigned_peptides: bool,
    ) -> Result<f64> {
        let higher = match identifier {
            None => map
                .features
                .iter()
                .find_map(|f| f.base.peptide_identifications.first())
                .is_some_and(|id| id.higher_score_better),
            Some(_) => {
                map.features
                    .first()
                    .and_then(|f| f.base.peptide_identifications.first())
                    .ok_or_else(|| {
                        bad("rocN needs an identification on the first consensus feature")
                    })?
                    .higher_score_better
            }
        };
        let mut pairs = Vec::new();
        for feature in &map.features {
            self.peptide_pairs(
                &feature.base.peptide_identifications,
                identifier,
                &mut pairs,
            )?;
        }
        if include_unassigned_peptides {
            self.peptide_pairs(
                &map.unassigned_peptide_identifications,
                identifier,
                &mut pairs,
            )?;
        }
        self.roc_n_from_pairs(pairs, higher, fp_cutoff)
    }

    fn protein_pairs(id: &ProteinIdentification) -> Result<Vec<(f64, f64)>> {
        if id.score_type != "Posterior Probability" {
            return Err(bad(
                "Proteins in ProteinIdentification do not have a posterior probability assigned. Please run an inference first.",
            ));
        }
        let mut pairs = id
            .hits
            .iter()
            .map(|hit| Ok((hit.score, td_label(&hit.metadata)?)))
            .collect::<Result<Vec<_>>>()?;
        if pairs.iter().any(|p| p.0.is_nan()) {
            return Err(bad("protein scores must not be NaN"));
        }
        sort_pairs_descending(&mut pairs);
        Ok(pairs)
    }

    /// `applyEvaluateProteinIDs(const ProteinIdentification&, ...)`: a quality
    /// score of posterior protein probabilities,
    /// `(1 - diff) * (1 - diff_weight) + rocN * diff_weight`, with `diff` from
    /// [`FalseDiscoveryRate::diff_estimated_empirical`] and `rocN` from
    /// [`FalseDiscoveryRate::roc_n_sorted`] over the hits sorted best first
    /// (source defaults: `pep_cutoff` 1, `fp_cutoff` 50, `diff_weight` 0.2).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] unless the score type is `Posterior
    /// Probability` (the source's `InvalidValue`) or for a NaN score;
    /// [`Error::MissingInformation`] for a hit without `target_decoy`.
    pub fn apply_evaluate_protein_ids(
        &self,
        id: &ProteinIdentification,
        pep_cutoff: f64,
        fp_cutoff: u32,
        diff_weight: f64,
    ) -> Result<f64> {
        let pairs = Self::protein_pairs(id)?;
        Ok(self.evaluate_sorted(&pairs, pep_cutoff, fp_cutoff, diff_weight))
    }

    /// `applyEvaluateProteinIDs(const std::vector<ProteinIdentification>&,
    /// ...)`: only the first run is evaluated (the source warns about the
    /// others), and, unlike the single-run overload, the combination is
    /// `diff * diff_weight + rocN * (1 - diff_weight)`, as the source
    /// computes it (a defect candidate, see `docs/FDR_SUPPORT.md`).
    ///
    /// # Errors
    ///
    /// As [`FalseDiscoveryRate::apply_evaluate_protein_ids`], and
    /// [`Error::InvalidValue`] for no runs (undefined in the source).
    pub fn apply_evaluate_protein_ids_first(
        &self,
        ids: &[ProteinIdentification],
        pep_cutoff: f64,
        fp_cutoff: u32,
        diff_weight: f64,
    ) -> Result<f64> {
        let first = ids
            .first()
            .ok_or_else(|| bad("no protein identification run to evaluate"))?;
        let pairs = Self::protein_pairs(first)?;
        let diff = self.diff_estimated_empirical(&pairs, pep_cutoff);
        let roc = self.roc_n_sorted(&pairs, fp_cutoff as usize);
        Ok(diff * diff_weight + roc * (1.0 - diff_weight))
    }

    /// `applyEvaluateProteinIDs(ScoreToTgtDecLabelPairs&, ...)`: sorts
    /// `scores_labels` best first in place, then evaluates as
    /// [`FalseDiscoveryRate::apply_evaluate_protein_ids`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a NaN score or label.
    pub fn apply_evaluate_score_labels(
        &self,
        scores_labels: &mut [(f64, f64)],
        pep_cutoff: f64,
        fp_cutoff: u32,
        diff_weight: f64,
    ) -> Result<f64> {
        if scores_labels.iter().any(|p| p.0.is_nan() || p.1.is_nan()) {
            return Err(bad("scores and labels must not be NaN"));
        }
        sort_pairs_descending(scores_labels);
        Ok(self.evaluate_sorted(scores_labels, pep_cutoff, fp_cutoff, diff_weight))
    }

    fn evaluate_sorted(
        &self,
        pairs: &[(f64, f64)],
        pep_cutoff: f64,
        fp_cutoff: u32,
        diff_weight: f64,
    ) -> f64 {
        let diff = self.diff_estimated_empirical(pairs, pep_cutoff);
        let auc = self.roc_n_sorted(pairs, fp_cutoff as usize);
        (1.0 - diff) * (1.0 - diff_weight) + auc * diff_weight
    }
}
