// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Resolve competing peptide annotations on features and repeated spectra.
//!
//! Source: ANALYSIS/ID/IDConflictResolverAlgorithm. Resolution preserves the
//! records it moves to the unassigned list and never changes feature geometry.
//! Mutations are committed only after the complete operation succeeds.

use crate::chemistry::AASequence;
use crate::identification::PeptideIdentification;
use crate::kernel::{ConsensusMap, FeatureMap};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ResolutionMethod {
    #[default]
    BestScore,
    /// Keep the winning record and matching sequences in other records.
    /// The winner keeps all its sorted hits, as in the source overload.
    KeepMatching,
    RankAggregation,
}

/// Counts are over the input. Inconsistency counts refer to peptidoform groups,
/// not individual identifications, matching the source counter behavior.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct UnresolvedIdentifications {
    pub removed: usize,
    pub multiply_identified_spectra: usize,
    pub without_spectrum_reference: usize,
    pub inconsistent_score_direction: usize,
    /// Native additional check: different score types cannot be compared.
    pub inconsistent_score_type: usize,
    pub example: String,
}

// The source orders length, terminal modifications and then residues. It uses
// modification pointers for residue ties; names/registry IDs give Rust a
// reproducible order without relying on allocation addresses. Full IDs distinguish
// anonymous mass tags and different specificities sharing a named registry ID.
// The final complete annotation key also preserves distinct caller chemistry
// whose names/accessions happen to coincide.
type ModificationKey = (
    String,
    Option<u32>,
    String,
    String,
    crate::chemistry::SequenceModification,
);
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct SequenceKey {
    length: usize,
    n_terminal: Option<ModificationKey>,
    residues: Vec<(u8, Option<ModificationKey>)>,
    c_terminal: Option<ModificationKey>,
}
impl SequenceKey {
    fn new(sequence: &AASequence) -> Self {
        let modification = |m: &crate::chemistry::SequenceModification| {
            (
                m.name().to_owned(),
                m.record_id(),
                m.full_id().to_owned(),
                m.known()
                    .map(|record| record.accession())
                    .unwrap_or_default(),
                m.clone(),
            )
        };
        Self {
            length: sequence.len(),
            n_terminal: sequence.n_terminal_modification().map(modification),
            residues: sequence
                .as_str()
                .bytes()
                .enumerate()
                .map(|(i, code)| {
                    (
                        code,
                        sequence
                            .residue_modification(i)
                            .expect("valid sequence index")
                            .map(modification),
                    )
                })
                .collect(),
            c_terminal: sequence.c_terminal_modification().map(modification),
        }
    }
}

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn validate(ids: &[PeptideIdentification]) -> Result<()> {
    ids.iter().try_for_each(PeptideIdentification::validate)
}
fn comparable(ids: &[PeptideIdentification]) -> Result<()> {
    if let Some(first) = ids.iter().find(|id| !id.hits.is_empty()) {
        for id in ids.iter().filter(|id| !id.hits.is_empty()) {
            if id.score_type != first.score_type
                || id.higher_score_better != first.higher_score_better
            {
                return Err(bad(
                    "competing identifications require one score type and direction",
                ));
            }
        }
    }
    Ok(())
}
fn annotate(id: &mut PeptideIdentification, feature_id: u64) {
    id.metadata
        .insert("feature_id".into(), feature_id.to_string().into());
}
fn best_record(ids: &[PeptideIdentification]) -> Option<usize> {
    let mut best: Option<usize> = None;
    for (i, id) in ids.iter().enumerate() {
        let Some(hit) = id.hits.first() else { continue };
        if best.is_none_or(|j| {
            let previous = ids[j].hits[0].score;
            if id.higher_score_better {
                hit.score > previous
            } else {
                hit.score < previous
            }
        }) {
            best = Some(i);
        }
    }
    best
}

fn resolve_inner(
    ids: &mut Vec<PeptideIdentification>,
    removed: &mut Vec<PeptideIdentification>,
    feature_id: u64,
    method: ResolutionMethod,
) -> Result<()> {
    if ids.is_empty() {
        return Ok(());
    }
    comparable(ids)?;
    for id in ids.iter_mut() {
        id.sort()?;
        if method != ResolutionMethod::KeepMatching {
            annotate(id, feature_id);
        }
        if method == ResolutionMethod::BestScore {
            id.hits.truncate(1);
        }
    }
    if method == ResolutionMethod::KeepMatching {
        // Ignore empty records when choosing the winner; the C++ overload
        // otherwise dereferences an empty hit vector for lower-better input.
        let Some(best) = best_record(ids) else {
            for id in ids.iter_mut().skip(1) {
                annotate(id, feature_id);
            }
            removed.extend(ids.drain(1..));
            return Ok(());
        };
        ids.swap(0, best);
        let sequence = ids[0].hits[0].sequence.clone();
        let mut keep = vec![ids[0].clone()];
        for mut id in ids.drain(1..) {
            if let Some(hit) = id.hits.iter().find(|hit| hit.sequence == sequence).cloned() {
                id.hits = vec![hit];
                keep.push(id);
            } else {
                annotate(&mut id, feature_id);
                removed.push(id);
            }
        }
        *ids = keep;
        return Ok(());
    }
    let best = if method == ResolutionMethod::RankAggregation {
        rank_winner(ids)?
    } else {
        (best_record(ids).unwrap_or(0), None)
    };
    let mut kept = ids[best.0].clone();
    if let Some(hit_index) = best.1 {
        kept.hits = vec![kept.hits[hit_index].clone()];
    } else {
        kept.hits.truncate(1);
    }
    for (i, mut id) in ids.drain(..).enumerate() {
        if i != best.0 {
            id.hits.truncate(1);
            removed.push(id);
        }
    }
    ids.push(kept);
    Ok(())
}

fn rank_winner(ids: &[PeptideIdentification]) -> Result<(usize, Option<usize>)> {
    let max_hits = ids.iter().map(|id| id.hits.len()).max().unwrap_or(0);
    if max_hits == 0 {
        return Ok((0, None));
    }
    let denominator = max_hits
        .checked_mul(ids.len())
        .ok_or_else(|| bad("rank population overflow"))?;
    let mut ranks: BTreeMap<SequenceKey, (usize, usize)> = BTreeMap::new();
    for id in ids {
        let mut seen = BTreeSet::new();
        for (rank, hit) in id.hits.iter().enumerate() {
            let key = SequenceKey::new(&hit.sequence);
            if seen.insert(key.clone()) {
                let entry = ranks.entry(key).or_default();
                entry.0 = entry
                    .0
                    .checked_add(rank)
                    .ok_or_else(|| bad("rank sum overflow"))?;
                entry.1 += 1;
            }
        }
    }
    let mut best_key = None;
    let mut best_aggregate = -1.0;
    for (key, (rank, count)) in ranks {
        let total = rank + (ids.len() - count) * max_hits; // bounded by denominator
        let aggregate = 1.0 - total as f64 / denominator as f64;
        if aggregate > best_aggregate {
            best_key = Some(key);
            best_aggregate = aggregate;
        }
    }
    let key = best_key.expect("at least one ranked sequence");
    let mut best: Option<(usize, usize)> = None;
    for (i, id) in ids.iter().enumerate() {
        if let Some(j) = id
            .hits
            .iter()
            .position(|hit| SequenceKey::new(&hit.sequence) == key)
        {
            if best.is_none_or(|(old_i, old_j)| {
                let previous = ids[old_i].hits[old_j].score;
                if id.higher_score_better {
                    id.hits[j].score > previous
                } else {
                    id.hits[j].score < previous
                }
            }) {
                best = Some((i, j));
            }
        }
    }
    Ok(best
        .map(|(i, j)| (i, Some(j)))
        .expect("winning sequence exists"))
}

/// Resolve one feature's IDs and append rejected IDs to `unassigned`.
/// Existing unassigned records are preserved verbatim by this direct adapter.
pub fn resolve_identifications(
    ids: &mut Vec<PeptideIdentification>,
    unassigned: &mut Vec<PeptideIdentification>,
    feature_id: u64,
    method: ResolutionMethod,
) -> Result<()> {
    validate(ids)?;
    validate(unassigned)?;
    let mut next = ids.clone();
    let mut removed = Vec::new();
    resolve_inner(&mut next, &mut removed, feature_id, method)?;
    *ids = next;
    unassigned.extend(removed);
    Ok(())
}

macro_rules! maps {
    ($map:ty, $resolve:ident, $between:ident) => {
        /// Resolve top-level features. Existing unassigned IDs receive the
        /// source `feature_id="not mapped"` annotation; subordinates are untouched.
        pub fn $resolve(map: &mut $map, method: ResolutionMethod) -> Result<()> {
            map.validate()?;
            let mut next = map.clone();
            for id in &mut next.unassigned_peptide_identifications {
                id.metadata.insert("feature_id".into(), "not mapped".into());
            }
            for feature in &mut next.features {
                let feature = &mut feature.base;
                feature.metadata.insert("feature_id".into(), feature.unique_id.to_string());
                resolve_inner(&mut feature.peptide_identifications, &mut next.unassigned_peptide_identifications, feature.unique_id, method)?;
            }
            *map = next;
            Ok(())
        }

        /// Keep each (feature charge, modified peptide) on its most intense
        /// feature. Equal intensity keeps the first; removed IDs are unassigned.
        pub fn $between(map: &mut $map) -> Result<()> {
            map.validate()?;
            let mut next = map.clone();
            let mut owners: BTreeMap<(i32, SequenceKey), usize> = BTreeMap::new();
            for i in 0..next.features.len() {
                let feature = &mut next.features[i].base;
                if feature.peptide_identifications.len() > 1 {
                    return Err(bad("between-feature resolution requires at most one identification per feature"));
                }
                let Some(id) = feature.peptide_identifications.first_mut() else { continue };
                id.sort()?;
                let Some(hit) = id.hits.first() else { continue };
                let key = (feature.charge, SequenceKey::new(&hit.sequence));
                if let Some(&old) = owners.get(&key) {
                    let obsolete = if next.features[i].intensity > next.features[old].intensity {
                        owners.insert(key, i);
                        old
                    } else { i };
                    next.unassigned_peptide_identifications.append(&mut next.features[obsolete].peptide_identifications);
                } else { owners.insert(key, i); }
            }
            *map = next;
            Ok(())
        }
    };
}
maps!(FeatureMap, resolve_feature_map, resolve_between_features);
maps!(
    ConsensusMap,
    resolve_consensus_map,
    resolve_between_consensus_features
);

/// Reduce duplicate (spectrum reference, stored top-hit peptidoform, charge)
/// claims within one run. Chimeric alternatives and inconsistent groups remain.
/// Hits are never sorted; survivors retain their original relative order.
pub fn reduce_to_one_per_spectrum(
    ids: &mut Vec<PeptideIdentification>,
) -> Result<UnresolvedIdentifications> {
    validate(ids)?;
    if ids
        .first()
        .is_some_and(|first| ids.iter().any(|id| id.identifier != first.identifier))
    {
        return Err(bad("spectrum reduction requires one identification run"));
    }
    let mut report = UnresolvedIdentifications::default();
    let mut groups: BTreeMap<(String, SequenceKey, i32), Vec<usize>> = BTreeMap::new();
    for (i, id) in ids.iter().enumerate() {
        let Some(hit) = id.hits.first() else { continue };
        let reference = id.spectrum_reference();
        if reference.is_empty() {
            report.without_spectrum_reference += 1;
            continue;
        }
        groups
            .entry((reference, SequenceKey::new(&hit.sequence), hit.charge))
            .or_default()
            .push(i);
    }
    let mut remove = vec![false; ids.len()];
    let mut spectra = BTreeMap::<String, usize>::new();
    for ((reference, _, charge), indices) in groups {
        *spectra.entry(reference.clone()).or_default() += 1;
        let first = &ids[indices[0]];
        let direction = indices
            .iter()
            .any(|&i| ids[i].higher_score_better != first.higher_score_better);
        let score_type = indices
            .iter()
            .any(|&i| ids[i].score_type != first.score_type);
        report.inconsistent_score_direction += usize::from(direction);
        report.inconsistent_score_type += usize::from(score_type);
        if indices.len() == 1 || direction || score_type {
            continue;
        }
        let mut best = indices[0];
        for &i in &indices[1..] {
            if if first.higher_score_better {
                ids[i].hits[0].score > ids[best].hits[0].score
            } else {
                ids[i].hits[0].score < ids[best].hits[0].score
            } {
                best = i;
            }
        }
        for &i in &indices {
            remove[i] = i != best;
        }
        report.removed += indices.len() - 1;
        if report.example.is_empty() {
            report.example = format!("{reference} / {} / charge {charge}", first.hits[0].sequence);
        }
    }
    report.multiply_identified_spectra = spectra.values().filter(|&&n| n > 1).count();
    let mut i = 0;
    ids.retain(|_| {
        let keep = !remove[i];
        i += 1;
        keep
    });
    Ok(report)
}
