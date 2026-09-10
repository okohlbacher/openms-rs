// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Native subset of OpenMS PROCESSING/ID/IDFilter.
//!
//! Filters validate all affected records before mutation. Retained records and
//! ties keep their input order unless score sorting is explicitly required.
//! Hit filters leave empty identifications and stored hit ranks intact. Cleanup
//! of identifications, evidence and protein groups is always explicit.

use crate::chemistry::AASequence;
use crate::comparison::Tolerance;
use crate::identification::{
    PeptideHit, PeptideIdentification, ProteinGroup, ProteinHit, ProteinIdentification,
};
use crate::kernel::{ConsensusMap, FeatureMap};
use crate::metadata::MetaInfo;
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatchAction {
    Keep,
    Remove,
}
impl MatchAction {
    fn keeps(self, matches: bool) -> bool {
        matches == (self == Self::Keep)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DuplicatePolicy {
    /// Compare the entire hit, including score, charge, evidence and metadata.
    #[default]
    Exact,
    /// Compare modified peptide sequences, ignoring every other hit field.
    Sequence,
}

/// Worst-case comparison budget for the source's full-record equality scan.
pub const MAX_EXACT_DUPLICATE_COMPARISONS: usize = 10_000_000;

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn finite(value: f64) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(bad("filter bounds must be finite"))
    }
}
fn bounds<T: PartialOrd>(min: T, max: Option<T>) -> Result<()> {
    if max.is_some_and(|max| max < min) {
        return Err(bad("filter upper bound is smaller than lower bound"));
    }
    Ok(())
}
fn validate_peptides(ids: &[PeptideIdentification]) -> Result<()> {
    ids.iter().try_for_each(PeptideIdentification::validate)
}
fn validate_proteins(ids: &[ProteinIdentification]) -> Result<()> {
    ids.iter().try_for_each(ProteinIdentification::validate)
}
fn order(a: f64, b: f64, higher_better: bool) -> Ordering {
    // Validation excludes NaN; partial_cmp treats signed zero as a tie.
    let cmp = a.partial_cmp(&b).expect("validated finite scores");
    if higher_better { cmp.reverse() } else { cmp }
}
fn has_string_value(metadata: &MetaInfo, key: &str, expected: &str) -> bool {
    metadata
        .get(key)
        .is_some_and(|value| value.unit().is_none() && value.as_str().ok() == Some(expected))
}
fn decoy(metadata: &MetaInfo) -> bool {
    // The source combines both annotations, even when they conflict.
    has_string_value(metadata, "target_decoy", "decoy")
        || has_string_value(metadata, "isDecoy", "true")
}

// Peptide and protein records have the same score/rank contract. Keeping the
// implementation here avoids a new public trait just for these two record types.
macro_rules! hit_filters {
    ($id:ty, $validate:ident, $score:ident, $top:ident, $rank:ident, $decoys:ident, $empty:ident) => {
        /// Keep scores at least as good as the inclusive cutoff, using each ID's direction.
        pub fn $score(ids: &mut [$id], cutoff: f64) -> Result<()> {
            finite(cutoff)?;
            $validate(ids)?;
            for id in ids {
                let higher = id.higher_score_better;
                id.hits.retain(|hit| {
                    if higher {
                        hit.score >= cutoff
                    } else {
                        hit.score <= cutoff
                    }
                });
            }
            Ok(())
        }
        /// Stable score sort followed by truncation; a tie may be split at N.
        pub fn $top(ids: &mut [$id], n: usize) -> Result<()> {
            $validate(ids)?;
            for id in ids {
                let higher = id.higher_score_better;
                id.hits.sort_by(|a, b| order(a.score, b.score, higher));
                id.hits.truncate(n);
            }
            Ok(())
        }
        /// Keep inclusive dense score ranks (1,1,2,...). Stored hit.rank is unchanged.
        /// None means no upper bound; a zero lower bound or reversed range is an error.
        pub fn $rank(ids: &mut [$id], min: usize, max: Option<usize>) -> Result<()> {
            bounds(min, max)?;
            if min == 0 {
                return Err(bad("score ranks start at one"));
            }
            $validate(ids)?;
            for id in ids {
                let higher = id.higher_score_better;
                id.hits.sort_by(|a, b| order(a.score, b.score, higher));
                let mut rank = 0;
                let mut previous = None;
                id.hits.retain(|hit| {
                    if previous != Some(hit.score) {
                        rank += 1;
                        previous = Some(hit.score);
                    }
                    rank >= min && max.is_none_or(|max| rank <= max)
                });
            }
            Ok(())
        }
        /// Remove exact string annotations target_decoy="decoy" OR isDecoy="true".
        pub fn $decoys(ids: &mut [$id]) -> Result<()> {
            $validate(ids)?;
            for id in ids {
                id.hits.retain(|hit| !decoy(&hit.metadata));
            }
            Ok(())
        }
        /// Remove records with no hits, even if they still contain other metadata.
        pub fn $empty(ids: &mut Vec<$id>) -> Result<()> {
            $validate(ids)?;
            ids.retain(|id| !id.hits.is_empty());
            Ok(())
        }
    };
}
hit_filters!(
    PeptideIdentification,
    validate_peptides,
    filter_peptides_by_score,
    keep_n_best_peptide_hits,
    filter_peptides_by_rank,
    remove_decoy_peptide_hits,
    remove_empty_peptide_identifications
);
hit_filters!(
    ProteinIdentification,
    validate_proteins,
    filter_proteins_by_score,
    keep_n_best_protein_hits,
    filter_proteins_by_rank,
    remove_decoy_protein_hits,
    remove_empty_protein_identifications
);

/// Keep all best-scoring ties, or remove the entire tied set when strict is true.
pub fn keep_best_peptide_hits(ids: &mut [PeptideIdentification], strict: bool) -> Result<()> {
    validate_peptides(ids)?;
    for id in ids {
        let higher = id.higher_score_better;
        id.hits.sort_by(|a, b| order(a.score, b.score, higher));
        if let Some(best) = id.hits.first().map(|hit| hit.score) {
            if strict && id.hits.get(1).is_some_and(|hit| hit.score == best) {
                id.hits.clear();
            } else {
                id.hits.retain(|hit| hit.score == best);
            }
        }
    }
    Ok(())
}

/// Rank spectra by their best hit, retaining all hits within selected spectra.
/// Score type and direction must agree across all records. Empty spectra sort last.
pub fn keep_n_best_spectra(ids: &mut Vec<PeptideIdentification>, n: usize) -> Result<()> {
    validate_peptides(ids)?;
    if let Some(first) = ids.first() {
        if ids.iter().any(|id| {
            id.score_type != first.score_type || id.higher_score_better != first.higher_score_better
        }) {
            return Err(bad(
                "spectrum ranking requires one score type and direction",
            ));
        }
    }
    for id in ids.iter_mut() {
        let higher = id.higher_score_better;
        id.hits.sort_by(|a, b| order(a.score, b.score, higher));
    }
    ids.sort_by(|a, b| match (a.hits.first(), b.hits.first()) {
        (Some(x), Some(y)) => order(x.score, y.score, a.higher_score_better),
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (None, None) => Ordering::Equal,
    });
    ids.truncate(n);
    Ok(())
}

fn retain_peptides(
    ids: &mut [PeptideIdentification],
    keep: impl Fn(&PeptideHit) -> bool,
) -> Result<()> {
    validate_peptides(ids)?;
    for id in ids {
        id.hits.retain(&keep);
    }
    Ok(())
}

/// Inclusive residue counts; modifications do not contribute to length.
pub fn filter_peptides_by_length(
    ids: &mut [PeptideIdentification],
    min: usize,
    max: Option<usize>,
) -> Result<()> {
    bounds(min, max)?;
    retain_peptides(ids, |hit| {
        hit.sequence.len() >= min && max.is_none_or(|max| hit.sequence.len() <= max)
    })
}

/// Inclusive signed charge range, including unknown charge zero when requested.
pub fn filter_peptides_by_charge(
    ids: &mut [PeptideIdentification],
    min: i32,
    max: Option<i32>,
) -> Result<()> {
    bounds(min, max)?;
    retain_peptides(ids, |hit| {
        hit.charge >= min && max.is_none_or(|max| hit.charge <= max)
    })
}

/// Match any evidence accession. This filters whole hits, leaving their evidence intact.
pub fn filter_peptides_by_accessions(
    ids: &mut [PeptideIdentification],
    accessions: &BTreeSet<String>,
    action: MatchAction,
) -> Result<()> {
    retain_peptides(ids, |hit| {
        action.keeps(
            hit.evidences.iter().any(|e| {
                !e.protein_accession.is_empty() && accessions.contains(&e.protein_accession)
            }),
        )
    })
}

pub fn filter_proteins_by_accessions(
    ids: &mut [ProteinIdentification],
    accessions: &BTreeSet<String>,
    action: MatchAction,
) -> Result<()> {
    validate_proteins(ids)?;
    for id in ids {
        id.hits
            .retain(|hit| action.keeps(accessions.contains(&hit.accession)));
    }
    Ok(())
}

fn matching_modification(sequence: &AASequence, names: &BTreeSet<String>) -> bool {
    if names.is_empty() {
        return sequence.is_modified();
    }
    sequence
        .n_terminal_modification()
        .is_some_and(|m| names.contains(m.full_id()))
        || sequence
            .c_terminal_modification()
            .is_some_and(|m| names.contains(m.full_id()))
        || (0..sequence.len()).any(|i| {
            sequence
                .residue_modification(i)
                .expect("index bounded by sequence length")
                .is_some_and(|m| names.contains(m.full_id()))
        })
}

/// Match exact full modification IDs, including both termini. Empty names match any modification.
/// Known modifications use registry full IDs; anonymous tags use the full ID
/// returned by the sequence annotation. This does not search registry mass ranges
/// or treat unresolved B/Z/X residues as modifications.
pub fn filter_peptides_by_modifications(
    ids: &mut [PeptideIdentification],
    names: &BTreeSet<String>,
    action: MatchAction,
) -> Result<()> {
    retain_peptides(ids, |hit| {
        action.keeps(matching_modification(&hit.sequence, names))
    })
}

/// Extract canonical modified strings, or unmodified sequences when requested.
pub fn extract_peptide_sequences(
    ids: &[PeptideIdentification],
    ignore_modifications: bool,
) -> Result<BTreeSet<String>> {
    validate_peptides(ids)?;
    Ok(ids
        .iter()
        .flat_map(|id| &id.hits)
        .map(|hit| {
            if ignore_modifications {
                hit.sequence.as_str().to_owned()
            } else {
                hit.sequence.to_string()
            }
        })
        .collect())
}

/// Match canonical strings from extract_peptide_sequences. Charge is ignored.
pub fn filter_peptides_by_sequences(
    ids: &mut [PeptideIdentification],
    sequences: &BTreeSet<String>,
    ignore_modifications: bool,
    action: MatchAction,
) -> Result<()> {
    retain_peptides(ids, |hit| {
        action.keeps(if ignore_modifications {
            sequences.contains(hit.sequence.as_str())
        } else {
            sequences.contains(&hit.sequence.to_string())
        })
    })
}

/// Keep the source's exact protein_references="unique" annotation, not an inferred count.
pub fn keep_unique_peptides_per_protein(ids: &mut [PeptideIdentification]) -> Result<()> {
    retain_peptides(ids, |hit| {
        has_string_value(&hit.metadata, "protein_references", "unique")
    })
}

/// Keep the first occurrence within each identification. No score sorting is performed.
pub fn remove_duplicate_peptide_hits(
    ids: &mut [PeptideIdentification],
    policy: DuplicatePolicy,
) -> Result<()> {
    validate_peptides(ids)?;
    if policy == DuplicatePolicy::Exact {
        let mut comparisons = 0_usize;
        for id in ids.iter() {
            let n = id.hits.len();
            let pairs = n
                .checked_mul(n.saturating_sub(1))
                .map(|v| v / 2)
                .ok_or_else(|| bad("duplicate comparison budget overflow"))?;
            comparisons = comparisons
                .checked_add(pairs)
                .ok_or_else(|| bad("duplicate comparison budget overflow"))?;
            if comparisons > MAX_EXACT_DUPLICATE_COMPARISONS {
                return Err(bad("exact duplicate filtering exceeds comparison budget"));
            }
        }
    }
    for id in ids {
        if policy == DuplicatePolicy::Sequence {
            let mut seen = BTreeSet::new();
            id.hits.retain(|hit| seen.insert(hit.sequence.clone()));
        } else {
            // ponytail: source full-record equality needs an O(n²) scan, bounded
            // above; add a canonical hash only if larger candidate sets need it.
            let mut unique = Vec::with_capacity(id.hits.len());
            for hit in std::mem::take(&mut id.hits) {
                if !unique.contains(&hit) {
                    unique.push(hit);
                }
            }
            id.hits = unique;
        }
    }
    Ok(())
}

/// Filter complete identifications by inclusive RT range. Missing RT does not match.
pub fn filter_peptides_by_rt(
    ids: &mut Vec<PeptideIdentification>,
    min: f64,
    max: f64,
) -> Result<()> {
    finite(min)?;
    finite(max)?;
    bounds(min, Some(max))?;
    validate_peptides(ids)?;
    ids.retain(|id| id.rt.is_some_and(|rt| rt >= min && rt <= max));
    Ok(())
}

/// Filter complete identifications by inclusive precursor m/z range. Missing m/z does not match.
pub fn filter_peptides_by_mz(
    ids: &mut Vec<PeptideIdentification>,
    min: f64,
    max: f64,
) -> Result<()> {
    finite(min)?;
    finite(max)?;
    bounds(min, Some(max))?;
    validate_peptides(ids)?;
    ids.retain(|id| id.mz.is_some_and(|mz| mz >= min && mz <= max));
    Ok(())
}

/// Keep hits within precursor mass tolerance (ppm is relative to observed m/z).
/// Unknown charge zero is treated as +1, matching OpenMS. Missing/nonpositive
/// observed m/z and charges whose calculated m/z is invalid are errors, atomically.
/// Unresolved residue masses likewise fail before any hit is removed. Absolute
/// mass tags with known monoisotopic masses can be filtered without a formula.
pub fn filter_peptides_by_mz_error(
    ids: &mut [PeptideIdentification],
    tolerance: Tolerance,
) -> Result<()> {
    let (Tolerance::Absolute(value) | Tolerance::Ppm(value)) = tolerance;
    finite(value)?;
    if value < 0.0 {
        return Err(bad("mass tolerance must be nonnegative"));
    }
    validate_peptides(ids)?;
    // Precompute fallible mass calculations before removing any hits.
    let masks: Vec<Vec<bool>> =
        ids.iter()
            .map(|id| {
                let mz = id.mz.filter(|mz| *mz > 0.0).ok_or_else(|| {
                    bad("precursor mass filtering requires positive observed m/z")
                })?;
                let delta = match tolerance {
                    Tolerance::Absolute(v) => v,
                    Tolerance::Ppm(v) => v * (mz / 1e6),
                };
                finite(delta)?;
                id.hits
                    .iter()
                    .map(|hit| {
                        let charge = if hit.charge == 0 { 1 } else { hit.charge };
                        Ok((mz - hit.sequence.mz(charge)?).abs() <= delta)
                    })
                    .collect()
            })
            .collect::<Result<_>>()?;
    for (id, mask) in ids.iter_mut().zip(masks) {
        let mut decisions = mask.into_iter();
        id.hits
            .retain(|_| decisions.next().expect("one decision per hit"));
    }
    Ok(())
}

type RunAccessions = BTreeMap<String, BTreeSet<String>>;
fn referenced_accessions<'a>(
    ids: impl Iterator<Item = &'a PeptideIdentification>,
) -> RunAccessions {
    let mut runs = RunAccessions::new();
    for id in ids {
        let accessions = runs.entry(id.identifier.clone()).or_default();
        for hit in &id.hits {
            accessions.extend(hit.protein_accessions().into_iter().map(str::to_owned));
        }
    }
    runs
}
fn protein_accessions(ids: &[ProteinIdentification]) -> RunAccessions {
    let mut runs = RunAccessions::new();
    for id in ids {
        runs.entry(id.identifier.clone())
            .or_default()
            .extend(id.hits.iter().map(|h| h.accession.clone()));
    }
    runs
}
fn retain_referenced_proteins(ids: &mut [ProteinIdentification], runs: &RunAccessions) {
    for id in ids {
        let accessions = runs.get(&id.identifier);
        id.hits
            .retain(|hit| accessions.is_some_and(|set| set.contains(&hit.accession)));
    }
}
fn retain_valid_references(
    ids: &mut [PeptideIdentification],
    runs: &RunAccessions,
    remove_unmatched: bool,
) {
    for id in ids {
        let accessions = runs.get(&id.identifier);
        for hit in &mut id.hits {
            hit.evidences
                .retain(|e| accessions.is_some_and(|set| set.contains(&e.protein_accession)));
        }
        if remove_unmatched {
            id.hits.retain(|hit| !hit.evidences.is_empty());
        }
    }
}

/// Keep proteins referenced by peptide hits in the same run identifier. Groups are unchanged.
pub fn remove_unreferenced_proteins(
    proteins: &mut [ProteinIdentification],
    peptides: &[PeptideIdentification],
) -> Result<()> {
    validate_proteins(proteins)?;
    validate_peptides(peptides)?;
    retain_referenced_proteins(proteins, &referenced_accessions(peptides.iter()));
    Ok(())
}

/// Keep evidence pointing to surviving proteins within the same run. Repeated
/// protein records with the same identifier contribute their union of accessions.
pub fn remove_dangling_protein_references(
    peptides: &mut [PeptideIdentification],
    proteins: &[ProteinIdentification],
    remove_peptides_without_reference: bool,
) -> Result<()> {
    validate_peptides(peptides)?;
    validate_proteins(proteins)?;
    retain_valid_references(
        peptides,
        &protein_accessions(proteins),
        remove_peptides_without_reference,
    );
    Ok(())
}

/// Remove absent group members and empty groups, preserving sample-indexed arrays.
/// Returns false only if a *surviving* group lost some members, as in OpenMS.
pub fn update_protein_groups(groups: &mut Vec<ProteinGroup>, hits: &[ProteinHit]) -> Result<bool> {
    groups.iter().try_for_each(ProteinGroup::validate)?;
    hits.iter().try_for_each(ProteinHit::validate)?;
    let accessions: BTreeSet<_> = hits.iter().map(|h| h.accession.as_str()).collect();
    let mut valid = true;
    for group in groups.iter_mut() {
        let previous = group.accessions.len();
        group.accessions.retain(|a| accessions.contains(a.as_str()));
        if !group.accessions.is_empty() && group.accessions.len() < previous {
            valid = false;
        }
    }
    groups.retain(|group| !group.accessions.is_empty());
    Ok(valid)
}

pub fn remove_ungrouped_proteins(
    groups: &[ProteinGroup],
    hits: &mut Vec<ProteinHit>,
) -> Result<()> {
    groups.iter().try_for_each(ProteinGroup::validate)?;
    hits.iter().try_for_each(ProteinHit::validate)?;
    let accessions: BTreeSet<_> = groups.iter().flat_map(|g| &g.accessions).collect();
    hits.retain(|h| accessions.contains(&h.accession));
    Ok(())
}

// Source map helpers visit top-level features and unassigned peptide IDs.
// Subordinate features and feature records themselves are never removed here.
macro_rules! map_filters {
    ($map:ty, $best:ident, $unreferenced:ident, $dangling:ident, $empty:ident) => {
        /// Keep N best hits in each top-level feature and unassigned peptide identification.
        pub fn $best(map: &mut $map, n: usize) -> Result<()> {
            map.validate()?;
            for feature in &mut map.features {
                keep_n_best_peptide_hits(&mut feature.peptide_identifications, n)?;
            }
            keep_n_best_peptide_hits(&mut map.unassigned_peptide_identifications, n)
        }
        /// Keep referenced proteins, optionally counting unassigned peptides too.
        pub fn $unreferenced(map: &mut $map, include_unassigned: bool) -> Result<()> {
            map.validate()?;
            let unassigned = if include_unassigned {
                &map.unassigned_peptide_identifications[..]
            } else {
                &[]
            };
            let runs = referenced_accessions(
                map.features
                    .iter()
                    .flat_map(|f| &f.peptide_identifications)
                    .chain(unassigned),
            );
            retain_referenced_proteins(&mut map.protein_identifications, &runs);
            Ok(())
        }
        /// Clean evidence on top-level features and unassigned peptides against stored protein runs.
        pub fn $dangling(map: &mut $map, remove_peptides_without_reference: bool) -> Result<()> {
            map.validate()?;
            let runs = protein_accessions(&map.protein_identifications);
            for feature in &mut map.features {
                retain_valid_references(
                    &mut feature.peptide_identifications,
                    &runs,
                    remove_peptides_without_reference,
                );
            }
            retain_valid_references(
                &mut map.unassigned_peptide_identifications,
                &runs,
                remove_peptides_without_reference,
            );
            Ok(())
        }
        /// Remove empty peptide identifications; all features and protein identifications remain.
        pub fn $empty(map: &mut $map) -> Result<()> {
            map.validate()?;
            for feature in &mut map.features {
                feature
                    .peptide_identifications
                    .retain(|id| !id.hits.is_empty());
            }
            map.unassigned_peptide_identifications
                .retain(|id| !id.hits.is_empty());
            Ok(())
        }
    };
}
map_filters!(
    FeatureMap,
    keep_n_best_hits_in_feature_map,
    remove_unreferenced_feature_proteins,
    remove_dangling_feature_references,
    remove_empty_feature_identifications
);
map_filters!(
    ConsensusMap,
    keep_n_best_hits_in_consensus_map,
    remove_unreferenced_consensus_proteins,
    remove_dangling_consensus_references,
    remove_empty_consensus_identifications
);
