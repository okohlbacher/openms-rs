// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `PROCESSING/ID/IDFilter.h` helpers and filters beyond the original subset
//! (core `bc9cc12`): `countHits`, `getBestHit`, `filterGroupsByScore`,
//! `filterPeptidesByRTPredictPValue`, `removePeptidesWithMatchingRegEx`, the
//! best-per-peptide annotation family, `extractUnassignedProteins` and the
//! reference-run overload of `removeDanglingProteinReferences`.

use super::predicates::{
    HasMatchingAccession, HasMaxMetaValue, HasMetaValue, Identification, IdentificationHit,
    PeptideIdentificationMap, move_matching_items, try_keep_matching_items,
};
use super::{RunAccessions, finite, retain_valid_references, validate_peptides};
use crate::concept::boost_regex::BoostRegex;
use crate::identification::{
    PeptideIdentification, ProteinGroup, ProteinHit, ProteinIdentification,
};
use crate::kernel::ConsensusMap;
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet, HashSet};

/// The meta value the best-per-peptide annotation writes (1 for the best hit
/// of its sequence, 0 for the others).
pub const BEST_PER_PEPTIDE: &str = "best_per_peptide";

/// Total number of hits in a list of peptide or protein identifications
/// (`countHits`).
pub fn count_hits<I: Identification>(ids: &[I]) -> usize {
    ids.iter().map(|id| id.hits().len()).sum()
}

/// The best-scoring hit across identifications (`getBestHit`); `None` when no
/// identification has a hit.
///
/// Source rules, kept: identifications without hits are skipped; the score
/// orientation is that of the *first* identification with hits; a later hit
/// replaces the current best only when strictly better, so the first of
/// several equal best hits wins; with `assume_sorted` only the first hit of
/// each identification is looked at.
///
/// # Errors
///
/// [`Error::InvalidValue`] ("Can't compare scores of different types") when two
/// identifications with hits have different score types, as the source throws
/// `Exception::InvalidValue`.
pub fn best_hit<I: Identification>(ids: &[I], assume_sorted: bool) -> Result<Option<I::Hit>> {
    let mut best: Option<(&I, &I::Hit)> = None;
    for id in ids {
        let Some(first) = id.hits().first() else {
            continue;
        };
        let (best_id, mut best_hit) = match best {
            None => (id, first),
            Some((best_id, best_hit)) => {
                if best_id.score_type() != id.score_type() {
                    return Err(Error::InvalidValue(format!(
                        "Can't compare scores of different types: {}/{}",
                        best_id.score_type(),
                        id.score_type()
                    )));
                }
                (best_id, best_hit)
            }
        };
        let higher_better = best_id.higher_score_better();
        for hit in id.hits() {
            if (higher_better && hit.score() > best_hit.score())
                || (!higher_better && hit.score() < best_hit.score())
            {
                best_hit = hit;
            }
            if assume_sorted {
                break;
            }
        }
        best = Some((best_id, best_hit));
    }
    Ok(best.map(|(_, hit)| hit.clone()))
}

/// Keeps the protein groups whose probability passes the threshold
/// (`filterGroupsByScore`).
///
/// The source's rule is kept although it is asymmetric: with `higher_better`
/// a group is kept only when its probability is **strictly** above
/// `threshold_score`, otherwise when it is at most `threshold_score`. A group
/// whose probability equals the threshold is therefore removed in the first
/// case and kept in the second; the hit filters are inclusive in both.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `threshold_score` is not finite or a group is
/// invalid; nothing is removed then.
pub fn filter_groups_by_score(
    groups: &mut Vec<ProteinGroup>,
    threshold_score: f64,
    higher_better: bool,
) -> Result<()> {
    finite(threshold_score)?;
    groups.iter().try_for_each(ProteinGroup::validate)?;
    groups.retain(|group| {
        !((higher_better && threshold_score >= group.probability)
            || (!higher_better && threshold_score < group.probability))
    });
    Ok(())
}

/// Filters peptide hits by the RTPredict p-value in `metavalue_key`
/// (`filterPeptidesByRTPredictPValue`), returning how many hits were removed
/// for lacking the meta value (the source's warning count).
///
/// A hit is kept when the meta value is present and at most
/// `1 - threshold`, the source's cutoff (its own comment asks "why?"). The
/// usual keys are `predicted_RT_p_value` and `predicted_RT_p_value_first_dim`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `threshold` is not finite, an identification
/// is invalid, or a hit's value is not numeric (the source then compares a
/// reinterpreted pointer); nothing is removed then.
pub fn filter_peptides_by_rt_predict_p_value(
    ids: &mut [PeptideIdentification],
    metavalue_key: &str,
    threshold: f64,
) -> Result<usize> {
    finite(threshold)?;
    validate_peptides(ids)?;
    let present = HasMetaValue::new(metavalue_key, MetaValue::default());
    let pvalue = HasMaxMetaValue {
        key: metavalue_key.to_owned(),
        value: 1.0 - threshold,
    };
    // Decide every hit first, so a non-numeric value leaves all unchanged.
    for id in ids.iter() {
        for hit in &id.hits {
            if present.matches(hit) {
                pvalue.matches(hit)?;
            }
        }
    }
    let mut missing = 0usize;
    for id in ids.iter_mut() {
        let before = id.hits.len();
        id.hits.retain(|hit| present.matches(hit));
        missing += before - id.hits.len();
        try_keep_matching_items(&mut id.hits, |hit| pvalue.matches(hit))?;
    }
    Ok(missing)
}

/// Removes the peptide hits whose unmodified sequence contains a match of
/// `regex` (`removePeptidesWithMatchingRegEx`).
///
/// The source compiles `std::regex` (ECMAScript grammar) and calls
/// `regex_search`. This compiles the crate's Boost-compatible expression
/// ([`BoostRegex`]); the two grammars agree on character classes, anchors,
/// alternation, groups, quantifiers and lookahead, and differ in constructs
/// only one of them has (e.g. lookbehind, possessive quantifiers, `\A`).
///
/// # Errors
///
/// [`Error::InvalidValue`] for an expression that does not compile (the source
/// throws `std::regex_error`) or an invalid identification; nothing is removed
/// then.
pub fn remove_peptides_with_matching_regex(
    ids: &mut [PeptideIdentification],
    regex: &str,
) -> Result<()> {
    validate_peptides(ids)?;
    let expression = BoostRegex::new(regex)?;
    let mut decisions = Vec::new();
    for id in ids.iter() {
        for hit in &id.hits {
            decisions.push(expression.is_search_match(hit.sequence.as_str().as_bytes())?);
        }
    }
    let mut decisions = decisions.into_iter();
    for id in ids {
        id.hits.retain(|_| !decisions.next().unwrap_or(false));
    }
    Ok(())
}

/// Position of a processed hit: list, identification, hit.
type HitPosition = (usize, usize, usize);

/// The best hit so far for one (run, sequence, charge) key.
struct Winner {
    position: HitPosition,
    score: f64,
}

/// Winners of the best-per-peptide annotation (the source's
/// `RunToSequenceToChargeToPepHitP`, with positions instead of pointers).
#[derive(Default)]
struct BestPerPeptideIndex {
    winners: BTreeMap<(String, String, i32), Winner>,
    processed: Vec<HitPosition>,
}

impl BestPerPeptideIndex {
    /// The source's `annotateBestPerPeptideWithData` for one identification,
    /// minus the meta-value writes, which [`Self::write`] applies afterwards:
    /// each processed hit ends with 1 exactly when it is its key's final
    /// winner, which is what the source's sequence of writes leaves.
    fn visit(
        &mut self,
        list: usize,
        index: usize,
        id: &mut PeptideIdentification,
        run: Option<&str>,
        ignore_mods: bool,
        ignore_charges: bool,
        nr_best_spectrum: usize,
    ) -> Result<()> {
        let higher_score_better = id.higher_score_better;
        // Source: make sure the first hit is the best one.
        id.sort()?;
        let limit = if nr_best_spectrum == 0 || id.hits.len() <= nr_best_spectrum {
            id.hits.len()
        } else {
            nr_best_spectrum
        };
        for (position, hit) in id.hits.iter().take(limit).enumerate() {
            let sequence = if ignore_mods {
                hit.sequence.as_str().to_owned()
            } else {
                hit.sequence.to_string()
            };
            let charge = if ignore_charges { 0 } else { hit.charge };
            let key = (run.unwrap_or("").to_owned(), sequence, charge);
            let here = (list, index, position);
            self.processed.push(here);
            match self.winners.get_mut(&key) {
                None => {
                    self.winners.insert(
                        key,
                        Winner {
                            position: here,
                            score: hit.score,
                        },
                    );
                }
                Some(winner) => {
                    if (higher_score_better && hit.score > winner.score)
                        || (!higher_score_better && hit.score < winner.score)
                    {
                        winner.position = here;
                        winner.score = hit.score;
                    }
                }
            }
        }
        Ok(())
    }

    /// Writes `best_per_peptide` (1 or 0) on every processed hit of list `list`.
    fn write(&self, list: usize, ids: &mut [PeptideIdentification]) {
        let winners: HashSet<HitPosition> = self.winners.values().map(|w| w.position).collect();
        for &(l, index, position) in &self.processed {
            if l != list {
                continue;
            }
            let value = i64::from(winners.contains(&(l, index, position)));
            ids[index].hits[position]
                .metadata
                .insert(BEST_PER_PEPTIDE.to_owned(), MetaValue::from(value));
        }
    }
}

fn annotate_lists(
    lists: &mut [&mut Vec<PeptideIdentification>],
    per_run: bool,
    ignore_mods: bool,
    ignore_charges: bool,
    nr_best_spectrum: usize,
) -> Result<()> {
    for list in lists.iter() {
        validate_peptides(list)?;
    }
    let mut index = BestPerPeptideIndex::default();
    for (l, list) in lists.iter_mut().enumerate() {
        for (i, id) in list.iter_mut().enumerate() {
            let run = per_run.then(|| id.identifier.clone());
            index.visit(
                l,
                i,
                id,
                run.as_deref(),
                ignore_mods,
                ignore_charges,
                nr_best_spectrum,
            )?;
        }
    }
    for (l, list) in lists.iter_mut().enumerate() {
        index.write(l, list);
    }
    Ok(())
}

fn keep_best(ids: &mut [PeptideIdentification]) {
    let best = HasMetaValue::new(BEST_PER_PEPTIDE, MetaValue::from(1i64));
    for id in ids {
        id.hits.retain(|hit| best.matches(hit));
    }
}

/// Annotates, over all identifications regardless of run, the best hit per
/// peptide sequence (`annotateBestPerPeptide`).
///
/// Each identification is sorted by score first (its hits are reordered),
/// then its first `nr_best_spectrum` hits (all when 0) are visited in order.
/// The key is the modified sequence (unmodified with `ignore_mods`) and the
/// charge (ignored with `ignore_charges`). The visited hit that ends as its
/// key's best gets `best_per_peptide` = 1, every other visited hit 0; hits
/// beyond `nr_best_spectrum` get no annotation. A later hit wins only with a
/// strictly better score in *its own* identification's orientation.
///
/// # Errors
///
/// [`Error::InvalidValue`] for an invalid identification, before anything is
/// sorted or annotated.
pub fn annotate_best_per_peptide(
    ids: &mut Vec<PeptideIdentification>,
    ignore_mods: bool,
    ignore_charges: bool,
    nr_best_spectrum: usize,
) -> Result<()> {
    annotate_lists(
        &mut [ids],
        false,
        ignore_mods,
        ignore_charges,
        nr_best_spectrum,
    )
}

/// Keeps only the best hit per peptide sequence (`keepBestPerPeptide`):
/// [`annotate_best_per_peptide`], then keeps the hits annotated 1, which
/// removes the hits beyond `nr_best_spectrum` too.
///
/// # Errors
///
/// As [`annotate_best_per_peptide`].
pub fn keep_best_per_peptide(
    ids: &mut Vec<PeptideIdentification>,
    ignore_mods: bool,
    ignore_charges: bool,
    nr_best_spectrum: usize,
) -> Result<()> {
    annotate_best_per_peptide(ids, ignore_mods, ignore_charges, nr_best_spectrum)?;
    keep_best(ids);
    Ok(())
}

/// [`annotate_best_per_peptide`] separately for each identification run
/// (`annotateBestPerPeptidePerRun`), the run being the identification's
/// [`PeptideIdentification::identifier`].
///
/// The source also takes the protein runs to pre-create their entries; a run
/// without a protein identification is still grouped by its identifier, so
/// the protein runs do not change the result and are not needed here.
///
/// # Errors
///
/// As [`annotate_best_per_peptide`].
pub fn annotate_best_per_peptide_per_run(
    ids: &mut Vec<PeptideIdentification>,
    ignore_mods: bool,
    ignore_charges: bool,
    nr_best_spectrum: usize,
) -> Result<()> {
    annotate_lists(
        &mut [ids],
        true,
        ignore_mods,
        ignore_charges,
        nr_best_spectrum,
    )
}

/// [`annotate_best_per_peptide_per_run`], then keeps the hits annotated 1
/// (`keepBestPerPeptidePerRun`).
///
/// # Errors
///
/// As [`annotate_best_per_peptide`].
pub fn keep_best_per_peptide_per_run(
    ids: &mut Vec<PeptideIdentification>,
    ignore_mods: bool,
    ignore_charges: bool,
    nr_best_spectrum: usize,
) -> Result<()> {
    annotate_best_per_peptide_per_run(ids, ignore_mods, ignore_charges, nr_best_spectrum)?;
    keep_best(ids);
    Ok(())
}

/// Collects a map's peptide identification lists, features first, then the
/// unassigned list, and runs `f` over them together.
fn with_map_lists(
    map: &mut impl PeptideIdentificationMap,
    f: impl FnOnce(&mut [&mut Vec<PeptideIdentification>]) -> Result<()>,
) -> Result<()> {
    let mut lists: Vec<Vec<PeptideIdentification>> = Vec::new();
    map.for_each_peptide_identification_list(&mut |ids| {
        lists.push(std::mem::take(ids));
        Ok(())
    })?;
    let mut refs: Vec<&mut Vec<PeptideIdentification>> = lists.iter_mut().collect();
    let result = f(&mut refs);
    // Put the lists back whether or not `f` failed; `f` validates before it
    // changes anything.
    let mut lists = lists.into_iter();
    map.for_each_peptide_identification_list(&mut |ids| {
        *ids = lists.next().unwrap_or_default();
        Ok(())
    })?;
    result
}

/// [`annotate_best_per_peptide_per_run`] across a whole map, one index shared
/// by every feature's identifications and the unassigned ones
/// (`annotateBestPerPeptidePerRun` for a map type).
///
/// # Errors
///
/// As [`annotate_best_per_peptide`]; the map is unchanged then.
pub fn annotate_best_per_peptide_per_run_in_map(
    map: &mut impl PeptideIdentificationMap,
    ignore_mods: bool,
    ignore_charges: bool,
    nr_best_spectrum: usize,
) -> Result<()> {
    with_map_lists(map, |lists| {
        annotate_lists(lists, true, ignore_mods, ignore_charges, nr_best_spectrum)
    })
}

/// [`annotate_best_per_peptide_per_run_in_map`], then keeps the hits annotated
/// 1 in every peptide identification of the map (`keepBestPerPeptidePerRun`
/// for a map type).
///
/// # Errors
///
/// As [`annotate_best_per_peptide`]; the map is unchanged then.
pub fn keep_best_per_peptide_per_run_in_map(
    map: &mut impl PeptideIdentificationMap,
    ignore_mods: bool,
    ignore_charges: bool,
    nr_best_spectrum: usize,
) -> Result<()> {
    with_map_lists(map, |lists| {
        annotate_lists(lists, true, ignore_mods, ignore_charges, nr_best_spectrum)?;
        for list in lists.iter_mut() {
            keep_best(list);
        }
        Ok(())
    })
}

/// Moves out of a consensus map's protein runs every protein hit that no
/// feature's peptide references, and returns them by run identifier
/// (`extractUnassignedProteins`).
///
/// Only peptides assigned to features count; unassigned peptides do not. A
/// run's hits are split with the source's `std::partition` order (see
/// [`move_matching_items`]), so both the kept and the returned hits can be
/// reordered exactly as in the source. Two protein runs with the same
/// identifier append to one result entry.
///
/// # Errors
///
/// [`Error::InvalidValue`] for an invalid map, before anything moves.
pub fn extract_unassigned_proteins(
    map: &mut ConsensusMap,
) -> Result<BTreeMap<String, Vec<ProteinHit>>> {
    map.validate()?;
    let mut referenced: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for feature in &map.features {
        for id in &feature.peptide_identifications {
            let accessions = referenced.entry(id.identifier.clone()).or_default();
            for hit in &id.hits {
                accessions.extend(hit.protein_accessions().into_iter().map(str::to_owned));
            }
        }
    }
    let empty = BTreeSet::new();
    let mut result: BTreeMap<String, Vec<ProteinHit>> = BTreeMap::new();
    for run in &mut map.protein_identifications {
        let accessions = referenced.get(&run.identifier).unwrap_or(&empty);
        let target = result.entry(run.identifier.clone()).or_default();
        let referenced_hit = HasMatchingAccession::new(accessions);
        move_matching_items(&mut run.hits, |hit| !referenced_hit.matches(hit), target);
    }
    Ok(result)
}

/// Removes, from every peptide hit of a consensus map (features and
/// unassigned), the evidences whose accession is not a hit of `reference_run`,
/// optionally removing hits left without evidence
/// (`removeDanglingProteinReferences(ConsensusMap&, const ProteinIdentification&, bool)`).
///
/// Unlike [`super::remove_dangling_consensus_references`], every peptide is
/// checked against the one reference run, whatever its own run identifier.
///
/// # Errors
///
/// [`Error::InvalidValue`] for an invalid map or reference run, before
/// anything is removed.
pub fn remove_dangling_consensus_references_to_run(
    map: &mut ConsensusMap,
    reference_run: &ProteinIdentification,
    remove_peptides_without_reference: bool,
) -> Result<()> {
    map.validate()?;
    reference_run.validate()?;
    let accessions: BTreeSet<String> = reference_run
        .hits
        .iter()
        .map(|hit| hit.accession.clone())
        .collect();
    // Every identification is looked up under its own run key, so give each
    // key present in the map the reference run's accessions.
    let mut runs = RunAccessions::new();
    for id in map
        .features
        .iter()
        .flat_map(|feature| &feature.peptide_identifications)
        .chain(&map.unassigned_peptide_identifications)
    {
        runs.entry(id.identifier.clone())
            .or_insert_with(|| accessions.clone());
    }
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
