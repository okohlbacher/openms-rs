// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The predicate functors of `PROCESSING/ID/IDFilter.h` (core `bc9cc12`) and
//! the higher-order filters that apply them.
//!
//! Each source functor is a small public struct whose `matches` method is the
//! source's `operator()`. The higher-order filters take any closure, so a
//! functor is passed as `|hit| predicate.matches(hit)`. Predicates that convert
//! a meta value to a number, or compute an m/z, can fail; they return
//! `Result<bool>` and are applied with the `try_` filters, which decide every
//! item before removing any.

use super::matching_modification;
use crate::format::fasta::FASTAEntry;
use crate::identification::{
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
};
use crate::kernel::{ConsensusMap, FeatureMap};
use crate::metadata::{MetaInfo, MetaValue, MetaValueData};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet, HashSet};

/// A hit the generic filters can test: [`PeptideHit`] or [`ProteinHit`]
/// (the source's `HitType` template parameter).
pub trait IdentificationHit: Clone {
    /// The hit's main score (`getScore`).
    fn score(&self) -> f64;
    /// The hit's meta values (`MetaInfoInterface`).
    fn metadata(&self) -> &MetaInfo;
    /// Whether any of the hit's protein accessions is in `accessions`
    /// (`HasMatchingAccession`): a peptide's non-empty evidence accessions, a
    /// protein's own accession.
    fn has_matching_accession(&self, accessions: &dyn AccessionSet) -> bool;
}

impl IdentificationHit for PeptideHit {
    fn score(&self) -> f64 {
        self.score
    }
    fn metadata(&self) -> &MetaInfo {
        &self.metadata
    }
    fn has_matching_accession(&self, accessions: &dyn AccessionSet) -> bool {
        // Source: `extractProteinAccessionsSet()`, which skips empty accessions.
        self.protein_accessions()
            .into_iter()
            .any(|accession| accessions.contains_accession(accession))
    }
}

impl IdentificationHit for ProteinHit {
    fn score(&self) -> f64 {
        self.score
    }
    fn metadata(&self) -> &MetaInfo {
        &self.metadata
    }
    fn has_matching_accession(&self, accessions: &dyn AccessionSet) -> bool {
        accessions.contains_accession(&self.accession)
    }
}

/// A peptide or protein identification: its hits and score settings (the
/// source's `IdentificationType` template parameter).
pub trait Identification {
    /// [`PeptideHit`] or [`ProteinHit`].
    type Hit: IdentificationHit;
    /// The hits (`getHits() const`).
    fn hits(&self) -> &[Self::Hit];
    /// The hits, for filtering (`getHits()`).
    fn hits_mut(&mut self) -> &mut Vec<Self::Hit>;
    /// The score type name (`getScoreType`).
    fn score_type(&self) -> &str;
    /// Score orientation (`isHigherScoreBetter`).
    fn higher_score_better(&self) -> bool;
}

impl Identification for PeptideIdentification {
    type Hit = PeptideHit;
    fn hits(&self) -> &[PeptideHit] {
        &self.hits
    }
    fn hits_mut(&mut self) -> &mut Vec<PeptideHit> {
        &mut self.hits
    }
    fn score_type(&self) -> &str {
        &self.score_type
    }
    fn higher_score_better(&self) -> bool {
        self.higher_score_better
    }
}

impl Identification for ProteinIdentification {
    type Hit = ProteinHit;
    fn hits(&self) -> &[ProteinHit] {
        &self.hits
    }
    fn hits_mut(&mut self) -> &mut Vec<ProteinHit> {
        &mut self.hits
    }
    fn score_type(&self) -> &str {
        &self.score_type
    }
    fn higher_score_better(&self) -> bool {
        self.higher_score_better
    }
}

/// A set of protein accessions: the source's `SetType`, `std::set` for
/// `HasMatchingAccession` and `std::unordered_set` for
/// `HasMatchingAccessionUnordered`.
pub trait AccessionSet {
    /// Whether `accession` is in the set (`count`/`contains`).
    fn contains_accession(&self, accession: &str) -> bool;
}
impl AccessionSet for BTreeSet<String> {
    fn contains_accession(&self, accession: &str) -> bool {
        self.contains(accession)
    }
}
impl AccessionSet for HashSet<String> {
    fn contains_accession(&self, accession: &str) -> bool {
        self.contains(accession)
    }
}

/// A feature or consensus map whose peptide identifications the map-level
/// filters visit: those of every top-level feature, then the unassigned ones
/// (the source's `MapType`; `for (auto& feat : map)` plus
/// `getUnassignedPeptideIdentifications()`).
pub trait PeptideIdentificationMap {
    /// Calls `visit` with each top-level feature's identification list, then
    /// with the unassigned list.
    fn for_each_peptide_identification_list(
        &mut self,
        visit: &mut dyn FnMut(&mut Vec<PeptideIdentification>) -> Result<()>,
    ) -> Result<()>;
    /// The map's protein identification runs (`getProteinIdentifications`).
    fn protein_identifications(&self) -> &[ProteinIdentification];
}

macro_rules! peptide_identification_map {
    ($map:ty) => {
        impl PeptideIdentificationMap for $map {
            fn for_each_peptide_identification_list(
                &mut self,
                visit: &mut dyn FnMut(&mut Vec<PeptideIdentification>) -> Result<()>,
            ) -> Result<()> {
                for feature in &mut self.features {
                    visit(&mut feature.peptide_identifications)?;
                }
                visit(&mut self.unassigned_peptide_identifications)
            }
            fn protein_identifications(&self) -> &[ProteinIdentification] {
                &self.protein_identifications
            }
        }
    };
}
peptide_identification_map!(FeatureMap);
peptide_identification_map!(ConsensusMap);

// ---------------------------------------------------------------- functors

/// Is the score of this hit at least as good as the given value
/// (`HasGoodScore`)? Inclusive, in the given orientation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HasGoodScore {
    /// The cutoff.
    pub score: f64,
    /// Whether higher scores are better.
    pub higher_score_better: bool,
}
impl HasGoodScore {
    /// `score >= cutoff` when higher is better, else `score <= cutoff`.
    pub fn matches(&self, hit: &impl IdentificationHit) -> bool {
        if self.higher_score_better {
            hit.score() >= self.score
        } else {
            hit.score() <= self.score
        }
    }
}

/// The source's `DataValue` equality: same type and unit; a double within
/// an absolute `1e-6` (strictly less); strings, integers and lists exactly.
fn data_value_eq(found: &MetaValue, expected: &MetaValue) -> bool {
    if found.unit() != expected.unit() {
        return false;
    }
    match (found.data(), expected.data()) {
        (MetaValueData::Float(a), MetaValueData::Float(b)) => (a - b).abs() < 1e-6,
        (a, b) => std::mem::discriminant(a) == std::mem::discriminant(b) && a == b,
    }
}

/// The source's `double(DataValue)` for the numeric meta-value predicates.
///
/// Source converts an integer, returns a double as stored and throws
/// `ConversionError` for an empty value; for a string or list it reads the
/// union's `double` member, which holds pointer bits, and compares that
/// garbage. The port refuses those types.
fn meta_number(value: &MetaValue, key: &str) -> Result<f64> {
    value.as_f64().map_err(|_| {
        Error::InvalidValue(format!(
            "meta value {key:?} is not numeric; the source reinterprets such a value's storage as a double"
        ))
    })
}

/// Is a meta value with the given key, and value, set on this hit
/// (`HasMetaValue`)? An empty `value` only checks that the key is set.
#[derive(Clone, Debug, PartialEq)]
pub struct HasMetaValue {
    /// The meta-value key.
    pub key: String,
    /// The required value; empty means "any value".
    pub value: MetaValue,
}
impl HasMetaValue {
    /// A predicate for `key` with `value` (empty for "any value").
    pub fn new(key: impl Into<String>, value: MetaValue) -> Self {
        Self {
            key: key.into(),
            value,
        }
    }
    /// False when the key is absent or holds an empty value; true when
    /// `value` is empty; otherwise the source's `DataValue` equality, which
    /// accepts doubles within `1e-6` and requires the same unit.
    pub fn matches(&self, hit: &impl IdentificationHit) -> bool {
        match hit.metadata().get(&self.key) {
            None => false,
            Some(found) if found.is_empty() => false,
            Some(_) if self.value.is_empty() => true,
            Some(found) => data_value_eq(found, &self.value),
        }
    }
}

/// Does a meta value of this hit have at most the given value
/// (`HasMaxMetaValue`)?
#[derive(Clone, Debug, PartialEq)]
pub struct HasMaxMetaValue {
    /// The meta-value key.
    pub key: String,
    /// The inclusive upper bound.
    pub value: f64,
}
impl HasMaxMetaValue {
    /// False when the key is absent; else whether the numeric value is at
    /// most [`Self::value`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the meta value is not numeric (see the
    /// module notes on the source's union read).
    pub fn matches(&self, hit: &impl IdentificationHit) -> Result<bool> {
        match hit.metadata().get(&self.key) {
            None => Ok(false),
            Some(found) if found.is_empty() => Ok(false),
            Some(found) => Ok(meta_number(found, &self.key)? <= self.value),
        }
    }
}

/// Does a meta value of this hit have at least the given value
/// (`HasMinMetaValue`)?
#[derive(Clone, Debug, PartialEq)]
pub struct HasMinMetaValue {
    /// The meta-value key.
    pub key: String,
    /// The inclusive lower bound.
    pub value: f64,
}
impl HasMinMetaValue {
    /// False when the key is absent; else whether the numeric value is at
    /// least [`Self::value`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the meta value is not numeric.
    pub fn matches(&self, hit: &impl IdentificationHit) -> Result<bool> {
        match hit.metadata().get(&self.key) {
            None => Ok(false),
            Some(found) if found.is_empty() => Ok(false),
            Some(found) => Ok(meta_number(found, &self.key)? >= self.value),
        }
    }
}

/// Is this a decoy hit (`HasDecoyAnnotation`)? True for the string meta value
/// `target_decoy` = `decoy` or `isDecoy` = `true`, compared as the source's
/// `DataValue`s (same type, no unit).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct HasDecoyAnnotation;
impl HasDecoyAnnotation {
    /// Either annotation marks the hit as a decoy, even when they conflict.
    pub fn matches(&self, hit: &impl IdentificationHit) -> bool {
        HasMetaValue::new("target_decoy", MetaValue::from("decoy")).matches(hit)
            || HasMetaValue::new("isDecoy", MetaValue::from("true")).matches(hit)
    }
}

/// Given a set of protein accessions, do any occur in the annotations of this
/// hit or evidence (`HasMatchingAccession`, `HasMatchingAccessionUnordered`)?
pub struct HasMatchingAccession<'a> {
    /// The accessions; a `BTreeSet` or `HashSet` of strings.
    pub accessions: &'a dyn AccessionSet,
}
impl<'a> HasMatchingAccession<'a> {
    /// A predicate over `accessions`.
    pub fn new(accessions: &'a dyn AccessionSet) -> Self {
        Self { accessions }
    }
    /// For a peptide hit: any non-empty evidence accession is in the set; for
    /// a protein hit: its accession is.
    pub fn matches(&self, hit: &impl IdentificationHit) -> bool {
        hit.has_matching_accession(self.accessions)
    }
    /// For a peptide evidence: its accession is in the set (the source's
    /// `operator()(const PeptideEvidence&)`, which does not skip an empty
    /// accession).
    pub fn matches_evidence(&self, evidence: &PeptideEvidence) -> bool {
        self.accessions
            .contains_accession(&evidence.protein_accession)
    }
}

/// An index of FASTA entries by identifier, to resolve peptide evidences
/// (`GetMatchingItems<PeptideEvidence, FASTAFile::FASTAEntry>`).
///
/// A later entry with the same identifier replaces an earlier one, as the
/// source's `items[key] = &entry` does.
#[derive(Clone, Debug, Default)]
pub struct GetMatchingItems<'a> {
    /// Entries by identifier (the source's `ItemMap`).
    pub items: BTreeMap<&'a str, &'a FASTAEntry>,
}
impl<'a> GetMatchingItems<'a> {
    /// Indexes `records` by identifier.
    pub fn new(records: &'a [FASTAEntry]) -> Self {
        let mut items = BTreeMap::new();
        for record in records {
            items.insert(record.identifier.as_str(), record);
        }
        Self { items }
    }
    /// The key of an entry (`getKey`): its identifier.
    pub fn key(entry: &FASTAEntry) -> &str {
        &entry.identifier
    }
    /// The key of an evidence (`getHitKey`): its protein accession.
    pub fn hit_key(evidence: &PeptideEvidence) -> &str {
        &evidence.protein_accession
    }
    /// Whether an entry exists for the evidence's accession (`exists`).
    pub fn exists(&self, evidence: &PeptideEvidence) -> bool {
        self.items.contains_key(evidence.protein_accession.as_str())
    }
    /// The entry for the evidence's accession (`getValue`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no entry has the accession, where the
    /// source throws `Exception::InvalidParameter`.
    pub fn value(&self, evidence: &PeptideEvidence) -> Result<&'a FASTAEntry> {
        self.items
            .get(evidence.protein_accession.as_str())
            .copied()
            .ok_or_else(|| {
                Error::InvalidValue(format!(
                    "Accession: '{}'. peptide evidence accession not in data",
                    evidence.protein_accession
                ))
            })
    }
}

/// Does the sequence of this peptide hit have at least the given length
/// (`HasMinPeptideLength`)? Residues are counted; modifications are not.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HasMinPeptideLength {
    /// The inclusive minimum.
    pub length: usize,
}
impl HasMinPeptideLength {
    /// `sequence length >= length`.
    pub fn matches(&self, hit: &PeptideHit) -> bool {
        hit.sequence.len() >= self.length
    }
}

/// Does the charge of this peptide hit have at least the given value
/// (`HasMinCharge`)?
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HasMinCharge {
    /// The inclusive minimum.
    pub charge: i32,
}
impl HasMinCharge {
    /// `charge >= self.charge`.
    pub fn matches(&self, hit: &PeptideHit) -> bool {
        hit.charge >= self.charge
    }
}

/// Is the m/z error of this peptide hit below the given value
/// (`HasLowMZError`)?
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HasLowMZError {
    /// Observed precursor m/z.
    pub precursor_mz: f64,
    /// Absolute tolerance in Th (a ppm tolerance is converted on construction).
    pub tolerance: f64,
}
impl HasLowMZError {
    /// A predicate for `precursor_mz`; a ppm `tolerance` is converted to Th
    /// relative to `precursor_mz`, as the source's constructor does.
    pub fn new(precursor_mz: f64, tolerance: f64, unit_ppm: bool) -> Self {
        let tolerance = if unit_ppm {
            tolerance * (precursor_mz / 1.0e6)
        } else {
            tolerance
        };
        Self {
            precursor_mz,
            tolerance,
        }
    }
    /// `|precursor m/z - peptide m/z| <= tolerance`, the peptide m/z at the
    /// hit's charge, charge 0 counting as 1.
    ///
    /// # Errors
    ///
    /// Those of [`crate::chemistry::AASequence::mz`], for a sequence without a mass.
    pub fn matches(&self, hit: &PeptideHit) -> Result<bool> {
        let charge = if hit.charge == 0 { 1 } else { hit.charge };
        Ok((self.precursor_mz - hit.sequence.mz(charge)?).abs() <= self.tolerance)
    }
}

/// Given a set of modifications (full IDs), do any occur in the sequence of
/// this peptide hit (`HasMatchingModification`)? An empty set matches any
/// modified sequence.
#[derive(Clone, Copy, Debug)]
pub struct HasMatchingModification<'a> {
    /// Full IDs such as `Oxidation (M)` or `Acetyl (N-term)`.
    pub modifications: &'a BTreeSet<String>,
}
impl HasMatchingModification<'_> {
    /// Residue and terminal modifications are checked by full ID.
    pub fn matches(&self, hit: &PeptideHit) -> bool {
        matching_modification(&hit.sequence, self.modifications)
    }
}

/// Is the sequence of this peptide hit among a set of sequences
/// (`HasMatchingSequence`)? With `ignore_mods`, the unmodified sequence is
/// compared; otherwise the source's modified string form.
#[derive(Clone, Copy, Debug)]
pub struct HasMatchingSequence<'a> {
    /// The sequences, in the form [`super::extract_peptide_sequences`]
    /// returns.
    pub sequences: &'a BTreeSet<String>,
    /// Compare unmodified sequences.
    pub ignore_mods: bool,
}
impl HasMatchingSequence<'_> {
    /// Whether the hit's sequence text is in the set.
    pub fn matches(&self, hit: &PeptideHit) -> bool {
        if self.ignore_mods {
            self.sequences.contains(hit.sequence.as_str())
        } else {
            self.sequences.contains(&hit.sequence.to_string())
        }
    }
}

/// Is the list of peptide evidences of this peptide hit empty
/// (`HasNoEvidence`)?
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct HasNoEvidence;
impl HasNoEvidence {
    /// Whether the hit has no evidence.
    pub fn matches(&self, hit: &PeptideHit) -> bool {
        hit.evidences.is_empty()
    }
}

/// Is the list of hits of this peptide or protein identification empty
/// (`HasNoHits`)?
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct HasNoHits;
impl HasNoHits {
    /// Whether the identification has no hits.
    pub fn matches(&self, id: &impl Identification) -> bool {
        id.hits().is_empty()
    }
}

/// Is the retention time of this peptide identification in the given range
/// (`HasRTInRange`)? Inclusive; an identification without RT does not match
/// (the source's unset RT is NaN, which compares false).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HasRTInRange {
    /// Inclusive lower bound.
    pub rt_min: f64,
    /// Inclusive upper bound.
    pub rt_max: f64,
}
impl HasRTInRange {
    /// `rt_min <= rt <= rt_max`.
    pub fn matches(&self, id: &PeptideIdentification) -> bool {
        id.rt
            .is_some_and(|rt| rt >= self.rt_min && rt <= self.rt_max)
    }
}

/// Is the precursor m/z of this peptide identification in the given range
/// (`HasMZInRange`)? Inclusive; an identification without m/z does not match.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HasMZInRange {
    /// Inclusive lower bound.
    pub mz_min: f64,
    /// Inclusive upper bound.
    pub mz_max: f64,
}
impl HasMZInRange {
    /// `mz_min <= mz <= mz_max`.
    pub fn matches(&self, id: &PeptideIdentification) -> bool {
        id.mz
            .is_some_and(|mz| mz >= self.mz_min && mz <= self.mz_max)
    }
}

// ------------------------------------------------------- higher-order filters

/// Removes the items that satisfy `predicate`, keeping the others in order
/// (`removeMatchingItems`).
pub fn remove_matching_items<T>(items: &mut Vec<T>, mut predicate: impl FnMut(&T) -> bool) {
    items.retain(|item| !predicate(item));
}

/// Keeps the items that satisfy `predicate`, in order (`keepMatchingItems`).
pub fn keep_matching_items<T>(items: &mut Vec<T>, mut predicate: impl FnMut(&T) -> bool) {
    items.retain(|item| predicate(item));
}

/// [`remove_matching_items`] for a fallible predicate: every item is decided
/// before any is removed, so an error leaves `items` unchanged.
///
/// # Errors
///
/// The first error `predicate` returns.
pub fn try_remove_matching_items<T>(
    items: &mut Vec<T>,
    predicate: impl FnMut(&T) -> Result<bool>,
) -> Result<()> {
    let decisions: Vec<bool> = items.iter().map(predicate).collect::<Result<_>>()?;
    let mut decisions = decisions.into_iter();
    items.retain(|_| !decisions.next().unwrap_or(false));
    Ok(())
}

/// [`keep_matching_items`] for a fallible predicate, atomically.
///
/// # Errors
///
/// The first error `predicate` returns.
pub fn try_keep_matching_items<T>(
    items: &mut Vec<T>,
    predicate: impl FnMut(&T) -> Result<bool>,
) -> Result<()> {
    let decisions: Vec<bool> = items.iter().map(predicate).collect::<Result<_>>()?;
    let mut decisions = decisions.into_iter();
    items.retain(|_| decisions.next().unwrap_or(true));
    Ok(())
}

/// Moves the items that satisfy `predicate` to the end of `target`
/// (`moveMatchingItems`).
///
/// Order follows the source, which partitions with `std::partition`: an
/// unstable two-pointer exchange (the algorithm libstdc++ and libc++ both use
/// for bidirectional iterators) that swaps the leftmost matching item with the
/// rightmost non-matching one. Both the kept and the moved items can therefore
/// change order; this reproduces the same order.
pub fn move_matching_items<T>(
    items: &mut Vec<T>,
    mut predicate: impl FnMut(&T) -> bool,
    target: &mut Vec<T>,
) {
    // std::partition(items, not_fn(predicate)): non-matching items first.
    let mut first = 0usize;
    let mut last = items.len();
    loop {
        while first != last && !predicate(&items[first]) {
            first += 1;
        }
        if first == last {
            break;
        }
        last -= 1;
        while first != last && predicate(&items[last]) {
            last -= 1;
        }
        if first == last {
            break;
        }
        items.swap(first, last);
        first += 1;
    }
    target.extend(items.drain(first..));
}

/// Removes, in every identification, the hits that satisfy `predicate`
/// (`removeMatchingItemsUnroll`).
pub fn remove_matching_items_unroll<I: Identification>(
    ids: &mut [I],
    mut predicate: impl FnMut(&I::Hit) -> bool,
) {
    for id in ids {
        remove_matching_items(id.hits_mut(), &mut predicate);
    }
}

/// Keeps, in every identification, the hits that satisfy `predicate`
/// (`keepMatchingItemsUnroll`).
pub fn keep_matching_items_unroll<I: Identification>(
    ids: &mut [I],
    mut predicate: impl FnMut(&I::Hit) -> bool,
) {
    for id in ids {
        keep_matching_items(id.hits_mut(), &mut predicate);
    }
}

/// Keeps the peptide hits that satisfy `predicate` in every peptide
/// identification of a map, assigned and unassigned
/// (`keepMatchingPeptideHits`).
///
/// # Errors
///
/// None today; the signature leaves room for a fallible map visit.
pub fn keep_matching_peptide_hits(
    map: &mut impl PeptideIdentificationMap,
    mut predicate: impl FnMut(&PeptideHit) -> bool,
) -> Result<()> {
    map.for_each_peptide_identification_list(&mut |ids| {
        keep_matching_items_unroll(ids, &mut predicate);
        Ok(())
    })
}

/// Removes the peptide hits that satisfy `predicate` from every peptide
/// identification of a map (`removeMatchingPeptideHits`).
///
/// # Errors
///
/// None today; the signature leaves room for a fallible map visit.
pub fn remove_matching_peptide_hits(
    map: &mut impl PeptideIdentificationMap,
    mut predicate: impl FnMut(&PeptideHit) -> bool,
) -> Result<()> {
    map.for_each_peptide_identification_list(&mut |ids| {
        remove_matching_items_unroll(ids, &mut predicate);
        Ok(())
    })
}

/// Removes the peptide identifications that satisfy `predicate` from every
/// feature and from the unassigned list of a map
/// (`removeMatchingPeptideIdentifications`). For a plain list use
/// [`remove_matching_items`], the source's `PeptideIdentificationList`
/// specialisation.
///
/// # Errors
///
/// None today; the signature leaves room for a fallible map visit.
pub fn remove_matching_peptide_identifications(
    map: &mut impl PeptideIdentificationMap,
    mut predicate: impl FnMut(&PeptideIdentification) -> bool,
) -> Result<()> {
    map.for_each_peptide_identification_list(&mut |ids| {
        remove_matching_items(ids, &mut predicate);
        Ok(())
    })
}

/// Keeps, in every hit, the peptide evidences `filter` accepts
/// (`FilterPeptideEvidences`; also the declared but undefined
/// `filterPeptideEvidences`).
pub fn filter_peptide_evidences(
    ids: &mut [PeptideIdentification],
    mut filter: impl FnMut(&PeptideEvidence) -> bool,
) {
    for id in ids {
        for hit in &mut id.hits {
            hit.evidences.retain(|evidence| filter(evidence));
        }
    }
}

/// [`filter_peptide_evidences`] for a fallible filter, atomically.
///
/// # Errors
///
/// The first error `filter` returns; nothing is removed then.
pub fn try_filter_peptide_evidences(
    ids: &mut [PeptideIdentification],
    mut filter: impl FnMut(&PeptideEvidence) -> Result<bool>,
) -> Result<()> {
    let mut decisions = Vec::new();
    for id in ids.iter() {
        for hit in &id.hits {
            for evidence in &hit.evidences {
                decisions.push(filter(evidence)?);
            }
        }
    }
    let mut decisions = decisions.into_iter();
    for id in ids {
        for hit in &mut id.hits {
            hit.evidences.retain(|_| decisions.next().unwrap_or(true));
        }
    }
    Ok(())
}
