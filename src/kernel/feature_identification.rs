// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Hendrik Weisser, Chris Bielow, Timo Sachsenberg, OpenMS Rust contributors $

//! Identification surface of the feature containers.
//!
//! Ports the annotation and identification members of `KERNEL/BaseFeature.h`,
//! `KERNEL/Feature.h` and `KERNEL/ConsensusFeature.h`: the annotation state of
//! attached peptide identifications, the score-ordered sort of those
//! identifications, the primary identified molecule and observation-match set
//! that reference an `IdentificationData` graph, the reference update after
//! that graph is copied, the subordinate-recursive unique-ID traversal, and the
//! experimental ratio list of a consensus feature.
//!
//! The types themselves live in
//! [`kernel::features`](crate::kernel::features); this module only adds
//! operations, exactly as
//! [`kernel::gap_closures`](crate::kernel::gap_closures) does, because
//! `kernel/features.rs` is frozen for this work package. Rust permits inherent
//! implementations for a crate-local type from any module of the same crate,
//! so callers see one type with one API.
//!
//! Identification records attach to features **by reference**. The source's
//! `FeatureMap` owns an `IdentificationData id_data_` member (`FeatureMap.h:294`);
//! here the graph is not `Clone` by design, so embedding it would remove `Clone`
//! and `PartialEq` from every map. Operations that must resolve a reference take
//! the graph, or the
//! [`ReferenceTranslator`](crate::identification::graph::ReferenceTranslator),
//! as a parameter instead. See
//! `docs/FEATURE_IDENTIFICATION_SUPPORT.md`.

use super::features::{BaseFeature, ConsensusFeature, Feature, Ratio};
use crate::identification::graph::{
    IdentificationData, IdentifiedMolecule, ObservationMatchId, ReferenceTranslator,
};
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::collections::BTreeSet;
use std::fmt;

/// Source metadata key recording the map a peptide identification came from.
const MAP_INDEX: &str = "map_index";

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn limit() -> Error {
    Error::InvalidValue("feature identification work limit exceeded".into())
}

/// State of the identifications attached to a feature.
///
/// Ports `BaseFeature::AnnotationState`. When one identification carries
/// several hits, only its best hit is considered, as in the source.
///
/// The source's trailing `SIZE_OF_ANNOTATIONSTATE` sentinel is not ported: it
/// exists only to size `NamesOfAnnotationState`, and [`AnnotationState::NAMES`]
/// is a fixed-length array whose length is the variant count.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AnnotationState {
    /// No identification, or none of the attached identifications has a hit.
    #[default]
    None,
    /// Exactly one identification with at least one hit.
    Single,
    /// Several identifications whose best hits all agree.
    MultipleSame,
    /// Several identifications whose best hits disagree, which usually means a
    /// bad mapping, as the source comment states.
    MultipleDivergent,
}

impl AnnotationState {
    /// Source `BaseFeature::NamesOfAnnotationState`, in variant order.
    pub const NAMES: [&'static str; 4] = [
        "no ID",
        "single ID",
        "multiple IDs (identical)",
        "multiple IDs (divergent)",
    ];

    /// The source name of this state, as written by featureXML and the TOPP tools.
    pub const fn name(self) -> &'static str {
        Self::NAMES[self as usize]
    }
}

impl fmt::Display for AnnotationState {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

impl BaseFeature {
    /// Attached peptide identifications a checked operation will process.
    pub const MAX_PEPTIDE_IDENTIFICATIONS: usize = 1_000_000;
    /// Peptide hits a checked operation will process across all identifications.
    pub const MAX_PEPTIDE_HITS: usize = 4_000_000;
    /// Observation matches one feature may reference.
    pub const MAX_ID_MATCHES: usize = 1_000_000;

    /// Copy this feature, recording `map_index` on every attached peptide
    /// identification.
    ///
    /// Ports `BaseFeature(const BaseFeature& rhs, UInt64 map_index)`, which the
    /// source uses to build a consensus feature from one map's feature so that
    /// the identifications remember where they came from. The metadata key is
    /// `map_index`, as in the source.
    ///
    /// [`ConsensusFeature::from_feature`] does *not* apply this, while the
    /// source's `ConsensusFeature(UInt64, const BaseFeature&)` constructor
    /// does; call this first when the stamp is wanted.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when more than
    /// [`BaseFeature::MAX_PEPTIDE_IDENTIFICATIONS`] identifications are
    /// attached, when one of them is invalid, or when `map_index` exceeds the
    /// signed metadata integer range. `self` is only read.
    pub fn clone_with_map_index(&self, map_index: u64) -> Result<Self> {
        if self.peptide_identifications.len() > Self::MAX_PEPTIDE_IDENTIFICATIONS {
            return Err(limit());
        }
        for identification in &self.peptide_identifications {
            identification.validate()?;
        }
        let value = MetaValue::try_from(map_index)?;
        let mut copy = self.clone();
        for identification in &mut copy.peptide_identifications {
            identification
                .metadata
                .insert(MAP_INDEX.into(), value.clone());
        }
        Ok(copy)
    }

    /// Sort the attached peptide identifications, best first.
    ///
    /// Ports `BaseFeature::sortPeptideIdentifications`. Each identification's
    /// hits are sorted by score first, then the identifications themselves are
    /// ordered by their best hit: descending score when
    /// `higher_score_better` is set, ascending otherwise. Identifications
    /// without hits sort last.
    ///
    /// The source's `@note` that identifications are assumed to share one score
    /// type is stated as a check here: the comparator reads
    /// `isHigherScoreBetter()` from its left operand only, so identifications
    /// that disagree give `std::sort` an asymmetric comparator and therefore
    /// undefined behaviour. This port rejects that input instead.
    ///
    /// Three further differences from the source are deliberate:
    ///
    /// * The source sorts an identification's hits inside the comparator, so an
    ///   element that never takes part in a comparison keeps unsorted hits. All
    ///   hits are sorted here, which is the evident intent.
    /// * The source comparator answers "less" for two empty identifications,
    ///   which is not a strict weak ordering. Empty identifications simply sort
    ///   last here and keep their relative order.
    /// * The order is stable, so identifications with equal best scores keep
    ///   their relative order; the source's `std::sort` may permute them.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when more than
    /// [`BaseFeature::MAX_PEPTIDE_IDENTIFICATIONS`] identifications or
    /// [`BaseFeature::MAX_PEPTIDE_HITS`] hits are attached, when an attached
    /// record is invalid (a non-finite score, for instance), or when the
    /// identifications disagree on `higher_score_better`. Every check runs
    /// before anything is reordered, so a failure leaves the feature unchanged.
    pub fn sort_peptide_identifications(&mut self) -> Result<()> {
        let identifications = &self.peptide_identifications;
        if identifications.len() > Self::MAX_PEPTIDE_IDENTIFICATIONS {
            return Err(limit());
        }
        let mut hits = 0usize;
        let mut direction: Option<bool> = None;
        for identification in identifications {
            identification.validate()?;
            hits = hits
                .checked_add(identification.hits.len())
                .ok_or_else(limit)?;
            if hits > Self::MAX_PEPTIDE_HITS {
                return Err(limit());
            }
            if identification.hits.is_empty() {
                continue;
            }
            match direction {
                Some(previous) if previous != identification.higher_score_better => {
                    return Err(invalid(
                        "peptide identifications disagree on score direction; \
                         the source comparator is undefined for mixed score types",
                    ));
                }
                _ => direction = Some(identification.higher_score_better),
            }
        }
        // Nothing below can fail, so the feature is either fully sorted or untouched.
        let higher_better = direction.unwrap_or(true);
        for identification in &mut self.peptide_identifications {
            identification
                .hits
                .sort_by(|a, b| score_order(a.score, b.score, higher_better));
        }
        self.peptide_identifications.sort_by(|a, b| {
            match (a.hits.first(), b.hits.first()) {
                (Some(left), Some(right)) => score_order(left.score, right.score, higher_better),
                // Empty identifications sort last, as the source places them.
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            }
        });
        Ok(())
    }

    /// State of the identifications attached to this feature.
    ///
    /// Ports `BaseFeature::getAnnotationState`. Observation matches take
    /// precedence: when [`BaseFeature::id_matches`] is non-empty only those are
    /// considered, and the attached peptide identifications are ignored
    /// entirely, exactly as the source dispatches.
    ///
    /// With matches, one match is [`AnnotationState::Single`] and several are
    /// [`AnnotationState::MultipleSame`] or
    /// [`AnnotationState::MultipleDivergent`] according to whether every match
    /// names the same identified molecule. That branch never reports
    /// [`AnnotationState::None`], because an empty match set takes the other
    /// branch. Every match is resolved even once divergence is established,
    /// where the source returns at the first difference; a stale reference late
    /// in the set is therefore still reported.
    ///
    /// Without matches, the peptide identifications decide: none, or none with
    /// hits, is [`AnnotationState::None`]; a single identification carrying at
    /// least one hit is [`AnnotationState::Single`]; otherwise the best hit of
    /// every identification that has one is compared by sequence. Note the
    /// source consequence, preserved here: two identifications of which only
    /// one has hits collect a single sequence and therefore report
    /// [`AnnotationState::MultipleSame`], not `Single`.
    ///
    /// # Arguments
    ///
    /// * `graph` — the `IdentificationData` the matches were registered in.
    ///   Pass `None` when the feature has no matches; the source reaches the
    ///   graph through the references themselves, which this port cannot do
    ///   because an ID is owner-tagged rather than an iterator.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when the feature has matches and
    /// `graph` is `None`, [`Error::InvalidValue`] when a match does not belong
    /// to `graph` or an attached peptide record is invalid, and the limit error
    /// of [`BaseFeature::sort_peptide_identifications`] when more records are
    /// attached than a checked operation processes. The feature is only read.
    pub fn annotation_state(&self, graph: Option<&IdentificationData>) -> Result<AnnotationState> {
        if !self.id_matches.is_empty() {
            if self.id_matches.len() > Self::MAX_ID_MATCHES {
                return Err(limit());
            }
            let graph = graph.ok_or_else(|| {
                Error::MissingInformation(
                    "feature references observation matches, but no identification graph was given"
                        .into(),
                )
            })?;
            if self.id_matches.len() == 1 {
                // Still resolved, so a stale or foreign reference is reported.
                for &id in &self.id_matches {
                    graph.observation_match(id)?;
                }
                return Ok(AnnotationState::Single);
            }
            let mut molecule: Option<IdentifiedMolecule> = None;
            let mut divergent = false;
            for &id in &self.id_matches {
                let current = graph.observation_match(id)?.identified_molecule;
                match molecule {
                    None => molecule = Some(current),
                    Some(first) if first != current => divergent = true,
                    Some(_) => {}
                }
            }
            return Ok(if divergent {
                AnnotationState::MultipleDivergent
            } else {
                AnnotationState::MultipleSame
            });
        }
        let identifications = &self.peptide_identifications;
        if identifications.len() > Self::MAX_PEPTIDE_IDENTIFICATIONS {
            return Err(limit());
        }
        if identifications.is_empty() {
            return Ok(AnnotationState::None);
        }
        if identifications.len() == 1 && !identifications[0].hits.is_empty() {
            return Ok(AnnotationState::Single);
        }
        let mut sequences: BTreeSet<String> = BTreeSet::new();
        let mut hits = 0usize;
        for identification in identifications {
            hits = hits
                .checked_add(identification.hits.len())
                .ok_or_else(limit)?;
            if hits > Self::MAX_PEPTIDE_HITS {
                return Err(limit());
            }
            // Source `getAnnotationState` sorts a copy and reads hit 0; the
            // native best-hit accessor is the same comparison without the copy.
            // Where scores tie, this keeps the first hit in stored order, while
            // the source's unsorted `std::sort` may report any tied hit.
            if let Some(hit) = identification.best_hit()? {
                sequences.insert(hit.sequence.to_string());
            }
        }
        Ok(match sequences.len() {
            0 => AnnotationState::None,
            1 => AnnotationState::MultipleSame,
            _ => AnnotationState::MultipleDivergent,
        })
    }

    /// Whether a primary identified molecule has been assigned.
    ///
    /// Ports `BaseFeature::hasPrimaryID`.
    pub const fn has_primary_id(&self) -> bool {
        self.primary_id.is_some()
    }

    /// The primary identified molecule (peptide, RNA or compound) of this feature.
    ///
    /// Ports `BaseFeature::getPrimaryID`. The returned ID is a reference into
    /// the `IdentificationData` graph it was registered in; resolve it there.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when no primary ID was assigned,
    /// matching the source `@throw Exception::MissingInformation`.
    pub fn primary_id(&self) -> Result<IdentifiedMolecule> {
        self.primary_id
            .ok_or_else(|| Error::MissingInformation("no primary ID assigned".into()))
    }

    /// Assign the primary identified molecule.
    ///
    /// Ports `BaseFeature::setPrimaryID`. Like the source, the ID is stored
    /// without consulting a graph; use
    /// [`BaseFeature::set_primary_id_checked`] to reject a foreign or stale
    /// reference at assignment time instead of at first use.
    pub fn set_primary_id(&mut self, id: impl Into<IdentifiedMolecule>) {
        self.primary_id = Some(id.into());
    }

    /// Assign the primary identified molecule after resolving it in `graph`.
    ///
    /// Native addition without a source counterpart: because an ID carries its
    /// owning graph generation, an assignment can be checked. The source cannot
    /// check, since its reference is a container iterator.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `id` does not belong to `graph` or
    /// its record has been removed. The feature is unchanged in that case.
    pub fn set_primary_id_checked(
        &mut self,
        graph: &IdentificationData,
        id: impl Into<IdentifiedMolecule>,
    ) -> Result<()> {
        let id = id.into();
        match id {
            IdentifiedMolecule::Peptide(value) => {
                graph.peptide(value)?;
            }
            IdentifiedMolecule::Compound(value) => {
                graph.compound(value)?;
            }
            IdentifiedMolecule::Oligo(value) => {
                graph.oligo(value)?;
            }
        }
        self.primary_id = Some(id);
        Ok(())
    }

    /// Remove any assigned primary identified molecule.
    ///
    /// Ports `BaseFeature::clearPrimaryID`. Attached observation matches and
    /// peptide identifications are not touched.
    pub const fn clear_primary_id(&mut self) {
        self.primary_id = None;
    }

    /// Reference one more observation match (for example a PSM) from this feature.
    ///
    /// Ports `BaseFeature::addIDMatch`. Re-adding a referenced match is a
    /// no-op, as for the source `std::set`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the feature already references
    /// [`BaseFeature::MAX_ID_MATCHES`] matches and `id` is not among them; the
    /// set is unchanged. The reference itself is not resolved here, matching
    /// the source; use [`BaseFeature::add_id_match_checked`] to resolve it.
    pub fn add_id_match(&mut self, id: ObservationMatchId) -> Result<()> {
        if self.id_matches.len() >= Self::MAX_ID_MATCHES && !self.id_matches.contains(&id) {
            return Err(limit());
        }
        self.id_matches.insert(id);
        Ok(())
    }

    /// Reference one more observation match after resolving it in `graph`.
    ///
    /// Native addition without a source counterpart, for the same reason as
    /// [`BaseFeature::set_primary_id_checked`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `id` does not belong to `graph`,
    /// its record has been removed, or the match limit is reached. The set is
    /// unchanged in each case.
    pub fn add_id_match_checked(
        &mut self,
        graph: &IdentificationData,
        id: ObservationMatchId,
    ) -> Result<()> {
        graph.observation_match(id)?;
        self.add_id_match(id)
    }

    /// Translate this feature's primary ID and observation matches into the
    /// graph the translator describes.
    ///
    /// Ports `BaseFeature::updateIDReferences`, needed after the
    /// `IdentificationData` holding the referenced records has been copied or
    /// merged into another graph. `translator` is the value returned by
    /// `IdentificationData::merge_from`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the translator has no entry for the
    /// primary ID or for one of the matches, and when more than
    /// [`BaseFeature::MAX_ID_MATCHES`] matches are referenced. The source has
    /// an `allow_missing` mode that keeps an untranslatable reference; that
    /// mode is not ported, because a reference the new graph does not own
    /// cannot be resolved later.
    ///
    /// Unlike the source, which swaps the match set out before translating and
    /// so leaves a partially translated set behind when a translation throws,
    /// every reference is translated into a temporary that replaces the stored
    /// values only once all of them succeed.
    pub fn update_id_references(&mut self, translator: &ReferenceTranslator) -> Result<()> {
        if self.id_matches.len() > Self::MAX_ID_MATCHES {
            return Err(limit());
        }
        let primary = self
            .primary_id
            .map(|id| translator.molecule(id))
            .transpose()?;
        let mut matches = BTreeSet::new();
        for &id in &self.id_matches {
            matches.insert(translator.observation_match(id)?);
        }
        self.primary_id = primary;
        self.id_matches = matches;
        Ok(())
    }
}

impl Feature {
    /// Features one checked traversal of a feature and its subordinates visits.
    pub const MAX_TRAVERSED_FEATURES: usize = 1_000_000;

    /// Apply `visit` to the unique ID of this feature and of every subordinate,
    /// returning the sum of its results.
    ///
    /// Ports the mutating `Feature::applyMemberFunction` overload, whose C++
    /// member-function pointer is a closure here. The source passes a
    /// `UniqueIdInterface` member such as `hasInvalidUniqueId` or
    /// `ensureUniqueId` and accumulates the returned `Size` values, so `u64`
    /// (which implements [`HasUniqueId`](crate::concept::HasUniqueId)) is the
    /// argument. Visiting order is the source's: this feature first, then each
    /// subordinate subtree in order, depth first.
    ///
    /// ```
    /// use openms::concept::HasUniqueId;
    /// use openms::kernel::features::Feature;
    ///
    /// let mut feature = Feature::default();
    /// feature.subordinates.push(Feature::default());
    /// // Source: `f.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId)`.
    /// let invalid = feature.for_each_unique_id(|id| usize::from(id.has_invalid_unique_id()))?;
    /// assert_eq!(invalid, 2);
    /// let cleared = feature.for_each_unique_id(|id| id.clear_unique_id())?;
    /// assert_eq!(cleared, 0);
    /// # Ok::<(), openms::Error>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the tree holds more than
    /// [`Feature::MAX_TRAVERSED_FEATURES`] features, when it is deeper than
    /// [`Feature::MAX_SUBORDINATE_DEPTH`], or when the accumulated count
    /// overflows. The limits are checked while walking, so `visit` may already
    /// have run on the features visited before the limit was reached; the
    /// source is unbounded and recursive.
    pub fn for_each_unique_id(
        &mut self,
        mut visit: impl FnMut(&mut u64) -> usize,
    ) -> Result<usize> {
        let mut total = 0usize;
        let mut visited = 0usize;
        // Reverse pushes so that popping yields the source's pre-order walk.
        let mut pending: Vec<(&mut Self, usize)> = vec![(self, 0)];
        while let Some((feature, depth)) = pending.pop() {
            check_traversal(&mut visited, depth)?;
            total = total
                .checked_add(visit(&mut feature.base.unique_id))
                .ok_or_else(limit)?;
            pending.extend(
                feature
                    .subordinates
                    .iter_mut()
                    .rev()
                    .map(|sub| (sub, depth + 1)),
            );
        }
        Ok(total)
    }

    /// Apply `visit` to the unique ID of this feature and of every subordinate
    /// without mutating them, returning the sum of its results.
    ///
    /// Ports the `const` `Feature::applyMemberFunction` overload, which exists
    /// in the source only because a member-function pointer is either `const`
    /// or not. Order, accumulation and limits match
    /// [`Feature::for_each_unique_id`].
    ///
    /// # Errors
    ///
    /// As [`Feature::for_each_unique_id`].
    pub fn count_unique_ids(&self, mut visit: impl FnMut(u64) -> usize) -> Result<usize> {
        let mut total = 0usize;
        let mut visited = 0usize;
        let mut pending: Vec<(&Self, usize)> = vec![(self, 0)];
        while let Some((feature, depth)) = pending.pop() {
            check_traversal(&mut visited, depth)?;
            total = total
                .checked_add(visit(feature.base.unique_id))
                .ok_or_else(limit)?;
            pending.extend(
                feature
                    .subordinates
                    .iter()
                    .rev()
                    .map(|sub| (sub, depth + 1)),
            );
        }
        Ok(total)
    }

    /// Translate the identification references of this feature and of every
    /// subordinate into the graph the translator describes.
    ///
    /// Ports `Feature::updateAllIDReferences`, which recurses into subordinate
    /// features after updating the feature itself.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`BaseFeature::update_id_references`], plus the
    /// traversal limits of [`Feature::for_each_unique_id`]. Every reference in
    /// the whole tree is translated in a first, read-only pass before any
    /// feature is changed, so a missing translation anywhere leaves the entire
    /// tree unchanged. The source updates as it descends and leaves a partly
    /// translated tree behind.
    pub fn update_all_id_references(&mut self, translator: &ReferenceTranslator) -> Result<()> {
        let mut visited = 0usize;
        let mut pending: Vec<(&Self, usize)> = vec![(self, 0)];
        while let Some((feature, depth)) = pending.pop() {
            check_traversal(&mut visited, depth)?;
            if feature.base.id_matches.len() > BaseFeature::MAX_ID_MATCHES {
                return Err(limit());
            }
            if let Some(id) = feature.base.primary_id {
                translator.molecule(id)?;
            }
            for &id in &feature.base.id_matches {
                translator.observation_match(id)?;
            }
            pending.extend(
                feature
                    .subordinates
                    .iter()
                    .rev()
                    .map(|sub| (sub, depth + 1)),
            );
        }
        let mut pending: Vec<&mut Self> = vec![self];
        while let Some(feature) = pending.pop() {
            feature.base.update_id_references(translator)?;
            pending.extend(feature.subordinates.iter_mut());
        }
        Ok(())
    }
}

impl Ratio {
    /// Check that the ratio value is finite.
    ///
    /// The source `ConsensusFeature::Ratio` validates nothing and its default
    /// constructor leaves `ratio_value_` uninitialised; the Rust
    /// [`Default`] is `0.0`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `ratio_value` is not finite.
    pub fn validate(&self) -> Result<()> {
        if !self.ratio_value.is_finite() {
            return Err(invalid("consensus feature ratio value must be finite"));
        }
        Ok(())
    }
}

impl ConsensusFeature {
    /// Ratios one consensus feature may store.
    pub const MAX_RATIOS: usize = 1_000_000;

    /// Attach one more quantification ratio.
    ///
    /// Ports `ConsensusFeature::addRatio`. The source `@note` still applies:
    /// ratios are experimental and the consensus feature handler ignores them.
    /// Duplicate ratios are appended, as by the source `push_back`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the ratio value is not finite or
    /// the feature already stores [`ConsensusFeature::MAX_RATIOS`] ratios. The
    /// stored list is unchanged in both cases; the source appends unchecked.
    pub fn add_ratio(&mut self, ratio: Ratio) -> Result<()> {
        ratio.validate()?;
        if self.ratios.len() >= Self::MAX_RATIOS {
            return Err(limit());
        }
        self.ratios.push(ratio);
        Ok(())
    }

    /// Replace the quantification ratios.
    ///
    /// Ports `ConsensusFeature::setRatios`, which takes a non-`const`
    /// `std::vector<Ratio>&` and therefore cannot be called with a temporary or
    /// a `const` vector; this takes the vector by value.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when more than
    /// [`ConsensusFeature::MAX_RATIOS`] ratios are given or any ratio value is
    /// not finite. All ratios are checked before the list is replaced, so a
    /// failure leaves the previous ratios in place.
    pub fn set_ratios(&mut self, ratios: Vec<Ratio>) -> Result<()> {
        if ratios.len() > Self::MAX_RATIOS {
            return Err(limit());
        }
        for ratio in &ratios {
            ratio.validate()?;
        }
        self.ratios = ratios;
        Ok(())
    }

    /// The attached quantification ratios, in insertion order.
    ///
    /// Ports both `ConsensusFeature::getRatios` overloads; the mutable one is
    /// the public [`ConsensusFeature::ratios`] field, which bypasses the checks
    /// of [`ConsensusFeature::set_ratios`] exactly as the source's mutable
    /// accessor does.
    pub fn ratios(&self) -> &[Ratio] {
        &self.ratios
    }
}

/// Count one visited feature against the traversal ceilings.
fn check_traversal(visited: &mut usize, depth: usize) -> Result<()> {
    if depth > Feature::MAX_SUBORDINATE_DEPTH {
        return Err(invalid("subordinate feature depth exceeds 128"));
    }
    *visited = visited.checked_add(1).ok_or_else(limit)?;
    if *visited > Feature::MAX_TRAVERSED_FEATURES {
        return Err(limit());
    }
    Ok(())
}

/// Better-first score order, as `PeptideIdentification::sort`.
fn score_order(left: f64, right: f64, higher_better: bool) -> std::cmp::Ordering {
    if higher_better {
        right
            .partial_cmp(&left)
            .unwrap_or(std::cmp::Ordering::Equal)
    } else {
        left.partial_cmp(&right)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}
