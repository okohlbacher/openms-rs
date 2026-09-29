// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Establish and validate canonical precursor/transition identifiers for
//! OpenSWATH libraries.
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathLibraryIDNormalizer.h` and its
//! translation unit. OpenSWATH algorithms use string identifiers in a
//! [`LightTargetedExperiment`], while the persistent PQP/OSWPQ schemas use
//! integer foreign keys. This module defines the application-level invariant
//! that the light experiment carries those canonical integer identifiers
//! encoded as decimal strings:
//!
//! - `LightCompound::id` is the canonical precursor id;
//! - `LightTransition::peptide_ref` references that canonical precursor id;
//! - `LightTransition::transition_name` is the canonical transition id.
//!
//! Arbitrary source identifiers (for example from TSV or TraML) are normalised
//! once at the OpenSWATH library-loading boundary
//! ([`OpenSwathLibraryIDNormalizer::normalize_source_ids`]). Database-backed
//! inputs (PQP/OSWPQ) already carry persistent numeric ids and are validated
//! ([`OpenSwathLibraryIDNormalizer::validate_canonical_ids`]) instead of
//! renumbered.
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_INFERENCE_SUPPORT.md).
//!
//! [`LightTargetedExperiment`]: crate::analysis::openswath::transition_experiment::LightTargetedExperiment

use crate::analysis::openswath::transition_experiment::{LightCompound, LightTargetedExperiment};
use crate::{Error, Result};
use std::collections::{HashMap, HashSet};

/// Library identifier normalisation and validation. Stateless; the source's
/// class has only static members.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OpenSwathLibraryIDNormalizer;

/// Source-to-canonical precursor provenance captured during library loading,
/// the source's `SourceIDMapping`.
///
/// The operational [`LightTargetedExperiment`] contains only canonical numeric
/// ids. This mapping retains the original precursor identifiers long enough
/// for callers that still need source-id semantics (for example configurable
/// target/decoy pairing) without reusing source identifiers as operational
/// foreign keys.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SourceIDMapping {
    /// Original/source precursor id to canonical operational precursor id.
    pub precursor_source_to_canonical: HashMap<String, String>,
    /// Canonical operational precursor id to original/source precursor id.
    pub precursor_canonical_to_source: HashMap<String, String>,
    /// Canonical operational transition id to original/source transition id.
    /// Source transition ids are not required to be unique, so only the
    /// canonical-to-source direction is retained.
    pub transition_canonical_to_source: HashMap<String, String>,
}

/// The source's file-static `hasConventionalDecoyPrefix_`: `DECOY`, `Decoy`
/// or `decoy` at the start.
fn has_conventional_decoy_prefix(id: &str) -> bool {
    id.starts_with("DECOY") || id.starts_with("Decoy") || id.starts_with("decoy")
}

/// The source's file-static `hasCanonicalIDFormat_`: a non-empty decimal
/// string without sign or leading zero (other than `"0"` itself) that parses
/// completely as a non-negative `i64`.
fn is_canonical_id(text: &str) -> bool {
    let bytes = text.as_bytes();
    if bytes.is_empty()
        || bytes[0] == b'-'
        || bytes[0] == b'+'
        || (bytes.len() > 1 && bytes[0] == b'0')
    {
        return false;
    }
    // `std::from_chars` accepts only digits after the sign checks above, and
    // fails on overflow; `str::parse` would also accept a leading `+`, which
    // is excluded above.
    bytes.iter().all(u8::is_ascii_digit) && text.parse::<i64>().is_ok_and(|v| v >= 0)
}

/// The source's file-static `validateCanonicalIDFormat_`.
fn validate_canonical_id_format(text: &str, entity: &str) -> Result<()> {
    if text.is_empty() {
        return Err(invalid_id(&format!("{entity} ID must not be empty"), text));
    }
    if !is_canonical_id(text) {
        return Err(invalid_id(
            &format!("{entity} ID must be a non-negative Int64 in canonical decimal form"),
            text,
        ));
    }
    Ok(())
}

/// The source's `Exception::InvalidValue(message, value)`.
fn invalid_id(message: &str, value: &str) -> Error {
    Error::InvalidValue(format!("{message}: '{value}'"))
}

impl OpenSwathLibraryIDNormalizer {
    /// Replace arbitrary source precursor/transition identifiers with
    /// deterministic canonical ids.
    ///
    /// Precursor ids are assigned by sorting all unique source compound ids
    /// lexicographically (byte order, as `std::sort` on `std::string`) and
    /// numbering them from zero, matching the persistent PQP convention.
    /// Compounds that are not referenced by any transition are then omitted
    /// from the operational experiment without compressing the remaining ids,
    /// so sparse precursor ids are expected; this keeps direct source loading
    /// consistent with source-to-PQP round-trips. Transition ids are assigned
    /// from their order in `exp.transitions`, also starting at zero. Existing
    /// transition `peptide_ref` values are rewritten to the new precursor ids.
    ///
    /// A transition that is not yet a decoy becomes one when its source
    /// precursor reference or its source transition id starts with `DECOY`,
    /// `Decoy` or `decoy`, while the source ids are still available.
    ///
    /// This is intended for source-oriented formats such as TSV and TraML. It
    /// must not be used to renumber libraries that already contain persistent
    /// canonical ids.
    ///
    /// The source rebuilds the experiment so that its internal
    /// compound-reference cache cannot keep source-id keys; the port's
    /// [`LightTargetedExperiment`] keeps no cache, and the rebuild is kept for
    /// atomicity: `exp` is replaced only after every check has passed, so an
    /// error leaves it unchanged, as in the source.
    ///
    /// # Returns
    ///
    /// Source-id provenance for the canonicalised experiment. The precursor
    /// maps cover every source compound, including compounds omitted from the
    /// operational `exp`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] (the source's
    /// `Exception::InvalidValue`) when a compound id is empty or duplicated,
    /// or a transition references an empty or unknown compound id.
    pub fn normalize_source_ids(exp: &mut LightTargetedExperiment) -> Result<SourceIDMapping> {
        let mut known: HashSet<&str> = HashSet::with_capacity(exp.compounds.len());
        for compound in &exp.compounds {
            if compound.id.is_empty() {
                return Err(invalid_id("Source compound ID must not be empty", ""));
            }
            if !known.insert(compound.id.as_str()) {
                return Err(invalid_id(
                    "Source compound IDs must be unique",
                    &compound.id,
                ));
            }
        }

        let mut source_ids: Vec<&str> = exp.compounds.iter().map(|c| c.id.as_str()).collect();
        source_ids.sort_unstable();
        let mut precursor_map: HashMap<String, String> = HashMap::with_capacity(source_ids.len());
        for (index, id) in source_ids.iter().enumerate() {
            precursor_map.insert((*id).to_string(), index.to_string());
        }

        let mut referenced: HashSet<&str> = HashSet::with_capacity(exp.transitions.len());
        for transition in &exp.transitions {
            let reference = transition.peptide_ref.as_str();
            if reference.is_empty() {
                return Err(invalid_id(
                    "Transition precursor reference must not be empty",
                    reference,
                ));
            }
            if !known.contains(reference) {
                return Err(invalid_id(
                    "Transition references an unknown source compound ID",
                    reference,
                ));
            }
            referenced.insert(reference);
        }

        let mut mapping = SourceIDMapping {
            precursor_canonical_to_source: precursor_map
                .iter()
                .map(|(source, canonical)| (canonical.clone(), source.clone()))
                .collect(),
            precursor_source_to_canonical: HashMap::new(),
            transition_canonical_to_source: HashMap::with_capacity(exp.transitions.len()),
        };

        let compounds: Vec<LightCompound> = exp
            .compounds
            .iter()
            .filter(|c| referenced.contains(c.id.as_str()))
            .map(|c| {
                let mut compound = c.clone();
                // Every source compound id is a key of `precursor_map`.
                if let Some(canonical) = precursor_map.get(&c.id) {
                    compound.id.clone_from(canonical);
                }
                compound
            })
            .collect();

        let mut transitions = exp.transitions.clone();
        for (index, transition) in transitions.iter_mut().enumerate() {
            // Validated above: every reference is a known source compound.
            let Some(canonical_ref) = precursor_map.get(&transition.peptide_ref) else {
                return Err(invalid_id(
                    "Transition references an unknown source compound ID",
                    &transition.peptide_ref,
                ));
            };
            if !transition.decoy()
                && (has_conventional_decoy_prefix(&transition.peptide_ref)
                    || has_conventional_decoy_prefix(&transition.transition_name))
            {
                transition.set_decoy(true);
            }
            let canonical_id = index.to_string();
            let source_transition_id =
                std::mem::replace(&mut transition.transition_name, canonical_id.clone());
            mapping
                .transition_canonical_to_source
                .insert(canonical_id, source_transition_id);
            transition.peptide_ref.clone_from(canonical_ref);
        }

        let normalized = LightTargetedExperiment {
            transitions,
            compounds,
            proteins: exp.proteins.clone(),
        };
        Self::validate_canonical_ids(&normalized)?;
        mapping.precursor_source_to_canonical = precursor_map;
        *exp = normalized;
        Ok(mapping)
    }

    /// Materialise decoy flags from a configurable source precursor prefix.
    ///
    /// Canonicalisation replaces source precursor references and transition
    /// names with numeric ids. This uses the provenance returned by
    /// [`Self::normalize_source_ids`] (or reconstructed by a persistent-format
    /// reader) to preserve prefix-based decoy semantics before filtering. The
    /// prefix on either the source precursor id or the source transition id
    /// sets the transition's decoy flag; existing decoy flags are retained. An
    /// empty `decoy_prefix` disables the materialisation.
    pub fn materialize_decoy_prefix(
        exp: &mut LightTargetedExperiment,
        source_ids: &SourceIDMapping,
        decoy_prefix: &str,
    ) {
        if decoy_prefix.is_empty() {
            return;
        }
        let decoy_precursors: HashSet<&str> = source_ids
            .precursor_source_to_canonical
            .iter()
            .filter(|(source, _)| source.starts_with(decoy_prefix))
            .map(|(_, canonical)| canonical.as_str())
            .collect();
        for transition in &mut exp.transitions {
            let source_transition_is_decoy = source_ids
                .transition_canonical_to_source
                .get(&transition.transition_name)
                .is_some_and(|source| source.starts_with(decoy_prefix));
            if decoy_precursors.contains(transition.peptide_ref.as_str())
                || source_transition_is_decoy
            {
                transition.set_decoy(true);
            }
        }
    }

    /// Resolve the canonical target precursor paired with a source-prefix
    /// decoy.
    ///
    /// `decoy_ref` is normally a canonical precursor id. The source precursor
    /// id is recovered from `source_ids`, `decoy_prefix` is removed there, and
    /// the matching target source id is mapped back to its canonical
    /// operational id. If no provenance mapping is available (both precursor
    /// maps empty), source-id input is supported as a compatibility path and
    /// the prefix is stripped from `decoy_ref` itself; that fallback is not
    /// used when a mapping is present, because a numeric-looking source id
    /// could otherwise collide with a canonical id.
    ///
    /// # Returns
    ///
    /// The matching target precursor id in the same operational domain as the
    /// loaded experiment, or `None` if the decoy cannot be paired (including
    /// for an empty `decoy_prefix`).
    pub fn canonical_target_for_decoy_precursor(
        decoy_ref: &str,
        source_ids: &SourceIDMapping,
        decoy_prefix: &str,
    ) -> Option<String> {
        if decoy_prefix.is_empty() {
            return None;
        }
        if let Some(source_decoy) = source_ids.precursor_canonical_to_source.get(decoy_ref) {
            let target_source = source_decoy.strip_prefix(decoy_prefix)?;
            return source_ids
                .precursor_source_to_canonical
                .get(target_source)
                .cloned();
        }
        if source_ids.precursor_canonical_to_source.is_empty()
            && source_ids.precursor_source_to_canonical.is_empty()
        {
            return decoy_ref.strip_prefix(decoy_prefix).map(str::to_string);
        }
        None
    }

    /// Report whether an experiment satisfies the canonical-id invariant:
    /// `true` exactly when [`Self::validate_canonical_ids`] would accept it.
    pub fn has_canonical_ids(exp: &LightTargetedExperiment) -> bool {
        if !Self::has_canonical_id_format(exp) {
            return false;
        }
        let mut precursors: HashSet<&str> = HashSet::with_capacity(exp.compounds.len());
        for compound in &exp.compounds {
            if !precursors.insert(compound.id.as_str()) {
                return false;
            }
        }
        let mut transitions: HashSet<&str> = HashSet::with_capacity(exp.transitions.len());
        for transition in &exp.transitions {
            if !transitions.insert(transition.transition_name.as_str())
                || !precursors.contains(transition.peptide_ref.as_str())
            {
                return false;
            }
        }
        true
    }

    /// Report whether all operational id fields use canonical decimal syntax.
    ///
    /// This checks only the textual representation: every compound id,
    /// transition id and transition precursor reference must be a
    /// non-negative `i64` in canonical decimal form. It deliberately does not
    /// check uniqueness or whether the precursor references resolve; writer
    /// compatibility paths use it to tell source-style identifiers from
    /// malformed canonical-looking input.
    pub fn has_canonical_id_format(exp: &LightTargetedExperiment) -> bool {
        exp.compounds.iter().all(|c| is_canonical_id(&c.id))
            && exp
                .transitions
                .iter()
                .all(|t| is_canonical_id(&t.transition_name) && is_canonical_id(&t.peptide_ref))
    }

    /// Validate that an experiment already uses canonical numeric ids.
    ///
    /// Ids must be unique, non-negative `i64` values written in canonical
    /// decimal form (for example `7`, not `007`, `+7` or `-0`). Ids may be
    /// sparse; filtering must not force renumbering. Every transition
    /// `peptide_ref` must exactly match an existing compound id.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] (the source's
    /// `Exception::InvalidValue`) for the first violated invariant, checking
    /// all compounds before the transitions, as the source does.
    pub fn validate_canonical_ids(exp: &LightTargetedExperiment) -> Result<()> {
        let mut precursors: HashSet<&str> = HashSet::with_capacity(exp.compounds.len());
        for compound in &exp.compounds {
            validate_canonical_id_format(&compound.id, "Precursor")?;
            if !precursors.insert(compound.id.as_str()) {
                return Err(invalid_id(
                    "Canonical precursor IDs must be unique",
                    &compound.id,
                ));
            }
        }
        let mut transitions: HashSet<&str> = HashSet::with_capacity(exp.transitions.len());
        for transition in &exp.transitions {
            validate_canonical_id_format(&transition.transition_name, "Transition")?;
            if !transitions.insert(transition.transition_name.as_str()) {
                return Err(invalid_id(
                    "Canonical transition IDs must be unique",
                    &transition.transition_name,
                ));
            }
            validate_canonical_id_format(
                &transition.peptide_ref,
                "Transition precursor reference",
            )?;
            if !precursors.contains(transition.peptide_ref.as_str()) {
                return Err(invalid_id(
                    "Transition references an unknown canonical precursor ID",
                    &transition.peptide_ref,
                ));
            }
        }
        Ok(())
    }
}
