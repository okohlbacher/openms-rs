// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Hendrik Weisser, OpenMS Rust contributors $

//! Sequence and parent-evidence conversions from IdentificationDataConverter.
//! Complete legacy run conversion, feature-map conversion and mzTab are separate APIs.

use super::{
    GraphWork, IdentificationData, MoleculeType, ParentMatch, ParentMatches, ParentSequence,
};
use crate::format::fasta::FASTAEntry;
use crate::identification::{FlankingResidue, PeptideEvidence, PeptideHit};
use crate::{Error, Result};
use std::mem::size_of;

/// Stateless bridges between sequence files, graph parents and legacy evidence.
pub struct IdentificationDataConverter;
impl IdentificationDataConverter {
    /// Import FASTA entries as parents in caller order, preserving text verbatim.
    /// The source defaults are Protein and an empty decoy pattern. A nonempty
    /// pattern matches any case-sensitive substring of the accession.
    /// Existing parents use normal registration merge and current-step rules.
    /// The entire import is atomic and shares the graph's operation limits.
    pub fn import_sequences(
        graph: &mut IdentificationData,
        entries: &[FASTAEntry],
        molecule_type: MoleculeType,
        decoy_pattern: &str,
    ) -> Result<()> {
        if entries.is_empty() {
            return Ok(());
        }
        graph.transaction(|staged| {
            staged.operation(|_, work| work.consume(add(entries.len(), decoy_pattern.len())?))?;
            for entry in entries {
                let is_decoy = staged.operation(|_, work| {
                    let bytes = add(
                        size_of::<ParentSequence>(),
                        add(
                            entry.identifier.len(),
                            add(entry.sequence.len(), entry.description.len())?,
                        )?,
                    )?;
                    // Charge owned text before cloning and substring comparisons before search.
                    work.copy(bytes)?;
                    if decoy_pattern.is_empty() {
                        Ok(false)
                    } else {
                        work.consume(mul(entry.identifier.len(), decoy_pattern.len())?)?;
                        Ok(entry.identifier.contains(decoy_pattern))
                    }
                })?;
                let mut parent = ParentSequence::new(entry.identifier.clone());
                parent.molecule_type = molecule_type;
                parent.sequence.clone_from(&entry.sequence);
                parent.description.clone_from(&entry.description);
                parent.is_decoy = is_decoy;
                staged.register_parent_sequence(parent)?;
            }
            Ok(())
        })
    }

    /// Append graph parent matches to a hit's existing evidences, then sort by
    /// accession, inclusive start/end and flanking characters. Duplicates remain.
    /// Empty flanks become X; nonempty flanks use only the first source byte.
    /// Unrepresentable legacy flanks, positions beyond signed 32-bit range and
    /// invalid intervals fail before mutation. Unrelated hit payload is untouched.
    pub fn export_parent_matches(
        graph: &IdentificationData,
        parent_matches: &ParentMatches,
        hit: &mut PeptideHit,
    ) -> Result<()> {
        let mut work = GraphWork::new(graph.limits());
        work.consume(add(hit.evidences.len(), parent_matches.len())?)?;
        let mut count = hit.evidences.len();
        let mut payload = 0usize;
        for evidence in &hit.evidences {
            evidence.validate()?;
            legacy_position(evidence.start)?;
            legacy_position(evidence.end)?;
            payload = add(payload, evidence_bytes(&evidence.protein_accession)?)?;
        }
        for (&id, matches) in parent_matches {
            let parent = graph.parent(id)?;
            work.consume(matches.len())?;
            count = add(count, matches.len())?;
            for parent_match in matches {
                checked_context(parent_match)?;
                payload = add(payload, evidence_bytes(&parent.accession)?)?;
            }
        }
        if count > graph.limits().max_edges {
            return Err(invalid(
                "exported peptide evidence count exceeds graph edge limit",
            ));
        }
        // Bound string comparisons, both traversals, all vector slots and cloned
        // accessions before allocating the replacement. Unstable sort needs no scratch.
        let depth = (usize::BITS - count.max(1).leading_zeros()) as usize + 3;
        work.consume(add(mul(payload, depth)?, mul(count, 2)?)?)?;
        work.copy(payload)?;
        let mut result = Vec::new();
        result
            .try_reserve_exact(count)
            .map_err(|_| invalid("cannot allocate exported peptide evidences"))?;
        result.extend(hit.evidences.iter().cloned());
        for (&id, matches) in parent_matches {
            let accession = &graph.parent(id)?.accession;
            for parent_match in matches {
                let (aa_before, aa_after) = checked_context(parent_match)?;
                result.push(PeptideEvidence {
                    protein_accession: accession.clone(),
                    start: parent_match.start_pos,
                    end: parent_match.end_pos,
                    aa_before,
                    aa_after,
                });
            }
        }
        result.sort_unstable();
        hit.evidences = result;
        Ok(())
    }
}

fn checked_context(value: &ParentMatch) -> Result<(FlankingResidue, FlankingResidue)> {
    legacy_position(value.start_pos)?;
    legacy_position(value.end_pos)?;
    if matches!((value.start_pos, value.end_pos), (Some(start), Some(end)) if start > end) {
        return Err(invalid("parent match start exceeds inclusive end"));
    }
    Ok((flank(&value.left_neighbor)?, flank(&value.right_neighbor)?))
}
fn flank(value: &str) -> Result<FlankingResidue> {
    match value.as_bytes().first().copied() {
        None => Ok(FlankingResidue::Unknown),
        Some(byte) => FlankingResidue::from_code(char::from(byte)).map_err(|_| {
            Error::Unsupported(
                "parent flank cannot be represented as a legacy amino-acid character".into(),
            )
        }),
    }
}
fn legacy_position(value: Option<usize>) -> Result<()> {
    if value.is_some_and(|position| position > i32::MAX as usize) {
        Err(Error::Unsupported(
            "parent position exceeds the source legacy signed 32-bit range".into(),
        ))
    } else {
        Ok(())
    }
}
fn evidence_bytes(accession: &str) -> Result<usize> {
    add(size_of::<PeptideEvidence>(), accession.len())
}
fn add(left: usize, right: usize) -> Result<usize> {
    left.checked_add(right)
        .ok_or_else(|| invalid("identification conversion size overflow"))
}
fn mul(left: usize, right: usize) -> Result<usize> {
    left.checked_mul(right)
        .ok_or_else(|| invalid("identification conversion work overflow"))
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
