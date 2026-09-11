// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Protein coverage by exact, overlapping unmodified peptide occurrences.
//! Source: OpenMS Core SDK 54a232f SequenceCoverage. Limits and native allocation
//! differences are documented in `docs/SEQUENCE_COVERAGE_SUPPORT.md`.

use super::AASequence;
use crate::{Error, Result};

pub const MAX_SEQUENCE_COVERAGE_RESIDUES: usize = 1_000_000;
pub const MAX_SEQUENCE_COVERAGE_PEPTIDES: usize = 1_000_000;
pub const MAX_SEQUENCE_COVERAGE_WORK: usize = 50_000_000;

/// Source coverage is a percentage, not a fraction. Annotations are ignored;
/// residue symbols, including B/Z/X, are matched literally without wildcards.
pub struct SequenceCoverage;

impl SequenceCoverage {
    /// Union all exact occurrences, including overlaps, of every nonempty
    /// peptide. Empty protein or empty peptide input returns zero immediately.
    /// Queries borrow all sequences and never inspect their chemical metadata.
    pub fn get_coverage(protein: &AASequence, peptides: &[AASequence]) -> Result<f64> {
        let mut work = MAX_SEQUENCE_COVERAGE_WORK;
        coverage(protein, peptides, &mut work)
    }
}

fn coverage(protein: &AASequence, peptides: &[AASequence], work: &mut usize) -> Result<f64> {
    if protein.is_empty() || peptides.is_empty() {
        return Ok(0.0);
    }
    let protein = protein.as_str().as_bytes();
    if protein.len() > MAX_SEQUENCE_COVERAGE_RESIDUES
        || peptides.len() > MAX_SEQUENCE_COVERAGE_PEPTIDES
    {
        return Err(limit());
    }
    // Precharge both peptide-list passes, coverage initialization and final sum.
    consume(work, 2 * (peptides.len() + protein.len()))?;
    let mut total_residues = 0usize;
    for peptide in peptides {
        let length = peptide.len();
        total_residues = total_residues
            .checked_add(length)
            .filter(|&n| n <= MAX_SEQUENCE_COVERAGE_RESIDUES)
            .ok_or_else(limit)?;
        if length == 0 || length > protein.len() {
            continue;
        }
        let windows = protein.len() - length + 1;
        // Byte-slice equality may return early, but every comparison is charged
        // its full possible length. This deliberate bound is data-independent.
        consume(work, windows.checked_mul(length + 1).ok_or_else(limit)?)?;
    }
    let mut covered = Vec::new();
    covered
        .try_reserve_exact(protein.len())
        .map_err(|_| limit())?;
    covered.resize(protein.len(), false);
    for peptide in peptides {
        let peptide = peptide.as_str().as_bytes();
        if peptide.is_empty() || peptide.len() > protein.len() {
            continue;
        }
        for (position, window) in protein.windows(peptide.len()).enumerate() {
            if window == peptide {
                // Repeated/overlapping writes are charged even when already true.
                consume(work, peptide.len())?;
                covered[position..position + peptide.len()].fill(true);
            }
        }
    }
    let covered_count = covered.into_iter().filter(|&value| value).count();
    Ok(covered_count as f64 * 100.0 / protein.len() as f64)
}
fn consume(work: &mut usize, amount: usize) -> Result<()> {
    *work = work.checked_sub(amount).ok_or_else(limit)?;
    Ok(())
}
fn limit() -> Error {
    Error::InvalidValue("sequence coverage resource limit exceeded".into())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn shared_work_covers_search_and_repeated_covered_span_writes() {
        let protein = AASequence::parse("AAAA").unwrap();
        let peptides = [AASequence::parse("AAA").unwrap()];
        // 2*(1+4) list/bitmap units, 2*(3+1) comparison units, 2*3 writes.
        let mut exact = 24;
        assert_eq!(coverage(&protein, &peptides, &mut exact).unwrap(), 100.0);
        assert_eq!(exact, 0);
        let mut late = 23;
        assert!(coverage(&protein, &peptides, &mut late).is_err());
        assert_eq!(late, 2); // First match written; second rejected before writes.
        assert_eq!(protein.as_str(), "AAAA");
        let mut early = 17;
        assert!(coverage(&protein, &peptides, &mut early).is_err());
        assert_eq!(early, 7); // Comparison upper bound fails before allocation.
    }
    #[test]
    fn empty_source_shortcuts_do_not_consume_work() {
        let empty = AASequence::default();
        let protein = AASequence::parse("A").unwrap();
        let mut work = 0;
        assert_eq!(
            coverage(&empty, std::slice::from_ref(&protein), &mut work).unwrap(),
            0.0
        );
        assert_eq!(coverage(&protein, &[], &mut work).unwrap(), 0.0);
        assert!(coverage(&protein, &[empty], &mut work).is_err());
    }
}
