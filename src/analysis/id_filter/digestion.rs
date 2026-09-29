// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The digestion predicates of `PROCESSING/ID/IDFilter.h` (core `bc9cc12`):
//! `PeptideDigestionFilter` (missed cleavages of a peptide hit) and
//! `DigestionFilter` (is a peptide evidence a digestion product of its
//! protein).

use super::predicates::{GetMatchingItems, try_filter_peptide_evidences};
use crate::Result;
use crate::chemistry::{AASequence, ProductValidation, ProteaseDigestion};
use crate::format::fasta::FASTAEntry;
use crate::identification::{PeptideEvidence, PeptideHit, PeptideIdentification};

/// Filters peptide hits by their number of missed cleavages
/// (`PeptideDigestionFilter`).
///
/// The source's `-1` "disabled" bound (`disabledValue()`) is `None`.
#[derive(Clone, Copy, Debug)]
pub struct PeptideDigestionFilter<'a> {
    /// The digestion whose enzyme counts the cleavage sites.
    pub digestion: &'a ProteaseDigestion,
    /// Inclusive minimum number of missed cleavages, if any.
    pub min_cleavages: Option<i32>,
    /// Inclusive maximum number of missed cleavages, if any.
    pub max_cleavages: Option<i32>,
}

impl<'a> PeptideDigestionFilter<'a> {
    /// A filter over `digestion` with inclusive bounds; `None` disables one.
    pub fn new(
        digestion: &'a ProteaseDigestion,
        min_cleavages: Option<i32>,
        max_cleavages: Option<i32>,
    ) -> Self {
        Self {
            digestion,
            min_cleavages,
            max_cleavages,
        }
    }

    /// The number of missed cleavages the source computes for an unmodified
    /// sequence (`countInternalCleavageSites`, cast to `Int` by
    /// `filterByMissedCleavages`).
    ///
    /// This is the source's token count, not
    /// [`ProteaseDigestion::count_internal_cleavage_sites`]: a peptide that
    /// starts at a lookahead-only site (Asp-N, Lys-N) counts one more missed
    /// cleavage.
    fn missed_cleavages(&self, sequence: &str) -> Result<i64> {
        self.digestion.enzyme.source_token_cleavage_count(sequence)
    }

    /// Whether the hit should be **removed** (the source's `operator()`, for
    /// use with `remove_if`): its missed cleavages exceed an enabled maximum
    /// or fall below an enabled minimum.
    ///
    /// # Errors
    ///
    /// Those of [`ProteaseDigestion::count_internal_cleavage_sites`] for the
    /// unmodified sequence (for example a sequence longer than the digestion
    /// ceiling).
    pub fn matches(&self, hit: &PeptideHit) -> Result<bool> {
        let missed = self.missed_cleavages(hit.sequence.as_str())?;
        let above = self
            .max_cleavages
            .is_some_and(|max| missed > i64::from(max));
        let below = self
            .min_cleavages
            .is_some_and(|min| missed < i64::from(min));
        Ok(above || below)
    }

    /// Removes the hits [`Self::matches`] flags, keeping the others in order
    /// (`filterPeptideSequences`).
    ///
    /// # Errors
    ///
    /// As [`Self::matches`]; nothing is removed then.
    pub fn filter_peptide_sequences(&self, hits: &mut Vec<PeptideHit>) -> Result<()> {
        super::predicates::try_remove_matching_items(hits, |hit| self.matches(hit))
    }
}

/// Keeps the peptide evidences that are digestion products of their protein
/// (`DigestionFilter`).
#[derive(Clone, Debug)]
pub struct DigestionFilter<'a> {
    /// Protein sequences by accession (`accession_resolver_`).
    pub accession_resolver: GetMatchingItems<'a>,
    /// The digestion defining valid products.
    pub digestion: &'a ProteaseDigestion,
    /// Accept products regardless of missed cleavages.
    pub ignore_missed_cleavages: bool,
    /// Accept a start at residue 1 or 2 of a protein beginning with M
    /// (the source's `methionine_cleavage`, passed as
    /// `allow_nterm_protein_cleavage`).
    pub methionine_cleavage: bool,
}

impl<'a> DigestionFilter<'a> {
    /// A filter over the protein `entries`, indexed by identifier.
    pub fn new(
        entries: &'a [FASTAEntry],
        digestion: &'a ProteaseDigestion,
        ignore_missed_cleavages: bool,
        methionine_cleavage: bool,
    ) -> Self {
        Self {
            accession_resolver: GetMatchingItems::new(entries),
            digestion,
            ignore_missed_cleavages,
            methionine_cleavage,
        }
    }

    /// Whether the evidence is **kept** (the source's `operator()`).
    ///
    /// An evidence without valid limits (unknown start or end, or end at the
    /// source's N-terminal position 0) is kept, as is one whose protein is not
    /// among the entries; the source logs a warning for both. Otherwise the
    /// source asks `isValidProduct(protein, start, end - start, ...)`. The
    /// evidence end is inclusive, so that length leaves the peptide's last
    /// residue out and the C-terminal check looks at the bond *before* it; the
    /// port keeps this, checking the range `start..end`.
    ///
    /// # Errors
    ///
    /// When the protein sequence does not parse as an [`AASequence`] (the
    /// source throws from `AASequence::fromString`) or the digestion rejects
    /// it.
    pub fn matches(&self, evidence: &PeptideEvidence) -> Result<bool> {
        let (Some(start), Some(end)) = (evidence.start, evidence.end) else {
            return Ok(true);
        };
        if end == 0 {
            return Ok(true);
        }
        if !self.accession_resolver.exists(evidence) {
            return Ok(true);
        }
        let entry = self.accession_resolver.value(evidence)?;
        let protein: AASequence = entry.sequence.parse()?;
        // Source `isValidProduct_` applies the methionine allowance before its
        // bounds and emptiness checks, so `start..end` with start <= 2 becomes
        // `0..end` even when it was empty; the native `is_valid_product`
        // rejects the empty range first, so the source order is applied here.
        let start = if self.methionine_cleavage && start <= 2 && protein.as_str().starts_with('M') {
            0
        } else {
            start
        };
        self.digestion.is_valid_product(
            &protein,
            start..end,
            ProductValidation {
                ignore_missed_cleavages: self.ignore_missed_cleavages,
                allow_nterm_protein_cleavage: self.methionine_cleavage,
                allow_random_asp_pro_cleavage: false,
            },
        )
    }

    /// Removes, from every hit, the evidences [`Self::matches`] rejects
    /// (`filterPeptideEvidences`).
    ///
    /// # Errors
    ///
    /// As [`Self::matches`]; nothing is removed then.
    pub fn filter_peptide_evidences(&self, ids: &mut [PeptideIdentification]) -> Result<()> {
        try_filter_peptide_evidences(ids, |evidence| self.matches(evidence))
    }
}
