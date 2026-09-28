// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Detected proteins and peptides as a proteomics QC metric (source
//! `QC/IdentificationSummary.h`): a summary of the identifications of an
//! idXML file. See `docs/QC_SUPPORT.md`.

use super::missed_cleavages::MissedCleavages;
use super::qc_base::{QCBase, Requires, Status};
use crate::identification::{PeptideIdentification, ProteinIdentification};
use crate::{Error, Result};
use std::collections::BTreeSet;

/// The number of unique identifications (by sequence or accession only) and
/// their significance threshold.
///
/// Source `IdentificationSummary::UniqueID`. The default is a count of 0 and
/// the threshold -1.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniqueId {
    /// The number of unique identifications.
    pub count: u32,
    /// The significance threshold when the score type is `FDR`, else -1.
    pub fdr_threshold: f32,
}

impl Default for UniqueId {
    fn default() -> Self {
        Self {
            count: 0,
            fdr_threshold: -1.0,
        }
    }
}

/// The summary values [`IdentificationSummary::compute`] calculates.
///
/// Source `IdentificationSummary::Result`; equality compares every field,
/// source `operator==`. A mean over nothing is NaN, as the source's division
/// by zero gives.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IdentificationSummaryResult {
    /// The total number of peptide-spectrum matches: peptide identifications
    /// that are not empty in the source's sense.
    pub peptide_spectrum_matches: u32,
    /// Unique peptides, by the unmodified sequence of each first hit.
    pub unique_peptides: UniqueId,
    /// Unique proteins, by accession.
    pub unique_proteins: UniqueId,
    /// The mean number of missed cleavages of the first hits.
    pub missed_cleavages_mean: f32,
    /// The mean score of all protein hits.
    pub protein_hit_scores_mean: f64,
    /// The mean length of the unique peptide sequences.
    pub peptide_length_mean: f64,
}

/// Detected proteins and peptides as a proteomics QC metric.
///
/// Source `IdentificationSummary`, a stateless metric that summarises the
/// identifications of an idXML file.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct IdentificationSummary;

impl IdentificationSummary {
    /// Source `getName`.
    pub const NAME: &'static str = "Summary of detected Proteins and Peptides from idXML file";

    /// A new metric, source `IdentificationSummary()`.
    pub fn new() -> Self {
        Self
    }

    /// Summarise an idXML file's identifications, source `compute`.
    ///
    /// - `peptide_spectrum_matches`: the peptide identifications that are not
    ///   empty in the source's sense (`PeptideIdentification::empty`: no
    ///   identifier, no hits, significance threshold 0, no score type and
    ///   higher-score-better set), so an identification without hits but with
    ///   a score type still counts;
    /// - `unique_peptides.count` and `peptide_length_mean`: the distinct
    ///   unmodified sequences of the first hits and their mean length;
    /// - `missed_cleavages_mean`: from [`MissedCleavages::compute`], the mean
    ///   number of missed cleavages over the first hits, with the enzyme of the
    ///   first protein identification's search. That call annotates each first
    ///   hit with the meta value `missed_cleavages`, which is why `pep_ids` is
    ///   mutable, as in the source;
    /// - `unique_proteins.count` and `protein_hit_scores_mean`: the distinct
    ///   accessions of all protein hits and the mean of all their scores;
    /// - the `fdr_threshold`s: the significance threshold of the first peptide
    ///   (protein) identification when its score type is exactly `FDR`, else
    ///   -1. The source assumes, as idXML files have it, that every
    ///   identification shares the first one's threshold.
    ///
    /// The missed-cleavage sum and count are accumulated in `u32` with
    /// wrapping, as the source's `UInt` arithmetic; this matters only for the
    /// 4294967295 missed cleavages [`MissedCleavages`] reports for an empty
    /// sequence.
    ///
    /// # Errors
    ///
    /// As [`MissedCleavages::compute`]: [`Error::MissingInformation`] when
    /// `prot_ids` is empty or names no enzyme. Also
    /// [`Error::MissingInformation`] when `pep_ids` is empty, where the source
    /// reads `pep_ids.front()` of an empty list (undefined behaviour).
    pub fn compute(
        &self,
        prot_ids: &[ProteinIdentification],
        pep_ids: &mut [PeptideIdentification],
    ) -> Result<IdentificationSummaryResult> {
        let mut result = IdentificationSummaryResult::default();
        let mut peptides = BTreeSet::new();
        for pep_id in pep_ids.iter() {
            if !source_empty(pep_id) {
                result.peptide_spectrum_matches = result.peptide_spectrum_matches.saturating_add(1);
                if let Some(hit) = pep_id.hits.first() {
                    peptides.insert(hit.sequence.as_str().to_owned());
                }
            }
        }
        let peptide_length_sum: u64 = peptides.iter().map(|p| p.len() as u64).sum();
        result.peptide_length_mean = peptide_length_sum as f64 / peptides.len() as f64;

        let mut mc = MissedCleavages::new();
        mc.compute(prot_ids, pep_ids)?;
        let (mut missed_cleavages, mut pep_count) = (0_u32, 0_u32);
        for map in mc.results() {
            for (&key, &val) in map {
                missed_cleavages = missed_cleavages.wrapping_add(key.wrapping_mul(val));
                pep_count = pep_count.wrapping_add(val);
            }
        }
        result.missed_cleavages_mean = (f64::from(missed_cleavages) / f64::from(pep_count)) as f32;

        let mut proteins = BTreeSet::new();
        let (mut protein_hit_scores_sum, mut protein_hit_count) = (0.0_f64, 0_u64);
        for prot_id in prot_ids {
            protein_hit_count += prot_id.hits.len() as u64;
            for hit in &prot_id.hits {
                proteins.insert(hit.accession.as_str());
                protein_hit_scores_sum += hit.score;
            }
        }
        result.protein_hit_scores_mean = protein_hit_scores_sum / protein_hit_count as f64;

        result.unique_peptides.count = saturating_count(peptides.len());
        result.unique_proteins.count = saturating_count(proteins.len());
        let first_peptide = pep_ids.first().ok_or_else(|| {
            Error::MissingInformation(
                "IdentificationSummary needs at least one peptide identification".into(),
            )
        })?;
        if first_peptide.score_type == "FDR" {
            result.unique_peptides.fdr_threshold = first_peptide.significance_threshold as f32;
        }
        // `mc.compute` has already refused an empty `prot_ids`.
        if let Some(first_protein) = prot_ids.first() {
            if first_protein.score_type == "FDR" {
                result.unique_proteins.fdr_threshold = first_protein.significance_threshold as f32;
            }
        }
        Ok(result)
    }
}

/// Source `PeptideIdentification::empty`, which is narrower than the crate's
/// [`PeptideIdentification::is_empty`] (no hits): it also requires an empty
/// identifier and score type, a zero significance threshold and
/// higher-score-better set.
fn source_empty(pep_id: &PeptideIdentification) -> bool {
    pep_id.identifier.is_empty()
        && pep_id.hits.is_empty()
        && pep_id.significance_threshold == 0.0
        && pep_id.score_type.is_empty()
        && pep_id.higher_score_better
}

fn saturating_count(n: usize) -> u32 {
    u32::try_from(n).unwrap_or(u32::MAX)
}

impl QCBase for IdentificationSummary {
    fn name(&self) -> &str {
        Self::NAME
    }
    /// An idXML file with protein identifications, source `requirements`.
    fn requirements(&self) -> Status {
        Status::from_flag(Requires::Id)
    }
}
