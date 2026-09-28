// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Missed cleavages as a QC metric for the QualityControl TOPP tool (source
//! `QC/MissedCleavages.h`). See `docs/QC_SUPPORT.md`.

use super::qc_base::{QCBase, Requires, Status, has_pep_id_feature_map};
use crate::chemistry::ProteaseDigestion;
use crate::identification::{PeptideIdentification, ProteinIdentification};
use crate::kernel::FeatureMap;
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::collections::BTreeMap;

/// Counts the missed cleavages of each peptide identification.
///
/// Source `MissedCleavages`. Each call to a `compute` method digests the
/// sequence of the first hit of every peptide identification with the
/// search's enzyme and no missed cleavages allowed; the number of fragments
/// minus one is that hit's number of missed cleavages. The first hit is
/// annotated with the meta value `missed_cleavages`, and one map
/// "missed cleavages -> occurrences" per call is appended to
/// [`MissedCleavages::results`].
///
/// The source class documentation also lists `FWHM` and `mass` annotations;
/// the source code writes neither, and neither does this port.
///
/// Warnings the source writes to `OPENMS_LOG_WARN` are collected in
/// [`MissedCleavages::log`], in order, without line ends.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MissedCleavages {
    mc_result: Vec<BTreeMap<u32, u32>>,
    /// The warning lines the source logs: a peptide identification without
    /// hits, an observed count above the search's maximum, an empty feature
    /// map.
    pub log: Vec<String>,
}

impl MissedCleavages {
    /// Source `getName`: `MissedCleavages`.
    pub const NAME: &'static str = "MissedCleavages";

    /// A new metric with no results, source `MissedCleavages()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// The results, one map "number of missed cleavages -> count" per call to
    /// a `compute` method, source `getResults`.
    pub fn results(&self) -> &[BTreeMap<u32, u32>] {
        &self.mc_result
    }

    /// Count the missed cleavages of every peptide identification in
    /// `pep_ids`, source `compute(prot_ids, pep_ids)`.
    ///
    /// The enzyme and the allowed maximum come from the search parameters of
    /// the first protein identification. The first hit of each peptide
    /// identification gains the meta value `missed_cleavages`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when `prot_ids` is empty or its
    /// first search names the enzyme `unknown_enzyme`, as the source's
    /// `Exception::MissingInformation`; the enzyme lookup's error when the name
    /// is not in the protease database (source: `ElementNotFound`); and the
    /// digestion's error for a sequence it cannot digest. On error nothing is
    /// appended to the results, though hits annotated before the failing one
    /// keep their annotation, as in the source.
    pub fn compute(
        &mut self,
        prot_ids: &[ProteinIdentification],
        pep_ids: &mut [PeptideIdentification],
    ) -> Result<()> {
        let first = prot_ids.first().ok_or_else(|| {
            Error::MissingInformation("Missing information in ProteinIdentifications.".into())
        })?;
        let (digestor, max_mc) = digestor(
            first,
            "No digestion enzyme in ID data detected. No computation possible.",
        )?;
        let mut result = BTreeMap::new();
        for pep_id in pep_ids.iter_mut() {
            self.count_peptide_identification(&digestor, &mut result, max_mc, pep_id)?;
        }
        self.mc_result.push(result);
        Ok(())
    }

    /// Count the missed cleavages of every peptide identification of `fmap`,
    /// source `compute(FeatureMap&)`: those of each top-level feature in
    /// order, then the unassigned ones (subordinate features are not
    /// visited), with the enzyme of the map's first protein identification.
    ///
    /// A map without any peptide identification appends an empty result, as
    /// does an empty map (with the source's warning, which its earlier check
    /// makes unreachable).
    ///
    /// # Errors
    ///
    /// As [`MissedCleavages::compute`], with the source's feature-map wording
    /// for a missing enzyme.
    pub fn compute_feature_map(&mut self, fmap: &mut FeatureMap) -> Result<()> {
        if !has_pep_id_feature_map(fmap) {
            self.mc_result.push(BTreeMap::new());
            return Ok(());
        }
        if fmap.features.is_empty() {
            self.log.push("FeatureXML is empty.".into());
            self.mc_result.push(BTreeMap::new());
            return Ok(());
        }
        let first = fmap.protein_identifications.first().ok_or_else(|| {
            Error::MissingInformation("Missing information in ProteinIdentifications.".into())
        })?;
        let (digestor, max_mc) = digestor(
            first,
            "No digestion enzyme in FeatureMap detected. No computation possible.",
        )?;
        let mut result = BTreeMap::new();
        for feature in &mut fmap.features {
            for pep_id in &mut feature.peptide_identifications {
                self.count_peptide_identification(&digestor, &mut result, max_mc, pep_id)?;
            }
        }
        for pep_id in &mut fmap.unassigned_peptide_identifications {
            self.count_peptide_identification(&digestor, &mut result, max_mc, pep_id)?;
        }
        self.mc_result.push(result);
        Ok(())
    }

    /// Source `get_missed_cleavages_from_peptide_identification_`.
    ///
    /// An empty sequence digests to no fragment at all, and the source's
    /// `UInt32(digest_output.size() - 1)` then wraps to 4294967295 missed
    /// cleavages; this port preserves that count so that its means agree with
    /// the source's.
    fn count_peptide_identification(
        &mut self,
        digestor: &ProteaseDigestion,
        result: &mut BTreeMap<u32, u32>,
        max_mc: u32,
        pep_id: &mut PeptideIdentification,
    ) -> Result<()> {
        let Some(hit) = pep_id.hits.first_mut() else {
            self.log.push(format!(
                "There is a Peptideidentification(RT: {}, MZ: {}) without PeptideHits.",
                coordinate(pep_id.rt),
                coordinate(pep_id.mz)
            ));
            return Ok(());
        };
        // The source digests into a vector and takes its size; counting the
        // same products allocates nothing.
        let fragments = digestor.peptide_count(&hit.sequence)?;
        let num_mc = (fragments as u64).wrapping_sub(1) as u32;
        if num_mc > max_mc {
            self.log.push(format!(
                "Observed number of missed cleavages: {num_mc} is greater than: {max_mc} the \
                 allowed maximum number of missed cleavages during MS2-Search in: {}",
                hit.sequence
            ));
        }
        let count = result.entry(num_mc).or_insert(0_u32);
        *count = count.saturating_add(1);
        hit.metadata
            .insert("missed_cleavages".into(), MetaValue::from(num_mc));
        Ok(())
    }
}

/// A digestor for the enzyme of `protein`'s search, allowing no missed
/// cleavage, and the search's allowed maximum.
fn digestor(protein: &ProteinIdentification, missing: &str) -> Result<(ProteaseDigestion, u32)> {
    let parameters = &protein.search_parameters;
    if parameters.digestion_enzyme == "unknown_enzyme" {
        return Err(Error::MissingInformation(missing.into()));
    }
    let mut digestor = ProteaseDigestion::default();
    digestor.set_enzyme(&parameters.digestion_enzyme)?;
    digestor.missed_cleavages = 0;
    Ok((digestor, parameters.missed_cleavages))
}

/// A peptide identification's RT or m/z for the log line: the stored value in
/// Rust's shortest form (the source's stream prints six significant digits),
/// or `nan` when unset, as the source prints its NaN default.
fn coordinate(value: Option<f64>) -> String {
    value.map_or_else(|| "nan".into(), |v| v.to_string())
}

impl QCBase for MissedCleavages {
    fn name(&self) -> &str {
        Self::NAME
    }
    /// Features with FDR-filtered peptide identifications, source
    /// `requirements`.
    fn requirements(&self) -> Status {
        Status::from_flag(Requires::PostFdrFeat)
    }
}
