// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `ANALYSIS/ID/IDMergerAlgorithm.h`: merge protein and peptide
//! identification runs into one run.

use crate::concept::constants::user_param::ID_MERGE_INDEX;
use crate::identification::{PeptideIdentification, ProteinHit, ProteinIdentification};
use crate::metadata::{MetaValue, MetaValueData};
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::{Error, Result};
use chrono::{Datelike, Timelike};
use std::collections::{BTreeMap, BTreeSet};

fn missing(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// Warning appended to every consistency warning, as the source's
/// `ProteinIdentification::peptideIDsMergeable`.
const MERGE_WARNING: &str = " You probably do not want to merge the results with this tool. For merging searches with different engines/settings please use ConsensusID or PercolatorAdapter to create a comparable score.";

/// The parameters of [`IDMergerAlgorithm`] (source handler
/// `IDMergerAlgorithm`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IDMergerParameters {
    /// `annotate_origin`: add an `id_merge_index` meta value to each peptide
    /// identification naming the run file it came from.
    pub annotate_origin: bool,
    /// `allow_disagreeing_settings`: merge runs whose search engine or search
    /// settings disagree anyway. Use at your own risk.
    pub allow_disagreeing_settings: bool,
}

impl Default for IDMergerParameters {
    fn default() -> Self {
        Self {
            annotate_origin: true,
            allow_disagreeing_settings: false,
        }
    }
}

impl IDMergerParameters {
    /// The default parameters, in source declaration order.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters in the layout of [`IDMergerParameters::defaults`],
    /// carrying the current values (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let mut param = Param::new();
        for (key, value, description) in [
            (
                "annotate_origin",
                self.annotate_origin,
                "If true, adds a map_index MetaValue to the PeptideIDs to annotate the IDRun they came from.",
            ),
            (
                "allow_disagreeing_settings",
                self.allow_disagreeing_settings,
                "Force merging of disagreeing runs. Use at your own risk.",
            ),
        ] {
            param.set_value(
                key,
                ParamValue::String(if value { "true" } else { "false" }.into()),
                description,
                &[],
            )?;
            param.set_valid_strings(key, &["true".to_owned(), "false".to_owned()])?;
        }
        Ok(param)
    }

    /// Build from parameters, as source `setParameters` followed by
    /// `updateMembers_`. Missing parameters take their defaults; unknown ones
    /// are ignored (see [`IDMergerParameters::from_param_with_warnings`]).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its restriction (source `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`IDMergerParameters::from_param`], also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`IDMergerParameters::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new("IDMergerAlgorithm")?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |param| {
            Ok(Self {
                annotate_origin: param.value("annotate_origin")?.to_bool()?,
                allow_disagreeing_settings: param.value("allow_disagreeing_settings")?.to_bool()?,
            })
        })
    }
}

/// Merges protein and peptide identification runs (source
/// `IDMergerAlgorithm`).
///
/// A new protein identification run is created, into which runs are
/// inserted: the protein hits of all inserted runs are united (by accession,
/// the first hit of an accession wins), the peptide identifications are
/// concatenated and re-pointed at the new run, the search engine and settings
/// of all runs are checked for consistency, and, with `annotate_origin`, every
/// peptide identification is annotated with the index of the spectrum file it
/// came from in the merged run's primary MS run paths.
///
/// Unlike the `IDMerger` tool's default mode, all peptide hits are kept, not
/// only the first per sequence.
///
/// # Notes
///
/// The source collects protein hits in an `std::unordered_set` keyed by
/// accession and returns them in that set's iteration order, which depends on
/// the standard library's hash function. This returns them in first-insertion
/// order, which is deterministic; the set of hits is the same.
///
/// The source logs a warning and does nothing when a call inserts no protein
/// runs; the warnings of each call are returned instead.
#[derive(Clone, Debug)]
pub struct IDMergerAlgorithm {
    /// The algorithm parameters.
    pub parameters: IDMergerParameters,
    prot_result: ProteinIdentification,
    pep_result: Vec<PeptideIdentification>,
    collected_protein_hits: Vec<ProteinHit>,
    collected_accessions: BTreeSet<String>,
    filled: bool,
    file_origin_to_idx: BTreeMap<String, usize>,
    id: String,
    fixed_identifier: bool,
}

impl Default for IDMergerAlgorithm {
    /// `IDMergerAlgorithm::new("merged", true)`, the source defaults.
    fn default() -> Self {
        Self::new("merged", true)
    }
}

impl IDMergerAlgorithm {
    /// A merger whose run identifier is `run_identifier`, followed by the
    /// local time as `%d-%m-%Y %H-%M-%S` when `add_time_stamp_to_id`, so that
    /// merged runs get unique identifiers.
    pub fn new(run_identifier: &str, add_time_stamp_to_id: bool) -> Self {
        let mut merger = Self {
            parameters: IDMergerParameters::default(),
            prot_result: ProteinIdentification::default(),
            pep_result: Vec::new(),
            collected_protein_hits: Vec::new(),
            collected_accessions: BTreeSet::new(),
            filled: false,
            file_origin_to_idx: BTreeMap::new(),
            id: run_identifier.to_owned(),
            fixed_identifier: !add_time_stamp_to_id,
        };
        merger.prot_result.identifier = merger.new_identifier();
        merger
    }

    /// The identifier the merged run currently carries.
    pub fn identifier(&self) -> &str {
        &self.prot_result.identifier
    }

    /// Source `getNewIdentifier_`.
    fn new_identifier(&self) -> String {
        if self.fixed_identifier {
            return self.id.clone();
        }
        let now = chrono::Local::now();
        format!(
            "{}{:02}-{:02}-{:04} {:02}-{:02}-{:02}",
            self.id,
            now.day(),
            now.month(),
            now.year(),
            now.hour(),
            now.minute(),
            now.second()
        )
    }

    /// Insert runs (source `insertRuns`, both the copying and the moving
    /// overload: pass clones to keep the input).
    ///
    /// Only peptide identifications of runs in `prots` are inserted; others are
    /// skipped silently, which merges the peptides of selected runs only. If
    /// several compatible runs come from the same spectrum file, the file
    /// appears once per run in the origins. Nothing happens when `prots` is
    /// empty.
    ///
    /// The first call checks the runs against the first of them and copies its
    /// search engine, version and parameters into the merged run; later calls
    /// check against the merged run. Mismatches produce warnings, which are
    /// returned; the source checks stop at the first mismatching run. The
    /// experiment type is assumed to be label-free.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] (source `Exception::MissingInformation`)
    /// when the settings disagree and `allow_disagreeing_settings` is off;
    /// when `annotate_origin` is on and a run has no primary MS run path; when
    /// an identification to annotate carries no `id_merge_index` although its
    /// run has several files; or when its `id_merge_index` is not a valid file
    /// index of its run. On error the merger is unchanged; the source may have
    /// inserted part of the input.
    pub fn insert_runs(
        &mut self,
        prots: Vec<ProteinIdentification>,
        peps: Vec<PeptideIdentification>,
    ) -> Result<Vec<String>> {
        if prots.is_empty() {
            return Ok(vec![
                "No ProteinIdentification(Runs) given. Skipping.".to_owned(),
            ]);
        }
        let mut next = self.clone();
        let mut warnings = Vec::new();
        if !next.filled {
            if prots.len() > 1 {
                next.check_run_consistency(&prots, &prots[0], "label-free", &mut warnings)?;
            }
            next.prot_result.search_engine = prots[0].search_engine.clone();
            next.prot_result.search_engine_version = prots[0].search_engine_version.clone();
            next.prot_result.search_parameters = prots[0].search_parameters.clone();
            next.filled = true;
        } else {
            let reference = next.prot_result.clone();
            next.check_run_consistency(&prots, &reference, "label-free", &mut warnings)?;
        }
        next.move_pep_ids_and_ref_proteins_to_result(peps, prots)?;
        *self = next;
        Ok(warnings)
    }

    /// Return the merged protein run and peptide identifications, and reset
    /// the merger for reuse with a new identifier (source
    /// `returnResultsAndClear`).
    ///
    /// The merged run's primary MS run paths are the spectrum files of all
    /// inserted runs, in the order first seen; its hits are the collected
    /// protein hits.
    pub fn return_results_and_clear(
        &mut self,
    ) -> (ProteinIdentification, Vec<PeptideIdentification>) {
        let mut origins = vec![String::new(); self.file_origin_to_idx.len()];
        for (path, &index) in &self.file_origin_to_idx {
            origins[index] = path.clone();
        }
        let mut prots = std::mem::take(&mut self.prot_result);
        prots.primary_ms_run_paths = origins;
        prots.hits.append(&mut self.collected_protein_hits);
        let peps = std::mem::take(&mut self.pep_result);
        self.prot_result.identifier = self.new_identifier();
        self.file_origin_to_idx.clear();
        self.collected_accessions.clear();
        self.filled = false;
        (prots, peps)
    }

    /// Source `checkOldRunConsistency_`: every run against `reference`.
    fn check_run_consistency(
        &self,
        runs: &[ProteinIdentification],
        reference: &ProteinIdentification,
        experiment_type: &str,
        warnings: &mut Vec<String>,
    ) -> Result<bool> {
        let mut ok = true;
        for run in runs {
            // The source's `ok = ok && mergeable(...)` stops checking, and
            // warning, after the first mismatch.
            if ok {
                ok = peptide_ids_mergeable(reference, run, experiment_type, warnings)?;
            }
        }
        if !ok && !self.parameters.allow_disagreeing_settings {
            return Err(missing(
                "Search settings are not matching across IdentificationRuns. See warnings. Aborting..",
            ));
        }
        Ok(ok)
    }

    /// Source `movePepIDsAndRefProteinsToResultFaster_`: merges peptides and
    /// proteins without checking that the proteins a peptide references exist.
    fn move_pep_ids_and_ref_proteins_to_result(
        &mut self,
        peps: Vec<PeptideIdentification>,
        prots: Vec<ProteinIdentification>,
    ) -> Result<()> {
        let annotate_origin = self.parameters.annotate_origin;
        let mut origin_files: Vec<Vec<String>> = Vec::with_capacity(prots.len());
        for run in &prots {
            if run.primary_ms_run_paths.is_empty() && annotate_origin {
                return Err(missing(format!(
                    "Annotation of origin requested during merge, but no origin present in run {}.",
                    run.identifier
                )));
            }
            for path in &run.primary_ms_run_paths {
                let next = self.file_origin_to_idx.len();
                self.file_origin_to_idx.entry(path.clone()).or_insert(next);
            }
            origin_files.push(run.primary_ms_run_paths.clone());
        }
        let mut run_id_to_run_idx = BTreeMap::new();
        for (index, run) in prots.iter().enumerate() {
            run_id_to_run_idx.insert(run.identifier.clone(), index);
        }
        self.update_and_move_pep_ids(peps, &run_id_to_run_idx, &origin_files, annotate_origin)?;
        for run in prots {
            for hit in run.hits {
                if self.collected_accessions.insert(hit.accession.clone()) {
                    self.collected_protein_hits.push(hit);
                }
            }
        }
        Ok(())
    }

    /// Source `updateAndMovePepIDs_`.
    fn update_and_move_pep_ids(
        &mut self,
        peps: Vec<PeptideIdentification>,
        run_id_to_run_idx: &BTreeMap<String, usize>,
        origin_files: &[Vec<String>],
        annotate_origin: bool,
    ) -> Result<()> {
        for mut pid in peps {
            let Some(&run_index) = run_id_to_run_idx.get(&pid.identifier) else {
                // An easy way to merge the peptides of certain runs only.
                continue;
            };
            let annotated = pid.metadata.contains_key(ID_MERGE_INDEX);
            if annotate_origin || annotated {
                let origins = &origin_files[run_index];
                let position = describe_position(&pid);
                let old_file_index = if annotated {
                    merge_index(&pid.metadata[ID_MERGE_INDEX], &position)?
                } else if origins.len() > 1 {
                    return Err(missing(format!(
                        "Trying to annotate new id_merge_index for PeptideIdentification ({position}) butno old id_merge_index present"
                    )));
                } else {
                    Some(0)
                };
                let file = old_file_index
                    .and_then(|index| origins.get(index))
                    .ok_or_else(|| {
                        missing(format!(
                            "Trying to annotate new id_merge_index for PeptideIdentification ({position}) but the index exceeds the number of files in the run."
                        ))
                    })?;
                let new_index = self.file_origin_to_idx[file];
                let new_index = i64::try_from(new_index)
                    .map_err(|_| missing("id_merge_index exceeds the integer meta value range"))?;
                pid.metadata
                    .insert(ID_MERGE_INDEX.into(), MetaValue::from(new_index));
            }
            pid.identifier = self.prot_result.identifier.clone();
            self.pep_result.push(pid);
        }
        Ok(())
    }
}

/// Source `(Size)(Int)getMetaValue(ID_MERGE_INDEX)`: an integer, or a float
/// truncated toward zero; a negative index is `None` (the source's unsigned
/// cast makes it too large for any run).
///
/// # Errors
///
/// A non-numeric value (source conversion error).
pub(super) fn merge_index(value: &MetaValue, position: &str) -> Result<Option<usize>> {
    let index = match value.data() {
        MetaValueData::Integer(v) => i64::from(i32::try_from(*v).unwrap_or(-1)),
        MetaValueData::Float(v) if v.is_finite() && v.abs() < 2_147_483_648.0 => v.trunc() as i64,
        _ => {
            return Err(missing(format!(
                "id_merge_index of PeptideIdentification ({position}) is not an integer"
            )));
        }
    };
    Ok(usize::try_from(index).ok())
}

/// `m/z, RT` as the source messages print them; unknown values print as the
/// source's NaN.
pub(super) fn describe_position(pid: &PeptideIdentification) -> String {
    let show = |v: Option<f64>| v.map_or_else(|| "nan".to_owned(), |v| v.to_string());
    format!("{}, {}", show(pid.mz), show(pid.rt))
}

/// Source `ProteinIdentification::peptideIDsMergeable`, with its warnings.
pub(super) fn peptide_ids_mergeable(
    reference: &ProteinIdentification,
    run: &ProteinIdentification,
    experiment_type: &str,
    warnings: &mut Vec<String>,
) -> Result<bool> {
    let mut ok = true;
    if run.search_engine != reference.search_engine
        || run.search_engine_version != reference.search_engine_version
    {
        ok = false;
        warnings.push(format!(
            "Search engine {}from IDRun {} does not match with the others.{MERGE_WARNING}",
            run.search_engine, run.identifier
        ));
    }
    if !reference
        .search_parameters
        .mergeable(&run.search_parameters, experiment_type)?
    {
        ok = false;
        warnings.push(format!(
            "Searchengine settings or modifications from IDRun {} do not match with the others.{MERGE_WARNING}",
            run.identifier
        ));
    }
    Ok(ok)
}
