// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `ANALYSIS/ID/ConsensusMapMergerAlgorithm.h`: merge the protein
//! identification runs of a consensus map, per condition or into one.

use super::id_merger::{describe_position, merge_index, peptide_ids_mergeable};
use crate::concept::constants::user_param::ID_MERGE_INDEX;
use crate::identification::{PeptideIdentification, ProteinHit, ProteinIdentification};
use crate::kernel::ConsensusMap;
use crate::metadata::{ExperimentalDesign, MetaValue};
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// Source `Exception::MissingInformation`.
fn missing(message: impl Into<String>) -> Error {
    Error::MissingInformation(message.into())
}

/// Merges the identification runs of a [`ConsensusMap`] (source
/// `ConsensusMapMergerAlgorithm`).
///
/// Protein hits are united by accession, the first hit of an accession
/// winning; peptide identifications are re-pointed at the new runs and their
/// `id_merge_index` shifted to index the new runs' primary MS run paths.
/// Protein groups are not carried over.
///
/// # Notes
///
/// The source collects protein hits in an unordered container keyed by
/// accession and emits them in its hash order; this emits them in
/// first-insertion order, which is deterministic. Every method returns the
/// warnings the source logs. On error the map is unchanged; the source may
/// have modified it partway.
///
/// The source's to-do notes: this could be folded into
/// [`IDMergerAlgorithm`](super::id_merger::IDMergerAlgorithm), and it is
/// untested for TMT/iTRAQ data, where one identification run covers several
/// multiplexed conditions and peptide identifications would need duplicating
/// per condition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ConsensusMapMergerAlgorithm {
    /// `annotate_origin`: declared by the source ("adds a map_index MetaValue
    /// to the PeptideIDs to annotate the IDRun they came from") but never read
    /// by it; kept for the parameter interface only.
    pub annotate_origin: bool,
    /// Reproduce the source's `id_merge_index` offsets in
    /// [`merge_protein_id_runs`](Self::merge_protein_id_runs).
    ///
    /// The source looks the offset of a peptide's old run inside its new
    /// condition run up under the peptide's *new* run identifier, which it has
    /// just assigned, so the lookup misses and every offset is 0: peptides of
    /// the second and later runs merged into a condition then point at the
    /// first run's files. By default the offset is looked up under the old run
    /// identifier, as intended; `true` selects the source result. Not a source
    /// parameter.
    pub source_merge_index_offsets: bool,
}

impl Default for ConsensusMapMergerAlgorithm {
    fn default() -> Self {
        Self {
            annotate_origin: true,
            source_merge_index_offsets: false,
        }
    }
}

impl ConsensusMapMergerAlgorithm {
    /// The default parameters, as the source declares them.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters in the layout of
    /// [`ConsensusMapMergerAlgorithm::defaults`], carrying the current values
    /// (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let mut param = Param::new();
        param.set_value(
            "annotate_origin",
            ParamValue::String(if self.annotate_origin { "true" } else { "false" }.into()),
            "If true, adds a map_index MetaValue to the PeptideIDs to annotate the IDRun they came from.",
            &[],
        )?;
        param.set_valid_strings("annotate_origin", &["true".to_owned(), "false".to_owned()])?;
        Ok(param)
    }

    /// Build from parameters, as source `setParameters`. Missing parameters
    /// take their defaults; unknown ones are ignored (see
    /// [`ConsensusMapMergerAlgorithm::from_param_with_warnings`]).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its restriction.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`ConsensusMapMergerAlgorithm::from_param`], also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`ConsensusMapMergerAlgorithm::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new("ConsensusMapMergerAlgorithm")?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |param| {
            Ok(Self {
                annotate_origin: param.value("annotate_origin")?.to_bool()?,
                ..Self::default()
            })
        })
    }

    /// Merge the map's identification runs (usually one per column) into one
    /// run per condition of `design`, i.e. per sample-section entry once
    /// replicate columns are removed (source
    /// `mergeProteinsAcrossFractionsAndReplicates`).
    ///
    /// Each column header is matched to a condition by its file name and
    /// channel (`channel_id` + 1, or 1 when absent); the mapping is then
    /// applied by [`merge_protein_id_runs`](Self::merge_protein_id_runs).
    /// Without a sample section that separates conditions, all runs count as
    /// replicates or fractions and are merged into one.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when a column header matches no
    /// condition (source `Exception::MissingInformation`),
    /// [`Error::InvalidValue`] when the design cannot produce its condition
    /// mapping, or as
    /// [`merge_protein_id_runs`](Self::merge_protein_id_runs).
    pub fn merge_proteins_across_fractions_and_replicates(
        &self,
        map: &mut ConsensusMap,
        design: &ExperimentalDesign,
    ) -> Result<Vec<String>> {
        let to_merge = design.condition_to_path_label_vector()?;
        let label_free = map.experiment_type == "label-free";
        let mut warnings = Vec::new();
        if !label_free {
            warnings.push("Merging untested for labelled experiments".to_owned());
        }
        let mut map_index_to_condition = BTreeMap::new();
        for (&column, header) in &map.column_headers {
            if !header.metadata.contains_key("channel_id") && !label_free {
                warnings.push(
                    "No channel id annotated in consensusXML. Assuming one channel.".to_owned(),
                );
            }
            let label = header.label_as_uint(&map.experiment_type)?;
            let found = to_merge.iter().position(|batch| {
                batch
                    .iter()
                    .any(|(path, lab)| *path == header.filename && *lab == label)
            });
            let Some(condition) = found else {
                return Err(missing(format!(
                    "ConsensusHeader entry ({}, {}) could not be matched to the given experimental design.",
                    header.filename, header.label
                )));
            };
            let column = u32::try_from(column).map_err(|_| bad("column index exceeds 32 bits"))?;
            let condition =
                u32::try_from(condition).map_err(|_| bad("condition index exceeds 32 bits"))?;
            map_index_to_condition.insert(column, condition);
        }
        warnings.extend(self.merge_protein_id_runs(map, &map_index_to_condition)?);
        Ok(warnings)
    }

    /// Merge the map's identification runs by a mapping from column index
    /// (the peptides' `map_index`) to new run index (source
    /// `mergeProteinIDRuns`).
    ///
    /// Each new run `i` collects the old runs whose primary MS run paths are
    /// all among the files of the columns mapped to it, is named
    /// `condition<i>`, takes the search engine and settings of its first old
    /// run and the concatenated primary MS run paths of all of them. Protein
    /// hits go to every new run their old run contributes to. Peptide
    /// identifications (assigned and unassigned) are re-pointed at their old
    /// run's new run (the last one, if several) and their `id_merge_index` is
    /// shifted by the position of their old run's files in the new run; see
    /// [`source_merge_index_offsets`](Self::source_merge_index_offsets).
    ///
    /// # Notes
    ///
    /// The source visits the old runs in the order of their identifiers (a
    /// `std::map`). Readers generate those identifiers with a random or
    /// process-global unique number, so for runs of equal search engine and
    /// date that order, and with it the order of a new run's primary MS run
    /// paths and the peptides' `id_merge_index` offsets, is not reproducible.
    /// This visits them in the order they are first matched: new runs in
    /// index order, old runs in map order.
    ///
    /// A run identifier seen again, either a duplicate or a run contributing
    /// to a second new run, is renamed by appending `_` and the byte `0x01`:
    /// the source appends a `Size` 1 to a `std::string`, which converts it to
    /// that character.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a mapped column has no header, the
    /// number of new runs is at least the number of columns (and not one:
    /// nothing would be merged), a new run index is not below the number of
    /// new runs for a run that contributes to it (the source reads out of
    /// bounds), a run renamed as a duplicate can no longer be found under its
    /// old identifier (the source dereferences an end iterator), or a peptide
    /// identification belongs to no merged run or has a negative or
    /// non-numeric `id_merge_index`.
    pub fn merge_protein_id_runs(
        &self,
        map: &mut ConsensusMap,
        map_index_to_new_run: &BTreeMap<u32, u32>,
    ) -> Result<Vec<String>> {
        let mut warnings = Vec::new();
        if map.experiment_type != "label-free" {
            warnings.push("Merging untested for labelled experiments".to_owned());
        }
        let mut new_indices: BTreeMap<u32, BTreeSet<String>> = BTreeMap::new();
        for (&column, &new_run) in map_index_to_new_run {
            let header = map
                .column_headers
                .get(&u64::from(column))
                .ok_or_else(|| bad(format!("no column header for map index {column}")))?;
            new_indices
                .entry(new_run)
                .or_default()
                .insert(header.filename.clone());
        }
        let new_size = new_indices.len();
        if new_size == 1 {
            warnings.push("Number of new protein ID runs is one. Consider using mergeAllProteinRuns for some additional speed.".to_owned());
        } else if new_size >= map.column_headers.len() {
            return Err(bad(format!(
                "Number of new protein runs after merging is bigger or equal to the original ones. Aborting. Nothing would be merged. ({new_size})"
            )));
        }

        let mut old_runs = map.protein_identifications.clone();
        let mut run_id_to_new_runs: BTreeMap<String, BTreeSet<usize>> = BTreeMap::new();
        // First-seen order of the run identifiers; see the Notes above.
        let mut run_order: Vec<String> = Vec::new();
        for (&new_run, merge_request) in &new_indices {
            for old in &mut old_runs {
                let content: BTreeSet<&str> = old
                    .primary_ms_run_paths
                    .iter()
                    .map(String::as_str)
                    .collect();
                if content.iter().all(|path| merge_request.contains(*path)) {
                    if new_run as usize >= new_size {
                        return Err(bad(format!(
                            "new protein run index {new_run} is not below the number of new runs {new_size}"
                        )));
                    }
                    if run_id_to_new_runs.contains_key(&old.identifier) {
                        warnings.push("Duplicate protein run ID found. Uniquifying it.".to_owned());
                        old.identifier.push_str("_\u{1}");
                    }
                    if !run_id_to_new_runs.contains_key(&old.identifier) {
                        run_order.push(old.identifier.clone());
                    }
                    run_id_to_new_runs
                        .entry(old.identifier.clone())
                        .or_default()
                        .insert(new_run as usize);
                }
            }
        }

        let mut new_runs = vec![ProteinIdentification::default(); new_size];
        let mut collected: Vec<(ProteinHit, BTreeSet<usize>)> = Vec::new();
        let mut collected_index: BTreeMap<String, usize> = BTreeMap::new();
        let mut offsets: BTreeMap<(String, usize), usize> = BTreeMap::new();
        for run_id in &run_order {
            let targets = &run_id_to_new_runs[run_id];
            let old = old_runs
                .iter_mut()
                .find(|run| &run.identifier == run_id)
                .ok_or_else(|| bad(format!("protein run {run_id} disappeared while merging")))?;
            for &target in targets {
                let new = &mut new_runs[target];
                if new.identifier.is_empty() {
                    new.search_engine = old.search_engine.clone();
                    new.search_engine_version = old.search_engine_version.clone();
                    new.search_parameters = old.search_parameters.clone();
                    new.primary_ms_run_paths = old.primary_ms_run_paths.clone();
                    new.identifier = format!("condition{target}");
                    offsets.entry((run_id.clone(), target)).or_insert(0);
                } else {
                    // The source only logs the mismatch warnings here.
                    peptide_ids_mergeable(old, new, &map.experiment_type, &mut warnings)?;
                    let offset = new.primary_ms_run_paths.len();
                    new.primary_ms_run_paths
                        .extend(old.primary_ms_run_paths.iter().cloned());
                    offsets.entry((run_id.clone(), target)).or_insert(offset);
                }
            }
            for hit in std::mem::take(&mut old.hits) {
                let slot = match collected_index.get(&hit.accession) {
                    Some(&slot) => slot,
                    None => {
                        collected_index.insert(hit.accession.clone(), collected.len());
                        collected.push((hit, BTreeSet::new()));
                        collected.len() - 1
                    }
                };
                collected[slot].1.extend(targets.iter().copied());
            }
        }
        for (hit, targets) in &collected {
            for &target in targets {
                new_runs[target].hits.push(hit.clone());
            }
        }

        let mut next = map.clone();
        let mut untested = false;
        let mut update = |pid: &mut PeptideIdentification| -> Result<()> {
            let targets = run_id_to_new_runs.get(&pid.identifier).ok_or_else(|| {
                bad(format!(
                    "peptide identification ({}) belongs to run '{}', which is not merged",
                    describe_position(pid),
                    pid.identifier
                ))
            })?;
            if targets.len() > 1 {
                untested = true;
            }
            let old_merge_index = match pid.metadata.get(ID_MERGE_INDEX) {
                Some(value) => merge_index(value, &describe_position(pid))?
                    .ok_or_else(|| bad("negative id_merge_index"))?,
                None => 0,
            };
            let old_run = pid.identifier.clone();
            for &target in targets {
                pid.identifier = new_runs[target].identifier.clone();
                let key = if self.source_merge_index_offsets {
                    pid.identifier.clone()
                } else {
                    old_run.clone()
                };
                let offset = offsets.get(&(key, target)).copied().unwrap_or(0);
                pid.metadata.insert(
                    ID_MERGE_INDEX.into(),
                    merge_index_value(old_merge_index + offset)?,
                );
            }
            Ok(())
        };
        for feature in &mut next.features {
            for pid in &mut feature.peptide_identifications {
                update(pid)?;
            }
        }
        for pid in &mut next.unassigned_peptide_identifications {
            update(pid)?;
        }
        if untested {
            warnings.push("Warning: Merging parts of IDRuns currently untested. If it is not a TMT/iTraq sample,something is wrong anyway.".to_owned());
        }
        next.protein_identifications = new_runs;
        *map = next;
        Ok(warnings)
    }

    /// Merge every identification run of the map into one run named `merged`
    /// (source `mergeAllIDRuns`).
    ///
    /// Proteins are inserted once; peptide identifications stay unfiltered, so
    /// the same peptide may occur in several of them afterwards. The new run
    /// takes the search engine and settings of the first run and the
    /// concatenated primary MS run paths of all; each peptide's
    /// `id_merge_index` is shifted by its old run's offset in that list (a
    /// peptide of an unknown run gets offset 0, as in the source). A map with
    /// one run is unchanged. When the first run carries an inference engine, a
    /// warning says that the inference must be redone.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] (source
    /// `Exception::MissingInformation`) when the runs' search engines or
    /// settings disagree (there is no force option) or a peptide of a
    /// multi-file run has no `id_merge_index`, and [`Error::InvalidValue`]
    /// when the map has no run (the source reads out of bounds) or an
    /// `id_merge_index` is negative or non-numeric.
    pub fn merge_all_id_runs(&self, map: &mut ConsensusMap) -> Result<Vec<String>> {
        let runs = &map.protein_identifications;
        if runs.len() == 1 {
            return Ok(Vec::new());
        }
        let Some(first) = runs.first() else {
            return Err(bad("the consensus map has no protein identification run"));
        };
        let mut warnings = Vec::new();
        let mut ok = true;
        for run in runs {
            if ok {
                ok = peptide_ids_mergeable(first, run, &map.experiment_type, &mut warnings)?;
            }
        }
        if !ok {
            return Err(missing(
                "Search settings are not matching across IdentificationRuns. See warnings. Aborting..",
            ));
        }
        let mut merged = ProteinIdentification {
            identifier: "merged".into(),
            search_engine: first.search_engine.clone(),
            search_engine_version: first.search_engine_version.clone(),
            search_parameters: first.search_parameters.clone(),
            ..ProteinIdentification::default()
        };
        if !first.inference_engine()?.is_empty() {
            warnings.push("Inference was already performed on the runs in this ConsensusXML. Merging their proteins, will invalidate correctness of the inference. You should redo it.".to_owned());
        }
        let mut offsets: BTreeMap<String, (usize, bool)> = BTreeMap::new();
        for run in runs {
            let offset = merged.primary_ms_run_paths.len();
            merged
                .primary_ms_run_paths
                .extend(run.primary_ms_run_paths.iter().cloned());
            offsets
                .entry(run.identifier.clone())
                .or_insert((offset, run.primary_ms_run_paths.len() > 1));
        }
        let mut accessions = BTreeSet::new();
        for run in runs {
            for hit in &run.hits {
                if accessions.insert(hit.accession.clone()) {
                    merged.hits.push(hit.clone());
                }
            }
        }
        let mut next = map.clone();
        let update = |pid: &mut PeptideIdentification| -> Result<()> {
            let (offset, multi) = offsets.get(&pid.identifier).copied().unwrap_or((0, false));
            pid.identifier = merged.identifier.clone();
            let old = match pid.metadata.get(ID_MERGE_INDEX) {
                Some(value) => merge_index(value, &describe_position(pid))?
                    .ok_or_else(|| bad("negative id_merge_index"))?,
                None if multi => {
                    return Err(missing("No id_merge_index value in a merged ID run."));
                }
                None => 0,
            };
            pid.metadata
                .insert(ID_MERGE_INDEX.into(), merge_index_value(old + offset)?);
            Ok(())
        };
        for feature in &mut next.features {
            for pid in &mut feature.peptide_identifications {
                update(pid)?;
            }
        }
        for pid in &mut next.unassigned_peptide_identifications {
            update(pid)?;
        }
        next.protein_identifications = vec![merged];
        *map = next;
        Ok(warnings)
    }
}

fn merge_index_value(index: usize) -> Result<MetaValue> {
    i64::try_from(index)
        .map(MetaValue::from)
        .map_err(|_| bad("id_merge_index exceeds the integer meta value range"))
}
