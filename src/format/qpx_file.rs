// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Peptide-spectrum matches as Arrow tables and QPX Parquet files.
//!
//! Ports `QPXFile` (`FORMAT/QPXFile.h` and `QPXFile.cpp`): the export of
//! `PeptideIdentification`s to two table layouts and the import of one of
//! them. See `docs/QPX_SUPPORT.md` for the API mapping.
//!
//! * [`QpxFile::export_to_arrow`](crate::format::qpx_file::QpxFile::export_to_arrow) writes the internal
//!   [`PSMSchema`](crate::format::arrow_schema_registry::PSMSchema) layout that the
//!   FeatureMap and ConsensusMap Parquet bundles use for lossless round
//!   trips, and [`QpxFile::import_from_arrow`](crate::format::qpx_file::QpxFile::import_from_arrow) reads it back.
//! * [`QpxFile::export_psms_to_qpx_arrow`](crate::format::qpx_file::QpxFile::export_psms_to_qpx_arrow) writes the QPX (Quantitative
//!   Proteomics eXchange, quantms) psm view,
//!   [`QPXPSMSchema`](crate::format::arrow_schema_registry::QPXPSMSchema), for cross-tool
//!   exchange; [`QpxFile::export_to_parquet`](crate::format::qpx_file::QpxFile::export_to_parquet),
//!   [`QpxFile::export_table_to_parquet`](crate::format::qpx_file::QpxFile::export_table_to_parquet) and
//!   [`QpxFile::export_to_parquet_streaming`](crate::format::qpx_file::QpxFile::export_to_parquet_streaming) write it as a QPX `psm_file`,
//!   whose values [`QpxValueValidation`](crate::format::qpx_value_validation::QpxValueValidation) checks before anything is written.
//!
//! The source marks the class experimental ("may change in future
//! versions").
//!
//! **Origin files.** A merged identification run (ProteomicsLFQ, IDMerger)
//! lists all its input files as primary MS run paths and points each
//! identification at one of them through its `id_merge_index` metadata. Both
//! exports resolve each identification's file through
//! [`IdentifierMSRunMapper`](crate::identification::IdentifierMSRunMapper) and refuse merged input without a usable index
//! (see [`QpxFile::require_resolvable_merge_indices`](crate::format::qpx_file::QpxFile::require_resolvable_merge_indices)) rather than labelling
//! every identification with the run's first file. The internal layout keeps
//! the full path (`reference_file_name`); the QPX view keeps only the file's
//! stem (`run_file_name`, a primary-key component), preferring a per-hit or
//! per-identification file reference
//! ([`qpx_psm_run_file_name`](crate::format::arrow_io_helpers::qpx_psm_run_file_name)).
//!
//! **Model mapping.** The source's `PeptideHit` keeps its rank as the `rank`
//! metadata value, written only when it is not 0; this crate's
//! [`PeptideHit`](crate::identification::PeptideHit) has a [`rank`](crate::identification::PeptideHit::rank) field, which the export
//! writes as that metadata entry (when not 0) and the import restores. A
//! spectrum reference is the identification's `spectrum_reference` metadata
//! in both.
//!
//! **Metadata order.** The source lists metadata in `MetaInfoRegistry` index
//! order, that is in the order the process first registered each name; this
//! crate's metadata has no registry and lists them in key order. This
//! decides the order of the `additional_scores`, `psm_metavalues` and
//! `spectrum_metavalues` entries of a row, not their content.
//!
//! **Warnings.** Where the source logs a warning and continues, the port
//! returns the warning lines ([`QpxOutput::warnings`](crate::format::qpx_file::QpxOutput::warnings) or the `Vec<String>`
//! of the import and file writers). The source's `false`/`nullptr` returns
//! with a logged error are [`Result`](crate::Result) errors here.
//!
//! **Evidence.** `tests/qpx_file.rs` compares every export cell by cell with
//! the tables the executed Linux x86_64 Release build produced for the same
//! identifications (oracle `F1c-qpx`, `tests/data/qpx/`), and carries every
//! section of the upstream class test `QPXFile_test.cpp`.

use crate::analysis::scores::{ScoreType, find_peptide_score};
use crate::chemistry::proforma::{ConversionPolicy, Peptidoform, WriteMode};
use crate::chemistry::{AASequence, ModificationsDB, PeptideFragmentType, SequenceModification};
use crate::concept::constants::user_param::{FRAGMENT_ANNOTATION_USERPARAM, RANK};
use crate::concept::fuzzy_string_comparator::format_g;
use crate::concept::parallel::{Threads, map_collect};
use crate::format::arrow_io_helpers::{
    self as helpers, MetaEntry, ParquetWriteConfig, QpxRunFileNameKeys, configured_properties,
    data_value_text, get_bool_value, get_double_value, get_int32_value, get_string_value, is_null,
    list_array, list_element, meta_value_type_token, metavalues_array, qpx_file_metadata,
    qpx_psm_run_file_name, qpx_run_file_name, qpx_scan_components, qpx_scan_format_of_all,
    read_meta_values, struct_array, write_table_to_parquet,
};
use crate::format::arrow_schema_registry::validation::{Mode, validate};
use crate::format::arrow_schema_registry::{PSMSchema, QPXPSMSchema};
use crate::format::parquet_file::{write_row_group, writer_options};
use crate::format::qpx_identity::{self, FeatureLinks};
use crate::format::qpx_value_validation::{QpxValueValidation, View};
use crate::identification::{
    FlankingResidue, IdentifierMSRunMapper, PeptideEvidence, PeptideHit, PeptideIdentification,
    ProteinIdentification,
};
use crate::metadata::{MetaInfo, MetaValue, MetaValueData};
use crate::{Error, Result};
use arrow_array::cast::AsArray;
use arrow_array::types::{Float64Type, Int32Type};
use arrow_array::{
    Array, ArrayRef, BooleanArray, Float32Array, Float64Array, Int16Array, Int32Array, Int64Array,
    RecordBatch, StringArray, new_null_array,
};
use arrow_schema::{DataType, Schema, SchemaRef};
use parquet::arrow::ArrowWriter;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::path::Path;
use std::sync::Arc;

/// Export and import of peptide-spectrum matches as Arrow tables and QPX
/// Parquet files.
///
/// Ports `OpenMS::QPXFile`, a class of static members; each is an
/// associated function here.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct QpxFile;

/// A result together with the warning lines the source logs while producing
/// it (`OPENMS_LOG_WARN`), in the order it logs them.
#[derive(Clone, Debug, PartialEq)]
pub struct QpxOutput<T> {
    /// The result.
    pub value: T,
    /// One line per warning.
    pub warnings: Vec<String>,
}

/// The streaming writer's knobs, the source's defaulted parameters of
/// `exportToParquetStreaming`.
#[derive(Clone, Debug, PartialEq)]
pub struct StreamingOptions {
    /// Export every hit of an identification rather than only its first
    /// (default `false`).
    pub export_all_psms: bool,
    /// Number of identifications materialised into one in-memory table at a
    /// time (default 1,000,000). This is the peak-memory knob only: with
    /// `export_all_psms` one identification can emit several rows, and the
    /// row-group size is `config.row_group_size`. 0 means the default.
    pub batch_size: usize,
    /// Parquet writing options; `row_group_size` is the maximum number of
    /// rows per Parquet row group.
    pub config: ParquetWriteConfig,
    /// Threads that build each batch in parallel: 1 = serial (default), 0 or
    /// less = every available core, N = that many. The file is the same
    /// whatever the count.
    pub n_threads: i32,
}

impl Default for StreamingOptions {
    fn default() -> Self {
        Self {
            export_all_psms: false,
            batch_size: DEFAULT_BATCH_SIZE,
            config: ParquetWriteConfig::default(),
            n_threads: 1,
        }
    }
}

/// The source's default `batch_size`.
const DEFAULT_BATCH_SIZE: usize = 1_000_000;

/// Largest row count Parquet C++ puts in one row group
/// (`max_row_group_length`, 1Mi rows); a larger `row_group_size` is capped.
const MAX_ROW_GROUP_ROWS: usize = 1 << 20;

/// Metadata keys of a hit with dedicated columns: neither an additional score
/// nor a `psm_metavalues` entry.
const EXCLUDED_HIT_KEYS: [&str; 8] = [
    "target_decoy",
    "predicted_RT",
    "predicted_rt",
    "ion_mobility",
    "IM",
    "scan",
    "reference_file_name",
    FRAGMENT_ANNOTATION_USERPARAM,
];

/// Metadata keys of an identification with dedicated columns.
const EXCLUDED_SPECTRUM_KEYS: [&str; 3] = ["spectrum_reference", "ion_mobility", "IM"];

/// `psm_metavalues` keys the import does not restore (their columns do).
const IMPORT_EXCLUDED_HIT_KEYS: [&str; 7] = [
    "target_decoy",
    "predicted_RT",
    "predicted_rt",
    "ion_mobility",
    "IM",
    "scan",
    "reference_file_name",
];

/// The warning both QPX exports log when identifications have no
/// resolvable origin file.
fn unattributable_warning(count: usize) -> String {
    format!(
        "QPXFile: {count} PSM(s) have no resolvable origin file, so their run_file_name is \
         empty -- a QPX primary-key component that must not be, and the key the psm view is \
         joined to feature and pg on. Protein inference commonly drops the path: set it with \
         ProteinIdentification::setPrimaryMSRunPath() before exporting."
    )
}

// ---------------------------------------------------------------------------
// Run mapping
// ---------------------------------------------------------------------------

/// Refuse identification runs whose identifiers cannot uniquely select an
/// MS-run path list.
fn require_unique_run_identifiers(proteins: &[ProteinIdentification]) -> Result<()> {
    let mut seen: HashSet<&str> = HashSet::with_capacity(proteins.len());
    for run in proteins {
        if !seen.insert(run.identifier.as_str()) {
            // The source's Exception::InvalidValue message, with the
            // identifier as its value.
            return Err(Error::InvalidValue(format!(
                "the value '{id}' was used but is not valid; QPXFile: several identification \
                 runs share the identifier '{id}'. Their MS run paths cannot be told apart, so a \
                 PSM could be attributed to the wrong run file.",
                id = run.identifier
            )));
        }
    }
    Ok(())
}

/// The run-identifier -> primary MS run paths mapping each PSM's origin file
/// is resolved with.
///
/// Duplicate path lists across runs make the mapper's reverse direction fail,
/// but its forward map (the only direction used here) is complete by then,
/// so that failure is a warning, as in the source. Duplicate identifiers are
/// refused first, because they overwrite the forward map itself.
fn build_run_mapper(
    proteins: &[ProteinIdentification],
    warnings: &mut Vec<String>,
) -> Result<IdentifierMSRunMapper> {
    require_unique_run_identifiers(proteins)?;
    let mut mapper = IdentifierMSRunMapper::new();
    if let Err(error) = mapper.create(proteins) {
        // A resource refusal happens before anything is mapped; only the
        // ambiguity failure leaves a usable forward map.
        if mapper.len() != proteins.len() {
            return Err(error);
        }
        let message = match &error {
            Error::InvalidValue(message) => message.clone(),
            other => other.to_string(),
        };
        warnings.push(format!(
            "QPXFile: could not build a complete MS-run mapping: {message}. Per-PSM \
             run_file_name resolution may be incomplete."
        ));
    }
    Ok(mapper)
}

/// One warning per QPX `run_file_name` that several distinct origin files
/// share after stripping path and extension, in name order.
fn stem_collision_warnings(mapper: &IdentifierMSRunMapper) -> Vec<String> {
    let mut stems: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
    for identifier in mapper.identifiers() {
        for path in mapper.ms_run_paths(identifier) {
            if path.is_empty() {
                continue;
            }
            stems
                .entry(qpx_run_file_name(path))
                .or_default()
                .insert(path.as_str());
        }
    }
    stems
        .into_iter()
        .filter(|(_, paths)| paths.len() >= 2)
        .map(|(stem, paths)| {
            format!(
                "QPXFile: several MS runs share the run_file_name '{stem}' after stripping path \
                 and extension: {}. They cannot be told apart in the exported QPX tables.",
                paths.into_iter().collect::<Vec<_>>().join(", ")
            )
        })
        .collect()
}

fn validate_merge_indices<'a>(
    mapper: &IdentifierMSRunMapper,
    peptides: impl Iterator<Item = &'a PeptideIdentification>,
) -> Result<()> {
    for (index, peptide) in peptides.enumerate() {
        mapper.validate_merge_index(peptide, index)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Per-hit values shared by both layouts
// ---------------------------------------------------------------------------

/// One modification position: (1-based position, one-letter residue);
/// position 0 is the N-terminus and `len + 1` the C-terminus, which carry no
/// residue.
type Position = (i32, Option<char>);

/// One modification entry: its name, UniMod accession and positions.
struct ModificationEntry {
    name: String,
    accession: Option<String>,
    positions: Vec<Position>,
}

/// The name a modification is grouped and written under: the source's
/// `getModificationName()`, which is the modification's id. A mass tag the
/// sequence parser could not resolve (`M[+123.456]`) is a modification the
/// source creates with only a full id, so its name is empty (measured).
fn modification_name(modification: &SequenceModification) -> String {
    match modification {
        SequenceModification::Known(known) => known.name().to_owned(),
        SequenceModification::MassTag(_) => String::new(),
    }
}

fn unimod_accession(modification: &SequenceModification) -> Option<String> {
    modification
        .record_id()
        .filter(|id| *id > 0)
        .map(|id| format!("UNIMOD:{id}"))
}

fn position(value: usize) -> Result<i32> {
    i32::try_from(value)
        .map_err(|_| Error::InvalidValue("QPXFile: a sequence position exceeds int32".into()))
}

/// The modifications of `sequence`, grouped by name in name order, each with
/// the UniMod accession of its first occurrence (`UNIMOD:<id>`, `None` when
/// the modification has no UniMod record) and every position it occurs at:
/// the N-terminus, residues in sequence order, then the C-terminus.
fn modification_entries(sequence: &AASequence) -> Result<Vec<ModificationEntry>> {
    let mut entries: BTreeMap<String, (Option<String>, Vec<Position>)> = BTreeMap::new();
    if let Some(modification) = sequence.n_terminal_modification() {
        // The source assigns rather than inserts here; the N-terminus is
        // always the first entry, so the two agree.
        entries.insert(
            modification_name(modification),
            (unimod_accession(modification), vec![(0, None)]),
        );
    }
    for (index, residue) in sequence.as_str().chars().enumerate() {
        if let Some(modification) = sequence.residue_modification(index)? {
            entries
                .entry(modification_name(modification))
                .or_insert_with(|| (unimod_accession(modification), Vec::new()))
                .1
                .push((position(index + 1)?, Some(residue)));
        }
    }
    if let Some(modification) = sequence.c_terminal_modification() {
        let at = position(sequence.len() + 1)?;
        entries
            .entry(modification_name(modification))
            .or_insert_with(|| (unimod_accession(modification), Vec::new()))
            .1
            .push((at, None));
    }
    Ok(entries
        .into_iter()
        .map(|(name, (accession, positions))| ModificationEntry {
            name,
            accession,
            positions,
        })
        .collect())
}

/// The m/z of a sequence at a non-zero charge, as `AASequence::getMZ(charge)`
/// computes it: the full ion's monoisotopic mass summed in the source's order
/// (charge protons, terminal modifications, residues, then water), divided by
/// the charge; 0 for an empty sequence, whose mass the source reports as 0.
///
/// The last bits can still differ from the source's: the sequence port sums
/// a modified residue as its internal formula's mass plus the modification's,
/// the source as the full residue's formula mass minus water (measured: 2 ULP
/// on `PEM(Oxidation)TIDER` at charge 2, see `docs/QPX_SUPPORT.md`).
fn sequence_mz(sequence: &AASequence, charge: i32) -> Result<f64> {
    Ok(sequence.mono_mass_for(PeptideFragmentType::Full, charge)? / f64::from(charge))
}

/// The ProForma canonical text of a sequence (`ProForma::toString(
/// fromAASequence(seq), CANONICAL)`).
fn canonical_peptidoform(sequence: &AASequence) -> Result<String> {
    Peptidoform::from_aa_sequence(sequence)?.to_text(WriteMode::Canonical)
}

/// A metadata value as a number, as the source's `DataValue` to `double`
/// conversion: an integer or a float.
///
/// The source converts a string or list value too, reading the bits of its
/// storage as a double (`DataValue::operator double` checks only for an
/// empty value); this refuses it.
fn numeric(value: &MetaValue, key: &str) -> Result<f64> {
    match value.data() {
        MetaValueData::Integer(number) => Ok(*number as f64),
        MetaValueData::Float(number) => Ok(*number),
        // The source's Exception::ConversionError message.
        MetaValueData::Empty => Err(Error::InvalidValue(
            "Could not convert DataValue of type 'Empty' to double".into(),
        )),
        _ => Err(Error::InvalidValue(format!(
            "QPXFile: metadata value '{key}' of type '{}' and value '{}' is not numeric",
            data_type_name(value),
            data_value_text(value)
        ))),
    }
}

/// The source's `DataValue::NamesOfDataType` entry of a value's type.
fn data_type_name(value: &MetaValue) -> &'static str {
    match value.data() {
        MetaValueData::String(_) => "String",
        MetaValueData::Integer(_) => "Int",
        MetaValueData::Float(_) => "Double",
        MetaValueData::StringList(_) => "StringList",
        MetaValueData::IntegerList(_) => "IntList",
        MetaValueData::FloatList(_) => "DoubleList",
        MetaValueData::Empty => "Empty",
    }
}

/// The first of `keys` present in `metadata`, as a number.
fn first_numeric(metadata: &MetaInfo, keys: &[&str]) -> Result<Option<f64>> {
    for key in keys {
        if let Some(value) = metadata.get(*key) {
            return numeric(value, key).map(Some);
        }
    }
    Ok(None)
}

/// The ion mobility of a hit: its own `ion_mobility` or `IM`, else its
/// identification's.
fn ion_mobility(hit: &PeptideHit, peptide: &PeptideIdentification) -> Result<Option<f64>> {
    match first_numeric(&hit.metadata, &["ion_mobility", "IM"])? {
        Some(value) => Ok(Some(value)),
        None => first_numeric(&peptide.metadata, &["ion_mobility", "IM"]),
    }
}

/// Whether the hit's `target_decoy` annotation starts with `decoy`, or
/// `None` when it has none.
fn decoy_annotation(hit: &PeptideHit) -> Option<bool> {
    hit.metadata
        .get("target_decoy")
        .map(|value| data_value_text(value).starts_with("decoy"))
}

/// The posterior error probability of one hit: the main score when the
/// identification's score type is a PEP, else the first PEP-named metadata
/// value of its first hit, looked up on this hit.
fn posterior_error_probability(
    hit: &PeptideHit,
    pep_score: Option<&crate::analysis::scores::ScoreSearchResult>,
) -> Result<Option<f64>> {
    match pep_score {
        Some(found) if found.is_main_score => Ok(Some(hit.score)),
        Some(found) => match hit.metadata.get(&found.name) {
            Some(value) => numeric(value, &found.name).map(Some),
            None => Ok(None),
        },
        None => Ok(None),
    }
}

/// One `additional_scores` entry.
struct AdditionalScore {
    name: String,
    value: f64,
    higher_better: Option<bool>,
}

/// The hit's metadata that the exports list: in key order, without the keys
/// that have dedicated columns, the rank field as the source's `rank` value
/// (when not 0).
fn hit_metadata(hit: &PeptideHit) -> Vec<(String, MetaValue)> {
    let mut entries: BTreeMap<String, MetaValue> = hit
        .metadata
        .iter()
        .filter(|(key, _)| !EXCLUDED_HIT_KEYS.contains(&key.as_str()))
        .map(|(key, value)| (key.clone(), value.clone()))
        .collect();
    if hit.rank != 0 {
        entries.insert(RANK.to_owned(), MetaValue::from(i64::from(hit.rank)));
    }
    entries.into_iter().collect()
}

/// Split the hit's listed metadata into known numeric scores and the rest.
///
/// A value is an additional score when it is an integer or a float and its
/// name, without a `_score` suffix, is one of the known score names
/// (`Scores::isKnownScoreType`); `higher_better` is the orientation of the
/// score category the name parses as (`IDScoreSwitcherAlgorithm::
/// toScoreTypeEnum`), null when it names none (e.g. `expect`).
fn split_scores(hit: &PeptideHit) -> (Vec<AdditionalScore>, Vec<MetaEntry>) {
    let mut scores = Vec::new();
    let mut rest = Vec::new();
    for (key, value) in hit_metadata(hit) {
        let number = match value.data() {
            MetaValueData::Integer(number) => Some(*number as f64),
            MetaValueData::Float(number) => Some(*number),
            _ => None,
        };
        match number {
            Some(number) if ScoreType::from_normalized_name(&key).is_some() => {
                scores.push(AdditionalScore {
                    higher_better: ScoreType::parse(&key).ok().map(ScoreType::higher_is_better),
                    name: key,
                    value: number,
                });
            }
            _ => rest.push([
                key,
                data_value_text(&value),
                meta_value_type_token(&value).to_owned(),
            ]),
        }
    }
    (scores, rest)
}

fn spectrum_metavalues(peptide: &PeptideIdentification) -> Vec<MetaEntry> {
    peptide
        .metadata
        .iter()
        .filter(|(key, _)| !EXCLUDED_SPECTRUM_KEYS.contains(&key.as_str()))
        .map(|(key, value)| {
            [
                key.clone(),
                data_value_text(value),
                meta_value_type_token(value).to_owned(),
            ]
        })
        .collect()
}

/// The rows an identification contributes: every hit or the first.
fn exported_hits(peptide: &PeptideIdentification, export_all_psms: bool) -> &[PeptideHit] {
    if export_all_psms {
        &peptide.hits
    } else {
        &peptide.hits[..peptide.hits.len().min(1)]
    }
}

fn int32(value: usize, what: &str) -> Result<i32> {
    i32::try_from(value)
        .map_err(|_| Error::InvalidValue(format!("QPXFile: {what} {value} exceeds int32")))
}

// ---------------------------------------------------------------------------
// Column assembly
// ---------------------------------------------------------------------------

fn utf8(values: Vec<Option<String>>) -> ArrayRef {
    Arc::new(StringArray::from(values))
}

fn field_type(schema: &Schema, name: &str) -> Result<DataType> {
    schema
        .field_with_name(name)
        .map(|field| field.data_type().clone())
        .map_err(|error| Error::InvalidValue(format!("QPXFile: {error}")))
}

/// A `list<primitive>` column: one row per entry (`None` null).
fn primitive_lists<T>(
    data_type: &DataType,
    rows: &[Option<Vec<T>>],
    values: impl FnOnce(Vec<T>) -> ArrayRef,
) -> Result<ArrayRef>
where
    T: Clone,
{
    let lengths: Vec<Option<usize>> = rows.iter().map(|row| row.as_ref().map(Vec::len)).collect();
    let flat: Vec<T> = rows.iter().flatten().flatten().cloned().collect();
    list_array(data_type, &lengths, values(flat))
}

fn additional_scores_array(
    data_type: &DataType,
    rows: &[Vec<AdditionalScore>],
) -> Result<ArrayRef> {
    let lengths: Vec<Option<usize>> = rows.iter().map(|row| Some(row.len())).collect();
    let flat: Vec<&AdditionalScore> = rows.iter().flatten().collect();
    let values = struct_array(
        list_element(data_type)?,
        vec![
            utf8(flat.iter().map(|s| Some(s.name.clone())).collect()),
            Arc::new(Float64Array::from(
                flat.iter().map(|s| s.value).collect::<Vec<_>>(),
            )),
            Arc::new(BooleanArray::from(
                flat.iter().map(|s| s.higher_better).collect::<Vec<_>>(),
            )),
        ],
    )?;
    list_array(data_type, &lengths, values)
}

/// A `modifications` column of either layout. `qpx` selects the QPX
/// position struct (`int32` position, amino acid, an empty `scores` list)
/// over the internal one (`"<residue>.<position>"`, `N-term.0`,
/// `C-term.<len+1>`; a null score).
fn modifications_array(
    data_type: &DataType,
    rows: &[Vec<ModificationEntry>],
    qpx: bool,
) -> Result<ArrayRef> {
    let entry_type = list_element(data_type)?;
    let DataType::Struct(entry_fields) = entry_type else {
        return Err(Error::InvalidValue(
            "QPXFile: modifications entries are not structs".into(),
        ));
    };
    let positions_type = entry_fields[2].data_type();
    let position_type = list_element(positions_type)?;
    let entries: Vec<&ModificationEntry> = rows.iter().flatten().collect();
    let positions: Vec<Position> = entries
        .iter()
        .flat_map(|entry| entry.positions.iter().copied())
        .collect();
    let position_children: Vec<ArrayRef> = if qpx {
        let DataType::Struct(position_fields) = position_type else {
            return Err(Error::InvalidValue(
                "QPXFile: modification positions are not structs".into(),
            ));
        };
        let scores_type = position_fields[2].data_type();
        let score_type = list_element(scores_type)?;
        let DataType::Struct(score_fields) = score_type else {
            return Err(Error::InvalidValue(
                "QPXFile: modification scores are not structs".into(),
            ));
        };
        let empty_scores = struct_array(
            score_type,
            score_fields
                .iter()
                .map(|field| arrow_array::new_empty_array(field.data_type()))
                .collect(),
        )?;
        vec![
            Arc::new(Int32Array::from(
                positions.iter().map(|(at, _)| *at).collect::<Vec<_>>(),
            )),
            utf8(
                positions
                    .iter()
                    .map(|(_, residue)| residue.map(String::from))
                    .collect(),
            ),
            list_array(scores_type, &vec![Some(0); positions.len()], empty_scores)?,
        ]
    } else {
        let text = |(at, residue): (i32, Option<char>)| match residue {
            Some(residue) => format!("{residue}.{at}"),
            None if at == 0 => "N-term.0".to_owned(),
            None => format!("C-term.{at}"),
        };
        vec![
            utf8(positions.iter().map(|at| Some(text(*at))).collect()),
            Arc::new(Float64Array::from(vec![None::<f64>; positions.len()])),
        ]
    };
    let position_values = struct_array(position_type, position_children)?;
    let positions_list = list_array(
        positions_type,
        &entries
            .iter()
            .map(|entry| Some(entry.positions.len()))
            .collect::<Vec<_>>(),
        position_values,
    )?;
    let entry_values = struct_array(
        entry_type,
        vec![
            utf8(
                entries
                    .iter()
                    .map(|entry| Some(entry.name.clone()))
                    .collect(),
            ),
            utf8(
                entries
                    .iter()
                    .map(|entry| entry.accession.clone())
                    .collect(),
            ),
            positions_list,
        ],
    )?;
    list_array(
        data_type,
        &rows.iter().map(|row| Some(row.len())).collect::<Vec<_>>(),
        entry_values,
    )
}

/// The four fragment-annotation list columns: one list per hit with
/// annotations, null for a hit without.
fn fragment_arrays(
    schema: &Schema,
    names: [&str; 4],
    hits: &[&PeptideHit],
) -> Result<[ArrayRef; 4]> {
    let rows: Vec<Option<&[crate::identification::PeakAnnotation]>> = hits
        .iter()
        .map(|hit| (!hit.peak_annotations.is_empty()).then_some(hit.peak_annotations.as_slice()))
        .collect();
    let lengths: Vec<Option<usize>> = rows.iter().map(|row| row.map(<[_]>::len)).collect();
    let flat: Vec<&crate::identification::PeakAnnotation> =
        rows.iter().flatten().flat_map(|row| row.iter()).collect();
    Ok([
        list_array(
            &field_type(schema, names[0])?,
            &lengths,
            Arc::new(Float32Array::from(
                flat.iter().map(|a| a.mz as f32).collect::<Vec<_>>(),
            )),
        )?,
        list_array(
            &field_type(schema, names[1])?,
            &lengths,
            Arc::new(Float32Array::from(
                flat.iter().map(|a| a.intensity as f32).collect::<Vec<_>>(),
            )),
        )?,
        list_array(
            &field_type(schema, names[2])?,
            &lengths,
            Arc::new(Int32Array::from(
                flat.iter().map(|a| a.charge).collect::<Vec<_>>(),
            )),
        )?,
        list_array(
            &field_type(schema, names[3])?,
            &lengths,
            utf8(flat.iter().map(|a| Some(a.annotation.clone())).collect()),
        )?,
    ])
}

/// Assemble a table and check it strictly against its registry schema, as
/// the source does before returning it.
fn finish_table(schema: SchemaRef, columns: Vec<ArrayRef>) -> Result<RecordBatch> {
    let table = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|error| Error::InvalidValue(format!("QPXFile: {error}")))?;
    let validation = validate(&table, &schema, Mode::Strict);
    if !validation.valid {
        return Err(Error::InvalidValue(format!(
            "QPXFile: Schema validation failed: {validation}"
        )));
    }
    Ok(table)
}

// ---------------------------------------------------------------------------
// The QPX psm view
// ---------------------------------------------------------------------------

/// The QPX psm table of `peptides`, and how many of its rows have no
/// resolvable origin file.
///
/// Every row is independent of the others, so any contiguous range of
/// identifications gives a self-contained table; the streaming writer relies
/// on that.
fn build_qpx_psm_table(
    mapper: &IdentifierMSRunMapper,
    peptides: &[&PeptideIdentification],
    export_all_psms: bool,
    feature_links: Option<&FeatureLinks>,
) -> Result<(RecordBatch, usize)> {
    let schema = QPXPSMSchema::schema();
    let keys = QpxRunFileNameKeys::new();
    let mut unattributable = 0usize;

    let mut psm_ids = Vec::new();
    let mut feature_ids: Vec<Option<i64>> = Vec::new();
    let mut sequences = Vec::new();
    let mut peptidoforms = Vec::new();
    let mut modifications = Vec::new();
    let mut charges = Vec::new();
    let mut peps = Vec::new();
    let mut decoys = Vec::new();
    let mut calculated = Vec::new();
    let mut observed = Vec::new();
    let mut mass_errors = Vec::new();
    let mut scores = Vec::new();
    let mut predicted = Vec::new();
    let mut run_names = Vec::new();
    let mut scans: Vec<Option<Vec<i32>>> = Vec::new();
    let mut rts = Vec::new();
    let mut mobilities = Vec::new();
    let mut missed = Vec::new();
    let mut accessions: Vec<Option<Vec<String>>> = Vec::new();
    let mut hits_out: Vec<&PeptideHit> = Vec::new();

    for peptide in peptides {
        if peptide.hits.is_empty() {
            continue;
        }
        let pep_score = find_peptide_score(peptide, ScoreType::PosteriorErrorProbability);
        let resolved = mapper.primary_ms_run_path(peptide)?;
        let spectrum_reference = peptide.spectrum_reference();
        let scan = qpx_scan_components(&spectrum_reference);
        for hit in exported_hits(peptide, export_all_psms) {
            let sequence = &hit.sequence;
            sequences.push(sequence.as_str().to_owned());
            let peptidoform = canonical_peptidoform(sequence)?;
            modifications.push(modification_entries(sequence)?);

            // Narrowed once: psm_id hashes the charge as persisted. The
            // source's static_cast<int16_t> wraps, and so does this.
            let charge = hit.charge as i16;
            charges.push(charge);
            peps.push(posterior_error_probability(hit, pep_score.as_ref())?);
            decoys.push(decoy_annotation(hit).unwrap_or(false));

            let calc_mz = if charge > 0 {
                sequence_mz(sequence, i32::from(charge))? as f32
            } else {
                0.0
            };
            calculated.push(calc_mz);
            let obs_mz = peptide
                .mz
                .map(|mz| mz as f32)
                .filter(|mz| !mz.is_nan())
                .unwrap_or(0.0);
            observed.push(obs_mz);
            mass_errors.push(if charge > 0 && obs_mz > 0.0 && calc_mz > 0.0 {
                // float arithmetic, then the source's float * 1e6 (double).
                let relative = (obs_mz - calc_mz) / calc_mz;
                Some((f64::from(relative) * 1e6) as f32)
            } else {
                None
            });

            let (hit_scores, _) = split_scores(hit);
            scores.push(hit_scores);
            predicted.push(
                first_numeric(&hit.metadata, &["predicted_RT", "predicted_rt"])?
                    .map(|value| value as f32),
            );

            let run_stem =
                qpx_psm_run_file_name(&hit.metadata, &peptide.metadata, &resolved, &keys);
            if run_stem.is_empty() {
                unattributable += 1;
            }
            let psm_id = qpx_identity::psm_id(&run_stem, &scan, &peptidoform, i64::from(charge));
            psm_ids.push(psm_id);
            feature_ids.push(feature_links.and_then(|links| links.get(&psm_id).copied()));
            run_names.push(run_stem);
            peptidoforms.push(peptidoform);
            scans.push(Some(scan.clone()));

            rts.push(peptide.rt.filter(|rt| !rt.is_nan()).map(|rt| rt as f32));
            mobilities.push(ion_mobility(hit, peptide)?.map(|value| value as f32));
            missed.push(match hit.metadata.get("missed_cleavages") {
                // The source's Exception::ConversionError message
                // (`DataValue::operator int`).
                Some(value) => Some(value.as_i64().map_err(|_| {
                    Error::InvalidValue(format!(
                        "Could not convert non-integer DataValue of type '{}' and value '{}' to \
                         int",
                        data_type_name(value),
                        data_value_text(value)
                    ))
                })? as i16),
                None => None,
            });
            accessions.push(Some(
                hit.evidences
                    .iter()
                    .map(|evidence| evidence.protein_accession.clone())
                    .collect(),
            ));
            hits_out.push(hit);
        }
    }

    let rows = sequences.len();
    let fragments = fragment_arrays(
        &schema,
        [
            QPXPSMSchema::MZ_ARRAY,
            QPXPSMSchema::INTENSITY_ARRAY,
            QPXPSMSchema::CHARGE_ARRAY,
            QPXPSMSchema::ION_TYPE_ARRAY,
        ],
        &hits_out,
    )?;
    let [mz_array, intensity_array, charge_array, ion_type_array] = fragments;
    let columns: Vec<ArrayRef> = vec![
        Arc::new(Int64Array::from(psm_ids)),
        utf8(sequences.into_iter().map(Some).collect()),
        utf8(peptidoforms.into_iter().map(Some).collect()),
        modifications_array(
            &field_type(&schema, QPXPSMSchema::MODIFICATIONS)?,
            &modifications,
            true,
        )?,
        Arc::new(Int16Array::from(charges)),
        Arc::new(Float64Array::from(peps)),
        Arc::new(BooleanArray::from(decoys)),
        Arc::new(Float32Array::from(calculated)),
        Arc::new(Float32Array::from(observed)),
        Arc::new(Float32Array::from(mass_errors)),
        additional_scores_array(
            &field_type(&schema, QPXPSMSchema::ADDITIONAL_SCORES)?,
            &scores,
        )?,
        Arc::new(Float32Array::from(predicted)),
        utf8(run_names.into_iter().map(Some).collect()),
        new_null_array(&field_type(&schema, QPXPSMSchema::CV_PARAMS)?, rows),
        primitive_lists(&field_type(&schema, QPXPSMSchema::SCAN)?, &scans, |flat| {
            Arc::new(Int32Array::from(flat))
        })?,
        Arc::new(Float32Array::from(rts)),
        Arc::new(Float32Array::from(mobilities)),
        Arc::new(Int16Array::from(missed)),
        primitive_lists(
            &field_type(&schema, QPXPSMSchema::PROTEIN_ACCESSIONS)?,
            &accessions,
            |flat| Arc::new(StringArray::from(flat)),
        )?,
        new_null_array(&field_type(&schema, QPXPSMSchema::CROSS_LINKS)?, rows),
        mz_array,
        intensity_array,
        charge_array,
        ion_type_array,
        new_null_array(
            &field_type(&schema, QPXPSMSchema::ION_MOBILITY_ARRAY)?,
            rows,
        ),
        Arc::new(Int64Array::from(feature_ids)),
    ];
    Ok((finish_table(schema, columns)?, unattributable))
}

/// The QPX `psm_file` metadata; `scan_format` is omitted when empty.
fn qpx_psm_metadata(
    config: &ParquetWriteConfig,
    scan_format: &str,
) -> Result<HashMap<String, String>> {
    let mut extra = BTreeMap::new();
    if !scan_format.is_empty() {
        extra.insert("scan_format".to_owned(), scan_format.to_owned());
    }
    Ok(qpx_file_metadata("psm_file", config, &extra)?.to_schema_metadata())
}

fn with_metadata(table: &RecordBatch, metadata: HashMap<String, String>) -> Result<RecordBatch> {
    let schema = Arc::new(Schema::new(table.schema().fields().clone()).with_metadata(metadata));
    RecordBatch::try_new(schema, table.columns().to_vec())
        .map_err(|error| Error::InvalidValue(format!("QPXFile: {error}")))
}

/// The QPX `scan_format` of the identifications' spectrum references, and
/// the warning when they mix conventions.
fn scan_format<'a>(
    peptides: impl Iterator<Item = &'a PeptideIdentification>,
    warnings: &mut Vec<String>,
) -> &'static str {
    let references: Vec<String> = peptides
        .map(PeptideIdentification::spectrum_reference)
        .filter(|reference| !reference.is_empty())
        .collect();
    let (format, warning) = qpx_scan_format_of_all(&references);
    warnings.extend(warning);
    format
}

impl QpxFile {
    /// Refuse identifications of a merged run that carry no usable
    /// `id_merge_index`.
    ///
    /// Ports `requireResolvableMergeIndices(protein_identifications,
    /// peptide_identifications)`. Without the index every identification of
    /// a merged run resolves to the run's first file, so `run_file_name` --
    /// a QPX primary-key component -- would be wrong rather than missing.
    /// Call it before an output file is opened: the source's streaming build
    /// could only report the failure after truncating its output. Runs with
    /// zero or one path, and identifications of unknown runs, pass.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when identification runs share an identifier;
    /// [`Error::MissingInformation`] when an identification of a merged run
    /// has no `id_merge_index`, a non-integer one or one out of range, with
    /// the message
    /// [`IdentifierMSRunMapper::validate_merge_index`] documents, `#<n>`
    /// counting from 0 in `peptides`.
    pub fn require_resolvable_merge_indices(
        proteins: &[ProteinIdentification],
        peptides: &[PeptideIdentification],
    ) -> Result<()> {
        let mapper = build_run_mapper(proteins, &mut Vec::new())?;
        validate_merge_indices(&mapper, peptides.iter())
    }

    /// [`require_resolvable_merge_indices`](Self::require_resolvable_merge_indices)
    /// over borrowed identifications.
    ///
    /// Ports the pointer overload, which avoids copying every identification
    /// when a caller holds them in another container. The source skips null
    /// pointers; a reference cannot be null.
    ///
    /// # Errors
    ///
    /// As the slice form.
    pub fn require_resolvable_merge_indices_of(
        proteins: &[ProteinIdentification],
        peptides: &[&PeptideIdentification],
    ) -> Result<()> {
        let mapper = build_run_mapper(proteins, &mut Vec::new())?;
        validate_merge_indices(&mapper, peptides.iter().copied())
    }

    /// Identifications as a [`PSMSchema`] table, for lossless round trips.
    ///
    /// Ports `exportToArrow(protein_identifications,
    /// peptide_identifications, export_all_psms)`. One row per exported hit
    /// (the first hit of each identification, or every hit with
    /// `export_all_psms`); identifications without hits add no row but are
    /// still counted by `peptide_identification_index`. Per row:
    ///
    /// * `sequence` (unmodified), `peptidoform` (ProForma canonical) and
    ///   `modifications`, grouped by name in name order, each with its
    ///   `UNIMOD:<id>` accession (null when there is none) and its positions
    ///   as `<residue>.<1-based position>`, `N-term.0` and
    ///   `C-term.<length + 1>`, each with a null score;
    /// * `precursor_charge`; `calculated_mz` (the sequence's m/z at that
    ///   charge, null for a charge <= 0); `observed_mz` and `rt` of the
    ///   identification (null when unset or NaN);
    /// * `posterior_error_probability`: the main score when the
    ///   identification's score type is a PEP, else the first PEP-named
    ///   metadata value of its first hit (looked up on each hit), else null;
    /// * `is_decoy`: whether `target_decoy` starts with `decoy`, null without
    ///   the annotation;
    /// * `additional_scores` and `psm_metavalues`: the hit's metadata without
    ///   the keys that have columns (`target_decoy`, `predicted_RT`,
    ///   `predicted_rt`, `ion_mobility`, `IM`, `scan`, `reference_file_name`,
    ///   `fragment_annotation`); known numeric scores go to
    ///   `additional_scores`, everything else, the rank included, to
    ///   `psm_metavalues` (value as `DataValue::toString`, a type token);
    /// * `protein_accessions`: every evidence with its flanking residues and
    ///   positions (null when unknown);
    /// * `predicted_rt` (`predicted_RT`, else `predicted_rt`),
    ///   `ion_mobility` (the hit's `ion_mobility` or `IM`, else the
    ///   identification's);
    /// * `reference_file_name`: the full path of the identification's origin
    ///   file, deliberately not stemmed (null when unresolved);
    ///   `cv_params`: null;
    /// * `spectrum_reference` and `scan` (its scan number; both null without
    ///   a reference, `scan` also when the reference has no number);
    /// * `score`, `score_type`, `higher_score_better`, `hit_index` (the hit's
    ///   position in its identification), `peptide_identification_index`,
    ///   `spectrum_metavalues` (without `spectrum_reference`,
    ///   `ion_mobility` and `IM`), `run_identifier`;
    /// * the fragment annotations as `mz_array`, `intensity_array`,
    ///   `charge_array` and `ion_type_array`, null without annotations.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when runs share an identifier, a metadata
    /// value of a numeric column is not numeric (the source reinterprets a
    /// string's storage as a double), a sequence has no ProForma form or no
    /// m/z, or more than `i32::MAX` identifications are given;
    /// [`Error::MissingInformation`] for a merged run without usable
    /// `id_merge_index` (see
    /// [`require_resolvable_merge_indices`](Self::require_resolvable_merge_indices)).
    pub fn export_to_arrow(
        proteins: &[ProteinIdentification],
        peptides: &[PeptideIdentification],
        export_all_psms: bool,
    ) -> Result<QpxOutput<RecordBatch>> {
        let mut warnings = Vec::new();
        let mapper = build_run_mapper(proteins, &mut warnings)?;
        validate_merge_indices(&mapper, peptides.iter())?;
        int32(peptides.len(), "identification count")?;
        let schema = PSMSchema::schema();

        let mut sequences = Vec::new();
        let mut peptidoforms = Vec::new();
        let mut modifications = Vec::new();
        let mut charges = Vec::new();
        let mut peps = Vec::new();
        let mut decoys = Vec::new();
        let mut calculated = Vec::new();
        let mut observed = Vec::new();
        let mut scores = Vec::new();
        let mut evidences: Vec<Vec<&PeptideEvidence>> = Vec::new();
        let mut predicted = Vec::new();
        let mut reference_files = Vec::new();
        let mut scan = Vec::new();
        let mut rts = Vec::new();
        let mut mobilities = Vec::new();
        let mut spectrum_references = Vec::new();
        let mut main_scores = Vec::new();
        let mut score_types = Vec::new();
        let mut higher_better = Vec::new();
        let mut hit_indices = Vec::new();
        let mut identification_indices = Vec::new();
        let mut psm_metavalues = Vec::new();
        let mut spectrum_metas = Vec::new();
        let mut run_identifiers = Vec::new();
        let mut hits_out: Vec<&PeptideHit> = Vec::new();

        for (p_id, peptide) in peptides.iter().enumerate() {
            if peptide.hits.is_empty() {
                continue;
            }
            let p_id = int32(p_id, "identification index")?;
            let pep_score = find_peptide_score(peptide, ScoreType::PosteriorErrorProbability);
            let resolved = mapper.primary_ms_run_path(peptide)?;
            let spectrum_reference = peptide.spectrum_reference();
            let spectrum_meta = spectrum_metavalues(peptide);
            for (hit_index, hit) in exported_hits(peptide, export_all_psms).iter().enumerate() {
                let sequence = &hit.sequence;
                sequences.push(Some(sequence.as_str().to_owned()));
                peptidoforms.push(Some(canonical_peptidoform(sequence)?));
                modifications.push(if sequence.is_modified() {
                    modification_entries(sequence)?
                } else {
                    Vec::new()
                });
                charges.push(hit.charge);
                peps.push(posterior_error_probability(hit, pep_score.as_ref())?);
                decoys.push(decoy_annotation(hit));
                calculated.push(if hit.charge > 0 {
                    Some(sequence_mz(sequence, hit.charge)?)
                } else {
                    None
                });
                observed.push(peptide.mz.filter(|mz| !mz.is_nan()));
                let (hit_scores, rest) = split_scores(hit);
                scores.push(hit_scores);
                psm_metavalues.push(Some(rest));
                evidences.push(hit.evidences.iter().collect());
                predicted.push(first_numeric(
                    &hit.metadata,
                    &["predicted_RT", "predicted_rt"],
                )?);
                reference_files.push((!resolved.is_empty()).then(|| resolved.to_string()));
                if spectrum_reference.is_empty() {
                    spectrum_references.push(None);
                    scan.push(None);
                } else {
                    spectrum_references.push(Some(spectrum_reference.clone()));
                    scan.push(qpx_scan_components(&spectrum_reference).first().copied());
                }
                rts.push(peptide.rt.filter(|rt| !rt.is_nan()));
                mobilities.push(ion_mobility(hit, peptide)?);
                main_scores.push(hit.score);
                score_types.push(Some(peptide.score_type.clone()));
                higher_better.push(peptide.higher_score_better);
                hit_indices.push(int32(hit_index, "hit index")?);
                identification_indices.push(p_id);
                spectrum_metas.push(Some(spectrum_meta.clone()));
                run_identifiers.push(Some(peptide.identifier.clone()));
                hits_out.push(hit);
            }
        }

        let rows = sequences.len();
        let accessions_type = field_type(&schema, PSMSchema::PROTEIN_ACCESSIONS)?;
        let flat: Vec<&PeptideEvidence> = evidences.iter().flatten().copied().collect();
        let flank = |residue: FlankingResidue| match residue {
            FlankingResidue::Unknown => None,
            other => Some(other.code().to_string()),
        };
        let position = |value: Option<usize>| -> Result<Option<i32>> {
            value.map(|at| int32(at, "evidence position")).transpose()
        };
        let starts = flat
            .iter()
            .map(|evidence| position(evidence.start))
            .collect::<Result<Vec<_>>>()?;
        let ends = flat
            .iter()
            .map(|evidence| position(evidence.end))
            .collect::<Result<Vec<_>>>()?;
        let accessions = list_array(
            &accessions_type,
            &evidences
                .iter()
                .map(|row| Some(row.len()))
                .collect::<Vec<_>>(),
            struct_array(
                list_element(&accessions_type)?,
                vec![
                    utf8(
                        flat.iter()
                            .map(|evidence| Some(evidence.protein_accession.clone()))
                            .collect(),
                    ),
                    utf8(
                        flat.iter()
                            .map(|evidence| flank(evidence.aa_before))
                            .collect(),
                    ),
                    utf8(
                        flat.iter()
                            .map(|evidence| flank(evidence.aa_after))
                            .collect(),
                    ),
                    Arc::new(Int32Array::from(starts)),
                    Arc::new(Int32Array::from(ends)),
                ],
            )?,
        )?;
        let [mz_array, intensity_array, charge_array, ion_type_array] = fragment_arrays(
            &schema,
            [
                PSMSchema::MZ_ARRAY,
                PSMSchema::INTENSITY_ARRAY,
                PSMSchema::CHARGE_ARRAY,
                PSMSchema::ION_TYPE_ARRAY,
            ],
            &hits_out,
        )?;
        let metavalues_type = field_type(&schema, PSMSchema::PSM_METAVALUES)?;
        let columns: Vec<ArrayRef> = vec![
            utf8(sequences),
            utf8(peptidoforms),
            modifications_array(
                &field_type(&schema, PSMSchema::MODIFICATIONS)?,
                &modifications,
                false,
            )?,
            Arc::new(Int32Array::from(charges)),
            Arc::new(Float64Array::from(peps)),
            Arc::new(BooleanArray::from(decoys)),
            Arc::new(Float64Array::from(calculated)),
            Arc::new(Float64Array::from(observed)),
            additional_scores_array(&field_type(&schema, PSMSchema::ADDITIONAL_SCORES)?, &scores)?,
            accessions,
            Arc::new(Float64Array::from(predicted)),
            utf8(reference_files),
            new_null_array(&DataType::Utf8, rows),
            Arc::new(Int32Array::from(scan)),
            Arc::new(Float64Array::from(rts)),
            Arc::new(Float64Array::from(mobilities)),
            utf8(spectrum_references),
            Arc::new(Float64Array::from(main_scores)),
            utf8(score_types),
            Arc::new(BooleanArray::from(higher_better)),
            Arc::new(Int32Array::from(hit_indices)),
            Arc::new(Int32Array::from(identification_indices)),
            metavalues_array(&metavalues_type, &psm_metavalues)?,
            metavalues_array(
                &field_type(&schema, PSMSchema::SPECTRUM_METAVALUES)?,
                &spectrum_metas,
            )?,
            utf8(run_identifiers),
            mz_array,
            intensity_array,
            charge_array,
            ion_type_array,
        ];
        Ok(QpxOutput {
            value: finish_table(schema, columns)?,
            warnings,
        })
    }

    /// Identifications as the QPX psm view ([`QPXPSMSchema`]).
    ///
    /// Ports `exportPSMsToQPXArrow(protein_identifications,
    /// peptide_identifications, export_all_psms, feature_links)`. Unlike
    /// [`export_to_arrow`](Self::export_to_arrow) this is the cross-tool
    /// exchange layout, not a lossless one. One row per exported hit, with:
    ///
    /// * `psm_id`, the opaque identity [`qpx_identity::psm_id`] derives from
    ///   the row's persisted `run_file_name`, `scan`, `peptidoform` and
    ///   `charge`; `feature_id`, looked up by `psm_id` in `feature_links`
    ///   (the psm-to-feature edge the consensus-map exporter collects), null
    ///   when the identification has no feature or no links are given;
    /// * `sequence`, `peptidoform` (ProForma canonical), `modifications`
    ///   (grouped by name, `int32` positions with 0 for the N-terminus and
    ///   `length + 1` for the C-terminus, the residue as `amino_acid` or null
    ///   at a terminus, an empty `scores` list), `charge` (`int16`, wrapping
    ///   as the source's cast);
    /// * `posterior_error_probability` as in the internal layout; `is_decoy`
    ///   (`false` without a `target_decoy` annotation);
    /// * `calculated_mz` (0 for a charge <= 0), `observed_mz` (0 when unset
    ///   or NaN), both `float32`, and `mass_error_ppm`, `(observed -
    ///   calculated) / calculated * 1e6`, null unless the charge and both m/z
    ///   are positive;
    /// * `additional_scores`: the known numeric scores only (other metadata
    ///   is not exported);
    /// * `predicted_rt`, `rt`, `ion_mobility` as `float32`,
    ///   `missed_cleavages` (`int16`, from an integer `missed_cleavages`
    ///   value), `protein_accessions` (every evidence's accession);
    /// * `run_file_name`: the stem of the hit's `reference_file_name` or
    ///   `run_file_name`, else the identification's `reference_file_name`,
    ///   else its resolved origin file ([`qpx_psm_run_file_name`]);
    /// * `scan`: the spectrum reference's scan number as a one-element list,
    ///   empty without one;
    /// * the fragment annotations as in the internal layout; `cv_params`,
    ///   `cross_links` and `ion_mobility_array` null.
    ///
    /// The returned warnings name `run_file_name`s that several distinct
    /// origin files share, and count the rows whose origin file could not be
    /// resolved: their `run_file_name` is empty, which [`QpxValueValidation`]
    /// refuses when the table is written.
    ///
    /// # Errors
    ///
    /// As [`export_to_arrow`](Self::export_to_arrow); also for a
    /// `missed_cleavages` value that is not an integer (the source throws
    /// `Exception::ConversionError`).
    pub fn export_psms_to_qpx_arrow(
        proteins: &[ProteinIdentification],
        peptides: &[PeptideIdentification],
        export_all_psms: bool,
        feature_links: Option<&FeatureLinks>,
    ) -> Result<QpxOutput<RecordBatch>> {
        let mut warnings = Vec::new();
        let mapper = build_run_mapper(proteins, &mut warnings)?;
        validate_merge_indices(&mapper, peptides.iter())?;
        warnings.extend(stem_collision_warnings(&mapper));
        let borrowed: Vec<&PeptideIdentification> = peptides.iter().collect();
        let (table, unattributable) =
            build_qpx_psm_table(&mapper, &borrowed, export_all_psms, feature_links)?;
        if unattributable > 0 {
            warnings.push(unattributable_warning(unattributable));
        }
        Ok(QpxOutput {
            value: table,
            warnings,
        })
    }

    /// Identifications as a QPX `psm_file` Parquet file.
    ///
    /// Ports `exportToParquet(protein_identifications,
    /// peptide_identifications, filename, export_all_psms, config,
    /// feature_links)`: [`export_psms_to_qpx_arrow`](Self::export_psms_to_qpx_arrow),
    /// then [`export_table_to_parquet`](Self::export_table_to_parquet) with the
    /// `scan_format` of the identifications' spectrum references. Returns the
    /// export's warnings, and the warning that the references mix scan
    /// conventions (the `scan_format` is then omitted).
    ///
    /// # Errors
    ///
    /// As the two functions it runs; nothing is written on error.
    pub fn export_to_parquet(
        proteins: &[ProteinIdentification],
        peptides: &[PeptideIdentification],
        filename: impl AsRef<Path>,
        export_all_psms: bool,
        config: &ParquetWriteConfig,
        feature_links: Option<&FeatureLinks>,
    ) -> Result<Vec<String>> {
        let QpxOutput {
            value: table,
            mut warnings,
        } = Self::export_psms_to_qpx_arrow(proteins, peptides, export_all_psms, feature_links)?;
        let format = scan_format(peptides.iter(), &mut warnings);
        Self::export_table_to_parquet(&table, filename, config, format)?;
        Ok(warnings)
    }

    /// Write a QPX psm table that is already built to a QPX `psm_file`.
    ///
    /// Ports the table overload of `exportToParquet(table, filename, config,
    /// scan_format)`. The table must follow [`QPXPSMSchema`] exactly (for
    /// example, from [`export_psms_to_qpx_arrow`](Self::export_psms_to_qpx_arrow)).
    /// Its schema metadata is replaced by the QPX file metadata
    /// ([`qpx_file_metadata`] for `psm_file`: `qpx_version`, `file_type`,
    /// a fresh `uuid` and `creation_date`, `compression_format`,
    /// `primary_key`, `identity_composite`), plus `scan_format` when it is
    /// not empty -- a built table no longer carries the native IDs, so the
    /// caller supplies it ([`qpx_scan_format_of_all`]);
    /// when empty the key is omitted rather than guessed. The table is then
    /// written by [`write_table_to_parquet`], which refuses it when it
    /// violates the QPX value contract.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the table does not match
    /// [`QPXPSMSchema`], when `config` selects LZ4 (QPX has no token for it),
    /// or when the value validation or the write fails; [`Error::Io`] when
    /// the file cannot be written. Nothing is written on error.
    pub fn export_table_to_parquet(
        table: &RecordBatch,
        filename: impl AsRef<Path>,
        config: &ParquetWriteConfig,
        scan_format: &str,
    ) -> Result<()> {
        let filename = filename.as_ref();
        let validation = validate(table, &QPXPSMSchema::schema(), Mode::Strict);
        if !validation.valid {
            return Err(Error::InvalidValue(format!(
                "QPXFile: table schema does not match QPXPSMSchema ({}): {validation}",
                filename.display()
            )));
        }
        let table = with_metadata(table, qpx_psm_metadata(config, scan_format)?)?;
        write_table_to_parquet(&table, filename, config)
    }

    /// Stream identifications to a QPX `psm_file` in batches, to cap peak
    /// memory.
    ///
    /// Ports `exportToParquetStreaming(protein_identifications,
    /// peptide_identification_ptrs, filename, export_all_psms, batch_size,
    /// config, n_threads, feature_links)`, with the four defaulted knobs in
    /// [`StreamingOptions`]. The output is equivalent to
    /// [`export_to_parquet`](Self::export_to_parquet) -- the same rows,
    /// schema and QPX metadata (one `uuid` and `creation_date` for the file)
    /// -- but it is built `options.batch_size` identifications at a time and
    /// written as several row groups: each batch is split into row groups of
    /// at most `options.config.row_group_size` rows (capped at 1,048,576),
    /// and a batch without rows is one empty row group, as Parquet C++'s
    /// `WriteTable` writes it. No identifications give a valid 0-row file.
    /// One [`QpxValueValidation`] checks every batch, so a key repeated in a
    /// later batch is still refused.
    ///
    /// `options.n_threads` builds each batch in parallel partitions whose
    /// rows keep their order. The source writes each partition as its own
    /// row groups, so its file's row-group boundaries (not its rows) depend
    /// on the thread count; this writes the batch as a whole, so the file is
    /// identical for every thread count and equal to the source's serial
    /// (default) output.
    ///
    /// Everything that can be refused from the input alone -- shared run
    /// identifiers, unresolvable merge indices, LZ4 -- is refused before the
    /// file is touched. The file is written to a temporary sibling and
    /// renamed into place when complete, so a batch refused after earlier
    /// ones were written leaves no file (the source removes its partial
    /// output) and an existing file at `filename` is untouched.
    ///
    /// Returns the warnings of [`export_psms_to_qpx_arrow`](Self::export_psms_to_qpx_arrow)
    /// (the unattributable count summed over all batches) and the scan
    /// convention warning.
    ///
    /// # Errors
    ///
    /// As [`export_psms_to_qpx_arrow`](Self::export_psms_to_qpx_arrow) and
    /// [`export_table_to_parquet`](Self::export_table_to_parquet); a batch
    /// failing [`QpxValueValidation`] is [`Error::InvalidValue`] naming it;
    /// a row group size <= 0 with rows to write is [`Error::InvalidValue`].
    pub fn export_to_parquet_streaming(
        proteins: &[ProteinIdentification],
        peptides: &[&PeptideIdentification],
        filename: impl AsRef<Path>,
        options: &StreamingOptions,
        feature_links: Option<&FeatureLinks>,
    ) -> Result<Vec<String>> {
        let filename = filename.as_ref();
        let batch_size = if options.batch_size == 0 {
            DEFAULT_BATCH_SIZE
        } else {
            options.batch_size
        };
        // The source asks OpenMP for every core whenever the count is not
        // positive.
        let threads = if options.n_threads <= 0 {
            Threads::all()
        } else {
            Threads::from_cli(i64::from(options.n_threads))
        };

        let mut warnings = Vec::new();
        let mapper = build_run_mapper(proteins, &mut warnings)?;
        validate_merge_indices(&mapper, peptides.iter().copied())?;
        warnings.extend(stem_collision_warnings(&mapper));
        let format = scan_format(peptides.iter().copied(), &mut warnings);
        let metadata = qpx_psm_metadata(&options.config, format)?;
        let schema = Arc::new(
            QPXPSMSchema::schema()
                .as_ref()
                .clone()
                .with_metadata(metadata),
        );
        let properties = configured_properties(
            &options.config,
            parquet::file::properties::WriterProperties::builder(),
        )?;
        let (bare, writer_options) = writer_options(&schema, properties, true);
        let rows_per_group = usize::try_from(options.config.row_group_size)
            .ok()
            .filter(|rows| *rows > 0)
            .map(|rows| rows.min(MAX_ROW_GROUP_ROWS));
        let failure = format!(
            "QPXFile: Failed to write PSM batch to {}",
            filename.display()
        );
        let fail = |error: parquet::errors::ParquetError| {
            Error::InvalidValue(format!("{failure}: {error}"))
        };

        let mut validator = QpxValueValidation::new(View::Psm);
        let mut unattributable = 0usize;
        let batches: Vec<&[&PeptideIdentification]> = if peptides.is_empty() {
            vec![&[]]
        } else {
            peptides.chunks(batch_size).collect()
        };

        crate::format::path_io::write_plain(filename, |out| {
            let writer =
                ArrowWriter::try_new_with_options(Vec::new(), bare.clone(), writer_options)
                    .map_err(fail)?;
            let (mut file, factory) = writer.into_serialized_writer().map_err(fail)?;
            let mut group_index = 0usize;
            for batch in batches {
                let partitions = partition(batch, threads.get());
                let built = map_collect(&partitions, threads, |part| {
                    build_qpx_psm_table(&mapper, part, options.export_all_psms, feature_links)
                });
                let mut tables = Vec::with_capacity(built.len());
                for part in built {
                    let (table, missing) = part?;
                    unattributable += missing;
                    tables.push(table);
                }
                let table = if tables.len() == 1 {
                    tables.swap_remove(0)
                } else {
                    arrow_select::concat::concat_batches(&QPXPSMSchema::schema(), &tables)
                        .map_err(|error| Error::InvalidValue(format!("{failure}: {error}")))?
                };
                let validation = validator.validate(&table);
                if !validation.valid {
                    return Err(Error::InvalidValue(format!(
                        "QPXFile: refusing invalid QPX PSM batch for {}: {validation}",
                        filename.display()
                    )));
                }
                let table = RecordBatch::try_new(bare.clone(), table.columns().to_vec())
                    .map_err(|error| Error::InvalidValue(format!("{failure}: {error}")))?;
                if table.num_rows() == 0 {
                    write_row_group(&mut file, &factory, group_index, &table, &failure)?;
                    group_index += 1;
                } else {
                    let Some(rows) = rows_per_group else {
                        return Err(Error::InvalidValue(format!(
                            "{failure}: Invalid: chunk size per row_group must be greater than 0"
                        )));
                    };
                    let mut offset = 0;
                    while offset < table.num_rows() {
                        let length = rows.min(table.num_rows() - offset);
                        write_row_group(
                            &mut file,
                            &factory,
                            group_index,
                            &table.slice(offset, length),
                            &failure,
                        )?;
                        group_index += 1;
                        offset += length;
                        out.write_all(&std::mem::take(file.inner_mut()))?;
                    }
                }
                out.write_all(&std::mem::take(file.inner_mut()))?;
            }
            let tail = file.into_inner().map_err(fail)?;
            out.write_all(&tail)?;
            Ok(())
        })?;
        if unattributable > 0 {
            warnings.push(unattributable_warning(unattributable));
        }
        Ok(warnings)
    }

    /// Append the identifications of a [`PSMSchema`] table.
    ///
    /// [`import_from_arrow_with_registry`](Self::import_from_arrow_with_registry)
    /// with a copy of the global modification database, which the
    /// peptidoforms' definitions may extend and which is then dropped.
    ///
    /// # Errors
    ///
    /// As the registry form.
    pub fn import_from_arrow(
        table: &RecordBatch,
        proteins: &mut Vec<ProteinIdentification>,
        peptides: &mut Vec<PeptideIdentification>,
    ) -> Result<Vec<String>> {
        let mut registry = ModificationsDB::global().clone();
        Self::import_from_arrow_with_registry(table, proteins, peptides, &mut registry)
    }

    /// Append the identifications of a [`PSMSchema`] table, resolving
    /// peptidoforms against `registry`.
    ///
    /// Ports `importFromArrow(table, protein_identifications,
    /// peptide_identifications)`. The table is checked in subset mode
    /// (unknown columns are ignored; a known one must have its type), even
    /// when it has no rows. Rows are grouped into identifications by
    /// `peptide_identification_index` -- by value, not adjacency -- and the
    /// identifications are appended in first-seen order. The first row of an
    /// identification sets its score type, `run_identifier`, score
    /// direction (the row's `higher_score_better`, else that of the known run
    /// with its identifier, else higher-is-better), RT, m/z, spectrum
    /// reference and `spectrum_metavalues`. Every row adds one hit:
    ///
    /// * the sequence from `peptidoform`, converted best-effort (a
    ///   peptidoform that cannot be represented completely adds a warning
    ///   listing what was skipped; one that cannot be parsed or converted
    ///   adds a warning and falls back to `sequence`);
    /// * a warning when the stored `calculated_mz` differs from the
    ///   reconstructed sequence's by more than 1e-4 Da;
    /// * charge, score, `is_decoy` as `target_decoy`, evidences,
    ///   `additional_scores` as metadata, `predicted_RT`, `ion_mobility`,
    ///   `scan`, `reference_file_name`, and `psm_metavalues` (without the
    ///   keys those columns restore); a `rank` entry becomes the hit's
    ///   [`rank`](PeptideHit::rank). `hit_index` is not consulted.
    ///
    /// Finally every `run_identifier` of an appended identification that no
    /// run in `proteins` has gets a shell run with that identifier, score
    /// type and direction, and the first non-empty `reference_file_name` of
    /// its rows as primary MS run path.
    ///
    /// Returns the warnings. `proteins`, `peptides` and `registry` are
    /// unchanged on error.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an incompatible schema, a column of the
    /// wrong type, a `sequence` that cannot be parsed (the source throws
    /// from `AASequence::fromString` too), a flanking residue that is not
    /// an uppercase letter or `X`/`[`/`]`, or a negative position other than
    /// -1 (unknown); [`Error::MissingInformation`] when
    /// `peptide_identification_index`, `precursor_charge`, `score` or
    /// `score_type` is missing.
    pub fn import_from_arrow_with_registry(
        table: &RecordBatch,
        proteins: &mut Vec<ProteinIdentification>,
        peptides: &mut Vec<PeptideIdentification>,
        registry: &mut ModificationsDB,
    ) -> Result<Vec<String>> {
        let mut warnings = Vec::new();
        let validation = validate(table, &PSMSchema::schema(), Mode::Subset);
        if !validation.valid {
            return Err(Error::InvalidValue(format!(
                "QPXFile::importFromArrow: Incompatible PSM schema: {validation}"
            )));
        }
        if table.num_rows() == 0 {
            return Ok(warnings);
        }
        let column = |name: &str| helpers::get_column(table, name, false);
        let (Some(p_ids), Some(charges), Some(score), Some(score_type)) = (
            column(PSMSchema::PEPTIDE_IDENTIFICATION_INDEX),
            column(PSMSchema::PRECURSOR_CHARGE),
            column(PSMSchema::SCORE),
            column(PSMSchema::SCORE_TYPE),
        ) else {
            return Err(Error::MissingInformation(
                "QPXFile::importFromArrow: Missing required PSM columns".into(),
            ));
        };
        let peptidoform_col = column(PSMSchema::PEPTIDOFORM);
        let sequence_col = column(PSMSchema::SEQUENCE);
        let calculated_col = column(PSMSchema::CALCULATED_MZ);
        let rt_col = column(PSMSchema::RT);
        let mz_col = column(PSMSchema::OBSERVED_MZ);
        let reference_col = column(PSMSchema::SPECTRUM_REFERENCE);
        let run_col = column(PSMSchema::RUN_IDENTIFIER);
        let decoy_col = column(PSMSchema::IS_DECOY);
        let accessions_col = column(PSMSchema::PROTEIN_ACCESSIONS);
        let scores_col = column(PSMSchema::ADDITIONAL_SCORES);
        let psm_meta_col = column(PSMSchema::PSM_METAVALUES);
        let spectrum_meta_col = column(PSMSchema::SPECTRUM_METAVALUES);
        let predicted_col = column(PSMSchema::PREDICTED_RT);
        let mobility_col = column(PSMSchema::ION_MOBILITY);
        let better_col = column(PSMSchema::HIGHER_SCORE_BETTER);
        let scan_col = column(PSMSchema::SCAN);
        let file_col = column(PSMSchema::REFERENCE_FILE_NAME);

        let mut better_by_run: HashMap<&str, bool> = HashMap::new();
        for run in proteins.iter() {
            better_by_run.insert(run.identifier.as_str(), run.higher_score_better);
        }
        let mut staged_registry = registry.clone();
        let mut added: Vec<PeptideIdentification> = Vec::new();
        let mut index_of: HashMap<i32, usize> = HashMap::new();
        let mut file_by_run: HashMap<String, String> = HashMap::new();

        for row in 0..table.num_rows() {
            let p_id = get_int32_value(Some(p_ids.as_ref()), row, -1)?;
            let at = match index_of.get(&p_id) {
                Some(at) => *at,
                None => {
                    let mut peptide = PeptideIdentification {
                        score_type: get_string_value(Some(score_type.as_ref()), row)?,
                        ..PeptideIdentification::default()
                    };
                    let run_present = !is_null(run_col.as_deref(), row);
                    if run_present {
                        peptide.identifier = get_string_value(run_col.as_deref(), row)?;
                    }
                    peptide.higher_score_better = if !is_null(better_col.as_deref(), row) {
                        get_bool_value(better_col.as_deref(), row, true)?
                    } else if run_present {
                        better_by_run
                            .get(peptide.identifier.as_str())
                            .copied()
                            .unwrap_or(true)
                    } else {
                        true
                    };
                    if !is_null(rt_col.as_deref(), row) {
                        peptide.rt = Some(get_double_value(rt_col.as_deref(), row, 0.0)?);
                    }
                    if !is_null(mz_col.as_deref(), row) {
                        peptide.mz = Some(get_double_value(mz_col.as_deref(), row, 0.0)?);
                    }
                    if !is_null(reference_col.as_deref(), row) {
                        peptide.set_spectrum_reference(get_string_value(
                            reference_col.as_deref(),
                            row,
                        )?);
                    }
                    read_meta_values(
                        spectrum_meta_col.as_deref(),
                        row,
                        &mut peptide.metadata,
                        &[],
                    )?;
                    index_of.insert(p_id, added.len());
                    added.push(peptide);
                    added.len() - 1
                }
            };

            let mut hit = PeptideHit {
                charge: get_int32_value(Some(charges.as_ref()), row, 0)?,
                ..PeptideHit::default()
            };
            let mut sequence_set = false;
            let mut peptidoform = String::new();
            if !is_null(peptidoform_col.as_deref(), row) {
                peptidoform = get_string_value(peptidoform_col.as_deref(), row)?;
                if !peptidoform.is_empty() {
                    match convert_peptidoform(&peptidoform, &mut staged_registry, &mut warnings) {
                        Ok(sequence) => {
                            hit.sequence = sequence;
                            sequence_set = true;
                        }
                        Err(message) => warnings.push(format!(
                            "QPXFile: failed to parse peptidoform '{peptidoform}' ({message}); \
                             falling back to the unmodified sequence column."
                        )),
                    }
                }
            }
            if !sequence_set && !is_null(sequence_col.as_deref(), row) {
                hit.sequence = AASequence::parse_with_registry(
                    &get_string_value(sequence_col.as_deref(), row)?,
                    &staged_registry,
                )?;
            }
            if !is_null(calculated_col.as_deref(), row)
                && hit.charge != 0
                && !hit.sequence.is_empty()
            {
                let stored = get_double_value(calculated_col.as_deref(), row, 0.0)?;
                // A sequence without known chemistry has no m/z to compare;
                // the source computes one regardless.
                if let Ok(reconstructed) = sequence_mz(&hit.sequence, hit.charge) {
                    let difference = (reconstructed - stored).abs();
                    if difference > 1e-4 {
                        warnings.push(format!(
                            "QPXFile: reconstructed sequence for peptidoform '{peptidoform}' has \
                             calculated m/z {}, but the row stores {} (difference {} Da). \
                             Unsupported modification data may have been lost.",
                            format_g(reconstructed),
                            format_g(stored),
                            format_g(difference)
                        ));
                    }
                }
            }
            hit.score = get_double_value(Some(score.as_ref()), row, 0.0)?;
            if !is_null(decoy_col.as_deref(), row) {
                let decoy = get_bool_value(decoy_col.as_deref(), row, false)?;
                hit.metadata.insert(
                    "target_decoy".into(),
                    if decoy { "decoy" } else { "target" }.into(),
                );
            }
            if let Some(column) = accessions_col.as_deref() {
                read_evidences(column, row, &mut hit)?;
            }
            if let Some(column) = scores_col.as_deref() {
                read_additional_scores(column, row, &mut hit.metadata)?;
            }
            if !is_null(predicted_col.as_deref(), row) {
                hit.metadata.insert(
                    "predicted_RT".into(),
                    MetaValue::source_float(get_double_value(predicted_col.as_deref(), row, 0.0)?),
                );
            }
            if !is_null(mobility_col.as_deref(), row) {
                hit.metadata.insert(
                    "ion_mobility".into(),
                    MetaValue::source_float(get_double_value(mobility_col.as_deref(), row, 0.0)?),
                );
            }
            if !is_null(scan_col.as_deref(), row) {
                hit.metadata.insert(
                    "scan".into(),
                    i64::from(get_int32_value(scan_col.as_deref(), row, 0)?).into(),
                );
            }
            if !is_null(file_col.as_deref(), row) {
                let file = get_string_value(file_col.as_deref(), row)?;
                if !file.is_empty() {
                    file_by_run
                        .entry(added[at].identifier.clone())
                        .or_insert_with(|| file.clone());
                }
                hit.metadata
                    .insert("reference_file_name".into(), file.into());
            }
            read_meta_values(
                psm_meta_col.as_deref(),
                row,
                &mut hit.metadata,
                &IMPORT_EXCLUDED_HIT_KEYS,
            )?;
            if let Some(rank) = hit
                .metadata
                .get(RANK)
                .and_then(|value| value.as_i64().ok())
                .and_then(|rank| u32::try_from(rank).ok())
            {
                hit.rank = rank;
                hit.metadata.remove(RANK);
            }
            added[at].hits.push(hit);
        }

        let mut known: HashSet<String> =
            proteins.iter().map(|run| run.identifier.clone()).collect();
        let mut shells = Vec::new();
        for peptide in &added {
            let identifier = &peptide.identifier;
            if !identifier.is_empty() && known.insert(identifier.clone()) {
                shells.push(ProteinIdentification {
                    identifier: identifier.clone(),
                    score_type: peptide.score_type.clone(),
                    higher_score_better: peptide.higher_score_better,
                    primary_ms_run_paths: file_by_run
                        .get(identifier)
                        .map(|file| vec![file.clone()])
                        .unwrap_or_default(),
                    ..ProteinIdentification::default()
                });
            }
        }
        peptides.extend(added);
        proteins.extend(shells);
        *registry = staged_registry;
        Ok(warnings)
    }
}

/// Contiguous partitions of `items`, at most `count` and never empty unless
/// `items` is: the first `len % count` get one more element, as the source
/// splits a batch across its OpenMP workers.
fn partition<T>(items: &[T], count: usize) -> Vec<&[T]> {
    let parts = count.min(items.len()).max(1);
    let base = items.len() / parts;
    let extra = items.len() % parts;
    let mut result = Vec::with_capacity(parts);
    let mut start = 0;
    for part in 0..parts {
        let length = base + usize::from(part < extra);
        result.push(&items[start..start + length]);
        start += length;
    }
    result
}

/// Parse and convert one peptidoform best-effort, adding the source's
/// warning when the conversion skips data; the error text on failure.
fn convert_peptidoform(
    text: &str,
    registry: &mut ModificationsDB,
    warnings: &mut Vec<String>,
) -> std::result::Result<AASequence, String> {
    let parsed = Peptidoform::parse(text).map_err(|error| error.to_string())?;
    let mut staged = registry.clone();
    let issues = parsed
        .aa_sequence_conversion_issues(&mut staged)
        .map_err(|error| error.to_string())?;
    // Logged before the conversion, as the source does, so a conversion that
    // then fails still reports what it could not represent.
    if !issues.value.is_empty() {
        let mut line = format!(
            "QPXFile: peptidoform '{text}' cannot be represented completely as an AASequence. \
             BEST_EFFORT conversion will skip unsupported modification data:"
        );
        for issue in &issues.value {
            line.push(' ');
            line.push_str(&issue.description);
            line.push(';');
        }
        warnings.push(line);
    }
    let sequence = parsed
        .to_aa_sequence_with_policy(ConversionPolicy::BestEffort, &mut staged)
        .map_err(|error| error.to_string())?;
    *registry = staged;
    Ok(sequence.value)
}

fn flanking(text: Option<&str>) -> Result<FlankingResidue> {
    match text.and_then(|text| text.chars().next()) {
        None => Ok(FlankingResidue::Unknown),
        Some(code) => FlankingResidue::from_code(code).map_err(|_| {
            Error::InvalidValue(format!(
                "QPXFile::importFromArrow: invalid flanking residue '{code}'"
            ))
        }),
    }
}

fn evidence_position(value: Option<i32>) -> Result<Option<usize>> {
    match value {
        None | Some(-1) => Ok(None),
        Some(at) => usize::try_from(at).map(Some).map_err(|_| {
            Error::InvalidValue(format!(
                "QPXFile::importFromArrow: invalid evidence position {at}"
            ))
        }),
    }
}

/// The evidences of one `protein_accessions` row.
fn read_evidences(column: &dyn Array, row: usize, hit: &mut PeptideHit) -> Result<()> {
    let Some((values, start, end)) = helpers::list_range(column, row)? else {
        return Ok(());
    };
    let entries = values.as_struct_opt().ok_or_else(|| {
        Error::InvalidValue(
            "QPXFile::importFromArrow: protein_accessions is not a list of structs".into(),
        )
    })?;
    let child = |name: &str| {
        entries.column_by_name(name).cloned().ok_or_else(|| {
            Error::InvalidValue(format!(
                "QPXFile::importFromArrow: protein_accessions has no '{name}'"
            ))
        })
    };
    let (accession, before, after, first, last) = (
        child("accession")?,
        child("aa_before")?,
        child("aa_after")?,
        child("start")?,
        child("end")?,
    );
    let text = |array: &ArrayRef, index: usize| -> Result<Option<String>> {
        if array.is_null(index) {
            Ok(None)
        } else {
            get_string_value(Some(array.as_ref()), index).map(Some)
        }
    };
    let int = |array: &ArrayRef, index: usize| -> Result<Option<i32>> {
        if array.is_null(index) {
            return Ok(None);
        }
        array
            .as_primitive_opt::<Int32Type>()
            .map(|values| Some(values.value(index)))
            .ok_or_else(|| {
                Error::InvalidValue(
                    "QPXFile::importFromArrow: evidence positions are not int32".into(),
                )
            })
    };
    for index in start..end {
        hit.evidences.push(PeptideEvidence {
            protein_accession: text(&accession, index)?.unwrap_or_default(),
            aa_before: flanking(text(&before, index)?.as_deref())?,
            aa_after: flanking(text(&after, index)?.as_deref())?,
            start: evidence_position(int(&first, index)?)?,
            end: evidence_position(int(&last, index)?)?,
        });
    }
    Ok(())
}

/// The `additional_scores` of one row as float metadata values.
fn read_additional_scores(column: &dyn Array, row: usize, metadata: &mut MetaInfo) -> Result<()> {
    let Some((values, start, end)) = helpers::list_range(column, row)? else {
        return Ok(());
    };
    let refuse = || {
        Error::InvalidValue(
            "QPXFile::importFromArrow: additional_scores is not list<struct<score_name: string, \
             score_value: double>>"
                .into(),
        )
    };
    let entries = values.as_struct_opt().ok_or_else(refuse)?;
    let names = entries
        .column_by_name("score_name")
        .and_then(|column| column.as_string_opt::<i32>())
        .ok_or_else(refuse)?;
    let scores = entries
        .column_by_name("score_value")
        .and_then(|column| column.as_primitive_opt::<Float64Type>())
        .ok_or_else(refuse)?;
    for index in start..end {
        // The source reads the value slot even when it is null.
        metadata.insert(
            names.value(index).to_owned(),
            MetaValue::source_float(scores.value(index)),
        );
    }
    Ok(())
}
