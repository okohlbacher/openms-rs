// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Read and write OpenSWATH Parquet assay libraries (`.oswpq`).
//!
//! Ports `ANALYSIS/OPENSWATH/TransitionParquetFile.h` and its translation
//! unit. The Parquet library format is a container with separate tables for
//! precursors and transitions; the reader materialises every row into a
//! [`LightTargetedExperiment`], the writer serialises one.
//!
//! The container is either a directory or a store-only ZIP archive of one:
//!
//! ```text
//! <library>.oswpq
//! └── library/
//!     ├── metadata.json
//!     ├── precursors.parquet
//!     └── transitions.parquet
//! ```
//!
//! `metadata.json` holds an OpenMS block with QC counts:
//!
//! ```json
//! {
//!   "openms": {
//!     "schema_version": 1,
//!     "generator": "TransitionParquetFile",
//!     "openms_version": "<version>",
//!     "build_time": "<build_time>",
//!     "tool": {"name": "OpenSwathWorkflow", "version": "<version>"},
//!     "counts": {
//!       "proteins": {"total": 0, "target": 0, "decoy": 0},
//!       "peptides": {"total": 0, "target": 0, "decoy": 0},
//!       "precursors": {"total": 0, "target": 0, "decoy": 0},
//!       "compounds": {"total": 0, "target": 0, "decoy": 0},
//!       "transitions": {"total": 0, "target": 0, "decoy": 0}
//!     },
//!     "fragment_type_counts": {
//!       "target": {"b": 0, "y": 0, "other": 0},
//!       "decoy": {"b": 0, "y": 0, "other": 0}
//!     },
//!     "charge_counts": {
//!       "precursor": {"target": {"2": 0, "3": 0}, "decoy": {"2": 0, "3": 0}},
//!       "transition": {"target": {"1": 0, "2": 0}, "decoy": {"1": 0, "2": 0}}
//!     }
//!   }
//! }
//! ```
//!
//! Required columns of `precursors.parquet` ([`OSWPrecursorSchema`]):
//! `precursor_id` (int64), `precursor_mz` (float64), `charge` (int32),
//! `library_rt` (float64). Optional: `library_drift_time` (float64),
//! `traml_id` (string), `decoy` (bool), `modified_sequence` (string),
//! `unmodified_sequence` (string), `protein_accessions` (string or
//! list<string>).
//!
//! Required columns of `transitions.parquet` ([`OSWTransitionSchema`]):
//! `transition_id` (int64), `precursor_id` (int64), `product_mz` (float64),
//! `charge` (int32), `type` (string), `ordinal` (int32), `detecting`,
//! `identifying`, `quantifying` (bool), `library_intensity` (float64),
//! `decoy` (bool). Optional: `traml_id` (string), `annotation` (string).
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_PARQUET_SUPPORT.md)
//! for the member mapping and the evidence.
//!
//! # Identity
//!
//! The persistent `precursor_id` and `transition_id` columns are the
//! operational identity of the experiment the reader returns: they become the
//! decimal strings of `LightCompound::id` and
//! `LightTransition::transition_name`, and a transition's `peptide_ref` is its
//! `precursor_id`. The `traml_id` columns are source provenance only; the
//! `*_with_source_ids` functions expose them as a [`SourceIDMapping`].
//!
//! # Differences from the source
//!
//! - **Compound order.** The source collects precursor rows in a
//!   `std::unordered_map<int64_t, …>` and emits compounds (and the proteins
//!   they introduce) in that map's iteration order, which the standard
//!   library's hash and bucket policy decide. The port emits them in
//!   precursor row order. The set of compounds and every field are the same;
//!   transitions keep row order in both.
//! - **Archive entries** are extracted to a temporary file
//!   ([`ZipArchiveFile::extract_entry_to_temp_file`]) and read from there, the
//!   source's fallback path; its first attempt, a random-access reader over a
//!   stored entry (`ZipRandomAccessFile`), reads the same bytes.
//! - **Unparseable sequences.** The source parses sequences with
//!   `AASequence::fromString` and catches only `Exception::InvalidValue` (an
//!   unknown modification): the modifications of such a compound are skipped
//!   and its `unmodified_sequence` is written empty. A `ParseError` (an
//!   unknown residue letter, an unterminated bracket) escapes and aborts the
//!   read or write. The port's [`AASequence`] does not tell the two apart, so
//!   every sequence it refuses is treated as the source treats an unknown
//!   modification.
//! - **`build_time`.** The source writes its compiler's `__DATE__, __TIME__`;
//!   this crate records no build time and writes the empty string, so the
//!   metadata of two writes of the same library are identical.
//! - **Failure.** The source resets the output experiment first and leaves it
//!   partially filled when it throws; the port leaves it (and the source-id
//!   mapping) empty on every error.
//!
//! The source contains no `#pragma omp`; the port is serial.
//!
//! [`LightTargetedExperiment`]: crate::analysis::openswath::transition_experiment::LightTargetedExperiment
//! [`OSWPrecursorSchema`]: crate::format::arrow_schema_registry::OSWPrecursorSchema
//! [`OSWTransitionSchema`]: crate::format::arrow_schema_registry::OSWTransitionSchema
//! [`SourceIDMapping`]: crate::analysis::openswath::openswath_library_id_normalizer::SourceIDMapping
//! [`ZipArchiveFile::extract_entry_to_temp_file`]: crate::format::zip_archive_file::ZipArchiveFile::extract_entry_to_temp_file
//! [`AASequence`]: crate::chemistry::AASequence

use crate::analysis::openswath::openswath_library_id_normalizer::{
    OpenSwathLibraryIDNormalizer, SourceIDMapping,
};
use crate::analysis::openswath::transition_experiment::{
    LightCompound, LightModification, LightProtein, LightTargetedExperiment, LightTransition,
};
use crate::chemistry::{AASequence, SequenceModification};
use crate::format::arrow_schema_registry::validation::{self, Mode};
use crate::format::arrow_schema_registry::{OSWPrecursorSchema, OSWTransitionSchema};
use crate::format::parquet_file::ParquetFile;
use crate::format::zip_archive_file::ZipArchiveFile;
use crate::format::container_io::{self as file, TempDir};
use crate::{Error, Result};
use arrow_array::{
    ArrayRef, BooleanArray, Float64Array, Int32Array, Int64Array, RecordBatch, StringArray,
};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Entry name of the precursor table inside a `.oswpq` container.
pub const PRECURSORS_ENTRY: &str = "library/precursors.parquet";
/// Entry name of the transition table inside a `.oswpq` container.
pub const TRANSITIONS_ENTRY: &str = "library/transitions.parquet";
/// Entry name of the library metadata inside a `.oswpq` container.
pub const METADATA_ENTRY: &str = "library/metadata.json";

/// Reader and writer of OpenSWATH Parquet assay libraries.
///
/// Stateless, as the source class is; the API is the two conversion
/// directions, each with a source-id-provenance variant.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct TransitionParquetFile;

/// Per-class QC counts of a written library (the source's file-static
/// `OpenMSLibraryStats`).
#[derive(Default)]
struct LibraryStats {
    proteins_total: usize,
    proteins_decoy: usize,
    peptides_total: usize,
    peptides_decoy: usize,
    precursors_total: usize,
    precursors_decoy: usize,
    compounds_total: usize,
    compounds_decoy: usize,
    transitions_total: usize,
    transitions_decoy: usize,
    fragment_b_target: usize,
    fragment_b_decoy: usize,
    fragment_y_target: usize,
    fragment_y_decoy: usize,
    fragment_other_target: usize,
    fragment_other_decoy: usize,
    precursor_charge_counts_target: BTreeMap<i32, usize>,
    precursor_charge_counts_decoy: BTreeMap<i32, usize>,
    transition_charge_counts_target: BTreeMap<i32, usize>,
    transition_charge_counts_decoy: BTreeMap<i32, usize>,
}

/// One precursor row as the reader keeps it (the source's `PrecursorInfo`).
struct PrecursorInfo {
    precursor_id: i64,
    precursor_mz: f64,
    library_rt: f64,
    drift_time: f64,
    charge: i32,
    modified_sequence: String,
    unmodified_sequence: String,
    protein_accessions: Vec<String>,
}

impl TransitionParquetFile {
    /// A reader/writer (the source's default constructor).
    pub fn new() -> Self {
        Self
    }

    /// Read a `.oswpq` library into `targeted_exp`.
    ///
    /// Ports the two-argument `convertParquetToTargetedExperiment`; see
    /// [`convert_parquet_to_targeted_experiment_with_source_ids`](Self::convert_parquet_to_targeted_experiment_with_source_ids).
    ///
    /// # Errors
    ///
    /// As [`convert_parquet_to_targeted_experiment_with_source_ids`](Self::convert_parquet_to_targeted_experiment_with_source_ids).
    pub fn convert_parquet_to_targeted_experiment(
        &self,
        oswpq_dir: impl AsRef<Path>,
        targeted_exp: &mut LightTargetedExperiment,
    ) -> Result<()> {
        self.convert_parquet_to_targeted_experiment_with_source_ids(oswpq_dir, targeted_exp, None)
    }

    /// Read a `.oswpq` library into `targeted_exp` and, optionally, return the
    /// source-id provenance.
    ///
    /// Ports the three-argument `convertParquetToTargetedExperiment`. Opens
    /// `library/precursors.parquet` and `library/transitions.parquet` from
    /// `oswpq_dir` (a ZIP archive or an extracted directory), validates each
    /// table against [`OSWPrecursorSchema`] / [`OSWTransitionSchema`] in
    /// subset mode (extra columns are tolerated) and materialises the rows.
    /// `targeted_exp` is reset first: pre-existing contents are discarded,
    /// not appended to.
    ///
    /// Every precursor row becomes one compound whose `id` is the decimal
    /// `precursor_id`; its `sequence` is `modified_sequence`, or
    /// `unmodified_sequence` when that is empty, and its modifications are
    /// parsed from that sequence. Every accession of `protein_accessions`
    /// (a `;`-separated string or a list) introduces one protein the first
    /// time it appears. Every transition row becomes one transition named by
    /// its decimal `transition_id`, whose `peptide_ref` is its
    /// `precursor_id` and whose precursor m/z and ion mobility are its
    /// precursor's; its fragment type is `type`, or `annotation` when `type`
    /// is empty. Missing optional values default as in the source: drift time
    /// -1, library RT and intensity 0, decoy false, detecting and quantifying
    /// true, identifying false, ordinal -1.
    ///
    /// When `source_ids` is given it is reset, and every non-empty precursor
    /// `traml_id` is mapped to and from its canonical `precursor_id`, every
    /// non-empty transition `traml_id` from its canonical `transition_id`.
    ///
    /// # Errors
    ///
    /// - [`Error::Io`] of kind `NotFound` when a required table is not in
    ///   `oswpq_dir`, where the source throws `Exception::FileNotFound` (its
    ///   header documents `Exception::MissingInformation`, which its
    ///   implementation cannot reach); [`Error::InvalidValue`] when the
    ///   archive cannot be opened or a table cannot be decoded.
    /// - [`Error::InvalidValue`] "Precursors table schema validation failed:
    ///   …" / "Transitions table schema validation failed: …", as the source.
    /// - [`Error::InvalidValue`] "Duplicate precursor_id in OSWPQ precursor
    ///   table", "Ambiguous OSWPQ precursor traml_id maps to multiple
    ///   precursor_id values" and "OSWPQ precursor_id has multiple traml_id
    ///   provenance values", as the source.
    /// - [`Error::MissingInformation`] for a null required value and
    ///   "Transition references unknown precursor_id", as the source.
    /// - Whatever [`OpenSwathLibraryIDNormalizer::validate_canonical_ids`]
    ///   reports for the result (for example a negative or duplicated
    ///   `transition_id`).
    ///
    /// [`OSWPrecursorSchema`]: crate::format::arrow_schema_registry::OSWPrecursorSchema
    /// [`OSWTransitionSchema`]: crate::format::arrow_schema_registry::OSWTransitionSchema
    pub fn convert_parquet_to_targeted_experiment_with_source_ids(
        &self,
        oswpq_dir: impl AsRef<Path>,
        targeted_exp: &mut LightTargetedExperiment,
        source_ids: Option<&mut SourceIDMapping>,
    ) -> Result<()> {
        *targeted_exp = LightTargetedExperiment::default();
        let mut ids = SourceIDMapping::default();
        let want_ids = source_ids.is_some();
        let result = read_library(oswpq_dir.as_ref(), want_ids, &mut ids);
        if let Some(out) = source_ids {
            *out = SourceIDMapping::default();
            if result.is_ok() {
                *out = ids;
            }
        }
        *targeted_exp = result?;
        Ok(())
    }

    /// Write `targeted_exp` as a `.oswpq` library.
    ///
    /// Ports the two-argument `convertLightTargetedExperimentToParquet`, the
    /// compatibility entry point: an experiment that already carries
    /// canonical ids ([`OpenSwathLibraryIDNormalizer::has_canonical_ids`]) is
    /// written as it is. One whose ids merely look canonical (decimal
    /// syntax, [`OpenSwathLibraryIDNormalizer::has_canonical_id_format`]) but
    /// break an invariant is refused rather than silently renumbered. Any
    /// other (source-style) experiment is copied, normalised once with
    /// [`OpenSwathLibraryIDNormalizer::normalize_source_ids`], and written
    /// with the resulting provenance, so the original ids land in the
    /// `traml_id` columns.
    ///
    /// # Errors
    ///
    /// [`OpenSwathLibraryIDNormalizer::validate_canonical_ids`]'s error for a
    /// malformed canonical-looking experiment, the normaliser's errors, and
    /// those of
    /// [`convert_light_targeted_experiment_to_parquet_with_source_ids`](Self::convert_light_targeted_experiment_to_parquet_with_source_ids).
    pub fn convert_light_targeted_experiment_to_parquet(
        &self,
        oswpq_path: impl AsRef<Path>,
        targeted_exp: &LightTargetedExperiment,
    ) -> Result<()> {
        if OpenSwathLibraryIDNormalizer::has_canonical_ids(targeted_exp) {
            return self.convert_light_targeted_experiment_to_parquet_with_source_ids(
                oswpq_path,
                targeted_exp,
                None,
            );
        }
        if OpenSwathLibraryIDNormalizer::has_canonical_id_format(targeted_exp) {
            // Canonical syntax with a broken invariant is malformed canonical
            // data, not source-style ids: report it instead of renumbering.
            OpenSwathLibraryIDNormalizer::validate_canonical_ids(targeted_exp)?;
        }
        let mut normalized = targeted_exp.clone();
        let source_ids = OpenSwathLibraryIDNormalizer::normalize_source_ids(&mut normalized)?;
        self.convert_light_targeted_experiment_to_parquet_with_source_ids(
            oswpq_path,
            &normalized,
            Some(&source_ids),
        )
    }

    /// Write a canonical `targeted_exp` as a `.oswpq` library, with optional
    /// source-id provenance.
    ///
    /// Ports the three-argument `convertLightTargetedExperimentToParquet`.
    /// Precursor and transition ids are preserved exactly, zero and sparse
    /// values included. When `oswpq_path` is an existing directory the tables
    /// are written below `<oswpq_path>/library/`; otherwise the layout is
    /// staged in a temporary directory, added entry by entry to the staging
    /// archive `<absolute oswpq_path>.tmp` (store-only), given its sidecar
    /// index ([`ZipArchiveFile::write_sidecar_index`]) and renamed over
    /// `oswpq_path`.
    ///
    /// One precursor row is written per compound, in compound order: its
    /// `precursor_mz` is that of the first transition referencing it, its
    /// `decoy` flag is set when any transition referencing it is a decoy, its
    /// `traml_id` is the source id from `source_ids` (else the canonical id),
    /// `modified_sequence` is the compound sequence and `unmodified_sequence`
    /// its residues, and `protein_accessions` the protein refs joined by `;`.
    /// One transition row is written per transition, in order, its
    /// `traml_id` again the source id when `source_ids` has one. Proteins
    /// only enter the QC counts (a protein id starting with `DECOY_` counts as
    /// a decoy).
    ///
    /// # Errors
    ///
    /// - [`OpenSwathLibraryIDNormalizer::validate_canonical_ids`]'s error when
    ///   the ids are not canonical, before anything is written.
    /// - [`Error::MissingInformation`] "No precursor_mz found for compound
    ///   '<id>'" for a compound no transition references, and "Transition
    ///   references unknown peptide_ref '<ref>'", as the source.
    /// - [`Error::InvalidValue`] when a built table fails its strict schema
    ///   validation or cannot be encoded, and the archive errors of
    ///   [`ZipArchiveFile`]; [`Error::Io`] when a file cannot be written (the
    ///   source's `Exception::FileNotWritable`).
    ///
    /// [`ZipArchiveFile`]: crate::format::zip_archive_file::ZipArchiveFile
    /// [`ZipArchiveFile::write_sidecar_index`]: crate::format::zip_archive_file::ZipArchiveFile::write_sidecar_index
    pub fn convert_light_targeted_experiment_to_parquet_with_source_ids(
        &self,
        oswpq_path: impl AsRef<Path>,
        targeted_exp: &LightTargetedExperiment,
        source_ids: Option<&SourceIDMapping>,
    ) -> Result<()> {
        OpenSwathLibraryIDNormalizer::validate_canonical_ids(targeted_exp)?;
        let oswpq_path = oswpq_path.as_ref();
        if file::is_directory(oswpq_path) {
            return write_library(oswpq_path, targeted_exp, source_ids);
        }
        let temp_dir = TempDir::new(false)?;
        let base_dir = temp_dir.path().join("oswpq_output");
        file::make_dir(&base_dir)?;
        write_library(&base_dir, targeted_exp, source_ids)?;

        let output_zip = file::absolute_path(oswpq_path)?;
        let mut staging_name = output_zip.clone().into_os_string();
        staging_name.push(".tmp");
        let staging_zip = PathBuf::from(staging_name);
        if file::exists(&staging_zip) {
            file::remove(&staging_zip)?;
        }
        // The staged tree holds exactly these three files. The source walks it
        // with a recursive directory iterator, whose order the file system
        // decides; on the Release oracle's ext4 that is creation order, which
        // is the order used here, so the archive is deterministic.
        for entry in [PRECURSORS_ENTRY, TRANSITIONS_ENTRY, METADATA_ENTRY] {
            ZipArchiveFile::add_or_replace_from_file(&staging_zip, entry, base_dir.join(entry))?;
        }
        ZipArchiveFile::write_sidecar_index(&staging_zip)?;
        if file::exists(&output_zip) {
            file::remove(&output_zip)?;
        }
        file::rename(&staging_zip, &output_zip, true)?;
        // The source's `TempDir` destructor removes the staging tree and
        // ignores a failure to do so; dropping the guard does the same.
        drop(temp_dir);
        Ok(())
    }
}

/// Read one table of the container, as the source's `open_table_from_entry`.
fn open_table(
    oswpq_dir: &Path,
    entry: &str,
    temp_dir: &mut Option<TempDir>,
) -> Result<RecordBatch> {
    let path = ZipArchiveFile::extract_entry_to_temp_file(oswpq_dir, entry, temp_dir)?;
    ParquetFile::read_table(path)
}

fn read_library(
    oswpq_dir: &Path,
    want_ids: bool,
    source_ids: &mut SourceIDMapping,
) -> Result<LightTargetedExperiment> {
    let mut temp_dir: Option<TempDir> = None;
    let precursors = open_table(oswpq_dir, PRECURSORS_ENTRY, &mut temp_dir)?;
    let transitions = open_table(oswpq_dir, TRANSITIONS_ENTRY, &mut temp_dir)?;

    let prec_validation = validation::validate(
        &precursors,
        OSWPrecursorSchema::schema().as_ref(),
        Mode::Subset,
    );
    if !prec_validation.valid {
        return Err(Error::InvalidValue(format!(
            "Precursors table schema validation failed: {prec_validation}"
        )));
    }
    let trans_validation = validation::validate(
        &transitions,
        OSWTransitionSchema::schema().as_ref(),
        Mode::Subset,
    );
    if !trans_validation.valid {
        return Err(Error::InvalidValue(format!(
            "Transitions table schema validation failed: {trans_validation}"
        )));
    }

    let precursor_id_col = ParquetFile::get_column(&precursors, OSWPrecursorSchema::PRECURSOR_ID)?;
    let precursor_mz_col = ParquetFile::get_column(&precursors, OSWPrecursorSchema::PRECURSOR_MZ)?;
    let charge_col = ParquetFile::get_column(&precursors, OSWPrecursorSchema::CHARGE)?;
    let library_rt_col = ParquetFile::get_column(&precursors, OSWPrecursorSchema::LIBRARY_RT)?;
    let drift_time_col =
        ParquetFile::get_optional_column(&precursors, OSWPrecursorSchema::LIBRARY_DRIFT_TIME);
    let decoy_col = ParquetFile::get_optional_column(&precursors, OSWPrecursorSchema::DECOY);
    let precursor_traml_id_col =
        ParquetFile::get_optional_column(&precursors, OSWPrecursorSchema::TRAML_ID);
    let modified_sequence_col =
        ParquetFile::get_optional_column(&precursors, OSWPrecursorSchema::MODIFIED_SEQUENCE);
    let unmodified_sequence_col =
        ParquetFile::get_optional_column(&precursors, OSWPrecursorSchema::UNMODIFIED_SEQUENCE);
    let protein_accessions_col =
        ParquetFile::get_optional_column(&precursors, OSWPrecursorSchema::PROTEIN_ACCESSIONS);

    // Row order, with an index by id; the source keeps only a hash map.
    let mut precursor_rows: Vec<PrecursorInfo> = Vec::with_capacity(precursors.num_rows());
    let mut precursor_index: HashMap<i64, usize> = HashMap::with_capacity(precursors.num_rows());
    for row in 0..precursors.num_rows() {
        let precursor_id = ParquetFile::get_int64(Some(precursor_id_col.as_ref()), row, 0, false)?;
        // Read in the source's order, so the first failing field is the same.
        let precursor_mz =
            ParquetFile::get_double(Some(precursor_mz_col.as_ref()), row, 0.0, false)?;
        // `static_cast<int>` of the int64 value, which wraps.
        let charge = ParquetFile::get_int64(Some(charge_col.as_ref()), row, 0, false)? as i32;
        let library_rt = ParquetFile::get_double(Some(library_rt_col.as_ref()), row, 0.0, true)?;
        let drift_time = ParquetFile::get_double(drift_time_col.as_deref(), row, -1.0, true)?;
        // The source stores `decoy` in the precursor record and never uses it.
        ParquetFile::get_bool(decoy_col.as_deref(), row, false, true)?;
        let info = PrecursorInfo {
            precursor_id,
            precursor_mz,
            charge,
            library_rt,
            drift_time,
            modified_sequence: ParquetFile::get_string(modified_sequence_col.as_deref(), row)?,
            unmodified_sequence: ParquetFile::get_string(unmodified_sequence_col.as_deref(), row)?,
            protein_accessions: ParquetFile::get_string_list(
                protein_accessions_col.as_deref(),
                row,
            )?,
        };
        if precursor_index
            .insert(precursor_id, precursor_rows.len())
            .is_some()
        {
            return Err(Error::InvalidValue(format!(
                "Duplicate precursor_id in OSWPQ precursor table: {precursor_id}"
            )));
        }
        precursor_rows.push(info);

        if want_ids {
            let source_id = ParquetFile::get_string(precursor_traml_id_col.as_deref(), row)?;
            if !source_id.is_empty() {
                let canonical_id = precursor_id.to_string();
                let forward = source_ids
                    .precursor_source_to_canonical
                    .entry(source_id.clone())
                    .or_insert_with(|| canonical_id.clone());
                if *forward != canonical_id {
                    return Err(Error::InvalidValue(format!(
                        "Ambiguous OSWPQ precursor traml_id maps to multiple precursor_id values: {source_id}"
                    )));
                }
                let reverse = source_ids
                    .precursor_canonical_to_source
                    .entry(canonical_id.clone())
                    .or_insert_with(|| source_id.clone());
                if *reverse != source_id {
                    return Err(Error::InvalidValue(format!(
                        "OSWPQ precursor_id has multiple traml_id provenance values: {canonical_id}"
                    )));
                }
            }
        }
    }

    let mut exp = LightTargetedExperiment::default();
    let mut seen_proteins: HashSet<String> = HashSet::new();
    exp.compounds.reserve(precursor_rows.len());
    for info in &precursor_rows {
        let mut compound = LightCompound {
            id: info.precursor_id.to_string(),
            drift_time: info.drift_time,
            rt: info.library_rt,
            charge: info.charge,
            sequence: if info.modified_sequence.is_empty() {
                info.unmodified_sequence.clone()
            } else {
                info.modified_sequence.clone()
            },
            protein_refs: info.protein_accessions.clone(),
            ..LightCompound::default()
        };
        if !compound.sequence.is_empty() {
            add_modifications(&mut compound);
        }
        exp.compounds.push(compound);
        for accession in &info.protein_accessions {
            if seen_proteins.insert(accession.clone()) {
                exp.proteins.push(LightProtein {
                    id: accession.clone(),
                    ..LightProtein::default()
                });
            }
        }
    }

    let transition_id_col =
        ParquetFile::get_column(&transitions, OSWTransitionSchema::TRANSITION_ID)?;
    let transition_precursor_id_col =
        ParquetFile::get_column(&transitions, OSWTransitionSchema::PRECURSOR_ID)?;
    let transition_traml_id_col =
        ParquetFile::get_optional_column(&transitions, OSWTransitionSchema::TRAML_ID);
    let product_mz_col = ParquetFile::get_column(&transitions, OSWTransitionSchema::PRODUCT_MZ)?;
    let fragment_charge_col = ParquetFile::get_column(&transitions, OSWTransitionSchema::CHARGE)?;
    let fragment_type_col = ParquetFile::get_column(&transitions, OSWTransitionSchema::TYPE)?;
    let fragment_annotation_col =
        ParquetFile::get_optional_column(&transitions, OSWTransitionSchema::ANNOTATION);
    let fragment_ordinal_col = ParquetFile::get_column(&transitions, OSWTransitionSchema::ORDINAL)?;
    let detecting_col = ParquetFile::get_column(&transitions, OSWTransitionSchema::DETECTING)?;
    let identifying_col = ParquetFile::get_column(&transitions, OSWTransitionSchema::IDENTIFYING)?;
    let quantifying_col = ParquetFile::get_column(&transitions, OSWTransitionSchema::QUANTIFYING)?;
    let transition_intensity_col =
        ParquetFile::get_column(&transitions, OSWTransitionSchema::LIBRARY_INTENSITY)?;
    let transition_decoy_col = ParquetFile::get_column(&transitions, OSWTransitionSchema::DECOY)?;

    exp.transitions.reserve(transitions.num_rows());
    for row in 0..transitions.num_rows() {
        let precursor_id =
            ParquetFile::get_int64(Some(transition_precursor_id_col.as_ref()), row, 0, false)?;
        let Some(&precursor_row) = precursor_index.get(&precursor_id) else {
            return Err(Error::MissingInformation(
                "Transition references unknown precursor_id".into(),
            ));
        };
        let precursor = &precursor_rows[precursor_row];
        let transition_id =
            ParquetFile::get_int64(Some(transition_id_col.as_ref()), row, 0, false)?;
        let transition_name = transition_id.to_string();
        if want_ids {
            let source_transition_id =
                ParquetFile::get_string(transition_traml_id_col.as_deref(), row)?;
            if !source_transition_id.is_empty() {
                source_ids
                    .transition_canonical_to_source
                    .entry(transition_name.clone())
                    .or_insert(source_transition_id);
            }
        }
        let fragment_type = ParquetFile::get_string(Some(fragment_type_col.as_ref()), row)?;
        let annotation = ParquetFile::get_string(fragment_annotation_col.as_deref(), row)?;

        let mut transition = LightTransition {
            transition_name,
            peptide_ref: precursor_id.to_string(),
            library_intensity: ParquetFile::get_double(
                Some(transition_intensity_col.as_ref()),
                row,
                0.0,
                true,
            )?,
            precursor_mz: precursor.precursor_mz,
            product_mz: ParquetFile::get_double(Some(product_mz_col.as_ref()), row, 0.0, false)?,
            precursor_im: precursor.drift_time,
            // `static_cast<int8_t>` / `static_cast<int16_t>`, which wrap.
            fragment_charge: ParquetFile::get_int64(
                Some(fragment_charge_col.as_ref()),
                row,
                0,
                true,
            )? as i8,
            fragment_nr: ParquetFile::get_int64(Some(fragment_ordinal_col.as_ref()), row, -1, true)?
                as i16,
            ..LightTransition::default()
        };
        transition.set_fragment_type(if fragment_type.is_empty() {
            &annotation
        } else {
            &fragment_type
        });
        transition.set_decoy(ParquetFile::get_bool(
            Some(transition_decoy_col.as_ref()),
            row,
            false,
            true,
        )?);
        transition.set_detecting_transition(ParquetFile::get_bool(
            Some(detecting_col.as_ref()),
            row,
            true,
            true,
        )?);
        transition.set_identifying_transition(ParquetFile::get_bool(
            Some(identifying_col.as_ref()),
            row,
            false,
            true,
        )?);
        transition.set_quantifying_transition(ParquetFile::get_bool(
            Some(quantifying_col.as_ref()),
            row,
            true,
            true,
        )?);
        exp.transitions.push(transition);
    }

    OpenSwathLibraryIDNormalizer::validate_canonical_ids(&exp)?;
    Ok(exp)
}

/// The UniMod record id of a sequence modification, -1 when it has none (a
/// mass tag), as `ResidueModification::getUniModRecordId`.
fn unimod_id(modification: &SequenceModification) -> i32 {
    modification
        .record_id()
        .and_then(|id| i32::try_from(id).ok())
        .unwrap_or(-1)
}

/// The source's file-static `addModifications_`: parse the compound's
/// sequence and record its terminal and residue modifications. A sequence
/// that does not parse contributes none.
fn add_modifications(compound: &mut LightCompound) {
    let Ok(sequence) = AASequence::parse(&compound.sequence) else {
        return;
    };
    if let Some(m) = sequence.n_terminal_modification() {
        compound.modifications.push(LightModification {
            location: -1,
            unimod_id: unimod_id(m),
        });
    }
    if let Some(m) = sequence.c_terminal_modification() {
        compound.modifications.push(LightModification {
            location: i32::try_from(sequence.len()).unwrap_or(i32::MAX),
            unimod_id: unimod_id(m),
        });
    }
    for index in 0..sequence.len() {
        if let Ok(Some(m)) = sequence.residue_modification(index) {
            compound.modifications.push(LightModification {
                location: i32::try_from(index).unwrap_or(i32::MAX),
                unimod_id: unimod_id(m),
            });
        }
    }
}

/// The residues of a sequence (`AASequence::toUnmodifiedString`), or the
/// empty string when it does not parse.
fn unmodified_sequence(sequence: &str) -> String {
    if sequence.is_empty() {
        return String::new();
    }
    AASequence::parse(sequence)
        .map(|s| s.as_str().to_owned())
        .unwrap_or_default()
}

fn batch(schema: arrow_schema::SchemaRef, columns: Vec<ArrayRef>) -> Result<RecordBatch> {
    RecordBatch::try_new(schema, columns)
        .map_err(|e| Error::InvalidValue(format!("Failed to build table: {e}")))
}

fn write_library(
    base_dir: &Path,
    exp: &LightTargetedExperiment,
    source_ids: Option<&SourceIDMapping>,
) -> Result<()> {
    let library_dir = base_dir.join("library");
    file::make_dir(&library_dir)?;

    // validate_canonical_ids guarantees canonical decimal, unique,
    // non-negative compound ids; parsing cannot fail, but is checked anyway.
    let mut compound_to_precursor: HashMap<&str, i64> = HashMap::with_capacity(exp.compounds.len());
    for compound in &exp.compounds {
        let id = compound.id.parse::<i64>().map_err(|_| {
            Error::InvalidValue(format!("Non-canonical precursor id '{}'", compound.id))
        })?;
        compound_to_precursor
            .entry(compound.id.as_str())
            .or_insert(id);
    }
    let mut compound_decoy: HashSet<&str> = HashSet::new();
    let mut precursor_mz: HashMap<&str, f64> = HashMap::new();
    for transition in &exp.transitions {
        if transition.decoy() {
            compound_decoy.insert(transition.peptide_ref.as_str());
        }
        precursor_mz
            .entry(transition.peptide_ref.as_str())
            .or_insert(transition.precursor_mz);
    }

    let mut stats = LibraryStats::default();
    let n = exp.compounds.len();
    let mut precursor_ids = Vec::with_capacity(n);
    let mut precursor_mzs = Vec::with_capacity(n);
    let mut charges = Vec::with_capacity(n);
    let mut library_rts = Vec::with_capacity(n);
    let mut drift_times = Vec::with_capacity(n);
    let mut decoys = Vec::with_capacity(n);
    let mut traml_ids = Vec::with_capacity(n);
    let mut modified_sequences = Vec::with_capacity(n);
    let mut unmodified_sequences = Vec::with_capacity(n);
    let mut protein_accessions = Vec::with_capacity(n);
    for compound in &exp.compounds {
        let is_decoy = compound_decoy.contains(compound.id.as_str());
        stats.compounds_total += 1;
        stats.precursors_total += 1;
        if compound.is_peptide() {
            stats.peptides_total += 1;
            if is_decoy {
                stats.peptides_decoy += 1;
            }
        }
        if is_decoy {
            stats.compounds_decoy += 1;
            stats.precursors_decoy += 1;
            *stats
                .precursor_charge_counts_decoy
                .entry(compound.charge)
                .or_insert(0) += 1;
        } else {
            *stats
                .precursor_charge_counts_target
                .entry(compound.charge)
                .or_insert(0) += 1;
        }
        precursor_ids.push(compound_to_precursor[compound.id.as_str()]);
        let Some(&mz) = precursor_mz.get(compound.id.as_str()) else {
            return Err(Error::MissingInformation(format!(
                "No precursor_mz found for compound '{}'",
                compound.id
            )));
        };
        precursor_mzs.push(mz);
        charges.push(compound.charge);
        library_rts.push(compound.rt);
        drift_times.push(compound.drift_time);
        decoys.push(is_decoy);
        let source = source_ids
            .and_then(|ids| ids.precursor_canonical_to_source.get(&compound.id))
            .unwrap_or(&compound.id);
        traml_ids.push(source.clone());
        modified_sequences.push(compound.sequence.clone());
        unmodified_sequences.push(unmodified_sequence(&compound.sequence));
        protein_accessions.push(compound.protein_refs.join(";"));
    }
    for protein in &exp.proteins {
        stats.proteins_total += 1;
        if protein.id.starts_with("DECOY_") {
            stats.proteins_decoy += 1;
        }
    }

    let precursors = batch(
        OSWPrecursorSchema::schema(),
        vec![
            Arc::new(Int64Array::from(precursor_ids)),
            Arc::new(Float64Array::from(precursor_mzs)),
            Arc::new(Int32Array::from(charges)),
            Arc::new(Float64Array::from(library_rts)),
            Arc::new(Float64Array::from(drift_times)),
            Arc::new(BooleanArray::from(decoys)),
            Arc::new(StringArray::from(traml_ids)),
            Arc::new(StringArray::from(modified_sequences)),
            Arc::new(StringArray::from(unmodified_sequences)),
            Arc::new(StringArray::from(protein_accessions)),
        ],
    )?;
    let prec_validation = validation::validate(
        &precursors,
        OSWPrecursorSchema::schema().as_ref(),
        Mode::Strict,
    );
    if !prec_validation.valid {
        return Err(Error::InvalidValue(format!(
            "Precursors table schema validation failed: {prec_validation}"
        )));
    }
    ParquetFile::write_table(&precursors, library_dir.join("precursors.parquet"))?;

    let n = exp.transitions.len();
    let mut transition_ids = Vec::with_capacity(n);
    let mut transition_precursor_ids = Vec::with_capacity(n);
    let mut transition_traml_ids = Vec::with_capacity(n);
    let mut product_mzs = Vec::with_capacity(n);
    let mut fragment_charges = Vec::with_capacity(n);
    let mut fragment_types = Vec::with_capacity(n);
    let mut annotations = Vec::with_capacity(n);
    let mut ordinals = Vec::with_capacity(n);
    let mut detecting = Vec::with_capacity(n);
    let mut identifying = Vec::with_capacity(n);
    let mut quantifying = Vec::with_capacity(n);
    let mut intensities = Vec::with_capacity(n);
    let mut transition_decoys = Vec::with_capacity(n);
    for transition in &exp.transitions {
        let Some(&precursor_ref) = compound_to_precursor.get(transition.peptide_ref.as_str())
        else {
            return Err(Error::MissingInformation(format!(
                "Transition references unknown peptide_ref '{}'",
                transition.peptide_ref
            )));
        };
        let transition_id = transition.transition_name.parse::<i64>().map_err(|_| {
            Error::InvalidValue(format!(
                "Non-canonical transition id '{}'",
                transition.transition_name
            ))
        })?;
        transition_ids.push(transition_id);
        transition_precursor_ids.push(precursor_ref);
        let source = source_ids
            .and_then(|ids| {
                ids.transition_canonical_to_source
                    .get(&transition.transition_name)
            })
            .unwrap_or(&transition.transition_name);
        transition_traml_ids.push(source.clone());
        product_mzs.push(transition.product_mz);
        fragment_charges.push(i32::from(transition.fragment_charge));
        fragment_types.push(transition.fragment_type_str());
        annotations.push(transition.annotation());
        ordinals.push(i32::from(transition.fragment_nr));
        detecting.push(transition.is_detecting_transition());
        identifying.push(transition.is_identifying_transition());
        quantifying.push(transition.is_quantifying_transition());
        intensities.push(transition.library_intensity);
        let is_decoy = transition.decoy();
        transition_decoys.push(is_decoy);

        stats.transitions_total += 1;
        let charge = i32::from(transition.fragment_charge);
        if is_decoy {
            stats.transitions_decoy += 1;
            *stats
                .transition_charge_counts_decoy
                .entry(charge)
                .or_insert(0) += 1;
        } else {
            *stats
                .transition_charge_counts_target
                .entry(charge)
                .or_insert(0) += 1;
        }
        match transition.fragment_type_str().to_ascii_lowercase().as_str() {
            "b" if is_decoy => stats.fragment_b_decoy += 1,
            "b" => stats.fragment_b_target += 1,
            "y" if is_decoy => stats.fragment_y_decoy += 1,
            "y" => stats.fragment_y_target += 1,
            _ if is_decoy => stats.fragment_other_decoy += 1,
            _ => stats.fragment_other_target += 1,
        }
    }

    let transitions = batch(
        OSWTransitionSchema::schema(),
        vec![
            Arc::new(Int64Array::from(transition_ids)),
            Arc::new(Int64Array::from(transition_precursor_ids)),
            Arc::new(StringArray::from(transition_traml_ids)),
            Arc::new(Float64Array::from(product_mzs)),
            Arc::new(Int32Array::from(fragment_charges)),
            Arc::new(StringArray::from(fragment_types)),
            Arc::new(StringArray::from(annotations)),
            Arc::new(Int32Array::from(ordinals)),
            Arc::new(BooleanArray::from(detecting)),
            Arc::new(BooleanArray::from(identifying)),
            Arc::new(BooleanArray::from(quantifying)),
            Arc::new(Float64Array::from(intensities)),
            Arc::new(BooleanArray::from(transition_decoys)),
        ],
    )?;
    let trans_validation = validation::validate(
        &transitions,
        OSWTransitionSchema::schema().as_ref(),
        Mode::Strict,
    );
    if !trans_validation.valid {
        return Err(Error::InvalidValue(format!(
            "Transitions table schema validation failed: {trans_validation}"
        )));
    }
    ParquetFile::write_table(&transitions, library_dir.join("transitions.parquet"))?;

    std::fs::write(library_dir.join("metadata.json"), library_metadata(&stats))?;
    Ok(())
}

/// `{"target":{…},"decoy":{…}}` with the charges in ascending order, as the
/// source's `jsonMapByClass_`.
fn json_map_by_class(target: &BTreeMap<i32, usize>, decoy: &BTreeMap<i32, usize>) -> String {
    let join = |map: &BTreeMap<i32, usize>| {
        map.iter()
            .map(|(charge, count)| format!("\"{charge}\":{count}"))
            .collect::<Vec<_>>()
            .join(",")
    };
    format!(
        "{{\"target\":{{{}}},\"decoy\":{{{}}}}}",
        join(target),
        join(decoy)
    )
}

/// The text of `library/metadata.json`, byte for byte the source's
/// `writeLibraryMetadata_` except for `build_time` (see the module notes).
fn library_metadata(stats: &LibraryStats) -> String {
    let version = ParquetFile::json_escape(crate::CORE_SDK_VERSION);
    let counts = |total: usize, decoy: usize| {
        format!(
            "{{\"total\": {total}, \"target\": {}, \"decoy\": {decoy}}}",
            total - decoy
        )
    };
    let mut out = String::new();
    out.push_str("{\n");
    out.push_str("  \"openms\": {\n");
    out.push_str("    \"schema_version\": 1,\n");
    out.push_str("    \"generator\": \"TransitionParquetFile\",\n");
    out.push_str(&format!("    \"openms_version\": \"{version}\",\n"));
    out.push_str("    \"build_time\": \"\",\n");
    out.push_str(&format!(
        "    \"tool\": {{\"name\": \"OpenSwathWorkflow\", \"version\": \"{version}\"}},\n"
    ));
    out.push_str("    \"counts\": {\n");
    out.push_str(&format!(
        "      \"proteins\": {},\n",
        counts(stats.proteins_total, stats.proteins_decoy)
    ));
    out.push_str(&format!(
        "      \"peptides\": {},\n",
        counts(stats.peptides_total, stats.peptides_decoy)
    ));
    out.push_str(&format!(
        "      \"precursors\": {},\n",
        counts(stats.precursors_total, stats.precursors_decoy)
    ));
    out.push_str(&format!(
        "      \"compounds\": {},\n",
        counts(stats.compounds_total, stats.compounds_decoy)
    ));
    out.push_str(&format!(
        "      \"transitions\": {}\n",
        counts(stats.transitions_total, stats.transitions_decoy)
    ));
    out.push_str("    },\n");
    out.push_str("    \"fragment_type_counts\": {\n");
    out.push_str(&format!(
        "      \"target\": {{\"b\": {}, \"y\": {}, \"other\": {}}},\n",
        stats.fragment_b_target, stats.fragment_y_target, stats.fragment_other_target
    ));
    out.push_str(&format!(
        "      \"decoy\": {{\"b\": {}, \"y\": {}, \"other\": {}}}\n",
        stats.fragment_b_decoy, stats.fragment_y_decoy, stats.fragment_other_decoy
    ));
    out.push_str("    },\n");
    out.push_str("    \"charge_counts\": {\n");
    out.push_str(&format!(
        "      \"precursor\": {},\n",
        json_map_by_class(
            &stats.precursor_charge_counts_target,
            &stats.precursor_charge_counts_decoy
        )
    ));
    out.push_str(&format!(
        "      \"transition\": {}\n",
        json_map_by_class(
            &stats.transition_charge_counts_target,
            &stats.transition_charge_counts_decoy
        )
    ));
    out.push_str("    }\n");
    out.push_str("  }\n");
    out.push_str("}\n");
    out
}
