// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Protein identifications as Arrow tables and Parquet files.
//!
//! Ports `FORMAT/ProteinIdentificationArrowIO.h` and
//! `ProteinIdentificationArrowIO.cpp`: protein hits, protein groups and
//! search parameters, each its own table (`proteins.parquet`,
//! `protein_groups.parquet`, `search_params.parquet`) with the columns of
//! [`ProteinSchema`](crate::format::arrow_schema_registry::ProteinSchema),
//! [`ProteinGroupSchema`](crate::format::arrow_schema_registry::ProteinGroupSchema) and
//! [`SearchParamsSchema`](crate::format::arrow_schema_registry::SearchParamsSchema), and
//! the run-identifier handling the FeatureMap and ConsensusMap Parquet
//! bundles share. See `docs/ARROW_PARQUET_SUPPORT.md`.
//!
//! **Model mapping.** The source's `ProteinIdentification` keeps several
//! things as metadata that this crate's
//! [`ProteinIdentification`](crate::identification::ProteinIdentification) keeps as
//! fields, and the reverse:
//!
//! * the primary MS run paths are the source's `spectra_data` (and
//!   `spectra_data_raw`) metadata; they are
//!   [`primary_ms_run_paths`](crate::identification::ProteinIdentification::primary_ms_run_paths)
//!   (and `raw_ms_run_paths`) here. The export writes them into the
//!   `metavalues` column under those names when non-empty, as the source's
//!   metadata would be; an empty list the source sets explicitly cannot be
//!   told from an unset one here.
//! * the inference engine and its version are `InferenceEngine` and
//!   `InferenceEngineVersion` in the search parameters' metadata, as in the
//!   source.
//! * the run date is the text [`date_time`](crate::identification::ProteinIdentification::date_time)
//!   holds, parsed as the source's `DateTime`.
//! * a digestion enzyme is its name; the import resolves it in
//!   [`ProteaseDB`](crate::chemistry::ProteaseDB) as the source does and stores the canonical name.
//!
//! **Warnings.** Where the source logs a warning and continues (an unknown
//! enzyme, an unknown modification, a hit whose run has no search
//! parameters), the import functions do the same and return the warning
//! lines.

use crate::chemistry::{ModificationsDB, ProteaseDB, SequenceModification};
use crate::comparison::Tolerance;
use crate::data_structures::datetime::DateTime;
use crate::format::arrow_io_helpers::{
    self as helpers, MetaEntry, ParquetWriteConfig, data_value_text, get_bool_value,
    get_double_value, get_int32_value, get_string_value, is_null, list_array, list_element,
    list_range, meta_entries, metavalues_array, read_meta_values, string_list_array,
    string_list_at, struct_array,
};
use crate::format::arrow_schema_registry::validation::{Mode, validate};
use crate::format::arrow_schema_registry::{ProteinGroupSchema, ProteinSchema, SearchParamsSchema};
use crate::format::parquet_file::ParquetFile;
use crate::identification::{
    EnzymeTermSpecificity, PeakMassType, PeptideIdentification, ProteinGroup, ProteinHit,
    ProteinIdentification, ProteinModification, TargetDecoyType,
};
use crate::kernel::DataArray;
use crate::metadata::MetaInfo;
use crate::{Error, Result};
use arrow_array::cast::AsArray;
use arrow_array::types::{
    Float64Type, Int32Type, Int64Type, TimestampMicrosecondType, TimestampMillisecondType,
    TimestampNanosecondType, TimestampSecondType,
};
use arrow_array::{
    Array, ArrayRef, BooleanArray, Float64Array, Int32Array, Int64Array, RecordBatch, StringArray,
    TimestampSecondArray,
};
use arrow_schema::{DataType, SchemaRef, TimeUnit};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;
use std::sync::Arc;

/// Import and export of protein identifications as Arrow tables and Parquet
/// files.
///
/// Ports `OpenMS::ProteinIdentificationArrowIO`, a class of static members;
/// each is an associated function here. The source's functions return
/// `nullptr` or `false` and log; these return a [`Result`].
#[derive(Clone, Copy, Debug)]
pub struct ProteinIdentificationArrowIO;

/// Metadata keys of a protein hit that have their own columns.
const EXCLUDED_HIT_KEYS: [&str; 2] = ["Description", "target_decoy"];

/// Metadata keys of a protein identification the import sets from their own
/// columns.
const EXCLUDED_RUN_KEYS: [&str; 4] = [
    "InferenceEngine",
    "InferenceEngineVersion",
    "spectra_data",
    "spectra_data_raw",
];

/// `Constants::UserParam::MODIFICATION_DEFINITIONS`.
const MODIFICATION_DEFINITIONS: &str = crate::format::modification_definitions::METADATA_KEY;

fn column_error(what: &str, error: arrow_schema::ArrowError) -> Error {
    Error::InvalidValue(format!(
        "ProteinIdentificationArrowIO: Failed to build the {what} table: {error}"
    ))
}

/// A table of `schema`, validated against it in strict mode as the source
/// validates what it built.
fn checked_table(schema: SchemaRef, columns: Vec<ArrayRef>, what: &str) -> Result<RecordBatch> {
    let table = RecordBatch::try_new(schema.clone(), columns).map_err(|e| column_error(what, e))?;
    let result = validate(&table, &schema, Mode::Strict);
    if !result.valid {
        return Err(Error::InvalidValue(format!(
            "Schema validation failed: {result}"
        )));
    }
    Ok(table)
}

fn optional_text(text: &str) -> Option<&str> {
    (!text.is_empty()).then_some(text)
}

/// `ProteinIdentification::hasInferenceEngineAsSearchEngine`.
fn inference_engine_is_search_engine(run: &ProteinIdentification) -> bool {
    let engine = run.search_engine.as_str();
    matches!(
        engine,
        "Fido" | "BayesianProteinInference" | "Epifany" | "ProteinInference"
    ) || (engine == "Percolator" && !run.indistinguishable_groups.is_empty())
}

/// `ProteinIdentification::getInferenceEngine`.
fn inference_engine(run: &ProteinIdentification) -> String {
    match run.search_parameters.metadata.get("InferenceEngine") {
        Some(value) => data_value_text(value),
        None if inference_engine_is_search_engine(run) => run.search_engine.clone(),
        None => String::new(),
    }
}

/// `ProteinIdentification::getInferenceEngineVersion`.
fn inference_engine_version(run: &ProteinIdentification) -> String {
    match run.search_parameters.metadata.get("InferenceEngineVersion") {
        Some(value) => data_value_text(value),
        None if !inference_engine(run).is_empty() => run.search_engine_version.clone(),
        None => String::new(),
    }
}

/// Days since 1970-01-01 of a proleptic Gregorian date.
fn days_from_civil(year: i64, month: i64, day: i64) -> i64 {
    let year = if month <= 2 { year - 1 } else { year };
    let era = year.div_euclid(400);
    let year_of_era = year - era * 400;
    let month_index = (month + 9) % 12;
    let day_of_year = (153 * month_index + 2) / 5 + day - 1;
    let day_of_era = year_of_era * 365 + year_of_era / 4 - year_of_era / 100 + day_of_year;
    era * 146_097 + day_of_era - 719_468
}

/// The date of a day count since 1970-01-01: year, month, day.
fn civil_from_days(days: i64) -> (i64, i64, i64) {
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let day_of_era = z.rem_euclid(146_097);
    let year_of_era =
        (day_of_era - day_of_era / 1460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let month_index = (5 * day_of_year + 2) / 153;
    let day = day_of_year - (153 * month_index + 2) / 5 + 1;
    let month = if month_index < 10 {
        month_index + 3
    } else {
        month_index - 9
    };
    let year = year_of_era + era * 400 + i64::from(month <= 2);
    (year, month, day)
}

/// The export's `date` value: the run date as seconds since the epoch, or
/// `None` when it is unset, does not parse, or lies before 1970 (the
/// source's `timegm` result below zero, which it writes as null).
fn epoch_seconds(date_time: Option<&str>) -> Option<i64> {
    let parsed = DateTime::parse(date_time?).ok()?;
    let text = parsed.format("yyyy-MM-ddThh:mm:ss").ok()?;
    if text.is_empty() {
        return None;
    }
    let number = |range: std::ops::Range<usize>| text.get(range)?.parse::<i64>().ok();
    let seconds = days_from_civil(number(0..4)?, number(5..7)?, number(8..10)?) * 86_400
        + number(11..13)? * 3600
        + number(14..16)? * 60
        + number(17..19)?;
    (seconds >= 0).then_some(seconds)
}

/// The import's run date: a timestamp value in `unit` as the source's
/// `DateTime` text, or `None` for a negative second count (which the source
/// skips) or a date `DateTime` cannot hold.
fn date_from_timestamp(raw: i64, unit: &TimeUnit) -> Option<String> {
    let seconds = match unit {
        TimeUnit::Second => raw,
        TimeUnit::Millisecond => raw.div_euclid(1_000),
        TimeUnit::Microsecond => raw.div_euclid(1_000_000),
        TimeUnit::Nanosecond => raw.div_euclid(1_000_000_000),
    };
    if seconds < 0 {
        return None;
    }
    let (year, month, day) = civil_from_days(seconds.div_euclid(86_400));
    let rest = seconds.rem_euclid(86_400);
    let mut date = DateTime::default();
    date.set_components(
        u32::try_from(month).ok()?,
        u32::try_from(day).ok()?,
        u32::try_from(year).ok()?,
        u32::try_from(rest / 3600).ok()?,
        u32::try_from(rest % 3600 / 60).ok()?,
        u32::try_from(rest % 60).ok()?,
    )
    .ok()?;
    date.format("yyyy-MM-ddThh:mm:ss").ok()
}

/// The index of the last run with each identifier, as the source's
/// `buildRunIdMap_`.
fn run_index(runs: &[ProteinIdentification]) -> HashMap<String, usize> {
    runs.iter()
        .enumerate()
        .map(|(index, run)| (run.identifier.clone(), index))
        .collect()
}

/// The run with `identifier`, created with only that identifier when there
/// is none.
fn run_for(
    identifier: String,
    runs: &mut Vec<ProteinIdentification>,
    index: &mut HashMap<String, usize>,
) -> (usize, bool) {
    if let Some(&at) = index.get(&identifier) {
        return (at, false);
    }
    let at = runs.len();
    runs.push(ProteinIdentification {
        identifier: identifier.clone(),
        ..ProteinIdentification::default()
    });
    index.insert(identifier, at);
    (at, true)
}

/// A column the import needs, with the source's error when it is missing.
fn required(table: &RecordBatch, name: &str, what: &str) -> Result<ArrayRef> {
    helpers::get_column(table, name, true).ok_or_else(|| {
        Error::MissingInformation(format!(
            "ProteinIdentificationArrowIO: Missing required column '{name}' in {what} table"
        ))
    })
}

fn optional(table: &RecordBatch, name: &str) -> Option<ArrayRef> {
    helpers::get_column(table, name, false)
}

/// Refuse a table whose known columns have the wrong type or nullability
/// (the source's subset-mode validation).
fn check_subset(table: &RecordBatch, schema: &SchemaRef) -> Result<()> {
    let result = validate(table, schema, Mode::Subset);
    if result.valid {
        Ok(())
    } else {
        Err(Error::InvalidValue(format!(
            "Incompatible schema: {result}"
        )))
    }
}

impl ProteinIdentificationArrowIO {
    // ==================== Export ====================

    /// Protein hits as an Arrow table, one row per hit of every run.
    ///
    /// Ports `exportProteinsToArrow(protein_identifications)`. Columns, in
    /// [`ProteinSchema`] order: `accession`, `score`, `rank`, `coverage`
    /// (null when unknown), `sequence` (null when empty), `description`
    /// (the `Description` metadata, null when unset or empty), `is_decoy`
    /// (null when the hit has no target/decoy annotation), `run_identifier`,
    /// `modifications` (position and full modification ID, null when none)
    /// and `metavalues` (every other metadata entry).
    ///
    /// A rank above `i32::MAX` wraps as the source's `static_cast<int32_t>`
    /// does, and the import's cast back restores it.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an invalid `target_decoy` annotation;
    /// [`Error::Unsupported`] for a modification that is a mass tag rather
    /// than a named modification (the source's `ProteinHit` holds
    /// `ResidueModification`s only).
    pub fn export_proteins_to_arrow(runs: &[ProteinIdentification]) -> Result<RecordBatch> {
        let rows = runs.iter().map(|run| run.hits.len()).sum::<usize>();
        let mut accession = Vec::with_capacity(rows);
        let mut score = Vec::with_capacity(rows);
        let mut rank = Vec::with_capacity(rows);
        let mut coverage = Vec::with_capacity(rows);
        let mut sequence = Vec::with_capacity(rows);
        let mut description = Vec::with_capacity(rows);
        let mut is_decoy = Vec::with_capacity(rows);
        let mut run_identifier = Vec::with_capacity(rows);
        let mut modification_lengths = Vec::with_capacity(rows);
        let mut positions = Vec::new();
        let mut names = Vec::new();
        let mut metavalues = Vec::with_capacity(rows);
        for run in runs {
            for hit in &run.hits {
                accession.push(hit.accession.clone());
                score.push(hit.score);
                // The source's static_cast<int32_t>(UInt).
                rank.push(hit.rank as i32);
                coverage.push(hit.coverage.filter(|c| *c >= 0.0));
                sequence.push(optional_text(&hit.sequence).map(str::to_owned));
                description.push(
                    hit.metadata
                        .get("Description")
                        .map(data_value_text)
                        .filter(|text| !text.is_empty()),
                );
                is_decoy.push(match hit.target_decoy_type()? {
                    TargetDecoyType::Unknown => None,
                    kind => Some(kind == TargetDecoyType::Decoy),
                });
                run_identifier.push(run.identifier.clone());
                if hit.modifications.is_empty() {
                    modification_lengths.push(None);
                } else {
                    modification_lengths.push(Some(hit.modifications.len()));
                    for modification in &hit.modifications {
                        let known = modification.modification.known().ok_or_else(|| {
                            Error::Unsupported(
                                "ProteinIdentificationArrowIO: a protein modification that is a \
                                 mass tag has no full modification ID to write"
                                    .into(),
                            )
                        })?;
                        positions.push(modification.position as i32);
                        names.push(known.full_id().to_owned());
                    }
                }
                metavalues.push(Some(meta_entries(&hit.metadata, &EXCLUDED_HIT_KEYS)));
            }
        }
        let schema = ProteinSchema::schema();
        let modifications_type = ProteinSchema::modifications_type();
        let modification_values = struct_array(
            list_element(&modifications_type)?,
            vec![
                Arc::new(Int32Array::from(positions)) as ArrayRef,
                Arc::new(StringArray::from(names)),
            ],
        )?;
        let columns: Vec<ArrayRef> = vec![
            Arc::new(StringArray::from(accession)),
            Arc::new(Float64Array::from(score)),
            Arc::new(Int32Array::from(rank)),
            Arc::new(Float64Array::from(coverage)),
            Arc::new(StringArray::from(sequence)),
            Arc::new(StringArray::from(description)),
            Arc::new(BooleanArray::from(is_decoy)),
            Arc::new(StringArray::from(run_identifier)),
            list_array(
                &modifications_type,
                &modification_lengths,
                modification_values,
            )?,
            metavalues_array(&ProteinSchema::metavalues_type(), &metavalues)?,
        ];
        checked_table(schema, columns, "proteins")
    }

    /// Protein hits to a Parquet file.
    ///
    /// Ports `exportProteinsToParquet(protein_identifications, filename,
    /// config)`: [`export_proteins_to_arrow`](Self::export_proteins_to_arrow)
    /// written with the file-level metadata of every table here (see
    /// [`write_table`](Self::write_table)).
    ///
    /// # Errors
    ///
    /// As the export, and as
    /// [`write_table_to_parquet`](helpers::write_table_to_parquet).
    pub fn export_proteins_to_parquet(
        runs: &[ProteinIdentification],
        filename: impl AsRef<Path>,
        config: &ParquetWriteConfig,
    ) -> Result<()> {
        let table = Self::export_proteins_to_arrow(runs)?;
        Self::write_table(&table, filename.as_ref(), "proteins", config)
    }

    /// Protein groups as an Arrow table: every run's protein groups, then
    /// its indistinguishable-protein groups.
    ///
    /// Ports `exportProteinGroupsToArrow(protein_identifications)`. Columns,
    /// in [`ProteinGroupSchema`] order: `group_type` (`protein_group` or
    /// `indistinguishable`), `probability`, `accessions`, `run_identifier`,
    /// `group_index` (0-based within its run and type), and `float_data`,
    /// `string_data`, `integer_data` (each data array's name and values,
    /// null when the group has none). Float values are widened to `double`
    /// and integers to `int64`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a column exceeds Arrow's 2^31-element
    /// list limit or a run has more than 2^31 groups of one type.
    pub fn export_protein_groups_to_arrow(runs: &[ProteinIdentification]) -> Result<RecordBatch> {
        let mut group_type = Vec::new();
        let mut probability = Vec::new();
        let mut accessions = Vec::new();
        let mut run_identifier = Vec::new();
        let mut group_index = Vec::new();
        let mut float_data = DataColumn::default();
        let mut string_data = DataColumn::default();
        let mut integer_data = DataColumn::default();
        let mut float_values = Vec::new();
        let mut string_values = Vec::new();
        let mut integer_values = Vec::new();
        for run in runs {
            for (kind, groups) in [
                ("protein_group", &run.protein_groups),
                ("indistinguishable", &run.indistinguishable_groups),
            ] {
                for (index, group) in groups.iter().enumerate() {
                    group_type.push(kind);
                    probability.push(group.probability);
                    accessions.push(Some(group.accessions.clone()));
                    run_identifier.push(run.identifier.clone());
                    group_index.push(i32::try_from(index).map_err(|_| {
                        Error::InvalidValue("more than 2^31 protein groups in one run".into())
                    })?);
                    float_data.push(&group.float_data_arrays, |array| {
                        float_values.extend(array.data.iter().map(|v| f64::from(*v)));
                    });
                    string_data.push(&group.string_data_arrays, |array| {
                        string_values.extend(array.data.iter().cloned());
                    });
                    integer_data.push(&group.integer_data_arrays, |array| {
                        integer_values.extend(array.data.iter().map(|v| i64::from(*v)));
                    });
                }
            }
        }
        let columns: Vec<ArrayRef> = vec![
            Arc::new(StringArray::from(group_type)),
            Arc::new(Float64Array::from(probability)),
            string_list_array(
                &ProteinGroupSchema::schema().field(2).data_type().clone(),
                &accessions,
            )?,
            Arc::new(StringArray::from(run_identifier)),
            Arc::new(Int32Array::from(group_index)),
            float_data.finish(
                &ProteinGroupSchema::float_data_type(),
                Arc::new(Float64Array::from(float_values)),
            )?,
            string_data.finish(
                &ProteinGroupSchema::string_data_type(),
                Arc::new(StringArray::from(string_values)),
            )?,
            integer_data.finish(
                &ProteinGroupSchema::integer_data_type(),
                Arc::new(Int64Array::from(integer_values)),
            )?,
        ];
        checked_table(ProteinGroupSchema::schema(), columns, "protein_groups")
    }

    /// Protein groups to a Parquet file.
    ///
    /// Ports `exportProteinGroupsToParquet(protein_identifications,
    /// filename, config)`.
    ///
    /// # Errors
    ///
    /// As the export and the write.
    pub fn export_protein_groups_to_parquet(
        runs: &[ProteinIdentification],
        filename: impl AsRef<Path>,
        config: &ParquetWriteConfig,
    ) -> Result<()> {
        let table = Self::export_protein_groups_to_arrow(runs)?;
        Self::write_table(&table, filename.as_ref(), "protein_groups", config)
    }

    /// Search parameters as an Arrow table, one row per run.
    ///
    /// Ports `exportSearchParamsToArrow(protein_identifications,
    /// definitions_by_run)`. Columns, in [`SearchParamsSchema`] order:
    /// `run_identifier`, `search_engine`, `search_engine_version`,
    /// `inference_engine` and `inference_engine_version` (as the source's
    /// getters resolve them, null when empty), `date` (seconds since the
    /// epoch, null when unset or before 1970), `score_type`,
    /// `higher_score_better`, `significance_threshold`, `db`, `db_version`,
    /// `taxonomy`, `charges` (null when empty), `mass_type` (`MONOISOTOPIC`
    /// or `AVERAGE`), the two tolerances with their ppm flags,
    /// `digestion_enzyme` (null when empty or `unknown_enzyme`),
    /// `enzyme_term_specificity` (`FULL`, `SEMI`, `NONE`, null when unknown),
    /// `missed_cleavages`, the fixed and variable modifications, the primary
    /// MS run paths, the run's `metavalues` and the search parameters'
    /// `sp_metavalues`.
    ///
    /// For a run listed in `definitions_by_run`, its
    /// `modification_definitions` entry is replaced by the given text, written
    /// last with type `string` (see
    /// [`modification_definitions`](crate::format::modification_definitions)).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a list column exceeds Arrow's limit.
    pub fn export_search_params_to_arrow(
        runs: &[ProteinIdentification],
        definitions_by_run: &BTreeMap<String, String>,
    ) -> Result<RecordBatch> {
        let n = runs.len();
        let mut run_identifier = Vec::with_capacity(n);
        let mut search_engine = Vec::with_capacity(n);
        let mut search_engine_version = Vec::with_capacity(n);
        let mut engine = Vec::with_capacity(n);
        let mut engine_version = Vec::with_capacity(n);
        let mut date = Vec::with_capacity(n);
        let mut score_type = Vec::with_capacity(n);
        let mut higher_score_better = Vec::with_capacity(n);
        let mut threshold = Vec::with_capacity(n);
        let mut db = Vec::with_capacity(n);
        let mut db_version = Vec::with_capacity(n);
        let mut taxonomy = Vec::with_capacity(n);
        let mut charges = Vec::with_capacity(n);
        let mut mass_type = Vec::with_capacity(n);
        let mut precursor = Vec::with_capacity(n);
        let mut precursor_ppm = Vec::with_capacity(n);
        let mut fragment = Vec::with_capacity(n);
        let mut fragment_ppm = Vec::with_capacity(n);
        let mut enzyme = Vec::with_capacity(n);
        let mut specificity = Vec::with_capacity(n);
        let mut missed = Vec::with_capacity(n);
        let mut fixed = Vec::with_capacity(n);
        let mut variable = Vec::with_capacity(n);
        let mut run_paths = Vec::with_capacity(n);
        let mut metavalues = Vec::with_capacity(n);
        let mut sp_metavalues = Vec::with_capacity(n);
        for run in runs {
            let sp = &run.search_parameters;
            run_identifier.push(run.identifier.clone());
            search_engine.push(run.search_engine.clone());
            search_engine_version
                .push(optional_text(&run.search_engine_version).map(str::to_owned));
            engine.push(Some(inference_engine(run)).filter(|t| !t.is_empty()));
            engine_version.push(Some(inference_engine_version(run)).filter(|t| !t.is_empty()));
            date.push(epoch_seconds(run.date_time.as_deref()));
            score_type.push(run.score_type.clone());
            higher_score_better.push(run.higher_score_better);
            threshold.push(Some(run.significance_threshold));
            db.push(optional_text(&sp.database).map(str::to_owned));
            db_version.push(optional_text(&sp.database_version).map(str::to_owned));
            taxonomy.push(optional_text(&sp.taxonomy).map(str::to_owned));
            charges.push(optional_text(&sp.charges).map(str::to_owned));
            mass_type.push(match sp.mass_type {
                PeakMassType::Monoisotopic => "MONOISOTOPIC",
                PeakMassType::Average => "AVERAGE",
            });
            let (value, ppm) = tolerance_parts(sp.precursor_tolerance);
            precursor.push(value);
            precursor_ppm.push(ppm);
            let (value, ppm) = tolerance_parts(sp.fragment_tolerance);
            fragment.push(value);
            fragment_ppm.push(ppm);
            enzyme.push(
                Some(sp.digestion_enzyme.clone())
                    .filter(|name| !name.is_empty() && name != "unknown_enzyme"),
            );
            specificity.push(match sp.enzyme_specificity {
                EnzymeTermSpecificity::Full => Some("FULL"),
                EnzymeTermSpecificity::Semi => Some("SEMI"),
                EnzymeTermSpecificity::None => Some("NONE"),
                EnzymeTermSpecificity::Unknown => None,
            });
            // The source's static_cast<int32_t>(UInt).
            missed.push(sp.missed_cleavages as i32);
            fixed.push(Some(sp.fixed_modifications.clone()));
            variable.push(Some(sp.variable_modifications.clone()));
            run_paths.push(Some(run.primary_ms_run_paths.clone()));
            metavalues.push(Some(run_meta_entries(run)));
            let mut entries;
            match definitions_by_run.get(&run.identifier) {
                None => entries = meta_entries(&sp.metadata, &[]),
                Some(definitions) => {
                    entries = meta_entries(&sp.metadata, &[MODIFICATION_DEFINITIONS]);
                    entries.push([
                        MODIFICATION_DEFINITIONS.to_owned(),
                        definitions.clone(),
                        "string".to_owned(),
                    ]);
                }
            }
            sp_metavalues.push(Some(entries));
        }
        let schema = SearchParamsSchema::schema();
        let list_type = |index: usize| schema.field(index).data_type().clone();
        let columns: Vec<ArrayRef> = vec![
            Arc::new(StringArray::from(run_identifier)),
            Arc::new(StringArray::from(search_engine)),
            Arc::new(StringArray::from(search_engine_version)),
            Arc::new(StringArray::from(engine)),
            Arc::new(StringArray::from(engine_version)),
            Arc::new(TimestampSecondArray::from(date)),
            Arc::new(StringArray::from(score_type)),
            Arc::new(BooleanArray::from(higher_score_better)),
            Arc::new(Float64Array::from(threshold)),
            Arc::new(StringArray::from(db)),
            Arc::new(StringArray::from(db_version)),
            Arc::new(StringArray::from(taxonomy)),
            Arc::new(StringArray::from(charges)),
            Arc::new(StringArray::from(mass_type)),
            Arc::new(Float64Array::from(precursor)),
            Arc::new(BooleanArray::from(precursor_ppm)),
            Arc::new(Float64Array::from(fragment)),
            Arc::new(BooleanArray::from(fragment_ppm)),
            Arc::new(StringArray::from(enzyme)),
            Arc::new(StringArray::from(specificity)),
            Arc::new(Int32Array::from(missed)),
            string_list_array(&list_type(21), &fixed)?,
            string_list_array(&list_type(22), &variable)?,
            string_list_array(&list_type(23), &run_paths)?,
            metavalues_array(&SearchParamsSchema::metavalues_type(), &metavalues)?,
            metavalues_array(&SearchParamsSchema::metavalues_type(), &sp_metavalues)?,
        ];
        checked_table(schema, columns, "search_params")
    }

    /// Search parameters to a Parquet file.
    ///
    /// Ports `exportSearchParamsToParquet(protein_identifications, filename,
    /// config, definitions_by_run)`.
    ///
    /// # Errors
    ///
    /// As the export and the write.
    pub fn export_search_params_to_parquet(
        runs: &[ProteinIdentification],
        filename: impl AsRef<Path>,
        config: &ParquetWriteConfig,
        definitions_by_run: &BTreeMap<String, String>,
    ) -> Result<()> {
        let table = Self::export_search_params_to_arrow(runs, definitions_by_run)?;
        Self::write_table(&table, filename.as_ref(), "search_params", config)
    }

    /// Write one of the three tables with the file-level metadata the source
    /// gives them.
    ///
    /// Ports the source's `writeArrowTableToParquet_`: the schema metadata
    /// becomes `creator` = `OpenMS`, `file_type` = `file_type`,
    /// `creation_date` (UTC, `yyyy-MM-ddThh:mm:ssZ`), `uuid` (a fresh
    /// version-4 UUID; the source draws it from `std::random_device`) and
    /// `software_provider` = `OpenMS`; the table is then written as
    /// [`write_table_to_parquet`](helpers::write_table_to_parquet) writes
    /// it, which is what the source's copy of that writer does.
    ///
    /// # Errors
    ///
    /// As [`write_table_to_parquet`](helpers::write_table_to_parquet).
    pub fn write_table(
        table: &RecordBatch,
        filename: &Path,
        file_type: &str,
        config: &ParquetWriteConfig,
    ) -> Result<()> {
        let metadata = HashMap::from([
            ("creator".to_owned(), "OpenMS".to_owned()),
            ("file_type".to_owned(), file_type.to_owned()),
            (
                "creation_date".to_owned(),
                DateTime::now_utc().format("yyyy-MM-ddThh:mm:ssZ")?,
            ),
            ("uuid".to_owned(), helpers::generate_uuid_v4()),
            ("software_provider".to_owned(), "OpenMS".to_owned()),
        ]);
        let schema = Arc::new(
            arrow_schema::Schema::new(table.schema().fields().clone()).with_metadata(metadata),
        );
        let table = RecordBatch::try_new(schema, table.columns().to_vec())
            .map_err(|e| column_error(file_type, e))?;
        helpers::write_table_to_parquet(&table, filename, config)
    }

    // ==================== Import ====================

    /// Append one run per row of a search-parameters table.
    ///
    /// Ports `importSearchParamsFromArrow(table, protein_identifications)`
    /// with `registry` as the modification database the source's
    /// `ModificationDefinitionIO::registerFrom` extends: each row's
    /// `modification_definitions` are registered in it. An empty table adds
    /// nothing. The table is checked in subset mode (unknown columns are
    /// ignored). Each row sets the identifier, search engine and version,
    /// score type, score direction, significance threshold (0 when null),
    /// run date (any timestamp unit; a date before 1970 is skipped), the
    /// primary MS run paths (only when the column exists), the search
    /// parameters, the inference engine and version (only when non-empty)
    /// and the run's metadata (without the keys set from columns).
    ///
    /// An enzyme [`ProteaseDB`] does not know, and an unknown
    /// `enzyme_term_specificity`, leave the defaults and add a warning line.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an incompatible schema or a column of the
    /// wrong type; [`Error::MissingInformation`] for a missing required
    /// column; the registration's errors for malformed modification
    /// definitions. `runs` is unchanged on error.
    pub fn import_search_params_from_arrow_with_registry(
        table: &RecordBatch,
        runs: &mut Vec<ProteinIdentification>,
        registry: &mut ModificationsDB,
    ) -> Result<Vec<String>> {
        let mut warnings = Vec::new();
        if table.num_rows() == 0 {
            return Ok(warnings);
        }
        check_subset(table, &SearchParamsSchema::schema())?;
        let what = "search_params";
        let run_id = required(table, SearchParamsSchema::RUN_IDENTIFIER, what)?;
        let engine = required(table, SearchParamsSchema::SEARCH_ENGINE, what)?;
        let score_type = required(table, SearchParamsSchema::SCORE_TYPE, what)?;
        let better = required(table, SearchParamsSchema::HIGHER_SCORE_BETTER, what)?;
        let mass_type = required(table, SearchParamsSchema::MASS_TYPE, what)?;
        let precursor = required(table, SearchParamsSchema::PRECURSOR_MASS_TOLERANCE, what)?;
        let precursor_ppm = required(
            table,
            SearchParamsSchema::PRECURSOR_MASS_TOLERANCE_PPM,
            what,
        )?;
        let fragment = required(table, SearchParamsSchema::FRAGMENT_MASS_TOLERANCE, what)?;
        let fragment_ppm = required(table, SearchParamsSchema::FRAGMENT_MASS_TOLERANCE_PPM, what)?;
        let missed = required(table, SearchParamsSchema::MISSED_CLEAVAGES, what)?;
        let version = optional(table, SearchParamsSchema::SEARCH_ENGINE_VERSION);
        let inference = optional(table, SearchParamsSchema::INFERENCE_ENGINE);
        let inference_version = optional(table, SearchParamsSchema::INFERENCE_ENGINE_VERSION);
        let date = optional(table, SearchParamsSchema::DATE);
        let threshold = optional(table, SearchParamsSchema::SIGNIFICANCE_THRESHOLD);
        let db = optional(table, SearchParamsSchema::DB);
        let db_version = optional(table, SearchParamsSchema::DB_VERSION);
        let taxonomy = optional(table, SearchParamsSchema::TAXONOMY);
        let charges = optional(table, SearchParamsSchema::CHARGES);
        let enzyme = optional(table, SearchParamsSchema::DIGESTION_ENZYME);
        let specificity = optional(table, SearchParamsSchema::ENZYME_TERM_SPECIFICITY);
        let fixed = optional(table, SearchParamsSchema::FIXED_MODIFICATIONS);
        let variable = optional(table, SearchParamsSchema::VARIABLE_MODIFICATIONS);
        let run_paths = optional(table, SearchParamsSchema::PRIMARY_MS_RUN_PATHS);
        let metavalues = optional(table, SearchParamsSchema::METAVALUES);
        let sp_metavalues = optional(table, SearchParamsSchema::SP_METAVALUES);
        let text = |array: &Option<ArrayRef>, row: usize| get_string_value(array.as_deref(), row);
        let mut staged = registry.clone();
        let mut added = Vec::with_capacity(table.num_rows());
        for row in 0..table.num_rows() {
            let mut run = ProteinIdentification {
                identifier: get_string_value(Some(run_id.as_ref()), row)?,
                search_engine: get_string_value(Some(engine.as_ref()), row)?,
                search_engine_version: text(&version, row)?,
                score_type: get_string_value(Some(score_type.as_ref()), row)?,
                higher_score_better: get_bool_value(Some(better.as_ref()), row, false)?,
                significance_threshold: get_double_value(threshold.as_deref(), row, 0.0)?,
                ..ProteinIdentification::default()
            };
            if let Some(dates) = date.as_deref() {
                if !is_null(Some(dates), row) {
                    run.date_time = timestamp_at(dates, row)?;
                }
            }
            if let Some(paths) = run_paths.as_deref() {
                run.primary_ms_run_paths = string_list_at(Some(paths), row)?;
            }
            let sp = &mut run.search_parameters;
            sp.database = text(&db, row)?;
            sp.database_version = text(&db_version, row)?;
            sp.taxonomy = text(&taxonomy, row)?;
            sp.charges = text(&charges, row)?;
            sp.mass_type = if get_string_value(Some(mass_type.as_ref()), row)? == "AVERAGE" {
                PeakMassType::Average
            } else {
                PeakMassType::Monoisotopic
            };
            sp.precursor_tolerance = tolerance(
                get_double_value(Some(precursor.as_ref()), row, 0.0)?,
                get_bool_value(Some(precursor_ppm.as_ref()), row, false)?,
            );
            sp.fragment_tolerance = tolerance(
                get_double_value(Some(fragment.as_ref()), row, 0.0)?,
                get_bool_value(Some(fragment_ppm.as_ref()), row, false)?,
            );
            // The source's static_cast<UInt>(int32).
            sp.missed_cleavages = get_int32_value(Some(missed.as_ref()), row, 0)? as u32;
            let enzyme_name = text(&enzyme, row)?;
            if !enzyme_name.is_empty() {
                match ProteaseDB::global().get_enzyme(&enzyme_name) {
                    Ok(known) => sp.digestion_enzyme = known.name().to_owned(),
                    Err(_) => warnings.push(format!(
                        "ProteinIdentificationArrowIO: Unknown enzyme '{enzyme_name}'"
                    )),
                }
            }
            match text(&specificity, row)?.as_str() {
                "FULL" => sp.enzyme_specificity = EnzymeTermSpecificity::Full,
                "SEMI" => sp.enzyme_specificity = EnzymeTermSpecificity::Semi,
                "NONE" => sp.enzyme_specificity = EnzymeTermSpecificity::None,
                "" => {}
                other => warnings.push(format!(
                    "ProteinIdentificationArrowIO: Unknown enzyme_term_specificity '{other}'"
                )),
            }
            sp.fixed_modifications = string_list_at(fixed.as_deref(), row)?;
            sp.variable_modifications = string_list_at(variable.as_deref(), row)?;
            read_meta_values(sp_metavalues.as_deref(), row, &mut sp.metadata, &[])?;
            crate::format::modification_definitions::register_search_parameters(sp, &mut staged)?;
            let engine_name = text(&inference, row)?;
            if !engine_name.is_empty() {
                sp.metadata
                    .insert("InferenceEngine".into(), engine_name.into());
                let engine_version = text(&inference_version, row)?;
                if !engine_version.is_empty() {
                    sp.metadata
                        .insert("InferenceEngineVersion".into(), engine_version.into());
                }
            }
            read_meta_values(
                metavalues.as_deref(),
                row,
                &mut run.metadata,
                &EXCLUDED_RUN_KEYS,
            )?;
            added.push(run);
        }
        runs.extend(added);
        *registry = staged;
        Ok(warnings)
    }

    /// [`import_search_params_from_arrow_with_registry`](Self::import_search_params_from_arrow_with_registry)
    /// with a copy of the global modification database, which the rows'
    /// `modification_definitions` extend and which is then dropped.
    ///
    /// # Errors
    ///
    /// As the registry form.
    pub fn import_search_params_from_arrow(
        table: &RecordBatch,
        runs: &mut Vec<ProteinIdentification>,
    ) -> Result<Vec<String>> {
        let mut registry = ModificationsDB::global().clone();
        Self::import_search_params_from_arrow_with_registry(table, runs, &mut registry)
    }

    /// Add protein hits to the runs they name.
    ///
    /// Ports `importProteinsFromArrow(table, protein_identifications)`. Each
    /// row becomes a hit of the last run with its `run_identifier`; a run
    /// that does not exist yet is created with only that identifier (and a
    /// warning line). `rank`, `coverage`, `sequence` and `description` are
    /// set when not null; `is_decoy` becomes the `target_decoy` metadata
    /// (`decoy` or `target`); each modification is resolved by its full ID in
    /// `registry`, and one it does not know is skipped with a warning line;
    /// the `metavalues` are read without `Description` and `target_decoy`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an incompatible schema or a column of the
    /// wrong type; [`Error::MissingInformation`] for a missing `accession`,
    /// `score` or `run_identifier`. `runs` is unchanged on error.
    pub fn import_proteins_from_arrow_with_registry(
        table: &RecordBatch,
        runs: &mut Vec<ProteinIdentification>,
        registry: &ModificationsDB,
    ) -> Result<Vec<String>> {
        let mut warnings = Vec::new();
        if table.num_rows() == 0 {
            return Ok(warnings);
        }
        check_subset(table, &ProteinSchema::schema())?;
        let what = "proteins";
        let accession = required(table, ProteinSchema::ACCESSION, what)?;
        let score = required(table, ProteinSchema::SCORE, what)?;
        let run_id = required(table, ProteinSchema::RUN_IDENTIFIER, what)?;
        let rank = optional(table, ProteinSchema::RANK);
        let coverage = optional(table, ProteinSchema::COVERAGE);
        let sequence = optional(table, ProteinSchema::SEQUENCE);
        let description = optional(table, ProteinSchema::DESCRIPTION);
        let is_decoy = optional(table, ProteinSchema::IS_DECOY);
        let modifications = optional(table, ProteinSchema::MODIFICATIONS);
        let metavalues = optional(table, ProteinSchema::METAVALUES);
        let mut staged = runs.clone();
        let mut index = run_index(&staged);
        for row in 0..table.num_rows() {
            let identifier = get_string_value(Some(run_id.as_ref()), row)?;
            let (at, created) = run_for(identifier.clone(), &mut staged, &mut index);
            if created {
                warnings.push(format!(
                    "ProteinIdentificationArrowIO: run_identifier '{identifier}' not found in \
                     search_params; creating minimal ProteinIdentification"
                ));
            }
            let mut hit = ProteinHit {
                accession: get_string_value(Some(accession.as_ref()), row)?,
                score: get_double_value(Some(score.as_ref()), row, 0.0)?,
                ..ProteinHit::default()
            };
            if !is_null(rank.as_deref(), row) {
                // The source's static_cast<UInt>(int32).
                hit.rank = get_int32_value(rank.as_deref(), row, 0)? as u32;
            }
            if !is_null(coverage.as_deref(), row) {
                hit.coverage = Some(get_double_value(coverage.as_deref(), row, 0.0)?);
            }
            if !is_null(sequence.as_deref(), row) {
                hit.sequence = get_string_value(sequence.as_deref(), row)?;
            }
            if !is_null(description.as_deref(), row) {
                hit.set_description(get_string_value(description.as_deref(), row)?);
            }
            if !is_null(is_decoy.as_deref(), row) {
                let decoy = get_bool_value(is_decoy.as_deref(), row, false)?;
                hit.metadata.insert(
                    "target_decoy".into(),
                    if decoy { "decoy" } else { "target" }.into(),
                );
            }
            if let Some(column) = modifications.as_deref() {
                read_modifications(column, row, registry, &mut hit, &mut warnings)?;
            }
            read_meta_values(
                metavalues.as_deref(),
                row,
                &mut hit.metadata,
                &EXCLUDED_HIT_KEYS,
            )?;
            staged[at].hits.push(hit);
        }
        *runs = staged;
        Ok(warnings)
    }

    /// [`import_proteins_from_arrow_with_registry`](Self::import_proteins_from_arrow_with_registry)
    /// with the global modification database.
    ///
    /// # Errors
    ///
    /// As the registry form.
    pub fn import_proteins_from_arrow(
        table: &RecordBatch,
        runs: &mut Vec<ProteinIdentification>,
    ) -> Result<Vec<String>> {
        Self::import_proteins_from_arrow_with_registry(table, runs, ModificationsDB::global())
    }

    /// Add protein groups to the runs they name.
    ///
    /// Ports `importProteinGroupsFromArrow(table, protein_identifications)`.
    /// Each row becomes a protein group of the last run with its
    /// `run_identifier` (created with only that identifier if missing, with
    /// no warning, as in the source), or an indistinguishable-protein group
    /// when `group_type` is `indistinguishable`. `group_index` is not
    /// consulted: the groups keep table order. Float data are narrowed to
    /// `float` and integer data to `Int` (wrapping, as the source's casts).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an incompatible schema or a column of the
    /// wrong type; [`Error::MissingInformation`] for a missing `group_type`,
    /// `probability`, `accessions` or `run_identifier`. `runs` is unchanged
    /// on error.
    pub fn import_protein_groups_from_arrow(
        table: &RecordBatch,
        runs: &mut Vec<ProteinIdentification>,
    ) -> Result<()> {
        if table.num_rows() == 0 {
            return Ok(());
        }
        check_subset(table, &ProteinGroupSchema::schema())?;
        let what = "protein_groups";
        let group_type = required(table, ProteinGroupSchema::GROUP_TYPE, what)?;
        let probability = required(table, ProteinGroupSchema::PROBABILITY, what)?;
        let accessions = required(table, ProteinGroupSchema::ACCESSIONS, what)?;
        let run_id = required(table, ProteinGroupSchema::RUN_IDENTIFIER, what)?;
        let float_data = optional(table, ProteinGroupSchema::FLOAT_DATA);
        let string_data = optional(table, ProteinGroupSchema::STRING_DATA);
        let integer_data = optional(table, ProteinGroupSchema::INTEGER_DATA);
        let mut staged = runs.clone();
        let mut index = run_index(&staged);
        for row in 0..table.num_rows() {
            let identifier = get_string_value(Some(run_id.as_ref()), row)?;
            let (at, _) = run_for(identifier, &mut staged, &mut index);
            let mut group = ProteinGroup {
                probability: get_double_value(Some(probability.as_ref()), row, 0.0)?,
                accessions: string_list_at(Some(accessions.as_ref()), row)?,
                ..ProteinGroup::default()
            };
            group.float_data_arrays = read_data_arrays(float_data.as_deref(), row, |values, i| {
                let doubles = values
                    .as_primitive_opt::<Float64Type>()
                    .ok_or_else(|| data_type_error("float_data", values.data_type()))?;
                // The source's static_cast<float>.
                Ok(doubles.value(i) as f32)
            })?;
            group.string_data_arrays =
                read_data_arrays(string_data.as_deref(), row, |values, i| {
                    let strings = values
                        .as_string_opt::<i32>()
                        .ok_or_else(|| data_type_error("string_data", values.data_type()))?;
                    Ok(if strings.is_null(i) {
                        String::new()
                    } else {
                        strings.value(i).to_owned()
                    })
                })?;
            group.integer_data_arrays =
                read_data_arrays(integer_data.as_deref(), row, |values, i| {
                    let integers = values
                        .as_primitive_opt::<Int64Type>()
                        .ok_or_else(|| data_type_error("integer_data", values.data_type()))?;
                    // The source's static_cast<Int>(int64_t).
                    Ok(integers.value(i) as i32)
                })?;
            if get_string_value(Some(group_type.as_ref()), row)? == "indistinguishable" {
                staged[at].indistinguishable_groups.push(group);
            } else {
                staged[at].protein_groups.push(group);
            }
        }
        *runs = staged;
        Ok(())
    }

    /// Read a Parquet file into a table as the source's `readParquetTable_`.
    fn read(filename: &Path) -> Result<RecordBatch> {
        ParquetFile::read_table(filename)
    }

    /// Search parameters from a Parquet file.
    ///
    /// Ports `importSearchParamsFromParquet(filename,
    /// protein_identifications)`.
    ///
    /// # Errors
    ///
    /// The read's errors and the import's.
    pub fn import_search_params_from_parquet(
        filename: impl AsRef<Path>,
        runs: &mut Vec<ProteinIdentification>,
    ) -> Result<Vec<String>> {
        Self::import_search_params_from_arrow(&Self::read(filename.as_ref())?, runs)
    }

    /// Protein hits from a Parquet file.
    ///
    /// Ports `importProteinsFromParquet(filename, protein_identifications)`.
    ///
    /// # Errors
    ///
    /// The read's errors and the import's.
    pub fn import_proteins_from_parquet(
        filename: impl AsRef<Path>,
        runs: &mut Vec<ProteinIdentification>,
    ) -> Result<Vec<String>> {
        Self::import_proteins_from_arrow(&Self::read(filename.as_ref())?, runs)
    }

    /// Protein groups from a Parquet file.
    ///
    /// Ports `importProteinGroupsFromParquet(filename,
    /// protein_identifications)`.
    ///
    /// # Errors
    ///
    /// The read's errors and the import's.
    pub fn import_protein_groups_from_parquet(
        filename: impl AsRef<Path>,
        runs: &mut Vec<ProteinIdentification>,
    ) -> Result<()> {
        Self::import_protein_groups_from_arrow(&Self::read(filename.as_ref())?, runs)
    }

    /// All three files, reconstructing the runs.
    ///
    /// Ports `importFromParquet(proteins_filename, protein_groups_filename,
    /// search_params_filename, protein_identifications)` with `registry` as
    /// the modification database: `runs` is cleared, then the search
    /// parameters create the runs (registering their modification
    /// definitions in `registry`), the proteins add hits (resolving
    /// modifications in `registry`) and the groups add groups. Returns the
    /// warning lines.
    ///
    /// # Errors
    ///
    /// Any read's or import's error; `runs` and `registry` are then left as
    /// the failing step found them (the source leaves its output partly
    /// filled as well).
    pub fn import_from_parquet_with_registry(
        proteins: impl AsRef<Path>,
        protein_groups: impl AsRef<Path>,
        search_params: impl AsRef<Path>,
        runs: &mut Vec<ProteinIdentification>,
        registry: &mut ModificationsDB,
    ) -> Result<Vec<String>> {
        runs.clear();
        let mut warnings = Self::import_search_params_from_arrow_with_registry(
            &Self::read(search_params.as_ref())?,
            runs,
            registry,
        )?;
        warnings.extend(Self::import_proteins_from_arrow_with_registry(
            &Self::read(proteins.as_ref())?,
            runs,
            registry,
        )?);
        Self::import_protein_groups_from_arrow(&Self::read(protein_groups.as_ref())?, runs)?;
        Ok(warnings)
    }

    /// [`import_from_parquet_with_registry`](Self::import_from_parquet_with_registry)
    /// with a copy of the global modification database.
    ///
    /// # Errors
    ///
    /// As the registry form.
    pub fn import_from_parquet(
        proteins: impl AsRef<Path>,
        protein_groups: impl AsRef<Path>,
        search_params: impl AsRef<Path>,
        runs: &mut Vec<ProteinIdentification>,
    ) -> Result<Vec<String>> {
        let mut registry = ModificationsDB::global().clone();
        Self::import_from_parquet_with_registry(
            proteins,
            protein_groups,
            search_params,
            runs,
            &mut registry,
        )
    }

    // ==================== Run identifiers ====================

    /// Give every run a fresh identifier, as the idXML reader does on load.
    ///
    /// Ports `synthesizeRunIdentifiers(protein_identifications)`: each run's
    /// identifier becomes `<search engine>_<date>_<unique id>`, with
    /// `unknown` for an empty search engine, `1900-01-01T00:00:00` for a
    /// missing or invalid date and a draw of the process-wide ID generator.
    /// Returns the map from each stored identifier to its new one; when
    /// several runs share a stored identifier, the map keeps the last one's
    /// (the source also logs a warning, returned here as the second value).
    pub fn synthesize_run_identifiers(
        runs: &mut [ProteinIdentification],
    ) -> (BTreeMap<String, String>, Vec<String>) {
        let mut rename = BTreeMap::new();
        let mut warnings = Vec::new();
        for run in runs.iter_mut() {
            let engine = if run.search_engine.is_empty() {
                "unknown".to_owned()
            } else {
                run.search_engine.clone()
            };
            let date = run
                .date_time
                .as_deref()
                .and_then(|text| DateTime::parse(text).ok())
                .filter(|date| date.is_valid())
                .and_then(|date| date.format("yyyy-MM-ddThh:mm:ss").ok())
                .unwrap_or_else(|| "1900-01-01T00:00:00".to_owned());
            let synthesized = format!("{engine}_{date}_{}", helpers::process_unique_id());
            let stored = std::mem::replace(&mut run.identifier, synthesized.clone());
            if rename.insert(stored.clone(), synthesized).is_some() {
                warnings.push(format!(
                    "ProteinIdentificationArrowIO: multiple ProtIDs share stored identifier \
                     '{stored}'; pep_id assignment ambiguity resolved by load order — regenerate \
                     input from a fixed-code build to eliminate."
                ));
            }
        }
        (rename, warnings)
    }

    /// Re-stamp peptide identifications with the new run identifiers.
    ///
    /// Ports `applyRunIdentifierRename(rename, pep_ids)`: an identification
    /// whose identifier is not in `rename` keeps it.
    pub fn apply_run_identifier_rename(
        rename: &BTreeMap<String, String>,
        peptides: &mut [PeptideIdentification],
    ) {
        for peptide in peptides {
            if let Some(new) = rename.get(&peptide.identifier) {
                peptide.identifier = new.clone();
            }
        }
    }

    /// Refuse runs whose identifiers are not unique.
    ///
    /// Ports `checkUniqueIdentifiers(protein_identifications)`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] with the source's message and the first
    /// duplicated identifier.
    pub fn check_unique_identifiers(runs: &[ProteinIdentification]) -> Result<()> {
        let mut seen = BTreeSet::new();
        for run in runs {
            if !seen.insert(run.identifier.as_str()) {
                return Err(Error::InvalidValue(format!(
                    "ProteinIdentification run identifiers are not unique. This can lead to loss \
                     of unique PeptideIdentification assignment. Duplicated Protein-ID is: {}",
                    run.identifier
                )));
            }
        }
        Ok(())
    }
}

/// A tolerance as the source's value and ppm flag.
fn tolerance_parts(tolerance: Tolerance) -> (f64, bool) {
    match tolerance {
        Tolerance::Absolute(value) => (value, false),
        Tolerance::Ppm(value) => (value, true),
    }
}

fn tolerance(value: f64, ppm: bool) -> Tolerance {
    if ppm {
        Tolerance::Ppm(value)
    } else {
        Tolerance::Absolute(value)
    }
}

/// The run's `metavalues` entries: its metadata, and the primary MS run
/// paths under the source's metadata names when set.
fn run_meta_entries(run: &ProteinIdentification) -> Vec<MetaEntry> {
    let mut meta = run.metadata.clone();
    for (key, paths) in [
        ("spectra_data", &run.primary_ms_run_paths),
        ("spectra_data_raw", &run.raw_ms_run_paths),
    ] {
        if !paths.is_empty() && !meta.contains_key(key) {
            meta.insert(key.into(), paths.clone().into());
        }
    }
    meta_entries(&meta, &[])
}

/// The run date at `row` of a timestamp column.
fn timestamp_at(array: &dyn Array, row: usize) -> Result<Option<String>> {
    let (raw, unit) = match array.data_type() {
        DataType::Timestamp(unit, _) => {
            let raw = match unit {
                TimeUnit::Second => array.as_primitive::<TimestampSecondType>().value(row),
                TimeUnit::Millisecond => {
                    array.as_primitive::<TimestampMillisecondType>().value(row)
                }
                TimeUnit::Microsecond => {
                    array.as_primitive::<TimestampMicrosecondType>().value(row)
                }
                TimeUnit::Nanosecond => array.as_primitive::<TimestampNanosecondType>().value(row),
            };
            (raw, *unit)
        }
        other => return Err(data_type_error("date", other)),
    };
    Ok(date_from_timestamp(raw, &unit))
}

fn data_type_error(column: &str, data_type: &DataType) -> Error {
    Error::InvalidValue(format!(
        "ProteinIdentificationArrowIO: column '{column}' has unexpected type {data_type}"
    ))
}

/// The modifications at `row`, resolved by full ID in `registry`.
fn read_modifications(
    column: &dyn Array,
    row: usize,
    registry: &ModificationsDB,
    hit: &mut ProteinHit,
    warnings: &mut Vec<String>,
) -> Result<()> {
    let Some((values, start, end)) = list_range(column, row)? else {
        return Ok(());
    };
    if start == end {
        return Ok(());
    }
    let entries = values
        .as_struct_opt()
        .filter(|s| s.num_columns() >= 2)
        .ok_or_else(|| data_type_error("modifications", column.data_type()))?;
    let positions = entries
        .column(0)
        .as_primitive_opt::<Int32Type>()
        .ok_or_else(|| data_type_error("modifications", column.data_type()))?;
    let names = entries
        .column(1)
        .as_string_opt::<i32>()
        .ok_or_else(|| data_type_error("modifications", column.data_type()))?;
    let mut modifications: Vec<ProteinModification> = Vec::with_capacity(end - start);
    for i in start..end {
        // The source's static_cast<Size>(int32).
        let position = positions.value(i) as usize;
        let name = if names.is_null(i) { "" } else { names.value(i) };
        match registry.get_modification_handle(name, None, None) {
            Ok(known) => modifications.push(ProteinModification {
                position,
                modification: SequenceModification::Known(known),
            }),
            Err(_) => warnings.push(format!(
                "ProteinIdentificationArrowIO: Could not find modification '{name}'"
            )),
        }
    }
    // The source collects them in a std::set<pair<Size, ResidueModification>>.
    modifications.sort_by(|a, b| {
        a.position.cmp(&b.position).then_with(|| {
            let id = |m: &ProteinModification| {
                m.modification
                    .known()
                    .map(|k| k.full_id().to_owned())
                    .unwrap_or_default()
            };
            id(a).cmp(&id(b))
        })
    });
    modifications.dedup_by(|a, b| {
        a.position == b.position
            && a.modification.known().map(|k| k.full_id())
                == b.modification.known().map(|k| k.full_id())
    });
    hit.modifications = modifications;
    Ok(())
}

/// The data arrays at `row` of a `list<struct<name, values: list<T>>>`
/// column, each value read by `value`.
fn read_data_arrays<T>(
    column: Option<&dyn Array>,
    row: usize,
    value: impl Fn(&ArrayRef, usize) -> Result<T>,
) -> Result<Vec<DataArray<T>>> {
    let Some(column) = column else {
        return Ok(Vec::new());
    };
    let Some((values, start, end)) = list_range(column, row)? else {
        return Ok(Vec::new());
    };
    let entries = values
        .as_struct_opt()
        .filter(|s| s.num_columns() >= 2)
        .ok_or_else(|| data_type_error("data arrays", column.data_type()))?;
    let names = entries
        .column(0)
        .as_string_opt::<i32>()
        .ok_or_else(|| data_type_error("data arrays", column.data_type()))?;
    let lists = entries.column(1);
    let mut arrays = Vec::with_capacity(end - start);
    for i in start..end {
        let mut array = DataArray::<T> {
            name: if names.is_null(i) {
                String::new()
            } else {
                names.value(i).to_owned()
            },
            data: Vec::new(),
            metadata: MetaInfo::new(),
            data_processing: Vec::new(),
        };
        if let Some((inner, from, to)) = list_range(lists.as_ref(), i)? {
            array.data.reserve(to - from);
            for j in from..to {
                array.data.push(value(&inner, j)?);
            }
        }
        arrays.push(array);
    }
    Ok(arrays)
}

/// Offsets of a `list<struct<name, values: list<T>>>` export column.
#[derive(Default)]
struct DataColumn {
    /// Arrays per row (`None` for a group with none).
    rows: Vec<Option<usize>>,
    /// Names of every array, in order.
    names: Vec<String>,
    /// Values per array.
    lengths: Vec<Option<usize>>,
}

impl DataColumn {
    fn push<T>(&mut self, arrays: &[DataArray<T>], mut values: impl FnMut(&DataArray<T>)) {
        if arrays.is_empty() {
            self.rows.push(None);
            return;
        }
        self.rows.push(Some(arrays.len()));
        for array in arrays {
            self.names.push(array.name.clone());
            self.lengths.push(Some(array.data.len()));
            values(array);
        }
    }

    fn finish(self, data_type: &DataType, values: ArrayRef) -> Result<ArrayRef> {
        let entry_type = list_element(data_type)?;
        let DataType::Struct(fields) = entry_type else {
            return Err(Error::InvalidValue(format!(
                "{entry_type} is not a struct type"
            )));
        };
        let inner = list_array(fields[1].data_type(), &self.lengths, values)?;
        let entries = struct_array(
            entry_type,
            vec![Arc::new(StringArray::from(self.names)) as ArrayRef, inner],
        )?;
        list_array(data_type, &self.rows, entries)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn civil_dates_round_trip() {
        for days in [-719_468, -1, 0, 1, 10_957, 19_675, 2_932_896] {
            let (y, m, d) = civil_from_days(days);
            assert_eq!(days_from_civil(y, m, d), days);
        }
        assert_eq!(civil_from_days(19_675), (2023, 11, 14));
        assert_eq!(
            epoch_seconds(Some("2023-11-14T22:13:20")),
            Some(1_700_000_000)
        );
        assert_eq!(epoch_seconds(Some("1969-12-31T23:59:59")), None);
        assert_eq!(
            date_from_timestamp(1_700_000_000_000, &TimeUnit::Millisecond).as_deref(),
            Some("2023-11-14T22:13:20")
        );
        assert_eq!(date_from_timestamp(-500, &TimeUnit::Millisecond), None);
    }
}
