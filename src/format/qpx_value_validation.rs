// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Stateful value and primary-key validator for QPX Arrow tables.
//!
//! Ports `QPXValueValidation` (`FORMAT/QPXValueValidation.h` and
//! `QPXValueValidation.cpp`). See `docs/QPX_SUPPORT.md` for the API mapping.
//!
//! Arrow schema validation checks declarations such as field type and
//! nullability. It does not reject an empty string in an identity-bearing
//! key, physical nulls under a non-nullable field, duplicate primary keys, or
//! invalid nested values. This validator supplies those write-time checks for
//! the OpenMS QPX psm, feature and pg views. QPX-permitted unmapped features
//! keep OpenMS's empty-peptidoform representation and are still distinguished
//! by their remaining key fields.
//!
//! A validator retains primary keys between [`QpxValueValidation::validate`](crate::format::qpx_value_validation::QpxValueValidation::validate)
//! calls. Streaming writers can therefore validate consecutive batches with
//! one instance and detect duplicates that cross a batch boundary. Use
//! [`QpxValueValidation::reset`](crate::format::qpx_value_validation::QpxValueValidation::reset) before starting a different output file.
//!
//! # The contract this validator enforces
//!
//! This is the single place the QPX input contract is written down; the tools
//! that offer `-out_qpx` (ProteinQuantifier, ProteomicsLFQ, IsobaricWorkflow,
//! ProSE) point here rather than restating it. A table that violates it is
//! refused outright: the export writes nothing, rather than emitting rows that
//! will not join.
//!
//! Consensus maps produced by `ProteomicsLFQ` and `IsobaricWorkflow` satisfy
//! the contract by construction, and they are the supported producers. A map
//! assembled by a different pipeline may not, and the requirements it can
//! fail are:
//!
//! * **Every identification needs a scan.** The `scan` primary-key component
//!   is derived from the identification's spectrum reference.
//!   Identifications carrying none -- protein-inference output, some ID
//!   converters, transferred (match-between-runs) identifications -- cannot
//!   be exported.
//! * **Every identification needs a resolvable origin file.** `run_file_name`
//!   is the key the psm view is joined to the feature and pg views on.
//!   Protein inference commonly drops the path; restore it on the
//!   `ProteinIdentification` (`primary_ms_run_paths`) before exporting.
//! * **Identifications must be unique per (peptidoform, charge, run, scan).**
//!   IDMapper maps one MS2 identification onto every consensus feature whose
//!   RT/m-z window contains it, so the same identification can end up on two
//!   features. Resolve conflicts first, e.g. with IDConflictResolver.
//! * **Every channel needs a QPX label.** The vocabulary is label-free, the
//!   TMT and iTRAQ plexes, and SILAC 2-/3-plex. Those labels are join keys
//!   against the run sidecar, so OpenMS refuses rather than guesses one. Maps
//!   from FeatureFinderMultiplex carrying `ICPL*`, `Leu3`, `"label N"` or
//!   `no_label` fall outside the vocabulary and cannot be exported.
//!
//! Feature rows are keyed on `(peptidoform, charge, run_file_name, rt,
//! observed_mz)`. `observed_mz` is part of that key for the sake of unmapped
//! features: QPX permits them, OpenMS writes them with an empty peptidoform,
//! and `rt` is narrowed to `float32` on write -- so without the m/z two
//! co-eluting unmapped features of one charge in one run would share a key and
//! the export would be refused. Rows agreeing on all five are genuinely
//! indistinguishable and are still refused.
//!
//! **Representation.** A source `std::shared_ptr<arrow::Table>` is an
//! [`arrow_array::RecordBatch`], one array per column, so the source's
//! "column has no chunks" and multi-chunk concatenation paths have no
//! counterpart: a table concatenated from batches
//! (`arrow_select::concat::concat_batches`) is one logical column here as
//! well. A `RecordBatch` reference cannot be null, so the source's "table is
//! null" refusal cannot occur. The source casts every column to the type its
//! view declares without checking (`std::static_pointer_cast`); the exact
//! schema check runs first in both, and a column that nevertheless has
//! another type is reported here as an error instead of being undefined
//! behaviour.
//!
//! The diagnostics are the source's, message for message, in the order the
//! source emits them.

use crate::format::arrow_io_helpers::qpx_is_canonical_intensity_label;
use crate::format::arrow_schema_registry::validation::{self, Mode};
use crate::format::arrow_schema_registry::{QPXFeatureSchema, QPXPSMSchema, QPXPgSchema};
use arrow_array::types::{Float32Type, Int16Type, Int32Type, Int64Type};
use arrow_array::{
    Array, ArrayRef, GenericListArray, PrimitiveArray, RecordBatch, StringArray, StructArray,
};
use arrow_schema::SchemaRef;
use std::collections::HashSet;
use std::fmt;

/// QPX table whose value contract is checked.
///
/// Ports `QPXValueValidation::View`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum View {
    /// The psm view (`QPXPSMSchema`); the source `PSM`.
    Psm,
    /// The feature view (`QPXFeatureSchema`); the source `FEATURE`.
    Feature,
    /// The protein-group view (`QPXPgSchema`); the source `PROTEIN_GROUP`.
    ProteinGroup,
}

impl View {
    /// The view's exact schema.
    fn schema(self) -> SchemaRef {
        match self {
            View::Psm => QPXPSMSchema::schema(),
            View::Feature => QPXFeatureSchema::schema(),
            View::ProteinGroup => QPXPgSchema::schema(),
        }
    }

    /// The name the diagnostics use: `psm`, `feature` or `pg`.
    fn name(self) -> &'static str {
        match self {
            View::Psm => "psm",
            View::Feature => "feature",
            View::ProteinGroup => "pg",
        }
    }

    /// The view's mandatory opaque identity column, which is also its declared
    /// primary key.
    fn id_column(self) -> &'static str {
        match self {
            View::Psm => QPXPSMSchema::PSM_ID,
            View::Feature => QPXFeatureSchema::FEATURE_ID,
            View::ProteinGroup => QPXPgSchema::PG_ID,
        }
    }
}

/// Result of validating one table or streaming batch.
///
/// Ports `QPXValueValidation::Result`. A default result is valid with no
/// errors, as the source's member initialisers make it. Its `Display` is the
/// source `toString()`: the diagnostics joined by `"; "`, so one log line
/// carries every refusal.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ValidationResult {
    /// Whether writing may proceed; the source `valid`, default `true`.
    pub valid: bool,
    /// One diagnostic per problem, in the order found; the source `errors`.
    pub errors: Vec<String>,
}

impl Default for ValidationResult {
    fn default() -> Self {
        Self {
            valid: true,
            errors: Vec::new(),
        }
    }
}

impl fmt::Display for ValidationResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.errors.join("; "))
    }
}

impl ValidationResult {
    fn add(&mut self, error: String) {
        self.valid = false;
        self.errors.push(error);
    }
}

/// Stateful value and primary-key validator for one QPX view.
///
/// Ports `QPXValueValidation`. The source is move-only (copy construction and
/// assignment are deleted); this type does not implement `Clone` either, so a
/// streaming writer cannot fork its accumulated key state by accident.
#[derive(Debug)]
pub struct QpxValueValidation {
    view: View,
    primary_keys: HashSet<String>,
    pg_run_keys: HashSet<String>,
    /// Identities seen in earlier batches: a streaming writer validates one
    /// batch at a time, so uniqueness has to hold across the whole file.
    ids: HashSet<i64>,
}

// ---------------------------------------------------------------------------
// Key encoding, byte for byte as the source builds it
// ---------------------------------------------------------------------------

fn append_string(key: &mut String, value: &str) {
    key.push_str(&value.len().to_string());
    key.push(':');
    key.push_str(value);
    key.push('|');
}

fn append_nullable_string(key: &mut String, value: Option<&str>) {
    match value {
        None => key.push_str("N|"),
        Some(value) => {
            key.push('S');
            append_string(key, value);
        }
    }
}

fn append_integer<I: fmt::Display>(key: &mut String, value: I) {
    key.push_str(&value.to_string());
    key.push('|');
}

fn append_float(key: &mut String, value: f32) {
    // QPX/Arrow equality does not distinguish +0 from -0; normalise before
    // encoding the bits so the in-memory comparison behaves the same.
    let value = if value == 0.0 { 0.0 } else { value };
    append_integer(key, value.to_bits());
}

/// A nullable float key component; a null is its own value, distinct from
/// any float.
fn append_nullable_float(key: &mut String, values: &PrimitiveArray<Float32Type>, row: usize) {
    if values.is_null(row) {
        key.push_str("N|");
        return;
    }
    key.push_str("F|");
    append_float(key, values.value(row));
}

/// The source's `std::isspace` over every byte in the C locale: an empty
/// string is blank too.
fn is_blank(value: &str) -> bool {
    value
        .bytes()
        .all(|b| matches!(b, b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r'))
}

// ---------------------------------------------------------------------------
// Typed column access
// ---------------------------------------------------------------------------

type ListArray = GenericListArray<i32>;

/// The named column, or an error. The source's `combinedColumn`: a missing
/// column is `missing required column 'x'`; the chunk handling has no
/// counterpart (see the module documentation).
fn column<'a>(
    table: &'a RecordBatch,
    name: &str,
    result: &mut ValidationResult,
) -> Option<&'a ArrayRef> {
    let column = table.column_by_name(name);
    if column.is_none() {
        result.add(format!("missing required column '{name}'"));
    }
    column
}

fn downcast<'a, T: 'static>(
    array: &'a dyn Array,
    name: &str,
    result: &mut ValidationResult,
) -> Option<&'a T> {
    let typed = array.as_any().downcast_ref::<T>();
    if typed.is_none() {
        result.add(format!(
            "column '{name}' has unexpected type {}",
            validation::cpp_type_string(array.data_type())
        ));
    }
    typed
}

fn typed_column<'a, T: 'static>(
    table: &'a RecordBatch,
    name: &str,
    result: &mut ValidationResult,
) -> Option<&'a T> {
    let array = column(table, name, result)?;
    downcast::<T>(array.as_ref(), name, result)
}

fn string_child<'a>(
    parent: &'a StructArray,
    index: usize,
    name: &str,
    result: &mut ValidationResult,
) -> Option<&'a StringArray> {
    downcast::<StringArray>(parent.column(index).as_ref(), name, result)
}

/// Offset range of one list row, as absolute indices into its values.
fn list_range(list: &ListArray, row: usize) -> std::ops::Range<usize> {
    let offsets = list.value_offsets();
    let start = usize::try_from(offsets[row]).unwrap_or(0);
    let end = usize::try_from(offsets[row + 1]).unwrap_or(start);
    start..end.max(start)
}

fn nonempty_string(
    values: &StringArray,
    row: usize,
    name: &str,
    result: &mut ValidationResult,
) -> bool {
    if values.is_null(row) {
        result.add(format!("row {row} has a null '{name}' primary-key value"));
        return false;
    }
    if is_blank(values.value(row)) {
        result.add(format!(
            "row {row} has an empty or blank '{name}' primary-key value"
        ));
        return false;
    }
    true
}

fn required_columns_have_values(
    table: &RecordBatch,
    expected: &SchemaRef,
    result: &mut ValidationResult,
) {
    for field in expected.fields() {
        if field.is_nullable() {
            continue;
        }
        if let Some(column) = table.column_by_name(field.name()) {
            let nulls = column.null_count();
            if nulls > 0 {
                result.add(format!(
                    "required column '{}' contains {nulls} physical null value(s)",
                    field.name()
                ));
            }
        }
    }
}

/// Validate a `list<struct{label, intensities: list<struct{intensity_name,
/// intensity_value}>}>` column; with `scalar_labels`, every entry's label
/// must equal its row's pg label.
fn validate_additional_intensities(
    table: &RecordBatch,
    column_name: &str,
    scalar_labels: Option<&StringArray>,
    result: &mut ValidationResult,
) {
    let Some(outer) = typed_column::<ListArray>(table, column_name, result) else {
        return;
    };
    let Some(entries) = downcast::<StructArray>(outer.values().as_ref(), column_name, result)
    else {
        return;
    };
    let Some(labels) = string_child(entries, 0, column_name, result) else {
        return;
    };
    let Some(intensity_lists) =
        downcast::<ListArray>(entries.column(1).as_ref(), column_name, result)
    else {
        return;
    };
    let Some(intensity_entries) =
        downcast::<StructArray>(intensity_lists.values().as_ref(), column_name, result)
    else {
        return;
    };
    let Some(names) = string_child(intensity_entries, 0, column_name, result) else {
        return;
    };
    let Some(values) = downcast::<PrimitiveArray<Float32Type>>(
        intensity_entries.column(1).as_ref(),
        column_name,
        result,
    ) else {
        return;
    };

    for row in 0..table.num_rows() {
        if outer.is_null(row) {
            continue;
        }
        let mut row_labels: HashSet<&str> = HashSet::new();
        for index in list_range(outer, row) {
            if entries.is_null(index) {
                result.add(format!("row {row} contains a null {column_name} entry"));
                continue;
            }

            if labels.is_null(index) || is_blank(labels.value(index)) {
                result.add(format!("row {row} contains an empty {column_name}[].label"));
            } else {
                let label = labels.value(index);
                if !qpx_is_canonical_intensity_label(label) {
                    result.add(format!(
                        "row {row} has non-canonical {column_name} label '{label}'"
                    ));
                }
                if !row_labels.insert(label) {
                    result.add(format!("row {row} repeats {column_name} label '{label}'"));
                }
                if let Some(scalar) = scalar_labels {
                    if scalar.is_null(row) || scalar.value(row) != label {
                        result.add(format!(
                            "row {row} has {column_name} label '{label}' that differs from its \
                             pg label"
                        ));
                    }
                }
            }

            if intensity_lists.is_null(index) {
                result.add(format!(
                    "row {row} has a null {column_name}[].intensities list"
                ));
                continue;
            }
            let mut row_names: HashSet<&str> = HashSet::new();
            for value_index in list_range(intensity_lists, index) {
                if intensity_entries.is_null(value_index) {
                    result.add(format!(
                        "row {row} contains a null {column_name}[].intensities entry"
                    ));
                    continue;
                }
                if names.is_null(value_index) || is_blank(names.value(value_index)) {
                    result.add(format!(
                        "row {row} contains an empty {column_name}[].intensities[].intensity_name"
                    ));
                } else if !row_names.insert(names.value(value_index)) {
                    result.add(format!(
                        "row {row} repeats {column_name} intensity name '{}'",
                        names.value(value_index)
                    ));
                }
                if values.is_null(value_index) || !values.value(value_index).is_finite() {
                    result.add(format!(
                        "row {row} has a null or non-finite {column_name} intensity value"
                    ));
                }
            }
        }
    }
}

impl QpxValueValidation {
    /// Construct a validator for one QPX view, with no keys accumulated.
    ///
    /// Ports `QPXValueValidation(View)`; `view` selects the schema, the
    /// values and the primary key that are checked.
    pub fn new(view: View) -> Self {
        Self {
            view,
            primary_keys: HashSet::new(),
            pg_run_keys: HashSet::new(),
            ids: HashSet::new(),
        }
    }

    /// The view this validator checks.
    ///
    /// Native accessor; the source keeps the view private.
    pub fn view(&self) -> View {
        self.view
    }

    /// Forget primary keys accumulated from earlier [`validate`](Self::validate)
    /// calls; the source `reset()`. Call it before starting a different output
    /// file.
    pub fn reset(&mut self) {
        self.primary_keys.clear();
        self.pg_run_keys.clear();
        self.ids.clear();
    }

    /// Validate one whole table or the next streaming batch.
    ///
    /// Ports `validate(table)`. Checks, in this order: the exact OpenMS QPX
    /// schema of the view ([`Mode::Strict`], each problem reported as
    /// `schema: ...`, and nothing else is checked when it fails); that the
    /// arrays are internally consistent (the source's `ValidateFull`);
    /// physical nulls in required columns; the uniqueness of the view's
    /// opaque identity column (`psm_id`, `feature_id`, `pg_id`), across every
    /// batch validated so far; and then the view's natural key and nested
    /// values:
    ///
    /// * **psm** -- non-empty, non-blank `peptidoform` and `run_file_name`, a
    ///   non-null `charge`, a non-empty `scan` list without null components,
    ///   and a unique `(peptidoform, charge, run_file_name, scan)`; the scan
    ///   list is sorted before comparison, as QPX canonicalises list-valued
    ///   key fields.
    /// * **feature** -- a non-empty `run_file_name`, non-null `peptidoform`
    ///   (it may be empty: an unmapped feature) and `charge`, finite `rt` and
    ///   `observed_mz` when set, a unique `(peptidoform, charge,
    ///   run_file_name, rt, observed_mz)`; each `intensities` entry non-null
    ///   with a non-blank, canonical label unique within the row and a finite
    ///   intensity; the same for `additional_intensities`, whose intensity
    ///   names must also be non-blank and unique per entry.
    /// * **pg** -- a non-blank `anchor_protein` (required, though not part of
    ///   the key); non-empty `pg_accessions` and `grouped_runs` without blank
    ///   members and without a run repeated in one row; `label` and
    ///   `intensity` set together or both null, the label non-blank and
    ///   canonical and the intensity finite; a unique `(pg_accessions as a
    ///   set, grouped_runs as a set, label)`, the composite `pg_id` is
    ///   derived from; runs disjoint within one `(pg_accessions, label)` so
    ///   no measurement is counted twice; and `additional_intensities` as for
    ///   the feature view, each label equal to the row's pg label.
    ///
    /// A 0-row table is valid: the streaming writers emit one for empty
    /// input.
    ///
    /// The accumulated key state advances only when the complete table is
    /// valid, so a rejected batch can be fixed and submitted again without
    /// first calling [`reset`](Self::reset).
    pub fn validate(&mut self, table: &RecordBatch) -> ValidationResult {
        let mut result = ValidationResult::default();
        let expected = self.view.schema();
        let schema_validation = validation::validate(table, &expected, Mode::Strict);
        if !schema_validation.valid {
            for error in schema_validation.errors {
                result.add(format!("schema: {error}"));
            }
            return result;
        }

        for column in table.columns() {
            if let Err(error) = column.to_data().validate_full() {
                result.add(format!("Arrow table is inconsistent: {error}"));
                return result;
            }
        }
        required_columns_have_values(table, &expected, &mut result);

        // The opaque identity column: QPX declares it non-nullable AND the
        // view's primary key, so uniqueness is checked here rather than left
        // to the reader. This runs alongside the natural-key check below: a
        // natural-key duplicate implies an id duplicate, but the reverse would
        // be a genuine hash collision worth naming separately.
        let mut new_ids: HashSet<i64> = HashSet::new();
        let id_column = self.view.id_column();
        let Some(ids) = typed_column::<PrimitiveArray<Int64Type>>(table, id_column, &mut result)
        else {
            return result;
        };
        for row in 0..table.num_rows() {
            if ids.is_null(row) {
                result.add(format!("row {row} has a null '{id_column}'"));
                continue;
            }
            let id = ids.value(row);
            if self.ids.contains(&id) || !new_ids.insert(id) {
                result.add(format!(
                    "row {row} repeats the QPX {} identity '{id_column}'",
                    self.view.name()
                ));
            }
        }

        let mut new_primary_keys: HashSet<String> = HashSet::new();
        let mut new_pg_run_keys: HashSet<String> = HashSet::new();
        match self.view {
            View::Psm => self.validate_psm(table, &mut new_primary_keys, &mut result),
            View::Feature => self.validate_feature(table, &mut new_primary_keys, &mut result),
            View::ProteinGroup => self.validate_pg(
                table,
                &mut new_primary_keys,
                &mut new_pg_run_keys,
                &mut result,
            ),
        }

        if result.valid {
            self.primary_keys.extend(new_primary_keys);
            self.pg_run_keys.extend(new_pg_run_keys);
            self.ids.extend(new_ids);
        }
        result
    }

    fn validate_psm(
        &self,
        table: &RecordBatch,
        new_primary_keys: &mut HashSet<String>,
        result: &mut ValidationResult,
    ) {
        let peptidoform = typed_column::<StringArray>(table, QPXPSMSchema::PEPTIDOFORM, result);
        let charge = typed_column::<PrimitiveArray<Int16Type>>(table, QPXPSMSchema::CHARGE, result);
        let run = typed_column::<StringArray>(table, QPXPSMSchema::RUN_FILE_NAME, result);
        let scan = typed_column::<ListArray>(table, QPXPSMSchema::SCAN, result);
        let (Some(peptidoform), Some(charge), Some(run), Some(scan)) =
            (peptidoform, charge, run, scan)
        else {
            return;
        };
        let Some(scan_values) = downcast::<PrimitiveArray<Int32Type>>(
            scan.values().as_ref(),
            QPXPSMSchema::SCAN,
            result,
        ) else {
            return;
        };

        for row in 0..table.num_rows() {
            let peptidoform_valid =
                nonempty_string(peptidoform, row, QPXPSMSchema::PEPTIDOFORM, result);
            let run_valid = nonempty_string(run, row, QPXPSMSchema::RUN_FILE_NAME, result);
            let mut charge_valid = true;
            if charge.is_null(row) {
                result.add(format!("row {row} has a null 'charge' primary-key value"));
                charge_valid = false;
            }
            let mut key_valid = peptidoform_valid && run_valid && charge_valid;
            let range = list_range(scan, row);
            if scan.is_null(row) || range.is_empty() {
                result.add(format!("row {row} has an empty 'scan' primary-key value"));
                key_valid = false;
            }

            let mut scans: Vec<i32> = Vec::new();
            if !scan.is_null(row) {
                scans.reserve(range.len());
                for index in range {
                    if scan_values.is_null(index) {
                        result.add(format!("row {row} has a null scan component"));
                        key_valid = false;
                    } else {
                        scans.push(scan_values.value(index));
                    }
                }
            }
            if !key_valid {
                continue;
            }

            // QPX canonicalises list-valued primary-key fields by sorting.
            scans.sort_unstable();
            let mut key = String::new();
            append_string(&mut key, peptidoform.value(row));
            append_integer(&mut key, charge.value(row));
            append_string(&mut key, run.value(row));
            append_integer(&mut key, scans.len());
            for value in &scans {
                append_integer(&mut key, value);
            }

            if self.primary_keys.contains(&key) || !new_primary_keys.insert(key) {
                result.add(format!(
                    "row {row} repeats the QPX psm primary key \
                     (peptidoform, charge, run_file_name, scan)"
                ));
            }
        }
    }

    fn validate_feature(
        &self,
        table: &RecordBatch,
        new_primary_keys: &mut HashSet<String>,
        result: &mut ValidationResult,
    ) {
        let peptidoform = typed_column::<StringArray>(table, QPXFeatureSchema::PEPTIDOFORM, result);
        let charge =
            typed_column::<PrimitiveArray<Int16Type>>(table, QPXFeatureSchema::CHARGE, result);
        let run = typed_column::<StringArray>(table, QPXFeatureSchema::RUN_FILE_NAME, result);
        let rt = typed_column::<PrimitiveArray<Float32Type>>(table, QPXFeatureSchema::RT, result);
        let observed_mz = typed_column::<PrimitiveArray<Float32Type>>(
            table,
            QPXFeatureSchema::OBSERVED_MZ,
            result,
        );
        let intensities = typed_column::<ListArray>(table, QPXFeatureSchema::INTENSITIES, result);
        let (
            Some(peptidoform),
            Some(charge),
            Some(run),
            Some(rt),
            Some(observed_mz),
            Some(intensities),
        ) = (peptidoform, charge, run, rt, observed_mz, intensities)
        else {
            return;
        };
        let name = QPXFeatureSchema::INTENSITIES;
        let Some(intensity_values) =
            downcast::<StructArray>(intensities.values().as_ref(), name, result)
        else {
            return;
        };
        let Some(labels) = string_child(intensity_values, 0, name, result) else {
            return;
        };
        let Some(values) = downcast::<PrimitiveArray<Float32Type>>(
            intensity_values.column(1).as_ref(),
            name,
            result,
        ) else {
            return;
        };

        for row in 0..table.num_rows() {
            // QPX permits unmapped features; OpenMS writes their unknown
            // peptidoform as an empty string. The run identity, unlike the
            // optional mapping, must always be present. 'observed_mz' is in
            // the key for those unmapped rows (see the module documentation):
            // do not drop it back out.
            let strings_valid = nonempty_string(run, row, QPXFeatureSchema::RUN_FILE_NAME, result);
            let peptidoform_valid = !peptidoform.is_null(row);
            if !peptidoform_valid {
                result.add(format!(
                    "row {row} has a null 'peptidoform' primary-key value"
                ));
            }
            let charge_valid = !charge.is_null(row);
            if !charge_valid {
                result.add(format!("row {row} has a null 'charge' primary-key value"));
            }
            let mut rt_valid = true;
            if !rt.is_null(row) && !rt.value(row).is_finite() {
                result.add(format!("row {row} has a non-finite 'rt' primary-key value"));
                rt_valid = false;
            }
            let mut observed_mz_valid = true;
            if !observed_mz.is_null(row) && !observed_mz.value(row).is_finite() {
                result.add(format!(
                    "row {row} has a non-finite 'observed_mz' primary-key value"
                ));
                observed_mz_valid = false;
            }
            if strings_valid && peptidoform_valid && charge_valid && rt_valid && observed_mz_valid {
                let mut key = String::new();
                append_string(&mut key, peptidoform.value(row));
                append_integer(&mut key, charge.value(row));
                append_string(&mut key, run.value(row));
                append_nullable_float(&mut key, rt, row);
                append_nullable_float(&mut key, observed_mz, row);
                if self.primary_keys.contains(&key) || !new_primary_keys.insert(key) {
                    result.add(format!(
                        "row {row} repeats the QPX feature primary key \
                         (peptidoform, charge, run_file_name, rt, observed_mz)"
                    ));
                }
            }

            if intensities.is_null(row) {
                continue;
            }
            let mut row_labels: HashSet<&str> = HashSet::new();
            for index in list_range(intensities, row) {
                if intensity_values.is_null(index) {
                    result.add(format!("row {row} contains a null intensities entry"));
                    continue;
                }
                if labels.is_null(index) || is_blank(labels.value(index)) {
                    result.add(format!("row {row} contains an empty intensities[].label"));
                    continue;
                }
                let label = labels.value(index);
                if !qpx_is_canonical_intensity_label(label) {
                    result.add(format!(
                        "row {row} has non-canonical intensity label '{label}'"
                    ));
                }
                if !row_labels.insert(label) {
                    result.add(format!("row {row} repeats intensity label '{label}'"));
                }
                if values.is_null(index) || !values.value(index).is_finite() {
                    result.add(format!(
                        "row {row} has a null or non-finite intensity for label '{label}'"
                    ));
                }
            }
        }
        validate_additional_intensities(
            table,
            QPXFeatureSchema::ADDITIONAL_INTENSITIES,
            None,
            result,
        );
    }

    fn validate_pg(
        &self,
        table: &RecordBatch,
        new_primary_keys: &mut HashSet<String>,
        new_pg_run_keys: &mut HashSet<String>,
        result: &mut ValidationResult,
    ) {
        let anchor = typed_column::<StringArray>(table, QPXPgSchema::ANCHOR_PROTEIN, result);
        let pg_accessions = typed_column::<ListArray>(table, QPXPgSchema::PG_ACCESSIONS, result);
        let grouped_runs = typed_column::<ListArray>(table, QPXPgSchema::GROUPED_RUNS, result);
        let label = typed_column::<StringArray>(table, QPXPgSchema::LABEL, result);
        let intensity =
            typed_column::<PrimitiveArray<Float32Type>>(table, QPXPgSchema::INTENSITY, result);
        let (Some(anchor), Some(pg_accessions), Some(grouped_runs), Some(label), Some(intensity)) =
            (anchor, pg_accessions, grouped_runs, label, intensity)
        else {
            return;
        };
        let Some(run_values) = downcast::<StringArray>(
            grouped_runs.values().as_ref(),
            QPXPgSchema::GROUPED_RUNS,
            result,
        ) else {
            return;
        };
        let Some(accession_values) = downcast::<StringArray>(
            pg_accessions.values().as_ref(),
            QPXPgSchema::PG_ACCESSIONS,
            result,
        ) else {
            return;
        };

        for row in 0..table.num_rows() {
            // anchor_protein is non-nullable in the schema, so it is still
            // required -- but it is descriptive, not part of the key. What
            // identifies the group is its membership.
            let mut key_valid = nonempty_string(anchor, row, QPXPgSchema::ANCHOR_PROTEIN, result);

            // The membership, canonicalised the way pg_id keys on it: a set.
            let mut accessions: Vec<&str> = Vec::new();
            let accession_range = list_range(pg_accessions, row);
            if pg_accessions.is_null(row) || accession_range.is_empty() {
                result.add(format!(
                    "row {row} has an empty 'pg_accessions' primary-key value"
                ));
                key_valid = false;
            } else {
                accessions.reserve(accession_range.len());
                for index in accession_range {
                    if accession_values.is_null(index) || is_blank(accession_values.value(index)) {
                        result.add(format!("row {row} has an empty member in pg_accessions"));
                        key_valid = false;
                        continue;
                    }
                    accessions.push(accession_values.value(index));
                }
                accessions.sort_unstable();
                accessions.dedup();
            }

            let mut runs: Vec<&str> = Vec::new();
            let mut row_runs: HashSet<&str> = HashSet::new();
            let run_range = list_range(grouped_runs, row);
            if grouped_runs.is_null(row) || run_range.is_empty() {
                result.add(format!(
                    "row {row} has an empty 'grouped_runs' primary-key value"
                ));
                key_valid = false;
            } else {
                runs.reserve(run_range.len());
                for index in run_range {
                    if run_values.is_null(index) || is_blank(run_values.value(index)) {
                        result.add(format!("row {row} has an empty member in grouped_runs"));
                        key_valid = false;
                        continue;
                    }
                    let run_name = run_values.value(index);
                    if !row_runs.insert(run_name) {
                        result.add(format!(
                            "row {row} repeats run '{run_name}' inside grouped_runs"
                        ));
                        key_valid = false;
                        continue;
                    }
                    runs.push(run_name);
                }
            }

            let has_label = !label.is_null(row);
            let has_intensity = !intensity.is_null(row);
            if has_label != has_intensity {
                result.add(format!(
                    "row {row} must set label and intensity together, or leave both null"
                ));
            }
            let mut label_value: Option<&str> = None;
            if has_label {
                let value = label.value(row);
                label_value = Some(value);
                if is_blank(value) {
                    result.add(format!("row {row} has an empty 'label' primary-key value"));
                    key_valid = false;
                } else if !qpx_is_canonical_intensity_label(value) {
                    result.add(format!("row {row} has non-canonical label '{value}'"));
                }
            }
            if has_intensity && !intensity.value(row).is_finite() {
                result.add(format!("row {row} has a non-finite intensity"));
            }

            runs.sort_unstable();
            if !key_valid {
                continue;
            }

            // Built from the same three values pg_id hashes, so this check and
            // the identity agree on what "the same pg row" means.
            let mut membership = String::new();
            append_integer(&mut membership, accessions.len());
            for accession in &accessions {
                append_string(&mut membership, accession);
            }

            let mut key = membership.clone();
            append_integer(&mut key, runs.len());
            for run_name in &runs {
                append_string(&mut key, run_name);
            }
            append_nullable_string(&mut key, label_value);
            if self.primary_keys.contains(&key) || !new_primary_keys.insert(key) {
                result.add(format!(
                    "row {row} repeats the QPX pg identity (pg_accessions, grouped_runs, label)"
                ));
            }

            // One raw file may contribute to at most one protein quantity of
            // the same group at this label; otherwise the same measurement is
            // double-counted.
            for run_name in &runs {
                let mut run_key = membership.clone();
                append_nullable_string(&mut run_key, label_value);
                append_string(&mut run_key, run_name);
                if self.pg_run_keys.contains(&run_key) || !new_pg_run_keys.insert(run_key) {
                    result.add(format!(
                        "row {row} reuses run '{run_name}' in another pg row with the same \
                         (pg_accessions, label)"
                    ));
                }
            }
        }
        validate_additional_intensities(
            table,
            QPXPgSchema::ADDITIONAL_INTENSITIES,
            Some(label),
            result,
        );
    }
}
