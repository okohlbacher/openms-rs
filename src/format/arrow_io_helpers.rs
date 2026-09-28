// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Helpers for writing, concatenating and reading Arrow tables as Parquet, and
//! the QPX naming helpers the Parquet exporters share.
//!
//! Ports the `ArrowIOHelpers` namespace (`FORMAT/ArrowIOHelpers.h` and
//! `ArrowIOHelpers.cpp`), together with the `ParquetWriteConfig` struct it
//! takes, which the source declares in `FORMAT/MSExperimentArrowExport.h`.
//! See `docs/ARROW_PARQUET_SUPPORT.md` for the full API mapping.
//!
//! **Tables and arrays** are represented as in [`crate::format::parquet_file`]:
//! a source `std::shared_ptr<arrow::Table>` is an [`arrow_array::RecordBatch`]
//! and a possibly-null `std::shared_ptr<arrow::Array>` an
//! `Option<&dyn Array>`. The source casts every array to the type a getter
//! expects without checking (`std::static_pointer_cast`), which is undefined
//! behaviour for any other type; every getter here checks the type and
//! returns [`Error::InvalidValue`] instead. A row outside the array reads as
//! null, which is what the header documents for [`get_string_value`] (the
//! source's implementation does not check the row at all).
//!
//! **Logging.** The source logs its failures through `OPENMS_LOG_ERROR` and
//! returns `false`, `nullptr` or an empty string. The write helpers here
//! return a [`Result`] carrying the same message; the label helpers keep the
//! source's empty-string contract, because an empty label is the documented
//! refusal their callers test for.
//!
//! **Evidence.** Every behaviour documented as "measured" was recorded from
//! the executed Linux x86_64 Release build (core `bc9cc12`, cli `c19e494`,
//! topp `174b576`, Arrow C++ 23.0.1) on ibminode06; `tests/arrow_parquet.rs`
//! transcribes those records.

use crate::data_structures::datetime::DateTime;
use crate::data_structures::list;
use crate::format::arrow_schema_registry::validation::types_equal;
use crate::format::parquet_file::{encode, publish, writer_properties};
use crate::format::percolator_infile::ScanNumberPattern;
use crate::kernel::features::{ColumnHeader, ConsensusMap};
use crate::metadata::{MetaInfo, MetaValue, MetaValueData};
use crate::{Error, Result};
use arrow_array::cast::AsArray;
use arrow_array::types::{Float32Type, Float64Type, Int32Type, Int64Type};
use arrow_array::{Array, ArrayRef, RecordBatch};
use arrow_schema::{DataType, Schema};
use parquet::arrow::arrow_writer::ArrowWriterOptions;
use parquet::basic::{Compression as ParquetCodec, GzipLevel, ZstdLevel};
use parquet::file::properties::EnabledStatistics;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;
use std::sync::{Arc, Mutex, OnceLock, PoisonError};

// ---------------------------------------------------------------------------
// ParquetWriteConfig
// ---------------------------------------------------------------------------

/// Compression algorithm of a Parquet file.
///
/// Ports `ParquetWriteConfig::Compression` (`MSExperimentArrowExport.h`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ParquetCompression {
    /// No compression (fastest write, largest files).
    None,
    /// Fast compression, moderate ratio (widely compatible).
    Snappy,
    /// Good compression, slower (universal compatibility).
    Gzip,
    /// Very fast compression, lower ratio. Written as Parquet's `LZ4_RAW`
    /// codec, as Arrow C++ writes `Compression::LZ4`.
    Lz4,
    /// Best ratio/speed trade-off (recommended, the default).
    #[default]
    Zstd,
}

/// Configuration for Parquet file writing.
///
/// Ports `ParquetWriteConfig` (`MSExperimentArrowExport.h`): compression, row
/// group size and the other Parquet-specific settings. The defaults are the
/// source's.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ParquetWriteConfig {
    /// Compression algorithm (default [`ParquetCompression::Zstd`]).
    pub compression: ParquetCompression,
    /// Compression level; only ZSTD and GZIP use it (default 3). See
    /// [`write_table_to_parquet`] for the levels each accepts.
    pub compression_level: i32,
    /// Maximum rows per row group (default 134,217,728).
    ///
    /// The source documents this as a size in **bytes** ("128MB row groups")
    /// but passes it to `parquet::arrow::WriteTable` as `chunk_size`, which is
    /// a row count, and Parquet C++ then caps it at 1,048,576 rows
    /// ([`crate::format::parquet_file::MAX_ROW_GROUP_LENGTH`]). This port does
    /// what the source does: rows, capped (measured: a 1,100,000-row table
    /// written with the default becomes row groups of 1,048,576 and 51,424
    /// rows).
    pub row_group_size: i64,
    /// Write column statistics (min/max/null count) per row group (default
    /// `true`).
    pub write_statistics: bool,
    /// Data page size in bytes (default 1,048,576).
    pub data_page_size: i64,
}

impl Default for ParquetWriteConfig {
    fn default() -> Self {
        Self {
            compression: ParquetCompression::Zstd,
            compression_level: 3,
            row_group_size: 128 * 1024 * 1024,
            write_statistics: true,
            data_page_size: 1024 * 1024,
        }
    }
}

// ---------------------------------------------------------------------------
// Key-value metadata
// ---------------------------------------------------------------------------

/// Ordered key-value metadata, as Arrow C++'s `arrow::KeyValueMetadata`.
///
/// The source's type keeps its pairs in insertion order and permits a key to
/// occur more than once; [`get`](Self::get) and [`find_key`](Self::find_key)
/// see the first occurrence. [`qpx_file_metadata`] returns one.
///
/// An arrow-rs [`Schema`] stores its metadata as a map, so
/// [`to_schema_metadata`](Self::to_schema_metadata) keeps the first value of
/// a repeated key and drops the others; the source writes every pair to the
/// Parquet footer.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct KeyValueMetadata {
    pairs: Vec<(String, String)>,
}

impl KeyValueMetadata {
    /// Metadata with no pairs.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append one pair, keeping any earlier pair with the same key.
    pub fn append(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.pairs.push((key.into(), value.into()));
    }

    /// Number of pairs, counting repeated keys.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    /// Whether there are no pairs.
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// The pairs, in insertion order.
    pub fn pairs(&self) -> &[(String, String)] {
        &self.pairs
    }

    /// The key at `index`, as Arrow C++'s `key(i)`.
    pub fn key(&self, index: usize) -> Option<&str> {
        self.pairs.get(index).map(|(key, _)| key.as_str())
    }

    /// The value at `index`, as Arrow C++'s `value(i)`.
    pub fn value(&self, index: usize) -> Option<&str> {
        self.pairs.get(index).map(|(_, value)| value.as_str())
    }

    /// Index of the first pair with `key`, as Arrow C++'s `FindKey` (which
    /// returns -1 where this returns `None`).
    pub fn find_key(&self, key: &str) -> Option<usize> {
        self.pairs.iter().position(|(k, _)| k == key)
    }

    /// Value of the first pair with `key`, as Arrow C++'s `Get`.
    pub fn get(&self, key: &str) -> Option<&str> {
        self.find_key(key).and_then(|index| self.value(index))
    }

    /// Whether any pair has `key`, as Arrow C++'s `Contains`.
    pub fn contains(&self, key: &str) -> bool {
        self.find_key(key).is_some()
    }

    /// The pairs as an arrow-rs schema metadata map; a repeated key keeps its
    /// first value.
    pub fn to_schema_metadata(&self) -> HashMap<String, String> {
        let mut map = HashMap::with_capacity(self.pairs.len());
        for (key, value) in &self.pairs {
            map.entry(key.clone()).or_insert_with(|| value.clone());
        }
        map
    }
}

// ---------------------------------------------------------------------------
// Identity and counts
// ---------------------------------------------------------------------------

fn uuid_generator() -> &'static Mutex<crate::concept::unique_id::UniqueIdGenerator> {
    static GENERATOR: OnceLock<Mutex<crate::concept::unique_id::UniqueIdGenerator>> =
        OnceLock::new();
    GENERATOR.get_or_init(|| Mutex::new(crate::concept::unique_id::UniqueIdGenerator::new()))
}

/// A lowercase hyphenated RFC 4122 version-4 UUID string, e.g.
/// `550e8400-e29b-41d4-a716-446655440000`.
///
/// Ports `generateUuidV4()`, which returns `UniqueIdGenerator::getUUID()` of
/// the process-wide generator. This draws from one process-wide
/// [`UniqueIdGenerator`](crate::concept::unique_id::UniqueIdGenerator) seeded
/// on first use as the source's singleton is, so two calls never repeat a
/// stream the way two freshly clock-seeded generators within the same
/// microsecond would. The identifiers are not cryptographic randomness.
pub fn generate_uuid_v4() -> String {
    uuid_generator()
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .get_uuid()
}

/// Number of rows of a table, or 0 when there is none.
///
/// Ports `tableRowCount(table)`; the source's null table is `None`.
pub fn table_row_count(table: Option<&RecordBatch>) -> usize {
    table.map_or(0, RecordBatch::num_rows)
}

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

/// The Parquet codec for `config`, with its level checked as the executed
/// source accepts it.
fn codec(config: &ParquetWriteConfig) -> Result<ParquetCodec> {
    let level = config.compression_level;
    Ok(match config.compression {
        ParquetCompression::None => ParquetCodec::UNCOMPRESSED,
        ParquetCompression::Snappy => ParquetCodec::SNAPPY,
        ParquetCompression::Lz4 => ParquetCodec::LZ4_RAW,
        ParquetCompression::Gzip => {
            // zlib accepts Z_DEFAULT_COMPRESSION (-1, level 6) and 0..=9.
            let level = match level {
                -1 => 6,
                0..=9 => level.unsigned_abs(),
                _ => {
                    return Err(Error::InvalidValue(format!(
                        "GZIP compression level {level} is outside -1..=9"
                    )));
                }
            };
            ParquetCodec::GZIP(GzipLevel::try_new(level).map_err(|error| {
                Error::InvalidValue(format!("GZIP compression level {level}: {error}"))
            })?)
        }
        ParquetCompression::Zstd => {
            // zstd reads 0 as its default level 3 and clamps above 22; its
            // negative "fast" levels have no arrow-rs equivalent, so they
            // become level 1, the fastest arrow-rs offers.
            let level = match level {
                0 => 3,
                i32::MIN..=-1 => 1,
                23..=i32::MAX => 22,
                _ => level,
            };
            ParquetCodec::ZSTD(ZstdLevel::try_new(level).map_err(|error| {
                Error::InvalidValue(format!("ZSTD compression level {level}: {error}"))
            })?)
        }
    })
}

/// Refuse a QPX view table this port cannot validate before it is written.
///
/// Ports the source's `validateQPXValuesForWrite`: a table whose schema
/// metadata has no `qpx_version` passes; one that has it but no `file_type`
/// is refused; a `psm_file`, `feature_file` or `pg_file` table is checked
/// with `QPXValueValidation`; every other `file_type` passes.
fn validate_qpx_values_for_write(table: &RecordBatch, filename: &Path) -> Result<()> {
    let schema = table.schema();
    let metadata = schema.metadata();
    if !metadata.contains_key("qpx_version") {
        return Ok(());
    }
    let Some(file_type) = metadata.get("file_type") else {
        return Err(Error::InvalidValue(format!(
            "ArrowIOHelpers: QPX table for {} has no file_type metadata.",
            filename.display()
        )));
    };
    match file_type.as_str() {
        "psm_file" | "feature_file" | "pg_file" => Err(Error::Unsupported(format!(
            "ArrowIOHelpers: refusing QPX {file_type} table for {}: QPX value validation \
             (FORMAT/QPXValueValidation.h) is not ported",
            filename.display()
        ))),
        _ => Ok(()),
    }
}

/// Write a table to a Parquet file.
///
/// Ports `writeTableToParquet(table, filename, config)`. The file stores the
/// Arrow schema (the source's `store_schema()`), so schema metadata, a
/// `large_string` column and the other Arrow-only types survive a round
/// trip; lists are written in the Parquet-compliant layout, element field
/// `element` (measured). `config` selects:
///
/// * the codec; LZ4 is written as Parquet's `LZ4_RAW`, as Arrow C++ does;
/// * the level, used for ZSTD and GZIP only. GZIP takes -1 (zlib's default,
///   6) and 0..=9 and refuses anything else (measured: -1, 0, 1 and 9
///   succeed, 10 fails). ZSTD accepts every level (measured: -1, 0, 1, 22 and
///   23 succeed): 0 means zstd's default 3 and a level above 22 is 22, as in
///   zstd. zstd's negative "fast" levels have no arrow-rs counterpart and are
///   written at level 1, so such a file decodes to the same table but its
///   compressed bytes differ from the source's;
/// * the row group size, in rows and capped at 1,048,576 (see
///   [`ParquetWriteConfig::row_group_size`]); a non-positive size refuses a
///   table with rows, as the source does;
/// * column statistics on or off;
/// * the data page size in bytes (0 is accepted, as measured).
///
/// A table whose schema metadata marks it as a QPX `psm_file`,
/// `feature_file` or `pg_file` view is refused: the source validates its
/// values with `QPXValueValidation`, which is not ported yet. A QPX table
/// without `file_type` is refused as in the source; any other QPX
/// `file_type` is written.
///
/// Native differences: the source opens (and truncates) the destination
/// before it encodes and removes the fragment when a later step fails, which
/// also destroys a file that existed there before; this encodes in memory,
/// writes a sibling temporary file and renames it into place, so a failure
/// leaves the destination as it was. The source's null table cannot be
/// passed.
///
/// # Errors
///
/// [`Error::InvalidValue`] for an invalid level, row group size or data page
/// size, a QPX table without `file_type`, or a Parquet encoding failure;
/// [`Error::Unsupported`] for a QPX view table; [`Error::Io`] when the file
/// cannot be written.
pub fn write_table_to_parquet(
    table: &RecordBatch,
    filename: impl AsRef<Path>,
    config: &ParquetWriteConfig,
) -> Result<()> {
    let filename = filename.as_ref();
    validate_qpx_values_for_write(table, filename)?;
    let codec = codec(config)?;
    let page_size = usize::try_from(config.data_page_size).map_err(|_| {
        Error::InvalidValue(format!(
            "Parquet data page size {} is negative",
            config.data_page_size
        ))
    })?;
    let statistics = if config.write_statistics {
        EnabledStatistics::Chunk
    } else {
        EnabledStatistics::None
    };
    let properties = writer_properties(table, config.row_group_size)?
        .set_compression(codec)
        .set_statistics_enabled(statistics)
        .set_data_page_size_limit(page_size)
        .set_coerce_types(true)
        .build();
    let options = ArrowWriterOptions::new().with_properties(properties);
    let bytes = encode(
        table,
        options,
        &format!("ArrowIOHelpers: Failed to write {}", filename.display()),
    )?;
    publish(filename, &bytes)
}

/// Whether two schemas are equal as Arrow C++'s
/// `Schema::Equals(other, check_metadata = false)` decides it: the same
/// field names, types (by [`types_equal`]) and nullability, in order.
fn schemas_equal(a: &Schema, b: &Schema) -> bool {
    a.fields().len() == b.fields().len()
        && a.fields().iter().zip(b.fields().iter()).all(|(x, y)| {
            x.name() == y.name()
                && x.is_nullable() == y.is_nullable()
                && types_equal(x.data_type(), y.data_type())
        })
}

/// Concatenate tables that share a schema and write the result to one
/// Parquet file.
///
/// Ports `concatenateAndWriteToParquet(tables, filename, config, metadata)`.
/// An empty `tables` writes nothing and succeeds. The tables must have equal
/// schemas as Arrow C++'s `ConcatenateTables` requires: the same field names,
/// types and nullability in the same order, metadata aside (measured: a type
/// mismatch and a nullability mismatch both fail). The result carries the
/// first table's schema metadata, or `metadata` when given (measured: the
/// merged file's keys are the first input's, or the replacement's). The file
/// is then written by [`write_table_to_parquet`].
///
/// # Errors
///
/// [`Error::InvalidValue`] when the schemas differ, and whatever
/// [`write_table_to_parquet`] returns.
pub fn concatenate_and_write_to_parquet(
    tables: &[RecordBatch],
    filename: impl AsRef<Path>,
    config: &ParquetWriteConfig,
    metadata: Option<&KeyValueMetadata>,
) -> Result<()> {
    let Some(first) = tables.first() else {
        return Ok(());
    };
    let filename = filename.as_ref();
    let schema = first.schema();
    for (index, table) in tables.iter().enumerate().skip(1) {
        if !schemas_equal(&schema, &table.schema()) {
            return Err(Error::InvalidValue(format!(
                "ArrowIOHelpers: Failed to concatenate tables for {}: schema at index {index} \
                 was different",
                filename.display()
            )));
        }
    }
    let combined = arrow_select::concat::concat_batches(&schema, tables).map_err(|error| {
        Error::InvalidValue(format!(
            "ArrowIOHelpers: Failed to concatenate tables for {}: {error}",
            filename.display()
        ))
    })?;
    let combined = match metadata {
        Some(metadata) => {
            let schema = Arc::new(
                Schema::new(schema.fields().clone()).with_metadata(metadata.to_schema_metadata()),
            );
            combined.with_schema(schema).map_err(|error| {
                Error::InvalidValue(format!(
                    "ArrowIOHelpers: Failed to replace schema metadata for {}: {error}",
                    filename.display()
                ))
            })?
        }
        None => combined,
    };
    write_table_to_parquet(&combined, filename, config)
}

// ---------------------------------------------------------------------------
// QPX file metadata and naming
// ---------------------------------------------------------------------------

/// Primary key and identity composite of the QPX views, transcribed from the
/// records of the executed source (`QPXIdentity.h` composites and the
/// `QPX*Schema` id column names).
fn qpx_identity(file_type: &str) -> Option<(&'static str, &'static str)> {
    use crate::format::arrow_schema_registry::{QPXFeatureSchema, QPXPSMSchema, QPXPgSchema};
    match file_type {
        "feature_file" => Some((
            QPXFeatureSchema::FEATURE_ID,
            "run_file_name,peptidoform,charge,rt,scan,observed_mz",
        )),
        "psm_file" => Some((
            QPXPSMSchema::PSM_ID,
            "run_file_name,scan,peptidoform,charge",
        )),
        "pg_file" => Some((QPXPgSchema::PG_ID, "pg_accessions,grouped_runs,label")),
        _ => None,
    }
}

/// The canonical QPX file-level key-value metadata.
///
/// Ports `qpxFileMetadata(file_type, config, extra)`. The pairs are, in
/// order: `qpx_version` (`1.1`), `file_type`, `creator` (`OpenMS`),
/// `software_provider` (`OpenMS <version>`), `creation_date` (UTC,
/// `yyyy-MM-ddThh:mm:ssZ`), `compression_format` (`none`, `snappy`, `gzip`
/// or `zstd`) and `uuid` (a fresh [`generate_uuid_v4`]); for `psm_file`,
/// `feature_file` and `pg_file` also `primary_key` and
/// `identity_composite`; then every `extra` pair in key order. An `extra` key
/// that repeats a canonical one is appended as a second pair, as in the
/// source (measured: an extra `uuid` gives ten pairs, the first `uuid` being
/// the generated one).
///
/// Build this once per output file: each call mints a new `uuid` and
/// `creation_date`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `config` selects LZ4, which QPX has no
/// `compression_format` token for (the source returns `nullptr`).
pub fn qpx_file_metadata(
    file_type: &str,
    config: &ParquetWriteConfig,
    extra: &BTreeMap<String, String>,
) -> Result<KeyValueMetadata> {
    let compression = match config.compression {
        ParquetCompression::None => "none",
        ParquetCompression::Snappy => "snappy",
        ParquetCompression::Gzip => "gzip",
        ParquetCompression::Zstd => "zstd",
        ParquetCompression::Lz4 => {
            return Err(Error::InvalidValue(
                "ArrowIOHelpers: QPX does not define a compression_format token for LZ4; use \
                 ZSTD (default), SNAPPY, GZIP or NONE for QPX output."
                    .into(),
            ));
        }
    };
    let mut metadata = KeyValueMetadata::new();
    metadata.append("qpx_version", "1.1");
    metadata.append("file_type", file_type);
    metadata.append("creator", "OpenMS");
    metadata.append(
        "software_provider",
        format!("OpenMS {}", crate::CORE_SDK_VERSION),
    );
    metadata.append(
        "creation_date",
        DateTime::now_utc().format("yyyy-MM-ddThh:mm:ssZ")?,
    );
    metadata.append("compression_format", compression);
    metadata.append("uuid", generate_uuid_v4());
    if let Some((primary_key, composite)) = qpx_identity(file_type) {
        metadata.append("primary_key", primary_key);
        metadata.append("identity_composite", composite);
    }
    for (key, value) in extra {
        metadata.append(key.clone(), value.clone());
    }
    Ok(metadata)
}

/// The prefixes `SpectrumNativeIDParser::isNativeID` recognises.
const NATIVE_ID_PREFIXES: [&str; 10] = [
    "scan=",
    "scanId=",
    "scanID=",
    "controllerType=",
    "function=",
    "sample=",
    "index=",
    "spectrum=",
    "file=",
    "frame=",
];

/// Classify one spectrum native ID as a QPX `scan_format` token.
///
/// Ports `qpxScanFormat(const std::string&)`: `index` for an ID starting with
/// `index=`, `scan` for any other ID `SpectrumNativeIDParser::isNativeID`
/// recognises (a prefix test, case-sensitive and without skipping
/// whitespace), and the empty string otherwise (measured: `Index=4`,
/// ` scan=3` and `merged=2` are unrecognised).
pub fn qpx_scan_format(native_id: &str) -> &'static str {
    if native_id.starts_with("index=") {
        "index"
    } else if NATIVE_ID_PREFIXES
        .iter()
        .any(|prefix| native_id.starts_with(prefix))
    {
        "scan"
    } else {
        ""
    }
}

/// One QPX `scan_format` token for a set of native IDs, and the warning the
/// source logs when they disagree.
///
/// Ports `qpxScanFormat(const std::vector<std::string>&)`. Unrecognised IDs
/// are ignored; the token is empty when none is recognised, or when two
/// conventions are mixed, in which case the second element holds the
/// source's `OPENMS_LOG_WARN` line.
pub fn qpx_scan_format_of_all<S: AsRef<str>>(native_ids: &[S]) -> (&'static str, Option<String>) {
    let mut resolved = "";
    for id in native_ids {
        let format = qpx_scan_format(id.as_ref());
        if format.is_empty() {
            continue;
        }
        if resolved.is_empty() {
            resolved = format;
        } else if resolved != format {
            return (
                "",
                Some(format!(
                    "ArrowIOHelpers: spectrum native IDs mix '{resolved}' and '{format}' \
                     conventions; omitting QPX scan_format."
                )),
            );
        }
    }
    (resolved, None)
}

/// Reduce an MS run path to the QPX `run_file_name` form: its stem.
///
/// Ports `qpxRunFileName(ms_run_path)`, which is `File::stemName`: both `/`
/// and `\` separate directories and known compound extensions are removed
/// whole (measured: `/data/proj/S1_Frontal_1.mzML` gives `S1_Frontal_1`,
/// `file.mzML.gz` gives `file`, `C:\x\y.raw` gives `y`, `.hidden` and `dir/`
/// give the empty string).
pub fn qpx_run_file_name(ms_run_path: &str) -> &str {
    crate::system::file::stem_name(ms_run_path)
}

/// The QPX `scan` column value for one spectrum reference.
///
/// Ports `qpxScanComponents(spectrum_reference)`: the scan number
/// `SpectrumNativeIDParser` extracts from the native ID, as a one-element
/// list, or an empty list when the reference is empty or has no scan number
/// that fits `Int32` (measured: `scan=99999999999` and `file=x` are empty,
/// `merged=2` and `abc12def34` fall back to their last digit run).
pub fn qpx_scan_components(spectrum_reference: &str) -> Vec<i32> {
    if spectrum_reference.is_empty() {
        return Vec::new();
    }
    ScanNumberPattern::for_native_id(spectrum_reference)
        .extract(spectrum_reference)
        .filter(|scan| *scan >= 0)
        .map_or_else(Vec::new, |scan| vec![scan])
}

/// The metadata keys [`qpx_psm_run_file_name`] consults.
///
/// Ports `QPXRunFileNameKeys`, which resolves `reference_file_name` and
/// `run_file_name` to `MetaInfoRegistry` indices once per export. This crate
/// keys metadata by name and has no registry, so the keys are the names and
/// always resolve. The source's keys do not resolve when they are built
/// before either name was registered anywhere in the process, in which case
/// every lookup silently falls back to the resolved run path (measured); this
/// port behaves as the source does after registration.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct QpxRunFileNameKeys {
    /// `reference_file_name`.
    pub reference_file_name: &'static str,
    /// `run_file_name`.
    pub run_file_name: &'static str,
}

impl Default for QpxRunFileNameKeys {
    fn default() -> Self {
        Self {
            reference_file_name: "reference_file_name",
            run_file_name: "run_file_name",
        }
    }
}

impl QpxRunFileNameKeys {
    /// The keys, as the source's constructor resolves them.
    pub fn new() -> Self {
        Self::default()
    }
}

/// A metadata value as the source's `DataValue::toString()` renders it: a
/// float with full precision (`StringUtils::toStr`), a list bracketed and
/// separated by `", "`.
fn data_value_text(value: &MetaValue) -> String {
    fn bracketed<I: Iterator<Item = String>>(parts: I) -> String {
        format!("[{}]", parts.collect::<Vec<_>>().join(", "))
    }
    match value.data() {
        MetaValueData::Empty => String::new(),
        MetaValueData::String(text) => text.clone(),
        MetaValueData::Integer(number) => number.to_string(),
        MetaValueData::Float(number) => crate::param::value::format_float(*number, true),
        MetaValueData::StringList(values) => bracketed(values.iter().cloned()),
        MetaValueData::IntegerList(values) => bracketed(values.iter().map(i64::to_string)),
        MetaValueData::FloatList(values) => bracketed(
            values
                .iter()
                .map(|number| crate::param::value::format_float(*number, true)),
        ),
    }
}

/// The QPX `run_file_name` the PSM view writes for one identification and
/// hit.
///
/// Ports `qpxPsmRunFileName(hit, identification, resolved_run_file, keys)`.
/// The first of these that exists supplies the path: the hit's
/// `reference_file_name`, the hit's `run_file_name`, the identification's
/// `reference_file_name` (its `run_file_name` is not consulted). The value is
/// rendered as `DataValue::toString` renders it, so a number or a list
/// becomes its text. When that is empty, `resolved_run_file` is used. The
/// result is reduced by [`qpx_run_file_name`] (measured, e.g. a hit
/// `reference_file_name` of `1.5` gives `1` and one of
/// `["/l/a.mzML", "b"]` gives `a`).
pub fn qpx_psm_run_file_name(
    hit: &MetaInfo,
    identification: &MetaInfo,
    resolved_run_file: &str,
    keys: &QpxRunFileNameKeys,
) -> String {
    let source = hit
        .get(keys.reference_file_name)
        .or_else(|| hit.get(keys.run_file_name))
        .or_else(|| identification.get(keys.reference_file_name));
    let mut run_file = source.map(data_value_text).unwrap_or_default();
    if run_file.is_empty() {
        run_file = resolved_run_file.to_owned();
    }
    qpx_run_file_name(&run_file).to_owned()
}

/// The run names two or more distinct paths collapse onto, each with the
/// source's warning line.
///
/// The detail behind [`qpx_warn_on_run_name_collisions`]: one entry per
/// shared stem, in stem order, holding the line the source writes to
/// `OPENMS_LOG_WARN`. Empty paths are ignored and a path listed twice is one
/// path.
pub fn qpx_run_name_collisions<S: AsRef<str>>(context: &str, ms_run_paths: &[S]) -> Vec<String> {
    let mut stems: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
    for path in ms_run_paths {
        let path = path.as_ref();
        if path.is_empty() {
            continue;
        }
        stems
            .entry(qpx_run_file_name(path))
            .or_default()
            .insert(path);
    }
    stems
        .into_iter()
        .filter(|(_, paths)| paths.len() >= 2)
        .map(|(stem, paths)| {
            format!(
                "{context}: several MS runs share the QPX run name '{stem}' after stripping \
                 path and extension: {}. They cannot be told apart in the exported QPX tables.",
                paths.into_iter().collect::<Vec<_>>().join(", ")
            )
        })
        .collect()
}

/// Whether every distinct MS run path yields a distinct QPX run name.
///
/// Ports `qpxWarnOnRunNameCollisions(context, ms_run_paths)`, which returns
/// this and logs one warning per collision; the warnings are
/// [`qpx_run_name_collisions`].
pub fn qpx_warn_on_run_name_collisions<S: AsRef<str>>(context: &str, ms_run_paths: &[S]) -> bool {
    qpx_run_name_collisions(context, ms_run_paths).is_empty()
}

// ---------------------------------------------------------------------------
// Intensity labels
// ---------------------------------------------------------------------------

/// The isobaric quantitation methods of `IsobaricQuantitationMethod`'s
/// registry (canonical name) with their channel names, in registry order.
///
/// Transcribed from the executed source: the `isobaric` records of the F1a
/// oracle, which enumerate `IsobaricQuantitationMethod::create(type)
/// ->getChannelInformation()` for every method type.
const ISOBARIC_METHODS: &[(&str, &[&str])] = &[
    ("tmt6plex", &["126", "127", "128", "129", "130", "131"]),
    (
        "tmt10plex",
        &[
            "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131",
        ],
    ),
    (
        "tmt11plex",
        &[
            "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C",
        ],
    ),
    (
        "tmt16plex",
        &[
            "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C",
            "132N", "132C", "133N", "133C", "134N",
        ],
    ),
    (
        "tmt18plex",
        &[
            "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C",
            "132N", "132C", "133N", "133C", "134N", "134C", "135N",
        ],
    ),
    (
        "tmt32plex",
        &[
            "126", "127N", "127C", "127D", "128N", "128C", "128ND", "128CD", "129N", "129C",
            "129ND", "129CD", "130N", "130C", "130ND", "130CD", "131N", "131C", "131ND", "131CD",
            "132N", "132C", "132ND", "132CD", "133N", "133C", "133ND", "133CD", "134N", "134ND",
            "134CD", "135ND",
        ],
    ),
    (
        "tmt35plex",
        &[
            "126", "127N", "127C", "127D", "128N", "128C", "128ND", "128CD", "129N", "129C",
            "129ND", "129CD", "130N", "130C", "130ND", "130CD", "131N", "131C", "131ND", "131CD",
            "132N", "132C", "132ND", "132CD", "133N", "133C", "133ND", "133CD", "134N", "134C",
            "134ND", "134CD", "135N", "135ND", "135CD",
        ],
    ),
    ("itraq4plex", &["114", "115", "116", "117"]),
    (
        "itraq8plex",
        &["113", "114", "115", "116", "117", "118", "119", "121"],
    ),
];

/// Whether `name` is a known isobaric method, as
/// `IsobaricQuantitationMethod::methodTypeFromName(name) != UNKNOWN`.
fn is_isobaric_method(name: &str) -> bool {
    ISOBARIC_METHODS.iter().any(|(method, _)| *method == name)
}

/// The non-isobaric part of the canonical vocabulary.
const CANONICAL_LABELS: [&str; 12] = [
    "LFQ",
    "SILAC light",
    "SILAC medium",
    "SILAC heavy",
    "MTRAQ0",
    "MTRAQ4",
    "MTRAQ8",
    "DIMETHYL0",
    "DIMETHYL2",
    "DIMETHYL4",
    "DIMETHYL6",
    "DIMETHYL8",
];

/// Whether a string is in the canonical SDRF/QPX intensity-label vocabulary.
///
/// Ports `qpxIsCanonicalIntensityLabel(label)`: `LFQ`, `SILAC light`,
/// `SILAC medium`, `SILAC heavy`, `MTRAQ0/4/8`, `DIMETHYL0/2/4/6/8`, and
/// `TMT<channel>` or `ITRAQ<channel>` for every channel of every TMT and
/// iTRAQ method OpenMS knows. Case-sensitive.
pub fn qpx_is_canonical_intensity_label(label: &str) -> bool {
    if CANONICAL_LABELS.contains(&label) {
        return true;
    }
    ISOBARIC_METHODS.iter().any(|(method, channels)| {
        let prefix = if method.starts_with("tmt") {
            "TMT"
        } else {
            "ITRAQ"
        };
        label
            .strip_prefix(prefix)
            .is_some_and(|channel| channels.contains(&channel))
    })
}

/// Whether the ASCII-lowercased label names a medium or heavy SILAC channel
/// by its amino-acid label.
fn medium_marker(lower: &str) -> bool {
    lower.contains("arg6") || lower.contains("lys4") || lower.contains("lys6")
}

fn heavy_marker(lower: &str) -> bool {
    lower.contains("arg10") || lower.contains("lys8")
}

/// QPX `intensities[].label` for one ConsensusMap column header.
///
/// Ports `qpxIntensityLabel(column_label, channel_name)`:
///
/// * an isobaric header resolves to `TMT<channel>` or `ITRAQ<channel>`. The
///   method is the part of `column_label` before its last `_`, and must be
///   one of the methods OpenMS knows; the channel is `channel_name`, or when
///   that is empty the part of `column_label` after the last `_`. The
///   channel itself is not checked (measured: `tmt6plex_abc` gives
///   `TMTabc`);
/// * with no `channel_name`, the label-free and plex vocabularies apply,
///   case-insensitively: empty, `label-free` and `lfq` give `LFQ`;
///   `no_label` and `light` give `SILAC light`; `medium` or an `Arg6`,
///   `Lys4` or `Lys6` label without `Arg10`/`Lys8` gives `SILAC medium`;
///   `heavy` or an `Arg10`/`Lys8` label gives `SILAC heavy`; `Dimethyl<n>`
///   gives `DIMETHYL<n>` when that is canonical; and an already canonical
///   label is kept;
/// * anything else gives the empty string, which the source logs as an
///   error and callers must refuse.
///
/// Prefer [`qpx_intensity_labels`] for a whole map: it can tell a two-plex
/// heavy channel from a three-plex medium one.
pub fn qpx_intensity_label(column_label: &str, channel_name: &str) -> String {
    let (method_name, suffix) = match column_label.rfind('_') {
        Some(at) => (&column_label[..at], Some(&column_label[at + 1..])),
        None => ("", None),
    };
    let method_known = is_isobaric_method(method_name);
    let mut channel = channel_name;
    if channel.is_empty() && method_known {
        channel = suffix.unwrap_or("");
    }
    if !channel.is_empty() && method_known {
        if method_name.starts_with("tmt") {
            return format!("TMT{channel}");
        }
        if method_name.starts_with("itraq") {
            return format!("ITRAQ{channel}");
        }
    }
    if !channel_name.is_empty() {
        return String::new();
    }
    if column_label.is_empty() || column_label == "label-free" {
        return "LFQ".into();
    }
    let lower = column_label.to_ascii_lowercase();
    if lower == "lfq" {
        return "LFQ".into();
    }
    if lower == "no_label" || lower == "silac light" || lower == "light" {
        return "SILAC light".into();
    }
    if lower == "silac medium"
        || lower == "medium"
        || (medium_marker(&lower) && !heavy_marker(&lower))
    {
        return "SILAC medium".into();
    }
    if lower == "silac heavy" || lower == "heavy" || heavy_marker(&lower) {
        return "SILAC heavy".into();
    }
    if lower.starts_with("dimethyl") {
        // The first eight bytes are the ASCII letters of "dimethyl", so the
        // cut is on a character boundary.
        let canonical = format!("DIMETHYL{}", column_label.get(8..).unwrap_or(""));
        if qpx_is_canonical_intensity_label(&canonical) {
            return canonical;
        }
    }
    if qpx_is_canonical_intensity_label(column_label) {
        return column_label.to_owned();
    }
    String::new()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SilacRole {
    Unknown,
    Light,
    Medium,
    Heavy,
}

fn silac_role(label: &str) -> SilacRole {
    let lower = label.to_ascii_lowercase();
    if lower.is_empty() || lower == "no_label" || lower == "silac light" || lower == "light" {
        return SilacRole::Light;
    }
    let medium = lower == "silac medium" || lower == "medium" || medium_marker(&lower);
    let heavy = lower == "silac heavy" || lower == "heavy" || heavy_marker(&lower);
    match (medium, heavy) {
        (true, false) => SilacRole::Medium,
        (false, true) => SilacRole::Heavy,
        _ => SilacRole::Unknown,
    }
}

fn has_silac_marker(label: &str) -> bool {
    matches!(silac_role(label), SilacRole::Medium | SilacRole::Heavy)
        || label.to_ascii_lowercase().contains("silac")
}

fn canonical_silac_label(role: SilacRole) -> &'static str {
    match role {
        SilacRole::Light => "SILAC light",
        SilacRole::Medium => "SILAC medium",
        SilacRole::Heavy => "SILAC heavy",
        SilacRole::Unknown => "",
    }
}

/// The header's `channel_name` metadata as text, or empty.
fn channel_name_of(header: &ColumnHeader) -> String {
    header
        .metadata
        .get("channel_name")
        .map(data_value_text)
        .unwrap_or_default()
}

/// Canonical QPX intensity labels for every column of a ConsensusMap.
///
/// Ports `qpxIntensityLabels(cmap)`. The headers are grouped by source file
/// (`ColumnHeader::filename`). A group with no isobaric channel
/// (`channel_name`) but a SILAC marker (a medium or heavy label, or any
/// label containing `silac`) is labelled by channel role: a two-plex must
/// have exactly one light channel and names the other `SILAC heavy`, even an
/// `Arg6` one; a three-plex must have one light, one medium and one heavy
/// channel. A group that is not such a plex has every header mapped to the
/// empty string (the source also logs an error). Every other group is
/// labelled header by header with [`qpx_intensity_label`].
///
/// Returns map index to label; an empty label means the header cannot be
/// represented, and exporters must refuse it.
pub fn qpx_intensity_labels(cmap: &ConsensusMap) -> BTreeMap<u64, String> {
    let mut by_source: BTreeMap<&str, Vec<(u64, &ColumnHeader)>> = BTreeMap::new();
    for (index, header) in &cmap.column_headers {
        by_source
            .entry(header.filename.as_str())
            .or_default()
            .push((*index, header));
    }
    let mut result = BTreeMap::new();
    for headers in by_source.values() {
        let has_isobaric_channel = headers
            .iter()
            .any(|(_, header)| !channel_name_of(header).is_empty());
        let has_marker = headers
            .iter()
            .any(|(_, header)| has_silac_marker(&header.label));
        if !has_isobaric_channel && has_marker {
            let roles: Vec<SilacRole> = headers
                .iter()
                .map(|(_, header)| silac_role(&header.label))
                .collect();
            let count = |role| roles.iter().filter(|r| **r == role).count();
            let supported = headers.len() == 2 || headers.len() == 3;
            let known = roles.iter().all(|role| *role != SilacRole::Unknown);
            let resolved: Option<Vec<SilacRole>> = if !supported || !known {
                None
            } else if headers.len() == 2 {
                (count(SilacRole::Light) == 1).then(|| {
                    roles
                        .iter()
                        .map(|role| {
                            if *role == SilacRole::Light {
                                SilacRole::Light
                            } else {
                                SilacRole::Heavy
                            }
                        })
                        .collect()
                })
            } else {
                (count(SilacRole::Light) == 1
                    && count(SilacRole::Medium) == 1
                    && count(SilacRole::Heavy) == 1)
                    .then(|| roles.clone())
            };
            match resolved {
                Some(roles) => {
                    for ((index, _), role) in headers.iter().zip(roles) {
                        result.insert(*index, canonical_silac_label(role).to_owned());
                    }
                }
                None => {
                    for (index, _) in headers {
                        result.insert(*index, String::new());
                    }
                }
            }
            continue;
        }
        for (index, header) in headers {
            result.insert(
                *index,
                qpx_intensity_label(&header.label, &channel_name_of(header)),
            );
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Read helpers
// ---------------------------------------------------------------------------

/// A named column of a table, or `None` when there is none.
///
/// Ports `getColumn(table, name, required)`. The source combines a
/// multi-chunk column into one array; a [`RecordBatch`] column is one array
/// already. `required` only decides whether the source logs the missing
/// column as an error; it does not change the result. When several columns
/// share the name the first is returned, where the source's
/// `GetColumnByName` returns none.
pub fn get_column(table: &RecordBatch, name: &str, required: bool) -> Option<ArrayRef> {
    let _ = required;
    crate::format::parquet_file::ParquetFile::get_optional_column(table, name)
}

/// The array when `row` holds a value: `None` for a missing array, a null
/// value or a row outside the array.
fn present(array: Option<&dyn Array>, row: usize) -> Option<&dyn Array> {
    array.filter(|array| row < array.len() && !array.is_null(row))
}

fn type_mismatch(getter: &str, expected: &str, actual: &DataType) -> Error {
    Error::InvalidValue(format!(
        "ArrowIOHelpers::{getter} expects a {expected} array, got {actual}"
    ))
}

/// Whether `array` has no value at `row`: a missing array, a null value, or
/// a row outside the array.
///
/// Ports `isNull(array, row)` ("or unset").
pub fn is_null(array: Option<&dyn Array>, row: usize) -> bool {
    present(array, row).is_none()
}

/// The string at `row`, or the empty string when there is none.
///
/// Ports `getStringValue(array, row)`, which documents the empty string for a
/// null value and a row outside the array.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the array is not `utf8` (the source's
/// unchecked cast).
pub fn get_string_value(array: Option<&dyn Array>, row: usize) -> Result<String> {
    let Some(array) = present(array, row) else {
        return Ok(String::new());
    };
    match array.data_type() {
        DataType::Utf8 => Ok(array.as_string::<i32>().value(row).to_owned()),
        other => Err(type_mismatch("getStringValue", "utf8", other)),
    }
}

/// The `double` at `row`, or `default_value` when there is none.
///
/// Ports `getDoubleValue(array, row, default_val)`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the array is not `double`.
pub fn get_double_value(array: Option<&dyn Array>, row: usize, default_value: f64) -> Result<f64> {
    let Some(array) = present(array, row) else {
        return Ok(default_value);
    };
    match array.data_type() {
        DataType::Float64 => Ok(array.as_primitive::<Float64Type>().value(row)),
        other => Err(type_mismatch("getDoubleValue", "double", other)),
    }
}

/// The `float` at `row`, or `default_value` when there is none.
///
/// Ports `getFloatValue(array, row, default_val)`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the array is not `float`.
pub fn get_float_value(array: Option<&dyn Array>, row: usize, default_value: f32) -> Result<f32> {
    let Some(array) = present(array, row) else {
        return Ok(default_value);
    };
    match array.data_type() {
        DataType::Float32 => Ok(array.as_primitive::<Float32Type>().value(row)),
        other => Err(type_mismatch("getFloatValue", "float", other)),
    }
}

/// The `int32` at `row`, or `default_value` when there is none.
///
/// Ports `getInt32Value(array, row, default_val)`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the array is not `int32`.
pub fn get_int32_value(array: Option<&dyn Array>, row: usize, default_value: i32) -> Result<i32> {
    let Some(array) = present(array, row) else {
        return Ok(default_value);
    };
    match array.data_type() {
        DataType::Int32 => Ok(array.as_primitive::<Int32Type>().value(row)),
        other => Err(type_mismatch("getInt32Value", "int32", other)),
    }
}

/// The `int64` at `row`, or `default_value` when there is none.
///
/// Ports `getInt64Value(array, row, default_val)`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the array is not `int64`.
pub fn get_int64_value(array: Option<&dyn Array>, row: usize, default_value: i64) -> Result<i64> {
    let Some(array) = present(array, row) else {
        return Ok(default_value);
    };
    match array.data_type() {
        DataType::Int64 => Ok(array.as_primitive::<Int64Type>().value(row)),
        other => Err(type_mismatch("getInt64Value", "int64", other)),
    }
}

/// The `bool` at `row`, or `default_value` when there is none.
///
/// Ports `getBoolValue(array, row, default_val)`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the array is not `bool`.
pub fn get_bool_value(array: Option<&dyn Array>, row: usize, default_value: bool) -> Result<bool> {
    let Some(array) = present(array, row) else {
        return Ok(default_value);
    };
    match array.data_type() {
        DataType::Boolean => Ok(array.as_boolean().value(row)),
        other => Err(type_mismatch("getBoolValue", "bool", other)),
    }
}

/// C locale `isspace`.
fn is_c_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// `std::stoll(text)`: leading whitespace, an optional sign and the longest
/// run of decimal digits, the rest ignored. `None` where `stoll` throws (no
/// digits, or out of range).
fn c_stoll(text: &str) -> Option<i64> {
    let bytes = text.as_bytes();
    let mut i = bytes.iter().take_while(|b| is_c_space(**b)).count();
    let negative = match bytes.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let digits = bytes[i..].iter().take_while(|b| b.is_ascii_digit()).count();
    if digits == 0 {
        return None;
    }
    // Accumulate negatively so that i64::MIN is reachable.
    let mut value: i64 = 0;
    for &digit in &bytes[i..i + digits] {
        value = value
            .checked_mul(10)?
            .checked_sub(i64::from(digit - b'0'))?;
    }
    if negative {
        Some(value)
    } else {
        value.checked_neg()
    }
}

/// Whether `text` starts with `prefix`, ASCII case-insensitively.
fn starts_with_ignore_case(text: &[u8], prefix: &str) -> bool {
    text.len() >= prefix.len() && text[..prefix.len()].eq_ignore_ascii_case(prefix.as_bytes())
}

/// Round `mantissa * 2^exponent` (plus a nonzero tail below it when
/// `sticky`) to a double, to nearest with ties to even. `None` where glibc's
/// `strtod` reports `ERANGE`: overflow, or an inexact result below the
/// smallest normal double (tininess after rounding, as on x86-64).
fn assemble_binary(mantissa: u64, exponent: i64, sticky: bool) -> Option<f64> {
    if mantissa == 0 {
        return Some(0.0);
    }
    let bits = 64 - i64::from(mantissa.leading_zeros());
    let top = exponent.saturating_add(bits - 1);
    if top > 1023 {
        return None;
    }
    let lsb = top.saturating_sub(52).max(-1074);
    let shift = lsb.saturating_sub(exponent);
    let wide = u128::from(mantissa);
    let (mut q, inexact) = if shift <= 0 {
        // At most 53 significant bits remain, so the shifted value fits.
        (wide << u32::try_from(-shift).unwrap_or(0).min(64), sticky)
    } else if shift >= 128 {
        (0, true)
    } else {
        let shift = u32::try_from(shift).unwrap_or(127);
        let q = wide >> shift;
        let rest = wide & ((1u128 << shift) - 1);
        let half = 1u128 << (shift - 1);
        let round_up = rest > half || (rest == half && (sticky || q & 1 == 1));
        (q + u128::from(round_up), rest != 0 || sticky)
    };
    let mut e = lsb;
    if q == 1u128 << 53 {
        q >>= 1;
        e += 1;
    }
    let q = u64::try_from(q).ok()?;
    let tiny = q < 1u64 << 52;
    if inexact && tiny {
        return None;
    }
    if tiny {
        return Some(f64::from_bits(q));
    }
    let biased = e + 52 + 1023;
    if biased >= 2047 {
        return None;
    }
    let biased = u64::try_from(biased).ok()?;
    Some(f64::from_bits((biased << 52) | (q & ((1u64 << 52) - 1))))
}

/// The hexadecimal float after `0x`: `None` when no hex digit follows (the
/// caller then reads the leading `0`), `Some(None)` for `ERANGE`.
fn c_hex_float(bytes: &[u8]) -> Option<Option<f64>> {
    let mut mantissa: u64 = 0;
    let mut exponent: i64 = 0;
    let mut sticky = false;
    let mut any = false;
    let mut dot = false;
    let mut i = 0;
    loop {
        match bytes.get(i) {
            Some(&c) if c.is_ascii_hexdigit() => {
                any = true;
                let digit = u64::from(char::from(c).to_digit(16).unwrap_or(0));
                if mantissa >> 60 == 0 {
                    mantissa = (mantissa << 4) | digit;
                    if dot {
                        exponent -= 4;
                    }
                } else {
                    sticky |= digit != 0;
                    if !dot {
                        exponent = exponent.saturating_add(4);
                    }
                }
            }
            Some(b'.') if !dot => dot = true,
            _ => break,
        }
        i += 1;
    }
    if !any {
        return None;
    }
    if matches!(bytes.get(i), Some(b'p' | b'P')) {
        let mut j = i + 1;
        let negative = match bytes.get(j) {
            Some(b'-') => {
                j += 1;
                true
            }
            Some(b'+') => {
                j += 1;
                false
            }
            _ => false,
        };
        let mut value: i64 = 0;
        let mut digits = 0;
        while let Some(digit) = bytes.get(j).filter(|b| b.is_ascii_digit()) {
            value = value
                .saturating_mul(10)
                .saturating_add(i64::from(digit - b'0'))
                .min(1 << 40);
            digits += 1;
            j += 1;
        }
        if digits > 0 {
            exponent = exponent.saturating_add(if negative { -value } else { value });
        }
    }
    Some(assemble_binary(mantissa, exponent, sticky))
}

/// Base-10^9 digits of `k * 5^1074`, most significant first, as a decimal
/// string: the exact digits of the subnormal double `k * 2^-1074` scaled by
/// 10^1074.
fn subnormal_digits(k: u64) -> String {
    const BASE: u64 = 1_000_000_000;
    let mut limbs: Vec<u64> = Vec::new();
    let mut rest = k;
    while rest > 0 {
        limbs.push(rest % BASE);
        rest /= BASE;
    }
    let mut multiply = |factor: u64| {
        let mut carry = 0u64;
        for limb in &mut limbs {
            let product = *limb * factor + carry;
            *limb = product % BASE;
            carry = product / BASE;
        }
        while carry > 0 {
            limbs.push(carry % BASE);
            carry /= BASE;
        }
    };
    // 1074 = 13 * 82 + 8; 5^13 = 1,220,703,125 keeps each product below 2^64.
    for _ in 0..82 {
        multiply(1_220_703_125);
    }
    multiply(390_625);
    let mut text = String::new();
    for (index, limb) in limbs.iter().rev().enumerate() {
        if index == 0 {
            text.push_str(&limb.to_string());
        } else {
            text.push_str(&format!("{limb:09}"));
        }
    }
    text
}

/// Strip leading and trailing zeros of a digit string, adjusting the decimal
/// exponent of its last digit.
fn normalized(digits: &str, exponent: i64) -> (&str, i64) {
    let digits = digits.trim_start_matches('0');
    let trimmed = digits.trim_end_matches('0');
    let dropped = i64::try_from(digits.len() - trimmed.len()).unwrap_or(i64::MAX);
    (trimmed, exponent.saturating_add(dropped))
}

/// `std::stod(text)` as glibc's `strtod` decides it: leading whitespace, an
/// optional sign, then the longest valid prefix — `inf`/`infinity`, `nan`
/// (with or without a parenthesised payload), a hexadecimal float after `0x`
/// or a decimal float — the rest ignored. `None` where `stod` throws: no
/// conversion, overflow, or an inexact result below the smallest normal
/// double (measured: `1e-310` and `0x1.8p-1074` throw, `0x1p-1074` and
/// `0e-400` do not).
fn c_stod(text: &str) -> Option<f64> {
    let bytes = text.as_bytes();
    let mut i = bytes.iter().take_while(|b| is_c_space(**b)).count();
    let negative = match bytes.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let sign = if negative { -1.0 } else { 1.0 };
    let rest = &bytes[i..];
    if starts_with_ignore_case(rest, "inf") {
        return Some(sign * f64::INFINITY);
    }
    if starts_with_ignore_case(rest, "nan") {
        return Some(if negative { -f64::NAN } else { f64::NAN });
    }
    if rest.len() >= 2 && rest[0] == b'0' && (rest[1] | 0x20) == b'x' {
        if let Some(value) = c_hex_float(&rest[2..]) {
            return value.map(|value| sign * value);
        }
    }
    // Decimal: digits, an optional point and digits, at least one digit, then
    // an exponent only when it has digits.
    let integer = rest.iter().take_while(|b| b.is_ascii_digit()).count();
    let mut end = integer;
    let mut fraction = 0;
    if rest.get(end) == Some(&b'.') {
        fraction = rest[end + 1..]
            .iter()
            .take_while(|b| b.is_ascii_digit())
            .count();
        if integer + fraction > 0 {
            end += 1 + fraction;
        }
    }
    if integer + fraction == 0 {
        return None;
    }
    let mantissa_end = end;
    let mut exponent: i64 = 0;
    if matches!(rest.get(end), Some(b'e' | b'E')) {
        let mut j = end + 1;
        let exponent_negative = match rest.get(j) {
            Some(b'-') => {
                j += 1;
                true
            }
            Some(b'+') => {
                j += 1;
                false
            }
            _ => false,
        };
        let digits = rest[j..].iter().take_while(|b| b.is_ascii_digit()).count();
        if digits > 0 {
            for &digit in &rest[j..j + digits] {
                exponent = exponent
                    .saturating_mul(10)
                    .saturating_add(i64::from(digit - b'0'))
                    .min(1 << 40);
            }
            if exponent_negative {
                exponent = -exponent;
            }
            end = j + digits;
        }
    }
    // Every byte in rest[..end] is ASCII.
    let literal = std::str::from_utf8(&rest[..end]).ok()?;
    let magnitude: f64 = literal.parse().ok()?;
    if magnitude.is_infinite() {
        return None;
    }
    let digits: String = rest[..mantissa_end]
        .iter()
        .filter(|b| b.is_ascii_digit())
        .map(|b| char::from(*b))
        .collect();
    if magnitude == 0.0 {
        return if digits.bytes().any(|b| b != b'0') {
            None
        } else {
            Some(sign * 0.0)
        };
    }
    if magnitude < f64::MIN_POSITIVE {
        let fraction = i64::try_from(fraction).unwrap_or(i64::MAX);
        let input = normalized(&digits, exponent.saturating_sub(fraction));
        let exact_digits = subnormal_digits(magnitude.to_bits());
        let exact = normalized(&exact_digits, -1074);
        if input != exact {
            return None;
        }
    }
    Some(sign * magnitude)
}

/// The string at `i`, the empty string for a null.
fn text_at(strings: &arrow_array::StringArray, i: usize) -> &str {
    if strings.is_null(i) {
        ""
    } else {
        strings.value(i)
    }
}

/// `s` without a surrounding `[` and `]`, as the source strips a list.
fn unbracketed(text: &str) -> &str {
    text.strip_prefix('[')
        .and_then(|inner| inner.strip_suffix(']'))
        .unwrap_or(text)
}

/// The typed value `readMetaValues` stores for one entry.
fn typed_meta_value(value: &str, value_type: &str) -> MetaValue {
    let parsed = match value_type {
        "int" => c_stoll(value).map(MetaValue::from),
        "double" | "float" => c_stod(value).map(MetaValue::source_float),
        "int_list" => list::create::<i32>(unbracketed(value), b',')
            .ok()
            .map(|values| MetaValue::from(values.into_iter().map(i64::from).collect::<Vec<_>>())),
        "double_list" => list::create::<f64>(unbracketed(value), b',')
            .ok()
            .map(MetaValue::source_float_list),
        "string_list" => list::create::<String>(unbracketed(value), b',')
            .ok()
            .map(|values| {
                MetaValue::from(
                    values
                        .iter()
                        .map(|item| list::trim(item).to_owned())
                        .collect::<Vec<_>>(),
                )
            }),
        _ => None,
    };
    parsed.unwrap_or_else(|| MetaValue::from(value))
}

/// Decode the metadata list at `row` of a `list<struct<name, value,
/// value_type>>` column into `target`.
///
/// Ports `readMetaValues(array, row, target, excluded_keys)`. The struct's
/// first three children are read by position as `name`, `value` and
/// `value_type` (a null reads as the empty string). Each entry not in
/// `excluded_keys` is stored under its name, a later entry replacing an
/// earlier one, typed by `value_type`:
///
/// * `int`: `std::stoll` — leading whitespace, a sign and the leading digit
///   run, the rest ignored (`12abc` and `1.9` give 12 and 1, `0x10` gives 0);
/// * `double` and `float`: `std::stod` as glibc parses it — decimal and
///   hexadecimal floats, `inf`, `infinity` and `nan` in any case, the rest
///   ignored; a value that overflows, or that underflows to an inexact
///   subnormal or zero, is out of range;
/// * `int_list`, `double_list`, `string_list`: an optional surrounding `[`
///   `]` removed, then `ListUtils::create` on `,`, string items trimmed of
///   space, tab, newline and carriage return (not vertical tab);
/// * any other type: the string.
///
/// A value that cannot be converted, or is out of range, is stored as the
/// string it was, as the source's `catch (...)` does. Non-finite floats are
/// stored as the source stores them; [`MetaValue::validate`] reports them.
///
/// A missing array, a null row and a row outside the array add nothing.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the array is not a `list` of a `struct` whose
/// first three children are `utf8` (the source's unchecked cast); `target`
/// is then unchanged.
pub fn read_meta_values(
    array: Option<&dyn Array>,
    row: usize,
    target: &mut MetaInfo,
    excluded_keys: &[&str],
) -> Result<()> {
    let Some(array) = array else {
        return Ok(());
    };
    let refuse = || {
        Error::InvalidValue(format!(
            "ArrowIOHelpers::readMetaValues expects list<struct<name: string, value: string, \
             value_type: string>>, got {}",
            array.data_type()
        ))
    };
    let list = array.as_list_opt::<i32>().ok_or_else(refuse)?;
    let entries = list.values().as_struct_opt().ok_or_else(refuse)?;
    if entries.num_columns() < 3 {
        return Err(refuse());
    }
    let child = |index: usize| {
        entries
            .column(index)
            .as_string_opt::<i32>()
            .ok_or_else(refuse)
    };
    let (names, values, types) = (child(0)?, child(1)?, child(2)?);
    if row >= list.len() || list.is_null(row) {
        return Ok(());
    }
    let offsets = list.value_offsets();
    let (start, end) = match (
        usize::try_from(offsets[row]),
        usize::try_from(offsets[row + 1]),
    ) {
        (Ok(start), Ok(end)) if start <= end && end <= entries.len() => (start, end),
        _ => return Err(refuse()),
    };
    for i in start..end {
        let name = text_at(names, i);
        if excluded_keys.contains(&name) {
            continue;
        }
        let value = typed_meta_value(text_at(values, i), text_at(types, i));
        target.insert(name.to_owned(), value);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stoll_follows_the_c_library() {
        assert_eq!(c_stoll("12"), Some(12));
        assert_eq!(c_stoll("\x0b\x0c 7"), Some(7));
        assert_eq!(c_stoll("- 1"), None);
        assert_eq!(c_stoll("-9223372036854775808"), Some(i64::MIN));
        assert_eq!(c_stoll("9223372036854775808"), None);
    }

    #[test]
    fn stod_reports_range_errors_as_glibc() {
        assert_eq!(c_stod("0x1p-1074"), Some(f64::from_bits(1)));
        assert_eq!(c_stod("0x1.8p-1074"), None);
        assert_eq!(c_stod("1e-310"), None);
        assert_eq!(c_stod("0e-400"), Some(0.0));
        assert_eq!(c_stod("0x"), Some(0.0));
        assert_eq!(c_stod("0x.8"), Some(0.5));
        assert_eq!(c_stod("0x1p1024"), None);
        assert_eq!(c_stod("0x1.fffffffffffffp1023"), Some(f64::MAX));
        // The exact decimal expansion of the smallest subnormal is exact.
        let exact = format!("{}e-1074", subnormal_digits(1));
        assert_eq!(c_stod(&exact), Some(f64::from_bits(1)));
    }
}
