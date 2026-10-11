// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Shared utilities for reading, writing and accessing Parquet-based formats.
//!
//! Ports `FORMAT/ParquetFile.h` and `ParquetFile.cpp`, the static helpers the
//! Parquet-backed readers and writers share (`TransitionParquetFile`,
//! `OpenSwathOSWParquetWriter`, `XICParquetFile` and their siblings): Arrow
//! builder helpers, reading and writing one table to one Parquet file, column
//! lookup and type-coercing value accessors. See
//! `docs/ARROW_PARQUET_SUPPORT.md` for the full API mapping.
//!
//! **Tables.** A source `std::shared_ptr<arrow::Table>` is an
//! [`arrow_array::RecordBatch`]: one array per column holding every row.
//! `readTable` combines chunks before returning anyway, so nothing the source
//! returns needs more; a caller holding several batches concatenates them
//! with `arrow_select::concat::concat_batches`. This removes one source
//! defect outright: `ParquetFile::getColumn` returns `chunk(0)` only, so a
//! multi-chunk table handed to it directly silently loses every later chunk
//! (measured on the Release build: 2 of 5 rows). A `RecordBatch` has exactly
//! one chunk per column.
//!
//! **Arrays and rows.** A source `std::shared_ptr<arrow::Array>` that may be
//! null is an `Option<&dyn Array>`; a row index is a `usize`. The source reads
//! a row outside the array without checking (undefined behaviour); every
//! accessor here refuses it with [`Error::InvalidValue`].
//!
//! The class header also documents zip/unzip of Parquet directory bundles;
//! the pinned header declares no such member, so there is nothing to port.

use crate::{Error, Result};
use arrow_array::builder::ArrayBuilder;
use arrow_array::cast::AsArray;
use arrow_array::types::{
    Int8Type, Int16Type, Int32Type, Int64Type, UInt8Type, UInt16Type, UInt32Type, UInt64Type,
};
use arrow_array::{Array, ArrayRef, OffsetSizeTrait, RecordBatch};
use arrow_schema::{ArrowError, DataType};
use parquet::arrow::ArrowWriter;
use parquet::arrow::arrow_reader::{
    ArrowReaderMetadata, ArrowReaderOptions, ParquetRecordBatchReaderBuilder,
};
use parquet::arrow::arrow_writer::ArrowWriterOptions;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use parquet::file::reader::{ChunkReader, FileReader};
use parquet::file::serialized_reader::SerializedFileReader;
use std::path::Path;

/// The source's default `row_group_size` for [`ParquetFile::write_table`]:
/// 262,144 rows.
pub const DEFAULT_ROW_GROUP_SIZE: i64 = 262_144;

/// Parquet C++'s `max_row_group_length` default, 1,048,576 rows.
///
/// `parquet::arrow::WriteTable` silently caps its `chunk_size` at this
/// writer property, which neither `ParquetFile` nor `ArrowIOHelpers` sets, so
/// no source row group ever exceeds it (measured: a 1,100,000-row table
/// written with a requested row group size of 134,217,728 becomes groups of
/// 1,048,576 and 51,424 rows). The writers here apply the same cap.
pub const MAX_ROW_GROUP_LENGTH: i64 = 1024 * 1024;

/// Resource ceilings for reading one Parquet file.
///
/// Native: the source reads whatever a file declares. A file is untrusted
/// input, and Parquet metadata declares row counts and uncompressed sizes that
/// drive allocation, so a reader checks them before it decodes anything.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Maximum total rows the footer may declare.
    pub max_rows: u64,
    /// Maximum total uncompressed bytes the footer's row groups may declare.
    pub max_uncompressed_bytes: u64,
}

impl Default for ReadLimits {
    /// 2^32 rows and 16 GiB uncompressed.
    fn default() -> Self {
        Self {
            max_rows: 1 << 32,
            max_uncompressed_bytes: 16 << 30,
        }
    }
}

/// Shared utilities for Parquet-backed file formats.
///
/// Ports `OpenMS::ParquetFile`, a class of static members; each is an
/// associated function here.
#[derive(Clone, Copy, Debug)]
pub struct ParquetFile;

impl ParquetFile {
    /// Turn a failed Arrow builder operation into this crate's error.
    ///
    /// Ports `appendOrThrow(status, column)`. The source takes the
    /// `arrow::Status` an append returned and throws
    /// `Exception::InvalidValue("Failed to append value for <column>")` when
    /// it is not OK; arrow-rs builders mostly append infallibly, and the
    /// fallible ones return a `Result`, which this maps the same way. On
    /// success the value passes through.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] naming the column and Arrow's message.
    pub fn append_or_error<T>(
        status: std::result::Result<T, ArrowError>,
        column: &str,
    ) -> Result<T> {
        status.map_err(|error| {
            Error::InvalidValue(format!("Failed to append value for {column}: {error}"))
        })
    }

    /// Finish an Arrow builder and return the resulting array.
    ///
    /// Ports the template `finishArray<BuilderT>(builder, name)`. arrow-rs's
    /// `finish` cannot report a failure, so the finished array is validated
    /// in full instead, and a failure is the source's
    /// `Exception::InvalidValue("Failed to finish array for <name>")`.
    ///
    /// arrow-rs's `StructBuilder` asserts that its children have equal
    /// lengths by panicking inside `finish`; that is a programming error in
    /// the caller, not something input can reach, and this crate's own
    /// writers build struct arrays with the fallible `StructArray::try_new`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the finished array is not valid.
    pub fn finish_array<B: ArrayBuilder + ?Sized>(builder: &mut B, name: &str) -> Result<ArrayRef> {
        let array = builder.finish();
        array.to_data().validate_full().map_err(|error| {
            Error::InvalidValue(format!("Failed to finish array for {name}: {error}"))
        })?;
        Ok(array)
    }

    /// Write a table to a Parquet file with the source's default row group
    /// size, [`DEFAULT_ROW_GROUP_SIZE`].
    ///
    /// Ports `writeTable(table, filename)`; see
    /// [`write_table_with_row_group_size`](Self::write_table_with_row_group_size).
    ///
    /// # Errors
    ///
    /// As [`write_table_with_row_group_size`](Self::write_table_with_row_group_size).
    pub fn write_table(table: &RecordBatch, filename: impl AsRef<Path>) -> Result<()> {
        Self::write_table_with_row_group_size(table, filename, DEFAULT_ROW_GROUP_SIZE)
    }

    /// Write a table to a Parquet file.
    ///
    /// Ports `writeTable(table, filename, row_group_size)`. The file is
    /// written as `parquet::arrow::WriteTable` writes it without properties:
    /// uncompressed, dictionary encoding and statistics on, and **no stored
    /// Arrow schema**, so a reader sees only the Parquet types — a
    /// `large_string` column reads back as `string`, a `timestamp[s]` as
    /// `timestamp[ms]` and a list's element field as `element` (all measured
    /// on the Release build; lists are written in the Parquet-compliant
    /// layout here too, so the same happens). `row_group_size` is in rows and
    /// capped at [`MAX_ROW_GROUP_LENGTH`], as in the source.
    ///
    /// Native differences: the source opens (and truncates) the destination
    /// first and removes the fragment when the write fails, which also
    /// destroys a file that existed there before; this writes to a sibling
    /// temporary file and renames it only on success, so a failure leaves the
    /// destination as it was. The row group size is checked before anything
    /// is written.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `row_group_size` is not positive and the
    /// table has rows (the source's `Invalid: chunk size per row_group must
    /// be greater than 0`), or when Parquet encoding fails; [`Error::Io`]
    /// when the file cannot be written (the source's
    /// `Exception::FileNotWritable`).
    pub fn write_table_with_row_group_size(
        table: &RecordBatch,
        filename: impl AsRef<Path>,
        row_group_size: i64,
    ) -> Result<()> {
        let properties = writer_properties(table, row_group_size)?
            .set_compression(Compression::UNCOMPRESSED)
            .set_coerce_types(true);
        let bytes = encode(table, properties, false, "Failed to write parquet table")?;
        publish(filename.as_ref(), &bytes)
    }

    /// Read a Parquet file into one table with default [`ReadLimits`].
    ///
    /// Ports `readTable(filename)`. The table holds every row in one array
    /// per column, as the source's combined chunks do, and the file is closed
    /// before this returns, so it can be removed or replaced at once. A file
    /// with no row groups yields a table with its schema and no rows (the
    /// source yields columns with no chunks, which its own `getColumn` then
    /// refuses).
    ///
    /// # Errors
    ///
    /// As [`read_table_with_limits`](Self::read_table_with_limits).
    pub fn read_table(filename: impl AsRef<Path>) -> Result<RecordBatch> {
        Self::read_table_with_limits(filename, &ReadLimits::default())
    }

    /// Read a Parquet file into one table, refusing it when its footer
    /// declares more than `limits` allow.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the file cannot be opened ("Failed to open
    /// parquet file", as the source), is not Parquet, cannot be decoded or
    /// exceeds `limits`.
    pub fn read_table_with_limits(
        filename: impl AsRef<Path>,
        limits: &ReadLimits,
    ) -> Result<RecordBatch> {
        let path = filename.as_ref();
        let file = std::fs::File::open(path).map_err(|error| {
            Error::InvalidValue(format!(
                "Failed to open parquet file {}: {error}",
                path.display()
            ))
        })?;
        Self::read_table_from_with_limits(file, limits)
    }

    /// Read Parquet data from an open reader into one table.
    ///
    /// Ports `readTable(std::shared_ptr<arrow::io::RandomAccessFile>)`, which
    /// reads from an already-open file (for example one inside an archive)
    /// and leaves it open. Pass anything implementing Parquet's
    /// [`ChunkReader`] — a `std::fs::File` (a `try_clone` keeps the caller's
    /// handle open) or an in-memory `bytes::Bytes`.
    ///
    /// # Errors
    ///
    /// As [`read_table_with_limits`](Self::read_table_with_limits), with the
    /// source's messages "Failed to create parquet reader from
    /// RandomAccessFile" and "Failed to read parquet table from
    /// RandomAccessFile".
    pub fn read_table_from<R: ChunkReader + 'static>(reader: R) -> Result<RecordBatch> {
        Self::read_table_from_with_limits(reader, &ReadLimits::default())
    }

    /// [`read_table_from`](Self::read_table_from) with explicit limits.
    ///
    /// # Errors
    ///
    /// As [`read_table_from`](Self::read_table_from).
    pub fn read_table_from_with_limits<R: ChunkReader + 'static>(
        reader: R,
        limits: &ReadLimits,
    ) -> Result<RecordBatch> {
        // arrow-rs reports malformed input as errors; a panic that escapes it
        // anyway becomes a typed error rather than taking the process down.
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            read_all(reader, limits, None)
        }))
        .map_err(|_| {
            Error::InvalidValue("Parquet reader panicked while decoding the file".into())
        })?
    }

    /// A required column by name.
    ///
    /// Ports `getColumn(table, name)`.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] `Missing required column '<name>'` when
    /// the table has no such column. The source's second error, a column with
    /// no chunks, cannot occur for a [`RecordBatch`]. When several columns
    /// share the name, the first is returned; the source's
    /// `Table::GetColumnByName` returns none in that case, so this port
    /// accepts a table the source refuses.
    pub fn get_column(table: &RecordBatch, name: &str) -> Result<ArrayRef> {
        Self::get_optional_column(table, name)
            .ok_or_else(|| Error::MissingInformation(format!("Missing required column '{name}'")))
    }

    /// An optional column by name, or `None` when the table has none.
    ///
    /// Ports `getOptionalColumn(table, name)`. Its "column has no chunks"
    /// error cannot occur for a [`RecordBatch`], so this cannot fail.
    pub fn get_optional_column(table: &RecordBatch, name: &str) -> Option<ArrayRef> {
        table
            .schema()
            .index_of(name)
            .ok()
            .map(|index| table.column(index).clone())
    }

    /// An integer value, coercing any signed or unsigned integer column.
    ///
    /// Ports `getInt64(array, row, default_value, allow_null)`. `None` (the
    /// source's null array pointer) yields `default_value` whatever
    /// `allow_null` says, as in the source; a null value yields it only when
    /// `allow_null` is set.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] "Required integer value is null" for a
    /// null value without `allow_null`; [`Error::InvalidValue`] for a `uint64`
    /// above `i64::MAX` ("Unsigned integer value too large to fit into
    /// int64_t"), for any non-integer column ("Unsupported integer column
    /// type") and for a row outside the array.
    pub fn get_int64(
        array: Option<&dyn Array>,
        row: usize,
        default_value: i64,
        allow_null: bool,
    ) -> Result<i64> {
        let Some(array) = array else {
            return Ok(default_value);
        };
        if checked_null(array, row)? {
            return if allow_null {
                Ok(default_value)
            } else {
                Err(Error::MissingInformation(
                    "Required integer value is null".into(),
                ))
            };
        }
        integer_value(array, row)?
            .ok_or_else(|| Error::InvalidValue("Unsupported integer column type".into()))
    }

    /// A floating-point value, coercing `float` and every integer column.
    ///
    /// Ports `getDouble(array, row, default_value, allow_null)`. An integer
    /// is read through [`get_int64`](Self::get_int64) and converted with
    /// Rust's `as f64`, which rounds to nearest like the source's
    /// `static_cast<double>` (so `i64::MAX` becomes 2^63); a `float` widens
    /// exactly.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] "Required numeric value is null";
    /// [`Error::InvalidValue`] for "Unsupported numeric column type", for the
    /// `uint64` overflow [`get_int64`](Self::get_int64) reports and for a row
    /// outside the array.
    pub fn get_double(
        array: Option<&dyn Array>,
        row: usize,
        default_value: f64,
        allow_null: bool,
    ) -> Result<f64> {
        let Some(array) = array else {
            return Ok(default_value);
        };
        if checked_null(array, row)? {
            return if allow_null {
                Ok(default_value)
            } else {
                Err(Error::MissingInformation(
                    "Required numeric value is null".into(),
                ))
            };
        }
        match array.data_type() {
            DataType::Float64 => Ok(array
                .as_primitive::<arrow_array::types::Float64Type>()
                .value(row)),
            DataType::Float32 => Ok(f64::from(
                array
                    .as_primitive::<arrow_array::types::Float32Type>()
                    .value(row),
            )),
            _ => match integer_value(array, row)? {
                Some(value) => Ok(value as f64),
                None => Err(Error::InvalidValue(
                    "Unsupported numeric column type".into(),
                )),
            },
        }
    }

    /// A boolean value, coercing every integer column (non-zero is `true`).
    ///
    /// Ports `getBool(array, row, default_value, allow_null)`.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] "Required boolean value is null";
    /// [`Error::InvalidValue`] for "Unsupported boolean column type", for the
    /// `uint64` overflow and for a row outside the array.
    pub fn get_bool(
        array: Option<&dyn Array>,
        row: usize,
        default_value: bool,
        allow_null: bool,
    ) -> Result<bool> {
        let Some(array) = array else {
            return Ok(default_value);
        };
        if checked_null(array, row)? {
            return if allow_null {
                Ok(default_value)
            } else {
                Err(Error::MissingInformation(
                    "Required boolean value is null".into(),
                ))
            };
        }
        if let DataType::Boolean = array.data_type() {
            return Ok(array.as_boolean().value(row));
        }
        match integer_value(array, row)? {
            Some(value) => Ok(value != 0),
            None => Err(Error::InvalidValue(
                "Unsupported boolean column type".into(),
            )),
        }
    }

    /// A string value from a `string` or `large_string` column; empty for
    /// `None` or a null value.
    ///
    /// Ports `getString(array, row)`. A null is checked before the type, so a
    /// null in a column of any type is `""`, as in the source.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] "Unsupported string column type" for any other
    /// column type, and for a row outside the array.
    pub fn get_string(array: Option<&dyn Array>, row: usize) -> Result<String> {
        let Some(array) = array else {
            return Ok(String::new());
        };
        if checked_null(array, row)? {
            return Ok(String::new());
        }
        match array.data_type() {
            DataType::Utf8 => Ok(array.as_string::<i32>().value(row).to_owned()),
            DataType::LargeUtf8 => Ok(array.as_string::<i64>().value(row).to_owned()),
            _ => Err(Error::InvalidValue("Unsupported string column type".into())),
        }
    }

    /// A list of strings.
    ///
    /// Ports `getStringList(array, row)`. A `string` or `large_string` value
    /// is split on `;` with every field kept: `"a;;b;"` is `["a", "", "b",
    /// ""]`, `";"` is `["", ""]` and `""` is empty (`StringUtils::split`, as
    /// measured). A `list` or `large_list` value yields its elements through
    /// [`get_string`](Self::get_string), a null element as `""`. `None`, a
    /// null value or an empty list yields an empty vector.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] "Unsupported list column type" for any other
    /// column type, "Unsupported string column type" for a list whose
    /// elements are not strings, and for a row outside the array.
    pub fn get_string_list(array: Option<&dyn Array>, row: usize) -> Result<Vec<String>> {
        let Some(array) = array else {
            return Ok(Vec::new());
        };
        if checked_null(array, row)? {
            return Ok(Vec::new());
        }
        match array.data_type() {
            DataType::Utf8 | DataType::LargeUtf8 => {
                let raw = Self::get_string(Some(array), row)?;
                if raw.is_empty() {
                    return Ok(Vec::new());
                }
                Ok(raw.split(';').map(str::to_owned).collect())
            }
            DataType::List(_) => list_strings(array.as_list::<i32>(), row),
            DataType::LargeList(_) => list_strings(array.as_list::<i64>(), row),
            _ => Err(Error::InvalidValue("Unsupported list column type".into())),
        }
    }

    /// Escape a string for embedding in a JSON string value.
    ///
    /// Ports `jsonEscape(input)`: backslash, double quote, newline, carriage
    /// return and tab get their two-character escapes; every other control
    /// character below U+0020 becomes `\u00XX` with lower-case hex. DEL
    /// (U+007F), `/` and non-ASCII characters pass through unchanged
    /// (measured on the Release build).
    pub fn json_escape(input: &str) -> String {
        const HEX: &[u8; 16] = b"0123456789abcdef";
        let mut out = String::with_capacity(input.len());
        for c in input.chars() {
            match c {
                '\\' => out.push_str("\\\\"),
                '"' => out.push_str("\\\""),
                '\n' => out.push_str("\\n"),
                '\r' => out.push_str("\\r"),
                '\t' => out.push_str("\\t"),
                c if u32::from(c) < 0x20 => {
                    let byte = u32::from(c) as usize;
                    out.push_str("\\u00");
                    out.push(char::from(HEX[byte >> 4]));
                    out.push(char::from(HEX[byte & 0x0f]));
                }
                c => out.push(c),
            }
        }
        out
    }

    /// The number of rows a Parquet file's footer declares; 0 when the file
    /// does not exist.
    ///
    /// Ports `rowCount(filename)`, which reads only the footer.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the file exists but is not Parquet. The
    /// source lets Parquet's own exception escape ("Parquet magic bytes not
    /// found in footer ..."), which is not an OpenMS exception.
    pub fn row_count(filename: impl AsRef<Path>) -> Result<i64> {
        let path = filename.as_ref();
        if !crate::system::file::exists(path) {
            return Ok(0);
        }
        let file = std::fs::File::open(path)?;
        let reader = SerializedFileReader::new(file).map_err(|error| {
            Error::InvalidValue(format!(
                "Failed to read parquet footer of {}: {error}",
                path.display()
            ))
        })?;
        Ok(reader.metadata().file_metadata().num_rows())
    }
}

/// Writer properties shared by [`ParquetFile`] and `ArrowIOHelpers`: the
/// row-group size check and cap. Everything else is left to the caller.
pub(crate) fn writer_properties(
    table: &RecordBatch,
    row_group_size: i64,
) -> Result<parquet::file::properties::WriterPropertiesBuilder> {
    let builder = WriterProperties::builder();
    if row_group_size <= 0 {
        if table.num_rows() > 0 {
            return Err(Error::InvalidValue(
                "Invalid: chunk size per row_group must be greater than 0".into(),
            ));
        }
        return Ok(builder);
    }
    let rows = row_group_size.min(MAX_ROW_GROUP_LENGTH);
    // `rows` is in 1..=2^20, so it fits every usize.
    Ok(builder.set_max_row_group_row_count(Some(usize::try_from(rows).unwrap_or(1 << 20))))
}

/// Encode one table as a complete Parquet file in memory.
///
/// Two things are done as Parquet C++'s `WriteTable` does them, which
/// arrow-rs does differently by itself:
///
/// * a `timestamp[s]` column, which Parquet cannot store, is written as
///   `timestamp[ms]` with its values multiplied by 1000 (see
///   [`coerce_second_timestamps`]); arrow-rs would store the raw seconds as a
///   plain `INT64`;
/// * a table with no rows is written with one empty row group whose column
///   chunks carry no statistics (measured on the Release build: `rows=0
///   row_groups=1 ... stats0=false`); arrow-rs writes no row group.
///
/// `store_schema` embeds the Arrow schema (Arrow C++'s `store_schema()`);
/// `failure` prefixes an encoding error's message.
pub(crate) fn encode(
    table: &RecordBatch,
    properties: parquet::file::properties::WriterPropertiesBuilder,
    store_schema: bool,
    failure: &str,
) -> Result<Vec<u8>> {
    let coerced = coerce_second_timestamps(table)?;
    let table = coerced.as_ref().unwrap_or(table);
    let fail =
        |error: parquet::errors::ParquetError| Error::InvalidValue(format!("{failure}: {error}"));
    let properties = if table.num_rows() == 0 {
        properties.set_statistics_enabled(parquet::file::properties::EnabledStatistics::None)
    } else {
        properties
    };
    let (bare, options) = writer_options(&table.schema(), properties, store_schema);
    let table = RecordBatch::try_new(bare.clone(), table.columns().to_vec())
        .map_err(|error| Error::InvalidValue(format!("{failure}: {error}")))?;
    let mut bytes = Vec::new();
    let mut writer = ArrowWriter::try_new_with_options(&mut bytes, bare, options).map_err(fail)?;
    if table.num_rows() == 0 {
        let (mut file, factory) = writer.into_serialized_writer().map_err(fail)?;
        write_row_group(&mut file, &factory, 0, &table, failure)?;
        file.close().map_err(fail)?;
    } else {
        writer.write(&table).map_err(fail)?;
        writer.close().map_err(fail)?;
    }
    Ok(bytes)
}

/// The stored schema and writer options of an Arrow-to-Parquet writer.
///
/// Schema metadata goes into the footer as key-value pairs in key order, as
/// Parquet C++ writes it (in insertion order there), and the stored Arrow
/// schema carries none: an arrow-rs schema keeps its metadata in a randomly
/// seeded hash map, so encoding it would make the file's bytes vary from run
/// to run. Both readers merge the footer pairs back into the schema
/// metadata. The returned schema is `schema` without metadata; the tables
/// written must use it. `store_schema` embeds it (Arrow C++'s
/// `store_schema()`).
pub(crate) fn writer_options(
    schema: &arrow_schema::Schema,
    properties: parquet::file::properties::WriterPropertiesBuilder,
    store_schema: bool,
) -> (arrow_schema::SchemaRef, ArrowWriterOptions) {
    let mut pairs: Vec<(&String, &String)> = schema.metadata().iter().collect();
    pairs.sort();
    let pairs: Vec<parquet::file::metadata::KeyValue> = pairs
        .into_iter()
        .map(|(key, value)| parquet::file::metadata::KeyValue {
            key: key.clone(),
            value: Some(value.clone()),
        })
        .collect();
    let bare = std::sync::Arc::new(arrow_schema::Schema::new(schema.fields().clone()));
    let mut properties = properties
        .set_key_value_metadata((!pairs.is_empty()).then_some(pairs))
        .build();
    if store_schema {
        parquet::arrow::add_encoded_arrow_schema_to_metadata(&bare, &mut properties);
    }
    let options = ArrowWriterOptions::new()
        .with_properties(properties)
        .with_skip_arrow_metadata(true);
    (bare, options)
}

/// Write `table` as the one row group `index` of `file`, whatever its row
/// count: a table with no rows becomes an empty row group, as Parquet C++'s
/// `WriteTable` writes one.
///
/// `table` must have the schema the writer was created with.
pub(crate) fn write_row_group<W: std::io::Write + Send>(
    file: &mut parquet::file::writer::SerializedFileWriter<W>,
    factory: &parquet::arrow::arrow_writer::ArrowRowGroupWriterFactory,
    index: usize,
    table: &RecordBatch,
    failure: &str,
) -> Result<()> {
    let fail =
        |error: parquet::errors::ParquetError| Error::InvalidValue(format!("{failure}: {error}"));
    let mut columns = factory.create_column_writers(index).map_err(fail)?;
    if table.num_rows() > 0 {
        let schema = table.schema();
        let mut writers = columns.iter_mut();
        for (field, array) in schema.fields().iter().zip(table.columns()) {
            for leaf in parquet::arrow::arrow_writer::compute_leaves(field, array).map_err(fail)? {
                let writer = writers.next().ok_or_else(|| {
                    Error::InvalidValue(format!("{failure}: more leaf columns than writers"))
                })?;
                writer.write(&leaf).map_err(fail)?;
            }
        }
    }
    let mut group = file.next_row_group().map_err(fail)?;
    for column in columns {
        column
            .close()
            .map_err(fail)?
            .append_to_row_group(&mut group)
            .map_err(fail)?;
    }
    group.close().map_err(fail)?;
    Ok(())
}

/// `table` with every `timestamp[s]` column, top level or inside a list,
/// large list, fixed-size list or struct, converted to `timestamp[ms]`, or
/// `None` when it has none.
///
/// Parquet has no seconds unit; Parquet C++ writes such a column as
/// milliseconds, and so reads it back as `timestamp[ms]` even when the Arrow
/// schema is stored (measured: `1700000000` s reads back as
/// `1700000000000` ms). A seconds column nested in any other type is left
/// as it is.
///
/// # Errors
///
/// [`Error::InvalidValue`] for a value whose millisecond count overflows
/// `i64`, which Arrow C++ refuses as an out-of-bounds timestamp too.
pub(crate) fn coerce_second_timestamps(table: &RecordBatch) -> Result<Option<RecordBatch>> {
    let mut changed = false;
    let mut columns = Vec::with_capacity(table.num_columns());
    for column in table.columns() {
        match coerce_array(column)? {
            Some(new) => {
                changed = true;
                columns.push(new);
            }
            None => columns.push(column.clone()),
        }
    }
    if !changed {
        return Ok(None);
    }
    let schema = table.schema();
    let fields: Vec<arrow_schema::Field> = schema
        .fields()
        .iter()
        .zip(&columns)
        .map(|(field, column)| {
            field
                .as_ref()
                .clone()
                .with_data_type(column.data_type().clone())
        })
        .collect();
    let schema = std::sync::Arc::new(
        arrow_schema::Schema::new(fields).with_metadata(schema.metadata().clone()),
    );
    RecordBatch::try_new(schema, columns)
        .map(Some)
        .map_err(|error| Error::InvalidValue(format!("Failed to convert timestamp[s]: {error}")))
}

fn converted_field(field: &arrow_schema::FieldRef, child: &ArrayRef) -> arrow_schema::FieldRef {
    std::sync::Arc::new(
        field
            .as_ref()
            .clone()
            .with_data_type(child.data_type().clone()),
    )
}

fn conversion_error(error: ArrowError) -> Error {
    Error::InvalidValue(format!("Failed to convert timestamp[s]: {error}"))
}

fn coerce_list<O: OffsetSizeTrait>(
    list: &arrow_array::GenericListArray<O>,
    field: &arrow_schema::FieldRef,
) -> Result<Option<ArrayRef>> {
    let Some(values) = coerce_array(list.values())? else {
        return Ok(None);
    };
    arrow_array::GenericListArray::<O>::try_new(
        converted_field(field, &values),
        list.offsets().clone(),
        values,
        list.nulls().cloned(),
    )
    .map(|array| Some(std::sync::Arc::new(array) as ArrayRef))
    .map_err(conversion_error)
}

fn coerce_array(array: &ArrayRef) -> Result<Option<ArrayRef>> {
    use arrow_array::types::TimestampSecondType;
    use arrow_array::{FixedSizeListArray, StructArray, TimestampMillisecondArray};
    use std::sync::Arc;
    match array.data_type() {
        DataType::Timestamp(arrow_schema::TimeUnit::Second, zone) => {
            let seconds = array.as_primitive::<TimestampSecondType>();
            let mut values = Vec::with_capacity(seconds.len());
            for index in 0..seconds.len() {
                if seconds.is_null(index) {
                    values.push(0);
                    continue;
                }
                let value = seconds.value(index);
                values.push(value.checked_mul(1000).ok_or_else(|| {
                    Error::InvalidValue(format!(
                        "Casting from timestamp[s] to timestamp[ms] would result in out of \
                         bounds timestamp: {value}"
                    ))
                })?);
            }
            let millis = TimestampMillisecondArray::new(values.into(), seconds.nulls().cloned())
                .with_timezone_opt(zone.clone());
            Ok(Some(Arc::new(millis)))
        }
        DataType::List(field) => coerce_list(array.as_list::<i32>(), field),
        DataType::LargeList(field) => coerce_list(array.as_list::<i64>(), field),
        DataType::FixedSizeList(field, size) => {
            let fixed = array.as_fixed_size_list();
            let Some(values) = coerce_array(fixed.values())? else {
                return Ok(None);
            };
            FixedSizeListArray::try_new(
                converted_field(field, &values),
                *size,
                values,
                fixed.nulls().cloned(),
            )
            .map(|array| Some(Arc::new(array) as ArrayRef))
            .map_err(conversion_error)
        }
        DataType::Struct(fields) => {
            let strukt = array.as_struct();
            let mut changed = false;
            let mut children = Vec::with_capacity(strukt.num_columns());
            for child in strukt.columns() {
                match coerce_array(child)? {
                    Some(new) => {
                        changed = true;
                        children.push(new);
                    }
                    None => children.push(child.clone()),
                }
            }
            if !changed {
                return Ok(None);
            }
            let fields: arrow_schema::Fields = fields
                .iter()
                .zip(&children)
                .map(|(field, child)| converted_field(field, child))
                .collect();
            StructArray::try_new(fields, children, strukt.nulls().cloned())
                .map(|array| Some(Arc::new(array) as ArrayRef))
                .map_err(conversion_error)
        }
        _ => Ok(None),
    }
}

/// The type with every `timestamp[s]`, at any list or struct depth, turned
/// into `timestamp[ms]`; `None` when it has none.
fn millisecond_type(data_type: &DataType) -> Option<DataType> {
    use std::sync::Arc;
    let field = |f: &arrow_schema::FieldRef| {
        millisecond_type(f.data_type()).map(|t| Arc::new(f.as_ref().clone().with_data_type(t)))
    };
    match data_type {
        DataType::Timestamp(arrow_schema::TimeUnit::Second, zone) => Some(DataType::Timestamp(
            arrow_schema::TimeUnit::Millisecond,
            zone.clone(),
        )),
        DataType::List(f) => field(f).map(DataType::List),
        DataType::LargeList(f) => field(f).map(DataType::LargeList),
        DataType::FixedSizeList(f, n) => field(f).map(|f| DataType::FixedSizeList(f, *n)),
        DataType::Struct(fields) => {
            let converted: Vec<Option<arrow_schema::FieldRef>> = fields.iter().map(field).collect();
            if converted.iter().all(Option::is_none) {
                return None;
            }
            Some(DataType::Struct(
                fields
                    .iter()
                    .zip(converted)
                    .map(|(original, new)| new.unwrap_or_else(|| original.clone()))
                    .collect(),
            ))
        }
        _ => None,
    }
}

/// `schema` with every `timestamp[s]` field turned into `timestamp[ms]`;
/// `None` when it has none.
fn millisecond_schema(schema: &arrow_schema::Schema) -> Option<arrow_schema::Schema> {
    let converted: Vec<Option<DataType>> = schema
        .fields()
        .iter()
        .map(|f| millisecond_type(f.data_type()))
        .collect();
    if converted.iter().all(Option::is_none) {
        return None;
    }
    let fields: Vec<arrow_schema::Field> = schema
        .fields()
        .iter()
        .zip(converted)
        .map(|(f, t)| match t {
            Some(t) => f.as_ref().clone().with_data_type(t),
            None => f.as_ref().clone(),
        })
        .collect();
    Some(arrow_schema::Schema::new(fields).with_metadata(schema.metadata().clone()))
}

/// Publish complete bytes at `path` atomically (a sibling temporary file,
/// renamed into place), whatever the file name's extension.
pub(crate) fn publish(path: &Path, bytes: &[u8]) -> Result<()> {
    crate::format::path_io::write_plain(path, |writer| {
        writer.write_all(bytes)?;
        Ok(())
    })
}

/// Read the named columns of a Parquet file into one table.
///
/// The crate-private counterpart of `ParquetReaderHelpers::readParquetTableColumns`
/// for one file: requested names absent from the file are skipped, and the
/// columns come back in file order. `None` for `columns` reads every column.
/// The footer is checked against `limits` before anything is decoded, and a
/// panic escaping arrow-rs becomes a typed error, as for
/// [`ParquetFile::read_table_from_with_limits`].
///
/// # Errors
///
/// [`Error::InvalidValue`] when the file cannot be opened or decoded or
/// exceeds `limits`; [`Error::MissingInformation`] when `columns` names no
/// column the file has.
pub(crate) fn read_file_columns(
    path: &Path,
    columns: Option<&[&str]>,
    limits: &ReadLimits,
) -> Result<RecordBatch> {
    let file = std::fs::File::open(path).map_err(|error| {
        Error::InvalidValue(format!(
            "Failed to open parquet file {}: {error}",
            path.display()
        ))
    })?;
    let projection = columns.map(|names| (names, path));
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        read_all(file, limits, projection)
    }))
    .map_err(|_| Error::InvalidValue("Parquet reader panicked while decoding the file".into()))?
}

/// The Arrow schema a Parquet file declares, without decoding any data.
///
/// The counterpart of `ParquetReaderHelpers::readParquetSchema`: the stored
/// Arrow schema when the file has one, else the one derived from the Parquet
/// types, with the same `timestamp[s]` to `timestamp[ms]` coercion the table
/// readers apply.
///
/// # Errors
///
/// [`Error::InvalidValue`] "Failed to read parquet schema" when the file
/// cannot be opened or its footer cannot be decoded.
pub(crate) fn read_file_schema(path: &Path) -> Result<arrow_schema::SchemaRef> {
    let failed = |detail: String| {
        Error::InvalidValue(format!(
            "Failed to read parquet schema {}: {detail}",
            path.display()
        ))
    };
    let file = std::fs::File::open(path).map_err(|error| failed(error.to_string()))?;
    let metadata = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        ArrowReaderMetadata::load(&file, ArrowReaderOptions::new())
    }))
    .map_err(|_| failed("the Parquet reader panicked".into()))?
    .map_err(|error| failed(error.to_string()))?;
    Ok(match millisecond_schema(metadata.schema()) {
        Some(schema) => std::sync::Arc::new(schema),
        None => metadata.schema().clone(),
    })
}

fn read_all<R: ChunkReader + 'static>(
    reader: R,
    limits: &ReadLimits,
    projection: Option<(&[&str], &Path)>,
) -> Result<RecordBatch> {
    let create = |error: parquet::errors::ParquetError| {
        Error::InvalidValue(format!(
            "Failed to create parquet reader from RandomAccessFile: {error}"
        ))
    };
    let mut metadata =
        ArrowReaderMetadata::load(&reader, ArrowReaderOptions::new()).map_err(create)?;
    // Parquet C++ reads a stored `timestamp[s]` column as `timestamp[ms]`;
    // arrow-rs would restore the seconds from the stored Arrow schema.
    if let Some(schema) = millisecond_schema(metadata.schema()) {
        metadata = ArrowReaderMetadata::try_new(
            metadata.metadata().clone(),
            ArrowReaderOptions::new().with_schema(std::sync::Arc::new(schema)),
        )
        .map_err(create)?;
    }
    let mut builder = ParquetRecordBatchReaderBuilder::new_with_metadata(reader, metadata);
    let mut schema = builder.schema().clone();
    if let Some((names, path)) = projection {
        // Arrow C++'s Schema::GetFieldIndex: a name found exactly once.
        let fields = builder.schema().fields();
        let mut indices: Vec<usize> = names
            .iter()
            .filter_map(|name| {
                let mut matches = fields
                    .iter()
                    .enumerate()
                    .filter(|(_, field)| field.name() == name);
                match (matches.next(), matches.next()) {
                    (Some((index, _)), None) => Some(index),
                    _ => None,
                }
            })
            .collect();
        if indices.is_empty() {
            return Err(Error::MissingInformation(format!(
                "No requested columns found in parquet file '{}'",
                path.display()
            )));
        }
        // Batches come back in file order, whatever order was requested.
        indices.sort_unstable();
        indices.dedup();
        schema = std::sync::Arc::new(schema.project(&indices).map_err(|error| {
            Error::InvalidValue(format!("Failed to read parquet table: {error}"))
        })?);
        let mask = parquet::arrow::ProjectionMask::roots(builder.parquet_schema(), indices);
        builder = builder.with_projection(mask);
    }
    let metadata = builder.metadata().clone();
    let rows = u64::try_from(metadata.file_metadata().num_rows())
        .map_err(|_| Error::InvalidValue("Parquet footer declares a negative row count".into()))?;
    if rows > limits.max_rows {
        return Err(Error::InvalidValue(format!(
            "Parquet file declares {rows} rows, above the limit of {}",
            limits.max_rows
        )));
    }
    let mut bytes: u64 = 0;
    for group in metadata.row_groups() {
        let size = u64::try_from(group.total_byte_size()).map_err(|_| {
            Error::InvalidValue("Parquet footer declares a negative row group size".into())
        })?;
        bytes = bytes.saturating_add(size);
    }
    if bytes > limits.max_uncompressed_bytes {
        return Err(Error::InvalidValue(format!(
            "Parquet file declares {bytes} uncompressed bytes, above the limit of {}",
            limits.max_uncompressed_bytes
        )));
    }
    let batch_size = usize::try_from(rows)
        .unwrap_or(usize::MAX)
        .clamp(1, 1 << 20);
    let batches = builder
        .with_batch_size(batch_size)
        .build()
        .map_err(|error| {
            Error::InvalidValue(format!(
                "Failed to read parquet table from RandomAccessFile: {error}"
            ))
        })?
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|error| {
            Error::InvalidValue(format!(
                "Failed to read parquet table from RandomAccessFile: {error}"
            ))
        })?;
    // The decoded batches carry no schema metadata; the table does.
    match batches.len() {
        0 => Ok(RecordBatch::new_empty(schema)),
        1 => match batches.into_iter().next() {
            Some(batch) => batch.with_schema(schema).map_err(|error| {
                Error::InvalidValue(format!(
                    "Failed to read parquet table from RandomAccessFile: {error}"
                ))
            }),
            None => Ok(RecordBatch::new_empty(schema)),
        },
        _ => arrow_select::concat::concat_batches(&schema, &batches).map_err(|error| {
            Error::InvalidValue(format!(
                "Failed to combine parquet chunks from RandomAccessFile: {error}"
            ))
        }),
    }
}

/// Whether `row` is null, refusing a row outside the array.
pub(crate) fn checked_null(array: &dyn Array, row: usize) -> Result<bool> {
    if row >= array.len() {
        return Err(Error::InvalidValue(format!(
            "row {row} is outside an array of {} values",
            array.len()
        )));
    }
    Ok(array.is_null(row))
}

/// The integer at `row` widened to `i64`, `None` for a non-integer column.
fn integer_value(array: &dyn Array, row: usize) -> Result<Option<i64>> {
    Ok(Some(match array.data_type() {
        DataType::Int64 => array.as_primitive::<Int64Type>().value(row),
        DataType::Int32 => i64::from(array.as_primitive::<Int32Type>().value(row)),
        DataType::Int16 => i64::from(array.as_primitive::<Int16Type>().value(row)),
        DataType::Int8 => i64::from(array.as_primitive::<Int8Type>().value(row)),
        DataType::UInt64 => {
            let value = array.as_primitive::<UInt64Type>().value(row);
            i64::try_from(value).map_err(|_| {
                Error::InvalidValue(format!(
                    "Unsigned integer value too large to fit into int64_t: {value}"
                ))
            })?
        }
        DataType::UInt32 => i64::from(array.as_primitive::<UInt32Type>().value(row)),
        DataType::UInt16 => i64::from(array.as_primitive::<UInt16Type>().value(row)),
        DataType::UInt8 => i64::from(array.as_primitive::<UInt8Type>().value(row)),
        _ => return Ok(None),
    }))
}

fn list_strings<O: OffsetSizeTrait>(
    list: &arrow_array::GenericListArray<O>,
    row: usize,
) -> Result<Vec<String>> {
    let offsets = list.value_offsets();
    let start = offsets[row].as_usize();
    let end = offsets[row + 1].as_usize();
    let values = list.values();
    (start..end)
        .map(|i| ParquetFile::get_string(Some(values.as_ref()), i))
        .collect()
}
