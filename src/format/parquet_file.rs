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
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
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
    pub fn append_or_error<T>(status: std::result::Result<T, ArrowError>, column: &str) -> Result<T> {
        status.map_err(|error| Error::InvalidValue(format!("Failed to append value for {column}: {error}")))
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
        array
            .to_data()
            .validate_full()
            .map_err(|error| Error::InvalidValue(format!("Failed to finish array for {name}: {error}")))?;
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
            .set_coerce_types(true)
            .build();
        let options = ArrowWriterOptions::new()
            .with_properties(properties)
            .with_skip_arrow_metadata(true);
        let bytes = encode(table, options)
            .map_err(|error| Error::InvalidValue(format!("Failed to write parquet table: {error}")))?;
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
    pub fn read_table_with_limits(filename: impl AsRef<Path>, limits: &ReadLimits) -> Result<RecordBatch> {
        let path = filename.as_ref();
        let file = std::fs::File::open(path)
            .map_err(|error| Error::InvalidValue(format!("Failed to open parquet file {}: {error}", path.display())))?;
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
    pub fn read_table_from_with_limits<R: ChunkReader + 'static>(reader: R, limits: &ReadLimits) -> Result<RecordBatch> {
        // arrow-rs reports malformed input as errors; a panic that escapes it
        // anyway becomes a typed error rather than taking the process down.
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| read_all(reader, limits)))
            .map_err(|_| Error::InvalidValue("Parquet reader panicked while decoding the file".into()))?
    }

    /// A required column by name.
    ///
    /// Ports `getColumn(table, name)`.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] "Missing required column '<name>'" when
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
        table.schema().index_of(name).ok().map(|index| table.column(index).clone())
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
    pub fn get_int64(array: Option<&dyn Array>, row: usize, default_value: i64, allow_null: bool) -> Result<i64> {
        let Some(array) = array else {
            return Ok(default_value);
        };
        if checked_null(array, row)? {
            return if allow_null {
                Ok(default_value)
            } else {
                Err(Error::MissingInformation("Required integer value is null".into()))
            };
        }
        integer_value(array, row)?.ok_or_else(|| Error::InvalidValue("Unsupported integer column type".into()))
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
    pub fn get_double(array: Option<&dyn Array>, row: usize, default_value: f64, allow_null: bool) -> Result<f64> {
        let Some(array) = array else {
            return Ok(default_value);
        };
        if checked_null(array, row)? {
            return if allow_null {
                Ok(default_value)
            } else {
                Err(Error::MissingInformation("Required numeric value is null".into()))
            };
        }
        match array.data_type() {
            DataType::Float64 => Ok(array.as_primitive::<arrow_array::types::Float64Type>().value(row)),
            DataType::Float32 => Ok(f64::from(array.as_primitive::<arrow_array::types::Float32Type>().value(row))),
            _ => match integer_value(array, row)? {
                Some(value) => Ok(value as f64),
                None => Err(Error::InvalidValue("Unsupported numeric column type".into())),
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
    pub fn get_bool(array: Option<&dyn Array>, row: usize, default_value: bool, allow_null: bool) -> Result<bool> {
        let Some(array) = array else {
            return Ok(default_value);
        };
        if checked_null(array, row)? {
            return if allow_null {
                Ok(default_value)
            } else {
                Err(Error::MissingInformation("Required boolean value is null".into()))
            };
        }
        if let DataType::Boolean = array.data_type() {
            return Ok(array.as_boolean().value(row));
        }
        match integer_value(array, row)? {
            Some(value) => Ok(value != 0),
            None => Err(Error::InvalidValue("Unsupported boolean column type".into())),
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
        let reader = SerializedFileReader::new(file)
            .map_err(|error| Error::InvalidValue(format!("Failed to read parquet footer of {}: {error}", path.display())))?;
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
    Ok(builder.set_max_row_group_size(usize::try_from(rows).unwrap_or(1 << 20)))
}

/// Encode one table as a complete Parquet file in memory.
pub(crate) fn encode(table: &RecordBatch, options: ArrowWriterOptions) -> std::result::Result<Vec<u8>, parquet::errors::ParquetError> {
    let mut bytes = Vec::new();
    let mut writer = ArrowWriter::try_new_with_options(&mut bytes, table.schema(), options)?;
    writer.write(table)?;
    writer.close()?;
    Ok(bytes)
}

/// Publish complete bytes at `path` atomically (a sibling temporary file,
/// renamed into place), whatever the file name's extension.
pub(crate) fn publish(path: &Path, bytes: &[u8]) -> Result<()> {
    crate::format::path_io::write_plain(path, |writer| {
        writer.write_all(bytes)?;
        Ok(())
    })
}

fn read_all<R: ChunkReader + 'static>(reader: R, limits: &ReadLimits) -> Result<RecordBatch> {
    let builder = ParquetRecordBatchReaderBuilder::try_new(reader)
        .map_err(|error| Error::InvalidValue(format!("Failed to create parquet reader from RandomAccessFile: {error}")))?;
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
        let size = u64::try_from(group.total_byte_size())
            .map_err(|_| Error::InvalidValue("Parquet footer declares a negative row group size".into()))?;
        bytes = bytes.saturating_add(size);
    }
    if bytes > limits.max_uncompressed_bytes {
        return Err(Error::InvalidValue(format!(
            "Parquet file declares {bytes} uncompressed bytes, above the limit of {}",
            limits.max_uncompressed_bytes
        )));
    }
    let schema = builder.schema().clone();
    let batch_size = usize::try_from(rows).unwrap_or(usize::MAX).clamp(1, 1 << 20);
    let batches = builder
        .with_batch_size(batch_size)
        .build()
        .map_err(|error| Error::InvalidValue(format!("Failed to read parquet table from RandomAccessFile: {error}")))?
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|error| Error::InvalidValue(format!("Failed to read parquet table from RandomAccessFile: {error}")))?;
    match batches.len() {
        0 => Ok(RecordBatch::new_empty(schema)),
        1 => Ok(batches.into_iter().next().unwrap_or_else(|| RecordBatch::new_empty(schema))),
        _ => arrow_select::concat::concat_batches(&schema, &batches)
            .map_err(|error| Error::InvalidValue(format!("Failed to combine parquet chunks from RandomAccessFile: {error}"))),
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

fn list_strings<O: OffsetSizeTrait>(list: &arrow_array::GenericListArray<O>, row: usize) -> Result<Vec<String>> {
    let offsets = list.value_offsets();
    let start = offsets[row].as_usize();
    let end = offsets[row + 1].as_usize();
    let values = list.values();
    (start..end).map(|i| ParquetFile::get_string(Some(values.as_ref()), i)).collect()
}
