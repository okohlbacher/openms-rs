// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Read SWATH/DIA windows and spectrum IDs from an sqMass database.
//!
//! This maps the installed `Internal::MzMLSqliteSwathHandler` API. It opens
//! each query read-only, bounds scanned records, and returns deterministic
//! ordering. See `docs/MZML_SQLITE_SWATH_SUPPORT.md` for source differences.

use super::sqlite_connector::{SqlOpenMode, SqliteConnector};
use crate::{Error, Result};
use rusqlite::{Connection, Row};
use std::collections::HashSet;
use std::path::{Path, PathBuf};

/// Isolation center and absolute boundaries of a SWATH window, in m/z.
///
/// This is the populated value subset of the source `OpenSwath::SwathMap`.
/// Its unpopulated `ms1=false`, ion-mobility limits `-1`, and absent spectrum
/// access pointer are not represented here. This is not a full SwathMap port.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SwathWindow {
    /// Precursor isolation center, in m/z.
    pub center: f64,
    /// Absolute lower boundary (center minus the lower offset), in m/z.
    pub lower: f64,
    /// Absolute upper boundary (center plus the upper offset), in m/z.
    pub upper: f64,
}

/// Native ceiling on source rows examined by one accessor.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SwathReadLimits {
    /// Maximum total rows across the tables scanned by one call.
    ///
    /// Defaults to 1,000,000. Window discovery counts SPECTRUM and PRECURSOR
    /// rows together; the other operations scan one table. Nonmatching rows
    /// count too. Zero accepts only empty scanned tables.
    pub max_records: usize,
}

impl Default for SwathReadLimits {
    fn default() -> Self {
        Self {
            max_records: 1_000_000,
        }
    }
}

/// Read-only accessor for SWATH windows and spectrum IDs in one sqMass file.
///
/// Construction performs no I/O. Every accessor opens and closes its own
/// connection. Returned IDs identify database records, not vector positions.
#[derive(Clone, Debug)]
pub struct MzMLSqliteSwathHandler {
    filename: PathBuf,
    limits: SwathReadLimits,
}

impl MzMLSqliteSwathHandler {
    /// Remember a database path, using the default record ceiling.
    ///
    /// The file is opened only when an accessor is called.
    pub fn new(filename: impl AsRef<Path>) -> Self {
        Self::with_limits(filename, SwathReadLimits::default())
    }

    /// Remember a database path and an explicit native record ceiling.
    ///
    /// This performs no I/O; a zero ceiling is valid for empty tables.
    pub fn with_limits(filename: impl AsRef<Path>, limits: SwathReadLimits) -> Self {
        Self {
            filename: filename.as_ref().to_owned(),
            limits,
        }
    }

    /// Read distinct `(center, lower, upper)` windows belonging to MS2 spectra.
    ///
    /// Boundaries are absolute m/z values. Equal centers with different bounds
    /// remain distinct, matching the source SQL despite its center-only docs.
    /// Results are sorted lexicographically by center, lower, then upper;
    /// the source promises only natural database order. A read transaction
    /// keeps both table scans in one snapshot.
    ///
    /// # Errors
    /// Returns [`Error::Io`] for open, schema, query, or SQL type errors.
    /// [`Error::InvalidValue`] reports negative spectrum IDs,
    /// nonfinite selected window values (including arithmetic overflow), a missing/nonordinary table, or
    /// exceeding the shared record ceiling. Required NULL coordinates are SQL
    /// type errors. No partial result is returned, and a missing file is not
    /// created.
    pub fn read_swath_windows(&self) -> Result<Vec<SwathWindow>> {
        let connector = SqliteConnector::open(&self.filename, SqlOpenMode::ReadOnly)?;
        let transaction = connector
            .connection()
            .unchecked_transaction()
            .map_err(sql_error)?;
        require_table(&transaction, "SPECTRUM")?;
        require_table(&transaction, "PRECURSOR")?;
        let mut remaining = self.limits.max_records;
        let mut ms2 = HashSet::new();
        scan(
            &transaction,
            "SELECT ID, MSLEVEL FROM main.SPECTRUM",
            &mut remaining,
            |row| {
                if row.get::<_, Option<i64>>(1).map_err(sql_error)? == Some(2) {
                    ms2.insert(record_id(row, 0)?);
                }
                Ok(())
            },
        )?;
        let mut windows = Vec::new();
        scan(
            &transaction,
            "SELECT SPECTRUM_ID, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER FROM main.PRECURSOR",
            &mut remaining,
            |row| {
                let Some(id) = row.get::<_, Option<i64>>(0).map_err(sql_error)? else {
                    return Ok(());
                };
                if !ms2.contains(&checked_id(id)?) {
                    return Ok(());
                }
                let center = finite(row.get(1).map_err(sql_error)?)?;
                let lower_offset = finite(row.get(2).map_err(sql_error)?)?;
                let upper_offset = finite(row.get(3).map_err(sql_error)?)?;
                windows.push(SwathWindow {
                    center,
                    lower: finite(center - lower_offset)?,
                    upper: finite(center + upper_offset)?,
                });
                Ok(())
            },
        )?;
        windows.sort_by(|a, b| {
            a.center
                .total_cmp(&b.center)
                .then(a.lower.total_cmp(&b.lower))
                .then(a.upper.total_cmp(&b.upper))
        });
        windows.dedup();
        Ok(windows)
    }

    /// Read all MS1 spectrum record IDs, sorted in ascending order.
    ///
    /// As in the source, MSLEVEL must equal 1; NULL levels do not match.
    /// The source's natural database ordering is replaced by explicit ID order.
    ///
    /// # Errors
    /// Returns [`Error::Io`] for open, schema, query, or SQL type errors, and
    /// [`Error::InvalidValue`] for an invalid matching ID, missing/nonordinary table, or record ceiling
    /// violation. Failed reads do not create files or return partial results.
    pub fn read_ms1_spectra(&self) -> Result<Vec<i64>> {
        let connector = SqliteConnector::open(&self.filename, SqlOpenMode::ReadOnly)?;
        let transaction = connector
            .connection()
            .unchecked_transaction()
            .map_err(sql_error)?;
        require_table(&transaction, "SPECTRUM")?;
        let mut remaining = self.limits.max_records;
        let mut indices = Vec::new();
        scan(
            &transaction,
            "SELECT ID, MSLEVEL FROM main.SPECTRUM",
            &mut remaining,
            |row| {
                if row.get::<_, Option<i64>>(1).map_err(sql_error)? == Some(1) {
                    indices.push(record_id(row, 0)?);
                }
                Ok(())
            },
        )?;
        indices.sort_unstable();
        Ok(indices)
    }

    /// Read spectrum IDs with precursor centers within `window.center ± 0.01`.
    ///
    /// Both endpoints are included. Only `center` is consulted: `lower` and
    /// `upper` are ignored, even when nonfinite. Like the source, this lookup
    /// does not filter by MS level or check membership in SPECTRUM, and repeated
    /// precursor rows produce repeated IDs. Results are sorted by ID. Rows
    /// belonging to chromatograms (NULL SPECTRUM_ID) are skipped rather than
    /// terminating the scan; NULL isolation targets do not match.
    ///
    /// # Errors
    /// Returns [`Error::InvalidValue`] for a nonfinite center/target, invalid
    /// matching spectrum ID, missing/nonordinary table, or record ceiling violation; [`Error::Io`] for
    /// open, schema, query, or SQL type errors. No partial result or missing
    /// input file is created.
    pub fn read_spectra_for_window(&self, window: &SwathWindow) -> Result<Vec<i64>> {
        let center = finite(window.center)?;
        let lower = finite(center - 0.01)?;
        let upper = finite(center + 0.01)?;
        let connector = SqliteConnector::open(&self.filename, SqlOpenMode::ReadOnly)?;
        let transaction = connector
            .connection()
            .unchecked_transaction()
            .map_err(sql_error)?;
        require_table(&transaction, "PRECURSOR")?;
        let mut remaining = self.limits.max_records;
        let mut indices = Vec::new();
        scan(
            &transaction,
            "SELECT SPECTRUM_ID, ISOLATION_TARGET FROM main.PRECURSOR",
            &mut remaining,
            |row| {
                let Some(id) = row.get::<_, Option<i64>>(0).map_err(sql_error)? else {
                    return Ok(());
                };
                let Some(target) = row.get::<_, Option<f64>>(1).map_err(sql_error)? else {
                    return Ok(());
                };
                let target = finite(target)?;
                if lower <= target && target <= upper {
                    indices.push(checked_id(id)?);
                }
                Ok(())
            },
        )?;
        indices.sort_unstable();
        Ok(indices)
    }
}

fn require_table(connection: &Connection, name: &str) -> Result<()> {
    let ordinary: bool = connection.query_row(
        "SELECT EXISTS(SELECT 1 FROM pragma_table_list() WHERE schema='main' AND name=?1 COLLATE NOCASE AND type='table')",
        [name], |row| row.get(0),
    ).map_err(sql_error)?;
    if !ordinary {
        return Err(Error::InvalidValue(format!(
            "sqMass requires an ordinary {name} table"
        )));
    }
    Ok(())
}

fn scan(
    connection: &Connection,
    query: &str,
    remaining: &mut usize,
    mut visit: impl FnMut(&Row<'_>) -> Result<()>,
) -> Result<()> {
    let mut statement = connection.prepare(query).map_err(sql_error)?;
    let mut rows = statement.query([]).map_err(sql_error)?;
    while let Some(row) = rows.next().map_err(sql_error)? {
        *remaining = remaining
            .checked_sub(1)
            .ok_or_else(|| Error::InvalidValue("sqMass SWATH record ceiling exceeded".into()))?;
        visit(row)?;
    }
    Ok(())
}

fn record_id(row: &Row<'_>, column: usize) -> Result<i64> {
    checked_id(row.get(column).map_err(sql_error)?)
}

fn checked_id(value: i64) -> Result<i64> {
    if value < 0 {
        Err(Error::InvalidValue("sqMass spectrum ID is negative".into()))
    } else {
        Ok(value)
    }
}

fn finite(value: f64) -> Result<f64> {
    if !value.is_finite() {
        return Err(Error::InvalidValue(
            "sqMass SWATH coordinate is not finite".into(),
        ));
    }
    // SQLite DISTINCT considers negative and positive zero equal.
    Ok(if value == 0.0 { 0.0 } else { value })
}

fn sql_error(error: rusqlite::Error) -> Error {
    Error::Io(std::io::Error::other(error))
}
