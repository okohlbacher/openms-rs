// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! SQLite connection helpers from `FORMAT/SqliteConnector.h`.
//!
//! The connection and prepared statements remain private and are closed by
//! ownership on every exit. SQL execution is caller-controlled: no runtime or
//! memory ceilings, implicit transactions, or automatic rollback are added.
//! See `docs/SQLITE_CONNECTOR_SUPPORT.md` for the source mapping and differences.

use crate::{Error, Result};
use rusqlite::fallible_iterator::FallibleIterator;
use rusqlite::{Batch, Connection, OpenFlags};
use std::path::Path;
use std::time::Duration;

/// How to open a SQLite database, matching the source's three distinct modes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SqlOpenMode {
    /// Open an existing database without permission to modify it.
    ReadOnly,
    /// Open an existing database for reading and writing; never create it.
    ReadWrite,
    /// Open for reading and writing, creating the database if necessary.
    #[default]
    ReadWriteOrCreate,
}

/// An owned SQLite connection with the public OpenMS connector operations.
///
/// This type is not clonable and exposes no native handle. Destruction closes
/// the connection and rolls back an outstanding transaction, as SQLite does.
#[derive(Debug)]
pub struct SqliteConnector {
    connection: Connection,
}

impl SqliteConnector {
    /// Borrow the owned connection for crate-internal typed format queries.
    pub(crate) fn connection(&self) -> &Connection {
        &self.connection
    }

    /// Open a database in [`SqlOpenMode::ReadWriteOrCreate`] mode.
    ///
    /// SQLite's `:memory:` filename creates a private in-memory database; an
    /// empty filename creates a private temporary database.
    ///
    /// # Errors
    /// Returns [`Error::Io`] if SQLite cannot open or initialize the connection.
    /// This replaces the source's `Exception::SqlOperationFailed`.
    pub fn new(filename: impl AsRef<Path>) -> Result<Self> {
        Self::open(filename, SqlOpenMode::ReadWriteOrCreate)
    }

    /// Open a database using the requested access mode.
    ///
    /// For ordinary filesystem paths, both non-creating modes require an
    /// existing file. The `:memory:` and empty filenames are passed through.
    /// SQLite's busy timeout
    /// is zero, matching the source: a lock conflict is reported immediately.
    ///
    /// # Errors
    /// As [`Self::new`]. A failed open leaves no owned database handle behind.
    pub fn open(filename: impl AsRef<Path>, mode: SqlOpenMode) -> Result<Self> {
        let flags = match mode {
            SqlOpenMode::ReadOnly => OpenFlags::SQLITE_OPEN_READ_ONLY,
            SqlOpenMode::ReadWrite => OpenFlags::SQLITE_OPEN_READ_WRITE,
            SqlOpenMode::ReadWriteOrCreate => {
                OpenFlags::SQLITE_OPEN_READ_WRITE | OpenFlags::SQLITE_OPEN_CREATE
            }
        };
        let connection = Connection::open_with_flags(filename, flags).map_err(sql_error)?;
        connection.busy_timeout(Duration::ZERO).map_err(sql_error)?;
        Ok(Self { connection })
    }

    /// Test whether the main database contains a table with this exact name.
    ///
    /// As in the source, views and temporary tables are excluded, and the name
    /// comparison is case-sensitive. The name is a bound value, not SQL syntax.
    ///
    /// # Errors
    /// Returns [`Error::Io`] on a SQLite error; a failed query is not `false`.
    pub fn table_exists(&self, table_name: &str) -> Result<bool> {
        self.connection
            .query_row(
                "SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type='table' AND name=?1)",
                [table_name],
                |row| row.get(0),
            )
            .map_err(sql_error)
    }

    /// Test whether SQLite reports a column with this exact name for the table.
    ///
    /// The source documents an existing-table precondition, but its
    /// `PRAGMA table_info` returns no rows for a missing table; this also returns
    /// `false`. Views are accepted by the pragma. Column-name comparison is
    /// case-sensitive; table lookup follows SQLite's identifier rules.
    ///
    /// # Errors
    /// Returns [`Error::InvalidValue`] for a NUL in the table identifier, or
    /// [`Error::Io`] for an underlying SQLite/query-conversion error.
    pub fn column_exists(&self, table_name: &str, column_name: &str) -> Result<bool> {
        let sql = format!("PRAGMA table_info({})", identifier(table_name)?);
        let mut statement = self.connection.prepare(&sql).map_err(sql_error)?;
        let mut rows = statement.query([]).map_err(sql_error)?;
        while let Some(row) = rows.next().map_err(sql_error)? {
            if row.get_ref(1).map_err(sql_error)?
                == rusqlite::types::ValueRef::Text(column_name.as_bytes())
            {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Count rows in the table identified by its literal name.
    ///
    /// The name is quoted as one identifier; it is not a SQL expression or a
    /// `schema.table` path. This corrects the source's unquoted interpolation.
    ///
    /// # Errors
    /// Returns [`Error::Io`] if the table is absent or SQLite cannot complete
    /// the query. The source header promises `SqlOperationFailed` for an absent
    /// table, but its prepare helper actually throws `IllegalArgument`.
    /// [`Error::InvalidValue`] reports NUL identifiers or an unrepresentable count.
    pub fn count_table_rows(&self, table_name: &str) -> Result<usize> {
        let sql = format!("SELECT count(*) FROM {}", identifier(table_name)?);
        let count: i64 = self
            .connection
            .query_row(&sql, [], |row| row.get(0))
            .map_err(sql_error)?;
        usize::try_from(count)
            .map_err(|_| Error::InvalidValue("SQLite row count does not fit usize".into()))
    }

    /// Execute SQL statements in order, discarding any result rows.
    ///
    /// Like the source's `sqlite3_exec`, this accepts batches, empty SQL and
    /// explicit transaction commands. **It is not atomic:** an error leaves
    /// earlier effects in place; an open transaction stays open until the caller
    /// commits, rolls back, or drops the connection. Every result is stepped to
    /// completion so errors after the first row are reported.
    ///
    /// # Errors
    /// Returns [`Error::InvalidValue`] for embedded NUL, which the source
    /// truncates at; this is rejected before executing any prefix.
    /// Returns [`Error::Io`] for SQL, lock, constraint or other SQLite errors,
    /// replacing `Exception::IllegalArgument`. No automatic rollback occurs.
    pub fn execute_statement(&self, sql: &str) -> Result<()> {
        check_sql(sql)?;
        let mut batch = Batch::new(&self.connection, sql);
        while let Some(mut statement) = batch.next().map_err(sql_error)? {
            let mut rows = statement.raw_query();
            while rows.next().map_err(sql_error)?.is_some() {}
        }
        Ok(())
    }

    /// Execute one prepared statement, binding each byte slice as a BLOB.
    ///
    /// Values bind to indexes `1..=data.len()`. Embedded NUL and non-UTF-8 bytes
    /// survive unchanged, and an empty slice is a zero-length BLOB, not NULL.
    /// Unbound parameters retain SQLite's NULL default, as in the source.
    ///
    /// Only the first SQL statement is prepared and executed, matching the
    /// source; **trailing SQL is ignored**, even if malformed. Use
    /// [`Self::execute_statement`] for batches without bound BLOB values.
    /// Statements must finish
    /// without returning a row, matching the source's `SQLITE_DONE` requirement.
    /// A statement producing rows may already have effects before it is rejected.
    ///
    /// # Errors
    /// Returns [`Error::InvalidValue`] for empty/comment-only SQL or embedded
    /// NUL (checked in the whole input before execution, rather than the
    /// source's truncation). [`Error::Io`] reports preparation, excessive bindings,
    /// constraint errors, or a returned row. Statements are finalized on errors.
    pub fn execute_bind_statement(&self, sql: &str, data: &[&[u8]]) -> Result<()> {
        check_sql(sql)?;
        let mut batch = Batch::new(&self.connection, sql);
        let mut statement = batch
            .next()
            .map_err(sql_error)?
            .ok_or_else(|| Error::InvalidValue("SQLite binding requires a statement".into()))?;
        for (index, value) in data.iter().enumerate() {
            statement
                .raw_bind_parameter(index + 1, *value)
                .map_err(sql_error)?;
        }
        // Step once without rusqlite's optional pre-execution `extra_check`:
        // the source checks the SQLite result after executing the statement.
        if statement.raw_query().next().map_err(sql_error)?.is_some() {
            return Err(sql_error(rusqlite::Error::ExecuteReturnedResults));
        }
        Ok(())
    }
}

fn check_sql(sql: &str) -> Result<()> {
    if sql.contains('\0') {
        return Err(Error::InvalidValue("SQLite SQL contains NUL".into()));
    }
    Ok(())
}

fn identifier(name: &str) -> Result<String> {
    if name.contains('\0') {
        return Err(Error::InvalidValue("SQLite identifier contains NUL".into()));
    }
    Ok(format!("\"{}\"", name.replace('"', "\"\"")))
}

fn sql_error(error: rusqlite::Error) -> Error {
    Error::Io(std::io::Error::other(error))
}
