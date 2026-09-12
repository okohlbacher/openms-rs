// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Public SqliteConnector class-test literals (source review, eight sections),
//! plus native SQL ownership, name, BLOB, error and transaction regressions.
//! These execute Rust/rusqlite; they do not execute the C++ SDK.
#![cfg(feature = "sqlite")]

use openms::Error;
use openms::format::sqlite_connector::{SqlOpenMode, SqliteConnector};
use openms::system::file::TempDir;
use std::path::PathBuf;

fn database() -> (TempDir, PathBuf) {
    let directory = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = directory.path().join("connector.sqlite");
    (directory, path)
}

fn io_error(error: Error) {
    match error {
        Error::Io(error) => assert!(error.get_ref().unwrap().is::<rusqlite::Error>()),
        other => panic!("expected preserved SQLite error, got {other:?}"),
    }
}

#[test]
fn construction_creates_a_usable_database() {
    let (_directory, path) = database();
    let connector = SqliteConnector::new(&path).unwrap();
    connector
        .execute_statement("CREATE TABLE IF NOT EXISTS liveness_check (ID INT)")
        .unwrap();
    assert!(connector.table_exists("liveness_check").unwrap());
    assert_eq!(SqlOpenMode::default(), SqlOpenMode::ReadWriteOrCreate);
}

#[test]
fn destruction_closes_the_connection() {
    let (_directory, path) = database();
    {
        let connector = SqliteConnector::new(&path).unwrap();
        connector
            .execute_statement("CREATE TABLE T (ID INT); BEGIN; INSERT INTO T VALUES(1)")
            .unwrap();
    }
    // Native strengthening of the source's destructor-only section.
    let connector = SqliteConnector::open(&path, SqlOpenMode::ReadWrite).unwrap();
    assert!(connector.table_exists("T").unwrap());
    assert_eq!(connector.count_table_rows("T").unwrap(), 0);
    // A leaked connection with the uncommitted INSERT could still permit a
    // reader to see zero rows, but would retain the writer lock.
    connector
        .execute_statement("INSERT INTO T VALUES(2)")
        .unwrap();
    assert_eq!(connector.count_table_rows("T").unwrap(), 1);
    drop(connector);
    std::fs::remove_file(path).unwrap();
}

#[test]
fn execute_statement_inserts_and_reports_malformed_sql() {
    let connector = SqliteConnector::new(":memory:").unwrap();
    connector
        .execute_statement("CREATE TABLE T (ID INT, NAME TEXT)")
        .unwrap();
    connector
        .execute_statement("INSERT INTO T (ID, NAME) VALUES (1, 'a')")
        .unwrap();
    assert!(connector.table_exists("T").unwrap());
    io_error(connector.execute_statement("THIS IS NOT SQL").unwrap_err());
}

#[test]
fn table_exists_matches_exact_main_table_names() {
    let connector = SqliteConnector::new(":memory:").unwrap();
    connector
        .execute_statement("CREATE TABLE T (ID INT, NAME TEXT)")
        .unwrap();
    assert!(connector.table_exists("T").unwrap());
    assert!(!connector.table_exists("DOES_NOT_EXIST").unwrap());
    // sqlite_master's exact name comparison excludes case aliases, views and
    // the separate temporary schema, matching the source query.
    assert!(!connector.table_exists("t").unwrap());
    connector
        .execute_statement("CREATE VIEW V AS SELECT * FROM T; CREATE TEMP TABLE TEMP_T(ID)")
        .unwrap();
    assert!(!connector.table_exists("V").unwrap());
    assert!(!connector.table_exists("TEMP_T").unwrap());
}

#[test]
fn column_exists_matches_exact_column_names() {
    let connector = SqliteConnector::new(":memory:").unwrap();
    connector
        .execute_statement("CREATE TABLE T (ID INT, NAME TEXT)")
        .unwrap();
    assert!(connector.column_exists("T", "ID").unwrap());
    assert!(connector.column_exists("T", "NAME").unwrap());
    assert!(!connector.column_exists("T", "MISSING").unwrap());
    assert!(!connector.column_exists("T", "name").unwrap());
    assert!(!connector.column_exists("MISSING", "ID").unwrap());
    assert!(connector.column_exists("t", "ID").unwrap());
    connector
        .execute_statement("CREATE VIEW V AS SELECT ID FROM T")
        .unwrap();
    assert!(connector.column_exists("V", "ID").unwrap());
}

#[test]
fn external_non_utf8_column_names_do_not_hide_valid_columns() {
    let (_directory, path) = database();
    let external = rusqlite::Connection::open(&path).unwrap();
    external
        .execute_batch("CREATE TABLE T(placeholder BLOB, ordinary INT); PRAGMA writable_schema=ON")
        .unwrap();
    // Build a schema that cannot be expressed as Rust SQL text, as an external
    // producer can. SQLite accepts byte-valued identifiers in its UTF-8 schema.
    external
        .execute(
            "UPDATE sqlite_master SET sql=CAST(?1 AS TEXT) WHERE name='T'",
            [b"CREATE TABLE T(\"\xff\" BLOB, ordinary INT)".as_slice()],
        )
        .unwrap();
    drop(external);
    let connector = SqliteConnector::open(&path, SqlOpenMode::ReadOnly).unwrap();
    assert!(connector.column_exists("T", "ordinary").unwrap());
    assert!(!connector.column_exists("T", "missing").unwrap());
}

#[test]
fn count_table_rows_preserves_the_source_zero_and_two_counts() {
    let connector = SqliteConnector::new(":memory:").unwrap();
    connector
        .execute_statement("CREATE TABLE T (ID INT, NAME TEXT)")
        .unwrap();
    assert_eq!(connector.count_table_rows("T").unwrap(), 0);
    connector
        .execute_statement(
            "INSERT INTO T (ID, NAME) VALUES (1, 'a'); INSERT INTO T (ID, NAME) VALUES (2, 'b')",
        )
        .unwrap();
    assert_eq!(connector.count_table_rows("T").unwrap(), 2);
    io_error(connector.count_table_rows("UNKNOWN").unwrap_err());
}

#[test]
fn execute_bind_statement_inserts_the_source_bound_value() {
    let (_directory, path) = database();
    let connector = SqliteConnector::new(&path).unwrap();
    connector
        .execute_statement("CREATE TABLE T (ID INT, NAME TEXT)")
        .unwrap();
    connector
        .execute_bind_statement("INSERT INTO T (ID, NAME) VALUES (1, ?1)", &[b"bound_value"])
        .unwrap();
    assert_eq!(connector.count_table_rows("T").unwrap(), 1);
    let observer = rusqlite::Connection::open(&path).unwrap();
    let stored: Vec<u8> = observer
        .query_row("SELECT NAME FROM T", [], |row| row.get(0))
        .unwrap();
    assert_eq!(stored, b"bound_value");
}

#[test]
fn open_modes_preserve_read_only_and_missing_file_behavior() {
    let (_directory, path) = database();
    for mode in [SqlOpenMode::ReadOnly, SqlOpenMode::ReadWrite] {
        io_error(SqliteConnector::open(&path, mode).unwrap_err());
        assert!(!path.exists());
    }
    {
        let connector = SqliteConnector::new(&path).unwrap();
        connector
            .execute_statement("CREATE TABLE RO(ID INT); INSERT INTO RO VALUES(1),(2),(3)")
            .unwrap();
    }
    let reader = SqliteConnector::open(&path, SqlOpenMode::ReadOnly).unwrap();
    assert!(reader.table_exists("RO").unwrap());
    assert_eq!(reader.count_table_rows("RO").unwrap(), 3);
    io_error(
        reader
            .execute_statement("INSERT INTO RO VALUES(4)")
            .unwrap_err(),
    );
    drop(reader);
    let writer = SqliteConnector::open(&path, SqlOpenMode::ReadWrite).unwrap();
    writer
        .execute_statement("INSERT INTO RO VALUES(4)")
        .unwrap();
    assert_eq!(writer.count_table_rows("RO").unwrap(), 4);
}

#[test]
fn blobs_preserve_empty_nul_and_non_utf8_bytes_and_missing_bindings() {
    let (_directory, path) = database();
    let connector = SqliteConnector::new(&path).unwrap();
    connector
        .execute_statement("CREATE TABLE T (ID INT, DATA BLOB, OPTIONAL BLOB)")
        .unwrap();
    connector
        .execute_bind_statement(
            "INSERT INTO T VALUES (1, ?1, ?2), (2, ?3, ?4)",
            &[&[], &[0, 255, 128, 39], b"text"],
        )
        .unwrap();
    let observer = rusqlite::Connection::open(&path).unwrap();
    let first: (String, Vec<u8>, Vec<u8>) = observer
        .query_row(
            "SELECT typeof(DATA),DATA,OPTIONAL FROM T WHERE ID=1",
            [],
            |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)),
        )
        .unwrap();
    assert_eq!(first, ("blob".into(), Vec::new(), vec![0, 255, 128, 39]));
    let second: (Vec<u8>, Option<Vec<u8>>) = observer
        .query_row("SELECT DATA,OPTIONAL FROM T WHERE ID=2", [], |row| {
            Ok((row.get(0)?, row.get(1)?))
        })
        .unwrap();
    assert_eq!(second, (b"text".to_vec(), None));
}

#[test]
fn quoted_identifiers_are_literal_and_cannot_change_the_query() {
    let connector = SqliteConnector::new(":memory:").unwrap();
    connector.execute_statement("CREATE TABLE \"odd'\"\";--\" (\"space name\" BLOB); INSERT INTO \"odd'\"\";--\" VALUES(1)").unwrap();
    assert!(connector.table_exists("odd'\";--").unwrap());
    assert!(connector.column_exists("odd'\";--", "space name").unwrap());
    assert_eq!(connector.count_table_rows("odd'\";--").unwrap(), 1);
    assert!(!connector.table_exists("' OR 1=1 --").unwrap());
    io_error(
        connector
            .count_table_rows("sqlite_master; DROP TABLE sqlite_master")
            .unwrap_err(),
    );
    assert!(matches!(
        connector.count_table_rows("x\0y"),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        connector.column_exists("x\0y", "a"),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn batch_failures_preserve_prior_effects_and_caller_transactions() {
    let connector = SqliteConnector::new(":memory:").unwrap();
    io_error(
        connector
            .execute_statement("CREATE TABLE T(ID); INSERT INTO T VALUES(1); NOT SQL")
            .unwrap_err(),
    );
    assert_eq!(connector.count_table_rows("T").unwrap(), 1);
    io_error(
        connector
            .execute_statement("BEGIN; INSERT INTO T VALUES(2); NOT SQL")
            .unwrap_err(),
    );
    assert_eq!(connector.count_table_rows("T").unwrap(), 2);
    connector.execute_statement("ROLLBACK").unwrap();
    assert_eq!(connector.count_table_rows("T").unwrap(), 1);
    connector
        .execute_statement("BEGIN; INSERT INTO T VALUES(3); COMMIT")
        .unwrap();
    assert_eq!(connector.count_table_rows("T").unwrap(), 2);
    connector
        .execute_statement("  -- only a comment\n")
        .unwrap();
    connector.execute_statement("").unwrap();
}

#[test]
fn batch_execution_observes_errors_after_the_first_result_row() {
    let connector = SqliteConnector::new(":memory:").unwrap();
    connector
        .execute_statement("CREATE TABLE T(ID); INSERT INTO T VALUES(1),(2)")
        .unwrap();
    io_error(connector.execute_statement("SELECT CASE ID WHEN 2 THEN abs(-9223372036854775808) ELSE ID END FROM T ORDER BY rowid; INSERT INTO T VALUES(3)").unwrap_err());
    assert_eq!(connector.count_table_rows("T").unwrap(), 2);
}

#[test]
fn failed_bindings_and_constraints_do_not_poison_the_connection() {
    let connector = SqliteConnector::new(":memory:").unwrap();
    connector
        .execute_statement("CREATE TABLE T(DATA BLOB UNIQUE)")
        .unwrap();
    io_error(
        connector
            .execute_bind_statement("INSERT INTO T VALUES(?1)", &[b"one", b"two"])
            .unwrap_err(),
    );
    assert_eq!(connector.count_table_rows("T").unwrap(), 0);
    connector
        .execute_bind_statement("INSERT INTO T VALUES(?1)", &[b"one"])
        .unwrap();
    io_error(
        connector
            .execute_bind_statement("INSERT INTO T VALUES(?1)", &[b"one"])
            .unwrap_err(),
    );
    io_error(
        connector
            .execute_bind_statement("SELECT ?1", &[b"one"])
            .unwrap_err(),
    );
    connector
        .execute_bind_statement("INSERT INTO T VALUES(?1)", &[b"two"])
        .unwrap();
    assert_eq!(connector.count_table_rows("T").unwrap(), 2);
    for sql in ["", " -- comment\n", "/* comment */"] {
        assert!(matches!(
            connector.execute_bind_statement(sql, &[]),
            Err(Error::InvalidValue(_))
        ));
    }
}

#[test]
fn binding_steps_transactions_and_returning_statements_once() {
    let connector = SqliteConnector::new(":memory:").unwrap();
    connector
        .execute_statement("CREATE TABLE T(DATA BLOB)")
        .unwrap();
    connector.execute_bind_statement("BEGIN", &[]).unwrap();
    let error = connector
        .execute_bind_statement("INSERT INTO T VALUES(?1) RETURNING DATA", &[b"one"])
        .unwrap_err();
    let Error::Io(error) = error else {
        panic!("expected SQLite error")
    };
    assert!(matches!(
        error.get_ref().unwrap().downcast_ref::<rusqlite::Error>(),
        Some(rusqlite::Error::ExecuteReturnedResults)
    ));
    assert_eq!(connector.count_table_rows("T").unwrap(), 1);
    connector.execute_bind_statement("COMMIT", &[]).unwrap();
    // Even with rusqlite/extra_check unified, evaluation reaches SQLite and
    // reports the overflow, rather than preemptively rejecting a SELECT.
    let error = connector
        .execute_bind_statement("SELECT abs(-9223372036854775808)", &[])
        .unwrap_err();
    let Error::Io(error) = error else {
        panic!("expected SQLite error")
    };
    assert!(matches!(
        error.get_ref().unwrap().downcast_ref::<rusqlite::Error>(),
        Some(rusqlite::Error::SqliteFailure(_, _))
    ));
}

#[test]
fn binding_ignores_the_tail_without_preparing_its_pragma() {
    let connector = SqliteConnector::new(":memory:").unwrap();
    connector.execute_statement("PRAGMA foreign_keys=ON; CREATE TABLE P(ID BLOB PRIMARY KEY); CREATE TABLE C(PID BLOB REFERENCES P(ID))").unwrap();
    connector
        .execute_bind_statement("INSERT INTO P VALUES(?1); THIS IS NOT SQL", &[b"parent"])
        .unwrap();
    connector
        .execute_bind_statement(
            "INSERT INTO C VALUES(?1); PRAGMA foreign_keys=OFF",
            &[b"parent"],
        )
        .unwrap();
    io_error(
        connector
            .execute_bind_statement("INSERT INTO C VALUES(?1)", &[b"missing"])
            .unwrap_err(),
    );
    assert_eq!(connector.count_table_rows("C").unwrap(), 1);
}

#[test]
fn nul_sql_is_rejected_before_running_any_prefix() {
    let connector = SqliteConnector::new(":memory:").unwrap();
    connector
        .execute_statement("CREATE TABLE T(DATA BLOB)")
        .unwrap();
    assert!(matches!(
        connector.execute_statement("INSERT INTO T VALUES(1);\0 INSERT INTO T VALUES(2)"),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        connector.execute_bind_statement("INSERT INTO T VALUES(?1);\0", &[b"one"]),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(connector.count_table_rows("T").unwrap(), 0);
}

#[test]
fn locked_query_errors_are_not_reported_as_absence() {
    let (_directory, path) = database();
    let locker = SqliteConnector::new(&path).unwrap();
    locker
        .execute_statement("CREATE TABLE T(ID); INSERT INTO T VALUES(1),(2),(3)")
        .unwrap();
    let reader = SqliteConnector::open(&path, SqlOpenMode::ReadOnly).unwrap();
    // Warm the schema as in the independently retained C++ lock probe, so
    // cached schema permits preparation and the query encounters the lock.
    assert!(reader.table_exists("T").unwrap());
    assert!(reader.column_exists("T", "ID").unwrap());
    assert_eq!(reader.count_table_rows("T").unwrap(), 3);
    locker.execute_statement("BEGIN EXCLUSIVE").unwrap();
    io_error(reader.table_exists("T").unwrap_err());
    io_error(reader.column_exists("T", "ID").unwrap_err());
    io_error(reader.count_table_rows("T").unwrap_err());
    locker.execute_statement("ROLLBACK").unwrap();
    assert!(reader.table_exists("T").unwrap());
    assert!(reader.column_exists("T", "ID").unwrap());
    assert_eq!(reader.count_table_rows("T").unwrap(), 3);
}

#[test]
fn sqlite_temporary_filename_is_supported() {
    let connector = SqliteConnector::new("").unwrap();
    connector.execute_statement("CREATE TABLE T(ID)").unwrap();
    assert!(connector.table_exists("T").unwrap());
}
