// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Source-reviewed SWATH class-test literals and retained mzML metadata projection.
//! Databases are built with rusqlite; no C++ execution or SqMassFile writer parity
//! is claimed by these tests.
#![cfg(feature = "sqlite")]

use openms::Error;
use openms::format::mzml_sqlite_swath_handler::{
    MzMLSqliteSwathHandler, SwathReadLimits, SwathWindow,
};
use openms::system::file::TempDir;
use rusqlite::Connection;
use std::path::PathBuf;

fn database() -> (TempDir, PathBuf, Connection) {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("swath.sqMass");
    let conn = Connection::open(&path).unwrap();
    conn.execute_batch("CREATE TABLE SPECTRUM(ID INT PRIMARY KEY NOT NULL, MSLEVEL INT); CREATE TABLE PRECURSOR(SPECTRUM_ID INT, CHROMATOGRAM_ID INT, ISOLATION_TARGET REAL, ISOLATION_LOWER REAL, ISOLATION_UPPER REAL)").unwrap();
    (dir, path, conn)
}

fn window(center: f64) -> SwathWindow {
    SwathWindow {
        center,
        lower: f64::NAN,
        upper: f64::INFINITY,
    }
}

fn source_fixture() -> (TempDir, PathBuf, Connection) {
    let (dir, path, conn) = database();
    // Full source-order ID/level/target projection from the hashed original
    // SwathFile.mzML. Every present lower/upper offset is exactly 12.5.
    // This bypasses SqMassFile writing and claims query-only source evidence.
    let rows = [
        (0, 1, None),
        (1, 2, Some(412.5)),
        (2, 2, Some(437.5)),
        (3, 2, Some(462.5)),
        (4, 2, Some(487.5)),
        (5, 2, Some(512.5)),
        (6, 1, None),
        (7, 2, Some(412.5)),
        (8, 2, Some(437.5)),
        (9, 2, Some(462.5)),
        (10, 2, Some(487.5)),
        (11, 2, Some(512.5)),
        (12, 1, None),
        (13, 2, Some(412.5)),
        (14, 2, Some(437.5)),
        (15, 2, Some(462.5)),
        (16, 2, Some(487.5)),
        (17, 2, Some(512.5)),
        (18, 1, None),
        (19, 2, Some(412.5)),
        (20, 2, Some(437.5)),
        (21, 2, Some(462.5)),
        (22, 2, Some(487.5)),
        (23, 2, Some(512.5)),
        (24, 1, None),
        (25, 2, Some(412.5)),
        (26, 2, Some(437.5)),
        (27, 2, Some(462.5)),
        (28, 2, Some(487.5)),
        (29, 2, Some(512.5)),
        (30, 1, None),
        (31, 2, Some(412.5)),
        (32, 2, Some(437.5)),
        (33, 2, Some(462.5)),
        (34, 2, Some(487.5)),
        (35, 2, Some(512.5)),
        (36, 1, None),
        (37, 2, Some(412.5)),
        (38, 2, Some(437.5)),
        (39, 2, Some(462.5)),
        (40, 2, Some(487.5)),
        (41, 2, Some(512.5)),
        (42, 1, None),
        (43, 2, Some(412.5)),
        (44, 2, Some(437.5)),
        (45, 2, Some(462.5)),
        (46, 2, Some(487.5)),
        (47, 2, Some(512.5)),
        (48, 1, None),
        (49, 2, Some(412.5)),
        (50, 2, Some(437.5)),
        (51, 2, Some(462.5)),
        (52, 2, Some(487.5)),
        (53, 2, Some(512.5)),
        (54, 1, None),
        (55, 2, Some(412.5)),
        (56, 2, Some(437.5)),
        (57, 2, Some(462.5)),
        (58, 2, Some(487.5)),
        (59, 2, Some(512.5)),
        (60, 1, None),
        (61, 2, Some(412.5)),
        (62, 2, Some(437.5)),
        (63, 2, Some(462.5)),
        (64, 2, Some(487.5)),
        (65, 2, Some(512.5)),
        (66, 1, None),
        (67, 2, Some(412.5)),
        (68, 2, Some(437.5)),
        (69, 2, Some(462.5)),
        (70, 2, Some(487.5)),
        (71, 2, Some(512.5)),
        (72, 1, None),
        (73, 2, Some(412.5)),
        (74, 2, Some(437.5)),
        (75, 2, Some(462.5)),
        (76, 2, Some(487.5)),
        (77, 2, Some(512.5)),
        (78, 1, None),
        (79, 2, Some(412.5)),
        (80, 2, Some(437.5)),
        (81, 2, Some(462.5)),
        (82, 2, Some(487.5)),
        (83, 2, Some(512.5)),
        (84, 1, None),
        (85, 2, Some(412.5)),
        (86, 2, Some(437.5)),
        (87, 2, Some(462.5)),
        (88, 2, Some(487.5)),
        (89, 2, Some(512.5)),
        (90, 1, None),
        (91, 2, Some(412.5)),
        (92, 2, Some(437.5)),
        (93, 2, Some(462.5)),
        (94, 2, Some(487.5)),
        (95, 2, Some(512.5)),
        (96, 1, None),
        (97, 2, Some(412.5)),
        (98, 2, Some(437.5)),
        (99, 2, Some(462.5)),
        (100, 2, Some(487.5)),
        (101, 2, Some(512.5)),
        (102, 1, None),
        (103, 2, Some(412.5)),
        (104, 2, Some(437.5)),
        (105, 2, Some(462.5)),
        (106, 2, Some(487.5)),
        (107, 2, Some(512.5)),
        (108, 1, None),
        (109, 2, Some(412.5)),
        (110, 2, Some(437.5)),
        (111, 2, Some(462.5)),
        (112, 2, Some(487.5)),
        (113, 2, Some(512.5)),
    ];
    assert_eq!(rows.len(), 114);
    for (id, level, center) in rows {
        conn.execute("INSERT INTO SPECTRUM VALUES (?1,?2)", (id, level))
            .unwrap();
        if let Some(center) = center {
            conn.execute(
                "INSERT INTO PRECURSOR VALUES (?1,NULL,?2,12.5,12.5)",
                (id, center),
            )
            .unwrap();
        }
    }
    (dir, path, conn)
}

#[test]
fn constructor_defers_io_and_source_defaults_map_to_value_subset() {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("missing.sqMass");
    let _handler = MzMLSqliteSwathHandler::new(&path);
    assert!(!path.exists());
    assert_eq!(
        SwathWindow::default(),
        SwathWindow {
            center: 0.0,
            lower: 0.0,
            upper: 0.0
        }
    );
    assert_eq!(SwathReadLimits::default().max_records, 1_000_000);
}

#[test]
fn destruction_and_queries_leave_no_connection_or_file_lock() {
    let (_dir, path, conn) = database();
    {
        let handler = MzMLSqliteSwathHandler::new(&path);
        assert!(handler.read_swath_windows().unwrap().is_empty());
        // All accessors release their transaction before returning, even while
        // the path-owning handler remains alive.
        conn.execute_batch("BEGIN EXCLUSIVE; COMMIT").unwrap();
    }
    drop(conn);
    std::fs::remove_file(path).unwrap();
}

#[test]
fn source_read_swath_windows_section() {
    let (_dir, path, _conn) = source_fixture();
    let maps = MzMLSqliteSwathHandler::new(path)
        .read_swath_windows()
        .unwrap();
    assert_eq!(maps.len(), 5);
    assert_eq!(
        maps[0],
        SwathWindow {
            center: 412.5,
            lower: 400.0,
            upper: 425.0
        }
    );
    assert_eq!((maps[1].lower, maps[1].upper), (425.0, 450.0));
    assert_eq!((maps[4].lower, maps[4].upper), (500.0, 525.0));
    // Source maps[0].ms1 == false has no field counterpart: this value type
    // represents discovered MS2 windows only (see support mapping).
}

#[test]
fn source_read_ms1_spectra_section() {
    let (_dir, path, _conn) = source_fixture();
    let ids = MzMLSqliteSwathHandler::new(path)
        .read_ms1_spectra()
        .unwrap();
    assert_eq!(ids.len(), 19);
    assert_eq!(ids[0], 0);
    assert_eq!(ids[18], 108);
    assert_eq!(ids, (0..19).map(|i| i * 6).collect::<Vec<_>>());
}

#[test]
fn source_read_spectra_for_window_section() {
    let (_dir, path, _conn) = source_fixture();
    let handler = MzMLSqliteSwathHandler::new(path);
    let maps = handler.read_swath_windows().unwrap();
    assert_eq!(maps.len(), 5);
    let first = handler.read_spectra_for_window(&maps[0]).unwrap();
    let second = handler.read_spectra_for_window(&maps[1]).unwrap();
    assert_eq!((first.len(), first[0], first[18]), (19, 1, 109));
    assert_eq!((second.len(), second[0], second[18]), (19, 2, 110));
}

#[test]
fn original_retained_sqmass_has_no_swath_windows_and_two_ms1_records() {
    // This is a retained C++-repository input, not a newly executed reference.
    let bytes = include_bytes!("data/sqlite_s1_source_review/SqliteMassFile_1.sqMass");
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("original.sqMass");
    std::fs::write(&path, bytes).unwrap();
    let handler = MzMLSqliteSwathHandler::new(path);
    assert!(handler.read_swath_windows().unwrap().is_empty());
    assert_eq!(handler.read_ms1_spectra().unwrap(), vec![0, 1]);
    // Its only precursor belongs to the TIC chromatogram (center 0).
    assert!(
        handler
            .read_spectra_for_window(&window(0.0))
            .unwrap()
            .is_empty()
    );
}

#[test]
fn discovery_deduplicates_tuples_and_sorts_independently_of_insertion_order() {
    let (_dir, path, conn) = database();
    conn.execute_batch("INSERT INTO SPECTRUM VALUES(8,2),(2,2),(1,1),(4,2); INSERT INTO PRECURSOR VALUES(8,NULL,437.5,12.5,12.5),(2,NULL,412.5,10,10),(2,NULL,412.5,12.5,12.5),(4,NULL,412.5,12.5,12.5),(1,NULL,99,1,1),(999,NULL,100,1,1)").unwrap();
    assert_eq!(
        MzMLSqliteSwathHandler::new(path)
            .read_swath_windows()
            .unwrap(),
        vec![
            SwathWindow {
                center: 412.5,
                lower: 400.0,
                upper: 425.0
            },
            SwathWindow {
                center: 412.5,
                lower: 402.5,
                upper: 422.5
            },
            SwathWindow {
                center: 437.5,
                lower: 425.0,
                upper: 450.0
            },
        ]
    );
}

#[test]
fn center_lookup_skips_null_rows_preserves_duplicates_and_ignores_other_fields() {
    let (_dir, path, conn) = database();
    conn.execute_batch("INSERT INTO SPECTRUM VALUES(9,1); INSERT INTO PRECURSOR VALUES(NULL,0,100,1,1),(9,NULL,100,1,1),(NULL,1,100,1,1),(4,NULL,100,1,1),(4,NULL,100,1,1),(7,NULL,NULL,1,1)").unwrap();
    assert_eq!(
        MzMLSqliteSwathHandler::new(path)
            .read_spectra_for_window(&window(100.0))
            .unwrap(),
        vec![4, 4, 9]
    );
}

#[test]
fn center_lookup_includes_both_exact_tolerance_boundaries() {
    let (_dir, path, conn) = database();
    let center = 100.0;
    for (id, target) in [
        (0, center - 0.01),
        (1, center + 0.01),
        (2, center - 0.01001),
        (3, center + 0.01001),
    ] {
        conn.execute("INSERT INTO PRECURSOR VALUES(?1,NULL,?2,0,0)", (id, target))
            .unwrap();
    }
    assert_eq!(
        MzMLSqliteSwathHandler::new(path)
            .read_spectra_for_window(&window(center))
            .unwrap(),
        vec![0, 1]
    );
}

#[test]
fn ids_are_checked_sorted_signed_64_bit_record_keys() {
    let (_dir, path, conn) = database();
    conn.execute_batch(
        "INSERT INTO SPECTRUM VALUES(4294967297,1),(7,1),(9223372036854775807,1),(0,NULL)",
    )
    .unwrap();
    let handler = MzMLSqliteSwathHandler::new(path);
    assert_eq!(
        handler.read_ms1_spectra().unwrap(),
        vec![7, 4294967297, i64::MAX]
    );
    conn.execute_batch("INSERT INTO SPECTRUM VALUES(-1,1)")
        .unwrap();
    assert!(matches!(
        handler.read_ms1_spectra(),
        Err(Error::InvalidValue(_))
    ));
    conn.execute_batch("INSERT INTO PRECURSOR VALUES(-1,NULL,10,0,0)")
        .unwrap();
    assert!(matches!(
        handler.read_spectra_for_window(&window(10.0)),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn nonfinite_input_and_selected_sql_values_are_rejected() {
    let (_dir, path, conn) = database();
    let handler = MzMLSqliteSwathHandler::new(path);
    for center in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            handler.read_spectra_for_window(&window(center)),
            Err(Error::InvalidValue(_))
        ));
    }
    conn.execute_batch("INSERT INTO SPECTRUM VALUES(0,2)")
        .unwrap();
    for (center, lo, hi) in [
        (f64::INFINITY, 0.0, 0.0),
        (0.0, f64::INFINITY, 0.0),
        (0.0, 0.0, f64::NEG_INFINITY),
        (f64::MAX, -f64::MAX, 0.0),
    ] {
        conn.execute_batch("DELETE FROM PRECURSOR").unwrap();
        conn.execute(
            "INSERT INTO PRECURSOR VALUES(0,NULL,?1,?2,?3)",
            (center, lo, hi),
        )
        .unwrap();
        assert!(matches!(
            handler.read_swath_windows(),
            Err(Error::InvalidValue(_))
        ));
    }
    conn.execute_batch("DELETE FROM PRECURSOR").unwrap();
    conn.execute(
        "INSERT INTO PRECURSOR VALUES(0,NULL,?1,0,0)",
        [f64::INFINITY],
    )
    .unwrap();
    assert!(matches!(
        handler.read_spectra_for_window(&window(10.0)),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn finite_negative_offsets_and_signed_zero_follow_source_arithmetic() {
    let (_dir, path, conn) = database();
    conn.execute_batch("INSERT INTO SPECTRUM VALUES(0,2); INSERT INTO PRECURSOR VALUES(0,NULL,10,-1,-2),(0,NULL,0,0,0),(0,NULL,-0.0,0,0)").unwrap();
    assert_eq!(
        MzMLSqliteSwathHandler::new(path)
            .read_swath_windows()
            .unwrap(),
        vec![
            SwathWindow::default(),
            SwathWindow {
                center: 10.0,
                lower: 11.0,
                upper: 8.0
            }
        ]
    );
}

#[test]
fn null_and_incorrect_sql_types_fail_without_returning_partial_results() {
    let (_dir, path, conn) = database();
    conn.execute_batch("INSERT INTO SPECTRUM VALUES(0,1),(1,1); INSERT INTO PRECURSOR VALUES(0,NULL,10,1,1),(1,NULL,'bad',1,1)").unwrap();
    let handler = MzMLSqliteSwathHandler::new(path);
    assert!(matches!(
        handler.read_spectra_for_window(&window(10.0)),
        Err(Error::Io(_))
    ));
    conn.execute_batch("UPDATE SPECTRUM SET MSLEVEL=2; UPDATE PRECURSOR SET ISOLATION_TARGET=10, ISOLATION_LOWER=NULL WHERE SPECTRUM_ID=1").unwrap();
    assert!(matches!(handler.read_swath_windows(), Err(Error::Io(_))));
    conn.execute_batch("UPDATE SPECTRUM SET MSLEVEL='bad' WHERE ID=1")
        .unwrap();
    assert!(matches!(handler.read_ms1_spectra(), Err(Error::Io(_))));
    conn.execute_batch("BEGIN EXCLUSIVE; COMMIT").unwrap();
}

#[test]
fn row_budget_counts_nonmatches_and_both_discovery_tables() {
    let (_dir, path, conn) = database();
    conn.execute_batch("INSERT INTO SPECTRUM VALUES(0,2),(1,1); INSERT INTO PRECURSOR VALUES(NULL,0,10,1,1),(0,NULL,10,1,1)").unwrap();
    let with = |n| MzMLSqliteSwathHandler::with_limits(&path, SwathReadLimits { max_records: n });
    assert!(matches!(
        with(1).read_ms1_spectra(),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(with(2).read_ms1_spectra().unwrap(), vec![1]);
    assert!(matches!(
        with(1).read_spectra_for_window(&window(10.0)),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(
        with(2).read_spectra_for_window(&window(10.0)).unwrap(),
        vec![0]
    );
    assert!(matches!(
        with(3).read_swath_windows(),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(
        with(4).read_swath_windows().unwrap(),
        vec![SwathWindow {
            center: 10.0,
            lower: 9.0,
            upper: 11.0
        }]
    );
    conn.execute_batch("DELETE FROM SPECTRUM; DELETE FROM PRECURSOR")
        .unwrap();
    assert!(with(0).read_ms1_spectra().unwrap().is_empty());
    assert!(with(0).read_swath_windows().unwrap().is_empty());
    assert!(
        with(0)
            .read_spectra_for_window(&window(10.0))
            .unwrap()
            .is_empty()
    );
}

#[test]
fn missing_paths_are_not_created_and_non_database_input_is_unchanged() {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("missing.sqMass");
    let handler = MzMLSqliteSwathHandler::new(&path);
    assert!(matches!(handler.read_ms1_spectra(), Err(Error::Io(_))));
    assert!(matches!(handler.read_swath_windows(), Err(Error::Io(_))));
    assert!(matches!(
        handler.read_spectra_for_window(&window(10.0)),
        Err(Error::Io(_))
    ));
    assert!(!path.exists());
    std::fs::write(&path, b"not a database").unwrap();
    assert!(matches!(handler.read_ms1_spectra(), Err(Error::Io(_))));
    assert_eq!(std::fs::read(&path).unwrap(), b"not a database");
}

#[test]
fn views_are_rejected_before_evaluation_and_missing_columns_report_errors() {
    let (_dir, path, conn) = database();
    conn.execute_batch("DROP TABLE SPECTRUM; CREATE VIEW SPECTRUM AS WITH RECURSIVE forever(x) AS (VALUES(0) UNION ALL SELECT x+1 FROM forever) SELECT x AS ID, 1 AS MSLEVEL FROM forever").unwrap();
    let handler = MzMLSqliteSwathHandler::new(path);
    assert!(matches!(
        handler.read_ms1_spectra(),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        handler.read_swath_windows(),
        Err(Error::InvalidValue(_))
    ));
    conn.execute_batch("DROP VIEW SPECTRUM; CREATE TABLE SPECTRUM(ID INT); DROP TABLE PRECURSOR; CREATE VIEW PRECURSOR AS SELECT 0 AS SPECTRUM_ID, 10 AS ISOLATION_TARGET").unwrap();
    assert!(matches!(handler.read_ms1_spectra(), Err(Error::Io(_))));
    assert!(matches!(
        handler.read_spectra_for_window(&window(10.0)),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn locked_reads_report_errors_and_recover_after_unlock() {
    let (_dir, path, conn) = database();
    let handler = MzMLSqliteSwathHandler::new(path);
    conn.execute_batch("BEGIN EXCLUSIVE").unwrap();
    assert!(matches!(handler.read_ms1_spectra(), Err(Error::Io(_))));
    assert!(matches!(handler.read_swath_windows(), Err(Error::Io(_))));
    assert!(matches!(
        handler.read_spectra_for_window(&window(10.0)),
        Err(Error::Io(_))
    ));
    conn.execute_batch("ROLLBACK").unwrap();
    assert!(handler.read_ms1_spectra().unwrap().is_empty());
}

#[test]
fn sql_step_error_after_a_valid_row_is_not_mistaken_for_end_of_results() {
    let (_dir, path, conn) = database();
    // Add the virtual generated column after inserting the source values.
    // An ordinary table can therefore fail during rows.next(): row 1 is
    // valid, then abs(i64::MIN) raises SQLITE_ERROR on the following row.
    conn.execute_batch("DROP TABLE SPECTRUM; CREATE TABLE SPECTRUM(ID INT); INSERT INTO SPECTRUM VALUES(1),(-9223372036854775808); ALTER TABLE SPECTRUM ADD COLUMN MSLEVEL INT GENERATED ALWAYS AS (abs(ID)) VIRTUAL").unwrap();
    let handler = MzMLSqliteSwathHandler::new(path);
    let error = handler.read_ms1_spectra().unwrap_err();
    let Error::Io(error) = error else {
        panic!("expected SQLite step error");
    };
    assert!(matches!(
        error.get_ref().unwrap().downcast_ref::<rusqlite::Error>(),
        Some(rusqlite::Error::SqliteFailure(_, _))
    ));
    conn.execute_batch("BEGIN EXCLUSIVE; ROLLBACK").unwrap();
}
