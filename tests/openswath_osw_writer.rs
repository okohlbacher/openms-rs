// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OpenSwathOSWWriter: the SQLite output.
//!
//! The first three tests transcribe sections of the upstream class test
//! `OpenSwathOSWWriter_test.cpp` at the core pin `bc9cc12` (evidence tier 3);
//! the values read back are those the upstream test queries. The remaining
//! tests check the schema against the statements of the source and the
//! transaction behaviour of `writeRows` (tier 4). Row *preparation* is
//! tested without a database in `tests/openswath_export.rs`.
#![cfg(feature = "sqlite")]

use openms::Error;
use openms::analysis::openswath::openswath_osw_writer::{
    FEATURE_COLUMNS, FEATURE_MS1_COLUMNS, FEATURE_MS2_COLUMNS, FEATURE_PRECURSOR_COLUMNS,
    FEATURE_TRANSITION_COLUMNS, OpenSwathOSWWriter, OswData, OswValue,
};
use openms::analysis::openswath::transition_experiment::LightCompound;
use openms::kernel::{Feature, FeatureMap};
use openms::metadata::MetaValue;
use rusqlite::Connection;
use std::path::PathBuf;

fn temp(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms_openswath_osw_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(name);
    let _ = std::fs::remove_file(&path);
    let _ = std::fs::remove_file(path.with_extension("osw-wal"));
    let _ = std::fs::remove_file(path.with_extension("osw-shm"));
    path
}

fn real(value: f64) -> MetaValue {
    MetaValue::try_from(value).unwrap()
}

fn query_i64(path: &PathBuf, sql: &str) -> i64 {
    let connection = Connection::open(path).unwrap();
    connection.query_row(sql, [], |row| row.get(0)).unwrap()
}

#[test]
fn run_id_is_stored_as_integer_not_blob() {
    // upstream regression: RUN.ID was a BLOB, so FEATURE.RUN_ID = RUN.ID
    // joins matched nothing
    let path = temp("run_id.osw");
    let large_run_id: u64 = 6_130_996_817_540_441_879;
    {
        let mut writer = OpenSwathOSWWriter::new(path.to_str().unwrap(), false);
        writer.write_header().unwrap();
        writer
            .add_run(large_run_id, "test_input_file.mzML")
            .unwrap();
        let rid = large_run_id & !(1_u64 << 63);
        let mut rows = OswData::default();
        rows.feature_rows.push([
            OswValue::parse(&(rid + 1).to_string()),
            OswValue::parse(&rid.to_string()),
            OswValue::parse("999"),
            OswValue::parse("100.0"),
            OswValue::parse("NULL"),
            OswValue::parse("100.0"),
            OswValue::parse("0.0"),
            OswValue::parse("90.0"),
            OswValue::parse("110.0"),
            OswValue::parse("NULL"),
            OswValue::parse("NULL"),
        ]);
        writer.write_rows(&rows).unwrap();
    }
    assert_eq!(
        query_i64(
            &path,
            "SELECT count(*) FROM RUN WHERE typeof(ID) = 'integer'"
        ),
        1
    );
    assert_eq!(
        query_i64(
            &path,
            "SELECT count(*) FROM FEATURE WHERE typeof(RUN_ID) = 'integer'"
        ),
        1
    );
    assert_eq!(
        query_i64(
            &path,
            "SELECT count(*) FROM RUN INNER JOIN FEATURE ON FEATURE.RUN_ID = RUN.ID"
        ),
        1
    );
    let rid = large_run_id & !(1_u64 << 63);
    assert_eq!(
        query_i64(&path, &format!("SELECT count(*) FROM RUN WHERE ID = {rid}")),
        1
    );
}

fn osw_feature() -> Feature {
    let mut feature = Feature::new(10.0, 0.0, 100.0);
    feature.unique_id = 0x1234;
    feature.metadata.insert("leftWidth".into(), real(9.0));
    feature.metadata.insert("rightWidth".into(), real(11.0));
    feature.metadata.insert("total_xic".into(), real(1000.0));
    feature
        .metadata
        .insert("peak_apices_sum".into(), real(100.0));
    feature
}

#[test]
fn prepare_line_and_write_lines_keep_masserror_by_ms2_position() {
    let path = temp("masserror.osw");
    {
        let mut writer = OpenSwathOSWWriter::new(path.to_str().unwrap(), false);
        writer.write_header().unwrap();
        writer.add_run(1, "test_input_file.mzML").unwrap();

        let mut feature = osw_feature();
        feature.metadata.insert(
            "masserror_ppm".into(),
            MetaValue::try_from(vec![12.5]).unwrap(),
        );
        let mut ms1 = Feature::new(0.0, 0.0, 10.0);
        ms1.metadata
            .insert("FeatureLevel".into(), MetaValue::from("MS1"));
        ms1.metadata
            .insert("native_id".into(), MetaValue::from("Precursor_i0"));
        ms1.metadata.insert("peak_apex_int".into(), real(5.0));
        let mut ms2 = Feature::new(0.0, 0.0, 50.0);
        ms2.metadata
            .insert("FeatureLevel".into(), MetaValue::from("MS2"));
        ms2.metadata
            .insert("native_id".into(), MetaValue::from(123_i64));
        ms2.metadata.insert("total_xic".into(), real(500.0));
        ms2.metadata.insert("peak_apex_int".into(), real(25.0));
        ms2.metadata.insert("peak_apex_position".into(), real(10.5));
        ms2.metadata.insert("width_at_50".into(), real(2.0));
        feature.subordinates.push(ms1);
        feature.subordinates.push(ms2);
        let map = FeatureMap::from_features(vec![feature]);
        let line = writer
            .prepare_line(&LightCompound::default(), None, &map, "77")
            .unwrap();
        writer.write_lines(&[line]).unwrap();
    }
    assert_eq!(
        query_i64(
            &path,
            "SELECT count(*) FROM FEATURE_TRANSITION WHERE TRANSITION_ID = 123 AND MASSERROR_PPM = 12.5"
        ),
        1
    );
    // the other tables received their rows too
    assert_eq!(query_i64(&path, "SELECT count(*) FROM FEATURE"), 1);
    assert_eq!(query_i64(&path, "SELECT count(*) FROM FEATURE_MS2"), 1);
    assert_eq!(
        query_i64(&path, "SELECT count(*) FROM FEATURE_PRECURSOR"),
        1
    );
    assert_eq!(query_i64(&path, "SELECT count(*) FROM FEATURE_MS1"), 0);
}

#[test]
fn uis_scoring_uses_the_available_transition_names() {
    let path = temp("uis.osw");
    {
        let mut writer = OpenSwathOSWWriter::new(path.to_str().unwrap(), true);
        writer.write_header().unwrap();
        writer.add_run(1, "test_input_file.mzML").unwrap();
        let mut feature = osw_feature();
        feature
            .metadata
            .insert("id_target_num_transitions".into(), MetaValue::from(2_i64));
        feature.metadata.insert(
            "id_target_transition_names".into(),
            MetaValue::from(vec![201_i64]),
        );
        for (name, value) in [
            ("id_target_area_intensity", 10.0),
            ("id_target_total_area_intensity", 20.0),
            ("id_target_apex_intensity", 30.0),
            ("id_target_peak_apex_position", 10.5),
            ("id_target_width_at_50", 2.0),
        ] {
            feature
                .metadata
                .insert(name.into(), MetaValue::try_from(vec![value]).unwrap());
        }
        let map = FeatureMap::from_features(vec![feature]);
        let line = writer
            .prepare_line(&LightCompound::default(), None, &map, "77")
            .unwrap();
        writer.write_lines(&[line]).unwrap();
    }
    assert_eq!(
        query_i64(&path, "SELECT count(*) FROM FEATURE_TRANSITION"),
        1
    );
    assert_eq!(
        query_i64(
            &path,
            "SELECT count(*) FROM FEATURE_TRANSITION WHERE TRANSITION_ID = 201 AND AREA_INTENSITY = 10.0"
        ),
        1
    );
}

fn column_names(path: &PathBuf, table: &str) -> Vec<String> {
    let connection = Connection::open(path).unwrap();
    let mut statement = connection
        .prepare(&format!("PRAGMA table_info({table})"))
        .unwrap();
    statement
        .query_map([], |row| row.get::<_, String>(1))
        .unwrap()
        .map(Result::unwrap)
        .collect()
}

#[test]
fn write_header_creates_the_source_schema() {
    let path = temp("schema.osw");
    let writer = OpenSwathOSWWriter::new(path.to_str().unwrap(), false);
    writer.write_header().unwrap();
    assert_eq!(column_names(&path, "RUN"), ["ID", "FILENAME"]);
    assert_eq!(
        column_names(&path, "FEATURE"),
        [
            "ID",
            "RUN_ID",
            "PRECURSOR_ID",
            "EXP_RT",
            "EXP_IM",
            "NORM_RT",
            "DELTA_RT",
            "LEFT_WIDTH",
            "RIGHT_WIDTH",
            "EXP_IM_LEFTWIDTH",
            "EXP_IM_RIGHTWIDTH",
        ]
    );
    // the column lists of the insert statements are columns of the tables
    // (the source declares FEATURE_MS1 with VAR_MI_SCORE before
    // VAR_IM_MS1_DELTA_SCORE and inserts the other way round, by name)
    for (table, columns) in [
        ("FEATURE", FEATURE_COLUMNS.to_vec()),
        ("FEATURE_MS1", FEATURE_MS1_COLUMNS.to_vec()),
        ("FEATURE_MS2", FEATURE_MS2_COLUMNS.to_vec()),
        ("FEATURE_PRECURSOR", FEATURE_PRECURSOR_COLUMNS.to_vec()),
        ("FEATURE_TRANSITION", FEATURE_TRANSITION_COLUMNS.to_vec()),
    ] {
        let declared = column_names(&path, table);
        assert_eq!(declared.len(), columns.len(), "{table}");
        let mut a = declared.clone();
        let mut b: Vec<String> = columns.iter().map(|c| c.to_string()).collect();
        a.sort();
        b.sort();
        assert_eq!(a, b, "{table}");
    }
    // creating the header twice fails: the tables exist
    assert!(matches!(writer.write_header(), Err(Error::Io(_))));
}

#[test]
fn write_rows_is_one_transaction() {
    let path = temp("transaction.osw");
    let mut writer = OpenSwathOSWWriter::new(path.to_str().unwrap(), false);
    writer.write_header().unwrap();
    writer.add_run(7, "a.mzML").unwrap();
    assert_eq!(writer.run_id(), 7);
    let feature_row = |id: i64| -> [OswValue; 11] {
        [
            OswValue::Int64(id),
            OswValue::Int64(7),
            OswValue::Int64(1),
            OswValue::Double(10.0),
            OswValue::Null,
            OswValue::Double(10.0),
            OswValue::Double(0.0),
            OswValue::Double(9.0),
            OswValue::Double(11.0),
            OswValue::Null,
            OswValue::Null,
        ]
    };
    let mut ok = OswData::default();
    ok.feature_rows.push(feature_row(1));
    ok.feature_rows.push(feature_row(2));
    writer.write_rows(&ok).unwrap();
    // an empty batch is a no-op
    writer.write_rows(&OswData::default()).unwrap();

    // a duplicate primary key rolls the whole batch back, the earlier rows of
    // the batch included
    let mut bad = OswData::default();
    bad.feature_rows.push(feature_row(3));
    bad.feature_rows.push(feature_row(1));
    assert!(matches!(writer.write_rows(&bad), Err(Error::Io(_))));
    drop(writer);
    assert_eq!(query_i64(&path, "SELECT count(*) FROM FEATURE"), 2);
    assert_eq!(
        query_i64(&path, "SELECT count(*) FROM FEATURE WHERE ID = 3"),
        0
    );
    // the write connection switched the file to WAL, as the source does
    let connection = Connection::open(&path).unwrap();
    let mode: String = connection
        .query_row("PRAGMA journal_mode", [], |row| row.get(0))
        .unwrap();
    assert_eq!(mode.to_lowercase(), "wal");
}

#[test]
fn nan_and_text_cells_are_bound_as_null_and_text() {
    let path = temp("cells.osw");
    let writer = OpenSwathOSWWriter::new(path.to_str().unwrap(), false);
    writer.write_header().unwrap();
    let mut rows = OswData::default();
    rows.feature_precursor_rows.push([
        OswValue::Int64(1),
        OswValue::Int64(0),
        OswValue::Double(f64::NAN), // NOT NULL column: refused
        OswValue::Double(1.0),
    ]);
    assert!(writer.write_rows(&rows).is_err());
    let mut rows = OswData::default();
    rows.feature_rows.push([
        OswValue::Int64(1),
        OswValue::Int64(1),
        OswValue::Text("PEPTIDE_2".into()),
        OswValue::Double(1.5),
        OswValue::Double(f64::NAN),
        OswValue::Double(1.5),
        OswValue::Double(0.0),
        OswValue::Double(1.0),
        OswValue::Double(2.0),
        OswValue::Null,
        OswValue::Null,
    ]);
    writer.write_rows(&rows).unwrap();
    drop(writer);
    assert_eq!(
        query_i64(
            &path,
            "SELECT count(*) FROM FEATURE WHERE typeof(PRECURSOR_ID) = 'text' AND EXP_IM IS NULL"
        ),
        1
    );
}

#[test]
fn inactive_writer_ignores_runs() {
    let mut writer = OpenSwathOSWWriter::new("", false);
    writer.add_run(5, "ignored.mzML").unwrap();
    assert_eq!(writer.run_id(), 0);
}
