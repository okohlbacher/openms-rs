// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `OpenSwathParquetExporter`.
//!
//! Expected values come from:
//!
//! - the class test `OpenSwathParquetExporter_test.cpp` at core `bc9cc12`,
//!   transcribed;
//! - the executed Release oracle `../oracle/O5c-export` (ibminode06,
//!   `openms4-release-bc9cc12-c19e494-174b576`, GCC 14.4 `-O2 -DNDEBUG`,
//!   Arrow C++ 23.0.1, run twice with identical records and files):
//!   `tests/data/openswath_parquet_exporter/` holds its record stream and
//!   every file the Release exporter wrote. The Rust exporter runs on the
//!   same tables (the driver's inputs, mirrored below) and every column is
//!   compared with the Release file: name, type, nullability, values and
//!   nulls.

#![cfg(feature = "arrow")]

use openms::Error;
use openms::analysis::openswath::openswath_export_data::{
    OpenSwathFeatureScoreRow, OpenSwathFeatureScoreTable, OpenSwathTransitionScoreRow,
    OpenSwathTransitionScoreTable,
};
use openms::analysis::openswath::openswath_parquet_exporter::OpenSwathParquetExporter;
use openms::format::parquet_file::ParquetFile;
use openms::system::file::TempDir;

use arrow_schema::Schema;
use std::path::Path;

const DATA: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/openswath_parquet_exporter"
);

fn tmp() -> TempDir {
    TempDir::new_in(std::env::temp_dir(), false).unwrap()
}

fn records() -> Vec<Vec<String>> {
    std::fs::read_to_string(Path::new(DATA).join("records.tsv"))
        .unwrap()
        .lines()
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect()
}

fn outcome(all: &[Vec<String>], case: &str) -> Vec<String> {
    all.iter()
        .find(|r| (r[0] == "OK" || r[0] == "E") && r[1] == case)
        .unwrap_or_else(|| panic!("no outcome for {case}"))
        .clone()
}

/// Every column of `got` equals the Release file's: schema (name, type,
/// nullability), row count, values and null bitmaps.
fn assert_matches_release(case: &str, got: &Path) {
    let all = records();
    assert_eq!(outcome(&all, case)[0], "OK", "{case}");
    let file = all
        .iter()
        .find(|r| r[0] == "F" && r[1] == case)
        .map(|r| r[2].clone())
        .unwrap();
    let want = ParquetFile::read_table(Path::new(DATA).join("files").join(file)).unwrap();
    let g = ParquetFile::read_table(got).unwrap();
    let fields = |s: &Schema| {
        s.fields()
            .iter()
            .map(|f| (f.name().clone(), f.data_type().clone(), f.is_nullable()))
            .collect::<Vec<_>>()
    };
    assert_eq!(fields(&g.schema()), fields(&want.schema()), "{case} schema");
    assert_eq!(g.num_rows(), want.num_rows(), "{case} rows");
    for (i, field) in want.schema().fields().iter().enumerate() {
        assert_eq!(
            g.column(i).to_data(),
            want.column(i).to_data(),
            "{case} column {}",
            field.name()
        );
    }
}

// ------------------------------------------------------------------ inputs

/// The class test's feature rows (also the oracle's `fs_class`).
fn class_feature_table() -> OpenSwathFeatureScoreTable {
    let r = OpenSwathFeatureScoreRow {
        protein_id: 1,
        peptide_id: 2,
        precursor_id: 3,
        protein_accession: "PROT1".into(),
        unmodified_sequence: "PEPTIDEK".into(),
        modified_sequence: "PEPTIDEK(UniMod:1)".into(),
        precursor_mz: 450.7,
        precursor_charge: 2,
        run_id: 99,
        feature_id: 7777,
        exp_rt: 120.5,
        norm_rt: 118.0,
        delta_rt: 2.5,
        left_width: 118.0,
        right_width: 123.0,
        ..OpenSwathFeatureScoreRow::default()
    };
    let mut r2 = r.clone();
    r2.feature_id = 8888;
    r2.modified_sequence = "ACDEFGHIK".into();
    OpenSwathFeatureScoreTable {
        rows: vec![r, r2],
        ..OpenSwathFeatureScoreTable::default()
    }
}

/// The class test's transition rows (also the oracle's `ts_class`).
fn class_transition_table() -> OpenSwathTransitionScoreTable {
    let r = OpenSwathTransitionScoreRow {
        precursor_id: 5,
        transition_id: 42,
        transition_traml_id: "tr42".into(),
        product_mz: 500.25,
        transition_charge: 1,
        transition_type: "y".into(),
        transition_ordinal: 3,
        annotation: "y3".into(),
        transition_detecting: true,
        transition_decoy: false,
        ..OpenSwathTransitionScoreRow::default()
    };
    let mut r2 = r.clone();
    r2.transition_id = 43;
    r2.annotation = "y4".into();
    r2.transition_decoy = true;
    OpenSwathTransitionScoreTable {
        rows: vec![r, r2],
        ..OpenSwathTransitionScoreTable::default()
    }
}

/// The driver's `fs_rich`: every optional member, dynamic columns, value
/// vectors shorter and longer than the column lists, NaN and infinity.
fn rich_feature_table() -> OpenSwathFeatureScoreTable {
    let nan = f64::NAN;
    let a = OpenSwathFeatureScoreRow {
        protein_id: -5,
        peptide_id: 0,
        ipf_peptide_id: Some(12),
        precursor_id: 1 << 40,
        protein_accession: "P1;P2".into(),
        precursor_traml_id: "tr\t\"x\"".into(),
        precursor_group_label: "grp".into(),
        precursor_mz: 1e-300,
        precursor_charge: -3,
        precursor_library_intensity: Some(0.5),
        precursor_library_rt: Some(nan),
        precursor_library_drift_time: Some(-0.0),
        gene_id: Some(7),
        gene_name: "GENE".into(),
        gene_decoy: Some(true),
        protein_decoy: true,
        peptide_decoy: false,
        precursor_decoy: true,
        run_id: -1,
        filename: "run.mzML".into(),
        feature_id: 42,
        exp_rt: nan,
        exp_im: Some(1.25),
        norm_rt: f64::INFINITY,
        delta_rt: -2.5,
        left_width: 1.0,
        right_width: 2.0,
        im_left_width: Some(0.9),
        im_right_width: Some(1.1),
        feature_ms1_values: vec![1.0, 2.0],
        feature_ms2_values: vec![0.1, nan, 0.3],
        score_ms1_score: Some(1.5),
        score_ms1_rank: Some(1),
        score_ms1_pvalue: Some(0.01),
        score_ms1_qvalue: Some(0.02),
        score_ms1_pep: Some(0.03),
        score_ms2_score: Some(2.5),
        score_ms2_peak_group_rank: Some(-2),
        score_ms2_pvalue: Some(0.04),
        score_ms2_qvalue: Some(0.05),
        score_ms2_pep: Some(0.06),
        score_ipf_precursor_peakgroup_pep: Some(0.07),
        score_ipf_pep: Some(0.08),
        score_ipf_qvalue: Some(0.09),
        score_peptide_global_score: Some(1.0),
        score_peptide_global_pvalue: Some(2.0),
        score_peptide_global_qvalue: Some(3.0),
        score_peptide_global_pep: Some(4.0),
        score_peptide_experiment_wide_score: Some(5.0),
        score_peptide_experiment_wide_pvalue: Some(6.0),
        score_peptide_experiment_wide_qvalue: Some(7.0),
        score_peptide_experiment_wide_pep: Some(8.0),
        score_peptide_run_specific_score: Some(9.0),
        score_peptide_run_specific_pvalue: Some(10.0),
        score_peptide_run_specific_qvalue: Some(11.0),
        score_peptide_run_specific_pep: Some(12.0),
        score_protein_global_score: Some(13.0),
        score_protein_global_pvalue: Some(14.0),
        score_protein_global_qvalue: Some(15.0),
        score_protein_global_pep: Some(16.0),
        score_protein_experiment_wide_score: Some(17.0),
        score_protein_experiment_wide_pvalue: Some(18.0),
        score_protein_experiment_wide_qvalue: Some(19.0),
        score_protein_experiment_wide_pep: Some(20.0),
        score_protein_run_specific_score: Some(21.0),
        score_protein_run_specific_pvalue: Some(22.0),
        score_protein_run_specific_qvalue: Some(23.0),
        score_protein_run_specific_pep: Some(24.0),
        score_gene_global_score: Some(25.0),
        score_gene_global_pvalue: Some(26.0),
        score_gene_global_qvalue: Some(27.0),
        score_gene_global_pep: Some(28.0),
        score_gene_experiment_wide_score: Some(29.0),
        score_gene_experiment_wide_pvalue: Some(30.0),
        score_gene_experiment_wide_qvalue: Some(31.0),
        score_gene_experiment_wide_pep: Some(32.0),
        score_gene_run_specific_score: Some(33.0),
        score_gene_run_specific_pvalue: Some(34.0),
        score_gene_run_specific_qvalue: Some(35.0),
        score_gene_run_specific_pep: Some(36.0),
        alignment_group_id: Some(100),
        alignment_reference_feature_id: Some(200),
        alignment_reference_rt: Some(300.5),
        alignment_pep: Some(0.5),
        alignment_qvalue: Some(0.25),
        from_alignment: true,
        ..OpenSwathFeatureScoreRow::default()
    };
    let b = OpenSwathFeatureScoreRow {
        feature_ms2_values: vec![9.0],
        ..OpenSwathFeatureScoreRow::default()
    };
    let mut c = a.clone();
    c.gene_decoy = Some(false);
    c.feature_ms1_values = Vec::new();
    c.feature_ms2_values = vec![1.0, 2.0, 3.0, 4.0];
    c.protein_accession = String::new();
    c.score_ms1_rank = None;
    c.alignment_group_id = None;
    OpenSwathFeatureScoreTable {
        feature_ms1_column_names: vec!["VAR_MS1_A".into(), "VAR_MS1_B".into()],
        feature_ms2_column_names: vec!["VAR_XCORR".into(), "VAR_SN".into(), "MAIN".into()],
        rows: vec![a, b, c],
    }
}

/// The driver's `ts_rich`.
fn rich_transition_table() -> OpenSwathTransitionScoreTable {
    let a = OpenSwathTransitionScoreRow {
        run_id: Some(3),
        ipf_peptide_id: Some(4),
        precursor_id: 5,
        transition_id: 6,
        transition_traml_id: String::new(),
        product_mz: f64::NAN,
        transition_charge: -1,
        transition_type: String::new(),
        transition_ordinal: 70000,
        annotation: "b7^2".into(),
        transition_detecting: false,
        transition_library_intensity: Some(12.5),
        transition_decoy: true,
        feature_id: Some(99),
        feature_transition_values: vec![Some(1.5), None, Some(f64::NAN)],
        score_transition_score: Some(0.5),
        score_transition_rank: Some(2),
        score_transition_pvalue: Some(0.1),
        score_transition_qvalue: Some(0.2),
        score_transition_pep: Some(0.3),
    };
    let b = OpenSwathTransitionScoreRow {
        feature_transition_values: vec![Some(2.0)],
        ..OpenSwathTransitionScoreRow::default()
    };
    OpenSwathTransitionScoreTable {
        feature_transition_column_names: vec!["VAR_A".into(), "VAR_B".into(), "VAR_C".into()],
        rows: vec![a, b],
    }
}

// ------------------------------------------------------------------ class test

/// `OpenSwathParquetExporter_test.cpp`, section `writeFeatureScores`.
#[test]
fn class_test_write_feature_scores() {
    let tmp = tmp();
    let fname = tmp.path().join("f.parquet");
    OpenSwathParquetExporter::write_feature_scores(&fname, &class_feature_table()).unwrap();
    let tbl = ParquetFile::read_table(&fname).unwrap();
    assert_eq!(tbl.num_rows(), 2);
    let col = |name: &str| ParquetFile::get_column(&tbl, name).unwrap();
    let fid = col("FEATURE_ID");
    assert_eq!(
        ParquetFile::get_int64(Some(fid.as_ref()), 0, -1, false).unwrap(),
        7777
    );
    assert_eq!(
        ParquetFile::get_int64(Some(fid.as_ref()), 1, -1, false).unwrap(),
        8888
    );
    assert_eq!(
        ParquetFile::get_string(Some(col("PROTEIN_ACCESSION").as_ref()), 0).unwrap(),
        "PROT1"
    );
    assert_eq!(
        ParquetFile::get_string(Some(col("MODIFIED_SEQUENCE").as_ref()), 0).unwrap(),
        "PEPTIDEK(UniMod:1)"
    );
    assert_eq!(
        ParquetFile::get_string(Some(col("MODIFIED_SEQUENCE").as_ref()), 1).unwrap(),
        "ACDEFGHIK"
    );
    // TEST_REAL_SIMILAR; the values round-trip exactly.
    assert_eq!(
        ParquetFile::get_double(Some(col("PRECURSOR_MZ").as_ref()), 0, -1.0, false).unwrap(),
        450.7
    );
    assert_eq!(
        ParquetFile::get_double(Some(col("EXP_RT").as_ref()), 0, -1.0, false).unwrap(),
        120.5
    );

    // An empty table writes a schema-only, 0-row file.
    let fname2 = tmp.path().join("f2.parquet");
    OpenSwathParquetExporter::write_feature_scores(&fname2, &OpenSwathFeatureScoreTable::default())
        .unwrap();
    assert_eq!(ParquetFile::row_count(&fname2).unwrap(), 0);
}

/// `OpenSwathParquetExporter_test.cpp`, section `writeTransitionScores`.
#[test]
fn class_test_write_transition_scores() {
    let tmp = tmp();
    let fname = tmp.path().join("t.parquet");
    OpenSwathParquetExporter::write_transition_scores(&fname, &class_transition_table()).unwrap();
    let tbl = ParquetFile::read_table(&fname).unwrap();
    assert_eq!(tbl.num_rows(), 2);
    let col = |name: &str| ParquetFile::get_column(&tbl, name).unwrap();
    let tid = col("TRANSITION_ID");
    assert_eq!(
        ParquetFile::get_int64(Some(tid.as_ref()), 0, -1, false).unwrap(),
        42
    );
    assert_eq!(
        ParquetFile::get_int64(Some(tid.as_ref()), 1, -1, false).unwrap(),
        43
    );
    assert_eq!(
        ParquetFile::get_double(Some(col("PRODUCT_MZ").as_ref()), 0, -1.0, false).unwrap(),
        500.25
    );
    assert_eq!(
        ParquetFile::get_string(Some(col("ANNOTATION").as_ref()), 1).unwrap(),
        "y4"
    );
    let dec = col("TRANSITION_DECOY");
    assert!(!ParquetFile::get_bool(Some(dec.as_ref()), 0, true, false).unwrap());
    assert!(ParquetFile::get_bool(Some(dec.as_ref()), 1, false, false).unwrap());

    let fname2 = tmp.path().join("t2.parquet");
    OpenSwathParquetExporter::write_transition_scores(
        &fname2,
        &OpenSwathTransitionScoreTable::default(),
    )
    .unwrap();
    assert_eq!(ParquetFile::row_count(&fname2).unwrap(), 0);
}

// ------------------------------------------------------------------ oracle

#[test]
fn oracle_feature_scores() {
    let tmp = tmp();
    let cols_only = OpenSwathFeatureScoreTable {
        feature_ms1_column_names: vec!["A".into()],
        feature_ms2_column_names: vec!["B".into()],
        rows: Vec::new(),
    };
    for (case, table) in [
        ("fs_class", class_feature_table()),
        ("fs_empty", OpenSwathFeatureScoreTable::default()),
        ("fs_rich", rich_feature_table()),
        ("fs_empty_cols", cols_only),
    ] {
        let path = tmp.path().join(format!("{case}.parquet"));
        OpenSwathParquetExporter::write_feature_scores(&path, &table).unwrap();
        assert_matches_release(case, &path);
    }
}

#[test]
fn oracle_transition_scores() {
    let tmp = tmp();
    for (case, table) in [
        ("ts_class", class_transition_table()),
        ("ts_empty", OpenSwathTransitionScoreTable::default()),
        ("ts_rich", rich_transition_table()),
    ] {
        let path = tmp.path().join(format!("{case}.parquet"));
        OpenSwathParquetExporter::write_transition_scores(&path, &table).unwrap();
        assert_matches_release(case, &path);
    }
}

/// The Release build throws `Exception::FileNotWritable` for a path in a
/// missing directory; the port returns [`Error::Io`] and writes nothing.
#[test]
fn oracle_unwritable_destination() {
    let record = outcome(&records(), "err_unwritable");
    assert_eq!(record[0], "E");
    assert_eq!(record[2], "OpenMS::Exception::FileNotWritable");
    assert!(record[3].contains("is not writable for the current user"));
    let tmp = tmp();
    let path = tmp.path().join("no/such/dir/x.parquet");
    let err = OpenSwathParquetExporter::write_feature_scores(
        &path,
        &OpenSwathFeatureScoreTable::default(),
    )
    .unwrap_err();
    assert!(matches!(err, Error::Io(_)), "{err:?}");
    assert!(!path.exists());
    let err = OpenSwathParquetExporter::write_transition_scores(
        &path,
        &OpenSwathTransitionScoreTable::default(),
    )
    .unwrap_err();
    assert!(matches!(err, Error::Io(_)), "{err:?}");
}
