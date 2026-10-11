// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OpenSwathExportData, OpenSwathExportConfig, OpenSwathResultsExporter,
//! OpenSwathMatrixExporter, and the row preparation of OpenSwathOSWWriter.
//!
//! Expected values are transcribed from the upstream class tests at the core
//! pin `bc9cc12` (evidence tier 3): `OpenSwathResultsExporter_test.cpp`,
//! `OpenSwathMatrixExporter_test.cpp`, the matrix section of
//! `OpenSwathExport_test.cpp` and the `prepareLine` sections of
//! `OpenSwathOSWWriter_test.cpp`. The cases marked "oracle" compare with
//! output of the executed C++ Release build, retained under
//! `tests/data/openswath_export/` (evidence tier 1; see
//! `tests/data/openswath_export_provenance.json`). Cases marked "native" check
//! behaviour the upstream tests do not reach.

use openms::Error;
use openms::analysis::openswath::openswath_export_config::{
    OpenSwathExportFilterConfig, OpenSwathMatrixExportConfig, OpenSwathParquetExportConfig,
    OpenSwathResultsExportConfig,
};
use openms::analysis::openswath::openswath_export_data::{
    OpenSwathExportFileFormat, OpenSwathExportRow, OpenSwathIPFExportMode, OpenSwathMatrixLevel,
    OpenSwathMatrixNormalization, OpenSwathQuantMatrix,
};
use openms::analysis::openswath::openswath_matrix_exporter::OpenSwathMatrixExporter;
use openms::analysis::openswath::openswath_osw_writer::{
    FEATURE_COLUMNS, FEATURE_MS1_COLUMNS, FEATURE_MS2_COLUMNS, FEATURE_PRECURSOR_COLUMNS,
    FEATURE_TRANSITION_COLUMNS, OpenSwathOSWWriter, OswData, OswValue,
};
use openms::analysis::openswath::openswath_results_exporter::{
    ExportColumn, OpenSwathResultsExporter, RESULTS_TSV_COLUMNS,
};
use openms::analysis::openswath::transition_experiment::LightCompound;
use openms::kernel::{Feature, FeatureMap};
use openms::metadata::MetaValue;
use std::path::PathBuf;

fn temp(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms_openswath_export_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

fn real(value: f64) -> MetaValue {
    MetaValue::try_from(value).unwrap()
}

fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

// --- enums and configuration ---------------------------------------------

#[test]
fn enum_names_are_the_source_to_string_values() {
    assert_eq!(OpenSwathExportFileFormat::Tsv.as_str(), "tsv");
    assert_eq!(OpenSwathExportFileFormat::Parquet.as_str(), "parquet");
    assert_eq!(OpenSwathIPFExportMode::Peptidoform.as_str(), "peptidoform");
    assert_eq!(OpenSwathIPFExportMode::Augmented.as_str(), "augmented");
    assert_eq!(OpenSwathIPFExportMode::Disable.as_str(), "disable");
    assert_eq!(OpenSwathMatrixLevel::Precursor.as_str(), "precursor");
    assert_eq!(OpenSwathMatrixLevel::Peptide.as_str(), "peptide");
    assert_eq!(OpenSwathMatrixLevel::Protein.as_str(), "protein");
    assert_eq!(OpenSwathMatrixLevel::Gene.as_str(), "gene");
    assert_eq!(OpenSwathMatrixNormalization::None.as_str(), "none");
    assert_eq!(OpenSwathMatrixNormalization::Median.as_str(), "median");
    assert_eq!(
        OpenSwathMatrixNormalization::MedianMedian.to_string(),
        "medianmedian"
    );
}

#[test]
fn configuration_defaults_are_those_of_the_header() {
    let f = OpenSwathExportFilterConfig::default();
    assert_eq!(f.ipf_mode, OpenSwathIPFExportMode::Peptidoform);
    assert!(f.transition_quantification);
    assert_eq!(f.max_transition_pep, 0.7);
    assert_eq!(f.ipf_max_peptidoform_pep, 0.4);
    assert_eq!(f.max_rs_peakgroup_qvalue, 0.05);
    assert!(f.peptide);
    assert_eq!(f.max_global_peptide_qvalue, 0.01);
    assert!(f.protein);
    assert_eq!(f.max_global_protein_qvalue, 0.01);
    assert!(!f.gene);
    assert_eq!(f.max_global_gene_qvalue, 0.01);
    assert!(!f.use_alignment);
    assert_eq!(f.max_alignment_pep, 0.7);
    assert!(f.exclude_decoys);

    let r = OpenSwathResultsExportConfig::default();
    assert_eq!(r.format, OpenSwathExportFileFormat::Tsv);
    assert_eq!(r.filters, f);

    let m = OpenSwathMatrixExportConfig::default();
    assert_eq!(m.level, OpenSwathMatrixLevel::Peptide);
    assert_eq!(m.normalization, OpenSwathMatrixNormalization::None);
    assert_eq!(m.top_n, 3);
    assert!(m.consistent_top);
    assert_eq!(m.format, OpenSwathExportFileFormat::Tsv);

    let p = OpenSwathParquetExportConfig::default();
    assert!(p.include_transition_data);
    assert_eq!(p.filters, f);
}

#[test]
fn row_defaults_use_minus_one_ids() {
    let row = OpenSwathExportRow::default();
    assert_eq!(row.run_id, -1);
    assert_eq!(row.feature_id, -1);
    assert_eq!(row.peptide_id, -1);
    assert_eq!(row.precursor_id, -1);
    assert_eq!(row.charge, 0);
    assert!(!row.decoy && !row.from_alignment);
    assert!(row.aggr_prec_peak_area.is_none());
}

// --- OpenSwathResultsExporter --------------------------------------------

fn results_row(feature_id: i64, tg: &str, seq: &str) -> OpenSwathExportRow {
    OpenSwathExportRow {
        run_id: 99,
        filename: "run.mzML".into(),
        run_name: "runA".into(),
        feature_id,
        transition_group_id: tg.into(),
        sequence: seq.into(),
        full_peptide_name: seq.into(),
        protein_name: "PROT1".into(),
        charge: 2,
        mz: 450.7,
        rt: 120.5,
        intensity: 1234.5,
        ..OpenSwathExportRow::default()
    }
}

#[test]
fn results_exporter_writes_the_tsv_of_the_upstream_test() {
    let mut rows = vec![
        results_row(7777, "TG1", "PEPTIDEK"),
        results_row(8888, "TG2", "SAMPLER"),
    ];
    rows[0].aggr_prec_peak_area = Some(321.5);
    rows[0].alignment_group_id = Some(1_234_567_890_123);
    rows[0].gene_name = "GENE1".into();
    rows[0].decoy = true;
    rows[0].from_alignment = true;

    let path = temp("results.tsv");
    OpenSwathResultsExporter::write(&path, &rows, &OpenSwathResultsExportConfig::default())
        .unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let r0: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let r1: Vec<&str> = lines.next().unwrap().split('\t').collect();
    assert!(lines.next().is_none());

    assert_eq!(header.len(), 59);
    assert_eq!(header, RESULTS_TSV_COLUMNS);
    assert_eq!(header[0], "run_id");
    assert_eq!(header[3], "feature_id");
    assert_eq!(header[5], "transition_group_id");
    assert_eq!(header[8], "Sequence");
    assert_eq!(header[14], "RT");
    assert_eq!(r0.len(), 59);
    assert_eq!(r0[0], "99");
    assert_eq!(r0[3], "7777");
    assert_eq!(r0[5], "TG1");
    assert_eq!(r0[8], "PEPTIDEK");
    assert_eq!(r0[14], "120.5");
    assert_eq!(r1[3], "8888");
    assert_eq!(r1[8], "SAMPLER");

    // the remaining set fields land where the header says
    let column = |name: &str| header.iter().position(|h| *h == name).unwrap();
    assert_eq!(r0[column("decoy")], "1");
    assert_eq!(r1[column("decoy")], "0");
    assert_eq!(r0[column("GeneName")], "GENE1");
    assert_eq!(r0[column("aggr_prec_Peak_Area")], "321.5");
    assert_eq!(r1[column("aggr_prec_Peak_Area")], "");
    assert_eq!(r0[column("alignment_group_id")], "1234567890123");
    assert_eq!(r0[column("from_alignment")], "1");
    assert_eq!(r0[column("Intensity")], "1234.5");
    assert_eq!(r0[column("mz")], "450.7");
}

#[test]
fn results_exporter_empty_export_is_a_header_only_tsv() {
    // OpenSwathExport_test: the header carries Sequence, no data line follows
    let path = temp("empty_results.tsv");
    OpenSwathResultsExporter::write(&path, &[], &OpenSwathResultsExportConfig::default()).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    assert_eq!(text.lines().count(), 1);
    assert!(text.lines().next().unwrap().contains("Sequence"));
    assert!(text.lines().next().unwrap().contains("m_score_gene_global"));
    assert!(text.ends_with('\n'));
}

#[test]
fn results_exporter_number_formats_follow_the_stream_state() {
    // native: doubles at %g precision 6 (os << double), optionals through
    // StringUtils::toStr (full precision)
    let mut row = results_row(1, "TG", "PEP");
    row.mz = 1_234_567.0;
    row.rt = 0.000_012_345_678;
    row.intensity = 100_000.0;
    row.d_score = -2.5;
    row.m_score = 1e-5;
    row.pep = Some(0.25);
    row.ms2_pep = Some(1e-5);
    row.exp_im = Some(36739.0);
    let text = OpenSwathResultsExporter::tsv_string(&[row]).unwrap();
    let fields: Vec<&str> = text.lines().nth(1).unwrap().split('\t').collect();
    let col = |n: &str| RESULTS_TSV_COLUMNS.iter().position(|h| *h == n).unwrap();
    assert_eq!(fields[col("mz")], "1.23457e+06");
    assert_eq!(fields[col("RT")], "1.23457e-05");
    assert_eq!(fields[col("Intensity")], "100000");
    assert_eq!(fields[col("d_score")], "-2.5");
    assert_eq!(fields[col("m_score")], "1e-05");
    assert_eq!(fields[col("pep")], "0.25");
    assert_eq!(fields[col("ms2_pep")], "1.0e-05");
    assert_eq!(fields[col("EXP_IM")], "3.6739e04");
}

#[test]
fn results_exporter_parquet_schema_is_59_typed_columns() {
    // OpenSwathResultsExporter_test: an empty Parquet export keeps the schema
    let empty = OpenSwathResultsExporter::table(&[]);
    assert_eq!(empty.num_rows(), 0);
    assert_eq!(empty.num_columns(), 59);
    assert_eq!(empty.columns[0].0, "run_id");
    assert_eq!(empty.columns[58].0, "ipf_peptidoform_m_score");

    let mut rows = vec![
        results_row(7777, "TG1", "PEPTIDEK"),
        results_row(8888, "TG2", "SAMPLER"),
    ];
    rows[0].aggr_prec_peak_area = Some(321.5);
    rows[0].alignment_group_id = Some(1_234_567_890_123);
    rows[0].gene_name = "GENE1".into();
    rows[0].decoy = true;
    rows[0].from_alignment = true;
    let table = OpenSwathResultsExporter::table(&rows);
    assert_eq!(table.num_rows(), 2);
    assert_eq!(
        table.column("feature_id"),
        Some(&ExportColumn::Int64(vec![Some(7777), Some(8888)]))
    );
    assert_eq!(
        table.column("Sequence"),
        Some(&ExportColumn::Utf8(vec![
            Some("PEPTIDEK".into()),
            Some("SAMPLER".into())
        ]))
    );
    assert_eq!(
        table.column("transition_group_id"),
        Some(&ExportColumn::Utf8(vec![
            Some("TG1".into()),
            Some("TG2".into())
        ]))
    );
    assert_eq!(
        table.column("aggr_prec_Peak_Area"),
        Some(&ExportColumn::Float64(vec![Some(321.5), None]))
    );
    assert_eq!(
        table.column("alignment_group_id"),
        Some(&ExportColumn::Int64(vec![Some(1_234_567_890_123), None]))
    );
    // an empty string is a null, as appendString_ does
    assert_eq!(
        table.column("GeneName"),
        Some(&ExportColumn::Utf8(vec![Some("GENE1".into()), None]))
    );
    assert_eq!(
        table.column("decoy"),
        Some(&ExportColumn::Boolean(vec![true, false]))
    );
    assert_eq!(
        table.column("from_alignment"),
        Some(&ExportColumn::Boolean(vec![true, false]))
    );
    // Charge is int64 although the row field is 32 bit
    assert_eq!(
        table.column("Charge"),
        Some(&ExportColumn::Int64(vec![Some(2), Some(2)]))
    );
}

#[test]
fn results_exporter_parquet_is_refused_until_parquet_file_is_ported() {
    let path = temp("results.parquet");
    let config = OpenSwathResultsExportConfig {
        format: OpenSwathExportFileFormat::Parquet,
        ..OpenSwathResultsExportConfig::default()
    };
    let error = OpenSwathResultsExporter::write(&path, &[], &config).unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)));
    assert!(!path.exists(), "a refused export must not create the file");
}

// --- OpenSwathMatrixExporter ---------------------------------------------

#[allow(clippy::too_many_arguments)]
fn matrix_row(
    run_id: i64,
    run_name: &str,
    tg: &str,
    seq: &str,
    prot: &str,
    charge: i32,
    intensity: f64,
    m_score: f64,
    feature_id: i64,
) -> OpenSwathExportRow {
    OpenSwathExportRow {
        run_id,
        run_name: run_name.into(),
        transition_group_id: tg.into(),
        sequence: seq.into(),
        full_peptide_name: seq.into(),
        protein_name: prot.into(),
        charge,
        intensity,
        m_score,
        feature_id,
        ..OpenSwathExportRow::default()
    }
}

fn precursor_config() -> OpenSwathMatrixExportConfig {
    OpenSwathMatrixExportConfig {
        level: OpenSwathMatrixLevel::Precursor,
        normalization: OpenSwathMatrixNormalization::None,
        ..OpenSwathMatrixExportConfig::default()
    }
}

#[test]
fn matrix_precursor_level_of_the_upstream_test() {
    // TG1 is seen in both runs; in runA it has two peak groups and the lower
    // m_score one (not the more intense one) wins. TG2 is only seen in runA.
    let rows = vec![
        matrix_row(1, "runA", "TG1", "PEPTIDEK", "PROT1", 2, 1000.0, 0.01, 11),
        matrix_row(1, "runA", "TG1", "PEPTIDEK", "PROT1", 2, 9999.0, 0.05, 12),
        matrix_row(2, "runB", "TG1", "PEPTIDEK", "PROT1", 2, 2000.0, 0.02, 21),
        matrix_row(1, "runA", "TG2", "SAMPLER", "PROT2", 3, 300.0, 0.03, 31),
    ];
    let m = OpenSwathMatrixExporter::build_matrix(&rows, &precursor_config()).unwrap();

    assert_eq!(m.identifier_column_names.len(), 5);
    assert_eq!(m.identifier_column_names[0], "transition_group_id");
    assert_eq!(m.identifier_column_names[4], "ProteinName");
    assert_eq!(m.sample_column_names, ["runA", "runB"]);
    assert_eq!(m.identifier_rows.len(), 2);
    assert_eq!(m.identifier_rows[0][0], "TG1");
    assert_eq!(m.identifier_rows[0][1], "PEPTIDEK");
    assert_eq!(m.identifier_rows[0][3], "2");
    assert_eq!(m.identifier_rows[0][4], "PROT1");
    assert_eq!(m.identifier_rows[1][0], "TG2");
    assert_eq!(m.identifier_rows[1][4], "PROT2");
    real_similar(m.values[0][0].unwrap(), 1000.0);
    real_similar(m.values[0][1].unwrap(), 2000.0);
    real_similar(m.values[1][0].unwrap(), 300.0);
    assert!(m.values[1][1].is_none());
}

#[test]
fn matrix_preconditions() {
    let config = precursor_config();
    assert!(matches!(
        OpenSwathMatrixExporter::build_matrix(&[], &config),
        Err(Error::InvalidValue(_))
    ));
    let rows = vec![matrix_row(
        1, "runA", "TG1", "PEPTIDEK", "PROT1", 2, 1000.0, 0.01, 11,
    )];
    let bad = OpenSwathMatrixExportConfig {
        level: OpenSwathMatrixLevel::Peptide,
        top_n: 0,
        ..OpenSwathMatrixExportConfig::default()
    };
    assert!(matches!(
        OpenSwathMatrixExporter::build_matrix(&rows, &bad),
        Err(Error::InvalidValue(_))
    ));
    // top_n 0 is fine at precursor level
    let ok = OpenSwathMatrixExportConfig {
        top_n: 0,
        ..precursor_config()
    };
    assert!(OpenSwathMatrixExporter::build_matrix(&rows, &ok).is_ok());
    // no protein annotation at protein level, no gene annotation at gene level
    let protein = OpenSwathMatrixExportConfig {
        level: OpenSwathMatrixLevel::Protein,
        ..OpenSwathMatrixExportConfig::default()
    };
    let unannotated = vec![matrix_row(
        1, "runA", "TG1", "PEPTIDEK", "", 2, 1000.0, 0.01, 11,
    )];
    assert!(OpenSwathMatrixExporter::build_matrix(&unannotated, &protein).is_err());
    let gene = OpenSwathMatrixExportConfig {
        level: OpenSwathMatrixLevel::Gene,
        ..OpenSwathMatrixExportConfig::default()
    };
    assert!(OpenSwathMatrixExporter::build_matrix(&rows, &gene).is_err());
}

#[test]
fn matrix_write_matrix_of_the_upstream_test() {
    let rows = vec![
        matrix_row(1, "runA", "TG1", "PEPTIDEK", "PROT1", 2, 1000.0, 0.01, 11),
        matrix_row(2, "runB", "TG1", "PEPTIDEK", "PROT1", 2, 2000.0, 0.02, 21),
    ];
    let config = precursor_config();
    let m = OpenSwathMatrixExporter::build_matrix(&rows, &config).unwrap();
    let path = temp("matrix.tsv");
    OpenSwathMatrixExporter::write_matrix(&path, &m, &config).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let mut lines = text.lines();
    let h: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let d: Vec<&str> = lines.next().unwrap().split('\t').collect();
    assert_eq!(h.len(), 7);
    assert_eq!(h[0], "transition_group_id");
    assert_eq!(h[5], "runA");
    assert_eq!(h[6], "runB");
    assert_eq!(d.len(), 7);
    assert_eq!(d[0], "TG1");
    assert_eq!(d[1], "PEPTIDEK");
    assert_eq!(d[3], "2");
    assert_eq!(d[4], "PROT1");
    assert_eq!(d[5], "1000");
    assert_eq!(d[6], "2000");

    // an empty matrix is rejected
    let error = OpenSwathMatrixExporter::write_matrix(
        temp("matrix_empty.tsv"),
        &OpenSwathQuantMatrix::default(),
        &config,
    );
    assert!(matches!(error, Err(Error::InvalidValue(_))));
}

#[test]
fn matrix_writer_rejects_malformed_shapes() {
    // OpenSwathExport_test: three malformed matrices are Preconditions and
    // nothing is written
    let path = temp("ignored.tsv");
    let _ = std::fs::remove_file(&path);
    let config = OpenSwathMatrixExportConfig::default();
    let mut matrix = OpenSwathQuantMatrix {
        identifier_column_names: vec!["Sequence".into()],
        sample_column_names: vec!["run1".into()],
        identifier_rows: vec![vec!["PEPTIDE".into()]],
        values: vec![],
    };
    assert!(matches!(
        OpenSwathMatrixExporter::write_matrix(&path, &matrix, &config),
        Err(Error::InvalidValue(_))
    ));
    matrix.values = vec![vec![Some(1.0)]];
    matrix.identifier_rows = vec![vec!["PEPTIDE".into(), "EXTRA".into()]];
    assert!(OpenSwathMatrixExporter::write_matrix(&path, &matrix, &config).is_err());
    matrix.identifier_rows = vec![vec!["PEPTIDE".into()]];
    matrix.values = vec![vec![Some(1.0), Some(2.0)]];
    assert!(OpenSwathMatrixExporter::write_matrix(&path, &matrix, &config).is_err());
    assert!(!path.exists());
    // the same checks guard the typed columns of the Parquet output
    assert!(OpenSwathMatrixExporter::table(&matrix).is_err());
}

#[test]
fn matrix_sample_names_are_made_unique() {
    // native: an empty run name becomes RUN_ID_<id>; a repeated label gets
    // _RUN<id> appended (the label of the run met first keeps its name)
    let rows = vec![
        matrix_row(5, "same", "TG", "PEP", "P", 2, 10.0, 0.01, 1),
        matrix_row(3, "same", "TG", "PEP", "P", 2, 20.0, 0.01, 2),
        matrix_row(9, "", "TG", "PEP", "P", 2, 30.0, 0.01, 3),
    ];
    let m = OpenSwathMatrixExporter::build_matrix(&rows, &precursor_config()).unwrap();
    // columns are ordered by run id: 3, 5, 9
    assert_eq!(m.sample_column_names, ["same_RUN3", "same", "RUN_ID_9"]);
    assert_eq!(m.values, [[Some(20.0), Some(10.0), Some(30.0)]]);
}

#[test]
fn matrix_peptide_protein_gene_summaries() {
    // native, hand-derived: two precursors of peptide A (charges 2 and 3) in
    // two runs, one of peptide B; proteins P1 (A, B) and P2 (C)
    let mut rows = Vec::new();
    let mut add =
        |run: i64, tg: &str, seq: &str, prot: &str, gene: &str, intensity: f64, id: i64| {
            let mut r = matrix_row(
                run,
                &format!("run{run}"),
                tg,
                seq,
                prot,
                2,
                intensity,
                0.01,
                id,
            );
            r.gene_name = gene.into();
            rows.push(r);
        };
    add(1, "A_2", "A", "P1", "G1", 100.0, 1);
    add(1, "A_3", "A", "P1", "G1", 300.0, 2);
    add(2, "A_2", "A", "P1", "G1", 200.0, 3);
    add(2, "A_3", "A", "P1", "G1", 400.0, 4);
    add(1, "B_2", "B", "P1", "G1", 50.0, 5);
    add(2, "B_2", "B", "P1", "G1", 150.0, 6);
    add(1, "C_2", "C", "P2", "G2", 10.0, 7);

    // peptide level, top 3 precursors consistent over runs: the mean of the
    // precursors of the peptide per run
    let peptide = OpenSwathMatrixExportConfig {
        level: OpenSwathMatrixLevel::Peptide,
        ..OpenSwathMatrixExportConfig::default()
    };
    let m = OpenSwathMatrixExporter::build_matrix(&rows, &peptide).unwrap();
    assert_eq!(m.identifier_column_names, ["Sequence", "FullPeptideName"]);
    assert_eq!(m.sample_column_names, ["run1", "run2"]);
    assert_eq!(m.identifier_rows[0], ["A", "A"]);
    assert_eq!(m.values[0], [Some(200.0), Some(300.0)]); // (100+300)/2, (200+400)/2
    assert_eq!(m.values[1], [Some(50.0), Some(150.0)]);
    assert_eq!(m.values[2], [Some(10.0), None]);

    // top 1 consistent: the precursor of highest median (A_3: 350) only
    let top1 = OpenSwathMatrixExportConfig {
        top_n: 1,
        ..peptide.clone()
    };
    let m = OpenSwathMatrixExporter::build_matrix(&rows, &top1).unwrap();
    assert_eq!(m.values[0], [Some(300.0), Some(400.0)]);

    // protein level: mean over the peptide summaries per run; P1 has A and B
    let protein = OpenSwathMatrixExportConfig {
        level: OpenSwathMatrixLevel::Protein,
        ..peptide.clone()
    };
    let m = OpenSwathMatrixExporter::build_matrix(&rows, &protein).unwrap();
    assert_eq!(m.identifier_column_names, ["ProteinName"]);
    assert_eq!(m.identifier_rows, [["P1"], ["P2"]]);
    assert_eq!(
        m.values[0],
        [Some((200.0 + 50.0) / 2.0), Some((300.0 + 150.0) / 2.0)]
    );
    assert_eq!(m.values[1], [Some(10.0), None]);

    // gene level
    let gene = OpenSwathMatrixExportConfig {
        level: OpenSwathMatrixLevel::Gene,
        ..peptide
    };
    let m = OpenSwathMatrixExporter::build_matrix(&rows, &gene).unwrap();
    assert_eq!(m.identifier_column_names, ["GeneName"]);
    assert_eq!(m.identifier_rows, [["G1"], ["G2"]]);
}

#[test]
fn matrix_normalization() {
    // native, hand-derived: precursor level, three precursors in two runs;
    // run1 values 10, 20, 30 (median 20), run2 values 2, 4, 6 (median 4)
    let mut rows = Vec::new();
    for (i, (a, b)) in [(10.0, 2.0), (20.0, 4.0), (30.0, 6.0)].iter().enumerate() {
        let tg = format!("TG{i}");
        rows.push(matrix_row(
            1,
            "r1",
            &tg,
            "P",
            "X",
            2,
            *a,
            0.01,
            10 + i as i64,
        ));
        rows.push(matrix_row(
            2,
            "r2",
            &tg,
            "P",
            "X",
            2,
            *b,
            0.01,
            20 + i as i64,
        ));
    }
    let median = OpenSwathMatrixExportConfig {
        normalization: OpenSwathMatrixNormalization::Median,
        ..precursor_config()
    };
    let m = OpenSwathMatrixExporter::build_matrix(&rows, &median).unwrap();
    assert_eq!(m.values[0], [Some(0.5), Some(0.5)]);
    assert_eq!(m.values[1], [Some(1.0), Some(1.0)]);
    assert_eq!(m.values[2], [Some(1.5), Some(1.5)]);
    // median of medians: (20 + 4) / 2 = 12
    let medmed = OpenSwathMatrixExportConfig {
        normalization: OpenSwathMatrixNormalization::MedianMedian,
        ..precursor_config()
    };
    let m = OpenSwathMatrixExporter::build_matrix(&rows, &medmed).unwrap();
    assert_eq!(m.values[0], [Some(6.0), Some(6.0)]);
    assert_eq!(m.values[1], [Some(12.0), Some(12.0)]);
    assert_eq!(m.values[2], [Some(18.0), Some(18.0)]);
}

// --- OpenSwathOSWWriter row preparation ----------------------------------

#[test]
fn osw_value_classification_follows_strtoll_and_strtod() {
    let cases: &[(&str, OswValue)] = &[
        ("", OswValue::Null),
        ("NULL", OswValue::Null),
        ("nan", OswValue::Null),
        ("-nan", OswValue::Null),
        ("NaN", OswValue::Null),
        ("-NaN", OswValue::Null),
        ("NAN", OswValue::Null),
        ("0", OswValue::Int64(0)),
        ("-17", OswValue::Int64(-17)),
        ("+5", OswValue::Int64(5)),
        (" 12 ", OswValue::Int64(12)),
        (
            "6130996817540441880",
            OswValue::Int64(6_130_996_817_540_441_880),
        ),
        ("100.0", OswValue::Double(100.0)),
        ("1e3", OswValue::Double(1000.0)),
        (".5", OswValue::Double(0.5)),
        ("99999999999999999999", OswValue::Double(1e20)),
        ("inf", OswValue::Double(f64::INFINITY)),
        ("1e999", OswValue::Text("1e999".into())), // ERANGE
        ("12abc", OswValue::Text("12abc".into())),
        ("NIAGGLEK_2", OswValue::Text("NIAGGLEK_2".into())),
        ("1 2", OswValue::Text("1 2".into())),
        ("Precursor_i0", OswValue::Text("Precursor_i0".into())),
    ];
    for (text, expected) in cases {
        assert_eq!(&OswValue::parse(text), expected, "text {text:?}");
    }
    // NaN as a real is a NULL, and `text` never parses
    assert_eq!(OswValue::real(f64::NAN), OswValue::Null);
    assert!(OswValue::Double(f64::NAN).is_null());
    assert_eq!(OswValue::text("5"), OswValue::Text("5".into()));
}

#[test]
fn osw_value_text_forms() {
    assert_eq!(OswValue::Null.to_text(), "NULL");
    assert_eq!(OswValue::Int64(-3).to_text(), "-3");
    assert_eq!(OswValue::Double(0.5).to_text(), "0.5");
    assert_eq!(OswValue::Double(36739.0).to_text(), "3.6739e04");
    assert_eq!(OswValue::Double(f64::NAN).to_text(), "NULL");
    assert_eq!(OswValue::Text("it's".into()).to_sql_literal(), "'it''s'");
    assert_eq!(OswValue::Int64(7).to_sql_literal(), "7");
    assert_eq!(OswValue::Null.to_sql_literal(), "NULL");
}

fn osw_feature() -> Feature {
    let mut feature = Feature::new(10.0, 0.0, 100.0);
    feature.unique_id = 0x8000_0000_0000_0123; // the sign bit is cleared
    feature.metadata.insert("leftWidth".into(), real(9.0));
    feature.metadata.insert("rightWidth".into(), real(11.0));
    feature.metadata.insert("total_xic".into(), real(1000.0));
    feature
        .metadata
        .insert("peak_apices_sum".into(), real(100.0));
    feature
}

fn subordinate(level: &str, intensity: f32, native_id: MetaValue) -> Feature {
    let mut sub = Feature::new(0.0, 0.0, intensity);
    sub.metadata
        .insert("FeatureLevel".into(), MetaValue::from(level));
    sub.metadata.insert("native_id".into(), native_id);
    sub
}

#[test]
fn prepare_rows_indexes_masserror_ppm_by_ms2_position() {
    // OpenSwathOSWWriter_test, "prepareLine indexes masserror_ppm by MS2
    // subordinate position": the MS1 subordinate does not consume an index
    let mut feature = osw_feature();
    feature.metadata.insert(
        "masserror_ppm".into(),
        MetaValue::try_from(vec![12.5]).unwrap(),
    );
    let mut ms1 = subordinate("MS1", 10.0, MetaValue::from("Precursor_i0"));
    ms1.metadata.insert("peak_apex_int".into(), real(5.0));
    let mut ms2 = subordinate("MS2", 50.0, MetaValue::from(123_i64));
    ms2.metadata.insert("total_xic".into(), real(500.0));
    ms2.metadata.insert("peak_apex_int".into(), real(25.0));
    ms2.metadata.insert("peak_apex_position".into(), real(10.5));
    ms2.metadata.insert("width_at_50".into(), real(2.0));
    feature.subordinates.push(ms1);
    feature.subordinates.push(ms2);
    let map = FeatureMap::from_features(vec![feature]);

    let writer = OpenSwathOSWWriter::new("", false);
    let rows = writer
        .prepare_rows(&LightCompound::default(), None, &map, "77")
        .unwrap();
    assert_eq!(rows.feature_rows.len(), 1);
    assert_eq!(rows.feature_ms2_rows.len(), 1);
    assert!(rows.feature_ms1_rows.is_empty(), "no var_ms1_ppm_diff");
    // FEATURE: ID with the sign bit cleared, RUN_ID, PRECURSOR_ID parsed from
    // the id text, EXP_RT, EXP_IM null, NORM_RT/DELTA_RT -1, widths
    let f = &rows.feature_rows[0];
    assert_eq!(f[0], OswValue::Int64(0x123));
    assert_eq!(f[1], OswValue::Int64(0));
    assert_eq!(f[2], OswValue::Int64(77));
    assert_eq!(f[3], OswValue::Double(10.0));
    assert_eq!(f[4], OswValue::Null);
    assert_eq!(f[5], OswValue::Double(-1.0));
    assert_eq!(f[6], OswValue::Double(-1.0));
    assert_eq!(f[7], OswValue::Double(9.0));
    assert_eq!(f[8], OswValue::Double(11.0));
    // FEATURE_PRECURSOR: the isotope is the text after Precursor_i
    assert_eq!(rows.feature_precursor_rows.len(), 1);
    let p = &rows.feature_precursor_rows[0];
    assert_eq!(p[1], OswValue::Int64(0));
    assert_eq!(p[2], OswValue::Double(10.0));
    assert_eq!(p[3], OswValue::Double(5.0));
    // FEATURE_TRANSITION: native_id 123 with MASSERROR_PPM 12.5
    assert_eq!(rows.feature_transition_rows.len(), 1);
    let t = &rows.feature_transition_rows[0];
    assert_eq!(t[1], OswValue::Int64(123));
    assert_eq!(t[2], OswValue::Double(50.0));
    assert_eq!(t[3], OswValue::Double(500.0));
    assert_eq!(t[4], OswValue::Double(10.5));
    assert_eq!(t[5], OswValue::Double(25.0));
    assert_eq!(t[6], OswValue::Double(2.0));
    assert_eq!(t[7], OswValue::Double(12.5));
    assert_eq!(t[8], OswValue::Null);
    assert_eq!(t[30], OswValue::Null);

    // the SQL text of prepareLine, in the source table order
    let sql = writer
        .prepare_line(&LightCompound::default(), None, &map, "77")
        .unwrap();
    let order: Vec<usize> = [
        "INSERT INTO FEATURE (",
        "INSERT INTO FEATURE_MS1 (",
        "INSERT INTO FEATURE_PRECURSOR (",
        "INSERT INTO FEATURE_MS2 (",
        "INSERT INTO FEATURE_TRANSITION (",
    ]
    .iter()
    .filter_map(|needle| sql.find(needle))
    .collect();
    assert_eq!(order.len(), 4, "FEATURE_MS1 has no row: {sql}");
    assert!(order.windows(2).all(|w| w[0] < w[1]));
    assert!(sql.contains("MASSERROR_PPM"));
    assert!(sql.contains("VALUES (291, 0, 77, 10.0, NULL, -1.0, -1.0, 9.0, 11.0, NULL, NULL); "));
}

#[test]
fn prepare_rows_uses_the_available_uis_transition_names() {
    // OpenSwathOSWWriter_test: id_target_num_transitions says 2 but only one
    // transition name is available: one row, with TRANSITION_ID 201 and
    // AREA_INTENSITY 10.0
    let mut feature = osw_feature();
    feature
        .metadata
        .insert("id_target_num_transitions".into(), MetaValue::from(2_i64));
    feature.metadata.insert(
        "id_target_transition_names".into(),
        MetaValue::from(vec![201_i64]),
    );
    let list = |name: &str, v: f64| (name.to_owned(), MetaValue::try_from(vec![v]).unwrap());
    for (name, value) in [
        list("id_target_area_intensity", 10.0),
        list("id_target_total_area_intensity", 20.0),
        list("id_target_apex_intensity", 30.0),
        list("id_target_peak_apex_position", 10.5),
        list("id_target_width_at_50", 2.0),
    ] {
        feature.metadata.insert(name, value);
    }
    let map = FeatureMap::from_features(vec![feature]);
    let writer = OpenSwathOSWWriter::new("", true);
    let rows = writer
        .prepare_rows(&LightCompound::default(), None, &map, "77")
        .unwrap();
    assert_eq!(rows.feature_transition_rows.len(), 1);
    let t = &rows.feature_transition_rows[0];
    assert_eq!(t[1], OswValue::Int64(201));
    assert_eq!(t[2], OswValue::Double(10.0));
    assert_eq!(t[3], OswValue::Double(20.0));
    assert_eq!(t[4], OswValue::Double(10.5));
    assert_eq!(t[5], OswValue::Double(30.0));
    assert_eq!(t[6], OswValue::Double(2.0));
    // the other id lists are absent: NULL cells
    assert_eq!(t[9], OswValue::Null);
}

#[test]
fn get_score_and_separate_score() {
    let mut feature = osw_feature();
    feature
        .metadata
        .insert("an_int".into(), MetaValue::from(7_i64));
    feature
        .metadata
        .insert("a_string".into(), MetaValue::from("0.25"));
    feature.metadata.insert(
        "a_list".into(),
        MetaValue::try_from(vec![1.0, 2.5]).unwrap(),
    );
    feature
        .metadata
        .insert("an_int_list".into(), MetaValue::from(vec![3_i64, 4]));
    let writer = OpenSwathOSWWriter::new("", false);
    assert_eq!(writer.get_score(&feature, "leftWidth"), "9.0");
    assert_eq!(writer.get_score(&feature, "an_int"), "7");
    assert_eq!(writer.get_score(&feature, "a_string"), "0.25");
    assert_eq!(writer.get_score(&feature, "missing"), "NULL");
    assert_eq!(
        writer.get_separate_score(&feature, "a_list"),
        ["1.0", "2.5"]
    );
    assert_eq!(
        writer.get_separate_score(&feature, "an_int_list"),
        ["3", "4"]
    );
    assert_eq!(writer.get_separate_score(&feature, "an_int"), ["7"]);
    assert!(writer.get_separate_score(&feature, "missing").is_empty());
}

#[test]
fn osw_data_buffer_operations() {
    let mut a = OswData::default();
    assert!(a.is_empty());
    a.reserve(4);
    a.feature_rows.push(Default::default());
    let mut b = OswData::default();
    b.feature_rows.push(Default::default());
    b.feature_transition_rows
        .push([const { OswValue::Null }; 43]);
    a.append(b);
    assert_eq!(a.feature_rows.len(), 2);
    assert_eq!(a.row_count(), 3);
    assert!(!a.is_empty());
    assert!(a.estimate_memory_usage() > 0);
    assert_eq!(FEATURE_COLUMNS.len(), 11);
    assert_eq!(FEATURE_MS1_COLUMNS.len(), 18);
    assert_eq!(FEATURE_MS2_COLUMNS.len(), 38);
    assert_eq!(FEATURE_PRECURSOR_COLUMNS.len(), 4);
    assert_eq!(FEATURE_TRANSITION_COLUMNS.len(), 43);
}

#[test]
fn osw_writer_state() {
    // OpenSwathOSWWriter_test: an empty file name is inactive
    let mut inactive = OpenSwathOSWWriter::new("", false);
    assert!(!inactive.is_active());
    assert!(OpenSwathOSWWriter::new("/tmp/some.osw", false).is_active());
    // setRunId clears the sign bit
    inactive.set_run_id(u64::MAX);
    assert_eq!(inactive.run_id(), u64::MAX >> 1);
}

#[cfg(not(feature = "sqlite"))]
#[test]
fn osw_io_without_sqlite_is_unsupported() {
    let mut writer = OpenSwathOSWWriter::new("/tmp/never_created.osw", false);
    assert!(matches!(writer.write_header(), Err(Error::Unsupported(_))));
    assert!(matches!(
        writer.add_run(1, "x.mzML"),
        Err(Error::Unsupported(_))
    ));
    assert!(matches!(
        writer.write_lines(&[]),
        Err(Error::Unsupported(_))
    ));
    // nothing to write is not an error, and an inactive writer ignores runs
    assert!(writer.write_rows(&OswData::default()).is_ok());
    assert!(OpenSwathOSWWriter::new("", false).add_run(1, "x").is_ok());
}
