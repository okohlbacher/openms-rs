// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Executed-differential tests of the OpenSWATH export classes (evidence
//! tier 1).
//!
//! The outputs under `tests/data/openswath_export/oracle_*.out` were produced
//! by the C++ Release build (`openms4-release-bc9cc12-c19e494-174b576`, run
//! twice on ibminode06, identical) from the driver `driver.cpp` and the case
//! file `cases.tsv` of the oracle directory `O5a-export`; hashes, compiler,
//! flags and host are in `tests/data/openswath_export_provenance.json`. Each
//! test replays the same cases through this crate and compares the text,
//! byte for byte. The case file is generated (`gen_cases.py`, fixed seed) and
//! committed: 34 matrix scenarios of the matrix exporter over precursor,
//! peptide, protein and gene level with every normalization, consistent and
//! per-run top-N selection, ties, label collisions and sparse runs; one
//! results scenario of 40 rows with every one of the 59 columns exercised; the
//! row preparation of the OSW writer on feature maps built in
//! [`osw_scenarios`]; and the planning functions of the scheduler over a grid.

use openms::analysis::openswath::openswath_export_config::{
    OpenSwathMatrixExportConfig, OpenSwathResultsExportConfig,
};
use openms::analysis::openswath::openswath_export_data::{
    OpenSwathExportRow, OpenSwathMatrixLevel, OpenSwathMatrixNormalization,
};
use openms::analysis::openswath::openswath_matrix_exporter::OpenSwathMatrixExporter;
use openms::analysis::openswath::openswath_osw_writer::{OpenSwathOSWWriter, OswData, OswValue};
use openms::analysis::openswath::openswath_results_exporter::OpenSwathResultsExporter;
use openms::analysis::openswath::openswath_workflow_scheduler::{
    ConcurrencyEstimate, OpenSwathWorkflowScheduler as Scheduler, SchedulerOptions,
};
use openms::analysis::openswath::swath_map::SwathMap;
use openms::analysis::openswath::transition_experiment::LightCompound;
use openms::kernel::{Feature, FeatureMap};
use openms::metadata::MetaValue;
use std::fmt::Write as _;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/openswath_export")
        .join(name)
}

fn read(name: &str) -> String {
    std::fs::read_to_string(data(name)).unwrap()
}

fn level_of(s: &str) -> OpenSwathMatrixLevel {
    match s {
        "precursor" => OpenSwathMatrixLevel::Precursor,
        "peptide" => OpenSwathMatrixLevel::Peptide,
        "protein" => OpenSwathMatrixLevel::Protein,
        _ => OpenSwathMatrixLevel::Gene,
    }
}

fn norm_of(s: &str) -> OpenSwathMatrixNormalization {
    match s {
        "none" => OpenSwathMatrixNormalization::None,
        "median" => OpenSwathMatrixNormalization::Median,
        _ => OpenSwathMatrixNormalization::MedianMedian,
    }
}

fn opt_f(s: &str) -> Option<f64> {
    (s != "-").then(|| s.parse().unwrap())
}

fn opt_i(s: &str) -> Option<i64> {
    (s != "-").then(|| s.parse().unwrap())
}

/// Replay the matrix scenarios; the output has the layout of the driver.
fn matrix_output() -> String {
    let mut out = String::new();
    let mut name = String::new();
    let mut config = OpenSwathMatrixExportConfig::default();
    let mut rows: Vec<OpenSwathExportRow> = Vec::new();
    for line in read("cases.tsv").lines() {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "M" => {
                name = f[1].to_owned();
                config = OpenSwathMatrixExportConfig {
                    level: level_of(f[2]),
                    normalization: norm_of(f[3]),
                    top_n: f[4].parse().unwrap(),
                    consistent_top: f[5] == "1",
                    ..OpenSwathMatrixExportConfig::default()
                };
                rows.clear();
            }
            "R" => rows.push(OpenSwathExportRow {
                run_id: f[1].parse().unwrap(),
                run_name: f[2].into(),
                transition_group_id: f[3].into(),
                sequence: f[4].into(),
                full_peptide_name: f[5].into(),
                protein_name: f[6].into(),
                gene_name: f[7].into(),
                charge: f[8].parse().unwrap(),
                intensity: f[9].parse().unwrap(),
                m_score: f[10].parse().unwrap(),
                feature_id: f[11].parse().unwrap(),
                ..OpenSwathExportRow::default()
            }),
            "E" if !name.is_empty() => {
                let _ = writeln!(out, "### {name}");
                match OpenSwathMatrixExporter::build_matrix(&rows, &config) {
                    Ok(matrix) => {
                        let mut bytes = Vec::new();
                        // writeMatrix rejects an empty matrix; build_matrix never makes one
                        OpenSwathMatrixExporter::write_tsv(&mut bytes, &matrix).unwrap();
                        out.push_str(&String::from_utf8(bytes).unwrap());
                    }
                    Err(error) => {
                        let _ = writeln!(out, "ERROR\t{error}");
                    }
                }
                name.clear();
            }
            _ => {}
        }
    }
    out
}

/// The expected output with the exception texts of C++ reduced to their
/// message: `ERROR<TAB>` lines keep only what follows the last `: `.
fn normalize_errors(text: &str) -> String {
    text.lines()
        .map(|l| {
            if let Some(rest) = l.strip_prefix("ERROR\t") {
                let message = rest.rsplit(": ").next().unwrap_or(rest);
                format!("ERROR\t{message}")
            } else {
                l.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

#[test]
fn matrix_exporter_matches_the_release_build() {
    let expected = normalize_errors(&read("oracle_matrix.out"));
    let actual = normalize_errors(&matrix_output());
    let expected_sections: Vec<&str> = expected.split("### ").collect();
    let actual_sections: Vec<&str> = actual.split("### ").collect();
    assert_eq!(expected_sections.len(), actual_sections.len());
    assert!(
        expected_sections.len() > 30,
        "{} sections",
        expected_sections.len()
    );
    for (e, a) in expected_sections.iter().zip(&actual_sections) {
        assert_eq!(a, e, "scenario {}", e.lines().next().unwrap_or(""));
    }
}

/// The 40 rows of the results scenario of `cases.tsv`.
fn results_rows() -> Vec<OpenSwathExportRow> {
    let mut rows = Vec::new();
    for line in read("cases.tsv").lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f[0] != "F" {
            continue;
        }
        let mut i = 1;
        let mut next = || {
            let value = f[i];
            i += 1;
            value
        };
        let mut r = OpenSwathExportRow {
            run_id: next().parse().unwrap(),
            filename: next().into(),
            run_name: next().into(),
            feature_id: next().parse().unwrap(),
            peptide_id: next().parse().unwrap(),
            transition_group_id: next().into(),
            precursor_id: next().parse().unwrap(),
            decoy: next() == "1",
            sequence: next().into(),
            full_peptide_name: next().into(),
            protein_name: next().into(),
            gene_name: next().into(),
            charge: next().parse().unwrap(),
            mz: next().parse().unwrap(),
            rt: next().parse().unwrap(),
            assay_rt: next().parse().unwrap(),
            delta_rt: next().parse().unwrap(),
            irt: next().parse().unwrap(),
            assay_irt: next().parse().unwrap(),
            delta_irt: next().parse().unwrap(),
            intensity: next().parse().unwrap(),
            aggr_prec_peak_area: opt_f(next()),
            aggr_prec_peak_apex: opt_f(next()),
            left_width: next().parse().unwrap(),
            right_width: next().parse().unwrap(),
            exp_im: opt_f(next()),
            im_left_width: opt_f(next()),
            im_right_width: opt_f(next()),
            ms1_pep: opt_f(next()),
            ms2_pep: opt_f(next()),
            precursor_pep: opt_f(next()),
            ipf_pep: opt_f(next()),
            peak_group_rank: next().parse().unwrap(),
            d_score: next().parse().unwrap(),
            m_score: next().parse().unwrap(),
            pep: opt_f(next()),
            ms2_m_score: opt_f(next()),
            ..OpenSwathExportRow::default()
        };
        r.peptide_run_specific_qvalue = opt_f(next());
        r.peptide_experiment_wide_qvalue = opt_f(next());
        r.peptide_global_qvalue = opt_f(next());
        r.protein_run_specific_qvalue = opt_f(next());
        r.protein_experiment_wide_qvalue = opt_f(next());
        r.protein_global_qvalue = opt_f(next());
        r.gene_run_specific_qvalue = opt_f(next());
        r.gene_experiment_wide_qvalue = opt_f(next());
        r.gene_global_qvalue = opt_f(next());
        r.alignment_group_id = opt_i(next());
        r.alignment_reference_feature_id = opt_i(next());
        r.alignment_reference_rt = opt_f(next());
        r.alignment_pep = opt_f(next());
        r.alignment_qvalue = opt_f(next());
        r.from_alignment = next() == "1";
        r.aggr_peak_area = next().into();
        r.aggr_peak_apex = next().into();
        r.aggr_fragment_annotation = next().into();
        r.ipf_full_peptide_name = next().into();
        r.ipf_precursor_peakgroup_pep = opt_f(next());
        r.ipf_peptidoform_pep = opt_f(next());
        r.ipf_peptidoform_m_score = opt_f(next());
        rows.push(r);
    }
    rows
}

#[test]
fn results_exporter_matches_the_release_build() {
    let rows = results_rows();
    assert_eq!(rows.len(), 40);
    let mut actual = String::from("### results_all_fields\n");
    actual.push_str(&OpenSwathResultsExporter::tsv_string(&rows).unwrap());
    assert_eq!(actual, read("oracle_results.out"));

    // and through `write`, which is what the C++ driver calls
    let path = std::env::temp_dir().join(format!("openms_o5a_results_{}.tsv", std::process::id()));
    OpenSwathResultsExporter::write(&path, &rows, &OpenSwathResultsExportConfig::default())
        .unwrap();
    let written = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    assert_eq!(
        format!("### results_all_fields\n{written}"),
        read("oracle_results.out")
    );
}

// --- The Parquet arm ----------------------------------------------------------
//
// `../oracle/integ-O-parquet` (ibminode06, Release install, run twice,
// identical) wrote `OpenSwathMatrixExporter::writeMatrix` and
// `OpenSwathResultsExporter::write` with format Parquet for every scenario of
// the case file above. The files are retained under
// `tests/data/openswath_export/parquet/`; each test writes the same scenario
// through this crate, reads both files back and compares the tables: schema
// (names, Arrow types, nullability) and every value.

#[cfg(feature = "arrow")]
fn parquet_fixture(name: &str) -> PathBuf {
    data(&format!("parquet/pq_{name}.parquet"))
}

#[cfg(feature = "arrow")]
fn assert_same_table(
    actual: &arrow_array::RecordBatch,
    expected: &arrow_array::RecordBatch,
    what: &str,
) {
    assert_eq!(actual.schema(), expected.schema(), "{what}: schema");
    assert_eq!(actual.num_rows(), expected.num_rows(), "{what}: rows");
    assert_eq!(
        actual.num_columns(),
        expected.num_columns(),
        "{what}: columns"
    );
    for (i, (a, e)) in actual.columns().iter().zip(expected.columns()).enumerate() {
        assert_eq!(
            a.as_ref(),
            e.as_ref(),
            "{what}: column {i} ({})",
            actual.schema().field(i).name()
        );
    }
}

#[cfg(feature = "arrow")]
#[test]
fn matrix_exporter_parquet_arm_matches_the_release_build() {
    use openms::analysis::openswath::openswath_export_data::OpenSwathExportFileFormat;
    use openms::format::parquet_file::ParquetFile;
    let mut name = String::new();
    let mut config = OpenSwathMatrixExportConfig::default();
    let mut rows: Vec<OpenSwathExportRow> = Vec::new();
    let mut scenarios = 0;
    let mut errors = 0;
    let dir = std::env::temp_dir().join(format!("openms_integ_o_pq_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let expected_status = read("oracle_matrix_parquet.out");
    for line in read("cases.tsv").lines() {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "M" => {
                name = f[1].to_owned();
                config = OpenSwathMatrixExportConfig {
                    level: level_of(f[2]),
                    normalization: norm_of(f[3]),
                    top_n: f[4].parse().unwrap(),
                    consistent_top: f[5] == "1",
                    format: OpenSwathExportFileFormat::Parquet,
                    ..OpenSwathMatrixExportConfig::default()
                };
                rows.clear();
            }
            "R" if !name.is_empty() => rows.push(OpenSwathExportRow {
                run_id: f[1].parse().unwrap(),
                run_name: f[2].into(),
                transition_group_id: f[3].into(),
                sequence: f[4].into(),
                full_peptide_name: f[5].into(),
                protein_name: f[6].into(),
                gene_name: f[7].into(),
                charge: f[8].parse().unwrap(),
                intensity: f[9].parse().unwrap(),
                m_score: f[10].parse().unwrap(),
                feature_id: f[11].parse().unwrap(),
                ..OpenSwathExportRow::default()
            }),
            "E" if !name.is_empty() => {
                scenarios += 1;
                let section = expected_status
                    .split("### ")
                    .find(|s| s.lines().next() == Some(name.as_str()))
                    .unwrap();
                let cpp_ok = section.lines().nth(1) == Some("OK");
                match OpenSwathMatrixExporter::build_matrix(&rows, &config) {
                    Ok(matrix) => {
                        assert!(
                            cpp_ok,
                            "{name}: the Release build failed, the port built a matrix"
                        );
                        let path = dir.join(format!("{name}.parquet"));
                        OpenSwathMatrixExporter::write_matrix(&path, &matrix, &config).unwrap();
                        let ours = ParquetFile::read_table(&path).unwrap();
                        let theirs = ParquetFile::read_table(parquet_fixture(&name)).unwrap();
                        assert_same_table(&ours, &theirs, &name);
                        // the typed columns the port builds are the same table
                        let built = OpenSwathMatrixExporter::table(&matrix)
                            .unwrap()
                            .to_record_batch()
                            .unwrap();
                        assert_same_table(&built, &theirs, &format!("{name} (in memory)"));
                    }
                    Err(_) => {
                        errors += 1;
                        assert!(
                            !cpp_ok,
                            "{name}: the port refused a matrix the Release build wrote"
                        );
                    }
                }
                name.clear();
            }
            _ => {}
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
    assert!(scenarios >= 34, "{scenarios} scenarios");
    assert!(scenarios > errors);
}

#[cfg(feature = "arrow")]
#[test]
fn results_exporter_parquet_arm_matches_the_release_build() {
    use openms::analysis::openswath::openswath_export_data::OpenSwathExportFileFormat;
    use openms::format::parquet_file::ParquetFile;
    let rows = results_rows();
    let config = OpenSwathResultsExportConfig {
        format: OpenSwathExportFileFormat::Parquet,
        ..OpenSwathResultsExportConfig::default()
    };
    let dir = std::env::temp_dir().join(format!("openms_integ_o_pqr_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("results.parquet");
    OpenSwathResultsExporter::write(&path, &rows, &config).unwrap();
    let ours = ParquetFile::read_table(&path).unwrap();
    let theirs = ParquetFile::read_table(parquet_fixture("results_all_fields")).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
    assert_eq!(theirs.num_rows(), 40);
    assert_eq!(theirs.num_columns(), 59);
    assert_same_table(&ours, &theirs, "results_all_fields");
}

// --- OSW row preparation --------------------------------------------------

fn dump_value(v: &OswValue) -> Option<String> {
    match v {
        OswValue::Null => None,
        OswValue::Int64(i) => Some(format!("I:{i}")),
        OswValue::Double(d) if d.is_nan() => None,
        OswValue::Double(d) => Some(format!("D:{:016x}", d.to_bits())),
        OswValue::Text(t) => Some(format!("T:{t}")),
    }
}

fn dump_value_with_null(v: &OswValue) -> String {
    dump_value(v).unwrap_or_else(|| "N".into())
}

fn dump_rows<const N: usize>(out: &mut String, table: &str, rows: &[[OswValue; N]]) {
    for (r, row) in rows.iter().enumerate() {
        let _ = write!(out, "{table}[{r}]");
        for (i, value) in row.iter().enumerate() {
            if let Some(text) = dump_value(value) {
                let _ = write!(out, " {i}={text}");
            }
        }
        out.push('\n');
    }
}

fn dump_data(out: &mut String, data: &OswData) {
    dump_rows(out, "FEATURE", &data.feature_rows);
    dump_rows(out, "FEATURE_MS1", &data.feature_ms1_rows);
    dump_rows(out, "FEATURE_MS2", &data.feature_ms2_rows);
    dump_rows(out, "FEATURE_PRECURSOR", &data.feature_precursor_rows);
    dump_rows(out, "FEATURE_TRANSITION", &data.feature_transition_rows);
}

fn real(v: f64) -> MetaValue {
    MetaValue::try_from(v).unwrap()
}

fn set_d(f: &mut Feature, name: &str, v: f64) {
    f.metadata.insert(name.into(), real(v));
}

fn base_feature(uid: u64, rt: f64, intensity: f64) -> Feature {
    let mut f = Feature::new(rt, 0.0, intensity as f32);
    f.unique_id = uid;
    f
}

fn subordinate(level: &str, intensity: f64, native_id: MetaValue) -> Feature {
    let mut s = Feature::new(0.0, 0.0, intensity as f32);
    s.metadata
        .insert("FeatureLevel".into(), MetaValue::from(level));
    s.metadata.insert("native_id".into(), native_id);
    s
}

fn dlist(v: &[f64]) -> MetaValue {
    MetaValue::try_from(v.to_vec()).unwrap()
}

fn ilist(v: &[i64]) -> MetaValue {
    MetaValue::from(v.to_vec())
}

fn slist(v: &[&str]) -> MetaValue {
    MetaValue::from(v.iter().map(|s| s.to_string()).collect::<Vec<_>>())
}

/// The scenarios of the driver's `runOsw`, built the same way.
fn osw_output() -> String {
    let pep = LightCompound::default();
    let mut out = String::new();

    // A
    out.push_str("### masserror_by_ms2_position\n");
    {
        let writer = OpenSwathOSWWriter::new("", false);
        let mut f = base_feature(0x8000_0000_0000_0123, 10.0, 100.0);
        set_d(&mut f, "leftWidth", 9.0);
        set_d(&mut f, "rightWidth", 11.0);
        set_d(&mut f, "total_xic", 1000.0);
        set_d(&mut f, "peak_apices_sum", 100.0);
        f.metadata.insert("masserror_ppm".into(), dlist(&[12.5]));
        let mut ms1 = subordinate("MS1", 10.0, MetaValue::from("Precursor_i0"));
        set_d(&mut ms1, "peak_apex_int", 5.0);
        let mut ms2 = subordinate("MS2", 50.0, MetaValue::from(123_i64));
        set_d(&mut ms2, "total_xic", 500.0);
        set_d(&mut ms2, "peak_apex_int", 25.0);
        set_d(&mut ms2, "peak_apex_position", 10.5);
        set_d(&mut ms2, "width_at_50", 2.0);
        f.subordinates.push(ms1);
        f.subordinates.push(ms2);
        let map = FeatureMap::from_features(vec![f]);
        dump_data(
            &mut out,
            &writer.prepare_rows(&pep, None, &map, "77").unwrap(),
        );
        out.push_str(&writer.prepare_line(&pep, None, &map, "77").unwrap());
        out.push('\n');
    }
    // B
    out.push_str("### uis_with_peak_shape\n");
    {
        let writer = OpenSwathOSWWriter::new("", true);
        let mut f = base_feature(7, 12.5, 250.0);
        set_d(&mut f, "leftWidth", 9.0);
        set_d(&mut f, "rightWidth", 15.0);
        let m = &mut f.metadata;
        m.insert("id_target_num_transitions".into(), MetaValue::from(3_i64));
        m.insert("id_target_transition_names".into(), ilist(&[11, 12, 13]));
        m.insert("id_target_area_intensity".into(), dlist(&[1.5, 2.5, 3.5]));
        m.insert(
            "id_target_total_area_intensity".into(),
            dlist(&[10.0, 20.0, 30.0]),
        );
        m.insert(
            "id_target_apex_intensity".into(),
            dlist(&[100.0, 200.0, 300.0]),
        );
        m.insert(
            "id_target_peak_apex_position".into(),
            dlist(&[12.4, 12.5, 12.6]),
        );
        m.insert("id_target_width_at_50".into(), dlist(&[1.0, 2.0, 3.0]));
        m.insert("id_target_total_mi".into(), dlist(&[0.1, 0.2, 0.3]));
        m.insert(
            "id_target_ind_massdev_score".into(),
            dlist(&[5.0, 6.0, 7.0]),
        );
        m.insert(
            "id_target_ind_start_position_at_5".into(),
            dlist(&[1.25, 1.5, 1.75]),
        );
        m.insert(
            "id_target_ind_end_position_at_5".into(),
            dlist(&[2.25, 2.5, 2.75]),
        );
        m.insert("id_target_ind_total_width".into(), dlist(&[3.0, 3.5, 4.0]));
        m.insert("id_decoy_num_transitions".into(), MetaValue::from(2_i64));
        m.insert("id_decoy_transition_names".into(), slist(&["401", "d_two"]));
        m.insert("id_decoy_area_intensity".into(), dlist(&[0.5, 0.25]));
        m.insert(
            "id_decoy_ind_start_position_at_5".into(),
            dlist(&[9.0, 8.0]),
        );
        m.insert("id_decoy_ind_im_drift".into(), dlist(&[0.9, 1.1]));
        let map = FeatureMap::from_features(vec![f]);
        dump_data(
            &mut out,
            &writer.prepare_rows(&pep, None, &map, "PEPTIDE_2").unwrap(),
        );
    }
    // C
    out.push_str("### uis_zero_start_and_fallback\n");
    {
        let writer = OpenSwathOSWWriter::new("", true);
        let mut f = base_feature(8, 20.0, 50.0);
        f.metadata
            .insert("id_target_num_transitions".into(), MetaValue::from(2_i64));
        f.metadata
            .insert("id_target_transition_names".into(), ilist(&[21, 22]));
        f.metadata
            .insert("id_target_ind_start_position_at_5".into(), ilist(&[0, 0]));
        f.metadata
            .insert("id_target_ind_end_position_at_5".into(), ilist(&[7, 7]));
        let mut g = base_feature(9, 21.0, 40.0);
        g.metadata
            .insert("id_target_num_transitions".into(), MetaValue::from(0_i64));
        g.metadata
            .insert("id_target_transition_names".into(), ilist(&[31]));
        let mut h = base_feature(10, 22.0, 30.0);
        let mut s = subordinate("MS2", 5.0, MetaValue::from("99"));
        set_d(&mut s, "total_xic", 42.0);
        s.metadata.insert("start_position_at_5".into(), real(3.0));
        s.metadata.insert("end_position_at_5".into(), real(4.0));
        h.subordinates.push(s);
        let map = FeatureMap::from_features(vec![f, g, h]);
        dump_data(
            &mut out,
            &writer.prepare_rows(&pep, None, &map, "5").unwrap(),
        );
    }
    // D
    out.push_str("### ms1_scores_and_precursor_rows\n");
    {
        let mut writer = OpenSwathOSWWriter::new("", false);
        writer.set_run_id(u64::MAX);
        let mut f = base_feature(1, 30.0, 1.0e6);
        set_d(&mut f, "var_ms1_ppm_diff", 1.5);
        set_d(&mut f, "ms1_area_intensity", 11.0);
        set_d(&mut f, "ms1_apex_intensity", 12.0);
        f.metadata.insert("norm_RT".into(), real(0.5));
        f.metadata.insert("delta_rt".into(), MetaValue::from(2_i64));
        set_d(&mut f, "im_drift", 1.01);
        set_d(&mut f, "im_drift_left", 1.0);
        set_d(&mut f, "im_drift_right", 1.02);
        set_d(&mut f, "var_xcorr_shape", 0.75);
        f.metadata
            .insert("var_library_corr".into(), MetaValue::from("0.5"));
        f.metadata
            .insert("var_log_sn_score".into(), MetaValue::from("not a number"));
        f.subordinates
            .push(subordinate("MS1", 0.0, MetaValue::from("7_Precursor_i0")));
        let mut p1 = subordinate("MS1", 3.0, MetaValue::from("7_Precursor_i2"));
        set_d(&mut p1, "peak_apex_int", 2.0);
        f.subordinates.push(p1);
        f.subordinates
            .push(subordinate("MS1", 4.0, MetaValue::from(55_i64)));
        f.subordinates
            .push(subordinate("MS2", 8.0, MetaValue::from("7_y3")));
        f.subordinates
            .push(subordinate("OTHER", 8.0, MetaValue::from("ignored")));
        let map = FeatureMap::from_features(vec![f]);
        dump_data(
            &mut out,
            &writer.prepare_rows(&pep, None, &map, "NIAGGLEK_2").unwrap(),
        );
    }
    // E
    out.push_str("### scores\n");
    {
        let writer = OpenSwathOSWWriter::new("", false);
        let mut f = base_feature(1, 1.0, 1.0);
        let m = &mut f.metadata;
        m.insert("an_int".into(), MetaValue::from(7_i64));
        m.insert("a_double".into(), real(0.1));
        m.insert("a_big".into(), real(123_456_789.125));
        m.insert("a_small".into(), real(1.0e-7));
        m.insert("a_string_num".into(), MetaValue::from("0.25"));
        m.insert("a_string_int".into(), MetaValue::from("12"));
        m.insert("a_string_text".into(), MetaValue::from("hello"));
        m.insert("a_double_list".into(), dlist(&[1.0, 2.5, 1e-7]));
        m.insert("an_int_list".into(), ilist(&[3, -4]));
        m.insert("a_string_list".into(), slist(&["5", "x", "NaN", "1.5"]));
        for n in [
            "an_int",
            "a_double",
            "a_big",
            "a_small",
            "a_string_num",
            "a_string_int",
            "a_string_text",
            "a_double_list",
            "an_int_list",
            "a_string_list",
            "missing",
        ] {
            let _ = writeln!(out, "getScore {n} = {}", writer.get_score(&f, n));
            let _ = write!(out, "getSeparateScore {n} =");
            for x in writer.get_separate_score(&f, n) {
                let _ = write!(out, " [{x}]");
            }
            out.push('\n');
        }
        for t in [
            "",
            "NULL",
            "nan",
            "-nan",
            "NaN",
            "-NaN",
            "NAN",
            "0",
            "-17",
            "+5",
            " 12 ",
            "0012",
            "100.0",
            "1e3",
            ".5",
            "5.",
            "99999999999999999999",
            "9223372036854775807",
            "9223372036854775808",
            "-9223372036854775808",
            "inf",
            "-inf",
            "Infinity",
            "1e999",
            "1e-999",
            "12abc",
            "NIAGGLEK_2",
            "1 2",
            "0x10",
            "0x1p3",
            "1,5",
            "1e5",
            "1E5",
            "+.5e-2",
            "  ",
            "\t7",
            "7\n",
            "--1",
            "nan(1)",
        ] {
            let _ = writeln!(
                out,
                "OSWValue [{t}] = {}",
                dump_value_with_null(&OswValue::parse(t))
            );
        }
    }
    out
}

#[test]
fn osw_row_preparation_matches_the_release_build() {
    let expected = read("oracle_osw.out");
    let actual = osw_output();
    let e: Vec<&str> = expected.lines().collect();
    let a: Vec<&str> = actual.lines().collect();
    for (i, (x, y)) in e.iter().zip(&a).enumerate() {
        assert_eq!(y, x, "line {} of the oracle output", i + 1);
    }
    assert_eq!(a.len(), e.len());
}

// --- scheduler ---------------------------------------------------------------

fn maps(n: usize) -> Vec<SwathMap> {
    let mut m = vec![SwathMap::new(0.0, 0.0, 0.0, true)];
    let mut lo = 400.0;
    for i in 0..n {
        if i == n / 2 {
            m.push(SwathMap::new(0.0, 0.0, 0.0, true));
        }
        m.push(SwathMap::new(lo, lo + 25.0, lo + 12.5, false));
        lo += 25.0;
    }
    m
}

fn scheduler_output() -> String {
    let mut out = String::new();
    for n in [0usize, 1, 7, 32] {
        for cap in [1usize, 3, 100] {
            for budget_mib in [0u64, 150, 250, 1000] {
                for avg_spectra in [0usize, 50, 2000] {
                    for avg_tr in [0usize, 10, 1000] {
                        let options = SchedulerOptions {
                            avg_transitions_per_swath: avg_tr,
                            ..SchedulerOptions::default()
                        };
                        let estimate = ConcurrencyEstimate {
                            max_concurrent_swaths: cap,
                            memory_budget_bytes: budget_mib * 1024 * 1024,
                            avg_spectra_per_swath: avg_spectra,
                            ..ConcurrencyEstimate::default()
                        };
                        let waves = Scheduler::plan_waves(&maps(n), &options, &estimate);
                        let _ = write!(
                            out,
                            "plan n={n} cap={cap} budget={budget_mib} spectra={avg_spectra} tr={avg_tr}:"
                        );
                        for w in &waves {
                            out.push_str(" [");
                            let idx: Vec<String> =
                                w.swath_indices.iter().map(|i| i.to_string()).collect();
                            out.push_str(&idx.join(","));
                            let _ = write!(out, "|{}]", w.estimated_bytes);
                        }
                        out.push('\n');
                    }
                }
            }
        }
    }
    let options = SchedulerOptions::default();
    for total in [0usize, 1, 999, 1000, 2000, 5000, 60000, 200_000, 1_000_000] {
        for active in [0usize, 1, 4, 50] {
            for threads in [0usize, 1, 4, 64] {
                for user in [-1i32, 0, 500, 50000] {
                    let _ = writeln!(
                        out,
                        "batch total={total} active={active} threads={threads} user={user} -> {}",
                        Scheduler::choose_inner_batch_size(total, active, threads, user, &options)
                    );
                }
            }
        }
    }
    let mut odd = SchedulerOptions {
        target_jobs_per_thread: 0,
        min_inner_batch_size: 0,
        max_inner_batch_size: 0,
        ..SchedulerOptions::default()
    };
    let _ = writeln!(
        out,
        "batch odd -> {}",
        Scheduler::choose_inner_batch_size(5000, 1, 4, -1, &odd)
    );
    odd.min_inner_batch_size = 300;
    odd.max_inner_batch_size = 100;
    let _ = writeln!(
        out,
        "batch inverted -> {}",
        Scheduler::choose_inner_batch_size(5000, 1, 4, -1, &odd)
    );
    out
}

#[test]
fn scheduler_planning_matches_the_release_build() {
    let expected = read("oracle_scheduler.out");
    let actual = scheduler_output();
    let e: Vec<&str> = expected.lines().collect();
    let a: Vec<&str> = actual.lines().collect();
    assert!(e.len() > 1000);
    for (i, (x, y)) in e.iter().zip(&a).enumerate() {
        assert_eq!(y, x, "line {} of the oracle output", i + 1);
    }
    assert_eq!(a.len(), e.len());
}
