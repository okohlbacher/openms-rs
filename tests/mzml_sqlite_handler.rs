// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source class-test mappings and native corruption/transaction regressions.
#![cfg(feature = "sqmass")]

use openms::format::mzml_sqlite_handler::{LossPolicy, MzMLSqliteHandler};
use openms::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D,
};
use openms::metadata::Polarity;
use openms::system::file::TempDir;
use std::path::{Path, PathBuf};

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/sqlite_s1_source_review/SqliteMassFile_1.sqMass")
}
fn original() -> MSExperiment {
    openms::format::mzml::load(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/mzml_sqlite_handler/MzMLSqliteHandler_1.mzML"),
    )
    .unwrap()
}
fn temporary() -> (TempDir, PathBuf) {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("test.sqMass");
    (dir, path)
}
fn sample() -> MSExperiment {
    let mut s = MSSpectrum {
        native_id: "scan=1".into(),
        rt: 1.5,
        peaks: vec![
            Peak1D {
                mz: 100.0,
                intensity: 50.0,
            },
            Peak1D {
                mz: 101.0,
                intensity: 60.0,
            },
        ],
        ..Default::default()
    };
    s.instrument_settings.polarity = Polarity::Positive;
    let c = MSChromatogram {
        native_id: "TIC".into(),
        peaks: vec![
            ChromatogramPeak {
                rt: 1.0,
                intensity: 50.0,
            },
            ChromatogramPeak {
                rt: 2.0,
                intensity: 60.0,
            },
        ],
        ..Default::default()
    };
    MSExperiment {
        spectra: vec![s],
        chromatograms: vec![c],
        ..Default::default()
    }
}
fn writer(path: &Path) -> MzMLSqliteHandler {
    let mut h = MzMLSqliteHandler::new(path, 12345);
    h.set_config(true, false, 0.0001, 500).unwrap();
    h.create_tables().unwrap();
    h
}
fn sql(path: &Path, statement: &str) {
    rusqlite::Connection::open(path)
        .unwrap()
        .execute_batch(statement)
        .unwrap();
}
fn copied_fixture() -> (TempDir, PathBuf) {
    let (d, p) = temporary();
    std::fs::copy(fixture(), &p).unwrap();
    (d, p)
}
fn similar(a: f64, b: f64, absolute: f64, relative: f64) -> bool {
    a.is_finite()
        && b.is_finite()
        && ((a - b).abs() <= absolute
            || (a != 0.0
                && b != 0.0
                && a.signum() == b.signum()
                && (a / b).abs().max((b / a).abs()) <= relative))
}
#[test]
fn comparison_policy_rejects_opposite_signs_except_within_absolute_tolerance() {
    assert!(!similar(1.0, -1.0, 1e-5, 1.001));
    assert!(similar(1e-7, -1e-7, 1e-5, 1.001));
    assert!(!similar(f64::INFINITY, f64::INFINITY, 1e-5, 1.001));
}

fn compare_data(actual: &MSExperiment, expected: &MSExperiment) {
    assert_eq!(actual.spectra.len(), expected.spectra.len());
    assert_eq!(actual.chromatograms.len(), expected.chromatograms.len());
    for (a, e) in actual.spectra.iter().zip(&expected.spectra) {
        assert_eq!(a.len(), e.len());
        for (a, e) in a.peaks.iter().zip(&e.peaks) {
            assert!(
                similar(a.mz, e.mz, 1e-5, 1.000001),
                "mz {} != {}",
                a.mz,
                e.mz
            );
            assert!(similar(a.intensity as f64, e.intensity as f64, 1e-4, 1.001));
        }
    }
    for (a, e) in actual.chromatograms.iter().zip(&expected.chromatograms) {
        assert_eq!(a.len(), e.len());
        for (a, e) in a.peaks.iter().zip(&e.peaks) {
            assert!(similar(a.rt, e.rt, 0.05, 1.000001));
            assert!(similar(a.intensity as f64, e.intensity as f64, 1e-4, 1.001));
        }
    }
}

#[test]
fn constructor_defaults_and_drop_do_not_create_a_file() {
    let (_d, p) = temporary();
    {
        let h = MzMLSqliteHandler::new(&p, u64::MAX);
        assert!(h.config().write_full_meta);
        assert!(h.config().use_lossy_compression);
        assert_eq!(h.config().sql_batch_size, 500);
    }
    assert!(!p.exists());
}
#[test]
fn retained_run_id_and_counts() {
    let h = MzMLSqliteHandler::new(fixture(), 0);
    assert_eq!(h.run_id().unwrap(), 12345);
    assert_eq!(h.nr_spectra().unwrap(), 2);
    assert_eq!(h.nr_chromatograms().unwrap(), 1);
}
#[test]
fn retained_experiment_metadata_and_all_numeric_values() {
    let h = MzMLSqliteHandler::new(fixture(), 0);
    let original = original();
    let meta = h.read_experiment(true).unwrap();
    assert_eq!(meta.spectra.len(), 2);
    assert_eq!(meta.chromatograms.len(), 1);
    assert_eq!(meta.sql_run_id, 12345);
    assert!(meta.spectra.iter().all(|s| s.is_empty()));
    assert!(meta.chromatograms.iter().all(|c| c.is_empty()));
    let mut source_settings = original.settings.clone();
    // C++ setSqlRunID adds a metadata transport key; native SQL ownership is dedicated.
    assert_eq!(
        source_settings
            .metadata
            .remove("sqMassRunID")
            .unwrap()
            .as_i64()
            .unwrap(),
        12345
    );
    assert_eq!(meta.settings, source_settings);
    let all = h.read_experiment(false).unwrap();
    compare_data(&all, &original);
    assert_eq!(all.sql_run_id, 12345);
    assert_eq!(all.settings, source_settings);
}
#[test]
fn selected_spectrum_class_test_literals() {
    let h = MzMLSqliteHandler::new(fixture(), 0);
    for (meta, n) in [(true, 0), (false, 19800)] {
        let s = h.read_spectra(&[1], meta).unwrap();
        assert_eq!(s.len(), 1);
        assert_eq!(s[0].len(), n);
        assert_eq!(s[0].rt, 0.4738);
    }
    let s = h.read_spectra(&[1, 0], false).unwrap();
    assert_eq!(s[0].len(), 19914);
    assert_eq!(s[1].len(), 19800);
    assert_eq!(s[0].rt, 0.2961);
    let meta = h.read_spectra(&[0, 1], true).unwrap();
    assert_eq!(meta.len(), 2);
    assert!(meta.iter().all(MSSpectrum::is_empty));
    assert_eq!((meta[0].rt, meta[1].rt), (0.2961, 0.4738));
    for ids in [&[0, 1, 2][..], &[5][..], &[-1][..], &[0, 0][..], &[][..]] {
        assert!(h.read_spectra(ids, false).is_err());
    }
}
#[test]
fn selected_chromatogram_class_test_literals() {
    let h = MzMLSqliteHandler::new(fixture(), 0);
    let c = h.read_chromatograms(&[0], true).unwrap();
    assert_eq!(c[0].native_id, "TIC");
    assert!(c[0].is_empty());
    for ids in [&[0, 1][..], &[5][..], &[-1][..], &[0, 0][..], &[][..]] {
        assert!(h.read_chromatograms(ids, false).is_err());
    }
    let (_d, p) = temporary();
    let mut h = writer(&p);
    h.set_loss_policy(LossPolicy::Source);
    let mut c = original().chromatograms;
    c.push(c[0].clone());
    c[1].native_id = "second".into();
    h.write_chromatograms(&c).unwrap();
    let read = h.read_chromatograms(&[1, 0], true).unwrap();
    assert_eq!(read[0].native_id, "TIC");
    assert_eq!(read[1].native_id, "second");
    assert_eq!(read.len(), 2);
    assert!(read.iter().all(MSChromatogram::is_empty));
    for (id, name) in [(0, "TIC"), (1, "second")] {
        let selected = h.read_chromatograms(&[id], true).unwrap();
        assert_eq!(selected.len(), 1);
        assert!(selected[0].is_empty());
        assert_eq!(selected[0].native_id, name);
    }
}
#[test]
fn retention_time_class_test_literals_and_stable_nearest() {
    let h = MzMLSqliteHandler::new(fixture(), 0);
    for (rt, delta, ids, expected) in [
        (0.4738, 0.1, vec![], vec![1]),
        (0.296, 0.1, vec![], vec![0]),
        (0.296, 1.1, vec![], vec![0, 1]),
        (0.296, 1.1, vec![1], vec![1]),
        (0.296, 1.1, vec![0], vec![0]),
        (0.0, 0.1, vec![], vec![]),
        (0.3, -0.1, vec![], vec![1]),
        (0.0, -0.1, vec![], vec![0]),
    ] {
        assert_eq!(h.spectra_indices_by_rt(rt, delta, &ids).unwrap(), expected);
    }
    let (_d, p) = copied_fixture();
    sql(
        &p,
        "UPDATE SPECTRUM SET RETENTION_TIME=CASE ID WHEN 0 THEN 9 ELSE 3 END",
    );
    let h = MzMLSqliteHandler::new(&p, 0);
    assert_eq!(h.spectra_indices_by_rt(2.0, 0.0, &[]).unwrap(), [1]);
    assert!(h.spectra_indices_by_rt(f64::NAN, 1.0, &[]).is_err());
}
#[test]
fn full_experiment_write_recreate_and_lossy_source_tolerances() {
    let (_d, p) = temporary();
    let input = original();
    let mut h = MzMLSqliteHandler::new(&p, 12345);
    assert!(h.write_experiment(&input).is_err());
    for _ in 0..2 {
        h.create_tables().unwrap();
        h.write_experiment(&input).unwrap();
        assert_eq!(h.nr_spectra().unwrap(), 2);
        compare_data(&h.read_experiment(false).unwrap(), &input);
        let mut settings = input.settings.clone();
        settings.metadata.remove("sqMassRunID");
        assert_eq!(h.read_experiment(true).unwrap().settings, settings);
    }
}
#[test]
fn run_level_loaded_file_path_is_bound_and_masked() {
    let (_d, p) = temporary();
    let mut h = writer(&p);
    h.set_run_id(u64::MAX);
    let mut e = sample();
    e.settings.document.loaded_file_path = "/tmp/o'brien'; DROP TABLE RUN;--/run.mzML".into();
    h.write_experiment(&e).unwrap();
    assert_eq!(h.run_id().unwrap(), i64::MAX as u64);
    assert_eq!(
        h.read_experiment(false)
            .unwrap()
            .settings
            .document
            .loaded_file_path,
        e.settings.document.loaded_file_path
    );
}
#[test]
fn repeated_spectrum_appends_preserve_source_literals() {
    let (_d, p) = temporary();
    let e = original();
    let mut h = writer(&p);
    h.set_loss_policy(LossPolicy::Source);
    h.set_config(false, true, 0.0001, 1).unwrap();
    for n in 1..=3 {
        h.write_spectra(&e.spectra).unwrap();
        assert_eq!(h.nr_spectra().unwrap(), n * 2);
    }
    h.write_run_level_information(&e, false).unwrap();
    let out = h.read_experiment(false).unwrap();
    assert_eq!(out.spectra.len(), 6);
    assert_eq!(
        out.spectra.iter().map(MSSpectrum::len).collect::<Vec<_>>(),
        [19914, 19800, 19914, 19800, 19914, 19800]
    );
    // Source helpers mutate global comparison state: these literals inherit
    // abs=.05 and relative=1.000001 from the earlier cmpDataRT call.
    assert!(similar(
        out.spectra[0].peaks[100].mz,
        204.817,
        0.05,
        1.000001
    ));
    assert!(similar(
        out.spectra[0].peaks[100].intensity as f64,
        3857.86,
        0.05,
        1.000001
    ));
    h.create_tables().unwrap();
    h.write_spectra(&e.spectra).unwrap();
    assert_eq!(h.read_spectra(&[0], false).unwrap()[0].len(), 19914);
}
#[test]
fn repeated_chromatogram_appends_lossless_and_lossy() {
    for lossy in [false, true] {
        let (_d, p) = temporary();
        let e = original();
        let mut h = writer(&p);
        h.set_loss_policy(LossPolicy::Source);
        h.set_config(false, lossy, 0.0001, 1).unwrap();
        for n in 1..=3 {
            h.write_chromatograms(&e.chromatograms).unwrap();
            assert_eq!(h.nr_chromatograms().unwrap(), n);
        }
        h.write_run_level_information(&e, false).unwrap();
        let out = h.read_experiment(false).unwrap();
        assert_eq!(out.chromatograms.len(), 3);
        for c in &out.chromatograms {
            assert_eq!(c.len(), 48);
            assert!(similar(c.peaks[20].rt, 0.200695, 0.05, 1.0002));
            assert!(similar(
                c.peaks[20].intensity as f64,
                147414.578125,
                0.05,
                1.0002
            ));
        }
        h.create_tables().unwrap();
        h.write_chromatograms(&e.chromatograms).unwrap();
        assert_eq!(h.nr_chromatograms().unwrap(), 1);
    }
}
#[test]
fn lossless_empty_arrays_and_empty_experiment_roundtrip() {
    let (_d, p) = temporary();
    let mut h = writer(&p);
    h.write_experiment(&MSExperiment::default()).unwrap();
    assert!(h.read_experiment(false).unwrap().spectra.is_empty());
    h.create_tables().unwrap();
    let mut e = sample();
    e.spectra[0].peaks.clear();
    e.chromatograms[0].peaks.clear();
    h.write_experiment(&e).unwrap();
    let out = h.read_experiment(false).unwrap();
    assert!(out.spectra[0].is_empty());
    assert!(out.chromatograms[0].is_empty());
}
#[test]
fn builder_can_write_snapshot_before_records() {
    let (_d, p) = temporary();
    let e = sample();
    let mut h = writer(&p);
    h.write_run_level_information(&e, true).unwrap();
    assert!(h.read_experiment(false).is_err());
    h.write_spectra(&e.spectra).unwrap();
    h.write_chromatograms(&e.chromatograms).unwrap();
    compare_data(&h.read_experiment(false).unwrap(), &e);
}
#[test]
fn gapped_ids_and_shuffled_data_use_record_identity() {
    let (_d, p) = copied_fixture();
    sql(
        &p,
        "UPDATE SPECTRUM SET ID=100 WHERE ID=0; UPDATE DATA SET SPECTRUM_ID=100 WHERE SPECTRUM_ID=0; UPDATE PRECURSOR SET SPECTRUM_ID=100 WHERE SPECTRUM_ID=0; UPDATE PRODUCT SET SPECTRUM_ID=100 WHERE SPECTRUM_ID=0; CREATE TEMP TABLE COPY AS SELECT * FROM DATA; DELETE FROM DATA; INSERT INTO DATA SELECT * FROM COPY ORDER BY SPECTRUM_ID DESC,DATA_TYPE DESC;",
    );
    let h = MzMLSqliteHandler::new(&p, 0);
    let s = h.read_spectra(&[100, 1], false).unwrap();
    assert_eq!(s[0].len(), 19800);
    assert_eq!(s[1].len(), 19914);
    assert!(h.read_experiment(false).is_err()); // snapshot order/identity no longer agrees
}
#[test]
fn every_public_write_rolls_back_and_reuses_counters() {
    let (_d, p) = temporary();
    let mut h = writer(&p);
    let mut e = sample();
    e.spectra.push(e.spectra[0].clone());
    e.spectra[1].native_id = "scan=2".into();
    sql(
        &p,
        "CREATE TRIGGER fail_late BEFORE INSERT ON SPECTRUM WHEN NEW.ID=1 BEGIN SELECT RAISE(ABORT,'injected'); END;",
    );
    assert!(h.write_experiment(&e).is_err());
    assert_eq!(h.nr_spectra().unwrap(), 0);
    assert_eq!(h.nr_chromatograms().unwrap(), 0);
    assert!(h.run_id().is_err());
    assert!(h.write_spectra(&e.spectra).is_err());
    assert_eq!(h.nr_spectra().unwrap(), 0);
    sql(&p, "DROP TRIGGER fail_late;");
    h.write_spectra(&e.spectra).unwrap();
    assert_eq!(h.read_spectra(&[0, 1], true).unwrap().len(), 2);
    h.create_tables().unwrap();
    sql(
        &p,
        "CREATE TRIGGER fail_data BEFORE INSERT ON DATA WHEN NEW.CHROMATOGRAM_ID=1 BEGIN SELECT RAISE(ABORT,'injected'); END;",
    );
    let c = vec![e.chromatograms[0].clone(); 2];
    assert!(h.write_chromatograms(&c).is_err());
    assert_eq!(h.nr_chromatograms().unwrap(), 0);
    sql(
        &p,
        "DROP TRIGGER fail_data; CREATE TRIGGER fail_extra BEFORE INSERT ON RUN_EXTRA BEGIN SELECT RAISE(ABORT,'injected'); END;",
    );
    assert!(h.write_run_level_information(&e, true).is_err());
    assert!(h.run_id().is_err());
}
#[test]
fn unsupported_or_corrupt_array_encodings_are_errors() {
    for code in [0, 2, 3, 4, 7, 8, -1] {
        let (_d, p) = copied_fixture();
        sql(
            &p,
            &format!("UPDATE DATA SET COMPRESSION={code} WHERE SPECTRUM_ID=0"),
        );
        assert!(
            MzMLSqliteHandler::new(&p, 0)
                .read_spectra(&[0], false)
                .is_err()
        );
    }
    for modification in [
        "UPDATE DATA SET DATA=X'789c' WHERE SPECTRUM_ID=0",
        "DELETE FROM DATA WHERE SPECTRUM_ID=0 AND DATA_TYPE=1",
        "UPDATE DATA SET DATA_TYPE=0 WHERE SPECTRUM_ID=0",
        "UPDATE DATA SET DATA=(SELECT DATA FROM DATA WHERE CHROMATOGRAM_ID=0 AND DATA_TYPE=1) WHERE SPECTRUM_ID=0 AND DATA_TYPE=1",
    ] {
        let (_d, p) = copied_fixture();
        sql(&p, modification);
        assert!(
            MzMLSqliteHandler::new(&p, 0)
                .read_spectra(&[0], false)
                .is_err(),
            "{modification}"
        );
    }
}
#[test]
fn invalid_ownership_run_and_activation_are_errors() {
    for modification in [
        "UPDATE DATA SET CHROMATOGRAM_ID=0 WHERE SPECTRUM_ID=0",
        "UPDATE DATA SET SPECTRUM_ID=999 WHERE SPECTRUM_ID=0",
        "UPDATE PRECURSOR SET ACTIVATION_METHOD=-2",
        "DELETE FROM RUN",
        "INSERT INTO RUN VALUES(5,'x','x')",
        "UPDATE RUN_EXTRA SET RUN_ID=5",
    ] {
        let (_d, p) = copied_fixture();
        sql(&p, modification);
        assert!(
            MzMLSqliteHandler::new(&p, 0)
                .read_experiment(false)
                .is_err(),
            "{modification}"
        );
    }
}
#[test]
fn missing_reads_do_not_create_files_and_existing_append_is_rejected() {
    let (_d, p) = temporary();
    let h = MzMLSqliteHandler::new(&p, 0);
    assert!(h.run_id().is_err());
    assert!(h.nr_spectra().is_err());
    assert!(h.nr_chromatograms().is_err());
    assert!(h.read_experiment(true).is_err());
    assert!(h.read_spectra(&[0], true).is_err());
    assert!(!p.exists());
    let mut uncreated = MzMLSqliteHandler::new(&p, 12345);
    assert!(uncreated.write_spectra(&sample().spectra).is_err());
    assert!(
        uncreated
            .write_chromatograms(&sample().chromatograms)
            .is_err()
    );
    assert!(!p.exists());
    let mut h = writer(&p);
    h.write_spectra(&sample().spectra).unwrap();
    let mut new = MzMLSqliteHandler::new(&p, 0);
    assert!(new.write_spectra(&sample().spectra).is_err());
}
#[test]
fn loss_policy_is_explicit_and_configuration_is_atomic() {
    let (_d, p) = temporary();
    let mut h = writer(&p);
    let old = h.config();
    assert!(h.set_config(false, true, f64::NAN, 500).is_err());
    assert!(h.set_config(false, true, 0.0001, 0).is_err());
    assert_eq!(h.config(), old);
    h.set_config(true, true, -1.0, 500).unwrap();
    h.write_experiment(&sample()).unwrap();
    compare_data(&h.read_experiment(false).unwrap(), &sample());
    h.create_tables().unwrap();
    let mut e = sample();
    e.spectra[0]
        .float_data_arrays
        .push(DataArray::new("aux", vec![1.0, 2.0]));
    assert!(h.write_experiment(&e).is_err());
    assert_eq!(h.nr_spectra().unwrap(), 0);
    h.set_loss_policy(LossPolicy::Source);
    h.write_experiment(&e).unwrap();
    assert!(
        h.read_experiment(false).unwrap().spectra[0]
            .float_data_arrays
            .is_empty()
    );
}
#[test]
fn limits_views_and_sidecars_fail_without_mutation() {
    let (_d, p) = copied_fixture();
    let mut h = MzMLSqliteHandler::new(&p, 0);
    h.limits.max_records = 1;
    assert!(h.read_experiment(false).is_err());
    h.limits.max_records = 100;
    h.limits.max_blob_bytes = 10;
    assert!(h.read_spectra(&[0], false).is_err());
    let (_d, p) = temporary();
    let mut h = writer(&p);
    sql(
        &p,
        "DROP TABLE DATA;CREATE VIEW DATA AS SELECT NULL AS SPECTRUM_ID,NULL AS CHROMATOGRAM_ID,1 AS COMPRESSION,0 AS DATA_TYPE,X'' AS DATA;",
    );
    assert!(h.nr_spectra().is_err());
    for suffix in ["-wal", "-shm", "-journal"] {
        let side = PathBuf::from(format!("{}{suffix}", p.display()));
        std::fs::write(&side, b"stale").unwrap();
        assert!(h.create_tables().is_err());
        assert!(side.exists());
        std::fs::remove_file(side).unwrap();
    }
}
#[test]
fn source_loss_policy_preserves_sql_text_nul_and_punctuation() {
    let (_d, p) = temporary();
    let mut h = writer(&p);
    h.set_loss_policy(LossPolicy::Source);
    let mut e = sample();
    e.spectra[0].native_id = "scan'\0;DROP TABLE DATA;--".into();
    h.write_spectra(&e.spectra).unwrap();
    assert_eq!(
        h.read_spectra(&[0], true).unwrap()[0].native_id,
        e.spectra[0].native_id
    );
}

#[test]
fn lossless_raw_blob_limit_is_checked_before_write() {
    let (_d, p) = temporary();
    let mut h = writer(&p);
    h.limits.max_blob_bytes = 16;
    let mut e = sample();
    e.spectra[0].peaks = vec![
        Peak1D {
            mz: 0.0,
            intensity: 0.0
        };
        3
    ];
    assert!(h.write_spectra(&e.spectra).is_err());
    assert_eq!(h.nr_spectra().unwrap(), 0);
}
#[test]
fn sampled_noise_metadata_grids_survive_full_snapshot() {
    let (_d, p) = temporary();
    let mut h = writer(&p);
    let mut e = sample();
    for (name, data) in [
        ("sampled noise m/z array", vec![100.0, 200.0]),
        ("sampled noise intensity array", vec![1.0, 2.0]),
        ("sampled noise baseline array", vec![0.5, 0.75]),
    ] {
        e.spectra[0].metadata.insert(
            name.into(),
            openms::metadata::MetaValue::try_from(data).unwrap(),
        );
    }
    h.write_experiment(&e).unwrap();
    let out = h.read_experiment(false).unwrap();
    assert_eq!(out.spectra[0].metadata, e.spectra[0].metadata);
}
#[test]
fn run_id_relationship_and_snapshot_acquisition_mismatch_are_rejected() {
    for statement in [
        "UPDATE SPECTRUM SET RUN_ID=999",
        "UPDATE CHROMATOGRAM SET RUN_ID=NULL",
        "UPDATE SPECTRUM SET MSLEVEL=5",
        "UPDATE SPECTRUM SET RETENTION_TIME=123.0",
    ] {
        let (_d, p) = copied_fixture();
        sql(&p, statement);
        assert!(MzMLSqliteHandler::new(&p, 0).read_experiment(true).is_err());
    }
}
#[test]
fn short_lossy_coordinate_arrays_use_readable_lossless_fallback() {
    for values in [
        vec![],
        vec![0.0],
        vec![0.0, 0.0],
        vec![100.0],
        vec![100.0, 101.0],
    ] {
        let (_d, p) = temporary();
        let mut h = writer(&p);
        h.set_config(true, true, 0.0001, 500).unwrap();
        let mut e = sample();
        e.spectra[0].peaks = values
            .iter()
            .map(|&mz| Peak1D { mz, intensity: 1.0 })
            .collect();
        e.chromatograms[0].peaks = values
            .iter()
            .map(|&rt| ChromatogramPeak { rt, intensity: 1.0 })
            .collect();
        h.write_experiment(&e).unwrap();
        let out = h.read_experiment(false).unwrap();
        compare_data(&out, &e);
        let db = rusqlite::Connection::open(&p).unwrap();
        let code: i64 = db
            .query_row(
                "SELECT COMPRESSION FROM DATA WHERE SPECTRUM_ID=0 AND DATA_TYPE=0",
                [],
                |r| r.get(0),
            )
            .unwrap();
        assert_eq!(code, 1);
    }
}

#[test]
fn rich_metadata_builder_matches_snapshot_and_rejects_disagreement() {
    let (_d, p) = temporary();
    let mut h = writer(&p);
    let mut e = sample();
    e.spectra[0].name = "descriptive scan".into();
    e.chromatograms[0].name = "descriptive trace".into();
    h.write_run_level_information(&e, true).unwrap();
    let mut wrong = e.spectra.clone();
    wrong[0].name = "different".into();
    assert!(h.write_spectra(&wrong).is_err());
    assert_eq!(h.nr_spectra().unwrap(), 0);
    h.write_spectra(&e.spectra).unwrap();
    h.write_chromatograms(&e.chromatograms).unwrap();
    let out = h.read_experiment(false).unwrap();
    assert_eq!(out.spectra[0].name, e.spectra[0].name);
    assert_eq!(out.chromatograms[0].name, e.chromatograms[0].name);
}
#[test]
fn reserved_legacy_run_key_cannot_conflict_with_native_owner() {
    let (_d, p) = temporary();
    let mut h = writer(&p);
    let mut e = sample();
    e.settings.metadata.insert(
        "sqMassRunID".into(),
        openms::metadata::MetaValue::from(999_i64),
    );
    e.settings.metadata.insert(
        "scientific annotation".into(),
        openms::metadata::MetaValue::from("keep"),
    );
    h.write_experiment(&e).unwrap();
    let out = h.read_experiment(false).unwrap();
    assert_eq!(out.sql_run_id, 12345);
    assert!(!out.settings.metadata.contains_key("sqMassRunID"));
    assert_eq!(
        out.settings.metadata["scientific annotation"],
        e.settings.metadata["scientific annotation"]
    );
}
