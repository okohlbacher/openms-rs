// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/sqmass_file.rs`, the port of `FORMAT/SqMassFile.h`.
//!
//! Every `START_SECTION` of the upstream `SqMassFile_test.cpp` (core `bc9cc12`)
//! is reproduced with its tolerances. The fixtures are the class-test files
//! `SqliteMassFile_1.sqMass` and `MzMLSqliteHandler_1.mzML`, which this
//! repository already holds (`tests/data/sqlite_s1_source_review/` and
//! `tests/data/mzml_sqlite_handler/`, byte-identical to core `bc9cc12`).
#![cfg(feature = "sqmass")]

use openms::format::mzml_sqlite_handler::LossPolicy;
use openms::format::sqmass_file::{SqMassConfig, SqMassFile};
use openms::interfaces::MSDataConsumer;
use openms::kernel::{MSChromatogram, MSExperiment, MSSpectrum, Precursor};
use openms::metadata::{ActivationMethod, DriftTimeUnit, ExperimentalSettings, MetaValue, Product};
use openms::system::file::TempDir;
use std::ops::ControlFlow;
use std::path::{Path, PathBuf};

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}
fn sqmass_1() -> PathBuf {
    root().join("tests/data/sqlite_s1_source_review/SqliteMassFile_1.sqMass")
}
fn mzml_1() -> MSExperiment {
    openms::format::mzml::load(
        root().join("tests/data/mzml_sqlite_handler/MzMLSqliteHandler_1.mzML"),
    )
    .unwrap()
}

/// The source keeps the sqMass run ID as the experiment's `sqMassRunID` meta
/// value, so its `ExperimentalSettings` comparison and `getSqlRunID` see it
/// there in both experiments. This crate's sqMass reader moves it into
/// `MSExperiment::sql_run_id`, and its mzML reader leaves it a meta value;
/// this moves it the same way in an mzML-loaded experiment, which is the
/// representation change and nothing more.
fn with_run_id_field(mut exp: MSExperiment) -> MSExperiment {
    if let Some(value) = exp.settings.metadata.remove("sqMassRunID") {
        match value.data() {
            openms::metadata::MetaValueData::Integer(id) => exp.sql_run_id = *id as u64,
            other => panic!("unexpected sqMassRunID {other:?}"),
        }
    }
    exp
}

/// The source's `TEST::isRealSimilar` under explicit tolerances: the absolute
/// criterion first, then the ratio.
fn similar(a: f64, b: f64, abs: f64, rel: f64) -> bool {
    if (a - b).abs() <= abs {
        return true;
    }
    if a == 0.0 || b == 0.0 {
        return false;
    }
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    ratio <= rel && ratio >= 1.0 / rel
}

fn cmp_intensity(a: &MSExperiment, b: &MSExperiment, abs: f64, rel: f64) {
    for (s1, s2) in a.spectra.iter().zip(&b.spectra) {
        assert_eq!(s1.peaks.len(), s2.peaks.len());
        for (p, q) in s1.peaks.iter().zip(&s2.peaks) {
            assert!(
                similar(f64::from(p.intensity), f64::from(q.intensity), abs, rel),
                "{} {}",
                p.intensity,
                q.intensity
            );
        }
    }
    for (c1, c2) in a.chromatograms.iter().zip(&b.chromatograms) {
        assert_eq!(c1.peaks.len(), c2.peaks.len());
        for (p, q) in c1.peaks.iter().zip(&c2.peaks) {
            assert!(
                similar(f64::from(p.intensity), f64::from(q.intensity), abs, rel),
                "{} {}",
                p.intensity,
                q.intensity
            );
        }
    }
}

fn cmp_mz(a: &MSExperiment, b: &MSExperiment, abs: f64, rel: f64) {
    for (s1, s2) in a.spectra.iter().zip(&b.spectra) {
        assert_eq!(s1.peaks.len(), s2.peaks.len());
        for (p, q) in s1.peaks.iter().zip(&s2.peaks) {
            assert!(similar(p.mz, q.mz, abs, rel), "{} {}", p.mz, q.mz);
        }
    }
}

fn cmp_rt(a: &MSExperiment, b: &MSExperiment, abs: f64, rel: f64) {
    for (c1, c2) in a.chromatograms.iter().zip(&b.chromatograms) {
        assert_eq!(c1.peaks.len(), c2.peaks.len());
        for (p, q) in c1.peaks.iter().zip(&c2.peaks) {
            assert!(similar(p.rt, q.rt, abs, rel), "{} {}", p.rt, q.rt);
        }
    }
}

// START_SECTION((SqMassFile())) / START_SECTION((~SqMassFile()))
#[test]
fn construction() {
    let file = SqMassFile::new();
    assert_eq!(file.config(), SqMassConfig::default());
    assert!(file.config().write_full_meta);
    assert!(!file.config().use_lossy_numpress);
    assert_eq!(file.config().linear_fp_mass_acc, -1.0);
}

// START_SECTION(void load(const std::string& filename, MapType& map))
#[test]
fn load_matches_the_source_section() {
    let exp = SqMassFile::new().load(sqmass_1()).unwrap();
    let exp2 = mzml_1();
    assert_eq!(exp.spectra.len(), exp2.spectra.len());
    assert_eq!(exp.chromatograms.len(), exp2.chromatograms.len());
    assert_eq!(exp.spectra.len(), 2);
    assert_eq!(exp.chromatograms.len(), 1);
    assert_ne!(exp.spectra[0], exp2.spectra[0]); // no exact duplicate

    cmp_intensity(&exp, &exp2, 1e-4, 1.001);
    cmp_mz(&exp, &exp2, 1e-5, 1.000001); // less than 1 ppm error for m/z
    cmp_rt(&exp, &exp2, 0.05, 1.000001); // max 0.05 seconds error in RT

    // mapping of experimental settings ...
    let exp2 = with_run_id_field(exp2);
    assert_eq!(exp.settings, exp2.settings);
    assert_eq!(exp.sql_run_id, exp2.sql_run_id);
    assert_eq!(exp.sql_run_id, 12345);
}

// START_SECTION(void store(const std::string& filename, MapType& map))
#[test]
fn store_matches_the_source_section() {
    let exp_orig = mzml_1();
    let config = SqMassConfig {
        use_lossy_numpress: false,
        linear_fp_mass_acc: -1.0,
        write_full_meta: false,
        // The section asserts that the experimental settings are lost.
        loss_policy: LossPolicy::Source,
    };
    let mut file = SqMassFile::new();
    file.set_config(config);
    let dir = TempDir::new(false).unwrap();
    let tmp_filename = dir.path().join("store.sqMass");
    file.store(&tmp_filename, &exp_orig).unwrap();
    let exp = file.load(&tmp_filename).unwrap();

    assert_eq!(exp.spectra.len(), exp_orig.spectra.len());
    assert_eq!(exp.chromatograms.len(), exp_orig.chromatograms.len());
    assert_eq!(exp.spectra.len(), 2);
    assert_eq!(exp.chromatograms.len(), 1);
    assert_ne!(exp.spectra[0], exp_orig.spectra[0]); // no exact duplicate

    let exp2 = exp_orig.clone();
    // since we specified no lossy compression, we expect high accuracy
    cmp_intensity(&exp, &exp2, 1e-8, 1.00000001);
    cmp_mz(&exp, &exp2, 1e-8, 1.00000001);
    cmp_rt(&exp, &exp2, 1e-8, 1.00000001);

    // no 1:1 mapping of experimental settings ...
    assert_ne!(exp.settings, exp2.settings);

    // Native: the default loss policy refuses the same store.
    let mut refusing = SqMassFile::new();
    refusing.set_config(SqMassConfig {
        loss_policy: LossPolicy::Reject,
        ..config
    });
    assert!(
        refusing
            .store(dir.path().join("refused.sqMass"), &exp_orig)
            .is_err()
    );
}

/// The precursor and product both lossy sections attach.
fn annotated() -> MSExperiment {
    let mut exp_orig = mzml_1();
    let mut p = Precursor::new(600.0, 4);
    p.activation_methods.insert(ActivationMethod::Bird);
    p.activation_energy = 500.0;
    p.isolation_window_upper_offset = 7.0;
    p.isolation_window_lower_offset = 14.0;
    p.drift_time = Some(0.5);
    p.drift_time_unit = DriftTimeUnit::Millisecond;
    p.cv_terms
        .metadata
        .insert("peptide_sequence".into(), MetaValue::from("PEPTIDEK"));
    exp_orig.spectra[0].precursors = vec![p.clone()];
    exp_orig.chromatograms[0].precursor = p;

    let pr = Product {
        mz: 300.0,
        isolation_window_upper_offset: 10.0,
        isolation_window_lower_offset: 15.0,
        ..Product::default()
    };
    exp_orig.spectra[0].products = vec![pr.clone()];
    exp_orig.chromatograms[0].product = pr;
    assert!(similar(
        exp_orig.spectra[0].precursors[0].activation_energy,
        500.0,
        1e-5,
        1.0 + 1e-5
    ));
    exp_orig
}

fn check_annotations(exp: &MSExperiment, exp_orig: &MSExperiment) {
    let s = &exp.spectra[0];
    let t = |a: f64, b: f64| assert!(similar(a, b, 1e-5, 1.0 + 1e-5), "{a} {b}");
    t(s.rt, exp_orig.spectra[0].rt);
    assert_eq!(s.native_id, exp_orig.spectra[0].native_id);
    assert_eq!(s.ms_level, exp_orig.spectra[0].ms_level);
    assert_eq!(
        s.instrument_settings.polarity,
        exp_orig.spectra[0].instrument_settings.polarity
    );
    assert_eq!(s.products.len(), 1);
    t(s.products[0].mz, 300.0);
    t(s.products[0].isolation_window_upper_offset, 10.0);
    t(s.products[0].isolation_window_lower_offset, 15.0);
    assert_eq!(s.precursors.len(), 1);
    let p = &s.precursors[0];
    t(p.activation_energy, 500.0);
    assert_eq!(p.charge, 4);
    t(p.drift_time.unwrap(), 0.5);
    t(p.mz, 600.0);
    t(p.isolation_window_upper_offset, 7.0);
    t(p.isolation_window_lower_offset, 14.0);
    assert_eq!(
        p.cv_terms.metadata["peptide_sequence"].to_string(),
        "PEPTIDEK"
    );

    let c: &MSChromatogram = &exp.chromatograms[0];
    assert_eq!(c.native_id, exp_orig.chromatograms[0].native_id);
    t(c.product.mz, 300.0);
    t(c.product.isolation_window_upper_offset, 10.0);
    t(c.product.isolation_window_lower_offset, 15.0);
    t(c.precursor.activation_energy, 500.0);
    assert_eq!(c.precursor.charge, 4);
    t(c.precursor.drift_time.unwrap(), 0.5);
    t(c.precursor.mz, 600.0);
    t(c.precursor.isolation_window_upper_offset, 7.0);
    t(c.precursor.isolation_window_lower_offset, 14.0);
    assert_eq!(
        c.precursor.cv_terms.metadata["peptide_sequence"].to_string(),
        "PEPTIDEK"
    );
}

// START_SECTION([EXTRA_LOSSY] void store(const std::string& filename, MapType& map))
#[test]
fn lossy_store_matches_the_source_section() {
    let exp_orig = annotated();
    let mut file = SqMassFile::new();
    file.set_config(SqMassConfig {
        use_lossy_numpress: true,
        linear_fp_mass_acc: 0.0001,
        write_full_meta: false,
        // The section asserts that the experimental settings are lost.
        loss_policy: LossPolicy::Source,
    });
    let dir = TempDir::new(false).unwrap();
    let tmp_filename = dir.path().join("lossy.sqMass");
    file.store(&tmp_filename, &exp_orig).unwrap();
    let exp = file.load(&tmp_filename).unwrap();

    assert_eq!(exp.spectra.len(), exp_orig.spectra.len());
    assert_eq!(exp.chromatograms.len(), exp_orig.chromatograms.len());
    assert_eq!(exp.spectra.len(), 2);
    assert_eq!(exp.chromatograms.len(), 1);
    assert_ne!(exp.spectra[0], exp_orig.spectra[0]);
    check_annotations(&exp, &exp_orig);

    // should not give 1:1 mapping of experimental settings ...
    assert_ne!(exp.settings, exp_orig.settings);
    cmp_intensity(&exp, &exp_orig, 1e-4, 1.001); // 0.1 % error for intensity
    cmp_mz(&exp, &exp_orig, 1e-5, 1.000001); // less than 1 ppm error for m/z
    cmp_rt(&exp, &exp_orig, 0.05, 1.00000001); // max 0.05 seconds error in RT
}

// START_SECTION([EXTRA_FULL_META] void store(const std::string& filename, MapType& map))
#[test]
fn full_meta_store_matches_the_source_section() {
    let exp_orig = annotated();
    let mut file = SqMassFile::new();
    file.set_config(SqMassConfig {
        use_lossy_numpress: true,
        linear_fp_mass_acc: 0.0001,
        write_full_meta: true,
        loss_policy: LossPolicy::Reject,
    });
    let dir = TempDir::new(false).unwrap();
    let tmp_filename = dir.path().join("full.sqMass");
    file.store(&tmp_filename, &exp_orig).unwrap();
    let exp = file.load(&tmp_filename).unwrap();

    assert_eq!(exp.spectra.len(), exp_orig.spectra.len());
    assert_eq!(exp.chromatograms.len(), exp_orig.chromatograms.len());
    assert_eq!(exp.spectra.len(), 2);
    assert_eq!(exp.chromatograms.len(), 1);
    assert_ne!(exp.spectra[0], exp_orig.spectra[0]);
    check_annotations(&exp, &exp_orig);

    // using full meta should give 1:1 mapping of experimental settings ...
    let exp2 = with_run_id_field(exp_orig.clone());
    let settings: &ExperimentalSettings = &exp.settings;
    assert_eq!(settings, &exp2.settings);
    assert_eq!(exp.sql_run_id, exp2.sql_run_id);
    assert_eq!(exp.sql_run_id, 12345);
    cmp_intensity(&exp, &exp_orig, 1e-4, 1.001);
    cmp_mz(&exp, &exp_orig, 1e-5, 1.000001);
    cmp_rt(&exp, &exp_orig, 0.05, 1.00000001);
}

// ----------------------------------------------------------------------------
// transform (no upstream section)
// ----------------------------------------------------------------------------

#[derive(Default)]
struct Collect {
    expected: Option<(usize, usize)>,
    settings: Option<ExperimentalSettings>,
    spectra: Vec<MSSpectrum>,
    chromatograms: Vec<MSChromatogram>,
    stop_after: Option<usize>,
}

impl MSDataConsumer for Collect {
    fn set_expected_size(&mut self, spectra: usize, chromatograms: usize) -> openms::Result<()> {
        self.expected = Some((spectra, chromatograms));
        Ok(())
    }
    fn set_experimental_settings(&mut self, settings: &ExperimentalSettings) -> openms::Result<()> {
        self.settings = Some(settings.clone());
        Ok(())
    }
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> openms::Result<ControlFlow<()>> {
        self.spectra.push(spectrum.clone());
        Ok(match self.stop_after {
            Some(n) if self.spectra.len() >= n => ControlFlow::Break(()),
            _ => ControlFlow::Continue(()),
        })
    }
    fn consume_chromatogram(
        &mut self,
        chromatogram: &mut MSChromatogram,
    ) -> openms::Result<ControlFlow<()>> {
        self.chromatograms.push(chromatogram.clone());
        Ok(ControlFlow::Continue(()))
    }
}

#[test]
fn transform_streams_what_load_reads() {
    let file = SqMassFile::new();
    let loaded = file.load(sqmass_1()).unwrap();
    let mut consumer = Collect::default();
    file.transform(sqmass_1(), &mut consumer, false, false)
        .unwrap();
    assert_eq!(consumer.expected, Some((2, 1)));
    assert_eq!(consumer.settings.as_ref(), Some(&loaded.settings));
    // The streamed records carry the SQL metadata only, as the source's
    // readSpectra/readChromatograms do; load also applies the RUN_EXTRA
    // snapshot. The data and native IDs agree.
    assert_eq!(consumer.spectra.len(), loaded.spectra.len());
    for (streamed, whole) in consumer.spectra.iter().zip(&loaded.spectra) {
        assert_eq!(streamed.native_id, whole.native_id);
        assert_eq!(streamed.rt, whole.rt);
        assert_eq!(streamed.peaks, whole.peaks);
    }
    assert_eq!(consumer.chromatograms.len(), loaded.chromatograms.len());
    for (streamed, whole) in consumer.chromatograms.iter().zip(&loaded.chromatograms) {
        assert_eq!(streamed.native_id, whole.native_id);
        assert_eq!(streamed.peaks, whole.peaks);
    }

    let mut early = Collect {
        stop_after: Some(1),
        ..Collect::default()
    };
    file.transform(sqmass_1(), &mut early, true, true).unwrap();
    assert_eq!(early.spectra.len(), 1);
    assert!(early.chromatograms.is_empty());
}

#[test]
fn transform_handles_a_multiple_of_the_batch_size() {
    // 500 spectra: the source issues an empty second batch and fails there.
    let mut exp = MSExperiment::new();
    for i in 0..500 {
        exp.spectra.push(MSSpectrum {
            native_id: format!("scan={i}"),
            rt: i as f64,
            ..MSSpectrum::default()
        });
    }
    let mut file = SqMassFile::new();
    file.set_config(SqMassConfig {
        write_full_meta: false,
        loss_policy: LossPolicy::Source,
        ..SqMassConfig::default()
    });
    let dir = TempDir::new(false).unwrap();
    let path = dir.path().join("batch.sqMass");
    file.store(&path, &exp).unwrap();
    let mut consumer = Collect::default();
    file.transform(&path, &mut consumer, false, false).unwrap();
    assert_eq!(consumer.spectra.len(), 500);
    assert_eq!(consumer.spectra[499].native_id, "scan=499");
}
