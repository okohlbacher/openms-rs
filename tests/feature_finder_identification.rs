// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FeatureFinderIdentificationAlgorithm` and `FeatureFinderAlgorithmMetaboIdent`:
//! the sections of the upstream class test, the parameters, and the behaviour
//! that does not need a retained file, on small synthetic MS1 data.
//!
//! Evidence. The three sections of
//! `FeatureFinderIdentificationAlgorithm_test.cpp` (core `bc9cc12`, tier 3) are
//! transcribed below; `FeatureFinderAlgorithmMetaboIdent` has no class test
//! upstream. The remaining tests are tier 4: values derived from the source's
//! documented rules (stated at each test), not from this port's output. The
//! executed comparisons against the C++ Release build are in
//! `tests/feature_finder_identification_oracle.rs` and
//! `tests/feature_finder_identification_topp.rs`.

use openms::Error;
use openms::analysis::featurefinder::feature_finder_identification::FeatureFinderIdentificationAlgorithm;
use openms::analysis::featurefinder::feature_finder_metabo_ident::{
    FeatureFinderAlgorithmMetaboIdent, FeatureFinderMetaboIdentCompound, MAX_LIBRARY_ENTRIES,
};
use openms::chemistry::AASequence;
use openms::concept::UniqueIdGenerator;
use openms::constants::C13C12_MASSDIFF_U;
use openms::identification::{PeptideHit, PeptideIdentification, ProteinIdentification};
use openms::kernel::{Feature, FeatureMap, MSExperiment, MSSpectrum, Peak1D};
use openms::metadata::MetaValue;
use openms::param::ParamValue;

// ---------------------------------------------------------------- data

/// Peptide used throughout: its doubly protonated ion.
const SEQUENCE: &str = "PEPTIDEK";

fn precursor_mz() -> f64 {
    AASequence::parse(SEQUENCE).unwrap().mz(2).unwrap()
}

/// One MS1 spectrum per `step` seconds in `0..=span`. The three isotopes of the
/// doubly charged ion at `mz` elute as a Gaussian of standard deviation `sigma`
/// around `apex`; `shift` moves all of them in m/z.
fn experiment(mz: f64, apex: f64, sigma: f64) -> MSExperiment {
    let ratios = [1.0f64, 0.55, 0.2];
    let mut exp = MSExperiment::new();
    let mut rt = 100.0;
    while rt <= 300.0 {
        let gauss = (-0.5 * ((rt - apex) / sigma).powi(2)).exp();
        let mut spectrum = MSSpectrum {
            rt,
            ms_level: 1,
            native_id: format!("scan={}", exp.spectra.len() + 1),
            ..MSSpectrum::default()
        };
        for (i, ratio) in ratios.iter().enumerate() {
            let isotope_mz = mz + C13C12_MASSDIFF_U * i as f64 / 2.0;
            let intensity = (20_000.0 * ratio * gauss) as f32 + 50.0;
            spectrum.peaks.push(Peak1D::new(isotope_mz, intensity));
        }
        // a peak elsewhere, so that the spectra are never empty
        spectrum.peaks.push(Peak1D::new(900.0, 80.0));
        exp.spectra.push(spectrum);
        rt += 2.0;
    }
    exp
}

fn identification(rt: f64, mz: f64, sequence: &str, charge: i32) -> PeptideIdentification {
    let hit = PeptideHit::new(0.01, 0, charge, AASequence::parse(sequence).unwrap()).unwrap();
    PeptideIdentification {
        identifier: "run".to_owned(),
        hits: vec![hit],
        score_type: "q-value".to_owned(),
        higher_score_better: false,
        rt: Some(rt),
        mz: Some(mz),
        ..PeptideIdentification::default()
    }
}

fn proteins() -> Vec<ProteinIdentification> {
    vec![ProteinIdentification {
        identifier: "run".to_owned(),
        ..ProteinIdentification::default()
    }]
}

fn parameters(extra: &[(&str, ParamValue)]) -> openms::param::Param {
    let mut param = FeatureFinderIdentificationAlgorithm::defaults().unwrap();
    param
        .set_value("extract:mz_window", ParamValue::Float(0.1), "", &[])
        .unwrap();
    param
        .set_value("detect:peak_width", ParamValue::Float(30.0), "", &[])
        .unwrap();
    param
        .set_value("model:type", ParamValue::String("none".into()), "", &[])
        .unwrap();
    for (key, value) in extra {
        param.set_value(key, value.clone(), "", &[]).unwrap();
    }
    param
}

fn ffid(extra: &[(&str, ParamValue)]) -> FeatureFinderIdentificationAlgorithm {
    let mut algorithm =
        FeatureFinderIdentificationAlgorithm::from_param(&parameters(extra)).unwrap();
    algorithm.set_unique_id_generator(UniqueIdGenerator::from_seed(7));
    algorithm
}

// ---------------------------------------------------------------- class test

/// `START_SECTION(FeatureFinderIdentificationAlgorithm())` and
/// `START_SECTION(~FeatureFinderIdentificationAlgorithm())`.
#[test]
fn class_test_constructor_and_destructor() {
    let algorithm = FeatureFinderIdentificationAlgorithm::new().unwrap();
    drop(algorithm);
}

/// `START_SECTION([EXTRA] summary statistics cannot wrap a missing-peptide
/// count)`: a feature with positive intensity, one identification with one hit,
/// and a peptide count of zero must report "0 peptides without features".
#[test]
fn class_test_summary_statistics_cannot_wrap_a_missing_peptide_count() {
    let hit = PeptideHit::new(0.0, 0, 0, AASequence::parse("PEPTIDEK").unwrap()).unwrap();
    let peptide = PeptideIdentification {
        hits: vec![hit],
        ..PeptideIdentification::default()
    };
    let mut feature = Feature::new(0.0, 0.0, 1.0);
    feature.peptide_identifications = vec![peptide];
    let mut features = FeatureMap::new();
    features.features.push(feature);

    let mut algorithm = FeatureFinderIdentificationAlgorithm::new().unwrap();
    algorithm.set_peptide_count(0);
    let statistics = algorithm.statistics(&features).unwrap();
    assert!(statistics.text.contains("0 peptides without features"));
    assert_eq!(statistics.missing, 0);
    assert_eq!(statistics.quantified, 1);
}

/// `START_SECTION([EXTRA] run() rejects MS data without spectra)`.
#[test]
fn class_test_run_rejects_ms_data_without_spectra() {
    let mut algorithm = FeatureFinderIdentificationAlgorithm::new().unwrap();
    assert!(algorithm.ms_data().spectra.is_empty());
    let mut peptide = identification(100.0, 500.0, SEQUENCE, 2);
    peptide.identifier = "run0".to_owned();
    let proteins = vec![ProteinIdentification {
        identifier: "run0".to_owned(),
        ..ProteinIdentification::default()
    }];
    let mut features = FeatureMap::new();
    let error = algorithm
        .run(
            vec![peptide],
            &proteins,
            &mut features,
            &FeatureMap::new(),
            "",
        )
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    assert!(features.features.is_empty());
}

/// An unquantified peptide is never counted as missing below zero, and the
/// quantified ones are distinct sequences.
#[test]
fn statistics_count_distinct_sequences_once() {
    let make = |sequence: &str| {
        let hit = PeptideHit::new(0.0, 0, 2, AASequence::parse(sequence).unwrap()).unwrap();
        let mut feature = Feature::new(0.0, 0.0, 1.0);
        feature.peptide_identifications = vec![PeptideIdentification {
            hits: vec![hit],
            ..PeptideIdentification::default()
        }];
        feature
    };
    let mut features = FeatureMap::new();
    features.features.push(make("PEPTIDEK"));
    features.features.push(make("PEPTIDEK"));
    features.features.push(make("AAAAK"));
    let mut algorithm = FeatureFinderIdentificationAlgorithm::new().unwrap();
    algorithm.set_peptide_count(5);
    let statistics = algorithm.statistics(&features).unwrap();
    assert_eq!(
        (
            statistics.identified,
            statistics.quantified,
            statistics.missing
        ),
        (5, 2, 3)
    );
}

// ---------------------------------------------------------------- parameters

#[test]
fn ffid_defaults_are_the_sources() {
    let defaults = FeatureFinderIdentificationAlgorithm::defaults().unwrap();
    let value = |key: &str| defaults.value(key).unwrap().clone();
    assert_eq!(value("extract:batch_size"), ParamValue::Integer(5000));
    assert_eq!(value("extract:mz_window"), ParamValue::Float(10.0));
    assert_eq!(value("extract:IM_window"), ParamValue::Float(0.06));
    assert_eq!(value("extract:n_isotopes"), ParamValue::Integer(2));
    assert_eq!(value("extract:isotope_pmin"), ParamValue::Float(0.0));
    assert_eq!(value("extract:rt_window"), ParamValue::Float(0.0));
    assert_eq!(value("detect:peak_width"), ParamValue::Float(60.0));
    assert_eq!(value("detect:min_peak_width"), ParamValue::Float(0.2));
    assert_eq!(value("detect:signal_to_noise"), ParamValue::Float(0.8));
    assert_eq!(value("detect:mapping_tolerance"), ParamValue::Float(0.0));
    assert_eq!(value("quantify_decoys"), ParamValue::String("false".into()));
    assert_eq!(value("min_psm_cutoff"), ParamValue::String("none".into()));
    assert_eq!(value("add_mass_offset_peptides"), ParamValue::Float(0.0));
    assert_eq!(value("seed_apex_rt_tolerance"), ParamValue::Float(5.0));
    assert_eq!(value("model:type"), ParamValue::String("symmetric".into()));
    assert_eq!(value("EMGScoring:max_iteration"), ParamValue::Integer(100));
    assert_eq!(
        value("EMGScoring:init_mom"),
        ParamValue::String("true".into())
    );
    assert_eq!(
        value("faims:merge_features"),
        ParamValue::String("true".into())
    );
    assert_eq!(value("debug"), ParamValue::Integer(0));
    // the model section is the elution model fitter's without `asymmetric`
    assert!(defaults.value("model:add_zeros").is_ok());
    assert!(defaults.value("model:check:min_area").is_ok());
    assert!(defaults.value("model:asymmetric").is_err());
}

#[test]
fn metabo_ident_defaults_are_the_sources() {
    let defaults = FeatureFinderAlgorithmMetaboIdent::defaults().unwrap();
    let value = |key: &str| defaults.value(key).unwrap().clone();
    assert_eq!(value("extract:mz_window"), ParamValue::Float(10.0));
    assert_eq!(value("extract:rt_window"), ParamValue::Float(0.0));
    assert_eq!(value("extract:im_window"), ParamValue::Float(0.0));
    assert_eq!(value("extract:n_isotopes"), ParamValue::Integer(2));
    assert_eq!(value("extract:isotope_pmin"), ParamValue::Float(0.0));
    assert_eq!(value("detect:peak_width"), ParamValue::Float(60.0));
    assert_eq!(value("detect:min_peak_width"), ParamValue::Float(0.2));
    assert_eq!(value("detect:signal_to_noise"), ParamValue::Float(0.8));
    assert_eq!(value("model:type"), ParamValue::String("symmetric".into()));
    assert!(defaults.value("model:asymmetric").is_err());
    assert_eq!(value("debug"), ParamValue::Integer(0));
    assert_eq!(
        value("faims:merge_features"),
        ParamValue::String("true".into())
    );
}

#[test]
fn parameters_are_validated_against_the_restrictions() {
    let mut param = FeatureFinderIdentificationAlgorithm::defaults().unwrap();
    param
        .set_value("extract:batch_size", ParamValue::Integer(0), "", &[])
        .unwrap();
    assert!(matches!(
        FeatureFinderIdentificationAlgorithm::from_param(&param),
        Err(Error::InvalidValue(_))
    ));
    let mut param = FeatureFinderIdentificationAlgorithm::defaults().unwrap();
    param
        .set_value("model:type", ParamValue::String("sideways".into()), "", &[])
        .unwrap();
    assert!(FeatureFinderIdentificationAlgorithm::from_param(&param).is_err());
    let mut param = FeatureFinderAlgorithmMetaboIdent::defaults().unwrap();
    param
        .set_value("detect:signal_to_noise", ParamValue::Float(0.0), "", &[])
        .unwrap();
    assert!(FeatureFinderAlgorithmMetaboIdent::from_param(&param).is_err());
}

#[test]
fn a_non_numeric_psm_cutoff_is_refused() {
    let mut param = FeatureFinderIdentificationAlgorithm::defaults().unwrap();
    param
        .set_value("min_psm_cutoff", ParamValue::String("high".into()), "", &[])
        .unwrap();
    assert!(FeatureFinderIdentificationAlgorithm::from_param(&param).is_err());
    param
        .set_value("min_psm_cutoff", ParamValue::String("0.05".into()), "", &[])
        .unwrap();
    assert!(FeatureFinderIdentificationAlgorithm::from_param(&param).is_ok());
}

#[test]
fn candidates_in_the_parquet_format_is_a_named_gap() {
    let mut algorithm = FeatureFinderIdentificationAlgorithm::from_param(&parameters(&[(
        "candidates_out",
        ParamValue::String("candidates.featureparquet".into()),
    )]))
    .unwrap();
    algorithm.set_ms_data(experiment(precursor_mz(), 200.0, 10.0));
    let mut features = FeatureMap::new();
    let error = algorithm
        .run(
            vec![identification(200.0, precursor_mz(), SEQUENCE, 2)],
            &proteins(),
            &mut features,
            &FeatureMap::new(),
            "",
        )
        .unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error:?}");
}

// ---------------------------------------------------------------- the run

#[test]
fn a_clean_peptide_is_found_assigned_and_annotated() {
    let mz = precursor_mz();
    let mut algorithm = ffid(&[]);
    algorithm.set_ms_data(experiment(mz, 200.0, 10.0));
    let mut features = FeatureMap::new();
    algorithm
        .run(
            vec![identification(200.0, mz, SEQUENCE, 2)],
            &proteins(),
            &mut features,
            &FeatureMap::new(),
            "scan.mzML",
        )
        .unwrap();

    // one feature for the one assay, found at the apex
    assert_eq!(features.features.len(), 1);
    let feature = &features.features[0];
    assert!((feature.rt - 200.0).abs() < 4.0, "apex at {}", feature.rt);
    assert_eq!(feature.charge, 2);
    assert!((feature.mz - mz).abs() < 1e-9);
    assert!(feature.intensity > 0.0);
    // classification: the one matching identification makes it "positive"
    let text = |key: &str| {
        feature
            .metadata
            .get(key)
            .unwrap_or_else(|| panic!("meta value {key}"))
            .as_str()
            .unwrap()
            .to_owned()
    };
    assert_eq!(text("feature_class"), "positive");
    assert_eq!(text("PeptideRef"), format!("{SEQUENCE}/2"));
    assert_eq!(
        feature
            .metadata
            .get("n_matching_ids")
            .unwrap()
            .as_i64()
            .unwrap(),
        1
    );
    assert_eq!(
        feature
            .metadata
            .get("n_total_ids")
            .unwrap()
            .as_i64()
            .unwrap(),
        1
    );
    // the identification and the proteins of the input
    assert_eq!(feature.peptide_identifications.len(), 1);
    assert_eq!(
        feature.peptide_identifications[0].hits[0]
            .sequence
            .to_string(),
        SEQUENCE
    );
    assert_eq!(features.protein_identifications.len(), 1);
    assert!(features.unassigned_peptide_identifications.is_empty());
    // one mass trace per isotope, annotated with its theoretical share
    assert_eq!(feature.subordinates.len(), 2);
    let mut total = 0.0;
    for sub in &feature.subordinates {
        total += sub
            .metadata
            .get("isotope_probability")
            .unwrap()
            .as_f64()
            .unwrap();
    }
    assert!(
        (total - 1.0).abs() < 1e-6,
        "isotope probabilities sum to {total}"
    );
    // hulls were written for every trace
    assert_eq!(feature.convex_hulls.len(), 2);
    // the run path falls back to the given file name
    assert_eq!(
        features.primary_ms_run_path().unwrap(),
        vec!["scan.mzML".to_owned()]
    );
    // the report says what the source logs
    let report = algorithm.report();
    assert_eq!(report.peptides_identified, 1);
    assert_eq!(report.peptides_quantified, 1);
    assert_eq!(report.peptides_missing, 0);
    assert_eq!(report.seeds_added, 0);
    // the library holds the assay and its two isotopes, with the protein stub
    let library = algorithm.library();
    assert_eq!(library.peptides().len(), 1);
    assert_eq!(library.transitions().len(), 2);
    assert_eq!(library.peptides()[0].id, format!("{SEQUENCE}/2"));
    assert_eq!(library.peptides()[0].protein_refs, vec!["not_available"]);
    // the spectra are released after the run, the chromatograms never kept
    assert!(algorithm.ms_data().spectra.is_empty());
    assert!(algorithm.chromatograms().chromatograms.is_empty());
}

#[test]
fn the_rt_window_follows_the_peak_width_and_the_mapping_tolerance() {
    // `extract:rt_window` 0: (2 * peak_width + tolerance) * 2, with a tolerance
    // below one relative to twice the peak width (`calculateRTWindow_`).
    let mz = precursor_mz();
    for (tolerance, expected) in [(0.0, 120.0), (0.5, 180.0), (10.0, 140.0)] {
        let mut algorithm = ffid(&[("detect:mapping_tolerance", ParamValue::Float(tolerance))]);
        algorithm.set_ms_data(experiment(mz, 200.0, 10.0));
        let mut features = FeatureMap::new();
        algorithm
            .run(
                vec![identification(200.0, mz, SEQUENCE, 2)],
                &proteins(),
                &mut features,
                &FeatureMap::new(),
                "",
            )
            .unwrap();
        // peak width 30 s: (60 + tolerance(abs or 0.5 * 60)) * 2
        assert_eq!(algorithm.rt_window(), expected, "tolerance {tolerance}");
    }
    // a set window is kept
    let mut algorithm = ffid(&[("extract:rt_window", ParamValue::Float(77.0))]);
    algorithm.set_ms_data(experiment(mz, 200.0, 10.0));
    let mut features = FeatureMap::new();
    algorithm
        .run(
            vec![identification(200.0, mz, SEQUENCE, 2)],
            &proteins(),
            &mut features,
            &FeatureMap::new(),
            "",
        )
        .unwrap();
    assert_eq!(algorithm.rt_window(), 77.0);
}

#[test]
fn a_relative_minimum_peak_width_becomes_seconds_in_place() {
    // `detect:min_peak_width` 0.2 is relative to the peak width (30 s): 6 s
    let mz = precursor_mz();
    let mut algorithm = ffid(&[]);
    assert_eq!(algorithm.min_peak_width(), 0.2);
    algorithm.set_ms_data(experiment(mz, 200.0, 10.0));
    let mut features = FeatureMap::new();
    algorithm
        .run(
            vec![identification(200.0, mz, SEQUENCE, 2)],
            &proteins(),
            &mut features,
            &FeatureMap::new(),
            "",
        )
        .unwrap();
    assert!((algorithm.min_peak_width() - 6.0).abs() < 1e-12);
}

#[test]
fn decoys_are_set_aside_unless_they_are_quantified() {
    let mz = precursor_mz();
    let mut decoy = identification(200.0, mz, SEQUENCE, 2);
    decoy.hits[0]
        .metadata
        .insert("target_decoy".to_owned(), MetaValue::from("decoy"));
    for (quantify, expect_feature) in [("false", false), ("true", true)] {
        let mut algorithm = ffid(&[("quantify_decoys", ParamValue::String(quantify.into()))]);
        algorithm.set_ms_data(experiment(mz, 200.0, 10.0));
        let mut features = FeatureMap::new();
        algorithm
            .run(
                vec![decoy.clone()],
                &proteins(),
                &mut features,
                &FeatureMap::new(),
                "",
            )
            .unwrap();
        assert_eq!(!features.features.is_empty(), expect_feature, "{quantify}");
        if !expect_feature {
            // the decoy is returned as an unassigned identification
            assert_eq!(features.unassigned_peptide_identifications.len(), 1);
        }
    }
}

#[test]
fn the_score_cutoff_depends_on_the_direction_of_the_score() {
    // q-values (lower is better): a cutoff of 0.05 keeps 0.01 and drops 0.5;
    // the dropped identification is unassigned
    let mz = precursor_mz();
    let kept = identification(200.0, mz, SEQUENCE, 2);
    let mut dropped = identification(200.0, mz, "AAAAAK", 2);
    dropped.hits[0].score = 0.5;
    let mut algorithm = ffid(&[("min_psm_cutoff", ParamValue::String("0.05".into()))]);
    algorithm.set_ms_data(experiment(mz, 200.0, 10.0));
    let mut features = FeatureMap::new();
    algorithm
        .run(
            vec![kept, dropped],
            &proteins(),
            &mut features,
            &FeatureMap::new(),
            "",
        )
        .unwrap();
    assert_eq!(features.features.len(), 1);
    assert_eq!(
        features
            .unassigned_peptide_identifications
            .iter()
            .filter(|id| id.hits[0].sequence.to_string() == "AAAAAK")
            .count(),
        1
    );
    assert_eq!(algorithm.report().peptides_identified, 1);
}

#[test]
fn identifications_without_hits_are_dropped_silently() {
    let mz = precursor_mz();
    let empty = PeptideIdentification {
        identifier: "run".to_owned(),
        rt: Some(200.0),
        mz: Some(mz),
        ..PeptideIdentification::default()
    };
    let mut algorithm = ffid(&[]);
    algorithm.set_ms_data(experiment(mz, 200.0, 10.0));
    let mut features = FeatureMap::new();
    algorithm
        .run(
            vec![empty, identification(200.0, mz, SEQUENCE, 2)],
            &proteins(),
            &mut features,
            &FeatureMap::new(),
            "",
        )
        .unwrap();
    assert_eq!(features.features.len(), 1);
    assert_eq!(algorithm.report().peptides_identified, 1);
}

#[test]
fn seeds_overlapping_an_identification_are_not_added() {
    // `addSeeds_`: a seed within half the seed window (30 s) in RT and within
    // the m/z tolerance of an identification - or of one of its first five
    // isotopes - is not a new target.
    let mz = precursor_mz();
    let seed = |rt: f64, seed_mz: f64| {
        let mut feature = Feature::new(rt, seed_mz, 1000.0);
        feature.charge = 2;
        feature.unique_id = 5000 + rt as u64;
        feature
    };
    let mut seeds = FeatureMap::new();
    seeds.features.push(seed(210.0, mz)); // overlaps in RT and m/z
    seeds
        .features
        .push(seed(210.0, mz + C13C12_MASSDIFF_U / 2.0)); // second isotope
    seeds.features.push(seed(290.0, mz)); // 90 s away
    seeds.features.push(seed(200.0, mz + 25.0)); // 25 m/z away
    let mut algorithm = ffid(&[]);
    algorithm.set_ms_data(experiment(mz, 200.0, 10.0));
    let mut features = FeatureMap::new();
    algorithm
        .run(
            vec![identification(200.0, mz, SEQUENCE, 2)],
            &proteins(),
            &mut features,
            &seeds,
            "",
        )
        .unwrap();
    assert_eq!(algorithm.report().seeds_added, 2);
    // seeds are assays of their own, without a sequence
    let library = algorithm.library();
    let seed_assays = library
        .peptides()
        .iter()
        .filter(|p| p.sequence.is_empty())
        .count();
    assert_eq!(seed_assays, 2);
    // their pseudo identifications never reach the output
    for feature in &features.features {
        for id in &feature.peptide_identifications {
            assert!(
                id.hits
                    .iter()
                    .all(|h| !h.sequence.as_str().starts_with("XXX"))
            );
        }
    }
    for id in &features.unassigned_peptide_identifications {
        assert!(
            id.hits
                .iter()
                .all(|h| !h.sequence.as_str().starts_with("XXX"))
        );
    }
}

#[test]
fn offset_peptides_double_the_targets_and_are_marked() {
    let mz = precursor_mz();
    let mut algorithm = ffid(&[("add_mass_offset_peptides", ParamValue::Float(5.0))]);
    algorithm.set_ms_data(experiment(mz, 200.0, 10.0));
    let mut features = FeatureMap::new();
    algorithm
        .run(
            vec![identification(200.0, mz, SEQUENCE, 2)],
            &proteins(),
            &mut features,
            &FeatureMap::new(),
            "",
        )
        .unwrap();
    assert_eq!(algorithm.report().offset_peptides_added, 1);
    // the offset peptide is an assay at m/z + 5 without sequence
    let offset = algorithm
        .library()
        .peptides()
        .iter()
        .find(|p| p.sequence.is_empty())
        .expect("the offset assay");
    assert_eq!(offset.rts[0].rt().unwrap(), 200.0 - 30.0);
    assert_eq!(offset.rts[1].rt().unwrap(), 200.0 + 30.0);
    let transitions: Vec<_> = algorithm
        .library()
        .transitions()
        .iter()
        .filter(|t| t.peptide_ref() == offset.id)
        .collect();
    assert!(!transitions.is_empty());
    assert!((transitions[0].precursor_mz() - (mz + 5.0)).abs() < 1e-9);
}

#[test]
fn a_run_is_deterministic() {
    let mz = precursor_mz();
    let run = || {
        let mut algorithm = ffid(&[]);
        algorithm.set_ms_data(experiment(mz, 200.0, 10.0));
        let mut features = FeatureMap::new();
        algorithm
            .run(
                vec![identification(200.0, mz, SEQUENCE, 2)],
                &proteins(),
                &mut features,
                &FeatureMap::new(),
                "",
            )
            .unwrap();
        features
    };
    assert_eq!(run(), run());
}

#[test]
fn candidates_are_reprocessed_without_the_spectra() {
    // a candidate without an identification is removed by the filter; with one
    // it stays (`runOnCandidates`)
    let mz = precursor_mz();
    let mut algorithm = ffid(&[]);
    algorithm.set_ms_data(experiment(mz, 200.0, 10.0));
    let mut features = FeatureMap::new();
    algorithm
        .run(
            vec![identification(200.0, mz, SEQUENCE, 2)],
            &proteins(),
            &mut features,
            &FeatureMap::new(),
            "",
        )
        .unwrap();
    assert_eq!(features.features.len(), 1);
    let mut candidates = features.clone();
    let mut orphan = candidates.features[0].clone();
    orphan.peptide_identifications.clear();
    orphan.rt += 1.0;
    candidates.features.push(orphan);
    let mut again = FeatureFinderIdentificationAlgorithm::from_param(&parameters(&[])).unwrap();
    again.run_on_candidates(&mut candidates).unwrap();
    assert_eq!(candidates.features.len(), 1);
}

// ---------------------------------------------------------------- MetaboIdent

fn glucose(rt: f64) -> FeatureFinderMetaboIdentCompound {
    FeatureFinderMetaboIdentCompound::new(
        "glucose",
        "C6H12O6",
        0.0,
        vec![1],
        vec![rt],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        "",
    )
}

/// The [M+H]+ ion of glucose, the formula of the compound `glucose`.
fn glucose_mz() -> f64 {
    181.070665 // monoisotopic 180.063388 + proton 1.007276
}

fn metabo() -> FeatureFinderAlgorithmMetaboIdent {
    let mut param = FeatureFinderAlgorithmMetaboIdent::defaults().unwrap();
    param
        .set_value("extract:mz_window", ParamValue::Float(0.1), "", &[])
        .unwrap();
    param
        .set_value("detect:peak_width", ParamValue::Float(30.0), "", &[])
        .unwrap();
    param
        .set_value("model:type", ParamValue::String("none".into()), "", &[])
        .unwrap();
    let mut algorithm = FeatureFinderAlgorithmMetaboIdent::from_param(&param).unwrap();
    algorithm.set_unique_id_generator(UniqueIdGenerator::from_seed(11));
    algorithm
}

#[test]
fn metabo_ident_rt_window_is_four_peak_widths() {
    assert_eq!(
        FeatureFinderAlgorithmMetaboIdent::new()
            .unwrap()
            .rt_window(),
        240.0
    );
    assert_eq!(metabo().rt_window(), 120.0);
}

#[test]
fn metabo_ident_finds_an_adduct_free_ion_and_annotates_the_target() {
    let mut algorithm = metabo();
    algorithm.set_ms_data(experiment(
        glucose_mz() - C13C12_MASSDIFF_U * 0.0,
        200.0,
        10.0,
    ));
    // the experiment has three isotopes spaced for charge two; the target is
    // singly charged, so only its monoisotopic trace matches and the other
    // isotopes are not extracted at their charge-one spacing
    let mut features = FeatureMap::new();
    algorithm
        .run(&[glucose(200.0)], &mut features, "glucose.mzML")
        .unwrap();
    assert_eq!(features.features.len(), 1);
    let feature = &features.features[0];
    assert_eq!(feature.charge, 1);
    let text = |key: &str| {
        feature
            .metadata
            .get(key)
            .unwrap_or_else(|| panic!("meta value {key}"))
            .as_str()
            .unwrap()
            .to_owned()
    };
    assert_eq!(text("label"), "glucose");
    assert_eq!(text("sum_formula"), "C6H12O6");
    assert_eq!(text("adduct"), "M+H;1+");
    // the id format of the Release build: `bad_adduct_m180.063385_z1_rt150.0_i0`
    assert_eq!(text("PeptideRef"), "glucose_m180.063385_z1_rt200.0");
    assert_eq!(
        feature
            .metadata
            .get("expected_rt")
            .unwrap()
            .as_f64()
            .unwrap(),
        200.0
    );
    assert_eq!(
        feature
            .metadata
            .get("num_of_masstraces")
            .unwrap()
            .as_i64()
            .unwrap(),
        feature.subordinates.len() as i64
    );
    // the retention time offset of the one data point
    let points = algorithm.transformations().data_points();
    assert_eq!(points.len(), 1);
    assert_eq!(points[0].x, 200.0);
    assert!((points[0].y - feature.rt).abs() < 1e-12);
    assert_eq!(algorithm.n_shared(), 0);
    assert!(features.unassigned_peptide_identifications.is_empty());
}

#[test]
fn metabo_ident_reports_targets_without_a_feature() {
    let mut algorithm = metabo();
    algorithm.set_ms_data(experiment(glucose_mz(), 200.0, 10.0));
    let mut features = FeatureMap::new();
    let other = FeatureFinderMetaboIdentCompound::new(
        "absent",
        "C10H15N5O10P2",
        0.0,
        vec![1],
        vec![150.0],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        "",
    );
    algorithm
        .run(&[glucose(200.0), other], &mut features, "")
        .unwrap();
    assert_eq!(features.features.len(), 1);
    let missing = &features.unassigned_peptide_identifications;
    assert_eq!(missing.len(), 1);
    assert_eq!(
        missing[0].metadata.get("label").unwrap().as_str().unwrap(),
        "absent"
    );
    assert_eq!(missing[0].rt, Some(150.0));
    assert_eq!(features.protein_identifications.len(), 1);
    assert_eq!(features.protein_identifications[0].identifier, "id");
}

#[test]
fn metabo_ident_skips_targets_that_violate_the_contract() {
    let mut algorithm = metabo();
    algorithm.set_ms_data(experiment(glucose_mz(), 200.0, 10.0));
    let mut features = FeatureMap::new();
    let no_mass = FeatureFinderMetaboIdentCompound::new(
        "no_mass",
        "",
        0.0,
        vec![1],
        vec![200.0],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        "",
    );
    let no_rt = FeatureFinderMetaboIdentCompound::new(
        "no_rt",
        "C6H12O6",
        0.0,
        vec![1],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        "",
    );
    let bad_ranges = FeatureFinderMetaboIdentCompound::new(
        "bad_ranges",
        "C6H12O6",
        0.0,
        vec![1],
        vec![200.0, 210.0],
        vec![1.0, 2.0, 3.0],
        Vec::new(),
        Vec::new(),
        "",
    );
    algorithm
        .run(
            &[no_mass, no_rt, bad_ranges, glucose(200.0)],
            &mut features,
            "",
        )
        .unwrap();
    assert_eq!(algorithm.report().targets_skipped, 3);
    assert_eq!(algorithm.library().compounds().len(), 1);
    assert_eq!(features.features.len(), 1);
}

#[test]
fn metabo_ident_refuses_an_oversized_library_before_building_it() {
    let mut algorithm = metabo();
    algorithm.set_ms_data(experiment(glucose_mz(), 200.0, 10.0));
    let charges: Vec<i32> = (1..=5000).collect();
    let rts: Vec<f64> = (0..5000).map(f64::from).collect();
    assert!(charges.len() * rts.len() > MAX_LIBRARY_ENTRIES);
    let target = FeatureFinderMetaboIdentCompound::new(
        "huge",
        "C6H12O6",
        0.0,
        charges,
        rts,
        Vec::new(),
        Vec::new(),
        Vec::new(),
        "",
    );
    let mut features = FeatureMap::new();
    let error = algorithm.run(&[target], &mut features, "").unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    assert!(algorithm.library().compounds().is_empty());
    assert!(features.features.is_empty());
}

#[test]
fn metabo_ident_without_ms1_spectra_leaves_the_features_unchanged() {
    let mut algorithm = metabo();
    let mut features = FeatureMap::new();
    algorithm
        .run(&[glucose(200.0)], &mut features, "empty.mzML")
        .unwrap();
    assert!(features.features.is_empty());
    assert_eq!(
        features.primary_ms_run_path().unwrap(),
        vec!["empty.mzML".to_owned()]
    );
    assert!(
        algorithm
            .report()
            .messages
            .iter()
            .any(|line| line.contains("No MS1 scans"))
    );
}

#[test]
fn metabo_ident_default_adducts_follow_the_charge() {
    // library entries of charge 1, 2 and -1: [M+H]+, [M+2H]2+ and [M-H]-
    let mut algorithm = metabo();
    algorithm.set_ms_data(experiment(glucose_mz(), 200.0, 10.0));
    let target = FeatureFinderMetaboIdentCompound::new(
        "multi",
        "C6H12O6",
        0.0,
        vec![1, 2, -1],
        vec![200.0],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        "",
    );
    let mut features = FeatureMap::new();
    algorithm.run(&[target], &mut features, "").unwrap();
    let library = algorithm.library();
    let adducts: Vec<String> = library
        .compounds()
        .iter()
        .map(|c| {
            c.cv_terms
                .metadata
                .get("adduct")
                .unwrap()
                .as_str()
                .unwrap()
                .to_owned()
        })
        .collect();
    assert_eq!(adducts, ["M+H;1+", "M+2H;2+", "M-H;1-"]);
    let neutral = 180.063_388_56;
    let proton = 1.007_276_466_771;
    let mzs: Vec<f64> = library
        .compounds()
        .iter()
        .map(|c| {
            library
                .transitions()
                .iter()
                .find(|t| t.compound_ref() == c.id)
                .unwrap()
                .precursor_mz()
        })
        .collect();
    // [M+H]+ and [M+2H]2+ lose or gain the electron mass the source accounts
    // for; the proton-mass formula agrees to well under a milli-Dalton
    assert!((mzs[0] - (neutral + proton)).abs() < 1e-3, "{}", mzs[0]);
    assert!(
        (mzs[1] - (neutral + 2.0 * proton) / 2.0).abs() < 1e-3,
        "{}",
        mzs[1]
    );
    assert!((mzs[2] - (neutral - proton)).abs() < 1e-3, "{}", mzs[2]);
}
