// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `MRMIonSeries_test.cpp` at core pin `bc9cc12`, plus native checks.
//! Literals are transcribed from that class test (evidence tier 3, source
//! review) unless a test says otherwise.

use openms::analysis::openswath::mrm_ion_series::{
    FRAGMENT_NEUTRAL_LOSS, IonSeries, IonSeriesOptions, MRMIonSeries, PRODUCT_ION_MZ_DELTA,
};
use openms::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition;
use openms::analysis::openswath::targeted_experiment_helper::Peptide;
use openms::chemistry::{AASequence, PeptideFragmentType};
use openms::metadata::MetaValue;

/// ClassTest `isRealSimilar` (`ClassTest.cpp:364-490`): absolute difference
/// at most 1e-5, or ratio at most 1 + 1e-5, with the source's zero cases.
fn real_similar(actual: f64, expected: f64) -> bool {
    if actual.is_nan() || expected.is_nan() {
        return false;
    }
    let small = (actual - expected).abs() <= 1e-5;
    if actual == 0.0 {
        return expected == 0.0 || small;
    }
    if expected == 0.0 {
        return small;
    }
    let mut ratio = actual / expected;
    if ratio < 0.0 {
        return small;
    }
    if ratio < 1.0 {
        ratio = 1.0 / ratio;
    }
    ratio <= 1.0 + 1e-5 || small
}

#[track_caller]
fn similar(actual: f64, expected: f64) {
    assert!(
        real_similar(actual, expected),
        "{actual} is not similar to {expected}"
    );
}

#[track_caller]
fn ion(series: &IonSeries, key: &str, expected: f64) {
    let actual = *series
        .get(key)
        .unwrap_or_else(|| panic!("{key} missing from {series:?}"));
    similar(actual, expected);
}

fn options(specific: bool, unspecific: bool) -> IonSeriesOptions {
    IonSeriesOptions::new(["b", "y"], [3, 2, 1]).with_losses(specific, unspecific)
}

fn series(sequence: &str, specific: bool, unspecific: bool) -> IonSeries {
    MRMIonSeries::get_ion_series(
        &AASequence::parse(sequence).unwrap(),
        3,
        &options(specific, unspecific),
    )
    .unwrap()
}

fn loss_value(tr: &ReactionMonitoringTransition) -> f64 {
    let term = &tr.product().interpretation_list()[0]
        .cv_terms
        .get(FRAGMENT_NEUTRAL_LOSS)
        .unwrap()[0];
    term.value.as_f64().unwrap()
}

// MRMIonSeries_test.cpp:39-146
#[test]
fn get_ion_series_standard_peptide() {
    let s = series("PEPTIDEK", false, false);
    assert_eq!(s.len(), 42);
    ion(&s, "b2^1", 227.10263491);
    ion(&s, "b2^2", 114.05495569);
    ion(&s, "b2^3", 76.37239595);
    ion(&s, "b5^1", 538.28714271);
    ion(&s, "b5^2", 269.64720959);
    ion(&s, "b5^3", 180.10056521);
    ion(&s, "b7^1", 782.35668109);
    ion(&s, "b7^2", 391.68197878);
    ion(&s, "b7^3", 261.45707801);
    ion(&s, "y1^1", 147.11280491);
    ion(&s, "y1^2", 74.06004069);
    ion(&s, "y1^3", 49.70911928);
    ion(&s, "y5^1", 605.31408687);
    ion(&s, "y5^2", 303.16068167);
    ion(&s, "y5^3", 202.44287993);
    ion(&s, "y7^1", 831.40944532);
    ion(&s, "y7^2", 416.20836089);
    ion(&s, "y7^3", 277.80799942);
}

#[test]
fn get_ion_series_with_neutral_losses() {
    let s = series("PEPTIDEK", true, true);
    assert_eq!(s.len(), 99);
    ion(&s, "b5^2", 269.64720959);
    ion(&s, "b5-H2O1^2", 260.64192709);
    // The source reads these through std::map::operator[], which inserts 0.
    for absent in ["b5-C1H4O1S1^2", "b5-H1O3P1^2", "b5-H3O4P1^2", "b5-C1O2^2"] {
        assert!(!s.contains_key(absent), "{absent}");
    }

    let s = series("ES(Phospho)", true, true);
    assert_eq!(s.len(), 12);
    ion(&s, "y1^1", 186.0162);
    ion(&s, "y1-H3O4P1^2", 44.5233);
    ion(&s, "y1^2", 93.5117);
    ion(&s, "y1^3", 62.6769);
    ion(&s, "b1-H2O1^1", 112.0393);
    ion(&s, "b1-H2O1^2", 56.5233);
    ion(&s, "y1-H3O4P1^1", 88.0393);
    ion(&s, "y1-H3O4P1^3", 30.018);
    ion(&s, "b1^2", 65.5285);
    ion(&s, "b1-H2O1^3", 38.018);
    ion(&s, "b1^1", 130.0498);
    ion(&s, "b1^3", 44.0214);

    let s = series("ES", true, true);
    ion(&s, "y1-H2O1^1", 88.0393);
    ion(&s, "y1^1", 106.0499);
    ion(&s, "y1-H2O1^2", 44.5233);
    ion(&s, "y1^2", 53.5286);
    ion(&s, "y1-H2O1^3", 30.0180);
    ion(&s, "y1^3", 36.0215);
    ion(&s, "b1-H2O1^1", 112.0393);
    ion(&s, "b1-H2O1^2", 56.5233);
    ion(&s, "b1^2", 65.5286);
    ion(&s, "b1-H2O1^3", 38.0180);
    ion(&s, "b1^1", 130.0499);
    ion(&s, "b1^3", 44.0215);

    let s = series("ES(Phospho)", true, false);
    assert_eq!(s.len(), 9);
    ion(&s, "y1^1", 186.0162);
    ion(&s, "y1-H3O4P1^2", 44.5233);
    ion(&s, "y1^2", 93.5117);
    ion(&s, "y1^3", 62.6769);
    ion(&s, "y1-H3O4P1^1", 88.0393);
    ion(&s, "y1-H3O4P1^3", 30.018);
    ion(&s, "b1^2", 65.5285);
    ion(&s, "b1^1", 130.0498);
    ion(&s, "b1^3", 44.0214);

    let s = series("ES(Phospho)", false, true);
    assert_eq!(s.len(), 9);
    ion(&s, "y1^1", 186.0162);
    ion(&s, "y1^2", 93.5117);
    ion(&s, "y1^3", 62.6769);
    ion(&s, "b1-H2O1^1", 112.0393);
    ion(&s, "b1-H2O1^2", 56.5233);
    ion(&s, "b1^2", 65.5285);
    ion(&s, "b1-H2O1^3", 38.018);
    ion(&s, "b1^1", 130.0498);
    ion(&s, "b1^3", 44.0214);
}

// MRMIonSeries_test.cpp:150-178
#[test]
fn annotate_ion() {
    let s = series("PEPTIDEK", false, false);
    let (name, mz) = MRMIonSeries::annotate_ion(&s, 202.44287993, 0.05);
    assert_eq!(name, "y5^3");
    similar(mz, 202.44287993);
    let (name, mz) = MRMIonSeries::annotate_ion(&s, 202.44287993, 0.0);
    assert_eq!(name, "unannotated");
    similar(mz, -1.0);
    let (name, mz) = MRMIonSeries::annotate_ion(&s, 202.4, 0.05);
    assert_eq!(name, "y5^3");
    similar(mz, 202.44287993);
}

// MRMIonSeries_test.cpp:182-201
#[test]
fn get_ion() {
    let s = series("PEPTIDEK", false, false);
    let (name, mz) = MRMIonSeries::get_ion(&s, "y5^3");
    assert_eq!(name, "y5^3");
    similar(mz, 202.44287993);
    assert_eq!(
        MRMIonSeries::get_ion(&s, "y9^3"),
        ("unannotated".to_owned(), -1.0)
    );
}

// MRMIonSeries_test.cpp:205-231
#[test]
fn annotate_transition_cv() {
    let mut tr = ReactionMonitoringTransition::new();
    let mut tr2 = ReactionMonitoringTransition::new();
    let mut tr3 = ReactionMonitoringTransition::new();
    MRMIonSeries::annotate_transition_cv(&mut tr, "y5^3").unwrap();
    MRMIonSeries::annotate_transition_cv(&mut tr2, "y5-H2O1^3").unwrap();
    MRMIonSeries::annotate_transition_cv(&mut tr3, "y5-18^3").unwrap();

    for t in [&tr, &tr2, &tr3] {
        assert_eq!(t.product().charge_state(), Some(3));
        let interpretation = &t.product().interpretation_list()[0];
        assert_eq!(interpretation.iontype, PeptideFragmentType::YIon);
        assert_eq!(interpretation.ordinal, 5);
    }
    assert!(
        !tr.product().interpretation_list()[0]
            .cv_terms
            .contains(FRAGMENT_NEUTRAL_LOSS)
    );
    similar(loss_value(&tr2), -18.0105650638);
    similar(loss_value(&tr3), -18.0);
    assert_eq!(
        tr.cv_terms.metadata.get("annotation"),
        Some(&MetaValue::from("y5^3"))
    );
}

// MRMIonSeries_test.cpp:235-277
#[test]
fn annotate_transition() {
    let mut peptide = Peptide::default();
    peptide.sequence = "PEPTIDEK".into();
    peptide.set_charge_state(3);

    let mut tr = ReactionMonitoringTransition::new();
    tr.set_product_mz(202.44);
    MRMIonSeries::annotate_transition(&mut tr, &peptide, 0.05, 0.05, true, &options(false, false))
        .unwrap();
    similar(tr.product_mz(), 202.442879934638);
    assert_eq!(tr.product().charge_state(), Some(3));
    let interpretation = &tr.product().interpretation_list()[0];
    assert_eq!(interpretation.iontype, PeptideFragmentType::YIon);
    assert_eq!(interpretation.ordinal, 5);
    assert!(!interpretation.cv_terms.contains(FRAGMENT_NEUTRAL_LOSS));

    let mut tr2 = ReactionMonitoringTransition::new();
    tr2.set_product_mz(196.44287993);
    MRMIonSeries::annotate_transition(&mut tr2, &peptide, 0.05, 0.05, true, &options(true, true))
        .unwrap();
    assert_eq!(tr2.product().charge_state(), Some(3));
    let interpretation = &tr2.product().interpretation_list()[0];
    assert_eq!(interpretation.iontype, PeptideFragmentType::YIon);
    assert_eq!(interpretation.ordinal, 5);
    assert_eq!(loss_value(&tr2) as i32, -18);

    let mut tr3 = ReactionMonitoringTransition::new();
    tr3.set_product_mz(202.44);
    MRMIonSeries::annotate_transition(
        &mut tr3,
        &peptide,
        0.05,
        0.05,
        false,
        &options(false, false),
    )
    .unwrap();
    similar(tr3.product_mz(), 202.44);
    assert!(!tr3.product().has_charge());
    let interpretation = &tr3.product().interpretation_list()[0];
    assert_eq!(interpretation.iontype, PeptideFragmentType::NonIdentified);
    assert!(!interpretation.cv_terms.contains(FRAGMENT_NEUTRAL_LOSS));
}

// Native: an annotated transition gains the m/z delta term, and a transition
// whose existing interpretation matches theory is kept without reannotation.
// The expected m/z is the class test's y5^3 value (tier 3).
#[test]
fn annotate_transition_without_reannotation_keeps_a_matching_interpretation() {
    let mut peptide = Peptide::default();
    peptide.sequence = "PEPTIDEK".into();
    peptide.set_charge_state(3);
    let mut tr = ReactionMonitoringTransition::new();
    MRMIonSeries::annotate_transition_cv(&mut tr, "y5^3").unwrap();
    tr.set_product_mz(202.44);
    tr.set_precursor_mz(AASequence::parse("PEPTIDEK").unwrap().mz(3).unwrap());
    MRMIonSeries::annotate_transition(&mut tr, &peptide, 0.05, 0.05, false, &options(false, false))
        .unwrap();
    similar(tr.product_mz(), 202.4429);
    let interpretation = &tr.product().interpretation_list()[0];
    assert_eq!(interpretation.iontype, PeptideFragmentType::YIon);
    let delta = interpretation.cv_terms.get(PRODUCT_ION_MZ_DELTA).unwrap()[0]
        .value
        .as_f64()
        .unwrap();
    similar(delta, 0.0029);
    assert_eq!(
        tr.cv_terms.metadata.get("annotation"),
        Some(&MetaValue::from("y5^3"))
    );
}

// Native error paths: invalid fragment type (source IllegalArgument), zero
// charge, entry ceiling, unparsable annotations. Each leaves input unchanged.
#[test]
fn error_paths() {
    let peptide = AASequence::parse("PEPTIDEK").unwrap();
    let bad_type = IonSeriesOptions::new(["q"], [1]);
    assert!(MRMIonSeries::get_ion_series(&peptide, 3, &bad_type).is_err());
    // The source raises only once the type is reached with an allowed charge.
    let skipped = IonSeriesOptions::new(["q"], [4]);
    assert!(
        MRMIonSeries::get_ion_series(&peptide, 3, &skipped)
            .unwrap()
            .is_empty()
    );
    let zero = IonSeriesOptions::new(["y"], [0]);
    assert!(MRMIonSeries::get_ion_series(&peptide, 3, &zero).is_err());
    let mut capped = IonSeriesOptions::new(["b", "y"], [1]);
    capped.max_ion_series_entries = 13;
    assert!(MRMIonSeries::get_ion_series(&peptide, 3, &capped).is_err());
    capped.max_ion_series_entries = 14;
    assert_eq!(
        MRMIonSeries::get_ion_series(&peptide, 3, &capped)
            .unwrap()
            .len(),
        14
    );
    assert!(
        MRMIonSeries::get_ion_series(&AASequence::default(), 3, &capped)
            .unwrap()
            .is_empty()
    );

    let mut tr = ReactionMonitoringTransition::new();
    let before = tr.clone();
    for bad in ["y", "", "y5^x", "y300^1", "y5-Q^2"] {
        assert!(
            MRMIonSeries::annotate_transition_cv(&mut tr, bad).is_err(),
            "{bad}"
        );
        assert_eq!(tr, before);
    }
    MRMIonSeries::annotate_transition_cv(&mut tr, "MS2_Precursor_i0").unwrap();
    assert_eq!(tr.product().interpretation_list()[0].ordinal, 0);
    MRMIonSeries::annotate_transition_cv(&mut tr, "q3^2").unwrap();
    assert_eq!(
        tr.product().interpretation_list()[0].iontype,
        PeptideFragmentType::NonIdentified
    );
}
