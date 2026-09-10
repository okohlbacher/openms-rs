// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Pinned source: TheoreticalSpectrumGenerator.cpp:1165–1261 and
// PrecursorPurity_test.cpp getPrefixAndSuffixIonsMZ consumer at revision7c029e8.

use openms::chemistry::theoretical::{
    MAX_THEORETICAL_PEAKS, MAX_THEORETICAL_RESIDUES, TheoreticalIonSeries as Ion,
    TheoreticalIsotopeModel as Isotopes, TheoreticalSpectrumGenerator as Generator,
};
use openms::chemistry::{AASequence, EmpiricalFormula, PROTON_MASS_U};

fn peptide(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}
fn formula(text: &str) -> f64 {
    EmpiricalFormula::parse(text).unwrap().mono_mass()
}
fn bits(values: &[f32]) -> Vec<u32> {
    values.iter().map(|value| value.to_bits()).collect()
}

#[test]
fn source_precursor_purity_consumer_goldens_and_full_length_counts() {
    let mut output = vec![2000., 0.25];
    let p = peptide("PEPTIDER");
    Generator::default()
        .append_mass_spectrum(&mut output, &p, 2)
        .unwrap();
    assert_eq!(output.len(), 2 + p.len() * 2 * 2);
    assert!(output.windows(2).all(|pair| pair[0] <= pair[1]));
    // Independently stated values in the pinned PrecursorPurity class test.
    for expected in [425.20308, 419.18849, 730.37299, 430.21143] {
        assert!(
            output
                .iter()
                .any(|&actual| (f64::from(actual) - expected).abs() < 0.0001),
            "{expected}"
        );
    }
    // Full-length y includes water; full-length b does not.
    for charge in [1, 2] {
        let y = p.mz(charge).unwrap() as f32;
        let b = ((p.mono_mass().unwrap() - formula("H2O") + f64::from(charge) * PROTON_MASS_U)
            / f64::from(charge)) as f32;
        assert!(output.contains(&y));
        assert!(output.contains(&b));
    }
}

#[test]
fn all_six_series_include_singleton_and_ignore_richer_spectrum_settings() {
    let mut g = Generator {
        ion_series: vec![
            Ion::ZPlusOne,
            Ion::Y,
            Ion::A,
            Ion::B,
            Ion::C,
            Ion::X,
            Ion::Z,
            Ion::B,
            Ion::ZPlusTwo,
        ],
        add_first_prefix_ion: false,
        add_losses: true,
        add_terminal_losses: true,
        add_precursor_peaks: true,
        add_all_precursor_charges: true,
        add_metainfo: true,
        sort_by_position: false,
        isotope_model: Isotopes::Coarse { max_peaks: 0 },
        relative_loss_intensity: f32::NAN,
        ..Default::default()
    };
    g.intensities.b = f32::NEG_INFINITY;
    let mut actual = vec![900., 0.];
    g.append_mass_spectrum(&mut actual, &peptide("A"), 1)
        .unwrap();
    // Independent one-residue elemental forms for ordinary a/b/c/x/y/z ions.
    let mut expected: Vec<f32> = ["C2H5N", "C3H5NO", "C3H8N2O", "C4H5NO3", "C3H7NO2", "C3H4O2"]
        .into_iter()
        .map(|text| (formula(text) + PROTON_MASS_U) as f32)
        .collect();
    expected.extend([900., 0.]);
    expected.sort_by(f32::total_cmp);
    assert_eq!(actual, expected);
}

#[test]
fn corresponding_terminal_deltas_apply_even_for_full_length_ions() {
    let p = peptide(".(Acetyl)A.(Amidated)");
    let before = p.clone();
    let mut output = Vec::new();
    Generator::default()
        .append_mass_spectrum(&mut output, &p, 1)
        .unwrap();
    let n_delta = p
        .n_terminal_modification()
        .unwrap()
        .diff_mono_mass()
        .unwrap();
    let c_delta = p
        .c_terminal_modification()
        .unwrap()
        .diff_mono_mass()
        .unwrap();
    let mut expected = vec![
        (PROTON_MASS_U + n_delta + formula("C3H5NO")) as f32,
        (PROTON_MASS_U + c_delta + formula("H2O") + formula("C3H5NO")) as f32,
    ];
    expected.sort_by(f32::total_cmp);
    assert_eq!(output, expected);
    assert_eq!(p, before);
    // Neither full-length ladder includes the opposite terminal modification.
    assert!(!output.contains(&(p.mz(1).unwrap() as f32)));
}

#[test]
fn mass_only_chemistry_and_independent_suffix_accumulation() {
    let p = peptide("X[1000000000000000000000000000000]A");
    assert!(p.formula().is_err());
    let mut output = Vec::new();
    Generator::default()
        .append_mass_spectrum(&mut output, &p, 1)
        .unwrap();
    assert_eq!(output.len(), 4);
    assert!(output.contains(&(peptide("A").mz(1).unwrap() as f32)));
    assert!(output.iter().all(|mass| mass.is_finite()));
}

#[test]
fn zero_empty_and_disabled_calls_sort_without_requesting_unused_chemistry() {
    let cases = [
        (Generator::default(), peptide("BZX"), 0),
        (Generator::default(), AASequence::default(), 3),
        (
            Generator {
                ion_series: vec![Ion::ZPlusOne, Ion::ZPlusTwo],
                ..Default::default()
            },
            peptide("BZX"),
            3,
        ),
    ];
    for (g, p, charge) in cases {
        let mut output = vec![9., -0., 2., 0.];
        g.append_mass_spectrum(&mut output, &p, charge).unwrap();
        assert_eq!(bits(&output), bits(&[-0., 0., 2., 9.]));
    }
}

#[test]
fn errors_preserve_existing_output_and_check_f32_range() {
    let g = Generator::default();
    for p in [peptide("BX"), peptide(&format!("X[1{}]", "0".repeat(100)))] {
        let mut output = vec![9., 2.];
        let before = output.clone();
        assert!(g.append_mass_spectrum(&mut output, &p, 1).is_err());
        assert_eq!(output, before);
    }
    // Charge two fits in f32, so the later charge-one failure happens after
    // valid additions have already been generated internally.
    let mut output = vec![9., 2.];
    let p = peptide(&format!("X[5{}]", "0".repeat(38)));
    assert!(g.append_mass_spectrum(&mut output, &p, 2).is_err());
    assert_eq!(output, [9., 2.]);
    let mut output = vec![9., 2.];
    let a_only = Generator {
        ion_series: vec![Ion::A],
        ..Default::default()
    };
    assert!(
        a_only
            .append_mass_spectrum(&mut output, &peptide("X[1]"), 1)
            .is_err()
    );
    assert_eq!(output, [9., 2.]);
    for invalid in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY, -1.] {
        let mut output = vec![9., invalid, 2.];
        let before = bits(&output);
        assert!(
            g.append_mass_spectrum(&mut output, &peptide("A"), 0)
                .is_err()
        );
        assert_eq!(bits(&output), before);
    }
}

#[test]
fn residue_and_combined_output_limits_preflight_generation() {
    let g = Generator {
        ion_series: vec![Ion::B],
        ..Default::default()
    };
    let at_limit = peptide(&"A".repeat(MAX_THEORETICAL_RESIDUES));
    let mut output = Vec::new();
    g.append_mass_spectrum(&mut output, &at_limit, 1).unwrap();
    assert_eq!(output.len(), MAX_THEORETICAL_RESIDUES);
    let mut output = vec![9., 2.];
    assert!(
        g.append_mass_spectrum(
            &mut output,
            &peptide(&"A".repeat(MAX_THEORETICAL_RESIDUES + 1)),
            1
        )
        .is_err()
    );
    assert_eq!(output, [9., 2.]);
    let mut output = vec![1.; MAX_THEORETICAL_PEAKS - 1];
    g.append_mass_spectrum(&mut output, &peptide("A"), 1)
        .unwrap();
    assert_eq!(output.len(), MAX_THEORETICAL_PEAKS);
    let before = output.clone();
    assert!(
        g.append_mass_spectrum(&mut output, &peptide("A"), 1)
            .is_err()
    );
    assert_eq!(output, before);
}
