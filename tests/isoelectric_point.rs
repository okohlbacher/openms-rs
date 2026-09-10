// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source: OpenMS4-core 7c029e8 IsoelectricPoint class tests and independent charge checks.

use openms::Error;
use openms::chemistry::isoelectric_point::MAX_PI_RESIDUES;
use openms::chemistry::{AASequence, IsoelectricPoint, ProteomicsPkaScale};

fn peptide(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}
fn model(scale: ProteomicsPkaScale) -> IsoelectricPoint {
    IsoelectricPoint {
        scale,
        ..Default::default()
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.15} != {expected:.15} (tolerance {tolerance})"
    );
}

#[test]
fn source_defaults_and_single_residue_pi_goldens() {
    let default = IsoelectricPoint::default();
    assert_eq!(default.scale, ProteomicsPkaScale::Lehninger);
    assert_eq!(default.tolerance, 1e-4);
    for (sequence, expected) in [("A", 6.015), ("K", 10.11), ("D", 2.995), ("U", 4.035)] {
        close(
            default.compute_pi(&peptide(sequence)).unwrap(),
            expected,
            5e-5,
        );
    }
    let bjellqvist = model(ProteomicsPkaScale::Bjellqvist);
    for (sequence, expected) in [("A", 5.57), ("P", 5.955)] {
        close(
            bjellqvist.compute_pi(&peptide(sequence)).unwrap(),
            expected,
            5e-5,
        );
    }
}

#[test]
fn source_emboss_literal_goldens_and_sillero_zero_crossing() {
    let emboss = IsoelectricPoint {
        scale: ProteomicsPkaScale::Emboss,
        tolerance: 1e-6,
    };
    // Four-decimal upstream expectations; retain that literal precision.
    for (sequence, expected) in [
        ("ACDEFGHIK", 5.4518),
        ("PEPTIDER", 3.9276),
        ("SAMPLER", 6.4101),
        ("YGGFMR", 9.3488),
        ("DEEEAANK", 3.7782),
    ] {
        close(
            emboss.compute_pi(&peptide(sequence)).unwrap(),
            expected,
            1e-4,
        );
    }
    let sequence = peptide("DEKRPEPTIDE");
    let sillero = model(ProteomicsPkaScale::Sillero);
    let pi = sillero.compute_pi(&sequence).unwrap();
    assert!(pi > 0.0 && pi < 14.0);
    assert!(sillero.compute_charge(&sequence, pi - 1e-4).unwrap() > 0.0);
    assert!(sillero.compute_charge(&sequence, pi + 1e-4).unwrap() < 0.0);
    close(sillero.compute_charge(&sequence, pi).unwrap(), 0.0, 1e-3);
    assert_ne!(
        pi,
        IsoelectricPoint::default().compute_pi(&sequence).unwrap()
    );
}

#[test]
fn terminal_charge_matches_independent_acid_base_equations() {
    let sequence = peptide("A");
    for (scale, nterm, cterm) in [
        (ProteomicsPkaScale::Lehninger, 9.69, 2.34),
        (ProteomicsPkaScale::Emboss, 8.6, 3.6),
        (ProteomicsPkaScale::Sillero, 8.2, 3.2),
        (ProteomicsPkaScale::Bjellqvist, 7.59, 3.55),
    ] {
        for ph in [-3.0, 0.0, 2.34, 7.0, 14.0, 17.0] {
            let expected =
                1.0 / (1.0 + 10.0_f64.powf(ph - nterm)) - 1.0 / (1.0 + 10.0_f64.powf(cterm - ph));
            assert_eq!(
                model(scale).compute_charge(&sequence, ph).unwrap(),
                expected
            );
        }
    }
    // A has its mathematical zero at the mean terminal pKa. The source stops
    // on pH width, not when an interior charge residual first becomes small.
    let width = 14.0_f64 / 2.0_f64.powi(18);
    let expected = (6.015 / width).floor() * width + width / 2.0;
    assert_eq!(
        IsoelectricPoint::default().compute_pi(&sequence).unwrap(),
        expected
    );
}

#[test]
fn named_and_anonymous_terminal_annotations_block_only_their_terminal_group() {
    let model = IsoelectricPoint::default();
    let plain = peptide("Y");
    let acetyl = peptide(".(Acetyl)Y");
    let amidated = peptide("Y.(Amidated)");
    for ph in [1.0, 7.0, 14.0] {
        let unmodified_charge = model.compute_charge(&plain, ph).unwrap();
        let n = 1.0 / (1.0 + 10.0_f64.powf(ph - 9.69));
        let c = -1.0 / (1.0 + 10.0_f64.powf(2.34 - ph));
        close(
            model.compute_charge(&acetyl, ph).unwrap(),
            unmodified_charge - n,
            2e-15,
        );
        close(
            model.compute_charge(&amidated, ph).unwrap(),
            unmodified_charge - c,
            2e-15,
        );
    }
    let mut tagged = peptide("A");
    tagged.set_n_terminal_mass_tag("+123.456789").unwrap();
    assert!(
        tagged
            .n_terminal_modification()
            .unwrap()
            .mass_tag()
            .is_some()
    );
    assert!(tagged.formula().is_err());
    close(
        model.compute_charge(&tagged, 7.0).unwrap(),
        -1.0 / (1.0 + 10.0_f64.powf(2.34 - 7.0)),
        0.0,
    );
    tagged.set_c_terminal_mass_tag("+234.567891").unwrap();
    for ph in [-f64::MAX, 0.0, 7.0, 14.0, f64::MAX] {
        assert_eq!(model.compute_charge(&tagged, ph).unwrap(), 0.0);
    }
    assert_eq!(model.compute_pi(&tagged).unwrap(), 0.0);
}

#[test]
fn parent_residue_sidechain_charge_ignores_named_and_mass_only_modifications() {
    let sequences = [
        (peptide("K"), peptide("K(Acetyl)")),
        (peptide("S"), peptide("S(Phospho)")),
        (peptide("C"), peptide("C(Carbamidomethyl)")),
        (peptide("X"), peptide("X[123.456789]")),
    ];
    for scale in [
        ProteomicsPkaScale::Lehninger,
        ProteomicsPkaScale::Emboss,
        ProteomicsPkaScale::Sillero,
        ProteomicsPkaScale::Bjellqvist,
    ] {
        let model = model(scale);
        for (plain, modified) in &sequences {
            for ph in [0.0, 5.73, 7.0, 14.0] {
                assert_eq!(
                    model.compute_charge(plain, ph).unwrap(),
                    model.compute_charge(modified, ph).unwrap()
                );
            }
            assert_eq!(
                model.compute_pi(plain).unwrap(),
                model.compute_pi(modified).unwrap()
            );
        }
    }
}

#[test]
fn ambiguous_residues_are_neutral_selenocysteine_is_acidic_and_pyrrolysine_errors() {
    for scale in [
        ProteomicsPkaScale::Lehninger,
        ProteomicsPkaScale::Emboss,
        ProteomicsPkaScale::Sillero,
        ProteomicsPkaScale::Bjellqvist,
    ] {
        let model = model(scale);
        // N has no ionizable sidechain and uses each scale's fallback terminal pKa.
        let neutral = peptide("N");
        for residue in ["B", "Z", "X", "J", "BZXJ"] {
            let ambiguous = peptide(residue);
            for ph in [0.0, 5.73, 7.0, 14.0] {
                assert_eq!(
                    model.compute_charge(&ambiguous, ph).unwrap(),
                    model.compute_charge(&neutral, ph).unwrap()
                );
            }
        }
        let mut sec = peptide("U");
        sec.set_n_terminal_mass_tag("+123.456789").unwrap();
        sec.set_c_terminal_mass_tag("+234.567891").unwrap();
        assert_eq!(model.compute_charge(&sec, 5.73).unwrap(), -0.5);
        assert!(matches!(
            model.compute_charge(&peptide("O"), 7.0),
            Err(Error::Unsupported(_))
        ));
        assert!(model.compute_pi(&peptide("O")).is_err());
    }
}

#[test]
fn charge_is_monotone_and_finite_extreme_ph_saturates_without_error() {
    let sequence = peptide("ACDEFGHIKLMNPQRSTVWYU");
    for scale in [
        ProteomicsPkaScale::Lehninger,
        ProteomicsPkaScale::Emboss,
        ProteomicsPkaScale::Sillero,
        ProteomicsPkaScale::Bjellqvist,
    ] {
        let model = model(scale);
        let mut last = model.compute_charge(&sequence, -f64::MAX).unwrap();
        assert_eq!(last, 4.0);
        for ph in [
            -1000.0,
            -1.0,
            0.0,
            1.0,
            4.0,
            7.0,
            11.0,
            14.0,
            20.0,
            1000.0,
            f64::MAX,
        ] {
            let charge = model.compute_charge(&sequence, ph).unwrap();
            assert!(charge <= last);
            last = charge;
        }
        assert_eq!(last, -6.0);
    }
}

#[test]
fn source_endpoint_priority_and_no_crossing_behavior() {
    let model = IsoelectricPoint::default();
    assert_eq!(model.compute_pi(&peptide(&"R".repeat(100))).unwrap(), 14.0);
    assert_eq!(model.compute_pi(&peptide(".(Acetyl)DDD")).unwrap(), 0.0);
    assert_eq!(model.compute_pi(&peptide("RRR.(Amidated)")).unwrap(), 14.0);
    // Both endpoint residuals fit this tolerance. The low endpoint wins before
    // the interval-width rule could otherwise return an interior midpoint.
    let large_tolerance = IsoelectricPoint {
        tolerance: 14.0,
        ..model
    };
    assert_eq!(large_tolerance.compute_pi(&peptide("A")).unwrap(), 0.0);
}

#[test]
fn invalid_arguments_and_unattainable_tolerance_are_checked() {
    let model = IsoelectricPoint::default();
    let sequence = peptide("A");
    assert!(model.compute_charge(&AASequence::default(), 7.0).is_err());
    assert!(model.compute_pi(&AASequence::default()).is_err());
    for ph in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(model.compute_charge(&sequence, ph).is_err());
    }
    for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let configured = IsoelectricPoint { tolerance, ..model };
        assert!(configured.compute_pi(&sequence).is_err());
        assert_eq!(
            configured.compute_charge(&sequence, 7.0).unwrap(),
            model.compute_charge(&sequence, 7.0).unwrap()
        );
    }
    let tiny_tolerance = IsoelectricPoint {
        tolerance: f64::from_bits(1),
        ..model
    };
    assert!(
        tiny_tolerance
            .compute_pi(&sequence)
            .unwrap_err()
            .to_string()
            .contains("cannot progress")
    );
}

#[test]
fn public_residue_limit_and_cumulative_work_limit_are_enforced() {
    let text = "X".repeat(MAX_PI_RESIDUES + 1);
    let model = IsoelectricPoint::default();
    {
        let too_long = peptide(&text);
        assert!(
            model
                .compute_charge(&too_long, 7.0)
                .unwrap_err()
                .to_string()
                .contains("residue limit")
        );
        assert!(
            model
                .compute_pi(&too_long)
                .unwrap_err()
                .to_string()
                .contains("residue limit")
        );
    }
    let at_limit = peptide(&text[..MAX_PI_RESIDUES]);
    assert_eq!(model.compute_charge(&at_limit, -f64::MAX).unwrap(), 1.0);
    // Neutral parent codes make each pass cheap, but all visits still count.
    // The budget fails before bisection can reach floating-point stagnation.
    let tiny_tolerance = IsoelectricPoint {
        tolerance: f64::from_bits(1),
        ..model
    };
    assert!(
        tiny_tolerance
            .compute_pi(&at_limit)
            .unwrap_err()
            .to_string()
            .contains("pI work limit")
    );
}
