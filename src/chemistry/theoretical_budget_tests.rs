// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::*;

#[test]
fn candidate_generation_shares_residue_and_loss_allowances() {
    let peptide = AASequence::parse("AG").unwrap();
    let generator = TheoreticalSpectrumGenerator::default();
    // For two residues/two series/one charge: n²*s + n*(s+charges+1) = 16.
    let mut work = TheoreticalGenerationWork {
        residue_remaining: 32,
        ..TheoreticalGenerationWork::default()
    };
    for _ in 0..2 {
        assert_eq!(
            generator
                .generate_with_work(&peptide, 1, 1, None, &mut work)
                .unwrap(),
            generator.generate(&peptide, 1, 1, None).unwrap()
        );
    }
    assert_eq!(work.residue_remaining, 0);
    assert!(
        generator
            .generate_with_work(&peptide, 1, 1, None, &mut work)
            .unwrap_err()
            .to_string()
            .contains("cumulative theoretical residue")
    );

    let generator = TheoreticalSpectrumGenerator {
        add_losses: true,
        ..generator
    };
    // Only the G suffix visits declarations: the source has no b1 by default.
    // The precharge reserves the default three declarations plus one visit.
    let mut work = TheoreticalGenerationWork {
        loss_remaining: 8,
        ..TheoreticalGenerationWork::default()
    };
    for _ in 0..2 {
        generator
            .generate_with_work(&peptide, 1, 1, None, &mut work)
            .unwrap();
    }
    assert_eq!(work.loss_remaining, 0);
    assert!(
        generator
            .generate_with_work(&peptide, 1, 1, None, &mut work)
            .unwrap_err()
            .to_string()
            .contains("cumulative theoretical loss")
    );
}

#[test]
fn candidate_generation_does_not_reset_fine_or_coarse_search_allowances() {
    let peptide = AASequence::parse("AG").unwrap();
    for (model, expected_error) in [
        (
            TheoreticalIsotopeModel::Fine {
                unexplained_probability: 0.01,
            },
            "fine isotope work",
        ),
        (
            TheoreticalIsotopeModel::Coarse { max_peaks: 3 },
            "isotope convolution exceeds work",
        ),
    ] {
        let generator = TheoreticalSpectrumGenerator {
            isotope_model: model,
            ..TheoreticalSpectrumGenerator::default()
        };
        let expected = generator.generate(&peptide, 1, 1, None).unwrap();
        let mut work = TheoreticalGenerationWork::default();
        // Precharge the allowance to keep this exhaustion regression small.
        work.fine
            .consume(super::super::fine_isotopes::MAX_FINE_WORK - 3_000)
            .unwrap();
        work.coarse
            .consume(super::super::isotopes::MAX_CONVOLUTION_PRODUCTS - 2_000)
            .unwrap();
        let mut successes = 0;
        let mut exhausted = false;
        for _ in 0..128 {
            match generator.generate_with_work(&peptide, 1, 1, None, &mut work) {
                Ok(spectrum) => {
                    assert_eq!(spectrum, expected);
                    successes += 1;
                }
                Err(error) => {
                    assert!(error.to_string().contains(expected_error), "{error}");
                    exhausted = true;
                    break;
                }
            }
        }
        assert!(successes > 0 && exhausted, "{model:?}, {successes} calls");
    }
}
