// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Golden values: OpenMS4-core 7c029e8 TheoreticalSpectrumGenerator_test.cpp.

use openms::chemistry::theoretical::{
    MAX_THEORETICAL_PEAKS, MAX_THEORETICAL_RESIDUES, TheoreticalIonSeries as Ion,
    TheoreticalIsotopeModel as Isotopes, TheoreticalSpectrumGenerator as Generator,
};
use openms::chemistry::{
    AASequence, C13C12_MASSDIFF_U, CoarseIsotopePatternGenerator, CoarseMassMode, EmpiricalFormula,
    PROTON_MASS_U,
};
use openms::kernel::{DataArray, MSSpectrum, Peak1D, Precursor, SpectrumType};
use std::collections::BTreeMap;

fn peptide(text: &str) -> AASequence {
    text.parse().unwrap()
}
fn formula(text: &str) -> EmpiricalFormula {
    text.parse().unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.12} != {expected:.12} (tol {tolerance})"
    );
}
fn names(spectrum: &MSSpectrum) -> &[String] {
    &spectrum
        .string_data_arrays
        .iter()
        .find(|a| a.name == "IonNames")
        .unwrap()
        .data
}
fn charges(spectrum: &MSSpectrum) -> &[i32] {
    &spectrum
        .integer_data_arrays
        .iter()
        .find(|a| a.name == "Charges")
        .unwrap()
        .data
}
fn named(spectrum: &MSSpectrum) -> BTreeMap<String, Peak1D> {
    names(spectrum)
        .iter()
        .cloned()
        .zip(spectrum.peaks.iter().copied())
        .collect()
}
fn golden_masses(spectrum: &MSSpectrum, mut expected: Vec<f64>) {
    expected.sort_by(f64::total_cmp);
    assert_eq!(spectrum.len(), expected.len());
    for (peak, mass) in spectrum.peaks.iter().zip(expected) {
        near(peak.mz, mass, 0.001);
    }
}

#[test]
fn default_by_spectrum_and_first_prefix_match_upstream() {
    let p = peptide("IFSQVGK");
    let mut g = Generator::default();
    let s = g.generate(&p, 1, 1, None).unwrap();
    golden_masses(
        &s,
        vec![
            147.113, 204.135, 261.160, 303.203, 348.192, 431.262, 476.251, 518.294, 575.319,
            632.341, 665.362,
        ],
    );
    assert!(s.peaks.iter().all(|p| p.intensity == 1.0));
    assert_eq!(s.ms_level, 2);
    assert_eq!(s.spectrum_type, SpectrumType::Centroid);
    assert_eq!(s.precursors.len(), 1);
    assert_eq!(s.precursors[0].charge, 2);
    near(s.precursors[0].mz, p.mz(2).unwrap(), 1e-12);
    assert!(s.string_data_arrays.is_empty());
    assert!(s.integer_data_arrays.is_empty());
    let doubled = g.generate(&p, 1, 2, None).unwrap();
    assert_eq!(doubled.len(), 22);
    assert_eq!(doubled.precursors[0].charge, 3);
    g.add_first_prefix_ion = true;
    g.add_metainfo = true;
    let first = g.generate(&p, 1, 1, Some(4)).unwrap();
    assert_eq!(first.len(), 12);
    near(named(&first)["b1+"].mz, 114.09139, 0.001);
    assert!(charges(&first).iter().all(|&c| c == 1));
    assert_eq!(first.precursors[0].charge, 4);
}

#[test]
fn all_six_series_and_precursors_match_upstream_mass_table() {
    let g = Generator {
        ion_series: vec![Ion::A, Ion::B, Ion::C, Ion::X, Ion::Y, Ion::Z],
        add_first_prefix_ion: true,
        add_precursor_peaks: true,
        add_metainfo: true,
        ..Generator::default()
    };
    let s = g.generate(&peptide("DFPLANGER"), 1, 1, None).unwrap();
    golden_masses(
        &s,
        vec![
            88.03990, 235.10831, 332.16108, 445.24514, 516.28225, 630.32518, 687.34664, 816.38924,
            116.03481, 263.10323, 360.15599, 473.24005, 544.27717, 658.32009, 715.34156, 844.38415,
            1000.48526, 133.06136, 280.12978, 377.18254, 490.26660, 561.30372, 675.34664,
            732.36811, 861.41070, 929.44815, 782.37973, 685.32697, 572.24291, 501.20579, 387.16287,
            330.14140, 201.09881, 1018.49583, 903.46888, 756.40047, 659.34771, 546.26364,
            475.22653, 361.18360, 304.16214, 175.11955, 1001.46928, 886.44233, 739.37392,
            642.32116, 529.23709, 458.19998, 344.15705, 287.13559, 158.09300,
        ],
    );
    for series in ["a", "b", "c", "x", "y", "z"] {
        for n in 1..=8 {
            assert!(named(&s).contains_key(&format!("{series}{n}+")));
        }
    }
}

#[test]
fn radical_z_variants_differ_by_neutral_hydrogen() {
    let p = peptide("PEPTIDE");
    let g = Generator {
        ion_series: vec![Ion::Z, Ion::ZPlusOne, Ion::ZPlusTwo],
        add_metainfo: true,
        ..Generator::default()
    };
    let s = g.generate(&p, 2, 2, None).unwrap();
    let rows = named(&s);
    for n in 1..p.len() {
        let z = rows[&format!("z{n}++")].mz;
        let expected = (p.suffix(n).unwrap().mono_mass().unwrap() - formula("NH3").mono_mass()
            + 2.0 * PROTON_MASS_U)
            / 2.0;
        near(z, expected, 1e-10);
        near(
            rows[&format!("z.{n}++")].mz - z,
            formula("H").mono_mass() / 2.0,
            1e-10,
        );
        near(
            rows[&format!("z'{n}++")].mz - z,
            formula("H2").mono_mass() / 2.0,
            1e-10,
        );
    }
}

#[test]
fn loss_annotations_intensities_and_charge_ranges_match_upstream() {
    let mut g = Generator {
        ion_series: vec![Ion::B, Ion::X],
        add_first_prefix_ion: true,
        add_losses: true,
        add_precursor_peaks: true,
        add_metainfo: true,
        ..Generator::default()
    };
    g.intensities.b = 0.6;
    g.intensities.x = 0.4;
    g.relative_loss_intensity = 0.25;
    g.intensities.precursor = 0.8;
    g.intensities.precursor_water_loss = 0.3;
    g.intensities.precursor_ammonia_loss = 0.2;
    let s = g.generate(&peptide("IFSQVGK"), 1, 3, None).unwrap();
    assert_eq!(s.len(), 84);
    for (charge, count) in [(1, 27), (2, 27), (3, 30)] {
        assert_eq!(charges(&s).iter().filter(|&&c| c == charge).count(), count);
    }
    let rows = named(&s);
    let losses = [
        "x1-H3N1", "x2-H3N1", "x3-H3N1", "b3-H2O1", "x4-H3N1", "b4-H2O1", "b4-H3N1", "x5-H2O1",
        "x5-H3N1", "b5-H2O1", "b5-H3N1", "b6-H2O1", "b6-H3N1", "x6-H2O1", "x6-H3N1",
    ];
    for loss in losses {
        near(
            f64::from(rows[&format!("{loss}++")].intensity),
            if loss.starts_with('b') { 0.15 } else { 0.1 },
            1e-7,
        );
    }
    near(f64::from(rows["b2+++"].intensity), 0.6, 1e-7);
    near(f64::from(rows["[M+3H]+++"].intensity), 0.8, 1e-7);
    near(f64::from(rows["[M+3H-H2O]+++"].intensity), 0.3, 1e-7);
    near(f64::from(rows["[M+3H-NH3]+++"].intensity), 0.2, 1e-7);
    assert!(!rows.contains_key("[M+H]+"));
    g.add_all_precursor_charges = true;
    let all = g.generate(&peptide("IFSQVGK"), 1, 3, Some(3)).unwrap();
    assert_eq!(all.len(), 90);
    assert!(named(&all).contains_key("[M+H]+"));
}

#[test]
fn terminal_losses_and_modified_residue_losses_are_distinct() {
    let mut g = Generator {
        add_first_prefix_ion: true,
        add_losses: true,
        add_metainfo: true,
        ..Generator::default()
    };
    let unmodified = g.generate(&peptide("ASA"), 1, 1, None).unwrap();
    assert!(named(&unmodified).contains_key("b2-H2O1+"));
    let phospho = g.generate(&peptide("AS(Phospho)A"), 1, 1, None).unwrap();
    let rows = named(&phospho);
    assert!(rows.contains_key("b2-H3O4P1+"));
    assert!(rows.contains_key("y2-H3O4P1+"));
    assert!(!rows.contains_key("b2-H2O1+"));
    near(
        rows["b2+"].mz - rows["b2-H3O4P1+"].mz,
        formula("H3PO4").mono_mass(),
        1e-10,
    );
    let plain = g.generate(&peptide("AAA"), 1, 1, None).unwrap();
    assert_eq!(plain.len(), 4);
    g.add_terminal_losses = true;
    let terminal = g.generate(&peptide("AAA"), 1, 1, None).unwrap();
    assert_eq!(terminal.len(), 10);
    assert!(named(&terminal).contains_key("b1-H2O1+"));
    assert!(named(&terminal).contains_key("y1-H3N1+"));
}

#[test]
fn isotope_clusters_and_loss_probabilities_match_upstream() {
    let p = peptide("ARRGH");
    let mut g = Generator {
        ion_series: vec![Ion::Y],
        isotope_model: Isotopes::Coarse { max_peaks: 2 },
        add_metainfo: true,
        ..Generator::default()
    };
    let s = g.generate(&p, 2, 2, None).unwrap();
    let base = [78.54206, 107.05279, 185.10335, 263.15390];
    golden_masses(
        &s,
        base.iter()
            .flat_map(|m| [*m, m + C13C12_MASSDIFF_U / 2.0])
            .collect(),
    );
    assert_eq!(names(&s), ["y1", "y1", "y2", "y2", "y3", "y3", "y4", "y4"]);
    for pair in s.peaks.chunks_exact(2) {
        near(
            f64::from(pair[0].intensity) + f64::from(pair[1].intensity),
            1.0,
            1e-7,
        );
    }
    g.add_losses = true;
    let s = g.generate(&p, 1, 2, None).unwrap();
    let base = [
        156.07675, 213.09821, 325.18569, 327.17753, 352.17278, 369.19932, 481.28680, 483.27864,
        508.27389, 525.30044,
    ];
    golden_masses(
        &s,
        base.iter()
            .flat_map(|m| {
                [
                    *m,
                    m + C13C12_MASSDIFF_U,
                    (m + PROTON_MASS_U) / 2.0,
                    (m + PROTON_MASS_U + C13C12_MASSDIFF_U) / 2.0,
                ]
            })
            .collect(),
    );
    near(f64::from(s.peaks[0].intensity), 0.927642, 1e-6);
    near(f64::from(s.peaks[1].intensity), 0.0723581, 1e-6);
    assert!(names(&s).iter().any(|n| n == "y3-C1H2N1O1++"));
    g.ion_series.clear();
    g.add_precursor_peaks = true;
    let precursors = g.generate(&p, 2, 2, None).unwrap();
    golden_masses(
        &precursors,
        [578.32698, 579.31100, 596.33755]
            .iter()
            .flat_map(|m| {
                [
                    (m + PROTON_MASS_U) / 2.0,
                    (m + PROTON_MASS_U + C13C12_MASSDIFF_U) / 2.0,
                ]
            })
            .collect(),
    );
}

#[test]
fn invalid_loss_formulas_are_skipped_in_both_modes() {
    let mut g = Generator {
        ion_series: vec![Ion::A],
        add_first_prefix_ion: true,
        add_losses: true,
        add_metainfo: true,
        ..Generator::default()
    };
    for model in [Isotopes::None, Isotopes::Coarse { max_peaks: 2 }] {
        g.isotope_model = model;
        let s = g.generate(&peptide("RDK"), 1, 1, None).unwrap();
        assert!(!names(&s).iter().any(|n| n == "a1-C1H2N1O1+"));
        assert!(names(&s).iter().any(|n| n == "a2-C1H2N1O1+"));
    }
    // The source's regression fixture covers negative CONH2 loss from a1(R).
    g.ion_series = vec![Ion::A, Ion::B, Ion::Y];
    let s = g.generate(&peptide("RDAGGPALKK"), 1, 1, None).unwrap();
    assert_eq!(s.len(), 212);
    g.add_first_prefix_ion = false;
    assert_eq!(
        g.generate(&peptide("RDAGGPALKK"), 1, 1, None)
            .unwrap()
            .len(),
        198
    );
}

#[test]
fn modified_peptides_retain_terminal_deltas_and_isotope_labels() {
    let p = peptide(".(Acetyl)ACDK(Label:13C(6)15N(2)).(Amidated)");
    let mut g = Generator {
        ion_series: vec![Ion::A, Ion::B, Ion::C, Ion::X, Ion::Y, Ion::Z],
        add_first_prefix_ion: true,
        add_metainfo: true,
        ..Generator::default()
    };
    let mono = g.generate(&p, 2, 2, None).unwrap();
    let rows = named(&mono);
    for &series in &g.ion_series {
        for n in 1..p.len() {
            let fragment = if series.is_prefix() {
                p.prefix(n)
            } else {
                p.suffix(n)
            }
            .unwrap();
            let expected = (fragment.mono_mass().unwrap()
                + series.formula_delta().mono_mass()
                + 2.0 * PROTON_MASS_U)
                / 2.0;
            near(
                rows[&format!("{}{n}++", series.label())].mz,
                expected,
                1e-10,
            );
        }
    }
    g.isotope_model = Isotopes::Coarse { max_peaks: 3 };
    let coarse = g.generate(&p, 2, 2, None).unwrap();
    let independently =
        CoarseIsotopePatternGenerator::new(Some(3), CoarseMassMode::Approximate).unwrap();
    for &series in &g.ion_series {
        for n in 1..p.len() {
            let fragment = if series.is_prefix() {
                p.prefix(n)
            } else {
                p.suffix(n)
            }
            .unwrap();
            let base = fragment
                .formula()
                .unwrap()
                .checked_add(&series.formula_delta())
                .unwrap();
            let expected = independently
                .run(&base.checked_add(&formula("H2")).unwrap())
                .unwrap();
            let selected: Vec<_> = coarse
                .peaks
                .iter()
                .zip(names(&coarse))
                .filter(|(_, name)| *name == &format!("{}{n}", series.label()))
                .map(|(peak, _)| peak)
                .collect();
            assert_eq!(selected.len(), 3);
            for (actual, expected) in selected.iter().zip(expected.peaks()) {
                near(actual.mz, expected.mass / 2.0, 1e-12);
                near(f64::from(actual.intensity), expected.probability, 1e-7);
            }
            // Coarse uses elemental formula and neutral H; mono keeps declared
            // terminal delta mass. This source distinction is observable at 1e-8.
            let delta = fragment.mono_mass().unwrap() - fragment.formula().unwrap().mono_mass();
            let expected_shift = formula("H").mono_mass() - PROTON_MASS_U - delta / 2.0;
            near(
                selected[0].mz - rows[&format!("{}{n}++", series.label())].mz,
                expected_shift,
                1e-10,
            );
        }
    }
}

#[test]
fn precursor_losses_are_independent_of_fragment_losses() {
    let p = peptide(".(Acetyl)PEP.(Amidated)");
    let mut g = Generator {
        ion_series: vec![],
        add_precursor_peaks: true,
        add_metainfo: true,
        ..Generator::default()
    };
    let s = g.generate(&p, 1, 2, Some(5)).unwrap();
    assert_eq!(s.len(), 3);
    assert!(charges(&s).iter().all(|&c| c == 2));
    assert_eq!(s.precursors[0].charge, 5);
    let rows = named(&s);
    near(rows["[M+2H]++"].mz, p.mz(2).unwrap(), 1e-12);
    near(
        rows["[M+2H-H2O]++"].mz,
        (p.formula()
            .unwrap()
            .checked_sub(&formula("H2O"))
            .unwrap()
            .mono_mass()
            + 2.0 * PROTON_MASS_U)
            / 2.0,
        1e-12,
    );
    g.add_all_precursor_charges = true;
    assert_eq!(g.generate(&p, 1, 2, None).unwrap().len(), 6);
}

#[test]
fn unsorted_emission_uses_source_series_order_and_sort_keeps_annotations() {
    let p = peptide("PEPTIDE");
    let g = Generator {
        ion_series: vec![Ion::A, Ion::Y, Ion::B],
        add_metainfo: true,
        sort_by_position: false,
        ..Generator::default()
    };
    let unsorted = g.generate(&p, 1, 2, None).unwrap();
    assert_eq!(&names(&unsorted)[..5], ["b2+", "b3+", "b4+", "b5+", "b6+"]);
    assert_eq!(
        &names(&unsorted)[5..11],
        ["y1+", "y2+", "y3+", "y4+", "y5+", "y6+"]
    );
    let mut expected = unsorted.clone();
    expected.sort_by_position().unwrap();
    let sorted = Generator {
        sort_by_position: true,
        ..g
    }
    .generate(&p, 1, 2, None)
    .unwrap();
    assert_eq!(sorted, expected);
}

#[test]
fn append_is_atomic_and_extends_named_annotations() {
    let g = Generator {
        add_metainfo: true,
        ..Generator::default()
    };
    let mut s = MSSpectrum {
        peaks: vec![Peak1D::new(5000.0, 9.0), Peak1D::new(50.0, 7.0)],
        rt: 7.5,
        name: "retained".into(),
        native_id: "scan=1".into(),
        precursors: vec![Precursor::new(999.0, 1)],
        ..MSSpectrum::default()
    };
    s.metadata.insert("source".into(), "retained".into());
    s.float_data_arrays
        .push(DataArray::new("placeholder", vec![]));
    g.append_to(&mut s, &peptide("PEPTIDE"), 1, 1, None)
        .unwrap();
    assert_eq!(s.len(), 13);
    assert_eq!(s.rt, 7.5);
    assert_eq!(s.name, "retained");
    assert_eq!(s.metadata["source"], "retained");
    assert_eq!(s.precursors.len(), 2);
    assert_eq!(names(&s)[0], "");
    assert_eq!(charges(&s)[0], 0);
    assert_eq!(names(&s)[12], "");
    let old = s.clone();
    Generator::default()
        .append_to(&mut s, &peptide("PEPTIDE"), 2, 2, None)
        .unwrap();
    assert_eq!(s.len(), 24);
    assert_eq!(names(&s).len(), 24);
    assert!(names(&s).iter().any(|n| n == "b2++"));
    s.validate().unwrap();
    let before = s.clone();
    g.append_to(&mut s, &peptide(""), 1, 1, None).unwrap();
    assert_eq!(s, before);
    assert!(g.append_to(&mut s, &peptide("PEP"), 1, 3, Some(1)).is_err());
    assert_eq!(s, before);
    let mut bad = old.clone();
    bad.float_data_arrays
        .push(DataArray::new("measured values", vec![0.0; bad.len()]));
    let before = bad.clone();
    assert!(g.append_to(&mut bad, &peptide("PEP"), 1, 1, None).is_err());
    assert_eq!(bad, before);
    bad.string_data_arrays
        .push(DataArray::new("IonNames", vec![]));
    let before = bad.clone();
    assert!(g.append_to(&mut bad, &peptide("PEP"), 1, 1, None).is_err());
    assert_eq!(bad, before);
    let mut bad = old;
    bad.integer_data_arrays[0].data.pop();
    let before = bad.clone();
    assert!(g.append_to(&mut bad, &peptide("PEP"), 1, 1, None).is_err());
    assert_eq!(bad, before);
}

#[test]
fn empty_monomer_and_invalid_configuration_boundaries() {
    let g = Generator::default();
    assert_eq!(
        g.generate(&peptide(""), 1, 1, None).unwrap(),
        MSSpectrum::default()
    );
    assert!(g.generate(&peptide("R"), 1, 1, None).unwrap().is_empty());
    assert_eq!(
        Generator {
            add_precursor_peaks: true,
            ..g.clone()
        }
        .generate(&peptide("R"), 1, 1, None)
        .unwrap()
        .len(),
        3
    );
    for ion in [Ion::C, Ion::X] {
        assert!(
            Generator {
                ion_series: vec![ion],
                ..g.clone()
            }
            .generate(&peptide("R"), 1, 1, None)
            .is_err()
        );
    }
    let p = peptide("PEPTIDE");
    for (min, max, precursor) in [(0, 1, None), (2, 1, None), (1, 2, Some(1)), (1, 1, Some(0))] {
        assert!(g.generate(&p, min, max, precursor).is_err());
    }
    for bad in [
        Generator {
            ion_series: vec![Ion::B, Ion::B],
            ..g.clone()
        },
        Generator {
            relative_loss_intensity: f32::NAN,
            ..g.clone()
        },
        Generator {
            relative_loss_intensity: 1.1,
            ..g.clone()
        },
        Generator {
            isotope_model: Isotopes::Coarse { max_peaks: 0 },
            ..g.clone()
        },
        Generator {
            isotope_model: Isotopes::Coarse {
                max_peaks: MAX_THEORETICAL_PEAKS + 1,
            },
            ..g.clone()
        },
        Generator {
            isotope_model: Isotopes::Coarse {
                max_peaks: MAX_THEORETICAL_PEAKS,
            },
            ..g.clone()
        },
        Generator {
            add_terminal_losses: true,
            ..g.clone()
        },
        Generator {
            add_terminal_losses: true,
            add_losses: true,
            isotope_model: Isotopes::Coarse { max_peaks: 2 },
            ..g.clone()
        },
    ] {
        assert!(bad.generate(&p, 1, 1, None).is_err());
    }
    assert!(
        g.generate(
            &peptide(&"A".repeat(MAX_THEORETICAL_RESIDUES + 1)),
            1,
            1,
            None
        )
        .is_err()
    );
    let maximum_charge = g.generate(&peptide("AA"), 255, 255, None).unwrap();
    assert_eq!(maximum_charge.precursors[0].charge, 256);
}
