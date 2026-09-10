// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Independent literal and analytical references from pinned OpenMS source.
//! See data/theoretical_extension_provenance.json. No C++ runtime is needed.

use openms::chemistry::{
    AASequence, TheoreticalIonSeries as Ion, TheoreticalIsotopeModel as Isotopes,
    TheoreticalSpectrumGenerator as Generator,
};
use openms::kernel::MSSpectrum;
use openms::metadata::ActivationMethod;
use std::collections::BTreeMap;

// ElementDB.cpp literals and Constants.h, independent of Rust chemistry getters.
const H: f64 = 1.007_825_031_9;
const C: f64 = 12.0;
const N: f64 = 14.003_074;
#[allow(clippy::excessive_precision)] // Keep the pinned source literal for audit.
const O: f64 = 15.994_915_000_000_001;
const PROTON: f64 = 1.007_276_466_771;
const WATER: f64 = 2.0 * H + O;
const CO: f64 = C + O;
const AMMONIA: f64 = 3.0 * H + N;
// ResidueDB declares the free forms C2H5NO2, C3H7NO2, C6H14N2O2.
const GLYCINE_INTERNAL: f64 = ((5.0 * H + 2.0 * C) + N + 2.0 * O) - WATER;
const ALANINE_INTERNAL: f64 = ((7.0 * H + 3.0 * C) + N + 2.0 * O) - WATER;
const LYSINE_INTERNAL: f64 = ((14.0 * H + 6.0 * C) + 2.0 * N + 2.0 * O) - WATER;

fn peptide(text: &str) -> AASequence {
    text.parse().unwrap()
}
fn internal_generator() -> Generator {
    Generator {
        ion_series: vec![],
        add_internal_fragments: true,
        add_metainfo: true,
        sort_by_position: false,
        ..Default::default()
    }
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
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.12} != {expected:.12} (tol {tolerance})"
    );
}

#[test]
fn internal_source_intervals_and_raw_labels_have_independent_analytical_masses() {
    let mut generator = internal_generator();
    generator.intensities.b = 2.;
    generator.intensities.a = 3.;
    let spectrum = generator.generate(&peptide("AGGGA"), 1, 2, None).unwrap();
    // l=1 is the only admitted start. l=2 would provide another GG interval,
    // but source l+3<n deliberately excludes it. Neither terminal A is retained.
    assert_eq!(
        names(&spectrum),
        [
            "GG", "GGG", "GG-CO", "GGG-CO", "GG", "GGG", "GG-CO", "GGG-CO"
        ]
    );
    assert_eq!(charges(&spectrum), [1, 1, 1, 1, 2, 2, 2, 2]);
    for (row, charge) in spectrum.peaks.chunks_exact(4).zip([1., 2.]) {
        for (peak, mass, intensity) in [
            (
                &row[0],
                (2. * GLYCINE_INTERNAL + charge * PROTON) / charge,
                2.,
            ),
            (
                &row[1],
                (3. * GLYCINE_INTERNAL + charge * PROTON) / charge,
                2.,
            ),
            (
                &row[2],
                (2. * GLYCINE_INTERNAL + charge * PROTON - CO) / charge,
                3.,
            ),
            (
                &row[3],
                (3. * GLYCINE_INTERNAL + charge * PROTON - CO) / charge,
                3.,
            ),
        ] {
            near(peak.mz, mass, 1e-9);
            assert_eq!(peak.intensity, intensity);
        }
    }
}

#[test]
fn internal_short_inputs_and_ten_residue_cap_match_source_loop_bounds() {
    let generator = internal_generator();
    for text in ["", "A", "AA", "AAA", "AAAA"] {
        assert!(
            generator
                .generate(&peptide(text), 1, 1, None)
                .unwrap()
                .is_empty(),
            "{text}"
        );
    }
    let spectrum = generator
        .generate(&peptide("AAAAAAAAAAAAAA"), 1, 1, None)
        .unwrap();
    assert_eq!(spectrum.len(), 124);
    let mut lengths = BTreeMap::<usize, usize>::new();
    for name in names(&spectrum) {
        let text = name.strip_suffix("-CO").unwrap_or(name);
        assert!(text.chars().all(|c| c == 'A'));
        *lengths.entry(text.len()).or_default() += 1;
    }
    // For each length k, min(10,13-k) starts, for each of b and a.
    assert_eq!(
        lengths,
        BTreeMap::from([
            (2, 20),
            (3, 20),
            (4, 18),
            (5, 16),
            (6, 14),
            (7, 12),
            (8, 10),
            (9, 8),
            (10, 6)
        ])
    );
    assert_eq!(
        generator
            .generate(&peptide("PEPTIDEK"), 1, 1, None)
            .unwrap()
            .len(),
        28
    );
}

#[test]
fn internal_loss_collection_skips_the_first_retained_residue_and_keeps_source_order() {
    let mut generator = internal_generator();
    generator.add_losses = true;
    generator.relative_loss_intensity = 0.25;
    generator.intensities.b = 2.;
    generator.intensities.a = 3.;
    let first = generator.generate(&peptide("AKGGA"), 1, 1, None).unwrap();
    assert_eq!(names(&first), ["KG", "KGG", "KG-CO", "KGG-CO"]);
    let second = generator.generate(&peptide("AGKGA"), 1, 1, None).unwrap();
    assert_eq!(
        names(&second),
        [
            "GK",
            "GKG",
            "GK-H3N1+",
            "GKG-H3N1+",
            "GK-CO",
            "GKG-CO",
            "GK-CO-H3N1+",
            "GKG-CO-H3N1+",
        ]
    );
    let base = PROTON + GLYCINE_INTERNAL + LYSINE_INTERNAL;
    for (index, expected, intensity) in [
        (0, base, 2.),
        (1, base + GLYCINE_INTERNAL, 2.),
        (2, base - AMMONIA, 0.5),
        (3, base + GLYCINE_INTERNAL - AMMONIA, 0.5),
        (4, base - CO, 3.),
        (5, base + GLYCINE_INTERNAL - CO, 3.),
        (6, base - CO - AMMONIA, 0.75),
        (7, base + GLYCINE_INTERNAL - CO - AMMONIA, 0.75),
    ] {
        near(second.peaks[index].mz, expected, 1e-9);
        assert_eq!(second.peaks[index].intensity, intensity);
    }
}

#[test]
fn internal_isotope_and_terminal_loss_flags_do_not_change_the_source_output() {
    let mut generator = internal_generator();
    generator.add_losses = true;
    let input = peptide("AGKGA");
    let base = generator.generate(&input, 1, 1, None).unwrap();
    generator.isotope_model = Isotopes::Coarse { max_peaks: 4 };
    let isotope_setting = generator.generate(&input, 1, 1, None).unwrap();
    assert_eq!(isotope_setting.peaks, base.peaks);
    assert_eq!(names(&isotope_setting), names(&base));
    generator.isotope_model = Isotopes::None;
    generator.add_terminal_losses = true;
    let terminal_setting = generator.generate(&input, 1, 1, None).unwrap();
    assert_eq!(terminal_setting.peaks, base.peaks);
    assert_eq!(names(&terminal_setting), names(&base));
}

#[test]
fn internal_terminal_modifications_are_excluded_and_mass_only_residues_are_retained() {
    let generator = internal_generator();
    let base = generator.generate(&peptide("AGGGA"), 1, 1, None).unwrap();
    let terminal = generator
        .generate(&peptide(".(Acetyl)AGGGA.(Amidated)"), 1, 1, None)
        .unwrap();
    assert_eq!(terminal.peaks, base.peaks);
    assert_eq!(names(&terminal), names(&base));
    let tagged = generator
        .generate(&peptide("AX[999]GGA"), 1, 1, None)
        .unwrap();
    assert_eq!(names(&tagged), ["XG", "XGG", "XG-CO", "XGG-CO"]);
    for (peak, expected) in tagged.peaks.iter().zip([
        PROTON + 999. + GLYCINE_INTERNAL,
        PROTON + 999. + 2. * GLYCINE_INTERNAL,
        PROTON + 999. + GLYCINE_INTERNAL - CO,
        PROTON + 999. + 2. * GLYCINE_INTERNAL - CO,
    ]) {
        near(peak.mz, expected, 1e-9);
    }
}

#[test]
fn source_output_groups_put_internal_between_terminal_and_precursor_peaks() {
    let generator = Generator {
        ion_series: vec![Ion::B],
        add_precursor_peaks: true,
        add_abundant_immonium_ions: true,
        ..internal_generator()
    };
    let spectrum = generator.generate(&peptide("PGGGP"), 2, 2, None).unwrap();
    assert_eq!(
        names(&spectrum),
        [
            "b2++",
            "b3++",
            "b4++",
            "GG",
            "GGG",
            "GG-CO",
            "GGG-CO",
            "[M+2H]++",
            "[M+2H-H2O]++",
            "[M+2H-NH3]++",
            "iP+",
        ]
    );
    assert_eq!(charges(&spectrum), [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1]);
}

#[test]
fn literal_source_immonium_and_activation_masses() {
    let generator = Generator {
        ion_series: vec![],
        add_abundant_immonium_ions: true,
        add_metainfo: true,
        ..Default::default()
    };
    let spectrum = generator.generate(&peptide("HFYLWCP"), 2, 3, None).unwrap();
    assert_eq!(
        names(&spectrum),
        ["iP+", "iC+", "iL/I+", "iH+", "iF+", "iY+", "iW+"]
    );
    let expected: [f64; 7] = [
        70.0656, 76.0221, 86.09698, 110.0718, 120.0813, 136.0762, 159.0922,
    ];
    for (peak, expected) in spectrum.peaks.iter().zip(expected) {
        assert_eq!(peak.mz.to_bits(), expected.to_bits());
        assert_eq!(peak.intensity, 1.);
    }
    assert_eq!(charges(&spectrum), [1; 7]);

    let cid =
        Generator::generate_for_activation(ActivationMethod::Cid, &peptide("HFYLWCP"), 1).unwrap();
    let cid_expected = [
        116.0706, 219.0797, 285.1346, 405.1591, 448.1979, 518.2431, 561.2819, 681.3064, 747.3613,
        828.3749, 850.3704,
    ];
    assert_eq!(cid.len(), cid_expected.len());
    for (peak, expected) in cid.peaks.iter().zip(cid_expected) {
        near(peak.mz, expected, 0.0001);
    }
    let ecd =
        Generator::generate_for_activation(ActivationMethod::Ecd, &peptide("HFYLWCP"), 1).unwrap();
    assert_eq!(ecd.len(), 17);
    // The pinned class test only gives the first eleven of seventeen ECD masses.
    for (peak, expected) in ecd.peaks.iter().zip([
        100.0518816,
        101.0597067,
        203.0610665,
        204.0688916,
        302.1611520,
        389.1403798,
        390.1482049,
        465.2244813,
        502.2244442,
        503.2322692,
        578.3085457,
    ]) {
        near(peak.mz, expected, 0.000001);
    }
}

#[test]
fn independent_float_helper_oracle_includes_first_and_full_length_ions_and_existing_output() {
    let generator = Generator::default();
    let mut actual = vec![999., 1.];
    generator
        .append_mass_spectrum(&mut actual, &peptide("AG"), 2)
        .unwrap();
    // Source emits b1,b2,y1,y2 for both charges, converts every result to f32,
    // then sorts together with the existing caller values. No generator call
    // supplies these expected masses.
    let mut expected = vec![999_f32, 1.];
    for charge in [2., 1.] {
        let b1 = charge * PROTON + ALANINE_INTERNAL;
        let b2 = b1 + GLYCINE_INTERNAL;
        let y1 = charge * PROTON + WATER + GLYCINE_INTERNAL;
        let y2 = y1 + ALANINE_INTERNAL;
        expected.extend([b1, b2, y1, y2].map(|mass| (mass / charge) as f32));
    }
    expected.sort_by(f32::total_cmp);
    assert_eq!(
        actual.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        expected.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
    );
}

#[test]
fn float_helper_boolean_series_and_sorted_noops_follow_its_separate_contract() {
    let baseline = Generator::default();
    let settings = Generator {
        ion_series: vec![Ion::B, Ion::B, Ion::Y, Ion::ZPlusOne, Ion::ZPlusTwo],
        add_losses: true,
        add_terminal_losses: true,
        add_internal_fragments: true,
        add_abundant_immonium_ions: true,
        add_precursor_peaks: true,
        isotope_model: Isotopes::Coarse { max_peaks: 3 },
        ..Default::default()
    };
    let mut expected = Vec::new();
    baseline
        .append_mass_spectrum(&mut expected, &peptide("AG"), 2)
        .unwrap();
    let mut actual = Vec::new();
    settings
        .append_mass_spectrum(&mut actual, &peptide("AG"), 2)
        .unwrap();
    assert_eq!(actual, expected);
    for (settings, input, charge) in [
        (baseline.clone(), peptide("X"), 0),
        (baseline, peptide(""), 2),
        (
            Generator {
                ion_series: vec![Ion::ZPlusOne],
                ..Default::default()
            },
            peptide("X"),
            2,
        ),
    ] {
        let mut output = vec![3., 1., 2.];
        settings
            .append_mass_spectrum(&mut output, &input, charge)
            .unwrap();
        assert_eq!(output, [1., 2., 3.]);
    }
}
