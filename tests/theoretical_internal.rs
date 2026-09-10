// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::theoretical::{MAX_THEORETICAL_PEAKS, MAX_THEORETICAL_RESIDUES};
use openms::chemistry::{
    AASequence, TheoreticalIonSeries as Ion, TheoreticalIsotopeModel as Isotopes,
    TheoreticalSpectrumGenerator as Generator,
};
use openms::kernel::{DataArray, MSSpectrum, Peak1D};

fn peptide(text: &str) -> AASequence {
    text.parse().unwrap()
}
fn generator() -> Generator {
    Generator {
        ion_series: vec![],
        add_internal_fragments: true,
        add_metainfo: true,
        sort_by_position: false,
        ..Default::default()
    }
}
fn names(s: &MSSpectrum) -> &[String] {
    &s.string_data_arrays[0].data
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a} != {b}");
}

#[test]
fn internal_boundaries_and_ten_residue_cap_are_independent_of_terminal_series() {
    let g = generator();
    for (n, count) in [
        (0, 0),
        (1, 0),
        (2, 0),
        (3, 0),
        (4, 0),
        (5, 4),
        (6, 10),
        (14, 124),
    ] {
        let s = g.generate(&peptide(&"G".repeat(n)), 1, 1, None).unwrap();
        assert_eq!(s.len(), count, "length {n}");
        if n > 0 {
            assert!(
                names(&s)
                    .iter()
                    .all(|name| (2..=10).contains(&name.trim_end_matches("-CO").len()))
            );
        }
    }
    let s = g.generate(&peptide("AGGGA"), 1, 2, None).unwrap();
    assert_eq!(
        names(&s),
        [
            "GG", "GGG", "GG-CO", "GGG-CO", "GG", "GGG", "GG-CO", "GGG-CO"
        ]
    );
    assert_eq!(s.integer_data_arrays[0].data, [1, 1, 1, 1, 2, 2, 2, 2]);
}

#[test]
fn internal_masses_ignore_terminal_modifications_and_keep_numeric_residue_tags() {
    let g = generator();
    let base = g.generate(&peptide("AGGGA"), 1, 1, None).unwrap();
    let terminal = g
        .generate(&peptide(".(Acetyl)AGGGA.(Amidated)"), 1, 1, None)
        .unwrap();
    assert_eq!(base.peaks, terminal.peaks);
    let tagged = g
        .generate(&peptide("AGG[+12.3456789]GA"), 1, 1, None)
        .unwrap();
    for (p, b) in tagged.peaks.iter().zip(&base.peaks) {
        near(p.mz - b.mz, 12.3456789);
    }
    assert_eq!(names(&base), names(&tagged));
    assert!(g.generate(&peptide("AGXGA"), 1, 1, None).is_err());
    assert!(
        g.generate(&peptide("AGX[111.23456789]GA"), 1, 1, None)
            .is_ok()
    );
}

#[test]
fn losses_skip_the_first_internal_residue_and_replace_modified_residue_losses() {
    let mut g = generator();
    g.add_losses = true;
    let first = g.generate(&peptide("AKGGA"), 1, 1, None).unwrap();
    assert_eq!(first.len(), 4);
    let second = g.generate(&peptide("AGKGA"), 1, 1, None).unwrap();
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
            "GKG-CO-H3N1+"
        ]
    );
    let shifted = g
        .generate(&peptide("AGK[+12.3456789]GA"), 1, 1, None)
        .unwrap();
    assert_eq!(shifted.len(), 4);
    let phospho = g.generate(&peptide("AAS(Phospho)GA"), 1, 1, None).unwrap();
    assert!(names(&phospho).iter().any(|n| n.contains("H3O4P1")));
    assert!(!names(&phospho).iter().any(|n| n.contains("H2O1")));
}

#[test]
fn internal_ions_use_series_intensities_and_stay_monoisotopic() {
    let mut g = generator();
    g.add_losses = true;
    g.intensities.b = 0.5;
    g.intensities.a = 0.25;
    g.relative_loss_intensity = 0.2;
    let p = peptide("AGKGA");
    let base = g.generate(&p, 2, 2, None).unwrap();
    assert_eq!(
        base.peaks.iter().map(|p| p.intensity).collect::<Vec<_>>(),
        [0.5, 0.5, 0.1, 0.1, 0.25, 0.25, 0.05, 0.05]
    );
    g.add_terminal_losses = true;
    assert_eq!(base, g.generate(&p, 2, 2, None).unwrap());
    g.add_terminal_losses = false;
    g.isotope_model = Isotopes::Coarse { max_peaks: 4 };
    assert_eq!(base, g.generate(&p, 2, 2, None).unwrap());
    let mass_only = peptide("AGK[+12.3456789]GA");
    let coarse = g.generate(&mass_only, 2, 2, None).unwrap();
    g.isotope_model = Isotopes::None;
    assert_eq!(coarse, g.generate(&mass_only, 2, 2, None).unwrap());
}

#[test]
fn appended_internal_peaks_keep_stable_annotations_and_fail_atomically_at_limits() {
    let mut g = generator();
    g.sort_by_position = true;
    let p = peptide("AGGGA");
    let mut s = MSSpectrum {
        peaks: vec![Peak1D::new(0.5, 2.0)],
        string_data_arrays: vec![DataArray::new("IonNames", vec!["observed".into()])],
        integer_data_arrays: vec![DataArray::new("Charges", vec![4])],
        ..Default::default()
    };
    g.append_to(&mut s, &p, 1, 1, None).unwrap();
    assert_eq!(s.len(), 5);
    assert_eq!(names(&s)[0], "observed");
    assert_eq!(s.integer_data_arrays[0].data[0], 4);
    s.validate().unwrap();
    let before = s.clone();
    assert!(
        g.append_to(
            &mut s,
            &peptide(&"G".repeat(MAX_THEORETICAL_RESIDUES)),
            1,
            2,
            None
        )
        .is_err()
    );
    assert_eq!(s, before);
    let mut full = MSSpectrum {
        peaks: vec![Peak1D::new(10.0, 1.0); MAX_THEORETICAL_PEAKS - 3],
        ..Default::default()
    };
    let before = full.clone();
    assert!(g.append_to(&mut full, &p, 1, 1, None).is_err());
    assert_eq!(full, before);
}

#[test]
fn output_groups_follow_terminal_then_internal_then_precursor_then_immonium_order() {
    let mut g = generator();
    g.ion_series = vec![Ion::B];
    g.add_precursor_peaks = true;
    g.add_abundant_immonium_ions = true;
    let s = g.generate(&peptide("AGGGP"), 1, 1, None).unwrap();
    assert_eq!(
        names(&s),
        [
            "b2+",
            "b3+",
            "b4+",
            "GG",
            "GGG",
            "GG-CO",
            "GGG-CO",
            "[M+H]+",
            "[M+H-H2O]+",
            "[M+H-NH3]+",
            "iP+"
        ]
    );
}

#[test]
fn custom_loss_work_and_template_storage_are_bounded_before_output_changes() {
    use openms::chemistry::{
        EmpiricalFormula, ModificationRecord, ModificationsDB, NeutralLoss, ResidueModification,
    };
    let loss = |text: &str| {
        let f: EmpiricalFormula = text.parse().unwrap();
        NeutralLoss::new(f.clone(), f.mono_mass(), f.average_mass()).unwrap()
    };
    let database = |neutral_losses| {
        ModificationsDB::from_records(vec![
            ResidueModification::from_record(ModificationRecord {
                name: "Bulk".into(),
                origin: Some('K'),
                neutral_losses,
                ..Default::default()
            })
            .unwrap(),
        ])
        .unwrap()
    };
    let repeated = database(vec![loss("H2O"); 10_000]);
    let p = AASequence::parse_with_registry("AGK(Bulk)GA", &repeated).unwrap();
    let g = Generator {
        add_losses: true,
        ..generator()
    };
    let mut output = MSSpectrum {
        peaks: vec![Peak1D::new(12.0, 3.0)],
        ..Default::default()
    };
    let before = output.clone();
    let error = g.append_to(&mut output, &p, 1, 255, None).unwrap_err();
    assert!(error.to_string().contains("declaration work"), "{error}");
    assert_eq!(output, before);
    let terminal = Generator {
        add_losses: true,
        ..Default::default()
    };
    let p = AASequence::parse_with_registry(&"K(Bulk)".repeat(50), &repeated).unwrap();
    let error = terminal.append_to(&mut output, &p, 1, 1, None).unwrap_err();
    assert!(error.to_string().contains("declaration work"), "{error}");
    assert_eq!(output, before);
    // Distinct, mostly impossible losses still consume storage before physical
    // filtering. One middle residue visits more than 100 terminal templates.
    let distinct = database((1..=1000).map(|n| loss(&format!("H{n}"))).collect());
    let p = AASequence::parse_with_registry(
        &format!("{}K(Bulk){}", "A".repeat(110), "A".repeat(109)),
        &distinct,
    )
    .unwrap();
    let error = terminal.append_to(&mut output, &p, 1, 1, None).unwrap_err();
    assert!(error.to_string().contains("template storage"), "{error}");
    assert_eq!(output, before);
}
