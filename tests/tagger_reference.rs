// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source literal expectations and independent control-flow/mass constructions.
// See data/tagger_provenance.json; no C++ execution generated these references.

use openms::chemistry::tagger::{Tagger, TaggerOptions};
use openms::chemistry::{
    AASequence, ModificationRecord, ModificationsDB, ResidueModification, TheoreticalIonSeries,
    TheoreticalSpectrumGenerator,
};
use openms::comparison::Tolerance;

fn options(minimum: usize, maximum: usize, tolerance: Tolerance) -> TaggerOptions {
    let mut options = TaggerOptions::new(minimum, tolerance);
    options.max_tag_length = maximum;
    options
}

fn source_spectrum(modified: bool) -> openms::MSSpectrum {
    let generator = TheoreticalSpectrumGenerator {
        ion_series: vec![
            TheoreticalIonSeries::A,
            TheoreticalIonSeries::B,
            TheoreticalIonSeries::Y,
        ],
        add_first_prefix_ion: true,
        add_losses: true,
        add_precursor_peaks: true,
        add_metainfo: false,
        ..Default::default()
    };
    let peptide = AASequence::parse(if modified {
        "PEPTID(Oxidation)ETESTTHISTAGGER"
    } else {
        "PEPTIDETESTTHISTAGGER"
    })
    .unwrap();
    let spectrum = generator
        .generate(&peptide, if modified { 2 } else { 1 }, 2, None)
        .unwrap();
    assert_eq!(spectrum.peaks.len(), if modified { 180 } else { 357 });
    spectrum
}

#[test]
fn all_six_source_counts_and_120_literal_membership_assertions() {
    let original = source_spectrum(false);
    let modified = source_spectrum(true);
    let contexts = [
        ("unmodified_q1", &original, 1, 1, 0),
        ("unmodified_q2", &original, 2, 2, 0),
        ("unmodified_q1_q2", &original, 1, 2, 0),
        ("modified_unaware", &modified, 1, 2, 0),
        ("modified_fixed", &modified, 2, 2, 1),
        ("modified_variable", &modified, 2, 2, 2),
    ];
    let mut assertions = 0;
    for (name, spectrum, min_charge, max_charge, modification) in contexts {
        let mut config = options(2, 5, Tolerance::Ppm(10.0));
        config.min_charge = min_charge;
        config.max_charge = max_charge;
        if modification == 1 {
            config.fixed_mods.push("Oxidation (D)".into());
        }
        if modification == 2 {
            config.variable_mods.push("Oxidation (D)".into());
        }
        let tagger = Tagger::new(config).unwrap();
        let tags = tagger.get_spectrum_tags(spectrum).unwrap();
        assert!(tags.windows(2).all(|pair| pair[0] < pair[1]));
        let expected = include_str!("data/tagger_source_counts.tsv")
            .lines()
            .skip(1)
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .find(|row| row[0] == name)
            .unwrap();
        assert_eq!(
            tags.len(),
            expected[1].parse::<usize>().unwrap(),
            "{name}, source line {}",
            expected[2]
        );
        for row in include_str!("data/tagger_source_membership.tsv")
            .lines()
            .skip(1)
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .filter(|row| row[0] == name)
        {
            let found = tags
                .binary_search_by(|tag| tag.as_str().cmp(row[1]))
                .is_ok();
            assert_eq!(
                found,
                row[2] == "true",
                "{name}: {}, source line {}",
                row[1],
                row[3]
            );
            assertions += 1;
        }
        let positions: Vec<_> = spectrum.peaks.iter().map(|peak| peak.mz).collect();
        assert_eq!(tagger.get_tags(&positions).unwrap(), tags);
    }
    assert_eq!(assertions, 120);
}

#[test]
fn literal_small_traces_preserve_source_da_and_ppm_claims() {
    // Original arithmetic expressions from Tagger_test.cpp, not rounded masses
    // extracted from a native generated spectrum.
    let trace = [
        150.0,
        150.0 + 97.0527,
        150.0 + 97.0527 + 129.0426,
        150.0 + 97.0527 + 129.0426 + 97.0527,
    ];
    let tagger = Tagger::new(options(2, 3, Tolerance::Absolute(0.02))).unwrap();
    let tags = tagger.get_tags(&trace).unwrap();
    for expected in ["PE", "EP", "PEP"] {
        assert!(tags.iter().any(|tag| tag == expected));
    }
    for (error, expected) in [(0.015, true), (0.025, false)] {
        let trace = [
            200.0,
            200.0 + 97.0527 + error,
            200.0 + 97.0527 + 129.0426 + error,
        ];
        let tagger = Tagger::new(options(2, 2, Tolerance::Absolute(0.02))).unwrap();
        assert_eq!(
            tagger
                .get_tags(&trace)
                .unwrap()
                .iter()
                .any(|tag| tag == "PE"),
            expected
        );
    }
    let compatibility = [200.0, 297.0527, 426.0953];
    let tagger = Tagger::new(options(2, 2, Tolerance::Ppm(20.0))).unwrap();
    assert_eq!(tagger.get_tags(&compatibility).unwrap(), ["PE"]);
    let positive = Tagger::new(options(2, 3, Tolerance::Absolute(0.02))).unwrap();
    let negative = Tagger::new(options(2, 3, Tolerance::Absolute(-0.02))).unwrap();
    assert_eq!(
        positive.get_tags(&trace).unwrap(),
        negative.get_tags(&trace).unwrap()
    );
}

// Source free-residue water, directly from the literal source element masses.
const WATER: f64 = 2.0 * 1.007_825_031_9 + 15.994_915;
fn mass_record(name: &str, origin: char, internal_mass: f64) -> ResidueModification {
    ResidueModification::from_record(ModificationRecord {
        name: name.into(),
        origin: Some(origin),
        diff_mono_mass: 1.0,
        mono_mass: WATER + internal_mass,
        ..Default::default()
    })
    .unwrap()
}
fn custom_tagger(
    records: Vec<ResidueModification>,
    names: &[&str],
    tolerance: f64,
    length: usize,
) -> Tagger {
    let registry = ModificationsDB::from_records(records).unwrap();
    let mut config = options(length, length, Tolerance::Absolute(tolerance));
    config.fixed_mods = names.iter().map(|name| (*name).to_owned()).collect();
    Tagger::with_registry(config, &registry).unwrap()
}

#[test]
fn strict_lower_boundary_nearest_tie_and_outside_successor_are_source_specific() {
    let build = |tolerance| {
        custom_tagger(
            vec![
                mass_record("Low", 'G', 32.0),
                mass_record("Mid", 'A', 34.0),
                mass_record("High", 'C', 36.0),
            ],
            &["Low (G)", "Mid (A)", "High (C)"],
            tolerance,
            1,
        )
    };
    // The first lower_bound candidate is 32 at the excluded lower edge.
    // Source returns immediately, even though the exact 34 match exists later.
    assert!(build(2.0).get_tags(&[0.0, 34.0]).unwrap().is_empty());
    assert_eq!(build(2.5).get_tags(&[0.0, 34.0]).unwrap(), ["A"]);
    assert_eq!(build(1.5).get_tags(&[0.0, 33.0]).unwrap(), ["G"]);
    // The next map entry is visited after increment but lies outside this
    // window; it cannot replace the closer valid lower mass.
    assert_eq!(build(0.5).get_tags(&[0.0, 32.25]).unwrap(), ["G"]);
    assert!(build(0.0).get_tags(&[0.0, 32.0]).unwrap().is_empty());
}

#[test]
fn fixed_mass_collisions_follow_input_order_and_l_expands_even_when_modified() {
    let build = |names: &[&str]| {
        custom_tagger(
            vec![
                mass_record("Low", 'G', 32.0),
                mass_record("Other", 'A', 32.0),
            ],
            names,
            0.01,
            1,
        )
    };
    assert_eq!(
        build(&["Low (G)", "Other (A)"])
            .get_tags(&[0.0, 32.0])
            .unwrap(),
        ["A"]
    );
    assert_eq!(
        build(&["Other (A)", "Low (G)"])
            .get_tags(&[0.0, 32.0])
            .unwrap(),
        ["G"]
    );
    let modified_l = custom_tagger(
        vec![mass_record("Shift", 'L', 64.0)],
        &["Shift (L)"],
        0.01,
        2,
    );
    assert_eq!(
        modified_l.get_tags(&[0.0, 64.0, 128.0]).unwrap(),
        ["II", "IL", "LI", "LL"]
    );
    let modified_i = custom_tagger(
        vec![mass_record("Shift", 'I', 64.0)],
        &["Shift (I)"],
        0.01,
        2,
    );
    assert_eq!(modified_i.get_tags(&[0.0, 64.0, 128.0]).unwrap(), ["II"]);
}

#[test]
fn source_zero_free_mass_for_bzx_remains_private_negative_gap_geometry() {
    for origin in ['B', 'Z', 'X'] {
        let record = ResidueModification::from_record(ModificationRecord {
            name: "Unchanged".into(),
            origin: Some(origin),
            ..Default::default()
        })
        .unwrap();
        let full_id = record.full_id().to_owned();
        let tagger = custom_tagger(vec![record], &[&full_id], 0.001, 2);
        let tags = tagger.get_tags(&[0.0, -WATER, -2.0 * WATER]).unwrap();
        assert_eq!(tags, [origin.to_string().repeat(2)]);
    }
}

#[test]
fn source_index_pruning_and_distinct_append_noops_are_retained() {
    let tagger = Tagger::new(options(2, 2, Tolerance::Absolute(0.02))).unwrap();
    assert!(
        tagger
            .get_tags(&[0.0, 400.0, 97.0527, 226.0953])
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        tagger.get_tags(&[0.0, 97.0527, 226.0953, 400.0]).unwrap(),
        ["PE"]
    );
    let mut tags = vec!["z".into(), "A".into(), "z".into()];
    let long = Tagger::new(options(3, 3, Tolerance::Ppm(10.0))).unwrap();
    let before = tags.clone();
    long.append_tags(&[0.0, 1.0], &mut tags).unwrap();
    assert_eq!(
        tags, before,
        "source min length > N skips output normalization"
    );
    long.append_tags(&[0.0, 1.0, 2.0], &mut tags).unwrap();
    assert_eq!(
        tags,
        ["A", "z"],
        "source min length == N still normalizes output"
    );
    let mut config = options(2, 2, Tolerance::Absolute(100.0));
    config.min_charge = 0;
    config.max_charge = 0;
    assert_eq!(
        Tagger::new(config)
            .unwrap()
            .get_tags(&[5.0, 5.0, -100.0])
            .unwrap(),
        ["GG"]
    );
    let zero_max = Tagger::new(options(0, 0, Tolerance::Absolute(1.0))).unwrap();
    let mut tags = vec!["z".into(), "a".into(), "z".into()];
    zero_max.append_tags(&[0.0, 32.0], &mut tags).unwrap();
    assert_eq!(tags, ["a", "z"]);
    let mut config = options(1, 1, Tolerance::Absolute(1.0));
    config.min_charge = 2;
    config.max_charge = 1;
    let mut inverted = Tagger::new(config).unwrap();
    let mut tags = vec!["b".into(), "a".into(), "b".into()];
    inverted.append_tags(&[0.0, 97.0527], &mut tags).unwrap();
    assert_eq!(tags, ["a", "b"]);
    inverted.set_max_charge(2).unwrap();
    assert_eq!(inverted.get_tags(&[0.0, 97.0527 / 2.0]).unwrap(), ["P"]);
}
