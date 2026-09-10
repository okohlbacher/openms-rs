// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::tagger::{MAX_TAGGER_PEAKS, MAX_TAGGER_TAGS};
use openms::chemistry::{
    EmpiricalFormula, ModificationRecord, ModificationsDB, ResidueModification, Tagger,
    TaggerOptions,
};
use openms::comparison::Tolerance;
use openms::kernel::DataArray;
use openms::{MSSpectrum, Peak1D};

fn pep() -> [f64; 4] {
    [150.0, 247.0527, 376.0953, 473.1480]
}
fn options() -> TaggerOptions {
    let mut options = TaggerOptions::new(2, Tolerance::Absolute(0.02));
    options.max_tag_length = 3;
    options
}

#[test]
fn defaults_negative_tolerance_and_checked_setter_keep_source_option_behavior() {
    let config = TaggerOptions::new(2, Tolerance::Ppm(-20.0));
    assert_eq!(
        (config.max_tag_length, config.min_charge, config.max_charge),
        (65_535, 1, 1)
    );
    assert!(config.fixed_mods.is_empty() && config.variable_mods.is_empty());
    let tagger = Tagger::new(config).unwrap();
    assert_eq!(tagger.options().tolerance, Tolerance::Ppm(20.0));
    let mut tagger = Tagger::new(options()).unwrap();
    let old = tagger.clone();
    assert!(tagger.set_max_charge(usize::MAX).is_err());
    assert_eq!(tagger, old);
    tagger.set_max_charge(0).unwrap();
    assert!(tagger.get_tags(&pep()).unwrap().is_empty());
    tagger.set_max_charge(1).unwrap();
    assert_eq!(tagger.get_tags(&pep()).unwrap(), ["EP", "PE", "PEP"]);
}

#[test]
fn atomic_append_moves_old_strings_and_sorts_deduplicates_both_sources() {
    let tagger = Tagger::new(options()).unwrap();
    let mut tags = vec![
        "é".to_owned(),
        "PE".to_owned(),
        String::new(),
        "PE".to_owned(),
        "A".to_owned(),
    ];
    let unicode_pointer = tags[0].as_ptr();
    let peptide_pointer = tags[1].as_ptr();
    tagger.append_tags(&pep(), &mut tags).unwrap();
    assert_eq!(tags, ["", "A", "EP", "PE", "PEP", "é"]);
    assert_eq!(tags[3].as_ptr(), peptide_pointer);
    assert_eq!(tags[5].as_ptr(), unicode_pointer);
}

#[test]
fn early_minimum_shortcut_and_equal_length_normalization_are_distinct() {
    let tagger = Tagger::new(options()).unwrap();
    let mut tags = vec!["z".into(), "a".into(), "a".into()];
    let saved = tags.clone();
    let pointer = tags.as_ptr();
    tagger.append_tags(&[f64::NAN], &mut tags).unwrap();
    assert_eq!(tags, saved);
    assert_eq!(tags.as_ptr(), pointer);
    tagger.append_tags(&[1.0, 2.0], &mut tags).unwrap();
    assert_eq!(tags, ["a", "z"]);
}

#[test]
fn zero_and_inverted_lengths_and_charges_are_not_normalized_away() {
    let mut config = TaggerOptions::new(0, Tolerance::Absolute(0.02));
    config.max_tag_length = 0;
    let tagger = Tagger::new(config).unwrap();
    assert!(tagger.get_tags(&pep()).unwrap().is_empty());
    let mut old = vec!["z".to_owned(), "a".to_owned(), "a".to_owned()];
    tagger.append_tags(&[], &mut old).unwrap();
    assert_eq!(old, ["a", "z"]);
    let mut config = options();
    config.min_tag_length = 3;
    config.max_tag_length = 1;
    assert!(
        Tagger::new(config)
            .unwrap()
            .get_tags(&pep())
            .unwrap()
            .is_empty()
    );
    let mut config = TaggerOptions::new(1, Tolerance::Absolute(100.0));
    config.min_charge = 0;
    config.max_charge = 0;
    // Zero charge still performs source mass matching against zero. A broad
    // absolute window chooses the lightest natural residue, glycine.
    assert_eq!(
        Tagger::new(config).unwrap().get_tags(&[0.0, 1.0]).unwrap(),
        ["G"]
    );
    let mut config = options();
    config.min_charge = usize::MAX;
    config.max_charge = 1;
    assert!(
        Tagger::new(config)
            .unwrap()
            .get_tags(&pep())
            .unwrap()
            .is_empty()
    );
}

#[test]
fn spectrum_entrypoints_ignore_unused_intensity_metadata_and_unaligned_arrays() {
    let tagger = Tagger::new(options()).unwrap();
    let mut spectrum = MSSpectrum::from_peaks(
        pep()
            .into_iter()
            .map(|mz| Peak1D::new(mz, f32::NAN))
            .collect(),
    );
    spectrum.rt = f64::NAN;
    spectrum.ms_level = 0;
    spectrum
        .float_data_arrays
        .push(DataArray::new("unused", vec![1.0]));
    assert_eq!(
        tagger.get_spectrum_tags(&spectrum).unwrap(),
        ["EP", "PE", "PEP"]
    );
    let mut tags = vec!["PEP".to_owned()];
    tagger.append_spectrum_tags(&spectrum, &mut tags).unwrap();
    assert_eq!(tags, ["EP", "PE", "PEP"]);
    assert!(spectrum.peaks.iter().all(|peak| peak.intensity.is_nan()));
    assert_eq!(spectrum.float_data_arrays[0].data.len(), 1);
}

#[test]
fn late_overflow_and_invalid_coordinates_leave_existing_output_unchanged() {
    let tagger = Tagger::new(TaggerOptions::new(1, Tolerance::Absolute(0.02))).unwrap();
    for mzs in [
        vec![f64::NAN, 100.0],
        vec![0.0, f64::INFINITY],
        vec![0.0, 97.0527, f64::MAX, -f64::MAX],
    ] {
        let mut tags = vec!["keep".to_owned(), "keep".to_owned()];
        let pointers: Vec<_> = tags.iter().map(|tag| tag.as_ptr()).collect();
        assert!(tagger.append_tags(&mzs, &mut tags).is_err());
        assert_eq!(tags, ["keep", "keep"]);
        assert_eq!(
            tags.iter().map(|tag| tag.as_ptr()).collect::<Vec<_>>(),
            pointers
        );
    }
    for tolerance in [Tolerance::Ppm(f64::NAN), Tolerance::Absolute(f64::INFINITY)] {
        assert!(Tagger::new(TaggerOptions::new(1, tolerance)).is_err());
    }
}

#[test]
fn caller_registry_lifetime_and_later_changes_do_not_change_resolved_tags() {
    let mut registry = ModificationsDB::from_records(vec![
        ResidueModification::from_record(ModificationRecord {
            name: "MeasuredA".into(),
            origin: Some('A'),
            diff_mono_mass: 1.0,
            mono_mass: 300.0,
            ..Default::default()
        })
        .unwrap(),
    ])
    .unwrap();
    let mut config = TaggerOptions::new(1, Tolerance::Absolute(0.001));
    config.fixed_mods.push("MeasuredA".into());
    let tagger = Tagger::with_registry(config, &registry).unwrap();
    registry
        .extend_records(vec![
            ResidueModification::from_record(ModificationRecord {
                name: "DifferentA".into(),
                origin: Some('A'),
                diff_mono_mass: 1.0,
                mono_mass: 400.0,
                ..Default::default()
            })
            .unwrap(),
        ])
        .unwrap();
    drop(registry);
    let mass = 300.0 - EmpiricalFormula::parse("H2O").unwrap().mono_mass();
    assert_eq!(tagger.get_tags(&[0.0, mass]).unwrap(), ["A"]);
    assert!(tagger.get_tags(&[0.0, 71.0371]).unwrap().is_empty());
    assert_eq!(tagger.clone().get_tags(&[0.0, mass]).unwrap(), ["A"]);
}

#[test]
fn actual_peak_and_existing_tag_caps_fail_before_large_search_or_copy() {
    let tagger = Tagger::new(TaggerOptions::new(0, Tolerance::Absolute(0.01))).unwrap();
    assert!(
        tagger
            .get_tags(&vec![0.0; MAX_TAGGER_PEAKS + 1])
            .unwrap_err()
            .to_string()
            .contains("peak count limit")
    );
    let mut tags = vec![String::new(); MAX_TAGGER_TAGS + 1];
    let pointer = tags.as_ptr();
    assert!(
        tagger
            .append_tags(&[], &mut tags)
            .unwrap_err()
            .to_string()
            .contains("tag count limit")
    );
    assert_eq!(tags.len(), MAX_TAGGER_TAGS + 1);
    assert_eq!(tags.as_ptr(), pointer);
}
