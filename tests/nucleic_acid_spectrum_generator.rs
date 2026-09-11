// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{
    EmpiricalFormula, NAFragmentType, NASequence, NucleicAcidSpectrumGenerator as Generator,
    PROTON_MASS_U, Ribonucleotide, RibonucleotideRecord,
};
use openms::kernel::{DataArray, MSSpectrum, Peak1D, SpectrumType};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

fn na(text: &str) -> NASequence {
    NASequence::parse(text).unwrap()
}
fn mass(text: &str) -> f64 {
    EmpiricalFormula::parse(text).unwrap().mono_mass()
}
fn off() -> Generator {
    Generator {
        add_b_ions: false,
        add_y_ions: false,
        ..Default::default()
    }
}
fn custom(records: &[(&str, &str, f64)]) -> NASequence {
    NASequence::from_records(
        records
            .iter()
            .map(|&(code, formula, mono_mass)| {
                Arc::new(
                    Ribonucleotide::from_record(RibonucleotideRecord {
                        code: code.into(),
                        formula: formula.parse().unwrap(),
                        mono_mass,
                        ..Default::default()
                    })
                    .unwrap(),
                )
            })
            .collect(),
    )
    .unwrap()
}
fn names(s: &MSSpectrum) -> &[String] {
    &s.string_data_arrays[0].data
}
fn charges(s: &MSSpectrum) -> &[i32] {
    &s.integer_data_arrays[0].data
}
fn key_set(values: &[i32]) -> BTreeSet<i32> {
    values.iter().copied().collect()
}

#[test]
fn defaults_and_declared_mass_provenance_leave_source_spectrum_settings() {
    let generator = Generator::default();
    let sequence = na("ACG");
    let spectrum = generator.generate(&sequence, -1, -1).unwrap();
    assert_eq!(spectrum.len(), 3); // b2, y1, y2; first-prefix option does not hide y1
    assert!(spectrum.string_data_arrays.is_empty());
    assert!(spectrum.integer_data_arrays.is_empty());
    assert_eq!(spectrum.ms_level, 1);
    assert_eq!(spectrum.spectrum_type, SpectrumType::Unknown);
    assert!(spectrum.precursors.is_empty());
    let a = sequence.residues()[0].mono_mass();
    let c = sequence.residues()[1].mono_mass();
    let g = sequence.residues()[2].mono_mass();
    let k = mass("H-1PO2");
    let mut expected = vec![
        (a + c + k) / -1.0 + PROTON_MASS_U,
        g / -1.0 + PROTON_MASS_U,
        (g + c + k) / -1.0 + PROTON_MASS_U,
    ];
    for mz in &mut expected {
        *mz = mz.abs();
    }
    expected.sort_by(f64::total_cmp);
    assert_eq!(
        spectrum.peaks.iter().map(|p| p.mz).collect::<Vec<_>>(),
        expected
    );
    assert_eq!(
        generator.clone().generate(&sequence, -1, -1).unwrap(),
        spectrum
    );
}

#[test]
fn single_charge_swaps_cutoff_and_source_precursor_selection() {
    let sequence = na("ACG");
    let generator = Generator {
        add_metainfo: true,
        add_precursor_peaks: true,
        ..Default::default()
    };
    assert_eq!(
        generator.generate(&sequence, -2, -1).unwrap(),
        generator.generate(&sequence, -1, -2).unwrap()
    );
    let clipped = generator.generate(&sequence, -1, -3).unwrap();
    assert!(!names(&clipped).iter().any(|n| n == "M"));
    assert_eq!(
        charges(&clipped).iter().copied().collect::<BTreeSet<_>>(),
        key_set(&[-1, -2])
    );
    let exact = generator.generate(&sequence, -1, -2).unwrap();
    assert_eq!(
        names(&exact).iter().filter(|n| n.as_str() == "M").count(),
        1
    );
    assert_eq!(
        charges(&exact)[names(&exact).iter().position(|n| n == "M").unwrap()],
        -2
    );
    assert!(generator.generate(&sequence, -1, 1).is_err());
    assert!(generator.generate(&sequence, i32::MIN, -1).is_err());
}

#[test]
fn zero_charge_errors_only_when_a_peak_is_divided_and_unused_options_stay_unused() {
    let sequence = na("AC");
    let empty = Generator {
        add_metainfo: true,
        a_intensity: f64::NAN,
        ..off()
    };
    let result = empty.generate(&sequence, 0, 0).unwrap();
    assert!(result.is_empty());
    assert!(names(&result).is_empty());
    assert!(Generator::default().generate(&sequence, 0, 1).is_err());
    let precursor = Generator {
        add_metainfo: true,
        add_precursor_peaks: true,
        ..off()
    };
    // Zero's precursor is excluded; max=-1 with min=0 still selects POSITIVE mode.
    let result = precursor.generate(&sequence, 0, -1).unwrap();
    assert_eq!(names(&result), &["M"]);
    assert_eq!(charges(&result), &[1]);
    assert!(precursor.generate(&sequence, 0, 0).is_err());
    assert!(
        precursor
            .generate(&NASequence::new(), 0, 0)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn multiple_source_skipped_key_rules_and_mixed_mode_are_preserved() {
    let sequence = na("ACG");
    let set = key_set(&[1, 2]);
    assert!(
        Generator::default()
            .generate_multiple(&sequence, &set, 3)
            .unwrap()
            .is_empty()
    );
    let generator = Generator {
        add_metainfo: true,
        ..Default::default()
    };
    let skipped = generator.generate_multiple(&sequence, &set, 3).unwrap();
    assert_eq!(skipped.len(), 2);
    assert!(
        skipped
            .values()
            .all(|s| s.is_empty() && s.integer_data_arrays.len() == 1)
    );
    let mixed = generator
        .generate_multiple(&sequence, &key_set(&[-2, 1]), 1)
        .unwrap();
    assert!(mixed[&1].is_empty());
    assert_eq!(
        charges(&mixed[&-2])
            .iter()
            .copied()
            .collect::<BTreeSet<_>>(),
        key_set(&[-1, -2])
    );
    let no_meta = Generator::default()
        .generate_multiple(&sequence, &key_set(&[-2, 1]), 1)
        .unwrap();
    assert_eq!(no_meta.keys().copied().collect::<Vec<_>>(), vec![-2]);
}

#[test]
fn multiple_final_precursor_uses_next_counter_and_positive_branch_has_no_abs() {
    let sequence = custom(&[("Q", "C-100", 100.0)]);
    let generator = Generator {
        add_metainfo: true,
        add_precursor_peaks: true,
        ..off()
    };
    assert!(generator.generate(&sequence, 1, 1).unwrap().is_empty());
    let positive = generator
        .generate_multiple(&sequence, &key_set(&[1, 3]), 1)
        .unwrap();
    let raw = sequence.mono_mass(NAFragmentType::Full, 0).unwrap();
    assert_eq!(charges(&positive[&1]), &[2]);
    assert_eq!(charges(&positive[&3]), &[4]);
    assert_eq!(positive[&1].peaks[0].mz, raw / 2.0 + PROTON_MASS_U);
    assert!(positive[&1].peaks[0].mz < 0.0);
    assert_eq!(positive[&3].len(), 1); // first final-only precursor was not inherited
    let negative = generator
        .generate_multiple(&sequence, &key_set(&[-1, -3]), 1)
        .unwrap();
    assert_eq!(charges(&negative[&-1]), &[-2]);
    assert_eq!(
        negative[&-1].peaks[0].mz,
        (raw / -2.0 + PROTON_MASS_U).abs()
    );
    let all = Generator {
        add_all_precursor_charges: true,
        ..generator
    };
    let multiple = all.generate_multiple(&sequence, &key_set(&[3]), 1).unwrap();
    assert_eq!(multiple[&3].len(), 3); // no charge<sequence-length guard in multi
    assert_eq!(
        charges(&multiple[&3])
            .iter()
            .copied()
            .collect::<BTreeSet<_>>(),
        key_set(&[1, 2, 3])
    );
}

#[test]
fn ambiguity_halves_f64_intensity_before_f32_and_keeps_both_labels() {
    let sequence = custom(&[("Q?", "C", 100.0), ("A", "C", 100.0)]);
    let generator = Generator {
        add_first_prefix_ion: true,
        add_a_minus_b_ions: true,
        add_metainfo: true,
        a_minus_b_intensity: 6.0e38,
        ..off()
    };
    let spectrum = generator.generate(&sequence, 1, 1).unwrap();
    assert_eq!(spectrum.len(), 2);
    assert_eq!(names(&spectrum), &["a1-B", "a1-B"]);
    assert!(
        spectrum
            .peaks
            .iter()
            .all(|p| p.intensity == (6.0e38_f64 * 0.5) as f32)
    );
    assert!(
        generator
            .generate(&custom(&[("Q?*", "C", 100.0), ("A", "C", 100.0)]), 1, 1)
            .is_err()
    );
    let negative = Generator {
        a_minus_b_intensity: -2.0,
        ..generator
    };
    assert!(
        negative
            .generate(&sequence, 1, 1)
            .unwrap()
            .peaks
            .iter()
            .all(|p| p.intensity == -1.0)
    );
}

#[test]
fn append_preserves_unrelated_metadata_and_first_array_names_with_aligned_padding() {
    let mut spectrum = MSSpectrum {
        rt: f64::NAN,
        ms_level: 0,
        native_id: "kept native id".into(),
        ..Default::default()
    };
    spectrum
        .metadata
        .insert("retained".into(), "retained string allocation".into());
    let metadata_pointer = spectrum.metadata["retained"].as_str().unwrap().as_ptr();
    spectrum.peaks = vec![Peak1D::new(900.0, 7.0), Peak1D::new(10.0, -3.0)];
    spectrum
        .integer_data_arrays
        .push(DataArray::new("original integer name", Vec::new()));
    spectrum.string_data_arrays.push(DataArray::new(
        "original label name",
        vec!["old-high".into(), "old-low".into()],
    ));
    let name_pointer = spectrum.string_data_arrays[0].name.as_ptr();
    spectrum
        .float_data_arrays
        .push(DataArray::new("placeholder", Vec::new()));
    let generator = Generator {
        add_metainfo: true,
        ..Default::default()
    };
    generator
        .append_to(&mut spectrum, &na("ACG"), 1, 1)
        .unwrap();
    assert!(spectrum.rt.is_nan());
    assert_eq!(spectrum.ms_level, 0);
    assert_eq!(
        spectrum.metadata["retained"].as_str().unwrap().as_ptr(),
        metadata_pointer
    );
    assert_eq!(spectrum.string_data_arrays[0].name.as_ptr(), name_pointer);
    assert_eq!(
        spectrum.integer_data_arrays[0].name,
        "original integer name"
    );
    assert_eq!(spectrum.string_data_arrays[0].name, "original label name");
    assert_eq!(names(&spectrum)[0], "old-low");
    assert_eq!(charges(&spectrum)[0], 0);
    assert!(spectrum.peaks.windows(2).all(|w| w[0].mz <= w[1].mz));
    assert_eq!(spectrum.peaks.len(), names(&spectrum).len());
    assert_eq!(spectrum.peaks.len(), charges(&spectrum).len());
    assert!(spectrum.float_data_arrays[0].data.is_empty());
}

#[test]
fn no_additions_still_sort_old_peaks_and_retain_empty_annotation_placeholders() {
    let mut spectrum = MSSpectrum {
        peaks: vec![Peak1D::new(2.0, 2.0), Peak1D::new(1.0, 1.0)],
        ..Default::default()
    };
    spectrum
        .float_data_arrays
        .push(DataArray::new("paired", vec![20.0, 10.0]));
    let generator = Generator {
        add_metainfo: true,
        ..off()
    };
    generator
        .append_to(&mut spectrum, &NASequence::new(), 1, 1)
        .unwrap();
    assert_eq!(spectrum.peaks[0].mz, 1.0);
    assert_eq!(spectrum.float_data_arrays[0].data, vec![10.0, 20.0]);
    assert!(names(&spectrum).is_empty());
    assert!(charges(&spectrum).is_empty());
}

#[test]
fn unsupported_arrays_and_late_mass_failure_leave_destination_unchanged() {
    let mut spectrum = MSSpectrum {
        peaks: vec![Peak1D::new(100.0, 2.0)],
        ..Default::default()
    };
    spectrum
        .float_data_arrays
        .push(DataArray::new("no new values", vec![3.0]));
    let old = spectrum.clone();
    let pointer = spectrum.peaks.as_ptr();
    assert!(
        Generator::default()
            .append_to(&mut spectrum, &na("ACG"), 1, 1)
            .is_err()
    );
    assert_eq!(spectrum, old);
    assert_eq!(spectrum.peaks.as_ptr(), pointer);
    spectrum.float_data_arrays.clear();
    let bad = custom(&[("Q", "C", 1.0), ("R", "C", f64::MAX), ("S", "C", f64::MAX)]);
    let old = spectrum.clone();
    let generator = Generator {
        add_first_prefix_ion: true,
        ..Default::default()
    };
    assert!(generator.append_to(&mut spectrum, &bad, 1, 1).is_err());
    assert_eq!(spectrum, old);
    assert_eq!(spectrum.peaks.as_ptr(), pointer);
}

#[test]
fn multiple_replacement_is_atomic_and_empty_charge_set_really_clears() {
    let mut output = BTreeMap::from([(
        7,
        MSSpectrum {
            native_id: "old".into(),
            ..Default::default()
        },
    )]);
    let old = output.clone();
    let final_only = Generator {
        add_precursor_peaks: true,
        ..off()
    };
    assert!(
        final_only
            .replace_multiple(&mut output, &NASequence::new(), &key_set(&[1]), 1)
            .is_err()
    );
    assert_eq!(output, old);
    assert!(
        final_only
            .replace_multiple(&mut output, &na("AC"), &key_set(&[i32::MAX]), 1)
            .is_err()
    );
    assert_eq!(output, old);
    final_only
        .replace_multiple(&mut output, &na("AC"), &BTreeSet::new(), 0)
        .unwrap();
    assert!(output.is_empty());
}

#[test]
fn multiple_cumulative_peak_count_and_empty_copy_work_are_bounded() {
    let all = Generator {
        add_precursor_peaks: true,
        add_all_precursor_charges: true,
        ..off()
    };
    let keys = (1..=500).collect();
    let error = all.generate_multiple(&na("AC"), &keys, 1).unwrap_err();
    assert!(error.to_string().contains("peak limit"), "{error}");
    let empty = off()
        .generate_multiple(&na("AC"), &key_set(&[i32::MAX - 1]), 1)
        .unwrap_err();
    assert!(empty.to_string().contains("work limit"), "{empty}");
}

#[test]
fn configured_multi_matches_single_without_precursors_with_stable_duplicate_rows() {
    let generator = Generator {
        add_metainfo: true,
        add_first_prefix_ion: true,
        add_a_ions: true,
        add_c_ions: true,
        add_d_ions: true,
        add_w_ions: true,
        add_x_ions: true,
        add_z_ions: true,
        add_a_minus_b_ions: true,
        ..Default::default()
    };
    let sequence = na("[m1A]UC[C*]AC[A*]Gp");
    let targets = key_set(&[-1, -3, -5]);
    let multiple = generator
        .generate_multiple(&sequence, &targets, -1)
        .unwrap();
    for charge in targets {
        assert_eq!(
            multiple[&charge],
            generator.generate(&sequence, -1, charge).unwrap()
        );
    }
}
