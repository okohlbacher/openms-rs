use openms::{
    MSSpectrum, Peak1D,
    chemistry::{
        AASequence, EmpiricalFormula, ModificationRecord, ModificationsDB, ProteinProteinCrossLink,
        ResidueModification, TermSpecificity, TheoreticalSpectrumGeneratorXLMS as Generator,
        XLMSOptions,
    },
    kernel::DataArray,
    metadata::DataProcessing,
};
use std::{collections::BTreeSet, sync::Arc};
fn peptide() -> AASequence {
    AASequence::parse("IFSQVGK").unwrap()
}
fn pair() -> ProteinProteinCrossLink {
    let mut p = ProteinProteinCrossLink::new(150.).unwrap();
    p.alpha = Some(Arc::new(peptide()));
    p.beta = Some(Arc::new(AASequence::parse("TESTPEP").unwrap()));
    p.cross_link_position = (3, 4);
    p
}
fn only(g: &mut Generator, types: &str) {
    g.options.add_a_ions = types.contains('a');
    g.options.add_b_ions = types.contains('b');
    g.options.add_c_ions = types.contains('c');
    g.options.add_x_ions = types.contains('x');
    g.options.add_y_ions = types.contains('y');
    g.options.add_z_ions = types.contains('z');
}
fn spectrum(g: &Generator, case: usize, alpha: bool, max: i32) -> MSSpectrum {
    let mut s = MSSpectrum::default();
    match case {
        0 => g.get_linear_ion_spectrum(&mut s, &peptide(), 3, alpha, max, 0),
        1 => g.get_xlink_ion_spectrum(&mut s, &peptide(), 3, 2000., alpha, 2, max, 0),
        _ => g.get_crosslink_ion_spectrum(&mut s, &pair(), alpha, 2, max),
    }
    .unwrap();
    s
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a} != {b}");
}
#[test]
fn all_three_source_literal_mass_arrays() {
    let mut g = Generator::default();
    for case in 0..3 {
        only(&mut g, if case == 0 { "aby" } else { "by" });
        let actual = spectrum(&g, case, true, if case == 0 { 2 } else { 3 });
        let expected: Vec<f64> = include_str!("data/xlms_source_masses.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
            .filter_map(|l| {
                let c: Vec<_> = l.split('\t').collect();
                (c[0].parse::<usize>().unwrap() == case).then(|| c[2].parse().unwrap())
            })
            .collect();
        assert_eq!(actual.len(), expected.len());
        for (p, x) in actual.peaks.iter().zip(expected) {
            assert!((p.mz - x).abs() < 0.001, "case{case}:{} vs{x}", p.mz);
        }
    }
}
#[test]
fn source_counts_all_series_losses_and_isotope_saturation() {
    let mut g = Generator::default();
    assert_eq!(spectrum(&g, 0, true, 3).len(), 27);
    only(&mut g, "abcxyz");
    assert_eq!(spectrum(&g, 0, true, 3).len(), 54);
    for case in 1..3 {
        only(&mut g, "by");
        g.options.add_losses = false;
        assert_eq!(spectrum(&g, case, true, 3).len(), 17);
        g.options.add_losses = true;
        assert_eq!(
            spectrum(&g, case, true, 3).len(),
            if case == 1 { 39 } else { 41 }
        );
        g.options.add_losses = false;
        assert_eq!(spectrum(&g, case, true, 4).len(), 24);
        only(&mut g, "abcxyz");
        assert_eq!(spectrum(&g, case, true, 4).len(), 60);
    }
    only(&mut g, "by");
    g.options.add_precursor_peaks = false;
    g.options.add_k_linked_ions = false;
    g.options.add_isotopes = true;
    for max_isotope in [1, 2, 3] {
        g.options.max_isotope = max_isotope;
        g.options.add_losses = true;
        assert_eq!(
            spectrum(&g, 0, true, 3).len(),
            if max_isotope < 2 { 30 } else { 48 }
        );
        g.options.add_losses = false;
        for case in 1..3 {
            assert_eq!(
                spectrum(&g, case, true, 5).len(),
                if max_isotope < 2 { 24 } else { 48 }
            );
        }
    }
}
#[test]
fn all_six_source_allowed_annotation_sets_and_charge_counts() {
    let mut g = Generator::default();
    only(&mut g, "bx");
    g.options.add_losses = true;
    for case in 0..3 {
        for alpha in [true, false] {
            let max = if case == 0 {
                3
            } else if alpha {
                5
            } else {
                4
            };
            let s = spectrum(&g, case, alpha, max);
            if alpha {
                assert_eq!(s.len(), [30, 75, 79][case]);
            }
            let allowed: BTreeSet<_> = include_str!("data/xlms_source_names.tsv")
                .lines()
                .filter(|l| !l.starts_with('#'))
                .filter_map(|l| {
                    let c: Vec<_> = l.split('\t').collect();
                    (c[0].parse::<usize>().unwrap() == case
                        && c[1] == if alpha { "1" } else { "0" })
                    .then_some(c[2])
                })
                .collect();
            assert_eq!(s.string_data_arrays[0].name, "IonNames");
            assert_eq!(s.integer_data_arrays[0].name, "charge");
            for name in &s.string_data_arrays[0].data {
                assert!(
                    allowed.contains(name.as_str()),
                    "case{case} {alpha}: {name}"
                );
            }
            if case == 0 {
                for z in 1..=3 {
                    assert_eq!(
                        s.integer_data_arrays[0]
                            .data
                            .iter()
                            .filter(|&&v| v == z)
                            .count(),
                        10
                    );
                }
            } else if !alpha {
                let expected = if case == 1 {
                    [0, 18, 18, 21, 0]
                } else {
                    [0, 19, 19, 22, 0]
                };
                for (offset, count) in expected.into_iter().enumerate() {
                    assert_eq!(
                        s.integer_data_arrays[0]
                            .data
                            .iter()
                            .filter(|&&z| z == offset as i32 + 1)
                            .count(),
                        count
                    );
                }
            }
        }
    }
}

#[test]
fn source_terminal_correction_and_numeric_tags_use_retained_mass() {
    let mut g = Generator::default();
    only(&mut g, "y");
    let p = AASequence::parse("AA.[+5.123456789]").unwrap();
    let mut s = MSSpectrum::default();
    g.get_linear_ion_spectrum(&mut s, &p, 0, true, 1, 0)
        .unwrap();
    let alanine = EmpiricalFormula::parse("C3H7NO2").unwrap().mono_mass();
    near(s.peaks[0].mz, 1.007276466771 + 5.123456789 + alanine);
    only(&mut g, "b");
    let p = AASequence::parse("A[+5.123456789]A").unwrap();
    s = MSSpectrum::default();
    g.get_linear_ion_spectrum(&mut s, &p, 1, true, 1, 0)
        .unwrap();
    near(
        s.peaks[0].mz,
        1.007276466771
            + (alanine + 5.123456789 - EmpiricalFormula::parse("H2O").unwrap().mono_mass()),
    );
}

#[test]
fn precursor_only_pair_uses_linker_once_and_ignores_record_annotations() {
    let mut g = Generator::default();
    only(&mut g, "");
    g.options.add_k_linked_ions = false;
    let mut p = pair();
    let mut a = MSSpectrum::default();
    g.get_crosslink_ion_spectrum(&mut a, &p, true, 1, 1)
        .unwrap();
    p.cross_linker_name = "annotation only".into();
    p.term_spec_alpha = TermSpecificity::NTerm;
    p.term_spec_beta = TermSpecificity::CTerm;
    p.precursor_correction = 99;
    let mut b = MSSpectrum::default();
    g.get_crosslink_ion_spectrum(&mut b, &p, false, 1, 1)
        .unwrap();
    assert_eq!(a, b);
    p.set_cross_linker_mass(151.).unwrap();
    let mut c = MSSpectrum::default();
    g.get_crosslink_ion_spectrum(&mut c, &p, true, 1, 1)
        .unwrap();
    for (x, y) in a.peaks.iter().zip(c.peaks) {
        near(y.mz - x.mz, 1.);
    }
}

#[test]
fn annotations_can_be_disabled_independently_and_negative_intensity_is_retained() {
    let mut g = Generator::default();
    only(&mut g, "b");
    g.options.b_intensity = -0.75;
    g.options.add_metainfo = false;
    let mut s = MSSpectrum::default();
    g.get_linear_ion_spectrum(&mut s, &peptide(), 3, true, 1, 0)
        .unwrap();
    assert!(s.string_data_arrays.is_empty());
    assert_eq!(s.integer_data_arrays.len(), 1);
    assert!(s.peaks.iter().all(|p| p.intensity == -0.75));
    g.options.add_metainfo = true;
    g.options.add_charges = false;
    s = MSSpectrum::default();
    g.get_linear_ion_spectrum(&mut s, &peptide(), 3, true, 1, 0)
        .unwrap();
    assert!(s.integer_data_arrays.is_empty());
    assert_eq!(s.string_data_arrays.len(), 1);
}

#[test]
fn dormant_pair_does_not_force_full_x_mass_and_length_endpoint_is_finite() {
    let mut g = Generator::default();
    only(&mut g, "");
    g.options.add_k_linked_ions = false;
    g.options.add_precursor_peaks = false;
    let mut p = ProteinProteinCrossLink::default();
    p.alpha = Some(Arc::new(AASequence::parse("X").unwrap()));
    p.cross_link_position = (-1, -1);
    let mut s = MSSpectrum::default();
    g.get_crosslink_ion_spectrum(&mut s, &p, true, 1, 1)
        .unwrap();
    assert!(s.is_empty());
    only(&mut g, "b");
    s = MSSpectrum::default();
    g.get_linear_ion_spectrum(&mut s, &AASequence::parse("AA").unwrap(), 2, true, 1, 0)
        .unwrap();
    assert_eq!(s.len(), 2);
    assert!(
        s.string_data_arrays[0]
            .data
            .contains(&"[alpha|ci$b2]".into())
    );
}
#[test]
fn source_smallest_loop_cases() {
    let mut g = Generator::default();
    only(&mut g, "by");
    g.options.add_precursor_peaks = false;
    g.options.add_k_linked_ions = false;
    for pos in [0, 1] {
        let p = AASequence::parse("HA").unwrap();
        let mut s = MSSpectrum::default();
        g.get_linear_ion_spectrum(&mut s, &p, pos, true, 1, 0)
            .unwrap();
        assert_eq!(s.len(), 1);
        s = MSSpectrum::default();
        g.get_xlink_ion_spectrum(&mut s, &p, pos, 2000., true, 1, 1, 0)
            .unwrap();
        assert_eq!(s.len(), 1);
    }
    let p = AASequence::parse("PEPTIDESAREWEIRD").unwrap();
    for (first, last, count) in [(1, 14, 2), (2, 14, 3), (2, 13, 4)] {
        let mut s = MSSpectrum::default();
        g.get_xlink_ion_spectrum(&mut s, &p, first, 2000., false, 1, 1, last)
            .unwrap();
        assert_eq!(s.len(), count);
        if last == 14 {
            s = MSSpectrum::default();
            g.get_linear_ion_spectrum(&mut s, &p, first, false, 1, last)
                .unwrap();
            assert_eq!(s.len(), count);
        }
    }
}
#[test]
fn configured_clone_and_two_inert_source_options() {
    let mut g = Generator::default();
    g.options.add_b_ions = false;
    g.options.a_intensity = 0.5;
    let before = spectrum(&g, 0, true, 2);
    assert_eq!(spectrum(&g.clone(), 0, true, 2), before);
    g.options.add_first_prefix_ion = false;
    g.options.add_abundant_immonium_ions = true;
    assert_eq!(spectrum(&g, 0, true, 2), before);
    let d = XLMSOptions::default();
    assert!(
        d.add_a_ions
            && d.add_b_ions
            && d.add_y_ions
            && d.add_metainfo
            && d.add_charges
            && d.add_precursor_peaks
            && d.add_k_linked_ions
    );
    assert!(!d.add_isotopes && !d.add_losses && !d.add_c_ions && !d.add_x_ions && !d.add_z_ions);
    assert_eq!(d.max_isotope, 2);
}
#[test]
fn cpp042_suffix_loss_source_value_is_distinct_from_chemical_expectation() {
    let mut g = Generator::default();
    only(&mut g, "y");
    g.options.add_losses = true;
    let mut s = MSSpectrum::default();
    g.get_linear_ion_spectrum(&mut s, &AASequence::parse("AS").unwrap(), 0, true, 2, 0)
        .unwrap();
    let row = |name: &str| {
        let i = (0..s.len())
            .find(|&i| {
                s.string_data_arrays[0].data[i] == name && s.integer_data_arrays[0].data[i] == 2
            })
            .unwrap();
        s.peaks[i].mz
    };
    let base = row("[alpha|ci$y1]");
    let loss = row("[alpha|ci$y1-H2O1]");
    near(base, 53.528573578421);
    near(loss, 17.7590042573105);
    let chemical = 44.523291046521;
    assert!((loss - chemical).abs() > 20.);
    near(
        loss,
        (base - EmpiricalFormula::parse("H2O").unwrap().mono_mass()) / 2.,
    );
}
#[test]
fn cpp043_precursor_companions_keep_undivided_source_mass() {
    let mut g = Generator::default();
    only(&mut g, "");
    g.options.add_k_linked_ions = false;
    g.options.add_isotopes = true;
    let mut s = MSSpectrum::default();
    g.get_xlink_ion_spectrum(&mut s, &AASequence::default(), 0, 1000., true, 2, 2, 0)
        .unwrap();
    assert_eq!(s.len(), 6);
    let xs: Vec<_> = (0..s.len())
        .filter(|&i| s.string_data_arrays[0].data[i] == "[M+H]")
        .map(|i| s.peaks[i].mz)
        .collect();
    near(xs[0], 501.007276466771);
    near(xs[1], 1002.516230352442);
    assert!((xs[1] - 501.508953885671).abs() > 500.);
    assert!(
        s.string_data_arrays[0]
            .data
            .iter()
            .any(|s| s == "[M+H]-NH3")
    );
}
#[test]
fn signed_charges_empty_ranges_and_consumed_zero_division() {
    let mut g = Generator::default();
    only(&mut g, "");
    g.options.add_k_linked_ions = false;
    let mut s = MSSpectrum::default();
    g.get_xlink_ion_spectrum(&mut s, &peptide(), 0, 1000., true, -2, -2, 0)
        .unwrap();
    assert_eq!(s.len(), 3);
    assert!(s.peaks.iter().all(|p| p.mz < 0.));
    assert!(s.integer_data_arrays[0].data.iter().all(|&z| z == -2));
    let mut inverted = MSSpectrum::default();
    g.get_xlink_ion_spectrum(&mut inverted, &peptide(), 0, 1000., true, 4, 2, 0)
        .unwrap();
    assert_eq!(inverted.len(), 3);
    let old = s.clone();
    assert!(
        g.get_xlink_ion_spectrum(&mut s, &peptide(), 0, 1000., true, 0, 0, 0)
            .is_err()
    );
    assert_eq!(s, old);
    g.options.add_precursor_peaks = false;
    g.get_xlink_ion_spectrum(&mut s, &peptide(), usize::MAX, 1000., true, 0, 0, 0)
        .unwrap();
    let mut empty = MSSpectrum::default();
    g.get_linear_ion_spectrum(&mut empty, &peptide(), 3, true, -1, 0)
        .unwrap();
    assert!(empty.is_empty());
    assert_eq!(empty.string_data_arrays.len(), 1);
}
#[test]
fn append_preserves_descriptions_and_moves_all_supported_arrays_stably() {
    let mut g = Generator::default();
    only(&mut g, "b");
    let p = AASequence::parse("AA").unwrap();
    let mut generated = MSSpectrum::default();
    g.get_linear_ion_spectrum(&mut generated, &p, 1, true, 1, 0)
        .unwrap();
    let mut s = MSSpectrum::from_peaks(vec![
        Peak1D::new(1000., 5.),
        Peak1D::new(generated.peaks[0].mz, 7.),
    ]);
    s.metadata.insert("keep".into(), "value".into());
    s.name = "caller".into();
    let history = Arc::new(DataProcessing::default());
    let mut names = DataArray::new(
        "first unrecognized name",
        vec!["old high".into(), "old tie".into()],
    );
    names.metadata.insert("keep".into(), 2_i64.into());
    names.data_processing.push(history.clone());
    s.string_data_arrays.push(names);
    s.integer_data_arrays
        .push(DataArray::new("first charge", vec![8, 9]));
    s.float_data_arrays
        .push(DataArray::new("empty placeholder", Vec::new()));
    g.get_linear_ion_spectrum(&mut s, &p, 1, true, 1, 0)
        .unwrap();
    assert_eq!(
        s.string_data_arrays[0].data,
        vec!["old tie", "[alpha|ci$b1]", "old high"]
    );
    assert_eq!(s.integer_data_arrays[0].data, vec![9, 1, 8]);
    assert!(Arc::ptr_eq(
        &s.string_data_arrays[0].data_processing[0],
        &history
    ));
    assert_eq!(s.string_data_arrays[0].metadata["keep"], 2_i64.into());
    assert_eq!(s.name, "caller");
    assert_eq!(s.ms_level, 1);
    assert!(s.precursors.is_empty());
    assert_eq!(s.metadata["keep"], "value");
}
#[test]
fn missing_old_annotations_are_padded_but_unrelated_data_reject_atomically() {
    let g = Generator::default();
    let mut s = MSSpectrum::from_peaks(vec![Peak1D::new(0., 1.)]);
    g.get_linear_ion_spectrum(&mut s, &peptide(), 3, true, 1, 0)
        .unwrap();
    assert_eq!(s.string_data_arrays[0].data[0], "");
    assert_eq!(s.integer_data_arrays[0].data[0], 0);
    for kind in 0..3 {
        let mut s = MSSpectrum::from_peaks(vec![Peak1D::new(0., 1.)]);
        match kind {
            0 => s
                .float_data_arrays
                .push(DataArray::new("unrelated", vec![1.])),
            1 => {
                s.integer_data_arrays.push(DataArray::new("first", vec![1]));
                s.integer_data_arrays
                    .push(DataArray::new("second", vec![2]));
            }
            _ => {
                s.string_data_arrays
                    .push(DataArray::new("first", vec!["x".into()]));
                s.string_data_arrays
                    .push(DataArray::new("second", vec!["y".into()]));
            }
        }
        let before = s.clone();
        assert!(
            g.get_linear_ion_spectrum(&mut s, &peptide(), 3, true, 1, 0)
                .is_err()
        );
        assert_eq!(s, before);
    }
}
#[test]
fn position_empty_loss_and_small_cx_boundaries_are_checked() {
    let mut g = Generator::default();
    for p in [
        AASequence::parse("AA").unwrap(),
        AASequence::parse("A").unwrap(),
    ] {
        let mut s = MSSpectrum::default();
        assert!(
            g.get_linear_ion_spectrum(&mut s, &p, p.len() + 1, true, 1, 0)
                .is_err()
        );
        assert!(s.is_empty());
    }
    g.options.add_losses = true;
    assert!(
        g.get_linear_ion_spectrum(
            &mut MSSpectrum::default(),
            &AASequence::default(),
            0,
            true,
            1,
            0
        )
        .is_err()
    );
    for code in ["O", "J", "B", "Z", "X"] {
        assert!(
            g.get_linear_ion_spectrum(
                &mut MSSpectrum::default(),
                &AASequence::parse(code).unwrap(),
                0,
                true,
                1,
                0
            )
            .is_err()
        );
    }
    g.options.add_losses = false;
    only(&mut g, "cx");
    assert!(
        g.get_linear_ion_spectrum(
            &mut MSSpectrum::default(),
            &AASequence::parse("A").unwrap(),
            0,
            true,
            1,
            0
        )
        .is_err()
    );
    let mut s = MSSpectrum::default();
    g.get_crosslink_ion_spectrum(&mut s, &ProteinProteinCrossLink::default(), true, 1, 1)
        .unwrap();
    assert_eq!(s, MSSpectrum::default());
}
#[test]
fn all_failure_stages_leave_the_target_unchanged() {
    let mut g = Generator::default();
    let target = MSSpectrum::from_peaks(vec![Peak1D::new(123., 4.)]);
    for mode in 0..5 {
        let mut s = target.clone();
        match mode {
            0 => g.limits.max_work = 1,
            1 => {
                g.limits = Default::default();
                g.limits.max_bytes = 100;
            }
            2 => {
                g.limits = Default::default();
                g.limits.max_peaks = 2;
            }
            3 => {
                g.limits = Default::default();
                g.options.b_intensity = f64::MAX;
            }
            _ => g.options.b_intensity = f64::NAN,
        };
        assert!(
            g.get_linear_ion_spectrum(&mut s, &peptide(), 3, true, 3, 0)
                .is_err()
        );
        assert_eq!(s, target);
    }
}
#[test]
fn known_modifications_and_terminal_deltas_follow_source_arithmetic() {
    let mut g = Generator::default();
    only(&mut g, "b");
    let p = AASequence::parse(".(Acetyl)M(Oxidation)A").unwrap();
    let mut s = MSSpectrum::default();
    g.get_linear_ion_spectrum(&mut s, &p, 1, true, 1, 0)
        .unwrap();
    let methionine = EmpiricalFormula::parse("C5H11NO2S").unwrap();
    let oxygen = EmpiricalFormula::parse("O").unwrap();
    let water = EmpiricalFormula::parse("H2O").unwrap().mono_mass();
    let expected = 1.007276466771
        + p.n_terminal_modification()
            .unwrap()
            .diff_mono_mass()
            .unwrap()
        + methionine.checked_add(&oxygen).unwrap().mono_mass()
        - water;
    near(s.peaks[0].mz, expected);
    let record = ResidueModification::from_record(ModificationRecord {
        name: "custom".into(),
        origin: Some('A'),
        term_specificity: TermSpecificity::Anywhere,
        diff_mono_mass: 5.,
        mono_mass: 200.,
        ..Default::default()
    })
    .unwrap();
    let db = ModificationsDB::from_records(vec![record]).unwrap();
    let p = AASequence::parse_with_registry("A(custom)A", &db).unwrap();
    s = MSSpectrum::default();
    g.get_linear_ion_spectrum(&mut s, &p, 1, true, 1, 0)
        .unwrap();
    near(s.peaks[0].mz, 1.007276466771 + 200. - water);
}

#[test]
fn linked_length_endpoint_skips_missing_losses_but_keeps_base_and_isotope() {
    let mut g = Generator::default();
    only(&mut g, "y");
    g.options.add_losses = true;
    g.options.add_isotopes = true;
    g.options.add_precursor_peaks = false;
    g.options.add_k_linked_ions = false;
    let p = AASequence::parse("KS").unwrap();
    let mut pair = ProteinProteinCrossLink::new(100.).unwrap();
    pair.alpha = Some(Arc::new(p.clone()));
    pair.beta = Some(Arc::new(p.clone()));
    pair.cross_link_position = (2, 2);
    for kind in 0..2 {
        let mut s = MSSpectrum::default();
        if kind == 0 {
            g.get_xlink_ion_spectrum(&mut s, &p, 2, 1000., true, 2, 2, 0)
                .unwrap();
        } else {
            g.get_crosslink_ion_spectrum(&mut s, &pair, true, 2, 2)
                .unwrap();
        }
        let names = &s.string_data_arrays[0].data;
        assert_eq!(
            names
                .iter()
                .filter(|n| n.as_str() == "[alpha|xi$y0]")
                .count(),
            2
        );
        assert!(!names.iter().any(|n| n.starts_with("[alpha|xi$y0-")));
    }
}
#[test]
fn empty_alpha_suppresses_beta_ladders_but_not_k_linked_or_precursor() {
    let g = Generator::default();
    let mut p = ProteinProteinCrossLink::new(100.).unwrap();
    p.alpha = Some(Arc::new(AASequence::default()));
    p.beta = Some(Arc::new(AASequence::parse("AKA").unwrap()));
    p.cross_link_position = (0, 1);
    let mut s = MSSpectrum::default();
    g.get_crosslink_ion_spectrum(&mut s, &p, false, 1, 1)
        .unwrap();
    assert_eq!(s.len(), 4);
    assert!(
        s.string_data_arrays[0]
            .data
            .iter()
            .all(|n| n.starts_with("[M+H]") || n == "[K-linked-alpha]")
    );
}
