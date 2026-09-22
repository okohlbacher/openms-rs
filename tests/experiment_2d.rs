// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::concept::constants::C13C12_MASSDIFF_U;
use openms::concept::constants::user_param::NUM_OF_MASSTRACES;
use openms::kernel::{BaseFeature, ConsensusFeature, Data2DLimits, Feature, Peak2D, RichPeak2D};
use openms::metadata::{MetaValue, Unit};
use openms::{Error, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};
fn scan(rt: f64, level: u32, points: &[(f64, f32)]) -> MSSpectrum {
    MSSpectrum {
        rt,
        ms_level: level,
        peaks: points.iter().map(|&(mz, i)| Peak1D::new(mz, i)).collect(),
        ..Default::default()
    }
}
fn old_map() -> MSExperiment {
    let mut m = MSExperiment {
        spectra: vec![scan(99., 2, &[(50., 60.)])],
        chromatograms: vec![MSChromatogram::new()],
        ..Default::default()
    };
    m.settings
        .metadata
        .insert("old".into(), "preserve previous ownership".into());
    m.spectra[0]
        .metadata
        .insert("annotation".into(), "owned".into());
    m
}
#[test]
fn source_ms1_export_literal_sequence_and_append() {
    let exp = MSExperiment {
        spectra: vec![
            scan(11.1, 1, &[(5., 47.11), (10., 48.11), (15., 48.11)]),
            scan(11.5, 2, &[(6., 48.11), (11., 48.11)]),
            scan(12.2, 1, &[(20., 48.11), (25., 48.11)]),
            scan(12.5, 2, &[(21., 48.11), (26., 48.11), (31., 48.11)]),
        ],
        ..Default::default()
    };
    let expected = vec![
        Peak2D::new(11.1, 5., 47.11),
        Peak2D::new(11.1, 10., 48.11),
        Peak2D::new(11.1, 15., 48.11),
        Peak2D::new(12.2, 20., 48.11),
        Peak2D::new(12.2, 25., 48.11),
    ];
    assert_eq!(exp.get_2d_data().unwrap(), expected);
    let mut output = vec![Peak2D::new(1., 2., 3.)];
    exp.append_2d_data(&mut output).unwrap();
    assert_eq!(output[0], Peak2D::new(1., 2., 3.));
    assert_eq!(output[1..], expected);
}
#[test]
fn source_plain_roundtrip_replaces_all_current_data_and_returns_old_ownership() {
    let input = [
        Peak2D::new(2., 3., 1.),
        Peak2D::new(5., 6., 4.),
        Peak2D::new(8.5, 9.5, 7.5),
    ];
    let mut exp = old_map();
    let before = exp.clone();
    let ptr = exp.settings.metadata["old"].as_str().unwrap().as_ptr();
    let peak_ptr = exp.spectra[0].peaks.as_ptr();
    let old = exp.set_2d_data(&input).unwrap();
    assert_eq!(old, before);
    assert_eq!(old.settings.metadata["old"].as_str().unwrap().as_ptr(), ptr);
    assert_eq!(old.spectra[0].peaks.as_ptr(), peak_ptr);
    assert!(exp.settings.metadata.is_empty());
    assert!(exp.chromatograms.is_empty());
    assert_eq!(exp.spectra.len(), 3);
    assert!(
        exp.spectra
            .iter()
            .all(|s| s.ms_level == 1 && s.metadata.is_empty())
    );
    assert_eq!(exp.get_2d_data().unwrap(), input);
    let replaced = exp.set_2d_data(&[]).unwrap();
    assert_eq!(replaced.get_2d_data().unwrap(), input);
    assert_eq!(exp, MSExperiment::new());
}
#[test]
fn source_rich_metadata_literals_missing_nan_and_duplicate_names() {
    let mut input = vec![
        RichPeak2D::new(2., 3., 1.),
        RichPeak2D::new(5., 6., 4.),
        RichPeak2D::new(8.5, 9.5, 7.5),
    ];
    input[0]
        .metadata
        .insert("meta1".into(), MetaValue::try_from(111.1).unwrap());
    input[2]
        .metadata
        .insert("meta3".into(), MetaValue::try_from(333.3).unwrap());
    input[0].unique_id = 99;
    let mut exp = MSExperiment::new();
    exp.set_2d_data_rich(&input, &["meta1".into(), "meta3".into(), "meta1".into()])
        .unwrap();
    assert_eq!(exp.spectra[0].float_data_arrays[0].data, [111.1f32]);
    assert_eq!(exp.spectra[2].float_data_arrays[1].data, [333.3f32]);
    assert_eq!(exp.spectra[0].float_data_arrays[2].data, [111.1f32]);
    assert!(exp.spectra[0].float_data_arrays[1].data[0].is_nan());
    assert!(
        exp.spectra[1]
            .float_data_arrays
            .iter()
            .all(|a| a.data[0].is_nan())
    );
    assert_eq!(
        exp.spectra[0]
            .float_data_arrays
            .iter()
            .map(|a| a.name.as_str())
            .collect::<Vec<_>>(),
        ["meta1", "meta3", "meta1"]
    );
    assert!(
        exp.spectra
            .iter()
            .flat_map(|s| &s.float_data_arrays)
            .all(|a| a.metadata.is_empty() && a.data_processing.is_empty())
    );
    assert_eq!(
        exp.get_2d_data().unwrap(),
        input.iter().map(|p| p.peak).collect::<Vec<_>>()
    );
}
#[test]
fn grouping_is_exact_rt_and_preserves_unsorted_mz_and_zero_sign() {
    let input = [
        Peak2D::new(-0., 5., 1.),
        Peak2D::new(0., 2., 2.),
        Peak2D::new(0.1, 9., 3.),
        Peak2D::new(0.1, 1., 4.),
    ];
    let mut exp = MSExperiment::new();
    exp.set_2d_data(&input).unwrap();
    assert_eq!(exp.spectra.len(), 2);
    assert_eq!(exp.spectra[0].rt.to_bits(), (-0f64).to_bits());
    assert_eq!(
        exp.spectra[0]
            .peaks
            .iter()
            .map(|p| p.mz)
            .collect::<Vec<_>>(),
        [5., 2.]
    );
    assert_eq!(
        exp.spectra[1]
            .peaks
            .iter()
            .map(|p| p.mz)
            .collect::<Vec<_>>(),
        [9., 1.]
    );
    let out = exp.get_2d_data().unwrap();
    assert_eq!(out, input);
    assert_eq!(out[1].rt().to_bits(), (-0f64).to_bits());
    // Unlike filtered area extraction, unfiltered export permits either ordering.
    exp.spectra.reverse();
    assert_eq!(exp.get_2d_data().unwrap()[0].rt(), 0.1);
}
#[test]
fn direct_integer_to_float_metadata_conversion_avoids_double_rounding() {
    let value = (1i64 << 62) + (1i64 << 38) + 1;
    let expected = f32::from_bits(((1u64 << 62) as f32).to_bits() + 1);
    assert_ne!(expected, (value as f64) as f32);
    let mut input = RichPeak2D::new(1., 2., 3.);
    input.metadata.insert("i".into(), value.into());
    input.metadata.insert(
        "f".into(),
        MetaValue::try_from(2.5)
            .unwrap()
            .with_unit(Unit::new("UO:0000010", "s", "UO").unwrap())
            .unwrap(),
    );
    let mut exp = MSExperiment::new();
    exp.set_2d_data_rich(&[input], &["i".into(), "f".into()])
        .unwrap();
    assert_eq!(exp.spectra[0].float_data_arrays[0].data, [expected]);
    assert_eq!(exp.spectra[0].float_data_arrays[1].data, [2.5]);
}
#[test]
fn late_import_and_export_errors_are_atomic_and_unused_payload_is_ignored() {
    let mut exp = old_map();
    let before = exp.clone();
    for input in [
        vec![Peak2D::new(2., 1., 1.), Peak2D::new(1., 1., 1.)],
        vec![Peak2D::new(1., 1., 1.), Peak2D::new(f64::NAN, 1., 1.)],
        vec![Peak2D::new(1., f64::INFINITY, 1.)],
    ] {
        assert!(exp.set_2d_data(&input).is_err());
        assert_eq!(exp, before);
    }
    let mut bad = RichPeak2D::new(2., 3., 4.);
    bad.metadata.insert("x".into(), "not numeric".into());
    assert!(
        exp.set_2d_data_rich(&[RichPeak2D::new(1., 1., 1.), bad.clone()], &["x".into()])
            .is_err()
    );
    assert_eq!(exp, before);
    bad.metadata.insert("x".into(), MetaValue::default());
    assert!(exp.set_2d_data_rich(&[bad.clone()], &["x".into()]).is_err());
    assert_eq!(exp, before);
    bad.metadata
        .insert("x".into(), MetaValue::try_from(f64::MAX).unwrap());
    assert!(exp.set_2d_data_rich(&[bad.clone()], &["x".into()]).is_err());
    assert_eq!(exp, before);
    exp.set_2d_data_rich(&[bad], &[]).unwrap(); // unrequested value not converted
    exp.spectra.push(scan(f64::NAN, 2, &[(f64::NAN, f32::NAN)]));
    assert_eq!(exp.get_2d_data().unwrap().len(), 1);
    exp.spectra.push(scan(f64::NAN, 1, &[])); // no point consumes its RT
    assert_eq!(exp.get_2d_data().unwrap().len(), 1);
    exp.spectra.push(scan(5., 1, &[(6., f32::NAN)]));
    let mut output = vec![Peak2D::new(9., 8., 7.)];
    let output_before = output.clone();
    assert!(exp.append_2d_data(&mut output).is_err());
    assert_eq!(output, output_before);
}
#[test]
fn shared_limits_bound_new_points_names_arrays_and_existing_export_prefix() {
    let input = [Peak2D::new(1., 2., 3.), Peak2D::new(2., 3., 4.)];
    for limits in [
        Data2DLimits {
            max_points: 1,
            ..Default::default()
        },
        Data2DLimits {
            max_spectra: 1,
            ..Default::default()
        },
        Data2DLimits {
            max_work: 1,
            ..Default::default()
        },
        Data2DLimits {
            max_bytes: 1,
            ..Default::default()
        },
    ] {
        let mut exp = old_map();
        let before = exp.clone();
        assert!(exp.set_2d_data_with_limits(&input, limits).is_err());
        assert_eq!(exp, before);
    }
    let rich = input.map(RichPeak2D::from);
    for limits in [
        Data2DLimits {
            max_arrays: 1,
            ..Default::default()
        },
        Data2DLimits {
            max_name_bytes: 1,
            ..Default::default()
        },
        Data2DLimits {
            max_work: 10,
            ..Default::default()
        },
    ] {
        let mut exp = old_map();
        let before = exp.clone();
        assert!(
            exp.set_2d_data_rich_with_limits(&rich, &["long".into()], limits)
                .is_err()
        );
        assert_eq!(exp, before);
    }
    let mut exp = MSExperiment::new();
    exp.set_2d_data(&input).unwrap();
    let mut out = vec![input[0]];
    assert!(
        exp.append_2d_data_with_limits(
            &mut out,
            Data2DLimits {
                max_points: 2,
                ..Default::default()
            }
        )
        .is_err()
    );
    assert_eq!(out, [input[0]]);
    let zero = Data2DLimits {
        max_points: 0,
        max_spectra: 0,
        max_arrays: 0,
        max_name_bytes: 0,
        max_work: 0,
        max_bytes: 0,
    };
    let old = exp
        .set_2d_data_rich_with_limits(&[], &["unused".into()], zero)
        .unwrap();
    assert_eq!(old.get_2d_data().unwrap(), input);
    assert!(exp.get_2d_data_with_limits(zero).unwrap().is_empty());
}
#[test]
fn bounded_random_sorted_roundtrips_and_unsorted_export_oracle() {
    let mut seed = 37u64;
    let mut next = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        (seed >> 32) as usize
    };
    for _ in 0..100 {
        let mut points = Vec::new();
        let mut rt = -2.;
        for _ in 0..next() % 30 {
            rt += (next() % 3) as f64 / 10.;
            points.push(Peak2D::new(
                rt,
                (next() % 200) as f64 / 10.,
                (next() % 50) as f32 - 20.,
            ));
        }
        let mut exp = MSExperiment::new();
        exp.set_2d_data(&points).unwrap();
        assert_eq!(exp.get_2d_data().unwrap(), points);
        exp.spectra.reverse();
        if let Some(scan) = exp.spectra.first_mut() {
            scan.ms_level = 2;
        }
        let expected: Vec<_> = exp
            .spectra
            .iter()
            .filter(|s| s.ms_level == 1)
            .flat_map(|s| {
                s.peaks
                    .iter()
                    .map(move |p| Peak2D::new(s.rt, p.mz, p.intensity))
            })
            .collect();
        assert_eq!(exp.get_2d_data().unwrap(), expected);
    }
    let mut exp = MSExperiment::new();
    assert!(matches!(
        exp.set_2d_data(&[Peak2D::new(1., 1., 1.), Peak2D::new(0., 1., 1.)]),
        Err(Error::UnsortedData)
    ));
}

// --------------------------------------------------------------------------
// set2DData<add_mass_traces = true>
// --------------------------------------------------------------------------

fn feature(rt: f64, mz: f64, intensity: f32) -> Feature {
    Feature::from(BaseFeature::new(rt, mz, intensity))
}
fn float(value: f64) -> MetaValue {
    MetaValue::try_from(value).unwrap()
}
fn traces(feature: &mut Feature, count: i64, intensities: &[MetaValue]) {
    feature
        .metadata
        .insert(NUM_OF_MASSTRACES.into(), count.into());
    for (i, value) in intensities.iter().enumerate() {
        feature
            .metadata
            .insert(format!("masstrace_intensity_{i}"), value.clone());
    }
}
/// Each spectrum as its retention time and its peaks' bit patterns.
fn peak_bits(exp: &MSExperiment) -> Vec<(f64, Vec<(u64, u32)>)> {
    exp.spectra
        .iter()
        .map(|s| {
            assert_eq!(s.ms_level, 1);
            assert!(s.float_data_arrays.is_empty());
            let peaks = s.peaks.iter().map(|p| (p.mz.to_bits(), p.intensity.to_bits()));
            (s.rt, peaks.collect())
        })
        .collect()
}
fn bits(mz: f64, intensity: f32) -> (u64, u32) {
    (mz.to_bits(), intensity.to_bits())
}
#[test]
fn source_set_2d_data_mass_traces_literals_and_precondition() {
    // START_SECTION((template <bool add_mass_traces, class Container> void set2DData(const Container& cont, const StringList& store_metadata_names = StringList())))
    // MSExperiment_test.cpp:359-393, the mass-trace half and the closing
    // precondition. The plain and rich halves (:310-357) are
    // `source_plain_roundtrip_replaces_all_current_data_and_returns_old_ownership`
    // and `source_rich_metadata_literals_missing_nan_and_duplicate_names`.
    let f1 = feature(8.5, 9.5, 7.5);
    let mut f2 = feature(18.5, 19.5, 17.5);
    traces(&mut f2, 2, &[float(11.0), float(12.0)]);
    f2.charge = 2;
    let mut exp = MSExperiment::new();
    exp.set_2d_data_mass_traces(&[f1, f2]).unwrap();
    // `get2DData(fm_out)` fills RT, m/z and intensity only, so those three are
    // the whole comparison with `fm2`: f1, then f2 once per trace, the m/z
    // advanced by `C13C12_MASSDIFF_U/2*i` and the intensity replaced.
    let mz0 = 19.5 + C13C12_MASSDIFF_U / 2.0 * 0.0;
    let mz1 = mz0 + C13C12_MASSDIFF_U / 2.0 * 1.0;
    assert_eq!(
        exp.get_2d_data().unwrap(),
        [
            Peak2D::new(8.5, 9.5, 7.5),
            Peak2D::new(18.5, mz0, 11.0),
            Peak2D::new(18.5, mz1, 12.0),
        ]
    );
    // The Release oracle's value, `0x1.4006dee6b5134p+4` (case mass_section11).
    assert_eq!(mz1.to_bits(), 0x4034_006d_ee6b_5134);
    assert_eq!(exp.spectra.len(), 2);

    // "test precondition (input sorted by RT)": `input.push_back(p1)`.
    let input = [
        Peak2D::new(2.0, 3.0, 1.0),
        Peak2D::new(5.0, 6.0, 4.0),
        Peak2D::new(8.5, 9.5, 7.5),
        Peak2D::new(2.0, 3.0, 1.0),
    ];
    assert!(matches!(exp.set_2d_data(&input), Err(Error::UnsortedData)));
}
#[test]
fn release_set_2d_data_mass_traces_empty_spectrum_negative_charge_and_integer_intensity() {
    // Release oracle (Linux x86_64 Release build, `../oracle/msexperiment-close`),
    // case mass_zero_negative_int_values. A retention time whose only feature
    // records zero traces still opens a spectrum; charge -2 spaces the traces
    // downward; an integer intensity converts; a -0.0 intensity keeps its sign.
    let mut zero = feature(1.0, 100.0, 5.0);
    traces(&mut zero, 0, &[]);
    let plain = feature(2.0, 200.0, 6.0);
    let mut negative = feature(2.0, 150.0, 7.0);
    traces(&mut negative, 3, &[1i64.into(), float(2.5), float(-0.0)]);
    negative.charge = -2;
    let mut triple = feature(3.0, 300.0, 8.0);
    traces(&mut triple, 2, &[float(4.0), float(0.0)]);
    triple.charge = 3;
    let mut exp = MSExperiment::new();
    exp.set_2d_data_mass_traces(&[zero, plain, negative, triple])
        .unwrap();
    assert_eq!(
        peak_bits(&exp),
        [
            (1.0, vec![]),
            (
                2.0,
                vec![
                    bits(200.0, 6.0),
                    bits(150.0, 1.0),
                    (0x4062_aff2_4232_95da, 2.5f32.to_bits()),
                    (0x4062_9fe4_8465_2bb3, (-0.0f32).to_bits()),
                ]
            ),
            (
                3.0,
                vec![bits(300.0, 4.0), (0x4072_c559_e9ef_2362, 0.0f32.to_bits())]
            ),
        ]
    );
}
#[test]
fn release_set_2d_data_mass_traces_errors_follow_source_order_and_replace_nothing() {
    // Release oracle cases mass_missing_trace, mass_count_huge,
    // mass_count_double, mass_count_negative, mass_count_string,
    // mass_intensity_empty_value, mass_unsorted, mass_error_order and
    // mass_intensity_overflow. The source has cleared the experiment and
    // added the peaks before the failing one; the port replaces nothing.
    let mut exp = old_map();
    let before = exp.clone();
    let mut run = |features: &[Feature]| exp.set_2d_data_mass_traces(features);
    let message = |result: openms::Result<MSExperiment>| match result {
        Err(Error::InvalidValue(message)) => message,
        other => panic!("expected an invalid value, got {other:?}"),
    };

    let mut missing = feature(1.0, 100.0, 5.0);
    traces(&mut missing, 2, &[float(1.0)]);
    assert_eq!(
        message(run(&[missing])),
        "Meta value 'masstrace_intensity_1' expected but not found in container."
    );
    // A huge count fails at the first missing name, not by running long.
    let mut huge = feature(1.0, 100.0, 5.0);
    traces(&mut huge, 1_000_000_000_000, &[float(1.0)]);
    assert_eq!(
        message(run(&[huge])),
        "Meta value 'masstrace_intensity_1' expected but not found in container."
    );
    // The source's DataValue conversions throw ConversionError here.
    for count in [float(2.0), (-1i64).into(), "2".into()] {
        let mut f = feature(1.0, 100.0, 5.0);
        f.metadata.insert(NUM_OF_MASSTRACES.into(), count);
        message(run(&[f]));
    }
    let mut empty = feature(1.0, 100.0, 5.0);
    traces(&mut empty, 1, &[MetaValue::default()]);
    message(run(&[empty]));

    assert!(matches!(
        run(&[feature(5.0, 100.0, 5.0), feature(3.0, 100.0, 5.0)]),
        Err(Error::UnsortedData)
    ));
    // The second feature's missing trace comes before the third's RT regression.
    let mut second = feature(7.0, 100.0, 5.0);
    traces(&mut second, 1, &[]);
    assert_eq!(
        message(run(&[
            feature(5.0, 100.0, 5.0),
            second,
            feature(1.0, 100.0, 5.0)
        ])),
        "Meta value 'masstrace_intensity_0' expected but not found in container."
    );
    // The source stores `float(1e40)`, infinity; the port's 2D import refuses
    // a nonfinite intensity, as it does for plain and rich points.
    let mut overflow = feature(1.0, 100.0, 5.0);
    traces(&mut overflow, 1, &[float(1e40)]);
    message(run(&[overflow]));

    assert_eq!(exp, before);
}
#[test]
fn release_set_2d_data_mass_traces_clears_settings_and_chromatograms() {
    // Release oracle cases mass_clears and mass_empty; the old map comes back.
    let mut exp = old_map();
    exp.settings.comment = "old".into();
    let old = exp
        .set_2d_data_mass_traces(&[feature(1.0, 100.0, 5.0)])
        .unwrap();
    assert_eq!(old.settings.comment, "old");
    assert_eq!(old.chromatograms.len(), 1);
    assert!(exp.settings.comment.is_empty() && exp.settings.metadata.is_empty());
    assert!(exp.chromatograms.is_empty());
    assert_eq!(exp.get_2d_data().unwrap(), [Peak2D::new(1.0, 100.0, 5.0)]);
    exp.set_2d_data_mass_traces::<Feature>(&[]).unwrap();
    assert_eq!(exp, MSExperiment::new());
}
#[test]
fn mass_trace_expansion_is_bounded_and_accepts_consensus_features() {
    // Native bounds: the expanded peaks count against `max_points`, the
    // retention times against `max_spectra`.
    let mut three = feature(1.0, 100.0, 5.0);
    traces(&mut three, 3, &[float(1.0), float(2.0), float(3.0)]);
    let mut exp = MSExperiment::new();
    let two_points = Data2DLimits {
        max_points: 2,
        ..Default::default()
    };
    assert!(matches!(
        exp.set_2d_data_mass_traces_with_limits(&[three.clone()], two_points),
        Err(Error::InvalidValue(_))
    ));
    let one_spectrum = Data2DLimits {
        max_spectra: 1,
        ..Default::default()
    };
    let two_rts = [feature(1.0, 1.0, 1.0), feature(2.0, 1.0, 1.0)];
    assert!(matches!(
        exp.set_2d_data_mass_traces_with_limits(&two_rts, one_spectrum),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(exp, MSExperiment::new());
    let three_points = Data2DLimits {
        max_points: 3,
        ..Default::default()
    };
    exp.set_2d_data_mass_traces_with_limits(&[three], three_points)
        .unwrap();
    assert_eq!(exp.spectra[0].peaks.len(), 3);

    // The source template takes any container of feature-like values.
    let mut consensus = ConsensusFeature::from(BaseFeature::new(4.0, 50.0, 2.0));
    consensus
        .metadata
        .insert(NUM_OF_MASSTRACES.into(), 1i64.into());
    consensus
        .metadata
        .insert("masstrace_intensity_0".into(), float(9.0));
    exp.set_2d_data_mass_traces(&[consensus]).unwrap();
    assert_eq!(exp.get_2d_data().unwrap(), [Peak2D::new(4.0, 50.0, 9.0)]);
}
