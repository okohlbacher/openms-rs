// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::kernel::{Data2DLimits, Peak2D, RichPeak2D};
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
    m.metadata
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
    let ptr = exp.metadata["old"].as_ptr();
    let peak_ptr = exp.spectra[0].peaks.as_ptr();
    let old = exp.set_2d_data(&input).unwrap();
    assert_eq!(old, before);
    assert_eq!(old.metadata["old"].as_ptr(), ptr);
    assert_eq!(old.spectra[0].peaks.as_ptr(), peak_ptr);
    assert!(exp.metadata.is_empty());
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
