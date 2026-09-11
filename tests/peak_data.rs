// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::kernel::{AreaBounds, FlatPeakData, PeakDataLimits, SpectrumPeakData};
use openms::{Error, MSExperiment, MSSpectrum, Peak1D};

fn scan(rt: f64, level: u32, points: &[(f64, f32)]) -> MSSpectrum {
    MSSpectrum {
        rt,
        ms_level: level,
        peaks: points.iter().map(|&(mz, i)| Peak1D::new(mz, i)).collect(),
        ..Default::default()
    }
}
fn experiment(scans: Vec<MSSpectrum>) -> MSExperiment {
    MSExperiment {
        spectra: scans,
        ..Default::default()
    }
}
fn bounds(rt: (f64, f64), mz: (f64, f64)) -> AreaBounds {
    AreaBounds::new(rt.0, rt.1, mz.0, mz.1).unwrap()
}
fn source() -> MSExperiment {
    experiment(vec![
        scan(1., 1, &[(100., 1000.), (200., 2000.)]),
        scan(2., 1, &[(150., 1500.), (250., 2500.)]),
        scan(3., 2, &[(175., 1750.)]),
    ])
}
fn flat(rt: &[f32], mz: &[f32], intensity: &[f32]) -> FlatPeakData {
    FlatPeakData {
        rt: rt.into(),
        mz: mz.into(),
        intensity: intensity.into(),
    }
}
fn rows(rt: &[f32], mz: &[&[f32]], intensity: &[&[f32]]) -> SpectrumPeakData {
    SpectrumPeakData {
        rt: rt.into(),
        mz: mz.iter().map(|row| row.to_vec()).collect(),
        intensity: intensity.iter().map(|row| row.to_vec()).collect(),
    }
}

#[test]
fn source_per_spectrum_literal_queries() {
    let exp = source();
    let cases = [
        (
            bounds((0., 4.), (0., 300.)),
            1,
            rows(
                &[1., 2.],
                &[&[100., 200.], &[150., 250.]],
                &[&[1000., 2000.], &[1500., 2500.]],
            ),
        ),
        (
            bounds((1.5, 2.5), (0., 300.)),
            1,
            rows(&[2.], &[&[150., 250.]], &[&[1500., 2500.]]),
        ),
        (
            bounds((0., 4.), (120., 180.)),
            1,
            rows(&[2.], &[&[150.]], &[&[1500.]]),
        ),
        (
            bounds((0., 4.), (0., 300.)),
            2,
            rows(&[3.], &[&[175.]], &[&[1750.]]),
        ),
        (bounds((5., 6.), (0., 300.)), 1, SpectrumPeakData::default()),
    ];
    for (b, level, expected) in cases {
        assert_eq!(
            exp.get_2d_peak_data_per_spectrum(b, level).unwrap(),
            expected
        );
    }
}

#[test]
fn source_flat_literal_queries() {
    let exp = source();
    let cases = [
        (
            bounds((0., 4.), (0., 300.)),
            flat(
                &[1., 1., 2., 2.],
                &[100., 200., 150., 250.],
                &[1000., 2000., 1500., 2500.],
            ),
        ),
        (
            bounds((1.5, 2.5), (0., 300.)),
            flat(&[2., 2.], &[150., 250.], &[1500., 2500.]),
        ),
        (
            bounds((0., 4.), (120., 180.)),
            flat(&[2.], &[150.], &[1500.]),
        ),
        (bounds((5., 6.), (0., 300.)), FlatPeakData::default()),
    ];
    for (b, expected) in cases {
        assert_eq!(exp.get_2d_peak_data(b, 1).unwrap(), expected);
    }
}

#[test]
fn source_grouping_merges_exact_rt_scans_and_splits_inexact_rt_per_peak() {
    let exp = experiment(vec![
        scan(0.1, 1, &[(1., 10.), (2., 20.)]),
        scan(0.1, 1, &[(3., 30.)]),
        scan(1., 1, &[(4., 40.)]),
        scan(1., 2, &[(5., 50.)]),
        scan(1., 1, &[]),
        scan(1., 1, &[(6., 60.), (7., 70.)]),
    ]);
    assert_eq!(
        exp.get_2d_peak_data_per_spectrum(AreaBounds::default(), 1)
            .unwrap(),
        rows(
            &[0.1, 0.1, 0.1, 1.],
            &[&[1.], &[2.], &[3.], &[4., 6., 7.]],
            &[&[10.], &[20.], &[30.], &[40., 60., 70.]]
        )
    );
    // Grouping compares raw RT, not rounded keys or spectrum identities.
    let rounded = f64::from(0.1f32);
    let exp = experiment(vec![
        scan(0.1, 1, &[(1., 1.)]),
        scan(rounded, 1, &[(2., 2.), (3., 3.)]),
    ]);
    assert_eq!(
        exp.get_2d_peak_data_per_spectrum(AreaBounds::default(), 1)
            .unwrap(),
        rows(&[0.1], &[&[1., 2., 3.]], &[&[1., 2., 3.]])
    );
}

#[test]
fn minus_one_requires_an_existing_row_and_preserves_its_rt() {
    let exp = experiment(vec![
        scan(-1., 1, &[(1., 2.), (2., 3.)]),
        scan(-1., 1, &[(3., 4.)]),
    ]);
    assert!(
        matches!(exp.get_2d_peak_data_per_spectrum(AreaBounds::default(),1), Err(Error::InvalidValue(s)) if s.contains("no existing output row"))
    );
    let mut target = rows(&[99.], &[&[100.]], &[&[101.]]);
    exp.append_2d_peak_data_per_spectrum(AreaBounds::default(), 1, &mut target)
        .unwrap();
    assert_eq!(
        target,
        rows(&[99.], &[&[100., 1., 2., 3.]], &[&[101., 2., 3., 4.]])
    );
    assert_eq!(
        exp.get_2d_peak_data(AreaBounds::default(), 1).unwrap(),
        flat(&[-1., -1., -1.], &[1., 2., 3.], &[2., 3., 4.])
    );
    let normal = experiment(vec![
        scan(-2., 1, &[(0., 1.)]),
        scan(-1., 1, &[(1., 2.)]),
        scan(0., 1, &[(2., 3.)]),
    ]);
    assert_eq!(
        normal
            .get_2d_peak_data_per_spectrum(AreaBounds::default(), 1)
            .unwrap(),
        rows(
            &[-2., -1., 0.],
            &[&[0.], &[1.], &[2.]],
            &[&[1.], &[2.], &[3.]]
        )
    );
}

#[test]
fn append_retains_prior_values_and_resets_source_grouping_each_call() {
    let exp = experiment(vec![scan(1., 1, &[(1., 2.)])]);
    let mut a = flat(&[9.], &[8.], &[7.]);
    let mut b = rows(&[1.], &[&[8.]], &[&[7.]]);
    for _ in 0..2 {
        exp.append_2d_peak_data(AreaBounds::default(), 1, &mut a)
            .unwrap();
        exp.append_2d_peak_data_per_spectrum(AreaBounds::default(), 1, &mut b)
            .unwrap();
    }
    assert_eq!(a, flat(&[9., 1., 1.], &[8., 1., 1.], &[7., 2., 2.]));
    assert_eq!(
        b,
        rows(
            &[1., 1., 1.],
            &[&[8.], &[1.], &[1.]],
            &[&[7.], &[2.], &[2.]]
        )
    );
    let b_before = b.clone();
    exp.append_2d_peak_data_per_spectrum(bounds((5., 6.), (0., 10.)), 1, &mut b)
        .unwrap();
    assert_eq!(b, b_before);
}

#[test]
fn filters_are_closed_before_f32_rounding_and_empty_scans_do_not_emit_rows() {
    let next_mz = f64::from_bits(100f64.to_bits() + 1);
    let next_rt = f64::from_bits(1f64.to_bits() + 1);
    let exp = experiment(vec![
        scan(1., 1, &[]),
        scan(1., 1, &[(100., -0.), (100., -1.), (next_mz, 2.)]),
        scan(next_rt, 1, &[(100., 3.)]),
    ]);
    let b = bounds((1., 1.), (100., 100.));
    let a = exp.get_2d_peak_data(b, 1).unwrap();
    assert_eq!(a, flat(&[1., 1.], &[100., 100.], &[-0., -1.]));
    assert_eq!(a.intensity[0].to_bits(), (-0f32).to_bits());
    assert_eq!(
        exp.get_2d_peak_data_per_spectrum(b, 1).unwrap(),
        rows(&[1.], &[&[100., 100.]], &[&[-0., -1.]])
    );
}

#[test]
fn source_size_level_conversion_is_explicit_and_zero_is_exact() {
    let exp = experiment(vec![
        scan(0., 0, &[(1., 1.)]),
        scan(1., 1, &[(2., 2.)]),
        scan(2., u32::MAX, &[(3., 3.)]),
    ]);
    assert_eq!(
        exp.get_2d_peak_data(AreaBounds::default(), 256).unwrap().mz,
        [1.]
    );
    assert_eq!(
        exp.get_2d_peak_data(AreaBounds::default(), 0).unwrap().mz,
        [1.]
    );
    assert_eq!(
        exp.get_2d_peak_data(AreaBounds::default(), 255).unwrap().mz,
        [3.]
    );
    if usize::BITS > 32 {
        let huge = usize::try_from((1u64 << 32) + 1).unwrap();
        assert_eq!(
            exp.get_2d_peak_data(AreaBounds::default(), huge)
                .unwrap()
                .mz,
            [2.]
        );
    }
}

#[test]
fn invalid_consumed_values_and_prefix_shapes_fail_atomically() {
    let mut exp = experiment(vec![scan(1., 1, &[(1., 1.), (f64::MAX, 2.)])]);
    let mut flat_target = flat(&[9.], &[8.], &[7.]);
    let flat_before = flat_target.clone();
    let mut row_target = rows(&[9.], &[&[8.]], &[&[7.]]);
    let row_before = row_target.clone();
    assert!(
        exp.append_2d_peak_data(AreaBounds::default(), 1, &mut flat_target)
            .is_err()
    );
    assert!(
        exp.append_2d_peak_data_per_spectrum(AreaBounds::default(), 1, &mut row_target)
            .is_err()
    );
    assert_eq!(flat_target, flat_before);
    assert_eq!(row_target, row_before);
    // An excluded peak's f32 conversion is not consumed.
    assert_eq!(
        exp.get_2d_peak_data(bounds((1., 1.), (1., 1.)), 1)
            .unwrap()
            .mz,
        [1.]
    );
    exp.spectra[0].peaks[1].mz = 2.;
    exp.spectra[0].peaks[1].intensity = f32::NAN;
    assert!(exp.get_2d_peak_data(AreaBounds::default(), 1).is_err());
    assert_eq!(
        exp.get_2d_peak_data(bounds((1., 1.), (1., 1.)), 1)
            .unwrap()
            .mz,
        [1.]
    );
    exp.spectra[0].peaks[1].intensity = 2.;
    flat_target.rt.clear();
    let before = flat_target.clone();
    assert!(
        exp.append_2d_peak_data(AreaBounds::default(), 1, &mut flat_target)
            .is_err()
    );
    assert_eq!(flat_target, before);
    row_target.intensity[0].clear();
    let before = row_target.clone();
    assert!(
        exp.append_2d_peak_data_per_spectrum(AreaBounds::default(), 1, &mut row_target)
            .is_err()
    );
    assert_eq!(row_target, before);
    // No selected points leave even unused malformed destinations unchanged.
    exp.append_2d_peak_data(bounds((10., 11.), (1., 2.)), 1, &mut flat_target)
        .unwrap();
    exp.append_2d_peak_data_per_spectrum(bounds((10., 11.), (1., 2.)), 1, &mut row_target)
        .unwrap();
    assert_eq!(row_target, before);
    exp.spectra.push(scan(2., 2, &[(2., 1.), (1., 1.)]));
    assert!(matches!(
        exp.get_2d_peak_data(bounds((1., 1.), (1., 1.)), 1),
        Err(Error::UnsortedData)
    ));
}

#[test]
fn cumulative_area_output_and_existing_prefix_limits() {
    let exp = experiment(vec![scan(1., 1, &[(1., 1.), (2., 2.)])]);
    for limits in [
        PeakDataLimits {
            max_peaks: 1,
            ..Default::default()
        },
        PeakDataLimits {
            max_output_points: 2,
            ..Default::default()
        },
        PeakDataLimits {
            max_work: 15,
            ..Default::default()
        },
        PeakDataLimits {
            max_bytes: 100,
            ..Default::default()
        },
    ] {
        let mut target = flat(&[9.], &[8.], &[7.]);
        let before = target.clone();
        assert!(
            exp.append_2d_peak_data_with_limits(AreaBounds::default(), 1, &mut target, limits)
                .is_err()
        );
        assert_eq!(target, before);
    }
    let mut target = flat(&[9.], &[8.], &[7.]);
    exp.append_2d_peak_data_with_limits(
        AreaBounds::default(),
        1,
        &mut target,
        PeakDataLimits {
            max_work: 16,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(target.mz, [8., 1., 2.]);
    let split = experiment(vec![scan(0.1, 1, &[(1., 1.), (2., 2.)])]);
    for limits in [
        PeakDataLimits {
            max_output_rows: 1,
            ..Default::default()
        },
        PeakDataLimits {
            max_bytes: 600,
            ..Default::default()
        },
    ] {
        let mut target = SpectrumPeakData::default();
        assert!(
            split
                .append_2d_peak_data_per_spectrum_with_limits(
                    AreaBounds::default(),
                    1,
                    &mut target,
                    limits
                )
                .is_err()
        );
        assert_eq!(target, SpectrumPeakData::default());
    }
    let zero = PeakDataLimits {
        max_spectra: 0,
        max_peaks: 0,
        max_output_points: 0,
        max_output_rows: 0,
        max_work: 0,
        max_bytes: 0,
    };
    assert_eq!(
        MSExperiment::new()
            .get_2d_peak_data_with_limits(AreaBounds::default(), 1, zero)
            .unwrap(),
        FlatPeakData::default()
    );
    assert_eq!(
        MSExperiment::new()
            .get_2d_peak_data_per_spectrum_with_limits(AreaBounds::default(), 1, zero)
            .unwrap(),
        SpectrumPeakData::default()
    );
}

#[test]
fn deterministic_queries_match_direct_enumeration_in_both_layouts() {
    let mut seed = 17u64;
    let mut next = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        (seed >> 32) as usize
    };
    for _ in 0..80 {
        let mut exp = MSExperiment::new();
        let mut rt = 0.;
        for _ in 0..next() % 8 {
            rt += (next() % 3) as f64 / 10.;
            let mut mz = 0.;
            let mut points = Vec::new();
            for i in 0..next() % 7 {
                mz += (next() % 3) as f64 / 10.;
                points.push((mz, i as f32 - 2.));
            }
            exp.spectra.push(scan(rt, (next() % 3) as u32, &points));
        }
        for _ in 0..8 {
            let lo = (next() % 5) as f64 / 10.;
            let hi = lo + (next() % 8) as f64 / 10.;
            let level = next() % 3;
            let b = bounds((0., 1.), (lo, hi));
            let mut expected = FlatPeakData::default();
            let mut expected_rows = SpectrumPeakData::default();
            let mut remembered = -1f32;
            for s in &exp.spectra {
                if s.ms_level as usize != level || !(0. ..=1.).contains(&s.rt) {
                    continue;
                }
                for p in &s.peaks {
                    if p.mz < lo || p.mz > hi {
                        continue;
                    }
                    expected.rt.push(s.rt as f32);
                    expected.mz.push(p.mz as f32);
                    expected.intensity.push(p.intensity);
                    if s.rt != f64::from(remembered) {
                        remembered = s.rt as f32;
                        expected_rows.rt.push(remembered);
                        expected_rows.mz.push(Vec::new());
                        expected_rows.intensity.push(Vec::new());
                    }
                    expected_rows.mz.last_mut().unwrap().push(p.mz as f32);
                    expected_rows
                        .intensity
                        .last_mut()
                        .unwrap()
                        .push(p.intensity);
                }
            }
            assert_eq!(exp.get_2d_peak_data(b, level).unwrap(), expected);
            assert_eq!(
                exp.get_2d_peak_data_per_spectrum(b, level).unwrap(),
                expected_rows
            );
        }
    }
}
