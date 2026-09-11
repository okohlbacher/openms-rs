use openms::Error;
use openms::kernel::{
    AreaBounds, AreaIter, AreaLimits, AreaOptions, DataArray, MSExperiment, MSSpectrum, MzRtRegion,
    NumericRange, Peak1D,
};

fn scan(rt: f64, level: u32, mzs: &[f64]) -> MSSpectrum {
    MSSpectrum {
        rt,
        ms_level: level,
        peaks: mzs.iter().map(|&mz| Peak1D::new(mz, mz as f32)).collect(),
        ..Default::default()
    }
}
fn source() -> MSExperiment {
    // AreaIterator_test.cpp literal five-scan fixture, including the empty scan.
    MSExperiment {
        spectra: vec![
            scan(2., 1, &[502., 510.]),
            scan(4., 1, &[504., 506.]),
            scan(6., 1, &[]),
            scan(8., 1, &[504.1, 506.1]),
            scan(10., 1, &[502.1, 510.1]),
        ],
        ..Default::default()
    }
}
fn options(rt: Option<(f64, f64)>, mz: Option<(f64, f64)>, level: u32) -> AreaOptions {
    AreaOptions::new(
        AreaBounds {
            rt: rt.map(|(min, max)| NumericRange { min, max }),
            mz: mz.map(|(min, max)| NumericRange { min, max }),
        },
        level,
    )
}
fn positions(exp: &MSExperiment, o: AreaOptions) -> Vec<(usize, usize)> {
    exp.area_iter(o)
        .unwrap()
        .map(|p| (p.spectrum_index, p.peak_index))
        .collect()
}

#[test]
fn source_area_fixture_quadrants_indices_rt_and_empty_scan() {
    let exp = source();
    let cases = [
        (
            (0., 15., 500., 520.),
            vec![502., 510., 504., 506., 504.1, 506.1, 502.1, 510.1],
        ),
        ((3., 9., 503., 509.), vec![504., 506., 504.1, 506.1]),
        ((0., 7., 505., 520.), vec![510., 506.]),
        ((5., 11., 505., 520.), vec![506.1, 510.1]),
        ((5., 11., 500., 505.), vec![504.1, 502.1]),
        ((0., 7., 500., 505.), vec![502., 504.]),
        ((5., 5.5, 500., 520.), vec![]),
        ((0., 15., 505., 505.5), vec![]),
    ];
    for ((rt0, rt1, mz0, mz1), expected) in cases {
        let observed: Vec<_> = exp
            .area_begin(rt0, rt1, mz0, mz1, 1)
            .unwrap()
            .map(|p| p.peak.mz)
            .collect();
        assert_eq!(observed, expected);
    }
    let points: Vec<_> = exp.area_iter(AreaOptions::default()).unwrap().collect();
    assert_eq!(
        points
            .iter()
            .map(|p| (p.spectrum_index, p.peak_index))
            .collect::<Vec<_>>(),
        [
            (0, 0),
            (0, 1),
            (1, 0),
            (1, 1),
            (3, 0),
            (3, 1),
            (4, 0),
            (4, 1)
        ]
    );
    assert_eq!(
        points.iter().map(|p| p.spectrum.rt).collect::<Vec<_>>(),
        [2., 2., 4., 4., 8., 8., 10., 10.]
    );
    for p in points {
        assert!(std::ptr::eq(
            p.peak,
            &exp.spectra[p.spectrum_index].peaks[p.peak_index]
        ));
        assert!(std::ptr::eq(p.spectrum, &exp.spectra[p.spectrum_index]));
    }
}

#[test]
fn source_experiment_fixture_and_mutable_endpoint_semantics() {
    let mut exp = MSExperiment {
        spectra: vec![
            scan(-1., 1, &[2.]),
            scan(1., 1, &[2., 3.]),
            scan(2., 1, &[10., 11., 12.]),
        ],
        ..Default::default()
    };
    let mut it = exp.area_begin_mut(0., 15., 0., 11., 1).unwrap();
    assert_eq!(it.len(), 4);
    let first = it.next().unwrap();
    assert_eq!(first.peak.mz, 2.);
    first.peak.mz = 4711.;
    // Earlier coordinates can change; the selected scan endpoint stays fixed.
    let second = it.next().unwrap();
    assert_eq!(second.peak.mz, 3.);
    second.peak.intensity = -3.;
    let rest: Vec<_> = it.map(|p| p.peak.mz).collect();
    assert_eq!(rest, [10., 11.]);
    assert_eq!(first.peak.mz, 4711.);
    assert_eq!(second.peak.intensity, -3.);
    assert_eq!(exp.spectra[2].peaks[2].mz, 12.);
    assert_eq!(exp.spectra[1].peaks[0].mz, 4711.);
    assert!(matches!(
        exp.area_iter(AreaOptions::default()),
        Err(Error::UnsortedData)
    ));
}

#[test]
fn inclusive_duplicates_zero_width_and_optional_dimensions() {
    let exp = MSExperiment {
        spectra: vec![
            scan(0., 1, &[1., 2., 2., 3.]),
            scan(0., 1, &[2.]),
            scan(1., 2, &[2.]),
            scan(2., 1, &[2.]),
        ],
        ..Default::default()
    };
    assert_eq!(
        positions(&exp, options(Some((0., 0.)), Some((2., 2.)), 1)),
        [(0, 1), (0, 2), (1, 0)]
    );
    assert_eq!(
        positions(&exp, options(None, Some((2., 2.)), 1)),
        [(0, 1), (0, 2), (1, 0), (3, 0)]
    );
    assert_eq!(positions(&exp, options(Some((0., 0.)), None, 1)).len(), 5);
    let region = MzRtRegion::new(2., 2., 0., 0.).unwrap();
    assert_eq!(
        positions(&exp, AreaOptions::new(region.into(), 1)),
        [(0, 1), (0, 2), (1, 0)]
    );
    let extremes = MSExperiment {
        spectra: vec![
            scan(-f64::MAX, 1, &[-f64::MAX, -0., 0., f64::MAX]),
            scan(f64::MAX, 1, &[f64::MAX]),
        ],
        ..Default::default()
    };
    assert_eq!(
        extremes.area_iter(AreaOptions::default()).unwrap().count(),
        5
    );
    assert_eq!(
        positions(&extremes, options(None, Some((0., 0.)), 1)),
        [(0, 1), (0, 2)]
    );
}

#[test]
fn source_level_conversion_is_explicit_and_native_levels_are_exact() {
    let exp = MSExperiment {
        spectra: vec![
            scan(0., 0, &[1.]),
            scan(1., 1, &[1.]),
            scan(2., 127, &[1.]),
            scan(3., 128, &[1.]),
            scan(4., 256, &[1.]),
            scan(5., u32::MAX, &[1.]),
            scan(6., u32::MAX - 127, &[1.]),
        ],
        ..Default::default()
    };
    assert_eq!(
        positions(&exp, AreaOptions::new(AreaBounds::default(), 0)),
        [(0, 0)]
    );
    assert_eq!(
        positions(&exp, AreaOptions::new(AreaBounds::default(), 128)),
        [(3, 0)]
    );
    assert_eq!(
        positions(&exp, AreaOptions::new(AreaBounds::default(), 256)),
        [(4, 0)]
    );
    assert_eq!(
        positions(
            &exp,
            AreaOptions::source_compatible(AreaBounds::default(), 256)
        ),
        [(0, 0)]
    );
    assert_eq!(
        positions(
            &exp,
            AreaOptions::source_compatible(AreaBounds::default(), 255)
        ),
        [(5, 0)]
    );
    assert_eq!(
        positions(
            &exp,
            AreaOptions::source_compatible(AreaBounds::default(), 128)
        ),
        [(6, 0)]
    );
    assert_eq!(
        exp.area_begin(-1., 10., 0., 2., 255)
            .unwrap()
            .next()
            .unwrap()
            .spectrum_index,
        5
    );
    assert!(positions(&exp, AreaOptions::new(AreaBounds::default(), 2)).is_empty());
}

#[test]
fn const_clones_have_independent_cursors_address_equality_and_fused_end() {
    let exp = source();
    let mut first = exp.area_iter(AreaOptions::default()).unwrap();
    let second = first.clone();
    assert_eq!(first, second);
    assert_eq!(first.peek().unwrap().peak.mz, 502.);
    assert_eq!(first.size_hint(), (8, Some(8)));
    first.next();
    assert_ne!(first, second);
    assert_eq!(first.len(), 7);
    let mut same_peak = exp.area_begin(0., 3., 509., 511., 1).unwrap();
    assert_eq!(same_peak, first);
    same_peak.next();
    assert_ne!(same_peak, first);
    let default = AreaIter::default();
    assert_eq!(same_peak, default);
    assert!(same_peak.peek().is_none());
    assert!(same_peak.next().is_none());
    assert!(same_peak.next().is_none());
    let other = source();
    assert_ne!(second, other.area_iter(AreaOptions::default()).unwrap());
    assert_eq!(second.clone().count(), 8);
    assert_eq!(second.len(), 8);
}

#[test]
fn mutable_rows_can_coexist_without_aliases_and_preserve_annotations() {
    let mut exp = source();
    exp.spectra[0]
        .string_data_arrays
        .push(DataArray::new("identity", vec!["a".into(), "b".into()]));
    let before = exp.spectra[0].string_data_arrays.clone();
    let mut iter = exp.area_iter_mut(options(Some((0., 4.)), None, 1)).unwrap();
    let a = iter.next().unwrap();
    let b = iter.next().unwrap();
    assert_ne!(a.peak as *mut _, b.peak as *mut _);
    a.peak.intensity = 11.;
    b.peak.intensity = 12.;
    let c = iter.next().unwrap();
    assert_eq!(
        (c.spectrum_index, c.peak_index, c.rt, c.ms_level),
        (1, 0, 4., 1)
    );
    c.peak.intensity = 13.;
    assert_eq!(iter.len(), 1);
    iter.next();
    assert_eq!(iter.len(), 0);
    assert!(iter.next().is_none());
    assert!(iter.next().is_none());
    drop(iter);
    assert_eq!(a.peak.intensity, 11.);
    assert_eq!(b.peak.intensity, 12.);
    assert_eq!(c.peak.intensity, 13.);
    assert_eq!(exp.spectra[0].string_data_arrays, before);
}

#[test]
fn invalid_boundaries_and_global_order_fail_before_mutable_borrow() {
    let mut exp = source();
    let before = exp.clone();
    for b in [
        AreaBounds {
            rt: Some(NumericRange { min: 2., max: 1. }),
            mz: None,
        },
        AreaBounds {
            rt: None,
            mz: Some(NumericRange {
                min: f64::NAN,
                max: 1.,
            }),
        },
        AreaBounds {
            rt: None,
            mz: Some(NumericRange {
                min: 0.,
                max: f64::INFINITY,
            }),
        },
    ] {
        assert!(exp.area_iter_mut(AreaOptions::new(b, 1)).is_err());
        assert_eq!(exp, before);
    }
    assert!(AreaBounds::new(1., 0., 0., 1.).is_err());
    assert!(AreaBounds::new(0., 1., 1., 0.).is_err());
    exp.spectra.swap(0, 1);
    assert!(matches!(
        exp.area_begin(100., 200., 0., 1., 3),
        Err(Error::UnsortedData)
    ));
    exp = before.clone();
    exp.spectra.last_mut().unwrap().peaks.swap(0, 1);
    assert!(matches!(
        exp.area_begin(0., 3., 0., 1., 3),
        Err(Error::UnsortedData)
    )); // even wrong-level/outside-area scans
    exp = before;
    exp.spectra[4].peaks[0].mz = f64::INFINITY;
    assert!(exp.area_iter(AreaOptions::default()).is_err());
    exp = source();
    exp.spectra[4].rt = f64::NAN;
    assert!(exp.area_iter(AreaOptions::default()).is_err());
}

#[test]
fn traversal_ignores_intensity_array_alignment_and_unrelated_records() {
    let mut exp = source();
    exp.spectra[0].peaks[0].intensity = f32::NAN;
    exp.spectra[0]
        .integer_data_arrays
        .push(DataArray::new("unused", vec![1]));
    exp.chromatograms
        .push(openms::MSChromatogram::from_peaks(vec![
            openms::ChromatogramPeak::new(f64::NAN, f32::INFINITY),
        ]));
    assert_eq!(exp.area_iter(AreaOptions::default()).unwrap().count(), 8);
    assert!(
        exp.area_iter(AreaOptions::default())
            .unwrap()
            .next()
            .unwrap()
            .peak
            .intensity
            .is_nan()
    );
}

#[test]
fn checked_resources_are_shared_and_empty_input_can_use_zero_limits() {
    let mut exp = source();
    let before = exp.clone();
    for limits in [
        AreaLimits {
            max_spectra: 4,
            ..Default::default()
        },
        AreaLimits {
            max_peaks: 7,
            ..Default::default()
        },
        AreaLimits {
            max_work: 12,
            ..Default::default()
        },
        // Validation and the interval plan succeed, but the final shared
        // traversal reservation does not. No mutable reference escapes.
        AreaLimits {
            max_work: 30,
            ..Default::default()
        },
        AreaLimits {
            max_bytes: 0,
            ..Default::default()
        },
    ] {
        assert!(
            matches!(exp.area_iter_mut_with_limits(AreaOptions::default(),limits),Err(Error::InvalidValue(s)) if s.contains("resource limit"))
        );
        assert_eq!(exp, before);
    }
    // Five validation visits + eight peak visits + five interval visits +
    // five future scan visits + eight future peak visits share one budget.
    assert_eq!(
        exp.area_iter_with_limits(
            AreaOptions::default(),
            AreaLimits {
                max_work: 31,
                ..Default::default()
            }
        )
        .unwrap()
        .count(),
        8
    );
    let empty = MSExperiment::new();
    let zero = AreaLimits {
        max_spectra: 0,
        max_peaks: 0,
        max_work: 0,
        max_bytes: 0,
    };
    assert_eq!(
        empty
            .area_iter_with_limits(AreaOptions::default(), zero)
            .unwrap()
            .count(),
        0
    );
    // Raw validation cannot be bypassed by choosing an empty region/MS level.
    assert!(
        exp.area_iter_with_limits(
            options(Some((100., 200.)), None, 9),
            AreaLimits {
                max_peaks: 7,
                ..Default::default()
            }
        )
        .is_err()
    );
}

#[test]
fn small_random_grids_match_independent_brute_enumeration_for_const_and_mutable() {
    let mut seed = 0x52a2_32fe_u64;
    let mut next = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (seed >> 32) as usize
    };
    for _ in 0..100 {
        let mut exp = MSExperiment::new();
        let mut rt = -3.;
        for _ in 0..next() % 8 {
            rt += (next() % 3) as f64;
            let mut mz = -2.;
            let mut peaks = Vec::new();
            for i in 0..next() % 7 {
                mz += (next() % 3) as f64;
                peaks.push(Peak1D::new(mz, i as f32));
            }
            exp.spectra.push(MSSpectrum {
                rt,
                ms_level: (next() % 4) as u32,
                peaks,
                ..Default::default()
            });
        }
        for _ in 0..8 {
            let rt0 = (next() % 12) as f64 - 4.;
            let rt1 = rt0 + (next() % 6) as f64;
            let mz0 = (next() % 12) as f64 - 3.;
            let mz1 = mz0 + (next() % 6) as f64;
            let level = (next() % 4) as u32;
            let bounds = AreaBounds::new(rt0, rt1, mz0, mz1).unwrap();
            let o = AreaOptions::new(bounds, level);
            let mut expected = Vec::new();
            for (si, s) in exp.spectra.iter().enumerate() {
                if s.ms_level == level && s.rt >= rt0 && s.rt <= rt1 {
                    for (pi, p) in s.peaks.iter().enumerate() {
                        if p.mz >= mz0 && p.mz <= mz1 {
                            expected.push((si, pi));
                        }
                    }
                }
            }
            assert_eq!(positions(&exp, o), expected);
            let mut edited = exp.clone();
            let actual: Vec<_> = edited
                .area_iter_mut(o)
                .unwrap()
                .map(|p| {
                    p.peak.intensity = -99.;
                    (p.spectrum_index, p.peak_index)
                })
                .collect();
            assert_eq!(actual, expected);
            for (si, s) in edited.spectra.iter().enumerate() {
                for (pi, p) in s.peaks.iter().enumerate() {
                    assert_eq!(
                        p.intensity,
                        if expected.contains(&(si, pi)) {
                            -99.
                        } else {
                            exp.spectra[si].peaks[pi].intensity
                        }
                    );
                }
            }
        }
    }
}
