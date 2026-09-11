// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::kernel::{MassTrace, MassTraceLimits, MassTraceQuantMethod as Quant, Peak2D, Point2D};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17} != {expected:.17}"
    );
}
fn trace(intensities: &[f32]) -> MassTrace {
    MassTrace::from_peaks(
        intensities
            .iter()
            .enumerate()
            .map(|(i, &v)| Peak2D::new(i as f64, 100. + i as f64, v))
            .collect(),
    )
    .unwrap()
}
fn source_trace() -> MassTrace {
    // Source input is f64, then narrowed by fillPeak to the f32 peak field.
    let rows: [(f64, f64, f64); 7] = [
        (152.22, 230.10223, 542.),
        (153.23, 230.10235, 542293.),
        (154.21, 230.10181, 18282393.),
        (155.24, 230.10229, 33329535.),
        (156.233, 230.10116, 17342933.),
        (157.24, 230.10198, 333291.),
        (158.238, 230.10254, 339.),
    ];
    MassTrace::from_peaks(
        rows.into_iter()
            .map(|(r, m, i)| Peak2D::new(r, m, i as f32))
            .collect(),
    )
    .unwrap()
}
const SMOOTH: [f64; 7] = [
    500., 540000., 18000000., 33000000., 17500000., 540000., 549223.,
];

#[test]
fn constructors_defaults_names_borrowed_access_and_complete_cached_clone() {
    let empty = MassTrace::new();
    assert!(empty.is_empty());
    assert_eq!(empty.len(), 0);
    assert_eq!(empty.label(), "");
    assert_eq!(empty.centroid_mz(), 0.);
    assert_eq!(empty.centroid_rt(), 0.);
    assert_eq!(empty.centroid_sd(), 0.);
    assert_eq!(empty.centroid_im(), 0.);
    assert!(!empty.contains_im_data());
    assert_eq!(empty.fwhm(), 0.);
    assert_eq!(empty.fwhm_borders(), (0, 0));
    assert_eq!(empty.quant_method(), Quant::Area);
    for (method, name) in Quant::ALL.into_iter().zip(Quant::NAMES) {
        assert_eq!(Quant::from_name(name), Some(method));
        assert_eq!(method.name(), name);
    }
    assert_eq!(Quant::from_name("somethingwrong"), None);
    assert_eq!(Quant::from_name("Area"), None);
    let mut t = source_trace();
    let copied = MassTrace::from_slice(t.peaks()).unwrap();
    assert_eq!(copied, t);
    assert_ne!(copied.peaks().as_ptr(), t.peaks().as_ptr());
    assert_eq!(t[1].rt(), 153.23);
    assert_eq!(t.get(99), None);
    assert_eq!(t.iter().next_back().unwrap().rt(), 158.238);
    t.get_mut(0).unwrap().set_mz(100.);
    t[0].set_mz(230.10223);
    for p in &mut t {
        p.intensity += 0.;
    }
    assert_eq!((&t).into_iter().count(), 7);
    t.set_label("TEST_TRACE").unwrap();
    t.set_centroid_sd(-2.).unwrap();
    t.set_centroid_im(-0.).unwrap();
    t.fwhm_mz_avg = 1.;
    t.fwhm_im_avg = 2.;
    t.set_smoothed_intensities(&SMOOTH).unwrap();
    t.estimate_fwhm(true).unwrap();
    t.update_weighted_mean_mz().unwrap();
    t.update_weighted_mean_rt().unwrap();
    let mut copy = t.clone();
    assert_eq!(copy, t);
    assert!(copy.contains_im_data());
    assert_eq!(copy.centroid_im().to_bits(), (-0f64).to_bits());
    copy.peaks_mut()[0].intensity = 123.;
    assert_ne!(copy[0], t[0]);
    assert_eq!(copy.centroid_mz(), t.centroid_mz()); // no implicit cache refresh
    assert_eq!(copy.smoothed_intensities(), SMOOTH);
}

#[test]
fn source_seven_peak_area_width_height_and_centroid_literals() {
    let mut t = source_trace();
    close(t.trace_length().unwrap(), 6.018, 1e-12);
    close(t.average_ms1_cycle_time().unwrap(), 6.018 / 6., 1e-12);
    close(t.compute_peak_area().unwrap(), 70303710.2575001, 1e-6);
    assert_eq!(t.compute_intensity_sum().unwrap(), 69831325.); // exact f32 source inputs
    t.update_weighted_mean_rt().unwrap();
    close(t.centroid_rt(), 155.214671250425, 1e-12);
    t.update_weighted_mean_mz().unwrap();
    close(t.centroid_mz(), 230.10188267861923, 1e-12);
    t.update_weighted_mz_sd().unwrap();
    close(t.centroid_sd(), 0.00046045716004770155, 1e-14);
    t.update_median_rt().unwrap();
    assert_eq!(t.centroid_rt(), 155.24);
    t.update_median_mz().unwrap();
    assert_eq!(t.centroid_mz(), 230.10223); // fourth of the seven sorted literal m/z values
    t.update_mean_mz().unwrap();
    close(t.centroid_mz(), 230.10205142857143, 1e-12);
    assert_eq!(t.find_max_by_int_peak(false).unwrap(), 3);
    assert_eq!(t.max_intensity(false).unwrap(), 33329536.);
    close(t.estimate_fwhm(false).unwrap(), 2.15250939241199, 1e-12);
    assert_eq!(t.fwhm_borders(), (1, 5));
    close(t.compute_fwhm_area().unwrap(), 69863097.2125001, 1e-6);
    t.set_smoothed_intensities(&SMOOTH).unwrap();
    close(
        t.compute_smoothed_peak_area().unwrap(),
        70303689.0475001,
        1e-6,
    );
    close(t.estimate_fwhm(true).unwrap(), 2.16656743986252, 1e-12);
    close(t.intensity(true).unwrap(), 69505990.0000001, 1e-6);
    close(t.intensity(false).unwrap(), 69863097.2125001, 1e-6);
    assert_eq!(t.find_max_by_int_peak(true).unwrap(), 3);
    assert_eq!(t.max_intensity(true).unwrap(), 33000000.);
    t.update_smoothed_max_rt().unwrap();
    assert_eq!(t.centroid_rt(), 155.24);
    t.update_smoothed_weighted_mean_rt().unwrap();
    close(t.centroid_rt(), 155.24680397032225, 1e-12);
}

#[test]
fn source_symmetric_asymmetric_and_open_flank_crossings() {
    for (values, width) in [
        ([10., 60., 100., 60., 10.], 2.4),
        ([10., 60., 100., 80., 40.], 2.95),
        ([10., 60., 100., 100., 100.], 3.2),
    ] {
        let mut t = trace(&values);
        close(t.estimate_fwhm(false).unwrap(), width, 1e-14);
        assert_eq!(t.fwhm_borders(), (0, 4));
        t.set_smoothed_intensities(&values.map(f64::from)).unwrap();
        close(t.estimate_fwhm(true).unwrap(), width, 1e-14);
    }
    let mut duplicates = trace(&[10., 60., 100., 60., 10.]);
    for p in duplicates.peaks_mut() {
        p.set_rt(1.);
    }
    assert_eq!(duplicates.estimate_fwhm(false).unwrap(), 0.);
    assert_eq!(duplicates.fwhm_borders(), (0, 4));
}

#[test]
fn endpoint_apex_clears_only_borders_and_failed_updates_leave_state_unchanged() {
    let mut t = trace(&[10., 60., 100., 60., 10.]);
    let width = t.estimate_fwhm(false).unwrap();
    t[0].intensity = 200.;
    assert_eq!(t.estimate_fwhm(false).unwrap(), 0.);
    assert_eq!(t.fwhm_borders(), (0, 0));
    assert_eq!(t.fwhm(), width);
    assert_eq!(t.intensity(false).unwrap(), 0.);
    t[0].intensity = 10.;
    t.estimate_fwhm(false).unwrap();
    t[1].set_rt(-1.);
    let old = t.clone();
    assert!(t.estimate_fwhm(false).is_err());
    assert_eq!(t, old);
    t[1].set_rt(1.);
    t[4].set_rt(f64::MAX);
    t[0].set_rt(-f64::MAX);
    let borders = t.fwhm_borders();
    let width = t.fwhm();
    assert!(t.estimate_fwhm(false).is_err());
    assert_eq!(t.fwhm_borders(), borders);
    assert_eq!(t.fwhm(), width);
}

#[test]
fn smoothed_area_mixes_raw_values_while_smoothed_centroids_use_positive_weights() {
    let mut t = MassTrace::from_peaks(vec![
        Peak2D::new(0., 1., 2.),
        Peak2D::new(2., 2., 6.),
        Peak2D::new(5., 3., 10.),
        Peak2D::new(6., 4., 14.),
    ])
    .unwrap();
    t.set_smoothed_intensities(&[100., 2., 4., 0.]).unwrap();
    assert_eq!(t.compute_smoothed_peak_area().unwrap(), 130.);
    t.set_smoothed_intensities(&[100., -3., 4., -1.]).unwrap();
    assert_eq!(t.compute_smoothed_peak_area().unwrap(), 24.);
    t.set_smoothed_intensities(&[-100., 2., 4., -1.]).unwrap();
    t.update_smoothed_weighted_mean_rt().unwrap();
    assert_eq!(t.centroid_rt(), 4.);
    t.update_smoothed_max_rt().unwrap();
    assert_eq!(t.centroid_rt(), 5.);
    t.update_weighted_mean_rt().unwrap();
    close(t.centroid_rt(), 258. / 56., 1e-14);
    t[0].intensity = f32::NAN; // right-rectangle centroid does not use first intensity
    t.update_weighted_mean_rt().unwrap();
    close(t.centroid_rt(), 258. / 56., 1e-14);
}

#[test]
fn empty_unset_negative_and_single_peak_rules_are_distinct() {
    let mut e = MassTrace::new();
    assert_eq!(e.compute_peak_area().unwrap(), 0.);
    assert_eq!(e.compute_intensity_sum().unwrap(), 0.);
    assert_eq!(e.trace_length().unwrap(), 0.);
    assert_eq!(e.average_ms1_cycle_time().unwrap(), 0.);
    assert_eq!(e.intensity(false).unwrap(), 0.);
    assert_eq!(e.intensity(true).unwrap(), 0.);
    assert_eq!(e.max_intensity(true).unwrap(), 0.);
    assert_eq!(e.max_intensity(false).unwrap(), 0.);
    assert!(e.compute_smoothed_peak_area().is_err());
    assert!(e.find_max_by_int_peak(false).is_err());
    assert!(e.find_max_by_int_peak(true).is_err());
    assert!(e.estimate_fwhm(false).is_err());
    assert!(e.update_mean_mz().is_err());
    assert!(e.update_median_mz().is_err());
    assert!(e.update_median_rt().is_err());
    assert!(e.update_weighted_mean_mz().is_err());
    assert!(e.update_weighted_mean_rt().is_err());
    assert!(e.update_weighted_mz_sd().is_err());
    assert!(e.update_smoothed_max_rt().is_err());
    assert!(e.update_smoothed_weighted_mean_rt().is_err());
    e.set_quant_method(Quant::Median);
    assert!(e.intensity(true).is_err());
    let mut one = MassTrace::from_peaks(vec![Peak2D::new(-0., -0., f32::NAN)]).unwrap();
    for update in [
        MassTrace::update_mean_mz,
        MassTrace::update_median_mz,
        MassTrace::update_weighted_mean_mz,
    ] {
        update(&mut one).unwrap();
        assert_eq!(one.centroid_mz().to_bits(), (-0f64).to_bits());
    }
    one.update_weighted_mean_rt().unwrap();
    assert_eq!(one.centroid_rt().to_bits(), (-0f64).to_bits());
    one.set_smoothed_intensities(&[-1.]).unwrap();
    one.update_smoothed_weighted_mean_rt().unwrap();
    one.update_smoothed_max_rt().unwrap();
    assert_eq!(one.centroid_rt().to_bits(), (-0f64).to_bits());
    let mut neg = trace(&[-4., -2., -3.]);
    assert_eq!(neg.find_max_by_int_peak(false).unwrap(), 1);
    assert_eq!(neg.max_intensity(false).unwrap(), 0.);
    assert!(neg.estimate_fwhm(false).is_err());
    neg.set_smoothed_intensities(&[-1., 0., -2.]).unwrap();
    assert!(neg.update_smoothed_max_rt().is_err());
    assert!(neg.update_smoothed_weighted_mean_rt().is_err());
}

#[test]
fn quantification_median_ignores_smoothing_and_preserves_raw_precision() {
    let mut t = source_trace();
    t.set_quant_method(Quant::Median);
    assert_eq!(t.intensity(false).unwrap(), 542293.);
    assert_eq!(t.intensity(true).unwrap(), 542293.);
    t.set_smoothed_intensities(&[1.; 7]).unwrap();
    assert_eq!(t.intensity(true).unwrap(), 542293.);
    t.set_quant_method(Quant::MaxHeight);
    assert_eq!(t.intensity(false).unwrap(), 33329536.);
    assert_eq!(t.intensity(true).unwrap(), 1.);
    let mut even = trace(&[9., 1., 7., 3.]);
    even.set_quant_method(Quant::Median);
    assert_eq!(even.intensity(true).unwrap(), 5.);
    even.update_median_rt().unwrap();
    assert_eq!(even.centroid_rt(), 1.5);
    even.update_median_mz().unwrap();
    assert_eq!(even.centroid_mz(), 101.5);
}

#[test]
fn source_missing_scan_area_invariance_and_signed_unsorted_integrals() {
    let rows = [
        542., 542293., 18282400., 33329500., 33329500., 33329500., 17342900., 333291., 339.,
    ];
    let mut full = trace(&rows);
    let mut reduced = MassTrace::from_peaks(
        full.peaks()
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != 4)
            .map(|(_, p)| *p)
            .collect(),
    )
    .unwrap();
    assert_eq!(
        full.compute_peak_area().unwrap(),
        reduced.compute_peak_area().unwrap()
    );
    full.estimate_fwhm(false).unwrap();
    reduced.estimate_fwhm(false).unwrap();
    assert_eq!(
        full.compute_fwhm_area().unwrap(),
        reduced.compute_fwhm_area().unwrap()
    );
    let reverse =
        MassTrace::from_peaks(vec![Peak2D::new(4., 1., 2.), Peak2D::new(1., 2., 4.)]).unwrap();
    assert_eq!(reverse.trace_length().unwrap(), 3.);
    assert_eq!(reverse.average_ms1_cycle_time().unwrap(), -3.);
    assert_eq!(reverse.compute_peak_area().unwrap(), -9.);
}

#[test]
fn hull_uses_rt_mz_envelope_and_ignores_intensity() {
    let mut t = source_trace();
    t[0].intensity = f32::NAN;
    let h = t.convex_hull().unwrap();
    assert!(h.encloses(Point2D::new(154.21, 230.10181)).unwrap());
    assert!(!h.encloses(Point2D::new(155.22, 230.10181)).unwrap());
    assert!(!h.encloses(Point2D::new(154.21, 229.10181)).unwrap());
    assert!(MassTrace::new().convex_hull().unwrap().is_empty());
    let h = MassTrace::from_peaks(vec![
        Peak2D::new(1., 3., 0.),
        Peak2D::new(1., 1., 0.),
        Peak2D::new(0., 2., 0.),
    ])
    .unwrap()
    .convex_hull()
    .unwrap();
    assert!(h.encloses(Point2D::new(1., 2.)).unwrap());
}

#[test]
fn setters_resource_and_late_numeric_errors_are_atomic() {
    let mut t = trace(&[1., 4., 1.]);
    t.set_label("old").unwrap();
    t.set_smoothed_intensities(&[2., 5., 2.]).unwrap();
    assert!(t.set_smoothed_intensities(&[1.]).is_err());
    assert!(t.set_smoothed_intensities(&[1., 2., f64::NAN]).is_err());
    assert_eq!(t.smoothed_intensities(), [2., 5., 2.]);
    assert!(t.set_centroid_im(f64::INFINITY).is_err());
    assert!(!t.contains_im_data());
    t.set_centroid_sd(2.).unwrap();
    assert!(t.set_centroid_sd(f64::NAN).is_err());
    assert_eq!(t.centroid_sd(), 2.);
    t.limits.max_bytes = 2;
    assert!(t.set_label("large").is_err());
    assert_eq!(t.label(), "old");
    assert!(t.set_smoothed_intensities(&[1., 2., 3.]).is_err());
    assert!(t.update_median_mz().is_err());
    assert!(t.convex_hull().is_err());
    t.limits = MassTraceLimits::default();
    t.update_mean_mz().unwrap();
    let old = t.centroid_mz();
    t[2].set_mz(f64::INFINITY);
    assert!(t.update_mean_mz().is_err());
    assert_eq!(t.centroid_mz(), old);
    t.limits.max_work = 0;
    let error = t.update_mean_mz().unwrap_err().to_string();
    assert!(error.contains("resource"));
    let limits = MassTraceLimits {
        max_peaks: 2,
        ..Default::default()
    };
    assert!(MassTrace::from_slice_with_limits(t.peaks(), limits).is_err());
    assert!(MassTrace::from_peaks_with_limits(t.peaks().to_vec(), limits).is_err());
    let limits = MassTraceLimits {
        max_bytes: 0,
        ..Default::default()
    };
    assert!(MassTrace::from_slice_with_limits(&[Peak2D::default()], limits).is_err());
    assert!(MassTrace::from_peaks_with_limits(vec![Peak2D::default()], limits).is_ok()); // ownership needs no new peak allocation
}

#[test]
fn weighted_epsilon_zero_variance_and_nonfinite_arithmetic_boundaries() {
    let mut t = trace(&[0., 0.]);
    assert!(t.update_weighted_mean_mz().is_err());
    assert!(t.update_weighted_mz_sd().is_err());
    assert!(t.update_weighted_mean_rt().is_err());
    t[0].intensity = f64::EPSILON as f32;
    t.update_weighted_mean_mz().unwrap();
    assert_eq!(t.centroid_mz(), 100.);
    t.update_weighted_mz_sd().unwrap();
    assert_eq!(t.centroid_sd(), 0.); // exp(2*ln(0))
    t[0].intensity = 1.;
    t[1].intensity = -1.;
    assert!(t.update_weighted_mean_mz().is_err());
    let mut tiny = trace(&[0., f32::MIN_POSITIVE]);
    tiny.update_weighted_mean_rt().unwrap();
    assert_eq!(tiny.centroid_rt(), 1.); // RT has no epsilon guard
    tiny[0].set_rt(1.);
    assert!(tiny.update_weighted_mean_rt().is_err());
    let huge = MassTrace::from_peaks(vec![
        Peak2D::new(-f64::MAX, 1., 1.),
        Peak2D::new(f64::MAX, 2., 1.),
    ])
    .unwrap();
    assert!(huge.compute_peak_area().is_err());
    assert!(huge.trace_length().is_err());
    assert!(huge.average_ms1_cycle_time().is_err());
}

#[test]
fn deterministic_small_grids_match_independent_weighted_and_integral_oracles() {
    let mut seed = 71u64;
    let mut next = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        (seed >> 32) as u32
    };
    for _ in 0..100 {
        let n = (next() % 12 + 2) as usize;
        let mut rt = -5.;
        let mut peaks = Vec::new();
        for _ in 0..n {
            rt += f64::from(next() % 3 + 1);
            peaks.push(Peak2D::new(
                rt,
                100. + f64::from(next() % 100) / 10.,
                (next() % 20 + 1) as f32,
            ));
        }
        let mut t = MassTrace::from_peaks(peaks).unwrap();
        let p = t.peaks();
        let weight: f64 = p.iter().map(|p| f64::from(p.intensity)).sum();
        let expected_mz = 100.
            + p.iter()
                .map(|p| f64::from(p.intensity) * (p.mz() - 100.))
                .sum::<f64>()
                / weight;
        let area = p
            .windows(2)
            .map(|w| {
                (w[1].rt() - w[0].rt())
                    * 0.5
                    * (f64::from(w[0].intensity) + f64::from(w[1].intensity))
            })
            .sum::<f64>();
        let rt_weights = p
            .windows(2)
            .map(|w| f64::from(w[1].intensity) * (w[1].rt() - w[0].rt()))
            .collect::<Vec<_>>();
        let total = rt_weights.iter().sum::<f64>();
        let expected_rt = rt_weights
            .iter()
            .zip(&p[1..])
            .map(|(w, p)| w / total * p.rt())
            .sum::<f64>();
        let expected_sd = (p
            .iter()
            .map(|p| f64::from(p.intensity) / weight * (p.mz() - expected_mz).powi(2))
            .sum::<f64>())
        .sqrt();
        close(t.compute_peak_area().unwrap(), area, 1e-10);
        t.update_weighted_mean_rt().unwrap();
        close(t.centroid_rt(), expected_rt, 1e-12);
        t.update_weighted_mean_mz().unwrap();
        close(t.centroid_mz(), expected_mz, 1e-12);
        t.update_weighted_mz_sd().unwrap();
        close(t.centroid_sd(), expected_sd, 1e-12);
    }
}
