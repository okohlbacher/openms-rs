// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::analysis::mass_trace_detection::{
    MassTraceDetection, MassTraceDetectionOptions, TraceTerminationCriterion,
};
use openms::concept::progress_logger::ProgressLogType;
use openms::kernel::{
    AreaIter, AreaOptions, DataArray, MSExperiment, MSSpectrum, MassTrace, MassTraceQuantMethod,
    Peak1D, Peak2D,
};

fn detector() -> MassTraceDetection {
    let mut d = MassTraceDetection::new();
    d.logger.set_log_type(ProgressLogType::None);
    d.options.min_trace_length = 0.;
    d.options.reestimate_mt_sd = false;
    d
}
fn scan(rt: f64, mz_int: &[(f64, f32)]) -> MSSpectrum {
    MSSpectrum {
        rt,
        peaks: mz_int.iter().map(|&(m, i)| Peak1D::new(m, i)).collect(),
        ..Default::default()
    }
}
fn input(rows: &[&[(f64, f32)]]) -> MSExperiment {
    MSExperiment {
        spectra: rows
            .iter()
            .enumerate()
            .map(|(i, p)| scan(i as f64, p))
            .collect(),
        ..Default::default()
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17} != {b:.17}");
}
fn source_input() -> MSExperiment {
    let mut result = MSExperiment::default();
    for line in include_str!("data/mass_trace_detection_source.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let v: Vec<_> = line.split('\t').collect();
        let index: usize = v[0].parse().unwrap();
        if index == result.spectra.len() {
            result.spectra.push(MSSpectrum {
                rt: v[1].parse().unwrap(),
                ms_level: v[2].parse().unwrap(),
                ..Default::default()
            });
        }
        if v[3] != "." {
            result.spectra[index]
                .peaks
                .push(Peak1D::new(v[3].parse().unwrap(), v[4].parse().unwrap()));
        }
    }
    result
}
fn add_arrays(input: &mut MSExperiment, ims: &[f32]) {
    for (i, scan) in input.spectra.iter_mut().enumerate() {
        scan.float_data_arrays = vec![
            DataArray::new("FWHM_ppm", vec![i as f32 + 2.; scan.len()]),
            DataArray::new("Ion Mobility", vec![ims[i]; scan.len()]),
            DataArray::new("IM Peak FWHM", vec![i as f32 + 8.; scan.len()]),
        ];
    }
}

#[test]
fn source_defaults_and_exact_fixture_three_traces_and_ms2_invariance() {
    let mut d = MassTraceDetection::new();
    assert_eq!(d.logger.log_type(), ProgressLogType::Cmd);
    assert_eq!(d.options, MassTraceDetectionOptions::default());
    d.logger.set_log_type(ProgressLogType::None);
    let mut input = source_input();
    assert_eq!(input.spectra.len(), 133);
    assert_eq!(d.run(&input, 0).unwrap().len(), 2);
    d.options.min_trace_length = 3.;
    let result = d.run(&input, 0).unwrap();
    assert_eq!(
        result.iter().map(MassTrace::len).collect::<Vec<_>>(),
        [86, 31, 16]
    );
    for (i, t) in result.iter().enumerate() {
        close(t.centroid_rt(), [348.667, 347.107, 346.888][i], 0.0005);
        close(
            t.centroid_mz(),
            [437.26675, 438.27241, 439.27594][i],
            0.0001,
        );
        close(
            t.compute_peak_area().unwrap(),
            [3381.72226139326, 664.763828332733, 109.490108620676][i],
            1e-8,
        );
        assert_eq!(t.label(), format!("T{}", i + 1));
    }
    let mut ignored = scan(f64::NAN, &[(f64::NAN, f32::NAN)]);
    ignored.ms_level = 2;
    input.spectra.splice(0..0, vec![ignored; 133]);
    assert_eq!(d.run(&input, 0).unwrap(), result);
    assert_eq!(d.run(&input, 2).unwrap(), result[..2]);
}

#[test]
fn strict_noise_apex_boundaries_and_reverse_encounter_ties() {
    let input = input(&[
        &[(100., 10.), (200., 30.), (300., 31.)],
        &[],
        &[(400., 31.)],
    ]);
    let mut d = detector();
    let traces = d.run(&input, 0).unwrap();
    assert_eq!(
        traces.iter().map(|t| t.centroid_mz()).collect::<Vec<_>>(),
        [400., 300.]
    );
    d.options.chrom_peak_snr = 1.;
    assert_eq!(
        d.run(&input, 0)
            .unwrap()
            .iter()
            .map(|t| t.centroid_mz())
            .collect::<Vec<_>>(),
        [400., 300., 200.]
    );
    d.options.noise_threshold_int = 31.;
    assert!(d.run(&input, 0).unwrap().is_empty());
}

#[test]
fn im_gate_chooses_nearest_mz_and_double_apex_centroid_with_fwhm_medians() {
    let mut input = input(&[&[(100., 40.)], &[(100., 100.)], &[(100., 60.)]]);
    add_arrays(&mut input, &[1., 1.1, 1.2]);
    let mut d = detector();
    d.options.ion_mobility_tolerance = 0.3;
    let traces = d.run(&input, 0).unwrap();
    assert_eq!(traces.len(), 1);
    let t = &traces[0];
    assert_eq!(t.len(), 3);
    close(
        t.centroid_im(),
        (40. * f64::from(1f32) + 200. * f64::from(1.1f32) + 60. * f64::from(1.2f32)) / 300.,
        1e-14,
    );
    assert_eq!(t.fwhm_mz_avg, 3.);
    assert_eq!(t.fwhm_im_avg, 9.);
    assert!(d.has_centroid_im() && d.has_fwhm_mz() && d.has_fwhm_im());
    for quant in MassTraceQuantMethod::ALL {
        d.options.quant_method = quant;
        assert_eq!(d.run(&input, 1).unwrap()[0].quant_method(), quant);
    }
    // The first in-window point has closer IM; the second has closer m/z.
    input.spectra[0].peaks = vec![Peak1D::new(99.999, 40.), Peak1D::new(100., 40.)];
    input.spectra[0].float_data_arrays[0].data = vec![2., 20.];
    input.spectra[0].float_data_arrays[1].data = vec![1.1, 1.2];
    input.spectra[0].float_data_arrays[2].data = vec![8., 80.];
    let t = d.run(&input, 1).unwrap().remove(0);
    assert_eq!(t[0].mz(), 100.);
    assert_eq!(t.fwhm_mz_avg, 4.); // [20,3,4]
}

#[test]
fn metadata_first_nonempty_array_scan_coherence_and_reuse_are_atomic() {
    let mut input = input(&[&[(100., 40.)], &[(100., 100.)], &[(100., 60.)]]);
    add_arrays(&mut input, &[1., 1., 1.]);
    let mut d = detector();
    d.run(&input, 0).unwrap();
    input.spectra[2].float_data_arrays.swap(0, 1);
    let mut old = vec![MassTrace::from_peaks(vec![Peak2D::new(7., 8., 9.)]).unwrap()];
    let before = old.clone();
    assert!(d.run_into(&input, &mut old, 0).is_err());
    assert_eq!(old, before);
    assert!(d.has_centroid_im());
    for s in &mut input.spectra {
        s.float_data_arrays.clear();
    }
    d.run(&input, 0).unwrap();
    assert!(!d.has_centroid_im() && !d.has_fwhm_mz() && !d.has_fwhm_im());
    // An unrelated first nonempty array vector prevents discovering later IM.
    input.spectra[0].float_data_arrays = vec![DataArray::new("unrelated", vec![])];
    input.spectra[1].float_data_arrays = vec![DataArray::new("Ion Mobility", vec![1.])];
    assert_eq!(d.run(&input, 0).unwrap()[0].len(), 3);
    assert!(!d.has_centroid_im());
    input.spectra[0].float_data_arrays = vec![DataArray::new("Ion Mobility", vec![])];
    assert!(d.run(&input, 0).is_err());
}

#[test]
fn ccs_warning_is_first_recognized_array_only_and_not_scientific_im() {
    let mut input = input(&[&[(100., 40.)], &[(100., 100.)], &[(100., 60.)]]);
    input.spectra[0].float_data_arrays = vec![DataArray::new("Ion Mobility MS:1002954", vec![])];
    let mut d = detector();
    d.run(&input, 0).unwrap();
    assert!(d.ccs_tolerance_warning());
    assert!(!d.has_centroid_im());
    d.options.ion_mobility_tolerance = 1.;
    d.run(&input, 0).unwrap();
    assert!(!d.ccs_tolerance_warning());
    d.options.ion_mobility_tolerance = 0.01;
    for name in [
        "mean ion mobility drift time array",
        "mean ion mobility array",
        "raw ion mobility array",
        "raw ion mobility drift time array",
        "deconvoluted ion mobility array",
        "deconvoluted ion mobility drift time array",
        "mean inverse reduced ion mobility array",
        "raw inverse reduced ion mobility array",
        "deconvoluted inverse reduced ion mobility array",
        "inverse reduced ion mobility vendor",
        "Ion Mobility MS:1003006 MS:1002954",
    ] {
        input.spectra[0]
            .float_data_arrays
            .insert(0, DataArray::new(name, vec![]));
        d.run(&input, 0).unwrap();
        assert!(!d.ccs_tolerance_warning(), "{name}");
        input.spectra[0].float_data_arrays.remove(0);
    }
    input.spectra[0]
        .float_data_arrays
        .insert(0, DataArray::new("ion mobility array", vec![]));
    d.run(&input, 0).unwrap();
    assert!(d.ccs_tolerance_warning()); // generic parent is not its own descendant
}

#[test]
fn outlier_greater_than_threshold_empty_scans_and_length_filters() {
    let mut d = detector();
    d.options.trace_termination_outliers = 0;
    let empty = input(&[&[(100., 40.)], &[], &[(100., 100.)], &[], &[(100., 60.)]]);
    assert_eq!(d.run(&empty, 1).unwrap()[0].len(), 3); // empties do not increment consecutive misses
    let misses = input(&[
        &[(100., 40.)],
        &[(200., 20.)],
        &[(100., 100.)],
        &[(200., 20.)],
        &[(100., 60.)],
    ]);
    assert_eq!(d.run(&misses, 1).unwrap()[0].len(), 1);
    d.options.trace_termination_outliers = 1;
    assert_eq!(d.run(&misses, 1).unwrap()[0].len(), 3); // one miss is allowed
    d.options.min_trace_length = 4.;
    assert_eq!(d.run(&misses, 0).unwrap().len(), 1);
    d.options.min_trace_length = 4.1;
    assert!(d.run(&misses, 0).unwrap().is_empty());
    d.options.min_trace_length = 0.;
    d.options.max_trace_length = 3.;
    assert!(d.run(&misses, 0).unwrap().is_empty());
    d.options.max_trace_length = -3.;
    assert_eq!(d.run(&misses, 0).unwrap().len(), 1);
}

#[test]
fn sample_rate_requires_six_direction_scans_and_final_quality_keeps_empty_scans() {
    let mut d = detector();
    d.options.trace_termination_criterion = TraceTerminationCriterion::SampleRate;
    d.options.min_sample_rate = 0.5;
    let mut e = input(&[
        &[],
        &[],
        &[],
        &[],
        &[],
        &[],
        &[(100., 100.)],
        &[],
        &[],
        &[],
        &[],
        &[],
        &[],
    ]);
    // Source empty scans don't become trailing misses; they still reduce quality.
    assert!(d.run(&e, 0).unwrap().is_empty());
    // Nonempty misses ARE subtracted from the final denominator.
    for s in &mut e.spectra {
        if s.is_empty() {
            s.peaks.push(Peak1D::new(200., 20.));
        }
    }
    assert_eq!(d.run(&e, 0).unwrap()[0].len(), 1);
    // At six scanned positions the direction is stopped; the seventh is unreachable.
    e.spectra.insert(0, scan(-1., &[(100., 40.)]));
    e.spectra.push(scan(13., &[(100., 40.)]));
    assert_eq!(d.run(&e, 1).unwrap()[0].len(), 1);
    d.options.trace_termination_criterion = TraceTerminationCriterion::Outlier;
    d.options.trace_termination_outliers = 6;
    d.options.min_sample_rate = 0.1;
    assert_eq!(d.run(&e, 1).unwrap()[0].len(), 3);
}

#[test]
fn nearest_visited_candidate_does_not_fall_back_and_im_negative_tolerance_is_empty() {
    let mut e = input(&[
        &[(100., 100.), (100.001, 50.)],
        &[(100., 200.), (100.003, 20.)],
        &[(100., 100.), (100.001, 50.)],
    ]);
    let mut d = detector();
    let t = d.run(&e, 0).unwrap();
    assert_eq!(t.iter().map(MassTrace::len).collect::<Vec<_>>(), [3, 2]);
    add_arrays(&mut e, &[1., 1., 1.]);
    d.options.ion_mobility_tolerance = -0.1;
    assert!(d.run(&e, 0).unwrap().iter().all(|t| t.len() == 1));
    d.options.ion_mobility_tolerance = 0.;
    assert_eq!(d.run(&e, 1).unwrap()[0].len(), 3);
}

#[test]
fn owned_output_replacement_does_not_inspect_prior_values() {
    let e = input(&[&[(100., 40.)], &[(100., 100.)], &[(100., 60.)]]);
    let mut d = detector();
    let mut old = vec![MassTrace::from_peaks(vec![Peak2D::new(f64::NAN, 0., 0.)]).unwrap()];
    let ptr = old.as_ptr();
    let returned = d.run_into(&e, &mut old, 0).unwrap();
    assert_eq!(returned.as_ptr(), ptr);
    assert!(returned[0][0].rt().is_nan());
    assert_eq!(old[0].len(), 3);
}

#[test]
fn area_success_atomic_failure_noop_and_source_sentinel_rt_grouping() {
    let e = input(&[&[(100., 40.)], &[(100., 100.)], &[(100., 60.)]]);
    let mut d = detector();
    let mut begin = e.area_iter(AreaOptions::default()).unwrap();
    let end = AreaIter::default();
    let mut output = vec![];
    assert!(d.run_area(&mut begin, &end, &mut output).unwrap().is_some());
    assert_eq!(begin, end);
    assert_eq!(output[0].len(), 3);
    d.options.mass_error_ppm = f64::NAN;
    d.limits.max_work = 0;
    output[0].peaks_mut()[0].intensity = f32::NAN;
    let pointer = output.as_ptr();
    assert!(d.run_area(&mut begin, &end, &mut output).unwrap().is_none());
    assert_eq!(pointer, output.as_ptr());
    d = detector();
    let mut begin = e.area_iter(AreaOptions::default()).unwrap();
    let mut middle = begin.clone();
    middle.next();
    let saved = begin.clone();
    assert!(d.run_area(&mut begin, &middle, &mut output).is_err());
    assert_eq!(begin, saved);
    assert_eq!(pointer, output.as_ptr());
    let mut e = e.clone();
    e.spectra.insert(0, scan(-1., &[(100., 1000.)]));
    let mut begin = e.area_iter(AreaOptions::default()).unwrap();
    d.run_area(&mut begin, &end, &mut output).unwrap();
    assert_eq!(output[0].len(), 3);
    assert_eq!(output[0][0].rt(), 0.); // sentinel group was discarded
    e.spectra[1].rt = -1.;
    let mut begin = e.area_iter(AreaOptions::default()).unwrap();
    let saved = begin.clone();
    assert!(d.run_area(&mut begin, &end, &mut output).is_err()); // only two reconstructed scans
    assert_eq!(begin, saved);
}

#[test]
fn checked_coordinates_order_counts_work_and_bytes_fail_atomically() {
    let good = input(&[&[(100., 40.)], &[(100., 100.)], &[(100., 60.)]]);
    let mut d = detector();
    let expected = d.run(&good, 0).unwrap();
    let mut output = expected.clone();
    for kind in 0..7 {
        let mut e = good.clone();
        d = detector();
        match kind {
            0 => e.spectra[0].peaks[0].mz = f64::NAN,
            1 => e.spectra[0].peaks.push(Peak1D::new(90., 40.)),
            2 => d.limits.max_spectra = 2,
            3 => d.limits.max_peaks = 2,
            4 => d.limits.max_work = 10,
            5 => d.limits.max_bytes = 1,
            _ => d.limits.max_traces = 0,
        }
        assert!(d.run_into(&e, &mut output, 0).is_err(), "{kind}");
        assert_eq!(output, expected);
    }
    // Defined finite negative configurations with re-estimation disabled remain usable.
    d = detector();
    d.options.mass_error_ppm = -1.;
    assert!(d.run(&good, 0).unwrap().iter().all(|t| t.len() == 1));
    d.options.mass_error_ppm = 0.;
    assert_eq!(d.run(&good, 0).unwrap()[0].len(), 3);
    // Zero m/z makes the source incremental ratio nonfinite: checked rejection.
    let zero = input(&[&[(0., 40.)], &[(0., 100.)], &[(0., 60.)]]);
    assert!(d.run(&zero, 0).is_err());
}

#[test]
fn independent_grid_enumeration_matches_complete_peak_membership_and_apex_order() {
    // Equal m/z groups are isolated by much more than the ppm tolerance. A
    // direct grouping oracle needs neither incremental means nor nearest search.
    let mut seed = 19u64;
    for _ in 0..100 {
        let mut e = MSExperiment::default();
        let mut groups: Vec<Vec<Peak2D>> = vec![vec![], vec![], vec![]];
        let mut apices: Vec<Option<(u32, usize, usize)>> = vec![None; 3];
        for i in 0..6 {
            let mut s = scan(i as f64 * 1.25, &[]);
            for g in 0..3 {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                let intensity = [0u32, 10, 20, 30, 40, 100][(seed >> 32) as usize % 6];
                if intensity == 0 {
                    continue;
                }
                let mz = (g + 1) as f64 * 100.;
                s.peaks.push(Peak1D::new(mz, intensity as f32));
                if intensity > 10 {
                    groups[g].push(Peak2D::new(s.rt, mz, intensity as f32));
                }
                if intensity > 30 {
                    apices[g] =
                        Some(apices[g].map_or((intensity, i, g), |old| old.max((intensity, i, g))));
                }
            }
            e.spectra.push(s);
        }
        let mut order: Vec<_> = apices
            .into_iter()
            .enumerate()
            .filter_map(|(g, a)| a.map(|a| (a, g)))
            .collect();
        order.sort_by(|a, b| b.cmp(a));
        let mut d = detector();
        d.options.min_sample_rate = 0.;
        d.options.trace_termination_outliers = 6;
        let observed = d.run(&e, 0).unwrap();
        assert_eq!(observed.len(), order.len());
        for (trace, (_, group)) in observed.iter().zip(order) {
            assert_eq!(trace.peaks(), groups[group]);
            close(trace.centroid_mz(), (group + 1) as f64 * 100., 1e-10);
            let expected_rt = if trace.len() == 1 {
                trace[0].rt()
            } else {
                let (mut numerator, mut denominator) = (0., 0.);
                for pair in groups[group].windows(2) {
                    let weight = f64::from(pair[1].intensity) * (pair[1].rt() - pair[0].rt());
                    numerator += weight * pair[1].rt();
                    denominator += weight;
                }
                numerator / denominator
            };
            close(trace.centroid_rt(), expected_rt, 1e-12);
        }
    }
}

#[test]
fn progress_failure_balances_nesting_and_preserves_last_successful_state() {
    use openms::concept::progress_logger::{
        ProgressBackend, ProgressLogger, ProgressNesting, ProgressTime,
    };
    use std::sync::{Arc, Mutex};
    struct Backend(Arc<Mutex<Vec<&'static str>>>);
    impl ProgressBackend for Backend {
        fn start_progress(&mut self, _: i64, _: i64, _: &str, _: usize) -> openms::Result<()> {
            self.0.lock().unwrap().push("start");
            Ok(())
        }
        fn set_progress(&mut self, _: i64, _: usize) -> openms::Result<()> {
            Ok(())
        }
        fn next_progress(&mut self) -> openms::Result<i64> {
            Ok(0)
        }
        fn end_progress(&mut self, _: usize, _: u64) -> openms::Result<()> {
            self.0.lock().unwrap().push("end");
            Err(openms::Error::InvalidValue("injected end failure".into()))
        }
    }
    let mut e = input(&[&[(100., 40.)], &[(100., 100.)], &[(100., 60.)]]);
    add_arrays(&mut e, &[1., 1., 1.]);
    let mut d = detector();
    let mut output = d.run(&e, 0).unwrap();
    let before = output.clone();
    let nesting = ProgressNesting::default();
    let calls = Arc::new(Mutex::new(vec![]));
    d.logger = ProgressLogger::with_clock_and_nesting(
        Arc::new(|| {
            Ok(ProgressTime {
                wall_second: 0,
                wall_seconds: 0.,
                cpu_seconds: None,
            })
        }),
        nesting.clone(),
    );
    d.logger.set_logger(Box::new(Backend(calls.clone())));
    for s in &mut e.spectra {
        s.float_data_arrays.clear();
    }
    assert!(d.run_into(&e, &mut output, 0).is_err());
    assert_eq!(output, before);
    assert!(d.has_centroid_im());
    assert_eq!(nesting.depth(), 0);
    assert_eq!(*calls.lock().unwrap(), ["start", "end"]);
}
