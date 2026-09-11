// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::Result;
use openms::analysis::elution_peak_detection::{
    ElutionPeakDetection, ElutionPeakDetectionOptions, ElutionPeakWidthFiltering,
};
use openms::analysis::mass_trace_detection::MassTraceDetection;
use openms::concept::progress_logger::{
    ProgressBackend, ProgressLogType, ProgressLogger, ProgressNesting, ProgressTime,
};
use openms::kernel::{MSExperiment, MSSpectrum, MassTrace, MassTraceQuantMethod, Peak1D, Peak2D};
use openms::processing::smoothing::SavitzkyGolayFilter;
use std::sync::{Arc, Mutex};

fn detector() -> ElutionPeakDetection {
    let mut d = ElutionPeakDetection::new();
    d.logger.set_log_type(ProgressLogType::None);
    d.options.width_filtering = ElutionPeakWidthFiltering::Off;
    d
}
fn trace(values: &[f32]) -> MassTrace {
    let mut t = MassTrace::from_peaks(
        values
            .iter()
            .enumerate()
            .map(|(i, &v)| Peak2D::new(i as f64, 400. + i as f64, v))
            .collect(),
    )
    .unwrap();
    t.set_label("native").unwrap();
    t
}
fn smoothed(values: &[f64]) -> MassTrace {
    let mut t = trace(&values.iter().map(|&v| v as f32).collect::<Vec<_>>());
    t.set_smoothed_intensities(values).unwrap();
    t
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17} != {expected:.17}"
    );
}
fn source_trace() -> MassTrace {
    let mut e = MSExperiment::default();
    for line in include_str!("data/elution_peak_detection_source.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let v: Vec<_> = line.split('\t').collect();
        let index: usize = v[0].parse().unwrap();
        if index == e.spectra.len() {
            e.spectra.push(MSSpectrum {
                rt: v[1].parse().unwrap(),
                ms_level: v[2].parse().unwrap(),
                ..Default::default()
            });
        }
        e.spectra[index]
            .peaks
            .push(Peak1D::new(v[3].parse().unwrap(), v[4].parse().unwrap()));
    }
    assert_eq!(e.spectra.len(), 333);
    let mut mtd = MassTraceDetection::new();
    mtd.logger.set_log_type(ProgressLogType::None);
    let mut traces = mtd.run(&e, 0).unwrap();
    assert_eq!(traces.len(), 1);
    assert_eq!(traces[0].label(), "T1");
    traces.remove(0)
}

#[test]
fn source_defaults_and_exact_projected_fixture_literals() {
    let native = ElutionPeakDetection::new();
    assert_eq!(native.options, ElutionPeakDetectionOptions::default());
    assert_eq!(native.logger.log_type(), ProgressLogType::Cmd);
    let mut d = detector();
    let mut inputs = vec![source_trace()];
    let split = d.detect_peaks_many(&mut inputs).unwrap();
    assert_eq!(split.len(), 3);
    assert_eq!(split[0].label(), "T1.1");
    assert_eq!(split[1].label(), "T1.2");
    // Third label follows the source expression; its upstream assertion is commented out.
    assert_eq!(split[2].label(), "T1.3");
    close(
        d.compute_mass_trace_noise(&inputs[0]).unwrap(),
        573.8585,
        573.8585 * 0.01,
    );
    for (i, t) in split.iter().enumerate() {
        let snr = [0.1907, 9.8855, 7.6432][i];
        let apex = [2.0427, 37.7893, 52.9933][i];
        close(d.compute_mass_trace_snr(t).unwrap(), snr, snr * 0.01);
        close(d.compute_apex_snr(t).unwrap(), apex, apex * 0.01);
    }
    let mut t = source_trace();
    d.smooth_data(&mut t, 20).unwrap();
    let e = d.find_local_extrema(&t, 10).unwrap();
    assert_eq!((e.maxima.len(), e.minima.len()), (4, 2));
    let e = d.find_local_extrema(&t, 70).unwrap();
    assert_eq!((e.maxima.len(), e.minima.len()), (2, 1));
}

#[test]
fn closed_form_noise_area_apex_and_empty_branches() {
    let d = detector();
    let mut t = trace(&[1., 4., 1.]);
    t.set_smoothed_intensities(&[1., 3., 1.]).unwrap();
    close(
        d.compute_mass_trace_noise(&t).unwrap(),
        (1.0_f64 / 3.).sqrt(),
        1e-15,
    );
    close(
        d.compute_mass_trace_snr(&t).unwrap(),
        2.5 * 3.0_f64.sqrt(),
        1e-14,
    );
    close(d.compute_apex_snr(&t).unwrap(), 3. * 3.0_f64.sqrt(), 1e-14);
    let empty = MassTrace::new();
    assert_eq!(d.compute_mass_trace_noise(&empty).unwrap(), 0.);
    assert_eq!(d.compute_mass_trace_snr(&empty).unwrap(), 0.);
    assert_eq!(d.compute_apex_snr(&empty).unwrap(), 0.);
    let unsmoothed = trace(&[1., 2., 1.]);
    assert_eq!(d.compute_mass_trace_noise(&unsmoothed).unwrap(), 0.);
    assert_eq!(d.compute_apex_snr(&unsmoothed).unwrap(), 0.);
    assert!(d.compute_mass_trace_snr(&unsmoothed).is_err());
    assert!(d.compute_mass_trace_snr(&smoothed(&[1.])).is_err());
}

#[test]
fn smoothing_short_even_odd_long_negative_and_f32_storage() {
    let d = detector();
    for raw in [vec![], vec![-2.], vec![-2., 3.]] {
        let mut t = trace(&raw);
        d.smooth_data(&mut t, i32::MAX).unwrap();
        assert_eq!(
            t.smoothed_intensities(),
            raw.iter().map(|&x| f64::from(x)).collect::<Vec<_>>()
        );
    }
    let raw = [-1., 2., 12., 3., -2., 1., 0.];
    let mut t = trace(&raw);
    d.smooth_data(&mut t, 4).unwrap();
    let a = t.smoothed_intensities().to_vec();
    d.smooth_data(&mut t, 5).unwrap();
    assert_eq!(a, t.smoothed_intensities());
    let positions: Vec<_> = (0..raw.len()).map(|x| x as f64).collect();
    let raw64: Vec<_> = raw.iter().map(|&x| f64::from(x)).collect();
    let fitted = SavitzkyGolayFilter::new(5, 2)
        .unwrap()
        .filter(&positions, &raw64)
        .unwrap();
    assert!(fitted.iter().zip(&a).any(|(&x, &y)| x != y));
    assert_eq!(
        a,
        fitted
            .iter()
            .map(|&x| f64::from(x as f32))
            .collect::<Vec<_>>()
    );
    assert!(a.iter().all(|&x| x >= 0.));
    d.smooth_data(&mut t, 9).unwrap();
    assert_eq!(t.smoothed_intensities(), raw64);
    let saved = t.clone();
    assert!(d.smooth_data(&mut t, 1024).is_err());
    assert_eq!(t, saved);
    d.smooth_data(&mut t, -1).unwrap();
    let negative = t.smoothed_intensities().to_vec();
    d.smooth_data(&mut t, 0).unwrap();
    assert_eq!(negative, t.smoothed_intensities());
}

#[test]
fn extrema_plateaus_exclusive_right_edge_and_separation_option() {
    let mut d = detector();
    assert!(
        d.find_local_extrema(&smoothed(&[]), 0)
            .unwrap()
            .maxima
            .is_empty()
    );
    assert_eq!(
        d.find_local_extrema(&smoothed(&[-1., -1.]), 100)
            .unwrap()
            .maxima,
        [0]
    );
    assert!(d.find_local_extrema(&trace(&[1., 2., 3.]), 2).is_err());
    // Neighborhood one excludes the successor, hence the ascending ramp seeds all indices.
    assert_eq!(
        d.find_local_extrema(&smoothed(&[1., 2., 3.]), 1)
            .unwrap()
            .maxima,
        [0, 1, 2]
    );
    assert_eq!(
        d.find_local_extrema(&smoothed(&[2., 2., 2., 2.]), 2)
            .unwrap()
            .maxima,
        [0, 2]
    );
    assert!(
        d.find_local_extrema(&smoothed(&[0., -1., 0.]), 0)
            .unwrap()
            .maxima
            .is_empty()
    );
    let t = smoothed(&[4., 1., 1., 4.]);
    let e = d.find_local_extrema(&t, 2).unwrap();
    assert_eq!(e.maxima, [0, 3]);
    assert_eq!(e.minima, [2]); // Equal valleys choose the right boundary.
    d.options.min_fwhm = 2.;
    assert_eq!(d.find_local_extrema(&t, 2).unwrap().minima, [2]);
    d.options.min_fwhm = 2.0001;
    assert!(d.find_local_extrema(&t, 2).unwrap().minima.is_empty());
    d.options.width_filtering = ElutionPeakWidthFiltering::Auto;
    assert!(d.find_local_extrema(&t, 2).unwrap().minima.is_empty());
    d.options.min_fwhm = -10.;
    assert_eq!(d.find_local_extrema(&t, 2).unwrap().minima, [2]);
    assert!(d.find_local_extrema(&t, usize::MAX).is_err());
}

#[test]
fn independent_exhaustive_small_grid_maxima() {
    let d = detector();
    // Independent enumeration of all 4^5 signals, sorted candidate membership
    // and overlapping-window suppression; minima are covered by literal cases.
    for encoding in 0..1024usize {
        let values: Vec<f64> = (0..5)
            .map(|i| ((encoding >> (2 * i)) & 3) as f64 - 1.)
            .collect();
        for neighbors in 0..=5 {
            let mut candidates: Vec<_> = (0..5)
                .filter(|&i| {
                    values[i] > 0.
                        && (0..5).all(|j| {
                            j < i.saturating_sub(neighbors)
                                || j >= (i + neighbors).min(5)
                                || values[j] <= values[i]
                        })
                })
                .collect();
            candidates.sort_by(|&a, &b| values[a].partial_cmp(&values[b]).unwrap().then(a.cmp(&b)));
            let mut taken = Vec::new();
            for candidate in candidates {
                if taken.iter().all(|&i: &usize| {
                    candidate < i.saturating_sub(neighbors) || candidate >= (i + neighbors).min(5)
                }) {
                    taken.push(candidate);
                }
            }
            taken.sort_unstable();
            assert_eq!(
                d.find_local_extrema(&smoothed(&values), neighbors)
                    .unwrap()
                    .maxima,
                taken,
                "{values:?}, {neighbors}"
            );
        }
    }
}

#[test]
fn successful_rejection_publishes_smoothing_and_single_fwhm_state() {
    let mut d = detector();
    d.options.chrom_fwhm = 10.; // Too-wide SG frame copies raw, single maximum.
    let mut t = trace(&[1., 5., 1.]);
    t.set_centroid_sd(7.).unwrap();
    t.update_mean_mz().unwrap();
    t.update_median_rt().unwrap();
    let mz = t.centroid_mz();
    d.options.width_filtering = ElutionPeakWidthFiltering::Fixed;
    d.options.min_fwhm = 5.;
    assert!(d.detect_peaks(&mut t).unwrap().is_empty());
    assert_eq!(t.smoothed_intensities(), [1., 5., 1.]);
    close(t.fwhm(), 1.25, 1e-15);
    assert_eq!(t.fwhm_borders(), (0, 2));
    assert_eq!(t.centroid_mz(), mz);
    assert_eq!(t.centroid_sd(), 7.);
    d.options.min_fwhm = -1.;
    let out = d.detect_peaks(&mut t).unwrap();
    assert_eq!(out, [t.clone()]);
    assert_eq!(t.label(), "native");
    let mut no_max = trace(&[-1., -2., -1.]);
    assert!(d.detect_peaks(&mut no_max).unwrap().is_empty());
    assert_eq!(no_max.smoothed_intensities(), [-1., -2., -1.]);
}

#[test]
fn source_split_snr_uses_original_and_preserves_state_and_nonoverlap() {
    let mut d = detector();
    let mut original = source_trace();
    original.set_centroid_im(2.5).unwrap();
    original.fwhm_mz_avg = 4.;
    original.fwhm_im_avg = 5.;
    original.set_quant_method(MassTraceQuantMethod::MaxHeight);
    let before = original.clone();
    let split = d.detect_peaks(&mut original).unwrap();
    assert_eq!(original.centroid_rt(), before.centroid_rt());
    assert_eq!(original.fwhm(), before.fwhm());
    assert_eq!(
        split
            .iter()
            .flat_map(|t| t.peaks().iter().copied())
            .collect::<Vec<_>>(),
        before.peaks()
    );
    for t in &split {
        assert_eq!(t.centroid_im(), 2.5);
        assert!(t.contains_im_data());
        assert_eq!(t.fwhm_mz_avg, 4.);
        assert_eq!(t.fwhm_im_avg, 5.);
        assert_eq!(t.quant_method(), MassTraceQuantMethod::MaxHeight);
    }
    let threshold = d.compute_apex_snr(&original).unwrap() / 2.;
    assert!(
        split
            .iter()
            .any(|t| d.compute_apex_snr(t).unwrap() < threshold)
    );
    d.options.masstrace_snr_filtering = true;
    d.options.chrom_peak_snr = threshold;
    let output = d.detect_peaks(&mut before.clone()).unwrap();
    assert_eq!(output, split);
    d.options.chrom_peak_snr = threshold * 2. + 1.;
    assert!(d.detect_peaks(&mut before.clone()).unwrap().is_empty());
}

#[test]
fn auto_is_explicit_and_width_quantiles_are_inclusive_stable() {
    let d = detector();
    for n in [0, 1, 19, 20, 21] {
        let mut traces: Vec<_> = (0..n)
            .map(|i| {
                let mut t = smoothed(&[0., 2., 0.]);
                let width = (n - i) as f64;
                for p in t.peaks_mut() {
                    p.position[0] *= width;
                }
                t.set_label(&i.to_string()).unwrap();
                t
            })
            .collect();
        let output = d.filter_by_peak_width(&mut traces).unwrap();
        let lower = (n as f64 * 0.05).floor() as usize;
        let upper = (n as f64 * 0.95).floor() as usize;
        let expected: Vec<_> = (0..n)
            .rev()
            .enumerate()
            .filter(|&(rank, _)| rank >= lower && rank <= upper)
            .map(|(_, i)| i.to_string())
            .collect();
        assert_eq!(
            output.iter().map(MassTrace::label).collect::<Vec<_>>(),
            expected
        );
        assert!(traces.iter().all(|t| t.fwhm() > 0.));
        if n == 20 {
            assert_eq!(output.len(), 19);
        }
    }
    let mut ties = vec![smoothed(&[0., 2., 0.]); 3];
    for (i, t) in ties.iter_mut().enumerate() {
        t.set_label(&i.to_string()).unwrap();
    }
    assert_eq!(
        d.filter_by_peak_width(&mut ties)
            .unwrap()
            .iter()
            .map(MassTrace::label)
            .collect::<Vec<_>>(),
        ["0", "1", "2"]
    );
    let mut d = detector();
    let mut off = source_trace();
    let a = d.detect_peaks(&mut off).unwrap();
    d.options.width_filtering = ElutionPeakWidthFiltering::Auto;
    let mut auto = source_trace();
    assert_eq!(d.detect_peaks(&mut auto).unwrap(), a);
    assert_eq!(off, auto);
}

#[test]
fn late_failure_preserves_all_inputs_output_and_unused_old_output_is_moved() {
    let mut d = detector();
    d.options.chrom_fwhm = 10.;
    let mut bad = trace(&[1., 5., 1.]);
    for p in bad.peaks_mut() {
        p.position[0] = 0.;
    }
    let mut batch = vec![trace(&[1., 5., 1.]), bad];
    let saved = batch.clone();
    let mut output = vec![trace(&[9.])];
    let old = output.clone();
    assert!(d.detect_peaks_many_into(&mut batch, &mut output).is_err());
    assert_eq!(batch, saved);
    assert_eq!(output, old);
    let mut widths = vec![smoothed(&[0., 2., 0.]), trace(&[1., 2., 1.])];
    let saved = widths.clone();
    assert!(
        d.filter_by_peak_width_into(&mut widths, &mut output)
            .is_err()
    );
    assert_eq!(widths, saved);
    assert_eq!(output, old);
    output = vec![trace(&vec![f32::NAN; 10000])];
    let old_pointer = output.as_ptr();
    d.limits.max_peaks = 3;
    let previous = d
        .detect_peaks_into(&mut trace(&[1., 5., 1.]), &mut output)
        .unwrap();
    assert_eq!(previous.as_ptr(), old_pointer);
    assert_eq!(previous[0].len(), 10000);
    assert_eq!(output.len(), 1);
}

#[test]
fn invalid_conversion_ordering_finite_results_and_cumulative_limits() {
    let mut d = detector();
    for values in [vec![], vec![1.]] {
        let mut t = trace(&values);
        assert!(d.detect_peaks(&mut t).is_err());
    }
    let mut t = trace(&[1., 3., 1.]);
    let saved = t.clone();
    for value in [-1., f64::MAX, f64::NAN, f64::INFINITY] {
        d.options.chrom_fwhm = value;
        assert!(d.detect_peaks(&mut t).is_err());
        assert_eq!(t, saved);
    }
    d.options.chrom_fwhm = 10.;
    t.peaks_mut().swap(0, 1);
    assert!(d.smooth_data(&mut t, 5).is_err());
    let mut t = saved.clone();
    let original_limits = d.limits;
    for bytes in [0, 1, 100] {
        d.limits.max_bytes = bytes;
        assert!(d.detect_peaks(&mut t).is_err());
        assert_eq!(t, saved);
    }
    d.limits = original_limits;
    d.limits.max_work = 100;
    assert!(
        d.detect_peaks_many(&mut [saved.clone(), saved.clone()])
            .is_err()
    );
    d.limits = original_limits;
    d.limits.max_peaks = 5;
    assert!(
        d.detect_peaks_many(&mut [saved.clone(), saved.clone()])
            .is_err()
    );
    d.limits = original_limits;
    d.limits.max_traces = 0;
    assert!(d.detect_peaks(&mut t).is_err());
    d.limits = original_limits;
    d.options.chrom_peak_snr = f64::NAN; // Unused settings do not poison smoothing.
    d.options.min_fwhm = f64::NAN;
    d.smooth_data(&mut t, 5).unwrap();
}

#[test]
fn batch_progress_events_and_single_no_events() {
    struct Recorder(Arc<Mutex<Vec<String>>>);
    impl ProgressBackend for Recorder {
        fn start_progress(&mut self, a: i64, b: i64, label: &str, _: usize) -> Result<()> {
            self.0
                .lock()
                .unwrap()
                .push(format!("start {a} {b} {label}"));
            Ok(())
        }
        fn set_progress(&mut self, v: i64, _: usize) -> Result<()> {
            self.0.lock().unwrap().push(format!("set {v}"));
            Ok(())
        }
        fn next_progress(&mut self) -> Result<i64> {
            unreachable!()
        }
        fn end_progress(&mut self, _: usize, _: u64) -> Result<()> {
            self.0.lock().unwrap().push("end".into());
            Ok(())
        }
    }
    let ticks = Arc::new(Mutex::new(0i64));
    let clock = {
        let ticks = ticks.clone();
        Arc::new(move || {
            let mut t = ticks.lock().unwrap();
            *t += 1;
            Ok(ProgressTime {
                wall_second: *t,
                wall_seconds: *t as f64,
                cpu_seconds: None,
            })
        })
    };
    let events = Arc::new(Mutex::new(Vec::new()));
    let mut d = detector();
    d.logger = ProgressLogger::with_clock_and_nesting(clock, ProgressNesting::default());
    d.logger.set_logger(Box::new(Recorder(events.clone())));
    d.options.chrom_fwhm = 10.;
    d.detect_peaks(&mut trace(&[1., 3., 1.])).unwrap();
    assert!(events.lock().unwrap().is_empty());
    d.detect_peaks_many(&mut [trace(&[1., 3., 1.]), trace(&[1., 3., 1.])])
        .unwrap();
    assert_eq!(
        *events.lock().unwrap(),
        ["start 0 2 elution peak detection", "set 0", "set 1", "end"]
    );
}

#[test]
fn zero_window_rejected_segment_labels_and_multiple_maxima_without_valleys() {
    let mut d = detector();
    d.options.chrom_fwhm = 0.;
    d.options.width_filtering = ElutionPeakWidthFiltering::Fixed;
    let mut t = trace(&[9., 1., 1., 1., 8., 1.]);
    let output = d.detect_peaks(&mut t).unwrap();
    assert_eq!(output.len(), 1);
    assert_eq!(output[0].label(), "native.2");
    assert_eq!(output[0].peaks(), &t.peaks()[3..]);
    d.options.width_filtering = ElutionPeakWidthFiltering::Off;
    let mut plateau = trace(&[1., 1., 1.]);
    let output = d.detect_peaks(&mut plateau).unwrap();
    assert_eq!(output.len(), 1);
    assert_eq!(output[0].label(), "native.1");
    assert_eq!(output[0].peaks(), plateau.peaks());
}

#[test]
fn rejected_width_still_evaluates_snr_and_batch_allowances_are_cumulative() {
    let mut d = detector();
    d.options.chrom_fwhm = 10.;
    d.options.width_filtering = ElutionPeakWidthFiltering::Fixed;
    d.options.min_fwhm = 100.;
    d.options.masstrace_snr_filtering = true;
    d.options.chrom_peak_snr = f64::NAN;
    let mut t = trace(&[1., 5., 1.]);
    let saved = t.clone();
    assert!(d.detect_peaks(&mut t).is_err());
    assert_eq!(t, saved);
    d = detector();
    d.options.chrom_fwhm = 10.;
    // Find each operation's actual acceptance boundary, then prove the same
    // allowance cannot silently reset while processing a second trace.
    for bytes in [false, true] {
        let mut low = 0;
        let mut high = 100_000;
        while low + 1 < high {
            let mid = (low + high) / 2;
            if bytes {
                d.limits.max_bytes = mid;
            } else {
                d.limits.max_work = mid;
            }
            if d.detect_peaks_many(&mut [saved.clone()]).is_ok() {
                high = mid;
            } else {
                low = mid;
            }
        }
        if bytes {
            d.limits.max_bytes = high;
        } else {
            d.limits.max_work = high;
        }
        assert!(d.detect_peaks_many(&mut [saved.clone()]).is_ok());
        let mut batch = [saved.clone(), saved.clone()];
        let previous = batch.clone();
        assert!(d.detect_peaks_many(&mut batch).is_err());
        assert_eq!(batch, previous);
        d.limits = Default::default();
    }
}

#[test]
fn progress_end_failure_precedes_scientific_publication_and_balances_nesting() {
    struct EndFailure;
    impl ProgressBackend for EndFailure {
        fn start_progress(&mut self, _: i64, _: i64, _: &str, _: usize) -> Result<()> {
            Ok(())
        }
        fn set_progress(&mut self, _: i64, _: usize) -> Result<()> {
            Ok(())
        }
        fn next_progress(&mut self) -> Result<i64> {
            unreachable!()
        }
        fn end_progress(&mut self, _: usize, _: u64) -> Result<()> {
            Err(openms::Error::Unsupported("end failed".into()))
        }
    }
    let nesting = ProgressNesting::default();
    let mut d = detector();
    d.logger = ProgressLogger::with_clock_and_nesting(
        Arc::new(|| {
            Ok(ProgressTime {
                wall_second: 1,
                wall_seconds: 1.,
                cpu_seconds: None,
            })
        }),
        nesting.clone(),
    );
    d.logger.set_logger(Box::new(EndFailure));
    d.options.chrom_fwhm = 10.;
    let mut input = [trace(&[1., 5., 1.])];
    let saved = input.clone();
    let mut output = vec![trace(&[7.])];
    let old = output.clone();
    assert!(d.detect_peaks_many_into(&mut input, &mut output).is_err());
    assert_eq!(input, saved);
    assert_eq!(output, old);
    assert_eq!(nesting.depth(), 0);
}
