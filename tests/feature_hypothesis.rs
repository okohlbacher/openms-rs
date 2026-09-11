// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::{
    analysis::feature_hypothesis::{
        FeatureHypothesis, FeatureHypothesisLimits, MetaboIsotopeMassWindow,
        hypothesis_score_greater, mass_trace_mz_less,
    },
    kernel::{MassTrace, MassTraceQuantMethod as Quant, Peak2D, Point2D, Precursor},
    metadata::{ChromatogramType, Product},
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17} != {expected:.17}"
    );
}
fn trace(label: &str, mz: f64, values: &[f32]) -> MassTrace {
    let mut trace = MassTrace::from_peaks(
        values
            .iter()
            .enumerate()
            .map(|(i, &v)| Peak2D::new(i as f64, mz, v))
            .collect(),
    )
    .unwrap();
    trace.set_label(label).unwrap();
    if !values.is_empty() {
        trace.update_median_rt().unwrap();
        trace.update_median_mz().unwrap();
    }
    trace
}
fn hypothesis<'a>(traces: impl IntoIterator<Item = &'a MassTrace>) -> FeatureHypothesis<'a> {
    let mut hypothesis = FeatureHypothesis::new();
    for trace in traces {
        hypothesis.add_mass_trace(trace).unwrap();
    }
    hypothesis
}

#[test]
fn empty_source_branches_and_native_checked_chromatogram_boundary() {
    let empty = FeatureHypothesis::new();
    assert!(empty.is_empty());
    assert_eq!(empty.len(), 0);
    assert!(empty.traces().is_empty());
    assert_eq!(empty.score(), 0.);
    assert_eq!(empty.charge(), 0);
    assert_eq!(empty.fwhm(), 0.);
    assert_eq!(empty.number_of_feature_points(), 0);
    assert_eq!(empty.label().unwrap(), "");
    assert!(empty.labels().unwrap().is_empty());
    assert!(empty.all_centroid_mz().unwrap().is_empty());
    assert!(empty.all_centroid_rt().unwrap().is_empty());
    assert!(empty.all_centroid_im().unwrap().is_empty());
    assert!(empty.isotope_distances().unwrap().is_empty());
    assert!(empty.convex_hulls().unwrap().is_empty());
    for smoothed in [false, true] {
        assert!(empty.all_intensities(smoothed).unwrap().is_empty());
        assert_eq!(empty.summed_feature_intensity(smoothed).unwrap(), 0.);
        assert_eq!(empty.max_intensity(smoothed).unwrap(), 0.);
        assert!(empty.monoisotopic_feature_intensity(smoothed).is_err());
    }
    assert!(empty.centroid_mz().is_err());
    assert!(empty.centroid_rt().is_err());
    assert!(empty.chromatograms(0).is_err());
}

#[test]
fn borrowed_membership_copy_and_source_label_separators() {
    let a = trace("", 102., &[1., 5., 1.]);
    let b = trace("T_2", 100., &[3.]);
    let mut h = hypothesis([&a, &b, &a]);
    h.set_score(-0.);
    h.set_charge(-17);
    assert_eq!(h.labels().unwrap(), ["", "T_2", ""]);
    assert_eq!(h.label().unwrap(), "_T_2_");
    assert_eq!(h.number_of_feature_points(), 7);
    for copy in [h.clone(), h.checked_clone().unwrap()] {
        assert_ne!(copy.traces().as_ptr(), h.traces().as_ptr());
        assert!(std::ptr::eq(copy.traces()[0], &a));
        assert!(std::ptr::eq(copy.traces()[2], &a));
        assert_eq!(copy.score().to_bits(), (-0f64).to_bits());
        assert_eq!(copy.charge(), -17);
        assert_eq!(copy.number_of_feature_points(), 7);
    }
    let empty_label = hypothesis([&a, &a, &a]);
    assert_eq!(empty_label.label().unwrap(), "__");
}

#[test]
fn first_cached_summary_is_independent_of_raw_positions_and_member_order() {
    let mut a = trace("a", 102., &[1., 5., 1.]);
    a.estimate_fwhm(false).unwrap();
    a.set_centroid_im(-0.).unwrap();
    a.peaks_mut()[1].set_mz(900.);
    a.peaks_mut()[1].set_rt(77.);
    let b = trace("b", 100., &[2.]);
    let h = hypothesis([&a, &b, &a]);
    assert_eq!(h.centroid_mz().unwrap(), 102.);
    assert_eq!(h.centroid_rt().unwrap(), 1.);
    close(h.fwhm(), 1.25, 0.);
    assert_eq!(h.all_centroid_mz().unwrap(), [102., 100., 102.]);
    assert_eq!(h.all_centroid_rt().unwrap(), [1., 0., 1.]);
    let im = h.all_centroid_im().unwrap();
    assert_eq!(im[0].to_bits(), (-0f64).to_bits());
    assert_eq!(im[1].to_bits(), 0f64.to_bits());
    assert_eq!(h.isotope_distances().unwrap(), [-2., 2.]);
    assert!(!mass_trace_mz_less(&a, &b));
    assert!(mass_trace_mz_less(&b, &a));
    assert!(!mass_trace_mz_less(&a, &a));
}

#[test]
fn all_quantification_modes_keep_source_raw_smoothed_and_apex_distinctions() {
    let mut area = trace("area", 100., &[1., 5., 1.]);
    area.estimate_fwhm(false).unwrap();
    area.set_smoothed_intensities(&[2., 10., 2.]).unwrap();
    let mut median = trace("median", 101., &[1., 6., 9., 2.]);
    median.set_quant_method(Quant::Median);
    median.set_smoothed_intensities(&[100.; 4]).unwrap();
    let mut height = trace("height", 102., &[4., 8., 3.]);
    height.set_quant_method(Quant::MaxHeight);
    height.set_smoothed_intensities(&[1., 2., 1.]).unwrap();
    let h = hypothesis([&area, &median, &height]);
    assert_eq!(h.all_intensities(false).unwrap(), [6., 4., 8.]);
    assert_eq!(h.all_intensities(true).unwrap(), [12., 4., 2.]);
    assert_eq!(h.monoisotopic_feature_intensity(false).unwrap(), 6.);
    assert_eq!(h.monoisotopic_feature_intensity(true).unwrap(), 12.);
    assert_eq!(h.summed_feature_intensity(false).unwrap(), 18.);
    assert_eq!(h.summed_feature_intensity(true).unwrap(), 18.);
    assert_eq!(h.max_intensity(false).unwrap(), 9.);
    assert_eq!(h.max_intensity(true).unwrap(), 100.);
    let negative = trace("negative", 1., &[-8., -1.]);
    assert_eq!(hypothesis([&negative]).max_intensity(false).unwrap(), 0.);
    assert_eq!(hypothesis([&negative]).max_intensity(true).unwrap(), 0.);
}

#[test]
fn source_mass_trace_literal_projection_through_hypothesis() {
    // Literal seven-peak class-test input, not generated by this implementation.
    let rows = include_str!("data/feature_hypothesis_source.tsv");
    let mut peaks = Vec::new();
    for line in rows.lines().filter(|line| !line.starts_with('#')) {
        let row: Vec<f64> = line
            .split('\t')
            .map(|value| value.parse().unwrap())
            .collect();
        peaks.push(Peak2D::new(row[0], row[1], row[2] as f32));
    }
    let mut t = MassTrace::from_peaks(peaks).unwrap();
    t.set_label("source_trace").unwrap();
    t.update_weighted_mean_rt().unwrap();
    t.update_weighted_mean_mz().unwrap();
    t.estimate_fwhm(false).unwrap();
    let h = hypothesis([&t, &t]);
    close(h.centroid_rt().unwrap(), 155.214671250425, 1e-12);
    close(h.centroid_mz().unwrap(), 230.10188267861923, 1e-12);
    close(h.fwhm(), 2.15250939241199, 1e-12);
    close(
        h.monoisotopic_feature_intensity(false).unwrap(),
        69863097.2125001,
        1e-6,
    );
    close(
        h.summed_feature_intensity(false).unwrap(),
        139726194.4250002,
        2e-6,
    );
    assert_eq!(h.max_intensity(false).unwrap(), 33329536.);
    assert_eq!(h.isotope_distances().unwrap(), [0.]);
    assert_eq!(h.number_of_feature_points(), 14);
    assert_eq!(h.label().unwrap(), "source_trace_source_trace");
}

#[test]
fn hulls_use_raw_scan_envelopes_and_preserve_one_output_per_reference() {
    let mut t = trace("h", 100., &[1., 2., 3., 4.]);
    t.peaks_mut().copy_from_slice(&[
        Peak2D::new(2., 105., f32::NAN),
        Peak2D::new(1., 103., 0.),
        Peak2D::new(1., 101., 0.),
        Peak2D::new(2., 102., 0.),
    ]);
    let empty = MassTrace::new();
    let hulls = hypothesis([&t, &empty, &t]).convex_hulls().unwrap();
    assert_eq!(hulls.len(), 3);
    assert!(hulls[1].is_empty());
    assert_eq!(hulls[0], hulls[2]);
    assert_eq!(
        hulls[0].hull_points(),
        [
            Point2D::new(1., 101.),
            Point2D::new(2., 102.),
            Point2D::new(2., 105.),
            Point2D::new(1., 103.),
        ]
    );
}

#[test]
fn chromatograms_use_raw_peaks_first_precursor_source_names_and_empty_defaults() {
    let mut a = trace("not_the_name", 500., &[1., 5., 1., 8.]);
    a.set_smoothed_intensities(&[100.; 4]).unwrap();
    a.peaks_mut().copy_from_slice(&[
        Peak2D::new(2., f64::NAN, 3.),
        Peak2D::new(-0., 900., 4.),
        Peak2D::new(0., 901., -0.),
        Peak2D::new(2., 902., -5.),
    ]);
    let b = trace("other", 501., &[7.]);
    let empty = MassTrace::new();
    let mut h = hypothesis([&a, &b, &empty]);
    h.set_charge(-3);
    let result = h.chromatograms(u64::MAX).unwrap();
    assert_eq!(result.len(), 3);
    for (index, chrom) in result.iter().enumerate() {
        assert_eq!(chrom.name, format!("18446744073709551615_{index}"));
        assert_eq!(chrom.native_id, chrom.name);
        assert_eq!(chrom.chromatogram_type, ChromatogramType::BasePeak);
        let mut expected_precursor = Precursor::new(500., -3);
        expected_precursor
            .cv_terms
            .metadata
            .insert("peptide_sequence".into(), "18446744073709551615".into());
        assert_eq!(chrom.precursor, expected_precursor);
        assert_eq!(chrom.product, Product::default());
        assert_eq!(chrom.instrument_settings, Default::default());
        assert_eq!(chrom.acquisition_info, Default::default());
        assert_eq!(chrom.source_file, Default::default());
        assert!(chrom.data_processing.is_empty());
        assert!(chrom.metadata.is_empty());
        assert!(chrom.float_data_arrays.is_empty());
        assert!(chrom.integer_data_arrays.is_empty());
        assert!(chrom.string_data_arrays.is_empty());
    }
    assert_eq!(
        result[0].peaks.iter().map(|p| p.rt).collect::<Vec<_>>(),
        [-0., 0., 2., 2.]
    );
    assert_eq!(
        result[0]
            .peaks
            .iter()
            .map(|p| p.intensity)
            .collect::<Vec<_>>(),
        [4., -0., 3., -5.]
    );
    assert_eq!(result[0].peaks[0].rt.to_bits(), (-0f64).to_bits());
    assert_eq!(result[0].peaks[1].intensity.to_bits(), (-0f32).to_bits());
    assert_eq!(result[1].precursor.mz, 500.);
    assert!(result[2].peaks.is_empty());
    assert_eq!(a[0].rt(), 2.); // source traces remain untouched
}

#[test]
fn score_predicates_retain_nan_infinity_and_signed_charge_storage() {
    let mut a = FeatureHypothesis::new();
    let mut b = FeatureHypothesis::new();
    a.set_score(f64::INFINITY);
    b.set_score(100.);
    assert!(hypothesis_score_greater(&a, &b));
    assert!(!hypothesis_score_greater(&b, &a));
    a.set_score(f64::NAN);
    assert!(a.score().is_nan());
    assert!(!hypothesis_score_greater(&a, &b));
    assert!(!hypothesis_score_greater(&b, &a));
    a.set_charge(i64::MIN);
    assert_eq!(a.charge(), i64::MIN);
    let range = MetaboIsotopeMassWindow {
        left_boundary: f64::INFINITY,
        right_boundary: f64::NEG_INFINITY,
    };
    assert!(range.left_boundary > range.right_boundary);
    let trace = MassTrace::new();
    let mut h = hypothesis([&trace]);
    for charge in [i64::MIN, i64::MAX] {
        h.set_charge(charge);
        assert_eq!(h.summed_feature_intensity(false).unwrap(), 0.);
        assert!(h.chromatograms(0).is_err());
    }
    for charge in [i32::MIN, i32::MAX] {
        h.set_charge(i64::from(charge));
        assert_eq!(h.chromatograms(0).unwrap()[0].precursor.charge, charge);
    }
}

#[test]
fn membership_limits_count_duplicates_and_fail_without_partial_addition() {
    let trace = trace("x", 1., &[1., 2.]);
    let mut h = FeatureHypothesis::with_limits(FeatureHypothesisLimits {
        max_traces: 3,
        max_peaks: 4,
        ..Default::default()
    });
    h.add_mass_trace(&trace).unwrap();
    h.add_mass_trace(&trace).unwrap();
    assert!(h.add_mass_trace(&trace).is_err());
    assert_eq!(h.len(), 2);
    assert_eq!(h.number_of_feature_points(), 4);
    h.limits.max_peaks = 6;
    h.limits.max_bytes = 0;
    assert!(h.add_mass_trace(&trace).is_err());
    assert_eq!(h.number_of_feature_points(), 4);
    h.limits.max_bytes = 1024;
    h.add_mass_trace(&trace).unwrap();
    assert!(h.add_mass_trace(&MassTrace::new()).is_err());
}

#[test]
fn label_clone_and_summaries_share_one_operation_budget() {
    let mut trace = trace("abcdef", 1., &[1., 2., 3., 4.]);
    trace.set_quant_method(Quant::Median);
    let mut h = hypothesis([&trace, &trace]);
    h.limits.max_bytes = 12;
    assert!(h.label().is_err()); // 6+1+6 bytes
    h.limits.max_bytes = 13;
    assert_eq!(h.label().unwrap(), "abcdef_abcdef");
    assert!(h.checked_clone().is_err()); // two borrowed pointer slots
    h.limits.max_bytes = usize::MAX;
    // Each four-value median costs 4 + 4*3*32 + 4 = 392 visits.
    h.limits.max_work = 2 + 392;
    assert_eq!(h.monoisotopic_feature_intensity(false).unwrap(), 2.5);
    assert!(h.summed_feature_intensity(false).is_err());
    h.limits.max_work = 2 + 2 * 392;
    assert_eq!(h.summed_feature_intensity(false).unwrap(), 5.);
}

#[test]
fn aggregate_geometry_and_chromatogram_bytes_do_not_reset_per_trace() {
    let t = trace("a", 1., &[1., 3., 1.]);
    let mut one = hypothesis([&t]);
    let mut two = hypothesis([&t, &t]);
    one.limits.max_bytes = 700;
    two.limits.max_bytes = 700;
    assert!(one.convex_hulls().is_ok());
    assert!(two.convex_hulls().is_err());
    // Discover a byte boundary from the declared public budgets, then verify
    // that identical repeated members consume the budget again.
    let mut minimum = None;
    for bytes in (1024..32_768).step_by(256) {
        one.limits.max_bytes = bytes;
        if one.chromatograms(1).is_ok() {
            minimum = Some(bytes);
            break;
        }
    }
    two.limits.max_bytes = minimum.unwrap();
    assert!(two.chromatograms(1).is_err());
}

#[test]
fn invalid_consumed_values_are_checked_without_touching_borrowed_traces() {
    let good = trace("good", 1., &[1., 3., 1.]);
    let mut invalid = trace("bad", 2., &[1., 2.]);
    invalid.peaks_mut()[1].set_rt(f64::NAN);
    let h = hypothesis([&good, &invalid]);
    assert!(h.convex_hulls().is_err());
    assert!(h.chromatograms(1).is_err());
    assert_eq!(h.labels().unwrap(), ["good", "bad"]);
    assert_eq!(h.centroid_mz().unwrap(), 1.);
    assert_eq!(h.number_of_feature_points(), 5);
    assert!(invalid[1].rt().is_nan());
    assert_eq!(good[0].intensity, 1.);
    let mut a = trace("a", -f64::MAX, &[1.]);
    let mut b = trace("b", f64::MAX, &[1.]);
    a.set_quant_method(Quant::MaxHeight);
    b.set_quant_method(Quant::MaxHeight);
    assert!(hypothesis([&a, &b]).isotope_distances().is_err());
}
