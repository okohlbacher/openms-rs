// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::{
    analysis::{
        elution_peak_detection::ElutionPeakDetection,
        feature_finding_metabo::{
            FeatureFindingMetabo, FeatureFindingMetaboOptions, IsotopeFilteringModel,
        },
        mass_trace_detection::MassTraceDetection,
    },
    concept::{progress_logger::ProgressLogType, unique_id::UniqueIdGenerator},
    kernel::{DataArray, FeatureMap, MSExperiment, MSSpectrum, MassTrace, Peak1D},
};
fn finder(options: FeatureFindingMetaboOptions) -> FeatureFindingMetabo {
    let mut f = FeatureFindingMetabo::with_options(options).unwrap();
    f.logger.set_log_type(ProgressLogType::None);
    f
}
fn source_input(text: &str) -> MSExperiment {
    let mut input = MSExperiment::default();
    for line in text.lines().filter(|s| !s.starts_with('#')) {
        let fields: Vec<_> = line.split('\t').collect();
        let index: usize = fields[0].parse().unwrap();
        if index == input.spectra.len() {
            let mut spectrum = MSSpectrum {
                rt: fields[1].parse().unwrap(),
                ms_level: fields[2].parse().unwrap(),
                ..Default::default()
            };
            if fields[5] != "." {
                spectrum
                    .float_data_arrays
                    .push(DataArray::new(fields[5], Vec::new()));
            }
            input.spectra.push(spectrum);
        }
        if fields[3] != "." {
            input.spectra[index].peaks.push(Peak1D::new(
                fields[3].parse().unwrap(),
                fields[4].parse().unwrap(),
            ));
            if fields[5] != "." {
                input.spectra[index].float_data_arrays[0]
                    .data
                    .push(fields[6].parse().unwrap());
            }
        }
    }
    input
}
fn source_traces(im: bool) -> Vec<MassTrace> {
    let input = source_input(if im {
        include_str!("data/feature_finding_metabo_im_source.tsv")
    } else {
        include_str!("data/feature_finding_metabo_source.tsv")
    });
    assert_eq!(input.spectra.len(), if im { 31 } else { 360 });
    assert_eq!(
        input.spectra.iter().map(|s| s.peaks.len()).sum::<usize>(),
        if im { 4873 } else { 11646 }
    );
    let mut detection = MassTraceDetection::new();
    detection.logger.set_log_type(ProgressLogType::None);
    if im {
        detection.options.mass_error_ppm = 10.0;
        detection.options.noise_threshold_int = 10.0;
        detection.options.ion_mobility_tolerance = 0.01;
    }
    let mut traces = detection.run(&input, 0).unwrap();
    assert!(!traces.is_empty());
    let mut split = ElutionPeakDetection::new();
    split.logger.set_log_type(ProgressLogType::None);
    split.detect_peaks_many(&mut traces).unwrap()
}
fn close(actual: f64, expected: f64) {
    // Exact source comparator contract: absolute1 OR multiplicative1.001.
    let ratio = actual / expected;
    let ratio_ok = actual != 0.0
        && expected != 0.0
        && actual.is_sign_positive() == expected.is_sign_positive()
        && ratio.max(1.0 / ratio) <= 1.001;
    assert!(
        (actual - expected).abs() <= 1.0 || ratio_ok,
        "{actual:.17} != {expected:.17}"
    );
}
fn expected_output(map: &FeatureMap) {
    use std::collections::BTreeMap;
    let mut expected_hulls = BTreeMap::<(usize, usize), Vec<(f64, f64)>>::new();
    for line in include_str!("data/feature_finding_metabo_expected.tsv")
        .lines()
        .filter(|s| !s.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        let index: usize = fields[1].parse().unwrap();
        let feature = &map.features[index];
        match fields[0] {
            "F" => {
                close(feature.rt, fields[2].parse().unwrap());
                close(feature.mz, fields[3].parse().unwrap());
                close(f64::from(feature.intensity), fields[4].parse().unwrap());
                close(f64::from(feature.quality), fields[5].parse().unwrap());
                assert_eq!(feature.charge, fields[6].parse::<i32>().unwrap());
            }
            "M" => {
                let value = &feature.metadata[fields[2]];
                match fields[3] {
                    "string" => assert_eq!(value.as_str().unwrap(), fields[4]),
                    "int" => assert_eq!(value.as_i64().unwrap(), fields[4].parse::<i64>().unwrap()),
                    "float" => close(value.as_f64().unwrap(), fields[4].parse().unwrap()),
                    "floatList" => {
                        let text = fields[4].trim_matches(['[', ']']);
                        let expected: Vec<f64> = if text.is_empty() {
                            Vec::new()
                        } else {
                            text.split(',').map(|v| v.trim().parse().unwrap()).collect()
                        };
                        let actual = value.as_float_list().unwrap();
                        assert_eq!(actual.len(), expected.len());
                        for (&a, b) in actual.iter().zip(expected) {
                            close(a, b);
                        }
                    }
                    _ => panic!("unknown source metadata type"),
                }
            }
            "H" => expected_hulls
                .entry((index, fields[2].parse().unwrap()))
                .or_default()
                .push((fields[3].parse().unwrap(), fields[4].parse().unwrap())),
            _ => panic!("unknown projection row"),
        }
    }
    for ((feature, hull), points) in expected_hulls {
        // Source FeatureXMLHandler compresses a copy before writing the golden.
        let mut native_hull = map.features[feature].convex_hulls[hull].clone();
        native_hull.compress();
        let actual = native_hull.hull_points();
        assert_eq!(actual.len(), points.len(), "feature{feature} hull{hull}");
        for (a, b) in actual.iter().zip(points) {
            close(a.rt, b.0);
            close(a.mz, b.1);
        }
    }
}
#[test]
fn exact_source_chain_counts_c13_kenar_elements_and_all_expected_feature_fields() {
    let traces = source_traces(false);
    let mut rng = UniqueIdGenerator::from_seed(123);
    for (c13, elements, count) in [(true, false, 83), (false, false, 81), (false, true, 80)] {
        let mut input = traces.clone();
        let mut f = finder(FeatureFindingMetaboOptions {
            mz_scoring_13c: c13,
            mz_scoring_by_elements: elements,
            report_convex_hulls: true,
            ..Default::default()
        });
        let output = f.run(&mut input, &mut rng).unwrap();
        assert_eq!(
            output.features.features.len(),
            count,
            "c13={c13}, elements={elements}"
        );
        if !c13 && !elements {
            expected_output(&output.features);
        }
        assert!(output.chromatograms.is_empty());
    }
}
#[test]
fn exact_source_im_chain_retains_im_feature_metadata() {
    let mut traces = source_traces(true);
    assert!(traces.iter().any(MassTrace::contains_im_data));
    let mut f = finder(FeatureFindingMetaboOptions {
        mz_scoring_13c: true,
        ..Default::default()
    });
    let output = f
        .run(&mut traces, &mut UniqueIdGenerator::from_seed(1))
        .unwrap();
    assert!(!output.features.features.is_empty());
    assert!(f.contains_im_data());
    assert!(
        output
            .features
            .features
            .iter()
            .all(|f| f.metadata.contains_key("masstrace_centroid_im"))
    );
}
#[test]
fn default_nineteen_options_and_effective_source_override() {
    let mut f = finder(FeatureFindingMetaboOptions::default());
    assert_eq!(f.options(), &FeatureFindingMetaboOptions::default());
    let mut options = f.options().clone();
    options.mz_scoring_by_elements = true;
    f.set_options(options.clone()).unwrap();
    let mut rng = UniqueIdGenerator::from_seed(1);
    let mut before = rng.clone();
    let output = f.run(&mut [], &mut rng).unwrap();
    assert!(output.features.features.is_empty());
    assert_eq!(output.features.unique_id, 0);
    assert_eq!(rng.get_unique_id(), before.get_unique_id());
    assert_eq!(
        f.effective_isotope_filtering_model(),
        IsotopeFilteringModel::None
    );
    assert_eq!(
        f.options().isotope_filtering_model,
        IsotopeFilteringModel::Metabolites5
    );
    assert!(output.diagnostics.isotope_filtering_disabled_by_elements);
    options.mz_scoring_by_elements = false;
    f.set_options(options).unwrap();
    assert_eq!(
        f.effective_isotope_filtering_model(),
        IsotopeFilteringModel::Metabolites5
    );
}

fn small_trace(label: &str, mz: f64, height: f32) -> MassTrace {
    use openms::kernel::{MassTraceQuantMethod, Peak2D};
    let mut t = MassTrace::from_peaks(vec![
        Peak2D::new(0., mz, 1.),
        Peak2D::new(1., mz, height),
        Peak2D::new(2., mz, 1.),
    ])
    .unwrap();
    t.set_label(label).unwrap();
    t.update_median_mz().unwrap();
    t.update_median_rt().unwrap();
    t.estimate_fwhm(false).unwrap();
    t.set_quant_method(MassTraceQuantMethod::MaxHeight);
    t.set_smoothed_intensities(&[2., f64::from(height) * 2., 2.])
        .unwrap();
    t
}
fn unfiltered() -> FeatureFindingMetaboOptions {
    FeatureFindingMetaboOptions {
        isotope_filtering_model: IsotopeFilteringModel::None,
        charge_lower_bound: 1,
        charge_upper_bound: 1,
        ..Default::default()
    }
}
#[test]
fn independent_isotope_pair_quantification_charge_metadata_and_raw_chromatograms() {
    let a = small_trace("mono", 100., 10.);
    let b = small_trace("iso", 101.001948, 5.);
    for (use_smooth, report_smooth, summed, expected) in [
        (false, true, false, 10.),
        (true, false, false, 10.),
        (true, true, false, 20.),
        (true, true, true, 30.),
    ] {
        let mut f = finder(FeatureFindingMetaboOptions {
            use_smoothed_intensities: use_smooth,
            report_smoothed_intensities: report_smooth,
            report_summed_intensities: summed,
            report_chromatograms: true,
            report_convex_hulls: true,
            ..unfiltered()
        });
        let mut traces = [b.clone(), a.clone()];
        let result = f
            .run(&mut traces, &mut UniqueIdGenerator::from_seed(1))
            .unwrap();
        assert_eq!(result.features.features.len(), 1);
        let feature = &result.features.features[0];
        assert_eq!(feature.intensity, expected);
        assert_eq!(feature.charge, 1);
        assert_eq!(feature.metadata["label"].as_str().unwrap(), "mono_iso");
        assert_eq!(feature.metadata["num_of_masstraces"].as_i64().unwrap(), 2);
        assert_eq!(
            feature.metadata["legal_isotope_pattern"].as_i64().unwrap(),
            -1
        );
        // Raw RT cosine uses [1,10,1] and [1,5,1], regardless of smoothing.
        let cosine = 52.0 / (102.0f64.sqrt() * 27.0f64.sqrt());
        let expected_score = 2.0 / 3.0 + cosine / 3.0;
        assert!((f64::from(feature.quality) - expected_score).abs() < 1e-6);
        assert_eq!(feature.convex_hulls.len(), 2);
        assert_eq!(result.chromatograms.len(), 1);
        assert_eq!(result.chromatograms[0].len(), 2);
        assert_eq!(result.chromatograms[0][0].peaks[1].intensity, 10.);
        assert_eq!(result.chromatograms[0][1].peaks[1].intensity, 5.);
        assert_eq!(result.chromatograms[0][1].precursor.mz, 100.);
        assert_eq!(traces[0].label(), "mono");
        assert_eq!(
            result.diagnostics.smoothed_reporting_disabled,
            !use_smooth && report_smooth
        );
    }
}
#[test]
fn source_acceptance_order_differs_from_final_map_order_and_rng_draws_follow_acceptance() {
    let options = FeatureFindingMetaboOptions {
        charge_lower_bound: 2,
        charge_upper_bound: 1,
        report_chromatograms: true,
        ..unfiltered()
    };
    let mut f = finder(options);
    let mut traces = [
        small_trace("high", 500., 30.),
        small_trace("low", 100., 10.),
        small_trace("mid", 300., 20.),
    ];
    let mut rng = UniqueIdGenerator::from_seed(47);
    let mut oracle = rng.clone();
    let (high, mid, low, map) = (
        oracle.get_unique_id(),
        oracle.get_unique_id(),
        oracle.get_unique_id(),
        oracle.get_unique_id(),
    );
    let result = f.run(&mut traces, &mut rng).unwrap();
    assert_eq!(
        result
            .features
            .features
            .iter()
            .map(|f| f.unique_id)
            .collect::<Vec<_>>(),
        [low, mid, high]
    );
    assert_eq!(result.features.unique_id, map);
    assert_eq!(rng.get_unique_id(), oracle.get_unique_id());
    assert_eq!(
        result
            .chromatograms
            .iter()
            .map(|g| g[0].precursor.mz)
            .collect::<Vec<_>>(),
        [500., 300., 100.]
    );
    assert_eq!(result.chromatograms[0][0].native_id, format!("{high}_0"));
    assert_eq!(
        traces
            .iter()
            .map(MassTrace::centroid_mz)
            .collect::<Vec<_>>(),
        [100., 300., 500.]
    );
}
#[test]
fn labels_define_exclusion_and_stable_ties_keep_serial_sorted_seed_order() {
    let mut f = finder(FeatureFindingMetaboOptions {
        charge_lower_bound: 2,
        charge_upper_bound: 1,
        ..unfiltered()
    });
    let mut traces = [
        small_trace("same", 300., 10.),
        small_trace("same", 100., 10.),
        small_trace("same", 200., 10.),
    ];
    let output = f
        .run(&mut traces, &mut UniqueIdGenerator::from_seed(1))
        .unwrap();
    assert_eq!(output.features.features.len(), 1);
    assert_eq!(output.features.features[0].mz, 100.);
    assert_eq!(output.features.features[0].charge, 0);
    let mut options = f.options().clone();
    options.remove_single_traces = true;
    f.set_options(options).unwrap();
    assert!(
        f.run(&mut traces, &mut UniqueIdGenerator::from_seed(1))
            .unwrap()
            .features
            .features
            .is_empty()
    );
}
#[test]
fn local_windows_zero_charge_and_element_alphabets_keep_source_defined_branches() {
    let original = [
        small_trace("a", 100., 10.),
        small_trace("b", 101.001948, 5.),
    ];
    for options in [
        FeatureFindingMetaboOptions {
            local_rt_range: -1.,
            enable_rt_filtering: false,
            ..unfiltered()
        },
        FeatureFindingMetaboOptions {
            local_mz_range: 0.,
            ..unfiltered()
        },
        FeatureFindingMetaboOptions {
            charge_lower_bound: 0,
            charge_upper_bound: 0,
            ..unfiltered()
        },
        FeatureFindingMetaboOptions {
            charge_lower_bound: 2,
            charge_upper_bound: 1,
            chrom_fwhm: f64::NAN,
            ..unfiltered()
        },
        FeatureFindingMetaboOptions {
            elements: "(13)C".into(),
            mz_scoring_by_elements: true,
            ..unfiltered()
        },
    ] {
        let mut f = finder(options);
        assert_eq!(
            f.run(&mut original.clone(), &mut UniqueIdGenerator::from_seed(1))
                .unwrap()
                .features
                .features
                .len(),
            2
        );
    }
    // The empty isotope window is still computed but unused by mean scoring.
    let mut f = finder(FeatureFindingMetaboOptions {
        elements: "".into(),
        ..unfiltered()
    });
    assert_eq!(
        f.run(&mut original.clone(), &mut UniqueIdGenerator::from_seed(1))
            .unwrap()
            .features
            .features
            .len(),
        1
    );
    let before = f.options().clone();
    let mut bad = before.clone();
    bad.elements = "NotAnElement".into();
    assert!(f.set_options(bad).is_err());
    assert_eq!(f.options(), &before);
    let mut bad = before.clone();
    bad.min_isotope_rt_overlap = 1.1;
    assert!(f.set_options(bad).is_err());
    assert_eq!(f.options(), &before);
}
#[test]
fn mixed_im_gates_candidates_and_flags_update_only_on_success() {
    let a = small_trace("a", 100., 10.);
    let mut b = small_trace("b", 101.001948, 5.);
    b.set_centroid_im(1.).unwrap();
    let mut f = finder(unfiltered());
    let output = f
        .run(
            &mut [a.clone(), b.clone()],
            &mut UniqueIdGenerator::from_seed(1),
        )
        .unwrap();
    assert_eq!(output.features.features.len(), 2);
    assert!(f.contains_im_data());
    assert_eq!(
        output.features.features[0].metadata["masstrace_centroid_im"]
            .as_float_list()
            .unwrap(),
        [0.]
    );
    assert_eq!(
        output.features.features[1].metadata["masstrace_centroid_im"]
            .as_float_list()
            .unwrap(),
        [1.]
    );
    let mut options = f.options().clone();
    options.local_im_range = 1.;
    f.set_options(options).unwrap();
    assert_eq!(
        f.run(&mut [a.clone(), b], &mut UniqueIdGenerator::from_seed(1))
            .unwrap()
            .features
            .features
            .len(),
        1
    );
    f.run(&mut [], &mut UniqueIdGenerator::from_seed(1))
        .unwrap();
    assert!(f.contains_im_data()); // source empty run leaves flag
    f.run(&mut [a], &mut UniqueIdGenerator::from_seed(1))
        .unwrap();
    assert!(!f.contains_im_data());
}
#[test]
fn complete_output_replacement_returns_previous_ownership_and_ignores_unused_old_payload() {
    use openms::kernel::{Feature, MSChromatogram};
    let mut f = finder(unfiltered());
    let mut map = FeatureMap {
        features: vec![Feature::new(f64::NAN, 0., f32::NAN)],
        identifier: "old".into(),
        loaded_file_path: "source".into(),
        ..Default::default()
    };
    let ptr = map.features.as_ptr();
    let mut chroms = vec![vec![MSChromatogram {
        native_id: "old".into(),
        ..Default::default()
    }]];
    let chrom_ptr = chroms.as_ptr();
    f.limits.max_bytes = 0;
    f.limits.max_work = 0;
    let (old, old_chroms, _) = f
        .run_into(
            &mut [],
            &mut UniqueIdGenerator::from_seed(1),
            &mut map,
            &mut chroms,
        )
        .unwrap();
    assert_eq!(old.features.as_ptr(), ptr);
    assert_eq!(old_chroms.as_ptr(), chrom_ptr);
    assert_eq!(old.identifier, "old");
    assert_eq!(old.loaded_file_path, "source");
    assert!(old.features[0].rt.is_nan());
    assert_eq!(map, FeatureMap::default());
    assert!(chroms.is_empty());
}
#[test]
fn late_report_failure_rolls_back_input_output_rng_and_effective_flags() {
    use openms::kernel::Feature;
    let good = small_trace("good", 100., 20.);
    let mut bad = small_trace("bad", 300., 10.);
    bad.peaks_mut()[1].intensity = f32::NAN;
    bad.set_centroid_im(1.).unwrap();
    let mut input = [bad, good];
    let before = input.clone();
    let mut f = finder(FeatureFindingMetaboOptions {
        charge_lower_bound: 2,
        charge_upper_bound: 1,
        report_smoothed_intensities: false,
        mz_scoring_by_elements: true,
        ..Default::default()
    });
    let mut output = FeatureMap {
        features: vec![Feature::new(1., 2., 3.)],
        ..Default::default()
    };
    let saved = output.clone();
    let mut chroms = Vec::new();
    let mut rng = UniqueIdGenerator::from_seed(97);
    let mut oracle = rng.clone();
    assert!(
        f.run_into(&mut input, &mut rng, &mut output, &mut chroms)
            .is_err()
    );
    assert_eq!(rng.get_unique_id(), oracle.get_unique_id());
    assert_eq!(output, saved);
    assert!(chroms.is_empty());
    assert_eq!(input[0].label(), before[0].label());
    assert!(input[0].peaks()[1].intensity.is_nan());
    assert_eq!(input[1], before[1]);
    assert!(!f.contains_im_data());
    assert_eq!(
        f.effective_isotope_filtering_model(),
        IsotopeFilteringModel::Metabolites5
    );
}
#[test]
fn zero_total_and_count_byte_work_limits_are_checked_atomically() {
    let mut original = [small_trace("b", 300., 10.), small_trace("a", 100., 20.)];
    let mut f = finder(unfiltered());
    for limits in [
        openms::analysis::feature_finding_metabo::FeatureFindingMetaboLimits {
            max_traces: 1,
            ..Default::default()
        },
        openms::analysis::feature_finding_metabo::FeatureFindingMetaboLimits {
            max_peaks: 5,
            ..Default::default()
        },
        openms::analysis::feature_finding_metabo::FeatureFindingMetaboLimits {
            max_hypotheses: 1,
            ..Default::default()
        },
        openms::analysis::feature_finding_metabo::FeatureFindingMetaboLimits {
            max_features: 1,
            ..Default::default()
        },
        openms::analysis::feature_finding_metabo::FeatureFindingMetaboLimits {
            max_bytes: 0,
            ..Default::default()
        },
        openms::analysis::feature_finding_metabo::FeatureFindingMetaboLimits {
            max_work: 0,
            ..Default::default()
        },
    ] {
        f.limits = limits;
        let before = original.clone();
        let mut rng = UniqueIdGenerator::from_seed(1);
        let mut copy = rng.clone();
        assert!(f.run(&mut original, &mut rng).is_err());
        assert_eq!(original, before);
        assert_eq!(rng.get_unique_id(), copy.get_unique_id());
    }
    f.limits = Default::default();
    let mut empty = MassTrace::new();
    empty.set_label("zero").unwrap();
    assert!(
        f.run(&mut [empty], &mut UniqueIdGenerator::from_seed(1))
            .is_err()
    );
}
#[test]
fn progress_end_error_prevents_all_scientific_publication() {
    use openms::concept::progress_logger::{
        ProgressBackend, ProgressLogger, ProgressNesting, ProgressTime,
    };
    use std::sync::Arc;
    struct Fail;
    impl ProgressBackend for Fail {
        fn start_progress(&mut self, _: i64, _: i64, _: &str, _: usize) -> openms::Result<()> {
            Ok(())
        }
        fn set_progress(&mut self, _: i64, _: usize) -> openms::Result<()> {
            Ok(())
        }
        fn next_progress(&mut self) -> openms::Result<i64> {
            unreachable!()
        }
        fn end_progress(&mut self, _: usize, _: u64) -> openms::Result<()> {
            Err(openms::Error::Unsupported("end".into()))
        }
    }
    let nesting = ProgressNesting::default();
    let mut f = finder(unfiltered());
    f.logger = ProgressLogger::with_clock_and_nesting(
        Arc::new(|| {
            Ok(ProgressTime {
                wall_second: 1,
                wall_seconds: 1.,
                cpu_seconds: None,
            })
        }),
        nesting.clone(),
    );
    f.logger.set_logger(Box::new(Fail));
    let mut input = [small_trace("b", 300., 10.), small_trace("a", 100., 20.)];
    let before = input.clone();
    let mut rng = UniqueIdGenerator::from_seed(1);
    let mut saved = rng.clone();
    assert!(f.run(&mut input, &mut rng).is_err());
    assert_eq!(input, before);
    assert_eq!(rng.get_unique_id(), saved.get_unique_id());
    assert_eq!(nesting.depth(), 0);
}

#[test]
fn both_model_choices_and_mass_cap_match_executed_libsvm_rows() {
    // External LIBSVM fixture rows: (600,.05,.03,.0125) and
    // (1000,.2,.12,.05) choose label1 for2%, label2 for5%.
    for (mz, ratios) in [(600., [0.05, 0.03, 0.0125]), (1100., [0.2, 0.12, 0.05])] {
        for (model, accepted) in [
            (IsotopeFilteringModel::Metabolites2, false),
            (IsotopeFilteringModel::Metabolites5, true),
        ] {
            let mut traces = vec![small_trace("mono", mz, 1000.)];
            for (index, ratio) in ratios.into_iter().enumerate() {
                let offset = index as f64 + 1.;
                traces.push(small_trace(
                    &format!("iso{index}"),
                    mz + 1.000857 * offset + 0.001091,
                    (1000. * ratio) as f32,
                ));
            }
            let mut f = finder(FeatureFindingMetaboOptions {
                isotope_filtering_model: model,
                use_smoothed_intensities: false,
                enable_rt_filtering: false,
                local_mz_range: 3.5,
                ..unfiltered()
            });
            let output = f
                .run(&mut traces, &mut UniqueIdGenerator::from_seed(3))
                .unwrap();
            let full = output
                .features
                .features
                .iter()
                .find(|feature| feature.metadata["num_of_masstraces"].as_i64().unwrap() == 4);
            assert_eq!(full.is_some(), accepted, "{mz} {model:?}");
            if let Some(feature) = full {
                assert_eq!(
                    feature.metadata["legal_isotope_pattern"].as_i64().unwrap(),
                    1
                );
            }
        }
    }
}

#[test]
fn cumulative_limits_cover_many_individually_valid_output_records() {
    for bytes in [false, true] {
        let mut f = finder(unfiltered());
        if bytes {
            f.limits.max_bytes = 100_000;
        } else {
            f.limits.max_work = 100_000;
        }
        let inputs: Vec<_> = (0..10)
            .map(|i| small_trace(&format!("trace{i}"), 100. + 100. * i as f64, 10.))
            .collect();
        for trace in &inputs {
            assert!(
                f.run(&mut [trace.clone()], &mut UniqueIdGenerator::from_seed(1))
                    .is_ok()
            );
        }
        let mut inputs = inputs;
        inputs.reverse();
        let before = inputs.clone();
        let mut rng = UniqueIdGenerator::from_seed(1);
        let mut expected_rng = rng.clone();
        assert!(f.run(&mut inputs, &mut rng).is_err());
        assert_eq!(inputs, before);
        assert_eq!(rng.get_unique_id(), expected_rng.get_unique_id());
    }
}

#[test]
fn peptide_candidate_uses_raw_previous_members_and_configured_new_intensity() {
    // source-expression peptide scores at candidate mass1000: [10,5] vs[10,10].
    for (smooth, peptide_score) in [(false, 0.9995951520171115), (true, 0.9572966255102681)] {
        let mut traces = [
            small_trace("mono", 1000. - 1.001948, 10.),
            small_trace("iso", 1000., 5.),
        ];
        let mut f = finder(FeatureFindingMetaboOptions {
            isotope_filtering_model: IsotopeFilteringModel::Peptides,
            use_smoothed_intensities: smooth,
            ..unfiltered()
        });
        let output = f
            .run(&mut traces, &mut UniqueIdGenerator::from_seed(5))
            .unwrap();
        assert_eq!(output.features.features.len(), 1);
        let feature = &output.features.features[0];
        let rt_cosine = 52. / (102.0f64.sqrt() * 27.0f64.sqrt());
        let expected_score = 2. / 3. + rt_cosine * peptide_score / 3.;
        assert!((f64::from(feature.quality) - expected_score).abs() < 1e-7);
        assert_eq!(
            feature.metadata["legal_isotope_pattern"].as_i64().unwrap(),
            -1
        );
        assert_eq!(feature.intensity, if smooth { 20. } else { 10. });
    }
}
