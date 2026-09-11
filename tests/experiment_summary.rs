use openms::kernel::{
    ChromatogramPeak, DataArray, ExperimentRanges, MSChromatogram, MSExperiment, MSSpectrum,
    NumericRange, Peak1D, SummaryLimits,
};
use openms::{Error, metadata::MetaValue};

fn spectrum(rt: f64, level: u32, intensities: &[f32]) -> MSSpectrum {
    MSSpectrum {
        rt,
        ms_level: level,
        peaks: intensities
            .iter()
            .enumerate()
            .map(|(i, &v)| Peak1D::new(i as f64, v))
            .collect(),
        ..Default::default()
    }
}
fn chrom(mz: f64, points: &[(f64, f32)]) -> MSChromatogram {
    let mut c = MSChromatogram::from_peaks(
        points
            .iter()
            .map(|&(rt, v)| ChromatogramPeak::new(rt, v))
            .collect(),
    );
    c.product.mz = mz;
    c
}
fn source_tic() -> MSExperiment {
    MSExperiment {
        spectra: vec![
            spectrum(0.0, 1, &[3.0, 5.0]),
            spectrum(2.0, 1, &[2.0]),
            spectrum(2.0, 2, &[0.5]),
            spectrum(5.0, 1, &[2.0, 3.0, 4.0]),
        ],
        ..Default::default()
    }
}
fn bounds(min: f64, max: f64) -> Option<NumericRange> {
    Some(NumericRange { min, max })
}
fn tiny() -> SummaryLimits {
    SummaryLimits {
        max_work: 0,
        max_bytes: 0,
        max_output_points: 0,
    }
}

#[test]
fn source_tic_literals_cover_negative_zero_and_two_positive_bins() {
    // Literal source MSExperiment_test.cpp calculateTIC fixture and assertions.
    let exp = source_tic();
    for bin in [-1.0, -0.0, 0.0] {
        assert_eq!(
            exp.calculate_tic_binned(bin, 1).unwrap().peaks,
            vec![
                ChromatogramPeak::new(0.0, 8.0),
                ChromatogramPeak::new(2.0, 2.0),
                ChromatogramPeak::new(5.0, 9.0)
            ]
        );
    }
    assert_eq!(
        exp.calculate_tic_binned(2.0, 1).unwrap().peaks,
        vec![
            ChromatogramPeak::new(0.0, 8.0),
            ChromatogramPeak::new(2.0, 2.0),
            ChromatogramPeak::new(4.0, 4.5),
            ChromatogramPeak::new(6.0, 4.5)
        ]
    );
    let wide = exp.calculate_tic_binned(6.0, 1).unwrap();
    let left_after_second = (8.0_f64 + 2.0 * 4.0 / 6.0) as f32;
    let right_after_second = (2.0_f64 * 2.0 / 6.0) as f32;
    assert_eq!(
        wide.peaks,
        vec![
            ChromatogramPeak::new(0.0, (f64::from(left_after_second) + 9.0 * 1.0 / 6.0) as f32),
            ChromatogramPeak::new(
                6.0,
                (f64::from(right_after_second) + 9.0 * 5.0 / 6.0) as f32
            ),
        ]
    );
    assert_eq!(
        exp.calculate_tic_binned(2.0, 2).unwrap().peaks,
        vec![ChromatogramPeak::new(2.0, 0.5)]
    );
    assert_eq!(exp.calculate_tic_binned(0.0, 0).unwrap().len(), 4);
    assert_eq!(
        exp.calculate_tic_binned(2.0, 0).unwrap().peaks[1].intensity,
        2.5
    );
    assert!(exp.calculate_tic_binned(2.0, 9).unwrap().is_empty());
}

#[test]
fn source_grid_ceil_addition_order_and_f32_scan_accumulation() {
    // ceil(q+1), not ceil(q)+1: q immediately above 1 rounds back to 2
    // after adding 1. Source therefore emits two points and clamps the tail.
    let above_one = f64::from_bits(1.0_f64.to_bits() + 1);
    let exp = MSExperiment {
        spectra: vec![
            spectrum(0.0, 1, &[16_777_216.0, 1.0, -16_777_216.0]),
            spectrum(above_one, 1, &[7.0]),
        ],
        ..Default::default()
    };
    assert_eq!(
        exp.calculate_tic_binned(1.0, 1).unwrap().peaks,
        vec![
            ChromatogramPeak::new(0.0, 0.0),
            ChromatogramPeak::new(1.0, 7.0)
        ]
    );
}

#[test]
fn empty_singleton_and_duplicate_rt_grids_preserve_source_boundaries() {
    let empty = MSExperiment::new();
    for bin in [-1.0, 0.0, 1.0] {
        assert!(
            empty
                .calculate_tic_binned_with_limits(bin, 1, tiny())
                .unwrap()
                .is_empty()
        );
    }
    let mut exp = MSExperiment {
        spectra: vec![
            spectrum(-1.0, 1, &[]),
            spectrum(-1.0, 1, &[2.0]),
            spectrum(-1.0, 1, &[-4.0]),
        ],
        ..Default::default()
    };
    assert_eq!(
        exp.calculate_tic_binned(0.5, 1).unwrap().peaks,
        vec![ChromatogramPeak::new(-1.0, -2.0)]
    );
    exp.spectra.truncate(1);
    assert_eq!(
        exp.calculate_tic_binned(100.0, 1).unwrap().peaks,
        vec![ChromatogramPeak::new(-1.0, 0.0)]
    );
}

#[test]
fn tic_checks_only_consumed_fields_and_positive_bin_sorting() {
    let mut exp = source_tic();
    exp.spectra[0].peaks[0].mz = f64::NAN; // TIC never reads m/z.
    exp.spectra[0]
        .string_data_arrays
        .push(DataArray::new("unaligned", vec!["kept".into()]));
    exp.spectra[2].rt = f64::NAN; // Excluded MS2.
    exp.spectra[2].peaks[0].intensity = f32::NAN;
    assert!(exp.calculate_tic_binned(2.0, 1).is_ok());
    exp.spectra.swap(0, 3);
    let unbinned = exp.calculate_tic_binned(-1.0, 1).unwrap();
    assert_eq!(
        unbinned.peaks.iter().map(|p| p.rt).collect::<Vec<_>>(),
        vec![5.0, 2.0, 0.0]
    );
    assert!(matches!(
        exp.calculate_tic_binned(1.0, 1),
        Err(Error::UnsortedData)
    ));
    assert!(exp.calculate_tic_binned(0.0, 0).is_err());
}

#[test]
fn tic_nonfinite_overflow_stagnation_and_cumulative_limits_are_checked() {
    let exp = source_tic();
    for bin in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
        assert!(exp.calculate_tic_binned(bin, 1).is_err());
    }
    for limits in [
        SummaryLimits {
            max_work: 4,
            ..Default::default()
        },
        SummaryLimits {
            max_bytes: 1,
            ..Default::default()
        },
        SummaryLimits {
            max_output_points: 6,
            ..Default::default()
        }, // 3 raw +4 raster points.
    ] {
        assert!(
            exp.calculate_tic_binned_with_limits(2.0, 1, limits)
                .is_err()
        );
    }
    let overflow = MSExperiment {
        spectra: vec![spectrum(0.0, 1, &[f32::MAX, f32::MAX])],
        ..Default::default()
    };
    assert!(overflow.calculate_tic_binned(0.0, 1).is_err());
    let mut grid = MSExperiment {
        spectra: vec![
            spectrum(1e20, 1, &[1.0]),
            spectrum(1e20 + 16384.0, 1, &[1.0]),
        ],
        ..Default::default()
    };
    assert!(grid.calculate_tic_binned(1.0, 1).is_err()); // First grid step rounds to same RT.
    grid.spectra[0].rt = -f64::MAX;
    grid.spectra[1].rt = f64::MAX;
    assert!(grid.calculate_tic_binned(f32::MIN_POSITIVE, 1).is_err());
    assert_eq!(source_tic(), exp); // Generation only borrows its input.
}

#[test]
fn source_chromatogram_and_combined_range_literals() {
    let exp = MSExperiment {
        spectra: vec![
            MSSpectrum {
                peaks: vec![Peak1D::new(100.0, 1000.0)],
                rt: 30.0,
                ..Default::default()
            },
            MSSpectrum {
                peaks: vec![Peak1D::new(200.0, 2000.0)],
                rt: 35.0,
                ms_level: 2,
                ..Default::default()
            },
        ],
        chromatograms: vec![
            chrom(305.0, &[(10.0, 500.0), (20.0, 1500.0)]),
            chrom(405.0, &[(15.0, 800.0), (25.0, 1800.0)]),
        ],
        ..Default::default()
    };
    assert_eq!(
        exp.chromatogram_ranges().unwrap(),
        ExperimentRanges {
            rt: bounds(10.0, 25.0),
            mz: bounds(305.0, 405.0),
            intensity: bounds(500.0, 1800.0)
        }
    );
    assert_eq!(
        exp.combined_ranges().unwrap(),
        ExperimentRanges {
            rt: bounds(10.0, 35.0),
            mz: bounds(100.0, 405.0),
            intensity: bounds(500.0, 2000.0)
        }
    );
}

#[test]
fn ranges_include_empty_anchors_and_ignore_unconsumed_metadata() {
    let mut exp = MSExperiment {
        spectra: vec![spectrum(-10.0, 1, &[])],
        chromatograms: vec![chrom(700.0, &[])],
        ..Default::default()
    };
    exp.chromatograms[0].product.isolation_window_lower_offset = f64::NAN;
    exp.chromatograms[0]
        .metadata
        .insert("product_mz".into(), "9999".into());
    exp.chromatograms[0]
        .product
        .cv_terms
        .metadata
        .insert("owned".into(), MetaValue::from("retained"));
    assert_eq!(
        exp.chromatogram_ranges().unwrap(),
        ExperimentRanges {
            mz: bounds(700.0, 700.0),
            ..Default::default()
        }
    );
    assert_eq!(
        exp.combined_ranges().unwrap(),
        ExperimentRanges {
            rt: bounds(-10.0, -10.0),
            mz: bounds(700.0, 700.0),
            intensity: None
        }
    );
    exp.chromatograms[0].peaks = vec![
        ChromatogramPeak::new(3.0, -9.0),
        ChromatogramPeak::new(1.0, -1.0),
    ];
    assert_eq!(exp.chromatogram_ranges().unwrap().rt, bounds(1.0, 3.0));
    assert_eq!(exp.combined_ranges().unwrap().intensity, bounds(-9.0, -1.0));
    assert!(
        exp.combined_ranges_with_limits(SummaryLimits {
            max_work: 3,
            ..Default::default()
        })
        .is_err()
    ); // Shared across both data kinds.
    exp.chromatograms[0].product.mz = f64::NAN;
    assert!(exp.chromatogram_ranges().is_err());
    assert_eq!(
        MSExperiment::new()
            .combined_ranges_with_limits(tiny())
            .unwrap(),
        ExperimentRanges::default()
    );
}

#[test]
fn source_chromatogram_sort_literal_and_stable_aligned_three_cycle() {
    // Source literal outer100/80 and inner0.3/0.2/0.1 ordering.
    let mut exp = MSExperiment {
        chromatograms: vec![
            chrom(100.0, &[(0.3, 10.0), (0.2, 10.2)]),
            chrom(80.0, &[(0.2, 10.2), (0.1, 10.4)]),
        ],
        ..Default::default()
    };
    exp.sort_chromatograms(false).unwrap();
    assert_eq!(exp.chromatograms[0].product.mz, 80.0);
    assert_eq!(exp.chromatograms[1].peaks[0].rt, 0.3);
    exp.sort_chromatograms(true).unwrap();
    assert_eq!(exp.chromatograms[0].peaks[0].rt, 0.1);
    assert_eq!(exp.chromatograms[1].peaks[0].rt, 0.2);

    let mut c = chrom(3.0, &[(2.0, 20.0), (3.0, 30.0), (1.0, 10.0), (2.0, 21.0)]);
    c.float_data_arrays
        .push(DataArray::new("f", vec![20.0, 30.0, 10.0, 21.0]));
    c.integer_data_arrays
        .push(DataArray::new("i", vec![20, 30, 10, 21]));
    c.string_data_arrays.push(DataArray::new(
        "s",
        vec!["b".into(), "c".into(), "a".into(), "b2".into()],
    ));
    c.string_data_arrays.push(DataArray::new("empty", vec![]));
    c.product
        .cv_terms
        .metadata
        .insert("label".into(), "identity".into());
    c.native_id = "original".into();
    let mut second = chrom(1.0, &[]);
    second.name = "first tie".into();
    let mut third = chrom(1.0, &[]);
    third.name = "second tie".into();
    let mut exp = MSExperiment {
        chromatograms: vec![c, second, third],
        ..Default::default()
    };
    exp.sort_chromatograms(true).unwrap();
    assert_eq!(exp.chromatograms[0].name, "first tie");
    assert_eq!(exp.chromatograms[1].name, "second tie");
    let c = &exp.chromatograms[2];
    assert_eq!(c.string_data_arrays[0].data, ["a", "b", "b2", "c"]);
    assert_eq!(c.float_data_arrays[0].data, [10.0, 20.0, 21.0, 30.0]);
    assert_eq!(c.integer_data_arrays[0].data, [10, 20, 21, 30]);
    assert!(c.string_data_arrays[1].data.is_empty());
    assert_eq!(
        c.product.cv_terms.metadata["label"],
        MetaValue::from("identity")
    );
    assert_eq!(c.native_id, "original");
}

#[test]
fn sort_failures_leave_every_chromatogram_and_annotation_unchanged() {
    let mut exp = MSExperiment {
        chromatograms: vec![
            chrom(2.0, &[(2.0, 1.0), (1.0, 2.0)]),
            chrom(1.0, &[(1.0, 1.0)]),
        ],
        ..Default::default()
    };
    for limits in [
        SummaryLimits {
            max_work: 1,
            ..Default::default()
        },
        SummaryLimits {
            max_bytes: 1,
            ..Default::default()
        },
        SummaryLimits {
            max_work: 30,
            ..Default::default()
        },
    ] {
        let original = exp.clone();
        assert!(exp.sort_chromatograms_with_limits(true, limits).is_err());
        assert_eq!(exp, original);
    }
    exp.chromatograms[1]
        .integer_data_arrays
        .push(DataArray::new("bad", vec![1, 2]));
    let before = exp.clone();
    assert!(exp.sort_chromatograms(true).is_err());
    assert_eq!(exp, before);
    // With RT sorting disabled, the malformed unused arrays are moved as-is.
    exp.sort_chromatograms(false).unwrap();
    assert_eq!(exp.chromatograms[0].integer_data_arrays[0].data, [1, 2]);
}

#[test]
fn sorting_moves_owned_string_payloads_without_cloning() {
    let mut c = chrom(1.0, &[(2.0, f32::NAN), (1.0, f32::INFINITY)]);
    let text = "x".repeat(10000);
    let pointer = text.as_ptr();
    c.string_data_arrays
        .push(DataArray::new("opaque", vec![text, "small".into()]));
    let mut exp = MSExperiment {
        chromatograms: vec![c],
        ..Default::default()
    };
    // This budget fits permutations, not the large string payload.
    exp.sort_chromatograms_with_limits(
        true,
        SummaryLimits {
            max_bytes: 1000,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        exp.chromatograms[0].string_data_arrays[0].data[1].as_ptr(),
        pointer
    );
    assert!(exp.chromatograms[0].peaks[0].intensity.is_infinite()); // Unconsumed field.
}

#[test]
fn total_count_and_exact_level_zero_predicates_use_the_source_domains() {
    let mut exp = MSExperiment {
        spectra: vec![
            spectrum(0.0, 1, &[1.0, -2.0]),
            spectrum(0.0, 2, &[-0.0]),
            spectrum(0.0, 3, &[]),
        ],
        chromatograms: vec![chrom(1.0, &[(0.0, 0.0), (1.0, 1.0)])],
        ..Default::default()
    };
    assert_eq!(exp.len(), 3);
    assert_eq!(exp.total_peak_count().unwrap(), 5);
    assert!(exp.contains_scan_of_level(3).unwrap());
    assert!(!exp.contains_scan_of_level(0).unwrap());
    assert!(!exp.contains_scan_of_level(usize::MAX).unwrap());
    assert!(!exp.has_zero_intensities(1).unwrap());
    assert!(exp.has_zero_intensities(2).unwrap());
    assert!(!exp.has_zero_intensities(3).unwrap());
    assert!(!exp.has_zero_intensities(0).unwrap());
    exp.spectra[1].ms_level = 0;
    assert!(exp.has_zero_intensities(0).unwrap());
    exp.spectra[0].peaks[0].intensity = f32::NAN;
    assert!(exp.has_zero_intensities(1).is_err());
    assert!(exp.has_zero_intensities(0).unwrap());
    assert_eq!(exp.total_peak_count().unwrap(), 5); // Counts do not validate scalars.
    assert!(exp.total_peak_count_with_limits(tiny()).is_err());
    assert!(exp.contains_scan_of_level_with_limits(1, tiny()).is_err());
    assert!(exp.has_zero_intensities_with_limits(1, tiny()).is_err());
}

#[test]
fn clear_arrays_is_source_spectrum_only_and_atomic_on_late_work_failure() {
    let mut exp = MSExperiment {
        spectra: vec![spectrum(0.0, 1, &[1.0]), spectrum(0.0, 2, &[2.0])],
        chromatograms: vec![chrom(1.0, &[(0.0, 3.0)])],
        ..Default::default()
    };
    exp.spectra[0]
        .integer_data_arrays
        .push(DataArray::new("empty", vec![]));
    exp.spectra[1]
        .float_data_arrays
        .push(DataArray::new("mismatched", vec![1.0, 2.0]));
    exp.spectra[1]
        .string_data_arrays
        .push(DataArray::new("s", vec!["owned".into()]));
    exp.spectra[1]
        .metadata
        .insert("ordinary".into(), "kept".into());
    exp.chromatograms[0]
        .integer_data_arrays
        .push(DataArray::new("retained", vec![8]));
    let before = exp.clone();
    assert!(
        exp.clear_meta_data_arrays_with_limits(SummaryLimits {
            max_work: 6,
            ..Default::default()
        })
        .is_err()
    );
    assert_eq!(exp, before);
    assert!(exp.clear_meta_data_arrays().unwrap());
    assert!(!exp.clear_meta_data_arrays().unwrap());
    assert_eq!(exp.total_peak_count().unwrap(), 3);
    assert_eq!(
        exp.spectra[1].metadata["ordinary"].as_str().unwrap(),
        "kept"
    );
    assert_eq!(exp.chromatograms, before.chromatograms);
    assert!(exp.spectra.iter().all(|s| s.float_data_arrays.is_empty()
        && s.integer_data_arrays.is_empty()
        && s.string_data_arrays.is_empty()));
}

#[test]
fn all_empty_operations_accept_zero_limits() {
    let mut exp = MSExperiment::new();
    assert_eq!(exp.total_peak_count_with_limits(tiny()).unwrap(), 0);
    assert!(!exp.contains_scan_of_level_with_limits(1, tiny()).unwrap());
    assert!(!exp.has_zero_intensities_with_limits(1, tiny()).unwrap());
    assert!(!exp.clear_meta_data_arrays_with_limits(tiny()).unwrap());
    exp.sort_chromatograms_with_limits(true, tiny()).unwrap();
    assert_eq!(
        exp.chromatogram_ranges_with_limits(tiny()).unwrap(),
        ExperimentRanges::default()
    );
}

#[test]
fn combined_ranges_preserve_source_spectra_first_signed_zero_endpoints() {
    let exp = MSExperiment {
        spectra: vec![MSSpectrum {
            rt: 0.0,
            peaks: vec![Peak1D::new(0.0, 0.0)],
            ..Default::default()
        }],
        chromatograms: vec![chrom(-0.0, &[(-0.0, -0.0)])],
        ..Default::default()
    };
    for range in [
        exp.combined_ranges().unwrap().rt,
        exp.combined_ranges().unwrap().mz,
        exp.combined_ranges().unwrap().intensity,
    ] {
        let range = range.unwrap();
        assert_eq!(range.min.to_bits(), 0.0_f64.to_bits());
        assert_eq!(range.max.to_bits(), 0.0_f64.to_bits());
    }
}
