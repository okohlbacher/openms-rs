// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::kernel::{AggregationLimits, DataArray, MzAggregation, MzRtRegion};
use openms::{Error, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, Result};

fn scan(rt: f64, level: u32, peaks: &[(f64, f32)]) -> MSSpectrum {
    MSSpectrum {
        rt,
        ms_level: level,
        peaks: peaks
            .iter()
            .map(|&(mz, intensity)| Peak1D::new(mz, intensity))
            .collect(),
        ..Default::default()
    }
}

// MSExperiment_test.cpp aggregate, aggregateFromMatrix and extractXICs cases.
fn source_experiment() -> MSExperiment {
    MSExperiment {
        spectra: vec![
            scan(1.0, 1, &[(100.0, 1000.0), (200.0, 2000.0), (300.0, 3000.0)]),
            scan(2.0, 2, &[(150.0, 1500.0), (250.0, 2500.0)]),
            scan(3.0, 1, &[(100.0, 1100.0), (200.0, 2100.0), (300.0, 3100.0)]),
            scan(4.0, 1, &[(100.0, 1200.0), (200.0, 2200.0), (300.0, 3200.0)]),
        ],
        ..Default::default()
    }
}

fn region(mz_min: f64, mz_max: f64, rt_min: f64, rt_max: f64) -> MzRtRegion {
    MzRtRegion::new(mz_min, mz_max, rt_min, rt_max).unwrap()
}

#[test]
fn source_custom_aggregate_and_ms_level_literals() -> Result<()> {
    let input = source_experiment();
    let windows = [
        region(90.0, 110.0, 0.0, 3.5),
        region(190.0, 210.0, 0.0, 5.0),
    ];
    let first = |peaks: &[Peak1D]| Ok(peaks.first().map_or(0.0, |p| f64::from(p.intensity)));
    assert_eq!(
        input.aggregate_with(&windows, 1, first)?,
        [vec![1000.0, 1100.0], vec![2000.0, 2100.0, 2200.0]]
    );
    assert_eq!(
        input.aggregate_with(&[region(140.0, 160.0, 1.5, 2.5)], 2, first)?,
        [vec![1500.0]]
    );
    assert_eq!(
        input.aggregate_with(&[region(90.0, 310.0, 0.0, 5.0)], 1, |p| MzAggregation::Mean
            .reduce(p))?,
        [vec![2000.0, 2100.0, 2200.0]]
    );
    assert!(input.aggregate(&windows, 0)?.is_empty());
    assert!(input.aggregate(&windows, 3)?.is_empty());
    Ok(())
}

#[test]
fn source_matrix_all_reducers() -> Result<()> {
    let input = source_experiment();
    assert_eq!(
        input.aggregate_from_matrix(
            &[[90.0, 110.0, 0.0, 3.5], [190.0, 210.0, 0.0, 5.0]],
            1,
            MzAggregation::Sum
        )?,
        [vec![1000.0, 1100.0], vec![2000.0, 2100.0, 2200.0]]
    );
    for (mode, expected) in [
        (MzAggregation::Sum, [6000.0, 6300.0, 6600.0]),
        (MzAggregation::Min, [1000.0, 1100.0, 1200.0]),
        (MzAggregation::Max, [3000.0, 3100.0, 3200.0]),
        (MzAggregation::Mean, [2000.0, 2100.0, 2200.0]),
    ] {
        assert_eq!(
            input.aggregate_from_matrix(&[[100.0, 300.0, 1.0, 4.0]], 1, mode)?,
            [expected.to_vec()]
        );
        let xics = input.extract_xics_from_matrix(&[[100.0, 300.0, 1.0, 4.0]], 1, mode)?;
        assert_eq!(
            xics[0]
                .peaks
                .iter()
                .map(|p| f64::from(p.intensity))
                .collect::<Vec<_>>(),
            expected
        );
        assert_eq!(xics[0].product.mz, 200.0);
    }
    for name in ["sum", "max", "min", "mean"] {
        assert!(name.parse::<MzAggregation>().is_ok());
    }
    for name in ["SUM", "", "median", " sum"] {
        assert!(name.parse::<MzAggregation>().is_err());
    }
    Ok(())
}

#[test]
fn source_xic_full_rt_product_and_default_metadata() -> Result<()> {
    let input = source_experiment();
    let xics = input.extract_xics(
        &[
            region(90.0, 110.0, 0.0, 5.0),
            region(190.0, 210.0, 0.0, 5.0),
        ],
        1,
    )?;
    assert_eq!(xics.len(), 2);
    assert_eq!(
        xics[0]
            .peaks
            .iter()
            .map(|p| (p.rt, p.intensity))
            .collect::<Vec<_>>(),
        [(1.0, 1000.0), (3.0, 1100.0), (4.0, 1200.0)]
    );
    assert_eq!(xics[0].product.mz, 100.0);
    assert_eq!(xics[1].product.mz, 200.0);
    assert_eq!(
        xics[1]
            .peaks
            .iter()
            .map(|p| p.intensity)
            .collect::<Vec<_>>(),
        [2000.0, 2100.0, 2200.0]
    );
    let ms2 = input.extract_xics(&[region(140.0, 160.0, 1.5, 2.5)], 2)?;
    assert_eq!(
        (
            ms2[0].peaks[0].rt,
            ms2[0].peaks[0].intensity,
            ms2[0].product.mz
        ),
        (2.0, 1500.0, 150.0)
    );
    let mut metadata_only = xics[0].clone();
    metadata_only.peaks.clear();
    let mut expected = MSChromatogram::new();
    expected.product.mz = 100.0;
    assert_eq!(metadata_only, expected);
    Ok(())
}

#[test]
fn default_and_matrix_accumulation_have_distinct_source_precision() -> Result<()> {
    let input = MSExperiment {
        spectra: vec![scan(
            1.0,
            1,
            &[(1.0, 16_777_216.0), (2.0, 1.0), (3.0, -16_777_216.0)],
        )],
        ..Default::default()
    };
    let window = [region(1.0, 3.0, 1.0, 1.0)];
    assert_eq!(input.aggregate(&window, 1)?, [vec![0.0]]);
    assert_eq!(
        input.aggregate_from_matrix(&[[1.0, 3.0, 1.0, 1.0]], 1, MzAggregation::Sum)?,
        [vec![1.0]]
    );
    assert_eq!(input.extract_xics(&window, 1)?[0].peaks[0].intensity, 0.0);
    assert_eq!(
        input.extract_xics_from_matrix(&[[1.0, 3.0, 1.0, 1.0]], 1, MzAggregation::Sum)?[0].peaks[0]
            .intensity,
        1.0
    );
    let third = input.extract_xics_with(&window, 1, |_| Ok(1.0 / 3.0))?[0].peaks[0].intensity;
    assert_eq!(third.to_bits(), (1.0_f32 / 3.0).to_bits());
    Ok(())
}

#[test]
fn inclusive_edges_duplicates_empty_slices_and_signed_values() -> Result<()> {
    let rt = 16_777_217.125;
    let input = MSExperiment {
        spectra: vec![
            scan(
                rt,
                1,
                &[(-1.0, -9.0), (0.0, -4.0), (0.0, -2.0), (1.0, -8.0)],
            ),
            scan(rt, 1, &[]),
        ],
        ..Default::default()
    };
    let rows = [
        [0.0, 0.0, rt, rt],
        [10.0, 11.0, rt, rt],
        [-1.0, 1.0, rt + 1.0, rt + 2.0],
    ];
    assert_eq!(
        input.aggregate_from_matrix(&rows, 1, MzAggregation::Max)?,
        [vec![-2.0, 0.0], vec![0.0, 0.0], vec![]]
    );
    assert_eq!(
        input.aggregate_from_matrix(&rows, 1, MzAggregation::Min)?,
        [vec![-4.0, 0.0], vec![0.0, 0.0], vec![]]
    );
    let xics = input.extract_xics_from_matrix(&rows, 1, MzAggregation::Sum)?;
    assert_eq!(xics[0].peaks.len(), 2);
    assert_eq!(xics[0].peaks[0].rt.to_bits(), rt.to_bits());
    assert_eq!(xics[0].peaks[1].rt.to_bits(), rt.to_bits());
    assert!(xics[2].is_empty());
    assert_eq!(xics[2].product.mz, 0.0);
    Ok(())
}

#[test]
fn custom_empty_slice_calls_and_deterministic_native_order() -> Result<()> {
    let input = source_experiment();
    let windows = [
        region(999.0, 999.0, 1.0, 3.0),
        region(100.0, 100.0, 1.0, 3.0),
    ];
    let mut seen = Vec::new();
    let output = input.aggregate_with(&windows, 1, |peaks| {
        seen.push(peaks.first().map(|p| p.intensity));
        Ok(seen.len() as f64)
    })?;
    assert_eq!(seen, [None, None, Some(1000.0), Some(1100.0)]);
    assert_eq!(output, [vec![1.0, 2.0], vec![3.0, 4.0]]);
    Ok(())
}

#[test]
fn source_empty_early_returns_and_matrix_constructor_validation() -> Result<()> {
    let input = MSExperiment {
        spectra: vec![scan(f64::NAN, 1, &[(f64::NAN, f32::NAN)])],
        ..Default::default()
    };
    let zero = AggregationLimits {
        max_work: 0,
        max_bytes: 0,
        max_output_points: 0,
    };
    assert!(
        input
            .aggregate_with_limits(&[], 1, |_| unreachable!(), &zero)?
            .is_empty()
    );
    assert!(
        input
            .extract_xics_with_limits(&[], 1, |_| unreachable!(), &zero)?
            .is_empty()
    );
    assert!(
        input
            .aggregate(&[region(0.0, 1.0, 0.0, 1.0)], 2)?
            .is_empty()
    );
    assert!(
        MSExperiment::new()
            .aggregate_from_matrix(&[[2.0, 1.0, 0.0, 1.0]], 2, MzAggregation::Sum)
            .is_err()
    );
    assert!(
        input
            .aggregate_from_matrix(&[], 1, MzAggregation::Sum)?
            .is_empty()
    );
    Ok(())
}

#[test]
fn validates_only_consumed_fields_and_relevant_sort_orders() -> Result<()> {
    let mut input = source_experiment();
    input.spectra[1].rt = f64::NAN; // Other MS level is not visited.
    input.spectra[1].peaks.reverse();
    input.spectra[0]
        .float_data_arrays
        .push(DataArray::new("unaligned", vec![0.0; 2]));
    input.spectra[0]
        .precursors
        .push(openms::Precursor::new(f64::NAN, 0));
    input.spectra[0].peaks[2].intensity = f32::NAN; // Outside selected m/z slice.
    let windows = [region(100.0, 100.0, 1.0, 1.0)];
    assert_eq!(input.aggregate(&windows, 1)?, [vec![1000.0]]);
    input.spectra[3].peaks.reverse(); // Outside selected RT interval.
    assert!(input.aggregate(&windows, 1).is_ok());
    input.spectra[0].peaks.reverse();
    assert!(matches!(
        input.aggregate(&windows, 1),
        Err(Error::UnsortedData)
    ));
    input.spectra[0].peaks.reverse();
    input.spectra.swap(0, 2);
    assert!(matches!(
        input.aggregate(&windows, 1),
        Err(Error::UnsortedData)
    ));
    Ok(())
}

#[test]
fn nonfinite_consumed_input_and_result_errors_leave_input_unchanged() -> Result<()> {
    let windows = [region(100.0, 300.0, 1.0, 4.0)];
    let input = source_experiment();
    let before = input.clone();
    let mut calls = 0;
    assert!(
        input
            .aggregate_with(&windows, 1, |_| {
                calls += 1;
                if calls == 3 {
                    Err(Error::InvalidValue("late callback failure".into()))
                } else {
                    Ok(1.0)
                }
            })
            .is_err()
    );
    assert_eq!(calls, 3);
    assert_eq!(input, before);
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(input.aggregate_with(&windows, 1, |_| Ok(value)).is_err());
    }
    assert!(
        input
            .extract_xics_with(&windows, 1, |_| Ok(f64::MAX))
            .is_err()
    );
    let mut malformed = source_experiment();
    malformed.spectra[0].peaks[0].intensity = f32::INFINITY;
    assert!(malformed.aggregate(&windows, 1).is_err());
    malformed.spectra[0].peaks[0] = Peak1D::new(f64::NAN, 1.0);
    assert!(malformed.aggregate(&windows, 1).is_err());
    assert_eq!(input, before);
    Ok(())
}

#[test]
fn f32_sum_and_product_midpoint_overflow_are_checked() -> Result<()> {
    let input = MSExperiment {
        spectra: vec![scan(0.0, 1, &[(1.0, f32::MAX), (2.0, f32::MAX)])],
        ..Default::default()
    };
    let window = [region(1.0, 2.0, 0.0, 0.0)];
    assert!(input.aggregate(&window, 1).is_err());
    assert_eq!(
        input.aggregate_from_matrix(&[[1.0, 2.0, 0.0, 0.0]], 1, MzAggregation::Sum)?,
        [vec![f64::from(f32::MAX) * 2.0]]
    );
    let far = [region(f64::MAX, f64::MAX, 10.0, 20.0)];
    assert_eq!(input.aggregate(&far, 1)?, [Vec::<f64>::new()]);
    assert!(input.extract_xics(&far, 1).is_err());
    Ok(())
}

#[test]
fn shared_limits_reject_before_output_expansion_or_callback() -> Result<()> {
    let input = source_experiment();
    let before = input.clone();
    let windows = [region(0.0, 500.0, 0.0, 5.0); 2];
    let points = AggregationLimits {
        max_output_points: 5,
        ..Default::default()
    };
    assert!(
        input
            .aggregate_with_limits(
                &windows,
                1,
                |_| panic!("output count must be preflighted"),
                &points
            )
            .unwrap_err()
            .to_string()
            .contains("output point limit")
    );
    let bytes = AggregationLimits {
        max_bytes: 0,
        ..Default::default()
    };
    assert!(
        input
            .aggregate_with_limits(
                &windows,
                1,
                |_| panic!("allocation must be preflighted"),
                &bytes
            )
            .unwrap_err()
            .to_string()
            .contains("byte limit")
    );
    let work = AggregationLimits {
        max_work: 3,
        ..Default::default()
    };
    assert!(
        input
            .aggregate_with_limits(
                &windows,
                1,
                |_| panic!("scan pass must be preflighted"),
                &work
            )
            .unwrap_err()
            .to_string()
            .contains("work limit")
    );
    let good = AggregationLimits {
        max_output_points: 6,
        ..Default::default()
    };
    assert_eq!(
        input
            .aggregate_with_limits(&windows, 1, |p| MzAggregation::Sum.reduce(p), &good)?
            .len(),
        2
    );
    assert_eq!(input, before);
    Ok(())
}

#[test]
fn overlapping_regions_share_work_instead_of_resetting() -> Result<()> {
    let input = MSExperiment {
        spectra: vec![scan(0.0, 1, &[(1.0, 1.0); 100])],
        ..Default::default()
    };
    let one = region(1.0, 1.0, 0.0, 0.0);
    let limits = AggregationLimits {
        max_work: 400,
        ..Default::default()
    };
    let sum = |p: &[Peak1D]| Ok(p.len() as f64);
    assert_eq!(
        input.aggregate_with_limits(&[one], 1, sum, &limits)?,
        [vec![100.0]]
    );
    let mut calls = 0;
    let error = input
        .aggregate_with_limits(
            &[one, one],
            1,
            |p| {
                calls += 1;
                sum(p)
            },
            &limits,
        )
        .unwrap_err();
    assert!(error.to_string().contains("work limit"));
    assert_eq!(calls, 1);
    assert!(
        input
            .extract_xics_with_limits(&[one, one], 1, sum, &limits)
            .is_err()
    );
    Ok(())
}

#[test]
fn product_metadata_survives_container_operations_and_clear_policy() -> Result<()> {
    let mut chromatogram = source_experiment()
        .extract_xics(&[region(90.0, 110.0, 0.0, 5.0)], 1)?
        .remove(0);
    chromatogram.product.isolation_window_lower_offset = 1.5;
    chromatogram.product.isolation_window_upper_offset = 2.5;
    let product = chromatogram.product.clone();
    chromatogram.sort_by_intensity(true)?;
    chromatogram.select(&[0, 1])?;
    chromatogram.retain_peaks(|p| p.intensity > 1000.0)?;
    assert_eq!(chromatogram.product, product);
    chromatogram.clear(false);
    assert_eq!(chromatogram.product, product);
    chromatogram.clear(true);
    assert_eq!(chromatogram, MSChromatogram::default());
    chromatogram.product.mz = f64::NAN;
    assert!(chromatogram.validate().is_err());
    Ok(())
}
