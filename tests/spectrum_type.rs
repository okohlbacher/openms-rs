// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::{
    MSSpectrum, Peak1D,
    kernel::{DataArray, SpectrumType, SpectrumTypeQueryLimits},
    metadata::{DataProcessing, ProcessingAction},
    processing::peak_picking::estimate_spectrum_type,
};
use std::{io::Cursor, sync::Arc};
fn signal(x: &[f64], y: &[f32]) -> MSSpectrum {
    MSSpectrum::from_peaks(
        x.iter()
            .zip(y)
            .map(|(&mz, &i)| Peak1D::new(mz, i))
            .collect(),
    )
}
fn history(action: ProcessingAction) -> Arc<DataProcessing> {
    Arc::new(DataProcessing {
        actions: [action].into(),
        ..Default::default()
    })
}
#[test]
fn all_literal_msspectrum_type_query_assertions() {
    let mut s = MSSpectrum::default();
    for query in [false, true] {
        assert_eq!(s.get_type(query).unwrap(), SpectrumType::Unknown);
    }
    s.spectrum_type = SpectrumType::Profile;
    for query in [false, true] {
        assert_eq!(s.get_type(query).unwrap(), SpectrumType::Profile);
    }
    let h = history(ProcessingAction::PeakPicking);
    s.data_processing.push(Arc::clone(&h));
    for query in [false, true] {
        assert_eq!(s.get_type(query).unwrap(), SpectrumType::Profile);
    }
    s.spectrum_type = SpectrumType::Unknown;
    for query in [false, true] {
        assert_eq!(s.get_type(query).unwrap(), SpectrumType::Centroid);
    }
    assert!(Arc::ptr_eq(&s.data_processing[0], &h));
    s.data_processing.clear();
    s.peaks = (1..=4)
        .map(|i| Peak1D::new(f64::from(i) * 100.0, 1.0))
        .collect();
    for query in [false, true] {
        assert_eq!(s.get_type(query).unwrap(), SpectrumType::Unknown);
    }
    s.peaks
        .extend([Peak1D::new(500.0, 1.0), Peak1D::new(600.0, 1.0)]);
    assert_eq!(s.get_type(false).unwrap(), SpectrumType::Unknown);
    assert_eq!(s.get_type(true).unwrap(), SpectrumType::Centroid);
    assert_eq!(s.spectrum_type, SpectrumType::Unknown);
}
#[test]
fn upstream_original_estimator_fixtures_keep_their_literal_classes() {
    for (input, expected) in [
        (
            include_str!("data/spectrum_type/PeakTypeEstimator_raw.dta"),
            SpectrumType::Profile,
        ),
        (
            include_str!("data/spectrum_type/PeakTypeEstimator_rawTOF.dta"),
            SpectrumType::Profile,
        ),
        (
            include_str!("data/spectrum_type/PeakTypeEstimator_peak.dta"),
            SpectrumType::Centroid,
        ),
    ] {
        let mut s = openms::format::dta::read(Cursor::new(input)).unwrap();
        assert_eq!(s.get_type(true).unwrap(), expected);
        assert_eq!(estimate_spectrum_type(&s).unwrap(), expected);
        s.peaks.truncate(4);
        assert_eq!(s.get_type(true).unwrap(), SpectrumType::Unknown);
    }
}
#[test]
fn finite_signed_unsorted_duplicate_shapes_have_source_shoulder_semantics() {
    let y = [0.0, 2.0, 3.0, 10.0, 3.0, 2.0, 0.0];
    for x in [[10.0; 7], [0.4, 0.1, 0.2, 0.3, 0.2, 0.5, 0.0]] {
        let s = signal(&x, &y);
        assert_eq!(s.get_type(true).unwrap(), SpectrumType::Profile);
        assert!(estimate_spectrum_type(&s).is_err());
    }
    let s = signal(
        &[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        &[-1.0, 2.0, 3.0, 10.0, 3.0, 2.0, -1.0],
    );
    assert_eq!(s.get_type(true).unwrap(), SpectrumType::Profile);
    assert!(estimate_spectrum_type(&s).is_err());
    for y in [[-1.0; 7], [0.0; 7]] {
        assert_eq!(
            signal(&[1.0; 7], &y).get_type(true).unwrap(),
            SpectrumType::Centroid
        );
    }
    assert_eq!(
        signal(&[f64::MAX; 7], &[1.0; 7]).get_type(true).unwrap(),
        SpectrumType::Centroid
    );
}
#[test]
fn strict_ratio_and_mz_window_boundaries_are_not_widened() {
    let x = [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3];
    assert_eq!(
        signal(&x, &[0.0, 0.1, 1.0, 10.0, 1.0, 0.1, 0.0])
            .get_type(true)
            .unwrap(),
        SpectrumType::Centroid
    );
    let y = [0.0, 2.0, 3.0, 10.0, 3.0, 2.0, 0.0];
    assert_eq!(
        signal(&[-2.0, -1.0, -0.1, 0.0, 0.1, 1.0, 2.0], &y)
            .get_type(true)
            .unwrap(),
        SpectrumType::Centroid
    );
    assert_eq!(
        signal(&[-2.0, -0.999, -0.1, 0.0, 0.1, 0.999, 2.0], &y)
            .get_type(true)
            .unwrap(),
        SpectrumType::Profile
    );
}
#[test]
fn source_three_of_four_is_centroid_four_of_five_is_profile() {
    // Equal maxima are encountered in coordinate order. Isolated center C
    // explains10; each two-sided center P explains20. Far-away unit noise
    // adjusts total TIC so the strict half-TIC stop occurs after4 or5 centers.
    for (profiles, noise, expected) in [
        (3, 50, SpectrumType::Centroid),
        (4, 80, SpectrumType::Profile),
    ] {
        let mut s = signal(&[0.0, 0.1, 0.2], &[0.0, 10.0, 0.0]);
        for i in 0..profiles {
            for (j, y) in [0.0, 2.0, 3.0, 10.0, 3.0, 2.0, 0.0].into_iter().enumerate() {
                s.peaks
                    .push(Peak1D::new((i + 1) as f64 * 10.0 + j as f64 * 0.1, y));
            }
        }
        for i in 0..noise {
            s.peaks.push(Peak1D::new(1000.0 + i as f64 * 2.0, 1.0));
        }
        assert_eq!(s.get_type(true).unwrap(), expected);
        assert_eq!(estimate_spectrum_type(&s).unwrap(), expected);
    }
}
#[test]
fn explicit_and_history_answers_ignore_unrelated_invalid_graphs() {
    let mut s = signal(&[f64::NAN; 6], &[f32::NAN; 6]);
    s.integer_data_arrays
        .push(DataArray::new("unaligned", vec![1]));
    let history = Arc::new(DataProcessing {
        actions: [ProcessingAction::PeakPicking].into(),
        metadata: [("huge".into(), "x".repeat(2_000_000).into())].into(),
        ..Default::default()
    });
    s.data_processing.push(Arc::clone(&history));
    s.spectrum_type = SpectrumType::Profile;
    let zero = SpectrumTypeQueryLimits {
        max_points: 0,
        max_work: 0,
        max_bytes: 0,
    };
    assert_eq!(
        s.get_type_with_limits(true, zero).unwrap(),
        SpectrumType::Profile
    );
    s.spectrum_type = SpectrumType::Unknown;
    assert_eq!(
        s.get_type_with_limits(
            true,
            SpectrumTypeQueryLimits {
                max_work: 13,
                ..zero
            }
        )
        .unwrap(),
        SpectrumType::Centroid
    );
    assert!(Arc::ptr_eq(&history, &s.data_processing[0]));
    s.data_processing.clear();
    assert_eq!(
        s.get_type_with_limits(false, zero).unwrap(),
        SpectrumType::Unknown
    );
    assert!(s.get_type(true).is_err());
    s.peaks.truncate(4);
    assert_eq!(
        s.get_type_with_limits(true, zero).unwrap(),
        SpectrumType::Unknown
    );
}
#[test]
fn history_and_estimation_limits_charge_only_consumed_paths_before_mutation() {
    let mut s = signal(&[1.0, 2.0, 3.0, 4.0, 5.0], &[1.0; 5]);
    let before = s.clone();
    for limit in [
        SpectrumTypeQueryLimits {
            max_points: 4,
            ..Default::default()
        },
        SpectrumTypeQueryLimits {
            max_work: 159,
            ..Default::default()
        },
        SpectrumTypeQueryLimits {
            max_bytes: 79,
            ..Default::default()
        },
    ] {
        assert!(s.get_type_with_limits(true, limit).is_err());
        assert_eq!(s, before);
    }
    assert_eq!(
        s.get_type_with_limits(
            true,
            SpectrumTypeQueryLimits {
                max_points: 5,
                max_work: 160,
                max_bytes: 80
            }
        )
        .unwrap(),
        SpectrumType::Centroid
    );
    let empty = Arc::new(DataProcessing::default());
    s.data_processing = vec![empty; 4];
    s.data_processing
        .push(history(ProcessingAction::PeakPicking));
    assert!(
        s.get_type_with_limits(
            false,
            SpectrumTypeQueryLimits {
                max_work: 16,
                ..Default::default()
            }
        )
        .is_err()
    );
    assert_eq!(
        s.get_type_with_limits(
            false,
            SpectrumTypeQueryLimits {
                max_work: 17,
                max_points: 0,
                max_bytes: 0
            }
        )
        .unwrap(),
        SpectrumType::Centroid
    );
}
