use openms::Error;
use openms::kernel::{DataArray, MobilityPeak1D, Mobilogram, MobilogramLimits};
use openms::metadata::{DataProcessing, DriftTimeUnit, MetaValue, MetaValueData, Unit};
use std::{
    collections::hash_map::DefaultHasher,
    hash::{Hash, Hasher},
    sync::Arc,
};

fn peak(m: f64, i: f32) -> MobilityPeak1D {
    MobilityPeak1D::new(m, i)
}
fn source() -> Mobilogram {
    // Literal source Mobilogram_test.cpp:945-967, with no C++ execution.
    Mobilogram::from_peaks(vec![
        peak(412.321, 29.),
        peak(412.824, 60.),
        peak(413.8, 34.),
        peak(414.301, 29.),
        peak(415.287, 37.),
        peak(416.293, 31.),
        peak(418.232, 31.),
        peak(419.113, 31.),
        peak(420.13, 201.),
        peak(423.269, 56.),
        peak(426.292, 34.),
        peak(427.28, 82.),
        peak(428.322, 87.),
        peak(430.269, 30.),
        peak(431.246, 29.),
        peak(432.289, 42.),
        peak(436.161, 32.),
        peak(437.219, 54.),
        peak(439.186, 40.),
        peak(440.27, 40.),
        peak(441.224, 23.),
    ])
}
fn annotated() -> Mobilogram {
    let mut m = Mobilogram::from_peaks(vec![
        peak(1., 10.),
        peak(2., 20.),
        peak(1., 20.),
        peak(2., 10.),
    ]);
    m.float_data_arrays
        .push(DataArray::new("f", vec![1., 2., 3., 4.]));
    m.integer_data_arrays
        .push(DataArray::new("i", vec![1, 2, 3, 4]));
    m.string_data_arrays.push(DataArray::new(
        "s",
        ["A", "B", "C", "D"].map(String::from).to_vec(),
    ));
    m.float_data_arrays.push(DataArray::new("empty", vec![]));
    m
}
fn hash(value: &MobilityPeak1D) -> u64 {
    let mut h = DefaultHasher::new();
    value.hash(&mut h);
    h.finish()
}

#[test]
fn peak_values_precision_display_hash_and_native_vector_operations() {
    assert_eq!(MobilityPeak1D::DIMENSION, 1);
    let mut p = MobilityPeak1D::default();
    assert_eq!(p, peak(0., 0.));
    p.mobility = 47.11;
    p.intensity = 29.;
    assert_eq!(format!("{p}"), "POS: 47.11 INT: 29");
    assert_eq!(format!("{p:.2}"), "POS: 47.11 INT: 29.00");
    assert_eq!(hash(&peak(-0., 0.)), hash(&peak(0., -0.)));
    assert_ne!(hash(&p), hash(&peak(47.12, 29.)));
    assert_ne!(hash(&p), hash(&peak(47.11, 30.)));
    let nan = peak(f64::NAN, 0.);
    assert_ne!(nan, nan);
    let mut m = Mobilogram::new();
    m.peaks.reserve(3);
    m.peaks.push(p);
    m.peaks.insert(0, peak(2., 1.));
    m.peaks.resize(3, MobilityPeak1D::default());
    assert_eq!(m.peaks.pop(), Some(peak(0., 0.)));
    m.peaks[1].mobility = 48.;
    assert_eq!(m.peaks.iter().next_back().unwrap().mobility, 48.);
    assert_eq!(m.peaks.remove(0), peak(2., 1.));
}

#[test]
fn source_defaults_units_ranges_and_stream_layout() {
    let mut m = Mobilogram::default();
    assert_eq!(m.rt, -1.);
    assert_eq!(m.drift_time_unit_as_str(), "<NONE>");
    assert_eq!(m.ranges().unwrap().mobility, None);
    m.rt = 0.451;
    m.drift_time_unit = DriftTimeUnit::Millisecond;
    assert_eq!(m.drift_time_unit_as_str(), "ms");
    m.peaks = vec![peak(2., 1.), peak(10., 2.), peak(2., 1.)];
    let r = m.ranges().unwrap();
    assert_eq!(
        (r.mobility.unwrap().min, r.mobility.unwrap().max),
        (2., 10.)
    );
    assert_eq!(
        (r.intensity.unwrap().min, r.intensity.unwrap().max),
        (1., 2.)
    );
    m.select(&[1]).unwrap();
    assert_eq!(m.ranges().unwrap().mobility.unwrap().min, 10.);
    assert_eq!(
        format!("{m}"),
        "-- MOBILOGRAM BEGIN --\nPOS: 10 INT: 2\n-- MOBILOGRAM END --\n"
    );
    assert_eq!(
        format!("{m:.1}"),
        "-- MOBILOGRAM BEGIN --\nPOS: 10.0 INT: 2.0\n-- MOBILOGRAM END --\n"
    );
}

#[test]
fn source_sort_ties_keep_every_parallel_array_aligned() {
    for (kind, expected) in [
        (0, vec![1, 3, 2, 4]),
        (1, vec![1, 4, 2, 3]),
        (2, vec![2, 3, 1, 4]),
    ] {
        let mut m = annotated();
        match kind {
            0 => m.sort_by_position().unwrap(),
            1 => m.sort_by_intensity(false).unwrap(),
            _ => m.sort_by_intensity(true).unwrap(),
        }
        assert_eq!(m.integer_data_arrays[0].data, expected);
        for (i, &identity) in expected.iter().enumerate() {
            assert_eq!(m.float_data_arrays[0].data[i], identity as f32);
            assert_eq!(
                m.string_data_arrays[0].data[i],
                char::from(b'A' + identity as u8 - 1).to_string()
            );
        }
        assert_eq!(m.float_data_arrays[0].name, "f");
        assert!(m.float_data_arrays[1].data.is_empty());
    }
}

#[test]
fn custom_index_sort_reads_original_arrays_and_checks_shapes_before_callback() {
    let mut m = annotated();
    m.sort_by(|m, a, b| m.string_data_arrays[0].data[a] > m.string_data_arrays[0].data[b])
        .unwrap();
    assert_eq!(m.integer_data_arrays[0].data, [4, 3, 2, 1]);
    assert!(
        m.is_sorted_by(
            |m, a, b| m.integer_data_arrays[0].data[a] > m.integer_data_arrays[0].data[b]
        )
        .unwrap()
    );
    m.float_data_arrays[0].data.pop();
    let before = m.clone();
    let mut calls = 0;
    assert!(
        m.sort_by(|_, _, _| {
            calls += 1;
            false
        })
        .is_err()
    );
    assert_eq!(calls, 0);
    assert_eq!(m, before);
}

#[test]
fn already_sorted_source_fast_path_does_not_use_malformed_arrays() {
    let mut m = Mobilogram::from_peaks(vec![peak(1., 1.), peak(2., 2.)]);
    m.integer_data_arrays
        .push(DataArray::new("unused", vec![9]));
    assert!(m.validate().is_err());
    m.sort_by_position().unwrap();
    m.sort_by_intensity(false).unwrap();
    assert!(m.sort_by_intensity(true).is_err());
    assert!(m.select(&[]).is_err());
}

#[test]
fn arbitrary_selection_is_atomic_and_moves_string_buffers_and_descriptions() {
    let mut m = annotated();
    let text = "x".repeat(100_000);
    m.string_data_arrays[0].data[3] = text;
    let pointer = m.string_data_arrays[0].data[3].as_ptr();
    m.string_data_arrays[0].metadata.insert(
        "unitful".into(),
        MetaValue::new(MetaValueData::Float(7.))
            .unwrap()
            .with_unit(Unit::new("UO:0000010", "second", "UO").unwrap())
            .unwrap(),
    );
    let processing = Arc::new(DataProcessing::default());
    m.string_data_arrays[0]
        .data_processing
        .push(processing.clone());
    let description = m.string_data_arrays[0].metadata.clone();
    m.select(&[3, 0, 2]).unwrap();
    assert_eq!(m.integer_data_arrays[0].data, [4, 1, 3]);
    assert_eq!(m.string_data_arrays[0].data[0].as_ptr(), pointer);
    assert_eq!(m.string_data_arrays[0].metadata, description);
    assert!(Arc::ptr_eq(
        &m.string_data_arrays[0].data_processing[0],
        &processing
    ));
    for indices in [&[0, 0][..], &[3][..], &[0, 1, 2, 0][..]] {
        let before = m.clone();
        assert!(m.select(indices).is_err());
        assert_eq!(m, before);
    }
    m.select(&[]).unwrap();
    assert!(m.is_empty());
    assert_eq!(m.string_data_arrays[0].metadata, description);
    assert!(m.string_data_arrays[0].data.is_empty());
}

#[test]
fn source_nearest_and_highest_literal_assertions() {
    let m = source();
    for (q, i) in [
        (400., 0),
        (500., 20),
        (412.4, 0),
        (441.224, 20),
        (426.29, 10),
        (426.3, 10),
        (427.2, 11),
        (427.3, 11),
    ] {
        assert_eq!(m.find_nearest(q).unwrap(), Some(i));
    }
    for (q, t, i) in [
        (400., 1., None),
        (500., 1., None),
        (412.4, 0.01, None),
        (412.4, 0.1, Some(0)),
        (441.3, 0.01, None),
        (441.3, 0.1, Some(20)),
        (427.3, 0.001, None),
    ] {
        assert_eq!(m.find_nearest_with_tolerance(q, t).unwrap(), i);
    }
    for (q, l, r, i) in [
        (427.3, 0.1, 0.001, Some(11)),
        (427.3, 0.001, 1.01, None),
        (427.3, 0.001, 1.1, Some(12)),
    ] {
        assert_eq!(m.find_nearest_in_window(q, l, r).unwrap(), i);
        assert_eq!(m.find_highest_in_window(q, l, r).unwrap(), i);
    }
    assert_eq!(m.find_highest_in_window(427.3, 9., 4.).unwrap(), Some(8));
    assert_eq!(
        m.find_highest_in_window(430.25, 1.9, 1.01).unwrap(),
        Some(13)
    );
    assert_eq!(m.base_peak_index().unwrap(), Some(8));
    assert_eq!(m.base_peak().unwrap().unwrap().intensity, 201.);
    assert_eq!(m.calculate_tic().unwrap(), 1032.);
}

#[test]
fn boundaries_duplicates_empty_and_one_sided_negative_tolerances() {
    let m = Mobilogram::from_peaks(vec![peak(1., 5.), peak(1., 7.), peak(3., 7.), peak(3., 9.)]);
    assert_eq!(m.find_nearest(1.).unwrap(), Some(0));
    assert_eq!(m.find_nearest(2.).unwrap(), Some(1));
    assert_eq!(m.find_nearest(4.).unwrap(), Some(3));
    assert_eq!(m.find_nearest_with_tolerance(2., 1.).unwrap(), Some(1));
    assert_eq!(m.find_nearest_in_window(2., 0.5, 1.).unwrap(), Some(2));
    assert_eq!(m.find_nearest_with_tolerance(2., -1.).unwrap(), None);
    let negative = Mobilogram::from_peaks(vec![peak(1., 1.), peak(2.5, 2.)]);
    assert_eq!(
        negative.find_nearest_in_window(2., -1., 1.).unwrap(),
        Some(1)
    ); // source checks only right bound
    assert!(negative.find_highest_in_window(2., -2., 1.).is_err());
    let empty = Mobilogram::new();
    assert_eq!(empty.find_nearest(1.).unwrap(), None);
    assert_eq!(empty.find_nearest_in_window(1., -2., -3.).unwrap(), None);
    assert_eq!(empty.find_highest_in_window(1., -2., -3.).unwrap(), None);
    assert_eq!(empty.base_peak().unwrap(), None);
    assert_eq!(empty.calculate_tic().unwrap(), 0.);
}

#[test]
fn lower_upper_bounds_and_subranges() {
    let m = Mobilogram::from_peaks(vec![
        peak(3., 0.),
        peak(1., 0.),
        peak(1., 0.),
        peak(2., 0.),
        peak(0., 0.),
    ]);
    assert!(matches!(m.mobility_begin(1.), Err(Error::UnsortedData)));
    assert_eq!(m.mobility_begin_in(1., 1..4).unwrap(), 1);
    assert_eq!(m.mobility_end_in(1., 1..4).unwrap(), 3);
    assert_eq!(m.mobility_begin_in(0., 1..4).unwrap(), 1);
    assert_eq!(m.mobility_end_in(5., 1..4).unwrap(), 4);
    assert_eq!(m.mobility_begin_in(1., 2..2).unwrap(), 2);
    assert!(m.mobility_end_in(1., 4..6).is_err());
    assert!(
        m.mobility_begin_in(1., std::ops::Range { start: 3, end: 2 })
            .is_err()
    );
}

#[test]
fn signed_intensities_first_maximum_and_f32_accumulation() {
    let mut m = Mobilogram::from_peaks(vec![peak(-2., -4.), peak(3., -2.), peak(1., -2.)]);
    assert_eq!(m.base_peak_index().unwrap(), Some(1));
    m.base_peak_mut().unwrap().unwrap().intensity = -1.;
    assert_eq!(m.peaks[1].intensity, -1.);
    assert_eq!(m.calculate_tic().unwrap(), -7.);
    m.peaks = vec![peak(0., 16_777_216.), peak(0., 1.), peak(0., 1.)];
    assert_eq!(m.calculate_tic().unwrap(), 16_777_216.);
    m.peaks = vec![peak(0., f32::MAX), peak(0., f32::MAX)];
    assert!(m.calculate_tic().is_err());
    m.peaks = vec![peak(-0., -0.), peak(0., 0.)];
    assert_eq!(
        m.ranges().unwrap().mobility.unwrap().min.to_bits(),
        (-0.0f64).to_bits()
    );
}

#[test]
fn full_equality_source_equality_partial_swap_and_clear() {
    let mut a = annotated();
    a.rt = 5.;
    a.drift_time_unit = DriftTimeUnit::Millisecond;
    let mut b = a.clone();
    b.string_data_arrays[0]
        .metadata
        .insert("x".into(), "different".into());
    assert_ne!(a, b);
    assert!(a.source_equal(&b));
    b.peaks = vec![peak(50., 5.)];
    b.rt = 10.;
    b.drift_time_unit = DriftTimeUnit::InverseReducedMobility;
    a.swap_peak_data(&mut b);
    assert_eq!(a.len(), 1);
    assert_eq!(a.rt, 10.);
    assert_eq!(a.string_data_arrays[0].data.len(), 4);
    assert!(a.validate().is_err());
    a.clear();
    assert!(a.is_empty());
    assert!(a.string_data_arrays.is_empty());
    assert_eq!(a.rt, 10.);
    assert_eq!(a.drift_time_unit, DriftTimeUnit::InverseReducedMobility);
}

#[test]
fn resource_errors_and_nonfinite_consumed_values_are_atomic() {
    let m = annotated();
    for limits in [
        MobilogramLimits {
            max_peaks: 3,
            ..Default::default()
        },
        MobilogramLimits {
            max_arrays: 0,
            ..Default::default()
        },
        MobilogramLimits {
            max_work: 1,
            ..Default::default()
        },
        MobilogramLimits {
            max_bytes: 0,
            ..Default::default()
        },
    ] {
        let mut candidate = m.clone();
        assert!(candidate.sort_by_position_with_limits(limits).is_err());
        assert_eq!(candidate, m);
        let mut candidate = m.clone();
        assert!(candidate.select_with_limits(&[3, 1], limits).is_err());
        assert_eq!(candidate, m);
    }
    let mut candidate = m.clone();
    let mut comparisons = 0;
    assert!(
        candidate
            .sort_by_with_limits(
                |m, a, b| {
                    comparisons += 1;
                    m.peaks[a].intensity < m.peaks[b].intensity
                },
                MobilogramLimits {
                    max_work: 23,
                    ..Default::default()
                }
            )
            .is_err()
    );
    assert!(comparisons > 0);
    assert_eq!(candidate, m);
    let mut bad = Mobilogram::from_peaks(vec![peak(f64::NAN, 1.)]);
    assert!(bad.find_nearest(1.).is_err());
    assert_eq!(bad.calculate_tic().unwrap(), 1.);
    bad.peaks = vec![peak(1., f32::NAN)];
    assert_eq!(bad.find_nearest(1.).unwrap(), Some(0));
    assert!(bad.base_peak().is_err());
    assert!(bad.find_nearest(f64::INFINITY).is_err());
}
