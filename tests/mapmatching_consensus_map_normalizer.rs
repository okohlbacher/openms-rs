// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `ConsensusMapNormalizerAlgorithmMedian_test.cpp`,
//! `ConsensusMapNormalizerAlgorithmQuantile_test.cpp` and
//! `ConsensusMapNormalizerAlgorithmThreshold_test.cpp` at source revision
//! bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`).
//!
//! The median sections are evidence tier 3: every expected value is a
//! transcribed class-test literal. The quantile and threshold class tests are
//! `TODO` and `NOT_TESTABLE` upstream, so their sections here are tier 4: the
//! expected values are derived by hand from the algorithms as the source
//! headers document them, with the arithmetic shown at each assertion. The
//! threshold algorithm also has a tier-1 differential against a retained C++
//! output in `tests/mapmatching_consensus_map_normalizer_topp.rs`. No expected
//! value was derived from this port's output.

use openms::analysis::mapmatching::consensus_map_normalizer_median::{
    ConsensusMapNormalizerAlgorithmMedian as Median, IdentificationFilter, NormalizationMethod,
};
use openms::analysis::mapmatching::consensus_map_normalizer_quantile::ConsensusMapNormalizerAlgorithmQuantile as Quantile;
use openms::analysis::mapmatching::consensus_map_normalizer_threshold::ConsensusMapNormalizerAlgorithmThreshold as Threshold;
use openms::chemistry::AASequence;
use openms::identification::{
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
};
use openms::kernel::features::{ColumnHeader, ConsensusFeature, ConsensusMap, FeatureHandle};

/// `TEST_REAL_SIMILAR` with the ClassTest default relative tolerance of 1e-5.
fn similar(actual: f64, expected: f64) {
    let tolerance = 1e-5 * expected.abs().max(actual.abs()).max(1e-5);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} is not similar to {expected}"
    );
}

/// Class-test `makeCF`: one consensus feature with one handle per
/// `(map index, intensity)` pair and a recognisable own intensity of 123456.
fn make_cf(handles: &[(u64, f32)], rt: f64, mz: f64, uid: &mut u64) -> ConsensusFeature {
    let mut cf = ConsensusFeature::new();
    cf.rt = rt;
    cf.mz = mz;
    cf.intensity = 123_456.0;
    for &(map_index, intensity) in handles {
        *uid += 1;
        cf.insert(FeatureHandle {
            map_index,
            unique_id: *uid,
            rt,
            mz,
            intensity,
            ..FeatureHandle::default()
        })
        .unwrap();
    }
    cf
}

/// Class-test `makeMap`: three maps, three features, per-map medians
/// 200 / 400 / 800.
fn make_map(header_sizes: &[usize]) -> ConsensusMap {
    let mut map = ConsensusMap::new();
    for (i, &size) in header_sizes.iter().enumerate() {
        map.column_headers.insert(
            i as u64,
            ColumnHeader {
                filename: format!("map{i}.mzML"),
                size,
                ..ColumnHeader::default()
            },
        );
    }
    let mut uid = 0;
    map.features.push(make_cf(
        &[(0, 100.0), (1, 200.0), (2, 400.0)],
        100.0,
        500.0,
        &mut uid,
    ));
    map.features.push(make_cf(
        &[(0, 300.0), (1, 600.0), (2, 1200.0)],
        200.0,
        600.0,
        &mut uid,
    ));
    map.features.push(make_cf(
        &[(0, 200.0), (1, 400.0), (2, 800.0)],
        300.0,
        700.0,
        &mut uid,
    ));
    map
}

/// Class-test `intensitiesOf`: all handle intensities of one map, in consensus
/// feature order.
fn intensities_of(map: &ConsensusMap, map_index: u64) -> Vec<f64> {
    map.features
        .iter()
        .flat_map(|cf| cf.handles().iter())
        .filter(|fh| fh.map_index == map_index)
        .map(|fh| f64::from(fh.intensity))
        .collect()
}

fn no_filter() -> IdentificationFilter {
    IdentificationFilter::new("", "").unwrap()
}

// ---------------------------------------------------------------------------
// ConsensusMapNormalizerAlgorithmMedian_test.cpp
// ---------------------------------------------------------------------------

#[test]
fn compute_medians_one_per_map() {
    let map = make_map(&[3, 9, 3]);
    let result = Median::compute_medians(&map, &no_filter()).unwrap();
    assert_eq!(result.medians.len(), 3);
    similar(result.medians[0], 200.0);
    similar(result.medians[1], 400.0);
    similar(result.medians[2], 800.0);
}

#[test]
fn compute_medians_reference_is_argmax_of_header_size() {
    let map = make_map(&[3, 9, 3]);
    assert_eq!(
        Median::compute_medians(&map, &no_filter())
            .unwrap()
            .reference_map,
        1
    );
}

#[test]
fn compute_medians_ties_keep_lowest_index() {
    let map = make_map(&[9, 9, 9]);
    assert_eq!(
        Median::compute_medians(&map, &no_filter())
            .unwrap()
            .reference_map,
        0
    );
}

/// OpenMS issue #9864, finding A3: with every stored size 0 the reference
/// collapses to map 0 whatever the data says. Pinned, as upstream pins it.
#[test]
fn compute_medians_degenerate_reference_collapses_to_zero() {
    let mut map = make_map(&[0, 0, 0]);
    let mut uid = 1000;
    map.features
        .push(make_cf(&[(2, 1600.0)], 400.0, 800.0, &mut uid));
    map.features
        .push(make_cf(&[(2, 3200.0)], 500.0, 900.0, &mut uid));
    assert_eq!(intensities_of(&map, 0).len(), 3);
    assert_eq!(intensities_of(&map, 2).len(), 5);
    assert_eq!(
        Median::compute_medians(&map, &no_filter())
            .unwrap()
            .reference_map,
        0
    );
}

#[test]
fn compute_medians_bails_out_when_a_map_has_no_features() {
    let mut map = make_map(&[3, 3, 3]);
    map.column_headers.insert(
        3,
        ColumnHeader {
            filename: "map3.mzML".into(),
            size: 0,
            ..ColumnHeader::default()
        },
    );
    let result = Median::compute_medians(&map, &no_filter()).unwrap();
    assert_eq!(result.reference_map, 0);
    assert_eq!(result.medians.len(), 4);
    for median in &result.medians {
        similar(*median, 1.0);
    }
    // Native addition: the bail-out is visible to the caller.
    assert!(!result.enough_features);
}

#[test]
fn compute_medians_missing_column_header_is_an_error() {
    let mut map = make_map(&[3, 3, 3]);
    map.column_headers.remove(&1);
    map.column_headers.insert(
        7,
        ColumnHeader {
            filename: "map7.mzML".into(),
            ..ColumnHeader::default()
        },
    );
    // Source: Exception::ElementNotFound.
    assert!(Median::compute_medians(&map, &no_filter()).is_err());
}

#[test]
fn normalize_maps_scale_to_reference_median() {
    let mut map = make_map(&[3, 9, 3]);
    Median::normalize_maps(&mut map, NormalizationMethod::Scale, &no_filter()).unwrap();
    for index in 0..3 {
        let m = intensities_of(&map, index);
        similar(m[0], 200.0);
        similar(m[1], 600.0);
        similar(m[2], 400.0);
    }
}

#[test]
fn normalize_maps_leaves_consensus_intensity_alone() {
    let mut map = make_map(&[3, 9, 3]);
    Median::normalize_maps(&mut map, NormalizationMethod::Scale, &no_filter()).unwrap();
    for cf in &map.features {
        similar(f64::from(cf.intensity), 123_456.0);
    }
}

#[test]
fn normalize_maps_scale_ratios_invariant_to_reference() {
    let mut map = make_map(&[9, 3, 3]);
    Median::normalize_maps(&mut map, NormalizationMethod::Scale, &no_filter()).unwrap();
    let m2 = intensities_of(&map, 2);
    similar(m2[0], 100.0);
    similar(m2[1], 300.0);
    similar(m2[2], 200.0);
    similar(m2[1] / m2[0], 3.0);
}

#[test]
fn normalize_maps_shift_towards_largest_median() {
    let mut map = make_map(&[3, 9, 3]);
    Median::normalize_maps(&mut map, NormalizationMethod::Shift, &no_filter()).unwrap();
    let m0 = intensities_of(&map, 0);
    let m2 = intensities_of(&map, 2);
    similar(m0[0], 700.0);
    similar(m0[1], 900.0);
    similar(m2[0], 400.0);
    similar(m2[1], 1200.0);
}

#[test]
fn normalize_maps_bail_out_is_a_no_op() {
    let mut map = make_map(&[3, 3, 3]);
    map.column_headers.insert(
        3,
        ColumnHeader {
            filename: "map3.mzML".into(),
            ..ColumnHeader::default()
        },
    );
    Median::normalize_maps(&mut map, NormalizationMethod::Scale, &no_filter()).unwrap();
    let m2 = intensities_of(&map, 2);
    similar(m2[0], 400.0);
    similar(m2[1], 1200.0);
    similar(m2[2], 800.0);
}

#[test]
fn passes_filters() {
    let mut map = make_map(&[3, 3, 3]);
    let passes = |map: &ConsensusMap, acc: &str, desc: &str| {
        IdentificationFilter::new(acc, desc)
            .unwrap()
            .passes(&map.features[0], map)
            .unwrap()
    };

    assert!(passes(&map, "", ""));
    assert!(!passes(&map, "ALBU", ""));

    let mut hit = PeptideHit {
        sequence: AASequence::parse("PEPTIDER").unwrap(),
        ..PeptideHit::default()
    };
    hit.evidences.push(PeptideEvidence {
        protein_accession: "P02769|ALBU_BOVIN".into(),
        ..PeptideEvidence::default()
    });
    let pid = PeptideIdentification {
        hits: vec![hit],
        ..PeptideIdentification::default()
    };
    map.features[0].peptide_identifications = vec![pid];

    assert!(passes(&map, "ALBU", ""));
    assert!(!passes(&map, "TRYP", ""));

    let mut protein_hit = ProteinHit {
        accession: "P02769|ALBU_BOVIN".into(),
        ..ProteinHit::default()
    };
    protein_hit.set_description(format!("Serum albumin {} bovine", "x".repeat(4096)));
    let protein_id = ProteinIdentification {
        hits: vec![protein_hit],
        ..ProteinIdentification::default()
    };
    map.protein_identifications = vec![protein_id];

    assert!(passes(&map, "ALBU", "^Serum albumin .* bovine$"));
    assert!(!passes(&map, "ALBU", "trypsin"));
    assert!(!passes(&map, "TRYP", "albumin"));
}

// Native checks of the median port (tier 4).

/// An expression that matches the empty string is treated as absent, as the
/// source's `regex_search("", ...)` test does: every feature passes, even an
/// unidentified one.
#[test]
fn filter_matching_empty_string_passes_everything() {
    let map = make_map(&[3, 3, 3]);
    let filter = IdentificationFilter::new("ALBU|", "x*").unwrap();
    assert!(filter.passes(&map.features[0], &map).unwrap());
}

#[test]
fn malformed_filter_is_an_error() {
    assert!(IdentificationFilter::new("(", "").is_err());
}

#[test]
fn handle_outside_column_headers_is_refused_without_mutation() {
    let mut map = make_map(&[3, 9, 3]);
    let mut uid = 50;
    map.features.push(make_cf(&[(5, 10.0)], 1.0, 1.0, &mut uid));
    let before = map.clone();
    assert!(Median::normalize_maps(&mut map, NormalizationMethod::Scale, &no_filter()).is_err());
    assert_eq!(map, before);
}

/// A zero median would divide by zero; the source stores the infinity, this
/// port refuses and leaves the map unchanged.
#[test]
fn zero_median_is_refused_without_mutation() {
    let mut map = make_map(&[3, 9, 3]);
    for cf in &mut map.features {
        let handles = cf
            .handles()
            .iter()
            .map(|h| FeatureHandle {
                intensity: if h.map_index == 0 { 0.0 } else { h.intensity },
                ..*h
            })
            .collect();
        cf.set_handles(handles).unwrap();
    }
    let before = map.clone();
    assert!(Median::normalize_maps(&mut map, NormalizationMethod::Scale, &no_filter()).is_err());
    assert_eq!(map, before);
}

/// The shift target search starts at the smallest positive double, so with
/// every median non-positive the target is map 0 (tier 4, from the source's
/// `numeric_limits<double>::min()` initialiser).
#[test]
fn shift_target_defaults_to_map_zero_for_non_positive_medians() {
    let mut map = ConsensusMap::new();
    for i in 0..2u64 {
        map.column_headers.insert(i, ColumnHeader::default());
    }
    let mut uid = 0;
    map.features
        .push(make_cf(&[(0, -10.0), (1, -2.0)], 1.0, 1.0, &mut uid));
    Median::normalize_maps(&mut map, NormalizationMethod::Shift, &no_filter()).unwrap();
    // Map 0 median -10 is the target although map 1's -2 is larger:
    // map 1 becomes -2 + (-10) - (-2) = -10.
    similar(intensities_of(&map, 0)[0], -10.0);
    similar(intensities_of(&map, 1)[0], -10.0);
}

// ---------------------------------------------------------------------------
// ConsensusMapNormalizerAlgorithmQuantile_test.cpp (all sections TODO
// upstream; tier 4)
// ---------------------------------------------------------------------------

/// Delta `(3 - 1) / (5 - 1) = 0.5`: positions 0, 0.5, 1, 1.5, 2.
#[test]
fn resample_interpolates_linearly() {
    let out = Quantile::resample(&[1.0, 2.0, 3.0], 5).unwrap();
    assert_eq!(out, vec![1.0, 1.5, 2.0, 2.5, 3.0]);
}

/// Down-sampling keeps the end points; delta `(5 - 1) / (3 - 1) = 2`.
#[test]
fn resample_down_keeps_end_points() {
    let out = Quantile::resample(&[1.0, 2.0, 3.0, 4.0, 10.0], 3).unwrap();
    assert_eq!(out, vec![1.0, 3.0, 10.0]);
}

/// Zero points clear; one point takes the back value, since front and back are
/// both written to slot 0 and back is written last.
#[test]
fn resample_zero_and_one_point() {
    assert!(Quantile::resample(&[1.0, 2.0], 0).unwrap().is_empty());
    assert_eq!(Quantile::resample(&[1.0, 2.0], 1).unwrap(), vec![2.0]);
    assert!(Quantile::resample(&[], 3).is_err());
}

#[test]
fn extract_intensity_vectors_in_feature_order() {
    let map = make_map(&[3, 3, 3]);
    let out = Quantile::extract_intensity_vectors(&map).unwrap();
    assert_eq!(out.len(), 3);
    assert_eq!(out[0], vec![100.0, 300.0, 200.0]);
    assert_eq!(out[1], vec![200.0, 600.0, 400.0]);
    assert_eq!(out[2], vec![400.0, 1200.0, 800.0]);
}

/// Sorted maps are [100,200,300], [200,400,600], [400,800,1200]; all have
/// three values, so resampling is the identity and the reference is the
/// per-rank mean [700/3, 1400/3, 700]. Every map's ranks are (low, high, mid)
/// in feature order, so every map becomes [700/3, 700, 1400/3].
#[test]
fn normalize_maps_quantile_equal_sizes() {
    let mut map = make_map(&[3, 3, 3]);
    Quantile::normalize_maps(&mut map).unwrap();
    for index in 0..3 {
        let m = intensities_of(&map, index);
        similar(m[0], 700.0 / 3.0);
        similar(m[1], 700.0);
        similar(m[2], 1400.0 / 3.0);
    }
}

/// Map 0 has two handles [10, 30], map 1 has three [5, 20, 15]. Largest is 3:
/// map 0 resamples to [10, 20, 30], map 1 sorts to [5, 15, 20]; the reference
/// is [7.5, 17.5, 25]. Back to two points map 0 gets [7.5, 25] in rank order;
/// map 1 keeps three points and gets 7.5 at its smallest (5), 17.5 at 15 and 25
/// at 20.
#[test]
fn normalize_maps_quantile_unequal_sizes() {
    let mut map = ConsensusMap::new();
    for i in 0..2u64 {
        map.column_headers.insert(i, ColumnHeader::default());
    }
    let mut uid = 0;
    map.features
        .push(make_cf(&[(0, 10.0), (1, 5.0)], 1.0, 1.0, &mut uid));
    map.features
        .push(make_cf(&[(0, 30.0), (1, 20.0)], 2.0, 2.0, &mut uid));
    map.features.push(make_cf(&[(1, 15.0)], 3.0, 3.0, &mut uid));
    Quantile::normalize_maps(&mut map).unwrap();
    let m0 = intensities_of(&map, 0);
    let m1 = intensities_of(&map, 1);
    similar(m0[0], 7.5);
    similar(m0[1], 25.0);
    similar(m1[0], 7.5);
    similar(m1[1], 25.0);
    similar(m1[2], 17.5);
}

/// A map with no handle while another has some: the source reads `front()` of
/// an empty vector. Refused here, without mutation.
#[test]
fn normalize_maps_quantile_empty_map_is_refused() {
    let mut map = make_map(&[3, 3, 3]);
    map.column_headers.insert(3, ColumnHeader::default());
    let before = map.clone();
    assert!(Quantile::normalize_maps(&mut map).is_err());
    assert_eq!(map, before);
}

#[test]
fn set_normalized_intensity_values_round_trip_and_shape_check() {
    let mut map = make_map(&[3, 3, 3]);
    let mut values = Quantile::extract_intensity_vectors(&map).unwrap();
    values[1][2] = 42.0;
    Quantile::set_normalized_intensity_values(&values, &mut map).unwrap();
    assert_eq!(intensities_of(&map, 1), vec![200.0, 600.0, 42.0]);
    values[0].pop();
    assert!(Quantile::set_normalized_intensity_values(&values, &mut map).is_err());
}

// ---------------------------------------------------------------------------
// ConsensusMapNormalizerAlgorithmThreshold_test.cpp (both sections
// NOT_TESTABLE upstream; tier 4 here, tier 1 in the TOPP differential)
// ---------------------------------------------------------------------------

/// Reference map 1 (size 9). Per feature the ratios reference / map are 2 for
/// map 0, 1 for map 1 and 0.5 for map 2. With threshold 0.4 the window
/// (0.4, 2.5) keeps them all.
#[test]
fn compute_correlation_mean_ratio_per_map() {
    let map = make_map(&[3, 9, 3]);
    let result = Threshold::compute_correlation(&map, 0.4, &no_filter()).unwrap();
    assert!(result.enough_features);
    assert_eq!(result.reference_map, 1);
    assert_eq!(result.passing_features, 3);
    similar(result.ratios[0], 2.0);
    similar(result.ratios[1], 1.0);
    similar(result.ratios[2], 0.5);
}

/// With threshold 0.5 the window is the open interval (0.5, 2): the ratios 2
/// of map 0 are excluded, so the source bails out with all 1.0.
#[test]
fn compute_correlation_open_window_bails_out() {
    let map = make_map(&[3, 9, 3]);
    let result = Threshold::compute_correlation(&map, 0.5, &no_filter()).unwrap();
    assert!(!result.enough_features);
    assert_eq!(result.ratios, vec![1.0, 1.0, 1.0]);
}

/// Applying the factors 2, 1, 0.5 brings every map to [200, 600, 400].
#[test]
fn threshold_normalize_maps_multiplies() {
    let mut map = make_map(&[3, 9, 3]);
    Threshold::normalize_maps(&mut map, &[2.0, 1.0, 0.5]).unwrap();
    for index in 0..3 {
        assert_eq!(intensities_of(&map, index), vec![200.0, 600.0, 400.0]);
    }
    let before = map.clone();
    assert!(Threshold::normalize_maps(&mut map, &[1.0, 1.0]).is_err());
    assert_eq!(map, before);
}

/// A zero-intensity handle counts as absent from the ratio.
#[test]
fn compute_correlation_ignores_zero_intensities() {
    let mut map = make_map(&[3, 9, 3]);
    let first = &mut map.features[0];
    let handles = first
        .handles()
        .iter()
        .map(|h| FeatureHandle {
            intensity: if h.map_index == 0 { 0.0 } else { h.intensity },
            ..*h
        })
        .collect();
    first.set_handles(handles).unwrap();
    let result = Threshold::compute_correlation(&map, 0.4, &no_filter()).unwrap();
    similar(result.ratios[0], 2.0);
}
