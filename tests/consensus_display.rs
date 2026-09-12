// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The two members of `KERNEL/ConsensusFeature.h` the feature-identification
//! package could not reach: `operator<<` and `Ratio::description_`.
//!
//! `ConsensusFeature_test.cpp` has no `START_SECTION` for either — all 39 of its
//! sections are ported in `tests/feature_identification.rs` — so the layout
//! expectations here are transcribed from the implementation,
//! `ConsensusFeature.cpp:391-418` (tier 3, source review, read off the
//! implementation rather than off a class-test literal). The description
//! ceilings and their atomicity are native invariants (tier 4): the source
//! member is an unchecked public vector.

use openms::kernel::Peak2D;
use openms::kernel::features::{
    ColumnHeader, ConsensusFeature, ConsensusMap, FeatureHandle, Ratio,
};
use openms::metadata::MetaValue;

/// The `cm_display` fixture of `tests/map_operations.rs`: one consensus feature
/// with one grouped handle and one meta value.
fn feature_with_one_handle() -> ConsensusFeature {
    let mut feature = ConsensusFeature::new();
    feature.base.rt = 1.5;
    feature.base.mz = 100.0;
    feature.base.intensity = 7.0;
    feature.base.quality = 0.25;
    feature
        .insert(FeatureHandle::from_peak(0, Peak2D::new(1.0, 99.0, 6.0), 3))
        .unwrap();
    feature.base.metadata.insert("note".into(), "x".into());
    feature
}

/// `ConsensusFeature.cpp:391-418`: banner, `Position: <rt> <mz>`, `Intensity`
/// and `Quality` without a colon, one five-line block per grouped handle,
/// `Meta information: ` and the closing banner whose line ends in a space.
#[test]
fn display_reproduces_the_source_stream_layout() {
    assert_eq!(
        feature_with_one_handle().to_string(),
        "---------- CONSENSUS ELEMENT BEGIN -----------------\n\
         Position: 1.5 100\n\
         Intensity 7\n\
         Quality 0.25\n\
         Grouped features: \n \
         - Map index: 0\n   \
         Feature id: 3\n   \
         RT: 1\n   \
         m/z: 99\n   \
         Intensity: 6\n\
         Meta information: \n   \
         note: x\n\
         ---------- CONSENSUS ELEMENT END ----------------- \n"
    );
}

/// An empty consensus feature still prints both banners, both label lines and
/// the zero-valued summary: the source loops write nothing but the surrounding
/// stream inserts are unconditional.
#[test]
fn display_of_an_empty_feature_is_banners_and_summary_only() {
    assert_eq!(
        ConsensusFeature::new().to_string(),
        "---------- CONSENSUS ELEMENT BEGIN -----------------\n\
         Position: 0 0\n\
         Intensity 0\n\
         Quality 0\n\
         Grouped features: \n\
         Meta information: \n\
         ---------- CONSENSUS ELEMENT END ----------------- \n"
    );
}

/// Handles print in handle order, which is the source `HandleSetType`'s
/// `IndexLess` order — map index, then unique ID — not insertion order. Meta
/// values print in key order.
#[test]
fn display_follows_handle_and_key_order() {
    let mut feature = ConsensusFeature::new();
    feature
        .insert(FeatureHandle::from_peak(7, Peak2D::new(2.0, 50.0, 1.0), 1))
        .unwrap();
    feature
        .insert(FeatureHandle::from_peak(1, Peak2D::new(3.0, 60.0, 2.0), 9))
        .unwrap();
    feature
        .insert(FeatureHandle::from_peak(1, Peak2D::new(4.0, 70.0, 3.0), 2))
        .unwrap();
    feature
        .base
        .metadata
        .insert("zeta".into(), MetaValue::from(1));
    feature
        .base
        .metadata
        .insert("alpha".into(), MetaValue::from(2));

    let printed = feature.to_string();
    let map_indices: Vec<&str> = printed
        .lines()
        .filter_map(|line| line.strip_prefix(" - Map index: "))
        .collect();
    assert_eq!(map_indices, vec!["1", "1", "7"]);
    let ids: Vec<&str> = printed
        .lines()
        .filter_map(|line| line.strip_prefix("   Feature id: "))
        .collect();
    assert_eq!(ids, vec!["2", "9", "1"]);
    let alpha = printed.find("   alpha: ").unwrap();
    let zeta = printed.find("   zeta: ").unwrap();
    assert!(alpha < zeta, "meta values must print in key order");
}

/// `ConsensusMap`'s own `Display` writes the same per-feature block through a
/// private helper in `src/kernel/map_operations.rs`, because that module could
/// not call an implementation that did not exist. The two must agree character
/// for character; if this fails, one of them has drifted.
#[test]
fn display_agrees_with_the_consensus_map_block() {
    let feature = feature_with_one_handle();
    let mut map = ConsensusMap::new();
    map.column_headers.insert(
        0,
        ColumnHeader {
            filename: "a.mzML".into(),
            label: "light".into(),
            size: 2,
            ..ColumnHeader::default()
        },
    );
    map.features.push(feature.clone());
    assert_eq!(
        map.to_string(),
        format!("Map 0: a.mzML - light - 2\n{feature}\n")
    );
}

/// The source `Ratio` default-constructs with an empty `description_` and
/// nothing in the SDK writes it.
#[test]
fn ratio_description_is_empty_by_default() {
    let ratio = Ratio::default();
    assert!(ratio.description().is_empty());
    assert!(ratio.description.is_empty());
    ratio.validate_description().unwrap();
}

/// Appending keeps insertion order and does not deduplicate: the source member
/// is a vector, so a repeated line is a repeated line.
#[test]
fn ratio_description_appends_in_order_and_keeps_duplicates() {
    let mut ratio = Ratio::default();
    ratio.add_description("numerator over denominator").unwrap();
    ratio.add_description("").unwrap();
    ratio.add_description("numerator over denominator").unwrap();
    assert_eq!(
        ratio.description(),
        [
            "numerator over denominator",
            "",
            "numerator over denominator"
        ]
    );
}

/// `set_description` replaces the whole list, and the public field stays
/// available for the unchecked assignment the source's public member allows.
#[test]
fn ratio_description_replaces_and_stays_publicly_assignable() {
    let mut ratio = Ratio::default();
    ratio.add_description("first").unwrap();
    ratio
        .set_description(vec!["second".into(), "third".into()])
        .unwrap();
    assert_eq!(ratio.description(), ["second", "third"]);
    ratio.description = vec!["straight to the field".into()];
    assert_eq!(ratio.description(), ["straight to the field"]);
}

/// Both ceilings reject before anything is stored, so a refused call leaves the
/// description byte-identical.
#[test]
fn ratio_description_ceilings_are_checked_and_atomic() {
    let mut ratio = Ratio::default();
    ratio.add_description("keep me").unwrap();
    let before = ratio.description.clone();

    let too_many = vec![String::new(); Ratio::MAX_DESCRIPTION_LINES + 1];
    assert!(ratio.set_description(too_many).is_err());
    assert_eq!(ratio.description, before);

    let too_long = vec!["x".repeat(Ratio::MAX_DESCRIPTION_BYTES + 1)];
    assert!(ratio.set_description(too_long).is_err());
    assert_eq!(ratio.description, before);

    // Exactly at each ceiling is accepted.
    let at_line_ceiling = vec![String::new(); Ratio::MAX_DESCRIPTION_LINES];
    ratio.set_description(at_line_ceiling).unwrap();
    assert_eq!(ratio.description().len(), Ratio::MAX_DESCRIPTION_LINES);
    // ... and one more line is then refused.
    assert!(ratio.add_description("").is_err());
    assert_eq!(ratio.description().len(), Ratio::MAX_DESCRIPTION_LINES);

    let at_byte_ceiling = vec!["y".repeat(Ratio::MAX_DESCRIPTION_BYTES)];
    ratio.set_description(at_byte_ceiling).unwrap();
    ratio.validate_description().unwrap();
    assert!(ratio.add_description("z").is_err());
    assert_eq!(ratio.description().len(), 1);
    assert_eq!(ratio.description()[0].len(), Ratio::MAX_DESCRIPTION_BYTES);
}

/// `Ratio::validate`, which `add_ratio` and `set_ratios` call, checks the ratio
/// value only. A description assigned straight to the field is therefore not
/// measured by the consensus feature, exactly as in the source, where nothing
/// validates a `Ratio` at all; `validate_description` is the explicit check.
#[test]
fn ratio_validate_ignores_the_description() {
    let mut ratio = Ratio {
        ratio_value: 2.0,
        numerator_ref: "heavy".into(),
        denominator_ref: "light".into(),
        description: vec![String::new(); Ratio::MAX_DESCRIPTION_LINES + 1],
    };
    ratio.validate().unwrap();
    assert!(ratio.validate_description().is_err());

    let mut feature = ConsensusFeature::new();
    feature.add_ratio(ratio.clone()).unwrap();
    assert_eq!(feature.ratios().len(), 1);
    assert!(feature.ratios()[0].validate_description().is_err());

    ratio.ratio_value = f64::NAN;
    assert!(ratio.validate().is_err());
}
