// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `BaseFeature_test.cpp`, `Feature_test.cpp` and
//! `ConsensusFeature_test.cpp` at Core SDK
//! `bc9cc12514c768385ce121d6ca4bb710fe1983c4`, plus native checks for the
//! identification surface added in `src/kernel/feature_identification.rs`.
//!
//! Every `START_SECTION` of the three class tests has one test function here;
//! `docs/FEATURE_IDENTIFICATION_SUPPORT.md` holds the section-to-test table.
//! All transcribed expectations are tier 3 (source review): the literals come
//! from the pinned class tests, no C++ was built or executed.

use openms::chemistry::{AASequence, PROTON_MASS_U};
use openms::concept::HasUniqueId;
use openms::identification::graph::{
    IdentificationData, IdentifiedCompound, IdentifiedPeptide, InputFile, Observation,
    ObservationMatch, ObservationMatchId, PeptideId, ProcessingSoftware, ProcessingStep,
};
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::kernel::feature_identification::AnnotationState;
use openms::kernel::features::{
    BaseFeature, ConsensusFeature, ConsensusMap, Feature, FeatureHandle, FeatureMap, Ratio,
};
use openms::kernel::geometry::{ConvexHull2D, Point2D};
use openms::kernel::{Peak2D, RichPeak2D};
use openms::metadata::MetaValue;

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a} != {b}");
}
fn close_f32(a: f32, b: f32) {
    assert!((a - b).abs() < 1e-4, "{a} != {b}");
}
/// The two-scan hulls the Feature class test reuses across its sections.
fn source_hulls() -> Vec<ConvexHull2D> {
    vec![
        ConvexHull2D::from_points(&[Point2D::new(1.0, 2.0), Point2D::new(3.0, 4.0)]).unwrap(),
        ConvexHull2D::from_points(&[Point2D::new(0.5, 0.0), Point2D::new(1.0, 1.0)]).unwrap(),
    ]
}
/// The same hulls after the class test adds one more point to each.
fn source_hulls_extended() -> Vec<ConvexHull2D> {
    let mut hulls = source_hulls();
    hulls[0].add_point(Point2D::new(3.0, 2.0)).unwrap();
    hulls[1].add_point(Point2D::new(2.0, 1.0)).unwrap();
    hulls
}
/// `tmp_feature`, `tmp_feature2` and `tmp_feature3` of `ConsensusFeature_test.cpp`.
fn consensus_fixture() -> [Feature; 3] {
    let make = |rt: f64, mz: f64, intensity: f32, id: u64| {
        let mut feature = Feature::new(rt, mz, intensity);
        feature.unique_id = id;
        feature
    };
    [
        make(1.0, 2.0, 200.0, 3),
        make(2.0, 3.0, 300.0, 5),
        make(3.0, 4.0, 400.0, 7),
    ]
}
/// The source's `BaseFeature p;`, which its sections then fill through setters.
/// `BaseFeature::new` with zero coordinates is exactly the default value.
fn base_feature() -> BaseFeature {
    BaseFeature::new(0.0, 0.0, 0.0)
}
/// The source's `Feature p;`, for the same reason.
fn plain_feature() -> Feature {
    Feature::new(0.0, 0.0, 0.0)
}
/// `BaseFeature::QualityLess`, which is `<` on the quality value.
fn quality_less(left: f32, right: f32) -> bool {
    left < right
}
/// `ConsensusFeature::SizeLess`, which is `<` on the number of handles.
fn size_less(left: usize, right: usize) -> bool {
    left < right
}
fn hit(sequence: &str, score: f64) -> PeptideHit {
    PeptideHit::new(score, 0, 0, AASequence::parse(sequence).unwrap()).unwrap()
}
fn identification(hits: Vec<PeptideHit>) -> PeptideIdentification {
    PeptideIdentification {
        hits,
        ..PeptideIdentification::default()
    }
}

// ---------------------------------------------------------------------------
// BaseFeature_test.cpp
// ---------------------------------------------------------------------------

/// `BaseFeature()` (line 32): default construction.
#[test]
fn bf_default_constructor() {
    let feature = BaseFeature::default();
    assert_eq!(feature.rt, 0.0);
    assert_eq!(feature.mz, 0.0);
    assert_eq!(feature.intensity, 0.0);
    assert_eq!(feature.quality, 0.0);
    assert_eq!(feature.charge, 0);
    assert_eq!(feature.width, 0.0);
    assert!(feature.peptide_identifications.is_empty());
    assert!(feature.primary_id.is_none());
    assert!(feature.id_matches.is_empty());
}

/// `~BaseFeature()` (line 39): the source deletes a heap instance. Rust drops
/// on scope exit; the observable property is that a clone is independent.
#[test]
fn bf_destructor() {
    let mut feature = BaseFeature::new(1.0, 2.0, 3.0);
    let copy = feature.clone();
    drop(feature.peptide_identifications);
    feature = BaseFeature::default();
    assert_eq!(feature, BaseFeature::default());
    assert_eq!(copy.rt, 1.0);
}

/// `QualityType getQuality() const` (line 45).
#[test]
fn bf_get_quality() {
    assert_eq!(BaseFeature::default().quality, 0.0);
}

/// `void setQuality(QualityType q)` (line 51).
#[test]
fn bf_set_quality() {
    let mut feature = base_feature();
    feature.quality = 123.456;
    close_f32(feature.quality, 123.456);
    feature.quality = -0.12345;
    close_f32(feature.quality, -0.12345);
    feature.quality = 0.0;
    assert_eq!(feature.quality, 0.0);
    // A signed quality is stored and remains valid; only non-finite fails.
    feature.validate().unwrap();
    feature.quality = f32::NAN;
    assert!(feature.validate().is_err());
}

/// `WidthType getWidth() const` (line 61).
#[test]
fn bf_get_width() {
    assert_eq!(BaseFeature::default().width, 0.0);
}

/// `void setWidth(WidthType fwhm)` (line 67). The source stores a negative
/// FWHM; the checked setter rejects it, so the section's -0.12345 is asserted
/// as an error plus the unchecked field assignment the port still allows.
#[test]
fn bf_set_width() {
    let mut feature = base_feature();
    feature.set_width(123.456).unwrap();
    close_f32(feature.width, 123.456);
    assert_eq!(
        feature.metadata["FWHM"].as_f64().unwrap(),
        f64::from(123.456f32)
    );
    let before = feature.clone();
    assert!(feature.set_width(-0.12345).is_err());
    assert_eq!(feature, before);
    feature.width = -0.12345;
    close_f32(feature.width, -0.12345);
    assert!(feature.validate().is_err());
    feature.set_width(0.0).unwrap();
    assert_eq!(feature.width, 0.0);
}

/// `[EXTRA] IntensityType getIntensity() const` (line 77).
#[test]
fn bf_get_intensity_const() {
    assert_eq!(BaseFeature::default().intensity, 0.0);
}

/// `[EXTRA] const PositionType& getPosition() const` (line 82).
#[test]
fn bf_get_position_const() {
    let feature = BaseFeature::default();
    assert_eq!(feature.position(), Point2D::new(0.0, 0.0));
}

/// `[EXTRA] IntensityType& getIntensity()` (line 88).
#[test]
fn bf_mutable_intensity() {
    let mut feature = base_feature();
    assert_eq!(feature.intensity, 0.0);
    feature.intensity = 123.456;
    close_f32(feature.intensity, 123.456);
    feature.intensity = -0.12345;
    close_f32(feature.intensity, -0.12345);
    feature.intensity = 0.0;
    assert_eq!(feature.intensity, 0.0);
}

/// `[EXTRA] PositionType& getPosition()` (line 99).
#[test]
fn bf_mutable_position() {
    let mut feature = base_feature();
    assert_eq!(feature.position(), Point2D::new(0.0, 0.0));
    feature.rt = 1.0;
    feature.mz = 2.0;
    let position = feature.position();
    assert_eq!(position, Point2D::new(1.0, 2.0));
}

/// `const ChargeType& getCharge() const` (line 113).
#[test]
fn bf_get_charge() {
    assert_eq!(BaseFeature::default().charge, 0);
}

/// `void setCharge(const ChargeType& ch)` (line 121).
#[test]
fn bf_set_charge() {
    let mut feature = base_feature();
    assert_eq!(feature.charge, 0);
    feature.charge = 17;
    assert_eq!(feature.charge, 17);
}

/// `BaseFeature(const BaseFeature& feature)` (line 130).
#[test]
fn bf_copy_constructor() {
    let mut feature = base_feature();
    feature.intensity = 123.456;
    feature.rt = 21.21;
    feature.mz = 22.22;
    feature
        .metadata
        .insert("cluster_id".into(), 4711_i64.into());
    feature.quality = 0.9;
    let copy = feature.clone();
    close_f32(copy.intensity, 123.456);
    close(copy.rt, 21.21);
    close(copy.mz, 22.22);
    assert_eq!(copy.metadata["cluster_id"].as_i64().unwrap(), 4711);
    close_f32(copy.quality, 0.9);
}

/// `BaseFeature(BaseFeature&& feature)` (line 154). The source asserts the move
/// constructor is `noexcept` so `std::vector` moves instead of copying; a Rust
/// move is a memcpy that cannot fail, so the value checks are what remains.
#[test]
fn bf_move_constructor() {
    let mut feature = base_feature();
    feature.intensity = 123.456;
    feature.rt = 21.21;
    feature.mz = 22.22;
    feature
        .metadata
        .insert("cluster_id".into(), 4711_i64.into());
    feature.quality = 0.9;
    let original = feature.clone();
    let moved = feature;
    close_f32(moved.intensity, 123.456);
    close(moved.rt, 21.21);
    close(moved.mz, 22.22);
    assert_eq!(moved.metadata["cluster_id"].as_i64().unwrap(), 4711);
    close_f32(moved.quality, 0.9);
    assert_eq!(moved, original);
}

/// `BaseFeature(const Peak2D& point)` (line 184). There is no converting
/// constructor in the port; the section's semantics are the field copy plus
/// zero-initialised quality, charge, width and identifications.
#[test]
fn bf_from_peak2d() {
    let point = Peak2D::new(1.23, 4.56, 7.89);
    let feature = BaseFeature {
        rt: point.rt(),
        mz: point.mz(),
        intensity: point.intensity,
        ..BaseFeature::default()
    };
    close(feature.rt, 1.23);
    close(feature.mz, 4.56);
    close_f32(feature.intensity, 7.89);
    assert_eq!(feature.quality, 0.0);
    assert_eq!(feature.charge, 0);
    assert_eq!(feature.width, 0.0);
    assert!(feature.peptide_identifications.is_empty());
}

/// `BaseFeature(const RichPeak2D& point)` (line 202). As above, with metadata.
#[test]
fn bf_from_rich_peak2d() {
    let mut point = RichPeak2D::new(1.23, 4.56, 7.89);
    point.metadata.insert("meta".into(), "test".into());
    let feature = BaseFeature {
        rt: point.peak.rt(),
        mz: point.peak.mz(),
        intensity: point.peak.intensity,
        unique_id: point.unique_id,
        metadata: point.metadata.clone(),
        ..BaseFeature::default()
    };
    close(feature.rt, 1.23);
    close(feature.mz, 4.56);
    close_f32(feature.intensity, 7.89);
    assert_eq!(feature.metadata["meta"].as_str().unwrap(), "test");
    assert_eq!(feature.quality, 0.0);
    assert_eq!(feature.charge, 0);
    assert_eq!(feature.width, 0.0);
    assert!(feature.peptide_identifications.is_empty());
}

/// `BaseFeature& operator=(const BaseFeature& rhs)` (line 222).
#[test]
fn bf_assignment_operator() {
    let mut feature = base_feature();
    feature.intensity = 123.456;
    feature.rt = 21.21;
    feature.mz = 22.22;
    feature.quality = 0.9;
    let mut copy = base_feature();
    assert_eq!(copy.intensity, 0.0);
    copy = feature.clone();
    close_f32(copy.intensity, 123.456);
    close(copy.rt, 21.21);
    close(copy.mz, 22.22);
    close_f32(copy.quality, 0.9);
}

/// `bool operator==(const BaseFeature& rhs) const` (line 244).
#[test]
fn bf_equality_operator() {
    let mut left = base_feature();
    let mut right = left.clone();
    assert_eq!(left, right);
    left.intensity = 5.0;
    left.quality = 0.9;
    assert_ne!(left, right);
    right.intensity = 5.0;
    right.quality = 0.9;
    assert_eq!(left, right);
    left.rt = 5.0;
    assert_ne!(left, right);
    right.rt = 5.0;
    assert_eq!(left, right);
    let peptides = vec![PeptideIdentification::default()];
    left.peptide_identifications = peptides.clone();
    assert_ne!(left, right);
    right.peptide_identifications = peptides;
    assert_eq!(left, right);
}

/// `bool operator!=(const BaseFeature& rhs) const` (line 270). Rust derives
/// `ne` from `eq`, so the section checks the same states through `assert_ne!`.
#[test]
fn bf_inequality_operator() {
    let mut left = base_feature();
    let mut right = left.clone();
    assert!(!(left != right));
    left.intensity = 5.0;
    assert!(left != right);
    right.intensity = 5.0;
    assert!(!(left != right));
    left.rt = 5.0;
    assert!(left != right);
    right.rt = 5.0;
    assert!(!(left != right));
    let peptides = vec![PeptideIdentification::default()];
    left.peptide_identifications = peptides.clone();
    assert!(left != right);
    right.peptide_identifications = peptides;
    assert!(!(left != right));
}

/// `[EXTRA] meta info with copy constructor` (line 292). The source addresses
/// metadata by registry index 2; the port has named keys only.
#[test]
fn bf_meta_info_with_copy_constructor() {
    let mut feature = base_feature();
    feature.metadata.insert("2".into(), "bla".into());
    let copy = feature.clone();
    assert_eq!(feature.metadata["2"].as_str().unwrap(), "bla");
    assert_eq!(copy.metadata["2"].as_str().unwrap(), "bla");
    feature.metadata.insert("2".into(), "bluff".into());
    assert_eq!(feature.metadata["2"].as_str().unwrap(), "bluff");
    assert_eq!(copy.metadata["2"].as_str().unwrap(), "bla");
}

/// `[EXTRA] meta info with assignment` (line 303).
#[test]
fn bf_meta_info_with_assignment() {
    let mut feature = base_feature();
    feature.metadata.insert("2".into(), "bla".into());
    let copy = feature.clone();
    assert_eq!(feature.metadata["2"].as_str().unwrap(), "bla");
    assert_eq!(copy.metadata["2"].as_str().unwrap(), "bla");
    feature.metadata.insert("2".into(), "bluff".into());
    assert_eq!(feature.metadata["2"].as_str().unwrap(), "bluff");
    assert_eq!(copy.metadata["2"].as_str().unwrap(), "bla");
}

/// `[BaseFeature::QualityLess] operator()(const BaseFeature&, const BaseFeature&)`
/// (line 314). The comparator is `<` on quality; the port's counterpart is the
/// quality ordering `FeatureMap::sort_by_quality` uses.
#[test]
fn bf_quality_less_feature_feature() {
    let mut first = base_feature();
    let mut second = base_feature();
    first.quality = 0.94;
    second.quality = 0.78;
    assert!(!quality_less(first.quality, second.quality));
    assert!(quality_less(second.quality, first.quality));
    let mut map = FeatureMap::from_features(vec![
        Feature::from(first.clone()),
        Feature::from(second.clone()),
    ]);
    map.sort_by_quality(false).unwrap();
    close_f32(map.features[0].quality, 0.78);
}

/// `[BaseFeature::QualityLess] operator()(const BaseFeature&, const QualityType&)`
/// (line 324).
#[test]
fn bf_quality_less_feature_value() {
    let mut first = base_feature();
    let mut second = base_feature();
    first.quality = 0.94;
    second.quality = 0.78;
    let right = first.quality;
    assert!(!quality_less(first.quality, right));
    assert!(quality_less(second.quality, right));
}

/// `[BaseFeature::QualityLess] operator()(const QualityType&, const BaseFeature&)`
/// (line 335).
#[test]
fn bf_quality_less_value_feature() {
    let mut first = base_feature();
    let mut second = base_feature();
    first.quality = 0.94;
    second.quality = 0.78;
    let left = second.quality;
    assert!(!quality_less(left, second.quality));
    assert!(quality_less(left, first.quality));
}

/// `[BaseFeature::QualityLess] operator()(const QualityType&, const QualityType&)`
/// (line 346).
#[test]
fn bf_quality_less_value_value() {
    let mut first = base_feature();
    let mut second = base_feature();
    first.quality = 0.94;
    second.quality = 0.78;
    let (left, right) = (first.quality, second.quality);
    assert!(!quality_less(left, right));
    assert!(quality_less(right, left));
}

/// `const PeptideIdentificationList& getPeptideIdentifications() const` (line 359).
#[test]
fn bf_get_peptide_identifications_const() {
    let feature = BaseFeature::default();
    assert_eq!(feature.peptide_identifications.len(), 0);
}

/// `void setPeptideIdentifications(const PeptideIdentificationList&)` (line 365).
#[test]
fn bf_set_peptide_identifications() {
    let mut feature = base_feature();
    feature.peptide_identifications = Vec::new();
    assert_eq!(feature.peptide_identifications.len(), 0);
    feature.peptide_identifications = vec![PeptideIdentification::default()];
    assert_eq!(feature.peptide_identifications.len(), 1);
}

/// `PeptideIdentificationList& getPeptideIdentifications()` (line 378).
#[test]
fn bf_get_peptide_identifications_mut() {
    let mut feature = base_feature();
    feature
        .peptide_identifications
        .resize(1, PeptideIdentification::default());
    assert_eq!(feature.peptide_identifications.len(), 1);
}

/// `AnnotationState getAnnotationState() const` (line 385).
#[test]
fn bf_get_annotation_state() {
    let mut feature = base_feature();
    assert_eq!(
        feature.annotation_state(None).unwrap(),
        AnnotationState::None
    );
    feature
        .peptide_identifications
        .resize(1, PeptideIdentification::default());
    assert_eq!(
        feature.annotation_state(None).unwrap(),
        AnnotationState::None
    );
    let first = hit("ABCDE", 0.0);
    feature.peptide_identifications[0].hits = vec![first.clone()];
    assert_eq!(
        feature.annotation_state(None).unwrap(),
        AnnotationState::Single
    );
    feature
        .peptide_identifications
        .resize(2, PeptideIdentification::default());
    feature.peptide_identifications[1].hits = vec![first];
    assert_eq!(
        feature.annotation_state(None).unwrap(),
        AnnotationState::MultipleSame
    );
    feature.peptide_identifications[1].hits = vec![hit("KRGH", 0.0)];
    assert_eq!(
        feature.annotation_state(None).unwrap(),
        AnnotationState::MultipleDivergent
    );
}

/// `sortPeptideIdentifications()` (line 409).
#[test]
fn bf_sort_peptide_identifications() {
    let mut feature = base_feature();
    feature.peptide_identifications = vec![
        identification(vec![hit("ABCDE", 0.8)]),
        identification(vec![hit("ABCDE", 0.5), hit("KRGH", 0.9)]),
        identification(Vec::new()),
    ];
    feature.sort_peptide_identifications().unwrap();
    assert_eq!(feature.peptide_identifications[0].hits[0].score, 0.9);
    assert!(feature.peptide_identifications[2].hits.is_empty());
    // The rest of the source order: the 0.8 identification follows the 0.9 one.
    assert_eq!(feature.peptide_identifications[1].hits[0].score, 0.8);
    // Hits inside each identification are sorted too, which the source only
    // achieves for elements its comparator happens to visit.
    assert_eq!(feature.peptide_identifications[0].hits[1].score, 0.5);
}

// ---------------------------------------------------------------------------
// Feature_test.cpp
// ---------------------------------------------------------------------------

/// `Feature()` (line 31).
#[test]
fn f_default_constructor() {
    let feature = Feature::default();
    assert_eq!(feature.quality, 0.0);
    assert_eq!(feature.quality_rt, 0.0);
    assert_eq!(feature.quality_mz, 0.0);
    assert!(feature.convex_hulls.is_empty());
    assert!(feature.subordinates.is_empty());
}

/// `~Feature()` (line 38): see `bf_destructor`; Rust drops the subordinate tree.
#[test]
fn f_destructor() {
    let mut feature = plain_feature();
    feature.subordinates.push(Feature::default());
    drop(feature);
    assert!(Feature::default().subordinates.is_empty());
}

/// `QualityType getOverallQuality() const` (line 44).
#[test]
fn f_get_overall_quality() {
    let mut feature = plain_feature();
    assert_eq!(feature.quality, 0.0);
    feature.quality = 123.456;
    close_f32(feature.quality, 123.456);
    feature.quality = -0.12345;
    close_f32(feature.quality, -0.12345);
    feature.quality = 0.0;
    assert_eq!(feature.quality, 0.0);
}

/// `void setOverallQuality(QualityType q)` (line 55).
#[test]
fn f_set_overall_quality() {
    let mut feature = plain_feature();
    feature.quality = 123.456;
    close_f32(feature.quality, 123.456);
    feature.quality = -0.12345;
    close_f32(feature.quality, -0.12345);
    feature.quality = 0.0;
    assert_eq!(feature.quality, 0.0);
}

/// `QualityType getQuality(Size index) const` (line 65). The port has named
/// `quality_rt`/`quality_mz` fields, so the section's
/// `TEST_PRECONDITION_VIOLATED(p.getQuality(10))` has no counterpart: an
/// out-of-range dimension index cannot be written.
#[test]
fn f_get_quality_by_index() {
    let mut feature = plain_feature();
    assert_eq!(feature.quality_rt, 0.0);
    feature.quality_rt = 123.456;
    close_f32(feature.quality_rt, 123.456);
    feature.quality_rt = -0.12345;
    close_f32(feature.quality_rt, -0.12345);
    feature.quality_rt = 0.0;
    assert_eq!(feature.quality_rt, 0.0);
    assert_eq!(feature.quality_mz, 0.0);
}

/// `void setQuality(Size index, QualityType q)` (line 78).
#[test]
fn f_set_quality_by_index() {
    let mut feature = plain_feature();
    feature.quality_mz = 123.456;
    close_f32(feature.quality_mz, 123.456);
    feature.quality_mz = -0.12345;
    close_f32(feature.quality_mz, -0.12345);
    feature.quality_mz = 0.0;
    assert_eq!(feature.quality_rt, 0.0);
    assert_eq!(feature.quality_mz, 0.0);
    // Non-finite dimension qualities are rejected by the checked operations.
    feature.quality_mz = f32::INFINITY;
    assert!(feature.validate().is_err());
}

/// `const vector<ConvexHull2D>& getConvexHulls() const` (line 97).
#[test]
fn f_get_convex_hulls_const() {
    assert_eq!(Feature::default().convex_hulls.len(), 0);
}

/// `vector<ConvexHull2D>& getConvexHulls()` (line 102).
#[test]
fn f_get_convex_hulls_mut() {
    let mut feature = plain_feature();
    feature.convex_hulls = source_hulls();
    assert_eq!(feature.convex_hulls.len(), 2);
    let first = feature.convex_hulls[0].hull_points();
    close(first[0].rt, 1.0);
    close(first[0].mz, 2.0);
    close(first[1].rt, 3.0);
    close(first[1].mz, 4.0);
    let second = feature.convex_hulls[1].hull_points();
    close(second[0].rt, 0.5);
    close(second[0].mz, 0.0);
    close(second[1].rt, 1.0);
    close(second[1].mz, 1.0);
}

/// `void setConvexHulls(const vector<ConvexHull2D>& hulls)` (line 116).
#[test]
fn f_set_convex_hulls() {
    let mut feature = plain_feature();
    feature.convex_hulls = source_hulls();
    assert_eq!(feature.convex_hulls.len(), 2);
    let first = feature.convex_hulls[0].hull_points();
    close(first[0].rt, 1.0);
    close(first[0].mz, 2.0);
    close(first[1].rt, 3.0);
    close(first[1].mz, 4.0);
    let second = feature.convex_hulls[1].hull_points();
    close(second[0].rt, 0.5);
    close(second[0].mz, 0.0);
    close(second[1].rt, 1.0);
    close(second[1].mz, 1.0);
}

/// `ConvexHull2D& getConvexHull() const` (line 130).
#[test]
fn f_get_convex_hull() {
    let mut feature = plain_feature();
    feature.convex_hulls = source_hulls();
    let bounds = feature.convex_hull().bounding_box().unwrap();
    close(bounds.min().rt, 0.5);
    close(bounds.min().mz, 0.0);
    close(bounds.max().rt, 3.0);
    close(bounds.max().mz, 4.0);
    let points = feature.convex_hull().hull_points();
    assert_eq!(points.len(), 4);
    close(points[0].rt, 0.5);
    close(points[0].mz, 0.0);
    close(points[1].rt, 3.0);
    close(points[1].mz, 0.0);
    close(points[2].rt, 3.0);
    close(points[2].mz, 4.0);
    close(points[3].rt, 0.5);
    close(points[3].mz, 4.0);
}

/// `bool encloses(double rt, double mz) const` (line 156).
#[test]
fn f_encloses() {
    let mut feature = plain_feature();
    assert_eq!(feature.convex_hull().bounding_box(), None);
    feature.convex_hulls = source_hulls_extended();
    for &(rt, mz, expected) in &[
        (0.0, 0.0, false),
        (1.0, 1.0, true),
        (2.0, 0.5, false),
        (2.0, 3.001, false),
        (2.0, 2.999, true),
        (2.0, 3.5, false),
        (4.0, 3.0, false),
        (1.5, 1.5, false),
        (2.0, 1.0, true),
        (0.5, 0.0, true),
        (3.0, 3.2, true),
    ] {
        assert_eq!(feature.encloses(rt, mz).unwrap(), expected, "{rt},{mz}");
    }
}

/// `Feature(const Feature& feature)` (line 174).
#[test]
fn f_copy_constructor() {
    let mut feature = plain_feature();
    feature.intensity = 123.456;
    feature.rt = 21.21;
    feature.mz = 22.22;
    feature
        .metadata
        .insert("cluster_id".into(), 4711_i64.into());
    feature.quality = 0.9;
    feature.quality_rt = 0.1;
    feature.quality_mz = 0.2;
    feature.convex_hulls = source_hulls_extended();
    let _ = feature.convex_hull();
    let copy = feature.clone();
    close_f32(copy.intensity, 123.456);
    close(copy.rt, 21.21);
    close(copy.mz, 22.22);
    assert_eq!(copy.metadata["cluster_id"].as_i64().unwrap(), 4711);
    close_f32(copy.quality, 0.9);
    close_f32(copy.quality_rt, 0.1);
    close_f32(copy.quality_mz, 0.2);
    assert_eq!(
        copy.convex_hull().hull_points().len(),
        feature.convex_hull().hull_points().len()
    );
    assert_eq!(copy.convex_hulls.len(), feature.convex_hulls.len());
}

/// `Feature(const Feature&& source)` (line 215). The source additionally
/// asserts the move constructor is `noexcept` and that the moved-from feature
/// has no hulls left; a Rust move consumes the value, so the moved-from
/// variable cannot be read at all, which is the stronger guarantee.
#[test]
fn f_move_constructor() {
    let mut feature = plain_feature();
    feature.intensity = 123.456;
    feature.rt = 21.21;
    feature.mz = 22.22;
    feature
        .metadata
        .insert("cluster_id".into(), 4711_i64.into());
    feature.quality = 0.9;
    feature.quality_rt = 0.1;
    feature.quality_mz = 0.2;
    feature.convex_hulls = source_hulls_extended();
    let _ = feature.convex_hull();
    let original = feature.clone();
    let moved = feature;
    close_f32(moved.intensity, 123.456);
    close(moved.rt, 21.21);
    close(moved.mz, 22.22);
    assert_eq!(moved.metadata["cluster_id"].as_i64().unwrap(), 4711);
    close_f32(moved.quality, 0.9);
    close_f32(moved.quality_rt, 0.1);
    close_f32(moved.quality_mz, 0.2);
    assert_eq!(
        moved.convex_hull().hull_points().len(),
        original.convex_hull().hull_points().len()
    );
    assert_eq!(moved.convex_hulls.len(), original.convex_hulls.len());
    assert_eq!(Feature::default().convex_hull().hull_points().len(), 0);
    assert_eq!(Feature::default().convex_hulls.len(), 0);
}

/// `Feature& operator=(const Feature& rhs)` (line 266). The source pre-computes
/// the overall hull on the destination to prove the recalculation flag is
/// copied; the port computes the hull on demand and caches nothing.
#[test]
fn f_assignment_operator() {
    let mut feature = plain_feature();
    feature.intensity = 123.456;
    feature.rt = 21.21;
    feature.mz = 22.22;
    feature.quality = 0.9;
    feature.quality_rt = 0.1;
    feature.quality_mz = 0.2;
    feature
        .metadata
        .insert("cluster_id".into(), 4712_i64.into());
    feature.convex_hulls = source_hulls_extended();
    let mut copy = plain_feature();
    let _ = copy.convex_hull();
    copy = feature.clone();
    close_f32(copy.intensity, 123.456);
    close(copy.rt, 21.21);
    close(copy.mz, 22.22);
    close_f32(copy.quality, 0.9);
    close_f32(copy.quality_rt, 0.1);
    close_f32(copy.quality_mz, 0.2);
    assert_eq!(
        copy.convex_hull().hull_points().len(),
        feature.convex_hull().hull_points().len()
    );
    assert_eq!(copy.convex_hulls.len(), feature.convex_hulls.len());
}

/// `bool operator==(const Feature& rhs) const` (line 304).
#[test]
fn f_equality_operator() {
    let mut left = plain_feature();
    let mut right = left.clone();
    assert_eq!(left, right);
    left.intensity = 5.0;
    left.quality = 0.9;
    left.quality_rt = 0.1;
    assert_ne!(left, right);
    right.intensity = 5.0;
    right.quality = 0.9;
    right.quality_rt = 0.1;
    assert_eq!(left, right);
    left.rt = 5.0;
    assert_ne!(left, right);
    right.rt = 5.0;
    assert_eq!(left, right);
}

/// `[EXTRA] Feature& operator!=(const Feature& rhs)` (line 324).
#[test]
fn f_inequality_operator() {
    let mut left = plain_feature();
    let mut right = left.clone();
    assert!(!(left != right));
    left.intensity = 5.0;
    assert!(left != right);
    right.intensity = 5.0;
    assert!(!(left != right));
    left.rt = 5.0;
    assert!(left != right);
    right.rt = 5.0;
    assert!(!(left != right));
}

/// `[EXTRA] meta info with copy constructor` (line 340).
#[test]
fn f_meta_info_with_copy_constructor() {
    let mut feature = plain_feature();
    feature.metadata.insert("2".into(), "bla".into());
    let copy = feature.clone();
    assert_eq!(feature.metadata["2"].as_str().unwrap(), "bla");
    assert_eq!(copy.metadata["2"].as_str().unwrap(), "bla");
    feature.metadata.insert("2".into(), "bluff".into());
    assert_eq!(feature.metadata["2"].as_str().unwrap(), "bluff");
    assert_eq!(copy.metadata["2"].as_str().unwrap(), "bla");
}

/// `[EXTRA] meta info with assignment` (line 353).
#[test]
fn f_meta_info_with_assignment() {
    let mut feature = plain_feature();
    feature.metadata.insert("2".into(), "bla".into());
    let copy = feature.clone();
    assert_eq!(feature.metadata["2"].as_str().unwrap(), "bla");
    assert_eq!(copy.metadata["2"].as_str().unwrap(), "bla");
    feature.metadata.insert("2".into(), "bluff".into());
    assert_eq!(feature.metadata["2"].as_str().unwrap(), "bluff");
    assert_eq!(copy.metadata["2"].as_str().unwrap(), "bla");
}

/// `std::vector<Feature>& getSubordinates()` (line 367): `NOT_TESTABLE` in the
/// source, which defers to the const section below.
#[test]
fn f_get_subordinates_mut() {
    let mut feature = plain_feature();
    feature.subordinates.push(Feature::new(1101.0, 1102.0, 0.0));
    assert_eq!(feature.subordinates.len(), 1);
}

/// `void setSubordinates(const std::vector<Feature>& rhs)` (line 374):
/// `NOT_TESTABLE` in the source; see the const section below.
#[test]
fn f_set_subordinates() {
    let mut feature = plain_feature();
    feature.subordinates = vec![Feature::new(1101.0, 1102.0, 0.0)];
    assert_eq!(feature.subordinates.len(), 1);
}

/// `const std::vector<Feature>& getSubordinates() const` (line 381).
#[test]
fn f_get_subordinates_const() {
    let mut first = plain_feature();
    first.rt = 1001.0;
    first.mz = 1002.0;
    first.charge = 1003;
    let original = first.clone();
    let sub = |rt: f64, mz: f64| {
        let mut feature = plain_feature();
        feature.rt = rt;
        feature.mz = mz;
        feature
    };
    assert!(first.subordinates.is_empty());
    first.subordinates.push(sub(1101.0, 1102.0));
    assert_eq!(first.subordinates.len(), 1);
    first.subordinates.push(sub(1201.0, 1202.0));
    assert_eq!(first.subordinates.len(), 2);
    first.subordinates.push(sub(1301.0, 1302.0));
    assert_eq!(first.subordinates.len(), 3);
    assert_eq!(first.rt, 1001.0);
    assert_eq!(first.subordinates[0].rt, 1101.0);
    assert_eq!(first.subordinates[1].rt, 1201.0);
    assert_eq!(first.subordinates[2].rt, 1301.0);
    assert_eq!(first.mz, 1002.0);
    assert_eq!(first.subordinates[0].mz, 1102.0);
    assert_eq!(first.subordinates[1].mz, 1202.0);
    assert_eq!(first.subordinates[2].mz, 1302.0);
    assert_ne!(first, original);
    let with_subordinates = first.clone();
    assert_eq!(with_subordinates, first);
    first.subordinates.clear();
    assert_eq!(first, original);
    let mut second = plain_feature();
    second.rt = 1001.0;
    second.mz = 1002.0;
    second.charge = 1003;
    assert!(!with_subordinates.subordinates.is_empty());
    second.subordinates = with_subordinates.subordinates.clone();
    assert_eq!(second, with_subordinates);
}

/// `template<typename Type> Size applyMemberFunction(Size (Type::*)())` (line 429).
#[test]
fn f_apply_member_function_mutable() {
    let mut feature = plain_feature();
    assert_eq!(
        feature
            .for_each_unique_id(|id| usize::from(id.has_invalid_unique_id()))
            .unwrap(),
        1
    );
    feature.unique_id = 42;
    assert_eq!(
        feature
            .for_each_unique_id(|id| usize::from(id.has_invalid_unique_id()))
            .unwrap(),
        0
    );
}

/// `template<typename Type> Size applyMemberFunction(Size (Type::*)() const) const`
/// (line 439).
#[test]
fn f_apply_member_function_const() {
    let mut feature = plain_feature();
    assert_eq!(
        feature
            .count_unique_ids(|id| usize::from(id.has_invalid_unique_id()))
            .unwrap(),
        1
    );
    feature.unique_id = 42;
    assert_eq!(
        feature
            .count_unique_ids(|id| usize::from(id.has_invalid_unique_id()))
            .unwrap(),
        0
    );
}

// ---------------------------------------------------------------------------
// ConsensusFeature_test.cpp
// ---------------------------------------------------------------------------

/// `ConsensusFeature()` (line 32).
#[test]
fn cf_default_constructor() {
    let consensus = ConsensusFeature::new();
    assert!(consensus.is_empty());
    assert!(consensus.ratios().is_empty());
    assert_eq!(consensus.base, BaseFeature::default());
}

/// `virtual ~ConsensusFeature()` (line 37): Rust drops handles and ratios.
#[test]
fn cf_destructor() {
    let mut consensus = ConsensusFeature::new();
    consensus.add_ratio(Ratio::default()).unwrap();
    drop(consensus);
    assert!(ConsensusFeature::new().ratios().is_empty());
}

/// `[ConsensusFeature::SizeLess] operator()(ConsensusFeature const&, ConsensusFeature const&)`
/// (line 60).
#[test]
fn cf_size_less_feature_feature() {
    let [first_feature, second_feature, third_feature] = consensus_fixture();
    let mut first = ConsensusFeature::from(first_feature.base.clone());
    first.insert(FeatureHandle::new(1, &first_feature)).unwrap();
    first.insert(FeatureHandle::new(2, &third_feature)).unwrap();
    let mut second = ConsensusFeature::from(second_feature.base.clone());
    second
        .insert(FeatureHandle::new(1, &second_feature))
        .unwrap();
    assert!(!size_less(first.len(), second.len()));
    assert!(size_less(second.len(), first.len()));
}

/// `[ConsensusFeature::SizeLess] operator()(ConsensusFeature const&, UInt64 const&)`
/// (line 74).
#[test]
fn cf_size_less_feature_value() {
    let [first_feature, second_feature, third_feature] = consensus_fixture();
    let mut first = ConsensusFeature::from(first_feature.base.clone());
    first.insert(FeatureHandle::new(1, &first_feature)).unwrap();
    first.insert(FeatureHandle::new(2, &third_feature)).unwrap();
    let mut second = ConsensusFeature::from(first_feature.base.clone());
    second
        .insert(FeatureHandle::new(1, &first_feature))
        .unwrap();
    second
        .insert(FeatureHandle::new(2, &second_feature))
        .unwrap();
    second
        .insert(FeatureHandle::new(3, &third_feature))
        .unwrap();
    let right = second.len();
    assert!(size_less(first.len(), right));
    assert!(!size_less(second.len(), right));
}

/// `[ConsensusFeature::SizeLess] operator()(UInt64 const&, ConsensusFeature const&)`
/// (line 92).
#[test]
fn cf_size_less_value_feature() {
    let [first_feature, second_feature, third_feature] = consensus_fixture();
    let mut first = ConsensusFeature::from(first_feature.base.clone());
    first.insert(FeatureHandle::new(1, &first_feature)).unwrap();
    first.insert(FeatureHandle::new(2, &third_feature)).unwrap();
    let mut second = ConsensusFeature::from(first_feature.base.clone());
    second
        .insert(FeatureHandle::new(1, &first_feature))
        .unwrap();
    second
        .insert(FeatureHandle::new(2, &second_feature))
        .unwrap();
    second
        .insert(FeatureHandle::new(3, &third_feature))
        .unwrap();
    let left = first.len();
    assert!(!size_less(left, first.len()));
    assert!(size_less(left, second.len()));
}

/// `[ConsensusFeature::SizeLess] operator()(const UInt64&, const UInt64&)` (line 110).
#[test]
fn cf_size_less_value_value() {
    let [first_feature, second_feature, third_feature] = consensus_fixture();
    let mut first = ConsensusFeature::from(first_feature.base.clone());
    first.insert(FeatureHandle::new(1, &first_feature)).unwrap();
    first.insert(FeatureHandle::new(2, &third_feature)).unwrap();
    let mut second = ConsensusFeature::from(first_feature.base.clone());
    second
        .insert(FeatureHandle::new(1, &first_feature))
        .unwrap();
    second
        .insert(FeatureHandle::new(2, &second_feature))
        .unwrap();
    second
        .insert(FeatureHandle::new(3, &third_feature))
        .unwrap();
    let (left, right) = (first.len(), second.len());
    assert!(size_less(left, right));
    assert!(!size_less(right, left));
    // The map-level counterpart sorts by descending size; the two consensus
    // features need distinct IDs because a map indexes them by unique ID.
    second.unique_id = 11;
    let mut map = ConsensusMap::from_features(vec![first, second]);
    map.sort_by_size().unwrap();
    assert_eq!(map.features[0].len(), 3);
}

/// `[ConsensusFeature::MapsLess] operator()(ConsensusFeature const&, ConsensusFeature const&)`
/// (line 128).
#[test]
fn cf_maps_less() {
    let [first_feature, second_feature, third_feature] = consensus_fixture();
    let mut first = ConsensusFeature::from(first_feature.base.clone());
    first.insert(FeatureHandle::new(1, &first_feature)).unwrap();
    first.insert(FeatureHandle::new(2, &third_feature)).unwrap();
    let mut second = ConsensusFeature::from(first_feature.base.clone());
    second
        .insert(FeatureHandle::new(3, &first_feature))
        .unwrap();
    second
        .insert(FeatureHandle::new(4, &second_feature))
        .unwrap();
    second
        .insert(FeatureHandle::new(5, &third_feature))
        .unwrap();
    let maps_less = |left: &ConsensusFeature, right: &ConsensusFeature| {
        left.handles()
            .iter()
            .map(FeatureHandle::key)
            .lt(right.handles().iter().map(FeatureHandle::key))
    };
    assert!(!maps_less(&first, &first));
    assert!(maps_less(&first, &second));
    assert!(!maps_less(&second, &first));
    assert!(!maps_less(&second, &second));
    second.unique_id = 11;
    let mut map = ConsensusMap::from_features(vec![second, first]);
    map.sort_by_maps().unwrap();
    assert_eq!(map.features[0].handles()[0].map_index, 1);
}

/// `ConsensusFeature& operator=(const ConsensusFeature& rhs)` (line 146).
#[test]
fn cf_assignment_operator() {
    let [feature, ..] = consensus_fixture();
    let mut consensus = ConsensusFeature::from(feature.base.clone());
    consensus.insert(FeatureHandle::new(1, &feature)).unwrap();
    let mut copy = ConsensusFeature::new();
    assert!(copy.is_empty());
    copy = consensus.clone();
    close(copy.rt, 1.0);
    close(copy.mz, 2.0);
    close_f32(copy.intensity, 200.0);
    assert_eq!(copy.handles()[0].map_index, 1);
    assert_eq!(copy.handles()[0].unique_id, 3);
    assert_eq!(copy.handles()[0].intensity, 200.0);
}

/// `ConsensusFeature(const ConsensusFeature& rhs)` (line 161).
#[test]
fn cf_copy_constructor() {
    let [feature, ..] = consensus_fixture();
    let mut consensus = ConsensusFeature::from(feature.base.clone());
    consensus.insert(FeatureHandle::new(1, &feature)).unwrap();
    let copy = consensus.clone();
    close(copy.rt, 1.0);
    close(copy.mz, 2.0);
    close_f32(copy.intensity, 200.0);
    assert_eq!(copy.handles()[0].map_index, 1);
    assert_eq!(copy.handles()[0].unique_id, 3);
    assert_eq!(copy.handles()[0].intensity, 200.0);
}

/// `ConsensusFeature(ConsensusFeature&& rhs)` (line 176). The source's
/// `noexcept` assertion has no Rust counterpart; a move cannot fail.
#[test]
fn cf_move_constructor() {
    let [feature, ..] = consensus_fixture();
    let mut consensus = ConsensusFeature::from(feature.base.clone());
    consensus.insert(FeatureHandle::new(1, &feature)).unwrap();
    let moved = consensus;
    close(moved.rt, 1.0);
    close(moved.mz, 2.0);
    close_f32(moved.intensity, 200.0);
    assert_eq!(moved.handles()[0].map_index, 1);
    assert_eq!(moved.handles()[0].unique_id, 3);
    assert_eq!(moved.handles()[0].intensity, 200.0);
}

/// `void insert(const HandleSetType& handle_set)` (line 200).
#[test]
fn cf_insert_handle_set() {
    let handles: Vec<FeatureHandle> = (0..3u64)
        .map(|i| FeatureHandle {
            map_index: i + 10,
            unique_id: i + 1000,
            rt: i as f64 * 77.7,
            ..FeatureHandle::default()
        })
        .collect();
    let mut consensus = ConsensusFeature::new();
    consensus.set_handles(handles).unwrap();
    assert_eq!(consensus.len(), 3);
    assert_eq!(consensus.handles().first().unwrap().map_index, 10);
    assert_eq!(consensus.handles().last().unwrap().map_index, 12);
}

/// `void insert(UInt64 map_index, const Peak2D& element, UInt64 element_index)`
/// (line 219).
#[test]
fn cf_insert_peak_with_element_index() {
    let mut consensus = ConsensusFeature::new();
    for i in 0..3u64 {
        let point = Peak2D::new(i as f64 * 77.7, 0.0, 0.0);
        consensus
            .insert(FeatureHandle::from_peak(10 - i, point, i + 1000))
            .unwrap();
        assert_eq!(consensus.len(), usize::try_from(i).unwrap() + 1);
        // The source reads begin(), which is the lowest map index, i.e. the
        // handle just inserted because map indices count down.
        let first = consensus.handles()[0];
        close(first.rt, i as f64 * 77.7);
        assert_eq!(first.map_index, 10 - i);
        assert_eq!(first.unique_id, i + 1000);
    }
}

/// `void insert(UInt64 map_index, const BaseFeature& element)` (line 234).
#[test]
fn cf_insert_map_index_and_base_feature() {
    let mut consensus = ConsensusFeature::new();
    for i in 0..3u64 {
        let mut element = base_feature();
        element.rt = i as f64 * 77.7;
        element.charge = 2 * i32::try_from(i).unwrap();
        element.unique_id = i + 1000;
        consensus
            .insert(FeatureHandle::new(10 - i, &element))
            .unwrap();
        assert_eq!(consensus.len(), usize::try_from(i).unwrap() + 1);
        let first = consensus.handles()[0];
        close(first.rt, i as f64 * 77.7);
        assert_eq!(first.charge, 2 * i32::try_from(i).unwrap());
        assert_eq!(first.map_index, 10 - i);
        assert_eq!(first.unique_id, i + 1000);
    }
}

/// `ConsensusFeature(const BaseFeature& feature)` (line 251).
#[test]
fn cf_from_base_feature() {
    let mut feature = base_feature();
    feature.charge = -17;
    feature.rt = 44324.6;
    feature.mz = 867.4;
    feature.unique_id = 23;
    feature
        .peptide_identifications
        .resize(1, PeptideIdentification::default());
    let consensus = ConsensusFeature::from(feature);
    assert_eq!(consensus.rt, 44324.6);
    assert_eq!(consensus.mz, 867.4);
    assert_eq!(consensus.charge, -17);
    assert_eq!(consensus.peptide_identifications.len(), 1);
    assert!(consensus.is_empty());
}

/// `ConsensusFeature(UInt64 map_index, const BaseFeature& element)` (line 268).
/// The source constructor stamps `map_index` onto each attached peptide
/// identification through `BaseFeature(const BaseFeature&, UInt64)`;
/// `ConsensusFeature::from_feature` does not, so the port applies
/// `clone_with_map_index` explicitly.
#[test]
fn cf_from_map_index_and_base_feature() {
    let mut feature = base_feature();
    feature.charge = -17;
    feature.rt = 44324.6;
    feature.mz = 867.4;
    feature.intensity = 1000.0;
    feature.unique_id = 23;
    feature
        .peptide_identifications
        .resize(1, PeptideIdentification::default());
    let consensus =
        ConsensusFeature::from_feature(99, &feature.clone_with_map_index(99).unwrap()).unwrap();
    assert_eq!(consensus.rt, 44324.6);
    assert_eq!(consensus.mz, 867.4);
    assert_eq!(consensus.charge, -17);
    assert_eq!(consensus.peptide_identifications.len(), 1);
    assert_eq!(
        consensus.peptide_identifications[0].metadata["map_index"]
            .as_i64()
            .unwrap(),
        99
    );
    assert_eq!(consensus.handles()[0].map_index, 99);
    assert_eq!(consensus.handles()[0].unique_id, 23);
    assert_eq!(consensus.handles()[0].intensity, 1000.0);
}

/// `[EXTRA] ConsensusFeature(UInt64 map_index, const Feature& element)` (line 290).
#[test]
fn cf_from_map_index_and_feature() {
    let [feature, ..] = consensus_fixture();
    let mut consensus = ConsensusFeature::from_feature(1, &feature).unwrap();
    consensus.unique_id = 3;
    close(consensus.rt, 1.0);
    close(consensus.mz, 2.0);
    close_f32(consensus.intensity, 200.0);
    assert_eq!(consensus.handles()[0].map_index, 1);
    assert_eq!(consensus.handles()[0].unique_id, 3);
    assert_eq!(consensus.handles()[0].intensity, 200.0);
}

/// `ConsensusFeature(UInt64 map_index, const Peak2D& element, UInt64 element_index)`
/// (line 303).
#[test]
fn cf_from_map_index_peak_and_element_index() {
    let point = Peak2D::new(0.0, 0.0, -17.0);
    let mut consensus = ConsensusFeature::from(BaseFeature {
        rt: point.rt(),
        mz: point.mz(),
        intensity: point.intensity,
        ..BaseFeature::default()
    });
    consensus
        .insert(FeatureHandle::from_peak(99, point, 23))
        .unwrap();
    assert_eq!(consensus.handles()[0].map_index, 99);
    assert_eq!(consensus.handles()[0].unique_id, 23);
    assert_eq!(consensus.handles()[0].intensity, -17.0);
}

/// `[EXTRA] ConsensusFeature(UInt64 map_index, const ConsensusFeature& element)`
/// (line 315).
#[test]
fn cf_from_map_index_and_consensus_feature() {
    let mut element = ConsensusFeature::new();
    element.unique_id = 23;
    element.intensity = -17.0;
    let consensus = ConsensusFeature::from_feature(99, &element.base).unwrap();
    assert_eq!(consensus.handles()[0].map_index, 99);
    assert_eq!(consensus.handles()[0].unique_id, 23);
    assert_eq!(consensus.handles()[0].intensity, -17.0);
}

/// `DRange<1> getIntensityRange() const` (line 328), including the all-negative
/// regression the source comment describes.
#[test]
fn cf_get_intensity_range() {
    let mut consensus = ConsensusFeature::new();
    let mut feature = plain_feature();
    feature.intensity = 0.0;
    feature.unique_id = 0;
    consensus.insert(FeatureHandle::new(0, &feature)).unwrap();
    feature.unique_id = 1;
    feature.intensity = 200.0;
    consensus.insert(FeatureHandle::new(0, &feature)).unwrap();
    let range = consensus.handle_ranges().intensity.unwrap();
    close(range.min, 0.0);
    close(range.max, 200.0);
    let mut negative = ConsensusFeature::new();
    let mut feature = plain_feature();
    feature.intensity = -50.0;
    feature.unique_id = 0;
    negative.insert(FeatureHandle::new(0, &feature)).unwrap();
    feature.unique_id = 1;
    feature.intensity = -10.0;
    negative.insert(FeatureHandle::new(0, &feature)).unwrap();
    let range = negative.handle_ranges().intensity.unwrap();
    close(range.min, -50.0);
    close(range.max, -10.0);
}

/// `DRange<2> getPositionRange() const` (line 358), including the all-negative
/// regression.
#[test]
fn cf_get_position_range() {
    let mut consensus = ConsensusFeature::new();
    let mut feature = plain_feature();
    feature.rt = 1.0;
    feature.mz = 500.0;
    feature.unique_id = 0;
    consensus.insert(FeatureHandle::new(0, &feature)).unwrap();
    feature.rt = 1000.0;
    feature.mz = 1500.0;
    feature.unique_id = 1;
    consensus.insert(FeatureHandle::new(0, &feature)).unwrap();
    let ranges = consensus.handle_ranges();
    close(ranges.rt.unwrap().min, 1.0);
    close(ranges.rt.unwrap().max, 1000.0);
    close(ranges.mz.unwrap().min, 500.0);
    close(ranges.mz.unwrap().max, 1500.0);
    let mut negative = ConsensusFeature::new();
    let mut feature = plain_feature();
    feature.rt = -1000.0;
    feature.mz = -1500.0;
    feature.unique_id = 0;
    negative.insert(FeatureHandle::new(0, &feature)).unwrap();
    feature.rt = -1.0;
    feature.mz = -500.0;
    feature.unique_id = 1;
    negative.insert(FeatureHandle::new(0, &feature)).unwrap();
    let ranges = negative.handle_ranges();
    close(ranges.rt.unwrap().min, -1000.0);
    close(ranges.rt.unwrap().max, -1.0);
    close(ranges.mz.unwrap().min, -1500.0);
    close(ranges.mz.unwrap().max, -500.0);
}

/// `const HandleSetType& getFeatures() const` (line 396).
#[test]
fn cf_get_features() {
    let [feature, ..] = consensus_fixture();
    let mut consensus = ConsensusFeature::new();
    consensus.insert(FeatureHandle::new(2, &feature)).unwrap();
    let copy = consensus.clone();
    let group = copy.handles();
    assert_eq!(group[0].map_index, 2);
    assert_eq!(group[0].unique_id, 3);
    assert_eq!(group[0].intensity, 200.0);
}

/// `std::vector<FeatureHandle> getFeatureList() const` (line 409).
#[test]
fn cf_get_feature_list() {
    let [feature, ..] = consensus_fixture();
    let mut consensus = ConsensusFeature::new();
    consensus.insert(FeatureHandle::new(2, &feature)).unwrap();
    let copy = consensus.clone();
    let group = copy.handles().to_vec();
    assert_eq!(group.len(), 1);
    assert_eq!(group[0].map_index, 2);
    assert_eq!(group[0].unique_id, 3);
    assert_eq!(group[0].intensity, 200.0);
}

/// `void insert(const ConsensusFeature& cf)` (line 424).
#[test]
fn cf_insert_consensus_feature() {
    let [feature, ..] = consensus_fixture();
    let mut source = ConsensusFeature::new();
    let mut first = FeatureHandle::new(2, &feature);
    first.unique_id = 3;
    let mut second = FeatureHandle::new(4, &feature);
    second.unique_id = 5;
    source.insert(first).unwrap();
    source.insert(second).unwrap();
    let mut target = ConsensusFeature::new();
    target.merge(&source).unwrap();
    let handles = target.handles();
    assert_eq!(handles[0].map_index, 2);
    assert_eq!(handles[0].unique_id, 3);
    assert_eq!(handles[0].intensity, 200.0);
    assert_eq!(handles[1].map_index, 4);
    assert_eq!(handles[1].unique_id, 5);
    assert_eq!(handles[1].intensity, 200.0);
    assert_eq!(handles.len(), 2);
}

/// `void insert(const FeatureHandle& handle)` (line 447).
#[test]
fn cf_insert_feature_handle() {
    let [feature, ..] = consensus_fixture();
    let mut consensus = ConsensusFeature::new();
    let mut first = FeatureHandle::new(2, &feature);
    first.unique_id = 3;
    let mut second = FeatureHandle::new(4, &feature);
    second.unique_id = 5;
    consensus.insert(first).unwrap();
    consensus.insert(second).unwrap();
    let handles = consensus.handles();
    assert_eq!(handles[0].map_index, 2);
    assert_eq!(handles[0].unique_id, 3);
    assert_eq!(handles[0].intensity, 200.0);
    assert_eq!(handles[1].map_index, 4);
    assert_eq!(handles[1].unique_id, 5);
    assert_eq!(handles[1].intensity, 200.0);
    assert_eq!(handles.len(), 2);
}

/// `void insert(UInt64 map_index, const BaseFeature& element)` (line 468), the
/// source's second section for the same overload.
#[test]
fn cf_insert_map_index_and_base_feature_again() {
    let [feature, ..] = consensus_fixture();
    let mut consensus = ConsensusFeature::new();
    consensus.insert(FeatureHandle::new(2, &feature)).unwrap();
    let handles = consensus.handles();
    assert_eq!(handles[0].map_index, 2);
    assert_eq!(handles[0].unique_id, 3);
    assert_eq!(handles[0].intensity, 200.0);
    assert_eq!(handles.len(), 1);
}

/// `void computeConsensus()` (line 481).
#[test]
fn cf_compute_consensus() {
    let [first, second, third] = consensus_fixture();
    let mut consensus = ConsensusFeature::new();
    consensus.insert(FeatureHandle::new(2, &first)).unwrap();
    consensus.compute_consensus().unwrap();
    close_f32(consensus.intensity, 200.0);
    close(consensus.rt, 1.0);
    close(consensus.mz, 2.0);
    consensus.insert(FeatureHandle::new(4, &second)).unwrap();
    consensus.compute_consensus().unwrap();
    close_f32(consensus.intensity, 250.0);
    close(consensus.rt, 1.5);
    close(consensus.mz, 2.5);
    consensus.insert(FeatureHandle::new(6, &third)).unwrap();
    consensus.compute_consensus().unwrap();
    close_f32(consensus.intensity, 300.0);
    close(consensus.rt, 2.0);
    close(consensus.mz, 3.0);
}

/// `void computeMonoisotopicConsensus()` (line 503).
#[test]
fn cf_compute_monoisotopic_consensus() {
    let [first, second, third] = consensus_fixture();
    let mut consensus = ConsensusFeature::new();
    consensus.insert(FeatureHandle::new(2, &first)).unwrap();
    consensus.compute_monoisotopic_consensus().unwrap();
    close_f32(consensus.intensity, 200.0);
    close(consensus.rt, 1.0);
    close(consensus.mz, 2.0);
    consensus.insert(FeatureHandle::new(4, &second)).unwrap();
    consensus.compute_monoisotopic_consensus().unwrap();
    close_f32(consensus.intensity, 250.0);
    close(consensus.rt, 1.5);
    close(consensus.mz, 2.0);
    consensus.insert(FeatureHandle::new(6, &third)).unwrap();
    consensus.compute_monoisotopic_consensus().unwrap();
    close_f32(consensus.intensity, 300.0);
    close(consensus.rt, 2.0);
    close(consensus.mz, 2.0);
}

/// `void computeDechargeConsensus(const FeatureMap&, bool)` (line 525).
///
/// The C++ section builds its m/z fixtures from the monoisotopic weight of the
/// hydrogen *atom* while `computeDechargeConsensus` subtracts
/// `Constants::PROTON_MASS_U`; its `TEST_REAL_SIMILAR` tolerance absorbs the
/// electron masses. This port uses the proton constant on both sides, so the
/// expected neutral masses are exact rather than approximately right.
#[test]
fn cf_compute_decharge_consensus() {
    let sodium = openms::chemistry::element("Na").unwrap().mono_mass();
    let mass = 1000.0;
    let (add1, add2, add3) = (0.5, 1.0, -0.5);
    let mz1 = (mass + add1 + 3.0 * PROTON_MASS_U) / 3.0;
    let mz2 = (mass + add2 + PROTON_MASS_U + 2.0 * sodium) / 3.0;
    let mz3 = (mass + add3 + 4.0 * PROTON_MASS_U + sodium) / 5.0;
    let mut map = FeatureMap::new();
    let mut consensus = ConsensusFeature::new();

    let mut first = Feature::new(100.0, mz1, 200.0);
    first.charge = 3;
    first.unique_id = 1;
    map.features.push(first.clone());
    consensus.insert(FeatureHandle::new(2, &first)).unwrap();
    consensus.compute_decharge_consensus(&map, false).unwrap();
    close_f32(consensus.intensity, 200.0);
    close(consensus.rt, 100.0);
    close(consensus.mz, mass + add1);

    let mut second = Feature::new(102.0, mz2, 400.0);
    second.charge = 3;
    second.unique_id = 2;
    second.metadata.insert(
        "dc_charge_adduct_mass".into(),
        MetaValue::try_from(2.0 * sodium + PROTON_MASS_U).unwrap(),
    );
    map.features.push(second.clone());
    consensus.insert(FeatureHandle::new(4, &second)).unwrap();
    consensus.compute_decharge_consensus(&map, true).unwrap();
    close_f32(consensus.intensity, 600.0);
    close(consensus.rt, 100.0 / 3.0 + 102.0 * 2.0 / 3.0);
    close(
        consensus.mz,
        (mass + add1) / 3.0 + (mass + add2) * 2.0 / 3.0,
    );
    consensus.compute_decharge_consensus(&map, false).unwrap();
    close_f32(consensus.intensity, 600.0);
    close(consensus.rt, 100.0 / 2.0 + 102.0 / 2.0);
    close(consensus.mz, (mass + add1) / 2.0 + (mass + add2) / 2.0);

    let mut third = Feature::new(101.0, mz3, 600.0);
    third.charge = 5;
    third.unique_id = 3;
    third.metadata.insert(
        "dc_charge_adduct_mass".into(),
        MetaValue::try_from(sodium + 4.0 * PROTON_MASS_U).unwrap(),
    );
    map.features.push(third.clone());
    consensus.insert(FeatureHandle::new(5, &third)).unwrap();
    consensus.compute_decharge_consensus(&map, true).unwrap();
    close_f32(consensus.intensity, 1200.0);
    close(consensus.rt, 100.0 / 6.0 + 102.0 / 3.0 + 101.0 / 2.0);
    close(
        consensus.mz,
        (mass + add1) / 6.0 + (mass + add2) / 3.0 + (mass + add3) / 2.0,
    );
    consensus.compute_decharge_consensus(&map, false).unwrap();
    close_f32(consensus.intensity, 1200.0);
    close(consensus.rt, 100.0 / 3.0 + 102.0 / 3.0 + 101.0 / 3.0);
    close(
        consensus.mz,
        (mass + add1) / 3.0 + (mass + add2) / 3.0 + (mass + add3) / 3.0,
    );
}

/// `Size size() const` (line 600).
#[test]
fn cf_size() {
    let [first_feature, second_feature, third_feature] = consensus_fixture();
    let mut first = ConsensusFeature::from(first_feature.base.clone());
    first.insert(FeatureHandle::new(1, &first_feature)).unwrap();
    first.insert(FeatureHandle::new(2, &third_feature)).unwrap();
    assert_eq!(first.len(), 2);
    let mut second = ConsensusFeature::new();
    assert_eq!(second.len(), 0);
    second
        .insert(FeatureHandle::new(1, &second_feature))
        .unwrap();
    assert_eq!(second.len(), 1);
}

/// `const_iterator begin() const` (line 614).
#[test]
fn cf_begin_const() {
    let [_, second_feature, _] = consensus_fixture();
    let mut consensus = ConsensusFeature::new();
    assert_eq!(consensus.handles().iter().next(), None);
    consensus
        .insert(FeatureHandle::new(1, &second_feature))
        .unwrap();
    assert_eq!(consensus.handles().iter().next().unwrap().unique_id, 5);
}

/// `iterator begin()` (line 626): the port exposes handles as an immutable
/// slice, so there is no mutable iterator; replacement goes through
/// `set_handles`, which re-checks identity uniqueness.
#[test]
fn cf_begin_mut() {
    let [_, second_feature, _] = consensus_fixture();
    let mut consensus = ConsensusFeature::new();
    assert!(consensus.handles().is_empty());
    consensus
        .insert(FeatureHandle::new(1, &second_feature))
        .unwrap();
    let mut handles = consensus.handles().to_vec();
    assert_eq!(handles[0].unique_id, 5);
    handles[0].rt = 42.0;
    consensus.set_handles(handles).unwrap();
    assert_eq!(consensus.handles()[0].rt, 42.0);
}

/// `const_iterator end() const` (line 635): `NOT_TESTABLE`, tested above.
#[test]
fn cf_end_const() {
    assert_eq!(ConsensusFeature::new().handles().iter().next(), None);
}

/// `iterator end()` (line 639): `NOT_TESTABLE`, tested above.
#[test]
fn cf_end_mut() {
    assert!(ConsensusFeature::new().handles().is_empty());
}

/// `const_reverse_iterator rbegin() const` (line 644).
#[test]
fn cf_rbegin_const() {
    let [_, second_feature, _] = consensus_fixture();
    let mut consensus = ConsensusFeature::new();
    assert_eq!(consensus.handles().iter().next_back(), None);
    consensus
        .insert(FeatureHandle::new(1, &second_feature))
        .unwrap();
    assert_eq!(consensus.handles().iter().next_back().unwrap().unique_id, 5);
}

/// `reverse_iterator rbegin()` (line 656).
#[test]
fn cf_rbegin_mut() {
    let [_, second_feature, _] = consensus_fixture();
    let mut consensus = ConsensusFeature::new();
    assert!(consensus.handles().last().is_none());
    consensus
        .insert(FeatureHandle::new(1, &second_feature))
        .unwrap();
    assert_eq!(consensus.handles().last().unwrap().unique_id, 5);
}

/// `const_reverse_iterator rend() const` (line 666): `NOT_TESTABLE`.
#[test]
fn cf_rend_const() {
    assert_eq!(ConsensusFeature::new().handles().iter().next_back(), None);
}

/// `reverse_iterator rend()` (line 670): `NOT_TESTABLE`.
#[test]
fn cf_rend_mut() {
    assert_eq!(ConsensusFeature::new().handles().len(), 0);
}

/// `void clear()` (line 673).
#[test]
fn cf_clear() {
    let [first_feature, second_feature, third_feature] = consensus_fixture();
    let mut first = ConsensusFeature::from(first_feature.base.clone());
    first.insert(FeatureHandle::new(1, &first_feature)).unwrap();
    first.insert(FeatureHandle::new(2, &third_feature)).unwrap();
    first.clear();
    assert_eq!(first.len(), 0);
    let mut second = ConsensusFeature::new();
    assert_eq!(second.len(), 0);
    second
        .insert(FeatureHandle::new(1, &second_feature))
        .unwrap();
    second.clear();
    assert_eq!(second.len(), 0);
}

/// `bool empty() const` (line 689).
#[test]
fn cf_empty() {
    let [first_feature, second_feature, third_feature] = consensus_fixture();
    let mut first = ConsensusFeature::from(first_feature.base.clone());
    first.insert(FeatureHandle::new(1, &first_feature)).unwrap();
    first.insert(FeatureHandle::new(2, &third_feature)).unwrap();
    assert!(!first.is_empty());
    first.clear();
    assert!(first.is_empty());
    let mut second = ConsensusFeature::new();
    assert_eq!(second.len(), 0);
    second
        .insert(FeatureHandle::new(1, &second_feature))
        .unwrap();
    assert!(!second.is_empty());
    second.clear();
    assert!(second.is_empty());
}

// ---------------------------------------------------------------------------
// Native checks for the identification surface. No C++ section covers these:
// the class tests do not exercise the primary ID, the observation matches, the
// reference update or the ratios at all.
// ---------------------------------------------------------------------------

struct Graph {
    graph: IdentificationData,
    peptide: PeptideId,
    first: ObservationMatchId,
    second: ObservationMatchId,
    third: ObservationMatchId,
}

/// Two matches on the same peptide plus one on a compound, so that the
/// same/divergent branches of `getAnnotationState` are both reachable.
fn graph_fixture() -> Graph {
    let mut graph = IdentificationData::new().unwrap();
    let file = graph
        .register_input_file(InputFile::new("test.mzML"))
        .unwrap();
    let software = graph
        .register_processing_software(ProcessingSoftware::new("Tool", "1.0"))
        .unwrap();
    graph
        .register_processing_step(ProcessingStep::new(software), None)
        .unwrap();
    let peptide = graph
        .register_identified_peptide(IdentifiedPeptide::new(
            AASequence::parse("PEPTIDE").unwrap(),
        ))
        .unwrap();
    let compound = graph
        .register_identified_compound(IdentifiedCompound::new("glucose"))
        .unwrap();
    let first_observation = graph
        .register_observation(Observation::new("spectrum_1", file))
        .unwrap();
    let second_observation = graph
        .register_observation(Observation::new("spectrum_2", file))
        .unwrap();
    let first = graph
        .register_observation_match(ObservationMatch::new(peptide, first_observation))
        .unwrap();
    let second = graph
        .register_observation_match(ObservationMatch::new(peptide, second_observation))
        .unwrap();
    let third = graph
        .register_observation_match(ObservationMatch::new(compound, first_observation))
        .unwrap();
    Graph {
        graph,
        peptide,
        first,
        second,
        third,
    }
}

#[test]
fn annotation_state_from_matches_ignores_legacy_identifications() {
    let fixture = graph_fixture();
    let mut feature = base_feature();
    feature.peptide_identifications = vec![identification(vec![hit("ABCDE", 0.5)])];
    assert_eq!(
        feature.annotation_state(Some(&fixture.graph)).unwrap(),
        AnnotationState::Single
    );
    feature.add_id_match(fixture.first).unwrap();
    assert_eq!(
        feature.annotation_state(Some(&fixture.graph)).unwrap(),
        AnnotationState::Single
    );
    // Without a graph the matches cannot be resolved at all.
    assert!(matches!(
        feature.annotation_state(None),
        Err(openms::Error::MissingInformation(_))
    ));
    feature.add_id_match(fixture.second).unwrap();
    assert_eq!(
        feature.annotation_state(Some(&fixture.graph)).unwrap(),
        AnnotationState::MultipleSame
    );
    feature.add_id_match(fixture.third).unwrap();
    assert_eq!(
        feature.annotation_state(Some(&fixture.graph)).unwrap(),
        AnnotationState::MultipleDivergent
    );
    // A match from another graph is rejected instead of silently resolving.
    let other = graph_fixture();
    let mut foreign = base_feature();
    foreign.add_id_match(other.first).unwrap();
    assert!(foreign.annotation_state(Some(&fixture.graph)).is_err());
}

#[test]
fn annotation_state_names_follow_the_source_array() {
    assert_eq!(AnnotationState::NAMES.len(), 4);
    assert_eq!(AnnotationState::None.name(), "no ID");
    assert_eq!(AnnotationState::Single.name(), "single ID");
    assert_eq!(
        AnnotationState::MultipleSame.name(),
        "multiple IDs (identical)"
    );
    assert_eq!(
        AnnotationState::MultipleDivergent.to_string(),
        "multiple IDs (divergent)"
    );
    assert_eq!(AnnotationState::default(), AnnotationState::None);
}

#[test]
fn annotation_state_counts_only_identifications_that_have_hits() {
    let mut feature = base_feature();
    feature.peptide_identifications = vec![
        identification(vec![hit("ABCDE", 0.5)]),
        identification(Vec::new()),
    ];
    // The source collects one sequence from two identifications and therefore
    // reports MULTIPLE_SAME rather than SINGLE; the port keeps that.
    assert_eq!(
        feature.annotation_state(None).unwrap(),
        AnnotationState::MultipleSame
    );
    feature.peptide_identifications = vec![identification(Vec::new()), identification(Vec::new())];
    assert_eq!(
        feature.annotation_state(None).unwrap(),
        AnnotationState::None
    );
    // The best hit decides, not the stored hit order.
    feature.peptide_identifications = vec![
        identification(vec![hit("ABCDE", 0.1), hit("KRGH", 0.9)]),
        identification(vec![hit("KRGH", 0.9)]),
    ];
    assert_eq!(
        feature.annotation_state(None).unwrap(),
        AnnotationState::MultipleSame
    );
    let mut lower_better = feature.clone();
    for identification in &mut lower_better.peptide_identifications {
        identification.higher_score_better = false;
    }
    assert_eq!(
        lower_better.annotation_state(None).unwrap(),
        AnnotationState::MultipleDivergent
    );
    // A non-finite score is reported instead of being compared.
    feature.peptide_identifications[0].hits[0].score = f64::NAN;
    assert!(feature.annotation_state(None).is_err());
}

#[test]
fn sort_peptide_identifications_is_checked_and_atomic() {
    let mut feature = base_feature();
    feature.peptide_identifications = vec![
        identification(vec![hit("ABCDE", 0.8)]),
        identification(vec![hit("KRGH", 0.9)]),
    ];
    // A hit constructor rejects a non-finite score, so the invalid state is
    // reachable only through the public field, as for the source's setters.
    feature.peptide_identifications[1].hits[0].score = f64::INFINITY;
    let before = feature.clone();
    assert!(feature.sort_peptide_identifications().is_err());
    assert_eq!(feature, before);
    feature.peptide_identifications[1].hits[0].score = 0.9;
    feature.peptide_identifications[1].higher_score_better = false;
    let before = feature.clone();
    assert!(feature.sort_peptide_identifications().is_err());
    assert_eq!(feature, before);
    feature.peptide_identifications[1].higher_score_better = true;
    feature.sort_peptide_identifications().unwrap();
    assert_eq!(feature.peptide_identifications[0].hits[0].score, 0.9);
    // Lower-is-better inverts the order of the identifications.
    for identification in &mut feature.peptide_identifications {
        identification.higher_score_better = false;
    }
    feature.sort_peptide_identifications().unwrap();
    assert_eq!(feature.peptide_identifications[0].hits[0].score, 0.8);
    // Equal best scores keep their relative order, which the source does not
    // guarantee because `std::sort` is not stable.
    feature.peptide_identifications = vec![
        identification(vec![hit("ABCDE", 0.5)]),
        identification(vec![hit("KRGH", 0.5)]),
    ];
    feature.sort_peptide_identifications().unwrap();
    assert_eq!(
        feature.peptide_identifications[0].hits[0]
            .sequence
            .to_string(),
        "ABCDE"
    );
    // Several empty identifications are fine here; the source comparator
    // answers "less" for that pair and breaks its strict weak ordering.
    feature.peptide_identifications = vec![
        identification(Vec::new()),
        identification(vec![hit("ABCDE", 0.5)]),
        identification(Vec::new()),
    ];
    feature.sort_peptide_identifications().unwrap();
    assert_eq!(feature.peptide_identifications[0].hits.len(), 1);
    assert!(feature.peptide_identifications[1].hits.is_empty());
    assert!(feature.peptide_identifications[2].hits.is_empty());
}

#[test]
fn primary_id_assignment_clearing_and_checking() {
    let fixture = graph_fixture();
    let mut feature = base_feature();
    assert!(!feature.has_primary_id());
    assert!(matches!(
        feature.primary_id(),
        Err(openms::Error::MissingInformation(_))
    ));
    feature.set_primary_id(fixture.peptide);
    assert!(feature.has_primary_id());
    assert_eq!(
        feature.primary_id().unwrap().peptide().unwrap(),
        fixture.peptide
    );
    feature.clear_primary_id();
    assert!(!feature.has_primary_id());
    assert!(feature.id_matches.is_empty());
    feature
        .set_primary_id_checked(&fixture.graph, fixture.peptide)
        .unwrap();
    assert!(feature.has_primary_id());
    // A reference from a different graph generation is refused.
    let other = graph_fixture();
    let mut foreign = base_feature();
    assert!(
        foreign
            .set_primary_id_checked(&fixture.graph, other.peptide)
            .is_err()
    );
    assert!(!foreign.has_primary_id());
    assert!(
        foreign
            .add_id_match_checked(&fixture.graph, other.first)
            .is_err()
    );
    assert!(foreign.id_matches.is_empty());
    foreign
        .add_id_match_checked(&fixture.graph, fixture.first)
        .unwrap();
    assert_eq!(foreign.id_matches.len(), 1);
    // Re-adding an existing match is a no-op, as for the source std::set.
    foreign.add_id_match(fixture.first).unwrap();
    assert_eq!(foreign.id_matches.len(), 1);
}

#[test]
fn update_id_references_translates_atomically() {
    let source = graph_fixture();
    let mut destination = IdentificationData::new().unwrap();
    let translator = destination.merge_from(&source.graph).unwrap();
    let mut feature = base_feature();
    feature.set_primary_id(source.peptide);
    feature.add_id_match(source.first).unwrap();
    feature.add_id_match(source.third).unwrap();
    feature.update_id_references(&translator).unwrap();
    let translated = feature.primary_id().unwrap();
    assert_eq!(
        translated,
        translator.molecule(source.peptide.into()).unwrap()
    );
    // The translated references resolve in the destination graph, the old ones
    // do not, and the annotation state is computable again.
    assert_eq!(
        feature.annotation_state(Some(&destination)).unwrap(),
        AnnotationState::MultipleDivergent
    );
    assert!(feature.annotation_state(Some(&source.graph)).is_err());
    // A second update has no translation for the already-translated IDs and
    // must leave the feature exactly as it was.
    let before = feature.clone();
    assert!(feature.update_id_references(&translator).is_err());
    assert_eq!(feature, before);
}

#[test]
fn update_all_id_references_covers_subordinates_or_changes_nothing() {
    let source = graph_fixture();
    let mut destination = IdentificationData::new().unwrap();
    let translator = destination.merge_from(&source.graph).unwrap();
    let mut feature = plain_feature();
    feature.set_primary_id(source.peptide);
    let mut child = plain_feature();
    child.add_id_match(source.second).unwrap();
    let mut grandchild = plain_feature();
    grandchild.set_primary_id(source.peptide);
    grandchild.add_id_match(source.first).unwrap();
    child.subordinates.push(grandchild);
    feature.subordinates.push(child);
    feature.update_all_id_references(&translator).unwrap();
    assert_eq!(
        feature.subordinates[0].subordinates[0]
            .primary_id()
            .unwrap(),
        translator.molecule(source.peptide.into()).unwrap()
    );
    assert_eq!(
        feature.subordinates[0]
            .annotation_state(Some(&destination))
            .unwrap(),
        AnnotationState::Single
    );
    // One untranslatable reference deep in the tree leaves everything alone,
    // unlike the source, which updates on the way down.
    let other = graph_fixture();
    let mut mixed = feature.clone();
    mixed.subordinates[0].subordinates[0]
        .add_id_match(other.first)
        .unwrap();
    let before = mixed.clone();
    assert!(mixed.update_all_id_references(&translator).is_err());
    assert_eq!(mixed, before);
}

#[test]
fn unique_id_traversal_is_pre_order_and_bounded() {
    let mut feature = plain_feature();
    feature.unique_id = 1;
    let mut child = plain_feature();
    child.unique_id = 2;
    child.subordinates.push(Feature::default());
    child.subordinates[0].unique_id = 3;
    feature.subordinates.push(child);
    let mut last = plain_feature();
    last.unique_id = 4;
    feature.subordinates.push(last);
    let mut order = Vec::new();
    let total = feature
        .for_each_unique_id(|id| {
            order.push(*id);
            usize::try_from(*id).unwrap()
        })
        .unwrap();
    assert_eq!(order, vec![1, 2, 3, 4]);
    assert_eq!(total, 10);
    let mut seen = Vec::new();
    feature
        .count_unique_ids(|id| {
            seen.push(id);
            0
        })
        .unwrap();
    assert_eq!(seen, order);
    // The mutable traversal can assign, matching `ensureUniqueId` use.
    feature
        .for_each_unique_id(|id| id.clear_unique_id())
        .unwrap();
    assert_eq!(
        feature.count_unique_ids(|id| usize::from(id == 0)).unwrap(),
        4
    );
    // Deeper than the checked subordinate depth is refused.
    let mut deep = plain_feature();
    for _ in 0..=Feature::MAX_SUBORDINATE_DEPTH {
        deep = Feature {
            subordinates: vec![deep],
            ..Feature::default()
        };
    }
    assert!(deep.count_unique_ids(|_| 0).is_err());
    assert!(deep.for_each_unique_id(|_| 0).is_err());
    // An accumulator overflow is an error, not a wrap.
    let mut pair = plain_feature();
    pair.subordinates.push(Feature::default());
    assert!(pair.for_each_unique_id(|_| usize::MAX).is_err());
}

#[test]
fn ratios_are_checked_and_replaced_atomically() {
    let mut consensus = ConsensusFeature::new();
    assert!(consensus.ratios().is_empty());
    let ratio = Ratio {
        ratio_value: 1.5,
        denominator_ref: "channel_0".into(),
        numerator_ref: "channel_1".into(),
        description: Vec::new(),
    };
    consensus.add_ratio(ratio.clone()).unwrap();
    consensus.add_ratio(ratio.clone()).unwrap();
    assert_eq!(consensus.ratios().len(), 2);
    assert_eq!(consensus.ratios()[1].numerator_ref, "channel_1");
    assert_eq!(consensus.ratios(), consensus.ratios.as_slice());
    let before = consensus.clone();
    assert!(
        consensus
            .add_ratio(Ratio {
                ratio_value: f64::NAN,
                ..Ratio::default()
            })
            .is_err()
    );
    assert_eq!(consensus, before);
    assert!(
        consensus
            .set_ratios(vec![
                ratio.clone(),
                Ratio {
                    ratio_value: f64::INFINITY,
                    ..Ratio::default()
                },
            ])
            .is_err()
    );
    assert_eq!(consensus, before);
    consensus.set_ratios(vec![ratio]).unwrap();
    assert_eq!(consensus.ratios().len(), 1);
    // The default ratio is a defined zero; the source leaves it uninitialised.
    assert_eq!(Ratio::default().ratio_value, 0.0);
    Ratio::default().validate().unwrap();
}

#[test]
fn clone_with_map_index_stamps_every_attached_identification() {
    let mut feature = BaseFeature::new(1.0, 2.0, 3.0);
    feature.peptide_identifications = vec![
        identification(vec![hit("ABCDE", 0.5)]),
        identification(Vec::new()),
    ];
    feature.metadata.insert("label".into(), "a".into());
    let stamped = feature.clone_with_map_index(7).unwrap();
    for identification in &stamped.peptide_identifications {
        assert_eq!(identification.metadata["map_index"].as_i64().unwrap(), 7);
    }
    assert_eq!(stamped.metadata["label"].as_str().unwrap(), "a");
    assert_eq!(stamped.rt, 1.0);
    // The original is untouched, and an invalid attached record is refused.
    assert!(
        !feature.peptide_identifications[0]
            .metadata
            .contains_key("map_index")
    );
    feature.peptide_identifications[0].hits[0].score = f64::NAN;
    assert!(feature.clone_with_map_index(7).is_err());
}
