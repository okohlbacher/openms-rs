// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `FeatureMap_test.cpp` (32 sections) and `ConsensusMap_test.cpp`
//! (39 sections) at Core SDK `bc9cc12514c768385ce121d6ca4bb710fe1983c4`, plus
//! native checks for the container operations added in
//! `src/kernel/map_operations.rs`.
//!
//! Every `START_SECTION` of both class tests is covered here;
//! `docs/MAP_OPERATIONS_SUPPORT.md` holds the section-to-test table. All
//! transcribed expectations are tier 3 (source review): the literals come from
//! the pinned class tests, no C++ was built or executed.

use openms::chemistry::AASequence;
use openms::concept::{HasUniqueId, UniqueIdGenerator};
use openms::format::FileType;
use openms::identification::{PeptideHit, PeptideIdentification, ProteinIdentification};
use openms::kernel::feature_identification::AnnotationState;
use openms::kernel::features::{
    ColumnHeader, ConsensusFeature, ConsensusMap, Feature, FeatureHandle, FeatureMap,
};
use openms::kernel::geometry::{ConvexHull2D, Point2D};
use openms::kernel::map_operations::{AnnotationStatistics, SplitMeta};
use openms::kernel::{MSExperiment, MSSpectrum, Peak1D, Peak2D};
use openms::metadata::{DataProcessing, MetaValue, Software};

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a} != {b}");
}
fn generator() -> UniqueIdGenerator {
    UniqueIdGenerator::from_seed(20_260_912)
}
fn hit(sequence: &str) -> PeptideHit {
    PeptideHit::new(0.0, 0, 0, AASequence::parse(sequence).unwrap()).unwrap()
}
fn identification(sequence: &str) -> PeptideIdentification {
    PeptideIdentification {
        hits: vec![hit(sequence)],
        ..PeptideIdentification::default()
    }
}

/// `feature1` .. `feature4` of `FeatureMap_test.cpp` (lines 52-88).
///
/// `feature1` carries one identification with one hit, `feature2` two with the
/// same hit, `feature3` two with different hits and `feature4` none but a
/// three-point convex hull.
fn source_features() -> [Feature; 4] {
    let mut feature1 = Feature::new(2.0, 3.0, 1.0);
    feature1.base.peptide_identifications = vec![identification("ABCDE")];

    let mut feature2 = Feature::new(0.0, 2.5, 0.5);
    feature2.base.peptide_identifications = vec![identification("ABCDE"), identification("ABCDE")];

    let mut feature3 = Feature::new(10.5, 0.0, 0.01);
    feature3.base.peptide_identifications = vec![identification("ABCDE"), identification("KRGH")];

    let mut feature4 = Feature::new(5.25, 1.5, 0.5);
    feature4.convex_hulls = vec![
        ConvexHull2D::from_points(&[
            Point2D::new(-1.0, 2.0),
            Point2D::new(4.0, 1.2),
            Point2D::new(5.0, 3.123),
        ])
        .unwrap(),
    ];

    [feature1, feature2, feature3, feature4]
}

/// The `map_const_1` fixture of `ConsensusMap_test.cpp` (lines 309-323).
fn map_const_1() -> ConsensusMap {
    let mut map = ConsensusMap::from_features(vec![ConsensusFeature::new(); 3]);
    map.metadata.insert("meta".into(), "value".into());
    map.identifier = "lsid".into();
    let header = map.column_headers.entry(0).or_default();
    header.filename = "blub".into();
    header.size = 47;
    header.label = "label".into();
    header.metadata.insert("meta".into(), "meta".into());
    map.data_processing.resize(1, DataProcessing::default());
    map.set_experiment_type("labeled_MS2").unwrap();
    map.protein_identifications
        .resize(1, ProteinIdentification::default());
    map.unassigned_peptide_identifications
        .resize(1, PeptideIdentification::default());
    map
}

/// The four `feature1` .. `feature4` values of `ConsensusMap_test.cpp`
/// (lines 97-115), which unlike the FeatureMap fixture carry no identifications.
fn consensus_source_features() -> [Feature; 4] {
    [
        Feature::new(2.0, 3.0, 1.0),
        Feature::new(0.0, 2.5, 0.5),
        Feature::new(10.5, 0.0, 0.01),
        Feature::new(5.25, 1.5, 0.5),
    ]
}

// ===========================================================================
// FeatureMap_test.cpp
// ===========================================================================

/// `FeatureMap()` (line 40) and `virtual ~FeatureMap()` (line 48).
///
/// The source allocates on the heap and checks the pointer, then deletes it in
/// the destructor section. Rust owns the value; the observable content is the
/// empty map with no ranges at all.
#[test]
fn fm_default_constructor_and_destructor() {
    let map = FeatureMap::new();
    assert_eq!(map.features.len(), 0);
    let ranges = map.ranges().unwrap();
    assert!(ranges.rt.is_none() && ranges.mz.is_none() && ranges.intensity.is_none());
    drop(map);
}

/// `const std::vector<ProteinIdentification>& getProteinIdentifications() const`
/// (line 90), the non-const overload (line 95) and
/// `setProteinIdentifications` (line 101).
#[test]
fn fm_protein_identification_accessors() {
    let mut map = FeatureMap::new();
    assert_eq!(map.protein_identifications.len(), 0);
    map.protein_identifications
        .resize(1, ProteinIdentification::default());
    assert_eq!(map.protein_identifications.len(), 1);
    let mut map = FeatureMap::new();
    map.protein_identifications = vec![ProteinIdentification::default(); 2];
    assert_eq!(map.protein_identifications.len(), 2);
}

/// `const PeptideIdentificationList& getUnassignedPeptideIdentifications() const`
/// (line 107), the non-const overload (line 112) and the setter (line 118).
#[test]
fn fm_unassigned_peptide_identification_accessors() {
    let mut map = FeatureMap::new();
    assert_eq!(map.unassigned_peptide_identifications.len(), 0);
    map.unassigned_peptide_identifications
        .resize(1, PeptideIdentification::default());
    assert_eq!(map.unassigned_peptide_identifications.len(), 1);
    let mut map = FeatureMap::new();
    map.unassigned_peptide_identifications = vec![PeptideIdentification::default(); 2];
    assert_eq!(map.unassigned_peptide_identifications.len(), 2);
}

/// `const std::vector<DataProcessing>& getDataProcessing() const` (line 124),
/// the non-const overload (line 129) and `setDataProcessing` (line 135).
#[test]
fn fm_data_processing_accessors() {
    let mut map = FeatureMap::new();
    assert_eq!(map.data_processing.len(), 0);
    map.data_processing.resize(1, DataProcessing::default());
    assert_eq!(map.data_processing.len(), 1);
    let mut map = FeatureMap::new();
    map.data_processing = vec![DataProcessing::default()];
    assert_eq!(map.data_processing.len(), 1);
}

/// `void updateRanges()` (line 143).
///
/// `FeatureMap::ranges` recomputes on demand, so calling it twice is the
/// source's "second time to check the initialization".
#[test]
fn fm_update_ranges() {
    let [feature1, feature2, feature3, feature4] = source_features();
    let mut map =
        FeatureMap::from_features(vec![feature1.clone(), feature2.clone(), feature3.clone()]);
    map.ranges().unwrap();
    let ranges = map.ranges().unwrap();
    close(ranges.intensity.unwrap().max, 1.0);
    close(ranges.intensity.unwrap().min, f64::from(0.01f32));
    close(ranges.rt.unwrap().max, 10.5);
    close(ranges.mz.unwrap().max, 3.0);
    close(ranges.rt.unwrap().min, 0.0);
    close(ranges.mz.unwrap().min, 0.0);

    map.features.push(feature4);
    let ranges = map.ranges().unwrap();
    close(ranges.intensity.unwrap().max, 1.0);
    close(ranges.intensity.unwrap().min, f64::from(0.01f32));
    close(ranges.rt.unwrap().max, 10.5);
    close(ranges.mz.unwrap().max, 3.123);
    close(ranges.rt.unwrap().min, -1.0);
    close(ranges.mz.unwrap().min, 0.0);
}

/// `FeatureMap(const FeatureMap& source)` (line 172), `operator=(const
/// FeatureMap&)` (line 195) and `operator=(FeatureMap&&)` (line 230).
///
/// The three source sections assert the same content; Rust's `Clone` covers the
/// first two and a move the third. The source's copy constructor additionally
/// merges the embedded `IdentificationData` and repoints the features; here the
/// graph is the caller's, so a clone carries the owner-tagged IDs unchanged.
#[test]
fn fm_copy_assign_and_move() {
    let [feature1, feature2, feature3, _] = source_features();
    let mut map1 = FeatureMap::from_features(vec![feature1, feature2, feature3]);
    map1.metadata.insert("meta".into(), "value".into());
    map1.identifier = "lsid".into();
    map1.data_processing.resize(1, DataProcessing::default());
    map1.protein_identifications
        .resize(1, ProteinIdentification::default());
    map1.unassigned_peptide_identifications
        .resize(1, PeptideIdentification::default());

    let check = |map: &FeatureMap| {
        assert_eq!(map.features.len(), 3);
        assert_eq!(map.metadata["meta"].as_str().unwrap(), "value");
        close(map.ranges().unwrap().intensity.unwrap().max, 1.0);
        assert_eq!(map.identifier, "lsid");
        assert_eq!(map.data_processing.len(), 1);
        assert_eq!(map.protein_identifications.len(), 1);
        assert_eq!(map.unassigned_peptide_identifications.len(), 1);
    };

    let map2 = map1.clone();
    check(&map2);
    let mut map2 = FeatureMap::new();
    assert!(map2.is_empty());
    map2 = map1.clone();
    check(&map2);

    // Assignment of an empty object.
    map2 = FeatureMap::new();
    assert_eq!(map2.features.len(), 0);
    assert!(map2.ranges().unwrap().rt.is_none());
    assert_eq!(map2.identifier, "");
    assert_eq!(map2.data_processing.len(), 0);
    assert_eq!(map2.protein_identifications.len(), 0);
    assert_eq!(map2.unassigned_peptide_identifications.len(), 0);

    let moved = map1;
    check(&moved);
}

/// `bool operator==(const FeatureMap&) const` (line 265) and
/// `bool operator!=(const FeatureMap&) const` (line 298).
///
/// The two source sections edit the same seven members. The final case is the
/// one place where this port differs: the source compares the cached range
/// information as well, and `clear(false)` deliberately leaves that cache
/// populated, so its two maps differ. Ranges are computed on demand here, so
/// there is no cache to differ and the maps compare equal.
#[test]
fn fm_equality_and_inequality() {
    let [feature1, feature2, _, _] = source_features();
    let empty = FeatureMap::new();

    let mut edit = FeatureMap::new();
    assert!(empty == edit);
    assert!(!(empty != edit));

    edit.identifier = "lsid".into();
    assert!(empty != edit);

    edit = empty.clone();
    edit.features.push(feature1.clone());
    assert!(empty != edit);

    edit = empty.clone();
    edit.data_processing.resize(1, DataProcessing::default());
    assert!(empty != edit);

    edit = empty.clone();
    edit.protein_identifications
        .resize(1, ProteinIdentification::default());
    assert!(edit != empty);
    edit = empty.clone();
    edit.protein_identifications
        .resize(10, ProteinIdentification::default());
    assert!(edit != empty);

    edit = empty.clone();
    edit.unassigned_peptide_identifications
        .resize(10, PeptideIdentification::default());
    assert!(empty != edit);

    edit = empty.clone();
    edit.features.push(feature1);
    edit.features.push(feature2);
    edit.ranges().unwrap();
    edit.clear(false);
    // Source expectation: `empty == edit` is false, because `clear(false)` keeps
    // the cached RangeManager the source's operator== compares.
    assert_eq!(empty, edit);
}

/// `FeatureMap operator+(const FeatureMap&) const` (line 331).
#[test]
fn fm_operator_plus() {
    let mut generator = generator();
    let mut m1 = FeatureMap::new();
    let m2 = FeatureMap::new();
    let mut m3 = FeatureMap::new();
    assert_eq!(m1.merged(&m2, &mut generator).unwrap(), m3);

    let mut f1 = Feature::default();
    f1.base.mz = 100.12;
    m1.features.push(f1);
    m3 = m1.clone();
    assert_eq!(m1.merged(&m2, &mut generator).unwrap(), m3);
}

/// `FeatureMap& operator+=(const FeatureMap&)` (line 344).
#[test]
fn fm_operator_plus_assign() {
    let mut generator = generator();
    let mut m1 = FeatureMap::new();
    let mut m2 = FeatureMap::new();
    let mut m3 = FeatureMap::new();

    // Adding empty maps has no effect.
    m1.append(&m2, &mut generator).unwrap();
    assert_eq!(m1, m3);

    let mut f1 = Feature::default();
    f1.base.mz = 100.12;
    m1.features.push(f1);
    m3 = m1.clone();
    m1.append(&m2, &mut generator).unwrap();
    assert_eq!(m1, m3);

    m1.identifier = "123".into();
    m1.data_processing.resize(1, DataProcessing::default());
    m1.protein_identifications
        .resize(1, ProteinIdentification::default());
    m1.unassigned_peptide_identifications
        .resize(1, PeptideIdentification::default());
    m1.unique_id.ensure_unique_id(&mut generator);

    m2.identifier = "321".into();
    m2.data_processing.resize(2, DataProcessing::default());
    m2.protein_identifications
        .resize(2, ProteinIdentification::default());
    m2.unassigned_peptide_identifications
        .resize(2, PeptideIdentification::default());
    m2.features.push(Feature::default());
    m2.features.push(Feature::default());

    let replaced = m1.append(&m2, &mut generator).unwrap();
    assert_eq!(m1.identifier, "");
    assert!(!m1.unique_id.has_valid_unique_id());
    assert_eq!(m1.data_processing.len(), 3);
    assert_eq!(m1.protein_identifications.len(), 3);
    assert_eq!(m1.unassigned_peptide_identifications.len(), 3);
    assert_eq!(m1.features.len(), 3);
    // Every feature is unassigned, so no unique ID collides and none is redrawn.
    assert_eq!(replaced, 0);
}

/// Native: `append` redraws colliding unique IDs and reports how many, and a
/// failure leaves the destination untouched.
#[test]
fn fm_append_resolves_unique_id_conflicts_atomically() {
    let mut generator = generator();
    let mut left = FeatureMap::from_features(vec![Feature::default()]);
    left.features[0].base.unique_id = 7;
    let mut right = FeatureMap::from_features(vec![Feature::default()]);
    right.features[0].base.unique_id = 7;

    let before = left.clone();
    let replaced = left.append(&right, &mut generator).unwrap();
    assert_eq!(replaced, 1);
    assert_eq!(left.features.len(), 2);
    assert_ne!(left.features[0].unique_id, left.features[1].unique_id);
    // Both IDs remain assigned, and the map now validates.
    left.validate().unwrap();
    assert_ne!(left, before);
}

/// `void sortByIntensity(bool reverse=false)` (line 383).
#[test]
fn fm_sort_by_intensity() {
    let mut map = FeatureMap::from_features(vec![
        Feature::new(0.0, 0.0, 10.0),
        Feature::new(0.0, 0.0, 5.0),
        Feature::new(0.0, 0.0, 3.0),
    ]);
    map.sort_by_intensity(false).unwrap();
    assert_eq!(map.features[0].intensity, 3.0);
    assert_eq!(map.features[1].intensity, 5.0);
    assert_eq!(map.features[2].intensity, 10.0);

    map.sort_by_intensity(true).unwrap();
    assert_eq!(map.features[0].intensity, 10.0);
    assert_eq!(map.features[1].intensity, 5.0);
    assert_eq!(map.features[2].intensity, 3.0);
}

/// `void sortByPosition()` (line 411).
#[test]
fn fm_sort_by_position() {
    let mut map = FeatureMap::from_features(vec![
        Feature::new(10.0, 0.0, 0.0),
        Feature::new(5.0, 0.0, 0.0),
        Feature::new(3.0, 0.0, 0.0),
    ]);
    map.sort_by_position().unwrap();
    assert_eq!(map.features[0].rt, 3.0);
    assert_eq!(map.features[1].rt, 5.0);
    assert_eq!(map.features[2].rt, 10.0);
}

/// `void sortByMZ()` (line 433).
#[test]
fn fm_sort_by_mz() {
    let mut map = FeatureMap::from_features(vec![
        Feature::new(10.0, 25.0, 0.0),
        Feature::new(5.0, 15.0, 0.0),
        Feature::new(3.0, 10.0, 0.0),
    ]);
    map.sort_by_mz().unwrap();
    assert_eq!(map.features[0].mz, 10.0);
    assert_eq!(map.features[1].mz, 15.0);
    assert_eq!(map.features[2].mz, 25.0);
}

/// `void sortByRT()` (line 458).
#[test]
fn fm_sort_by_rt() {
    let mut map = FeatureMap::from_features(vec![
        Feature::new(10.0, 25.0, 0.0),
        Feature::new(5.0, 15.0, 0.0),
        Feature::new(3.0, 10.0, 0.0),
    ]);
    map.sort_by_rt().unwrap();
    assert_eq!(map.features[0].rt, 3.0);
    assert_eq!(map.features[1].rt, 5.0);
    assert_eq!(map.features[2].rt, 10.0);
}

/// `void swap(FeatureMap& from)` (line 483).
#[test]
fn fm_swap() {
    let [feature1, feature2, _, _] = source_features();
    let mut map1 = FeatureMap::from_features(vec![feature1, feature2]);
    map1.identifier = "stupid comment".into();
    map1.data_processing.resize(1, DataProcessing::default());
    map1.protein_identifications
        .resize(1, ProteinIdentification::default());
    map1.unassigned_peptide_identifications
        .resize(1, PeptideIdentification::default());
    // Native addition: the source does not swap the MetaInfoInterface.
    map1.metadata.insert("keep".into(), "left".into());
    let mut map2 = FeatureMap::new();

    map1.swap(&mut map2);

    assert_eq!(map1.identifier, "");
    assert_eq!(map1.features.len(), 0);
    assert!(map1.ranges().unwrap().rt.is_none());
    assert_eq!(map1.data_processing.len(), 0);
    assert_eq!(map1.protein_identifications.len(), 0);
    assert_eq!(map1.unassigned_peptide_identifications.len(), 0);

    assert_eq!(map2.identifier, "stupid comment");
    assert_eq!(map2.features.len(), 2);
    close(map2.ranges().unwrap().intensity.unwrap().min, 0.5);
    assert_eq!(map2.data_processing.len(), 1);
    assert_eq!(map2.protein_identifications.len(), 1);
    assert_eq!(map2.unassigned_peptide_identifications.len(), 1);

    // Source quirk: meta values stay with their map.
    assert_eq!(map1.metadata["keep"].as_str().unwrap(), "left");
    assert!(!map2.metadata.contains_key("keep"));
}

/// `void swapFeaturesOnly(FeatureMap& from)` (line 512).
#[test]
fn fm_swap_features_only() {
    let [feature1, feature2, _, _] = source_features();
    let mut map1 = FeatureMap::from_features(vec![feature1, feature2]);
    map1.identifier = "stupid comment".into();
    map1.data_processing.resize(1, DataProcessing::default());
    map1.protein_identifications
        .resize(1, ProteinIdentification::default());
    map1.unassigned_peptide_identifications
        .resize(1, PeptideIdentification::default());
    let mut map2 = FeatureMap::new();

    map1.swap_features_only(&mut map2);

    assert_eq!(map1.identifier, "stupid comment");
    assert_eq!(map1.features.len(), 0);
    assert!(map1.ranges().unwrap().rt.is_none());
    assert_eq!(map1.data_processing.len(), 1);
    assert_eq!(map1.protein_identifications.len(), 1);
    assert_eq!(map1.unassigned_peptide_identifications.len(), 1);

    assert_eq!(map2.identifier, "");
    assert_eq!(map2.features.len(), 2);
    close(map2.ranges().unwrap().intensity.unwrap().min, 0.5);
    assert_eq!(map2.data_processing.len(), 0);
    assert_eq!(map2.protein_identifications.len(), 0);
    assert_eq!(map2.unassigned_peptide_identifications.len(), 0);
}

/// `void sortByOverallQuality(bool reverse=false)` (line 541).
#[test]
fn fm_sort_by_overall_quality() {
    let mut map = FeatureMap::new();
    for (rt, quality) in [(1.0, 10.0f32), (2.0, 30.0), (3.0, 20.0)] {
        let mut feature = Feature::new(rt, rt, 0.0);
        feature.base.quality = quality;
        map.features.push(feature);
    }
    map.sort_by_quality(false).unwrap();
    assert_eq!(map.features[0].rt, 1.0);
    assert_eq!(map.features[1].rt, 3.0);
    assert_eq!(map.features[2].rt, 2.0);
    assert_eq!(map.features[0].quality, 10.0);
    assert_eq!(map.features[1].quality, 20.0);
    assert_eq!(map.features[2].quality, 30.0);

    map.sort_by_quality(true).unwrap();
    assert_eq!(map.features[0].rt, 2.0);
    assert_eq!(map.features[1].rt, 3.0);
    assert_eq!(map.features[2].rt, 1.0);
    assert_eq!(map.features[0].quality, 30.0);
    assert_eq!(map.features[1].quality, 20.0);
    assert_eq!(map.features[2].quality, 10.0);
}

/// `void clear(bool clear_meta_data=true)` (line 583).
#[test]
fn fm_clear() {
    let [feature1, feature2, _, _] = source_features();
    let mut map = FeatureMap::from_features(vec![feature1, feature2]);
    map.identifier = "stupid comment".into();
    map.data_processing.resize(1, DataProcessing::default());
    map.protein_identifications
        .resize(1, ProteinIdentification::default());
    map.unassigned_peptide_identifications
        .resize(1, PeptideIdentification::default());
    map.ranges().unwrap();

    map.clear(false);
    assert_eq!(map.features.len(), 0);
    assert!(map != FeatureMap::new());
    assert!(map.is_empty());

    map.clear(true);
    assert_eq!(map.features.len(), 0);
    assert_eq!(map, FeatureMap::new());
    assert!(map.is_empty());
}

/// `[EXTRA] void uniqueIdToIndex()` (line 605).
///
/// The source shuffles with `Math::RandomShuffler{0}`; the shuffle itself is
/// not the subject, so a fixed rotation is used here. The final part of the
/// section builds a container with seven elements, four of which carry a valid
/// unique ID but only three distinct ones, and expects
/// `updateUniqueIdToIndex()` to throw `Exception::Postcondition`; the port's
/// index lookup reports that duplication as an error instead.
#[test]
fn fm_unique_id_to_index() {
    let mut generator = generator();
    let mut map = FeatureMap::new();
    let mut pairs: Vec<(usize, u64)> = Vec::new();
    for i in 0..4usize {
        let mut feature = Feature::new(i as f64 * 100.0, 23.9, 0.0);
        feature.base.unique_id.assign_new_unique_id(&mut generator);
        pairs.push((i, feature.unique_id));
        map.features.push(feature);
    }
    for &(index, id) in &pairs {
        assert_eq!(map.unique_id_to_index(id).unwrap(), Some(index));
    }

    // "shuffling ..." - a deterministic rotation stands in for the source's
    // portable_random_shuffle; the invariant checked is the same.
    map.features.rotate_left(1);
    pairs.rotate_left(1);
    for (i, pair) in pairs.iter().enumerate() {
        let id = map.features[i].unique_id;
        assert_eq!(map.unique_id_to_index(id).unwrap(), Some(i));
        assert_eq!(
            map.features[map.unique_id_to_index(pair.1).unwrap().unwrap()].unique_id,
            pair.1
        );
    }

    let mut extra = Feature::new(98_765_421.0, 23.9, 0.0);
    extra.base.unique_id.assign_new_unique_id(&mut generator);
    // Source returns Size(-1) for an unknown id; the port returns None.
    assert_eq!(map.unique_id_to_index(extra.unique_id).unwrap(), None);
    map.features.push(extra.clone());
    assert_eq!(
        map.unique_id_to_index(extra.unique_id).unwrap(),
        Some(map.features.len() - 1)
    );

    map.features.push(Feature::default());
    map.features.push(extra);
    map.features.push(Feature::default());
    map.features.push(Feature::default());
    assert_eq!(map.features.len(), 9);
    map.features.remove(1);
    map.features.remove(2);
    assert_eq!(map.features.len(), 7);
    let valid = map
        .features
        .iter()
        .filter(|feature| feature.unique_id.has_valid_unique_id())
        .count();
    let distinct: std::collections::BTreeSet<u64> = map
        .features
        .iter()
        .map(|feature| feature.unique_id)
        .filter(|id| *id != 0)
        .collect();
    // "RandomAccessContainer has size()==7, num_valid_unique_id==4,
    //  uniqueid_to_index_.size()==3"
    assert_eq!((map.features.len(), valid, distinct.len()), (7, 4, 3));
    assert!(
        map.unique_id_to_index(distinct.iter().next().copied().unwrap())
            .is_err()
    );
}

/// `template <typename Type> Size applyMemberFunction(Size(Type::*)())`
/// (line 654) and the const overload (line 673).
#[test]
fn fm_apply_member_function() {
    let mut generator = generator();
    let mut map = FeatureMap::from_features(vec![Feature::default(), Feature::default()]);
    map.features[1].subordinates.push(Feature::default());

    let invalid = |map: &FeatureMap| {
        map.count_unique_ids(|id| usize::from(id.has_invalid_unique_id()))
            .unwrap()
    };
    let valid = |map: &FeatureMap| {
        map.count_unique_ids(|id| usize::from(id.has_valid_unique_id()))
            .unwrap()
    };

    assert_eq!(invalid(&map), 4);
    map.unique_id.assign_new_unique_id(&mut generator);
    assert_eq!(invalid(&map), 3);
    map.for_each_unique_id(|id| id.assign_new_unique_id(&mut generator))
        .unwrap();
    assert_eq!(valid(&map), 4);
    assert_eq!(invalid(&map), 0);
    map.features[0].base.unique_id.clear_unique_id();
    assert_eq!(valid(&map), 3);
    assert_eq!(invalid(&map), 1);

    // The mutating overload accumulates identically.
    assert_eq!(
        map.for_each_unique_id(|id| usize::from(id.has_invalid_unique_id()))
            .unwrap(),
        1
    );
}

/// `AnnotationStatistics getAnnotationStatistics() const` (line 693).
#[test]
fn fm_annotation_statistics() {
    let [feature1, feature2, feature3, feature4] = source_features();
    let mut map = FeatureMap::new();
    let mut expected = AnnotationStatistics::new();
    assert_eq!(map.annotation_statistics(None).unwrap(), expected);

    map.features.push(feature1); // single hit
    expected.add(AnnotationState::Single).unwrap();
    assert_eq!(map.annotation_statistics(None).unwrap(), expected);

    map.features.push(feature4.clone()); // single hit + no hit
    expected.add(AnnotationState::None).unwrap();
    assert_eq!(map.annotation_statistics(None).unwrap(), expected);

    map.features.push(feature4); // single hit + 2x no hit
    expected.add(AnnotationState::None).unwrap();
    assert_eq!(map.annotation_statistics(None).unwrap(), expected);

    map.features.push(feature2); // + multi-hit (same)
    expected.add(AnnotationState::MultipleSame).unwrap();
    assert_eq!(map.annotation_statistics(None).unwrap(), expected);

    map.features.push(feature3); // + multi (divergent)
    expected.add(AnnotationState::MultipleDivergent).unwrap();
    assert_eq!(map.annotation_statistics(None).unwrap(), expected);

    assert_eq!(expected.count(AnnotationState::None), 2);
    assert_eq!(expected.states(), &[2, 1, 1, 1]);
}

/// Native: the `operator<<` layout of `AnnotationStatistics`, printed by the
/// source section above through `std::cout << res`.
#[test]
fn fm_annotation_statistics_display() {
    let mut stats = AnnotationStatistics::new();
    stats += AnnotationState::Single;
    stats += AnnotationState::Single;
    assert_eq!(
        stats.to_string(),
        "Feature annotation with identifications:\n    no ID: 0\n    single ID: 2\n    \
         multiple IDs (identical): 0\n    multiple IDs (divergent): 0\n\n"
    );
}

/// `[EXTRA] ExposedVector Ctor` (line 732): `FeatureMap fm(10)` and
/// `FeatureMap fm2(10, f4)`.
#[test]
fn fm_exposed_vector_constructors() {
    let map = FeatureMap::from_features(vec![Feature::default(); 10]);
    assert_eq!(map.features.len(), 10);

    let mut f4 = Feature::default();
    f4.base.rt = 5.25;
    let map2 = FeatureMap::from_features(vec![f4; 10]);
    assert_eq!(map2.features.len(), 10);
    assert_eq!(map2.features[6].rt, 5.25);
}

/// Native: the `operator<<` layout of `FeatureMap`.
#[test]
fn fm_display() {
    let mut feature = Feature::new(2.0, 3.0, 1.0);
    feature.base.quality = 0.5;
    feature.base.charge = 2;
    feature.base.unique_id = 42;
    let map = FeatureMap::from_features(vec![feature]);
    assert_eq!(
        map.to_string(),
        "# -- DFEATUREMAP BEGIN --\n# POS \tINTENS\tOVALLQ\tCHARGE\tUniqueID\n\
         2 3\t1\t0.5\t2\t42\n# -- DFEATUREMAP END --\n"
    );
}

/// Native: protein-identification lookup by identifier on both containers.
#[test]
fn find_protein_identification_returns_the_first_match() {
    let mut map = FeatureMap::new();
    assert!(map.find_protein_identification("run").is_none());
    for name in ["other", "run", "run"] {
        map.protein_identifications.push(ProteinIdentification {
            identifier: name.into(),
            search_engine: name.into(),
            ..ProteinIdentification::default()
        });
    }
    assert_eq!(
        map.find_protein_identification("run").unwrap().identifier,
        "run"
    );
    map.find_protein_identification_mut("run")
        .unwrap()
        .score_type = "q-value".into();
    assert_eq!(map.protein_identifications[1].score_type, "q-value");
    assert_eq!(map.protein_identifications[2].score_type, "");

    let mut consensus = ConsensusMap::new();
    assert!(consensus.find_protein_identification("run").is_none());
    consensus.protein_identifications = map.protein_identifications.clone();
    assert_eq!(
        consensus
            .find_protein_identification("run")
            .unwrap()
            .score_type,
        "q-value"
    );
    consensus
        .find_protein_identification_mut("other")
        .unwrap()
        .score_type = "e-value".into();
    assert_eq!(consensus.protein_identifications[0].score_type, "e-value");
}

/// Native: `setPrimaryMSRunPath` / `getPrimaryMSRunPath` on a feature map,
/// including the `UNKNOWN` placeholder the source pushes when nothing is
/// annotated.
#[test]
fn fm_primary_ms_run_path() {
    let mut map = FeatureMap::new();
    assert_eq!(
        map.primary_ms_run_path().unwrap(),
        vec!["UNKNOWN".to_string()]
    );

    map.set_primary_ms_run_path(&[]).unwrap();
    assert_eq!(
        map.primary_ms_run_path().unwrap(),
        vec!["UNKNOWN".to_string()]
    );
    assert!(map.metadata.contains_key("spectra_data"));

    map.set_primary_ms_run_path(&["a.mzML".into(), "b.txt".into()])
        .unwrap();
    assert_eq!(
        map.primary_ms_run_path().unwrap(),
        vec!["a.mzML".to_string(), "b.txt".to_string()]
    );

    // A non-list annotation is a checked error rather than a silent cast.
    map.metadata
        .insert("spectra_data".into(), MetaValue::from(3_i64));
    assert!(map.primary_ms_run_path().is_err());
}

/// Native: the experiment overload falls back to the given list whenever the
/// experiment does not name exactly one existing mzML file.
#[test]
fn primary_ms_run_path_from_experiment_falls_back() {
    let experiment = MSExperiment::new();
    let fallback = vec!["fallback.mzML".to_string()];

    let mut map = FeatureMap::new();
    map.set_primary_ms_run_path_from_experiment(&fallback, &experiment)
        .unwrap();
    assert_eq!(map.primary_ms_run_path().unwrap(), fallback);

    let mut consensus = ConsensusMap::new();
    consensus
        .set_primary_ms_run_path_from_experiment(&fallback, &experiment)
        .unwrap();
    assert_eq!(consensus.primary_ms_run_path(), fallback);
}

// ===========================================================================
// ConsensusMap_test.cpp
// ===========================================================================

/// `ConsensusMap()` (line 34) and `~ConsensusMap()` (line 40).
#[test]
fn cm_default_constructor_and_destructor() {
    let map = ConsensusMap::new();
    assert!(map.metadata.is_empty());
    assert_eq!(map.features.len(), 0);
    assert_eq!(map.experiment_type, "label-free");
    drop(map);
}

/// `getProteinIdentifications() const` (line 44), non-const (line 49) and
/// `setProteinIdentifications` (line 55).
///
/// The source sections build a `FeatureMap` although they live in the
/// ConsensusMap test; both containers are exercised here.
#[test]
fn cm_protein_identification_accessors() {
    let mut map = FeatureMap::new();
    assert_eq!(map.protein_identifications.len(), 0);
    map.protein_identifications
        .resize(1, ProteinIdentification::default());
    assert_eq!(map.protein_identifications.len(), 1);
    let mut map = FeatureMap::new();
    map.protein_identifications = vec![ProteinIdentification::default(); 2];
    assert_eq!(map.protein_identifications.len(), 2);

    let mut consensus = ConsensusMap::new();
    assert_eq!(consensus.protein_identifications.len(), 0);
    consensus.protein_identifications = vec![ProteinIdentification::default(); 2];
    assert_eq!(consensus.protein_identifications.len(), 2);
}

/// `getUnassignedPeptideIdentifications() const` (line 61), non-const
/// (line 66) and the setter (line 72), again written against `FeatureMap`.
#[test]
fn cm_unassigned_peptide_identification_accessors() {
    let mut map = FeatureMap::new();
    assert_eq!(map.unassigned_peptide_identifications.len(), 0);
    map.unassigned_peptide_identifications
        .resize(1, PeptideIdentification::default());
    assert_eq!(map.unassigned_peptide_identifications.len(), 1);
    let mut map = FeatureMap::new();
    map.unassigned_peptide_identifications = vec![PeptideIdentification::default(); 2];
    assert_eq!(map.unassigned_peptide_identifications.len(), 2);

    let mut consensus = ConsensusMap::new();
    consensus.unassigned_peptide_identifications = vec![PeptideIdentification::default(); 2];
    assert_eq!(consensus.unassigned_peptide_identifications.len(), 2);
}

/// `getDataProcessing() const` (line 78), non-const (line 83) and
/// `setDataProcessing` (line 89).
#[test]
fn cm_data_processing_accessors() {
    let mut map = ConsensusMap::new();
    assert_eq!(map.data_processing.len(), 0);
    map.data_processing.resize(1, DataProcessing::default());
    assert_eq!(map.data_processing.len(), 1);
    let mut map = ConsensusMap::new();
    map.data_processing = vec![DataProcessing::default()];
    assert_eq!(map.data_processing.len(), 1);
}

/// `void updateRanges()` (line 117).
#[test]
fn cm_update_ranges() {
    let [mut feature1, mut feature2, mut feature3, mut feature4] = consensus_source_features();
    feature1.base.unique_id = 1;
    feature2.base.unique_id = 2;
    feature3.base.unique_id = 3;
    feature4.base.unique_id = 4;

    let mut map = ConsensusMap::new();
    let mut f = ConsensusFeature::new();
    f.base.intensity = 1.0;
    f.base.rt = 2.0;
    f.base.mz = 3.0;
    f.insert(FeatureHandle::new(1, &feature1.base)).unwrap();
    map.features.push(f.clone());

    for _ in 0..2 {
        let ranges = map.ranges().unwrap();
        close(ranges.intensity.unwrap().min, 1.0);
        close(ranges.intensity.unwrap().max, 1.0);
        close(ranges.rt.unwrap().max, 2.0);
        close(ranges.mz.unwrap().max, 3.0);
        close(ranges.rt.unwrap().min, 2.0);
        close(ranges.mz.unwrap().min, 3.0);
    }

    f.insert(FeatureHandle::new(1, &feature2.base)).unwrap();
    map.features.push(f.clone());
    let ranges = map.ranges().unwrap();
    close(ranges.intensity.unwrap().min, 0.5);
    close(ranges.intensity.unwrap().max, 1.0);
    close(ranges.rt.unwrap().max, 2.0);
    close(ranges.mz.unwrap().max, 3.0);
    close(ranges.rt.unwrap().min, 0.0);
    close(ranges.mz.unwrap().min, 2.5);

    f.insert(FeatureHandle::new(1, &feature3.base)).unwrap();
    f.insert(FeatureHandle::new(1, &feature4.base)).unwrap();
    map.features.push(f);
    let ranges = map.ranges().unwrap();
    close(ranges.intensity.unwrap().min, f64::from(0.01f32));
    close(ranges.intensity.unwrap().max, 1.0);
    close(ranges.rt.unwrap().max, 10.5);
    close(ranges.mz.unwrap().max, 3.0);
    close(ranges.rt.unwrap().min, 0.0);
    close(ranges.mz.unwrap().min, 0.0);
}

/// `ConsensusMap& appendRows(const ConsensusMap& rhs)` (line 171).
#[test]
fn cm_append_rows() {
    let mut generator = generator();
    let mut m1 = ConsensusMap::new();
    let mut m2 = ConsensusMap::new();
    let mut m3 = ConsensusMap::new();

    m1.append_rows(&m2, &mut generator).unwrap();
    assert_eq!(m1, m3);

    let mut f1 = ConsensusFeature::new();
    f1.base.mz = 100.12;
    m1.features.push(f1);
    m3 = m1.clone();
    m1.append_rows(&m2, &mut generator).unwrap();
    assert_eq!(m1, m3);

    m1.identifier = "123".into();
    m1.data_processing.resize(1, DataProcessing::default());
    m1.protein_identifications
        .resize(1, ProteinIdentification::default());
    m1.unassigned_peptide_identifications
        .resize(1, PeptideIdentification::default());
    m1.unique_id.ensure_unique_id(&mut generator);
    m1.column_headers.entry(0).or_default().filename = "m1".into();

    m2.identifier = "321".into();
    m2.data_processing.resize(2, DataProcessing::default());
    m2.protein_identifications
        .resize(2, ProteinIdentification::default());
    m2.unassigned_peptide_identifications
        .resize(2, PeptideIdentification::default());
    m2.features.push(ConsensusFeature::new());
    m2.features.push(ConsensusFeature::new());
    m2.column_headers.entry(1).or_default().filename = "m2".into();

    m1.append_rows(&m2, &mut generator).unwrap();
    assert_eq!(m1.identifier, "");
    assert!(!m1.unique_id.has_valid_unique_id());
    assert_eq!(m1.data_processing.len(), 3);
    assert_eq!(m1.protein_identifications.len(), 3);
    assert_eq!(m1.unassigned_peptide_identifications.len(), 3);
    assert_eq!(m1.features.len(), 3);
    assert_eq!(m1.column_headers.len(), 2);
    // Source quirk: the positional zip renames only the first merged header.
    assert_eq!(m1.column_headers[&0].filename, "mergedConsensusXMLFile");
    assert_eq!(m1.column_headers[&1].filename, "m2");
}

/// `ConsensusMap& appendColumns(const ConsensusMap& rhs)` (line 213).
#[test]
fn cm_append_columns() {
    let mut generator = generator();
    let mut m1 = ConsensusMap::new();
    let mut m2 = ConsensusMap::new();

    // Test1: adding an empty map has no effect.
    m1.append_columns(&ConsensusMap::new(), &mut generator)
        .unwrap();
    assert_eq!(m1, ConsensusMap::new());

    let f1 = Feature::new(1.0, 1.0, 1.0);
    let mut cf1 = ConsensusFeature::new();
    cf1.insert(FeatureHandle::from_peak(
        0,
        Peak2D::new(f1.rt, f1.mz, f1.intensity),
        0,
    ))
    .unwrap();
    cf1.base.mz = 100.12;
    m1.features.push(cf1);

    // Test2: adding an empty map to a map with content.
    let old_m1 = m1.clone();
    m1.append_columns(&ConsensusMap::new(), &mut generator)
        .unwrap();
    assert_eq!(m1, old_m1);

    m1.identifier = "123".into();
    m1.data_processing.resize(1, DataProcessing::default());
    m1.protein_identifications
        .resize(1, ProteinIdentification::default());
    m1.unassigned_peptide_identifications
        .resize(1, PeptideIdentification::default());
    m1.unique_id.ensure_unique_id(&mut generator);
    m1.column_headers.entry(0).or_default().filename = "m1".into();

    m2.identifier = "321".into();
    m2.data_processing.resize(2, DataProcessing::default());
    m2.protein_identifications
        .resize(2, ProteinIdentification::default());
    m2.unassigned_peptide_identifications
        .resize(2, PeptideIdentification::default());
    m2.column_headers.entry(0).or_default().filename = "m2_1".into();
    m2.column_headers.entry(1).or_default().filename = "m2_2".into();

    let f2 = Feature::new(2.0, 2.0, 2.0);
    let f3 = Feature::new(3.0, 3.0, 3.0);
    let mut cf2 = ConsensusFeature::new();
    cf2.insert(FeatureHandle::from_peak(
        0,
        Peak2D::new(f2.rt, f2.mz, f2.intensity),
        0,
    ))
    .unwrap();
    cf2.insert(FeatureHandle::from_peak(
        1,
        Peak2D::new(f3.rt, f3.mz, f3.intensity),
        1,
    ))
    .unwrap();
    m2.features.push(cf2);

    m1.append_columns(&m2, &mut generator).unwrap();

    assert_eq!(m1.identifier, "");
    assert!(!m1.unique_id.has_valid_unique_id());
    assert_eq!(m1.data_processing.len(), 3);
    assert_eq!(m1.protein_identifications.len(), 3);
    assert_eq!(m1.unassigned_peptide_identifications.len(), 3);
    assert_eq!(m1.features.len(), 2);
    assert_eq!(m1.column_headers.len(), 3);
    assert_eq!(m1.column_headers[&0].filename, "m1");
    assert_eq!(m1.column_headers[&1].filename, "m2_1");
    assert_eq!(m1.column_headers[&2].filename, "m2_2");

    let first = m1.features[0].handles();
    assert_eq!(first[0].intensity, 1.0);
    assert_eq!(first[0].unique_id, 0);
    assert_eq!(first[0].map_index, 0);

    let second = m1.features[1].handles();
    assert_eq!(second[0].intensity, 2.0);
    assert_eq!(second[0].unique_id, 0);
    assert_eq!(second[0].map_index, 1);
    assert_eq!(second[1].intensity, 3.0);
    assert_eq!(second[1].unique_id, 1);
    assert_eq!(second[1].map_index, 2);
}

/// Native: `appendColumns` shifts the `map_index` meta value of assigned and
/// unassigned peptide identifications alike.
#[test]
fn cm_append_columns_shifts_map_index_annotations() {
    let mut generator = generator();
    let mut left = ConsensusMap::new();
    left.column_headers.entry(0).or_default().filename = "left.mzML".into();
    left.column_headers.entry(1).or_default().filename = "left2.mzML".into();

    let mut right = ConsensusMap::new();
    right.column_headers.entry(0).or_default().filename = "right.mzML".into();
    let mut assigned = identification("AAA");
    assigned.metadata.insert("map_index".into(), 0_i64.into());
    let mut unassigned = identification("KKK");
    unassigned.metadata.insert("map_index".into(), 0_i64.into());
    let mut feature = ConsensusFeature::new();
    feature.base.peptide_identifications = vec![assigned];
    right.features.push(feature);
    right.unassigned_peptide_identifications = vec![unassigned];

    left.append_columns(&right, &mut generator).unwrap();
    assert_eq!(
        left.features[0].base.peptide_identifications[0].metadata["map_index"]
            .as_i64()
            .unwrap(),
        2
    );
    assert_eq!(
        left.unassigned_peptide_identifications[0].metadata["map_index"]
            .as_i64()
            .unwrap(),
        2
    );
    assert_eq!(left.column_headers.len(), 3);
}

/// `ConsensusMap& operator=(const ConsensusMap&)` (line 327),
/// `ConsensusMap(const ConsensusMap&)` (line 356),
/// `ConsensusMap& operator=(ConsensusMap&&)` (line 372) and
/// `ConsensusMap(ConsensusMap&&)` (line 402).
#[test]
fn cm_copy_assign_and_move() {
    let source = map_const_1();
    let check = |map: &ConsensusMap| {
        assert_eq!(map.features.len(), 3);
        assert_eq!(map.identifier, "lsid");
        assert_eq!(map.metadata["meta"].as_str().unwrap(), "value");
        assert_eq!(map.column_headers[&0].filename, "blub");
        assert_eq!(map.column_headers[&0].label, "label");
        assert_eq!(map.column_headers[&0].size, 47);
        assert_eq!(
            map.column_headers[&0].metadata["meta"].as_str().unwrap(),
            "meta"
        );
        assert_eq!(map.experiment_type, "labeled_MS2");
        assert_eq!(map.data_processing.len(), 1);
        assert_eq!(map.protein_identifications.len(), 1);
        assert_eq!(map.unassigned_peptide_identifications.len(), 1);
    };

    let mut map2 = ConsensusMap::new();
    assert!(map2.is_empty());
    map2 = source.clone();
    check(&map2);
    // Copy construction.
    check(&source.clone());
    // Move construction and move assignment.
    let moved_source = source.clone();
    let map_moved = moved_source;
    check(&map_moved);
    let map3 = source.clone();
    check(&map3);

    // Assignment of an empty object.
    map2 = ConsensusMap::new();
    assert_eq!(map2.identifier, "");
    assert_eq!(map2.column_headers.len(), 0);
    assert_eq!(map2.experiment_type, "label-free");
    assert_eq!(map2.data_processing.len(), 0);
    assert_eq!(map2.protein_identifications.len(), 0);
    assert_eq!(map2.unassigned_peptide_identifications.len(), 0);
}

/// `ConsensusMap(size_type n)` (line 419).
#[test]
fn cm_size_constructor() {
    let map = ConsensusMap::with_size(5).unwrap();
    assert_eq!(map.features.len(), 5);
    assert_eq!(map.experiment_type, "label-free");
    assert!(ConsensusMap::with_size(ConsensusMap::MAX_ITEMS + 1).is_err());
}

/// `[ConsensusMap::ColumnHeader] ColumnHeader()` (line 427).
#[test]
fn cm_column_header_default_constructor() {
    let header = ColumnHeader::default();
    assert_eq!(header.filename, "");
    assert_eq!(header.label, "");
    assert_eq!(header.size, 0);
    // The source default is UniqueIdInterface::INVALID, which is zero.
    assert_eq!(header.unique_id, 0);
    assert!(header.metadata.is_empty());
}

/// `const ColumnHeaders& getColumnHeaders() const` (line 434) and the
/// non-const overload (line 440).
#[test]
fn cm_column_header_accessors() {
    let mut map = ConsensusMap::new();
    assert_eq!(map.column_headers.len(), 0);
    map.column_headers.entry(0).or_default().filename = "blub".into();
    assert_eq!(map.column_headers[&0].filename, "blub");
}

/// `const std::string& getExperimentType() const` (line 447) and
/// `void setExperimentType(const std::string&)` (line 452).
#[test]
fn cm_experiment_type() {
    let mut map = ConsensusMap::new();
    assert_eq!(map.experiment_type, "label-free");
    map.set_experiment_type("labeled_MS2").unwrap();
    assert_eq!(map.experiment_type, "labeled_MS2");
    map.set_experiment_type("labeled_MS1").unwrap();
    map.set_experiment_type("label-free").unwrap();
    // The source throws Exception::IllegalArgument for anything else.
    assert!(map.set_experiment_type("SILAC").is_err());
    assert_eq!(map.experiment_type, "label-free");
}

/// `void swap(ConsensusMap& from)` (line 458).
#[test]
fn cm_swap() {
    let mut map1 = ConsensusMap::new();
    let mut f = ConsensusFeature::new();
    f.insert(FeatureHandle::new(1, &Feature::default().base))
        .unwrap();
    map1.features.push(f);
    let header = map1.column_headers.entry(1).or_default();
    header.filename = "bla".into();
    header.size = 5;
    map1.identifier = "LSID".into();
    map1.set_experiment_type("labeled_MS2").unwrap();
    map1.data_processing.resize(1, DataProcessing::default());
    map1.protein_identifications
        .resize(1, ProteinIdentification::default());
    map1.unassigned_peptide_identifications
        .resize(1, PeptideIdentification::default());
    map1.metadata.insert("keep".into(), "left".into());

    let mut map2 = ConsensusMap::new();
    map1.swap(&mut map2);

    assert_eq!(map1.features.len(), 0);
    assert_eq!(map1.column_headers.len(), 0);
    assert_eq!(map1.identifier, "");
    assert_eq!(map1.data_processing.len(), 0);
    assert_eq!(map1.protein_identifications.len(), 0);
    assert_eq!(map1.unassigned_peptide_identifications.len(), 0);

    assert_eq!(map2.features.len(), 1);
    assert_eq!(map2.column_headers.len(), 1);
    assert_eq!(map2.identifier, "LSID");
    assert_eq!(map2.experiment_type, "labeled_MS2");
    assert_eq!(map2.data_processing.len(), 1);
    assert_eq!(map2.protein_identifications.len(), 1);
    assert_eq!(map2.unassigned_peptide_identifications.len(), 1);

    // Source quirk: meta values stay with their map.
    assert_eq!(map1.metadata["keep"].as_str().unwrap(), "left");
    assert!(!map2.metadata.contains_key("keep"));
}

/// `bool operator==(const ConsensusMap&) const` (line 489) and
/// `bool operator!=(const ConsensusMap&) const` (line 537).
///
/// As for `FeatureMap`, the last case diverges: the source's `operator==`
/// compares the cached ranges that `clear(false)` leaves in place, and this
/// port has no such cache.
#[test]
fn cm_equality_and_inequality() {
    let [feature1, feature2, _, _] = consensus_source_features();
    let empty = ConsensusMap::new();
    let mut edit = ConsensusMap::new();
    assert!(empty == edit);
    assert!(!(empty != edit));

    edit.identifier = "lsid".into();
    assert!(empty != edit);

    edit = empty.clone();
    edit.features
        .push(ConsensusFeature::from(feature1.base.clone()));
    assert!(empty != edit);

    edit = empty.clone();
    edit.data_processing.resize(1, DataProcessing::default());
    assert!(empty != edit);

    edit = empty.clone();
    edit.metadata
        .insert("bla".into(), MetaValue::try_from(4.1_f64).unwrap());
    assert!(empty != edit);

    edit = empty.clone();
    edit.column_headers.entry(0).or_default().filename = "bla".into();
    assert!(empty != edit);

    edit = empty.clone();
    edit.set_experiment_type("labeled_MS2").unwrap();
    assert!(empty != edit);

    edit = empty.clone();
    edit.protein_identifications
        .resize(10, ProteinIdentification::default());
    assert!(empty != edit);

    edit = empty.clone();
    edit.unassigned_peptide_identifications
        .resize(10, PeptideIdentification::default());
    assert!(empty != edit);

    edit = empty.clone();
    edit.features
        .push(ConsensusFeature::from(feature1.base.clone()));
    edit.features
        .push(ConsensusFeature::from(feature2.base.clone()));
    edit.ranges().unwrap();
    edit.clear(false);
    // Source expectation: `empty == edit` is false, for the cached-range reason.
    assert_eq!(empty, edit);
}

/// `void sortByIntensity(bool reverse=false)` (line 582),
/// `sortByRT()` (588), `sortByMZ()` (594), `sortByPosition()` (600),
/// `sortByQuality(bool)` (606), `sortBySize()` (612) and `sortByMaps()` (618).
///
/// All seven are `NOT_TESTABLE` upstream ("tested within TOPP TextExporter");
/// they are ported as real tests rather than mapped away.
#[test]
fn cm_sorting() {
    let make = |rt: f64, mz: f64, intensity: f32, quality: f32| {
        let mut feature = ConsensusFeature::new();
        feature.base.rt = rt;
        feature.base.mz = mz;
        feature.base.intensity = intensity;
        feature.base.quality = quality;
        feature
    };
    let mut map = ConsensusMap::from_features(vec![
        make(3.0, 30.0, 3.0, 30.0),
        make(1.0, 10.0, 1.0, 10.0),
        make(2.0, 20.0, 2.0, 20.0),
    ]);

    map.sort_by_intensity(false).unwrap();
    assert_eq!(
        map.features.iter().map(|f| f.intensity).collect::<Vec<_>>(),
        vec![1.0, 2.0, 3.0]
    );
    map.sort_by_intensity(true).unwrap();
    assert_eq!(
        map.features.iter().map(|f| f.intensity).collect::<Vec<_>>(),
        vec![3.0, 2.0, 1.0]
    );
    map.sort_by_rt().unwrap();
    assert_eq!(
        map.features.iter().map(|f| f.rt).collect::<Vec<_>>(),
        vec![1.0, 2.0, 3.0]
    );
    map.sort_by_intensity(true).unwrap();
    map.sort_by_mz().unwrap();
    assert_eq!(
        map.features.iter().map(|f| f.mz).collect::<Vec<_>>(),
        vec![10.0, 20.0, 30.0]
    );
    map.sort_by_intensity(true).unwrap();
    map.sort_by_position().unwrap();
    assert_eq!(
        map.features.iter().map(|f| f.rt).collect::<Vec<_>>(),
        vec![1.0, 2.0, 3.0]
    );
    map.sort_by_quality(false).unwrap();
    assert_eq!(
        map.features.iter().map(|f| f.quality).collect::<Vec<_>>(),
        vec![10.0, 20.0, 30.0]
    );
    map.sort_by_quality(true).unwrap();
    assert_eq!(
        map.features.iter().map(|f| f.quality).collect::<Vec<_>>(),
        vec![30.0, 20.0, 10.0]
    );

    // sortBySize orders by decreasing number of handles; sortByMaps orders by
    // the lexicographic handle identity list.
    let mut sized = ConsensusMap::new();
    let mut one = ConsensusFeature::new();
    one.insert(FeatureHandle::from_peak(5, Peak2D::new(0.0, 0.0, 0.0), 1))
        .unwrap();
    let mut two = ConsensusFeature::new();
    two.insert(FeatureHandle::from_peak(1, Peak2D::new(0.0, 0.0, 0.0), 1))
        .unwrap();
    two.insert(FeatureHandle::from_peak(2, Peak2D::new(0.0, 0.0, 0.0), 1))
        .unwrap();
    sized.features.push(one);
    sized.features.push(two);
    sized.sort_by_size().unwrap();
    assert_eq!(sized.features[0].len(), 2);
    sized.sort_by_maps().unwrap();
    assert_eq!(sized.features[0].handles()[0].map_index, 1);
}

/// `void sortPeptideIdentificationsByMapIndex()` (line 624), upstream
/// `NOT_TESTABLE` ("tested within TOPP IDMapper").
#[test]
fn cm_sort_peptide_identifications_by_map_index() {
    let annotated = |sequence: &str, index: Option<i64>| {
        let mut record = identification(sequence);
        if let Some(index) = index {
            record.metadata.insert("map_index".into(), index.into());
        }
        record
    };
    let mut feature = ConsensusFeature::new();
    feature.base.peptide_identifications = vec![
        annotated("AAA", Some(3)),
        annotated("CCC", None),
        annotated("DDD", Some(1)),
        annotated("EEE", None),
        annotated("BBB", Some(1)),
    ];
    let mut map = ConsensusMap::from_features(vec![feature]);
    map.sort_peptide_identifications_by_map_index().unwrap();
    let order: Vec<String> = map.features[0]
        .base
        .peptide_identifications
        .iter()
        .map(|record| record.hits[0].sequence.to_string())
        .collect();
    // Annotated first, ascending and stable; unannotated last in input order.
    assert_eq!(order, vec!["DDD", "BBB", "AAA", "CCC", "EEE"]);

    // A non-integer map index is a checked error and changes nothing.
    let before = map.clone();
    map.features[0].base.peptide_identifications[0]
        .metadata
        .insert("map_index".into(), "first".into());
    assert!(map.sort_peptide_identifications_by_map_index().is_err());
    assert_eq!(
        map.features[0].base.peptide_identifications[1].hits[0]
            .sequence
            .to_string(),
        before.features[0].base.peptide_identifications[1].hits[0]
            .sequence
            .to_string()
    );
}

/// `void clear(bool clear_meta_data = true)` (line 630).
#[test]
fn cm_clear() {
    let mut map = ConsensusMap::new();
    let mut f = ConsensusFeature::new();
    f.insert(FeatureHandle::new(1, &Feature::default().base))
        .unwrap();
    map.features.push(f);
    let header = map.column_headers.entry(1).or_default();
    header.filename = "bla".into();
    header.size = 5;
    map.identifier = "LSID".into();
    map.set_experiment_type("labeled_MS2").unwrap();
    map.data_processing.resize(1, DataProcessing::default());
    map.protein_identifications
        .resize(1, ProteinIdentification::default());
    map.unassigned_peptide_identifications
        .resize(1, PeptideIdentification::default());

    map.clear(false);
    assert_eq!(map.features.len(), 0);
    assert!(map != ConsensusMap::new());
    assert!(map.is_empty());

    map.clear(true);
    assert_eq!(map.features.len(), 0);
    assert_eq!(map, ConsensusMap::new());
    assert!(map.is_empty());
}

/// `template <typename Type> Size applyMemberFunction(Size(Type::*)())`
/// (line 656) and the const overload (line 675).
#[test]
fn cm_apply_member_function() {
    let mut generator = generator();
    let mut map = ConsensusMap::from_features(vec![ConsensusFeature::new(); 3]);

    let invalid = |map: &ConsensusMap| {
        map.count_unique_ids(|id| usize::from(id.has_invalid_unique_id()))
            .unwrap()
    };
    let valid = |map: &ConsensusMap| {
        map.count_unique_ids(|id| usize::from(id.has_valid_unique_id()))
            .unwrap()
    };

    assert_eq!(invalid(&map), 4);
    map.unique_id.assign_new_unique_id(&mut generator);
    assert_eq!(invalid(&map), 3);
    map.for_each_unique_id(|id| id.assign_new_unique_id(&mut generator))
        .unwrap();
    assert_eq!(valid(&map), 4);
    assert_eq!(invalid(&map), 0);
    map.features[0].base.unique_id.clear_unique_id();
    assert_eq!(valid(&map), 3);
    assert_eq!(invalid(&map), 1);
    assert_eq!(
        map.for_each_unique_id(|id| usize::from(id.has_invalid_unique_id()))
            .unwrap(),
        1
    );
}

/// The `split` fixture of `ConsensusMap_test.cpp` (lines 697-730).
fn split_fixture() -> ConsensusMap {
    let mut map = ConsensusMap::new();
    map.column_headers.entry(0).or_default().filename = "file.FeatureXML".into();
    map.column_headers.entry(1).or_default().filename = "file2.FeatureXML".into();

    let mut cf1 = ConsensusFeature::new();
    cf1.insert(FeatureHandle::from_peak(
        0,
        Peak2D::new(10.0, 433.33, 100_000.0),
        0,
    ))
    .unwrap();
    cf1.insert(FeatureHandle::from_peak(
        1,
        Peak2D::new(11.0, 434.33, 200_000.0),
        0,
    ))
    .unwrap();
    let mut id1 = PeptideIdentification {
        rt: Some(10.0),
        hits: vec![PeptideHit::new(0.1, 1, 3, AASequence::parse("AAA").unwrap()).unwrap()],
        ..PeptideIdentification::default()
    };
    id1.metadata.insert("map_index".into(), 0_i64.into());
    cf1.base.peptide_identifications.push(id1);
    cf1.base
        .metadata
        .insert("test".into(), "some information".into());
    map.features.push(cf1);

    let mut cf2 = ConsensusFeature::new();
    cf2.insert(FeatureHandle::from_peak(
        0,
        Peak2D::new(20.0, 433.33, 300_000.0),
        0,
    ))
    .unwrap();
    cf2.insert(FeatureHandle::from_peak(
        1,
        Peak2D::new(21.0, 433.33, 400_000.0),
        0,
    ))
    .unwrap();
    let mut id2 = PeptideIdentification {
        rt: Some(20.0),
        hits: vec![PeptideHit::new(0.1, 1, 3, AASequence::parse("WWW").unwrap()).unwrap()],
        ..PeptideIdentification::default()
    };
    id2.metadata.insert("map_index".into(), 1_i64.into());
    cf2.base.peptide_identifications.push(id2);
    map.features.push(cf2);

    let mut uid1 = PeptideIdentification {
        hits: vec![PeptideHit::new(0.1, 1, 3, AASequence::parse("LLL").unwrap()).unwrap()],
        ..PeptideIdentification::default()
    };
    uid1.metadata.insert("map_index".into(), 0_i64.into());
    let mut uid2 = PeptideIdentification {
        hits: vec![PeptideHit::new(0.1, 1, 3, AASequence::parse("KKK").unwrap()).unwrap()],
        ..PeptideIdentification::default()
    };
    uid2.metadata.insert("map_index".into(), 1_i64.into());
    map.unassigned_peptide_identifications.push(uid1);
    map.unassigned_peptide_identifications.push(uid2);
    map
}

/// A `DataProcessing` record naming `IsobaricAnalyzer`, as the source section
/// builds before its second case.
fn isobaric_processing() -> DataProcessing {
    let mut processing = DataProcessing {
        software: Software {
            name: "IsobaricAnalyzer".into(),
            ..Software::default()
        },
        ..DataProcessing::default()
    };
    processing
        .actions
        .insert(openms::metadata::ProcessingAction::Quantitation);
    processing
}

/// `void split(std::vector<FeatureMap>&, SplitMeta mode) const` (line 695).
#[test]
fn cm_split() {
    let mut map = split_fixture();
    let sequence = |record: &PeptideIdentification| record.hits[0].sequence.to_string();

    // Non-isobaric data.
    let maps = map.split(SplitMeta::Discard).unwrap();
    assert_eq!(maps.len(), 2);
    assert_eq!(maps[0].features.len(), 2);
    assert_eq!(maps[1].features.len(), 2);
    assert_eq!(maps[0].features[0].rt, 10.0);
    assert_eq!(maps[0].features[0].intensity, 100_000.0);
    assert_eq!(
        sequence(&maps[0].features[0].base.peptide_identifications[0]),
        "AAA"
    );
    assert!(!maps[0].features[0].base.metadata.contains_key("test"));
    assert_eq!(maps[0].features[1].rt, 20.0);
    assert_eq!(maps[0].features[1].intensity, 300_000.0);
    assert!(maps[0].features[1].base.peptide_identifications.is_empty());
    assert_eq!(maps[1].features[0].rt, 11.0);
    assert_eq!(maps[1].features[0].intensity, 200_000.0);
    assert!(maps[1].features[0].base.peptide_identifications.is_empty());
    assert_eq!(maps[1].features[1].rt, 21.0);
    assert_eq!(maps[1].features[1].intensity, 400_000.0);
    assert_eq!(
        sequence(&maps[1].features[1].base.peptide_identifications[0]),
        "WWW"
    );
    assert_eq!(
        sequence(&maps[0].unassigned_peptide_identifications[0]),
        "LLL"
    );
    assert_eq!(
        sequence(&maps[1].unassigned_peptide_identifications[0]),
        "KKK"
    );

    // Isobaric data.
    map.data_processing.push(isobaric_processing());
    let maps = map.split(SplitMeta::Discard).unwrap();
    assert_eq!(maps.len(), 2);
    assert_eq!(maps[0].features.len(), 2);
    assert_eq!(maps[1].features.len(), 2);
    assert_eq!(
        sequence(&maps[0].features[0].base.peptide_identifications[0]),
        "AAA"
    );
    assert_eq!(
        sequence(&maps[0].features[1].base.peptide_identifications[0]),
        "WWW"
    );
    assert_eq!(
        sequence(&maps[0].unassigned_peptide_identifications[0]),
        "LLL"
    );
    assert_eq!(
        sequence(&maps[0].unassigned_peptide_identifications[1]),
        "KKK"
    );
    assert!(maps[1].features[0].base.peptide_identifications.is_empty());
    assert!(maps[1].features[1].base.peptide_identifications.is_empty());

    // Meta-value modes.
    let maps = map.split(SplitMeta::CopyFirst).unwrap();
    assert_eq!(maps.len(), 2);
    assert_eq!(maps[0].features.len(), 2);
    assert_eq!(maps[1].features.len(), 2);
    assert!(maps[0].features[0].base.metadata.contains_key("test"));
    assert_eq!(
        maps[0].features[0].base.metadata["test"].as_str().unwrap(),
        "some information"
    );
    assert!(!maps[1].features[0].base.metadata.contains_key("test"));

    let maps = map.split(SplitMeta::CopyAll).unwrap();
    assert_eq!(maps.len(), 2);
    assert_eq!(maps[0].features.len(), 2);
    assert_eq!(maps[1].features.len(), 2);
    assert!(maps[0].features[0].base.metadata.contains_key("test"));
    assert_eq!(
        maps[0].features[0].base.metadata["test"].as_str().unwrap(),
        "some information"
    );
    assert!(maps[1].features[0].base.metadata.contains_key("test"));
    assert_eq!(
        maps[1].features[0].base.metadata["test"].as_str().unwrap(),
        "some information"
    );
}

/// Native: `CopyFirst` tests the smallest *handle* map index, which the source
/// captures before routing the identifications, so a feature that exists at
/// index 0 only because an identification was routed there does not satisfy it.
#[test]
fn cm_split_copy_first_uses_the_handle_map_index() {
    let mut map = split_fixture();
    // Move the first consensus feature's handles to columns 1 and 2 ...
    map.column_headers.entry(2).or_default().filename = "file3.FeatureXML".into();
    let handles: Vec<FeatureHandle> = map.features[0]
        .handles()
        .iter()
        .map(|handle| FeatureHandle {
            map_index: handle.map_index + 1,
            ..*handle
        })
        .collect();
    map.features[0].set_handles(handles).unwrap();
    // ... while its identification still routes a feature into column 0.
    let maps = map.split(SplitMeta::Discard).unwrap();
    assert_eq!(maps[0].features.len(), 2);
    // Discard succeeds, CopyFirst does not.
    assert!(map.split(SplitMeta::CopyFirst).is_err());
    // CopyAll is unaffected by the handle indices.
    map.split(SplitMeta::CopyAll).unwrap();
}

/// Native: `split` carries the data processing to every produced map, drops the
/// handle unique IDs and keeps the handle charge and width.
#[test]
fn cm_split_carries_processing_and_handle_fields() {
    let mut map = split_fixture();
    map.data_processing.push(DataProcessing::default());
    let mut handles = map.features[0].handles().to_vec();
    handles[0].charge = 2;
    handles[0].width = 0.75;
    handles[0].unique_id = 77;
    map.features[0].set_handles(handles).unwrap();

    let maps = map.split(SplitMeta::Discard).unwrap();
    assert_eq!(maps[0].data_processing.len(), 1);
    assert_eq!(maps[1].data_processing.len(), 1);
    assert_eq!(maps[0].features[0].charge, 2);
    assert_eq!(maps[0].features[0].width, 0.75);
    // The source slices the handle to its Peak2D base, so the ID is not carried.
    assert_eq!(maps[0].features[0].unique_id, 0);
}

/// Native: the error paths of `split`, all of which the source leaves to
/// undefined behaviour or an exception.
#[test]
fn cm_split_error_paths() {
    // A non-isobaric identification without map_index.
    let mut map = split_fixture();
    map.unassigned_peptide_identifications[0]
        .metadata
        .remove("map_index");
    assert!(matches!(
        map.split(SplitMeta::Discard),
        Err(openms::Error::MissingInformation(_))
    ));

    // A map index that does not name a column.
    let mut map = split_fixture();
    map.unassigned_peptide_identifications[0]
        .metadata
        .insert("map_index".into(), 9_i64.into());
    assert!(map.split(SplitMeta::Discard).is_err());

    // Isobaric data whose smallest map index is not zero.
    let mut map = split_fixture();
    map.data_processing.push(isobaric_processing());
    let handles: Vec<FeatureHandle> = map.features[0]
        .handles()
        .iter()
        .map(|handle| FeatureHandle {
            map_index: handle.map_index + 1,
            ..*handle
        })
        .collect();
    map.features[0].set_handles(handles).unwrap();
    assert!(map.split(SplitMeta::Discard).is_err());

    // CopyFirst without a feature at map index 0.
    let mut map = split_fixture();
    let handles: Vec<FeatureHandle> = map.features[0]
        .handles()
        .iter()
        .map(|handle| FeatureHandle {
            map_index: handle.map_index + 1,
            ..*handle
        })
        .collect();
    map.features[0].set_handles(handles).unwrap();
    map.features[0].base.peptide_identifications[0]
        .metadata
        .insert("map_index".into(), 1_i64.into());
    assert!(map.split(SplitMeta::CopyFirst).is_err());
}

/// Native: `ConsensusMap::setPrimaryMSRunPath` writes into the column headers
/// and preserves the source's positional, default-inserting behaviour.
#[test]
fn cm_primary_ms_run_path() {
    let mut map = ConsensusMap::new();
    assert!(map.primary_ms_run_path().is_empty());

    // No headers yet: the paths create columns 0..n-1.
    map.set_primary_ms_run_path(&["a.mzML".into(), "b.mzML".into()])
        .unwrap();
    assert_eq!(
        map.primary_ms_run_path(),
        vec!["a.mzML".to_string(), "b.mzML".to_string()]
    );

    // A mismatched count is refused and changes nothing.
    let before = map.clone();
    assert!(map.set_primary_ms_run_path(&["only.mzML".into()]).is_err());
    assert_eq!(map, before);

    // An empty list renames every column to UNKNOWN.
    map.set_primary_ms_run_path(&[]).unwrap();
    assert_eq!(
        map.primary_ms_run_path(),
        vec!["UNKNOWN".to_string(), "UNKNOWN".to_string()]
    );

    // Source quirk: writing by position into headers keyed otherwise adds
    // columns rather than renaming the existing ones.
    let mut sparse = ConsensusMap::new();
    sparse.column_headers.entry(5).or_default().filename = "five.mzML".into();
    sparse
        .set_primary_ms_run_path(&["zero.mzML".into()])
        .unwrap();
    assert_eq!(sparse.column_headers.len(), 2);
    assert_eq!(sparse.column_headers[&0].filename, "zero.mzML");
    assert_eq!(sparse.column_headers[&5].filename, "five.mzML");
}

/// Native: the `operator<<` layouts of `ConsensusMap` and the consensus feature
/// blocks it prints.
#[test]
fn cm_display() {
    let mut map = ConsensusMap::new();
    let header = map.column_headers.entry(0).or_default();
    header.filename = "a.mzML".into();
    header.label = "light".into();
    header.size = 2;
    let mut feature = ConsensusFeature::new();
    feature.base.rt = 1.5;
    feature.base.mz = 100.0;
    feature.base.intensity = 7.0;
    feature.base.quality = 0.25;
    feature
        .insert(FeatureHandle::from_peak(0, Peak2D::new(1.0, 99.0, 6.0), 3))
        .unwrap();
    feature.base.metadata.insert("note".into(), "x".into());
    map.features.push(feature);

    assert_eq!(
        map.to_string(),
        "Map 0: a.mzML - light - 2\n\
         ---------- CONSENSUS ELEMENT BEGIN -----------------\n\
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
         ---------- CONSENSUS ELEMENT END ----------------- \n\n"
    );
}

// ===========================================================================
// Identification-graph surface (by reference, not embedded)
// ===========================================================================

/// Native: `getUnassignedIDMatches` on both containers, and the reference
/// translation both maps expose instead of the source's embedded merge.
#[test]
fn unassigned_id_matches_and_reference_translation() {
    use openms::identification::graph::{
        IdentificationData, IdentifiedPeptide, InputFile, Observation, ObservationMatch,
        ProcessingSoftware, ProcessingStep,
    };

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

    let mut feature = Feature::default();
    feature.base.add_id_match(first).unwrap();
    let mut map = FeatureMap::from_features(vec![feature]);
    let unassigned = map.unassigned_id_matches(&graph).unwrap();
    assert_eq!(unassigned.len(), 1);
    assert!(unassigned.contains(&second));

    let mut consensus_feature = ConsensusFeature::new();
    consensus_feature.base.add_id_match(second).unwrap();
    let consensus = ConsensusMap::from_features(vec![consensus_feature]);
    let unassigned = consensus.unassigned_id_matches(&graph).unwrap();
    assert_eq!(unassigned.len(), 1);
    assert!(unassigned.contains(&first));

    // Copying the graph renumbers the references; both maps repoint atomically.
    let (copy, translator) = graph.try_clone_with_translation().unwrap();
    map.update_id_references(&translator).unwrap();
    assert_eq!(map.unassigned_id_matches(&copy).unwrap().len(), 1);

    let mut broken =
        ConsensusMap::from_features(vec![ConsensusFeature::new(), ConsensusFeature::new()]);
    broken.features[0].base.add_id_match(first).unwrap();
    // The second feature's reference belongs to the already translated graph,
    // so it has no entry in this translator and the whole map stays untouched.
    let before = broken.clone();
    broken.features[1]
        .base
        .add_id_match(
            *copy
                .observation_matches()
                .map(|(id, _)| id)
                .collect::<Vec<_>>()
                .first()
                .unwrap(),
        )
        .unwrap();
    let before_second = broken.clone();
    assert!(broken.update_id_references(&translator).is_err());
    assert_eq!(broken, before_second);
    assert_ne!(before, before_second);
}

// ===========================================================================
// Bounded work
// ===========================================================================

/// Native: the container ceilings are checked before anything is allocated.
#[test]
fn ceilings_are_declared_and_checked() {
    assert_eq!(FeatureMap::MAX_ITEMS, 10_000_000);
    assert_eq!(ConsensusMap::MAX_ITEMS, 10_000_000);
    assert_eq!(ConsensusMap::MAX_COLUMNS, 1_000_000);
    assert!(ConsensusMap::with_size(ConsensusMap::MAX_ITEMS + 1).is_err());

    let mut map = FeatureMap::new();
    assert!(map.set_primary_ms_run_path(&[]).is_ok());
    // A spectrum-free experiment still yields the fallback, not a panic.
    let experiment = MSExperiment {
        spectra: vec![MSSpectrum {
            peaks: vec![Peak1D::new(100.0, 1.0)],
            ..MSSpectrum::default()
        }],
        ..MSExperiment::default()
    };
    map.set_primary_ms_run_path_from_experiment(&["x.mzML".into()], &experiment)
        .unwrap();
    assert_eq!(
        map.primary_ms_run_path().unwrap(),
        vec!["x.mzML".to_string()]
    );
}

/// Native: the loaded-file identity moves with `swap` and is reset by the two
/// append operations, as the source's `DocumentIdentifier` assignment does.
#[test]
fn document_identity_moves_and_resets() {
    let mut generator = generator();
    let mut left = FeatureMap::new();
    left.identifier = "left".into();
    left.loaded_file_path = "left.featureXML".into();
    left.loaded_file_type = FileType::FeatureXml;
    let mut right = FeatureMap::new();
    right.identifier = "right".into();

    let mut empty = FeatureMap::new();
    left.swap(&mut empty);
    assert_eq!(empty.loaded_file_type, FileType::FeatureXml);
    assert_eq!(left.loaded_file_type, FileType::Unknown);
    empty.swap(&mut left);

    left.append(&right, &mut generator).unwrap();
    assert_eq!(left.identifier, "");
    assert_eq!(left.loaded_file_path, "");
    assert_eq!(left.loaded_file_type, FileType::Unknown);

    let mut consensus = ConsensusMap::new();
    consensus.identifier = "consensus".into();
    consensus.loaded_file_type = FileType::ConsensusXml;
    consensus
        .append_rows(&ConsensusMap::new(), &mut generator)
        .unwrap();
    assert_eq!(consensus.identifier, "");
    assert_eq!(consensus.loaded_file_type, FileType::Unknown);
}
