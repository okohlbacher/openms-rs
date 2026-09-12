// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `MRMFeature_test.cpp` (16 sections) and `MRMTransitionGroup_test.cpp`
//! (28 sections), plus native checks for the invariants the source only asserts
//! in debug builds. Every literal below is transcribed from the pinned class
//! tests at `bc9cc12514c768385ce121d6ca4bb710fe1983c4` (evidence tier 3: source
//! review). No C++ was built or executed.
//!
//! The C++ tests instantiate the group as
//! `MRMTransitionGroup<MSChromatogram, ReactionMonitoringTransition>`.
//! `ReactionMonitoringTransition` is not ported, so the transition is
//! `SimpleTransition` here; the C++ `setLibraryIntensity` / `setNativeID` calls
//! map onto its fields one for one, and the `setMetaValue("detecting_transition",
//! ...)` calls map onto its typed `detecting` flag. Neither `subset` nor
//! `getLibraryIntensity` reads that flag, so the asserted values are unaffected.

use openms::Error;
use openms::kernel::mrm::{
    DuplicateKeyPolicy, MRMFeature, MRMTransitionGroup, OpenSwathIndScores, OpenSwathScores,
    SimpleTransition,
};
use openms::kernel::{Feature, MSChromatogram};
use openms::metadata::MetaValue;

type Group = MRMTransitionGroup<MSChromatogram, SimpleTransition>;

fn feature_with_dummy() -> Feature {
    let mut feature = Feature::default();
    feature
        .base
        .metadata
        .insert("dummy".into(), MetaValue::from(1i64));
    feature
}

fn chromatogram(native_id: &str) -> MSChromatogram {
    MSChromatogram {
        native_id: native_id.into(),
        ..MSChromatogram::default()
    }
}

fn chromatogram_with_value() -> MSChromatogram {
    let mut chrom = MSChromatogram::default();
    chrom
        .metadata
        .insert("some_value".into(), MetaValue::from(1i64));
    chrom
}

fn ids(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| (*value).to_owned()).collect()
}

// ---------------------------------------------------------------------------
// MRMFeature_test.cpp
// ---------------------------------------------------------------------------

/// START_SECTION(MRMFeature())
#[test]
fn mrm_feature_default_constructor() {
    let feature = MRMFeature::new();
    assert_eq!(feature, MRMFeature::default());
    assert!(feature.features().is_empty());
    assert!(feature.precursor_features().is_empty());
    assert_eq!(feature.scores(), &OpenSwathScores::default());
}

/// START_SECTION(~MRMFeature())
///
/// The source section only deletes a heap-allocated feature. Rust drops the
/// owned lists, maps and metadata without an explicit destructor; the check is
/// that a populated value can be dropped and that nothing it owned is shared.
#[test]
fn mrm_feature_destructor() {
    let mut feature = MRMFeature::new();
    feature
        .add_feature(feature_with_dummy(), "chromatogram1")
        .unwrap();
    feature
        .add_precursor_feature(feature_with_dummy(), "precursor1")
        .unwrap();
    let boxed = Box::new(feature);
    drop(boxed);
}

/// START_SECTION(MRMFeature(const MRMFeature &rhs))
#[test]
fn mrm_feature_copy_constructor() {
    let mut tmp = MRMFeature::new();
    tmp.base.intensity = 100.0;
    tmp.add_score("testscore", 200.0).unwrap();

    let tmp2 = tmp.clone();

    assert_eq!(tmp2.score("testscore"), Some(200.0));
    assert_eq!(tmp2.base.intensity, 100.0);
}

/// START_SECTION((MRMFeature(const MRMFeature&& source)))
///
/// The source asserts that the move constructor is `noexcept`, so that
/// `std::vector` moves instead of copying. Rust moves are unconditional and
/// cannot unwind; the observable equivalent is that a move transfers the whole
/// state unchanged.
#[test]
fn mrm_feature_move_constructor() {
    let mut tmp = MRMFeature::new();
    tmp.base.intensity = 100.0;
    tmp.add_score("testscore", 200.0).unwrap();
    tmp.add_feature(feature_with_dummy(), "chromatogram1")
        .unwrap();
    let expected = tmp.clone();

    let moved = tmp;

    assert_eq!(moved, expected);
    assert_eq!(moved.score("testscore"), Some(200.0));
    assert_eq!(moved.features().len(), 1);
}

/// START_SECTION(MRMFeature& operator=(const MRMFeature &rhs))
#[test]
fn mrm_feature_assignment_operator() {
    let mut tmp = MRMFeature::new();
    tmp.base.intensity = 100.0;
    tmp.add_score("testscore", 200.0).unwrap();

    let mut tmp2 = MRMFeature::new();
    assert_eq!(tmp2.score("testscore"), None);
    tmp2 = tmp.clone();

    assert_eq!(tmp2.score("testscore"), Some(200.0));
    assert_eq!(tmp2.base.intensity, 100.0);
}

/// START_SECTION(const PGScoresType & getScores() const)
///
/// NOT_TESTABLE in the source, which says it is covered by set/add score. The
/// accessor is exercised here directly: the record round-trips through
/// `set_scores` and the mutable accessor.
#[test]
fn mrm_feature_get_scores() {
    let mut feature = MRMFeature::new();
    assert_eq!(feature.scores().library_sangle, 0.0);
    assert_eq!(feature.scores().ms1_mi_score, -1.0);

    feature.scores_mut().library_sangle = 99.0;
    assert_eq!(feature.scores().library_sangle, 99.0);
}

/// START_SECTION(double getScore(const std::string & score_name))
///
/// NOT_TESTABLE in the source, and the member no longer exists in the pinned
/// header; the native `score` accessor replaces it by reading the metadata that
/// `addScore` writes.
#[test]
fn mrm_feature_get_score() {
    let mut feature = MRMFeature::new();
    feature.add_score("score1", 1.0).unwrap();
    assert_eq!(feature.score("score1"), Some(1.0));
    assert_eq!(feature.score("absent"), None);
}

/// START_SECTION(Feature & getFeature(std::string key))
#[test]
fn mrm_feature_get_feature() {
    let mut mrmfeature = MRMFeature::new();
    let f1 = feature_with_dummy();
    mrmfeature.add_feature(f1.clone(), "chromatogram1").unwrap();
    mrmfeature.add_feature(f1, "chromatogram2").unwrap();
    assert_eq!(
        mrmfeature
            .feature("chromatogram1")
            .unwrap()
            .base
            .metadata
            .get("dummy")
            .unwrap()
            .as_i64()
            .unwrap(),
        1
    );
}

/// START_SECTION(void setScores(const PGScoresType & scores))
#[test]
fn mrm_feature_set_scores() {
    let mut mrmfeature = MRMFeature::new();
    let scores = OpenSwathScores {
        library_sangle: 99.0,
        ..OpenSwathScores::default()
    };
    mrmfeature.set_scores(scores);

    assert_eq!(scores.library_sangle, mrmfeature.scores().library_sangle);
}

/// START_SECTION(void addScore(const std::string & score_name, double score))
#[test]
fn mrm_feature_add_score() {
    let mut mrmfeature = MRMFeature::new();
    mrmfeature.add_score("score1", 1.0).unwrap();
    mrmfeature.add_score("score2", 2.0).unwrap();
    assert_eq!(mrmfeature.score("score1"), Some(1.0));
    assert_eq!(mrmfeature.score("score2"), Some(2.0));
}

/// START_SECTION(void addFeature(Feature & feature, const std::string & key))
///
/// NOT_TESTABLE in the source ("tested in getFeature"). Exercised here for the
/// list and key bookkeeping the source leaves implicit.
#[test]
fn mrm_feature_add_feature() {
    let mut mrmfeature = MRMFeature::new();
    mrmfeature
        .add_feature(feature_with_dummy(), "chromatogram1")
        .unwrap();
    assert!(mrmfeature.has_feature("chromatogram1"));
    assert!(!mrmfeature.has_feature("chromatogram2"));
    assert_eq!(mrmfeature.features().len(), 1);
}

/// START_SECTION(const std::vector<Feature> & getFeatures() const)
#[test]
fn mrm_feature_get_features() {
    let mut mrmfeature = MRMFeature::new();
    let f1 = feature_with_dummy();
    mrmfeature.add_feature(f1.clone(), "chromatogram1").unwrap();
    mrmfeature.add_feature(f1, "chromatogram2").unwrap();
    assert_eq!(mrmfeature.features().len(), 2);
}

/// START_SECTION(void getFeatureIDs(std::vector<std::string> & result) const)
#[test]
fn mrm_feature_get_feature_ids() {
    let mut mrmfeature = MRMFeature::new();
    let f1 = feature_with_dummy();
    mrmfeature.add_feature(f1.clone(), "chromatogram1").unwrap();
    mrmfeature.add_feature(f1, "chromatogram2").unwrap();
    let result: Vec<&str> = mrmfeature.feature_ids().collect();
    assert_eq!(result.len(), 2);
    assert_eq!(result[0], "chromatogram1");
    assert_eq!(result[1], "chromatogram2");
}

/// START_SECTION(void addPrecursorFeature(Feature & feature, const std::string & key))
#[test]
fn mrm_feature_add_precursor_feature() {
    // Initially, there should be no feature present
    let mut mrmfeature = MRMFeature::new();
    assert_eq!(mrmfeature.precursor_feature_ids().count(), 0);

    // After adding a feature, there should be one feature present
    let f1 = Feature::default();
    mrmfeature
        .add_precursor_feature(f1, "precursor_chromatogram1")
        .unwrap();
    assert_eq!(mrmfeature.precursor_feature_ids().count(), 1);
}

/// START_SECTION(void getPrecursorFeatureIDs(std::vector<std::string> & result) const)
#[test]
fn mrm_feature_get_precursor_feature_ids() {
    let mut mrmfeature = MRMFeature::new();
    let f1 = feature_with_dummy();
    mrmfeature
        .add_precursor_feature(f1.clone(), "chromatogram1")
        .unwrap();
    mrmfeature
        .add_precursor_feature(f1, "chromatogram2")
        .unwrap();
    let result: Vec<&str> = mrmfeature.precursor_feature_ids().collect();
    assert_eq!(result.len(), 2);
    assert_eq!(result[0], "chromatogram1");
    assert_eq!(result[1], "chromatogram2");
}

/// START_SECTION(Feature & getPrecursorFeature(std::string key))
#[test]
fn mrm_feature_get_precursor_feature() {
    let mut mrmfeature = MRMFeature::new();
    let f1 = feature_with_dummy();
    mrmfeature
        .add_precursor_feature(f1.clone(), "chromatogram1")
        .unwrap();
    mrmfeature
        .add_precursor_feature(f1, "chromatogram2")
        .unwrap();
    assert_eq!(
        mrmfeature
            .precursor_feature("chromatogram1")
            .unwrap()
            .base
            .metadata
            .get("dummy")
            .unwrap()
            .as_i64()
            .unwrap(),
        1
    );
}

// ---------------------------------------------------------------------------
// MRMTransitionGroup_test.cpp
// ---------------------------------------------------------------------------

/// START_SECTION(MRMTransitionGroup())
#[test]
fn group_default_constructor() {
    let group = Group::new();
    assert_eq!(group.size(), 0);
    assert!(group.is_empty());
    assert_eq!(group.transition_group_id(), "");
    assert!(group.transitions().is_empty());
    assert!(group.features().is_empty());
}

/// START_SECTION(~MRMTransitionGroup())
#[test]
fn group_destructor() {
    let mut group = Group::new();
    group
        .add_chromatogram(chromatogram("dummy1"), "dummy1")
        .unwrap();
    let boxed = Box::new(group);
    drop(boxed);
}

/// START_SECTION(MRMTransitionGroup(const MRMTransitionGroup &rhs))
#[test]
fn group_copy_constructor() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_chromatogram(MSChromatogram::default(), "dummy1")
        .unwrap();
    mrmtrgroup
        .add_chromatogram(MSChromatogram::default(), "dummy2")
        .unwrap();

    let tmp = mrmtrgroup.clone();
    assert_eq!(mrmtrgroup.size(), tmp.size());
}

/// START_SECTION(MRMTransitionGroup& operator=(const MRMTransitionGroup &rhs))
#[test]
fn group_assignment_operator() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_chromatogram(MSChromatogram::default(), "dummy1")
        .unwrap();
    mrmtrgroup
        .add_chromatogram(MSChromatogram::default(), "dummy2")
        .unwrap();

    let mut tmp = Group::new();
    assert_eq!(tmp.size(), 0);
    tmp = mrmtrgroup.clone();
    assert_eq!(mrmtrgroup.size(), tmp.size());
}

/// START_SECTION(Size size() const)
#[test]
fn group_size() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_chromatogram(MSChromatogram::default(), "dummy1")
        .unwrap();
    assert_eq!(mrmtrgroup.size(), 1);
    mrmtrgroup
        .add_chromatogram(MSChromatogram::default(), "dummy2")
        .unwrap();
    assert_eq!(mrmtrgroup.size(), 2);
}

/// START_SECTION(const std::string & getTransitionGroupID() const)
#[test]
fn group_get_transition_group_id() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup.set_transition_group_id("some_id");
    assert_eq!(mrmtrgroup.transition_group_id(), "some_id");
}

/// START_SECTION(void setTransitionGroupID(const std::string & tr_gr_id))
///
/// NOT_TESTABLE in the source ("tested above"). Exercised here for the
/// overwrite the source leaves implicit.
#[test]
fn group_set_transition_group_id() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup.set_transition_group_id("some_id");
    mrmtrgroup.set_transition_group_id("other_id");
    assert_eq!(mrmtrgroup.transition_group_id(), "other_id");
}

/// START_SECTION(std::vector<TransitionType>& getTransitionsMuteable())
#[test]
fn group_get_transitions_muteable() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_transition(SimpleTransition::default(), "dummy1")
        .unwrap();
    mrmtrgroup
        .add_transition(SimpleTransition::default(), "dummy2")
        .unwrap();
    assert_eq!(mrmtrgroup.transitions_mut().len(), 2);

    mrmtrgroup.transitions_mut()[1].library_intensity = 7.0;
    assert_eq!(mrmtrgroup.transitions()[1].library_intensity, 7.0);
}

/// START_SECTION(void addTransition(const TransitionType &transition, std::string key))
///
/// NOT_TESTABLE in the source ("tested above"). Exercised here for the
/// duplicate-key rejection the source states as a throw.
#[test]
fn group_add_transition() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_transition(SimpleTransition::default(), "dummy1")
        .unwrap();
    assert_eq!(mrmtrgroup.transitions().len(), 1);
    let error = mrmtrgroup
        .add_transition(SimpleTransition::default(), "dummy1")
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)));
    assert_eq!(mrmtrgroup.transitions().len(), 1);
}

/// START_SECTION(const TransitionType& getTransition(std::string key))
#[test]
fn group_get_transition() {
    let mut mrmtrgroup = Group::new();
    let trans1 = SimpleTransition::new("", 42.0);
    mrmtrgroup.add_transition(trans1, "dummy1").unwrap();
    assert_eq!(
        mrmtrgroup.transition("dummy1").unwrap().library_intensity,
        42.0
    );
}

/// START_SECTION(const std::vector<TransitionType>& getTransitions() const)
#[test]
fn group_get_transitions() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_transition(SimpleTransition::new("", 42.0), "dummy1")
        .unwrap();
    mrmtrgroup
        .add_transition(SimpleTransition::new("", -2.0), "dummy2")
        .unwrap();
    assert_eq!(mrmtrgroup.transitions()[0].library_intensity, 42.0);
    assert_eq!(mrmtrgroup.transitions()[1].library_intensity, -2.0);
}

/// START_SECTION(bool hasTransition(std::string key))
#[test]
fn group_has_transition() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_transition(SimpleTransition::default(), "dummy1")
        .unwrap();
    assert!(mrmtrgroup.has_transition("dummy1"));
    assert!(!mrmtrgroup.has_transition("dummy2"));
}

/// START_SECTION(const std::vector<SpectrumType>& getChromatograms() const)
#[test]
fn group_get_chromatograms_const() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_chromatogram(MSChromatogram::default(), "dummy1")
        .unwrap();
    mrmtrgroup
        .add_chromatogram(MSChromatogram::default(), "dummy2")
        .unwrap();
    assert_eq!(mrmtrgroup.chromatograms().len(), 2);
}

/// START_SECTION(std::vector<SpectrumType>& getChromatograms())
#[test]
fn group_get_chromatograms_mutable() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_chromatogram(MSChromatogram::default(), "dummy1")
        .unwrap();
    mrmtrgroup
        .add_chromatogram(MSChromatogram::default(), "dummy2")
        .unwrap();
    assert_eq!(mrmtrgroup.chromatograms_mut().len(), 2);
}

/// START_SECTION(void addChromatogram(SpectrumType &chromatogram, std::string key))
///
/// NOT_TESTABLE in the source ("tested above"). Exercised here for the
/// duplicate-key rejection the source states as a throw.
#[test]
fn group_add_chromatogram() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_chromatogram(MSChromatogram::default(), "dummy1")
        .unwrap();
    let error = mrmtrgroup
        .add_chromatogram(MSChromatogram::default(), "dummy1")
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)));
    assert_eq!(mrmtrgroup.size(), 1);
}

/// START_SECTION(SpectrumType& getChromatogram(std::string key))
#[test]
fn group_get_chromatogram() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_chromatogram(chromatogram_with_value(), "dummy1")
        .unwrap();
    assert_eq!(
        mrmtrgroup
            .chromatogram("dummy1")
            .unwrap()
            .metadata
            .get("some_value")
            .unwrap()
            .as_i64()
            .unwrap(),
        1
    );
}

/// START_SECTION(bool hasChromatogram(std::string key))
#[test]
fn group_has_chromatogram() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_chromatogram(MSChromatogram::default(), "dummy1")
        .unwrap();
    assert!(mrmtrgroup.has_chromatogram("dummy1"));
    assert!(!mrmtrgroup.has_chromatogram("dummy2"));
}

/// START_SECTION(void addPrecusorChromatogram(SpectrumType &chromatogram, std::string key))
///
/// NOT_TESTABLE in the source ("tested below"). Exercised here for the
/// duplicate-key rejection and for the independence of the two key maps.
#[test]
fn group_add_precursor_chromatogram() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_precursor_chromatogram(MSChromatogram::default(), "dummy1")
        .unwrap();
    // The same key is free in the fragment-ion map.
    mrmtrgroup
        .add_chromatogram(MSChromatogram::default(), "dummy1")
        .unwrap();
    let error = mrmtrgroup
        .add_precursor_chromatogram(MSChromatogram::default(), "dummy1")
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)));
    assert_eq!(mrmtrgroup.precursor_chromatograms().len(), 1);
}

/// START_SECTION(SpectrumType& getPrecursorChromatogram(std::string key))
#[test]
fn group_get_precursor_chromatogram() {
    let mut mrmtrgroup = Group::new();
    let chrom1 = chromatogram_with_value();
    mrmtrgroup
        .add_precursor_chromatogram(chrom1.clone(), "dummy1")
        .unwrap();
    assert_eq!(
        mrmtrgroup
            .precursor_chromatogram("dummy1")
            .unwrap()
            .metadata
            .get("some_value")
            .unwrap()
            .as_i64()
            .unwrap(),
        1
    );

    // Add a few feature chromatograms and then add a precursor chromatogram ->
    // it should still work
    mrmtrgroup
        .add_chromatogram(chrom1.clone(), "feature1")
        .unwrap();
    mrmtrgroup
        .add_chromatogram(chrom1.clone(), "feature2")
        .unwrap();
    mrmtrgroup
        .add_chromatogram(chrom1.clone(), "feature3")
        .unwrap();
    mrmtrgroup
        .add_precursor_chromatogram(chrom1, "dummy2")
        .unwrap();
    assert_eq!(
        mrmtrgroup
            .precursor_chromatogram("dummy2")
            .unwrap()
            .metadata
            .get("some_value")
            .unwrap()
            .as_i64()
            .unwrap(),
        1
    );
}

/// START_SECTION(bool hasPrecursorChromatogram(std::string key))
#[test]
fn group_has_precursor_chromatogram() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_precursor_chromatogram(MSChromatogram::default(), "dummy1")
        .unwrap();
    assert!(mrmtrgroup.has_precursor_chromatogram("dummy1"));
    assert!(!mrmtrgroup.has_precursor_chromatogram("dummy2"));
}

/// START_SECTION(const std::vector<MRMFeature> & getFeatures() const)
#[test]
fn group_get_features() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup.add_feature(MRMFeature::default()).unwrap();
    mrmtrgroup.add_feature(MRMFeature::default()).unwrap();
    assert_eq!(mrmtrgroup.features().len(), 2);
}

/// START_SECTION(std::vector<MRMFeature> & getFeaturesMuteable())
#[test]
fn group_get_features_muteable() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup.add_feature(MRMFeature::default()).unwrap();
    mrmtrgroup.add_feature(MRMFeature::default()).unwrap();
    assert_eq!(mrmtrgroup.features_mut().len(), 2);

    mrmtrgroup.features_mut()[0].base.intensity = 5.0;
    assert_eq!(mrmtrgroup.features()[0].base.intensity, 5.0);
}

/// START_SECTION(void addFeature(MRMFeature & feature))
///
/// NOT_TESTABLE in the source ("tested above"). Exercised here for insertion
/// order, which `getBestFeature` depends on.
#[test]
fn group_add_feature() {
    let mut mrmtrgroup = Group::new();
    let mut first = MRMFeature::default();
    first.base.intensity = 1.0;
    let mut second = MRMFeature::default();
    second.base.intensity = 2.0;
    mrmtrgroup.add_feature(first).unwrap();
    mrmtrgroup.add_feature(second).unwrap();
    assert_eq!(mrmtrgroup.features()[0].base.intensity, 1.0);
    assert_eq!(mrmtrgroup.features()[1].base.intensity, 2.0);
}

/// START_SECTION(void getLibraryIntensity(std::vector<double> & result) const)
#[test]
fn group_get_library_intensity() {
    let mut mrmtrgroup = Group::new();
    mrmtrgroup
        .add_transition(SimpleTransition::new("", 3.0), "dummy1")
        .unwrap();
    mrmtrgroup
        .add_transition(SimpleTransition::new("", -2.0), "dummy2")
        .unwrap();
    let result = mrmtrgroup.library_intensity().unwrap();
    assert_eq!(result.len(), 2);
    assert_eq!(result[0], 3.0);
    assert_eq!(result[1], 0.0);
}

/// START_SECTION(MRMTransitionGroup subset(std::vector<std::string> tr_ids))
#[test]
fn group_subset() {
    let mut new_trans1 = SimpleTransition::new("new_trans1", 3.0);
    new_trans1.detecting = true;
    let mut new_trans2 = SimpleTransition::new("new_trans2", -2.0);
    new_trans2.detecting = false;
    let mut mrmtrgroup = Group::new();
    mrmtrgroup.add_transition(new_trans1, "new_trans1").unwrap();
    mrmtrgroup.add_transition(new_trans2, "new_trans2").unwrap();
    let transition_ids = ids(&["new_trans1"]);

    let mrmtrgroupsub = mrmtrgroup.subset(&transition_ids).unwrap();
    let result = mrmtrgroupsub.library_intensity().unwrap();
    assert_eq!(result.len(), 1);
    assert_eq!(result[0], 3.0);
}

/// START_SECTION(inline bool isInternallyConsistent() const)
#[test]
fn group_is_internally_consistent() {
    let mrmtrgroup = Group::new();
    assert!(mrmtrgroup.is_internally_consistent());
}

/// START_SECTION(inline bool chromatogramIdsMatch() const)
#[test]
fn group_chromatogram_ids_match() {
    {
        let mut mrmtrgroup = Group::new();
        let c = chromatogram("test");
        mrmtrgroup.add_chromatogram(c.clone(), "test").unwrap();

        assert!(mrmtrgroup.chromatogram_ids_match());
        mrmtrgroup.add_chromatogram(c, "test2").unwrap();
        assert!(!mrmtrgroup.chromatogram_ids_match());
    }

    {
        let mut mrmtrgroup = Group::new();
        let c = chromatogram("test");
        mrmtrgroup
            .add_precursor_chromatogram(c.clone(), "test")
            .unwrap();

        assert!(mrmtrgroup.chromatogram_ids_match());
        mrmtrgroup.add_precursor_chromatogram(c, "test2").unwrap();
        assert!(!mrmtrgroup.chromatogram_ids_match());
    }
}

/// START_SECTION(MRMTransitionGroup subsetDependent(std::vector<std::string> tr_ids))
///
/// The source section is named for `subsetDependent` but calls `subset`; the
/// port keeps the call the section actually makes, so the transcribed literals
/// stay meaningful. `subset_dependent` itself is covered by
/// `subset_dependent_keeps_whole_features_and_drops_precursors` below.
#[test]
fn group_subset_dependent_section() {
    let mut new_trans1 = SimpleTransition::new("new_trans1", 3.0);
    new_trans1.detecting = true;
    let mut new_trans2 = SimpleTransition::new("new_trans2", -2.0);
    new_trans2.detecting = false;
    let mut mrmtrgroup = Group::new();
    mrmtrgroup.add_transition(new_trans1, "new_trans1").unwrap();
    mrmtrgroup.add_transition(new_trans2, "new_trans2").unwrap();
    let transition_ids = ids(&["new_trans1", "new_trans2"]);

    let mrmtrgroupsub = mrmtrgroup.subset(&transition_ids).unwrap();
    let result = mrmtrgroupsub.library_intensity().unwrap();
    assert_eq!(result.len(), 2);
    assert_eq!(result[0], 3.0);
    assert_eq!(result[1], 0.0);
}

// ---------------------------------------------------------------------------
// Native checks: the invariants the source only asserts in debug builds, the
// error paths its preconditions compile out, and the port's own ceilings.
// ---------------------------------------------------------------------------

/// `isInternallyConsistent` reports the three source conditions rather than
/// always returning true, which is what the release-mode C++ does.
#[test]
fn internal_consistency_detects_every_source_condition() {
    // Unequal list lengths.
    let mut group = Group::new();
    group
        .add_transition(SimpleTransition::new("a", 1.0), "a")
        .unwrap();
    assert!(!group.is_internally_consistent());
    group.add_chromatogram(chromatogram("a"), "a").unwrap();
    assert!(group.is_internally_consistent());

    // Equal lengths, but a chromatogram key that names no transition.
    let mut skewed = Group::new();
    skewed
        .add_transition(SimpleTransition::new("a", 1.0), "a")
        .unwrap();
    skewed.add_chromatogram(chromatogram("b"), "b").unwrap();
    assert_eq!(skewed.transitions().len(), skewed.chromatograms().len());
    assert!(!skewed.is_internally_consistent());
}

/// A precursor chromatogram whose key matches while a fragment one does not is
/// still reported as a mismatch: the fragment map is checked first.
#[test]
fn chromatogram_ids_match_checks_both_maps() {
    let mut group = Group::new();
    group
        .add_chromatogram(chromatogram("wrong"), "key")
        .unwrap();
    group
        .add_precursor_chromatogram(chromatogram("p"), "p")
        .unwrap();
    assert!(!group.chromatogram_ids_match());
}

/// Lookups the source guards only with `OPENMS_PRECONDITION` become checked
/// errors here instead of returning the first element of the list.
#[test]
fn unknown_keys_are_missing_information() {
    let mut group = Group::new();
    group
        .add_transition(SimpleTransition::new("a", 1.0), "a")
        .unwrap();
    group.add_chromatogram(chromatogram("a"), "a").unwrap();
    group
        .add_precursor_chromatogram(chromatogram("p"), "p")
        .unwrap();

    assert!(matches!(
        group.transition("zzz").unwrap_err(),
        Error::MissingInformation(_)
    ));
    assert!(matches!(
        group.chromatogram("zzz").unwrap_err(),
        Error::MissingInformation(_)
    ));
    assert!(matches!(
        group.precursor_chromatogram("zzz").unwrap_err(),
        Error::MissingInformation(_)
    ));

    let mut feature = MRMFeature::new();
    feature.add_feature(feature_with_dummy(), "a").unwrap();
    assert!(matches!(
        feature.feature("zzz").unwrap_err(),
        Error::MissingInformation(_)
    ));
    assert!(matches!(
        feature.precursor_feature("zzz").unwrap_err(),
        Error::MissingInformation(_)
    ));
}

/// The source `addFeature` overwrite strands the previously keyed feature; the
/// native default refuses, the explicit policy reproduces it.
#[test]
fn duplicate_feature_keys_reject_by_default() {
    let mut feature = MRMFeature::new();
    feature
        .add_feature(Feature::new(1.0, 2.0, 3.0), "k")
        .unwrap();
    let error = feature
        .add_feature(Feature::new(4.0, 5.0, 6.0), "k")
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)));
    assert_eq!(feature.features().len(), 1);
    assert_eq!(feature.feature("k").unwrap().base.rt, 1.0);

    let mut source_like = MRMFeature::new();
    source_like
        .add_feature_with(
            Feature::new(1.0, 2.0, 3.0),
            "k",
            DuplicateKeyPolicy::SourceOverwrite,
        )
        .unwrap();
    source_like
        .add_feature_with(
            Feature::new(4.0, 5.0, 6.0),
            "k",
            DuplicateKeyPolicy::SourceOverwrite,
        )
        .unwrap();
    // Both features are stored, only the second is reachable by key.
    assert_eq!(source_like.features().len(), 2);
    assert_eq!(source_like.feature_ids().count(), 1);
    assert_eq!(source_like.feature("k").unwrap().base.rt, 4.0);

    let mut precursors = MRMFeature::new();
    precursors
        .add_precursor_feature(Feature::default(), "k")
        .unwrap();
    assert!(matches!(
        precursors
            .add_precursor_feature(Feature::default(), "k")
            .unwrap_err(),
        Error::InvalidValue(_)
    ));
    precursors
        .add_precursor_feature_with(
            Feature::new(9.0, 9.0, 9.0),
            "k",
            DuplicateKeyPolicy::SourceOverwrite,
        )
        .unwrap();
    assert_eq!(precursors.precursor_features().len(), 2);
    assert_eq!(precursors.precursor_feature("k").unwrap().base.rt, 9.0);
}

/// A non-finite score is refused rather than stored, unlike the source.
#[test]
fn non_finite_scores_are_refused() {
    let mut feature = MRMFeature::new();
    assert!(matches!(
        feature.add_score("bad", f64::NAN).unwrap_err(),
        Error::InvalidValue(_)
    ));
    assert!(feature.base.metadata.is_empty());

    let scores = OpenSwathIndScores {
        ind_area_intensity: vec![1.0, f64::INFINITY],
        ..OpenSwathIndScores::default()
    };
    assert!(matches!(
        feature.id_scores_as_meta_value(false, &scores).unwrap_err(),
        Error::InvalidValue(_)
    ));
    // The whole block is committed at once, so nothing was written.
    assert!(feature.base.metadata.is_empty());
}

/// `IDScoresAsMetaValue` writes the whole 42-key block under the prefix chosen
/// by `decoy`, including empty lists.
#[test]
fn id_scores_as_meta_value_writes_the_whole_block() {
    let scores = OpenSwathIndScores {
        ind_num_transitions: 2,
        ind_transition_names: vec!["y7".into(), "b3".into()],
        ind_area_intensity: vec![10.0, 20.0],
        ind_fwhm: vec![1.5, 2.5],
        ind_points_across_half_height: vec![7.0],
        ..OpenSwathIndScores::default()
    };

    let mut target = MRMFeature::new();
    target.id_scores_as_meta_value(false, &scores).unwrap();
    assert_eq!(target.base.metadata.len(), 42);
    for suffix in OpenSwathIndScores::KEY_SUFFIXES {
        assert!(
            target
                .base
                .metadata
                .contains_key(&format!("id_target_{suffix}")),
            "missing id_target_{suffix}"
        );
    }
    assert_eq!(
        target
            .base
            .metadata
            .get("id_target_num_transitions")
            .unwrap()
            .as_i64()
            .unwrap(),
        2
    );
    assert_eq!(
        target
            .base
            .metadata
            .get("id_target_transition_names")
            .unwrap()
            .as_string_list()
            .unwrap(),
        ["y7".to_owned(), "b3".to_owned()]
    );
    assert_eq!(
        target
            .base
            .metadata
            .get("id_target_area_intensity")
            .unwrap()
            .as_float_list()
            .unwrap(),
        [10.0, 20.0]
    );
    // `ind_fwhm` is written under `width_at_50`, not under its field name.
    assert_eq!(
        target
            .base
            .metadata
            .get("id_target_width_at_50")
            .unwrap()
            .as_float_list()
            .unwrap(),
        [1.5, 2.5]
    );
    // An untouched list still gets its key, with no values.
    assert!(
        target
            .base
            .metadata
            .get("id_target_ind_total_width")
            .unwrap()
            .as_float_list()
            .unwrap()
            .is_empty()
    );

    let mut decoy = MRMFeature::new();
    decoy.id_scores_as_meta_value(true, &scores).unwrap();
    assert_eq!(decoy.base.metadata.len(), 42);
    assert!(decoy.base.metadata.contains_key("id_decoy_width_at_50"));
    assert!(!decoy.base.metadata.contains_key("id_target_width_at_50"));
}

/// The best feature is the first of the maxima, by `Feature::quality`.
#[test]
fn best_feature_takes_the_first_maximum() {
    let mut group = Group::new();
    assert!(matches!(
        group.best_feature().unwrap_err(),
        Error::MissingInformation(_)
    ));

    for (index, quality) in [0.5_f32, 0.9, 0.9, 0.2].iter().enumerate() {
        let mut feature = MRMFeature::default();
        feature.base.quality = *quality;
        feature.base.rt = index as f64;
        group.add_feature(feature).unwrap();
    }
    assert_eq!(group.best_feature().unwrap().base.rt, 1.0);
}

/// `subset` carries every precursor chromatogram, re-keyed by its own native
/// ID, and rebuilds each peak group from intensity, RT and metadata only.
#[test]
fn subset_rebuilds_features_and_rekeys_precursors() {
    let mut group = Group::new();
    group.set_transition_group_id("group_1");
    group
        .add_transition(SimpleTransition::new("t1", 5.0), "t1")
        .unwrap();
    group
        .add_transition(SimpleTransition::new("t2", 6.0), "t2")
        .unwrap();
    group.add_chromatogram(chromatogram("t1"), "t1").unwrap();
    group.add_chromatogram(chromatogram("t2"), "t2").unwrap();
    // Stored under a key that is not the chromatogram's native ID.
    group
        .add_precursor_chromatogram(chromatogram("ms1"), "stored_under_this")
        .unwrap();

    let mut feature = MRMFeature::default();
    feature.base.rt = 33.0;
    feature.base.intensity = 44.0;
    feature.base.quality = 0.75;
    feature.base.charge = 2;
    feature
        .base
        .metadata
        .insert("label".into(), MetaValue::from("kept".to_owned()));
    feature
        .add_feature(Feature::new(1.0, 1.0, 1.0), "t1")
        .unwrap();
    feature
        .add_feature(Feature::new(2.0, 2.0, 2.0), "t2")
        .unwrap();
    feature
        .add_precursor_feature(Feature::new(3.0, 3.0, 3.0), "ms1")
        .unwrap();
    group.add_feature(feature).unwrap();

    let subset = group.subset(&ids(&["t1"])).unwrap();

    assert_eq!(subset.transition_group_id(), "group_1");
    assert_eq!(subset.transitions().len(), 1);
    assert_eq!(subset.size(), 1);
    assert_eq!(subset.precursor_chromatograms().len(), 1);
    assert!(subset.has_precursor_chromatogram("ms1"));
    assert!(!subset.has_precursor_chromatogram("stored_under_this"));

    let rebuilt = &subset.features()[0];
    assert_eq!(rebuilt.base.rt, 33.0);
    assert_eq!(rebuilt.base.intensity, 44.0);
    assert_eq!(
        rebuilt
            .base
            .metadata
            .get("label")
            .unwrap()
            .as_str()
            .unwrap(),
        "kept"
    );
    // Quality and charge are not carried by the source's rebuild.
    assert_eq!(rebuilt.base.quality, 0.0);
    assert_eq!(rebuilt.base.charge, 0);
    assert_eq!(rebuilt.features().len(), 1);
    assert!(rebuilt.has_feature("t1"));
    assert!(!rebuilt.has_feature("t2"));
    assert_eq!(rebuilt.precursor_features().len(), 1);
    assert_eq!(rebuilt.precursor_feature("ms1").unwrap().base.rt, 3.0);
}

/// A transition registered under a key other than its native ID is dropped by
/// `subset`, while a chromatogram stored under that native ID is still carried.
#[test]
fn subset_selects_by_native_id_not_by_registration_key() {
    let mut group = Group::new();
    group
        .add_transition(SimpleTransition::new("native", 5.0), "other_key")
        .unwrap();
    group
        .add_chromatogram(chromatogram("native"), "native")
        .unwrap();

    let subset = group.subset(&ids(&["native"])).unwrap();
    assert!(subset.transitions().is_empty());
    assert_eq!(subset.size(), 1);
}

/// `subset` needs a per-transition feature for every selected transition; the
/// source throws `std::out_of_range` there.
#[test]
fn subset_requires_a_feature_per_selected_transition() {
    let mut group = Group::new();
    group
        .add_transition(SimpleTransition::new("t1", 5.0), "t1")
        .unwrap();
    group.add_chromatogram(chromatogram("t1"), "t1").unwrap();
    group.add_feature(MRMFeature::default()).unwrap();

    assert!(matches!(
        group.subset(&ids(&["t1"])).unwrap_err(),
        Error::MissingInformation(_)
    ));
}

/// `subsetDependent` copies the peak groups whole, keeps no precursor
/// chromatogram, and requires a chromatogram for each selected transition.
#[test]
fn subset_dependent_keeps_whole_features_and_drops_precursors() {
    let mut group = Group::new();
    group.set_transition_group_id("group_1");
    group
        .add_transition(SimpleTransition::new("t1", 5.0), "t1")
        .unwrap();
    group
        .add_transition(SimpleTransition::new("t2", 6.0), "t2")
        .unwrap();
    group.add_chromatogram(chromatogram("t1"), "t1").unwrap();
    group.add_chromatogram(chromatogram("t2"), "t2").unwrap();
    group
        .add_precursor_chromatogram(chromatogram("ms1"), "ms1")
        .unwrap();

    let mut feature = MRMFeature::default();
    feature.base.quality = 0.75;
    feature
        .add_feature(Feature::new(1.0, 1.0, 1.0), "t1")
        .unwrap();
    feature
        .add_feature(Feature::new(2.0, 2.0, 2.0), "t2")
        .unwrap();
    group.add_feature(feature).unwrap();

    let subset = group.subset_dependent(&ids(&["t1"])).unwrap();
    assert_eq!(subset.transition_group_id(), "group_1");
    assert_eq!(subset.transitions().len(), 1);
    assert_eq!(subset.size(), 1);
    assert!(subset.precursor_chromatograms().is_empty());
    // The feature is copied whole: the dropped transition's feature survives.
    let copied = &subset.features()[0];
    assert_eq!(copied.base.quality, 0.75);
    assert_eq!(copied.features().len(), 2);
    assert!(copied.has_feature("t2"));

    // Without a chromatogram under the native ID, the source indexes `at()`.
    let mut bare = Group::new();
    bare.add_transition(SimpleTransition::new("t1", 5.0), "t1")
        .unwrap();
    assert!(matches!(
        bare.subset_dependent(&ids(&["t1"])).unwrap_err(),
        Error::MissingInformation(_)
    ));
}

/// Both subset operations build into a temporary, so a rejected call leaves the
/// receiver untouched.
#[test]
fn subset_leaves_the_receiver_unchanged_on_error() {
    let mut group = Group::new();
    group
        .add_transition(SimpleTransition::new("t1", 5.0), "t1")
        .unwrap();
    group.add_chromatogram(chromatogram("t1"), "t1").unwrap();
    group.add_feature(MRMFeature::default()).unwrap();
    let before = group.clone();

    assert!(group.subset(&ids(&["t1"])).is_err());
    assert_eq!(group, before);
}

/// The subset identifier ceiling is checked before anything is allocated.
#[test]
fn subset_refuses_too_many_identifiers() {
    let group = Group::new();
    let too_many = vec![String::new(); Group::MAX_ITEMS + 1];
    assert!(matches!(
        group.subset(&too_many).unwrap_err(),
        Error::InvalidValue(_)
    ));
    assert!(matches!(
        group.subset_dependent(&too_many).unwrap_err(),
        Error::InvalidValue(_)
    ));
}

/// The score-block ceiling is checked before any metadata is written.
#[test]
fn id_scores_refuse_oversized_input() {
    let scores = OpenSwathIndScores {
        ind_area_intensity: vec![0.0; MRMFeature::MAX_SCORE_VALUES + 1],
        ..OpenSwathIndScores::default()
    };
    let mut feature = MRMFeature::new();
    assert!(matches!(
        feature.id_scores_as_meta_value(false, &scores).unwrap_err(),
        Error::InvalidValue(_)
    ));
    assert!(feature.base.metadata.is_empty());
}

/// `MRMFeature` dereferences to its `Feature`, so the base kernel API is
/// reachable without naming the field, as `MRMFeature : public Feature` is in
/// the source.
#[test]
fn mrm_feature_derefs_to_feature() {
    let mut feature = MRMFeature::from(Feature::new(10.0, 500.0, 7.0));
    assert_eq!(feature.rt, 10.0);
    assert_eq!(feature.mz, 500.0);
    assert_eq!(feature.intensity, 7.0);
    feature.quality_rt = 0.5;
    assert_eq!(feature.feature.quality_rt, 0.5);
    feature.validate().unwrap();
}

/// `MSSpectrum` is accepted as the raw-data type, as the source header requires.
#[test]
fn a_group_of_spectra_tracks_native_ids() {
    use openms::kernel::MSSpectrum;
    let mut group: MRMTransitionGroup<MSSpectrum, SimpleTransition> = MRMTransitionGroup::new();
    let spectrum = MSSpectrum {
        native_id: "scan=1".into(),
        ..MSSpectrum::default()
    };
    group.add_chromatogram(spectrum, "scan=1").unwrap();
    assert!(group.chromatogram_ids_match());
    assert_eq!(group.size(), 1);
}
