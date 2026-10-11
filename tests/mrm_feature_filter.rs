// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! MRMFeatureFilter.
//!
//! Every section of the upstream class test at the core pin `bc9cc12`
//! (`src/tests/class_tests/openms/source/MRMFeatureFilter_test.cpp`, evidence
//! tier 3) is transcribed here, mechanically for the sections that build
//! feature maps and QC criteria (a line-for-line translation of each C++
//! statement; no expected value was edited), by hand for the helper
//! sections. The section named in each comment is the upstream section.
//!
//! One convention of the upstream framework matters: `TEST_EQUAL(a, b)`
//! compares `a == T1(b)` (`ClassTest.h` `testEqual`), so an expected value is
//! converted to the type of the tested value, truncating a `double` expected
//! for an `Int`. The transcription keeps that: `n_detecting_u, 43.30127019`
//! asserts `43`.
//!
//! Tests marked "native" cover behaviour the upstream test does not reach
//! (tier 4: independently derived values or invariants).
#![allow(
    unused_mut,
    unused_assignments,
    clippy::needless_range_loop,
    clippy::excessive_precision,
    clippy::approx_constant,
    clippy::identity_op,
    clippy::too_many_lines
)]

use openms::Error;
use openms::analysis::openswath::mrm_feature_filter::{
    FlagOrFilter, MRMFeatureFilter, QuantitationBounds,
};
use openms::analysis::openswath::mrm_feature_qc::{ComponentGroupQCs, ComponentQCs, MRMFeatureQC};
use openms::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition;
use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
use openms::kernel::{Feature, FeatureMap};
use openms::metadata::MetaValue;
use openms::param::{Param, ParamValue};

/// `TEST_REAL_SIMILAR` with tolerances stricter than the upstream defaults.
fn real_similar(a: f64, b: f64) {
    if a == b {
        return;
    }
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

fn with_meta(pairs: &[(&str, MetaValue)]) -> Feature {
    let mut f = Feature::default();
    for (k, v) in pairs {
        f.metadata.insert((*k).to_owned(), v.clone());
    }
    f
}

fn float(v: f64) -> MetaValue {
    MetaValue::try_from(v).unwrap()
}

#[test]
fn constructor_defaults() {
    // MRMFeatureFilter_test.cpp: MRMFeatureFilter() and its defaults.
    let filter = MRMFeatureFilter::new();
    assert_eq!(filter.flag_or_filter, FlagOrFilter::Flag);
    let defaults = MRMFeatureFilter::defaults().unwrap();
    assert_eq!(
        defaults.value("flag_or_filter").unwrap(),
        &ParamValue::String("flag".into())
    );
    assert!(defaults.has_tag("flag_or_filter", "advanced").unwrap());
    assert_eq!(
        defaults.valid_strings("flag_or_filter").unwrap(),
        ["flag".to_owned(), "filter".to_owned()]
    );
}

#[test]
fn range_helpers() {
    // MRMFeatureFilter_test.cpp: the commented-out checkRange, updateRange,
    // setRange and initRange sections (templates, not instantiated upstream).
    // Expected values follow the header comments and the source bodies.
    let f = MRMFeatureFilter::default();
    assert!(f.check_range(5.0, 4.0, 6.0));
    assert!(f.check_range(6.0, 4.0, 6.0));
    assert!(!f.check_range(3.0, 4.0, 6.0));
    assert!(f.check_range(5, 4, 6));
    let (mut l, mut u) = (4.0, 6.0);
    f.update_range(5.0, &mut l, &mut u);
    assert_eq!((l, u), (4.0, 6.0));
    f.update_range(7.0, &mut l, &mut u);
    assert_eq!((l, u), (4.0, 7.0));
    f.update_range(3.0, &mut l, &mut u);
    assert_eq!((l, u), (3.0, 7.0));
    f.set_range(5.0, &mut l, &mut u);
    assert_eq!((l, u), (0.0, 5.0));
    f.set_range(-2.0, &mut l, &mut u);
    assert_eq!((l, u), (-2.0, 0.0));
    let (mut a, mut b) = (1, 2);
    f.init_range(6, &mut a, &mut b);
    assert_eq!((a, b), (6, 6));
}

#[test]
fn section_3_calculate_ion_ratio() {
    // START_SECTION(double calculateIonRatio(...))
    let mrmff = MRMFeatureFilter::default();
    let feature_name = "peak_apex_int";
    let inf = f64::INFINITY;
    let mut component_1 = Feature::default();
    let mut component_2 = Feature::default();
    component_1.metadata.insert(feature_name.into(), float(5.0));
    component_1
        .metadata
        .insert("native_id".into(), "component1".into());
    component_2.metadata.insert(feature_name.into(), float(5.0));
    component_2
        .metadata
        .insert("native_id".into(), "component2".into());
    real_similar(
        mrmff
            .calculate_ion_ratio(&component_1, &component_2, feature_name)
            .unwrap(),
        1.0,
    );
    component_2.metadata.insert(feature_name.into(), float(0.0));
    real_similar(
        mrmff
            .calculate_ion_ratio(&component_1, &component_2, feature_name)
            .unwrap(),
        inf,
    );
    let mut component_3 = Feature::default();
    let mut component_4 = Feature::default();
    component_3.metadata.insert("peak_area".into(), float(5.0));
    component_3
        .metadata
        .insert("native_id".into(), "component3".into());
    component_4.metadata.insert("peak_area".into(), float(5.0));
    component_4
        .metadata
        .insert("native_id".into(), "component4".into());
    real_similar(
        mrmff
            .calculate_ion_ratio(&component_1, &component_4, feature_name)
            .unwrap(),
        5.0,
    );
    real_similar(
        mrmff
            .calculate_ion_ratio(&component_3, &component_4, feature_name)
            .unwrap(),
        0.0,
    );
    // feature_name == "intensity"
    let mut component_5 = Feature::default();
    let mut component_6 = Feature::default();
    let mut component_7 = Feature::default();
    let component_8 = Feature::default();
    component_5
        .metadata
        .insert("native_id".into(), "component5".into());
    component_6
        .metadata
        .insert("native_id".into(), "component6".into());
    component_5.intensity = 3.0;
    component_6.intensity = 4.0;
    let ratio = |a: &Feature, b: &Feature| mrmff.calculate_ion_ratio(a, b, "intensity").unwrap();
    real_similar(ratio(&component_5, &component_6), 0.75);
    real_similar(ratio(&component_6, &component_5), 1.33333333333333);
    component_7
        .metadata
        .insert("native_id".into(), "component7".into());
    real_similar(ratio(&component_5, &component_7), inf);
    real_similar(ratio(&component_5, &component_8), 3.0);
}

#[test]
fn section_4_check_meta_value() {
    // START_SECTION(bool checkMetaValue(...))
    let mrmff = MRMFeatureFilter::default();
    let feature_name = "peak_apex_int";
    let mut component_1 = with_meta(&[
        (feature_name, float(5.0)),
        ("native_id", "component1".into()),
    ]);
    let (l, u) = (4.0, 6.0);
    assert_eq!(
        mrmff
            .check_meta_value(&component_1, feature_name, l, u)
            .unwrap(),
        (true, true)
    ); // pass case
    component_1.metadata.insert(feature_name.into(), float(6.0));
    assert_eq!(
        mrmff
            .check_meta_value(&component_1, feature_name, l, u)
            .unwrap(),
        (true, true)
    ); // edge pass case
    component_1.metadata.insert(feature_name.into(), float(3.0));
    assert_eq!(
        mrmff
            .check_meta_value(&component_1, feature_name, l, u)
            .unwrap(),
        (false, true)
    ); // fail case
    assert_eq!(
        mrmff
            .check_meta_value(&component_1, "peak_area", l, u)
            .unwrap(),
        (true, false)
    ); // not found case
}

#[test]
fn section_5_update_meta_value() {
    // START_SECTION(void updateMetaValue(...))
    let mrmff = MRMFeatureFilter::default();
    let feature_name = "peak_apex_int";
    let mut component_1 = with_meta(&[
        (feature_name, float(5.0)),
        ("native_id", "component1".into()),
    ]);
    let (mut l, mut u) = (4.0, 6.0);
    let exists = mrmff
        .update_meta_value(&component_1, feature_name, &mut l, &mut u)
        .unwrap();
    assert_eq!((l, u, exists), (4.0, 6.0, true));
    component_1.metadata.insert(feature_name.into(), float(7.0));
    let exists = mrmff
        .update_meta_value(&component_1, feature_name, &mut l, &mut u)
        .unwrap();
    assert_eq!((l, u, exists), (4.0, 7.0, true));
    component_1.metadata.insert(feature_name.into(), float(3.0));
    let exists = mrmff
        .update_meta_value(&component_1, feature_name, &mut l, &mut u)
        .unwrap();
    assert_eq!((l, u, exists), (3.0, 7.0, true));
    let exists = mrmff
        .update_meta_value(&component_1, "peak_area", &mut l, &mut u)
        .unwrap();
    assert_eq!((l, u, exists), (3.0, 7.0, false));
}

#[test]
fn section_6_set_meta_value() {
    // START_SECTION(void setMetaValue(...))
    let mrmff = MRMFeatureFilter::default();
    let feature_name = "peak_apex_int";
    let mut component_1 = with_meta(&[
        (feature_name, float(5.0)),
        ("native_id", "component1".into()),
    ]);
    let (mut l, mut u) = (4.0, 6.0);
    let exists = mrmff
        .set_meta_value(&component_1, feature_name, &mut l, &mut u)
        .unwrap();
    assert_eq!((l, u, exists), (0.0, 5.0, true));
    component_1.metadata.insert(feature_name.into(), float(7.0));
    l = 4.0;
    u = 6.0;
    let exists = mrmff
        .set_meta_value(&component_1, feature_name, &mut l, &mut u)
        .unwrap();
    assert_eq!((l, u, exists), (0.0, 7.0, true));
    component_1
        .metadata
        .insert(feature_name.into(), float(-1.0));
    let exists = mrmff
        .set_meta_value(&component_1, feature_name, &mut l, &mut u)
        .unwrap();
    assert_eq!((l, u, exists), (-1.0, 0.0, true));
    let exists = mrmff
        .set_meta_value(&component_1, "peak_area", &mut l, &mut u)
        .unwrap();
    assert_eq!((l, u, exists), (-1.0, 0.0, false));
}

#[test]
fn section_7_init_meta_value() {
    // START_SECTION(void initMetaValue(...))
    let mrmff = MRMFeatureFilter::default();
    let feature_name = "peak_apex_int";
    let component_1 = with_meta(&[
        (feature_name, float(5.0)),
        ("native_id", "component1".into()),
    ]);
    let (mut l, mut u) = (4.0, 6.0);
    let exists = mrmff
        .init_meta_value(&component_1, feature_name, &mut l, &mut u)
        .unwrap();
    assert_eq!((l, u, exists), (5.0, 5.0, true));
    l = 4.0;
    u = 6.0;
    let exists = mrmff
        .init_meta_value(&component_1, "peak_area", &mut l, &mut u)
        .unwrap();
    assert_eq!((l, u, exists), (4.0, 6.0, false));
}

#[test]
fn section_8_count_labels_and_transition_types() {
    // START_SECTION(countLabelsAndTransitionTypes(...))
    let mrmff = MRMFeatureFilter::default();
    let mut component_1 = Feature::default();
    let mut subordinates = Vec::new();
    for (id, label) in [
        ("component1.1.Heavy", "Heavy"),
        ("component1.1.Light", "Light"),
        ("component1.2.Light", "Light"),
    ] {
        subordinates.push(with_meta(&[
            ("native_id", id.into()),
            ("LabelType", label.into()),
        ]));
    }
    component_1.subordinates = subordinates;

    let mut transitions = TargetedExperiment::new();
    for (id, group, quantifying) in [
        ("component1.1.Heavy", "component_group1", true),
        ("component1.1.Light", "component_group1", true),
        ("component1.2.Light", "component_group1", false),
        ("component2.1.Heavy", "component_group2", true),
        ("component2.1.Light", "component_group2", true),
    ] {
        let mut t = ReactionMonitoringTransition::new();
        t.set_native_id(id);
        t.set_peptide_ref(group);
        t.set_detecting_transition(true);
        t.set_identifying_transition(false);
        t.set_quantifying_transition(quantifying);
        transitions.add_transition(t);
    }
    let counts = mrmff.count_labels_and_transition_types(&component_1, &transitions);
    assert_eq!(counts.n_heavy, 1);
    assert_eq!(counts.n_light, 2);
    assert_eq!(counts.n_quantifying, 2);
    assert_eq!(counts.n_identifying, 0);
    assert_eq!(counts.n_detecting, 3);
    assert_eq!(counts.n_transitions, 3);
}

#[test]
fn native_unmatched_subordinate_counts_as_a_default_transition() {
    // native (tier 4): a subordinate whose native ID is not in the experiment
    // is counted with the flags of a default-constructed transition
    // (detecting and quantifying, not identifying), as the source's local
    // `ReactionMonitoringTransition transition;` is.
    let mrmff = MRMFeatureFilter::default();
    let mut group = Feature::default();
    group.subordinates.push(with_meta(&[
        ("native_id", "unknown".into()),
        ("LabelType", "Heavy".into()),
    ]));
    group.subordinates.push(Feature::default());
    let counts = mrmff.count_labels_and_transition_types(&group, &TargetedExperiment::new());
    assert_eq!(counts.n_heavy, 1);
    assert_eq!(counts.n_quantifying, 2);
    assert_eq!(counts.n_detecting, 2);
    assert_eq!(counts.n_identifying, 0);
    assert_eq!(counts.n_transitions, 2);
}

#[test]
fn native_parameters_round_trip() {
    // native (tier 4): flag_or_filter is restricted to its two valid strings.
    let mut params = Param::new();
    params
        .set_value(
            "flag_or_filter",
            ParamValue::String("filter".into()),
            "",
            &[],
        )
        .unwrap();
    let f = MRMFeatureFilter::from_param(&params).unwrap();
    assert_eq!(f.flag_or_filter, FlagOrFilter::Filter);
    assert_eq!(
        f.to_param().unwrap().value("flag_or_filter").unwrap(),
        &ParamValue::String("filter".into())
    );
    params
        .set_value(
            "flag_or_filter",
            ParamValue::String("remove".into()),
            "",
            &[],
        )
        .unwrap();
    assert!(matches!(
        MRMFeatureFilter::from_param(&params),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn native_filter_mode_keeps_the_map_level_members() {
    // native (tier 4): filter mode replaces the features only. The source
    // assigns a fresh FeatureMap and so also drops the identifier.
    let mrmff = MRMFeatureFilter {
        flag_or_filter: FlagOrFilter::Filter,
    };
    let mut sub = with_meta(&[("native_id", "c1".into())]);
    sub.intensity = 10.0;
    let mut feature = with_meta(&[("PeptideRef", "g1".into())]);
    feature.subordinates.push(sub);
    let mut map = FeatureMap::from_features(vec![feature]);
    map.identifier = "map-id".into();
    let mut qc = MRMFeatureQC::default();
    qc.component_group_qcs.push(ComponentGroupQCs {
        component_group_name: "g1".into(),
        ..ComponentGroupQCs::default()
    });
    mrmff
        .filter_feature_map(&mut map, &qc, &TargetedExperiment::new())
        .unwrap();
    assert_eq!(map.len(), 1);
    assert_eq!(map.identifier, "map-id");
}

#[test]
fn native_counts_divide_like_the_source_or_refuse() {
    // native (tier 4): the mean of the Int counts divides an Int by the number
    // of samples; the source crashes on zero samples (integer division by
    // zero), this port refuses.
    let mrmff = MRMFeatureFilter::default();
    let mut template = MRMFeatureQC::default();
    template
        .component_group_qcs
        .push(ComponentGroupQCs::default());
    let mut mean = MRMFeatureQC::default();
    let err = mrmff
        .calculate_filter_values_mean(&mut mean, &[], &template)
        .unwrap_err();
    assert!(matches!(err, Error::InvalidValue(_)));
    // without component group QCs only doubles are divided: 0 / 0 is NaN
    let mut only_components = MRMFeatureQC::default();
    only_components.component_qcs.push(ComponentQCs::default());
    mrmff
        .calculate_filter_values_mean(&mut mean, &[], &only_components)
        .unwrap();
    assert!(mean.component_qcs[0].intensity_l.is_nan());
    // an Int count is divided as `unsigned long`: 7 / 2 = 3
    let mut a = template.clone();
    let mut b = template.clone();
    a.component_group_qcs[0].n_heavy_l = 3;
    b.component_group_qcs[0].n_heavy_l = 4;
    mrmff
        .calculate_filter_values_mean(&mut mean, &[a, b], &template)
        .unwrap();
    assert_eq!(mean.component_group_qcs[0].n_heavy_l, 3);
}

struct Method(&'static str, f64, f64);
impl QuantitationBounds for Method {
    fn component_name(&self) -> &str {
        self.0
    }
    fn lloq(&self) -> f64 {
        self.1
    }
    fn uloq(&self) -> f64 {
        self.2
    }
}

#[test]
fn native_transfer_lloq_and_uloq() {
    // native (tier 4): MRMFeatureFilter.cpp TransferLLOQAndULOQToCalculated-
    // ConcentrationBounds; the upstream test has no section for it.
    let mrmff = MRMFeatureFilter::default();
    let mut qc = MRMFeatureQC::default();
    let mut a = ComponentQCs {
        component_name: "a".into(),
        ..ComponentQCs::default()
    };
    a.meta_value_qc
        .insert("calculated_concentration".into(), (0.0, 1e12));
    let b = ComponentQCs {
        component_name: "b".into(),
        ..ComponentQCs::default()
    };
    qc.component_qcs = vec![a, b];
    // zero LLOQ and ULOQ is skipped, even for a component without the bound
    mrmff
        .transfer_lloq_and_uloq_to_calculated_concentration_bounds(
            &[Method("b", 0.0, 0.0)],
            &mut qc,
        )
        .unwrap();
    mrmff
        .transfer_lloq_and_uloq_to_calculated_concentration_bounds(
            &[Method("a", 1.5, 9.5)],
            &mut qc,
        )
        .unwrap();
    assert_eq!(
        qc.component_qcs[0].meta_value_qc["calculated_concentration"],
        (1.5, 9.5)
    );
    // a matching component with no such bound is the source's out_of_range
    assert!(
        mrmff
            .transfer_lloq_and_uloq_to_calculated_concentration_bounds(
                &[Method("b", 1.0, 2.0)],
                &mut qc,
            )
            .is_err()
    );
}

// START_SECTION(void FilterFeatureMap(FeatureMap& features, MRMFeatureQC& filter_criteria, const TargetedExperiment & transitions))

#[test]
fn section_9_filter_feature_map_flag_and_filter() {
    {
        let mut mrmff = MRMFeatureFilter::default();
        let mut components = FeatureMap::default();
        let mut component_1 = Feature::default();
        let mut subordinate = Feature::default();
        let mut subordinates: Vec<Feature> = Vec::new();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(500_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component2.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(1000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component2.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(1000_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group2"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        let mut transitions = TargetedExperiment::new();
        let mut transition = ReactionMonitoringTransition::new();
        transition.set_native_id("component1.1.Heavy");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component1.1.Light");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component1.2.Light");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(false);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component2.1.Heavy");
        transition.set_peptide_ref("component_group2");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component2.1.Light");
        transition.set_peptide_ref("component_group2");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        let mut qc_criteria = MRMFeatureQC::default();
        let mut cgqcs = ComponentGroupQCs::default();
        let mut cqcs = ComponentQCs::default();
        let lbub: (f64, f64) = (501.0, 4e6);
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.n_heavy_l = 1;
        cgqcs.n_heavy_u = 1;
        cgqcs.n_light_l = 2;
        cgqcs.n_light_u = 2;
        cgqcs.n_detecting_l = 2;
        cgqcs.n_detecting_u = 3;
        cgqcs.n_quantifying_l = 2;
        cgqcs.n_quantifying_u = 2;
        cgqcs.n_identifying_l = 0;
        cgqcs.n_identifying_u = 3;
        cgqcs.n_transitions_l = 2;
        cgqcs.n_transitions_u = 3;
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.5;
        cgqcs.ion_ratio_u = 10.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc.insert("peak_apex_int".into(), lbub);
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc.insert("peak_apex_int".into(), lbub);
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc.insert("peak_apex_int".into(), lbub);
        qc_criteria.component_group_qcs.push(cgqcs.clone());
        qc_criteria.component_qcs.push(cqcs.clone());
        cgqcs.component_group_name = "component_group2".to_owned();
        cgqcs.n_heavy_l = 1;
        cgqcs.n_heavy_u = 1;
        cgqcs.n_light_l = 2;
        cgqcs.n_light_u = 2;
        cgqcs.n_detecting_l = 2;
        cgqcs.n_detecting_u = 3;
        cgqcs.n_quantifying_l = 2;
        cgqcs.n_quantifying_u = 2;
        cgqcs.n_identifying_l = 0;
        cgqcs.n_identifying_u = 3;
        cgqcs.n_transitions_l = 3;
        cgqcs.n_transitions_u = 3;
        cgqcs.ion_ratio_pair_name_1 = "component2.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component2.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.5;
        cgqcs.ion_ratio_u = 10.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component2.1.Heavy".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc.insert("peak_apex_int".into(), lbub);
        cqcs.component_name = "component2.1.Light".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc.insert("peak_apex_int".into(), lbub);
        qc_criteria.component_group_qcs.push(cgqcs.clone());
        qc_criteria.component_qcs.push(cqcs.clone());
        let mut params = Param::new();
        params
            .set_value("flag_or_filter", ParamValue::String("flag".into()), "", &[])
            .unwrap();
        mrmff = MRMFeatureFilter::from_param(&params).unwrap();
        mrmff
            .filter_feature_map(&mut components, &qc_criteria, &transitions)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_message"]
                .as_string_list()
                .unwrap()[0],
            "peak_apex_int"
        );
        assert_eq!(
            components.features[1].metadata["QC_transition_group_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[1].metadata["QC_transition_group_message"]
                .as_string_list()
                .unwrap()
                .len(),
            2
        );
        assert_eq!(
            components.features[1].metadata["QC_transition_group_message"]
                .as_string_list()
                .unwrap()[0],
            "n_light"
        );
        assert_eq!(
            components.features[1].metadata["QC_transition_group_message"]
                .as_string_list()
                .unwrap()[1],
            "n_transitions"
        );
        assert_eq!(
            components.features[1].subordinates[0].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[1].subordinates[1].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[2].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            0.75,
        );
        real_similar(
            components.features[1].metadata["QC_transition_group_score"]
                .as_f64()
                .unwrap(),
            0.777777777777778,
        );
        real_similar(
            components.features[1].subordinates[0].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[1].subordinates[1].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        params
            .set_value(
                "flag_or_filter",
                ParamValue::String("filter".into()),
                "",
                &[],
            )
            .unwrap();
        mrmff = MRMFeatureFilter::from_param(&params).unwrap();
        mrmff
            .filter_feature_map(&mut components, &qc_criteria, &transitions)
            .unwrap();
        assert_eq!(components.len(), 1);
        assert_eq!(components.features[0].subordinates.len(), 2);
    }
}

// START_SECTION(void FilterFeatureMap(FeatureMap& features, MRMFeatureQC& filter_criteria, const TargetedExperiment & transitions))

#[test]
fn section_10_filter_feature_map_individual_checks() {
    {
        let mut mrmff = MRMFeatureFilter::default();
        let mut components = FeatureMap::default();
        let mut component_1 = Feature::default();
        let mut subordinate = Feature::default();
        let mut subordinates: Vec<Feature> = Vec::new();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        let mut transitions = TargetedExperiment::new();
        let mut transition = ReactionMonitoringTransition::new();
        transition.set_native_id("component1.1.Heavy");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component1.1.Light");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component1.2.Light");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(false);
        transitions.add_transition(transition.clone());
        let mut qc_criteria = MRMFeatureQC::default();
        let mut cgqcs = ComponentGroupQCs::default();
        let mut cqcs = ComponentQCs::default();
        let lbub: (f64, f64) = (500.0, 4e6);
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.n_heavy_l = 1;
        cgqcs.n_heavy_u = 1;
        cgqcs.n_light_l = 1;
        cgqcs.n_light_u = 2;
        cgqcs.n_detecting_l = 2;
        cgqcs.n_detecting_u = 3;
        cgqcs.n_quantifying_l = 2;
        cgqcs.n_quantifying_u = 2;
        cgqcs.n_identifying_l = 0;
        cgqcs.n_identifying_u = 3;
        cgqcs.n_transitions_l = 3;
        cgqcs.n_transitions_u = 3;
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.5;
        cgqcs.ion_ratio_u = 2.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc.insert("peak_apex_int".into(), lbub);
        cqcs.meta_value_qc.insert("peak_area".into(), lbub);
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc.insert("peak_apex_int".into(), lbub);
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc.insert("peak_apex_int".into(), lbub);
        qc_criteria.component_group_qcs.push(cgqcs.clone());
        qc_criteria.component_qcs.push(cqcs.clone());
        let mut params = Param::new();
        params
            .set_value("flag_or_filter", ParamValue::String("flag".into()), "", &[])
            .unwrap();
        mrmff = MRMFeatureFilter::from_param(&params).unwrap();
        mrmff
            .filter_feature_map(&mut components, &qc_criteria, &transitions)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[2].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        components.features.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 6.0;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        mrmff
            .filter_feature_map(&mut components, &qc_criteria, &transitions)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_message"]
                .as_string_list()
                .unwrap()[0],
            "retention_time"
        );
        components.features.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 0.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        mrmff
            .filter_feature_map(&mut components, &qc_criteria, &transitions)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_message"]
                .as_string_list()
                .unwrap()[0],
            "intensity"
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[2].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            0.75,
        );
        components.features.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 0.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        mrmff
            .filter_feature_map(&mut components, &qc_criteria, &transitions)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_message"]
                .as_string_list()
                .unwrap()[0],
            "overall_quality"
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[2].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            0.75,
        );
        components.features.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(500_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(400_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        mrmff
            .filter_feature_map(&mut components, &qc_criteria, &transitions)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_message"]
                .as_string_list()
                .unwrap()[0],
            "peak_apex_int"
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[2].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            0.75,
        );
        components.features.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        mrmff
            .filter_feature_map(&mut components, &qc_criteria, &transitions)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].metadata["QC_transition_group_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].metadata["QC_transition_group_message"]
                .as_string_list()
                .unwrap()[0],
            "n_heavy"
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_score"]
                .as_f64()
                .unwrap(),
            0.892857142857143,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[2].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        components.features.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from(""));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from(""));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        mrmff
            .filter_feature_map(&mut components, &qc_criteria, &transitions)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].metadata["QC_transition_group_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].metadata["QC_transition_group_message"]
                .as_string_list()
                .unwrap()[0],
            "n_light"
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_score"]
                .as_f64()
                .unwrap(),
            0.892857142857143,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[2].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        components.features.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        mrmff
            .filter_feature_map(&mut components, &qc_criteria, &transitions)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].metadata["QC_transition_group_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].metadata["QC_transition_group_message"]
                .as_string_list()
                .unwrap()[0],
            "n_transitions"
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_score"]
                .as_f64()
                .unwrap(),
            0.888888888888889,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        components.features.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(500_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        mrmff
            .filter_feature_map(&mut components, &qc_criteria, &transitions)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].metadata["QC_transition_group_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].metadata["QC_transition_group_message"]
                .as_string_list()
                .unwrap()[0],
            "ion_ratio_pair[component1.1.Light/component1.2.Light]"
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_pass"],
            MetaValue::from(1_i32)
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_score"]
                .as_f64()
                .unwrap(),
            0.964285714285714,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[2].metadata["QC_transition_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        components.features.clear();
    }
}

// START_SECTION(void EstimateDefaultMRMFeatureQCValues(const std::vector<FeatureMap>& samples, MRMFeatureQC& filter_template, const TargetedExperiment&

#[test]
fn section_11_estimate_default_mrm_feature_qc_values() {
    {
        let mut mrmff = MRMFeatureFilter::default();
        let mut samples: Vec<FeatureMap> = Vec::new();
        let mut components = FeatureMap::default();
        let mut component_1 = Feature::default();
        let mut subordinate = Feature::default();
        let mut subordinates: Vec<Feature> = Vec::new();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(500_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.rt = 2.5;
        component_1.intensity = 15000.0_f32;
        component_1.quality = 300.0_f32;
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        samples.push(components.clone());
        components.features.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 3.0;
        subordinate.intensity = 1000.0_f32;
        subordinate.quality = 200.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(1000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 3.0;
        subordinate.intensity = 1000.0_f32;
        subordinate.quality = 400.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(2000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 3.1;
        subordinate.intensity = 2000.0_f32;
        subordinate.quality = 300.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(800_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.rt = 3.0;
        component_1.intensity = 5000.0_f32;
        component_1.quality = 200.0_f32;
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        samples.push(components.clone());
        components.features.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.0;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.0;
        subordinate.intensity = 4000.0_f32;
        subordinate.quality = 150.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(6000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.0;
        subordinate.intensity = 1000.0_f32;
        subordinate.quality = 300.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(100_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.rt = 2.0;
        component_1.intensity = 15000.0_f32;
        component_1.quality = 500.0_f32;
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        samples.push(components.clone());
        components.features.clear();
        let mut transitions = TargetedExperiment::new();
        let mut transition = ReactionMonitoringTransition::new();
        transition.set_native_id("component1.1.Heavy");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component1.1.Light");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component1.2.Light");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(false);
        transitions.add_transition(transition.clone());
        let mut qc_criteria1 = MRMFeatureQC::default();
        let mut qc_criteria2 = MRMFeatureQC::default();
        let mut cgqcs = ComponentGroupQCs::default();
        let mut cqcs = ComponentQCs::default();
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.n_heavy_l = 0;
        cgqcs.n_heavy_u = 0;
        cgqcs.n_light_l = 0;
        cgqcs.n_light_u = 0;
        cgqcs.n_detecting_l = 0;
        cgqcs.n_detecting_u = 0;
        cgqcs.n_quantifying_l = 0;
        cgqcs.n_quantifying_u = 0;
        cgqcs.n_identifying_l = 0;
        cgqcs.n_identifying_u = 0;
        cgqcs.n_transitions_l = 0;
        cgqcs.n_transitions_u = 0;
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 0.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 0.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 0.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 0.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 0.0));
        cqcs.meta_value_qc.insert("peak_area".into(), (0.0, 0.0));
        qc_criteria1.component_group_qcs.push(cgqcs.clone());
        qc_criteria1.component_qcs.push(cqcs.clone());
        qc_criteria2 = qc_criteria1.clone();
        mrmff
            .estimate_default_mrm_feature_qc_values(
                &samples,
                &mut qc_criteria1,
                &transitions,
                false,
            )
            .unwrap();
        assert_eq!(
            qc_criteria1.component_group_qcs[0].component_group_name,
            "component_group1"
        );
        assert_eq!(qc_criteria1.component_group_qcs[0].n_heavy_l, 0);
        assert_eq!(qc_criteria1.component_group_qcs[0].n_heavy_u, 1);
        assert_eq!(qc_criteria1.component_group_qcs[0].n_light_l, 0);
        assert_eq!(qc_criteria1.component_group_qcs[0].n_light_u, 2);
        assert_eq!(qc_criteria1.component_group_qcs[0].n_detecting_l, 0);
        assert_eq!(qc_criteria1.component_group_qcs[0].n_detecting_u, 3);
        assert_eq!(qc_criteria1.component_group_qcs[0].n_quantifying_l, 0);
        assert_eq!(qc_criteria1.component_group_qcs[0].n_quantifying_u, 2);
        assert_eq!(qc_criteria1.component_group_qcs[0].n_identifying_l, 0);
        assert_eq!(qc_criteria1.component_group_qcs[0].n_identifying_u, 0);
        assert_eq!(qc_criteria1.component_group_qcs[0].n_transitions_l, 0);
        assert_eq!(qc_criteria1.component_group_qcs[0].n_transitions_u, 3);
        assert_eq!(
            qc_criteria1.component_group_qcs[0].ion_ratio_pair_name_1,
            "component1.1.Light"
        );
        assert_eq!(
            qc_criteria1.component_group_qcs[0].ion_ratio_pair_name_2,
            "component1.2.Light"
        );
        real_similar(qc_criteria1.component_group_qcs[0].ion_ratio_l, 0.0);
        real_similar(qc_criteria1.component_group_qcs[0].ion_ratio_u, 60.0);
        assert_eq!(
            qc_criteria1.component_group_qcs[0].ion_ratio_feature_name,
            "peak_apex_int"
        );
        assert_eq!(
            qc_criteria1.component_qcs[0].component_name,
            "component1.1.Heavy"
        );
        real_similar(qc_criteria1.component_qcs[0].retention_time_l, 0.0);
        real_similar(qc_criteria1.component_qcs[0].retention_time_u, 3.0);
        real_similar(qc_criteria1.component_qcs[0].intensity_l, 0.0);
        real_similar(qc_criteria1.component_qcs[0].intensity_u, 5000.0);
        real_similar(qc_criteria1.component_qcs[0].overall_quality_l, 0.0);
        real_similar(qc_criteria1.component_qcs[0].overall_quality_u, 200.0);
        real_similar(
            qc_criteria1.component_qcs[0].meta_value_qc["peak_apex_int"].0,
            0.0,
        );
        real_similar(
            qc_criteria1.component_qcs[0].meta_value_qc["peak_apex_int"].1,
            5000.0,
        );
        real_similar(
            qc_criteria1.component_qcs[0].meta_value_qc["peak_area"].0,
            0.0,
        );
        real_similar(
            qc_criteria1.component_qcs[0].meta_value_qc["peak_area"].1,
            0.0,
        );
        mrmff
            .estimate_default_mrm_feature_qc_values(&samples, &mut qc_criteria2, &transitions, true)
            .unwrap();
        assert_eq!(
            qc_criteria2.component_group_qcs[0].component_group_name,
            "component_group1"
        );
        assert_eq!(qc_criteria2.component_group_qcs[0].n_heavy_l, 1);
        assert_eq!(qc_criteria2.component_group_qcs[0].n_heavy_u, 1);
        assert_eq!(qc_criteria2.component_group_qcs[0].n_light_l, 2);
        assert_eq!(qc_criteria2.component_group_qcs[0].n_light_u, 2);
        assert_eq!(qc_criteria2.component_group_qcs[0].n_detecting_l, 3);
        assert_eq!(qc_criteria2.component_group_qcs[0].n_detecting_u, 3);
        assert_eq!(qc_criteria2.component_group_qcs[0].n_quantifying_l, 2);
        assert_eq!(qc_criteria2.component_group_qcs[0].n_quantifying_u, 2);
        assert_eq!(qc_criteria2.component_group_qcs[0].n_identifying_l, 0);
        assert_eq!(qc_criteria2.component_group_qcs[0].n_identifying_u, 0);
        assert_eq!(qc_criteria2.component_group_qcs[0].n_transitions_l, 3);
        assert_eq!(qc_criteria2.component_group_qcs[0].n_transitions_u, 3);
        assert_eq!(
            qc_criteria2.component_group_qcs[0].ion_ratio_pair_name_1,
            "component1.1.Light"
        );
        assert_eq!(
            qc_criteria2.component_group_qcs[0].ion_ratio_pair_name_2,
            "component1.2.Light"
        );
        real_similar(qc_criteria2.component_group_qcs[0].ion_ratio_l, 2.5);
        real_similar(qc_criteria2.component_group_qcs[0].ion_ratio_u, 60.0);
        assert_eq!(
            qc_criteria2.component_group_qcs[0].ion_ratio_feature_name,
            "peak_apex_int"
        );
        assert_eq!(
            qc_criteria2.component_qcs[0].component_name,
            "component1.1.Heavy"
        );
        real_similar(qc_criteria2.component_qcs[0].retention_time_l, 2.0);
        real_similar(qc_criteria2.component_qcs[0].retention_time_u, 3.0);
        real_similar(qc_criteria2.component_qcs[0].intensity_l, 1000.0);
        real_similar(qc_criteria2.component_qcs[0].intensity_u, 5000.0);
        real_similar(qc_criteria2.component_qcs[0].overall_quality_l, 100.0);
        real_similar(qc_criteria2.component_qcs[0].overall_quality_u, 200.0);
        real_similar(
            qc_criteria2.component_qcs[0].meta_value_qc["peak_apex_int"].0,
            1000.0,
        );
        real_similar(
            qc_criteria2.component_qcs[0].meta_value_qc["peak_apex_int"].1,
            5000.0,
        );
        real_similar(
            qc_criteria2.component_qcs[0].meta_value_qc["peak_area"].0,
            0.0,
        );
        real_similar(
            qc_criteria2.component_qcs[0].meta_value_qc["peak_area"].1,
            0.0,
        );
    }
}

// START_SECTION(void EstimatePercRSD(const std::vector<FeatureMap>& samples, MRMFeatureQC& filter_template, const TargetedExperiment& transitions))

#[test]
fn section_12_estimate_perc_rsd() {
    {
        let mut mrmff = MRMFeatureFilter::default();
        let mut samples: Vec<FeatureMap> = Vec::new();
        let mut components = FeatureMap::default();
        let mut component_1 = Feature::default();
        let mut subordinate = Feature::default();
        let mut subordinates: Vec<Feature> = Vec::new();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(500_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component2.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(1000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component2.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(1000_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group2"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        samples.push(components.clone());
        samples.push(components.clone());
        samples.push(components.clone());
        let mut transitions = TargetedExperiment::new();
        let mut transition = ReactionMonitoringTransition::new();
        transition.set_native_id("component1.1.Heavy");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component1.1.Light");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component1.2.Light");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(false);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component2.1.Heavy");
        transition.set_peptide_ref("component_group2");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component2.1.Light");
        transition.set_peptide_ref("component_group2");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        let mut qc_criteria = MRMFeatureQC::default();
        let mut cgqcs = ComponentGroupQCs::default();
        let mut cqcs = ComponentQCs::default();
        let lbub: (f64, f64) = (500.0, 4e6);
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.n_heavy_l = 1;
        cgqcs.n_heavy_u = 1;
        cgqcs.n_light_l = 1;
        cgqcs.n_light_u = 2;
        cgqcs.n_detecting_l = 2;
        cgqcs.n_detecting_u = 3;
        cgqcs.n_quantifying_l = 2;
        cgqcs.n_quantifying_u = 2;
        cgqcs.n_identifying_l = 0;
        cgqcs.n_identifying_u = 3;
        cgqcs.n_transitions_l = 3;
        cgqcs.n_transitions_u = 3;
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.5;
        cgqcs.ion_ratio_u = 2.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc.insert("peak_apex_int".into(), lbub);
        cqcs.meta_value_qc.insert("peak_area".into(), lbub);
        qc_criteria.component_qcs.push(cqcs.clone());
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc.insert("peak_apex_int".into(), lbub);
        cqcs.meta_value_qc.insert("peak_area".into(), lbub);
        qc_criteria.component_qcs.push(cqcs.clone());
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc.insert("peak_apex_int".into(), lbub);
        cqcs.meta_value_qc.insert("peak_area".into(), lbub);
        qc_criteria.component_group_qcs.push(cgqcs.clone());
        qc_criteria.component_qcs.push(cqcs.clone());
        let mut filter_zeros = qc_criteria.clone();
        mrmff
            .estimate_perc_rsd(&samples, &mut filter_zeros, &transitions)
            .unwrap();
        assert_eq!(
            filter_zeros.component_group_qcs[0].component_group_name,
            "component_group1"
        );
        assert_eq!(filter_zeros.component_group_qcs[0].n_heavy_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_heavy_u, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_light_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_light_u, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_detecting_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_detecting_u, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_quantifying_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_quantifying_u, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_identifying_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_identifying_u, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_transitions_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_transitions_u, 0);
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_pair_name_1,
            "component1.1.Light"
        );
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_pair_name_2,
            "component1.2.Light"
        );
        real_similar(filter_zeros.component_group_qcs[0].ion_ratio_l, 0.0);
        real_similar(filter_zeros.component_group_qcs[0].ion_ratio_u, 0.0);
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_feature_name,
            "peak_apex_int"
        );
        assert_eq!(
            filter_zeros.component_qcs[0].component_name,
            "component1.1.Heavy"
        );
        real_similar(filter_zeros.component_qcs[0].retention_time_l, 0.0);
        real_similar(filter_zeros.component_qcs[0].retention_time_u, 0.0);
        real_similar(filter_zeros.component_qcs[0].intensity_l, 0.0);
        real_similar(filter_zeros.component_qcs[0].intensity_u, 0.0);
        real_similar(filter_zeros.component_qcs[0].overall_quality_l, 0.0);
        real_similar(filter_zeros.component_qcs[0].overall_quality_u, 0.0);
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_apex_int"].0,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_apex_int"].1,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_area"].0,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_area"].1,
            0.0,
        );
        assert_eq!(
            filter_zeros.component_qcs[1].component_name,
            "component1.1.Light"
        );
        real_similar(filter_zeros.component_qcs[1].retention_time_l, 0.0);
        real_similar(filter_zeros.component_qcs[1].retention_time_u, 0.0);
        real_similar(filter_zeros.component_qcs[1].intensity_l, 0.0);
        real_similar(filter_zeros.component_qcs[1].intensity_u, 0.0);
        real_similar(filter_zeros.component_qcs[1].overall_quality_l, 0.0);
        real_similar(filter_zeros.component_qcs[1].overall_quality_u, 0.0);
        real_similar(
            filter_zeros.component_qcs[1].meta_value_qc["peak_apex_int"].0,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[1].meta_value_qc["peak_apex_int"].1,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[1].meta_value_qc["peak_area"].0,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[1].meta_value_qc["peak_area"].1,
            0.0,
        );
        assert_eq!(
            filter_zeros.component_qcs[2].component_name,
            "component1.2.Light"
        );
        real_similar(filter_zeros.component_qcs[2].retention_time_l, 0.0);
        real_similar(filter_zeros.component_qcs[2].retention_time_u, 0.0);
        real_similar(filter_zeros.component_qcs[2].intensity_l, 0.0);
        real_similar(filter_zeros.component_qcs[2].intensity_u, 0.0);
        real_similar(filter_zeros.component_qcs[2].overall_quality_l, 0.0);
        real_similar(filter_zeros.component_qcs[2].overall_quality_u, 0.0);
        real_similar(
            filter_zeros.component_qcs[2].meta_value_qc["peak_apex_int"].0,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[2].meta_value_qc["peak_apex_int"].1,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[2].meta_value_qc["peak_area"].0,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[2].meta_value_qc["peak_area"].1,
            0.0,
        );
    }
}

// START_SECTION(void EstimateBackgroundInterferences(const std::vector<FeatureMap>& samples, MRMFeatureQC& filter_template, const TargetedExperiment& tr

#[test]
fn section_13_estimate_background_interferences() {
    {
        let mut mrmff = MRMFeatureFilter::default();
        let mut samples: Vec<FeatureMap> = Vec::new();
        let mut components = FeatureMap::default();
        let mut component_1 = Feature::default();
        let mut subordinate = Feature::default();
        let mut subordinates: Vec<Feature> = Vec::new();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(500_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component2.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(1000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component2.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(1000_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group2"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        samples.push(components.clone());
        samples.push(components.clone());
        samples.push(components.clone());
        let mut transitions = TargetedExperiment::new();
        let mut transition = ReactionMonitoringTransition::new();
        transition.set_native_id("component1.1.Heavy");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component1.1.Light");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component1.2.Light");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(false);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component2.1.Heavy");
        transition.set_peptide_ref("component_group2");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component2.1.Light");
        transition.set_peptide_ref("component_group2");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        let mut qc_criteria = MRMFeatureQC::default();
        let mut cgqcs = ComponentGroupQCs::default();
        let mut cqcs = ComponentQCs::default();
        let lbub: (f64, f64) = (500.0, 4e6);
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.n_heavy_l = 1;
        cgqcs.n_heavy_u = 1;
        cgqcs.n_light_l = 1;
        cgqcs.n_light_u = 2;
        cgqcs.n_detecting_l = 2;
        cgqcs.n_detecting_u = 3;
        cgqcs.n_quantifying_l = 2;
        cgqcs.n_quantifying_u = 2;
        cgqcs.n_identifying_l = 0;
        cgqcs.n_identifying_u = 3;
        cgqcs.n_transitions_l = 3;
        cgqcs.n_transitions_u = 3;
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.5;
        cgqcs.ion_ratio_u = 2.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc.insert("peak_apex_int".into(), lbub);
        cqcs.meta_value_qc.insert("peak_area".into(), lbub);
        qc_criteria.component_qcs.push(cqcs.clone());
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc.insert("peak_apex_int".into(), lbub);
        qc_criteria.component_qcs.push(cqcs.clone());
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc.insert("peak_apex_int".into(), lbub);
        qc_criteria.component_group_qcs.push(cgqcs.clone());
        qc_criteria.component_qcs.push(cqcs.clone());
        let mut filter_zeros = qc_criteria.clone();
        mrmff
            .estimate_background_interferences(&samples, &mut filter_zeros, &transitions)
            .unwrap();
        assert_eq!(
            filter_zeros.component_group_qcs[0].component_group_name,
            "component_group1"
        );
        assert_eq!(filter_zeros.component_group_qcs[0].n_heavy_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_heavy_u, 1);
        assert_eq!(filter_zeros.component_group_qcs[0].n_light_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_light_u, 2);
        assert_eq!(filter_zeros.component_group_qcs[0].n_detecting_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_detecting_u, 3);
        assert_eq!(filter_zeros.component_group_qcs[0].n_quantifying_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_quantifying_u, 2);
        assert_eq!(filter_zeros.component_group_qcs[0].n_identifying_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_identifying_u, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_transitions_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_transitions_u, 3);
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_pair_name_1,
            "component1.1.Light"
        );
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_pair_name_2,
            "component1.2.Light"
        );
        real_similar(filter_zeros.component_group_qcs[0].ion_ratio_l, 0.0);
        real_similar(filter_zeros.component_group_qcs[0].ion_ratio_u, 10.0);
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_feature_name,
            "peak_apex_int"
        );
        assert_eq!(
            filter_zeros.component_qcs[0].component_name,
            "component1.1.Heavy"
        );
        real_similar(filter_zeros.component_qcs[0].retention_time_l, 0.0);
        real_similar(filter_zeros.component_qcs[0].retention_time_u, 2.5);
        real_similar(filter_zeros.component_qcs[0].intensity_l, 0.0);
        real_similar(filter_zeros.component_qcs[0].intensity_u, 5000.0);
        real_similar(filter_zeros.component_qcs[0].overall_quality_l, 0.0);
        real_similar(filter_zeros.component_qcs[0].overall_quality_u, 100.0);
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_apex_int"].0,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_apex_int"].1,
            5000.0,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_area"].0,
            500.0,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_area"].1,
            4e6,
        );
        assert_eq!(
            filter_zeros.component_qcs[1].component_name,
            "component1.1.Light"
        );
        real_similar(filter_zeros.component_qcs[1].retention_time_l, 0.0);
        real_similar(filter_zeros.component_qcs[1].retention_time_u, 2.5);
        real_similar(filter_zeros.component_qcs[1].intensity_l, 0.0);
        real_similar(filter_zeros.component_qcs[1].intensity_u, 5000.0);
        real_similar(filter_zeros.component_qcs[1].overall_quality_l, 0.0);
        real_similar(filter_zeros.component_qcs[1].overall_quality_u, 100.0);
        real_similar(
            filter_zeros.component_qcs[1].meta_value_qc["peak_apex_int"].0,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[1].meta_value_qc["peak_apex_int"].1,
            5000.0,
        );
        real_similar(
            filter_zeros.component_qcs[1].meta_value_qc["peak_area"].0,
            500.0,
        );
        real_similar(
            filter_zeros.component_qcs[1].meta_value_qc["peak_area"].1,
            4e6,
        );
        assert_eq!(
            filter_zeros.component_qcs[2].component_name,
            "component1.2.Light"
        );
        real_similar(filter_zeros.component_qcs[2].retention_time_l, 0.0);
        real_similar(filter_zeros.component_qcs[2].retention_time_u, 2.5);
        real_similar(filter_zeros.component_qcs[2].intensity_l, 0.0);
        real_similar(filter_zeros.component_qcs[2].intensity_u, 5000.0);
        real_similar(filter_zeros.component_qcs[2].overall_quality_l, 0.0);
        real_similar(filter_zeros.component_qcs[2].overall_quality_u, 100.0);
        real_similar(
            filter_zeros.component_qcs[2].meta_value_qc["peak_apex_int"].0,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[2].meta_value_qc["peak_apex_int"].1,
            500.0,
        );
        real_similar(
            filter_zeros.component_qcs[2].meta_value_qc["peak_area"].0,
            500.0,
        );
        real_similar(
            filter_zeros.component_qcs[2].meta_value_qc["peak_area"].1,
            4e6,
        );
    }
}

// START_SECTION(void accumulateFilterValues(std::vector<MRMFeatureQC>& filter_values, const std::vector<FeatureMap>& samples, const MRMFeatureQC& filter

#[test]
fn section_14_accumulate_filter_values() {
    {
        let mut mrmff = MRMFeatureFilter::default();
        let mut samples: Vec<FeatureMap> = Vec::new();
        let mut components = FeatureMap::default();
        let mut component_1 = Feature::default();
        let mut subordinate = Feature::default();
        let mut subordinates: Vec<Feature> = Vec::new();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(500_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.rt = 2.5;
        component_1.intensity = 15000.0_f32;
        component_1.quality = 300.0_f32;
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        samples.push(components.clone());
        components.features.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 3.0;
        subordinate.intensity = 1000.0_f32;
        subordinate.quality = 200.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(1000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 3.0;
        subordinate.intensity = 1000.0_f32;
        subordinate.quality = 400.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(2000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 3.1;
        subordinate.intensity = 2000.0_f32;
        subordinate.quality = 300.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(800_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.rt = 3.0;
        component_1.intensity = 5000.0_f32;
        component_1.quality = 200.0_f32;
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        samples.push(components.clone());
        components.features.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.0;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.0;
        subordinate.intensity = 4000.0_f32;
        subordinate.quality = 150.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(6000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.0;
        subordinate.intensity = 1000.0_f32;
        subordinate.quality = 300.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(100_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.rt = 2.0;
        component_1.intensity = 15000.0_f32;
        component_1.quality = 500.0_f32;
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        samples.push(components.clone());
        components.features.clear();
        let mut transitions = TargetedExperiment::new();
        let mut transition = ReactionMonitoringTransition::new();
        transition.set_native_id("component1.1.Heavy");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component1.1.Light");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(true);
        transitions.add_transition(transition.clone());
        transition.set_native_id("component1.2.Light");
        transition.set_peptide_ref("component_group1");
        transition.set_detecting_transition(true);
        transition.set_identifying_transition(false);
        transition.set_quantifying_transition(false);
        transitions.add_transition(transition.clone());
        let mut filter_values: Vec<MRMFeatureQC> = Vec::new();
        let mut qc_criteria1 = MRMFeatureQC::default();
        let mut qc_criteria2 = MRMFeatureQC::default();
        let mut qc_criteria3 = MRMFeatureQC::default();
        let mut cgqcs = ComponentGroupQCs::default();
        let mut cqcs = ComponentQCs::default();
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.n_heavy_l = 0;
        cgqcs.n_heavy_u = 1;
        cgqcs.n_light_l = 0;
        cgqcs.n_light_u = 2;
        cgqcs.n_detecting_l = 0;
        cgqcs.n_detecting_u = 3;
        cgqcs.n_quantifying_l = 0;
        cgqcs.n_quantifying_u = 2;
        cgqcs.n_identifying_l = 0;
        cgqcs.n_identifying_u = 0;
        cgqcs.n_transitions_l = 0;
        cgqcs.n_transitions_u = 3;
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 10.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 2.5;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 5000.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 100.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 5000.0));
        cqcs.meta_value_qc.insert("peak_area".into(), (500.0, 4e6));
        qc_criteria1.component_group_qcs.push(cgqcs.clone());
        qc_criteria1.component_qcs.push(cqcs.clone());
        filter_values.push(qc_criteria1.clone());
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.n_heavy_l = 0;
        cgqcs.n_heavy_u = 1;
        cgqcs.n_light_l = 0;
        cgqcs.n_light_u = 2;
        cgqcs.n_detecting_l = 0;
        cgqcs.n_detecting_u = 3;
        cgqcs.n_quantifying_l = 0;
        cgqcs.n_quantifying_u = 2;
        cgqcs.n_identifying_l = 0;
        cgqcs.n_identifying_u = 0;
        cgqcs.n_transitions_l = 0;
        cgqcs.n_transitions_u = 3;
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 2.5;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 1000.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 200.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 1000.0));
        cqcs.meta_value_qc.insert("peak_area".into(), (500.0, 4e6));
        qc_criteria2.component_group_qcs.push(cgqcs.clone());
        qc_criteria2.component_qcs.push(cqcs.clone());
        filter_values.push(qc_criteria2.clone());
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.n_heavy_l = 0;
        cgqcs.n_heavy_u = 1;
        cgqcs.n_light_l = 0;
        cgqcs.n_light_u = 2;
        cgqcs.n_detecting_l = 0;
        cgqcs.n_detecting_u = 3;
        cgqcs.n_quantifying_l = 0;
        cgqcs.n_quantifying_u = 2;
        cgqcs.n_identifying_l = 0;
        cgqcs.n_identifying_u = 0;
        cgqcs.n_transitions_l = 0;
        cgqcs.n_transitions_u = 3;
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 60.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 2.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 5000.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 100.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 5000.0));
        cqcs.meta_value_qc.insert("peak_area".into(), (500.0, 4e6));
        qc_criteria3.component_group_qcs.push(cgqcs.clone());
        qc_criteria3.component_qcs.push(cqcs.clone());
        filter_values.push(qc_criteria3.clone());
        let mut filter_values_test: Vec<MRMFeatureQC> = Vec::new();
        mrmff
            .accumulate_filter_values(
                &mut filter_values_test,
                &samples,
                &qc_criteria1,
                &transitions,
            )
            .unwrap();
        for i in 0..filter_values.len() {
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].component_group_name,
                filter_values[i].component_group_qcs[0].component_group_name
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].n_heavy_l,
                filter_values[i].component_group_qcs[0].n_heavy_l
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].n_heavy_u,
                filter_values[i].component_group_qcs[0].n_heavy_u
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].n_light_l,
                filter_values[i].component_group_qcs[0].n_light_l
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].n_light_u,
                filter_values[i].component_group_qcs[0].n_light_u
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].n_detecting_l,
                filter_values[i].component_group_qcs[0].n_detecting_l
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].n_detecting_u,
                filter_values[i].component_group_qcs[0].n_detecting_u
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].n_quantifying_l,
                filter_values[i].component_group_qcs[0].n_quantifying_l
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].n_quantifying_u,
                filter_values[i].component_group_qcs[0].n_quantifying_u
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].n_identifying_l,
                filter_values[i].component_group_qcs[0].n_identifying_l
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].n_identifying_u,
                filter_values[i].component_group_qcs[0].n_identifying_u
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].n_transitions_l,
                filter_values[i].component_group_qcs[0].n_transitions_l
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].n_transitions_u,
                filter_values[i].component_group_qcs[0].n_transitions_u
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].ion_ratio_pair_name_1,
                filter_values[i].component_group_qcs[0].ion_ratio_pair_name_1
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].ion_ratio_pair_name_2,
                filter_values[i].component_group_qcs[0].ion_ratio_pair_name_2
            );
            real_similar(
                filter_values_test[i].component_group_qcs[0].ion_ratio_l,
                filter_values[i].component_group_qcs[0].ion_ratio_l,
            );
            real_similar(
                filter_values_test[i].component_group_qcs[0].ion_ratio_u,
                filter_values[i].component_group_qcs[0].ion_ratio_u,
            );
            assert_eq!(
                filter_values_test[i].component_group_qcs[0].ion_ratio_feature_name,
                filter_values[i].component_group_qcs[0].ion_ratio_feature_name
            );
            assert_eq!(
                filter_values_test[i].component_qcs[0].component_name,
                filter_values[i].component_qcs[0].component_name
            );
            real_similar(
                filter_values_test[i].component_qcs[0].retention_time_l,
                filter_values[i].component_qcs[0].retention_time_l,
            );
            real_similar(
                filter_values_test[i].component_qcs[0].retention_time_u,
                filter_values[i].component_qcs[0].retention_time_u,
            );
            real_similar(
                filter_values_test[i].component_qcs[0].intensity_l,
                filter_values[i].component_qcs[0].intensity_l,
            );
            real_similar(
                filter_values_test[i].component_qcs[0].intensity_u,
                filter_values[i].component_qcs[0].intensity_u,
            );
            real_similar(
                filter_values_test[i].component_qcs[0].overall_quality_l,
                filter_values[i].component_qcs[0].overall_quality_l,
            );
            real_similar(
                filter_values_test[i].component_qcs[0].overall_quality_u,
                filter_values[i].component_qcs[0].overall_quality_u,
            );
            real_similar(
                filter_values_test[i].component_qcs[0].meta_value_qc["peak_apex_int"].0,
                filter_values[i].component_qcs[0].meta_value_qc["peak_apex_int"].0,
            );
            real_similar(
                filter_values_test[i].component_qcs[0].meta_value_qc["peak_apex_int"].1,
                filter_values[i].component_qcs[0].meta_value_qc["peak_apex_int"].1,
            );
            real_similar(
                filter_values_test[i].component_qcs[0].meta_value_qc["peak_area"].0,
                filter_values[i].component_qcs[0].meta_value_qc["peak_area"].0,
            );
            real_similar(
                filter_values_test[i].component_qcs[0].meta_value_qc["peak_area"].1,
                filter_values[i].component_qcs[0].meta_value_qc["peak_area"].1,
            );
        }
    }
}

// START_SECTION(void zeroFilterValues(MRMFeatureQC& filter_zeros, const MRMFeatureQC& filter_template) const)

#[test]
fn section_15_zero_filter_values() {
    {
        let mut mrmff = MRMFeatureFilter::default();
        let mut qc_criteria = MRMFeatureQC::default();
        let mut cgqcs = ComponentGroupQCs::default();
        let mut cqcs = ComponentQCs::default();
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.n_heavy_l = 1;
        cgqcs.n_heavy_u = 1;
        cgqcs.n_light_l = 1;
        cgqcs.n_light_u = 2;
        cgqcs.n_detecting_l = 2;
        cgqcs.n_detecting_u = 3;
        cgqcs.n_quantifying_l = 2;
        cgqcs.n_quantifying_u = 2;
        cgqcs.n_identifying_l = 0;
        cgqcs.n_identifying_u = 3;
        cgqcs.n_transitions_l = 3;
        cgqcs.n_transitions_u = 3;
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.5;
        cgqcs.ion_ratio_u = 2.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4e6;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (500.0, 4e6));
        cqcs.meta_value_qc.insert("peak_area".into(), (500.0, 4e6));
        qc_criteria.component_group_qcs.push(cgqcs.clone());
        qc_criteria.component_qcs.push(cqcs.clone());
        let mut filter_zeros = MRMFeatureQC::default();
        mrmff.zero_filter_values(&mut filter_zeros, &qc_criteria);
        assert_eq!(
            filter_zeros.component_group_qcs[0].component_group_name,
            "component_group1"
        );
        assert_eq!(filter_zeros.component_group_qcs[0].n_heavy_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_heavy_u, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_light_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_light_u, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_detecting_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_detecting_u, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_quantifying_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_quantifying_u, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_identifying_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_identifying_u, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_transitions_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_transitions_u, 0);
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_pair_name_1,
            "component1.1.Light"
        );
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_pair_name_2,
            "component1.2.Light"
        );
        real_similar(filter_zeros.component_group_qcs[0].ion_ratio_l, 0.0);
        real_similar(filter_zeros.component_group_qcs[0].ion_ratio_u, 0.0);
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_feature_name,
            "peak_apex_int"
        );
        assert_eq!(
            filter_zeros.component_qcs[0].component_name,
            "component1.1.Heavy"
        );
        real_similar(filter_zeros.component_qcs[0].retention_time_l, 0.0);
        real_similar(filter_zeros.component_qcs[0].retention_time_u, 0.0);
        real_similar(filter_zeros.component_qcs[0].intensity_l, 0.0);
        real_similar(filter_zeros.component_qcs[0].intensity_u, 0.0);
        real_similar(filter_zeros.component_qcs[0].overall_quality_l, 0.0);
        real_similar(filter_zeros.component_qcs[0].overall_quality_u, 0.0);
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_apex_int"].0,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_apex_int"].1,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_area"].0,
            0.0,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_area"].1,
            0.0,
        );
    }
}

// START_SECTION(void calculateFilterValuesMean(MRMFeatureQC& filter_mean, const std::vector<MRMFeatureQC>& filter_values, const MRMFeatureQC& filter_tem

#[test]
fn section_16_calculate_filter_values_mean_var_and_perc_rsd() {
    {
        let mut mrmff = MRMFeatureFilter::default();
        let mut filter_values: Vec<MRMFeatureQC> = Vec::new();
        let mut qc_criteria1 = MRMFeatureQC::default();
        let mut qc_criteria2 = MRMFeatureQC::default();
        let mut qc_criteria3 = MRMFeatureQC::default();
        let mut cgqcs = ComponentGroupQCs::default();
        let mut cqcs = ComponentQCs::default();
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.retention_time_l = 2.0;
        cgqcs.retention_time_u = 3.0;
        cgqcs.n_heavy_l = 0;
        cgqcs.n_heavy_u = 1;
        cgqcs.n_light_l = 1;
        cgqcs.n_light_u = 2;
        cgqcs.n_detecting_l = 2;
        cgqcs.n_detecting_u = 3;
        cgqcs.n_quantifying_l = 2;
        cgqcs.n_quantifying_u = 2;
        cgqcs.n_identifying_l = 0;
        cgqcs.n_identifying_u = 3;
        cgqcs.n_transitions_l = 3;
        cgqcs.n_transitions_u = 3;
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.5;
        cgqcs.ion_ratio_u = 2.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 2.0;
        cqcs.retention_time_u = 3.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 4.00E+06;
        cqcs.overall_quality_l = 100.0;
        cqcs.overall_quality_u = 500.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (500.0, 4e6));
        cqcs.meta_value_qc.insert("peak_area".into(), (500.0, 4e6));
        qc_criteria1.component_group_qcs.push(cgqcs.clone());
        qc_criteria1.component_qcs.push(cqcs.clone());
        filter_values.push(qc_criteria1.clone());
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.retention_time_l = 1.0;
        cgqcs.retention_time_u = 2.0;
        cgqcs.n_heavy_l = 1;
        cgqcs.n_heavy_u = 1;
        cgqcs.n_light_l = 2;
        cgqcs.n_light_u = 3;
        cgqcs.n_detecting_l = 2;
        cgqcs.n_detecting_u = 6;
        cgqcs.n_quantifying_l = 2;
        cgqcs.n_quantifying_u = 2;
        cgqcs.n_identifying_l = 0;
        cgqcs.n_identifying_u = 3;
        cgqcs.n_transitions_l = 1;
        cgqcs.n_transitions_u = 2;
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.25;
        cgqcs.ion_ratio_u = 3.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 1.0;
        cqcs.retention_time_u = 2.0;
        cqcs.intensity_l = 400.0;
        cqcs.intensity_u = 5.00E+05;
        cqcs.overall_quality_l = 50.0;
        cqcs.overall_quality_u = 700.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (400.0, 5e5));
        cqcs.meta_value_qc.insert("peak_area".into(), (400.0, 5e5));
        qc_criteria2.component_group_qcs.push(cgqcs.clone());
        qc_criteria2.component_qcs.push(cqcs.clone());
        filter_values.push(qc_criteria2.clone());
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.retention_time_l = 1.0;
        cgqcs.retention_time_u = 4.0;
        cgqcs.n_heavy_l = 1;
        cgqcs.n_heavy_u = 2;
        cgqcs.n_light_l = 1;
        cgqcs.n_light_u = 2;
        cgqcs.n_detecting_l = 2;
        cgqcs.n_detecting_u = 3;
        cgqcs.n_quantifying_l = 2;
        cgqcs.n_quantifying_u = 4;
        cgqcs.n_identifying_l = 1;
        cgqcs.n_identifying_u = 3;
        cgqcs.n_transitions_l = 0;
        cgqcs.n_transitions_u = 4;
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.4;
        cgqcs.ion_ratio_u = 2.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 1.0;
        cqcs.retention_time_u = 4.0;
        cqcs.intensity_l = 500.0;
        cqcs.intensity_u = 3.00E+06;
        cqcs.overall_quality_l = 10.0;
        cqcs.overall_quality_u = 600.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (500.0, 3e6));
        cqcs.meta_value_qc.insert("peak_area".into(), (500.0, 3e6));
        qc_criteria3.component_group_qcs.push(cgqcs.clone());
        qc_criteria3.component_qcs.push(cqcs.clone());
        filter_values.push(qc_criteria3.clone());
        let mut filter_zeros = MRMFeatureQC::default();
        mrmff
            .calculate_filter_values_mean(&mut filter_zeros, &filter_values, &qc_criteria1)
            .unwrap();
        assert_eq!(
            filter_zeros.component_group_qcs[0].component_group_name,
            "component_group1"
        );
        real_similar(
            filter_zeros.component_group_qcs[0].retention_time_l,
            1.333333333,
        );
        real_similar(filter_zeros.component_group_qcs[0].retention_time_u, 3.0);
        assert_eq!(
            filter_zeros.component_group_qcs[0].n_heavy_l,
            0.666666667_f64 as i32
        );
        assert_eq!(
            filter_zeros.component_group_qcs[0].n_heavy_u,
            1.333333333_f64 as i32
        );
        assert_eq!(
            filter_zeros.component_group_qcs[0].n_light_l,
            1.333333333_f64 as i32
        );
        assert_eq!(
            filter_zeros.component_group_qcs[0].n_light_u,
            2.333333333_f64 as i32
        );
        assert_eq!(filter_zeros.component_group_qcs[0].n_detecting_l, 2);
        assert_eq!(filter_zeros.component_group_qcs[0].n_detecting_u, 4);
        assert_eq!(filter_zeros.component_group_qcs[0].n_quantifying_l, 2);
        assert_eq!(
            filter_zeros.component_group_qcs[0].n_quantifying_u,
            2.666666667_f64 as i32
        );
        assert_eq!(
            filter_zeros.component_group_qcs[0].n_identifying_l,
            0.333333333_f64 as i32
        );
        assert_eq!(filter_zeros.component_group_qcs[0].n_identifying_u, 3);
        assert_eq!(
            filter_zeros.component_group_qcs[0].n_transitions_l,
            1.333333333_f64 as i32
        );
        assert_eq!(filter_zeros.component_group_qcs[0].n_transitions_u, 3);
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_pair_name_1,
            "component1.1.Light"
        );
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_pair_name_2,
            "component1.2.Light"
        );
        real_similar(filter_zeros.component_group_qcs[0].ion_ratio_l, 0.383333333);
        real_similar(filter_zeros.component_group_qcs[0].ion_ratio_u, 2.333333333);
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_feature_name,
            "peak_apex_int"
        );
        assert_eq!(
            filter_zeros.component_qcs[0].component_name,
            "component1.1.Heavy"
        );
        real_similar(filter_zeros.component_qcs[0].retention_time_l, 1.333333333);
        real_similar(filter_zeros.component_qcs[0].retention_time_u, 3.0);
        real_similar(filter_zeros.component_qcs[0].intensity_l, 466.6666667);
        real_similar(filter_zeros.component_qcs[0].intensity_u, 2500000.0);
        real_similar(filter_zeros.component_qcs[0].overall_quality_l, 53.33333333);
        real_similar(filter_zeros.component_qcs[0].overall_quality_u, 600.0);
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_apex_int"].0,
            466.6666667,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_apex_int"].1,
            2500000.0,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_area"].0,
            466.6666667,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_area"].1,
            2500000.0,
        );
        let mut filter_means = filter_zeros.clone();
        mrmff
            .calculate_filter_values_var(
                &mut filter_zeros,
                &filter_values,
                &filter_means,
                &qc_criteria1,
            )
            .unwrap();
        assert_eq!(
            filter_zeros.component_group_qcs[0].component_group_name,
            "component_group1"
        );
        real_similar(
            filter_zeros.component_group_qcs[0].retention_time_l,
            0.333333333,
        );
        real_similar(filter_zeros.component_group_qcs[0].retention_time_u, 1.0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_heavy_l, 1);
        assert_eq!(filter_zeros.component_group_qcs[0].n_heavy_u, 0);
        assert_eq!(
            filter_zeros.component_group_qcs[0].n_light_l,
            0.333333333_f64 as i32
        );
        assert_eq!(
            filter_zeros.component_group_qcs[0].n_light_u,
            0.333333333_f64 as i32
        );
        assert_eq!(filter_zeros.component_group_qcs[0].n_detecting_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_detecting_u, 3);
        assert_eq!(filter_zeros.component_group_qcs[0].n_quantifying_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_quantifying_u, 2);
        assert_eq!(
            filter_zeros.component_group_qcs[0].n_identifying_l,
            0.333333333_f64 as i32
        );
        assert_eq!(filter_zeros.component_group_qcs[0].n_identifying_u, 0);
        assert_eq!(
            filter_zeros.component_group_qcs[0].n_transitions_l,
            2.333333333_f64 as i32
        );
        assert_eq!(filter_zeros.component_group_qcs[0].n_transitions_u, 1);
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_pair_name_1,
            "component1.1.Light"
        );
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_pair_name_2,
            "component1.2.Light"
        );
        real_similar(filter_zeros.component_group_qcs[0].ion_ratio_l, 0.015833333);
        real_similar(filter_zeros.component_group_qcs[0].ion_ratio_u, 0.333333333);
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_feature_name,
            "peak_apex_int"
        );
        assert_eq!(
            filter_zeros.component_qcs[0].component_name,
            "component1.1.Heavy"
        );
        real_similar(filter_zeros.component_qcs[0].retention_time_l, 0.333333333);
        real_similar(filter_zeros.component_qcs[0].retention_time_u, 1.0);
        real_similar(filter_zeros.component_qcs[0].intensity_l, 3333.333333);
        real_similar(filter_zeros.component_qcs[0].intensity_u, 3.25e12);
        real_similar(filter_zeros.component_qcs[0].overall_quality_l, 2033.333333);
        real_similar(filter_zeros.component_qcs[0].overall_quality_u, 10000.0);
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_apex_int"].0,
            3333.333333,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_apex_int"].1,
            3.25e12,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_area"].0,
            3333.333333,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_area"].1,
            3.25e12,
        );
        let mut filter_vars = filter_zeros.clone();
        mrmff
            .calculate_filter_values_perc_rsd(&mut filter_zeros, &filter_means, &filter_vars)
            .unwrap();
        assert_eq!(
            filter_zeros.component_group_qcs[0].component_group_name,
            "component_group1"
        );
        real_similar(
            filter_zeros.component_group_qcs[0].retention_time_l,
            43.30127019,
        );
        real_similar(
            filter_zeros.component_group_qcs[0].retention_time_u,
            33.33333333,
        );
        assert_eq!(filter_zeros.component_group_qcs[0].n_heavy_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_heavy_u, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_light_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_light_u, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_detecting_l, 0);
        assert_eq!(
            filter_zeros.component_group_qcs[0].n_detecting_u,
            43.30127019_f64 as i32
        );
        assert_eq!(filter_zeros.component_group_qcs[0].n_quantifying_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_quantifying_u, 70);
        assert_eq!(filter_zeros.component_group_qcs[0].n_identifying_l, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_identifying_u, 0);
        assert_eq!(filter_zeros.component_group_qcs[0].n_transitions_l, 141);
        assert_eq!(
            filter_zeros.component_group_qcs[0].n_transitions_u,
            33.33333333_f64 as i32
        );
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_pair_name_1,
            "component1.1.Light"
        );
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_pair_name_2,
            "component1.2.Light"
        );
        real_similar(filter_zeros.component_group_qcs[0].ion_ratio_l, 32.82536711);
        real_similar(filter_zeros.component_group_qcs[0].ion_ratio_u, 24.74358297);
        assert_eq!(
            filter_zeros.component_group_qcs[0].ion_ratio_feature_name,
            "peak_apex_int"
        );
        assert_eq!(
            filter_zeros.component_qcs[0].component_name,
            "component1.1.Heavy"
        );
        real_similar(filter_zeros.component_qcs[0].retention_time_l, 43.30127019);
        real_similar(filter_zeros.component_qcs[0].retention_time_u, 33.33333333);
        real_similar(filter_zeros.component_qcs[0].intensity_l, 12.37179148);
        real_similar(filter_zeros.component_qcs[0].intensity_u, 72.11102551);
        real_similar(filter_zeros.component_qcs[0].overall_quality_l, 84.54843287);
        real_similar(filter_zeros.component_qcs[0].overall_quality_u, 16.66666667);
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_apex_int"].0,
            12.37179148,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_apex_int"].1,
            72.11102551,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_area"].0,
            12.37179148,
        );
        real_similar(
            filter_zeros.component_qcs[0].meta_value_qc["peak_area"].1,
            72.11102551,
        );
    }
}

// START_SECTION(void FilterFeatureMapPercRSD(FeatureMap& features, MRMFeatureQC& filter_criteria, const MRMFeatureQC & filter_values))

#[test]
fn section_17_filter_feature_map_perc_rsd_flag_and_filter() {
    {
        let mut mrmff = MRMFeatureFilter::default();
        let mut components = FeatureMap::default();
        let mut component_1 = Feature::default();
        let mut subordinate = Feature::default();
        let mut subordinates: Vec<Feature> = Vec::new();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(500_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component2.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(1000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component2.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(1000_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group2"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        let mut qc_criteria = MRMFeatureQC::default();
        let mut qc_rsd_values = MRMFeatureQC::default();
        let mut cgqcs = ComponentGroupQCs::default();
        let mut cqcs = ComponentQCs::default();
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 30.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        qc_criteria.component_group_qcs.push(cgqcs.clone());
        qc_criteria.component_qcs.push(cqcs.clone());
        cgqcs.component_group_name = "component_group2".to_owned();
        cgqcs.ion_ratio_pair_name_1 = "component2.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component2.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 30.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component2.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        cqcs.component_name = "component2.1.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        qc_criteria.component_group_qcs.push(cgqcs.clone());
        qc_criteria.component_qcs.push(cqcs.clone());
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 30.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 20.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 20.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 10.0));
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 10.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 10.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 20.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 10.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 100.0));
        qc_rsd_values.component_group_qcs.push(cgqcs.clone());
        qc_rsd_values.component_qcs.push(cqcs.clone());
        cgqcs.component_group_name = "component_group2".to_owned();
        cgqcs.ion_ratio_pair_name_1 = "component2.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component2.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 40.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component2.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        cqcs.component_name = "component2.1.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 10.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        qc_rsd_values.component_group_qcs.push(cgqcs.clone());
        qc_rsd_values.component_qcs.push(cqcs.clone());
        let mut params = Param::new();
        params
            .set_value("flag_or_filter", ParamValue::String("flag".into()), "", &[])
            .unwrap();
        mrmff = MRMFeatureFilter::from_param(&params).unwrap();
        mrmff
            .filter_feature_map_perc_rsd(&mut components, &qc_criteria, &qc_rsd_values)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_message"]
                .as_string_list()
                .unwrap()[0],
            "peak_apex_int"
        );
        assert_eq!(
            components.features[1].metadata["QC_transition_group_%RSD_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[1].metadata["QC_transition_group_%RSD_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[1].metadata["QC_transition_group_%RSD_message"]
                .as_string_list()
                .unwrap()[0],
            "ion_ratio_pair[component2.1.Light/component2.2.Light]"
        );
        assert_eq!(
            components.features[1].subordinates[0].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[1].subordinates[1].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            0.75,
        );
        real_similar(
            components.features[1].metadata["QC_transition_group_%RSD_score"]
                .as_f64()
                .unwrap(),
            0.75,
        );
        real_similar(
            components.features[1].subordinates[0].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[1].subordinates[1].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        params
            .set_value(
                "flag_or_filter",
                ParamValue::String("filter".into()),
                "",
                &[],
            )
            .unwrap();
        mrmff = MRMFeatureFilter::from_param(&params).unwrap();
        mrmff
            .filter_feature_map_perc_rsd(&mut components, &qc_criteria, &qc_rsd_values)
            .unwrap();
        assert_eq!(components.len(), 1);
        assert_eq!(components.features[0].subordinates.len(), 2);
    }
}

// START_SECTION(void FilterFeatureMapPercRSD(FeatureMap& features, MRMFeatureQC& filter_criteria, const MRMFeatureQC & filter_values))

#[test]
fn section_18_filter_feature_map_perc_rsd_individual_checks() {
    {
        let mut mrmff = MRMFeatureFilter::default();
        let mut components = FeatureMap::default();
        let mut component_1 = Feature::default();
        let mut subordinate = Feature::default();
        let mut subordinates: Vec<Feature> = Vec::new();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        let mut qc_criteria = MRMFeatureQC::default();
        let mut qc_rsd_values = MRMFeatureQC::default();
        let mut cgqcs = ComponentGroupQCs::default();
        let mut cqcs = ComponentQCs::default();
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 30.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        qc_criteria.component_group_qcs.push(cgqcs.clone());
        qc_criteria.component_qcs.push(cqcs.clone());
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 30.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 20.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 20.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 10.0));
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 10.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 10.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 20.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 10.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 1.0));
        qc_rsd_values.component_group_qcs.push(cgqcs.clone());
        qc_rsd_values.component_qcs.push(cqcs.clone());
        let mut params = Param::new();
        params
            .set_value("flag_or_filter", ParamValue::String("flag".into()), "", &[])
            .unwrap();
        mrmff = MRMFeatureFilter::from_param(&params).unwrap();
        mrmff
            .filter_feature_map_perc_rsd(&mut components, &qc_criteria, &qc_rsd_values)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        qc_rsd_values.component_group_qcs.clear();
        qc_rsd_values.component_qcs.clear();
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 30.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 20.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 20.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 10.0));
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 10.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 10.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 20.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 80.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 10.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 1.0));
        qc_rsd_values.component_group_qcs.push(cgqcs.clone());
        qc_rsd_values.component_qcs.push(cqcs.clone());
        mrmff
            .filter_feature_map_perc_rsd(&mut components, &qc_criteria, &qc_rsd_values)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_message"]
                .as_string_list()
                .unwrap()[0],
            "retention_time"
        );
        qc_rsd_values.component_group_qcs.clear();
        qc_rsd_values.component_qcs.clear();
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 30.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 20.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 20.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 10.0));
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 10.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 10.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 20.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 100.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 1.0));
        qc_rsd_values.component_group_qcs.push(cgqcs.clone());
        qc_rsd_values.component_qcs.push(cqcs.clone());
        mrmff
            .filter_feature_map_perc_rsd(&mut components, &qc_criteria, &qc_rsd_values)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_message"]
                .as_string_list()
                .unwrap()[0],
            "intensity"
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            0.75,
        );
        qc_rsd_values.component_group_qcs.clear();
        qc_rsd_values.component_qcs.clear();
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 30.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 20.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 20.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 10.0));
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 10.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 10.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 20.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 10.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 300.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 1.0));
        qc_rsd_values.component_group_qcs.push(cgqcs.clone());
        qc_rsd_values.component_qcs.push(cqcs.clone());
        mrmff
            .filter_feature_map_perc_rsd(&mut components, &qc_criteria, &qc_rsd_values)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_message"]
                .as_string_list()
                .unwrap()[0],
            "overall_quality"
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            0.75,
        );
        qc_rsd_values.component_group_qcs.clear();
        qc_rsd_values.component_qcs.clear();
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 30.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 20.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 20.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 10.0));
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 10.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 10.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 20.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 10.0));
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 10.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 100.0));
        qc_rsd_values.component_group_qcs.push(cgqcs.clone());
        qc_rsd_values.component_qcs.push(cqcs.clone());
        mrmff
            .filter_feature_map_perc_rsd(&mut components, &qc_criteria, &qc_rsd_values)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_message"]
                .as_string_list()
                .unwrap()[0],
            "peak_apex_int"
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            0.75,
        );
        qc_rsd_values.component_group_qcs.clear();
        qc_rsd_values.component_qcs.clear();
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.ion_ratio_pair_name_1 = "component1.1.Light".to_owned();
        cgqcs.ion_ratio_pair_name_2 = "component1.2.Light".to_owned();
        cgqcs.ion_ratio_l = 0.0;
        cgqcs.ion_ratio_u = 50.0;
        cgqcs.ion_ratio_feature_name = "peak_apex_int".to_owned();
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 20.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 20.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 10.0));
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 10.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 10.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 20.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 30.0));
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.retention_time_l = 0.0;
        cqcs.retention_time_u = 30.0;
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 10.0;
        cqcs.overall_quality_l = 0.0;
        cqcs.overall_quality_u = 30.0;
        cqcs.meta_value_qc
            .insert("peak_apex_int".into(), (0.0, 1.0));
        qc_rsd_values.component_group_qcs.push(cgqcs.clone());
        qc_rsd_values.component_qcs.push(cqcs.clone());
        mrmff
            .filter_feature_map_perc_rsd(&mut components, &qc_criteria, &qc_rsd_values)
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_%RSD_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[0].metadata["QC_transition_group_%RSD_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[0].metadata["QC_transition_group_%RSD_message"]
                .as_string_list()
                .unwrap()[0],
            "ion_ratio_pair[component1.1.Light/component1.2.Light]"
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_pass"],
            MetaValue::from(1_i32)
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_%RSD_score"]
                .as_f64()
                .unwrap(),
            0.75,
        );
        real_similar(
            components.features[0].subordinates[0].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[1].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(
            components.features[0].subordinates[2].metadata["QC_transition_%RSD_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
    }
}

// START_SECTION(void FilterFeatureMapBackgroundInterference(FeatureMap& features, MRMFeatureQC& filter_criteria, const MRMFeatureQC & filter_values))

#[test]
fn section_19_filter_feature_map_background_interference() {
    {
        let mut mrmff = MRMFeatureFilter::default();
        let mut components = FeatureMap::default();
        let mut component_1 = Feature::default();
        let mut subordinate = Feature::default();
        let mut subordinates: Vec<Feature> = Vec::new();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(5000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component1.2.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(500_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group1"));
        component_1.intensity = 5000.0_f32;
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component2.1.Heavy"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Heavy"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(1000_i32));
        subordinates.push(subordinate.clone());
        subordinate
            .metadata
            .insert("native_id".into(), MetaValue::from("component2.1.Light"));
        subordinate.rt = 2.5;
        subordinate.intensity = 5000.0_f32;
        subordinate.quality = 100.0_f32;
        subordinate
            .metadata
            .insert("LabelType".into(), MetaValue::from("Light"));
        subordinate
            .metadata
            .insert("peak_apex_int".into(), MetaValue::from(1000_i32));
        subordinates.push(subordinate.clone());
        component_1
            .metadata
            .insert("PeptideRef".into(), MetaValue::from("component_group2"));
        component_1.intensity = 5000.0_f32;
        component_1.subordinates = subordinates.clone();
        components.features.push(component_1.clone());
        subordinates.clear();
        let mut qc_criteria = MRMFeatureQC::default();
        let mut qc_background_values = MRMFeatureQC::default();
        let mut cgqcs = ComponentGroupQCs::default();
        let mut cqcs = ComponentQCs::default();
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.intensity_l = 0.0;
        cgqcs.intensity_u = 30.0;
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        qc_criteria.component_group_qcs.push(cgqcs.clone());
        qc_criteria.component_qcs.push(cqcs.clone());
        cgqcs.component_group_name = "component_group2".to_owned();
        cgqcs.intensity_l = 0.0;
        cgqcs.intensity_u = 30.0;
        cqcs.component_name = "component2.1.Heavy".to_owned();
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        cqcs.component_name = "component2.1.Light".to_owned();
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 30.0;
        qc_criteria.component_group_qcs.push(cgqcs.clone());
        qc_criteria.component_qcs.push(cqcs.clone());
        cgqcs.component_group_name = "component_group1".to_owned();
        cgqcs.intensity_l = 0.0;
        cgqcs.intensity_u = 0.0;
        cqcs.component_name = "component1.1.Heavy".to_owned();
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 0.0;
        cqcs.component_name = "component1.1.Light".to_owned();
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 0.0;
        cqcs.component_name = "component1.2.Light".to_owned();
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 10000.0;
        qc_background_values.component_group_qcs.push(cgqcs.clone());
        qc_background_values.component_qcs.push(cqcs.clone());
        cgqcs.component_group_name = "component_group2".to_owned();
        cgqcs.intensity_l = 0.0;
        cgqcs.intensity_u = 10000.0;
        cqcs.component_name = "component2.1.Heavy".to_owned();
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 0.0;
        cqcs.component_name = "component2.1.Light".to_owned();
        cqcs.intensity_l = 0.0;
        cqcs.intensity_u = 0.0;
        qc_background_values.component_group_qcs.push(cgqcs.clone());
        qc_background_values.component_qcs.push(cqcs.clone());
        let mut params = Param::new();
        params
            .set_value("flag_or_filter", ParamValue::String("flag".into()), "", &[])
            .unwrap();
        mrmff = MRMFeatureFilter::from_param(&params).unwrap();
        mrmff
            .filter_feature_map_background_interference(
                &mut components,
                &qc_criteria,
                &qc_background_values,
            )
            .unwrap();
        assert_eq!(
            components.features[0].metadata["QC_transition_group_%BackgroundInterference_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[0].metadata["QC_transition_%BackgroundInterference_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[1].metadata["QC_transition_%BackgroundInterference_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[0].subordinates[2].metadata["QC_transition_%BackgroundInterference_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(components.features[0].subordinates[2].metadata["QC_transition_%BackgroundInterference_message"].as_string_list().unwrap().len(), 1);
        assert_eq!(components.features[0].subordinates[2].metadata["QC_transition_%BackgroundInterference_message"].as_string_list().unwrap()[0], "intensity");
        assert_eq!(
            components.features[1].metadata["QC_transition_group_%BackgroundInterference_pass"],
            MetaValue::from(0_i32)
        );
        assert_eq!(
            components.features[1].metadata["QC_transition_group_%BackgroundInterference_message"]
                .as_string_list()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            components.features[1].metadata["QC_transition_group_%BackgroundInterference_message"]
                .as_string_list()
                .unwrap()[0],
            "intensity"
        );
        assert_eq!(
            components.features[1].subordinates[0].metadata["QC_transition_%BackgroundInterference_pass"],
            MetaValue::from(1_i32)
        );
        assert_eq!(
            components.features[1].subordinates[1].metadata["QC_transition_%BackgroundInterference_pass"],
            MetaValue::from(1_i32)
        );
        real_similar(
            components.features[0].metadata["QC_transition_group_%BackgroundInterference_score"]
                .as_f64()
                .unwrap(),
            1.0,
        );
        real_similar(components.features[0].subordinates[0].metadata["QC_transition_%BackgroundInterference_score"].as_f64().unwrap(), 1.0);
        real_similar(components.features[0].subordinates[1].metadata["QC_transition_%BackgroundInterference_score"].as_f64().unwrap(), 1.0);
        real_similar(components.features[0].subordinates[2].metadata["QC_transition_%BackgroundInterference_score"].as_f64().unwrap(), 0.0);
        real_similar(
            components.features[1].metadata["QC_transition_group_%BackgroundInterference_score"]
                .as_f64()
                .unwrap(),
            0.0,
        );
        real_similar(components.features[1].subordinates[0].metadata["QC_transition_%BackgroundInterference_score"].as_f64().unwrap(), 1.0);
        real_similar(components.features[1].subordinates[1].metadata["QC_transition_%BackgroundInterference_score"].as_f64().unwrap(), 1.0);
        params
            .set_value(
                "flag_or_filter",
                ParamValue::String("filter".into()),
                "",
                &[],
            )
            .unwrap();
        mrmff = MRMFeatureFilter::from_param(&params).unwrap();
        mrmff
            .filter_feature_map_background_interference(
                &mut components,
                &qc_criteria,
                &qc_background_values,
            )
            .unwrap();
        assert_eq!(components.len(), 1);
        assert_eq!(components.features[0].subordinates.len(), 2);
    }
}

// START_SECTION(double calculateRTDifference(Feature & component_1, Feature & component_2))

#[test]
fn section_20_calculate_rt_difference() {
    {
        let mut mrmff = MRMFeatureFilter::default();
        let mut f1 = Feature::default();
        let mut f2 = Feature::default();
        f1.rt = 100.0;
        f2.rt = 105.5;
        let mut rt_diff = mrmff.calculate_rt_difference(&f1, &f2);
        real_similar(rt_diff, 5.5);
        rt_diff = mrmff.calculate_rt_difference(&f2, &f1);
        real_similar(rt_diff, 5.5);
    }
}

// START_SECTION(double calculateResolution(Feature & component_1, Feature & component_2))

#[test]
fn section_21_calculate_resolution() {
    {
        let mut mrmff = MRMFeatureFilter::default();
        let mut f1 = Feature::default();
        let mut f2 = Feature::default();
        f1.rt = 100.0;
        f1.metadata
            .insert("width_at_50".into(), MetaValue::try_from(2.0_f64).unwrap());
        f2.rt = 106.0;
        f2.metadata
            .insert("width_at_50".into(), MetaValue::try_from(2.0_f64).unwrap());
        let mut resolution = mrmff.calculate_resolution(&f1, &f2).unwrap();
        real_similar(resolution, 1.7647);
        let mut f3 = Feature::default();
        let mut f4 = Feature::default();
        f3.rt = 100.0;
        f4.rt = 106.0;
        resolution = mrmff.calculate_resolution(&f3, &f4).unwrap();
        real_similar(resolution, 0.0);
    }
}
