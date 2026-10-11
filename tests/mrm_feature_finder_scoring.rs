// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! MRMFeatureFinderScoring.
//!
//! Every section of the upstream class test at the core pin `bc9cc12`
//! (`MRMFeatureFinderScoring_test.cpp`, evidence tier 3) is transcribed, with
//! its tolerances (`TOLERANCE_ABSOLUTE(0.1)`, and `TOLERANCE_RELATIVE(1.001)`
//! for `rt_score`; a pair is similar when either bound holds).
//!
//! The class test loads `OpenSwath_generic_input.TraML` and
//! `OpenSwath_identification_input.TraML` with `TraMLFile`, which is not
//! ported (FORMAT/TraMLFile.h is a backlog item). Their content is built here
//! as a [`TargetedExperiment`] by `generic_transitions()` and
//! `identification_transitions()`, element for element from the two files
//! shipped in `tests/data`, and converted with `convert_targeted_exp`; the mzML
//! is read with the crate's reader.
#![cfg(feature = "mzml")]

use openms::Error;
use openms::analysis::openswath::data_access_helper::convert_targeted_exp;
use openms::analysis::openswath::mrm_feature_finder_scoring::{
    MRMFeatureFinderScoring, MRMTransitionGroupType, RUN_IDENTIFIER, TransitionGroupMapType,
};
use openms::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition;
use openms::analysis::openswath::spectrum_access_openms::SimpleOpenMSSpectraFactory;
use openms::analysis::openswath::swath_map::SwathMap;
use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
use openms::analysis::openswath::targeted_experiment_helper::{
    Peptide, Protein, RTType, RetentionTime,
};
use openms::analysis::transformations::TransformationDescription;
use openms::concept::unique_id::UniqueIdGenerator;
use openms::kernel::{FeatureMap, MSExperiment};
use openms::param::ParamValue;
use std::path::PathBuf;
use std::sync::Arc;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// `TEST_REAL_SIMILAR` under the given `TOLERANCE_ABSOLUTE` and
/// `TOLERANCE_RELATIVE` (ratio bound): similar when either bound holds.
fn similar(a: f64, b: f64, absolute: f64, relative: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= absolute || (ratio - 1.0).abs() <= relative || a == b,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

fn near(a: f64, b: f64) {
    similar(a, b, 0.1, 1e-5);
}

fn peptide(id: &str, sequence: &str, normalized_rt: f64) -> Peptide {
    let mut p = Peptide::default();
    p.id = id.into();
    p.sequence = sequence.into();
    p.protein_refs.push(String::new());
    let mut rt = RetentionTime::default();
    rt.set_rt(normalized_rt);
    rt.retention_time_type = RTType::Normalized;
    p.rts.push(rt);
    p
}

fn transition(
    id: &str,
    peptide_ref: &str,
    precursor_mz: f64,
    product_mz: f64,
    intensity: f64,
    flags: Option<(bool, bool, bool)>,
) -> ReactionMonitoringTransition {
    let mut t = ReactionMonitoringTransition::new();
    t.set_native_id(id);
    t.set_peptide_ref(peptide_ref);
    t.set_precursor_mz(precursor_mz);
    t.set_product_mz(product_mz);
    t.set_library_intensity(intensity);
    if let Some((detecting, identifying, quantifying)) = flags {
        t.set_detecting_transition(detecting);
        t.set_identifying_transition(identifying);
        t.set_quantifying_transition(quantifying);
    }
    t
}

fn proteins_and_peptides() -> TargetedExperiment {
    let mut exp = TargetedExperiment::new();
    exp.add_protein(Protein::default());
    exp.add_peptide(peptide("tr_gr1", "PEPTIDEA", 0.44));
    exp.add_peptide(peptide("tr_gr2", "PEPTIDEB", 0.2));
    exp
}

/// `OpenSwath_generic_input.TraML`.
fn generic_transitions() -> TargetedExperiment {
    let mut exp = proteins_and_peptides();
    for t in [
        transition("tr1", "tr_gr1", 500.0, 628.435, 1.0, None),
        transition("tr2", "tr_gr1", 500.0, 654.38, 2.0, None),
        transition("tr3", "tr_gr2", 501.0, 618.31, 10000.0, None),
        transition("tr4", "tr_gr2", 501.0, 628.435, 2000.0, None),
        transition("tr5", "tr_gr2", 501.0, 651.3, 4300.0, None),
    ] {
        exp.add_transition(t);
    }
    exp
}

/// `OpenSwath_identification_input.TraML`: the same assays with the
/// detecting, identifying and quantifying flags, and `tr2` once more in
/// `tr_gr2`, as an identifying but not detecting transition.
fn identification_transitions() -> TargetedExperiment {
    let mut exp = proteins_and_peptides();
    let on = Some((true, false, true));
    for t in [
        transition("tr1", "tr_gr1", 500.0, 628.435, 1.0, on),
        transition("tr2", "tr_gr1", 500.0, 654.38, 2.0, on),
        transition("tr3", "tr_gr2", 501.0, 618.31, 10000.0, on),
        transition("tr4", "tr_gr2", 501.0, 628.435, 2000.0, on),
        transition(
            "tr5",
            "tr_gr2",
            501.0,
            651.3,
            4300.0,
            Some((true, true, true)),
        ),
        transition(
            "tr2",
            "tr_gr2",
            500.0,
            654.38,
            2.0,
            Some((false, true, true)),
        ),
    ] {
        exp.add_transition(t);
    }
    exp
}

fn load_chromatograms() -> Arc<MSExperiment> {
    let bytes = std::fs::read(data("OpenSwath_generic_input.mzML")).unwrap();
    Arc::new(openms::format::mzml::read(bytes.as_slice()).unwrap())
}

#[test]
fn constructor_and_defaults() {
    // MRMFeatureFinderScoring() and the source's parameter defaults.
    let ff = MRMFeatureFinderScoring::new().unwrap();
    assert_eq!(ff.stop_report_after_feature, -1);
    assert_eq!(ff.rt_extraction_window, -1.0);
    assert_eq!(ff.rt_normalization_factor, 1.0);
    assert!(ff.strict);
    assert_eq!(ff.enzyme, "Trypsin");
    assert_eq!(ff.sn_bin_count, 30);
    let defaults = MRMFeatureFinderScoring::defaults().unwrap();
    assert_eq!(
        defaults.value("Scores:use_total_mi_score").unwrap(),
        &ParamValue::String("false".into())
    );
    assert!(
        defaults
            .exists("TransitionGroupPicker:PeakPickerChromatogram:sn_win_len")
            .unwrap()
    );
    assert!(defaults.exists("DIAScoring:dia_extraction_window").unwrap());
}

fn pick(
    ff: &mut MRMFeatureFinderScoring,
    transitions: &TargetedExperiment,
) -> (FeatureMap, TransitionGroupMapType) {
    let exp = load_chromatograms();
    let light = convert_targeted_exp(transitions).unwrap();
    let chromatogram_ptr = SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(exp).unwrap();
    let swath_ptr =
        SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(Arc::new(MSExperiment::default()))
            .unwrap();
    let swath_maps = [SwathMap {
        sptr: Some(swath_ptr),
        ..SwathMap::default()
    }];
    let mut feature_file = FeatureMap::default();
    let mut transition_group_map = TransitionGroupMapType::new();
    ff.set_unique_id_generator(UniqueIdGenerator::from_seed(42));
    ff.pick_experiment_light(
        &chromatogram_ptr,
        &mut feature_file,
        &light,
        &TransformationDescription::default(),
        &swath_maps,
        &mut transition_group_map,
    )
    .unwrap();
    (feature_file, transition_group_map)
}

#[test]
fn pick_experiment_default_and_legacy_parameters() {
    // START_SECTION(void pickExperiment(...)), first section.
    let mut ff = MRMFeatureFinderScoring::new().unwrap();
    let (feature_file, transition_group_map) = pick(&mut ff, &generic_transitions());

    // the number of groups found
    assert_eq!(transition_group_map.len(), 2);

    //// scores of the first group
    let group = &transition_group_map["tr_gr1"];
    assert_eq!(group.size(), 2);
    assert_eq!(group.features().len(), 1);
    let feature = &group.features()[0];
    near(feature.feature.base.rt, 3119.092);
    near(f64::from(feature.feature.base.intensity), 3599.98);
    let meta = &feature.feature.base.metadata;
    near(meta["leftWidth"].as_f64().unwrap(), 3092.85009765625);
    near(meta["rightWidth"].as_f64().unwrap(), 3151.10009765625);
    near(meta["total_xic"].as_f64().unwrap(), 3680.16);
    // feature scores
    near(meta["var_xcorr_coelution"].as_f64().unwrap(), 0.0);
    near(meta["var_xcorr_shape"].as_f64().unwrap(), 0.9981834605);
    near(meta["var_library_rmsd"].as_f64().unwrap(), 0.108663236);
    near(meta["var_library_corr"].as_f64().unwrap(), 1.0);
    near(
        meta["var_elution_model_fit_score"].as_f64().unwrap(),
        0.9854,
    );
    near(meta["var_intensity_score"].as_f64().unwrap(), 0.971);
    near(meta["sn_ratio"].as_f64().unwrap(), 65.2167155356372);
    near(meta["var_log_sn_score"].as_f64().unwrap(), 4.17771580932057);
    similar(meta["rt_score"].as_f64().unwrap(), 3118.651968, 0.1, 1e-3);

    //// scores of the second group
    let group = &transition_group_map["tr_gr2"];
    assert_eq!(group.size(), 3);
    assert_eq!(group.features().len(), 2);
    assert_eq!(feature_file.len(), 3);
    let feature = &group.features()[0];
    near(feature.feature.base.rt, 3119.092);
    near(f64::from(feature.feature.base.intensity), 1049.65);
    let meta = &feature.feature.base.metadata;
    near(meta["leftWidth"].as_f64().unwrap(), 3092.85009765625);
    near(meta["rightWidth"].as_f64().unwrap(), 3147.67993164062);
    near(meta["total_xic"].as_f64().unwrap(), 1610.27);
    near(meta["var_xcorr_coelution"].as_f64().unwrap(), 7.0);
    near(meta["var_xcorr_shape"].as_f64().unwrap(), 0.7245);
    near(meta["var_library_rmsd"].as_f64().unwrap(), 0.43566);
    near(meta["var_library_corr"].as_f64().unwrap(), -0.784);
    near(meta["var_elution_model_fit_score"].as_f64().unwrap(), 0.902);
    near(meta["var_intensity_score"].as_f64().unwrap(), 0.642);
    near(meta["sn_ratio"].as_f64().unwrap(), 22.9654496908248);
    near(meta["var_log_sn_score"].as_f64().unwrap(), 3.13399089915908);

    // the protein run the finder writes into the output
    assert_eq!(feature_file.protein_identifications.len(), 1);
    assert_eq!(
        feature_file.protein_identifications[0].identifier,
        RUN_IDENTIFIER
    );

    // test legacy parameters
    let mut picker_param = MRMFeatureFinderScoring::defaults().unwrap();
    picker_param
        .set_value(
            "TransitionGroupPicker:PeakPickerChromatogram:method",
            ParamValue::String("legacy".into()),
            "",
            &[],
        )
        .unwrap();
    picker_param
        .set_value(
            "TransitionGroupPicker:PeakPickerChromatogram:peak_width",
            ParamValue::Float(40.0),
            "",
            &[],
        )
        .unwrap();
    ff.set_parameters(&picker_param).unwrap();
    let (feature_file, transition_group_map) = pick(&mut ff, &generic_transitions());

    let group = &transition_group_map["tr_gr1"];
    assert_eq!(group.size(), 2);
    assert_eq!(group.features().len(), 1);
    let feature = &group.features()[0];
    near(feature.feature.base.rt, 3119.092);
    near(f64::from(feature.feature.base.intensity), 3574.23);
    let meta = &feature.feature.base.metadata;
    near(meta["leftWidth"].as_f64().unwrap(), 3099.69995117188);
    near(meta["rightWidth"].as_f64().unwrap(), 3147.68);
    near(meta["total_xic"].as_f64().unwrap(), 3680.16);

    // the more accurate (interpolated) S/N estimate resolves a small nearby
    // peak as its own feature, so this group yields 3 features
    let group = &transition_group_map["tr_gr2"];
    assert_eq!(group.size(), 3);
    assert_eq!(group.features().len(), 3);
    assert_eq!(feature_file.len(), 4);
    let feature = &group.features()[0];
    near(feature.feature.base.rt, 3119.092);
    near(f64::from(feature.feature.base.intensity), 1034.55);
    let meta = &feature.feature.base.metadata;
    near(meta["leftWidth"].as_f64().unwrap(), 3099.7);
    near(meta["rightWidth"].as_f64().unwrap(), 3144.25);
    near(meta["total_xic"].as_f64().unwrap(), 1610.27);
    near(meta["var_xcorr_coelution"].as_f64().unwrap(), 2.265);
}

#[test]
fn pick_experiment_with_identification_scores() {
    // START_SECTION(void pickExperiment(...)), second section: UIS scores.
    let mut param = MRMFeatureFinderScoring::defaults().unwrap();
    param
        .set_value(
            "Scores:use_uis_scores",
            ParamValue::String("true".into()),
            "",
            &[],
        )
        .unwrap();
    param
        .set_value(
            "TransitionGroupPicker:PeakPickerChromatogram:method",
            ParamValue::String("legacy".into()),
            "",
            &[],
        )
        .unwrap();
    param
        .set_value(
            "TransitionGroupPicker:PeakPickerChromatogram:peak_width",
            ParamValue::Float(40.0),
            "",
            &[],
        )
        .unwrap();
    let mut ff = MRMFeatureFinderScoring::from_param(&param).unwrap();
    let (feature_file, transition_group_map) = pick(&mut ff, &identification_transitions());

    assert_eq!(transition_group_map.len(), 2);
    let group = &transition_group_map["tr_gr2"];
    assert_eq!(group.size(), 3);
    assert_eq!(group.features().len(), 3);
    assert_eq!(feature_file.len(), 4);
    let feature = &group.features()[0];
    near(feature.feature.base.rt, 3119.092);
    near(f64::from(feature.feature.base.intensity), 1034.55);
    let meta = &feature.feature.base.metadata;
    near(meta["leftWidth"].as_f64().unwrap(), 3099.7);
    near(meta["rightWidth"].as_f64().unwrap(), 3144.25);
    near(meta["total_xic"].as_f64().unwrap(), 1610.27);
    // feature scores
    near(meta["var_xcorr_coelution"].as_f64().unwrap(), 2.265);
    near(meta["var_xcorr_shape"].as_f64().unwrap(), 0.7245);
    near(meta["var_library_rmsd"].as_f64().unwrap(), 0.43566);
    near(meta["var_library_corr"].as_f64().unwrap(), -0.784);
    near(meta["var_elution_model_fit_score"].as_f64().unwrap(), 0.902);
    near(
        meta["var_intensity_score"].as_f64().unwrap(),
        2.36342573991536,
    );
    near(meta["sn_ratio"].as_f64().unwrap(), 22.9654496908248);
    near(meta["var_log_sn_score"].as_f64().unwrap(), 3.13399089915908);
    // feature identification scores
    let names = meta["id_target_transition_names"].as_string_list().unwrap();
    assert_eq!(names[0], "tr5");
    assert_eq!(names[1], "tr2");
    let list = |key: &str| meta[key].as_float_list().unwrap().to_vec();
    near(list("id_target_ind_log_intensity")[0], 5.03352);
    near(list("id_target_ind_log_intensity")[1], 7.92704);
    near(meta["id_target_num_transitions"].as_f64().unwrap(), 2.0);
    near(list("id_target_ind_xcorr_coelution")[0], 1.0);
    near(list("id_target_ind_xcorr_coelution")[1], 1.66667);
    near(list("id_target_ind_xcorr_shape")[0], 0.68631);
    near(list("id_target_ind_xcorr_shape")[1], 0.690494);
    near(list("id_target_ind_log_sn_score")[0], 1.16692);
    near(list("id_target_ind_log_sn_score")[1], 4.18437274576997);
    assert_eq!(list("id_target_ind_isotope_correlation").len(), 0);
    assert_eq!(list("id_target_ind_isotope_overlap").len(), 0);
    assert_eq!(list("id_target_ind_massdev_score").len(), 0);
}

#[test]
fn map_experiment_to_transition_list() {
    // START_SECTION(void mapExperimentToTransitionList(...))
    let mut ff = MRMFeatureFinderScoring::new().unwrap();
    let light = convert_targeted_exp(&generic_transitions()).unwrap();
    let chromatogram_ptr =
        SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(load_chromatograms()).unwrap();
    ff.prepare_protein_peptide_maps(&light);
    let mut transition_group_map = TransitionGroupMapType::new();
    ff.map_experiment_to_transition_list(
        &chromatogram_ptr,
        &light,
        &mut transition_group_map,
        TransformationDescription::default(),
        -1.0,
    )
    .unwrap();
    assert_eq!(transition_group_map.len(), 2);

    // the first group
    let group = &transition_group_map["tr_gr1"];
    assert_eq!(group.size(), 2);
    assert_eq!(group.transitions().len(), 2);
    assert_eq!(group.chromatograms().len(), 2);
    assert!(group.has_chromatogram("tr1"));
    assert!(group.has_chromatogram("tr2"));
    assert_eq!(group.chromatogram("tr2").unwrap().native_id, "tr2");
    assert_eq!(group.transition("tr2").unwrap().native_id(), "tr2");

    // the second group
    let group: &MRMTransitionGroupType = &transition_group_map["tr_gr2"];
    assert_eq!(group.size(), 3);
    assert_eq!(group.transitions().len(), 3);
    assert_eq!(group.chromatograms().len(), 3);
    assert!(group.has_chromatogram("tr3"));
    assert!(group.has_chromatogram("tr4"));
    assert!(group.has_chromatogram("tr5"));
    assert_eq!(group.chromatogram("tr5").unwrap().native_id, "tr5");
    assert_eq!(group.transition("tr5").unwrap().native_id(), "tr5");
}

#[test]
fn native_strict_mapping_and_window() {
    // tier 4: a transition without chromatogram is an error when strict and
    // skipped otherwise; an RT window cuts the chromatograms and needs the
    // compound map (MRMFeatureFinderScoring.cpp mapExperimentToTransitionList).
    let mut ff = MRMFeatureFinderScoring::new().unwrap();
    let mut exp = generic_transitions();
    exp.add_transition(transition("missing", "tr_gr1", 500.0, 700.0, 5.0, None));
    let light = convert_targeted_exp(&exp).unwrap();
    let ptr = SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(load_chromatograms()).unwrap();
    let mut map = TransitionGroupMapType::new();
    let err = ff
        .map_experiment_to_transition_list(
            &ptr,
            &light,
            &mut map,
            TransformationDescription::default(),
            -1.0,
        )
        .unwrap_err();
    assert!(
        matches!(err, Error::InvalidValue(m) if m.contains("does not have a corresponding chromatogram"))
    );
    ff.set_strict_flag(false);
    let mut map = TransitionGroupMapType::new();
    ff.map_experiment_to_transition_list(
        &ptr,
        &light,
        &mut map,
        TransformationDescription::default(),
        -1.0,
    )
    .unwrap();
    assert_eq!(map["tr_gr1"].size(), 2);

    // an RT window around the expected RT keeps only nearby points: the
    // normalised RT 0.44 maps to itself under the identity trafo, so the
    // window [-9.56, 10.44] holds none of the points near RT 3100
    ff.prepare_protein_peptide_maps(&light);
    let mut windowed = TransitionGroupMapType::new();
    ff.map_experiment_to_transition_list(
        &ptr,
        &light,
        &mut windowed,
        TransformationDescription::default(),
        10.0,
    )
    .unwrap(); // not strict: empty traces are kept
    assert!(windowed["tr_gr1"].chromatograms()[0].peaks.is_empty());
    ff.set_strict_flag(true);
    let complete = convert_targeted_exp(&generic_transitions()).unwrap();
    let mut windowed = TransitionGroupMapType::new();
    let err = ff
        .map_experiment_to_transition_list(
            &ptr,
            &complete,
            &mut windowed,
            TransformationDescription::default(),
            10.0,
        )
        .unwrap_err();
    assert!(matches!(err, Error::InvalidValue(m) if m.contains("Could not find any points")));
}

#[test]
fn native_score_peakgroups_requires_the_maps() {
    // tier 4: scorePeakgroups' preconditions.
    let ff = MRMFeatureFinderScoring::new().unwrap();
    let mut group = MRMTransitionGroupType::new();
    let mut out = FeatureMap::default();
    let err = ff
        .score_peakgroups(
            &mut group,
            &TransformationDescription::default(),
            &[],
            &mut out,
            false,
            None,
        )
        .unwrap_err();
    assert!(matches!(err, Error::InvalidValue(m) if m.contains("Peptide reference map is empty")));
    let mut ff = ff;
    ff.prepare_protein_peptide_maps(&convert_targeted_exp(&generic_transitions()).unwrap());
    let err = ff
        .score_peakgroups(
            &mut group,
            &TransformationDescription::default(),
            &[],
            &mut out,
            false,
            None,
        )
        .unwrap_err();
    assert!(matches!(err, Error::InvalidValue(m) if m.contains("Transition group id is empty")));
}
