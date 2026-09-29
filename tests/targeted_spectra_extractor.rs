// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `ANALYSIS/OPENSWATH/TargetedSpectraExtractor`.
//!
//! Expected values are the literals of `TargetedSpectraExtractor_test.cpp` at
//! core `bc9cc12` unless a test says otherwise; see
//! `docs/TARGETED_SPECTRA_EXTRACTOR_SUPPORT.md` for the sections that need
//! unported readers (`TransitionTSVFile`, `MSPGenericFile`) or engines
//! (`AccurateMassSearchEngine`).

// The class test's literals are transcribed verbatim, including digits beyond
// f64 precision.
#![allow(clippy::excessive_precision)]

use openms::analysis::openswath::targeted_spectra_extractor::{
    BinnedSpectrumComparator, Comparator, TargetedSpectraExtractor,
};
use openms::concept::unique_id::UniqueIdGenerator;
use openms::kernel::features::{BaseFeature, Feature, FeatureMap};
use openms::kernel::{MSExperiment, MSSpectrum, Peak1D, Precursor};
use openms::metadata::MetaValue;
use openms::param::{Param, ParamValue};
use std::collections::BTreeMap;

fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

fn set(param: &mut Param, key: &str, value: ParamValue) {
    let description = param.description(key).unwrap_or("").to_owned();
    param.set_value(key, value, &description, &[]).unwrap();
}

/// The parameters most class-test sections set.
fn class_test_extractor(extra: &[(&str, ParamValue)]) -> TargetedSpectraExtractor {
    let mut params = TargetedSpectraExtractor::new().unwrap().to_param().unwrap();
    set(
        &mut params,
        "GaussFilter:gaussian_width",
        ParamValue::Float(0.25),
    );
    set(&mut params, "peak_height_min", ParamValue::Float(15000.0));
    set(&mut params, "peak_height_max", ParamValue::Float(110000.0));
    set(&mut params, "fwhm_threshold", ParamValue::Float(0.23));
    for (key, value) in extra {
        set(&mut params, key, value.clone());
    }
    TargetedSpectraExtractor::from_param(&params).unwrap()
}

fn meta_str(feature: &Feature, key: &str) -> String {
    feature.metadata[key].as_str().unwrap().to_owned()
}

fn meta_f64(feature: &Feature, key: &str) -> f64 {
    feature.metadata[key].as_f64().unwrap()
}

#[test]
fn parameters() {
    let tse = TargetedSpectraExtractor::new().unwrap();
    let params = tse.to_param().unwrap();
    let float = |key: &str| params.value(key).unwrap().to_f64().unwrap();
    let text = |key: &str| params.value(key).unwrap().as_str().unwrap().to_owned();
    let int = |key: &str| params.value(key).unwrap().to_i64().unwrap();
    assert_eq!(float("rt_window"), 30.0);
    assert_eq!(float("min_select_score"), 0.7);
    assert_eq!(float("mz_tolerance"), 0.1);
    assert_eq!(text("mz_unit_is_Da"), "true");
    assert_eq!(int("SavitzkyGolayFilter:frame_length"), 15);
    assert_eq!(int("SavitzkyGolayFilter:polynomial_order"), 3);
    assert_eq!(float("GaussFilter:gaussian_width"), 0.2);
    assert_eq!(text("use_gauss"), "true");
    assert_eq!(float("PeakPickerHiRes:signal_to_noise"), 1.0);
    assert_eq!(float("peak_height_min"), 0.0);
    assert_eq!(float("peak_height_max"), f64::MAX);
    assert_eq!(float("fwhm_threshold"), 0.0);
    assert_eq!(float("tic_weight"), 1.0);
    assert_eq!(float("fwhm_weight"), 1.0);
    assert_eq!(float("snr_weight"), 1.0);
    assert_eq!(int("top_matches_to_report"), 5);
    assert_eq!(float("min_match_score"), 0.8);
    // Source: the constructor overrides AccurateMassSearchEngine's own default.
    assert_eq!(
        text("AccurateMassSearchEngine:keep_unidentified_masses"),
        "false"
    );

    // getDefaultParameters: the class's own entries only.
    let defaults = TargetedSpectraExtractor::default_parameters().unwrap();
    assert_eq!(defaults.value("rt_window").unwrap().to_f64().unwrap(), 30.0);
    assert_eq!(
        defaults
            .value("top_matches_to_report")
            .unwrap()
            .to_i64()
            .unwrap(),
        5
    );
    assert_eq!(
        defaults.value("min_match_score").unwrap().to_f64().unwrap(),
        0.8
    );
    assert!(!defaults.exists("GaussFilter:gaussian_width").unwrap());

    // Round trip, and a restriction violation is refused.
    let again = TargetedSpectraExtractor::from_param(&params).unwrap();
    assert!(again.to_param().unwrap().source_equal(&params).unwrap());
    let mut bad = params.clone();
    set(&mut bad, "min_match_score", ParamValue::Float(1.5));
    assert!(TargetedSpectraExtractor::from_param(&bad).is_err());
}

const MZ: [f64; 73] = [
    61.92, 68.88, 71.4, 79.56, 84.6, 84.72, 84.84, 84.96, 85.08, 85.2, 85.32, 85.44, 85.68, 85.8,
    85.92, 86.04, 86.16, 86.28, 86.4, 87.72, 87.96, 88.08, 90.36, 94.44, 99.84, 100.8, 101.04,
    101.88, 102.0, 102.96, 110.16, 110.88, 111.0, 111.12, 111.24, 111.84, 111.96, 112.08, 112.2,
    112.32, 112.44, 112.56, 112.68, 114.0, 128.16, 128.4, 128.88, 129.0, 129.12, 129.84, 129.96,
    130.08, 130.2, 130.32, 130.44, 130.56, 132.12, 138.0, 139.08, 140.16, 144.12, 146.04, 146.16,
    156.0, 156.12, 156.36, 173.76, 174.0, 174.12, 174.24, 174.36, 174.6, 175.08,
];
const INTENSITY: [f64; 73] = [
    6705.41660838088,
    1676.35415209522,
    1676.35415209522,
    1676.35415209522,
    3352.70830419044,
    5029.06245628566,
    8381.7707604761,
    53643.332867047,
    51966.9787149518,
    6705.41660838088,
    8381.7707604761,
    1676.35415209522,
    11734.4790646665,
    25145.3122814283,
    68730.520235904,
    112315.72819038,
    6705.41660838088,
    6705.41660838088,
    3352.70830419044,
    1676.35415209522,
    1676.35415209522,
    1676.35415209522,
    3352.70830419044,
    1676.35415209522,
    1676.35415209522,
    1676.35415209522,
    5029.06245628566,
    3352.70830419044,
    3352.70830419044,
    3352.70830419044,
    1676.35415209522,
    5029.06245628566,
    3352.70830419044,
    5029.06245628566,
    3352.70830419044,
    5029.06245628566,
    18439.8956730474,
    20116.2498251426,
    5029.06245628566,
    1676.35415209522,
    1676.35415209522,
    3352.70830419044,
    3352.70830419044,
    3352.70830419044,
    6705.41660838088,
    1676.35415209522,
    3352.70830419044,
    3352.70830419044,
    6705.41660838088,
    5029.06245628566,
    10058.1249125713,
    31850.7288898092,
    10058.1249125713,
    1676.35415209522,
    1676.35415209522,
    3352.70830419044,
    1676.35415209522,
    1676.35415209522,
    1676.35415209522,
    3352.70830419044,
    1676.35415209522,
    3352.70830419044,
    1676.35415209522,
    1676.35415209522,
    5029.06245628566,
    1676.35415209522,
    1676.35415209522,
    1676.35415209522,
    6705.41660838088,
    11734.4790646665,
    6705.41660838088,
    1676.35415209522,
    1676.35415209522,
];

/// The class test's C-aconitate product-ion spectrum (QTRAP 5500).
fn class_spectrum() -> MSSpectrum {
    let mut spectrum = MSSpectrum::default();
    for (mz, intensity) in MZ.iter().zip(INTENSITY) {
        spectrum.peaks.push(Peak1D::new(*mz, intensity as f32));
    }
    spectrum
}

#[test]
fn pick_spectrum() {
    let mut params = TargetedSpectraExtractor::new().unwrap().to_param().unwrap();
    set(
        &mut params,
        "GaussFilter:gaussian_width",
        ParamValue::Float(0.25),
    );
    set(&mut params, "peak_height_min", ParamValue::Float(0.0));
    set(&mut params, "peak_height_max", ParamValue::Float(200000.0));
    set(&mut params, "fwhm_threshold", ParamValue::Float(0.0));
    let tse = TargetedSpectraExtractor::from_param(&params).unwrap();
    let spectrum = class_spectrum();

    let picked = tse.pick_spectrum(&spectrum).unwrap();
    assert_ne!(spectrum.len(), picked.len());
    assert_eq!(picked.len(), 6);
    let expected = [
        (85.014, 60754.7),
        (86.0196, 116036.0),
        (112.033, 21941.9),
        (129.396, 10575.5),
        (130.081, 31838.1),
        (174.24, 11731.3),
    ];
    for (peak, (mz, intensity)) in picked.peaks.iter().zip(expected) {
        real_similar(peak.mz, mz);
        real_similar(f64::from(peak.intensity), intensity);
    }

    set(&mut params, "peak_height_min", ParamValue::Float(15000.0));
    set(&mut params, "peak_height_max", ParamValue::Float(110000.0));
    let tse = TargetedSpectraExtractor::from_param(&params).unwrap();
    let picked = tse.pick_spectrum(&spectrum).unwrap();
    // With the new filters on peaks' heights, fewer peaks get picked.
    assert_eq!(picked.len(), 3);
    for (peak, (mz, intensity)) in
        picked
            .peaks
            .iter()
            .zip([(85.014, 60754.7), (112.033, 21941.9), (130.081, 31838.1)])
    {
        real_similar(peak.mz, mz);
        real_similar(f64::from(peak.intensity), intensity);
    }

    set(&mut params, "fwhm_threshold", ParamValue::Float(0.23));
    let tse = TargetedSpectraExtractor::from_param(&params).unwrap();
    let picked = tse.pick_spectrum(&spectrum).unwrap();
    // Filtering also on FWHM, even fewer peaks get picked.
    assert_eq!(picked.len(), 2);
    for (peak, (mz, intensity)) in picked
        .peaks
        .iter()
        .zip([(85.014, 60754.7), (112.033, 21941.9)])
    {
        real_similar(peak.mz, mz);
        real_similar(f64::from(peak.intensity), intensity);
    }

    let mut unordered = MSSpectrum::default();
    unordered.peaks.push(Peak1D::new(10.0, 100.0));
    unordered.peaks.push(Peak1D::new(9.0, 100.0));
    let error = tse.pick_spectrum(&unordered).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("Spectrum must be sorted by position")
    );
}

fn feature(rt: f64, mz: f64, intensity: f32) -> Feature {
    Feature::from(BaseFeature::new(rt, mz, intensity))
}

#[test]
fn merge_features() {
    let mut ids = UniqueIdGenerator::from_seed(1);
    let mut features = FeatureMap::new();
    for (identifier, peptide_ref, intensity, mz, rt) in [
        (vec!["ident1"], "PeptideRef1", 1.0, 10.0, 100.0),
        (vec!["ident1", "ident2"], "PeptideRef1", 2.0, 20.0, 100.0),
        (vec!["ident3"], "PeptideRef3", 3.0, 30.0, 300.0),
    ] {
        let mut f = feature(rt, mz, intensity);
        f.unique_id = ids.get_unique_id();
        let list: Vec<String> = identifier.iter().map(|s| (*s).to_owned()).collect();
        f.metadata
            .insert("identifier".into(), MetaValue::from(list));
        f.metadata
            .insert("PeptideRef".into(), MetaValue::from(peptide_ref));
        features.features.push(f);
    }
    let tse = TargetedSpectraExtractor::new().unwrap();
    let mut merged = FeatureMap::new();
    tse.merge_features(&features, &mut merged, &mut ids)
        .unwrap();

    assert_eq!(merged.features.len(), 2);
    let merged_f1 = &merged.features[0];
    assert_eq!(meta_str(merged_f1, "PeptideRef"), "PeptideRef1");
    real_similar(merged_f1.mz, 16.6667);
    real_similar(merged_f1.rt, 100.0);
    assert_eq!(merged_f1.subordinates.len(), 2);
    assert_eq!(
        merged_f1.subordinates[0].metadata["identifier"]
            .as_string_list()
            .unwrap(),
        ["ident1"]
    );
    real_similar(merged_f1.subordinates[0].mz, 10.0);
    real_similar(merged_f1.subordinates[0].rt, 100.0);
    assert_eq!(
        merged_f1.subordinates[1].metadata["identifier"]
            .as_string_list()
            .unwrap(),
        ["ident1", "ident2"]
    );
    real_similar(merged_f1.subordinates[1].mz, 20.0);
    real_similar(merged_f1.subordinates[1].rt, 100.0);

    let merged_f2 = &merged.features[1];
    assert_eq!(meta_str(merged_f2, "PeptideRef"), "PeptideRef3");
    real_similar(merged_f2.mz, 30.0);
    real_similar(merged_f2.rt, 300.0);
    assert_eq!(merged_f2.subordinates.len(), 1);
    assert_eq!(
        merged_f2.subordinates[0].metadata["identifier"]
            .as_string_list()
            .unwrap(),
        ["ident3"]
    );
    real_similar(merged_f2.subordinates[0].mz, 30.0);
    real_similar(merged_f2.subordinates[0].rt, 300.0);
}

#[test]
fn annotate_spectra_with_ms1_features() {
    let mut ids = UniqueIdGenerator::from_seed(7);
    let mut f1 = feature(100.0, 10.0, 1.0);
    f1.metadata.insert(
        "identifier".into(),
        MetaValue::from(vec!["ident1".to_owned()]),
    );
    let mut sub1 = feature(110.0, 9.0, 2.0);
    sub1.metadata
        .insert("PeptideRef".into(), MetaValue::from("ident1"));
    // This one will not be used in the annotation.
    let mut sub2 = feature(210.0, 29.0, 2.0);
    sub2.metadata
        .insert("PeptideRef".into(), MetaValue::from("ident1"));
    f1.subordinates = vec![sub1, sub2];
    let mut ms1 = FeatureMap::new();
    ms1.features.push(f1);

    let spectrum = MSSpectrum {
        ms_level: 2,
        rt: 100.0,
        ..MSSpectrum::default()
    };

    let tse = TargetedSpectraExtractor::new().unwrap();
    let (annotated, ms2) = tse
        .annotate_spectra_with_ms1_features(&[spectrum], &ms1, &mut ids)
        .unwrap();
    assert_eq!(ms2.features.len(), 1);
    real_similar(ms2.features[0].mz, 0.0);
    real_similar(ms2.features[0].rt, 100.0);
    assert_eq!(meta_str(&ms2.features[0], "PeptideRef"), "ident1");
    assert!(ms2.features[0].subordinates.is_empty());
    assert_eq!(annotated.len(), 1);
    assert_eq!(annotated[0].name, "ident1");
    real_similar(annotated[0].rt, 100.0);
}

#[test]
fn construct_transitions_list() {
    let mut ms1_f1 = feature(100.0, 10.0, 1.0);
    ms1_f1
        .metadata
        .insert("PeptideRef".into(), MetaValue::from("ident1"));
    let mut ms1 = FeatureMap::new();
    ms1.features.push(ms1_f1);

    let mut ms2_f1 = feature(100.0, 10.0, 1.0);
    let mut sub1 = feature(110.0, 9.0, 2.0);
    sub1.metadata
        .insert("PeptideRef".into(), MetaValue::from("ident1"));
    sub1.metadata
        .insert("native_id".into(), MetaValue::from("ms2_f1_sub1"));
    let mut sub2 = feature(210.0, 29.0, 2.0);
    sub2.metadata
        .insert("PeptideRef".into(), MetaValue::from("ident1"));
    sub2.metadata
        .insert("native_id".into(), MetaValue::from("ms2_f1_sub2"));
    ms2_f1.subordinates = vec![sub1, sub2];
    let mut ms2 = FeatureMap::new();
    ms2.features.push(ms2_f1);

    let tse = TargetedSpectraExtractor::new().unwrap();
    let mut t_exp = openms::analysis::openswath::targeted_experiment::TargetedExperiment::new();
    tse.construct_transitions_list(&ms1, &ms2, &mut t_exp)
        .unwrap();
    assert_eq!(t_exp.transitions().len(), 1);
    let transition = &t_exp.transitions()[0];
    assert_eq!(transition.peptide_ref(), "ident1");
    assert_eq!(
        transition.cv_terms.metadata["PeptideRef"].as_str().unwrap(),
        "ident1"
    );
    assert_eq!(
        transition.cv_terms.metadata["native_id"].as_str().unwrap(),
        "ms2_f1_sub1"
    );
    // Native: the transition's native id, as the source streams it.
    assert_eq!(transition.native_id(), "ident1_ms2_f1_sub1_110");
}

#[test]
fn store_spectra_msp_preparation_removes_peaks() {
    // storeSpectraMSP's class test: the in-place peak removal before writing.
    // (The MSP writer itself is not ported; see the support document.)
    let mut spectrum = MSSpectrum {
        ms_level: 2,
        rt: 100.0,
        name: "spectr1".into(),
        ..MSSpectrum::default()
    };
    spectrum.precursors = vec![Precursor {
        mz: 100.0,
        ..Precursor::default()
    }];
    spectrum.peaks.push(Peak1D::new(100.0, 10.0));
    // The source pushes the second peak before setting it: (0, 0).
    spectrum.peaks.push(Peak1D::new(0.0, 0.0));
    let mut experiment = MSExperiment::default();
    experiment.spectra.push(spectrum);

    let tse = TargetedSpectraExtractor::new().unwrap();
    assert!(!tse.deisotoping.use_deisotoper);
    tse.remove_ms2_spectra_peaks(&mut experiment).unwrap();
    assert_eq!(experiment.spectra.len(), 1);
    assert_eq!(experiment.spectra[0].len(), 1);
    assert_eq!(experiment.spectra[0].peaks[0].intensity, 10.0);
    assert_eq!(experiment.spectra[0].peaks[0].mz, 100.0);

    // Native: an MS2 spectrum without a precursor is refused, unchanged.
    let mut no_precursor = MSExperiment::default();
    let ms2 = MSSpectrum {
        ms_level: 2,
        peaks: vec![Peak1D::new(5.0, 5.0)],
        ..MSSpectrum::default()
    };
    no_precursor.spectra.push(ms2);
    let before = no_precursor.clone();
    assert!(tse.remove_ms2_spectra_peaks(&mut no_precursor).is_err());
    assert_eq!(no_precursor, before);
}

/// A comparator with fixed scores, to check the ordering and bookkeeping of
/// the matching functions (independently derived; the class test's matching
/// sections need the MoNA MSP library, whose reader is not ported).
struct FixedScores {
    library: Vec<MSSpectrum>,
    scores: Vec<f64>,
}

impl Comparator for FixedScores {
    fn generate_scores(
        &self,
        _spectrum: &MSSpectrum,
        min_score: f64,
    ) -> openms::Result<Vec<(usize, f64)>> {
        Ok(self
            .scores
            .iter()
            .copied()
            .enumerate()
            .filter(|(_, s)| *s >= min_score)
            .collect())
    }
    fn init(
        &mut self,
        library: &[MSSpectrum],
        _options: &BTreeMap<String, MetaValue>,
    ) -> openms::Result<()> {
        self.library = library.to_vec();
        Ok(())
    }
    fn library(&self) -> &[MSSpectrum] {
        &self.library
    }
}

#[test]
fn matching_orders_filters_and_annotates() {
    let named = |name: &str| MSSpectrum {
        name: name.into(),
        ..MSSpectrum::default()
    };
    let mut library = vec![named("a"), named("b"), named("c"), named("d")];
    library[2]
        .metadata
        .insert("Comments".into(), MetaValue::from("\"accession=X\""));
    let mut cmp = FixedScores {
        library: Vec::new(),
        scores: vec![0.85, 0.5, 0.95, 0.9],
    };
    cmp.init(&library, &BTreeMap::new()).unwrap();

    let mut params = TargetedSpectraExtractor::new().unwrap().to_param().unwrap();
    set(&mut params, "top_matches_to_report", ParamValue::Integer(2));
    let tse = TargetedSpectraExtractor::from_param(&params).unwrap();
    let matches = tse.match_spectrum(&MSSpectrum::default(), &cmp).unwrap();
    // min_match_score 0.8 drops "b"; the best two, descending.
    assert_eq!(matches.len(), 2);
    assert_eq!(matches[0].spectrum.name, "c");
    assert_eq!(matches[0].score, 0.95);
    assert_eq!(matches[1].spectrum.name, "d");

    let mut features = FeatureMap::new();
    features.features.push(feature(1.0, 2.0, 0.0));
    let unmatched = tse
        .targeted_matching(&[MSSpectrum::default()], &cmp, &mut features)
        .unwrap();
    assert!(unmatched.is_empty());
    assert_eq!(
        meta_str(&features.features[0], "spectral_library_name"),
        "c"
    );
    assert_eq!(
        meta_f64(&features.features[0], "spectral_library_score"),
        0.95
    );
    assert_eq!(
        meta_str(&features.features[0], "spectral_library_comments"),
        "\"accession=X\""
    );

    let none = FixedScores {
        library: library.clone(),
        scores: vec![0.1, 0.2, 0.3, 0.4],
    };
    let unmatched = tse
        .targeted_matching(&[MSSpectrum::default()], &none, &mut features)
        .unwrap();
    assert_eq!(unmatched, [0]);
    assert_eq!(meta_str(&features.features[0], "spectral_library_name"), "");
    assert_eq!(
        meta_f64(&features.features[0], "spectral_library_score"),
        0.0
    );
    assert!(tse.targeted_matching(&[], &none, &mut features).is_err());
}

#[test]
fn binned_spectrum_comparator_scores_identity_as_one() {
    // Independently derived: the contrast angle of a spectrum with itself is 1.
    let mut spectrum = MSSpectrum::default();
    for (mz, intensity) in [(100.0, 10.0), (150.5, 20.0), (201.2, 5.0)] {
        spectrum.peaks.push(Peak1D::new(mz, intensity));
    }
    let mut cmp = BinnedSpectrumComparator::new().unwrap();
    let options = BTreeMap::from([
        ("bin_size".to_owned(), MetaValue::try_from(1.0).unwrap()),
        ("peak_spread".to_owned(), MetaValue::from(0_i64)),
        ("bin_offset".to_owned(), MetaValue::try_from(0.4).unwrap()),
    ]);
    cmp.init(std::slice::from_ref(&spectrum), &options).unwrap();
    assert_eq!(cmp.bin_size(), 1.0);
    assert_eq!(cmp.bin_offset(), 0.4);
    let scores = cmp.generate_scores(&spectrum, 0.5).unwrap();
    assert_eq!(scores.len(), 1);
    real_similar(scores[0].1, 1.0);
    let bad = BTreeMap::from([("bin_size".to_owned(), MetaValue::from("wide"))]);
    assert!(cmp.init(&[], &bad).is_err());
    assert_eq!(cmp.library().len(), 1);
}

#[cfg(feature = "mzml")]
mod class_test_data {
    use super::*;
    use openms::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition;
    use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
    use openms::analysis::openswath::targeted_experiment_helper::{Peptide, RTUnit, RetentionTime};

    fn experiment() -> MSExperiment {
        openms::format::mzml::read(
            include_bytes!(
                "data/targeted_spectra_extractor/TargetedSpectraExtractor_13C1_spectra0to100.mzML"
            )
            .as_slice(),
        )
        .unwrap()
    }

    /// The target list the class test loads with `TransitionTSVFile`
    /// (`retentionTimeInterpretation = minutes`), built here for exactly the
    /// fields `TargetedSpectraExtractor` reads: one transition per row, named
    /// `transition_name`, referencing the peptide `transition_group_id` with
    /// the row's precursor and product m/z; one peptide per group, created at
    /// its first row, with the row's `RetentionTime` in minutes. The reader is
    /// not ported; the annotation counts and names below, which the class test
    /// asserts, confirm the construction.
    fn targeted_experiment() -> TargetedExperiment {
        let text = include_str!(
            "data/targeted_spectra_extractor/TargetedSpectraExtractor_13CFlux_TraML.csv"
        );
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap().split(',').collect();
        let column = |name: &str| header.iter().position(|h| *h == name).unwrap();
        let (group, name, rt, precursor, product) = (
            column("transition_group_id"),
            column("transition_name"),
            column("RetentionTime"),
            column("PrecursorMz"),
            column("ProductMz"),
        );
        let mut exp = TargetedExperiment::new();
        for line in lines.filter(|l| !l.trim().is_empty()) {
            let fields: Vec<&str> = line.split(',').collect();
            let mut transition = ReactionMonitoringTransition::new();
            transition.set_name(fields[name]);
            transition.set_peptide_ref(fields[group]);
            transition.set_precursor_mz(fields[precursor].parse().unwrap());
            transition.set_product_mz(fields[product].parse().unwrap());
            exp.add_transition(transition);
            if !exp.has_peptide(fields[group]) {
                let mut peptide = Peptide::default();
                peptide.id = fields[group].to_owned();
                let mut retention_time = RetentionTime::default();
                retention_time.set_rt(fields[rt].parse().unwrap());
                retention_time.retention_time_unit = RTUnit::Minute;
                peptide.rts.push(retention_time);
                exp.add_peptide(peptide);
            }
        }
        exp
    }

    #[test]
    fn annotate_spectra() {
        let tse = class_test_extractor(&[]);
        let experiment = experiment();
        assert_eq!(experiment.spectra.len(), 100);
        let targets = targeted_experiment();
        let (annotated, features) = tse.annotate_spectra(&experiment.spectra, &targets).unwrap();
        assert_eq!(annotated.len(), 30);
        assert_eq!(annotated.len(), features.features.len());
        for (index, name, size) in [
            (0, "met-L.met-L_m0-0", 121),
            (4, "met-L.met-L_m1-0", 135),
            (8, "asp-L.asp-L_m0-0", 55),
            (12, "asp-L.asp-L_m1-0", 389),
            (16, "asp-L.asp-L_m2-1", 143),
            (20, "glu-L.glu-L_m5-5", 82),
            (24, "glu-L.glu-L_m2-2", 94),
            (29, "skm.skm_m4-4", 552),
        ] {
            assert_eq!(annotated[index].name, name);
            assert_eq!(annotated[index].len(), size);
        }
        for (index, name, rt, mz) in [
            (0, "met-L.met-L_m0-0", 80.22100000002, 148.052001953125),
            (4, "met-L.met-L_m1-0", 87.927, 149.054992675781),
            (8, "asp-L.asp-L_m0-0", 126.37699999998, 132.029998779297),
            (12, "asp-L.asp-L_m1-0", 131.73100000002, 133.033004760742),
            (16, "asp-L.asp-L_m2-1", 138.29599999998, 134.035995483398),
            (20, "glu-L.glu-L_m5-5", 141.70399999998, 151.061996459961),
            (24, "glu-L.glu-L_m2-2", 148.473, 148.052001953125),
            (29, "skm.skm_m4-4", 166.95400000002, 177.057998657227),
        ] {
            assert_eq!(meta_str(&features.features[index], "transition_name"), name);
            real_similar(features.features[index].rt, rt);
            real_similar(features.features[index].mz, mz);
        }

        let only = tse
            .annotate_spectra_only(&experiment.spectra, &targets)
            .unwrap();
        assert_eq!(only.len(), 30);
        for (index, name, size) in [
            (0, "met-L.met-L_m0-0", 121),
            (4, "met-L.met-L_m1-0", 135),
            (20, "glu-L.glu-L_m5-5", 82),
            (24, "glu-L.glu-L_m2-2", 94),
            (29, "skm.skm_m4-4", 552),
        ] {
            assert_eq!(only[index].name, name);
            assert_eq!(only[index].len(), size);
        }
    }

    /// Annotate, pick and drop the spectra whose picked spectrum is empty, as
    /// the class test's scoring and selection sections do.
    fn annotated_and_picked(
        tse: &TargetedSpectraExtractor,
    ) -> (Vec<MSSpectrum>, Vec<MSSpectrum>, FeatureMap) {
        let experiment = experiment();
        let (annotated, features) = tse
            .annotate_spectra(&experiment.spectra, &targeted_experiment())
            .unwrap();
        let mut kept = (Vec::new(), Vec::new(), FeatureMap::new());
        for (i, spectrum) in annotated.iter().enumerate() {
            let picked = tse.pick_spectrum(spectrum).unwrap();
            if !picked.is_empty() {
                kept.0.push(spectrum.clone());
                kept.1.push(picked);
                kept.2.features.push(features.features[i].clone());
            }
        }
        kept
    }

    const SCORES: [(usize, &str, f64, f64, f64, f64); 6] = [
        (
            0,
            "met-L.met-L_m0-0",
            13.095,
            5.3508939743042,
            3.96267318725586,
            3.7814,
        ),
        (4, "asp-L.asp-L_m1-0", 10.6598, 6.49946, 2.65215, 1.50815),
        (8, "asp-L.asp-L_m2-1", 16.3603, 5.52142, 3.44492, 7.39392),
        (11, "asp-L.asp-L_m2-2", 14.6208, 5.48532, 3.79353, 5.34195),
        (15, "glu-L.glu-L_m1-1", 12.3466, 5.49089, 3.53584, 3.3199),
        (19, "skm.skm_m4-4", 10.1824, 6.60354, 2.02877, 1.55012),
    ];

    #[test]
    fn score_spectra() {
        let tse = class_test_extractor(&[]);
        let (annotated, picked, mut features) = annotated_and_picked(&tse);
        assert_eq!(annotated.len(), 20);
        assert_eq!(picked.len(), features.features.len());
        let scored = tse
            .score_spectra(&annotated, &picked, &mut features)
            .unwrap();
        assert_eq!(scored.len(), 20);
        for (index, name, score, tic, fwhm, snr) in SCORES {
            assert_eq!(scored[index].name, name);
            let arrays = &scored[index].float_data_arrays;
            real_similar(f64::from(arrays[1].data[0]), score);
            real_similar(f64::from(arrays[2].data[0]), tic);
            real_similar(f64::from(arrays[3].data[0]), fwhm);
            real_similar(f64::from(arrays[4].data[0]), snr);
        }
        for (index, name, score, tic, fwhm, snr, avg_fwhm) in [
            (
                0,
                "met-L.met-L_m0-0",
                13.095,
                5.3508939743042,
                3.96267318725586,
                3.7814010941664,
                0.252354895075162,
            ),
            (
                4,
                "asp-L.asp-L_m1-0",
                10.6598,
                6.49945796336373,
                2.65214624318674,
                1.50815156166727,
                0.377053114084097,
            ),
            (
                8,
                "asp-L.asp-L_m2-1",
                16.3603,
                5.52141560620828,
                3.44491858720322,
                7.39392424854829,
                0.290282621979713,
            ),
            (
                11,
                "asp-L.asp-L_m2-2",
                14.6208,
                5.48531541983726,
                3.79353289768299,
                5.34195250209406,
                0.263606518507004,
            ),
            (
                15,
                "glu-L.glu-L_m1-1",
                12.3466,
                5.49089446225569,
                3.53583924309525,
                3.31989908680269,
                0.282818287611008,
            ),
            (
                19,
                "skm.skm_m4-4",
                10.1824,
                6.60354130105922,
                2.02876571323501,
                1.55011855538165,
                0.492910538400923,
            ),
        ] {
            let feature = &features.features[index];
            assert_eq!(meta_str(feature, "transition_name"), name);
            real_similar(f64::from(feature.intensity), score);
            real_similar(meta_f64(feature, "log10_total_tic"), tic);
            real_similar(meta_f64(feature, "inverse_avgFWHM"), fwhm);
            real_similar(meta_f64(feature, "avgSNR"), snr);
            real_similar(meta_f64(feature, "avgFWHM"), avg_fwhm);
        }
        features.features.pop();
        assert!(
            tse.score_spectra(&annotated, &picked, &mut features)
                .is_err()
        );

        let only = tse.score_spectra_only(&annotated, &picked).unwrap();
        assert_eq!(only.len(), 20);
        for (index, name, score, tic, fwhm, snr) in SCORES {
            assert_eq!(only[index].name, name);
            let arrays = &only[index].float_data_arrays;
            real_similar(f64::from(arrays[1].data[0]), score);
            real_similar(f64::from(arrays[2].data[0]), tic);
            real_similar(f64::from(arrays[3].data[0]), fwhm);
            real_similar(f64::from(arrays[4].data[0]), snr);
        }
    }

    fn find<'a>(spectra: &'a [MSSpectrum], name: &str) -> &'a MSSpectrum {
        spectra.iter().find(|s| s.name == name).unwrap()
    }

    #[test]
    fn select_spectra() {
        let min_select_score = 15.0;
        let tse =
            class_test_extractor(&[("min_select_score", ParamValue::Float(min_select_score))]);
        let (annotated, picked, mut features) = annotated_and_picked(&tse);
        let scored = tse
            .score_spectra(&annotated, &picked, &mut features)
            .unwrap();
        let (selected, selected_features) = tse.select_spectra(&scored, &features).unwrap();
        // Note: "met-L.met-L_m0-0" no longer reaches min_select_score with the
        // more accurate (interpolated) S/N estimate.
        assert_eq!(selected.len(), 2);
        assert_eq!(selected.len(), selected_features.features.len());
        for (spectrum, feature) in selected.iter().zip(&selected_features.features) {
            assert_ne!(spectrum.name, "");
            assert_eq!(spectrum.name, meta_str(feature, "transition_name"));
            assert_eq!(spectrum.float_data_arrays[1].data[0], feature.intensity);
            assert!(f64::from(spectrum.float_data_arrays[1].data[0]) >= min_select_score);
        }
        real_similar(
            f64::from(find(&selected, "asp-L.asp-L_m2-1").float_data_arrays[1].data[0]),
            16.3603,
        );
        real_similar(
            f64::from(find(&selected, "asp-L.asp-L_m2-2").float_data_arrays[1].data[0]),
            16.3603,
        );
        features.features.pop();
        assert!(tse.select_spectra(&scored, &features).is_err());

        let scored_only = tse.score_spectra_only(&annotated, &picked).unwrap();
        let selected = tse.select_spectra_only(&scored_only).unwrap();
        assert_eq!(selected.len(), 2);
        for spectrum in &selected {
            assert_ne!(spectrum.name, "");
            assert!(f64::from(spectrum.float_data_arrays[1].data[0]) >= min_select_score);
        }
        real_similar(
            f64::from(find(&selected, "asp-L.asp-L_m2-1").float_data_arrays[1].data[0]),
            16.3603,
        );
        real_similar(
            f64::from(find(&selected, "asp-L.asp-L_m2-2").float_data_arrays[1].data[0]),
            16.3603,
        );
    }

    #[test]
    fn extract_spectra() {
        let tse = class_test_extractor(&[("min_select_score", ParamValue::Float(15.0))]);
        let experiment = experiment();
        let targets = targeted_experiment();
        let (extracted, features) = tse.extract_spectra(&experiment, &targets).unwrap();
        assert_eq!(extracted.len(), features.features.len());
        // Note: "met-L.met-L_m0-0" no longer reaches min_select_score with the
        // more accurate (interpolated) S/N estimate.
        assert_eq!(extracted.len(), 2);
        real_similar(
            f64::from(find(&extracted, "asp-L.asp-L_m2-1").float_data_arrays[1].data[0]),
            16.3603,
        );
        real_similar(
            f64::from(find(&extracted, "asp-L.asp-L_m2-2").float_data_arrays[1].data[0]),
            16.3603,
        );

        let only = tse.extract_spectra_only(&experiment, &targets).unwrap();
        assert_eq!(only.len(), 2);
        real_similar(
            f64::from(find(&only, "asp-L.asp-L_m2-1").float_data_arrays[1].data[0]),
            16.3603,
        );
        real_similar(
            f64::from(find(&only, "asp-L.asp-L_m2-2").float_data_arrays[1].data[0]),
            16.3603,
        );
    }
}
