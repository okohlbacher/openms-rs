// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! QC metrics (`QC/QCBase.h`, `SpectrumCount.h`, `TIC.h`, `FeatureSummary.h`,
//! `MissedCleavages.h`, `IdentificationSummary.h`) and their flag sets
//! (`DATASTRUCTURES/FlagSet.h`).
//!
//! Tier 3: every section of the class tests `FlagSet_test.cpp`,
//! `QCBase_test.cpp`, `SpectrumCount_test.cpp`, `TIC_test.cpp`,
//! `FeatureSummary_test.cpp`, `MissedCleavages_test.cpp` and
//! `IdentificationSummary_test.cpp` (core bc9cc12,
//! `src/tests/class_tests/openms/source/`), with their literals. The metrics
//! are also exercised end to end, tier 1, by `tests/mzqc_file.rs`.
//! `TEST_REAL_SIMILAR` is ClassTest's relative 1e-5 / absolute 1e-5 check.

use openms::analysis::qc::qc_base::{
    NAMES_OF_REQUIRES, NAMES_OF_TOLERANCE_UNIT, has_pep_id_consensus_map, has_pep_id_feature_map,
    is_labeled_experiment,
};
use openms::analysis::qc::{
    FeatureSummary, FeatureSummaryResult, IdentificationSummary, IdentificationSummaryResult,
    MissedCleavages, QCBase, Requires, SpectraMap, SpectrumCount, Status, Tic, TicResult,
    ToleranceUnit, UniqueId,
};
use openms::chemistry::AASequence;
use openms::data_structures::{FlagEnum, FlagSet};
use openms::format::mztab::MzTabMetaData;
use openms::identification::{
    PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
};
use openms::kernel::{ConsensusMap, Feature, FeatureMap, MSExperiment, MSSpectrum, Peak1D};
use openms::metadata::{DataProcessing, MetaValue};

/// ClassTest `TEST_REAL_SIMILAR` with its default tolerances.
fn similar(actual: f64, expected: f64) {
    let ok = (actual - expected).abs() <= 1e-5
        || (expected != 0.0 && (actual / expected).max(expected / actual) <= 1.0 + 1e-5);
    assert!(ok, "{actual} is not similar to {expected}");
}

// ---------------------------------------------------------------- FlagSet

#[derive(Clone, Copy, Debug)]
enum EnumTest {
    Red,
    Green,
    Blue,
}
impl FlagEnum for EnumTest {
    fn bit(self) -> u32 {
        self as u32
    }
}
type Fst = FlagSet<EnumTest>;

#[derive(Clone, Copy, Debug)]
struct TooLarge;
impl FlagEnum for TooLarge {
    fn bit(self) -> u32 {
        64
    }
}

#[test]
fn flag_set_constructors_and_equality() {
    // FlagSet(), explicit FlagSet(const ENUM&), copy, assignment, ==.
    assert_ne!(Fst::from_flag(EnumTest::Red), Fst::new());
    assert_eq!(Fst::from_flag(EnumTest::Red).value(), 1);
    assert_eq!(Fst::from(EnumTest::Green).value(), 2);
    let gg = Fst::from_flag(EnumTest::Green);
    let target = gg;
    assert_eq!(target, gg);
    assert_ne!(gg, Fst::default());
    // The source's TEST_PRECONDITION_VIOLATED for a value of 64: here such a
    // value contributes no bit rather than aliasing another flag.
    assert!(FlagSet::from_flag(TooLarge).is_empty());
}

#[test]
fn flag_set_and() {
    let gg = Fst::from_flag(EnumTest::Green);
    let rr = Fst::from_flag(EnumTest::Red);
    assert_eq!(gg & EnumTest::Red, Fst::new());
    assert_eq!(gg & EnumTest::Green, gg);
    assert_eq!(gg & rr, Fst::new());
    assert_eq!(gg & gg, gg);
    let mut empty = gg;
    empty &= EnumTest::Red;
    assert_eq!(empty, Fst::new());
    let mut just_green = gg;
    just_green &= EnumTest::Green;
    assert_eq!(just_green, gg);
    let mut empty = gg;
    empty &= rr;
    assert_eq!(empty, Fst::new());
    let mut just_green = gg;
    just_green &= gg;
    assert_eq!(just_green, gg);
}

#[test]
fn flag_set_or_and_plus() {
    let gg = Fst::from_flag(EnumTest::Green);
    assert_eq!((gg | EnumTest::Red).value(), 3);
    assert_eq!(gg | EnumTest::Green, gg);
    let empty_or_red = Fst::new() | EnumTest::Red;
    assert_eq!(empty_or_red, Fst::from_flag(EnumTest::Red));
    assert_eq!((empty_or_red | Fst::from_flag(EnumTest::Green)).value(), 3);
    let mut green_or_red = gg;
    green_or_red |= EnumTest::Red;
    assert_eq!(green_or_red.value(), 3);
    let mut red_or_green = empty_or_red;
    red_or_green |= Fst::from_flag(EnumTest::Green);
    assert_eq!(red_or_green.value(), 3);

    assert_eq!((gg + EnumTest::Red).value(), 3);
    assert_eq!(gg + EnumTest::Green, gg);
    assert_eq!(Fst::new() + EnumTest::Red, Fst::from_flag(EnumTest::Red));
    assert_eq!((Fst::from_flag(EnumTest::Red) + gg).value(), 3);
    let mut sum = gg;
    sum += EnumTest::Red;
    assert_eq!(sum.value(), 3);
    let mut sum = Fst::new();
    sum += Fst::from_flag(EnumTest::Red);
    sum += Fst::from_flag(EnumTest::Green);
    assert_eq!(sum.value(), 3);
}

#[test]
fn flag_set_minus() {
    let red = Fst::from_flag(EnumTest::Red);
    assert_eq!(Fst::new() - red, Fst::new());
    let red_or_green = red + EnumTest::Green;
    assert_eq!(red_or_green - Fst::from_flag(EnumTest::Blue), red_or_green);
    assert_eq!(red_or_green - EnumTest::Green, red);
    assert_eq!(Fst::new() - EnumTest::Red, Fst::new());
    assert_eq!(red_or_green - EnumTest::Blue, red_or_green);
    let mut only_red = red_or_green;
    only_red -= Fst::from_flag(EnumTest::Green);
    assert_eq!(only_red, red);
    let mut only_red = red_or_green;
    only_red -= EnumTest::Green;
    assert_eq!(only_red, red);
    let mut unchanged = red_or_green;
    unchanged -= EnumTest::Blue;
    assert_eq!(unchanged, red_or_green);
}

#[test]
fn flag_set_super_sets_and_empty() {
    let empty = Fst::new() - EnumTest::Red;
    assert!(Fst::new().is_super_set_of(&empty));
    assert!(empty.is_super_set_of(&Fst::new()));
    let red_or_green = Fst::from_flag(EnumTest::Red) + EnumTest::Green;
    let all = red_or_green + EnumTest::Blue;
    assert!(all.is_super_set_of(&red_or_green));
    assert!(all.is_super_set_of(&all));
    assert!(all.is_super_set_of(&empty));
    assert!(!red_or_green.is_super_set_of(&all));
    assert!(!empty.is_super_set_of(&all));
    assert!(!Fst::from_flag(EnumTest::Red).is_super_set_of(&all));

    assert!(!Fst::new().is_super_set_of_flag(EnumTest::Red));
    assert!(all.is_super_set_of_flag(EnumTest::Red));
    assert!(all.is_super_set_of_flag(EnumTest::Blue));
    assert!(all.is_super_set_of_flag(EnumTest::Green));

    assert!(Fst::new().is_empty());
    assert!(!Fst::from_flag(EnumTest::Red).is_empty());
    assert!(!red_or_green.is_empty());
}

// ----------------------------------------------------------------- QCBase

fn spectrum(native_id: &str) -> MSSpectrum {
    MSSpectrum {
        native_id: native_id.into(),
        ..MSSpectrum::default()
    }
}

#[test]
fn spectra_map_indexes_native_ids() {
    assert!(SpectraMap::new().is_empty());
    let mut exp = MSExperiment::new();
    exp.spectra = vec![
        spectrum("XTandem::0"),
        spectrum("XTandem::1"),
        spectrum("XTandem::2"),
    ];
    let mut map = SpectraMap::new();
    map.calculate_map(&exp);
    assert_eq!(map.len(), 3);
    assert_eq!(map.at("XTandem::0").unwrap(), 0);
    assert_eq!(map.at("XTandem::1").unwrap(), 1);
    assert_eq!(map.at("XTandem::2").unwrap(), 2);
    assert!(map.at("XTandem::15").is_err());
    assert_eq!(SpectraMap::from_experiment(&exp).len(), 3);
    assert!(!map.is_empty());
    map.clear();
    assert!(map.is_empty());
}

#[test]
fn names_of_requires_has_one_entry_per_value() {
    assert_eq!(NAMES_OF_REQUIRES.len(), Requires::SIZE_OF_REQUIRES);
    assert_eq!(NAMES_OF_REQUIRES[Requires::Nothing as usize], "fail");
    assert_eq!(
        NAMES_OF_REQUIRES[Requires::TrafoAlign as usize],
        "trafoAlign.trafoXML"
    );
    assert_eq!(NAMES_OF_REQUIRES[Requires::Id as usize], "id.idXML");
    assert_eq!(Requires::Id.name(), "id.idXML");
    assert_eq!(NAMES_OF_TOLERANCE_UNIT, ["auto", "ppm", "da"]);
    assert_eq!(ToleranceUnit::Ppm.name(), "ppm");
}

#[test]
fn is_labeled_experiment_needs_exactly_isobaric_analyzer() {
    let mut map = ConsensusMap::default();
    assert!(!is_labeled_experiment(&map));
    map.data_processing.push(DataProcessing::default());
    for name in [
        "",
        "FeatureFinder",
        "isobaricanalyzer",
        "IsobaricAnalyzer extra",
    ] {
        map.data_processing[0].software.name = name.into();
        assert!(!is_labeled_experiment(&map), "{name}");
    }
    let mut processing = DataProcessing::default();
    processing.software.name = "IsobaricAnalyzer".into();
    map.data_processing.push(processing);
    map.data_processing.push(DataProcessing::default());
    assert!(is_labeled_experiment(&map));
    assert_eq!(map.data_processing.len(), 3);
}

/// Independently derived from `QCBase::isRunnable`: a status must hold every
/// required bit, and each missing one yields the source's warning line.
#[test]
fn is_runnable_checks_the_required_inputs() {
    let tic = Tic::new();
    let status = Status::new();
    assert!(!tic.is_runnable(&status));
    assert_eq!(
        tic.runnable_warnings(&status),
        ["Note: Metric 'TIC' cannot run because input data 'raw.mzML' is missing!"]
    );
    let status = Status::from_flag(Requires::RawMzml) | Requires::Id;
    assert!(tic.is_runnable(&status));
    assert!(tic.runnable_warnings(&status).is_empty());
}

#[test]
fn has_pep_id_looks_at_features_and_unassigned_identifications() {
    let mut fmap = FeatureMap::default();
    fmap.features.push(Feature::default());
    assert!(!has_pep_id_feature_map(&fmap));
    fmap.features[0]
        .peptide_identifications
        .push(PeptideIdentification::default());
    assert!(has_pep_id_feature_map(&fmap));
    let mut cmap = ConsensusMap::default();
    assert!(!has_pep_id_consensus_map(&cmap));
    cmap.unassigned_peptide_identifications
        .push(PeptideIdentification::default());
    assert!(has_pep_id_consensus_map(&cmap));
}

// ---------------------------------------------------------- SpectrumCount

fn level(ms_level: u32) -> MSSpectrum {
    MSSpectrum {
        ms_level,
        ..MSSpectrum::default()
    }
}

#[test]
fn spectrum_count_counts_per_level() {
    let sc = SpectrumCount::new();
    let mut exp = MSExperiment::new();
    exp.spectra.extend((0..3).map(|_| level(1)));
    exp.spectra.extend((0..5).map(|_| level(2)));
    exp.spectra.push(level(3));
    let counts = sc.compute(&exp);
    assert_eq!(counts.len(), 3);
    assert_eq!(counts[&1], 3);
    assert_eq!(counts[&2], 5);
    assert_eq!(counts[&3], 1);
    assert!(sc.compute(&MSExperiment::new()).is_empty());
    assert_eq!(sc.name(), "SpectrumCount");
    assert_eq!(sc.requirements(), Status::from_flag(Requires::RawMzml));
}

// -------------------------------------------------------------------- TIC

fn tic_experiment() -> MSExperiment {
    let mut exp = MSExperiment::new();
    for (rt, peaks) in [
        (0.5, [(100.0, 30.0), (150.0, 70.0)]),
        (1.0, [(200.0, 500.0), (250.0, 2000.0)]),
        (1.5, [(300.0, 5.0), (350.0, 15.0)]),
        (2.0, [(400.0, 1000.0), (450.0, 3000.0)]),
    ] {
        let mut spec = level(1);
        spec.rt = rt;
        spec.peaks = peaks
            .iter()
            .map(|&(mz, intensity)| Peak1D { mz, intensity })
            .collect();
        exp.spectra.push(spec);
    }
    exp
}

#[test]
fn tic_name_and_requirements() {
    let tic = Tic::new();
    assert_eq!(tic.name(), "TIC");
    assert_eq!(tic.requirements(), Status::from_flag(Requires::RawMzml));
}

#[test]
fn tic_computes_trace_area_and_stability_counters() {
    let tic = Tic::new();
    let result = tic.compute(&tic_experiment(), 0.0, 1).unwrap();
    assert_eq!(result.intensities, [100, 2500, 20, 4000]);
    assert_eq!(result.retention_times.len(), 4);
    assert_eq!(result.relative_intensities.len(), 4);
    for (actual, expected) in result
        .relative_intensities
        .iter()
        .zip([2.5, 62.5, 0.5, 100.0])
    {
        similar(f64::from(*actual), expected);
    }
    for (actual, expected) in result.retention_times.iter().zip([0.5, 1.0, 1.5, 2.0]) {
        similar(f64::from(*actual), expected);
    }
    assert_eq!(result.area, 6620);
    assert_eq!(result.jump, 2);
    assert_eq!(result.fall, 1);
    // An empty experiment still yields an empty result.
    assert_eq!(
        tic.compute(&MSExperiment::new(), 0.0, 1).unwrap(),
        TicResult::default()
    );
}

#[test]
fn tic_result_equality_includes_relative_intensities() {
    let mut a = TicResult::default();
    let mut b = TicResult::default();
    a.relative_intensities = vec![1.0];
    b.relative_intensities = vec![2.0];
    assert_ne!(a, b);
    b.relative_intensities = vec![1.0];
    assert_eq!(a, b);
}

/// Independently derived from `TIC::addMetaDataMetricsToMzTab`: a parameter
/// per non-empty result, named `TIC_<1-based index>`, keyed by the map's size,
/// with RTs as `StringUtils::toStr(float, false)` writes them (three
/// fractional digits, trailing zeros trimmed to one) and integer intensities.
#[test]
fn tic_adds_its_traces_to_mztab_metadata() {
    let tic = Tic::new();
    let result = tic.compute(&tic_experiment(), 0.0, 1).unwrap();
    let mut meta = MzTabMetaData::default();
    tic.add_meta_data_metrics_to_mz_tab(&mut meta, &[TicResult::default(), result]);
    assert_eq!(meta.custom.len(), 1);
    let parameter = &meta.custom[&0];
    assert_eq!(parameter.cv_label(), "total ion current");
    assert_eq!(parameter.accession(), "MS:1000285");
    assert_eq!(parameter.name(), "TIC_2");
    assert_eq!(
        parameter.value(),
        "[0.5, 100, 1.0, 2500, 1.5, 20, 2.0, 4000]"
    );
}

/// Independently derived: binning a trace whose points already sit on the
/// grid keeps them, and a non-positive bin size disables binning.
#[test]
fn tic_bins_only_for_a_positive_bin_size() {
    let tic = Tic::new();
    let exact = tic.compute(&tic_experiment(), 0.0, 1).unwrap();
    assert_eq!(tic.compute(&tic_experiment(), -1.0, 1).unwrap(), exact);
    let binned = tic.compute(&tic_experiment(), 0.5, 1).unwrap();
    assert_eq!(binned.intensities, exact.intensities);
    assert_eq!(binned.retention_times, exact.retention_times);
}

// --------------------------------------------------------- FeatureSummary

fn feature_with_deviation(deviation: Option<f64>) -> Feature {
    let mut feature = Feature::default();
    if let Some(value) = deviation {
        feature
            .metadata
            .insert("rt_deviation".into(), MetaValue::try_from(value).unwrap());
    }
    feature
}

#[test]
fn feature_summary_counts_features_and_averages_rt_deviation() {
    let fs = FeatureSummary::new();
    let mut fm = FeatureMap::default();
    fm.features = vec![
        feature_with_deviation(Some(10.0)),
        feature_with_deviation(Some(20.0)),
        feature_with_deviation(None),
    ];
    let r = fs.compute(&fm).unwrap();
    assert_eq!(r.feature_count, 3);
    similar(f64::from(r.rt_shift_mean), 15.0);

    let r0 = fs.compute(&FeatureMap::default()).unwrap();
    assert_eq!(r0.feature_count, 0);
    similar(f64::from(r0.rt_shift_mean), 0.0);

    let mut no_dev = FeatureMap::default();
    no_dev.features = vec![Feature::default(), Feature::default()];
    let rn = fs.compute(&no_dev).unwrap();
    assert_eq!(rn.feature_count, 2);
    similar(f64::from(rn.rt_shift_mean), 0.0);

    assert_eq!(fs.name(), "Summary of features from featureXML file");
    assert_eq!(fs.requirements(), Status::from_flag(Requires::PreFdrFeat));
}

#[test]
fn feature_summary_result_equality() {
    let mut a = FeatureSummaryResult::default();
    let mut b = FeatureSummaryResult::default();
    assert_eq!(a, b);
    a.feature_count = 5;
    assert_ne!(a, b);
    b.feature_count = 5;
    assert_eq!(a, b);
    a.rt_shift_mean = 3.5;
    assert_ne!(a, b);
    b.rt_shift_mean = 3.5;
    assert_eq!(a, b);
}

/// Native: a non-numeric `rt_deviation` is an error where the source's
/// `DataValue` conversion throws.
#[test]
fn feature_summary_refuses_a_non_numeric_deviation() {
    let mut fm = FeatureMap::default();
    let mut feature = Feature::default();
    feature
        .metadata
        .insert("rt_deviation".into(), MetaValue::from("ten"));
    fm.features.push(feature);
    assert!(FeatureSummary::new().compute(&fm).is_err());
}

// -------------------------------------------------------- MissedCleavages

fn pep_id(sequences: &[&str]) -> PeptideIdentification {
    PeptideIdentification {
        hits: sequences
            .iter()
            .map(|s| PeptideHit {
                sequence: AASequence::parse(s).unwrap(),
                ..PeptideHit::default()
            })
            .collect(),
        ..PeptideIdentification::default()
    }
}

fn protein_ids(enzyme: &str, missed_cleavages: u32) -> Vec<ProteinIdentification> {
    let mut prot = ProteinIdentification::default();
    prot.search_parameters.digestion_enzyme = enzyme.into();
    prot.search_parameters.missed_cleavages = missed_cleavages;
    vec![prot]
}

fn feature_with(pep_ids: Vec<PeptideIdentification>) -> Feature {
    let mut feature = Feature::default();
    feature.peptide_identifications = pep_ids;
    feature
        .metadata
        .insert("FWHM".into(), MetaValue::try_from(32.21).unwrap());
    feature
}

fn annotated(pep_id: &PeptideIdentification) -> Option<i64> {
    pep_id
        .hits
        .first()
        .and_then(|hit| hit.metadata.get("missed_cleavages"))
        .map(|value| value.as_i64().unwrap())
}

#[test]
fn missed_cleavages_on_a_feature_map() {
    let pep_ids = vec![
        pep_id(&["AAAAAAAAAAAAAK"]),
        pep_id(&["AAAAKAAAAAR"]),
        pep_id(&[]),
    ];
    let mut feature_map = FeatureMap::default();
    feature_map.features = vec![
        feature_with(Vec::new()),
        feature_with(pep_ids.clone()),
        feature_with(Vec::new()),
    ];
    feature_map.unassigned_peptide_identifications =
        vec![pep_id(&["AAAAKAAAAAR"]), pep_id(&["AAAAKAAAAAR"])];
    // ProteaseDB::getEnzyme("trypsin"), whose name is "Trypsin".
    feature_map.protein_identifications = protein_ids("Trypsin", 2);

    let mut mc = MissedCleavages::new();
    mc.compute_feature_map(&mut feature_map).unwrap();
    let result = &mc.results()[0];
    assert_eq!(result.len(), 2);
    assert_eq!(result[&0], 1);
    assert_eq!(result[&1], 3);

    // The first hits carry the annotation, visited as
    // applyFunctionOnPeptideIDs visits them: features, then unassigned.
    let frequ: Vec<i64> = feature_map
        .features
        .iter()
        .flat_map(|f| f.peptide_identifications.iter())
        .chain(&feature_map.unassigned_peptide_identifications)
        .filter_map(annotated)
        .collect();
    assert_eq!(frequ, [0, 1, 1, 1]);

    // An empty feature map has an empty result.
    let mut mc_empty = MissedCleavages::new();
    mc_empty
        .compute_feature_map(&mut FeatureMap::default())
        .unwrap();
    assert!(mc_empty.results()[0].is_empty());

    // No protein identifications.
    let mut no_protein = FeatureMap::default();
    no_protein.features.push(feature_with(pep_ids.clone()));
    let error = MissedCleavages::new()
        .compute_feature_map(&mut no_protein)
        .unwrap_err();
    assert_eq!(
        error.to_string(),
        "missing information: Missing information in ProteinIdentifications."
    );

    // No enzyme: a default ProteinIdentification names "unknown_enzyme".
    let mut no_enzyme = FeatureMap::default();
    no_enzyme.features.push(feature_with(pep_ids));
    no_enzyme.protein_identifications = vec![ProteinIdentification::default()];
    let error = MissedCleavages::new()
        .compute_feature_map(&mut no_enzyme)
        .unwrap_err();
    assert_eq!(
        error.to_string(),
        "missing information: No digestion enzyme in FeatureMap detected. No computation possible."
    );

    // More missed cleavages than the search allowed: counted, and warned.
    let mut feature_map_3 = FeatureMap::default();
    feature_map_3
        .features
        .push(feature_with(vec![pep_id(&["AAAKAAARAAAKAAAR"])]));
    feature_map_3.protein_identifications = protein_ids("trypsin", 2);
    let mut mc_3 = MissedCleavages::new();
    mc_3.compute_feature_map(&mut feature_map_3).unwrap();
    let result_3 = &mc_3.results()[0];
    assert_eq!(result_3.len(), 1);
    assert_eq!(result_3[&3], 1);
    assert_eq!(mc_3.log.len(), 1);
    assert!(mc_3.log[0].starts_with("Observed number of missed cleavages: 3 is greater than: 2"));

    let mc = MissedCleavages::new();
    assert_eq!(mc.name(), "MissedCleavages");
    assert_eq!(mc.requirements(), Status::from_flag(Requires::PostFdrFeat));
}

// -------------------------------------------------- IdentificationSummary

fn scored_pep(sequence: &str, score_type: &str, significance: f32) -> PeptideIdentification {
    PeptideIdentification {
        score_type: score_type.into(),
        significance_threshold: f64::from(significance),
        ..pep_id(&[sequence])
    }
}

fn scored_proteins(score_type: &str, significance: f32) -> Vec<ProteinIdentification> {
    let mut prot_ids = protein_ids("Trypsin", 4);
    prot_ids[0].score_type = score_type.into();
    prot_ids[0].significance_threshold = f64::from(significance);
    prot_ids[0].hits = vec![
        ProteinHit {
            accession: "PROT1".into(),
            score: 0.9,
            ..ProteinHit::default()
        },
        ProteinHit {
            accession: "PROT2".into(),
            score: 0.5,
            ..ProteinHit::default()
        },
    ];
    prot_ids
}

#[test]
fn identification_summary_with_fdr_scores() {
    let prot_ids = scored_proteins("FDR", 0.01);
    let mut pep_ids = vec![
        scored_pep("AAAAAAAAAAAAAK", "FDR", 0.05),
        scored_pep("AAAAKAAAAAR", "FDR", 0.05),
        scored_pep("AAAKAAARAAAKAAAR", "FDR", 0.05),
    ];
    let r = IdentificationSummary::new()
        .compute(&prot_ids, &mut pep_ids)
        .unwrap();
    assert_eq!(r.peptide_spectrum_matches, 3);
    assert_eq!(r.unique_peptides.count, 3);
    assert_eq!(r.unique_proteins.count, 2);
    similar(f64::from(r.unique_peptides.fdr_threshold), 0.05);
    similar(f64::from(r.unique_proteins.fdr_threshold), 0.01);
    similar(f64::from(r.missed_cleavages_mean), 4.0 / 3.0);
    similar(r.protein_hit_scores_mean, 0.7);
    similar(r.peptide_length_mean, 41.0 / 3.0);
}

#[test]
fn identification_summary_without_fdr_scores() {
    let prot_ids = scored_proteins("Posterior Error Probability", 0.0);
    let mut pep_ids = vec![
        scored_pep("AAAAAAAAAAAAAK", "Posterior Error Probability", 0.0),
        scored_pep("AAAAKAAAAAR", "Posterior Error Probability", 0.0),
    ];
    let r = IdentificationSummary::new()
        .compute(&prot_ids, &mut pep_ids)
        .unwrap();
    assert_eq!(r.peptide_spectrum_matches, 2);
    assert_eq!(r.unique_peptides.count, 2);
    similar(f64::from(r.unique_peptides.fdr_threshold), -1.0);
    similar(f64::from(r.unique_proteins.fdr_threshold), -1.0);
}

#[test]
fn identification_summary_defaults_equality_and_identity() {
    let u = UniqueId::default();
    assert_eq!(u.count, 0);
    similar(f64::from(u.fdr_threshold), -1.0);

    let mut a = IdentificationSummaryResult::default();
    let mut b = IdentificationSummaryResult::default();
    assert_eq!(a, b);
    a.peptide_spectrum_matches = 7;
    assert_ne!(a, b);
    b.peptide_spectrum_matches = 7;
    assert_eq!(a, b);
    a.unique_peptides.count = 3;
    assert_ne!(a, b);
    b.unique_peptides.count = 3;
    assert_eq!(a, b);
    a.protein_hit_scores_mean = 0.7;
    assert_ne!(a, b);
    b.protein_hit_scores_mean = 0.7;
    assert_eq!(a, b);

    let is = IdentificationSummary::new();
    assert_eq!(
        is.name(),
        "Summary of detected Proteins and Peptides from idXML file"
    );
    assert_eq!(is.requirements(), Status::from_flag(Requires::Id));
}

/// Independently derived from `PeptideIdentification::empty`: an
/// identification without hits still counts as a PSM once it has a score
/// type, but contributes no peptide.
#[test]
fn identification_summary_counts_psms_by_the_source_emptiness() {
    let prot_ids = scored_proteins("FDR", 0.01);
    let mut pep_ids = vec![
        scored_pep("AAAAAAAAAAAAAK", "FDR", 0.05),
        PeptideIdentification::default(),
        PeptideIdentification {
            score_type: "FDR".into(),
            ..PeptideIdentification::default()
        },
    ];
    let r = IdentificationSummary::new()
        .compute(&prot_ids, &mut pep_ids)
        .unwrap();
    assert_eq!(r.peptide_spectrum_matches, 2);
    assert_eq!(r.unique_peptides.count, 1);
    assert!(
        IdentificationSummary::new()
            .compute(&prot_ids, &mut [])
            .is_err()
    );
}
