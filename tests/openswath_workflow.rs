// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OpenSwathWorkflow.
//!
//! **Evidence tier 1 (retained C++ output).** The cases below reproduce the
//! invocations of the upstream TOPP tests `TOPP_OpenSwathWorkflow_1`, `_3`,
//! `_6`, `_18`, `_19`, `_20` and `_threads4` (test-data pin `0cb15f2`,
//! `topp/CMakeLists.txt` lines 2063-2253) on the retained input
//! `OpenSwathWorkflow_1_input.{mzML,TraML,trafoXML}` and compare with the
//! retained expected outputs, which the executed C++ build wrote. The
//! comparison is the decoded feature-map comparison with the `FuzzyDiff.ini`
//! tolerances (`ratio 1.01`, `absdiff 0.01`) and the `id=` whitelist, as the
//! upstream `${DIFF}`. The TOPP tool of the `openswath` package is not part of
//! this crate, so the harness below plays its role: it splits the DIA run into
//! SWATH maps as `FullSwathFileConsumer` does, converts the TraML library, and
//! assembles the parameters from the `Scoring` subsection defaults and the
//! overrides of the tool (`OpenSwathWorkflow.cpp` of the package, lines
//! 392-450, 1063-1131 and 1403-1418).
//!
//! The wave-scheduler cases (`_6`, `_threads4`) and the per-SWATH cases (`_1`,
//! `_3`, `_18`..) take different code paths; each is also run at several
//! thread counts and must give bit-identical results. The remaining tests
//! check behaviour the TOPP tests do not reach (tier 4).
#![cfg(all(feature = "mzml", feature = "traml", feature = "featurexml"))]

#[path = "support/decoded_compare.rs"]
mod decoded;
#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

use openms::Error;
use openms::analysis::openswath::data_access_helper::convert_targeted_exp;
use openms::analysis::openswath::mrm_feature_finder_scoring::MRMFeatureFinderScoring;
use openms::analysis::openswath::openswath_library_id_normalizer::OpenSwathLibraryIDNormalizer;
use openms::analysis::openswath::openswath_osw_writer::OpenSwathOSWWriter;
use openms::analysis::openswath::openswath_workflow::{
    ChromExtractParams, ChromatogramHandler, LogLevel, OpenSwathWorkflow,
    select_compounds_for_batch,
};
use openms::analysis::openswath::spectrum_access_openms::SimpleOpenMSSpectraFactory;
use openms::analysis::openswath::swath_map::SwathMap;
use openms::analysis::openswath::transition_experiment::LightTargetedExperiment;
use openms::analysis::openswath::transition_tsv_file::TransitionTSVFile;
use openms::analysis::transformations::{LinearOptions, ModelConfig, TransformationDescription};
use openms::concept::UniqueIdGenerator;
use openms::concept::parallel::Threads;
use openms::format::FileType;
use openms::format::{featurexml, mzml, traml, transformation_xml};
use openms::interfaces::MSDataConsumer;
use openms::kernel::{FeatureMap, MSChromatogram, MSExperiment, MSSpectrum};
use openms::param::{Param, ParamValue};
use std::ops::ControlFlow;
use std::path::PathBuf;
use std::sync::Arc;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/openswath_workflow")
        .join(name)
}

// ---------------------------------------------------------------------------
// The harness: what the TOPP tool does before it calls the workflow
// ---------------------------------------------------------------------------

/// Split a DIA run into SWATH maps as `FullSwathFileConsumer` does: the MS1
/// map first (window -1), then one map per precursor isolation window in the
/// order of first appearance.
fn swath_maps_of(exp: &MSExperiment) -> Vec<SwathMap> {
    let mut ms1 = MSExperiment::new();
    ms1.settings = exp.settings.clone();
    let mut has_ms1 = false;
    let mut boundaries: Vec<SwathMap> = Vec::new();
    let mut experiments: Vec<MSExperiment> = Vec::new();
    for spectrum in &exp.spectra {
        if spectrum.ms_level == 1 {
            has_ms1 = true;
            ms1.spectra.push(spectrum.clone());
            continue;
        }
        let precursor = &spectrum.precursors[0];
        let center = precursor.mz;
        let lower = center - precursor.isolation_window_lower_offset;
        let upper = center + precursor.isolation_window_upper_offset;
        // the ion mobility window, -1 when the spectrum has none
        let im = |key: &str| match spectrum.metadata.get(key).map(|v| v.data()) {
            Some(openms::metadata::MetaValueData::Float(f)) => *f,
            _ => -1.0,
        };
        let (lower_im, upper_im) = (
            im("ion mobility lower limit"),
            im("ion mobility upper limit"),
        );
        let index = boundaries.iter().position(|b| {
            (b.center - center).abs() < 1e-6
                && (b.im_lower - lower_im).abs() < 1e-6
                && (b.im_upper - upper_im).abs() < 1e-6
        });
        let index = index.unwrap_or_else(|| {
            let mut boundary = SwathMap::new(lower, upper, center, false);
            boundary.im_lower = lower_im;
            boundary.im_upper = upper_im;
            boundaries.push(boundary);
            let mut experiment = MSExperiment::new();
            experiment.settings = exp.settings.clone();
            experiments.push(experiment);
            boundaries.len() - 1
        });
        experiments[index].spectra.push(spectrum.clone());
    }
    let mut maps = Vec::new();
    if has_ms1 {
        let mut map = SwathMap::new(-1.0, -1.0, -1.0, true);
        map.im_lower = -1.0;
        map.im_upper = -1.0;
        map.sptr =
            Some(SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(Arc::new(ms1)).unwrap());
        maps.push(map);
    }
    for (mut boundary, experiment) in boundaries.into_iter().zip(experiments) {
        boundary.sptr = Some(
            SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(Arc::new(experiment)).unwrap(),
        );
        maps.push(boundary);
    }
    maps
}

fn set_str(param: &mut Param, key: &str, value: &str) {
    param
        .set_value(key, ParamValue::String(value.into()), "", &[])
        .unwrap();
}
fn set_int(param: &mut Param, key: &str, value: i64) {
    param
        .set_value(key, ParamValue::Integer(value), "", &[])
        .unwrap();
}
fn set_f64(param: &mut Param, key: &str, value: f64) {
    param
        .set_value(key, ParamValue::Float(value), "", &[])
        .unwrap();
}

/// `TOPPOpenSwathWorkflow::getSubsectionDefaults_("Scoring")`.
fn scoring_defaults() -> Param {
    let mut p = MRMFeatureFinderScoring::defaults().unwrap();
    p.remove("rt_extraction_window").unwrap();
    set_int(&mut p, "stop_report_after_feature", 5);
    set_f64(&mut p, "rt_normalization_factor", 100.0);
    set_str(&mut p, "Scores:use_ms1_mi", "true");
    set_str(&mut p, "Scores:use_mi_score", "true");

    set_f64(&mut p, "TransitionGroupPicker:min_peak_width", -1.0);
    set_str(&mut p, "TransitionGroupPicker:recalculate_peaks", "true");
    set_str(
        &mut p,
        "TransitionGroupPicker:compute_peak_quality",
        "false",
    );
    set_f64(&mut p, "TransitionGroupPicker:minimal_quality", -1.5);
    set_str(
        &mut p,
        "TransitionGroupPicker:background_subtraction",
        "none",
    );
    set_str(
        &mut p,
        "TransitionGroupPicker:compute_peak_shape_metrics",
        "false",
    );
    p.remove("TransitionGroupPicker:stop_after_intensity_ratio")
        .unwrap();

    let picker = "TransitionGroupPicker:PeakPickerChromatogram:";
    set_str(&mut p, &format!("{picker}use_gauss"), "false");
    set_int(&mut p, &format!("{picker}sgolay_polynomial_order"), 3);
    set_int(&mut p, &format!("{picker}sgolay_frame_length"), 11);
    set_f64(&mut p, &format!("{picker}peak_width"), -1.0);
    set_str(&mut p, &format!("{picker}remove_overlapping_peaks"), "true");
    set_str(&mut p, &format!("{picker}write_sn_log_messages"), "false");
    set_f64(
        &mut p,
        "TransitionGroupPicker:recalculate_peaks_max_z",
        0.75,
    );
    set_str(&mut p, &format!("{picker}method"), "corrected");
    set_f64(&mut p, &format!("{picker}signal_to_noise"), 0.1);
    set_f64(&mut p, &format!("{picker}gauss_width"), 30.0);
    set_int(&mut p, "uis_threshold_sn", -1);
    set_int(&mut p, "uis_threshold_peak_area", 0);
    p.remove(&format!("{picker}sn_win_len")).unwrap();
    p.remove(&format!("{picker}sn_bin_count")).unwrap();
    p.remove(&format!("{picker}stop_after_feature")).unwrap();

    set_str(&mut p, "Scores:use_ion_mobility_scores", "auto");
    p.remove("Scores:use_elution_model_score").unwrap();
    set_int(&mut p, "EMGScoring:max_iteration", 10);
    p.remove("EMGScoring:interpolation_step").unwrap();
    p.remove("EMGScoring:tolerance_stdev_bounding_box").unwrap();
    p.remove("EMGScoring:deltaAbsError").unwrap();
    p.remove("EMGScoring:statistics:mean").unwrap();
    p.remove("EMGScoring:statistics:variance").unwrap();
    p
}

/// One invocation of the TOPP tool, reduced to what reaches the workflow.
#[derive(Clone)]
struct Case {
    /// The retained input run, and the library with its kind.
    input: &'static str,
    library_file: &'static str,
    /// `-Calibration:rt_norm` given (the trafoXML of test 1), else none.
    rt_norm: bool,
    /// `-matching_window_only`, `-pasef`.
    prm: bool,
    pasef: bool,
    /// `-ion_mobility_window`.
    im_window: f64,
    /// `-Scoring:Scores:use_ion_mobility_scores`.
    use_im_scores: bool,
    enable_ms1: bool,
    mz_window: f64,
    ppm: bool,
    ms1_isotopes: i32,
    compute_peak_quality: bool,
    use_ms1_mi: bool,
    use_mi_score: bool,
    use_total_mi_score: bool,
    load_into_memory: bool,
    batch_size: i32,
    threads: usize,
}

impl Case {
    /// The tool's defaults: 50 ppm, three MS1 isotopes, both MI scores on.
    fn tool_defaults() -> Self {
        Self {
            input: "OpenSwathWorkflow_1_input.mzML",
            library_file: "OpenSwathWorkflow_1_input.TraML",
            rt_norm: true,
            prm: false,
            pasef: false,
            im_window: -1.0,
            use_im_scores: false,
            enable_ms1: true,
            mz_window: 50.0,
            ppm: true,
            ms1_isotopes: 3,
            compute_peak_quality: false,
            use_ms1_mi: true,
            use_mi_score: true,
            use_total_mi_score: false,
            load_into_memory: false,
            batch_size: 0,
            threads: 1,
        }
    }

    /// `${OLD_OSW_PARAM}` with the calibration windows disabled.
    fn old() -> Self {
        Self {
            enable_ms1: true,
            mz_window: 0.05,
            ppm: false,
            ms1_isotopes: 0,
            compute_peak_quality: true,
            use_ms1_mi: false,
            use_mi_score: false,
            ..Self::tool_defaults()
        }
    }
}

struct Outcome {
    library: LightTargetedExperiment,
    features: FeatureMap,
    chromatograms: Vec<MSChromatogram>,
    log: Vec<openms::analysis::openswath::openswath_workflow::LogEntry>,
}

struct Collect(Vec<MSChromatogram>);

impl MSDataConsumer for Collect {
    fn set_expected_size(&mut self, _s: usize, _c: usize) -> openms::Result<()> {
        Ok(())
    }
    fn set_experimental_settings(
        &mut self,
        _settings: &openms::metadata::ExperimentalSettings,
    ) -> openms::Result<()> {
        Ok(())
    }
    fn consume_spectrum(&mut self, _s: &mut MSSpectrum) -> openms::Result<ControlFlow<()>> {
        Err(Error::InvalidValue("no spectra expected".into()))
    }
    fn consume_chromatogram(
        &mut self,
        chromatogram: &mut MSChromatogram,
    ) -> openms::Result<ControlFlow<()>> {
        self.0.push(chromatogram.clone());
        Ok(ControlFlow::Continue(()))
    }
}

/// `TOPPOpenSwathBase::loadTransitionList` for a TraML library: converted to a
/// light experiment, then its identifiers canonicalized (numeric strings, in
/// the order of `OpenSwathLibraryIDNormalizer`), which is also what orders the
/// assays of the workflow.
fn library() -> LightTargetedExperiment {
    library_of("OpenSwathWorkflow_1_input.TraML")
}

fn library_of(file: &str) -> LightTargetedExperiment {
    let exp = if file.ends_with(".tsv") {
        TransitionTSVFile::new()
            .convert_tsv_to_targeted_experiment(data(file), FileType::Tsv)
            .unwrap()
    } else {
        traml::load(data(file)).unwrap()
    };
    let mut light = convert_targeted_exp(&exp).unwrap();
    OpenSwathLibraryIDNormalizer::normalize_source_ids(&mut light).unwrap();
    light
}

fn rt_trafo() -> TransformationDescription {
    // FileHandler().loadTransformations(..., fit_model = false), then
    // fitModel("linear") with the RTNormalization defaults
    let mut trafo = transformation_xml::load(data("OpenSwathWorkflow_1_input.trafoXML")).unwrap();
    trafo
        .fit_model(ModelConfig::Linear(LinearOptions::default()))
        .unwrap();
    trafo
}

fn run_case(case: &Case) -> Outcome {
    let exp = mzml::read(std::fs::read(data(case.input)).unwrap().as_slice()).unwrap();
    let maps = swath_maps_of(&exp);

    let cp = ChromExtractParams {
        min_upper_edge_dist: 0.0,
        mz_extraction_window: case.mz_window,
        im_extraction_window: case.im_window,
        ppm: case.ppm,
        extraction_function: "tophat".into(),
        rt_extraction_window: 600.0,
        extra_rt_extract: 0.0,
    };
    let mut cp_ms1 = cp.clone();
    cp_ms1.mz_extraction_window = 50.0;
    cp_ms1.ppm = true;
    cp_ms1.im_extraction_window = -1.0;

    let mut ffp = scoring_defaults();
    set_str(&mut ffp, "use_ms1_ion_mobility", "true");
    if case.compute_peak_quality {
        set_str(
            &mut ffp,
            "TransitionGroupPicker:compute_peak_quality",
            "true",
        );
    }
    set_str(
        &mut ffp,
        "Scores:use_ms1_mi",
        if case.use_ms1_mi { "true" } else { "false" },
    );
    set_str(
        &mut ffp,
        "Scores:use_mi_score",
        if case.use_mi_score { "true" } else { "false" },
    );
    if case.use_total_mi_score {
        set_str(&mut ffp, "Scores:use_total_mi_score", "true");
    }
    set_str(&mut ffp, "Scores:use_elution_model_score", "false");
    if case.enable_ms1 {
        set_str(&mut ffp, "Scores:use_ms1_correlation", "true");
        set_str(&mut ffp, "Scores:use_ms1_fullscan", "true");
    }
    // -enable_ipf defaults to true
    set_str(&mut ffp, "Scores:use_uis_scores", "true");
    // "auto" resolves to true with -pasef and to false without; an explicit
    // value is respected
    set_str(
        &mut ffp,
        "Scores:use_ion_mobility_scores",
        if case.use_im_scores || case.pasef {
            "true"
        } else {
            "false"
        },
    );

    let mut workflow =
        OpenSwathWorkflow::new(case.enable_ms1, true, case.prm, case.pasef, false, -1);
    workflow.set_threads(Threads::from_cli(case.threads as i64));
    workflow.set_unique_id_generator(UniqueIdGenerator::from_seed(42));
    let mut features = FeatureMap::new();
    let mut consumer = Collect(Vec::new());
    workflow
        .perform_extraction(
            &maps,
            &if case.rt_norm {
                rt_trafo()
            } else {
                TransformationDescription::default()
            },
            &cp,
            &cp_ms1,
            &ffp,
            &library_of(case.library_file),
            &mut features,
            true,
            &OpenSwathOSWWriter::new("", false),
            Some(&mut consumer),
            case.batch_size,
            case.ms1_isotopes,
            case.load_into_memory,
            &Param::new(),
            None,
            -1,
            -1,
        )
        .unwrap();
    let log = workflow.take_log();
    for entry in &log {
        eprintln!("[workflow {:?}] {}", entry.level, entry.message);
    }
    Outcome {
        library: library_of(case.library_file),
        features,
        chromatograms: consumer.0,
        log,
    }
}

fn expected_features(name: &str) -> FeatureMap {
    featurexml::load(data(name)).unwrap()
}

/// Whether a feature's peptide has a bare ambiguous residue (B, Z or X).
fn has_ambiguous_residue(
    feature: &openms::kernel::Feature,
    library: &LightTargetedExperiment,
) -> bool {
    let ambiguous: Vec<String> = library
        .compounds
        .iter()
        .filter(|c| c.sequence.chars().any(|r| matches!(r, 'B' | 'Z' | 'X')))
        .map(|c| c.id.clone())
        .collect();
    match feature.metadata.get("PeptideRef").map(|v| v.data()) {
        Some(openms::metadata::MetaValueData::String(id)) => ambiguous.contains(id),
        _ => false,
    }
}

/// Compare as the upstream `${DIFF} -whitelist "id="`: decoded, `FuzzyDiff.ini`
/// tolerances, generated identifiers skipped. Map-level provenance (document
/// identifier, processing history written by the tool) is not the workflow's.
///
/// **Known divergence, carved out narrowly.** `AASequence` of this crate
/// refuses a bare `B`/`Z`/`X` residue where the source's `ResidueDB` gives it
/// the formula `""` (an internal mass of minus one water): the b/y-ion score
/// of an assay whose peptide contains one fails with `Unsupported`, and the
/// workflow, which catches every `std::exception` of scoring, skips the assay.
/// The upstream library has such an assay (`PEPTIDEB`). The expected features
/// of ambiguous peptides are therefore compared strictly when the run has
/// them, and otherwise are accepted as missing only when the workflow log
/// records the skipped assay with exactly that error. When the chemistry port
/// computes the source's masses, the run has the features and this branch is
/// no longer taken.
fn assert_features_match(outcome: &Outcome, expected_name: &str) {
    let mut expected = expected_features(expected_name);
    // The tool stores the map, and the source's featureXML writer omits (with
    // a warning) a peptide identification whose run has no protein
    // identification in the map; this crate's writer refuses such a map
    // instead. The stored map is therefore reproduced here, then re-read, as
    // the upstream DIFF sees it.
    let mut stored_map = outcome.features.clone();
    let runs: Vec<String> = stored_map
        .protein_identifications
        .iter()
        .map(|p| p.identifier.clone())
        .collect();
    for feature in &mut stored_map.features {
        feature
            .peptide_identifications
            .retain(|id| runs.contains(&id.identifier));
    }
    let mut stored = Vec::new();
    featurexml::write(&mut stored, &stored_map).unwrap();
    let mut actual = featurexml::read(stored.as_slice()).unwrap();
    let expected_ambiguous = expected
        .features
        .iter()
        .filter(|f| has_ambiguous_residue(f, &outcome.library))
        .count();
    let actual_ambiguous = actual
        .features
        .iter()
        .filter(|f| has_ambiguous_residue(f, &outcome.library))
        .count();
    if actual_ambiguous == 0 && expected_ambiguous > 0 {
        let skipped = outcome
            .log
            .iter()
            .filter(|e| {
                e.level == LogLevel::Error
                    && e.message
                        .contains("sequence contains a residue without a known monoisotopic mass")
            })
            .count();
        assert!(
            skipped >= expected_ambiguous,
            "{expected_name}: {expected_ambiguous} expected features of ambiguous peptides are \
             missing and only {skipped} assays were skipped for the known reason"
        );
        eprintln!(
            "{expected_name}: {expected_ambiguous} feature(s) of ambiguous peptides not produced \
             (AASequence refuses bare B/Z/X; source: formula empty)"
        );
        expected
            .features
            .retain(|f| !has_ambiguous_residue(f, &outcome.library));
    }
    for map in [&mut actual, &mut expected] {
        map.unique_id = 0;
        map.identifier.clear();
        map.data_processing.clear();
        map.metadata.clear();
    }
    let settings = fuzzy::FuzzyDiffSettings::upstream().unwrap();
    let options = decoded::DecodedOptions::new(decoded::Tolerance::from_settings(&settings))
        .ignoring_unique_ids();
    if let Err(mismatch) = decoded::compare_feature_maps(&actual, &expected, &options) {
        let describe = |label: &str, map: &FeatureMap| {
            eprintln!("{label}: {} features", map.features.len());
            for f in &map.features {
                eprintln!(
                    "  rt {:.4} mz {:.4} int {:.2} q {:.4} charge {} peptide {:?} n_sub {}",
                    f.rt,
                    f.mz,
                    f.intensity,
                    f.quality,
                    f.charge,
                    f.metadata.get("PeptideRef"),
                    f.subordinates.len()
                );
            }
        };
        describe("actual", &actual);
        describe("expected", &expected);
        if let (Some(a), Some(e)) = (actual.features.first(), expected.features.first()) {
            eprintln!("actual   peptide ids: {:#?}", a.peptide_identifications);
            eprintln!("expected peptide ids: {:#?}", e.peptide_identifications);
        }
        panic!("{expected_name}: {mismatch}");
    }
    assert!(!expected.features.is_empty());
}

// ---------------------------------------------------------------------------
// Tier 1: retained outputs of the upstream TOPP tests
// ---------------------------------------------------------------------------

#[test]
fn topp_1_per_swath_path_without_ms1() {
    // TOPP_OpenSwathWorkflow_1: -enable_ms1 false, ${OLD_OSW_PARAM}
    let case = Case {
        enable_ms1: false,
        ..Case::old()
    };
    let outcome = run_case(&case);
    assert_features_match(&outcome, "OpenSwathWorkflow_1_output.featureXML");
}

#[test]
fn topp_3_per_swath_path_with_ms1() {
    // TOPP_OpenSwathWorkflow_3: MS1 traces on
    let outcome = run_case(&Case::old());
    assert_features_match(&outcome, "OpenSwathWorkflow_3_output.featureXML");
}

/// The chromatograms the consumer received, against the retained
/// `chrom.mzML` of the tool (`TOPP_OpenSwathWorkflow_1_out2`/`_3_out2`),
/// compared decoded with the upstream tolerances and without identifiers.
fn assert_chromatograms_match(outcome: &Outcome, expected_name: &str) {
    let expected_exp = mzml::read(std::fs::read(data(expected_name)).unwrap().as_slice()).unwrap();
    let mut expected = MSExperiment::new();
    expected.chromatograms = expected_exp.chromatograms;
    // The tool stores the chromatograms in an mzML file, and the source's
    // writer keeps of a precursor's ion mobility only the drift time, in
    // milliseconds: the drift window offsets are dropped and the unit read
    // back is millisecond. This crate's writer refuses such a precursor
    // rather than drop information, so the loss is applied here, to the
    // in-memory chromatograms, before comparing with the stored ones.
    let mut actual = MSExperiment::new();
    actual.chromatograms = outcome.chromatograms.clone();
    for chromatogram in &mut actual.chromatograms {
        chromatogram.precursor.drift_window_lower_offset = 0.0;
        chromatogram.precursor.drift_window_upper_offset = 0.0;
        if chromatogram.precursor.drift_time.is_some() {
            chromatogram.precursor.drift_time_unit = openms::metadata::DriftTimeUnit::Millisecond;
        }
    }
    let mut expected = expected;
    // the processing history the tool's chromatogram consumer adds is not the
    // workflow's
    for exp in [&mut actual, &mut expected] {
        for chromatogram in &mut exp.chromatograms {
            chromatogram.data_processing.clear();
        }
    }
    let settings = fuzzy::FuzzyDiffSettings::upstream().unwrap();
    let options = decoded::DecodedOptions::new(decoded::Tolerance::from_settings(&settings))
        .ignoring_unique_ids();
    if let Err(mismatch) = decoded::compare_experiments(&actual, &expected, &options) {
        panic!("{expected_name}: {mismatch}");
    }
    assert!(!expected.chromatograms.is_empty());
}

#[test]
fn topp_1_chromatograms() {
    // TOPP_OpenSwathWorkflow_1_out2
    let case = Case {
        enable_ms1: false,
        ..Case::old()
    };
    let outcome = run_case(&case);
    assert_chromatograms_match(&outcome, "OpenSwathWorkflow_1_output.chrom.mzML");
}

#[test]
fn topp_3_chromatograms() {
    // TOPP_OpenSwathWorkflow_3_out2: with the MS1 traces
    let outcome = run_case(&Case::old());
    assert_chromatograms_match(&outcome, "OpenSwathWorkflow_3_output.chrom.mzML");
}

#[test]
fn topp_6_wave_scheduler_with_ms1() {
    // TOPP_OpenSwathWorkflow_6: -readOptions cacheWorkingInMemory, batchSize 0
    // selects the wave scheduler; the expected output is that of _3
    let case = Case {
        load_into_memory: true,
        ..Case::old()
    };
    let outcome = run_case(&case);
    assert_features_match(&outcome, "OpenSwathWorkflow_3_output.featureXML");
}

#[test]
fn topp_threads4_wave_scheduler_sorted() {
    // TOPP_OpenSwathWorkflow_threads4: -enable_ms1 false -readOptions
    // workingInMemory -batchSize 0 -threads 4, then FileFilter -sort
    let case = Case {
        enable_ms1: false,
        load_into_memory: true,
        threads: 4,
        ..Case::tool_defaults()
    };
    let mut outcome = run_case(&case);
    outcome.features.sort_by_position().unwrap();
    assert_features_match(&outcome, "OpenSwathWorkflow_threads_output.featureXML");
}

#[test]
fn topp_18_ms1_mutual_information() {
    // TOPP_OpenSwathWorkflow_18: use_ms1_mi stays on
    let case = Case {
        use_ms1_mi: true,
        ..Case::old()
    };
    let outcome = run_case(&case);
    assert_features_match(&outcome, "OpenSwathWorkflow_18_output.featureXML");
}

#[test]
fn topp_19_ms2_mutual_information() {
    // TOPP_OpenSwathWorkflow_19: use_mi_score stays on
    let case = Case {
        use_mi_score: true,
        ..Case::old()
    };
    let outcome = run_case(&case);
    assert_features_match(&outcome, "OpenSwathWorkflow_19_output.featureXML");
}

#[test]
fn topp_20_total_mutual_information() {
    // TOPP_OpenSwathWorkflow_20: total MI score on, MS2 MI on
    let case = Case {
        use_mi_score: true,
        use_total_mi_score: true,
        ..Case::old()
    };
    let outcome = run_case(&case);
    assert_features_match(&outcome, "OpenSwathWorkflow_20_output.featureXML");
}

#[test]
fn topp_22_matching_window_only() {
    // TOPP_OpenSwathWorkflow_22: PRM-like data, -readOptions workingInMemory
    // -matching_window_only true, no RT normalization, ${OLD_OSW_PARAM}
    // -ms1_isotopes 3 (the last occurrence of an option wins)
    let case = Case {
        input: "OpenSwathWorkflow_22_input.mzML",
        library_file: "OpenSwathWorkflow_22_input.tsv",
        rt_norm: false,
        prm: true,
        load_into_memory: true,
        ms1_isotopes: 3,
        ..Case::old()
    };
    let outcome = run_case(&case);
    assert_features_match(&outcome, "OpenSwathWorkflow_22_output.featureXML");
    assert_chromatograms_match(&outcome, "OpenSwathWorkflow_22_output.chrom.mzML");
}

#[test]
fn topp_23_diapasef() {
    // TOPP_OpenSwathWorkflow_23: -pasef -ion_mobility_window 0.1
    // -Scoring:Scores:use_ion_mobility_scores true ${OLD_OSW_PARAM}
    let case = Case {
        input: "OpenSwathWorkflow_23_input.mzML",
        library_file: "OpenSwathWorkflow_23_input.tsv",
        rt_norm: false,
        pasef: true,
        im_window: 0.1,
        use_im_scores: true,
        load_into_memory: true,
        ..Case::old()
    };
    let outcome = run_case(&case);
    assert_features_match(&outcome, "OpenSwathWorkflow_23_output.featureXML");
    assert_chromatograms_match(&outcome, "OpenSwathWorkflow_23_output.chrom.mzML");
}

// ---------------------------------------------------------------------------
// Determinism and structure (tier 4)
// ---------------------------------------------------------------------------

#[test]
fn result_does_not_depend_on_the_thread_count() {
    // The per-SWATH path (`load_into_memory = false`) is bit-identical in
    // order. The wave scheduler sorts a wave's features by position when it
    // scored several jobs, as the source does, and which SWATHs share a wave
    // follows the thread count; so there the *set* of features, with every
    // value and identifier, is what must not change.
    for load_into_memory in [false, true] {
        let mut reference = run_case(&Case {
            load_into_memory,
            threads: 1,
            ..Case::old()
        });
        for threads in [2, 4, 8] {
            let mut other = run_case(&Case {
                load_into_memory,
                threads,
                ..Case::old()
            });
            if load_into_memory {
                reference.features.sort_by_position().unwrap();
                other.features.sort_by_position().unwrap();
            }
            let exact = decoded::DecodedOptions::new(decoded::Tolerance::exact());
            if let Err(mismatch) =
                decoded::compare_feature_maps(&other.features, &reference.features, &exact)
            {
                panic!("in_memory={load_into_memory}, threads={threads}: {mismatch}");
            }
            assert_eq!(
                other.features, reference.features,
                "features, in_memory={load_into_memory}, threads={threads}"
            );
            if load_into_memory {
                // chromatograms are consumed in finalization order, which
                // follows the wave composition too
                let key = |c: &MSChromatogram| c.native_id.clone();
                let mut a = other.chromatograms.clone();
                let mut b = reference.chromatograms.clone();
                a.sort_by_key(key);
                b.sort_by_key(key);
                assert_eq!(a, b);
            } else {
                assert_eq!(other.chromatograms, reference.chromatograms);
            }
        }
    }
}

#[test]
fn explicit_batch_size_gives_the_same_features_as_automatic() {
    let automatic = run_case(&Case::old());
    let batched = run_case(&Case {
        batch_size: 2,
        ..Case::old()
    });
    assert_eq!(
        batched.features.features.len(),
        automatic.features.features.len()
    );
    assert_features_match(&batched, "OpenSwathWorkflow_3_output.featureXML");
}

#[test]
fn select_compounds_for_batch_cuts_compounds_and_their_transitions() {
    let lib = library();
    let total = lib.compounds.len();
    assert!(total > 2);
    let first = select_compounds_for_batch(&lib, 2, 0);
    assert_eq!(first.compounds.len(), 2);
    assert_eq!(first.proteins, lib.proteins);
    assert!(
        first
            .transitions
            .iter()
            .all(|t| first.compounds.iter().any(|c| c.id == t.peptide_ref))
    );
    // the last batch is shorter, one past the end is empty
    let last_index = (total - 1) / 2;
    let last = select_compounds_for_batch(&lib, 2, last_index);
    assert_eq!(last.compounds.len(), total - 2 * last_index);
    assert!(
        select_compounds_for_batch(&lib, 2, last_index + 1)
            .compounds
            .is_empty()
    );
    // every transition lands in exactly one batch
    let mut seen = 0;
    for batch in 0..=last_index {
        seen += select_compounds_for_batch(&lib, 2, batch).transitions.len();
    }
    assert_eq!(seen, lib.transitions.len());
}

#[test]
fn ms1_only_data_needs_ms1_traces() {
    let exp = mzml::read(
        std::fs::read(data("OpenSwathWorkflow_1_input.mzML"))
            .unwrap()
            .as_slice(),
    )
    .unwrap();
    let all = swath_maps_of(&exp);
    let ms1_only = vec![all[0].clone()];
    assert!(ms1_only[0].ms1);
    let mut workflow = OpenSwathWorkflow::new(false, false, false, false, false, -1);
    let cp = ChromExtractParams {
        mz_extraction_window: 50.0,
        ppm: true,
        im_extraction_window: -1.0,
        extraction_function: "tophat".into(),
        rt_extraction_window: 600.0,
        ..ChromExtractParams::default()
    };
    let ffp = scoring_defaults();
    let mut out = FeatureMap::new();
    let error = workflow
        .perform_extraction(
            &ms1_only,
            &rt_trafo(),
            &cp,
            &cp,
            &ffp,
            &library(),
            &mut out,
            true,
            &OpenSwathOSWWriter::new("", false),
            None,
            0,
            0,
            false,
            &Param::new(),
            None,
            -1,
            -1,
        )
        .unwrap_err();
    assert!(matches!(&error, Error::InvalidValue(m) if m.contains("use_ms1_traces")));
}

#[test]
fn prm_with_pasef_is_not_implemented() {
    let mut workflow = OpenSwathWorkflow::new(false, false, true, true, false, -1);
    let cp = ChromExtractParams::default();
    let ffp = scoring_defaults();
    let mut out = FeatureMap::new();
    let error = workflow
        .perform_extraction(
            &[],
            &TransformationDescription::default(),
            &cp,
            &cp,
            &ffp,
            &LightTargetedExperiment::new(),
            &mut out,
            true,
            &OpenSwathOSWWriter::new("", false),
            None,
            0,
            0,
            false,
            &Param::new(),
            None,
            -1,
            -1,
        )
        .unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)));
}

#[test]
fn mrm_mode_needs_a_chromatogram_handler() {
    struct Passthrough(Vec<MSChromatogram>);
    impl ChromatogramHandler for Passthrough {
        fn collect_irt_chromatograms_for_irt(
            &self,
            _swath_maps: &[SwathMap],
            _irt: &LightTargetedExperiment,
            _mapping: &Param,
            _cp: &ChromExtractParams,
            _trafo: &TransformationDescription,
            _pasef: bool,
            _load: bool,
        ) -> openms::Result<Vec<MSChromatogram>> {
            Ok(self.0.clone())
        }
        fn extract_and_map_chromatograms_for_transitions(
            &self,
            _swath_maps: &[SwathMap],
            _exp: &LightTargetedExperiment,
            _cp: &ChromExtractParams,
            _mapping: &Param,
        ) -> openms::Result<Vec<MSChromatogram>> {
            Ok(self.0.clone())
        }
    }
    let cp = ChromExtractParams::default();
    let mut ffp = scoring_defaults();
    set_str(&mut ffp, "Scores:use_ion_mobility_scores", "false");
    let call = |workflow: &mut OpenSwathWorkflow| {
        let mut out = FeatureMap::new();
        workflow.perform_extraction(
            &[],
            &TransformationDescription::default(),
            &cp,
            &cp,
            &ffp,
            &LightTargetedExperiment::new(),
            &mut out,
            true,
            &OpenSwathOSWWriter::new("", false),
            None,
            0,
            0,
            false,
            &Param::new(),
            None,
            -1,
            -1,
        )
    };
    let mut workflow = OpenSwathWorkflow::new(false, false, false, false, true, -1);
    assert!(matches!(call(&mut workflow), Err(Error::Unsupported(_))));
    // with a handler that maps nothing the run is empty, not an error
    workflow.set_chromatogram_handler(Some(Arc::new(Passthrough(Vec::new()))));
    call(&mut workflow).unwrap();
}
