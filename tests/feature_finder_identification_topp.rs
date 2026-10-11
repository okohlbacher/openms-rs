// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FeatureFinderIdentificationAlgorithm` and `FeatureFinderAlgorithmMetaboIdent`
//! against the retained outputs of the upstream TOPP tests (evidence tier 1).
//!
//! The upstream suite pins both algorithms only through TOPP tools:
//!
//! * `TOPP_FeatureFinderIdentification_1` (`-model:type none`), `_3`
//!   (`-model:type symmetric`) and `_5` (`-extract:batch_size 10` with
//!   `-candidates_out`), `topp` test data `0cb15f2`, `CMakeLists.txt:434-459`;
//!   `_4` is disabled upstream and is not asserted.
//! * `TOPP_FeatureFinderIdentification_summary_count`, which expects the text
//!   `22 peptides identified` on the tool's standard output.
//! * `TOPP_FeatureFinderMetaboIdent_1`, `CMakeLists.txt:4267-4270`.
//!
//! Each test performs what the tool does around the algorithm (load the MS1
//! spectra and the identifications or targets, run, store) and compares the
//! stored file with the retained expectation under the upstream rule of
//! `FuzzyDiff` (ratio `1.01` or absolute difference `0.01`), decoded rather than
//! line by line, the `spectra_data` meta value skipped as the upstream
//! `-whitelist "spectra_data"` does and generated identifiers skipped as its
//! `-whitelist "feature id"` does. The data processing records are the tool
//! wrapper's, not the algorithm's, and are not compared.
//!
//! The tight comparison against a run of the C++ Release build is
//! `tests/feature_finder_identification_oracle.rs`.
#![cfg(all(feature = "mzml", feature = "featurexml", feature = "idxml"))]

#[path = "support/decoded_compare.rs"]
mod decoded;
#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

use openms::analysis::featurefinder::feature_finder_identification::FeatureFinderIdentificationAlgorithm;
use openms::analysis::featurefinder::feature_finder_metabo_ident::{
    FeatureFinderAlgorithmMetaboIdent, FeatureFinderMetaboIdentCompound,
};
use openms::concept::UniqueIdGenerator;
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::format::peak_options::PeakFileOptions;
use openms::kernel::{FeatureMap, MSExperiment};
use openms::param::ParamValue;
use openms::system::file::TempDir;
use std::fs;
use std::path::{Path, PathBuf};

/// The seed of the unique ID generator under `-test`.
const TEST_MODE_UNIQUE_ID_SEED: u64 = 19_991_231_235_959;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/feature_finder_identification")
        .join(name)
}

fn text(path: impl AsRef<Path>) -> String {
    path.as_ref().to_string_lossy().into_owned()
}

/// The MS1 spectra of an mzML file, as the tools load them
/// (`mzml.getOptions().addMSLevel(1)`).
fn load_ms1(name: &str) -> MSExperiment {
    let mut options = PeakFileOptions::default();
    options.add_ms_level(1).unwrap();
    FileHandler::load_experiment_with_options(data(name), &[FileType::MzMl], &options).unwrap()
}

fn expected(name: &str) -> FeatureMap {
    FileHandler::load_feature_map(data(name), &[FileType::FeatureXml]).unwrap()
}

/// A feature map as the tool leaves it in its output file.
fn through_file(map: &FeatureMap, work: &TempDir, name: &str) -> FeatureMap {
    let path = work.path().join(name);
    FileHandler::store_feature_map(&path, map, Some(FileType::FeatureXml)).unwrap();
    FileHandler::load_feature_map(&path, &[FileType::FeatureXml]).unwrap()
}

/// The upstream rule: `FuzzyDiff` ratio 1.01 / absolute 0.01, ids skipped; the
/// run path and the processing records are not the algorithm's.
fn assert_matches(actual: &FeatureMap, expected: &FeatureMap, what: &str) {
    assert_matches_with(actual, expected, what, true);
}

/// As [`assert_matches`]; `skip_run_path` false compares `spectra_data` too.
fn assert_matches_with(
    actual: &FeatureMap,
    expected: &FeatureMap,
    what: &str,
    skip_run_path: bool,
) {
    let mut actual = actual.clone();
    let mut expected = expected.clone();
    for map in [&mut actual, &mut expected] {
        if skip_run_path {
            map.metadata.remove("spectra_data");
        }
        map.data_processing.clear();
    }
    let options =
        decoded::DecodedOptions::new(decoded::Tolerance::new(1.01, 0.01)).ignoring_unique_ids();
    if let Err(mismatch) = decoded::compare_feature_maps(&actual, &expected, &options) {
        panic!("{what}: {mismatch}");
    }
}

fn ffid_parameters(overrides: &[(&str, ParamValue)]) -> openms::param::Param {
    let mut param = FeatureFinderIdentificationAlgorithm::defaults().unwrap();
    // the tool does not register `candidates_out` of the algorithm's own
    // section; the case sets it itself
    for (key, value) in overrides {
        param.set_value(key, value.clone(), "", &[]).unwrap();
    }
    param
}

struct FfidRun {
    features: FeatureMap,
    algorithm: FeatureFinderIdentificationAlgorithm,
}

/// What `FeatureFinderIdentification` does for `-in/-id/-out`.
fn run_ffid(overrides: &[(&str, ParamValue)]) -> FfidRun {
    let param = ffid_parameters(overrides);
    let mut algorithm = FeatureFinderIdentificationAlgorithm::from_param(&param).unwrap();
    algorithm.set_unique_id_generator(UniqueIdGenerator::from_seed(TEST_MODE_UNIQUE_ID_SEED));
    algorithm.set_ms_data(load_ms1("FeatureFinderIdentification_1_input.mzML.gz"));
    let document = FileHandler::load_identifications(
        data("FeatureFinderIdentification_1_input.idXML"),
        &[FileType::IdXml],
    )
    .unwrap();
    let mut features = FeatureMap::new();
    algorithm
        .run(
            document.peptide_identifications,
            &document.protein_identifications,
            &mut features,
            &FeatureMap::new(),
            "FeatureFinderIdentification_1_input.mzML",
        )
        .unwrap();
    FfidRun {
        features,
        algorithm,
    }
}

fn float(value: f64) -> ParamValue {
    ParamValue::Float(value)
}

fn string(value: &str) -> ParamValue {
    ParamValue::String(value.to_owned())
}

#[test]
fn topp_1_model_none_matches_the_retained_output() {
    // -extract:mz_window 0.1 -detect:peak_width 60 -model:type none
    let run = run_ffid(&[
        ("extract:mz_window", float(0.1)),
        ("detect:peak_width", float(60.0)),
        ("model:type", string("none")),
    ]);
    let work = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let actual = through_file(&run.features, &work, "FeatureFinderIdentification_1.featureXML");
    assert_matches(
        &actual,
        &expected("FeatureFinderIdentification_1_output.featureXML"),
        "TOPP_FeatureFinderIdentification_1",
    );
}

#[test]
fn topp_summary_count_reports_22_identified_peptides() {
    // TOPP_FeatureFinderIdentification_summary_count:
    // PASS_REGULAR_EXPRESSION "22 peptides identified"
    let run = run_ffid(&[
        ("extract:mz_window", float(0.1)),
        ("detect:peak_width", float(60.0)),
        ("model:type", string("none")),
    ]);
    assert_eq!(run.algorithm.report().peptides_identified, 22);
    assert!(
        run.algorithm
            .report()
            .messages
            .iter()
            .any(|line| line.contains("22 peptides identified"))
    );
}

#[test]
fn topp_3_symmetric_model_matches_the_retained_output() {
    // -extract:mz_window 0.1 -detect:peak_width 60 -model:type symmetric
    let run = run_ffid(&[
        ("extract:mz_window", float(0.1)),
        ("detect:peak_width", float(60.0)),
        ("model:type", string("symmetric")),
    ]);
    let work = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let actual = through_file(&run.features, &work, "FeatureFinderIdentification_3.featureXML");
    assert_matches(
        &actual,
        &expected("FeatureFinderIdentification_3_output.featureXML"),
        "TOPP_FeatureFinderIdentification_3",
    );
}

#[test]
fn topp_5_batches_and_candidates_match_the_retained_outputs() {
    // -extract:mz_window 0.1 -extract:batch_size 10 -detect:peak_width 60
    // -model:type none -candidates_out <file>
    let work = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let candidates = work.path().join("candidates.featureXML");
    let run = run_ffid(&[
        ("extract:mz_window", float(0.1)),
        ("extract:batch_size", ParamValue::Integer(10)),
        ("detect:peak_width", float(60.0)),
        ("model:type", string("none")),
        ("candidates_out", string(&text(&candidates))),
    ]);
    // the final output equals the unbatched one (TOPP_5_out1 compares with the
    // output of _1)
    let actual = through_file(&run.features, &work, "FeatureFinderIdentification_5.featureXML");
    assert_matches(
        &actual,
        &expected("FeatureFinderIdentification_1_output.featureXML"),
        "TOPP_FeatureFinderIdentification_5_out1",
    );
    // and the candidates file equals the retained candidates (TOPP_5_out2)
    let written = FileHandler::load_feature_map(&candidates, &[FileType::FeatureXml]).unwrap();
    assert_matches(
        &written,
        &expected("FeatureFinderIdentification_5_candidates.featureXML"),
        "TOPP_FeatureFinderIdentification_5_out2",
    );
}

#[test]
fn candidates_of_the_retained_run_reprocess_to_the_retained_output() {
    // `-candidates_in`: only the filtering and the model fit run. The retained
    // candidates of TOPP_5 came from a run with `-model:type none`, so the
    // result is the retained output of TOPP_1.
    let mut features = expected("FeatureFinderIdentification_5_candidates.featureXML");
    let param = ffid_parameters(&[("model:type", string("none"))]);
    let mut algorithm = FeatureFinderIdentificationAlgorithm::from_param(&param).unwrap();
    algorithm.run_on_candidates(&mut features).unwrap();
    assert_matches(
        &features,
        &expected("FeatureFinderIdentification_1_output.featureXML"),
        "runOnCandidates",
    );
    assert_eq!(algorithm.report().peptides_identified, 22);
}

/// The targets of a TOPP table, as `readTargets_` parses a line.
fn read_targets(path: &Path) -> Vec<FeatureFinderMetaboIdentCompound> {
    let content = fs::read_to_string(path).unwrap();
    let mut lines = content.lines();
    let header = lines.next().unwrap();
    assert!(header.starts_with(
        "CompoundName\tSumFormula\tMass\tCharge\tRetentionTime\tRetentionTimeRange\tIsoDistribution"
    ));
    let list = |text: &str| -> Vec<f64> {
        text.split(',')
            .filter(|part| !part.trim().is_empty())
            .map(|part| part.trim().parse::<f64>().unwrap())
            .collect()
    };
    let mut names = std::collections::BTreeSet::new();
    let mut targets = Vec::new();
    for line in lines {
        if line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 7 || parts[0].is_empty() || !names.insert(parts[0].to_owned()) {
            continue;
        }
        targets.push(FeatureFinderMetaboIdentCompound::new(
            parts[0],
            parts[1],
            parts[2].parse::<f64>().unwrap(),
            list(parts[3]).into_iter().map(|c| c as i32).collect(),
            list(parts[4]),
            list(parts[5]),
            list(parts[6]),
            Vec::new(),
            "",
        ));
    }
    targets
}

#[test]
fn metabo_ident_1_matches_the_retained_output() {
    // -extract:mz_window 5 -extract:rt_window 20 -detect:peak_width 3
    let mut param = FeatureFinderAlgorithmMetaboIdent::defaults().unwrap();
    param
        .set_value("extract:mz_window", float(5.0), "", &[])
        .unwrap();
    param
        .set_value("extract:rt_window", float(20.0), "", &[])
        .unwrap();
    param
        .set_value("detect:peak_width", float(3.0), "", &[])
        .unwrap();
    let mut algorithm = FeatureFinderAlgorithmMetaboIdent::from_param(&param).unwrap();
    algorithm.set_unique_id_generator(UniqueIdGenerator::from_seed(TEST_MODE_UNIQUE_ID_SEED));
    algorithm.set_ms_data(load_ms1("FeatureFinderMetaboIdent_1_input.mzML.gz"));
    let targets = read_targets(&data("FeatureFinderMetaboIdent_1_input.tsv"));
    let mut features = FeatureMap::new();
    algorithm
        .run(&targets, &mut features, "FeatureFinderMetaboIdent_1_input.mzML")
        .unwrap();
    // -test: the run path without directory
    features
        .set_primary_ms_run_path(&["file://FeatureFinderMetaboIdent_1_input.mzML".to_owned()])
        .unwrap();
    let work = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let actual = through_file(&features, &work, "FeatureFinderMetaboIdent_1.featureXML");
    assert_matches_with(
        &actual,
        &expected("FeatureFinderMetaboIdent_1_output.featureXML"),
        "TOPP_FeatureFinderMetaboIdent_1",
        false,
    );
}
