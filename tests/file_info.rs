// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A4-FILEINFO-CORE: the FileInfo library (`FORMAT/FileInfo.h`), peak-file and
//! featureXML branches with `-m`, `-p` and `-s`, text and TSV reports.
//!
//! Evidence, in order of strength (see `tests/data/file_info_provenance.json`):
//! - tier 1, executed differential: the text and TSV reports of the product-SDK
//!   FileInfo (Debug, core `4fdec46`) from the C1 oracle
//!   (`../oracle/topp-early-bundle`: FileInfo_1/2/3/9 with `-out_tsv`, the
//!   empty featureXML and mzML, the retained FeatureFinderCentroided_1 output)
//!   and from this package's oracle (`../oracle/file-info-core`: 24 more cases),
//!   compared byte for byte with only the `File name` lines normalised;
//! - tier 1, retained upstream outputs: FileInfo_1, _2, _3 and _9 at test-data
//!   `0cb15f2`, compared as registered, with FuzzyDiff (`FuzzyDiff.ini`, ratio
//!   1.01, absdiff 0.01) and the whitelist `File name`;
//! - tier 3, the seven `FileInfo_test.cpp` sections at core `bc9cc12` that this
//!   package covers, transcribed; the consensusXML and FASTA sections are
//!   pending for their branches;
//! - tier 4, the refusal of every flag and branch not ported, and the error
//!   classes of unknown, unloadable and corrupt inputs.

#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

use openms::Error;
use openms::format::FileType;
use openms::format::file_info::model::{FileInfoResult, Options};
use openms::format::file_info::report::FileInfo;
use std::path::{Path, PathBuf};

fn data(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(relative)
}

fn read_text(relative: &str) -> String {
    std::fs::read_to_string(data(relative)).unwrap_or_else(|e| panic!("{relative}: {e}"))
}

fn mps() -> Options {
    Options {
        meta: true,
        processing: true,
        statistics: true,
        ..Options::default()
    }
}

/// Replace the value of the one `File name: ` text line and the one
/// `general: file name` TSV line, which embed the path the report was run on.
fn normalise_file_name(report: &str) -> String {
    report
        .split_inclusive('\n')
        .map(|line| {
            if line.starts_with("File name: ") {
                "File name: <input>\n".to_owned()
            } else if line.starts_with("general: file name\t") {
                "general: file name\t<input>\n".to_owned()
            } else {
                line.to_owned()
            }
        })
        .collect()
}

/// The first differing line, for a readable failure.
fn first_difference(actual: &str, expected: &str) -> String {
    let mut expected_lines = expected.split_inclusive('\n');
    for (index, line) in actual.split_inclusive('\n').enumerate() {
        match expected_lines.next() {
            Some(other) if other == line => {}
            other => {
                return format!(
                    "line {}: actual {:?}, expected {:?}",
                    index + 1,
                    line,
                    other
                );
            }
        }
    }
    match expected_lines.next() {
        Some(line) => format!("actual ends early; expected next {line:?}"),
        None => "no line differs".to_owned(),
    }
}

fn assert_report(actual: &str, expected_file: &str) {
    let expected = read_text(expected_file);
    let actual = normalise_file_name(actual);
    let expected = normalise_file_name(&expected);
    assert!(
        actual == expected,
        "{expected_file}: {}",
        first_difference(&actual, &expected)
    );
}

/// Run FileInfo and compare text and TSV with the oracle's pair of files
/// `tests/data/file_info/expected/<case>.{txt,tsv}`.
fn check_oracle(input: &str, options: &Options, case: &str) -> FileInfoResult {
    let result = FileInfo::new()
        .run(data(input), options)
        .unwrap_or_else(|e| panic!("{case}: {e}"));
    assert_report(&result.text, &format!("file_info/expected/{case}.txt"));
    assert_report(&result.tsv, &format!("file_info/expected/{case}.tsv"));
    assert_eq!(FileInfo::to_text(&result), result.text);
    assert_eq!(FileInfo::to_tsv(&result), result.tsv);
    result
}

fn dta_forced() -> Options {
    Options {
        forced_type: FileType::Dta,
        ..Options::default()
    }
}

// ---------------------------------------------------------------------------
// Tier 1: C1 oracle (../oracle/topp-early-bundle, run1)
// ---------------------------------------------------------------------------

#[test]
fn c1_file_info_1_dta_forced_type() {
    check_oracle(
        "file_info/inputs/FileInfo_1_input.dta",
        &dta_forced(),
        "FileInfo_1_tsv",
    );
}

#[test]
fn c1_file_info_2_dta2d() {
    check_oracle(
        "file_info/inputs/FileInfo_2_input.dta2d",
        &Options::default(),
        "FileInfo_2_tsv",
    );
}

#[cfg(feature = "featurexml")]
#[test]
fn c1_file_info_3_featurexml_mps() {
    check_oracle(
        "file_info/inputs/FileInfo_3_input.featureXML",
        &mps(),
        "FileInfo_3_tsv",
    );
}

// The original FileInfo_9 input does not load in the strict mzML reader, for
// three reasons outside this package, each needing a source-compatibility read
// option (D10). In the order the reader meets them:
// - it repeats the spectrum-level userParam 'name' (lines 159-160, 206-207 and
//   276-277); C++ MetaInfo overwrites the first value, the reader refuses with
//   Error::Unsupported("duplicate userParam name name");
// - its six m/z and intensity arrays carry dataProcessingRef
//   "XcaliburProcessing"; C++ keeps it on the array, the reader refuses with
//   Error::Unsupported("primary-array processing has no independent native
//   owner");
// - its 'charge array' is stored as 64-bit float; C++ loads it, the reader
//   refuses with Error::Unsupported("canonical auxiliary array binary type").
#[cfg(feature = "mzml")]
#[test]
#[ignore = "mzML reader gaps outside A4: duplicate userParam 'name', processing on primary arrays, float charge array"]
fn c1_file_info_9_mzml_mps() {
    check_oracle(
        "file_info/inputs/FileInfo_9_input.mzML",
        &mps(),
        "FileInfo_9_tsv",
    );
}

#[cfg(feature = "mzml")]
#[test]
fn c1_file_info_9_strict_reader_mps() {
    // The derived input deletes the three repeated userParams and the six
    // primary-array processing references, and stores the charge array's
    // integral values as 32-bit integers. The C++ reports on the derived and
    // the original input are identical apart from the file name
    // (../oracle/file-info-core case fileinfo_9_strict_reader_mps against C1
    // FileInfo_9_tsv), so this compares with both.
    check_oracle(
        "file_info/inputs/FileInfo_9_strict_reader.mzML",
        &mps(),
        "FileInfo_9_tsv",
    );
    check_oracle(
        "file_info/inputs/FileInfo_9_strict_reader.mzML",
        &mps(),
        "fileinfo_9_strict_reader_mps",
    );
}

#[cfg(feature = "featurexml")]
#[test]
fn c1_empty_featurexml_mps() {
    check_oracle(
        "file_info/inputs/empty.featureXML",
        &mps(),
        "FileInfo_empty_featureXML",
    );
}

#[cfg(feature = "mzml")]
#[test]
#[ignore = "mzML reader gap outside A4: the dangling spectrumList defaultDataProcessingRef \
            'dp_sp_0' is refused (unresolved dataProcessingRef); P2's source-compatible option (D10)"]
fn c1_empty_mzml_mps() {
    check_oracle("file_info/inputs/empty.mzML", &mps(), "FileInfo_empty_mzML");
}

#[cfg(feature = "mzml")]
#[test]
fn c1_empty_mzml_with_resolved_reference_mps() {
    // The derived input names the declared 'pwizconversion' instead of the
    // dangling 'dp_sp_0'. The C++ reports on both are identical apart from the
    // file name, so this compares with both.
    check_oracle(
        "file_info/inputs/empty_resolved_ref.mzML",
        &mps(),
        "FileInfo_empty_mzML",
    );
    check_oracle(
        "file_info/inputs/empty_resolved_ref.mzML",
        &mps(),
        "empty_resolved_ref_mps",
    );
}

#[cfg(feature = "featurexml")]
#[test]
fn c1_feature_finder_centroided_1_output_mps() {
    check_oracle(
        "mzml_mobility/FeatureFinderCentroided_1_1_output.featureXML",
        &mps(),
        "FileInfo_on_FFC_1_retained_output",
    );
}

// ---------------------------------------------------------------------------
// Tier 1: A4 oracle (../oracle/file-info-core, run1)
// ---------------------------------------------------------------------------

#[cfg(feature = "featurexml")]
#[test]
fn a4_class_featurexml_all_flags() {
    check_oracle(
        "file_info/inputs/AccurateMassSearchEngine_input1.featureXML",
        &mps(),
        "class_featurexml_mps",
    );
}

#[cfg(feature = "featurexml")]
#[test]
fn a4_class_featurexml_default_flags() {
    check_oracle(
        "file_info/inputs/AccurateMassSearchEngine_input1.featureXML",
        &Options::default(),
        "class_featurexml_default",
    );
}

#[cfg(feature = "featurexml")]
#[test]
fn a4_class_featurexml_statistics_only() {
    let options = Options {
        statistics: true,
        ..Options::default()
    };
    check_oracle(
        "file_info/inputs/AccurateMassSearchEngine_input1.featureXML",
        &options,
        "class_featurexml_s",
    );
}

#[cfg(feature = "featurexml")]
#[test]
fn a4_forced_featurexml_on_tmp_name() {
    let dir = openms::system::file::TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let copy = dir.path().join("test_forced_type.tmp");
    std::fs::copy(
        data("file_info/inputs/AccurateMassSearchEngine_input1.featureXML"),
        &copy,
    )
    .unwrap();
    let options = Options {
        forced_type: FileType::FeatureXml,
        ..mps()
    };
    let result = FileInfo::new().run(&copy, &options).unwrap();
    assert_report(&result.text, "file_info/expected/class_forced_type_tmp.txt");
    assert_report(&result.tsv, "file_info/expected/class_forced_type_tmp.tsv");
    assert_eq!(result.meta.file_type, FileType::FeatureXml);
    assert!(result.feature.is_some());
}

#[cfg(feature = "mzml")]
#[test]
fn a4_class_mzml_minimal_all_flags() {
    check_oracle(
        "file_info/inputs/MzMLFile_2_minimal.mzML",
        &mps(),
        "class_mzml_minimal_mps",
    );
}

#[cfg(feature = "featurexml")]
#[test]
fn a4_file_info_3_default_flags() {
    check_oracle(
        "file_info/inputs/FileInfo_3_input.featureXML",
        &Options::default(),
        "fileinfo_3_default",
    );
}

#[cfg(feature = "mzml")]
#[test]
#[ignore = "mzML reader gaps outside A4: duplicate userParam 'name', processing on primary arrays, float charge array"]
fn a4_file_info_9_default_flags() {
    check_oracle(
        "file_info/inputs/FileInfo_9_input.mzML",
        &Options::default(),
        "fileinfo_9_default",
    );
}

#[cfg(feature = "mzml")]
#[test]
fn a4_file_info_9_strict_reader_default_flags() {
    check_oracle(
        "file_info/inputs/FileInfo_9_strict_reader.mzML",
        &Options::default(),
        "fileinfo_9_default",
    );
    check_oracle(
        "file_info/inputs/FileInfo_9_strict_reader.mzML",
        &Options::default(),
        "fileinfo_9_strict_reader_default",
    );
}

#[cfg(feature = "mzml")]
#[test]
#[ignore = "mzML reader gap outside A4: 'charge array' stored as 64-bit float is refused \
            (canonical auxiliary array binary type); C++ converts it"]
fn a4_indexed_file_info_12_all_flags() {
    check_oracle(
        "file_info/inputs/FileInfo_12_input.mzML",
        &mps(),
        "fileinfo_12_mps",
    );
}

#[cfg(feature = "mzml")]
#[test]
fn a4_faims_test_data_all_flags() {
    let result = check_oracle(
        "mzml_mobility/FAIMS_test_data.mzML",
        &mps(),
        "faims_test_data_mps",
    );
    let peak = result.peak.as_ref().unwrap();
    assert_eq!(peak.faims_cvs, [-65.0]);
    assert!(result.warnings.is_empty());
}

#[cfg(feature = "mzml")]
#[test]
fn a4_faims_interleaved_all_flags() {
    check_oracle(
        "mzml_mobility/FAIMS_CV-60C_V-45_Interleaved.mzML",
        &mps(),
        "faims_interleaved_mps",
    );
}

#[cfg(feature = "mzml")]
#[test]
fn a4_im_faims_test_all_flags() {
    check_oracle(
        "faims_helper/IM_FAIMS_test.mzML",
        &mps(),
        "im_faims_test_mps",
    );
}

#[cfg(feature = "mzml")]
#[test]
fn a4_feature_finder_centroided_input_all_flags() {
    check_oracle(
        "mzml_mobility/FeatureFinderCentroided_1_input.mzML",
        &mps(),
        "ffc_input_mps",
    );
}

#[cfg(feature = "mzml")]
#[test]
fn a4_srm_chromatograms_all_flags() {
    check_oracle(
        "mzml_numpress_source_original.mzML",
        &mps(),
        "srm_chromatograms_mps",
    );
}

#[cfg(feature = "mzml")]
#[test]
#[ignore = "mzML reader gap outside A4 (A3 request 5): C++ copies the selected-ion drift time \
            8.1 onto the MS2 spectrum (MzMLHandler.cpp:1871-1875), so its ion-mobility ranges \
            end at 8.10; the Rust reader does not"]
fn a4_mzml_file_1_all_flags() {
    check_oracle("mzml_validator/MzMLFile_1.mzML", &mps(), "mzmlfile_1_mps");
}

#[cfg(feature = "mzml")]
#[test]
fn a4_mzml_file_1_without_selected_ion_drift_time_all_flags() {
    // The derived input deletes the selected-ion drift time; the scan-level
    // drift time 7.1 of the MS1 spectrum stays.
    let result = check_oracle(
        "file_info/inputs/MzMLFile_1_no_selected_ion_drift.mzML",
        &mps(),
        "mzmlfile_1_no_selected_ion_drift_mps",
    );
    let combined = result.ranges.combined.mobility.unwrap();
    assert_eq!((combined.min, combined.max), (7.1, 7.1));
    assert!(result.ranges.per_ms_level[&2].mobility.is_none());
}

#[cfg(feature = "mzml")]
#[test]
fn a4_precursor_purity_input_all_flags() {
    check_oracle(
        "precursor_purity_input.mzML",
        &mps(),
        "precursor_purity_mps",
    );
}

#[test]
fn a4_profile_dta_all_flags() {
    check_oracle(
        "spectrum_type/PeakTypeEstimator_raw.dta",
        &mps(),
        "pte_raw_dta_mps",
    );
}

#[test]
fn a4_centroid_dta_statistics_only() {
    let options = Options {
        statistics: true,
        ..Options::default()
    };
    check_oracle(
        "spectrum_type/PeakTypeEstimator_peak.dta",
        &options,
        "pte_peak_dta_s",
    );
}

#[test]
fn a4_dta2d_header_variants_all_flags() {
    for (index, case) in [(1, "dta2d_1_mps"), (2, "dta2d_2_mps"), (3, "dta2d_3_mps")] {
        check_oracle(
            &format!("text_peak_lists/DTA2DFile_test_{index}.dta2d"),
            &mps(),
            case,
        );
    }
}

// ---------------------------------------------------------------------------
// Tier 1: retained upstream outputs through FuzzyDiff, as registered
// (test-data 0cb15f2 topp/CMakeLists.txt:881-904)
// ---------------------------------------------------------------------------

fn fuzzy_diff_against_retained(result: &FileInfoResult, retained: &str) {
    let settings =
        fuzzy::FuzzyDiffSettings::load_ini(&data("fuzzy_string_comparator/FuzzyDiff.ini"))
            .unwrap()
            .with_whitelist(&["File name"]);
    assert_eq!(settings.ratio, 1.01);
    assert_eq!(settings.absdiff, 0.01);
    let expected = std::fs::read(data(retained)).unwrap();
    if let Err(log) = settings.compare_bytes(result.text.as_bytes(), &expected) {
        panic!("{retained}: FuzzyDiff failed:\n{log}");
    }
}

#[test]
fn retained_file_info_1_passes_fuzzy_diff() {
    // TOPP_FileInfo_1: -in FileInfo_1_input.dta -in_type dta -no_progress
    let result = FileInfo::new()
        .run(data("file_info/inputs/FileInfo_1_input.dta"), &dta_forced())
        .unwrap();
    fuzzy_diff_against_retained(&result, "file_info/retained/FileInfo_1_output.txt");
}

#[test]
fn retained_file_info_2_passes_fuzzy_diff() {
    // TOPP_FileInfo_2: -in FileInfo_2_input.dta2d -no_progress
    let result = FileInfo::new()
        .run_default(data("file_info/inputs/FileInfo_2_input.dta2d"))
        .unwrap();
    fuzzy_diff_against_retained(&result, "file_info/retained/FileInfo_2_output.txt");
}

#[cfg(feature = "featurexml")]
#[test]
fn retained_file_info_3_passes_fuzzy_diff() {
    // TOPP_FileInfo_3: -in FileInfo_3_input.featureXML -m -s -p -no_progress.
    // The retained file writes six spaces after `intensity:`; the source at the
    // pin writes one, which only FuzzyDiff's whitespace rule accepts.
    let result = FileInfo::new()
        .run(data("file_info/inputs/FileInfo_3_input.featureXML"), &mps())
        .unwrap();
    assert!(result.text.contains("  intensity: 1376.00 .. 6712.00\n"));
    let retained = read_text("file_info/retained/FileInfo_3_output.txt");
    assert!(retained.contains("  intensity:      1376.00 .. 6712.00\n"));
    fuzzy_diff_against_retained(&result, "file_info/retained/FileInfo_3_output.txt");
}

#[cfg(feature = "mzml")]
#[test]
fn retained_file_info_9_passes_fuzzy_diff() {
    // TOPP_FileInfo_9: -in FileInfo_9_input.mzML -m -p -s -no_progress. Run on
    // the derived single-name input, whose C++ report equals the original's
    // apart from the file name (see c1_file_info_9_strict_reader_mps).
    let result = FileInfo::new()
        .run(
            data("file_info/inputs/FileInfo_9_strict_reader.mzML"),
            &mps(),
        )
        .unwrap();
    fuzzy_diff_against_retained(&result, "file_info/retained/FileInfo_9_output.txt");
}

// ---------------------------------------------------------------------------
// Tier 3: FileInfo_test.cpp (core bc9cc12), section by section
// ---------------------------------------------------------------------------

#[test]
fn class_test_constructor_and_destructor() {
    // START_SECTION((FileInfo())) and ((~FileInfo())): construction and drop.
    let file_info = FileInfo::new();
    assert_eq!(file_info, FileInfo);
    let boxed = Box::new(file_info);
    drop(boxed);
}

#[cfg(feature = "featurexml")]
#[test]
fn class_test_run_featurexml() {
    // START_SECTION((Result run(...)) - featureXML)
    let result = FileInfo::new()
        .run_all(data(
            "file_info/inputs/AccurateMassSearchEngine_input1.featureXML",
        ))
        .unwrap();
    assert_eq!(result.meta.file_type, FileType::FeatureXml);
    assert_eq!(result.meta.file_type_name, "featureXML");
    assert!(result.feature.is_some());
    let feature = result.feature.as_ref().unwrap();
    assert!(!feature.is_consensus);
    assert!(feature.num_features > 0);
    assert!(!result.ranges.is_experiment);
    assert!(result.text.contains("-- General information --"));
    assert!(result.tsv.contains("general: number of features"));
}

#[cfg(feature = "consensusxml")]
#[test]
fn class_test_run_consensusxml_is_pending() {
    // START_SECTION((Result run ...) - consensusXML) expects a consensus
    // FeatureInfo. The consensusXML branch is not ported yet, so the section is
    // pending; until then the run refuses explicitly instead of returning a
    // result without it.
    let error = FileInfo::new()
        .run_all("AccurateMassSearchEngine_input1.consensusXML")
        .unwrap_err();
    assert!(
        matches!(&error, Error::Unsupported(message) if message.contains("consensusXML branch")),
        "{error}"
    );
}

#[cfg(feature = "mzml")]
#[test]
fn class_test_run_mzml_peaks() {
    // START_SECTION((Result run ...) - mzML peaks)
    let result = FileInfo::new()
        .run_all(data("file_info/inputs/MzMLFile_2_minimal.mzML"))
        .unwrap();
    assert_eq!(result.meta.file_type, FileType::MzMl);
    assert!(result.peak.is_some());
    assert!(result.ranges.is_experiment);
    let peak = result.peak.as_ref().unwrap();
    let sum: u64 = peak.spectra_per_ms_level.values().sum();
    assert_eq!(sum, peak.num_spectra);
}

#[test]
fn class_test_run_fasta_is_pending() {
    // START_SECTION((Result run ...) - FASTA) expects a FastaInfo; the FASTA
    // branch is not ported yet. The refusal comes before any file access.
    let error = FileInfo::new().run_all("FASTAFile_test.fasta").unwrap_err();
    assert!(
        matches!(&error, Error::Unsupported(message) if message.contains("fasta branch")),
        "{error}"
    );
}

#[cfg(feature = "featurexml")]
#[test]
fn class_test_to_text_and_to_tsv() {
    // START_SECTION((static std::string toText(const Result& r)) and toTSV)
    let result = FileInfo::new()
        .run_all(data(
            "file_info/inputs/AccurateMassSearchEngine_input1.featureXML",
        ))
        .unwrap();
    assert_eq!(FileInfo::to_text(&result), result.text);
    assert_eq!(FileInfo::to_tsv(&result), result.tsv);
    assert!(!result.text.is_empty());
    // The options overloads ignore their options, as in the source.
    let other = Options {
        statistics: false,
        ..Options::default()
    };
    assert_eq!(FileInfo::to_text_with_options(&result, &other), result.text);
    assert_eq!(FileInfo::to_tsv_with_options(&result, &other), result.tsv);
}

#[cfg(feature = "featurexml")]
#[test]
fn class_test_options_gating() {
    // START_SECTION((Options gating))
    let input = data("file_info/inputs/AccurateMassSearchEngine_input1.featureXML");
    let file_info = FileInfo::new();
    let default = file_info.run_default(&input).unwrap();
    assert!(!default.text.contains("-- Statistics --"));
    let options = Options {
        statistics: true,
        ..Options::default()
    };
    let stats = file_info.run(&input, &options).unwrap();
    assert!(stats.text.contains("-- Statistics --"));
}

#[cfg(feature = "featurexml")]
#[test]
fn class_test_forced_type_selects_the_parse_branch() {
    // START_SECTION((forced type selects the parse branch for an unrecognized extension))
    let dir = openms::system::file::TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let tmp_path = dir.path().join("test_forced_type.tmp");
    std::fs::copy(
        data("file_info/inputs/AccurateMassSearchEngine_input1.featureXML"),
        &tmp_path,
    )
    .unwrap();
    let options = Options {
        forced_type: FileType::FeatureXml,
        ..Options::default()
    };
    let result = FileInfo::new().run(&tmp_path, &options).unwrap();
    assert_eq!(result.meta.file_type, FileType::FeatureXml);
    assert!(result.feature.is_some());
}

// ---------------------------------------------------------------------------
// Structured result, checked against the executed reports
// ---------------------------------------------------------------------------

#[cfg(feature = "mzml")]
#[test]
fn peak_info_of_file_info_9_matches_its_report() {
    let result = FileInfo::new()
        .run(
            data("file_info/inputs/FileInfo_9_strict_reader.mzML"),
            &mps(),
        )
        .unwrap();
    let peak = result.peak.as_ref().unwrap();
    assert_eq!(peak.instrument_name, "LCQ Deca");
    assert_eq!(
        peak.mass_analyzers,
        [
            ("Quadrupole ion trap".to_owned(), 0.0),
            ("Linear ion trap".to_owned(), 0.0)
        ]
    );
    assert_eq!(peak.ms_levels, [1, 2]);
    assert_eq!(peak.total_peaks, 40);
    assert_eq!(peak.num_spectra, 4);
    assert_eq!(
        peak.spectra_per_ms_level.iter().collect::<Vec<_>>(),
        [(&1, &3), (&2, &1)]
    );
    assert_eq!(peak.peak_type_per_ms_level[&1], "Centroid (Centroid)");
    assert_eq!(peak.peak_type_per_ms_level[&2], "Centroid (Unknown)");
    assert_eq!(
        peak.activation_methods_flat(),
        [(2, "Collision-induced dissociation".to_owned(), 1)]
    );
    assert_eq!(
        peak.precursor_charges.iter().collect::<Vec<_>>(),
        [(&2, &1)]
    );
    let mut names: Vec<&str> = peak
        .float_arrays
        .keys()
        .chain(peak.int_arrays.keys())
        .chain(peak.string_arrays.keys())
        .map(String::as_str)
        .collect();
    names.sort_unstable();
    assert_eq!(names, ["charge array", "signal to noise array"]);
    assert!(peak.faims_cvs.is_empty());
    assert_eq!(peak.num_chromatograms, 0);
    assert_eq!(peak.num_chrom_peaks, 0);
    assert!(peak.chromatogram_types.is_empty());

    let ranges = &result.ranges;
    assert!(ranges.is_experiment);
    assert!(ranges.combined.has_mobility);
    assert!(ranges.combined.mobility.is_none());
    assert!(!ranges.chromatograms.has_mobility);
    assert!(ranges.chromatograms.rt.is_none());
    assert_eq!(ranges.per_ms_level.keys().collect::<Vec<_>>(), [&1, &2]);
    let level_2 = ranges.per_ms_level[&2];
    assert_eq!(level_2.rt.map(|r| (r.min, r.max)), Some((5.2, 5.2)));
    assert_eq!(level_2.intensity.map(|r| (r.min, r.max)), Some((2.0, 20.0)));

    assert_eq!(result.processing.len(), 1);
    let step = &result.processing[0];
    assert_eq!(step.software_name, "ProteoWizard software");
    assert_eq!(step.software_version, "1.0");
    assert_eq!(step.completion_time, "0000-00-00 00:00:00");
    assert_eq!(step.actions, ["Conversion to mzML format"]);
    // Declared by the source, never filled by its run.
    assert!(result.experiment_meta.is_none());
    assert!(result.statistics.is_empty());
    assert!(result.feature.is_none());
}

#[cfg(feature = "mzml")]
#[test]
fn processing_steps_keep_order_and_action_enum_order() {
    let result = FileInfo::new()
        .run(data("mzml_validator/MzMLFile_1.mzML"), &mps())
        .unwrap();
    let steps: Vec<(&str, &str, Vec<&str>)> = result
        .processing
        .iter()
        .map(|step| {
            (
                step.software_name.as_str(),
                step.completion_time.as_str(),
                step.actions.iter().map(String::as_str).collect(),
            )
        })
        .collect();
    assert_eq!(
        steps,
        [
            (
                "Xcalibur",
                "2001-02-03 04:05:00",
                vec!["Charge deconvolution", "Deisotoping"]
            ),
            (
                "ProteoWizard software",
                "2001-02-03 04:10:00",
                vec!["Conversion to mzML format"]
            ),
        ]
    );
    // Without -p the steps are not collected, as in the source.
    let plain = FileInfo::new()
        .run_default(data("mzml_validator/MzMLFile_1.mzML"))
        .unwrap();
    assert!(plain.processing.is_empty());
    let peak = plain.peak.as_ref().unwrap();
    assert_eq!(
        peak.chromatogram_types.iter().collect::<Vec<_>>(),
        [(&"total ion current chromatogram".to_owned(), &2)]
    );
    assert_eq!((peak.num_chromatograms, peak.num_chrom_peaks), (2, 25));
}

#[cfg(feature = "featurexml")]
#[test]
fn feature_info_of_file_info_3_matches_its_report() {
    let result = FileInfo::new()
        .run(data("file_info/inputs/FileInfo_3_input.featureXML"), &mps())
        .unwrap();
    let feature = result.feature.as_ref().unwrap();
    assert!(!feature.is_consensus);
    assert_eq!(feature.num_features, 10);
    assert_eq!(feature.tic, 36739.0);
    assert_eq!(
        feature.charges.iter().collect::<Vec<_>>(),
        [(&2, &4), (&3, &6)]
    );
    assert_eq!(
        feature.ids_per_element.iter().collect::<Vec<_>>(),
        [(&0, &7), (&1, &2), (&8, &1)]
    );
    assert_eq!((feature.assigned_ids, feature.unassigned_ids), (10, 2));
    assert!(feature.size_distribution.is_empty());
    assert!(feature.map_columns.is_empty());
    let ranges = &result.ranges;
    assert!(!ranges.is_experiment);
    assert!(!ranges.combined.has_mobility);
    assert!(ranges.combined.mobility.is_none());
    assert!(ranges.per_ms_level.is_empty());
    assert_eq!(
        ranges.combined.intensity.map(|r| (r.min, r.max)),
        Some((1376.0, 6712.0))
    );
    assert!(result.peak.is_none());
    assert_eq!(result.processing.len(), 1);
    assert_eq!(result.processing[0].actions, ["Data filtering"]);
}

// ---------------------------------------------------------------------------
// Tier 4: unknown, refused and failing inputs
// ---------------------------------------------------------------------------

#[test]
fn unknown_type_returns_an_empty_result() {
    let input = data("mzml_mobility/FileInfo_8_input.notype");
    let result = FileInfo::new().run(&input, &mps()).unwrap();
    assert_eq!(result.meta.file_type, FileType::Unknown);
    assert_eq!(result.meta.file_type_name, "unknown");
    assert_eq!(result.meta.file_name, input.to_str().unwrap());
    assert!(result.text.is_empty());
    assert!(result.tsv.is_empty());
    assert!(result.peak.is_none() && result.feature.is_none());
    assert_eq!(result.ranges, Default::default());
}

#[test]
fn unknown_type_wins_over_every_refused_flag() {
    // The source returns before report_ for an unknown type, so no flag is looked at.
    let options = Options {
        validate: true,
        check_index: true,
        detailed: true,
        check_corrupt: true,
        ..Options::default()
    };
    let result = FileInfo::new()
        .run(data("mzml_mobility/FileInfo_8_input.notype"), &options)
        .unwrap();
    assert!(result.text.is_empty());
}

fn unsupported_message(input: &str, options: &Options) -> String {
    match FileInfo::new().run(input, options) {
        Err(Error::Unsupported(message)) => message,
        other => panic!("{input}: expected Error::Unsupported, got {other:?}"),
    }
}

#[test]
fn validation_and_index_check_are_refused_for_every_type() {
    // The file names do not exist: the refusal precedes any file access.
    for input in ["missing.mzML", "missing.featureXML", "missing.dta"] {
        let validate = Options {
            validate: true,
            ..mps()
        };
        assert!(unsupported_message(input, &validate).contains("(-v)"));
        let index = Options {
            check_index: true,
            ..mps()
        };
        assert!(unsupported_message(input, &index).contains("(-i)"));
    }
}

#[test]
fn detailed_listing_and_corrupt_check_are_refused_for_peak_files() {
    for input in ["missing.mzML", "missing.dta", "missing.dta2d"] {
        let detailed = Options {
            detailed: true,
            ..Options::default()
        };
        assert!(unsupported_message(input, &detailed).contains("(-d)"));
        let corrupt = Options {
            check_corrupt: true,
            ..Options::default()
        };
        assert!(unsupported_message(input, &corrupt).contains("(-c)"));
    }
}

#[cfg(feature = "featurexml")]
#[test]
fn detailed_listing_and_corrupt_check_do_not_change_a_featurexml_report() {
    // FileInfo.cpp reads options.detailed and options.check_corrupt only in the
    // peak-file branch, so the featureXML report is the default one.
    let options = Options {
        detailed: true,
        check_corrupt: true,
        ..Options::default()
    };
    check_oracle(
        "file_info/inputs/AccurateMassSearchEngine_input1.featureXML",
        &options,
        "class_featurexml_default",
    );
}

#[test]
fn unported_branches_are_refused_by_name() {
    for (input, branch) in [
        ("missing.consensusXML", "consensusXML branch"),
        ("missing.idXML", "idXML branch"),
        ("missing.mzid", "mzid branch"),
        ("missing.fasta", "fasta branch"),
        ("missing.pepXML", "pepXML branch"),
        ("missing.mzTab", "mzTab branch"),
        ("missing.trafoXML", "trafoXML branch"),
        ("missing.pqp", "pqp branch"),
        ("missing.mzXML", "peak-file branch for mzXML"),
        ("missing.mzData", "peak-file branch for mzData"),
        ("missing.mgf", "peak-file branch for mgf"),
        ("missing.ms2", "peak-file branch for ms2"),
        ("missing.sqMass", "peak-file branch for sqMass"),
        ("missing.msp", "peak-file branch for msp"),
    ] {
        let message = unsupported_message(input, &Options::default());
        assert!(message.contains(branch), "{input}: {message}");
    }
    let forced = Options {
        forced_type: FileType::Xmass,
        ..Options::default()
    };
    assert!(unsupported_message("missing", &forced).contains("peak-file branch for fid"));
}

#[test]
fn a_type_the_source_cannot_load_is_a_parse_error() {
    let dir = openms::system::file::TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let text = dir.path().join("notes.txt");
    std::fs::write(&text, "plain text\n").unwrap();
    // FileHandler::loadExperiment's default branch throws ParseError.
    let error = FileInfo::new().run(&text, &Options::default()).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.ends_with("type is not supported for loading experiments")),
        "{error}"
    );
    // A forced type the detected type contradicts is refused as not allowed.
    let forced = Options {
        forced_type: FileType::Tsv,
        ..Options::default()
    };
    let error = FileInfo::new().run(&text, &forced).unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(message) if message.contains("not allowed")),
        "{error}"
    );
    // imzML is refused as an imaging format.
    let imzml = dir.path().join("image.imzML");
    std::fs::write(&imzml, "").unwrap();
    let error = FileInfo::new()
        .run(&imzml, &Options::default())
        .unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(message) if message.contains("imaging format")),
        "{error}"
    );
}

/// Tripwire for the ignored oracle cases: the three original inputs still fail
/// in the mzML reader, before FileInfo sees any data. When a reader change makes
/// one load, this test fails; remove that `#[ignore]` and its line here.
#[cfg(feature = "mzml")]
#[test]
fn reader_gaps_behind_the_ignored_cases_are_still_present() {
    for (input, expected) in [
        (
            "file_info/inputs/FileInfo_9_input.mzML",
            "duplicate userParam name name",
        ),
        (
            "file_info/inputs/empty.mzML",
            "unresolved dataProcessingRef",
        ),
        (
            "file_info/inputs/FileInfo_12_input.mzML",
            "canonical auxiliary array binary type",
        ),
    ] {
        let error = FileInfo::new().run(data(input), &mps()).unwrap_err();
        assert!(error.to_string().contains(expected), "{input}: {error}");
    }
    // MzMLFile_1: the selected-ion drift time is not copied onto the spectrum.
    let result = FileInfo::new()
        .run(data("mzml_validator/MzMLFile_1.mzML"), &mps())
        .unwrap();
    assert!(result.ranges.per_ms_level[&2].mobility.is_none());
}

#[cfg(feature = "mzml")]
#[test]
fn a_detected_type_other_than_the_forced_one_is_refused() {
    let options = Options {
        forced_type: FileType::Dta,
        ..Options::default()
    };
    let error = FileInfo::new()
        .run(data("file_info/inputs/FileInfo_9_input.mzML"), &options)
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
}

#[test]
fn missing_files_are_io_errors() {
    let dir = openms::system::file::TempDir::new_in(std::env::temp_dir(), false).unwrap();
    // An unknown extension needs the content, which cannot be read.
    let error = FileInfo::new()
        .run(dir.path().join("absent.notype"), &Options::default())
        .unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error}");
    // A known extension selects the branch; the loader cannot open the file.
    let error = FileInfo::new()
        .run(dir.path().join("absent.dta"), &Options::default())
        .unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error}");
}

#[cfg(feature = "featurexml")]
#[test]
fn a_truncated_featurexml_is_an_error_without_a_report() {
    let dir = openms::system::file::TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("corrupt.featureXML");
    let bytes = std::fs::read(data("file_info/inputs/FileInfo_3_input.featureXML")).unwrap();
    // The C1 oracle's derivation: the first 3000 bytes; C++ exits 3 (ParseError).
    std::fs::write(&path, &bytes[..3000]).unwrap();
    let error = FileInfo::new().run(&path, &mps()).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error}");
}

#[cfg(feature = "mzml")]
#[test]
fn a_truncated_mzml_is_an_error_without_a_report() {
    let dir = openms::system::file::TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("corrupt.mzML");
    let bytes = std::fs::read(data("file_info/inputs/FileInfo_9_input.mzML")).unwrap();
    std::fs::write(&path, &bytes[..3000]).unwrap();
    let error = FileInfo::new().run(&path, &mps()).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error}");
}

#[cfg(unix)]
#[test]
fn a_non_utf8_file_name_is_refused() {
    use std::os::unix::ffi::OsStrExt;
    let name = std::ffi::OsStr::from_bytes(b"bad\xff.dta");
    let error = FileInfo::new().run(name, &Options::default()).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
}

#[test]
fn file_name_is_printed_exactly_as_given() {
    let relative = "tests/data/file_info/inputs/FileInfo_2_input.dta2d";
    // Integration tests run in the package root, where this relative path resolves.
    let result = FileInfo::new().run_default(relative).unwrap();
    assert!(result.text.starts_with(&format!(
        "\n-- General information --\n\nFile name: {relative}\n"
    )));
    assert!(result.tsv.starts_with(&format!(
        "general: file name\t{relative}\ngeneral: file type\tdta2d\n"
    )));
    assert_eq!(result.meta.file_name, relative);
}

#[test]
fn options_and_result_defaults_follow_the_source_initialisers() {
    let options = Options::default();
    assert_eq!(options.forced_type, FileType::Unknown);
    assert!(!(options.meta || options.processing || options.statistics));
    assert!(
        !(options.detailed || options.check_corrupt || options.validate || options.check_index)
    );
    assert_eq!(
        options.log_type,
        openms::concept::progress_logger::ProgressLogType::None
    );
    let result = FileInfoResult::default();
    assert_eq!(result.meta.file_type, FileType::Unknown);
    assert!(result.validation.supported);
    assert!(!result.validation.performed && !result.validation.valid);
    assert_eq!(FileInfo::MAX_STATISTICS_VALUES, 1 << 27);
}
