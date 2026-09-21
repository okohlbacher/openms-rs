// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A8-FILEINFO, the cheap half: the FileInfo peak-file branch on mzXML, mzData,
//! MGF and MS2 (`FORMAT/FileInfo.cpp:1532-1965` and its `-m`, `-p`, `-s`, `-d`
//! and `-c` arms, core `bc9cc12`), loaded through the reader the source's
//! `FileHandler::loadExperiment` names for each type (`FileHandler.cpp:886-937`).
//!
//! Evidence (see `tests/data/file_info_a8_provenance.json` and
//! `docs/FILE_INFO_A8_SUPPORT.md`):
//!
//! - tier 1, executed differential: `../oracle/a8-fileinfo`, run twice on
//!   ibminode06 and reproduced, against the **Release** build
//!   `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`. The
//!   FileInfo tool ran 68 cases; `driver/fileinfo_driver.cpp` ran 14 more
//!   through `OpenMS::FileInfo::run` itself, because the tool's `-in` refuses
//!   the `ms2` extension and so never reaches `MS2File::load`. Of the 65 that
//!   wrote a report, 60 are compared here byte for byte, text and TSV, with
//!   only the input path normalised (the `File name: ` and `general: file name`
//!   lines, and the `-i` failure line). The other five: `x4_bare` wrote its
//!   report to standard output, which `tests/topp_file_info.rs` compares, and
//!   four are the refusals below;
//! - tier 1, retained upstream definition: TOPP_FileInfo_4, _5 and _6
//!   (`topp/CMakeLists.txt:890-898`, test-data `0cb15f2`), with the
//!   registrations' own flags (cases `x4`, `d5`, `d6`); `tests/topp_file_info.rs`
//!   runs the three through the tool against the retained outputs;
//! - tier 4, explicit refusal where the port declines what the Release build
//!   accepts: an mzData spectrum whose scan window the kernel's range validation
//!   refuses, and an MGF `MSLEVEL=0`, both recorded with the Release build's
//!   report in `tests/data/file_info_a8/expected`.
//!
//! mzXML and mzData need the `mzml` feature, as their readers do; MGF and MS2
//! need none.

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

/// A fixture this group added, under `tests/data/file_info_a8/inputs`.
fn a8(name: &str) -> String {
    format!("file_info_a8/inputs/{name}")
}

fn read_text(path: &Path) -> String {
    std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
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
            other => return format!("line {}: actual {line:?}, expected {other:?}", index + 1),
        }
    }
    match expected_lines.next() {
        Some(line) => format!("actual ends early; expected next {line:?}"),
        None => "no line differs".to_owned(),
    }
}

fn assert_report(actual: &str, expected_file: &Path, label: &str) {
    let expected = normalise_file_name(&read_text(expected_file));
    let actual = normalise_file_name(actual);
    assert!(
        actual == expected,
        "{label}: {}",
        first_difference(&actual, &expected)
    );
}

/// The options of an oracle case: the flags `m`, `p`, `s`, `d`, `c` and `i`,
/// and a forced type.
fn options(flags: &str, forced_type: FileType) -> Options {
    let mut options = Options {
        forced_type,
        ..Options::default()
    };
    for flag in flags.chars() {
        match flag {
            'm' => options.meta = true,
            'p' => options.processing = true,
            's' => options.statistics = true,
            'd' => options.detailed = true,
            'c' => options.check_corrupt = true,
            'i' => options.check_index = true,
            other => panic!("unknown flag {other}"),
        }
    }
    options
}

/// One oracle case: its id, the fixture below `tests/data`, the flags and the
/// forced type.
type Case = (&'static str, String, &'static str, FileType);

fn run(case: &Case) -> Result<FileInfoResult, Error> {
    let (_, input, flags, forced) = case;
    FileInfo::new().run(data(input), &options(flags, *forced))
}

/// The oracle's spelling of the input path, `<FIX>/<file name>`, in place of
/// the path this run read. The fixtures keep their oracle file names, so only
/// the directory differs; the `-i` report names the file in a second line.
fn oracle_path(report: &str, input: &Path) -> String {
    let name = input
        .file_name()
        .and_then(|name| name.to_str())
        .expect("name");
    report.replace(
        input.to_str().expect("UTF-8 path"),
        &format!("<FIX>/{name}"),
    )
}

/// Run `case` and compare both reports with the Release build's.
fn check(case: &Case) -> FileInfoResult {
    let id = case.0;
    let result = run(case).unwrap_or_else(|e| panic!("{id}: {e}"));
    let input = data(&case.1);
    assert_report(
        &oracle_path(&result.text, &input),
        &data(&format!("file_info_a8/expected/{id}.txt")),
        &format!("{id} text"),
    );
    assert_report(
        &oracle_path(&result.tsv, &input),
        &data(&format!("file_info_a8/expected/{id}.tsv")),
        &format!("{id} tsv"),
    );
    assert_eq!(FileInfo::to_text(&result), result.text);
    assert_eq!(FileInfo::to_tsv(&result), result.tsv);
    result
}

fn check_all(cases: &[Case]) {
    for case in cases {
        check(case);
    }
}

const ALL: &str = "mps";
const EVERYTHING: &str = "mpsdc";
const NONE: FileType = FileType::Unknown;

// ---------------------------------------------------------------------------
// mzXML (MzXMLFile with the handler's default PeakFileOptions)
// ---------------------------------------------------------------------------

#[cfg(feature = "mzml")]
fn mzxml_cases() -> Vec<Case> {
    let fi4 = "topp_file_info/inputs/FileInfo_4_input.mzXML".to_owned();
    vec![
        // TOPP_FileInfo_4's own flags.
        ("x4", fi4.clone(), "m", NONE),
        ("x4_all", fi4.clone(), EVERYTHING, NONE),
        // The class itself, on the same input (driver).
        ("lib_x4_all", fi4, EVERYTHING, NONE),
        ("x_mzxml1", "MzXMLFile_1.mzXML".into(), "", NONE),
        ("x_mzxml1_all", "MzXMLFile_1.mzXML".into(), ALL, NONE),
        ("x_mzxml1_dc", "MzXMLFile_1.mzXML".into(), "dc", NONE),
        (
            "x_mzxml1c_all",
            "MzXMLFile_1_compressed.mzXML".into(),
            EVERYTHING,
            NONE,
        ),
        (
            "x_mzxml2_all",
            "MzXMLFile_2_minimal.mzXML".into(),
            EVERYTHING,
            NONE,
        ),
        (
            "x_mzxml3_all",
            "MzXMLFile_3_64bit.mzXML".into(),
            EVERYTHING,
            NONE,
        ),
        (
            "x_inspect1_all",
            a8("InspectOutfile_test_1.mzXML"),
            EVERYTHING,
            NONE,
        ),
        (
            "x_fc4_all",
            a8("FileConverter_4_input.mzXML"),
            EVERYTHING,
            NONE,
        ),
        (
            "x_fc6_all",
            a8("FileConverter_6_output.mzXML"),
            EVERYTHING,
            NONE,
        ),
        (
            "x_fc26_all",
            a8("FileConverter_26_output.mzXML"),
            EVERYTHING,
            NONE,
        ),
        (
            "x_osw2_all",
            a8("OpenSwathWorkflow_2_input.mzXML"),
            EVERYTHING,
            NONE,
        ),
        (
            "x_spectrast_all",
            a8("spectra_spectrast.mzXML"),
            EVERYTHING,
            NONE,
        ),
        ("x_edges", a8("a8_mzxml_edges.mzXML"), "", NONE),
        ("x_edges_all", a8("a8_mzxml_edges.mzXML"), ALL, NONE),
        ("x_edges_dc", a8("a8_mzxml_edges.mzXML"), "dc", NONE),
        (
            "x_no_scans_all",
            a8("a8_mzxml_no_scans.mzXML"),
            EVERYTHING,
            NONE,
        ),
    ]
}

/// Every mzXML case of the oracle that wrote a report, byte for byte.
#[cfg(feature = "mzml")]
#[test]
fn mzxml_reports_match_the_release_build() {
    check_all(&mzxml_cases());
}

/// TOPP_FileInfo_4 (`topp/CMakeLists.txt:890-892`): the structured result the
/// report is written from, with the numbers of the Release build's report.
#[cfg(feature = "mzml")]
#[test]
fn mzxml_upstream_4_structured_result() {
    let cases = mzxml_cases();
    let result = check(&cases[0]);
    assert_eq!(result.meta.file_type, FileType::MzXml);
    let peak = result.peak.expect("the peak branch fills peak");
    assert_eq!(peak.num_spectra, 20);
    assert_eq!(peak.total_peaks, 6864);
    assert_eq!(peak.ms_levels, vec![1, 2]);
    assert_eq!(peak.spectra_per_ms_level[&1], 14);
    assert_eq!(peak.spectra_per_ms_level[&2], 6);
    assert_eq!(peak.precursor_charges[&0], 6);
    assert_eq!(
        peak.mass_analyzers,
        vec![("Quadrupole ion trap".to_owned(), 0.0)]
    );
}

/// The generated mzXML: `activationMethod` short names are looked up, an
/// unknown one (`FOO`) is dropped, a scan with `msLevel="0"` is read as MS1
/// with a warning (`MzXMLHandler.cpp:247-252`), and every precursor of a scan
/// counts its activation methods while only the first counts its charge.
#[cfg(feature = "mzml")]
#[test]
fn mzxml_activation_and_ms_level_zero() {
    let cases = mzxml_cases();
    let edges = cases
        .iter()
        .find(|case| case.0 == "x_edges")
        .expect("x_edges");
    let peak = check(edges).peak.expect("peak");
    assert_eq!(peak.ms_levels, vec![1, 2, 3]);
    assert_eq!(peak.spectra_per_ms_level[&1], 3);
    let methods: Vec<_> = peak.activation_methods.keys().cloned().collect();
    assert_eq!(
        methods,
        vec![
            (2, "Collision-induced dissociation".to_owned()),
            (2, "Electron transfer dissociation".to_owned()),
            (2, "beam-type collision-induced dissociation".to_owned()),
        ]
    );
    assert_eq!(peak.precursor_charges.get(&-1), Some(&1));
    assert_eq!(peak.precursor_charges.get(&0), Some(&1));
    assert_eq!(peak.precursor_charges.get(&2), Some(&2));
}

/// `-i` on mzXML through the class: the index check reads the file as indexed
/// mzML, finds no index and ends the report there, without its trailing
/// newlines. The tool refuses `-i` on anything but mzML before the class runs
/// (oracle case `x4_i`); the class does not (driver case `lib_x4_i`).
#[cfg(feature = "mzml")]
#[test]
fn mzxml_index_check_through_the_class() {
    let case: Case = (
        "lib_x4_i",
        "topp_file_info/inputs/FileInfo_4_input.mzXML".into(),
        "i",
        NONE,
    );
    let result = check(&case);
    assert!(result.validation.index_checked);
    assert!(!result.validation.index_valid);
    assert!(result.peak.is_none());
}

// ---------------------------------------------------------------------------
// mzData (MzDataFile with the handler's default PeakFileOptions)
// ---------------------------------------------------------------------------

#[cfg(feature = "mzml")]
fn mzdata_cases() -> Vec<Case> {
    let fi5 = "mzml_mobility/FileInfo_5_input.mzDat".to_owned();
    let fi6 = "topp_file_info/inputs/FileInfo_6_input.mzData".to_owned();
    vec![
        // TOPP_FileInfo_5's own flags, `-in_type mzData` included.
        ("d5", fi5.clone(), "ms", FileType::MzData),
        ("d5_all", fi5.clone(), EVERYTHING, FileType::MzData),
        ("lib_d5_all", fi5.clone(), EVERYTHING, FileType::MzData),
        // `.mzDat` is no mzData extension: the type is found by content.
        ("d5_detected", fi5, "", NONE),
        // TOPP_FileInfo_6's own flags.
        ("d6", fi6.clone(), "ds", NONE),
        ("d6_all", fi6, EVERYTHING, NONE),
        (
            "d_mzdata3_all",
            "MzDataFile_3_minimal.mzData".into(),
            EVERYTHING,
            NONE,
        ),
        (
            "d_mzdata4_all",
            "MzDataFile_4_64bit.mzData".into(),
            EVERYTHING,
            NONE,
        ),
        (
            "d_fc8_all",
            a8("FileConverter_8_output.mzData"),
            EVERYTHING,
            NONE,
        ),
        ("d_pepnovo_all", a8("PepNovo.mzData"), EVERYTHING, NONE),
        (
            "d_picked_all",
            a8("PickedPeakTestData.mzData"),
            EVERYTHING,
            NONE,
        ),
        (
            "d_extender_all",
            a8("SimpleExtender_test.mzData"),
            EVERYTHING,
            NONE,
        ),
        ("d_edges", a8("a8_mzdata_edges.mzData"), "", NONE),
        ("d_edges_all", a8("a8_mzdata_edges.mzData"), ALL, NONE),
        ("d_edges_dc", a8("a8_mzdata_edges.mzData"), "dc", NONE),
        (
            "d_no_spectra_all",
            a8("a8_mzdata_no_spectra.mzData"),
            EVERYTHING,
            NONE,
        ),
    ]
}

/// Every mzData case of the oracle that wrote a report and that this port
/// runs, byte for byte; `MzDataFile_1.mzData` is the refusal below.
#[cfg(feature = "mzml")]
#[test]
fn mzdata_reports_match_the_release_build() {
    check_all(&mzdata_cases());
}

/// TOPP_FileInfo_5 and _6 (`topp/CMakeLists.txt:893-898`): the structured
/// result, with the numbers of the Release build's reports.
#[cfg(feature = "mzml")]
#[test]
fn mzdata_upstream_5_and_6_structured_result() {
    let cases = mzdata_cases();
    let five = check(&cases[0]);
    assert_eq!(five.meta.file_type, FileType::MzData);
    let peak = five.peak.expect("peak");
    assert_eq!(peak.num_spectra, 10);
    assert_eq!(peak.total_peaks, 3149);
    assert_eq!(peak.instrument_name, "MS-Instrument");
    assert_eq!(
        peak.mass_analyzers,
        vec![
            ("Quadrupole ion trap".to_owned(), 22.33),
            ("Quadrupole".to_owned(), 12.3)
        ]
    );
    assert!(peak.precursor_charges.is_empty());

    let six = cases.iter().find(|case| case.0 == "d6").expect("d6");
    let peak = check(six).peak.expect("peak");
    assert_eq!(peak.num_spectra, 2);
    assert_eq!(peak.total_peaks, 9);
    assert_eq!(peak.float_arrays.len(), 7);
    assert_eq!(peak.float_arrays["SignalToNoise"], 2);
}

/// The generated mzData: a doubled `ChargeState` resets the charge to zero
/// (`MzDataHandler.cpp:1194-1205`), an activation method outside the handler's
/// vocabulary (`ETD`) is read as the first entry, CID, and an invalid
/// `spectrumType` leaves the type unknown, which peak picking in the data
/// processing then reports as centroid.
#[cfg(feature = "mzml")]
#[test]
fn mzdata_charge_activation_and_spectrum_type() {
    let case: Case = ("d_edges", a8("a8_mzdata_edges.mzData"), "", NONE);
    let peak = check(&case).peak.expect("peak");
    assert_eq!(peak.precursor_charges.get(&0), Some(&1));
    assert_eq!(peak.precursor_charges.get(&1), Some(&1));
    assert_eq!(
        peak.activation_methods
            .get(&(2, "Collision-induced dissociation".to_owned())),
        Some(&2)
    );
    assert_eq!(peak.peak_type_per_ms_level[&1], "Profile (Profile)");
    assert_eq!(peak.peak_type_per_ms_level[&2], "Centroid (Unknown)");
}

/// Refused, and recorded: `MzDataFile_1.mzData` holds a spectrum with
/// `mzRangeStart="110"` and no `mzRangeStop`. The mzData reader keeps the
/// source's scan window `[110, 0]` (`MzDataHandler.cpp:392-396`), and the
/// kernel's `MSSpectrum::range_manager` validates the whole spectrum, scan
/// windows included, before it reads the peaks, so the range computation
/// refuses it. The Release build, whose `updateRanges` looks at the peaks only,
/// prints the three reports kept as `expected/d_mzdata1*.txt` and `.tsv`; the
/// refusal belongs to the kernel, not to this branch or the reader.
#[cfg(feature = "mzml")]
#[test]
fn mzdata_inverted_scan_window_is_refused_by_the_kernel_range_validation() {
    for (id, flags) in [
        ("d_mzdata1", ""),
        ("d_mzdata1_all", ALL),
        ("d_mzdata1_dc", "dc"),
    ] {
        let case: Case = (id, "MzDataFile_1.mzData".into(), flags, NONE);
        match run(&case) {
            Err(Error::InvalidValue(message)) => {
                assert_eq!(message, "scan window begin exceeds end", "{id}");
            }
            other => panic!("{id}: {other:?}"),
        }
        assert!(data(&format!("file_info_a8/expected/{id}.txt")).is_file());
    }
}

// ---------------------------------------------------------------------------
// MGF (MascotGenericFile::load, with the source's carry-over)
// ---------------------------------------------------------------------------

fn mgf_cases() -> Vec<Case> {
    let gnps = "MascotGenericFile_GNPS.mgf".to_owned();
    vec![
        ("g_gnps", gnps.clone(), "", NONE),
        ("g_gnps_all", gnps.clone(), ALL, NONE),
        ("g_gnps_dc", gnps, "dc", NONE),
        (
            "g_content1_all",
            a8("FileHandler_MGFbyContent1.mgf"),
            EVERYTHING,
            NONE,
        ),
        (
            "g_content2_all",
            a8("FileHandler_MGFbyContent2.mgf"),
            EVERYTHING,
            NONE,
        ),
        (
            "g_idfc32_all",
            a8("IDFileConverter_32_output.mgf"),
            EVERYTHING,
            NONE,
        ),
        (
            "g_gnpsexport1_all",
            a8("GNPSExport_1_out.mgf"),
            EVERYTHING,
            NONE,
        ),
        ("g_carry", a8("a8_mgf_carry.mgf"), "", NONE),
        ("g_carry_all", a8("a8_mgf_carry.mgf"), ALL, NONE),
        ("g_carry_dc", a8("a8_mgf_carry.mgf"), "dc", NONE),
        ("lib_g_carry_all", a8("a8_mgf_carry.mgf"), ALL, NONE),
        (
            "g_title_min_all",
            a8("a8_mgf_title_min.mgf"),
            EVERYTHING,
            NONE,
        ),
        ("g_merged_all", a8("a8_mgf_merged.mgf"), EVERYTHING, NONE),
        ("g_columns_all", a8("a8_mgf_columns.mgf"), EVERYTHING, NONE),
        ("g_no_ions_all", a8("a8_mgf_no_ions.mgf"), EVERYTHING, NONE),
        ("g_mslevel_all", a8("a8_mgf_mslevel.mgf"), EVERYTHING, NONE),
        (
            "g_mslevel_negative_all",
            a8("a8_mgf_mslevel_negative.mgf"),
            EVERYTHING,
            NONE,
        ),
    ]
}

/// Every MGF case of the oracle that wrote a report and that this port runs,
/// byte for byte; `MSLEVEL=0` is the refusal below.
#[test]
fn mgf_reports_match_the_release_build() {
    check_all(&mgf_cases());
}

/// The source reader keeps one spectrum across blocks
/// (`MascotGenericFile.h:89-104`, `:141-157`): a block without `CHARGE=`,
/// `PEPMASS=`, `RTINSECONDS=` or `MSLEVEL=` inherits the previous block's, so
/// the four blocks count charges 2, 2, 3, 3 and MS levels 2, 2, 1, 1. Read
/// with the native default, a fresh record per block, the charges would be
/// 2, 0, 3, 0.
#[test]
fn mgf_blocks_inherit_the_previous_blocks_values() {
    let case: Case = ("g_carry", a8("a8_mgf_carry.mgf"), "", NONE);
    let peak = check(&case).peak.expect("peak");
    assert_eq!(peak.precursor_charges.get(&2), Some(&2));
    assert_eq!(peak.precursor_charges.get(&3), Some(&2));
    assert_eq!(peak.precursor_charges.get(&0), None);
    assert_eq!(peak.spectra_per_ms_level[&1], 2);
    assert_eq!(peak.spectra_per_ms_level[&2], 2);
    // MGF is centroided by definition (MascotGenericFile.h:95).
    assert_eq!(peak.peak_type_per_ms_level[&2], "Centroid (Centroid)");
}

/// `MSLEVEL=-1`: `std::stoi` gives `-1` and `setMSLevel` stores it in a
/// `UInt`, so the report prints MS level `4294967295`, and the result's `Int`
/// key is `static_cast<Int>` of that, `-1` (`FORMAT/FileInfo.cpp:1645-1648`).
#[test]
fn mgf_negative_ms_level_wraps_as_in_the_source() {
    let case: Case = (
        "g_mslevel_negative_all",
        a8("a8_mgf_mslevel_negative.mgf"),
        EVERYTHING,
        NONE,
    );
    let result = check(&case);
    assert!(result.text.contains("MS levels: 4294967295\n"));
    let peak = result.peak.expect("peak");
    assert_eq!(peak.ms_levels, vec![-1]);
    assert_eq!(peak.spectra_per_ms_level[&-1], 1);
}

/// Refused, and recorded: `MSLEVEL=0` loads as MS level 0, as in the source,
/// but the kernel's range computation validates each spectrum and refuses an
/// MS level of 0 outside the three optical scan modes. The Release build
/// prints the report kept as `expected/g_mslevel_zero_all.txt`, its `-c` block
/// flagging the level as an error.
#[test]
fn mgf_ms_level_zero_is_refused_by_the_kernel_range_validation() {
    let case: Case = (
        "g_mslevel_zero_all",
        a8("a8_mgf_mslevel_zero.mgf"),
        EVERYTHING,
        NONE,
    );
    match run(&case) {
        Err(Error::InvalidValue(message)) => {
            assert_eq!(message, "spectrum MS level must be positive");
        }
        other => panic!("{other:?}"),
    }
    assert!(data("file_info_a8/expected/g_mslevel_zero_all.txt").is_file());
}

// ---------------------------------------------------------------------------
// MS2 (MS2File::load), through the class only
// ---------------------------------------------------------------------------

fn ms2_cases() -> Vec<Case> {
    let test = "text_peak_lists/MS2File_test_spectra.ms2".to_owned();
    vec![
        ("lib_m_test", test.clone(), "", NONE),
        ("lib_m_test_all", test.clone(), ALL, NONE),
        ("lib_m_test_dc", test, "dc", NONE),
        ("lib_m_edges", a8("a8_ms2_edges.ms2"), "", NONE),
        ("lib_m_edges_all", a8("a8_ms2_edges.ms2"), ALL, NONE),
        ("lib_m_edges_dc", a8("a8_ms2_edges.ms2"), "dc", NONE),
        (
            "lib_m_no_scans_all",
            a8("a8_ms2_no_scans.ms2"),
            EVERYTHING,
            NONE,
        ),
    ]
}

/// Every MS2 case of the driver, byte for byte. The FileInfo tool's `-in`
/// refuses the `ms2` extension before the class runs, in the Release build as
/// here (oracle cases `m_*`, exit 6), so these reports come from
/// `OpenMS::FileInfo::run` directly.
#[test]
fn ms2_reports_match_the_release_build_class() {
    check_all(&ms2_cases());
}

/// `MS2File::load` drops a peak line before the first `S` line, keeps an `S`
/// line without peaks as an empty spectrum, never reads a charge (`Z` lines
/// are skipped) and never sets a retention time, which stays at the default
/// `-1`.
#[test]
fn ms2_scans_carry_no_charge_and_no_retention_time() {
    let case: Case = ("lib_m_edges", a8("a8_ms2_edges.ms2"), "", NONE);
    let result = check(&case);
    let peak = result.peak.expect("peak");
    assert_eq!(peak.num_spectra, 3);
    assert_eq!(peak.total_peaks, 14);
    assert_eq!(peak.precursor_charges.get(&0), Some(&3));
    let rt = result.ranges.spectra_overall.rt.expect("rt range");
    assert_eq!((rt.min, rt.max), (-1.0, -1.0));
}

// ---------------------------------------------------------------------------
// What the Release build refuses too
// ---------------------------------------------------------------------------

/// Malformed MGF and MS2 input: each reader refuses it with a parse error, as
/// the source's does (oracle `g_truncated` and `g_one_value` exit 3 with
/// `ParseError`; driver `lib_m_bad_s` and `lib_m_bad_peak` exit 3). `CHARGE=2-`
/// is a parse error here too; the source's `StringUtils::toInt32` throws a
/// `ConversionError` there, which the tool reports as an unexpected internal
/// error, exit 8 (oracle `g_charge_minus`).
#[test]
fn malformed_mgf_and_ms2_are_parse_errors() {
    for (id, input) in [
        ("g_truncated", a8("a8_mgf_truncated.mgf")),
        ("g_one_value", a8("a8_mgf_one_value.mgf")),
        ("g_charge_minus", a8("a8_mgf_charge_minus.mgf")),
        ("lib_m_bad_s", a8("a8_ms2_bad_s.ms2")),
        ("lib_m_bad_peak", a8("a8_ms2_bad_peak.ms2")),
    ] {
        let case: Case = (id, input, "", NONE);
        match run(&case) {
            Err(Error::Parse { .. }) => {}
            other => panic!("{id}: {other:?}"),
        }
    }
}

/// A forced type the loader's own detection contradicts: the source throws
/// `ParseError` (`FileHandler.cpp:858-864`; oracle `g_forced_on_mzxml` and
/// driver `lib_m_forced_mgf`, exit 3). This port refuses with
/// [`Error::InvalidValue`], as for every other peak type
/// (`docs/TOPP_FILE_INFO_SUPPORT.md`, native difference 3).
#[test]
fn a_forced_type_the_loader_detects_otherwise_is_refused() {
    for (id, input) in [
        (
            "g_forced_on_mzxml",
            "topp_file_info/inputs/FileInfo_4_input.mzXML".to_owned(),
        ),
        (
            "lib_m_forced_mgf",
            "text_peak_lists/MS2File_test_spectra.ms2".to_owned(),
        ),
    ] {
        let case: Case = (id, input, "", FileType::Mgf);
        match run(&case) {
            Err(Error::InvalidValue(message)) => {
                assert!(
                    message.ends_with("is not an allowed input format"),
                    "{id}: {message}"
                );
            }
            other => panic!("{id}: {other:?}"),
        }
    }
}

// ---------------------------------------------------------------------------
// What stays refused: the peak types without a native reader
// ---------------------------------------------------------------------------

/// sqMass, XMass (`fid`) and MSP have no reader here that fills an
/// `MSExperiment`, so their branch is still refused before the file is read,
/// with the branch named. The four types this group wired are no longer.
#[test]
fn peak_types_without_a_reader_stay_refused() {
    let input = data("text_peak_lists/MS2File_test_spectra.ms2");
    for (forced, name) in [
        (FileType::SqMass, "sqMass"),
        (FileType::Xmass, "fid"),
        (FileType::Msp, "msp"),
    ] {
        match FileInfo::new().run(&input, &options("", forced)) {
            Err(Error::Unsupported(message)) => assert_eq!(
                message,
                format!("FileInfo peak-file branch for {name} input is not ported")
            ),
            other => panic!("{name}: {other:?}"),
        }
    }
    for forced in [
        FileType::MzXml,
        FileType::MzData,
        FileType::Mgf,
        FileType::Ms2,
    ] {
        let outcome = FileInfo::new().run(&input, &options("", forced));
        assert!(
            !matches!(outcome, Err(Error::Unsupported(ref message)) if message.contains("is not ported")),
            "{forced:?}: {outcome:?}"
        );
    }
}
