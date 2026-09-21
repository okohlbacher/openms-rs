// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A8-FILEINFO, the cheap half: the FileInfo peak-file branch on mzXML, mzData,
//! MGF and MS2 (`FORMAT/FileInfo.cpp:1532-1965` and its `-m`, `-p`, `-s`, `-d`
//! and `-c` arms, core `bc9cc12`), loaded through the reader the source's
//! `FileHandler::loadExperiment` names for each type (`FORMAT/FileHandler.cpp:886-937`).
//!
//! Evidence (see `tests/data/file_info_a8_provenance.json` and
//! `docs/FILE_INFO_A8_SUPPORT.md`):
//!
//! - tier 1, executed differential: `../oracle/a8-fileinfo`, run twice on
//!   ibminode06 and reproduced, against the **Release** build
//!   `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`. The
//!   FileInfo tool ran 70 cases; `driver/fileinfo_driver.cpp` ran 14 more
//!   through `OpenMS::FileInfo::run` itself, because the tool's `-in` refuses
//!   the `ms2` extension and so never reaches `MS2File::load`. Of the 67 that
//!   wrote a report, 61 are compared here byte for byte, text and TSV, with
//!   only the input path normalised (the `File name: ` and `general: file name`
//!   lines, and the `-i` failure line). The other six: `x4_bare` wrote its
//!   report to standard output, which `tests/topp_file_info.rs` compares, and
//!   five are the refusals below;
//! - tier 1, retained upstream definition: TOPP_FileInfo_4, _5 and _6
//!   (`topp/CMakeLists.txt:890-898`, test-data `0cb15f2`), with the
//!   registrations' own flags (cases `x4`, `d5`, `d6`); `tests/topp_file_info.rs`
//!   runs the three through the tool against the retained outputs;
//! - tier 4, explicit refusal where the port declines what the Release build
//!   accepts: an mzData spectrum whose scan window the kernel's range validation
//!   refuses, an MGF `MSLEVEL=0`, and a gzip-compressed MGF, each recorded with
//!   the Release build's report in `tests/data/file_info_a8/expected`;
//! - tier 1, executed differential, truncated input (the repair of verifier
//!   finding F2): `../oracle/a8-truncated`, 42 tool cases, 13 class-driver
//!   cases and 8 `MzXMLFile::load` cases on the same Release install, twice,
//!   reproduced, every one checked in the last section against
//!   `tests/data/file_info_a8/truncated`.
//!
//! mzXML and mzData need the `mzml` feature, as their readers do; MGF and MS2
//! need none.

#[cfg(all(feature = "mzml", feature = "paramxml", feature = "featurexml"))]
#[path = "support/took_line.rs"]
mod took_line;

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
    check_in(case, "file_info_a8/expected")
}

/// As [`check`], with the Release reports in `expected`, below `tests/data`.
fn check_in(case: &Case, expected: &str) -> FileInfoResult {
    let id = case.0;
    let result = run(case).unwrap_or_else(|e| panic!("{id}: {e}"));
    let input = data(&case.1);
    assert_report(
        &oracle_path(&result.text, &input),
        &data(&format!("{expected}/{id}.txt")),
        &format!("{id} text"),
    );
    assert_report(
        &oracle_path(&result.tsv, &input),
        &data(&format!("{expected}/{id}.tsv")),
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
        // Compressed: the XML readers decompress, in the source as here.
        ("x_edges_gz", a8("a8_mzxml_edges_gz.mzXML.gz"), "", NONE),
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
/// with a warning (`FORMAT/HANDLERS/MzXMLHandler.cpp:247-252`), and every precursor of a scan
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
/// (`FORMAT/HANDLERS/MzDataHandler.cpp:1194-1204`), an activation method outside the handler's
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
/// source's scan window `[110, 0]` (`FORMAT/HANDLERS/MzDataHandler.cpp:392-395`), and the
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
/// (`FORMAT/MascotGenericFile.h:89-104`, `:141-157`): a block without `CHARGE=`,
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
    // MGF is centroided by definition (FORMAT/MascotGenericFile.h:95).
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

/// `mascot_generic::ReadOptions::source_ms_level`, the switch this group added
/// to the MGF reader: off, the default, `MSLEVEL=0` and `MSLEVEL=-1` are parse
/// errors; on, each is stored as the source's `setMSLevel(std::stoi(...))`
/// stores it (`FORMAT/MascotGenericFile.h:325-331`), `0` and `4294967295`, which are
/// the MS levels the Release build's reports print.
#[test]
fn mgf_source_ms_level_stores_what_std_stoi_returns() {
    use openms::format::mascot_generic::{ReadOptions, read_with_options};
    for (line, level, report) in [
        ("MSLEVEL=0", 0, "g_mslevel_zero_all"),
        ("MSLEVEL=-1", u32::MAX, "g_mslevel_negative_all"),
    ] {
        let block = format!("BEGIN IONS\n{line}\nPEPMASS=300\n100 1\nEND IONS\n");
        let strict = read_with_options(block.as_bytes(), &ReadOptions::default());
        assert!(
            matches!(strict, Err(Error::Parse { .. })),
            "{line}: {strict:?}"
        );
        let source = ReadOptions {
            source_ms_level: true,
            ..ReadOptions::default()
        };
        let experiment = read_with_options(block.as_bytes(), &source).expect(line);
        assert_eq!(experiment.spectra[0].ms_level, level, "{line}");
        let expected = read_text(&data(&format!("file_info_a8/expected/{report}.txt")));
        assert!(
            expected.contains(&format!("\nMS levels: {level}\n")),
            "{report}"
        );
    }
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

/// Refused, and recorded: a gzip-compressed MGF. `FileHandler::getType` strips
/// the `.gz` and answers MGF in both builds, and `MascotGenericFile::load`
/// reads with a plain `std::ifstream`, so the source sees the compressed bytes
/// as text, finds no `BEGIN IONS` line and reports an empty map
/// (`expected/g_gzipped_all.txt`: 0 spectra, exit 0). This reader does not
/// decompress either, and refuses bytes that are not text instead of reporting
/// a file it did not read. The XML readers decompress in both builds
/// (`x_edges_gz`, compared above).
#[test]
fn a_compressed_mgf_is_refused_where_the_source_reports_an_empty_map() {
    let case: Case = (
        "g_gzipped_all",
        a8("a8_mgf_gzipped.mgf.gz"),
        EVERYTHING,
        NONE,
    );
    match run(&case) {
        Err(Error::Io(error)) => {
            assert_eq!(error.kind(), std::io::ErrorKind::InvalidData, "{error}");
        }
        other => panic!("{other:?}"),
    }
    let expected = read_text(&data("file_info_a8/expected/g_gzipped_all.txt"));
    assert!(expected.contains("\nNumber of spectra: 0\n"), "{expected}");
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
/// `ParseError` (`FORMAT/FileHandler.cpp:858-864`; oracle `g_forced_on_mzxml` and
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

// ---------------------------------------------------------------------------
// Truncated input: the repair of verifier finding F2
// (../oracle/a8-truncated, tests/data/file_info_a8/truncated)
// ---------------------------------------------------------------------------

/// The truncation oracle's evidence, below `tests/data`: the fixtures under
/// `inputs`, byte prefixes of five A8 inputs, and under `expected` the
/// Release build's reports of every case that completed and one row per case
/// in `release_outcomes.tsv`.
const TRUNCATED: &str = "file_info_a8/truncated";
const TRUNCATED_EXPECTED: &str = "file_info_a8/truncated/expected";

/// What the Release build did on one case of the truncation oracle.
// Read only by the `mzml`-gated truncation test; unused on a build without it.
#[cfg_attr(not(feature = "mzml"), allow(dead_code))]
struct ReleaseOutcome {
    /// `tool` (the FileInfo executable), `driver` (`OpenMS::FileInfo::run`)
    /// or `mzxml_driver` (`MzXMLFile::load`).
    program: String,
    exit_code: i32,
    stdout_bytes: usize,
    /// The files the case left in its working directory, with their sizes.
    cwd_files: Vec<(String, usize)>,
    /// The last line on stderr, or on stdout for a load `mzxml_driver`
    /// completed, with the fixture directory spelled `<FIX>`.
    last_line: String,
}

/// Every row of `release_outcomes.tsv`, in the oracle's order.
fn release_outcomes() -> Vec<(String, ReleaseOutcome)> {
    let text = read_text(&data(&format!("{TRUNCATED_EXPECTED}/release_outcomes.tsv")));
    let mut lines = text.lines();
    assert_eq!(
        lines.next(),
        Some("id\tprogram\texit_code\tstdout_bytes\tcwd_files\tlast_line")
    );
    lines
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(fields.len(), 6, "{line}");
            let cwd_files = if fields[4] == "-" {
                Vec::new()
            } else {
                fields[4]
                    .split(';')
                    .map(|entry| {
                        let (name, bytes) = entry.rsplit_once(':').expect(entry);
                        (name.to_owned(), bytes.parse().expect(entry))
                    })
                    .collect()
            };
            let outcome = ReleaseOutcome {
                program: fields[1].to_owned(),
                exit_code: fields[2].parse().expect(line),
                stdout_bytes: fields[3].parse().expect(line),
                cwd_files,
                last_line: fields[5].to_owned(),
            };
            (fields[0].to_owned(), outcome)
        })
        .collect()
}

fn release_outcome(id: &str) -> ReleaseOutcome {
    release_outcomes()
        .into_iter()
        .find(|(case, _)| case == id)
        .map(|(_, outcome)| outcome)
        .unwrap_or_else(|| panic!("no Release outcome for {id}"))
}

/// The fixture a case ran on: the oracle names each case after its fixture's
/// stem, with `_all` for the run with every flag and `lib_` for the class
/// driver.
fn truncated_fixture(id: &str) -> String {
    let stem = id.strip_prefix("lib_").unwrap_or(id);
    let stem = stem.strip_suffix("_all").unwrap_or(stem);
    let inputs = data(&format!("{TRUNCATED}/inputs"));
    let mut names = std::fs::read_dir(&inputs)
        .unwrap_or_else(|e| panic!("{}: {e}", inputs.display()))
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .filter(|name| name.split_once('.').map(|(s, _)| s) == Some(stem));
    let name = names
        .next()
        .unwrap_or_else(|| panic!("no fixture for {id}"));
    assert!(names.next().is_none(), "{id}");
    format!("{TRUNCATED}/inputs/{name}")
}

/// The clause Xerces reports for a document that ends with an element open,
/// `input ended before all started tags were ended; last tag started is
/// '<tag>'`, out of a Release message; `None` for another parse error.
// Read only by the `mzml`-gated truncation test; unused on a build without it.
#[cfg_attr(not(feature = "mzml"), allow(dead_code))]
fn open_element_clause(message: &str) -> Option<&str> {
    const LEAD: &str = "input ended before all started tags were ended; last tag started is '";
    let start = message.find(LEAD)?;
    let tag_start = start + LEAD.len();
    let tag_end = tag_start + message[tag_start..].find('\'')?;
    Some(&message[start..=tag_end])
}

/// The input line a Release parse error names: `line (<n>)` in `MS2File`'s
/// messages (`FORMAT/MS2File.h:113`, `:144`), `line #<n>` in `MascotGenericFile`'s.
fn release_line(message: &str) -> usize {
    let digits = |rest: &str| -> usize {
        let end = rest
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(rest.len());
        rest[..end].parse().expect(message)
    };
    if let Some((_, rest)) = message.split_once("line (") {
        return digits(rest);
    }
    let (_, rest) = message
        .split_once("line #")
        .unwrap_or_else(|| panic!("no line number in {message:?}"));
    digits(rest)
}

/// Every truncated mzXML and mzData document is a parse error of the class,
/// as `MzXMLFile::load` and `MzDataFile::load` throw `ParseError` for each in
/// the Release build (FileInfo exit 3 on every `t_mzxml_*` and `t_mzdata_*`
/// tool case, `Parse Error` from the class driver on three). The truncations: after a complete record, inside a nested scan, in a
/// base64 payload, inside a start tag, before the first record, after the
/// record list, without the closing root tag, inside it, inside an mzXML scan
/// index, and a complete gzip member holding a cut document.
///
/// Before this repair, seven of the ten mzXML documents loaded, and FileInfo
/// printed a full report with exit 0. Where Xerces reports the element left
/// open, the port's message is its clause verbatim, without the line and
/// column Xerces adds. The mzData reader already refused every one of them.
#[cfg(feature = "mzml")]
#[test]
fn truncated_xml_peak_files_are_parse_errors_as_in_the_release_build() {
    let (mut refused, mut clauses) = (0, 0);
    for (id, release) in release_outcomes() {
        let stem = id.strip_prefix("lib_").unwrap_or(&id);
        let mzxml = stem.starts_with("t_mzxml_");
        if !mzxml && !stem.starts_with("t_mzdata_") {
            continue;
        }
        // The tool's bare runs and the class driver's (the `_all` runs repeat
        // the bare ones with flags the refusal never reaches).
        let lead = match release.program.as_str() {
            "tool" if !id.ends_with("_all") => "Error: Unable to read file (While loading '",
            "driver" => "Parse Error: While loading '",
            _ => continue,
        };
        assert_eq!(release.exit_code, 3, "{id}");
        assert!(
            release.last_line.starts_with(lead),
            "{id}: {}",
            release.last_line
        );
        if stem == "t_mzxml_gz_cut" {
            // A known gap outside the reader: see a_cut_gzip_member_is_refused.
            continue;
        }
        let input = data(&truncated_fixture(&id));
        match FileInfo::new().run(&input, &options("", NONE)) {
            Err(Error::Parse { message, .. }) => {
                let clause = open_element_clause(&release.last_line);
                if let Some(clause) = clause.filter(|_| mzxml) {
                    assert_eq!(message, clause, "{id}");
                    clauses += 1;
                }
            }
            other => panic!("{id}: {other:?}"),
        }
        refused += 1;
    }
    assert_eq!((refused, clauses), (18, 10));
}

/// `MzXMLFile::load` itself (driver `mzxml_driver.cpp`): a full load of a
/// document that ends with an element open throws `ParseError`; a
/// metadata-only load throws `EndParsingSoftly` at the first `<scan>`
/// (`FORMAT/HANDLERS/MzXMLHandler.cpp:242-245`), which `XMLFile::parse_` swallows
/// (`FORMAT/XMLFile.cpp:104-108`), so a document cut after that point loads with no
/// spectrum, and one cut before it (`t_mzxml_header`) throws.
#[cfg(feature = "mzml")]
#[test]
fn a_metadata_only_load_ends_before_the_truncation_as_in_the_release_build() {
    use openms::format::PeakFileOptions;
    use openms::format::mzxml;
    let mut seen = 0;
    for (id, release) in release_outcomes() {
        if release.program != "mzxml_driver" {
            continue;
        }
        seen += 1;
        let (mode, stem) = id
            .strip_prefix("mzxml_")
            .and_then(|rest| rest.split_once('_'))
            .expect(&id);
        let input = if stem == "MzXMLFile_1" {
            data("MzXMLFile_1.mzXML")
        } else {
            data(&truncated_fixture(stem))
        };
        let mut peaks = PeakFileOptions::default();
        peaks.metadata_only = mode == "meta";
        let options = mzxml::ReadOptions {
            peaks,
            ..mzxml::ReadOptions::default()
        };
        match (
            mzxml::load_with_options(&input, &options),
            release.exit_code,
        ) {
            (Ok(map), 0) => assert_eq!(
                format!(
                    "spectra {} instrument_name {}",
                    map.spectra.len(),
                    map.settings.instrument.name
                ),
                release.last_line,
                "{id}"
            ),
            (Err(Error::Parse { message, .. }), 3) => assert_eq!(
                Some(message.as_str()),
                open_element_clause(&release.last_line),
                "{id}"
            ),
            (other, code) => panic!("{id}: Release exit {code}, port {other:?}"),
        }
    }
    assert_eq!(seen, 8);
}

/// **Known gap, outside the readers:** a gzip member cut short. The Release
/// build's type sniffing reads at most 8191 decompressed bytes and takes what
/// the stream gives (`FORMAT/FileHandler.cpp:356-362`), and Xerces then parses the
/// decompressed prefix and throws `ParseError` (tool exit 3, class driver
/// `Parse Error`). Here `DocumentIdentifier::set_loaded_file_type` reads its
/// 64 KiB preview with `read_to_end`, which the gzip decoder fails with an
/// I/O error before the reader runs, so the class returns [`Error::Io`] and
/// the tool exits 8. The file is refused either way; only the exit code
/// differs. It is recorded in `docs/FILE_INFO_A8_SUPPORT.md` and in the
/// provenance manifest's `known_gaps`.
#[cfg(feature = "mzml")]
#[test]
fn a_cut_gzip_member_is_refused() {
    assert_eq!(release_outcome("t_mzxml_gz_cut").exit_code, 3);
    let class = release_outcome("lib_t_mzxml_gz_cut");
    assert_eq!(class.exit_code, 3);
    assert!(
        class.last_line.starts_with("Parse Error: "),
        "{}",
        class.last_line
    );
    let input = data(&truncated_fixture("t_mzxml_gz_cut"));
    let outcome = FileInfo::new().run(&input, &options("", NONE));
    assert!(outcome.is_err(), "{outcome:?}");
}

/// MGF has no closing element to miss, and what a cut does depends on where
/// it falls, in the Release build as here. Between two blocks, or in a
/// block's header before its first peak line, the file loads with the blocks
/// before the cut: the header loop runs to the end of the file, the outer
/// loop's `getline` fails too, and `getNextSpectrum_` returns false
/// (`FORMAT/MascotGenericFile.h:159-168`, `:399`). Both reports are the Release
/// build's byte for byte. Inside the peak list or its `END IONS` line, the cut
/// is a parse error on the line the Release build names (tool exit 3, and
/// `Parse Error` from the class driver on one).
#[test]
fn truncated_mgf_ends_as_in_the_release_build() {
    for id in ["t_mgf_between_blocks_all", "t_mgf_in_header_all"] {
        assert_eq!(release_outcome(id).exit_code, 0, "{id}");
        check_in(
            &(id, truncated_fixture(id), EVERYTHING, NONE),
            TRUNCATED_EXPECTED,
        );
    }
    for id in [
        "t_mgf_mid_peak",
        "t_mgf_end_cut",
        "t_mgf_end_word",
        "lib_t_mgf_end_cut",
    ] {
        let release = release_outcome(id);
        assert_eq!(release.exit_code, 3, "{id}");
        let input = data(&truncated_fixture(id));
        match FileInfo::new().run(&input, &options("", NONE)) {
            Err(Error::Parse { line, .. }) => {
                assert_eq!(line, release_line(&release.last_line), "{id}");
            }
            other => panic!("{id}: {other:?}"),
        }
    }
}

/// MS2 has no closing record either (driver cases, as the tool's `-in`
/// refuses the extension). A cut on a line boundary, or inside a number that
/// leaves both of a peak line's values, loads what is there, and both class
/// reports are the Release build's byte for byte; a cut that leaves an `S` line
/// three values or a peak line one is a parse error on the line the Release
/// build names (`FORMAT/MS2File.h:113`, `:144`).
#[test]
fn truncated_ms2_ends_as_in_the_release_build_class() {
    for (id, flags) in [
        ("lib_t_ms2_line_boundary", ""),
        ("lib_t_ms2_line_boundary_all", EVERYTHING),
        ("lib_t_ms2_mid_number", ""),
        ("lib_t_ms2_mid_number_all", EVERYTHING),
    ] {
        assert_eq!(release_outcome(id).exit_code, 0, "{id}");
        check_in(
            &(id, truncated_fixture(id), flags, NONE),
            TRUNCATED_EXPECTED,
        );
    }
    for id in ["lib_t_ms2_mid_s", "lib_t_ms2_mid_peak"] {
        let release = release_outcome(id);
        assert_eq!(release.exit_code, 3, "{id}");
        let input = data(&truncated_fixture(id));
        match FileInfo::new().run(&input, &options("", NONE)) {
            Err(Error::Parse { line, .. }) => {
                assert_eq!(line, release_line(&release.last_line), "{id}");
            }
            other => panic!("{id}: {other:?}"),
        }
    }
}

/// The FileInfo executable on every truncated mzXML, mzData and MGF input, as
/// the Release build ran it: bare, and with `-m -p -s -d -c -out -out_tsv`.
/// The exit code is the Release build's (except the known gap of
/// [`a_cut_gzip_member_is_refused`], which is refused all the same). A refused
/// file prints nothing on the output stream and `Error: Unable to read file
/// (` on the error stream, naming the element Xerces names for mzXML, and
/// leaves `-out` and `-out_tsv` as empty as the Release build does. A completed
/// run's reports, on the output stream or in the two files, are the Release
/// build's byte for byte.
#[cfg(all(feature = "mzml", feature = "paramxml", feature = "featurexml"))]
#[test]
fn the_tool_ends_every_truncated_input_as_the_release_build_does() {
    use openms::cli::tools::FileInfo as FileInfoTool;
    use openms::cli::{ExitCode, run_with};
    let dir = openms::system::file::TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let mut seen = 0;
    for (id, release) in release_outcomes() {
        if release.program != "tool" {
            continue;
        }
        seen += 1;
        let input = data(&truncated_fixture(&id));
        let input_text = input.to_str().unwrap().to_owned();
        let out = dir.path().join(format!("{id}.tmp.txt"));
        let out_tsv = dir.path().join(format!("{id}.tmp.tsv"));
        let mut args: Vec<String> = ["FileInfo", "-test", "-in", &input_text, "-no_progress"]
            .map(str::to_owned)
            .to_vec();
        if id.ends_with("_all") {
            args.extend(["-m", "-p", "-s", "-d", "-c"].map(str::to_owned));
            args.extend(["-out".to_owned(), out.to_str().unwrap().to_owned()]);
            args.extend(["-out_tsv".to_owned(), out_tsv.to_str().unwrap().to_owned()]);
        }
        let (mut stdout, mut stderr) = (Vec::new(), Vec::new());
        let code = run_with::<FileInfoTool>(&args, &mut stdout, &mut stderr);
        let (stdout, stderr) = (
            String::from_utf8(stdout).unwrap(),
            String::from_utf8(stderr).unwrap(),
        );
        let gap = id.starts_with("t_mzxml_gz_cut");
        if gap {
            assert_ne!(code, ExitCode::ExecutionOk, "{id}");
        } else {
            assert_eq!(code.as_i32(), release.exit_code, "{id}: {stderr}");
        }
        if release.exit_code != 0 {
            assert_eq!(release.stdout_bytes, 0, "{id}");
            assert!(stdout.is_empty(), "{id}: {stdout}");
            assert!(
                release
                    .last_line
                    .starts_with("Error: Unable to read file ("),
                "{id}"
            );
            if !gap {
                assert!(
                    stderr.starts_with("Error: Unable to read file ("),
                    "{id}: {stderr}"
                );
            }
            let clause = open_element_clause(&release.last_line);
            if let Some(clause) = clause.filter(|_| id.starts_with("t_mzxml_")) {
                assert!(stderr.contains(clause), "{id}: {stderr}");
            }
            for (name, bytes) in &release.cwd_files {
                let written = std::fs::metadata(dir.path().join(name))
                    .unwrap_or_else(|e| panic!("{id}: {name}: {e}"))
                    .len();
                assert_eq!(written, *bytes as u64, "{id}: {name}");
            }
        } else if id.ends_with("_all") {
            // The reports go to the two files, and only the timing line to
            // the output stream, which is all the Release build printed there.
            assert!(release.stdout_bytes > 0, "{id}");
            let (rest, took) = took_line::split_took_line("FileInfo", &stdout);
            assert!(rest.is_empty() && took.is_some(), "{id}: {stdout}");
            assert_eq!(release.cwd_files.len(), 2, "{id}");
            for (path, suffix) in [(&out, "txt"), (&out_tsv, "tsv")] {
                assert_report(
                    &oracle_path(&read_text(path), &input),
                    &data(&format!("{TRUNCATED_EXPECTED}/{id}.{suffix}")),
                    &format!("{id} {suffix}"),
                );
            }
        } else {
            // The report on the output stream, then the timing line, which
            // differs from run to run and is taken off both.
            let expected = read_text(&data(&format!("{TRUNCATED_EXPECTED}/{id}.stdout.txt")));
            let expected = expected
                .strip_suffix("FileInfo took <masked>\n")
                .unwrap_or_else(|| panic!("{id}: no timing line"));
            let (actual, took) = took_line::split_took_line("FileInfo", &stdout);
            assert!(took.is_some(), "{id}: {stdout}");
            let (actual, expected) = (
                normalise_file_name(&oracle_path(&actual, &input)),
                normalise_file_name(expected),
            );
            assert!(
                actual == expected,
                "{id} stdout: {}",
                first_difference(&actual, &expected)
            );
        }
    }
    assert_eq!(seen, 42);
}
