// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A6-FILEINFO: the FileInfo `-i`, `-d` and `-c` checks
//! (`FORMAT/FileInfo.cpp:827-846`, `:1779-1795`, `:1799-1848`, `:1851-1964`,
//! core `bc9cc12`; `OpenMS4-topp/src/FileInfo.cpp:118-148`, topp `174b576`).
//!
//! Evidence, in order of strength (see
//! `tests/data/file_info_checks_provenance.json` and
//! `docs/FILE_INFO_CHECKS_SUPPORT.md`):
//! - tier 1, executed differential: 59 cases of `../oracle/a6-fileinfo` run
//!   against the **Release** C++ FileInfo of
//!   `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576` on
//!   ibminode06, twice and reproduced. 38 of them have their `-out` and
//!   `-out_tsv` reports compared here byte for byte, one (`i_truncated_index`)
//!   its text and one (`i_window_below`) its report against the padded file's,
//!   which is native difference 11 itself; only the three lines that embed the
//!   input path are normalised — `File name: `, `general: file name` and the
//!   `-i` failure line — and the rest of the cases carry exit codes, counts and
//!   refusals;
//! - tier 1, retained upstream definition: TOPP_FileInfo_11
//!   (`topp/CMakeLists.txt:908-909`, test-data `0cb15f2`, `WILL_FAIL 1`) and
//!   TOPP_FileInfo_19 (`:928-930`) reproduced, the first including its exit
//!   code. TOPP_FileInfo_12 (`:910-911`) is not: its input's `charge array` is
//!   stored as 64-bit float, which the strict mzML reader refuses, so only its
//!   index is checked here;
//! - tier 4, the refusal of the one place the source's behaviour is undefined
//!   (an empty SRM chromatogram), the infinity and every NaN coordinate that
//!   are *not* refused, and the result fields the source leaves at their
//!   defaults. A NaN in either of the two `std::sort` calls was refused until
//!   lead decision D16; it is now reproduced, and the unit tests of
//!   `src/format/file_info/checks.rs` pin the permutation libstdc++ leaves.
//!
//! Four oracle cases load the original `FileInfo_9_input.mzML`, which the strict
//! mzML reader refuses for three reasons outside this package (see
//! `tests/file_info.rs`); `FileInfo_9_strict_reader.mzML`, the input A4 derived
//! for them, carries those cases here instead. The `-i` cases on that same file
//! need no substitute, because the check returns before the file is loaded.
//!
//! Every tool case runs in its own temporary directory, because the tests run in
//! parallel.

#![cfg(feature = "mzml")]

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

fn input(name: &str) -> PathBuf {
    data(&format!("file_info/inputs/{name}"))
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
            } else if line.starts_with("Could not detect a valid index for the mzML file ") {
                "Could not detect a valid index for the mzML file <input>\n".to_owned()
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

/// Run the library on `path` and compare both reports with the Release C++
/// output the oracle retained for `case`.
fn check(path: &Path, options: &Options, case: &str) -> FileInfoResult {
    let result = FileInfo::new()
        .run(path, options)
        .unwrap_or_else(|e| panic!("{case}: {e}"));
    assert_report(
        &result.text,
        &data(&format!("file_info_checks/expected/{case}.txt")),
        &format!("{case} text"),
    );
    assert_report(
        &result.tsv,
        &data(&format!("file_info_checks/expected/{case}.tsv")),
        &format!("{case} tsv"),
    );
    assert_eq!(FileInfo::to_text(&result), result.text);
    assert_eq!(FileInfo::to_tsv(&result), result.tsv);
    // FileInfo.h:206-217 declares both aggregates; report_ never fills either.
    assert_eq!(result.corruption, Default::default(), "{case} corruption");
    assert_eq!(result.detail, Default::default(), "{case} detail");
    result
}

/// The oracle ran the FileInfo tool, which reads a dangling mzML header
/// reference the way the source loader does (decision D10,
/// `src/cli/tools/file_info.rs`); the library default is strict. Every case
/// here therefore sets it, so the two sides load the same files.
fn flags(detailed: bool, check_corrupt: bool, check_index: bool) -> Options {
    Options {
        detailed,
        check_corrupt,
        check_index,
        source_dangling_references: true,
        ..Options::default()
    }
}

fn index_only() -> Options {
    flags(false, false, true)
}

fn detailed_only() -> Options {
    flags(true, false, false)
}

fn corrupt_only() -> Options {
    flags(false, true, false)
}

// ---------------------------------------------------------------------------
// -i, the indexed-mzML check (FORMAT/FileInfo.cpp:827-846)
// ---------------------------------------------------------------------------

/// A valid index, then the content of the file.
///
/// TOPP_FileInfo_12's own input cannot carry this case: its `charge array` is
/// stored as 64-bit float, which the strict mzML reader refuses (`A4`,
/// `tests/file_info.rs`), so the core `IndexedmzMLFile_1` fixture does. The
/// index of TOPP_FileInfo_12's input is still checked, in
/// [`index_check_reads_the_upstream_test_12_input`].
#[test]
fn index_valid_then_content() {
    let result = check(
        &data("indexed_mzml/IndexedmzMLFile_1.mzML"),
        &index_only(),
        "i_valid_indexed1",
    );
    assert!(result.validation.index_checked);
    assert!(result.validation.index_valid);
    assert_eq!(result.validation.indexed_spectra, 2);
    assert_eq!(result.validation.indexed_chromatograms, 1);
    // The rest of the ValidationInfo belongs to -v, which does not run.
    assert!(!result.validation.performed);
    assert!(result.validation.supported);
    assert!(!result.validation.valid);
}

/// TOPP_FileInfo_12's own input reaches the index check and passes it, but its
/// content cannot be read here: its `charge array` is stored as 64-bit float,
/// which the strict mzML reader refuses (`A4`, `tests/file_info.rs`). The index
/// itself is this package's part, and a unit test of
/// `src/format/file_info/checks.rs` asserts the three spectra and no
/// chromatogram the C++ reports for it.
#[test]
fn the_upstream_test_12_input_passes_the_index_but_not_the_reader() {
    let error = FileInfo::new()
        .run(input("FileInfo_12_input.mzML"), &index_only())
        .unwrap_err();
    assert!(
        matches!(&error, Error::Unsupported(message)
            if message == "canonical auxiliary array binary type"),
        "{error:?}"
    );
}

/// TOPP_FileInfo_11's input: no index. The source returns from `report_` there,
/// so the report ends with the failure text — no content and not even the two
/// trailing newlines.
#[test]
fn index_invalid_ends_the_report() {
    let result = check(
        &input("FileInfo_9_input.mzML"),
        &index_only(),
        "i_invalid_11",
    );
    assert!(result.validation.index_checked);
    assert!(!result.validation.index_valid);
    assert_eq!(result.validation.indexed_spectra, 0);
    assert_eq!(result.validation.indexed_chromatograms, 0);
    assert!(result.peak.is_none(), "no content was computed");
    assert!(!result.text.ends_with("\n\n"), "no trailing blank line");
    assert!(
        result
            .text
            .ends_with("Either the index is not present or is not correct.\n")
    );
}

/// The upstream TOPP_FileInfo_11 and _12 inputs are the same bytes; this port
/// keeps one copy, so the two cases must agree on it.
#[test]
fn the_upstream_11_and_9_inputs_are_one_file() {
    let bytes = std::fs::read(input("FileInfo_9_input.mzML")).unwrap();
    assert_eq!(bytes.len(), 22523);
}

/// A failed index check skips every other flag, because it returns first.
#[test]
fn index_invalid_skips_every_other_flag() {
    let options = Options {
        meta: true,
        processing: true,
        statistics: true,
        ..flags(true, true, true)
    };
    let result = check(
        &input("FileInfo_9_input.mzML"),
        &options,
        "i_invalid_11_all",
    );
    assert!(result.processing.is_empty());
    assert!(result.peak.is_none());
}

#[test]
fn index_valid_then_meta_processing_statistics() {
    let options = Options {
        meta: true,
        processing: true,
        statistics: true,
        ..index_only()
    };
    check(
        &data("indexed_mzml/IndexedmzMLFile_1.mzML"),
        &options,
        "i_valid_indexed1_mps",
    );
}

#[test]
fn index_valid_then_detailed_and_corrupt() {
    check(
        &data("indexed_mzml/IndexedmzMLFile_1.mzML"),
        &flags(true, true, true),
        "i_valid_indexed1_dc",
    );
}

/// An mzML with no index at all: the same failure as a truncated one.
#[test]
fn index_missing_on_an_empty_mzml() {
    check(&input("empty.mzML"), &index_only(), "i_invalid_empty");
}

#[test]
fn index_missing_on_the_derived_file_info_9() {
    check(
        &input("FileInfo_9_strict_reader.mzML"),
        &index_only(),
        "i_invalid_9_strict",
    );
}

/// An index section cut in half: the footer offset is gone with it, so the
/// source's `findIndexListOffset` fails and the check reports the same failure.
/// The C++ also prints a `IndexedMzMLDecoder::findIndexListOffset Error:` dump
/// of the searched bytes on stderr; this port has no such stream and writes
/// nothing there, which the report does not show either way.
#[test]
fn index_truncated_is_not_valid() {
    let directory = openms::system::file::TempDir::new(false).unwrap();
    let source = std::fs::read(input("FileInfo_12_input.mzML")).unwrap();
    let start = source
        .windows(11)
        .position(|window| window == b"<indexList ")
        .expect("the fixture is an indexed mzML");
    let path = directory.path().join("truncated.mzML");
    std::fs::write(&path, &source[..start + 200]).unwrap();

    let result = FileInfo::new().run(&path, &index_only()).unwrap();
    assert!(result.validation.index_checked);
    assert!(!result.validation.index_valid);
    assert_report(
        &result.text,
        &data("file_info_checks/expected/i_truncated_index.txt"),
        "i_truncated_index",
    );
}

/// The control for the two cases below: an indexed mzML of 1217 bytes, whose
/// footer is inside the 1023-byte window the source searches and whose
/// `<index name="spectrum">` is followed immediately by a newline, so the
/// source's DOM walk starts on a text node and loses no offset. Both
/// implementations parse it, and both reports agree byte for byte.
#[test]
fn index_inside_the_footer_window_agrees() {
    let result = check(
        &input("index_window_above.mzML"),
        &index_only(),
        "i_window_above",
    );
    assert!(result.validation.index_valid);
    assert_eq!(result.validation.indexed_spectra, 1);
    assert_eq!(result.validation.indexed_chromatograms, 0);
}

/// Native difference 11: the same file with its padding removed, 967 bytes.
///
/// `findIndexListOffset` seeks `-1023` from the end (`IndexedMzMLDecoder.cpp:165-168`),
/// which fails on a shorter file, so its regex searches a heap buffer nothing
/// wrote. The Release build therefore reports no index and exits
/// `ILLEGAL_PARAMETERS`, and its stderr dump of those bytes differs from run to
/// run — the oracle records both hashes under `indeterminate_stderr` and masks
/// the dump alone. There is nothing defined to reproduce, so this port reads
/// `min(length, 1023)` bytes instead and finds the index: its report for this
/// file is the C++ report for the padded one, line for line.
///
/// The difference belongs to `src/format/indexed_mzml.rs`, which every index
/// reader in the crate shares, not to this package.
#[test]
fn index_below_the_footer_window_diverges_from_the_source() {
    let cpp = read_text(&data("file_info_checks/expected/i_window_below.txt"));
    assert!(
        cpp.contains("Could not detect a valid index for the mzML file"),
        "the Release build reports no index below the window: {cpp}"
    );

    let result = FileInfo::new()
        .run(input("index_window_below.mzML"), &index_only())
        .expect("the port finds the index the source's failed seek hides");
    assert!(result.validation.index_valid);
    assert_eq!(result.validation.indexed_spectra, 1);
    assert_report(
        &result.text,
        &data("file_info_checks/expected/i_window_above.txt"),
        "i_window_below text against the padded file's C++ report",
    );
    assert_report(
        &result.tsv,
        &data("file_info_checks/expected/i_window_above.tsv"),
        "i_window_below tsv against the padded file's C++ report",
    );
}

/// Native difference 12: an index whose first child is an `<offset>`, i.e. one
/// with no whitespace immediately after the opening `<index ...>` tag.
///
/// `domParseIndexedEnd_` walks the children of each `<index>` as
/// `iter = getFirstChild(); while (iter != lastChild) { iter = getNextSibling(); ... }`
/// (`IndexedMzMLDecoder.cpp:280-282` sets `iter`, and the walk itself is
/// `:290-293`), advancing before it reads, so the first child is never looked
/// at. A text node — any whitespace — in that one position takes the place and
/// the walk loses nothing, which is why every index an OpenMS writer produces
/// parses; whitespace between the offsets or before `</index>` does not save
/// the first one. The Release build counts one spectrum in this two-offset
/// file; this port counts both.
///
/// The difference belongs to `src/format/indexed_mzml.rs`, not to this package.
#[test]
fn an_unspaced_index_keeps_the_offset_the_source_skips() {
    let cpp = read_text(&data("file_info_checks/expected/i_offsets_unspaced.txt"));
    assert!(
        cpp.contains("Found a valid indexed mzML XML File with 1 spectra and 0 chromatograms.\n"),
        "the Release build loses the first of the two offsets: {cpp}"
    );

    let result = FileInfo::new()
        .run(input("index_offsets_unspaced.mzML"), &index_only())
        .expect("the index parses");
    assert!(result.validation.index_valid);
    assert_eq!(result.validation.indexed_spectra, 2);
    assert_eq!(result.validation.indexed_chromatograms, 0);
    assert!(
        result
            .text
            .contains("Found a valid indexed mzML XML File with 2 spectra and 0 chromatograms.\n"),
        "{}",
        result.text
    );
}

/// The check reads the file itself, so a missing one is an I/O error, as the
/// source's `Exception::FileNotFound` is. The FileInfo tool never reaches this:
/// `registerInputFile_` refuses a missing `-in` first, exit 1.
#[test]
fn index_check_on_a_missing_file_is_an_io_error() {
    let directory = openms::system::file::TempDir::new(false).unwrap();
    let error = FileInfo::new()
        .run(directory.path().join("absent.mzML"), &index_only())
        .unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error:?}");
}

/// The source applies the check to any type; only the tool restricts `-i` to
/// mzML. A DTA has no footer, so the check fails and the report ends there,
/// before the branch this port does not refuse.
#[test]
fn index_check_runs_before_the_branch_is_judged() {
    let options = Options {
        forced_type: FileType::Dta,
        ..index_only()
    };
    let result = FileInfo::new()
        .run(input("FileInfo_1_input.dta"), &options)
        .unwrap();
    assert!(result.validation.index_checked);
    assert!(!result.validation.index_valid);
    assert!(result.text.contains("Could not detect a valid index"));
}

/// An unparsable index on a branch this port does not run is still reported in
/// full, because the source checks the index before it reaches the content.
/// `FileType` carries every type whatever the feature set, so the branch this
/// forces is refused in any build.
#[test]
fn index_check_precedes_an_unported_branch() {
    let options = Options {
        forced_type: FileType::ConsensusXml,
        ..index_only()
    };
    let result = FileInfo::new()
        .run(input("empty.mzML"), &options)
        .expect("the index check ends the report before the branch");
    assert!(!result.validation.index_valid);
}

// ---------------------------------------------------------------------------
// -d, the detailed listing (FORMAT/FileInfo.cpp:1779-1795 and :1799-1848)
// ---------------------------------------------------------------------------

/// TOPP_FileInfo_19 verbatim.
#[test]
fn detailed_faims_upstream_test_19() {
    check(
        &data("mzml_mobility/FAIMS_test_data.mzML"),
        &detailed_only(),
        "d_faims_19",
    );
}

/// SRM spectra become chromatograms on load, so the transition listing is all
/// there is: the experiment holds no spectrum and the spectrum listing is
/// suppressed.
#[test]
fn detailed_srm_transitions() {
    let result = check(&input("srm_spectra.mzML"), &detailed_only(), "d_srm");
    let text = &result.text;
    assert!(text.contains("\n -- Detailed chromatogram listing -- \n"));
    assert!(!text.contains("-- Detailed spectrum listing --"));
}

#[test]
fn detailed_srm_mixed_with_ordinary_spectra() {
    let result = check(
        &input("srm_spectra_mixed.mzML"),
        &detailed_only(),
        "d_srm_mixed",
    );
    assert!(
        result
            .text
            .contains(" -- Detailed chromatogram listing -- ")
    );
    assert!(result.text.contains("-- Detailed spectrum listing --"));
}

#[test]
fn detailed_spectra_and_a_chromatogram() {
    check(
        &data("indexed_mzml/IndexedmzMLFile_1.mzML"),
        &detailed_only(),
        "d_indexed1",
    );
}

/// The listing belongs to the content, so `-m`, `-p` and `-s` follow it.
#[test]
fn detailed_precedes_meta_processing_statistics() {
    let options = Options {
        meta: true,
        processing: true,
        statistics: true,
        ..detailed_only()
    };
    check(
        &data("indexed_mzml/IndexedmzMLFile_1.mzML"),
        &options,
        "d_indexed1_mps",
    );
}

#[test]
fn detailed_on_the_derived_file_info_9() {
    check(
        &input("FileInfo_9_strict_reader.mzML"),
        &detailed_only(),
        "d_9_strict",
    );
}

/// No spectrum and no chromatogram: neither block is written.
#[test]
fn detailed_on_an_empty_experiment() {
    let result = check(&input("empty.mzML"), &detailed_only(), "d_empty");
    assert!(!result.text.contains("Detailed"));
}

#[test]
fn detailed_on_an_mzml_without_a_run() {
    check(
        &input("MzMLFile_2_minimal.mzML"),
        &detailed_only(),
        "d_minimal",
    );
}

#[test]
fn detailed_on_a_dta() {
    let options = Options {
        forced_type: FileType::Dta,
        ..detailed_only()
    };
    check(&input("FileInfo_1_input.dta"), &options, "d_dta");
}

#[test]
fn detailed_on_a_dta2d() {
    check(
        &input("FileInfo_2_input.dta2d"),
        &detailed_only(),
        "d_dta2d",
    );
}

/// A drift time unit adds the `IM:` line; without one there is none.
#[test]
fn detailed_writes_the_ion_mobility_line() {
    let result = check(
        &data("faims_helper/IM_FAIMS_test.mzML"),
        &detailed_only(),
        "d_im_faims",
    );
    assert!(result.text.contains("\n  IM:         -55 FAIMS_CV\n"));
}

/// An empty spectrum writes no m/z extent and leaves the `m/z:` line open, so
/// the next line continues it. The source does this and the port keeps it.
#[test]
fn detailed_leaves_an_empty_spectrum_line_open() {
    let result = check(
        &input("corrupt_scans.mzML"),
        &detailed_only(),
        "d_corrupt_scans",
    );
    assert!(result.text.contains("  m/z:        Precursors:  0\n"));
    assert!(result.text.contains("  mslevel:    0\n"));
}

// ---------------------------------------------------------------------------
// -c, the corrupt-data check (FORMAT/FileInfo.cpp:1851-1964)
// ---------------------------------------------------------------------------

/// A duplicate m/z and a negative intensity. The mzML reader sorts a spectrum's
/// peaks on load (`MzMLHandler.cpp:218-221`), so the messages name the sorted
/// positions and the unsorted-peaks line cannot be reached through mzML.
#[test]
fn corrupt_duplicate_mz_and_negative_intensity() {
    let result = check(&input("corrupt_peaks.mzML"), &corrupt_only(), "c_peaks");
    assert!(
        result
            .text
            .contains("Warning: Negative peak intensity peak (RT: 10.5 MZ: 200 intensity: -2.5)\n")
    );
    assert!(
        result
            .text
            .contains("Error: Duplicate peak m/z 200 in spectrum (RT: 10.5)\n")
    );
}

/// Descending retention times, an MS-level-0 scan, an empty scan and two MS1
/// scans at the same retention time.
#[test]
fn corrupt_scan_level_problems() {
    let result = check(&input("corrupt_scans.mzML"), &corrupt_only(), "c_scans");
    for line in [
        "Error: Spectrum retention times are not sorted in ascending order\n",
        "Error: MS-level 0 in spectrum (RT: 30)\n",
        "Warning: No peaks in spectrum (RT: 40)\n",
        "Error: Duplicate spectrum retention time: 10\n",
    ] {
        assert!(result.text.contains(line), "missing {line:?}");
    }
}

/// The duplicate-data-array-name line of `-c` cannot be reached through this
/// port's mzML reader: `Record::check_array_kind` (`src/format/mzml.rs:1038`)
/// refuses a repeated auxiliary array name with
/// `Error::Parse("duplicate auxiliary array name")`, where the C++ reader loads
/// the file and leaves it to `-c` to report. The oracle records what the C++
/// writes for `corrupt_arrays.mzML` and `corrupt_arrays_mixed.mzML`, including
/// the name shared by a float and an integer array; the port's own rendering of
/// that line is covered by the unit tests of
/// `src/format/file_info/checks.rs`, on an experiment built in memory.
#[test]
fn a_duplicate_data_array_name_is_refused_by_the_reader() {
    for name in ["corrupt_arrays.mzML", "corrupt_arrays_mixed.mzML"] {
        let error = FileInfo::new()
            .run(input(name), &corrupt_only())
            .unwrap_err();
        assert!(
            matches!(&error, Error::Parse { message, .. }
                if message == "duplicate auxiliary array name"),
            "{name}: {error:?}"
        );
    }
}

/// Clean data still writes the header, and nothing after it.
#[test]
fn corrupt_check_on_clean_data() {
    let result = check(&input("clean_pair.mzML"), &corrupt_only(), "c_clean");
    let block = result
        .text
        .split_once("-- Checking for corrupt data --\n")
        .expect("the header is always written")
        .1;
    assert_eq!(block, "\n\n\n");
}

#[test]
fn corrupt_check_on_the_derived_file_info_9() {
    let result = check(
        &input("FileInfo_9_strict_reader.mzML"),
        &corrupt_only(),
        "c_9_strict",
    );
    assert!(
        result
            .text
            .contains("Warning: No peaks in spectrum (RT: 5.4)\n")
    );
}

#[test]
fn corrupt_check_on_an_indexed_mzml() {
    check(
        &data("indexed_mzml/IndexedmzMLFile_1.mzML"),
        &corrupt_only(),
        "c_indexed1",
    );
}

#[test]
fn corrupt_check_on_an_empty_experiment() {
    check(&input("empty.mzML"), &corrupt_only(), "c_empty");
}

#[test]
fn corrupt_check_on_a_dta() {
    let options = Options {
        forced_type: FileType::Dta,
        ..corrupt_only()
    };
    check(&input("FileInfo_1_input.dta"), &options, "c_dta");
}

/// A DTA2D whose scans are not in retention-time order, with one repeated.
#[test]
fn corrupt_check_on_a_dta2d() {
    let result = check(&input("FileInfo_2_input.dta2d"), &corrupt_only(), "c_dta2d");
    assert!(
        result
            .text
            .contains("Error: Spectrum retention times are not sorted in ascending order\n")
    );
    assert!(
        result
            .text
            .contains("Error: Duplicate spectrum retention time: 2\n")
    );
}

/// DTA and DTA2D are not sorted on load, so they are the only way to reach the
/// unsorted-peaks line.
#[test]
fn corrupt_unsorted_peaks_in_a_dta() {
    let options = Options {
        forced_type: FileType::Dta,
        ..corrupt_only()
    };
    let result = check(&input("unsorted_peaks.dta"), &options, "c_unsorted_dta");
    assert!(result.text.contains(
        "Error: Peak m/z positions are not sorted in ascending order in spectrum (RT: -1)\n"
    ));
}

#[test]
fn corrupt_unsorted_peaks_in_a_dta2d() {
    let result = check(
        &input("unsorted_peaks.dta2d"),
        &corrupt_only(),
        "c_unsorted_dta2d",
    );
    assert!(result.text.contains(
        "Error: Peak m/z positions are not sorted in ascending order in spectrum (RT: 1)\n"
    ));
}

/// Nothing is left to check once the SRM spectra have become chromatograms.
#[test]
fn corrupt_check_after_srm_conversion() {
    check(&input("srm_spectra.mzML"), &corrupt_only(), "c_srm");
}

#[test]
fn corrupt_check_on_the_faims_file() {
    check(
        &data("mzml_mobility/FAIMS_test_data.mzML"),
        &corrupt_only(),
        "c_faims",
    );
}

// ---------------------------------------------------------------------------
// -d and -c together, and the branches that ignore them
// ---------------------------------------------------------------------------

#[test]
fn detailed_precedes_the_corrupt_check() {
    let result = check(
        &input("FileInfo_9_strict_reader.mzML"),
        &flags(true, true, false),
        "dc_9_strict",
    );
    let listing = result.text.find("-- Detailed spectrum listing --").unwrap();
    let corrupt = result.text.find("-- Checking for corrupt data --").unwrap();
    assert!(listing < corrupt);
}

#[test]
fn every_flag_in_the_source_section_order() {
    let options = Options {
        meta: true,
        processing: true,
        statistics: true,
        ..flags(true, true, false)
    };
    check(
        &input("FileInfo_9_strict_reader.mzML"),
        &options,
        "dc_all_9_strict",
    );
}

#[test]
fn every_flag_on_srm_chromatograms() {
    let options = Options {
        meta: true,
        processing: true,
        statistics: true,
        ..flags(true, true, false)
    };
    check(&input("srm_spectra.mzML"), &options, "dc_all_srm");
}

#[test]
fn detailed_and_corrupt_on_a_dta2d() {
    check(
        &input("unsorted_peaks.dta2d"),
        &flags(true, true, false),
        "dc_unsorted_dta2d",
    );
}

#[test]
fn detailed_and_corrupt_on_degenerate_scans() {
    check(
        &input("corrupt_scans.mzML"),
        &flags(true, true, false),
        "dc_corrupt_scans",
    );
}

/// The source guards both flags inside the peak-file branch, so a featureXML
/// map ignores them and its report is the one it would have written without.
#[cfg(feature = "featurexml")]
#[test]
fn detailed_and_corrupt_are_ignored_on_featurexml() {
    let result = check(
        &input("empty.featureXML"),
        &flags(true, true, false),
        "dc_featurexml",
    );
    assert!(!result.text.contains("Detailed"));
    assert!(!result.text.contains("Checking for corrupt data"));
    let plain = FileInfo::new()
        .run(input("empty.featureXML"), &Options::default())
        .unwrap();
    assert_eq!(result.text, plain.text);
    assert_eq!(result.tsv, plain.tsv);
}

// ---------------------------------------------------------------------------
// Tier 4: where the source's behaviour is undefined, and what stays untouched
// ---------------------------------------------------------------------------

/// The one place the source's behaviour is undefined — an empty SRM
/// chromatogram — is refused; an infinity is not, and since lead decision D16
/// neither is a NaN in either `std::sort`, whose libstdc++ permutation the port
/// reproduces instead. All of those need an experiment no loader produces, so
/// they are unit tests of `src/format/file_info/checks.rs` instead of cases
/// here.
///
/// `-v` is the one flag of this block that is still refused.
#[test]
fn validation_is_still_refused() {
    let options = Options {
        validate: true,
        ..Options::default()
    };
    let error = FileInfo::new()
        .run(input("empty.mzML"), &options)
        .unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error:?}");
}

// ---------------------------------------------------------------------------
// The tool: the exit codes of -i (OpenMS4-topp/src/FileInfo.cpp:118-148)
// ---------------------------------------------------------------------------

#[cfg(all(feature = "paramxml", feature = "featurexml"))]
#[path = "support/took_line.rs"]
mod took_line;

/// Run the tool as `FileInfo args...`, as `tests/topp_file_info.rs` does, and
/// return its standard output without the closing `FileInfo took …` line,
/// which every one of these runs of the Release oracle ends with
/// (`i_invalid_11_bare`, `i_valid_indexed1`, `i_on_dta`), so it is required.
#[cfg(all(feature = "paramxml", feature = "featurexml"))]
fn run_tool(args: &[&str]) -> (openms::cli::ExitCode, String, String) {
    use openms::cli::{Tool, run_with};
    let arguments: Vec<String> = std::iter::once(openms::cli::tools::FileInfo::NAME.to_owned())
        .chain(args.iter().map(|a| (*a).to_owned()))
        .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<openms::cli::tools::FileInfo>(&arguments, &mut out, &mut err);
    let out = String::from_utf8(out).expect("UTF-8 output stream");
    let (out, took) = took_line::split_took_line("FileInfo", &out);
    assert!(took.is_some(), "no closing line: {out}");
    (
        code,
        out,
        String::from_utf8(err).expect("UTF-8 error stream"),
    )
}

/// TOPP_FileInfo_11 verbatim (`topp/CMakeLists.txt:908-909`, `WILL_FAIL 1`):
/// the mzML has no index, so the tool exits `ILLEGAL_PARAMETERS`, and the
/// Release oracle's `i_invalid_11_bare` exits 6 for the same input.
#[cfg(all(feature = "paramxml", feature = "featurexml"))]
#[test]
fn upstream_test_11_exits_illegal_parameters() {
    let path = input("FileInfo_9_input.mzML");
    let (code, out, err) =
        run_tool(&["-test", "-in", path.to_str().unwrap(), "-i", "-no_progress"]);
    assert_eq!(code, openms::cli::ExitCode::IllegalParameters, "{err}");
    assert!(
        out.contains("Could not detect a valid index for the mzML file"),
        "{out}"
    );
    assert!(
        out.ends_with("Either the index is not present or is not correct.\n"),
        "{out}"
    );
}

/// A valid index lets the run finish, so the tool exits 0, as the Release
/// oracle's `i_valid_indexed1` does. TOPP_FileInfo_12's own input
/// (`:910-911`, `WILL_FAIL 0`) cannot reach this: the strict mzML reader
/// refuses its float `charge array`.
#[cfg(all(feature = "paramxml", feature = "featurexml"))]
#[test]
fn a_valid_index_exits_zero() {
    let path = data("indexed_mzml/IndexedmzMLFile_1.mzML");
    let (code, _out, err) =
        run_tool(&["-test", "-in", path.to_str().unwrap(), "-i", "-no_progress"]);
    assert_eq!(code, openms::cli::ExitCode::ExecutionOk, "{err}");
}

/// The tool refuses `-i` on anything but mzML before the library runs, so the
/// index of a DTA is never looked at.
#[cfg(all(feature = "paramxml", feature = "featurexml"))]
#[test]
fn index_check_on_a_non_mzml_file_is_refused_by_the_tool() {
    let path = input("FileInfo_1_input.dta");
    let (code, out, err) = run_tool(&[
        "-test",
        "-in",
        path.to_str().unwrap(),
        "-in_type",
        "dta",
        "-i",
        "-no_progress",
    ]);
    assert_eq!(code, openms::cli::ExitCode::IllegalParameters);
    assert!(out.is_empty(), "{out}");
    assert!(
        err.starts_with("Error: Can only validate indices for mzML files\n"),
        "{err}"
    );
}
