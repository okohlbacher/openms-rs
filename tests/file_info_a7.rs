// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A7-FILEINFO: the FileInfo consensusXML, identification and FASTA branches
//! (`FORMAT/FileInfo.cpp:853-1076`, `:1146-1311`, `:1312-1470` and their `-m`,
//! `-p` and `-s` arms at `:1985-2004`, `:2101-2114` and `:2257-2379`, core
//! `bc9cc12`; `OpenMS4-topp/src/FileInfo.cpp`, topp `174b576`).
//!
//! Evidence, in order of strength (see
//! `tests/data/file_info_a7_provenance.json` and
//! `docs/FILE_INFO_A7_SUPPORT.md`):
//!
//! - tier 1, executed differential: 60 cases of `../oracle/a7-fileinfo` run
//!   against the **Release** C++ FileInfo of
//!   `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576` on
//!   ibminode06, twice and reproduced. 40 of them have their `-out` and
//!   `-out_tsv` reports compared here byte for byte; only the two lines that
//!   embed the input path are normalised, `File name: ` and
//!   `general: file name`;
//! - tier 1, retained upstream definition: TOPP_FileInfo_7, _10, _13, _17, _18
//!   and _20 (`topp/CMakeLists.txt:899-901`, `:905-907`, `:912`, `:922-927`,
//!   `:931-933`,
//!   test-data `0cb15f2`) reproduced with their own flags. TOPP_FileInfo_14 and
//!   _15 pass `-v`, which is not ported, so their inputs are exercised without
//!   it;
//! - tier 1, executed probe: `oracle/a7-fileinfo/scripts/probe_std_hash.cpp`,
//!   which pins libstdc++'s `std::hash<std::string>` — the FASTA duplicate
//!   detection is sensitive to it, because the source overwrites each hash
//!   bucket instead of appending to it;
//! - tier 4, the refusal of the three places the source's behaviour is an
//!   out-of-bounds `std::vector` access (lead decision D1), each of which the
//!   Release build answers with a segmentation fault or with a silently wrong
//!   number.
//!
//! Two oracle cases have no differential here: `FileFilter_25_input.idXML` and
//! `FalseDiscoveryRate_5_input.idXML` carry more modified peptide hits than the
//! shared identification-XML reader's document-wide work budget allows, which
//! is a limit of that reader and not of this branch.

#[cfg(any(feature = "consensusxml", feature = "idxml"))]
use openms::Error;
#[cfg(feature = "idxml")]
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

/// Run the library on `name` and compare both reports with the Release C++
/// output the oracle retained for `case`.
fn check(name: &str, options: &Options, case: &str) -> FileInfoResult {
    let result = FileInfo::new()
        .run(input(name), options)
        .unwrap_or_else(|e| panic!("{case}: {e}"));
    assert_report(
        &result.text,
        &data(&format!("file_info_a7/expected/{case}.txt")),
        &format!("{case} text"),
    );
    assert_report(
        &result.tsv,
        &data(&format!("file_info_a7/expected/{case}.tsv")),
        &format!("{case} tsv"),
    );
    assert_eq!(FileInfo::to_text(&result), result.text);
    assert_eq!(FileInfo::to_tsv(&result), result.tsv);
    // FileInfo.h:206-217 declares both aggregates; none of these branches fills
    // either, because -d and -c live inside the peak-file arm alone.
    assert_eq!(result.corruption, Default::default(), "{case} corruption");
    assert_eq!(result.detail, Default::default(), "{case} detail");
    assert_eq!(result.validation, Default::default(), "{case} validation");
    result
}

fn bare() -> Options {
    Options::default()
}

fn all_flags() -> Options {
    Options {
        meta: true,
        processing: true,
        statistics: true,
        ..Options::default()
    }
}

// ---------------------------------------------------------------------------
// FASTA
// ---------------------------------------------------------------------------

/// TOPP_FileInfo_17 (`topp/CMakeLists.txt:922-924`): eleven amino-acid
/// sequences with two duplicate headers and one duplicate sequence.
#[test]
fn fasta_upstream_17() {
    let result = check("FileInfo_17_input.fasta", &bare(), "f17");
    let fasta = result.fasta.expect("the FASTA branch fills fasta");
    assert!(!fasta.is_nucleic_acid);
    assert_eq!(fasta.num_sequences, 11);
    assert_eq!(fasta.total_residues, 1933);
    assert_eq!(fasta.dup_headers, 2);
    assert_eq!(fasta.dup_sequences, 1);
    assert_eq!(fasta.seq_with_ambiguous, 0);
    assert_eq!(fasta.ambiguity_counts["BZX"], 0);
    assert_eq!(fasta.ambiguity_counts["BZXJ"], 0);
    // FileInfo.cpp:1055-1057: three or more sequences use SummaryStatistics.
    assert_eq!(fasta.length_stats.count, 11);
    assert_eq!(fasta.length_stats.min, 174.0);
    assert_eq!(fasta.length_stats.max, 177.0);
    // FileInfo.cpp:1045: the print loop's operator[] inserts no zero-count key
    // for an amino-acid file, and the structured map skips zero counts anyway.
    assert!(fasta.residue_counts.values().all(|&count| count != 0));
    assert_eq!(fasta.residue_counts[&b'*'], 2);
    assert_eq!(fasta.residue_counts[&b'-'], 5);
}

/// The duplicate warnings the source writes with `OPENMS_LOG_WARN` rather than
/// into the report, in the order it writes them. `crab_chick` at index 7
/// matches index 2, five entries earlier, which is what the hash bucket has to
/// carry.
#[test]
fn fasta_duplicate_warnings_match_the_source_log() {
    let result = FileInfo::new()
        .run(input("FileInfo_17_input.fasta"), &bare())
        .expect("FileInfo_17 loads");
    assert_eq!(
        result.warnings,
        vec![
            "Warning: Duplicate header, #7, ID: crab_chick = #2, ID: crab_chick".to_owned(),
            "Warning: Duplicate header, #10, ID: crab_squac = #9, ID: crab_squac".to_owned(),
            "Warning: Duplicate sequence, #10, ID: crab_squac == #9, ID: crab_squac".to_owned(),
        ]
    );
}

/// Three identical entries. The source assigns each hash bucket a one-element
/// vector instead of appending to it (`FileInfo.cpp:931`, `:949`), so entry #1
/// is compared with #0 and entry #2 with #1, and both count.
#[test]
fn fasta_three_identical_entries_count_two_duplicates() {
    let result = check("a7_fasta_dup3.fasta", &bare(), "f_dup3");
    let fasta = result.fasta.expect("fasta");
    assert_eq!(fasta.dup_headers, 2);
    assert_eq!(fasta.dup_sequences, 2);
    assert_eq!(
        result.warnings,
        vec![
            "Warning: Duplicate header, #1, ID: DUP = #0, ID: DUP".to_owned(),
            "Warning: Duplicate sequence, #1, ID: DUP == #0, ID: DUP".to_owned(),
            "Warning: Duplicate header, #2, ID: DUP = #1, ID: DUP".to_owned(),
            "Warning: Duplicate sequence, #2, ID: DUP == #1, ID: DUP".to_owned(),
        ]
    );
}

/// `FASTAEntry::headerMatches` compares identifier *and* description, so two
/// entries that share only the identifier are not duplicate headers, while
/// their identical sequences are duplicate sequences.
#[test]
fn fasta_header_match_uses_the_description_too() {
    let result = check(
        "a7_fasta_dup_header_only.fasta",
        &bare(),
        "f_dup_header_only",
    );
    let fasta = result.fasta.expect("fasta");
    assert_eq!(fasta.dup_headers, 0);
    assert_eq!(fasta.dup_sequences, 1);
}

/// TOPP_FileInfo_18 (`topp/CMakeLists.txt:925-927`) and TOPP_FileInfo_20
/// (`:931-933`).
#[test]
fn fasta_upstream_18_and_20() {
    check("FileInfo_18_input.fasta", &bare(), "f18");
    check("FileInfo_20_input.fasta", &bare(), "f20");
}

/// One and two sequences: below three, `FileInfo.cpp:986-1006` never asks for
/// quartiles and prints the minimum and the maximum in their place, and
/// `:1051-1064` fills the structured statistics field by field.
#[test]
fn fasta_below_three_sequences_falls_back_to_the_extremes() {
    let single = check("a7_fasta_single.fasta", &bare(), "f_single")
        .fasta
        .expect("fasta");
    assert_eq!(single.length_stats.count, 1);
    assert_eq!(single.length_stats.lowerq, single.length_stats.min);
    assert_eq!(single.length_stats.upperq, single.length_stats.max);
    assert_eq!(single.length_stats.variance, 0.0);

    let two = check("a7_fasta_two.fasta", &bare(), "f_two")
        .fasta
        .expect("fasta");
    assert_eq!(two.length_stats.count, 2);
    assert_eq!(two.length_stats.lowerq, 40.0);
    assert_eq!(two.length_stats.upperq, 60.0);
    assert_eq!(two.length_stats.median, 50.0);

    let one_residue = check("a7_fasta_one_residue.fasta", &bare(), "f_one_residue")
        .fasta
        .expect("fasta");
    assert_eq!(one_residue.total_residues, 1);
    assert_eq!(one_residue.length_stats.count, 1);
}

/// Three sequences: the first case with real quartiles.
#[test]
fn fasta_three_sequences_use_the_quantiles() {
    let fasta = check("a7_fasta_three.fasta", &bare(), "f_three")
        .fasta
        .expect("fasta");
    assert_eq!(fasta.length_stats.count, 3);
    assert_eq!(fasta.length_stats.min, 20.0);
    assert_eq!(fasta.length_stats.max, 60.0);
    assert_eq!(fasta.length_stats.median, 40.0);
}

/// Every sequence byte inside the IUPAC nucleotide alphabet selects the
/// nucleotide labels and the `N` / IUPAC totals (`FileInfo.cpp:900-912`).
#[test]
fn fasta_nucleic_acid_detection_and_ambiguity_buckets() {
    let clean = check("a7_fasta_nucleic.fasta", &bare(), "f_nucleic")
        .fasta
        .expect("fasta");
    assert!(clean.is_nucleic_acid);
    assert_eq!(clean.seq_with_ambiguous, 0);
    assert_eq!(clean.ambiguity_counts["N"], 0);
    assert_eq!(clean.ambiguity_counts["IUPAC"], 0);

    let ambiguous = check("a7_fasta_nucleic_ambig.fasta", &bare(), "f_nucleic_ambig")
        .fasta
        .expect("fasta");
    assert!(ambiguous.is_nucleic_acid);
    // Two of the three sequences carry at least one ambiguity code.
    assert_eq!(ambiguous.seq_with_ambiguous, 2);
    assert_eq!(ambiguous.ambiguity_counts["N"], 4);
    assert_eq!(ambiguous.ambiguity_counts["IUPAC"], 24);
    // FileInfo.cpp:1040-1041 skips the zero-count keys that the verbatim print
    // code may have operator[]-inserted for 'N' and 'n'.
    assert!(ambiguous.residue_counts.values().all(|&count| count != 0));
}

/// One byte outside the nucleotide alphabet makes the whole file amino acid,
/// and then `(B/Z/X)` and `(B/Z/X/J)` differ.
#[test]
fn fasta_amino_acid_ambiguity_buckets_differ_on_j() {
    let fasta = check("a7_fasta_ambig_aa.fasta", &bare(), "f_ambig_aa")
        .fasta
        .expect("fasta");
    assert!(!fasta.is_nucleic_acid);
    assert_eq!(fasta.seq_with_ambiguous, 2);
    assert_eq!(fasta.ambiguity_counts["BZX"], 6);
    assert_eq!(fasta.ambiguity_counts["BZXJ"], 8);
}

/// `std::map<char, int>` keeps upper and lower case apart, and the ASCII order
/// of the printed table is the signed-char order the source iterates in.
#[test]
fn fasta_residue_table_separates_case() {
    let fasta = check("a7_fasta_lowercase.fasta", &bare(), "f_lowercase")
        .fasta
        .expect("fasta");
    assert_eq!(fasta.residue_counts[&b'M'], 1);
    assert_eq!(fasta.residue_counts[&b'm'], 1);
    assert_eq!(fasta.residue_counts[&b'L'], 2);
    assert_eq!(fasta.residue_counts[&b'l'], 2);
}

/// `-m`, `-p` and `-s` have empty FASTA arms (`FileInfo.cpp:2001-2004`,
/// `:2112-2114`, `:2377-2379`), so `-m` and `-s` add only their titles and `-p`
/// the no-information line.
#[test]
fn fasta_flag_sections_are_titles_only() {
    let result = check("FileInfo_17_input.fasta", &all_flags(), "f17_all");
    assert!(result.processing.is_empty());
    assert!(result.statistics.is_empty());
    assert!(result.experiment_meta.is_none());
    assert!(result.text.contains("-- Meta information --\n\n\n"));
    assert!(
        result
            .text
            .contains("No information about data processing available!")
    );
    assert!(result.text.ends_with("-- Statistics --\n\n\n\n"));
}

// ---------------------------------------------------------------------------
// consensusXML
// ---------------------------------------------------------------------------

/// TOPP_FileInfo_7 (`topp/CMakeLists.txt:899-901`) with its own `-s -m -p`.
///
/// The upstream reference output records the source defect this reproduces:
/// five consensus features, `Intensities ... num. of values: 5` and
/// `Qualities ... num. of values: 10`, because `FileInfo.cpp:2263-2266`
/// declares `vector<double> qualities(size)` and then appends to it.
#[cfg(feature = "consensusxml")]
#[test]
fn consensus_upstream_7() {
    let result = check("FileInfo_7_input.consensusXML", &bare(), "c7");
    let feature = result
        .feature
        .expect("the consensusXML branch fills feature");
    assert!(feature.is_consensus);
    assert_eq!(feature.num_features, 5);
    assert_eq!(feature.size_distribution[&2], 5);
    assert_eq!(feature.assigned_ids, 0);
    assert_eq!(feature.unassigned_ids, 0);
    assert_eq!(feature.map_columns.len(), 2);
    assert_eq!(feature.map_columns[0].identifier, "0");
    assert_eq!(feature.map_columns[0].label, "light");
    assert_eq!(feature.map_columns[1].label, "heavy");
    assert_eq!(feature.map_columns[0].size, 16);
    // A consensus map has no mobility dimension.
    assert!(!result.ranges.is_experiment);
    assert!(!result.ranges.combined.has_mobility);
    assert!(result.ranges.combined.mobility.is_none());
}

#[cfg(feature = "consensusxml")]
#[test]
fn consensus_upstream_7_with_all_flags() {
    let result = check("FileInfo_7_input.consensusXML", &all_flags(), "c7_all");
    assert_eq!(result.processing.len(), 1);
    assert_eq!(result.processing[0].software_name, "FileFilter");
    // The source's own defect, asserted on the rendered numbers rather than
    // only on the recorded text: eleven blocks, the second of which summarises
    // twice as many values as there are consensus features.
    assert!(
        result
            .text
            .contains("Intensities of consensus features:\n  num. of values: 5\n")
    );
    assert!(
        result
            .text
            .contains("Qualities of consensus features:\n  num. of values: 10\n")
    );
    // FileInfo.cpp:2257-2372 writes none of the -s block to the TSV report.
    assert!(!result.tsv.contains("statistics: "));
}

/// Each of `-m`, `-p` and `-s` alone, so that the section order and the
/// consensusXML arm of each is pinned separately. `-m` has no TSV twin here,
/// unlike the featureXML arm.
#[cfg(feature = "consensusxml")]
#[test]
fn consensus_flag_sections_individually() {
    let meta = check(
        "FileInfo_7_input.consensusXML",
        &Options {
            meta: true,
            ..Options::default()
        },
        "c7_m",
    );
    assert!(meta.text.contains("Document ID: cons\n"));
    assert!(!meta.tsv.contains("document ID"));

    check(
        "FileInfo_7_input.consensusXML",
        &Options {
            processing: true,
            ..Options::default()
        },
        "c7_p",
    );
    check(
        "FileInfo_7_input.consensusXML",
        &Options {
            statistics: true,
            ..Options::default()
        },
        "c7_s",
    );
}

/// TOPP_FileInfo_13 (`topp/CMakeLists.txt:912`, "empty file should not
/// crash"): no consensus feature, so the histogram, the peptide rows, the
/// totals and the ranges are all skipped and two lines replace them, while the
/// column headers and the identification counts are still written.
#[cfg(feature = "consensusxml")]
#[test]
fn consensus_empty_map() {
    let result = check("FileInfo_13_input.consensusXML", &bare(), "c13_out");
    let feature = result.feature.expect("feature");
    assert_eq!(feature.num_features, 0);
    assert!(feature.size_distribution.is_empty());
    assert_eq!(feature.map_columns.len(), 2);
    assert!(result.text.contains("Number of consensus features: 0\n"));
    assert!(
        result
            .text
            .contains("No consensus features found, map is empty!\n")
    );
    assert!(result.text.contains("File descriptions:\n"));
    // The ranges block belongs to the non-empty arm.
    assert!(!result.text.contains("Ranges:"));
    assert!(!result.tsv.contains("general: ranges: "));

    // Eleven all-zero statistics blocks over no values at all.
    let all = check("FileInfo_13_input.consensusXML", &all_flags(), "c13_all");
    assert_eq!(all.processing.len(), 4);
    assert_eq!(all.text.matches("  num. of values: 0\n").count(), 11);
}

/// Consensus features with peptide identifications: the with-at-least-one-ID
/// columns and the aggregated peptide rows of `FileInfo.cpp:1193-1209`.
#[cfg(feature = "consensusxml")]
#[test]
fn consensus_identification_columns_and_peptide_rows() {
    let result = check("a7_cons_ids.consensusXML", &bare(), "c_ids");
    let feature = result.feature.expect("feature");
    assert_eq!(feature.num_features, 4);
    assert_eq!(feature.size_distribution[&2], 3);
    assert_eq!(feature.size_distribution[&1], 1);
    assert_eq!(feature.assigned_ids, 3);
    assert!(
        result.text.contains(
            "  of size 2: 3\t (features: 6 )\t with at least one ID: 2\t (features: 4 )\n"
        )
    );
    // One peptide in two maps with four sub-features, one in a single map.
    assert!(result.text.contains(
        "  peptides (with different mod. and charge) observed in 2 maps: 1\t (features: 4 )\n"
    ));
    assert!(result.text.contains(
        "  peptides (with different mod. and charge) observed in 1 maps: 1\t (features: 1 )\n"
    ));
    assert!(
        result
            .text
            .contains("  total consensus features:    4  with at least one ID: 3\n")
    );
    assert!(
        result
            .text
            .contains("  total features:              7  with at least one ID:  5\n")
    );
    check("a7_cons_ids.consensusXML", &all_flags(), "c_ids_all");
}

/// A consensus feature of size 12 makes `field_width` two characters, which
/// right-aligns the size column and widens the spacer on the total-features
/// line. The width is the source's `largest / 10 + 1`, not a digit count.
#[cfg(feature = "consensusxml")]
#[test]
fn consensus_size_column_width() {
    let result = check("a7_cons_wide.consensusXML", &bare(), "c_wide");
    assert!(result.text.contains("  of size 12: 1\t"));
    assert!(result.text.contains("  of size  2: 1\t"));
    assert!(result.text.contains("  of size  1: 1\t"));
    assert!(result.text.contains(
        "  peptides (with different mod. and charge) observed in 12 maps: 1\t (features: 14 )\n"
    ));
    assert!(
        result
            .text
            .contains("  total features:              15  with at least one ID:   14\n")
    );
}

/// A single consensus feature with a single sub-feature, and a map with column
/// headers but no consensus feature.
#[cfg(feature = "consensusxml")]
#[test]
fn consensus_degenerate_maps() {
    let one = check("a7_cons_one.consensusXML", &all_flags(), "c_one");
    assert_eq!(one.feature.expect("feature").num_features, 1);
    let headers = check(
        "a7_cons_headers_only.consensusXML",
        &all_flags(),
        "c_headers_only",
    );
    let feature = headers.feature.expect("feature");
    assert_eq!(feature.num_features, 0);
    assert_eq!(feature.map_columns.len(), 2);
}

/// Upstream consensus maps whose map ids do run from zero, so the occurrence
/// vector is indexed in bounds throughout.
#[cfg(feature = "consensusxml")]
#[test]
fn consensus_upstream_maps_with_zero_based_ids() {
    check("Epifany_2_input.consensusXML", &all_flags(), "c_epifany2");
    check("ConsensusXMLFile_1.consensusXML", &all_flags(), "c_cxml1");
}

/// A sub-feature whose map index is outside the column headers is only reached
/// when its consensus feature also carries an identification, so this file —
/// the same map index without one — is reported normally, exactly as the
/// Release build reports it.
#[cfg(feature = "consensusxml")]
#[test]
fn consensus_out_of_range_map_index_without_an_identification_is_reported() {
    check(
        "a7_cons_mapindex_high_noid.consensusXML",
        &all_flags(),
        "c_mapindex_high_noid",
    );
}

/// D1: `FileInfo.cpp:1176-1183` sizes the occurrence vector from the number of
/// column headers and then indexes it with `FeatureHandle::getMapIndex()`,
/// which is the file's `map=` id rather than a position. Every file below puts
/// that index outside the vector, and the port refuses instead of reproducing
/// what the adjacent heap happened to hold.
///
/// The Release build answers `a7_cons_no_headers.consensusXML` with a
/// segmentation fault; it answers the other three with an exit code of zero and
/// a wrong peptide row.
#[cfg(feature = "consensusxml")]
#[test]
fn consensus_out_of_bounds_map_index_is_refused() {
    for name in [
        "a7_cons_no_headers.consensusXML",
        "a7_cons_mapindex_high.consensusXML",
        "a7_cons_ids_one_based.consensusXML",
        // An upstream fixture: ConsensusID_3 numbers its two maps 1 and 2, so
        // the Release build writes past a two-element vector and reports
        // `observed in 1 maps` for peptides that are in both.
        "ConsensusID_3_input.consensusXML",
    ] {
        let error = FileInfo::new()
            .run(input(name), &bare())
            .expect_err(&format!("{name} must be refused"));
        let Error::InvalidValue(message) = &error else {
            panic!("{name}: expected InvalidValue, got {error}");
        };
        assert!(
            message.contains("column header") && message.contains("out of bounds"),
            "{name}: {message}"
        );
    }
}

// ---------------------------------------------------------------------------
// idXML and mzIdentML
// ---------------------------------------------------------------------------

/// TOPP_FileInfo_10 (`topp/CMakeLists.txt:905-907`): one run with no hits.
#[cfg(feature = "idxml")]
#[test]
fn identifications_upstream_10() {
    let result = check("FileInfo_10_input.idXML", &bare(), "id10");
    let ident = result.ident.expect("the branch fills ident");
    assert_eq!(ident.num_runs, 1);
    assert_eq!(ident.protein_hits, 0);
    assert_eq!(ident.matched_spectra, 0);
    assert_eq!(ident.peptide_hits, 0);
    assert_eq!(
        ident.search_engines,
        vec!["Unknown (version: 0)".to_owned()]
    );
    // FileInfo.cpp:1396-1399 pushes a single zero so that the mean is defined.
    assert_eq!(ident.avg_peptide_length, 0.0);
    assert_eq!(ident.psms_per_spectrum, 0.0);
    assert!(ident.modification_counts.is_empty());
    // A hit-less run still prints no percentage after the modified-top-hit line.
    assert!(result.text.contains("  modified top-hits:          0/0\n"));
}

/// `-m` on idXML prints the document identifier and writes a TSV line that,
/// uniquely, has no trailing newline (`FileInfo.cpp:1994-1995`). `-p` and `-s`
/// have an idXML arm that contributes nothing beyond their titles.
#[cfg(feature = "idxml")]
#[test]
fn identifications_flag_sections() {
    let result = check("FileInfo_10_input.idXML", &all_flags(), "id10_all");
    assert!(result.processing.is_empty());
    assert!(result.statistics.is_empty());
    assert!(result.tsv.ends_with("meta: document ID\t"));
    assert!(
        result
            .text
            .contains("No information about data processing available!")
    );
    assert!(result.text.ends_with("-- Statistics --\n\n\n\n"));
}

/// `PSMs / spectrum` is an integer division in the text while the structured
/// field carries the real ratio (`FileInfo.cpp:1416`, `:1464`).
#[cfg(feature = "idxml")]
#[test]
fn identifications_psms_per_spectrum_is_truncated_in_the_text_only() {
    let result = check("a7_id_hits.idXML", &all_flags(), "id_hits");
    let ident = result.ident.expect("ident");
    assert_eq!(ident.matched_spectra, 3);
    assert_eq!(ident.peptide_hits, 7);
    assert!(
        result
            .text
            .contains("  PSMs / spectrum (ignoring unidentified spectra):    2\n")
    );
    assert!((ident.psms_per_spectrum - 7.0 / 3.0).abs() < 1e-12);
    // (8 + 8 + 8 + 9 + 8 + 6 + 14) / 7 = 8.714..., rounded to 9 in the text.
    assert!(result.text.contains("(avg. length: 9)"));
    assert!((ident.avg_peptide_length - 61.0 / 7.0).abs() < 1e-12);
}

/// A terminal modification is counted under `getId()` and a residue
/// modification under `getFullId()` (`FileInfo.cpp:1353-1372`), and the
/// modification line ends without a newline.
#[cfg(feature = "idxml")]
#[test]
fn identifications_modification_counts_use_two_identities() {
    let result = check("a7_id_mods.idXML", &all_flags(), "id_mods");
    let ident = result.ident.expect("ident");
    assert_eq!(ident.modified_tophits, 4);
    assert_eq!(ident.matched_spectra, 5);
    assert_eq!(ident.modification_counts["Oxidation (M)"], 2);
    assert_eq!(ident.modification_counts["Carbamidomethyl (C)"], 1);
    // Terminal modifications keep the bare identifier, with no origin suffix.
    assert_eq!(ident.modification_counts["Dimethyl"], 1);
    assert_eq!(ident.modification_counts["Amidated"], 1);
    assert!(result.text.contains(
        "  Modification count (top-hits only): Amidated 1, Carbamidomethyl (C) 1, \
         Dimethyl 1, Oxidation (M) 2\n-- Meta information --"
    ));
    // FileInfo.cpp:1418: the percentage goes through StringUtils::toStr.
    assert!(
        result
            .text
            .contains("  modified top-hits:          4/5 (80.0%)\n")
    );
}

/// Upstream identification files: several runs and protein hits, and
/// N-terminal modifications.
#[cfg(feature = "idxml")]
#[test]
fn identifications_upstream_files() {
    let several_runs = check("ConsensusID_1_input.idXML", &all_flags(), "id_cid1")
        .ident
        .expect("ident");
    assert_eq!(several_runs.num_runs, 3);
    assert_eq!(several_runs.protein_hits, 3);
    check("IDFileConverter_10_output.idXML", &all_flags(), "id_conv10");
}

/// One run with no hits and no identifications at all: `proteins[0]` exists, so
/// the branch runs to the end.
#[cfg(feature = "idxml")]
#[test]
fn identifications_run_without_hits() {
    let ident = check("a7_id_empty_run.idXML", &all_flags(), "id_empty_run")
        .ident
        .expect("ident");
    assert_eq!(ident.num_runs, 1);
    assert_eq!(ident.db_name, "a7.fasta");
    assert_eq!(ident.db_version, "7");
    assert_eq!(ident.taxonomy, "9606");
}

/// mzIdentML content, and the fall-through that `-m`, `-p` and `-s` give it:
/// none of the three has an mzIdentML arm, so all three land in the peak-file
/// arm and report the `MSExperiment` this branch never loaded.
#[cfg(feature = "idxml")]
#[test]
fn identifications_mzidentml_and_its_peak_file_fall_through() {
    let plain = check("FileInfo_14_input.mzid", &bare(), "mzid14");
    let ident = plain.ident.expect("ident");
    assert_eq!(ident.num_runs, 4);
    assert_eq!(ident.protein_hits, 67);
    assert_eq!(ident.non_redundant_protein_hits, 44);
    assert_eq!(ident.matched_spectra, 69);

    let all = check("FileInfo_14_input.mzid", &all_flags(), "mzid14_all");
    // The peak-file metadata layout, every field empty.
    assert!(
        all.text
            .contains("Document ID:        \nDate:               0000-00-00 00:00:00\n")
    );
    assert!(all.text.contains("\nSample:\n  name:             \n"));
    assert!(all.text.contains("\n  detector(s):      \n\n"));
    assert!(all.tsv.contains("\ndate\t0000-00-00 00:00:00\n"));
    // The peak-file statistics arm over no values.
    assert!(all.text.contains("Intensities:\n  num. of values: 0\n"));
    // An mzIdentML input still has no document identifier of its own here.
    assert!(all.experiment_meta.is_none());

    check("FileInfo_15_input.mzid", &bare(), "mzid15");
    check("FileInfo_15_input.mzid", &all_flags(), "mzid15_all");
}

/// D1: `FileInfo.cpp:1336-1341` reads `id_data.proteins[0]` before it has
/// established that a run exists, while the structured block at `:1451` guards
/// the same access. The Release build segmentation-faults.
///
/// The idXML reader refuses the file first, for its own reason; the refusal is
/// what matters, and the identification branch carries the same guard for the
/// inputs that do reach it.
#[cfg(feature = "idxml")]
#[test]
fn identifications_without_a_run_are_refused() {
    let error = FileInfo::new()
        .run(input("a7_id_no_runs.idXML"), &bare())
        .expect_err("a file with no IdentificationRun must be refused");
    assert!(
        matches!(error, Error::Parse { .. } | Error::InvalidValue(_)),
        "{error}"
    );
}

/// D1: `FileInfo.cpp:1354` reads `getHits()[0]` behind a guard that tests for a
/// default-constructed `PeptideIdentification` rather than for an empty hit
/// list. `IdXMLFile::load` always fills the identifier, so the guard never
/// protects the read for a loaded file, whatever the score type is; the Release
/// build segmentation-faults on both files below.
#[cfg(feature = "idxml")]
#[test]
fn identifications_with_a_hitless_identification_are_refused() {
    for name in ["a7_id_empty_hitlist.idXML", "a7_id_empty_hitlist_ok.idXML"] {
        let error = FileInfo::new()
            .run(input(name), &bare())
            .expect_err(&format!("{name} must be refused"));
        let Error::InvalidValue(message) = &error else {
            panic!("{name}: expected InvalidValue, got {error}");
        };
        assert!(
            message.contains("carries no hit") && message.contains("FileInfo.cpp:1347-1354"),
            "{name}: {message}"
        );
    }
}

/// The type the branch is entered with selects the loader, and a mismatch is
/// refused, as `FileHandler::loadIdentifications` refuses one (Release exit 3).
#[cfg(feature = "idxml")]
#[test]
fn identifications_forced_type_must_match_the_content() {
    let error = FileInfo::new()
        .run(
            input("FileInfo_10_input.idXML"),
            &Options {
                forced_type: FileType::MzIdentMl,
                ..Options::default()
            },
        )
        .expect_err("an idXML forced to mzIdentML must be refused");
    let Error::InvalidValue(message) = &error else {
        panic!("expected InvalidValue, got {error}");
    };
    assert!(message.contains("idXML"), "{message}");
}
