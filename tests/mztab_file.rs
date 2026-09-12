// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! MzTab file adapter: `FORMAT/MzTabFile.h`.
//!
//! Every one of the five `START_SECTION`s of `MzTabFile_test.cpp` is ported
//! below; the section each test belongs to is named in a banner comment.
//! Expected values marked "MzTabFile_test.cpp" are transcribed from that class
//! test (tier 3, source review); values marked "fixture" are read out of the
//! unmodified upstream reference files copied into `tests/data/` (also tier 3,
//! because they are input data the upstream suite round-trips rather than
//! output a C++ build produced here). The resource ceilings, the hostile
//! headers, the error-variant choices and every test of a divergence from the
//! source are independently derived (tier 4). No C++ was built or executed.

use openms::Error;
use openms::format::mztab::{
    MzTab, MzTabCell, MzTabDouble, MzTabInteger, MzTabMetaData, MzTabModification,
    MzTabModificationList, MzTabNucleicAcidSectionRow, MzTabOSMSectionRow,
    MzTabOligonucleotideSectionRow, MzTabOptionalColumnEntry, MzTabPSMSectionRow, MzTabParameter,
    MzTabParameterList, MzTabPeptideSectionRow, MzTabProteinSectionRow,
    MzTabSmallMoleculeSectionRow, MzTabString,
};
use openms::format::mztab_file::{
    MzTabFile, SectionLayout, extract_bracket_index, extract_index_pairs_from_brackets,
    optional_column_cells,
};
use openms::system::file::TempDir;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn mzstring(value: &str) -> MzTabString {
    MzTabString::from_text(value)
}

fn parameter(value: &str) -> MzTabParameter {
    MzTabParameter::parse(value).expect("parameter cell parses")
}

fn parameter_list(value: &str) -> MzTabParameterList {
    let mut list = MzTabParameterList::default();
    list.read_cell(value).expect("parameter list cell parses");
    list
}

/// The five reference documents `MzTabFile_test.cpp` round-trips, in its order.
const REFERENCE_FILES: [&str; 5] = [
    "MzTabFile_SILAC.mzTab",
    "MzTabFile_SILAC2.mzTab",
    "MzTabFile_labelfree.mzTab",
    "MzTabFile_iTRAQ.mzTab",
    "MzTabFile_Cytidine.mzTab",
];

// ---------------------------------------------------------------------------
// The upstream comparison, reproduced
// ---------------------------------------------------------------------------

/// `MzTabFile_test.cpp`'s store section compares two files by loading both into
/// a `TextFile`, sorting the lines, removing every space, and handing the result
/// to `TEST_FILE_SIMILAR`, i.e. `FuzzyStringComparator` at its default
/// tolerances (relative 1.0, absolute 0.0). That comparator skips any line that
/// is empty or whitespace only (`FuzzyStringComparator.cpp:843`) and compares
/// numbers numerically wherever a number begins on both sides, so `46` and
/// `46.0` are equal and so are `5035500000` and `5.0355e09`, while any other
/// difference fails.
fn similar_lines(left: &str, right: &str) -> Result<(), String> {
    let prepare = |text: &str| {
        let mut lines: Vec<&str> = text.lines().collect();
        lines.sort_unstable();
        lines
            .into_iter()
            .filter(|line| !line.trim().is_empty())
            .map(|line| line.replace(' ', ""))
            .collect::<Vec<String>>()
    };
    let left_lines = prepare(left);
    let right_lines = prepare(right);
    if left_lines.len() != right_lines.len() {
        return Err(format!(
            "line counts differ: {} vs {}",
            left_lines.len(),
            right_lines.len()
        ));
    }
    for (index, (a, b)) in left_lines.iter().zip(right_lines.iter()).enumerate() {
        if !similar_line(a, b) {
            return Err(format!("sorted line {index} differs:\n  {a}\n  {b}"));
        }
    }
    Ok(())
}

/// One line of the comparison above: walk both strings together, and wherever a
/// numeric literal begins on both sides consume the whole literal from each and
/// compare the values; otherwise compare one character.
fn similar_line(left: &str, right: &str) -> bool {
    let a: Vec<char> = left.chars().collect();
    let b: Vec<char> = right.chars().collect();
    let mut i = 0usize;
    let mut j = 0usize;
    loop {
        match (a.get(i), b.get(j)) {
            (None, None) => return true,
            (Some(x), Some(y)) => match (number_at(&a, i), number_at(&b, j)) {
                (Some((left_value, next_i)), Some((right_value, next_j))) => {
                    if left_value != right_value {
                        return false;
                    }
                    i = next_i;
                    j = next_j;
                }
                _ => {
                    if x != y {
                        return false;
                    }
                    i += 1;
                    j += 1;
                }
            },
            _ => return false,
        }
    }
}

/// The numeric literal starting at `index`, with its value and its end.
fn number_at(chars: &[char], index: usize) -> Option<(f64, usize)> {
    let digits = |from: usize| {
        let mut end = from;
        while chars.get(end).is_some_and(char::is_ascii_digit) {
            end += 1;
        }
        end
    };
    let mut end = index;
    if matches!(chars.get(end), Some('+' | '-')) {
        end += 1;
    }
    let integer_start = end;
    end = digits(end);
    let mut significant = end > integer_start;
    if matches!(chars.get(end), Some('.')) {
        let fraction = digits(end + 1);
        if fraction > end + 1 || significant {
            end = fraction;
            significant = true;
        }
    }
    if !significant {
        return None;
    }
    if matches!(chars.get(end), Some('e' | 'E')) {
        let mut exponent = end + 1;
        if matches!(chars.get(exponent), Some('+' | '-')) {
            exponent += 1;
        }
        let exponent_digits = digits(exponent);
        if exponent_digits > exponent {
            end = exponent_digits;
        }
    }
    let text: String = chars.get(index..end)?.iter().collect();
    text.parse::<f64>().ok().map(|value| (value, end))
}

/// Number of non-empty sections, which is how many blank separator lines a
/// fresh write introduces.
fn section_count(document: &MzTab) -> usize {
    usize::from(!document.protein_data.is_empty())
        + usize::from(!document.peptide_data.is_empty())
        + usize::from(!document.psm_data.is_empty())
        + usize::from(!document.small_molecule_data.is_empty())
        + usize::from(!document.nucleic_acid_data.is_empty())
        + usize::from(!document.oligonucleotide_data.is_empty())
        + usize::from(!document.osm_data.is_empty())
}

/// Write `document`, read it back, and assert the round trip is a fixed point:
/// writing the reload reproduces the same bytes, and reading those reproduces
/// the same document. The reload is returned so the caller can assert against
/// it.
///
/// The first write canonicalises a document built in memory, in two ways the
/// format makes unavoidable. A section's header declares one column set for
/// every row, so a member of an indexed family that only some rows carry
/// becomes a `null` cell — and therefore a null entry — in the rest; and the
/// blank line the writer puts before each section is recorded by the reader as
/// an empty row, whose count this asserts exactly.
fn round_trip(adapter: &MzTabFile, document: &MzTab) -> MzTab {
    let rendered = adapter.write_to_string(document).expect("renders");
    let reloaded = adapter.load_str(&rendered).expect("reloads");
    assert_eq!(
        reloaded.empty_rows.len(),
        section_count(document),
        "one blank separator line per section"
    );
    let again = adapter.write_to_string(&reloaded).expect("renders again");
    assert_eq!(again, rendered, "the round trip is a fixed point");
    assert_eq!(adapter.load_str(&again).expect("reloads again"), reloaded);
    reloaded
}

/// `document` with only the blank separator lines of `reloaded` added, for the
/// cases where nothing else is canonicalised.
fn with_blanks_of(document: &MzTab, reloaded: &MzTab) -> MzTab {
    let mut expected = document.clone();
    expected.empty_rows = reloaded.empty_rows.clone();
    expected
}

// ---------------------------------------------------------------------------
// MzTabFile_test.cpp START_SECTION(MzTabFile())
// ---------------------------------------------------------------------------

#[test]
fn default_construction_disables_every_optional_column() {
    // The section asserts only that `new MzTabFile()` is not the null pointer.
    // The observable equivalent is that the adapter exists and that all sixteen
    // flags its constructor lists are false (MzTabFile.cpp:45-62).
    let adapter = MzTabFile::new();
    assert_eq!(adapter, MzTabFile::default());
    assert!(!adapter.store_protein_reliability);
    assert!(!adapter.store_peptide_reliability);
    assert!(!adapter.store_psm_reliability);
    assert!(!adapter.store_small_molecule_reliability);
    assert!(!adapter.store_protein_uri);
    assert!(!adapter.store_peptide_uri);
    assert!(!adapter.store_psm_uri);
    assert!(!adapter.store_small_molecule_uri);
    assert!(!adapter.store_protein_go_terms);
    assert!(!adapter.store_nucleic_acid_reliability);
    assert!(!adapter.store_oligonucleotide_reliability);
    assert!(!adapter.store_osm_reliability);
    assert!(!adapter.store_nucleic_acid_uri);
    assert!(!adapter.store_oligonucleotide_uri);
    assert!(!adapter.store_osm_uri);
    assert!(!adapter.store_nucleic_acid_go_terms);
    // The seven flags the header keeps protected without a setter are reachable
    // here, which is the whole point of `lossless`.
    let lossless = MzTabFile::lossless();
    assert!(lossless.store_nucleic_acid_go_terms);
    assert!(lossless.store_osm_uri);
    assert_ne!(lossless, adapter);
}

// ---------------------------------------------------------------------------
// MzTabFile_test.cpp START_SECTION(~MzTabFile())
// ---------------------------------------------------------------------------

#[test]
fn dropping_the_adapter_releases_nothing() {
    // The section is `delete ptr;` with no assertion: it exists so the
    // destructor runs under the leak checker. `MzTabFile::~MzTabFile() = default`
    // (MzTabFile.cpp:65) and the Rust adapter is sixteen `bool`s with no owned
    // resource, so there is nothing to release and no `Drop` impl to test. The
    // observable equivalent is that the type is `Copy`, which a type owning a
    // resource could not be, and that dropping one leaves an independent copy
    // untouched.
    assert!(
        !std::mem::needs_drop::<MzTabFile>(),
        "the adapter owns no resource, so there is no destructor to test"
    );
    let adapter = MzTabFile::lossless();
    let copy = adapter;
    assert!(copy.store_protein_go_terms);
    assert_eq!(copy, adapter);
}

// ---------------------------------------------------------------------------
// MzTabFile_test.cpp
// START_SECTION(void load(const std::string& filename, MzTab& mzTab))
// ---------------------------------------------------------------------------

#[test]
fn load_reads_the_silac_reference_document() {
    // The section only calls `MzTabFile().load(...)` on MzTabFile_SILAC.mzTab
    // and asserts nothing. Every literal below is read out of that fixture.
    let document = MzTabFile::new()
        .load(data("MzTabFile_SILAC.mzTab"))
        .expect("SILAC reference document loads");
    let meta = &document.meta_data;
    assert_eq!(meta.mz_tab_version.get(), "1.0.0");
    assert_eq!(meta.mz_tab_mode.get(), "Complete");
    assert_eq!(meta.mz_tab_type.get(), "Quantification");
    assert_eq!(meta.ms_run.len(), 6);
    assert_eq!(meta.assay.len(), 12);
    assert_eq!(meta.study_variable.len(), 2);
    assert_eq!(meta.protein_search_engine_score.len(), 1);
    assert_eq!(
        meta.protein_search_engine_score[&1].to_cell_string(),
        "[MS, MS:1002338, Andromeda:score, ]"
    );
    // fixture: MTD study_variable[1]-assay_refs
    //          assay[1], assay[2], assay[3], assay[4], assay[5], assay[6]
    assert_eq!(meta.study_variable[&1].assay_refs, vec![1, 2, 3, 4, 5, 6]);
    assert_eq!(
        meta.study_variable[&2].assay_refs,
        vec![7, 8, 9, 10, 11, 12]
    );
    assert_eq!(meta.study_variable[&1].sample_refs, vec![1]);
    assert_eq!(meta.study_variable[&1].description.get(), "DB");
    assert_eq!(meta.study_variable[&2].description.get(), "HT");

    assert_eq!(document.protein_data.len(), 57);
    assert_eq!(document.peptide_data.len(), 80);
    assert_eq!(document.psm_data.len(), 946);
    assert!(document.small_molecule_data.is_empty());
    assert_eq!(document.comment_rows.len(), 2);
    assert_eq!(document.empty_rows.len(), 4);

    // The PRH declares one score type over six runs, and no
    // num_peptides_distinct/unique columns at all.
    let protein = &document.protein_data[0];
    assert_eq!(protein.search_engine_score_ms_run.len(), 1);
    assert_eq!(protein.search_engine_score_ms_run[&1].len(), 6);
    assert_eq!(protein.num_psms_ms_run.len(), 6);
    assert!(protein.num_peptides_distinct_ms_run.is_empty());
    assert!(protein.num_peptides_unique_ms_run.is_empty());
    assert_eq!(protein.protein_abundance_assay.len(), 12);
    assert_eq!(protein.protein_abundance_study_variable.len(), 2);
    // Four optional columns, discovered from the header in name order.
    assert_eq!(protein.opt.len(), 4);
    assert_eq!(
        protein.opt[0].name,
        "opt_study_variable[1]_ratio_heavy_to_light"
    );
    // The peptide section declares two of them.
    assert_eq!(document.peptide_data[0].opt.len(), 2);
    // fixture: the last PEP cell of the first peptide row is the text NaN in an
    // optional column, which is a string cell and not the NaN numeric state.
    let last_peptide_option = document.peptide_data[0]
        .opt
        .last()
        .expect("peptide optional column");
    assert_eq!(last_peptide_option.value.get(), "NaN");
}

#[test]
fn load_reads_every_reference_document() {
    // Section counts per file, from the tag census of each fixture.
    let expected = [
        ("MzTabFile_SILAC.mzTab", 57, 80, 946, 0),
        ("MzTabFile_SILAC2.mzTab", 5, 0, 30, 0),
        ("MzTabFile_labelfree.mzTab", 5, 0, 58, 0),
        ("MzTabFile_iTRAQ.mzTab", 5, 0, 36, 0),
        ("MzTabFile_Cytidine.mzTab", 0, 0, 0, 1),
    ];
    for (name, proteins, peptides, psms, small_molecules) in expected {
        let document = MzTabFile::new()
            .load(data(name))
            .unwrap_or_else(|error| panic!("{name} loads: {error}"));
        assert_eq!(document.protein_data.len(), proteins, "{name} PRT rows");
        assert_eq!(document.peptide_data.len(), peptides, "{name} PEP rows");
        assert_eq!(document.psm_data.len(), psms, "{name} PSM rows");
        assert_eq!(
            document.small_molecule_data.len(),
            small_molecules,
            "{name} SML rows"
        );
    }
}

#[test]
fn load_reads_the_cytidine_small_molecule_document() {
    // Cytidine is the only reference file with an SML section, and the only one
    // with CRLF line endings and an `mzTab-version` of `1.0 rc5`.
    let document = MzTabFile::new()
        .load(data("MzTabFile_Cytidine.mzTab"))
        .expect("Cytidine reference document loads");
    let meta = &document.meta_data;
    assert_eq!(meta.mz_tab_version.get(), "1.0 rc5");
    assert_eq!(meta.mz_tab_mode.get(), "Summary");
    assert_eq!(meta.mz_tab_id.get(), "Cytidine");
    assert_eq!(meta.description.get(), "LC-MS/MS Reference Standard");
    assert_eq!(meta.sample_processing[&1].get().len(), 5);
    assert_eq!(meta.instrument[&1].analyzer.len(), 1);
    assert_eq!(meta.software.len(), 2);
    assert_eq!(meta.software[&1].setting[&1].get(), "Peak Picking MS1");
    assert_eq!(meta.contact[&1].email.get(), "beiken@ebi.ac.uk");
    assert_eq!(meta.variable_mod.len(), 2);
    assert_eq!(meta.smallmolecule_search_engine_score.len(), 1);
    assert!(meta.assay.is_empty());
    assert!(meta.study_variable.is_empty());

    let row = &document.small_molecule_data[0];
    assert_eq!(row.identifier.get()[0].get(), "CHEBI:17562");
    assert_eq!(row.chemical_formula.get(), "C9H13N3O5");
    assert_eq!(row.description.get(), "Cytidine");
    assert_eq!(row.exp_mass_to_charge.get().unwrap(), 244.0928);
    // fixture: calc_mass_to_charge is the text `null`, taxid and species too.
    assert!(row.calc_mass_to_charge.is_null());
    assert!(row.taxid.is_null());
    assert!(row.species.is_null());
    assert_eq!(row.charge.get().unwrap(), 1);
    assert_eq!(row.retention_time.get()[0].get().unwrap(), 193.25);
    assert_eq!(row.database.get(), "ChEBI");
    assert_eq!(row.database_version.get(), "109");
    assert!(row.spectra_ref.is_null());
    assert_eq!(row.best_search_engine_score[&1].get().unwrap(), 977.0);
    assert_eq!(row.modifications.get(), "CHEMMOD:2M+H,CHEMMOD:M-C5H8O4");
    assert!(row.smallmolecule_abundance_assay.is_empty());
}

// ---------------------------------------------------------------------------
// MzTabFile_test.cpp
// START_SECTION(void store(const std::string& filename, MzTab& mzTab))
// ---------------------------------------------------------------------------

#[test]
fn store_reproduces_every_reference_document() {
    // The section loads each of the five reference files, stores it, and
    // compares the stored text with the original after sorting the lines and
    // removing spaces. `similar_lines` reproduces that comparison; see its
    // documentation for the tolerances.
    let adapter = MzTabFile::new();
    for name in REFERENCE_FILES {
        let document = MzTabFile::new()
            .load(data(name))
            .unwrap_or_else(|error| panic!("{name} loads: {error}"));
        let stored = adapter
            .write_to_string(&document)
            .unwrap_or_else(|error| panic!("{name} stores: {error}"));
        let original = std::fs::read_to_string(data(name)).expect("fixture readable");
        if let Err(report) = similar_lines(&stored, &original) {
            panic!("{name} did not round-trip: {report}");
        }
    }
}

#[test]
fn store_round_trips_every_reference_document_through_a_file() {
    // The same five documents, but written to and read back from disk, and
    // compared as models rather than as text. Model equality is the stronger
    // statement: it covers the null/NaN/Inf state of every numeric cell, the
    // optional columns present on each row, and the recorded comment and
    // blank-line positions, none of which a sorted space-stripped text
    // comparison can see.
    let adapter = MzTabFile::new();
    let directory = TempDir::new_in(std::env::temp_dir(), false).expect("temp dir");
    for name in REFERENCE_FILES {
        let document = adapter
            .load(data(name))
            .unwrap_or_else(|error| panic!("{name} loads: {error}"));
        let path = directory.path().join(name);
        adapter
            .store(&path, &document)
            .unwrap_or_else(|error| panic!("{name} stores: {error}"));
        let reloaded = adapter
            .load(&path)
            .unwrap_or_else(|error| panic!("{name} reloads: {error}"));
        assert_eq!(reloaded.meta_data, document.meta_data, "{name} metadata");
        assert_eq!(reloaded.protein_data, document.protein_data, "{name} PRT");
        assert_eq!(reloaded.peptide_data, document.peptide_data, "{name} PEP");
        assert_eq!(reloaded.psm_data, document.psm_data, "{name} PSM");
        assert_eq!(
            reloaded.small_molecule_data, document.small_molecule_data,
            "{name} SML"
        );
        assert_eq!(reloaded.empty_rows, document.empty_rows, "{name} blanks");
        assert_eq!(
            reloaded.comment_rows, document.comment_rows,
            "{name} comments"
        );
        assert_eq!(reloaded, document, "{name} document");
    }
}

#[test]
fn store_restores_comments_and_blank_lines_in_place() {
    // SILAC2 opens with two COM lines, carries four more inside its sections,
    // and has two blank lines. Both the count and the recorded line numbers
    // survive a round trip, which the source's doubled blank insertion
    // (MzTabFile.cpp:3328) would not allow.
    let adapter = MzTabFile::new();
    let document = adapter
        .load(data("MzTabFile_SILAC2.mzTab"))
        .expect("SILAC2 loads");
    assert_eq!(document.comment_rows.len(), 6);
    assert_eq!(document.empty_rows.len(), 2);
    assert_eq!(document.comment_rows.keys().copied().next(), Some(0));

    let lines = adapter
        .document_lines(&document)
        .expect("SILAC2 renders to lines");
    let original = std::fs::read_to_string(data("MzTabFile_SILAC2.mzTab")).expect("readable");
    assert_eq!(
        lines.len(),
        original.lines().count(),
        "one output line per input line"
    );
    for (&index, comment) in &document.comment_rows {
        assert_eq!(&lines[index], comment, "comment stays on line {index}");
    }
    for &index in &document.empty_rows {
        assert!(lines[index].is_empty(), "line {index} stays blank");
    }
}

// ---------------------------------------------------------------------------
// MzTabFile_test.cpp
// START_SECTION(generateMzTabPSMSectionRow_(const MzTabPSMSectionRow& row,
//               const vector<std::string>& optional_columns) const)
// ---------------------------------------------------------------------------

#[test]
fn psm_row_fills_requested_optional_columns_and_nulls_the_rest() {
    // Every literal in this test is transcribed from MzTabFile_test.cpp:105-160.
    let mut row = MzTabPSMSectionRow::default();
    row.sequence.read_cell("NDYKAPPQPAPGK").unwrap();
    row.psm_id.read_cell("38").unwrap();
    row.accession.read_cell("IPI:B1").unwrap();
    row.unique.read_cell("1").unwrap();
    row.database.read_cell("null").unwrap();
    row.database_version.read_cell("null").unwrap();
    row.search_engine.read_cell("[, , Percolator, ]").unwrap();
    // The class test writes `row.search_engine_score[0]`, which is
    // std::map::operator[] and creates key zero.
    let mut score = MzTabDouble::default();
    score.read_cell("51.9678841193106").unwrap();
    row.search_engine_score.insert(0, score);

    for (name, value) in [
        ("Percolator_score", "0.359083"),
        ("Percolator_qvalue", "0.00649874"),
        ("Percolator_PEP", "0.0420992"),
        ("search_engine_sequence", "NDYKAPPQPAPGK"),
    ] {
        row.opt
            .push(MzTabOptionalColumnEntry::new(name, mzstring(value)));
    }

    let optional_columns: Vec<String> = [
        "Percolator_score",
        "Percolator_qvalue",
        "EMPTY",
        "Percolator_PEP",
        "search_engine_sequence",
        "AScore_1",
    ]
    .iter()
    .map(|name| (*name).to_owned())
    .collect();
    let layout = SectionLayout {
        search_engine_score: [0].into_iter().collect(),
        optional_columns,
        ..SectionLayout::default()
    };

    let rendered = MzTabFile::new()
        .psm_row(&row, &layout)
        .expect("PSM row renders");
    let cells: Vec<&str> = rendered.split('\t').collect();
    let last = cells.len();
    // TEST_EQUAL(substrings[substrings.size() - 1], "null")
    assert_eq!(cells[last - 1], "null");
    // TEST_EQUAL(substrings[substrings.size() - 2], "NDYKAPPQPAPGK")
    assert_eq!(cells[last - 2], "NDYKAPPQPAPGK");
    // TEST_EQUAL(substrings[substrings.size() - 3], "0.0420992")
    assert_eq!(cells[last - 3], "0.0420992");
    // TEST_EQUAL(substrings[substrings.size() - 4], "null")
    assert_eq!(cells[last - 4], "null");

    // The rest of the row, which the class test does not inspect.
    assert_eq!(cells[0], "PSM");
    assert_eq!(cells[1], "NDYKAPPQPAPGK");
    assert_eq!(cells[2], "38");
    assert_eq!(cells[3], "IPI:B1");
    assert_eq!(cells[4], "1");
    assert_eq!(cells[5], "null");
    assert_eq!(cells[6], "null");
    assert_eq!(cells[7], "[, , Percolator, ]");
    // 51.9678841193106 renders with the source's fifteen fractional digits.
    assert_eq!(cells[8], "51.967884119310597");
    assert_eq!(cells.len(), 25);
    // The header declares the same 24 columns, which is what makes the source's
    // postcondition unreachable here.
    let header = MzTabFile::new()
        .psm_header(&layout)
        .expect("PSM header renders");
    assert_eq!(header.split('\t').count(), cells.len());
}

#[test]
fn optional_column_cells_take_the_first_match_and_null_a_miss() {
    // Source addOptionalColumnsToSectionRow_ breaks on the first matching name,
    // so a duplicated entry name is shadowed by the earlier one.
    let entries = vec![
        MzTabOptionalColumnEntry::new("opt_global_a", mzstring("first")),
        MzTabOptionalColumnEntry::new("opt_global_a", mzstring("second")),
        MzTabOptionalColumnEntry::new("opt_global_b", MzTabString::null()),
    ];
    let names: Vec<String> = ["opt_global_a", "opt_global_b", "opt_global_missing"]
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
    let cells = optional_column_cells(&names, &entries).expect("cells render");
    assert_eq!(cells, vec!["first", "null", "null"]);
}

// ---------------------------------------------------------------------------
// Bracket indices
// ---------------------------------------------------------------------------

#[test]
fn bracket_index_extraction_follows_the_source_and_refuses_zero() {
    assert_eq!(extract_bracket_index("assay[3]", "assay[").unwrap(), 3);
    assert_eq!(
        extract_bracket_index("sample_processing[12]", "sample_processing[").unwrap(),
        12
    );
    // The source trims and accepts a leading '+', because StringUtils::toInt32
    // does.
    assert_eq!(extract_bracket_index("uri[ +7 ]", "uri[").unwrap(), 7);
    // Native strictness: the source casts a signed Int to Size, so these become
    // key 0 and key usize::MAX and are written back as column names the format
    // does not define.
    assert!(matches!(
        extract_bracket_index("assay[0]", "assay["),
        Err(Error::Parse { .. })
    ));
    assert!(matches!(
        extract_bracket_index("assay[-1]", "assay["),
        Err(Error::Parse { .. })
    ));
    // StringUtils::toInt32 throws on non-numeric text; so does this.
    assert!(matches!(
        extract_bracket_index("colunit", "colunit["),
        Err(Error::Parse { .. })
    ));
    assert!(matches!(
        extract_bracket_index("assay[2000000]", "assay["),
        Err(Error::Parse { .. })
    ));
}

#[test]
fn index_pair_extraction_needs_two_bracketed_groups() {
    assert_eq!(
        extract_index_pairs_from_brackets("search_engine_score[1]_ms_run[2]").unwrap(),
        (1, 2)
    );
    assert_eq!(
        extract_index_pairs_from_brackets("search_engine_score[11]_ms_run[203]").unwrap(),
        (11, 203)
    );
    // The source's regex leaves the missing member at zero; this refuses.
    assert!(matches!(
        extract_index_pairs_from_brackets("search_engine_score[1]"),
        Err(Error::Parse { .. })
    ));
    assert!(matches!(
        extract_index_pairs_from_brackets("search_engine_score[]_ms_run[]"),
        Err(Error::Parse { .. })
    ));
}

// ---------------------------------------------------------------------------
// Optional column discovery
// ---------------------------------------------------------------------------

const MINIMAL_METADATA: &str = "MTD\tmzTab-version\t1.0.0\n\
     MTD\tmzTab-mode\tSummary\n\
     MTD\tmzTab-type\tIdentification\n\
     MTD\tdescription\td\n";

#[test]
fn psm_optional_columns_are_discovered_from_the_header() {
    // Source MzTabFile.cpp:1287 tests `cells[i] == "opt_"` instead of the `opt_`
    // prefix the other three sections use, so no PSM optional column is ever
    // discovered and every `opt_` cell of a PSM section is dropped on load.
    let text = format!(
        "{MINIMAL_METADATA}\n\
         PSH\tsequence\tPSM_ID\topt_global_score\topt_global_note\n\
         PSM\tPEPTIDE\t1\t0.5\thello\n"
    );
    let document = MzTabFile::new().load_str(&text).expect("document loads");
    let row = &document.psm_data[0];
    assert_eq!(row.sequence.get(), "PEPTIDE");
    assert_eq!(row.psm_id.get().unwrap(), 1);
    assert_eq!(row.opt.len(), 2);
    assert_eq!(row.opt[0].name, "opt_global_note");
    assert_eq!(row.opt[0].value.get(), "hello");
    assert_eq!(row.opt[1].name, "opt_global_score");
    assert_eq!(row.opt[1].value.get(), "0.5");
}

#[test]
fn optional_columns_tolerate_any_header_order_and_unknown_names() {
    // The required columns may appear in any order, an unknown column name is
    // ignored, and a required column the header omits leaves its field null.
    let text = format!(
        "{MINIMAL_METADATA}\n\
         PSH\topt_global_b\tcharge\tsequence\tunknown_column\topt_global_a\tPSM_ID\n\
         PSM\tsecond\t2\tPEPTIDE\tignored\tfirst\t7\n"
    );
    let document = MzTabFile::new().load_str(&text).expect("document loads");
    let row = &document.psm_data[0];
    assert_eq!(row.sequence.get(), "PEPTIDE");
    assert_eq!(row.charge.get().unwrap(), 2);
    assert_eq!(row.psm_id.get().unwrap(), 7);
    // Columns the header never declared stay at their defaults rather than
    // picking up cells[0], which is what the source's zero sentinel does.
    assert!(row.accession.is_null());
    assert!(row.exp_mass_to_charge.is_null());
    assert!(row.retention_time.is_null());
    // Optional columns come back in name order, as the source's std::map does.
    assert_eq!(row.opt.len(), 2);
    assert_eq!(row.opt[0].name, "opt_global_a");
    assert_eq!(row.opt[0].value.get(), "first");
    assert_eq!(row.opt[1].name, "opt_global_b");
    assert_eq!(row.opt[1].value.get(), "second");
}

#[test]
fn a_column_absent_from_one_row_is_written_as_null_not_dropped() {
    // Two rows declaring different optional columns produce one header with the
    // union of both, and the row missing a column gets a null cell. Reading
    // that back gives both rows all the columns, which is the only shape a
    // tab-separated section can express.
    let mut first = MzTabPSMSectionRow::default();
    first.sequence.set("AAA");
    first.psm_id.set(1);
    first
        .opt
        .push(MzTabOptionalColumnEntry::new("opt_global_a", mzstring("1")));
    let mut second = MzTabPSMSectionRow::default();
    second.sequence.set("BBB");
    second.psm_id.set(2);
    second
        .opt
        .push(MzTabOptionalColumnEntry::new("opt_global_b", mzstring("2")));

    let mut document = MzTab::default();
    document.meta_data.mz_tab_mode.set("Summary");
    document.meta_data.mz_tab_type.set("Identification");
    document.psm_data = vec![first, second];

    let adapter = MzTabFile::new();
    let rendered = adapter.write_to_string(&document).expect("renders");
    assert!(rendered.contains("opt_global_a\topt_global_b"));
    let reloaded = adapter.load_str(&rendered).expect("reloads");
    assert_eq!(reloaded.psm_data[0].opt.len(), 2);
    assert_eq!(reloaded.psm_data[0].opt[1].name, "opt_global_b");
    assert!(reloaded.psm_data[0].opt[1].value.is_null());
    assert_eq!(reloaded.psm_data[1].opt[0].name, "opt_global_a");
    assert!(reloaded.psm_data[1].opt[0].value.is_null());
    // A second write is now stable, which is what makes the round trip a fixed
    // point rather than only reversible once.
    round_trip(&adapter, &reloaded);
}

// ---------------------------------------------------------------------------
// The null / NaN / Inf distinction
// ---------------------------------------------------------------------------

#[test]
fn null_nan_and_inf_survive_a_round_trip_in_every_numeric_column() {
    let mut row = MzTabPSMSectionRow::default();
    row.sequence.set("PEPTIDE");
    row.psm_id.set(1);
    row.exp_mass_to_charge = MzTabDouble::nan();
    row.calc_mass_to_charge = MzTabDouble::inf();
    row.charge = MzTabInteger::nan();
    row.search_engine_score.insert(1, MzTabDouble::null());
    row.search_engine_score.insert(2, MzTabDouble::nan());
    row.search_engine_score.insert(3, MzTabDouble::inf());
    row.search_engine_score.insert(4, MzTabDouble::new(-0.5));
    row.opt.push(MzTabOptionalColumnEntry::new(
        "opt_global_absent",
        MzTabString::null(),
    ));

    let mut document = MzTab::default();
    document.meta_data.mz_tab_mode.set("Summary");
    document.meta_data.mz_tab_type.set("Identification");
    document.psm_data = vec![row];

    let adapter = MzTabFile::new();
    let rendered = adapter.write_to_string(&document).expect("renders");
    assert!(rendered.contains("\tNaN\t"));
    assert!(rendered.contains("\tInf\t"));
    assert!(rendered.contains("\tnull\t"));
    let reloaded = round_trip(&adapter, &document);
    let back = &reloaded.psm_data[0];
    assert!(back.exp_mass_to_charge.is_nan());
    assert!(back.calc_mass_to_charge.is_inf());
    assert!(back.charge.is_nan());
    assert!(back.search_engine_score[&1].is_null());
    assert!(back.search_engine_score[&2].is_nan());
    assert!(back.search_engine_score[&3].is_inf());
    assert_eq!(back.search_engine_score[&4].get().unwrap(), -0.5);
    assert!(back.opt[0].value.is_null());
    // Nothing here is canonicalised: every score index the layout declares is
    // carried by the one row, so the only difference is the separator line.
    assert_eq!(reloaded, with_blanks_of(&document, &reloaded));
}

// ---------------------------------------------------------------------------
// Every section type in one document
// ---------------------------------------------------------------------------

fn every_section_document() -> MzTab {
    let mut document = MzTab::default();
    let meta = &mut document.meta_data;
    meta.mz_tab_mode.set("Complete");
    meta.mz_tab_type.set("Quantification");
    meta.mz_tab_id.set("EVERY-SECTION");
    meta.title.set("every section type");
    meta.description.set("one row per section");
    meta.sample_processing
        .insert(1, parameter_list("[MS, MS:1000544, Conversion to mzML, ]"));
    meta.protein_search_engine_score
        .insert(1, parameter("[MS, MS:1001171, Mascot:score, ]"));
    meta.peptide_search_engine_score
        .insert(1, parameter("[MS, MS:1001171, Mascot:score, ]"));
    meta.psm_search_engine_score
        .insert(1, parameter("[MS, MS:1001171, Mascot:score, ]"));
    meta.psm_search_engine_score
        .insert(2, parameter("[MS, MS:1002355, PSM-level FDR, ]"));
    meta.smallmolecule_search_engine_score.insert(
        1,
        parameter("[MS, MS:1001153, search engine specific score, ]"),
    );
    meta.nucleic_acid_search_engine_score
        .insert(1, parameter("[, , NASearch, ]"));
    meta.oligonucleotide_search_engine_score
        .insert(1, parameter("[, , OligoSearch, ]"));
    meta.osm_search_engine_score
        .insert(1, parameter("[, , OsmSearch, ]"));
    meta.instrument.entry(1).or_default().name = parameter("[MS, MS:1000483, model, LTQ]");
    meta.instrument.entry(1).or_default().analyzer.insert(
        1,
        parameter("[MS, MS:1000443, Mass Analyzer Type, Orbitrap]"),
    );
    meta.software.entry(1).or_default().software = parameter("[MS, MS:1001583, MaxQuant, ]");
    meta.software
        .entry(1)
        .or_default()
        .setting
        .insert(1, mzstring("Peak Picking MS1"));
    meta.false_discovery_rate = parameter_list("[MS, MS:1001364, pep:global FDR, 0.01]");
    meta.publication.insert(1, mzstring("pubmed:12345"));
    meta.contact.entry(1).or_default().name = mzstring("A Scientist");
    meta.contact.entry(1).or_default().email = mzstring("a@example.org");
    meta.uri.insert(1, mzstring("http://example.org/study"));
    meta.fixed_mod.entry(1).or_default().modification =
        parameter("[UNIMOD, UNIMOD:4, Carbamidomethyl, ]");
    meta.fixed_mod.entry(1).or_default().site = mzstring("C");
    meta.fixed_mod.entry(1).or_default().position = mzstring("Anywhere");
    meta.variable_mod.entry(1).or_default().modification =
        parameter("[UNIMOD, UNIMOD:35, Oxidation, ]");
    meta.quantification_method = parameter("[MS, MS:1001835, SILAC, ]");
    meta.protein_quantification_unit = parameter("[PRIDE, PRIDE:0000393, Relative, ]");
    meta.peptide_quantification_unit = parameter("[PRIDE, PRIDE:0000393, Relative, ]");
    meta.small_molecule_quantification_unit = parameter("[PRIDE, PRIDE:0000393, Relative, ]");
    for run in 1..=2usize {
        let entry = meta.ms_run.entry(run).or_default();
        entry.location.set(&format!("file:///run{run}.mzML"));
        entry.format = parameter("[MS, MS:1000584, mzML file, ]");
        entry.id_format = parameter("[MS, MS:1000768, Thermo nativeID format, ]");
        entry.fragmentation_method = parameter_list("[MS, MS:1000133, CID, ]");
    }
    meta.custom.insert(1, parameter("[, , custom note, x]"));
    let sample = meta.sample.entry(1).or_default();
    sample.description.set("one sample");
    sample
        .species
        .insert(1, parameter("[NEWT, 9606, Homo sapiens, ]"));
    sample
        .tissue
        .insert(1, parameter("[BTO, BTO:0000759, liver, ]"));
    sample
        .cell_type
        .insert(1, parameter("[CL, CL:0000182, hepatocyte, ]"));
    sample
        .disease
        .insert(1, parameter("[DOID, DOID:684, carcinoma, ]"));
    sample.custom.insert(1, parameter("[, , sample note, y]"));
    for assay in 1..=2usize {
        let entry = meta.assay.entry(assay).or_default();
        entry.quantification_reagent = parameter("[PRIDE, PRIDE:0000326, SILAC light, ]");
        entry.sample_ref.set("sample[1]");
        entry.ms_run_ref = vec![i32::try_from(assay).unwrap()];
        let modification = entry.quantification_mod.entry(1).or_default();
        modification.modification = parameter("[UNIMOD, UNIMOD:259, Label, ]");
        modification.site = mzstring("L");
        modification.position = mzstring("Anywhere");
    }
    let variable = meta.study_variable.entry(1).or_default();
    variable.assay_refs = vec![1, 2];
    variable.sample_refs = vec![1];
    variable.description.set("only variable");
    let cv = meta.cv.entry(1).or_default();
    cv.label.set("MS");
    cv.full_name.set("PSI-MS controlled vocabulary");
    cv.version.set("4.1.0");
    cv.url.set("https://example.org/psi-ms.obo");
    meta.colunit_protein
        .push("protein_coverage=[UO, UO:0000187, percent, ]".to_owned());
    meta.colunit_peptide
        .push("retention_time=[UO, UO:0000031, minute, ]".to_owned());
    meta.colunit_psm
        .push("retention_time=[UO, UO:0000010, second, ]".to_owned());
    meta.colunit_small_molecule
        .push("retention_time=[UO, UO:0000031, minute, ]".to_owned());

    let mut protein = MzTabProteinSectionRow::default();
    protein.accession.set("P63017");
    protein.description.set("Heat shock cognate 71 kDa protein");
    protein.taxid.set(10090);
    protein.species.set("Mus musculus");
    protein.database.set("UniProtKB");
    protein.database_version.set("2013_08");
    protein.search_engine = parameter_list("[MS, MS:1001207, Mascot, ]");
    protein
        .best_search_engine_score
        .insert(1, MzTabDouble::new(46.0));
    protein
        .search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(1, MzTabDouble::new(46.0));
    protein
        .search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(2, MzTabDouble::null());
    protein.reliability.set(1);
    protein.num_psms_ms_run.insert(1, 3.into());
    protein.num_peptides_distinct_ms_run.insert(1, 2.into());
    protein.num_peptides_unique_ms_run.insert(1, 1.into());
    protein
        .ambiguity_members
        .set(vec![mzstring("Q340U4"), mzstring("Q5K0U2")]);
    let mut modification = MzTabModification::null();
    modification.set_modification_identifier(mzstring("UNIMOD:35"));
    modification.set_positions_and_parameters(vec![(12, MzTabParameter::null())]);
    protein.modifications.set(vec![modification]);
    protein.uri.set("http://example.org/P63017");
    protein.go_terms.set(vec![mzstring("GO:0005524")]);
    protein.coverage.set(0.34);
    protein
        .protein_abundance_assay
        .insert(1, MzTabDouble::new(44.5));
    protein
        .protein_abundance_study_variable
        .insert(1, MzTabDouble::new(36.25));
    protein
        .protein_abundance_stdev_study_variable
        .insert(1, MzTabDouble::new(9.25));
    protein
        .protein_abundance_std_error_study_variable
        .insert(1, MzTabDouble::nan());
    protein.opt.push(MzTabOptionalColumnEntry::new(
        "opt_global_ratio",
        mzstring("1.5"),
    ));
    document.protein_data.push(protein);

    let mut peptide = MzTabPeptideSectionRow::default();
    peptide.sequence.set("QTQTFTTYSDNQPGVL");
    peptide.accession.set("P63017");
    peptide.unique.set(true);
    peptide.database.set("UniProtKB");
    peptide.database_version.set("2013_08");
    peptide.search_engine = parameter_list("[MS, MS:1001207, Mascot, ]");
    peptide
        .best_search_engine_score
        .insert(1, MzTabDouble::new(46.0));
    peptide
        .search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(1, MzTabDouble::new(46.0));
    peptide.reliability.set(2);
    peptide.retention_time.set(vec![MzTabDouble::new(1336.5)]);
    peptide
        .retention_time_window
        .set(vec![MzTabDouble::new(1330.0), MzTabDouble::new(1340.0)]);
    peptide.charge.set(3);
    peptide.mass_to_charge.set(600.25);
    peptide.uri.set("http://example.org/peptide/1");
    peptide
        .spectra_ref
        .read_cell("ms_run[1]:scan=1296")
        .unwrap();
    peptide
        .peptide_abundance_assay
        .insert(1, MzTabDouble::new(85.0));
    // Only the value map is filled: the source's lock-step triple would then
    // write no study-variable columns at all.
    peptide
        .peptide_abundance_study_variable
        .insert(1, MzTabDouble::new(85.0));
    peptide.opt.push(MzTabOptionalColumnEntry::new(
        "opt_global_ratio",
        mzstring("2.5"),
    ));
    document.peptide_data.push(peptide);

    let mut psm = MzTabPSMSectionRow::default();
    psm.sequence.set("QTQTFTTYSDNQPGVL");
    psm.psm_id.set(1);
    psm.accession.set("P63017");
    psm.unique.set(true);
    psm.database.set("UniProtKB");
    psm.database_version.set("2013_08");
    psm.search_engine = parameter_list("[MS, MS:1001207, Mascot, ]");
    psm.search_engine_score.insert(1, MzTabDouble::new(46.0));
    psm.search_engine_score.insert(2, MzTabDouble::new(0.01));
    psm.reliability.set(3);
    psm.retention_time.set(vec![MzTabDouble::new(1336.5)]);
    psm.charge.set(3);
    psm.exp_mass_to_charge.set(600.25);
    psm.calc_mass_to_charge = MzTabDouble::nan();
    psm.uri.set("http://example.org/psm/1");
    psm.spectra_ref.read_cell("ms_run[1]:scan=1296").unwrap();
    psm.pre.set("K");
    psm.post.set("I");
    psm.start.set("424");
    psm.end.set("439");
    document.psm_data.push(psm);

    let mut small_molecule = MzTabSmallMoleculeSectionRow::default();
    small_molecule.identifier.set(vec![mzstring("CHEBI:17562")]);
    small_molecule.chemical_formula.set("C9H13N3O5");
    small_molecule
        .smiles
        .set("Nc1ccn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)n1");
    small_molecule.inchi_key.set("UHDGCWIWMRVCDJ-XVFCMESISA-N");
    small_molecule.description.set("Cytidine");
    small_molecule.exp_mass_to_charge.set(244.0);
    small_molecule.calc_mass_to_charge = MzTabDouble::null();
    small_molecule.charge.set(1);
    small_molecule
        .retention_time
        .set(vec![MzTabDouble::new(193.25)]);
    small_molecule.taxid.set(9606);
    small_molecule.species.set("Homo sapiens");
    small_molecule.database.set("ChEBI");
    small_molecule.database_version.set("109");
    small_molecule.reliability.set(1);
    small_molecule.uri.set("http://example.org/chebi/17562");
    small_molecule
        .spectra_ref
        .read_cell("ms_run[2]:scan=977")
        .unwrap();
    small_molecule.search_engine = parameter_list("[MS, MS:1001083, ms-ms search, MassBank]");
    small_molecule
        .best_search_engine_score
        .insert(1, MzTabDouble::new(977.0));
    small_molecule
        .search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(2, MzTabDouble::new(977.0));
    small_molecule.modifications.set("CHEMMOD:2M+H");
    // The source's SML row writer emits no assay cells at all, so this document
    // cannot be written by it.
    small_molecule
        .smallmolecule_abundance_assay
        .insert(1, MzTabDouble::new(12.5));
    small_molecule
        .smallmolecule_abundance_study_variable
        .insert(1, MzTabDouble::new(12.5));
    small_molecule
        .smallmolecule_abundance_stdev_study_variable
        .insert(1, MzTabDouble::new(0.5));
    small_molecule
        .smallmolecule_abundance_std_error_study_variable
        .insert(1, MzTabDouble::inf());
    small_molecule.opt.push(MzTabOptionalColumnEntry::new(
        "opt_global_adduct",
        mzstring("M+H"),
    ));
    document.small_molecule_data.push(small_molecule);

    let mut nucleic_acid = MzTabNucleicAcidSectionRow::default();
    nucleic_acid.accession.set("NA-1");
    nucleic_acid.description.set("a nucleic acid");
    nucleic_acid.taxid.set(9606);
    nucleic_acid.species.set("Homo sapiens");
    nucleic_acid.database.set("RNAcentral");
    nucleic_acid.database_version.set("14");
    nucleic_acid.search_engine = parameter_list("[, , NASearch, ]");
    nucleic_acid
        .best_search_engine_score
        .insert(1, MzTabDouble::new(12.5));
    nucleic_acid
        .search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(1, MzTabDouble::new(12.5));
    nucleic_acid.reliability.set(1);
    // Keyed from one, which the source's zero-based header cannot express.
    nucleic_acid.num_osms_ms_run.insert(1, 4.into());
    nucleic_acid.num_oligos_distinct_ms_run.insert(1, 3.into());
    nucleic_acid.num_oligos_unique_ms_run.insert(1, 2.into());
    nucleic_acid.ambiguity_members.set(vec![mzstring("NA-2")]);
    nucleic_acid.modifications = MzTabModificationList::null();
    nucleic_acid.uri.set("http://example.org/na/1");
    nucleic_acid.go_terms.set(vec![mzstring("GO:0003723")]);
    nucleic_acid.coverage.set(0.5);
    nucleic_acid.opt.push(MzTabOptionalColumnEntry::new(
        "opt_global_na",
        mzstring("yes"),
    ));
    document.nucleic_acid_data.push(nucleic_acid);

    let mut oligonucleotide = MzTabOligonucleotideSectionRow::default();
    oligonucleotide.sequence.set("AUGCAU");
    oligonucleotide.accession.set("NA-1");
    oligonucleotide.unique.set(false);
    oligonucleotide.search_engine = parameter_list("[, , OligoSearch, ]");
    oligonucleotide
        .best_search_engine_score
        .insert(1, MzTabDouble::new(8.5));
    oligonucleotide
        .search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(1, MzTabDouble::new(8.5));
    oligonucleotide.reliability.set(2);
    oligonucleotide
        .retention_time
        .set(vec![MzTabDouble::new(64.0)]);
    oligonucleotide
        .retention_time_window
        .set(vec![MzTabDouble::new(60.0), MzTabDouble::new(68.0)]);
    oligonucleotide.uri.set("http://example.org/oli/1");
    oligonucleotide.pre.set("A");
    oligonucleotide.post.set("U");
    oligonucleotide.start.set(3);
    oligonucleotide.end.set(8);
    oligonucleotide.opt.push(MzTabOptionalColumnEntry::new(
        "opt_global_oli",
        mzstring("no"),
    ));
    document.oligonucleotide_data.push(oligonucleotide);

    let mut osm = MzTabOSMSectionRow::default();
    osm.sequence.set("AUGCAU");
    osm.search_engine = parameter_list("[, , OsmSearch, ]");
    osm.search_engine_score.insert(1, MzTabDouble::new(8.5));
    osm.reliability.set(3);
    osm.retention_time.set(vec![MzTabDouble::new(64.0)]);
    osm.charge.set(2);
    osm.exp_mass_to_charge.set(500.5);
    osm.calc_mass_to_charge = MzTabDouble::inf();
    osm.uri.set("http://example.org/osm/1");
    osm.spectra_ref.read_cell("ms_run[2]:scan=42").unwrap();
    osm.opt.push(MzTabOptionalColumnEntry::new(
        "opt_global_osm",
        mzstring("maybe"),
    ));
    document.osm_data.push(osm);

    document
}

#[test]
fn every_section_type_and_metadata_key_round_trips() {
    let document = every_section_document();
    let adapter = MzTabFile::lossless();
    let rendered = adapter.write_to_string(&document).expect("renders");
    for tag in [
        "PRH\t", "PRT\t", "PEH\t", "PEP\t", "PSH\t", "PSM\t", "SMH\t", "SML\t", "NUH\t", "NUC\t",
        "OLH\t", "OLI\t", "OSH\t", "OSM\t",
    ] {
        assert!(rendered.contains(tag), "output is missing a {tag} line");
    }
    let reloaded = round_trip(&adapter, &document);
    // The metadata is reproduced key for key: every one of the fifty-odd keys
    // the writer knows is exercised by this document.
    assert_eq!(reloaded.meta_data, document.meta_data, "metadata");
    assert_eq!(reloaded.protein_data.len(), 1, "PRT");
    assert_eq!(reloaded.peptide_data.len(), 1, "PEP");
    assert_eq!(reloaded.psm_data, document.psm_data, "PSM");
    assert_eq!(reloaded.nucleic_acid_data.len(), 1, "NUC");
    assert_eq!(reloaded.oligonucleotide_data.len(), 1, "OLI");
    assert_eq!(reloaded.osm_data, document.osm_data, "OSM");
    assert_eq!(reloaded.empty_rows.len(), 7, "one blank line per section");
    // Every section's rows come back with their own fields; the exact
    // differences a write introduces are asserted one by one in
    // `the_canonicalisation_a_write_performs_is_exactly_the_declared_columns`.
    assert_eq!(reloaded.protein_data[0].accession.get(), "P63017");
    assert_eq!(reloaded.protein_data[0].coverage.get().unwrap(), 0.34);
    assert_eq!(
        reloaded.protein_data[0]
            .modifications
            .to_cell_string()
            .unwrap(),
        "12-UNIMOD:35"
    );
    assert_eq!(reloaded.peptide_data[0].sequence.get(), "QTQTFTTYSDNQPGVL");
    assert_eq!(
        reloaded.small_molecule_data[0].identifier.get()[0].get(),
        "CHEBI:17562"
    );
    assert!(
        reloaded.small_molecule_data[0].smallmolecule_abundance_std_error_study_variable[&1]
            .is_inf()
    );
    assert_eq!(
        reloaded.nucleic_acid_data[0].num_osms_ms_run[&1]
            .get()
            .unwrap(),
        4
    );
    assert_eq!(reloaded.oligonucleotide_data[0].end.get().unwrap(), 8);
    assert!(reloaded.osm_data[0].calc_mass_to_charge.is_inf());
}

#[test]
fn the_canonicalisation_a_write_performs_is_exactly_the_declared_columns() {
    // The document of `every_section_document` declares `assay[2]`,
    // `study_variable[1]` and, in Complete mode, `ms_run[1]` and `ms_run[2]`,
    // and some of its rows carry no value for some of those. The header must
    // declare a column for each, so the reload sees a null entry where the
    // original had no entry at all. Nothing else changes, which this asserts by
    // listing every addition.
    let document = every_section_document();
    let adapter = MzTabFile::lossless();
    let reloaded = round_trip(&adapter, &document);

    let mut expected = with_blanks_of(&document, &reloaded);
    let null = MzTabDouble::null();
    expected.protein_data[0]
        .protein_abundance_assay
        .insert(2, null);
    expected.peptide_data[0]
        .peptide_abundance_assay
        .insert(2, null);
    expected.peptide_data[0]
        .peptide_abundance_stdev_study_variable
        .insert(1, null);
    expected.peptide_data[0]
        .peptide_abundance_std_error_study_variable
        .insert(1, null);
    expected.peptide_data[0]
        .search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(2, null);
    expected.small_molecule_data[0]
        .smallmolecule_abundance_assay
        .insert(2, null);
    expected.small_molecule_data[0]
        .search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(1, null);
    expected.nucleic_acid_data[0]
        .search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(2, null);
    expected.oligonucleotide_data[0]
        .search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(2, null);
    assert_eq!(reloaded, expected);
}

#[test]
fn small_molecule_assay_cells_are_written() {
    // Source MzTabFile.cpp:2548: the SML row writer jumps from `modifications`
    // to the study-variable triple and never emits
    // `smallmolecule_abundance_assay[n]`, so any document with an assay trips
    // the header/content postcondition.
    let document = every_section_document();
    let adapter = MzTabFile::lossless();
    let layout =
        SectionLayout::for_small_molecule(&document.small_molecule_data, &document.meta_data)
            .expect("layout");
    assert_eq!(layout.abundance_assay.len(), 2);
    let header = adapter.small_molecule_header(&layout).expect("header");
    let row = adapter
        .small_molecule_row(&document.small_molecule_data[0], &layout)
        .expect("row");
    assert!(header.contains("smallmolecule_abundance_assay[1]"));
    assert!(header.contains("smallmolecule_abundance_assay[2]"));
    assert_eq!(header.split('\t').count(), row.split('\t').count());
    // Assay 2 has no value in the row, so it is a null cell rather than a
    // missing column.
    let names: Vec<&str> = header.split('\t').collect();
    let cells: Vec<&str> = row.split('\t').collect();
    let first = names
        .iter()
        .position(|name| *name == "smallmolecule_abundance_assay[1]")
        .expect("assay column");
    assert_eq!(cells[first], "12.5");
    assert_eq!(cells[first + 1], "null");
}

#[test]
fn peptide_study_variable_triple_survives_a_missing_stdev_map() {
    // Source MzTabFile.cpp:2367 advances the value, stdev and std_error
    // iterators together and stops at the first exhausted one.
    let document = every_section_document();
    let peptide = &document.peptide_data[0];
    assert!(peptide.peptide_abundance_stdev_study_variable.is_empty());
    let adapter = MzTabFile::lossless();
    let layout =
        SectionLayout::for_peptide(&document.peptide_data, &document.meta_data).expect("layout");
    let header = adapter.peptide_header(&layout).expect("header");
    let row = adapter.peptide_row(peptide, &layout).expect("row");
    assert!(header.contains("peptide_abundance_stdev_study_variable[1]"));
    assert_eq!(header.split('\t').count(), row.split('\t').count());
    let names: Vec<&str> = header.split('\t').collect();
    let cells: Vec<&str> = row.split('\t').collect();
    let value = names
        .iter()
        .position(|name| *name == "peptide_abundance_study_variable[1]")
        .expect("study variable column");
    assert_eq!(cells[value], "85.0");
    assert_eq!(cells[value + 1], "null");
    assert_eq!(cells[value + 2], "null");
}

#[test]
fn nucleic_acid_count_columns_are_numbered_from_the_data() {
    // Source MzTabFile.cpp:2603 writes num_osms_ms_run[0] while its row writer
    // emits the value keyed one.
    let document = every_section_document();
    let adapter = MzTabFile::lossless();
    let layout = SectionLayout::for_nucleic_acid(&document.nucleic_acid_data, &document.meta_data)
        .expect("layout");
    let header = adapter.nucleic_acid_header(&layout).expect("header");
    assert!(header.contains("num_osms_ms_run[1]"));
    assert!(!header.contains("num_osms_ms_run[0]"));
    assert!(header.contains("num_oligos_distinct_ms_run[1]"));
    assert!(header.contains("num_oligos_unique_ms_run[1]"));
    assert!(header.contains("sequence_coverage"));
    let row = adapter
        .nucleic_acid_row(&document.nucleic_acid_data[0], &layout)
        .expect("row");
    assert_eq!(header.split('\t').count(), row.split('\t').count());
}

#[test]
fn uri_columns_are_read_into_their_own_section() {
    // Source MzTabFile.cpp:1089 and :1265 assign the peptide and PSM `uri`
    // column index to `protein_uri_index`, so a peptide or PSM `uri` cell is
    // never read and, worse, a later PRT row reads its protein `uri` from the
    // peptide section's column number.
    let document = every_section_document();
    let adapter = MzTabFile::lossless();
    let rendered = adapter.write_to_string(&document).expect("renders");
    let reloaded = adapter.load_str(&rendered).expect("reloads");
    assert_eq!(
        reloaded.protein_data[0].uri.get(),
        "http://example.org/P63017"
    );
    assert_eq!(
        reloaded.peptide_data[0].uri.get(),
        "http://example.org/peptide/1"
    );
    assert_eq!(reloaded.psm_data[0].uri.get(), "http://example.org/psm/1");
    assert_eq!(
        reloaded.small_molecule_data[0].uri.get(),
        "http://example.org/chebi/17562"
    );
    assert_eq!(
        reloaded.nucleic_acid_data[0].uri.get(),
        "http://example.org/na/1"
    );
    assert_eq!(
        reloaded.oligonucleotide_data[0].uri.get(),
        "http://example.org/oli/1"
    );
    assert_eq!(reloaded.osm_data[0].uri.get(), "http://example.org/osm/1");
}

#[test]
fn clearing_an_optional_column_flag_drops_its_cells() {
    // Documented loss: the `reliability`, `uri` and `go_terms` columns are
    // optional in the format, so the default adapter cannot round-trip a row
    // that carries them.
    let document = every_section_document();
    let default_adapter = MzTabFile::new();
    let rendered = default_adapter.write_to_string(&document).expect("renders");
    assert!(!rendered.contains("\treliability\t"));
    assert!(!rendered.contains("\tgo_terms\t"));
    let reloaded = default_adapter.load_str(&rendered).expect("reloads");
    assert!(reloaded.protein_data[0].reliability.is_null());
    assert!(reloaded.protein_data[0].uri.is_null());
    assert!(reloaded.protein_data[0].go_terms.is_null());
    assert_ne!(reloaded, document);
    // Everything else still matches, so the loss is confined to those columns.
    assert_eq!(reloaded.meta_data, document.meta_data);
    // assay[2] is declared by the metadata and carried by no row, so the
    // header declares its column and the reload sees a null entry for it: a
    // tab-separated section cannot distinguish "absent" from "null".
    assert_eq!(
        reloaded.protein_data[0].protein_abundance_assay[&1]
            .get()
            .unwrap(),
        44.5
    );
    assert!(reloaded.protein_data[0].protein_abundance_assay[&2].is_null());
    assert_eq!(reloaded.protein_data[0].protein_abundance_assay.len(), 2);
}

// ---------------------------------------------------------------------------
// Multiple score types over multiple MS runs
// ---------------------------------------------------------------------------

#[test]
fn two_score_types_over_two_runs_keep_header_and_values_aligned() {
    // Source MzTabFile.cpp:2037 writes the header runs-outer, scores-inner and
    // :2119 writes the values scores-outer, runs-inner. The counts agree, so the
    // postcondition stays silent and every value after the first lands under the
    // wrong header. This port uses one order on both sides.
    let mut row = MzTabProteinSectionRow::default();
    row.accession.set("P1");
    for score in 1..=2usize {
        for run in 1..=2usize {
            let value = MzTabDouble::new(f64::from(
                u32::try_from(score * 10 + run).expect("small value"),
            ));
            row.search_engine_score_ms_run
                .entry(score)
                .or_default()
                .insert(run, value);
        }
    }
    let mut document = MzTab::default();
    document.meta_data.mz_tab_mode.set("Complete");
    document.meta_data.mz_tab_type.set("Identification");
    document
        .meta_data
        .ms_run
        .entry(1)
        .or_default()
        .location
        .set("a.mzML");
    document
        .meta_data
        .ms_run
        .entry(2)
        .or_default()
        .location
        .set("b.mzML");
    document.protein_data.push(row.clone());

    let adapter = MzTabFile::new();
    let layout =
        SectionLayout::for_protein(&document.protein_data, &document.meta_data).expect("layout");
    let header = adapter.protein_header(&layout).expect("header");
    let rendered = adapter.protein_row(&row, &layout).expect("row");
    let names: Vec<&str> = header.split('\t').collect();
    let cells: Vec<&str> = rendered.split('\t').collect();
    assert_eq!(names.len(), cells.len());
    for (score, run, expected) in [
        (1usize, 1usize, "11.0"),
        (1, 2, "12.0"),
        (2, 1, "21.0"),
        (2, 2, "22.0"),
    ] {
        let column = format!("search_engine_score[{score}]_ms_run[{run}]");
        let position = names
            .iter()
            .position(|name| *name == column)
            .unwrap_or_else(|| panic!("{column} is in the header"));
        assert_eq!(cells[position], expected, "{column}");
    }
    // And the document round-trips, which is what the alignment buys.
    let reloaded = round_trip(&adapter, &document);
    assert_eq!(reloaded, with_blanks_of(&document, &reloaded));
}

#[test]
fn complete_mode_declares_a_score_column_for_every_declared_ms_run() {
    // Source: store() sets search_ms_runs to the metadata's ms_run count in
    // Complete mode (MzTabFile.cpp:3168) even when no row carries a score for
    // one of them.
    let mut row = MzTabPeptideSectionRow::default();
    row.sequence.set("PEPTIDE");
    row.search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(1, MzTabDouble::new(5.0));
    let mut document = MzTab::default();
    document.meta_data.mz_tab_mode.set("Complete");
    document.meta_data.mz_tab_type.set("Identification");
    for run in 1..=3usize {
        document
            .meta_data
            .ms_run
            .entry(run)
            .or_default()
            .location
            .set(&format!("run{run}.mzML"));
    }
    document.peptide_data.push(row);
    let layout =
        SectionLayout::for_peptide(&document.peptide_data, &document.meta_data).expect("layout");
    assert_eq!(layout.ms_runs, [1, 2, 3].into_iter().collect());
    let header = MzTabFile::new().peptide_header(&layout).expect("header");
    assert!(header.contains("search_engine_score[1]_ms_run[3]"));

    // In Summary mode only the runs the rows carry are declared.
    let mut summary = document.clone();
    summary.meta_data.mz_tab_mode.set("Summary");
    let summary_layout =
        SectionLayout::for_peptide(&summary.peptide_data, &summary.meta_data).expect("layout");
    assert_eq!(summary_layout.ms_runs, [1].into_iter().collect());
}

// ---------------------------------------------------------------------------
// colunit
// ---------------------------------------------------------------------------

#[test]
fn colunit_keys_round_trip() {
    // Source MzTabFile.cpp:1984 writes `MTD\tcolunit-protein<value>` with no
    // separator, and :685 feeds the literal key `colunit` to
    // StringUtils::toInt32 while looking for a bracketed index, so any colunit
    // line makes its own reader throw a ConversionError.
    let mut document = MzTab::default();
    document.meta_data.mz_tab_mode.set("Summary");
    document.meta_data.mz_tab_type.set("Identification");
    document
        .meta_data
        .colunit_protein
        .push("protein_coverage=[UO, UO:0000187, percent, ]".to_owned());
    document
        .meta_data
        .colunit_psm
        .push("retention_time=[UO, UO:0000010, second, ]".to_owned());

    let adapter = MzTabFile::new();
    let rendered = adapter.write_to_string(&document).expect("renders");
    assert!(rendered.contains("MTD\tcolunit-protein\tprotein_coverage="));
    assert!(rendered.contains("MTD\tcolunit-psm\tretention_time="));
    assert_eq!(adapter.load_str(&rendered).expect("reloads"), document);

    // The source writes the PSM key in upper case; reading accepts either.
    let upper =
        format!("{MINIMAL_METADATA}MTD\tcolunit-PSM\tretention_time=[UO, UO:0000010, second, ]\n");
    let loaded = adapter.load_str(&upper).expect("upper-case key loads");
    assert_eq!(loaded.meta_data.colunit_psm.len(), 1);
}

// ---------------------------------------------------------------------------
// Diagnostics
// ---------------------------------------------------------------------------

#[test]
fn load_reporting_returns_the_mandatory_column_diagnostics() {
    // Source MzTabFile.cpp:832-934 prints each of these to std::cout, once per
    // data row, and the caller never sees them.
    let text = "MTD\tmzTab-version\t1.0.0\n\
         MTD\tmzTab-mode\tComplete\n\
         MTD\tmzTab-type\tQuantification\n\
         MTD\tdescription\td\n\
         MTD\tms_run[1]-location\ta.mzML\n\
         \n\
         PRH\taccession\tdescription\ttaxid\n\
         PRT\tP1\tone\t9606\n\
         PRT\tP2\ttwo\t9606\n";
    let (document, diagnostics) = MzTabFile::new()
        .load_reporting(text.as_bytes())
        .expect("document loads despite the missing columns");
    assert_eq!(document.protein_data.len(), 2);
    assert!(!diagnostics.is_empty());
    let joined = diagnostics.join("\n");
    for expected in [
        "mandatory protein species column missing",
        "mandatory protein best_search_engine_score[1-n] column missing",
        "mandatory protein ambiguity_members column missing",
        "mandatory protein_abundance_study_variable column(s) missing",
    ] {
        assert!(joined.contains(expected), "missing diagnostic: {expected}");
    }
    // Reported once per section, not once per row.
    assert_eq!(
        diagnostics
            .iter()
            .filter(|line| line.contains("mandatory protein species column missing"))
            .count(),
        1
    );
}

#[test]
fn a_row_shorter_than_its_header_is_reported_and_leaves_cells_unset() {
    let text = format!(
        "{MINIMAL_METADATA}\n\
         PSH\tsequence\tPSM_ID\tcharge\tspectra_ref\n\
         PSM\tPEPTIDE\t1\n"
    );
    let (document, diagnostics) = MzTabFile::new()
        .load_reporting(text.as_bytes())
        .expect("short row loads");
    let row = &document.psm_data[0];
    assert_eq!(row.sequence.get(), "PEPTIDE");
    assert_eq!(row.psm_id.get().unwrap(), 1);
    assert!(row.charge.is_null());
    assert!(row.spectra_ref.is_null());
    assert!(diagnostics.iter().any(|line| line.contains("shorter than")));
}

// ---------------------------------------------------------------------------
// Non-ASCII input and output
// ---------------------------------------------------------------------------

#[test]
fn non_ascii_content_round_trips_byte_for_byte() {
    let adapter = MzTabFile::new();
    let document = adapter
        .load(data("MzTabFile_unicode.mzTab"))
        .expect("unicode fixture loads");
    assert_eq!(document.meta_data.mz_tab_id.get(), "ユニコード試験");
    assert_eq!(
        document.meta_data.description.get(),
        "ペプチド同定の要約 — em dash and Ω"
    );
    assert_eq!(
        document.meta_data.ms_run[&1].location.get(),
        "file:///データ/試料.mzML"
    );
    assert_eq!(document.meta_data.colunit_psm.len(), 1);
    assert_eq!(document.psm_data.len(), 2);
    assert_eq!(document.psm_data[0].opt[0].name, "opt_global_備考");
    assert_eq!(document.psm_data[0].opt[0].value.get(), "日本語の注記");
    assert!(document.psm_data[0].calc_mass_to_charge.is_nan());
    // A sequence of two-byte characters: the second row's first cell after the
    // tag starts mid-way through the three-byte window a section tag occupies.
    assert_eq!(document.psm_data[1].sequence.get(), "ÄÖÜQRS");
    assert!(document.psm_data[1].search_engine_score[&1].is_inf());
    assert!(document.psm_data[1].calc_mass_to_charge.is_null());
    assert!(document.psm_data[1].opt[0].value.is_null());

    // Every value in this fixture has an exact binary representation and the
    // metadata is already in the writer's order, so the output is byte-exact.
    let rendered = adapter.write_to_string(&document).expect("renders");
    let original = std::fs::read_to_string(data("MzTabFile_unicode.mzTab")).expect("readable");
    assert_eq!(rendered, original);
    assert_eq!(adapter.load_str(&rendered).expect("reloads"), document);
}

#[test]
fn a_line_beginning_mid_character_is_not_mistaken_for_a_section_tag() {
    // `Ä` and `Ö` are two bytes each, so byte three of this line falls inside a
    // character. The source takes the first three bytes unconditionally; this
    // port asks for a character-boundary slice and treats the line as an
    // unknown section, exactly as the source's garbage prefix does.
    let text = format!("{MINIMAL_METADATA}ÄÖÜ\tone\ttwo\n");
    let document = MzTabFile::new().load_str(&text).expect("no panic");
    assert!(document.psm_data.is_empty());
    assert!(document.protein_data.is_empty());
    assert_eq!(document.meta_data.description.get(), "d");
}

#[test]
fn a_non_ascii_output_path_is_accepted_and_read_back() {
    let directory = TempDir::new_in(std::env::temp_dir(), false).expect("temp dir");
    let adapter = MzTabFile::new();
    let document = adapter
        .load(data("MzTabFile_unicode.mzTab"))
        .expect("fixture loads");
    for name in ["日本語.mzTab", "Ärger.mzTab", "résumé.tsv"] {
        let path = directory.path().join(name);
        adapter
            .store(&path, &document)
            .unwrap_or_else(|error| panic!("{name} stores: {error}"));
        let reloaded = adapter
            .load(&path)
            .unwrap_or_else(|error| panic!("{name} reloads: {error}"));
        assert_eq!(reloaded, document, "{name}");
    }
}

// ---------------------------------------------------------------------------
// Refusals and bounds
// ---------------------------------------------------------------------------

#[test]
fn store_rejects_an_extension_that_is_neither_mztab_nor_tsv() {
    let directory = TempDir::new_in(std::env::temp_dir(), false).expect("temp dir");
    let adapter = MzTabFile::new();
    let document = MzTab::default();
    // The source accepts mzTab, tsv and an unrecognised extension.
    for name in ["out.mzTab", "out.tsv", "out.unknown-suffix", "plain"] {
        adapter
            .store(directory.path().join(name), &document)
            .unwrap_or_else(|error| panic!("{name} should be accepted: {error}"));
    }
    for name in ["out.mzML", "out.csv", "out.idXML"] {
        let path = directory.path().join(name);
        let error = adapter
            .store(&path, &document)
            .expect_err("extension is rejected");
        assert!(matches!(error, Error::InvalidValue(_)), "{name}: {error}");
        assert!(!path.exists(), "{name} must not be created");
    }
}

#[test]
fn a_failed_store_leaves_the_previous_file_intact() {
    // The whole document is rendered and validated before the destination is
    // touched, and the bytes are published by renaming a temporary file, so a
    // rejected document cannot truncate a good output. The source opens the
    // destination with ios::trunc first.
    let directory = TempDir::new_in(std::env::temp_dir(), false).expect("temp dir");
    let path = directory.path().join("out.mzTab");
    let adapter = MzTabFile::new();
    let good = every_section_document();
    adapter.store(&path, &good).expect("first store succeeds");
    let before = std::fs::read_to_string(&path).expect("readable");

    // A modification cell with positions but no identifier cannot be rendered.
    let mut broken = good.clone();
    let mut modification = MzTabModification::null();
    modification.set_positions_and_parameters(vec![(3, MzTabParameter::null())]);
    broken.protein_data[0].modifications.set(vec![modification]);
    let error = adapter
        .store(&path, &broken)
        .expect_err("unrenderable modification is refused");
    assert!(matches!(error, Error::MissingInformation(_)), "{error}");
    assert_eq!(
        std::fs::read_to_string(&path).expect("readable"),
        before,
        "the previous output must survive"
    );

    // A cell carrying a tab is refused for the same reason.
    let mut tabbed = good.clone();
    tabbed.protein_data[0].description.set("a\tb");
    let error = adapter
        .store(&path, &tabbed)
        .expect_err("embedded tab is refused");
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    assert_eq!(std::fs::read_to_string(&path).expect("readable"), before);
}

#[test]
fn a_data_line_with_fewer_than_three_cells_is_a_parse_error() {
    // Source MzTabFile.cpp:249 throws Exception::ParseError.
    let text = format!("{MINIMAL_METADATA}PRT\tP1\n");
    let error = MzTabFile::new()
        .load_str(&text)
        .expect_err("two-cell data line is refused");
    match error {
        Error::Parse { line, ref message } => {
            assert_eq!(line, 5);
            assert!(message.contains("tabulator"), "{message}");
        }
        other => panic!("expected a parse error, got {other}"),
    }
}

#[test]
fn short_and_blank_lines_are_recorded_rather_than_parsed() {
    // Source MzTabFile.cpp:229: a line whose trimmed length is under three
    // bytes is recorded as an empty row, so a one- or two-character line never
    // reaches the section dispatch.
    let text = "MTD\tmzTab-version\t1.0.0\nab\n\n   \nMTD\tdescription\td\n";
    let document = MzTabFile::new().load_str(text).expect("loads");
    assert_eq!(document.empty_rows, vec![1, 2, 3]);
    assert_eq!(document.meta_data.description.get(), "d");
}

#[test]
fn comment_lines_are_preserved_verbatim_and_keyed_by_line() {
    let text = "COM\tfirst comment\nMTD\tmzTab-version\t1.0.0\nCOM\tsecond\nMTD\tdescription\td\n";
    let document = MzTabFile::new().load_str(text).expect("loads");
    assert_eq!(document.comment_rows.len(), 2);
    assert_eq!(document.comment_rows[&0], "COM\tfirst comment");
    assert_eq!(document.comment_rows[&2], "COM\tsecond");
}

#[test]
fn a_bad_index_in_a_metadata_key_is_a_parse_error() {
    let adapter = MzTabFile::new();
    for line in [
        "MTD\tassay[0]-sample_ref\tsample[1]\n",
        "MTD\tassay[9999999]-sample_ref\tsample[1]\n",
        "MTD\tms_run[x]-location\ta.mzML\n",
        "MTD\tsample_processing[0]\t[MS, MS:1000544, x, ]\n",
    ] {
        let text = format!("{MINIMAL_METADATA}{line}");
        let error = adapter
            .load_str(&text)
            .expect_err(&format!("{line:?} carries an invalid index"));
        assert!(matches!(error, Error::Parse { .. }), "{line:?}: {error}");
    }
}

#[test]
fn a_bad_index_in_a_section_header_is_a_parse_error() {
    let adapter = MzTabFile::new();
    for column in [
        "best_search_engine_score[0]",
        "best_search_engine_score[nope]",
        "search_engine_score[1]",
    ] {
        let text = format!("{MINIMAL_METADATA}\nPRH\taccession\t{column}\tmodifications\n");
        let error = adapter
            .load_str(&text)
            .expect_err(&format!("{column} carries an invalid index"));
        assert!(matches!(error, Error::Parse { .. }), "{column}: {error}");
    }
    // A valid pair is accepted.
    let text = format!(
        "{MINIMAL_METADATA}\nPRH\taccession\tsearch_engine_score[1]_ms_run[2]\tmodifications\n\
         PRT\tP1\t7.5\tnull\n"
    );
    let document = adapter.load_str(&text).expect("valid pair loads");
    assert_eq!(
        document.protein_data[0].search_engine_score_ms_run[&1][&2]
            .get()
            .unwrap(),
        7.5
    );
}

#[test]
fn a_cell_whose_text_its_column_rejects_is_a_parse_error() {
    let text = format!(
        "{MINIMAL_METADATA}\nPSH\tsequence\tPSM_ID\tcharge\nPSM\tPEPTIDE\t1\tnot-a-number\n"
    );
    let error = MzTabFile::new()
        .load_str(&text)
        .expect_err("a non-numeric charge is refused");
    assert!(matches!(error, Error::Parse { .. }), "{error}");
}

#[test]
fn section_layout_refuses_a_column_explosion() {
    // MAX_COLUMNS bounds the header a document can ask for, checked before any
    // row is rendered.
    let mut row = MzTabProteinSectionRow::default();
    row.accession.set("P1");
    for index in 1..=1000usize {
        row.best_search_engine_score
            .insert(index, MzTabDouble::new(1.0));
    }
    let meta = MzTabMetaData::default();
    let ok = SectionLayout::for_protein(std::slice::from_ref(&row), &meta)
        .expect("a thousand columns is within the limit");
    assert_eq!(ok.best_search_engine_score.len(), 1000);

    let mut wide = MzTabProteinSectionRow::default();
    wide.accession.set("P2");
    for index in 1..=200_001usize {
        wide.num_psms_ms_run.insert(index, 1.into());
    }
    let error = SectionLayout::for_protein(std::slice::from_ref(&wide), &meta)
        .expect_err("beyond MAX_COLUMNS the layout is refused");
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
}

#[test]
fn a_second_section_header_replaces_the_first() {
    // The source accumulates header columns into the same maps, so a document
    // with two PRH lines keeps the first one's column numbers and reads the
    // second section's rows through them.
    let text = format!(
        "{MINIMAL_METADATA}\n\
         PRH\taccession\tdescription\tmodifications\n\
         PRT\tP1\tfirst\tnull\n\
         \n\
         PRH\tmodifications\tdescription\taccession\n\
         PRT\tnull\tsecond\tP2\n"
    );
    let document = MzTabFile::new().load_str(&text).expect("loads");
    assert_eq!(document.protein_data.len(), 2);
    assert_eq!(document.protein_data[0].accession.get(), "P1");
    assert_eq!(document.protein_data[0].description.get(), "first");
    assert_eq!(document.protein_data[1].accession.get(), "P2");
    assert_eq!(document.protein_data[1].description.get(), "second");
}

#[test]
fn an_unknown_section_tag_is_ignored() {
    let text = format!("{MINIMAL_METADATA}XYZ\tone\ttwo\nMTD\ttitle\tkept\n");
    let document = MzTabFile::new().load_str(&text).expect("loads");
    assert_eq!(document.meta_data.title.get(), "kept");
    assert!(document.protein_data.is_empty());
}

#[test]
fn an_unknown_metadata_key_is_ignored() {
    let text = format!("{MINIMAL_METADATA}MTD\tnot_a_real_key[1]-suffix\tvalue\n");
    let document = MzTabFile::new().load_str(&text).expect("loads");
    assert_eq!(document.meta_data.description.get(), "d");
}

#[test]
fn an_empty_document_still_carries_the_three_mandatory_metadata_lines() {
    let document = MzTab::default();
    let lines = MzTabFile::new().document_lines(&document).expect("renders");
    assert_eq!(
        lines,
        vec![
            "MTD\tmzTab-version\t1.0.0".to_owned(),
            "MTD\tmzTab-mode\tnull".to_owned(),
            "MTD\tmzTab-type\tnull".to_owned(),
            "MTD\tdescription\tnull".to_owned(),
        ]
    );
    let reloaded = MzTabFile::new()
        .load_str(
            &MzTabFile::new()
                .write_to_string(&document)
                .expect("renders"),
        )
        .expect("reloads");
    assert_eq!(reloaded, document);
}
