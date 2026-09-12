// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! MzTab data model: `FORMAT/MzTabBase.h` and `FORMAT/MzTab.h`.
//!
//! Expected values marked "MzTab_test.cpp" are transcribed from that class
//! test (tier 3, source review). Everything else is derived from the two
//! headers and `MzTab.cpp`/`MzTabBase.cpp` by source review, or is a native
//! invariant or resource bound (tier 4). No C++ was built or executed.

use openms::Error;
use openms::chemistry::ModificationsDB;
use openms::format::mztab::{
    MAX_CELL_ITEMS, MzTab, MzTabBoolean, MzTabCVMetaData, MzTabCell, MzTabCellState,
    MzTabContactMetaData, MzTabDouble, MzTabDoubleList, MzTabInstrumentMetaData, MzTabInteger,
    MzTabIntegerList, MzTabModification, MzTabModificationList, MzTabNucleicAcidSectionRow,
    MzTabOSMSectionRow, MzTabOligonucleotideSectionRow, MzTabOptionalColumnEntry,
    MzTabPSMSectionRow, MzTabParameter, MzTabParameterList, MzTabPeptideSectionRow,
    MzTabProteinSectionRow, MzTabSampleMetaData, MzTabSmallMoleculeSectionRow,
    MzTabSoftwareMetaData, MzTabSpectraRef, MzTabString, MzTabStringList,
    add_meta_info_to_optional_columns, fixed_modification_metadata, modification_metadata,
    modification_metadata_with, optional_column_names, parse_cell, variable_modification_metadata,
};
use openms::identification::{FlankingResidue, PeptideEvidence};
use openms::metadata::{MetaInfo, MetaValue, MetaValueData};
use std::collections::{BTreeMap, BTreeSet};

fn text(value: &str) -> MzTabString {
    MzTabString::from_text(value)
}

// ---------------------------------------------------------------------------
// MzTab_test.cpp START_SECTION(MzTab())
// ---------------------------------------------------------------------------

#[test]
fn document_default_construction_declares_version_1_0_0() {
    // MzTab_test.cpp asserts only that `new MzTab()` is not the null pointer.
    // The equivalent observable in Rust is that a default document exists and
    // carries the metadata its constructor sets: MzTabMetaData() runs
    // mz_tab_version.fromCellString("1.0.0") (MzTab.cpp:301).
    let document = MzTab::default();
    assert_eq!(document.meta_data.mz_tab_version.get(), "1.0.0");
    assert_eq!(document.meta_data.mz_tab_version.to_cell_string(), "1.0.0");
    assert!(document.protein_data.is_empty());
    assert!(document.peptide_data.is_empty());
    assert!(document.psm_data.is_empty());
    assert!(document.small_molecule_data.is_empty());
    assert!(document.nucleic_acid_data.is_empty());
    assert!(document.oligonucleotide_data.is_empty());
    assert!(document.osm_data.is_empty());
    assert!(document.empty_rows.is_empty());
    assert!(document.comment_rows.is_empty());
    // Every other metadata field starts null or empty.
    assert!(document.meta_data.mz_tab_mode.is_null());
    assert!(document.meta_data.title.is_null());
    assert!(document.meta_data.quantification_method.is_null());
    assert!(document.meta_data.false_discovery_rate.is_null());
    assert!(document.meta_data.ms_run.is_empty());
    assert!(document.meta_data.colunit_protein.is_empty());
}

// ---------------------------------------------------------------------------
// MzTab_test.cpp START_SECTION(~MzTab())
// ---------------------------------------------------------------------------

#[test]
fn document_destruction_releases_a_populated_document() {
    // The upstream section has no assertion macro at all: it calls `delete
    // ptr` and relies on the test binary not crashing. Rust's Drop is
    // compiler-generated, so the equivalent is that a populated document and
    // every owned cell inside it can be dropped while an independent clone
    // stays valid.
    let mut document = MzTab::default();
    let mut row = MzTabPSMSectionRow::default();
    row.sequence.set("NDYKAPPQPAPGK");
    row.opt.push(MzTabOptionalColumnEntry::new(
        "Percolator_score",
        text("0.359083"),
    ));
    document.psm_data.push(row);
    let survivor = document.clone();
    drop(document);
    assert_eq!(survivor.psm_data.len(), 1);
    assert_eq!(survivor.psm_data[0].sequence.get(), "NDYKAPPQPAPGK");
    assert_eq!(survivor.psm_data[0].opt[0].value.get(), "0.359083");
}

// ---------------------------------------------------------------------------
// MzTab_test.cpp START_SECTION(std::vector<std::string> getPSMOptionalColumnNames() const)
// ---------------------------------------------------------------------------

#[test]
fn psm_optional_column_names_from_the_upstream_two_row_fixture() {
    // Transcribed from MzTab_test.cpp lines 46-129. The upstream fixture reuses
    // one `row` for both rows and never clears `opt_`, so the second row
    // carries nine optional entries and the union has five distinct names.
    let mut document = MzTab::default();
    let mut rows: Vec<MzTabPSMSectionRow> = Vec::new();
    let mut row = MzTabPSMSectionRow::default();

    // row 1
    row.sequence.from_cell_string("NDYKAPPQPAPGK");
    row.psm_id.from_cell_string("38").unwrap();
    row.accession.from_cell_string("IPI:B1");
    row.unique.from_cell_string("1").unwrap();
    row.database.from_cell_string("null");
    row.database_version.from_cell_string("null");
    row.search_engine
        .from_cell_string("[, , Percolator, ]")
        .unwrap();
    row.search_engine_score
        .entry(0)
        .or_default()
        .from_cell_string("51.9678841193106")
        .unwrap();
    for (name, value) in [
        ("Percolator_score", "0.359083"),
        ("Percolator_qvalue", "0.00649874"),
        ("Percolator_PEP", "0.0420992"),
        ("search_engine_sequence", "NDYKAPPQPAPGK"),
    ] {
        row.opt
            .push(MzTabOptionalColumnEntry::new(name, text(value)));
    }
    rows.push(row.clone());

    // row 2 — the upstream fixture keeps mutating the same object
    row.sequence.from_cell_string("IRRS(Phospho)SFSSK");
    row.psm_id.from_cell_string("39").unwrap();
    row.accession.from_cell_string("IPI:IPI00009899.4");
    row.unique.from_cell_string("0").unwrap();
    row.database.from_cell_string("null");
    row.database_version.from_cell_string("null");
    row.search_engine
        .from_cell_string("[, , Percolator, ]")
        .unwrap();
    row.search_engine_score
        .entry(0)
        .or_default()
        .from_cell_string("9.55915773892318")
        .unwrap();
    for (name, value) in [
        ("Percolator_score", "0.157068"),
        ("Percolator_qvalue", "0.00774619"),
        ("Percolator_PEP", "0.0779777"),
        ("search_engine_sequence", "IRRSSFS(Phospho)SK"),
        ("AScore_1", "3.64384830671351"),
    ] {
        row.opt
            .push(MzTabOptionalColumnEntry::new(name, text(value)));
    }
    rows.push(row);

    document.set_psm_section_rows(rows);
    let optional_columns = document.psm_optional_column_names().unwrap();

    // The two upstream assertions.
    assert_eq!(document.psm_section_rows().len(), 2);
    assert_eq!(optional_columns.len(), 5);

    // Column order is first-occurrence order, which the upstream test does not
    // check but the source's vector-and-linear-find deliberately preserves.
    assert_eq!(
        optional_columns,
        vec![
            "Percolator_score".to_owned(),
            "Percolator_qvalue".to_owned(),
            "Percolator_PEP".to_owned(),
            "search_engine_sequence".to_owned(),
            "AScore_1".to_owned(),
        ]
    );
    assert_eq!(document.psm_section_rows()[0].opt.len(), 4);
    assert_eq!(document.psm_section_rows()[1].opt.len(), 9);

    // Every cell in the fixture parsed as the source parses it.
    let first = &document.psm_section_rows()[0];
    assert_eq!(first.sequence.get(), "NDYKAPPQPAPGK");
    assert_eq!(first.psm_id.get().unwrap(), 38);
    assert_eq!(first.accession.get(), "IPI:B1");
    assert_eq!(first.unique.as_bool(), Some(true));
    assert!(first.database.is_null());
    assert!(first.database_version.is_null());
    assert_eq!(first.search_engine.get().len(), 1);
    assert_eq!(first.search_engine.get()[0].name(), "Percolator");
    assert_eq!(
        first.search_engine.to_cell_string(),
        "[, , Percolator, ]",
        "a Param whose parts contain no \", \" round-trips byte for byte"
    );
    assert_eq!(
        first.search_engine_score[&0].get().unwrap(),
        51.9678841193106
    );
    let second = &document.psm_section_rows()[1];
    assert_eq!(second.unique.as_bool(), Some(false));
    assert_eq!(second.psm_id.get().unwrap(), 39);
    assert_eq!(
        second.search_engine_score[&0].get().unwrap(),
        9.55915773892318
    );

    // getNumberOfPSMs counts distinct PSM_IDs, not rows.
    assert_eq!(document.number_of_psms().unwrap(), 2);
}

// ---------------------------------------------------------------------------
// MzTab_test.cpp START_SECTION(static void addMetaInfoToOptionalColumns(...))
// ---------------------------------------------------------------------------

#[test]
fn add_meta_info_to_optional_columns_matches_the_upstream_seven_assertions() {
    // MzTab_test.cpp lines 132-151: keys have spaces replaced by underscores,
    // values are taken as they are, and a key missing from the metadata still
    // produces a column whose value is null.
    let keys: BTreeSet<String> = ["FWHM", "with space", "ppm_errors"]
        .iter()
        .map(|key| (*key).to_owned())
        .collect();
    let mut meta = MetaInfo::new();
    meta.insert("FWHM".to_owned(), MetaValue::try_from(34.5).unwrap());
    meta.insert(
        "ppm_errors".to_owned(),
        MetaValue::try_from(vec![0.5, 1.4, -2.0, 0.1]).unwrap(),
    );
    let mut opt: Vec<MzTabOptionalColumnEntry> = Vec::new();
    add_meta_info_to_optional_columns(&keys, &mut opt, "global", &meta).unwrap();

    assert_eq!(opt.len(), 3);
    assert_eq!(opt[0].name, "opt_global_FWHM");
    assert_eq!(opt[1].name, "opt_global_ppm_errors");
    assert_eq!(opt[2].name, "opt_global_with_space");
    assert_eq!(opt[0].value.to_cell_string(), "34.5");
    assert_eq!(opt[1].value.to_cell_string(), "[0.5, 1.4, -2.0, 0.1]");
    assert_eq!(opt[2].value.to_cell_string(), "null");
}

#[test]
fn add_meta_info_to_optional_columns_appends_and_renders_every_value_kind() {
    let mut keys = BTreeSet::new();
    for key in ["empty", "count", "counts", "names", "note"] {
        keys.insert(key.to_owned());
    }
    let mut meta = MetaInfo::new();
    meta.insert(
        "empty".to_owned(),
        MetaValue::new(MetaValueData::Empty).unwrap(),
    );
    meta.insert("count".to_owned(), MetaValue::from(-7_i64));
    meta.insert("counts".to_owned(), MetaValue::from(vec![1_i64, -2, 3]));
    meta.insert(
        "names".to_owned(),
        MetaValue::from(vec!["a".to_owned(), "b".to_owned()]),
    );
    meta.insert("note".to_owned(), MetaValue::from("kept, verbatim"));

    let mut opt = vec![MzTabOptionalColumnEntry::new(
        "opt_global_existing",
        text("x"),
    )];
    add_meta_info_to_optional_columns(&keys, &mut opt, "global", &meta).unwrap();
    let rendered: Vec<(String, String)> = opt
        .iter()
        .map(|entry| (entry.name.clone(), entry.value.to_cell_string()))
        .collect();
    assert_eq!(
        rendered,
        vec![
            ("opt_global_existing".to_owned(), "x".to_owned()),
            ("opt_global_count".to_owned(), "-7".to_owned()),
            ("opt_global_counts".to_owned(), "[1, -2, 3]".to_owned()),
            // An EMPTY DataValue stringifies to "", which MzTabString reads as
            // null: a present-but-empty meta value is indistinguishable from an
            // absent one in the rendered cell.
            ("opt_global_empty".to_owned(), "null".to_owned()),
            ("opt_global_names".to_owned(), "[a, b]".to_owned()),
            ("opt_global_note".to_owned(), "kept, verbatim".to_owned()),
        ]
    );
}

#[test]
fn add_meta_info_to_optional_columns_refuses_past_the_column_ceiling() {
    let keys: BTreeSet<String> = (0..3).map(|i| format!("k{i}")).collect();
    let mut opt: Vec<MzTabOptionalColumnEntry> = (0..MzTab::MAX_OPTIONAL_COLUMNS)
        .map(|i| MzTabOptionalColumnEntry::new(format!("opt_global_{i}"), MzTabString::default()))
        .collect();
    let before = opt.len();
    let error = add_meta_info_to_optional_columns(&keys, &mut opt, "global", &MetaInfo::new());
    assert!(matches!(error, Err(Error::InvalidValue(_))));
    assert_eq!(opt.len(), before, "the vector is unchanged on refusal");
}

// ---------------------------------------------------------------------------
// MzTab_test.cpp START_SECTION([EXTRA] MzTabBoolean setNull / isNull polarity)
// ---------------------------------------------------------------------------

#[test]
fn boolean_set_null_polarity() {
    // MzTab_test.cpp lines 195-206, all five assertions.
    let mut cell = MzTabBoolean::new(true);
    // The upstream macros are TEST_EQUAL(b.isNull(), false/true); clippy
    // rejects comparing against a bool literal, so the same three states are
    // asserted with assert!.
    assert!(!cell.is_null());
    cell.set_null(true);
    assert!(cell.is_null());
    assert_eq!(cell.to_cell_string(), "null");
    cell.set_null(false);
    assert!(!cell.is_null());
    // set_null(false) stores 0, not the previous `true`.
    assert_eq!(cell.value(), 0);
    assert_eq!(cell.to_cell_string(), "0");
}

// ---------------------------------------------------------------------------
// Cell vocabulary: MzTabDouble
// ---------------------------------------------------------------------------

#[test]
fn double_default_is_null_not_zero() {
    let cell = MzTabDouble::default();
    assert!(cell.is_null());
    assert_eq!(cell.state(), MzTabCellState::Null);
    assert_eq!(cell.to_cell_string(), "null");
    assert_eq!(cell.raw_value(), 0.0);
    assert!(matches!(cell.get(), Err(Error::MissingInformation(_))));
}

#[test]
fn double_states_render_and_parse() {
    let cases = [
        ("null", MzTabCellState::Null, "null"),
        ("NULL", MzTabCellState::Null, "null"),
        ("  null  ", MzTabCellState::Null, "null"),
        ("NaN", MzTabCellState::NaN, "NaN"),
        ("nan", MzTabCellState::NaN, "NaN"),
        ("Inf", MzTabCellState::Inf, "Inf"),
        ("INF", MzTabCellState::Inf, "Inf"),
    ];
    for (input, state, rendered) in cases {
        let cell: MzTabDouble = parse_cell(input).unwrap();
        assert_eq!(cell.state(), state, "{input:?}");
        assert_eq!(cell.to_cell_string(), rendered, "{input:?}");
        assert!(cell.get().is_err(), "{input:?} carries no value");
    }
    assert!(MzTabDouble::nan().is_nan());
    assert!(MzTabDouble::inf().is_inf());
    assert!(MzTabDouble::null().is_null());
}

#[test]
fn double_value_rendering_follows_the_source_number_convention() {
    // Fifteen fractional digits, trailing zeros trimmed to at least one.
    assert_eq!(MzTabDouble::new(34.5).to_cell_string(), "34.5");
    assert_eq!(MzTabDouble::new(-2.0).to_cell_string(), "-2.0");
    assert_eq!(MzTabDouble::new(0.0).to_cell_string(), "0.0");
    assert_eq!(MzTabDouble::new(0.1).to_cell_string(), "0.1");
    assert_eq!(MzTabDouble::new(-18.010565).to_cell_string(), "-18.010565");
    // Fifteen FRACTIONAL digits is not fifteen significant digits, so the
    // shortest form is not what the source writes.
    assert_eq!(
        MzTabDouble::new(51.9678841193106).to_cell_string(),
        "51.967884119310597"
    );
    // Scientific outside [1e-2, 1e4), with a '+'-free two-digit exponent.
    assert_eq!(MzTabDouble::new(1e5).to_cell_string(), "1.0e05");
    assert_eq!(MzTabDouble::new(0.001).to_cell_string(), "1.0e-03");
    // A non-finite value in the Default state renders the number, not the
    // state keyword: lowercase "inf" against the "Inf" state.
    assert_eq!(MzTabDouble::new(f64::INFINITY).to_cell_string(), "inf");
    assert_eq!(MzTabDouble::new(f64::NEG_INFINITY).to_cell_string(), "-inf");
    assert_eq!(MzTabDouble::new(f64::NAN).to_cell_string(), "NaN");
}

#[test]
fn double_round_trips_through_its_own_text() {
    for value in [0.0_f64, 34.5, -18.010565, 51.9678841193106, 1e5, 0.001] {
        let rendered = MzTabDouble::new(value).to_cell_string();
        let parsed: MzTabDouble = parse_cell(&rendered).unwrap();
        assert_eq!(parsed.get().unwrap(), value, "{rendered}");
    }
}

#[test]
fn double_parse_rejects_text_that_is_not_a_number() {
    for input in ["", "abc", "1.5x", "--1", "0x10", "1,5", "日本語"] {
        let result: Result<MzTabDouble, Error> = parse_cell(input);
        assert!(
            matches!(result, Err(Error::Parse { .. })),
            "{input:?} must be a conversion error"
        );
    }
    // The source's std::from_chars reports result_out_of_range for an
    // overflowing literal; Rust's parse would silently yield infinity.
    let result: Result<MzTabDouble, Error> = parse_cell("1e999");
    assert!(matches!(result, Err(Error::Parse { .. })));
    // An explicit signed infinity spelling is accepted, as toDouble accepts it.
    let cell: MzTabDouble = parse_cell("-inf").unwrap();
    assert_eq!(cell.state(), MzTabCellState::Default);
    assert_eq!(cell.get().unwrap(), f64::NEG_INFINITY);
}

#[test]
fn double_source_comparison_ignores_the_cell_state() {
    let null = MzTabDouble::null();
    let zero = MzTabDouble::new(0.0);
    // operator== compares only the stored number, so these are "equal".
    assert!(null.source_equal(&zero));
    // Rust equality is stricter and separates the two states.
    assert_ne!(null, zero);
    assert!(null.source_less(&MzTabDouble::new(1.0)));
    assert!(!MzTabDouble::new(1.0).source_less(&null));
    // A NaN-valued Default cell is neither less than nor equal to anything.
    let nan = MzTabDouble::new(f64::NAN);
    assert!(!nan.source_less(&zero));
    assert!(!nan.source_equal(&nan));
}

#[test]
fn double_set_null_false_publishes_the_stored_value() {
    let mut cell = MzTabDouble::default();
    cell.set_null(false);
    assert_eq!(cell.get().unwrap(), 0.0);
    cell.set(7.25);
    cell.set_nan();
    assert!(cell.get().is_err());
    cell.set_null(false);
    assert_eq!(
        cell.get().unwrap(),
        7.25,
        "the state changes, not the value"
    );
}

// ---------------------------------------------------------------------------
// Cell vocabulary: MzTabInteger
// ---------------------------------------------------------------------------

#[test]
fn integer_default_is_null_and_states_round_trip() {
    let cell = MzTabInteger::default();
    assert!(cell.is_null());
    assert_eq!(cell.to_cell_string(), "null");
    assert!(matches!(cell.get(), Err(Error::MissingInformation(_))));
    assert_eq!(MzTabInteger::new(-5).to_cell_string(), "-5");
    assert_eq!(MzTabInteger::new(0).to_cell_string(), "0");
    assert_eq!(MzTabInteger::nan().to_cell_string(), "NaN");
    assert_eq!(MzTabInteger::inf().to_cell_string(), "Inf");
    let parsed: MzTabInteger = parse_cell("NaN").unwrap();
    assert!(parsed.is_nan());
    assert!(parsed.get().is_err());
}

#[test]
fn integer_accepts_a_float_spelling_because_external_files_write_one() {
    // MzTabBase.cpp:673 comment: "some mzTab files from external sources
    // contain floating point numbers in integer columns".
    let cell: MzTabInteger = parse_cell("4.0").unwrap();
    assert_eq!(cell.get().unwrap(), 4);
    assert_eq!(cell.to_cell_string(), "4");
    let cell: MzTabInteger = parse_cell("-0.0").unwrap();
    assert_eq!(cell.get().unwrap(), 0);
    // A fractional value is still a conversion error.
    let result: Result<MzTabInteger, Error> = parse_cell("4.5");
    assert!(matches!(result, Err(Error::Parse { .. })));
}

#[test]
fn integer_refuses_a_value_outside_the_32_bit_range() {
    for input in ["3000000000", "-3000000000", "1e30"] {
        let result: Result<MzTabInteger, Error> = parse_cell(input);
        assert!(
            matches!(result, Err(Error::Parse { .. })),
            "{input:?} is outside Int"
        );
    }
    let cell: MzTabInteger = parse_cell("2147483647").unwrap();
    assert_eq!(cell.get().unwrap(), i32::MAX);
    let cell: MzTabInteger = parse_cell("-2147483648").unwrap();
    assert_eq!(cell.get().unwrap(), i32::MIN);
}

// ---------------------------------------------------------------------------
// Cell vocabulary: MzTabBoolean
// ---------------------------------------------------------------------------

#[test]
fn boolean_accepts_only_the_two_exact_digits() {
    let mut cell = MzTabBoolean::default();
    assert!(cell.is_null());
    assert_eq!(cell.value(), -1);
    assert_eq!(cell.as_bool(), None);
    cell.from_cell_string("1").unwrap();
    assert_eq!(cell.as_bool(), Some(true));
    cell.from_cell_string("0").unwrap();
    assert_eq!(cell.as_bool(), Some(false));
    cell.from_cell_string("NULL").unwrap();
    assert!(cell.is_null());
    // The source trims for the "null" test but compares the two digits against
    // the untrimmed text, so a padded digit is a conversion error.
    for input in [" 1", "1 ", "true", "2", "-1", ""] {
        let mut cell = MzTabBoolean::default();
        assert!(
            matches!(cell.from_cell_string(input), Err(Error::Parse { .. })),
            "{input:?} must be a conversion error"
        );
    }
}

// ---------------------------------------------------------------------------
// Cell vocabulary: MzTabString
// ---------------------------------------------------------------------------

#[test]
fn string_trims_and_treats_the_literal_null_as_null() {
    assert_eq!(MzTabString::from_text("  hello  ").get(), "hello");
    assert_eq!(
        MzTabString::from_text("\t x \r\n").get(),
        "x",
        "trimming covers space, tab, CR and LF"
    );
    assert!(MzTabString::from_text("null").is_null());
    assert!(MzTabString::from_text("NULL").is_null());
    assert!(MzTabString::from_text(" Null ").is_null());
    assert!(MzTabString::from_text("").is_null());
    assert_eq!(MzTabString::default().to_cell_string(), "null");
    assert_eq!(MzTabString::from_text("nullable").get(), "nullable");
    let mut cell = MzTabString::from_text("kept");
    cell.set_null(false);
    assert_eq!(cell.get(), "kept", "set_null(false) is ignored");
    cell.set_null(true);
    assert!(cell.is_null());
}

#[test]
fn string_keeps_non_ascii_text_intact() {
    // The source trims four ASCII bytes only, so U+00A0 NO-BREAK SPACE stays.
    let cell = MzTabString::from_text(" 日本語/dir ");
    assert_eq!(cell.get(), "日本語/dir");
    assert_eq!(cell.to_cell_string(), "日本語/dir");
    let padded = MzTabString::from_text("\u{a0}日本語\u{a0}");
    assert_eq!(padded.get(), "\u{a0}日本語\u{a0}");
    // A fullwidth spelling of "null" is not the null keyword.
    let fullwidth = MzTabString::from_text("ＮＵＬＬ");
    assert!(!fullwidth.is_null());
    assert_eq!(fullwidth.get(), "ＮＵＬＬ");
    // Round trip through a cell.
    let mut parsed = MzTabString::default();
    parsed.read_cell(&cell.to_cell_string()).unwrap();
    assert_eq!(parsed, cell);
}

// ---------------------------------------------------------------------------
// Cell vocabulary: MzTabParameter
// ---------------------------------------------------------------------------

#[test]
fn parameter_parses_the_four_bracketed_fields() {
    let parameter =
        MzTabParameter::parse("[MS, MS:1002453, No fixed modifications searched, ]").unwrap();
    assert_eq!(parameter.cv_label(), "MS");
    assert_eq!(parameter.accession(), "MS:1002453");
    assert_eq!(parameter.name(), "No fixed modifications searched");
    assert_eq!(parameter.value(), "");
    assert!(!parameter.is_null());
    assert_eq!(
        parameter.to_cell_string(),
        "[MS, MS:1002453, No fixed modifications searched, ]"
    );
    // An all-empty parameter is null.
    let empty = MzTabParameter::parse("[, , , ]").unwrap();
    assert!(empty.is_null());
    assert_eq!(empty.to_cell_string(), "null");
    // The literal null.
    let mut parameter = MzTabParameter::from_parts("MS", "MS:1", "n", "v");
    parameter.from_cell_string("null").unwrap();
    assert!(parameter.is_null());
}

#[test]
fn parameter_quotes_only_a_comma_followed_by_a_space() {
    let quoted = MzTabParameter::from_parts("MS", "MS:1", "a, b", "c, d");
    assert_eq!(quoted.to_cell_string(), "[MS, MS:1, \"a, b\", \"c, d\"]");
    let reparsed = MzTabParameter::parse(&quoted.to_cell_string()).unwrap();
    assert_eq!(reparsed, quoted);
    // A bare comma is NOT quoted, and the cell it produces has five fields and
    // cannot be read back. This is a source defect, preserved.
    let unquoted = MzTabParameter::from_parts("MS", "MS:1", "a,b", "v");
    assert_eq!(unquoted.to_cell_string(), "[MS, MS:1, a,b, v]");
    let result = MzTabParameter::parse(&unquoted.to_cell_string());
    assert!(matches!(result, Err(Error::Parse { .. })));
}

#[test]
fn parameter_scanner_drops_brackets_even_inside_quotes() {
    // MzTabBase.cpp:393 skips '[' and ']' unconditionally, so the source's own
    // "quoted brackets" example loses them.
    let parameter = MzTabParameter::parse("[MS, MS:1, \"a, [b]\", v]").unwrap();
    assert_eq!(parameter.name(), "a, b");
    assert_eq!(parameter.value(), "v");
}

#[test]
fn parameter_requires_exactly_four_fields() {
    for input in ["[MS, MS:1, n]", "[MS, MS:1, n, v, extra]", "", "[]"] {
        let result = MzTabParameter::parse(input);
        assert!(
            matches!(result, Err(Error::Parse { .. })),
            "{input:?} must be a conversion error"
        );
    }
}

#[test]
fn parameter_keeps_non_ascii_names() {
    let parameter = MzTabParameter::parse("[MS, MS:1, 日本語, ]").unwrap();
    assert_eq!(parameter.name(), "日本語");
    assert_eq!(parameter.to_cell_string(), "[MS, MS:1, 日本語, ]");
}

// ---------------------------------------------------------------------------
// Cell vocabulary: the four list types
// ---------------------------------------------------------------------------

#[test]
fn parameter_list_joins_with_a_pipe_and_refuses_a_null_member() {
    let mut list = MzTabParameterList::default();
    assert!(list.is_null());
    assert_eq!(list.to_cell_string(), "null");
    list.from_cell_string("[MS, MS:1, a, ]|[MS, MS:2, b, ]")
        .unwrap();
    assert_eq!(list.get().len(), 2);
    assert_eq!(list.get()[1].accession(), "MS:2");
    assert_eq!(list.to_cell_string(), "[MS, MS:1, a, ]|[MS, MS:2, b, ]");
    // A member spelled "null" is a conversion error (MzTabBase.cpp:68).
    let mut list = MzTabParameterList::default();
    let result = list.from_cell_string("[MS, MS:1, a, ]|null");
    assert!(matches!(result, Err(Error::Parse { .. })));
    assert!(list.is_null(), "a refused parse leaves the list untouched");
}

#[test]
fn string_list_uses_its_configured_separator() {
    let mut list = MzTabStringList::default();
    assert_eq!(list.separator(), '|');
    list.from_cell_string("a|b|c").unwrap();
    assert_eq!(list.get().len(), 3);
    assert_eq!(list.to_cell_string(), "a|b|c");
    let mut list = MzTabStringList::default();
    list.set_separator(',');
    list.from_cell_string("GO:0001,GO:0002").unwrap();
    assert_eq!(list.get().len(), 2);
    assert_eq!(list.get()[0].get(), "GO:0001");
    assert_eq!(list.to_cell_string(), "GO:0001,GO:0002");
    // An empty subject produces no entries at all, so the list stays null.
    let mut list = MzTabStringList::default();
    list.from_cell_string("").unwrap();
    assert!(list.is_null());
    // Multi-byte entries survive the split.
    let mut list = MzTabStringList::default();
    list.from_cell_string("日|本").unwrap();
    assert_eq!(list.get().len(), 2);
    assert_eq!(list.get()[0].get(), "日");
    assert_eq!(list.get()[1].get(), "本");
}

#[test]
fn integer_list_joins_with_a_comma_and_keeps_member_states() {
    let mut list = MzTabIntegerList::default();
    list.from_cell_string("1,null,3").unwrap();
    assert_eq!(list.get().len(), 3);
    assert!(list.get()[1].is_null());
    assert_eq!(list.to_cell_string(), "1,null,3");
    let mut list = MzTabIntegerList::default();
    list.set(vec![MzTabInteger::new(-1), MzTabInteger::inf()]);
    assert_eq!(list.to_cell_string(), "-1,Inf");
    let result = MzTabIntegerList::default().from_cell_string("1,x");
    assert!(matches!(result, Err(Error::Parse { .. })));
}

#[test]
fn double_list_joins_with_a_pipe_and_keeps_member_states() {
    let mut list = MzTabDoubleList::default();
    list.from_cell_string("1.5|NaN|null").unwrap();
    assert_eq!(list.get().len(), 3);
    assert_eq!(list.get()[0].get().unwrap(), 1.5);
    assert!(list.get()[1].is_nan());
    assert!(list.get()[2].is_null());
    assert_eq!(list.to_cell_string(), "1.5|NaN|null");
    let mut list = MzTabDoubleList::default();
    assert_eq!(list.to_cell_string(), "null");
    list.from_cell_string("null").unwrap();
    assert!(list.is_null());
}

#[test]
fn a_second_parse_replaces_rather_than_appends() {
    // The source's list parsers push onto the existing vector without clearing
    // it, so parsing twice concatenates. This port replaces.
    let mut list = MzTabDoubleList::default();
    list.from_cell_string("1.0|2.0").unwrap();
    list.from_cell_string("3.0").unwrap();
    assert_eq!(list.get().len(), 1);
    assert_eq!(list.get()[0].get().unwrap(), 3.0);

    let mut list = MzTabIntegerList::default();
    list.from_cell_string("1,2").unwrap();
    list.from_cell_string("3").unwrap();
    assert_eq!(list.get().len(), 1);

    let mut list = MzTabStringList::default();
    list.from_cell_string("a|b").unwrap();
    list.from_cell_string("c").unwrap();
    assert_eq!(list.get().len(), 1);

    let mut list = MzTabParameterList::default();
    list.from_cell_string("[MS, MS:1, a, ]").unwrap();
    list.from_cell_string("[MS, MS:2, b, ]").unwrap();
    assert_eq!(list.get().len(), 1);
    assert_eq!(list.get()[0].accession(), "MS:2");

    let mut list = MzTabModificationList::default();
    list.from_cell_string("UNIMOD:35").unwrap();
    list.from_cell_string("UNIMOD:4").unwrap();
    assert_eq!(list.get().len(), 1);
}

#[test]
fn list_parsers_refuse_past_the_entry_ceiling() {
    let oversized = "1|".repeat(MAX_CELL_ITEMS + 1);
    let mut list = MzTabDoubleList::default();
    assert!(matches!(
        list.from_cell_string(&oversized),
        Err(Error::InvalidValue(_))
    ));
    assert!(list.is_null(), "the list is untouched on refusal");
    let mut list = MzTabStringList::default();
    assert!(matches!(
        list.from_cell_string(&oversized),
        Err(Error::InvalidValue(_))
    ));
    // Exactly at the ceiling is accepted.
    let at_ceiling = vec!["1"; MAX_CELL_ITEMS].join("|");
    let mut list = MzTabDoubleList::default();
    list.from_cell_string(&at_ceiling).unwrap();
    assert_eq!(list.get().len(), MAX_CELL_ITEMS);
}

// ---------------------------------------------------------------------------
// Cell vocabulary: MzTabSpectraRef
// ---------------------------------------------------------------------------

#[test]
fn spectra_ref_renders_and_parses_its_two_parts() {
    let reference = MzTabSpectraRef::new(3, "index=5").unwrap();
    assert_eq!(reference.ms_file(), 3);
    assert_eq!(reference.spec_ref(), "index=5");
    assert_eq!(reference.resolved(), Some((3, "index=5")));
    assert_eq!(reference.to_cell_string(), "ms_run[3]:index=5");
    let parsed: MzTabSpectraRef = parse_cell("ms_run[3]:index=5").unwrap();
    assert_eq!(parsed, reference);
    // A Thermo-style native identifier has no colon and parses.
    let parsed: MzTabSpectraRef =
        parse_cell("ms_run[2]:controllerType=0 controllerNumber=1 scan=17").unwrap();
    assert_eq!(parsed.ms_file(), 2);
    assert_eq!(
        parsed.spec_ref(),
        "controllerType=0 controllerNumber=1 scan=17"
    );
}

#[test]
fn spectra_ref_is_null_until_both_parts_are_set() {
    let mut reference = MzTabSpectraRef::default();
    assert!(reference.is_null());
    assert_eq!(reference.ms_file(), 0);
    assert_eq!(reference.spec_ref(), "");
    assert_eq!(reference.resolved(), None);
    assert_eq!(reference.to_cell_string(), "null");
    reference.set_ms_file(1).unwrap();
    assert!(reference.is_null(), "still null without a reference text");
    reference.set_spec_ref_file("scan=1").unwrap();
    assert!(!reference.is_null());
    reference.set_null(true);
    assert!(reference.is_null());
    reference.set_null(false);
    assert!(reference.is_null(), "set_null(false) is ignored");
}

#[test]
fn spectra_ref_refuses_the_states_the_source_only_asserts() {
    let mut reference = MzTabSpectraRef::default();
    // setMSFile(0) asserts, then silently does nothing in a release build.
    assert!(matches!(
        reference.set_ms_file(0),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        reference.set_spec_ref(""),
        Err(Error::InvalidValue(_))
    ));
    // A native identifier containing a colon splits into three fields.
    let result: Result<MzTabSpectraRef, Error> = parse_cell("ms_run[1]:file=a:1");
    assert!(matches!(result, Err(Error::Parse { .. })));
    // The source casts a negative index to Size; this rejects it.
    let result: Result<MzTabSpectraRef, Error> = parse_cell("ms_run[-1]:scan=1");
    assert!(matches!(result, Err(Error::Parse { .. })));
    let result: Result<MzTabSpectraRef, Error> = parse_cell("ms_run[x]:scan=1");
    assert!(matches!(result, Err(Error::Parse { .. })));
    // A zero index parses and leaves the reference null, as the source.
    let parsed: MzTabSpectraRef = parse_cell("ms_run[0]:scan=1").unwrap();
    assert!(parsed.is_null());
    assert_eq!(parsed.to_cell_string(), "null");
    let parsed: MzTabSpectraRef = parse_cell("null").unwrap();
    assert!(parsed.is_null());
}

/// The source reads the run index with `StringUtils::toInt32`, which parses
/// into an `Int32` — a token outside that range is a `ConversionError`
/// (`StringUtils.cpp:150-156`, pinned by `StringUtils_test.cpp:272`) — and
/// which advances past exactly one `+` (`StringUtils.cpp:148`), so a second
/// sign still reaches `std::from_chars` and fails.
#[test]
fn spectra_ref_run_index_is_an_int32_with_at_most_one_leading_plus() {
    for input in [
        "ms_run[2147483648]:scan=1",
        "ms_run[3000000000]:scan=1",
        "ms_run[++5]:scan=1",
        "ms_run[+-5]:scan=1",
    ] {
        let result: Result<MzTabSpectraRef, Error> = parse_cell(input);
        assert!(
            matches!(result, Err(Error::Parse { .. })),
            "{input:?} must be a conversion error"
        );
    }
    // One '+' is skipped, and the largest Int32 is still a legal index.
    let parsed: MzTabSpectraRef = parse_cell("ms_run[+5]:scan=1").unwrap();
    assert_eq!(parsed.ms_file(), 5);
    let parsed: MzTabSpectraRef = parse_cell("ms_run[2147483647]:scan=1").unwrap();
    assert_eq!(parsed.ms_file(), 2_147_483_647);
}

// ---------------------------------------------------------------------------
// Cell vocabulary: MzTabModification and MzTabModificationList
// ---------------------------------------------------------------------------

#[test]
fn modification_without_positions_is_a_bare_identifier() {
    let mut modification = MzTabModification::default();
    assert!(modification.is_null());
    assert_eq!(modification.to_cell_string().unwrap(), "null");
    modification.from_cell_string("UNIMOD:35").unwrap();
    assert!(modification.positions_and_parameters().is_empty());
    assert_eq!(modification.mod_or_subst_identifier().get(), "UNIMOD:35");
    assert_eq!(modification.to_cell_string().unwrap(), "UNIMOD:35");
    modification.set_null(true);
    assert!(modification.is_null());
}

#[test]
fn modification_positions_render_as_decimal_digits() {
    // Numeric string append in the source uses StringUtils.h's decimal
    // overload. Keep decimal spelling and parse/write round-trip coverage.
    let mut modification = MzTabModification::default();
    modification.from_cell_string("3-UNIMOD:35").unwrap();
    assert_eq!(
        modification.positions_and_parameters().len(),
        1,
        "one position"
    );
    assert_eq!(modification.positions_and_parameters()[0].0, 3);
    assert!(modification.positions_and_parameters()[0].1.is_null());
    let rendered = modification.to_cell_string().unwrap();
    assert_eq!(rendered, "3-UNIMOD:35");
    assert!(
        !rendered.bytes().any(|b| b < 0x20),
        "no control byte reaches the cell"
    );
    let round_trip: MzTabModification = parse_cell(&rendered).unwrap();
    assert_eq!(round_trip, modification);
}

#[test]
fn modification_positions_carry_optional_parameters() {
    let cell = "3|4[MS, MS:1001876, modification probability, 0.8]-UNIMOD:35";
    let modification: MzTabModification = parse_cell(cell).unwrap();
    let positions = modification.positions_and_parameters();
    assert_eq!(positions.len(), 2);
    assert_eq!(positions[0].0, 3);
    assert!(positions[0].1.is_null());
    assert_eq!(positions[1].0, 4);
    assert_eq!(positions[1].1.accession(), "MS:1001876");
    assert_eq!(positions[1].1.value(), "0.8");
    assert_eq!(modification.to_cell_string().unwrap(), cell);
}

#[test]
fn modification_rendering_requires_an_identifier() {
    let mut modification = MzTabModification::default();
    modification.set_positions_and_parameters(vec![(7, MzTabParameter::default())]);
    assert!(matches!(
        modification.to_cell_string(),
        Err(Error::MissingInformation(_))
    ));
    modification.set_modification_identifier(text("CHEMMOD:15.994915"));
    assert_eq!(
        modification.to_cell_string().unwrap(),
        "7-CHEMMOD:15.994915"
    );
}

#[test]
fn modification_cannot_read_back_a_negative_chemmod_delta() {
    // getModificationIdentifier_ (MzTab.cpp:1528) writes "CHEMMOD:" plus the
    // delta mono mass, which is negative for a loss. fromCellString splits on
    // '-' and requires exactly two fields, so neither spelling survives.
    let result: Result<MzTabModification, Error> = parse_cell("CHEMMOD:-18.010565");
    assert!(
        matches!(result, Err(Error::Parse { .. })),
        "the position half is 'CHEMMOD:' and is not an integer"
    );
    let result: Result<MzTabModification, Error> = parse_cell("8-CHEMMOD:-18.010565");
    assert!(
        matches!(result, Err(Error::Parse { .. })),
        "three '-'-separated fields"
    );
    // A positive delta does read back, because it carries no '-' of its own.
    let modification: MzTabModification = parse_cell("8-CHEMMOD:18.010565").unwrap();
    assert_eq!(
        modification.mod_or_subst_identifier().get(),
        "CHEMMOD:18.010565"
    );
    assert_eq!(modification.positions_and_parameters()[0].0, 8);
    assert_eq!(
        modification.to_cell_string().unwrap(),
        "8-CHEMMOD:18.010565"
    );
    // The negative spelling that does parse is misread: the identifier keeps
    // only the digits after the '-' and the CV prefix becomes a position.
    let mut modification = MzTabModification::default();
    let misread = modification.from_cell_string("1-CHEMMOD:18.010565");
    assert!(misread.is_ok());
    assert_eq!(
        modification.mod_or_subst_identifier().get(),
        "CHEMMOD:18.010565"
    );
}

#[test]
fn modification_position_must_be_a_non_negative_integer() {
    for input in ["x-UNIMOD:35", "-1|2-UNIMOD:35", "日本語-UNIMOD:35"] {
        let result: Result<MzTabModification, Error> = parse_cell(input);
        assert!(
            matches!(result, Err(Error::Parse { .. })),
            "{input:?} must be a conversion error"
        );
    }
    // Position 0 is the N-terminus and is legal.
    let modification: MzTabModification = parse_cell("0-UNIMOD:1").unwrap();
    assert_eq!(modification.positions_and_parameters()[0].0, 0);
}

/// Positions go through `StringUtils::toInt32` in the source, so the same
/// 32-bit range and single-`+` rules hold as for a `spectra_ref` index.
#[test]
fn modification_position_is_an_int32_with_at_most_one_leading_plus() {
    for input in [
        "2147483648-UNIMOD:35",
        "3000000000-UNIMOD:35",
        "++5-UNIMOD:35",
        "1|++2-UNIMOD:35",
    ] {
        let result: Result<MzTabModification, Error> = parse_cell(input);
        assert!(
            matches!(result, Err(Error::Parse { .. })),
            "{input:?} must be a conversion error"
        );
    }
    let modification: MzTabModification = parse_cell("+5-UNIMOD:35").unwrap();
    assert_eq!(modification.positions_and_parameters()[0].0, 5);
    let modification: MzTabModification = parse_cell("2147483647-UNIMOD:35").unwrap();
    assert_eq!(
        modification.positions_and_parameters()[0].0,
        2_147_483_647_usize
    );
}

#[test]
fn modification_list_splits_on_commas_outside_brackets() {
    let mut list = MzTabModificationList::default();
    assert!(list.is_null());
    assert_eq!(list.to_cell_string().unwrap(), "null");
    list.from_cell_string("UNIMOD:35,UNIMOD:4").unwrap();
    assert_eq!(list.get().len(), 2);
    assert_eq!(list.to_cell_string().unwrap(), "UNIMOD:35,UNIMOD:4");
    // A parameter's own commas do not split the list.
    let mut list = MzTabModificationList::default();
    list.from_cell_string("3|4[a,b,,v]-mod:123").unwrap();
    assert_eq!(list.get().len(), 1);
    let positions = list.get()[0].positions_and_parameters();
    assert_eq!(positions.len(), 2);
    assert_eq!(positions[1].1.cv_label(), "a");
    assert_eq!(positions[1].1.value(), "v");
    assert_eq!(list.get()[0].mod_or_subst_identifier().get(), "mod:123");
}

#[test]
fn modification_list_splits_inside_a_quoted_parameter_too() {
    // MzTabBase-style protection requires being outside quotes AND inside a
    // bracket (MzTab.cpp:263), so the source's own worked example is split at
    // the comma inside the quoted text, contrary to the comment above the loop.
    let cell = "3|4[a,b,,v]|8[,,\"blabla, [bla]\",v],1|2|3[a,b,,v]-mod:123";
    let mut list = MzTabModificationList::default();
    list.from_cell_string(cell).unwrap();
    assert_eq!(list.get().len(), 3, "not the two the comment implies");
    assert_eq!(
        list.get()[0].mod_or_subst_identifier().get(),
        "3|4[a,b,,v]|8[,,\"blabla"
    );
    assert_eq!(list.get()[1].mod_or_subst_identifier().get(), "[bla]\",v]");
    let last = &list.get()[2];
    assert_eq!(last.positions_and_parameters().len(), 3);
    assert_eq!(last.mod_or_subst_identifier().get(), "mod:123");
}

#[test]
fn modification_list_rendering_fails_on_the_first_bad_entry() {
    let mut incomplete = MzTabModification::default();
    incomplete.set_positions_and_parameters(vec![(1, MzTabParameter::default())]);
    let mut list = MzTabModificationList::default();
    list.set(vec![parse_cell("UNIMOD:35").unwrap(), incomplete]);
    assert!(matches!(
        list.to_cell_string(),
        Err(Error::MissingInformation(_))
    ));
    assert!(matches!(
        list.write_cell(),
        Err(Error::MissingInformation(_))
    ));
}

// ---------------------------------------------------------------------------
// The MzTabCell trait
// ---------------------------------------------------------------------------

#[test]
fn every_cell_type_drives_through_the_shared_trait() {
    fn round_trip<T: Default + MzTabCell + PartialEq + std::fmt::Debug>(rendered: &str) {
        let mut cell = T::default();
        cell.read_cell(rendered).unwrap();
        assert_eq!(cell.write_cell().unwrap(), rendered);
        let parsed: T = parse_cell(rendered).unwrap();
        assert_eq!(parsed, cell);
    }
    round_trip::<MzTabDouble>("1.5");
    round_trip::<MzTabDouble>("NaN");
    round_trip::<MzTabInteger>("-4");
    round_trip::<MzTabBoolean>("1");
    round_trip::<MzTabString>("text");
    round_trip::<MzTabParameter>("[MS, MS:1, a, b]");
    round_trip::<MzTabParameterList>("[MS, MS:1, a, b]|[MS, MS:2, c, d]");
    round_trip::<MzTabStringList>("a|b");
    round_trip::<MzTabIntegerList>("1,2");
    round_trip::<MzTabDoubleList>("1.5|2.5");
    round_trip::<MzTabSpectraRef>("ms_run[1]:scan=1");
    round_trip::<MzTabModification>("2-UNIMOD:35");
    round_trip::<MzTabModificationList>("UNIMOD:35,UNIMOD:4");

    // set_null through the trait, for a type whose null is an empty payload.
    let mut list: MzTabStringList = parse_cell("a|b").unwrap();
    MzTabCell::set_null(&mut list, true);
    assert!(MzTabCell::is_null(&list));
    MzTabCell::set_null(&mut list, false);
    assert!(
        MzTabCell::is_null(&list),
        "there is no value to restore, as the source"
    );
    // And for a type that carries an explicit state.
    let mut cell: MzTabDouble = parse_cell("1.5").unwrap();
    MzTabCell::set_null(&mut cell, true);
    assert!(MzTabCell::is_null(&cell));
    MzTabCell::set_null(&mut cell, false);
    assert!(!MzTabCell::is_null(&cell));
}

#[test]
fn cell_state_spellings() {
    assert_eq!(MzTabCellState::default(), MzTabCellState::Null);
    assert_eq!(MzTabCellState::Null.as_str(), "null");
    assert_eq!(MzTabCellState::NaN.as_str(), "NaN");
    assert_eq!(MzTabCellState::Inf.as_str(), "Inf");
    assert_eq!(MzTabCellState::Default.as_str(), "");
}

// ---------------------------------------------------------------------------
// Section rows
// ---------------------------------------------------------------------------

#[test]
fn protein_row_default_uses_comma_separated_lists_and_nucleic_acid_does_not() {
    // MzTab.cpp:294 sets both separators to ',' "because '|' can be used for
    // go terms and protein accessions". MzTabNucleicAcidSectionRow declares no
    // constructor, so it keeps '|'.
    let protein = MzTabProteinSectionRow::default();
    assert_eq!(protein.go_terms.separator(), ',');
    assert_eq!(protein.ambiguity_members.separator(), ',');
    let nucleic_acid = MzTabNucleicAcidSectionRow::default();
    assert_eq!(nucleic_acid.go_terms.separator(), '|');
    assert_eq!(nucleic_acid.ambiguity_members.separator(), '|');

    let mut protein = MzTabProteinSectionRow::default();
    protein
        .ambiguity_members
        .from_cell_string("P1|A,P2")
        .unwrap();
    assert_eq!(protein.ambiguity_members.get().len(), 2);
    assert_eq!(protein.ambiguity_members.get()[0].get(), "P1|A");
}

#[test]
fn row_order_follows_the_source_comparators() {
    let mut first = MzTabProteinSectionRow::default();
    first.accession.set("P2");
    let mut second = MzTabProteinSectionRow::default();
    second.accession.set("P10");
    // Text order, not numeric: "P10" < "P2".
    assert_eq!(first.row_order(&second), std::cmp::Ordering::Greater);

    let mut a = MzTabPeptideSectionRow::default();
    a.sequence.set("AAA");
    a.accession.set("P2");
    let mut b = MzTabPeptideSectionRow::default();
    b.sequence.set("AAA");
    b.accession.set("P1");
    assert_eq!(a.row_order(&b), std::cmp::Ordering::Greater);

    let mut a = MzTabPSMSectionRow::default();
    a.sequence.set("PEPTIDE");
    a.spectra_ref = MzTabSpectraRef::new(1, "scan=2").unwrap();
    a.psm_id.set(99);
    let mut b = MzTabPSMSectionRow::default();
    b.sequence.set("PEPTIDE");
    b.spectra_ref = MzTabSpectraRef::new(1, "scan=1").unwrap();
    b.psm_id.set(1);
    assert_eq!(
        a.row_order(&b),
        std::cmp::Ordering::Greater,
        "PSM_ID takes no part in the order"
    );

    let mut a = MzTabNucleicAcidSectionRow::default();
    a.accession.set("N1");
    let mut b = MzTabNucleicAcidSectionRow::default();
    b.accession.set("N2");
    assert_eq!(a.row_order(&b), std::cmp::Ordering::Less);

    let mut a = MzTabOSMSectionRow::default();
    a.sequence.set("ACGU");
    a.spectra_ref = MzTabSpectraRef::new(1, "scan=1").unwrap();
    let mut b = MzTabOSMSectionRow::default();
    b.sequence.set("ACGU");
    b.spectra_ref = MzTabSpectraRef::new(2, "scan=1").unwrap();
    assert_eq!(a.row_order(&b), std::cmp::Ordering::Less);
}

#[test]
fn oligonucleotide_row_order_reports_an_unset_position() {
    let mut a = MzTabOligonucleotideSectionRow::default();
    a.sequence.set("ACGU");
    a.start.set(1);
    a.end.set(4);
    let mut b = a.clone();
    b.start.set(2);
    assert_eq!(a.row_order(&b).unwrap(), std::cmp::Ordering::Less);
    // The source comparator calls MzTabInteger::get(), which throws for a null
    // cell; a default row therefore cannot be ordered.
    let unset = MzTabOligonucleotideSectionRow::default();
    assert!(matches!(
        a.row_order(&unset),
        Err(Error::MissingInformation(_))
    ));
}

#[test]
fn psm_row_pep_evidence_fills_five_cells_from_the_evidence_list() {
    let mut row = MzTabPSMSectionRow::default();
    row.accession.set("stale");
    let evidences = vec![
        PeptideEvidence {
            protein_accession: "P1".to_owned(),
            start: Some(9),
            end: Some(20),
            aa_before: FlankingResidue::Residue('K'),
            aa_after: FlankingResidue::Residue('R'),
        },
        PeptideEvidence {
            protein_accession: "P2".to_owned(),
            start: None,
            end: None,
            aa_before: FlankingResidue::NTerminus,
            aa_after: FlankingResidue::Unknown,
        },
    ];
    row.add_pep_evidence_to_rows(&evidences).unwrap();
    // MzTab counts from one, PeptideEvidence from zero.
    assert_eq!(row.start.get(), "10,null");
    assert_eq!(row.end.get(), "21,null");
    assert_eq!(row.pre.get(), "K,-");
    assert_eq!(row.post.get(), "R,null");
    assert_eq!(row.accession.get(), "P1,P2");
}

#[test]
fn psm_row_pep_evidence_clears_four_cells_but_keeps_the_accession() {
    let mut row = MzTabPSMSectionRow::default();
    row.accession.set("P9");
    row.pre.set("K");
    row.post.set("R");
    row.start.set("1");
    row.end.set("9");
    row.add_pep_evidence_to_rows(&[]).unwrap();
    assert!(row.pre.is_null());
    assert!(row.post.is_null());
    assert!(row.start.is_null());
    assert!(row.end.is_null());
    assert_eq!(
        row.accession.get(),
        "P9",
        "the source's early return does not touch the accession"
    );
    // A single unknown-everything evidence renders four "null" texts, which
    // MzTabString reads back as the null cell.
    let mut row = MzTabPSMSectionRow::default();
    row.add_pep_evidence_to_rows(&[PeptideEvidence::default()])
        .unwrap();
    assert!(row.pre.is_null());
    assert!(row.start.is_null());
    assert!(row.accession.is_null());
}

#[test]
fn optional_column_names_deduplicates_across_every_section_kind() {
    let mut document = MzTab::default();
    let mut protein = MzTabProteinSectionRow::default();
    protein
        .opt
        .push(MzTabOptionalColumnEntry::new("opt_global_b", text("1")));
    protein
        .opt
        .push(MzTabOptionalColumnEntry::new("opt_global_a", text("2")));
    let mut protein_two = MzTabProteinSectionRow::default();
    protein_two
        .opt
        .push(MzTabOptionalColumnEntry::new("opt_global_a", text("3")));
    protein_two
        .opt
        .push(MzTabOptionalColumnEntry::new("opt_global_c", text("4")));
    document.protein_data = vec![protein, protein_two];
    assert_eq!(
        document.protein_optional_column_names().unwrap(),
        vec![
            "opt_global_b".to_owned(),
            "opt_global_a".to_owned(),
            "opt_global_c".to_owned()
        ],
        "first-occurrence order, not sorted"
    );

    let mut peptide = MzTabPeptideSectionRow::default();
    peptide
        .opt
        .push(MzTabOptionalColumnEntry::new("opt_global_p", text("1")));
    document.peptide_data = vec![peptide];
    let mut small = MzTabSmallMoleculeSectionRow::default();
    small
        .opt
        .push(MzTabOptionalColumnEntry::new("opt_global_s", text("1")));
    document.small_molecule_data = vec![small];
    let mut nucleic = MzTabNucleicAcidSectionRow::default();
    nucleic
        .opt
        .push(MzTabOptionalColumnEntry::new("opt_global_n", text("1")));
    document.nucleic_acid_data = vec![nucleic];
    let mut oligo = MzTabOligonucleotideSectionRow::default();
    oligo
        .opt
        .push(MzTabOptionalColumnEntry::new("opt_global_o", text("1")));
    document.oligonucleotide_data = vec![oligo];
    let mut osm = MzTabOSMSectionRow::default();
    osm.opt
        .push(MzTabOptionalColumnEntry::new("opt_global_m", text("1")));
    document.osm_data = vec![osm];

    assert_eq!(
        document.peptide_optional_column_names().unwrap(),
        vec!["opt_global_p".to_owned()]
    );
    assert_eq!(
        document.small_molecule_optional_column_names().unwrap(),
        vec!["opt_global_s".to_owned()]
    );
    assert_eq!(
        document.nucleic_acid_optional_column_names().unwrap(),
        vec!["opt_global_n".to_owned()]
    );
    assert_eq!(
        document.oligonucleotide_optional_column_names().unwrap(),
        vec!["opt_global_o".to_owned()]
    );
    assert_eq!(
        document.osm_optional_column_names().unwrap(),
        vec!["opt_global_m".to_owned()]
    );
    assert!(document.psm_optional_column_names().unwrap().is_empty());
    // The free function works for any row type through the trait.
    assert_eq!(
        optional_column_names(&document.protein_data).unwrap().len(),
        3
    );
}

#[test]
fn optional_column_names_refuses_past_the_column_ceiling() {
    let mut row = MzTabPSMSectionRow::default();
    for index in 0..=MzTab::MAX_OPTIONAL_COLUMNS {
        row.opt.push(MzTabOptionalColumnEntry::new(
            format!("opt_global_{index}"),
            MzTabString::default(),
        ));
    }
    let rows = vec![row];
    assert!(matches!(
        optional_column_names(&rows),
        Err(Error::InvalidValue(_))
    ));
}

// ---------------------------------------------------------------------------
// The document
// ---------------------------------------------------------------------------

#[test]
fn number_of_psms_counts_distinct_ids_and_reports_an_unset_one() {
    let mut document = MzTab::default();
    for id in [1, 1, 2, 7] {
        let mut row = MzTabPSMSectionRow::default();
        row.psm_id.set(id);
        document.psm_data.push(row);
    }
    assert_eq!(document.psm_data.len(), 4);
    assert_eq!(document.number_of_psms().unwrap(), 3);
    // The source's @note says it relies on PSM_ID being set; MzTabInteger::get
    // throws otherwise.
    document.psm_data.push(MzTabPSMSectionRow::default());
    assert!(matches!(
        document.number_of_psms(),
        Err(Error::MissingInformation(_))
    ));
    assert_eq!(MzTab::default().number_of_psms().unwrap(), 0);
}

#[test]
fn document_accessors_and_setters_reach_every_section() {
    let mut document = MzTab::default();
    let mut meta = document.meta_data().clone();
    meta.title.set("a title");
    meta.mz_tab_mode.set("Complete");
    meta.mz_tab_type.set("Identification");
    document.set_meta_data(meta);
    assert_eq!(document.meta_data().title.get(), "a title");

    document.set_protein_section_rows(vec![MzTabProteinSectionRow::default()]);
    document.set_peptide_section_rows(vec![MzTabPeptideSectionRow::default(); 2]);
    document.set_psm_section_rows(vec![MzTabPSMSectionRow::default(); 3]);
    document.set_small_molecule_section_rows(vec![MzTabSmallMoleculeSectionRow::default(); 4]);
    document.set_nucleic_acid_section_rows(vec![MzTabNucleicAcidSectionRow::default(); 5]);
    document.set_oligonucleotide_section_rows(vec![MzTabOligonucleotideSectionRow::default(); 6]);
    document.set_osm_section_rows(vec![MzTabOSMSectionRow::default(); 7]);
    document.set_empty_rows(vec![4, 11]);
    let mut comments = BTreeMap::new();
    comments.insert(2_usize, "a comment".to_owned());
    document.set_comment_rows(comments);

    assert_eq!(document.protein_section_rows().len(), 1);
    assert_eq!(document.peptide_section_rows().len(), 2);
    assert_eq!(document.psm_section_rows().len(), 3);
    assert_eq!(document.small_molecule_section_rows().len(), 4);
    assert_eq!(document.nucleic_acid_section_rows().len(), 5);
    assert_eq!(document.oligonucleotide_section_rows().len(), 6);
    assert_eq!(document.osm_section_rows().len(), 7);
    assert_eq!(document.empty_rows(), &[4, 11]);
    assert_eq!(document.comment_rows()[&2], "a comment");
}

#[test]
fn metadata_records_every_indexed_key_in_index_order() {
    let mut meta = MzTab::default().meta_data;
    meta.ms_run.insert(
        2,
        openms::format::mztab::MzTabMSRunMetaData {
            location: text("file_b.mzML"),
            ..Default::default()
        },
    );
    meta.ms_run.insert(
        1,
        openms::format::mztab::MzTabMSRunMetaData {
            location: text("file_a.mzML"),
            ..Default::default()
        },
    );
    let locations: Vec<&str> = meta.ms_run.values().map(|run| run.location.get()).collect();
    assert_eq!(locations, vec!["file_a.mzML", "file_b.mzML"]);

    let mut software = MzTabSoftwareMetaData {
        software: MzTabParameter::parse("[MS, MS:1000752, TOPP software, ]").unwrap(),
        ..Default::default()
    };
    software.setting.insert(1, text("db=target.fasta"));
    meta.software.insert(1, software);
    assert_eq!(meta.software[&1].setting[&1].get(), "db=target.fasta");

    let mut sample = MzTabSampleMetaData::default();
    sample.description.set("a sample");
    sample
        .species
        .insert(1, MzTabParameter::parse("[NEWT, 9606, human, ]").unwrap());
    sample.tissue.insert(
        1,
        MzTabParameter::parse("[BTO, BTO:0000759, liver, ]").unwrap(),
    );
    sample.cell_type.insert(
        1,
        MzTabParameter::parse("[CL, CL:0000182, hepatocyte, ]").unwrap(),
    );
    sample.disease.insert(
        1,
        MzTabParameter::parse("[DOID, DOID:684, carcinoma, ]").unwrap(),
    );
    sample.custom.insert(
        1,
        MzTabParameter::parse("[, , description, free text]").unwrap(),
    );
    meta.sample.insert(1, sample);
    assert_eq!(meta.sample[&1].species[&1].name(), "human");
    assert_eq!(meta.sample[&1].custom[&1].value(), "free text");

    let mut instrument = MzTabInstrumentMetaData {
        name: MzTabParameter::parse("[MS, MS:1000449, LTQ Orbitrap, ]").unwrap(),
        source: MzTabParameter::parse("[MS, MS:1000073, ESI, ]").unwrap(),
        detector: MzTabParameter::parse("[MS, MS:1000253, EM, ]").unwrap(),
        ..Default::default()
    };
    instrument
        .analyzer
        .insert(1, MzTabParameter::parse("[MS, MS:1000291, LIT, ]").unwrap());
    meta.instrument.insert(1, instrument);
    assert_eq!(meta.instrument[&1].analyzer[&1].accession(), "MS:1000291");

    meta.contact.insert(
        1,
        MzTabContactMetaData {
            name: text("A Scientist"),
            affiliation: text("A University"),
            email: text("a@example.org"),
        },
    );
    assert_eq!(meta.contact[&1].email.get(), "a@example.org");

    meta.cv.insert(
        1,
        MzTabCVMetaData {
            label: text("MS"),
            full_name: text("PSI-MS"),
            version: text("4.1.0"),
            url: text("https://example.org/psi-ms.obo"),
        },
    );
    assert_eq!(meta.cv[&1].version.get(), "4.1.0");

    meta.assay.insert(
        1,
        openms::format::mztab::MzTabAssayMetaData {
            sample_ref: text("sample[1]"),
            ms_run_ref: vec![1, 2],
            ..Default::default()
        },
    );
    assert_eq!(meta.assay[&1].ms_run_ref, vec![1, 2]);
    meta.study_variable.insert(
        1,
        openms::format::mztab::MzTabStudyVariableMetaData {
            assay_refs: vec![1],
            sample_refs: vec![1],
            description: text("control"),
        },
    );
    assert_eq!(meta.study_variable[&1].assay_refs, vec![1]);

    meta.colunit_protein
        .push("protein_abundance_assay[1]=[UO, UO:0000187, percent, ]".to_owned());
    assert_eq!(meta.colunit_protein.len(), 1);
    meta.publication.insert(1, text("pubmed:1234567"));
    meta.uri.insert(1, text("https://example.org/run"));
    meta.custom
        .insert(1, MzTabParameter::parse("[, , operator, me]").unwrap());
    assert_eq!(meta.publication[&1].get(), "pubmed:1234567");
    assert_eq!(meta.uri[&1].get(), "https://example.org/run");
    assert_eq!(meta.custom[&1].value(), "me");

    for scores in [
        &mut meta.protein_search_engine_score,
        &mut meta.peptide_search_engine_score,
        &mut meta.psm_search_engine_score,
        &mut meta.smallmolecule_search_engine_score,
        &mut meta.nucleic_acid_search_engine_score,
        &mut meta.oligonucleotide_search_engine_score,
        &mut meta.osm_search_engine_score,
    ] {
        scores.insert(
            1,
            MzTabParameter::parse("[MS, MS:1002015, Percolator q-value, ]").unwrap(),
        );
    }
    assert_eq!(meta.osm_search_engine_score[&1].accession(), "MS:1002015");

    meta.sample_processing.insert(1, {
        let mut list = MzTabParameterList::default();
        list.from_cell_string("[MS, MS:1000544, Conversion to mzML, ]")
            .unwrap();
        list
    });
    assert_eq!(meta.sample_processing[&1].get().len(), 1);
    meta.false_discovery_rate
        .from_cell_string("[MS, MS:1001364, pep:global FDR, 0.01]")
        .unwrap();
    assert_eq!(meta.false_discovery_rate.get()[0].value(), "0.01");
    meta.quantification_method =
        MzTabParameter::parse("[MS, MS:1001834, LC-MS label-free, ]").unwrap();
    meta.protein_quantification_unit =
        MzTabParameter::parse("[PRIDE, PRIDE:0000393, Relative, ]").unwrap();
    meta.peptide_quantification_unit = meta.protein_quantification_unit.clone();
    meta.small_molecule_quantification_unit = meta.protein_quantification_unit.clone();
    assert_eq!(meta.quantification_method.accession(), "MS:1001834");
    assert_eq!(meta.small_molecule_quantification_unit.name(), "Relative");
}

#[test]
fn every_row_kind_carries_its_full_column_set() {
    let mut protein = MzTabProteinSectionRow::default();
    protein.accession.set("P1");
    protein.description.set("a protein");
    protein.taxid.set(9606);
    protein.species.set("Homo sapiens");
    protein.database.set("target.fasta");
    protein.database_version.set("2026-01");
    protein
        .search_engine
        .from_cell_string("[MS, MS:1001476, Percolator, ]")
        .unwrap();
    protein
        .best_search_engine_score
        .insert(1, MzTabDouble::new(0.99));
    protein
        .search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(2, MzTabDouble::new(0.98));
    protein.reliability.set(1);
    protein.num_psms_ms_run.insert(1, MzTabInteger::new(12));
    protein
        .num_peptides_distinct_ms_run
        .insert(1, MzTabInteger::new(5));
    protein
        .num_peptides_unique_ms_run
        .insert(1, MzTabInteger::new(4));
    protein.modifications.from_cell_string("UNIMOD:35").unwrap();
    protein.uri.set("https://example.org/P1");
    protein.coverage.set(0.42);
    protein
        .protein_abundance_assay
        .insert(1, MzTabDouble::new(10.0));
    protein
        .protein_abundance_study_variable
        .insert(1, MzTabDouble::default());
    protein
        .protein_abundance_stdev_study_variable
        .insert(1, MzTabDouble::new(0.5));
    protein
        .protein_abundance_std_error_study_variable
        .insert(1, MzTabDouble::new(0.25));
    assert_eq!(
        protein.search_engine_score_ms_run[&1][&2].get().unwrap(),
        0.98
    );
    assert!(protein.protein_abundance_study_variable[&1].is_null());
    assert_eq!(protein.coverage.to_cell_string(), "0.42");
    assert_eq!(protein.taxid.to_cell_string(), "9606");

    let mut peptide = MzTabPeptideSectionRow::default();
    peptide.sequence.set("PEPTIDER");
    peptide.unique.set(true);
    peptide.retention_time.from_cell_string("100.5").unwrap();
    peptide
        .retention_time_window
        .from_cell_string("99.0|101.0")
        .unwrap();
    peptide.charge.set(2);
    peptide.mass_to_charge.set(500.25);
    peptide.spectra_ref = MzTabSpectraRef::new(1, "scan=3").unwrap();
    peptide
        .peptide_abundance_assay
        .insert(1, MzTabDouble::new(5.0));
    peptide
        .peptide_abundance_study_variable
        .insert(1, MzTabDouble::default());
    peptide
        .peptide_abundance_stdev_study_variable
        .insert(1, MzTabDouble::default());
    peptide
        .peptide_abundance_std_error_study_variable
        .insert(1, MzTabDouble::default());
    assert_eq!(peptide.retention_time_window.get().len(), 2);
    assert_eq!(peptide.mass_to_charge.to_cell_string(), "500.25");

    let mut small = MzTabSmallMoleculeSectionRow::default();
    small
        .identifier
        .from_cell_string("CID:5793|HMDB:HMDB0000122")
        .unwrap();
    small.chemical_formula.set("C6H12O6");
    small.smiles.set("OCC1OC(O)C(O)C(O)C1O");
    small.inchi_key.set("WQZGKKKJIJFFOK-GASJEMHNSA-N");
    small.description.set("glucose");
    small.exp_mass_to_charge.set(181.0707);
    small.calc_mass_to_charge.set(181.0707);
    small.charge.set(1);
    small.retention_time.from_cell_string("42.0").unwrap();
    small.taxid.set(9606);
    small.species.set("Homo sapiens");
    small.database.set("HMDB");
    small.database_version.set("5.0");
    small.reliability.set(2);
    small.uri.set("https://example.org/glucose");
    small.spectra_ref = MzTabSpectraRef::new(1, "scan=9").unwrap();
    small
        .search_engine
        .from_cell_string("[MS, MS:1002889, AccurateMassSearch, ]")
        .unwrap();
    small
        .best_search_engine_score
        .insert(1, MzTabDouble::new(0.8));
    small
        .search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(1, MzTabDouble::new(0.8));
    // The SML section's `modifications` is a plain string cell, not a list.
    small.modifications.set("some annotation");
    small
        .smallmolecule_abundance_assay
        .insert(1, MzTabDouble::new(3.0));
    small
        .smallmolecule_abundance_study_variable
        .insert(1, MzTabDouble::new(3.0));
    small
        .smallmolecule_abundance_stdev_study_variable
        .insert(1, MzTabDouble::default());
    small
        .smallmolecule_abundance_std_error_study_variable
        .insert(1, MzTabDouble::default());
    assert_eq!(small.identifier.get().len(), 2);
    assert_eq!(small.modifications.get(), "some annotation");

    let mut nucleic = MzTabNucleicAcidSectionRow::default();
    nucleic.accession.set("N1");
    nucleic.description.set("a nucleic acid");
    nucleic.taxid.set(9606);
    nucleic.species.set("Homo sapiens");
    nucleic.database.set("rna.fasta");
    nucleic.database_version.set("1");
    nucleic
        .search_engine
        .from_cell_string("[MS, MS:1, NASE, ]")
        .unwrap();
    nucleic
        .best_search_engine_score
        .insert(1, MzTabDouble::new(1.0));
    nucleic
        .search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(1, MzTabDouble::new(1.0));
    nucleic.reliability.set(1);
    nucleic.num_osms_ms_run.insert(1, MzTabInteger::new(3));
    nucleic
        .num_oligos_distinct_ms_run
        .insert(1, MzTabInteger::new(2));
    nucleic
        .num_oligos_unique_ms_run
        .insert(1, MzTabInteger::new(1));
    nucleic.modifications.from_cell_string("UNIMOD:1").unwrap();
    nucleic.uri.set("https://example.org/N1");
    nucleic.coverage.set(0.5);
    assert_eq!(nucleic.num_osms_ms_run[&1].get().unwrap(), 3);

    let mut oligo = MzTabOligonucleotideSectionRow::default();
    oligo.sequence.set("ACGU");
    oligo.accession.set("N1");
    oligo.unique.set(false);
    oligo
        .search_engine
        .from_cell_string("[MS, MS:1, NASE, ]")
        .unwrap();
    oligo
        .best_search_engine_score
        .insert(1, MzTabDouble::new(0.5));
    oligo
        .search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(1, MzTabDouble::new(0.5));
    oligo.reliability.set(3);
    oligo.modifications.from_cell_string("null").unwrap();
    oligo.retention_time.from_cell_string("10.0").unwrap();
    oligo
        .retention_time_window
        .from_cell_string("9.0|11.0")
        .unwrap();
    oligo.uri.set("https://example.org/oli");
    oligo.pre.set("G");
    oligo.post.set("C");
    oligo.start.set(3);
    oligo.end.set(6);
    assert_eq!(oligo.end.get().unwrap(), 6);
    assert!(oligo.modifications.is_null());

    let mut osm = MzTabOSMSectionRow::default();
    osm.sequence.set("ACGU");
    osm.search_engine
        .from_cell_string("[MS, MS:1, NASE, ]")
        .unwrap();
    osm.search_engine_score.insert(1, MzTabDouble::new(0.1));
    osm.reliability.set(1);
    osm.modifications.from_cell_string("UNIMOD:1").unwrap();
    osm.retention_time.from_cell_string("10.0").unwrap();
    osm.charge.set(-2);
    osm.exp_mass_to_charge.set(600.5);
    osm.calc_mass_to_charge.set(600.5);
    osm.uri.set("https://example.org/osm");
    osm.spectra_ref = MzTabSpectraRef::new(1, "scan=4").unwrap();
    assert_eq!(osm.charge.to_cell_string(), "-2");
    assert_eq!(osm.search_engine_score[&1].get().unwrap(), 0.1);
}

// ---------------------------------------------------------------------------
// Modification metadata generators
// ---------------------------------------------------------------------------

#[test]
fn modification_metadata_reports_the_unimod_accession_upper_cased() {
    let report = modification_metadata(&["Oxidation (M)".to_owned()]).unwrap();
    assert!(report.skipped.is_empty());
    assert_eq!(report.metadata.len(), 1);
    let entry = &report.metadata[&1];
    assert_eq!(entry.modification.cv_label(), "UNIMOD");
    assert_eq!(entry.modification.accession(), "UNIMOD:35");
    assert_eq!(entry.modification.name(), "Oxidation");
    assert_eq!(entry.site.get(), "M");
    assert_eq!(entry.position.get(), "Anywhere");
    assert_eq!(
        entry.modification.to_cell_string(),
        "[UNIMOD, UNIMOD:35, Oxidation, ]"
    );
}

#[test]
fn modification_metadata_skips_an_unknown_name_and_leaves_an_index_gap() {
    // The source catches every exception, logs "Skipping unknown residue
    // modification" and still increments its index, so the map has a hole.
    let names = vec![
        "No Such Modification".to_owned(),
        "Oxidation (M)".to_owned(),
    ];
    let report = modification_metadata(&names).unwrap();
    assert_eq!(report.skipped, vec!["No Such Modification".to_owned()]);
    assert_eq!(report.metadata.len(), 1);
    assert!(!report.metadata.contains_key(&1));
    assert_eq!(report.metadata[&2].modification.name(), "Oxidation");
    // A caller-owned registry replaces the source's global singleton.
    let report = modification_metadata_with(ModificationsDB::global(), &names).unwrap();
    assert_eq!(report.metadata.len(), 1);
    // Terminal specificity spellings.
    let report = modification_metadata(&["Acetyl (Protein N-term)".to_owned()]).unwrap();
    assert_eq!(report.metadata[&1].position.get(), "Protein N-term");
    assert_eq!(report.metadata[&1].site.get(), "X");
}

#[test]
fn empty_modification_lists_produce_the_specified_placeholders() {
    let fixed = fixed_modification_metadata(&[]).unwrap();
    assert_eq!(fixed.metadata.len(), 1);
    assert_eq!(fixed.metadata[&1].modification.cv_label(), "MS");
    assert_eq!(fixed.metadata[&1].modification.accession(), "MS:1002453");
    assert_eq!(
        fixed.metadata[&1].modification.name(),
        "No fixed modifications searched"
    );
    assert_eq!(
        fixed.metadata[&1].modification.to_cell_string(),
        "[MS, MS:1002453, No fixed modifications searched, ]"
    );
    assert!(fixed.metadata[&1].site.is_null());

    let variable = variable_modification_metadata(&[]).unwrap();
    assert_eq!(variable.metadata[&1].modification.accession(), "MS:1002454");
    assert_eq!(
        variable.metadata[&1].modification.name(),
        "No variable modifications searched"
    );

    // A non-empty list delegates to the generator.
    let variable = variable_modification_metadata(&["Oxidation (M)".to_owned()]).unwrap();
    assert_eq!(variable.metadata[&1].modification.accession(), "UNIMOD:35");
    let fixed = fixed_modification_metadata(&["Oxidation (M)".to_owned()]).unwrap();
    assert_eq!(fixed.metadata[&1].modification.accession(), "UNIMOD:35");
}

#[test]
fn modification_metadata_refuses_an_oversized_name_list() {
    let names: Vec<String> = (0..=openms::format::mztab::MAX_MODIFICATION_NAMES)
        .map(|i| format!("mod{i}"))
        .collect();
    assert!(matches!(
        modification_metadata(&names),
        Err(Error::InvalidValue(_))
    ));
}
