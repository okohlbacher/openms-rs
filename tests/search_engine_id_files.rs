// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Class tests for the search-engine identification adapters.
//!
//! Every expected value is transcribed from the upstream class test named at
//! each case, at the core pin (tier 2), and run against the unmodified
//! upstream fixtures described in `tests/data/search_engine_files_provenance.json`.
//! Cases the upstream tests do not cover are derived from the source's control
//! flow and say so.

use openms::Error;
use openms::chemistry::AASequence;
use openms::format::omssa_csv_file::OMSSACSVFile;
use openms::format::sequest_outfile::SequestOutfile;
use openms::identification::{FlankingResidue, PeptideHit, ProteinIdentification};
use openms::metadata::MetaValueData;
use std::collections::{BTreeMap, BTreeSet};
use std::io::Cursor;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

#[test]
fn omssa_csv_load_matches_the_class_test() {
    // OMSSACSVFile_test.cpp, section load: one identification (the protein
    // identification the source leaves untouched is not part of the Rust API).
    let ids = OMSSACSVFile::new()
        .load(data("OMSSACSVFile_test_1.csv"))
        .unwrap();
    assert_eq!(ids.len(), 1);
    // Derived from OMSSACSVFile.cpp: the three lines of spectrum 0 are three
    // hits with score type OMSSA, the P-value as score and charge 1.
    let id = &ids[0];
    assert_eq!(id.score_type, "OMSSA");
    assert_eq!(id.hits.len(), 3);
    for hit in &id.hits {
        assert_eq!(hit.sequence, AASequence::parse("MSHHWGYGK").unwrap());
        assert_eq!(hit.score, 1.3086e-08);
        assert_eq!(hit.charge, 1);
    }
}

#[test]
fn omssa_csv_refuses_an_unbalanced_quote() {
    // OMSSACSVFile_test.cpp: TEST_EXCEPTION(ParseError) on the unbalanced-quote
    // regression fixture (15 fields, 5 quoted commas).
    match OMSSACSVFile::new().load(data("OMSSACSVFile_test_unbalanced_quote.csv")) {
        Err(Error::Parse { line, message }) => {
            assert_eq!(line, 2);
            assert_eq!(message, "number of columns should be 14 in line 1");
        }
        other => panic!("expected a parse error, got {other:?}"),
    }
}

#[test]
fn omssa_csv_reads_shifted_columns_after_a_quoted_defline() {
    // Derived from OMSSACSVFile.cpp: two quoted commas shift the charge and
    // P-value two fields right; a new spectrum number starts a new
    // identification.
    let text = "Spectrum number, ...\n\
                0,,PEPTIDE,1,2,0,3,1,7,plain,,2,3,0.5\n\
                1,,PEPTIDER,1,2,0,3,1,8,\"a, b, c\",,3,4,0.25\n";
    let ids = OMSSACSVFile::new().load_reader(Cursor::new(text)).unwrap();
    assert_eq!(ids.len(), 2);
    assert_eq!(ids[0].hits[0].charge, 2);
    assert_eq!(ids[0].hits[0].score, 0.5);
    assert_eq!(ids[1].hits[0].charge, 3);
    assert_eq!(ids[1].hits[0].score, 0.25);
}

#[test]
fn omssa_csv_requires_the_spectrum_header() {
    // Derived from OMSSACSVFile.cpp: hasPrefix(line, "Spectrum") on the first
    // line, an empty file included.
    for text in ["", "number, ...\n"] {
        match OMSSACSVFile::new().load_reader(Cursor::new(text)) {
            Err(Error::Parse { line, message }) => {
                assert_eq!(line, 1);
                assert_eq!(message, "first line does not contain a description");
            }
            other => panic!("expected a parse error, got {other:?}"),
        }
    }
}

#[test]
fn omssa_csv_missing_file_is_an_io_error() {
    // Derived from OMSSACSVFile.cpp: FileNotFound when the stream fails.
    let error = OMSSACSVFile::new()
        .load(data("OMSSACSVFile_test_doesnotexist.csv"))
        .unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error:?}");
}

/// ClassTest `isRealSimilar` (`ClassTest.cpp:364`, testframework): absolute
/// difference at most 1e-5, or ratio at most 1 + 1e-5, with the source's zero
/// cases.
fn real_similar(actual: f64, expected: f64) -> bool {
    if actual.is_nan() || expected.is_nan() {
        return false;
    }
    let small = (actual - expected).abs() <= 1e-5;
    if actual == 0.0 {
        return expected == 0.0 || small;
    }
    if expected == 0.0 {
        return small;
    }
    let mut ratio = actual / expected;
    if ratio < 0.0 {
        return small;
    }
    if ratio < 1.0 {
        ratio = 1.0 / ratio;
    }
    ratio <= 1.0 + 1e-5 || small
}

fn float_meta(hit: &PeptideHit, name: &str) -> f64 {
    match hit.metadata.get(name).map(|value| value.data()) {
        Some(MetaValueData::Float(value)) => *value,
        other => panic!("meta value {name} is not a double: {other:?}"),
    }
}

fn string_meta<'a>(hit: &'a PeptideHit, name: &str) -> &'a str {
    match hit.metadata.get(name).map(|value| value.data()) {
        Some(MetaValueData::String(value)) => value,
        other => panic!("meta value {name} is not a string: {other:?}"),
    }
}

fn int_meta(hit: &PeptideHit, name: &str) -> i64 {
    match hit.metadata.get(name).map(|value| value.data()) {
        Some(MetaValueData::Integer(value)) => *value,
        other => panic!("meta value {name} is not an integer: {other:?}"),
    }
}

/// One expected Sequest hit of `SequestOutfile.out`, from the class test.
struct ExpectedSequestHit {
    score: f64,
    sequence: &'static str,
    before: FlankingResidue,
    after: FlankingResidue,
    rank: u32,
    rank_sp: &'static str,
    mh: f64,
    delta_cn: f64,
    xcorr: f64,
    sp: f64,
    sf: f64,
    ions: &'static str,
}

fn check_sequest_hit(hit: &PeptideHit, expected: &ExpectedSequestHit) {
    assert!(real_similar(hit.score, expected.score), "{}", hit.score);
    assert_eq!(hit.sequence.to_string(), expected.sequence);
    assert_eq!(hit.evidences[0].aa_before, expected.before);
    assert_eq!(hit.evidences[0].aa_after, expected.after);
    assert_eq!(hit.rank, expected.rank);
    assert_eq!(hit.charge, 3);
    assert_eq!(string_meta(hit, "RankSp"), expected.rank_sp);
    assert_eq!(int_meta(hit, "SequestId"), 0);
    assert!(real_similar(float_meta(hit, "MH"), expected.mh));
    assert!(real_similar(float_meta(hit, "DeltCn"), expected.delta_cn));
    assert!(real_similar(float_meta(hit, "XCorr"), expected.xcorr));
    assert!(real_similar(float_meta(hit, "Sp"), expected.sp));
    assert!(real_similar(float_meta(hit, "Sf"), expected.sf));
    assert_eq!(string_meta(hit, "Ions"), expected.ions);
}

// SequestOutfile_test.cpp writes the Sequest terminus marker as the raw
// character '-'; the Rust evidence stores it as the C-terminus marker.
const EXPECTED_SEQUEST_HITS: [ExpectedSequestHit; 4] = [
    ExpectedSequestHit {
        score: 0.05,
        sequence: "ETQAWSIATILETLYDL",
        before: FlankingResidue::Residue('C'),
        after: FlankingResidue::CTerminus,
        rank: 1,
        rank_sp: "1/80",
        mh: 1967.0013,
        delta_cn: 0.0,
        xcorr: 1.5789,
        sp: 310.3,
        sf: 0.05,
        ions: "18/64",
    },
    ExpectedSequestHit {
        score: 0.04,
        sequence: "QVLNPLLVLIFIPLFDL",
        before: FlankingResidue::Residue('M'),
        after: FlankingResidue::Residue('V'),
        rank: 2,
        rank_sp: "2/85",
        mh: 1967.1985,
        delta_cn: 0.0390,
        xcorr: 1.5173,
        sp: 308.3,
        sf: 0.04,
        ions: "19/64",
    },
    ExpectedSequestHit {
        score: 0.02,
        sequence: "WVELGPSVLAGVGVMVLLI",
        before: FlankingResidue::Residue('L'),
        after: FlankingResidue::Residue('P'),
        rank: 3,
        rank_sp: "3/117",
        mh: 1968.1244,
        delta_cn: 0.0501,
        xcorr: 1.4998,
        sp: 292.4,
        sf: 0.02,
        ions: "17/72",
    },
    ExpectedSequestHit {
        score: 0.14,
        sequence: "FDEITAMTGDGVNDAPALK",
        before: FlankingResidue::Residue('S'),
        after: FlankingResidue::Residue('K'),
        rank: 4,
        rank_sp: "4/1",
        mh: 1964.9275,
        delta_cn: 0.0627,
        xcorr: 1.4799,
        sp: 530.9,
        sf: 0.14,
        ions: "24/72",
    },
];

fn accessions(hit: &PeptideHit) -> Vec<&str> {
    hit.evidences
        .iter()
        .map(|evidence| evidence.protein_accession.as_str())
        .collect()
}

#[test]
fn sequest_load_exceptions_match_the_class_test() {
    // SequestOutfile_test.cpp, section load, "test exceptions".
    let file = SequestOutfile::new();
    let mut protein = ProteinIdentification::default();
    let error = file
        .load(data("a"), &mut protein, 0.01, &[], false)
        .unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error:?}");
    match file.load(data("SequestOutfile.out1"), &mut protein, 0.01, &[], false) {
        Err(Error::Parse { line, message }) => {
            assert_eq!(line, 16);
            assert_eq!(
                message,
                "Wrong number of columns in line 16! (11 present, should be 12)"
            );
        }
        other => panic!("expected a parse error, got {other:?}"),
    }
    for threshold in [2.0, -1.0] {
        assert!(matches!(
            file.load("", &mut protein, threshold, &[], false),
            Err(Error::InvalidValue(_))
        ));
    }
    assert_eq!(protein, ProteinIdentification::default());
}

#[test]
fn sequest_load_matches_the_class_test() {
    // SequestOutfile_test.cpp, section load, "test the actual program".
    let file = SequestOutfile::new();
    let mut protein = ProteinIdentification::default();
    assert_eq!(
        file.load(data("SequestOutfile2.out"), &mut protein, 1.0, &[], false)
            .unwrap(),
        None
    );

    let id = file
        .load(data("SequestOutfile.out"), &mut protein, 1.0, &[], false)
        .unwrap()
        .expect("one identification");
    assert_eq!(id.hits.len(), 4);
    assert_eq!(id.score_type, "SEQUEST");
    assert_eq!(id.identifier, "TurboSEQUEST_2004-03-16");
    assert!(real_similar(id.significance_threshold, 1.0));
    for (hit, expected) in id.hits.iter().zip(&EXPECTED_SEQUEST_HITS) {
        check_sequest_hit(hit, expected);
    }
    let set: BTreeSet<&str> = id.hits[0].protein_accessions();
    assert_eq!(
        set.into_iter().collect::<Vec<_>>(),
        ["2136928", "L10605", "P35574"]
    );
    assert_eq!(accessions(&id.hits[1]), ["P46029", "U32507"]);
    assert_eq!(accessions(&id.hits[2]), ["e148876"]);
    assert_eq!(
        accessions(&id.hits[3]),
        [
            "P20647", "P04192", "67962", "67961", "109166", "224621", "X02814", "J04703"
        ]
    );

    let id = file
        .load(
            data("SequestOutfile.out"),
            &mut protein,
            0.01,
            &[0.001, 0.01, 0.05, 0.5],
            false,
        )
        .unwrap()
        .expect("one identification");
    assert_eq!(id.score_type, "SEQUEST");
    assert_eq!(id.identifier, "TurboSEQUEST_2004-03-16");
    assert_eq!(id.hits.len(), 2);
    assert!(real_similar(id.significance_threshold, 0.01));
    for (hit, expected) in id.hits.iter().zip(&EXPECTED_SEQUEST_HITS) {
        check_sequest_hit(hit, expected);
    }
    assert_eq!(accessions(&id.hits[0]), ["P35574", "2136928", "L10605"]);
    assert_eq!(accessions(&id.hits[1]), ["P46029", "U32507"]);

    assert_eq!(protein.search_engine, "TurboSEQUEST");
    assert_eq!(protein.search_engine_version, "v.27 (rev. 12)");
    assert_eq!(protein.identifier, "TurboSEQUEST_2004-03-16");
}

#[test]
fn sequest_get_columns_matches_the_class_test() {
    // SequestOutfile_test.cpp, section getColumns.
    let file = SequestOutfile::new();
    let mut columns = vec![
        "1.",
        "1/80",
        "0",
        "1967.0013",
        "0.0000",
        "1.5789",
        "310.3",
        "0.05",
        "0",
        "18/64",
        "gi|544379|sp|P35574|GDE RABIT+2",
        "C.ETQAWSIATILETLYDL.-",
    ];
    assert_eq!(file.get_columns("", 12, 10), None);
    let line = "  1.   1/80          0 1967.0013  0.0000  1.5789   310.3 0.05    0 18/64 gi|544379|sp|P35574|GDE RABIT   +2   C.ETQAWSIATILETLYDL.-";
    assert_eq!(file.get_columns(line, 12, 10).unwrap(), columns);
    let line = "  1.   1/80          0 1967.0013  0.0000  1.5789   310.3 0.05    0 18/64 gi|544379|sp|P35574|GDE RABIT+2   C.ETQAWSIATILETLYDL.-";
    assert_eq!(file.get_columns(line, 12, 10).unwrap(), columns);
    let line = "  1.   1/80          0 1967.0013  0.0000  1.5789   310.3 0.05    0 18/64 gi|544379|sp|P35574|GDE RABIT   +X   C.ETQAWSIATILETLYDL.-";
    columns[10] = "gi|544379|sp|P35574|GDE RABIT +X";
    assert_eq!(file.get_columns(line, 12, 10).unwrap(), columns);
}

#[test]
fn sequest_get_sequences_matches_the_class_test() {
    // SequestOutfile_test.cpp, section getSequences; the source appends to
    // the caller's lists across the two calls, which the test does here.
    let file = SequestOutfile::new();
    assert!(matches!(
        file.get_sequences(data("a"), &BTreeMap::new()),
        Err(Error::Io(_))
    ));
    let ac_position_map: BTreeMap<String, usize> =
        [("P02666", 0), ("Q9CQV8", 1), ("Q5EEQ7", 2), ("P68509", 3)]
            .into_iter()
            .map(|(accession, position)| (accession.to_owned(), position))
            .collect();
    let sequences = [
        "MKVLILACLVALALARELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPPFLQPEVMGVSKVKEAMAPKHKEMPFPKYPVEPFTESQSLTLTDVENLHLPLPLLQSWMHQPHQPLPPTVMFPPQSVLSLSQSKVLPVPQKAVPYPQRDMPIQAFLLYQEPVLGPVRGPFPIIV",
        "TMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLILNATQAESKVFYLKMKGDYFRYLSEVASGENKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGDAGEGEN",
        "SAPPSLLVLYFGKKELRAMKVLILACLVALALARELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPP",
        "GDREQLLQRARLAEQAERYDDMASAMKAVTELNEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVCNDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEISKEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQLLRDNLTLWTSDQQDEEAGEGN",
    ];
    let first = file
        .get_sequences(data("Sequest_test.fasta"), &ac_position_map)
        .unwrap();
    assert_eq!(first.found.len(), 2);
    assert_eq!(first.not_found.len(), 2);
    assert_eq!(first.found[0], ("P68509".to_owned(), 3));
    assert_eq!(first.sequences[0], sequences[3]);
    assert_eq!(first.found[1], ("Q9CQV8".to_owned(), 1));
    assert_eq!(first.sequences[1], sequences[1]);

    let second = file
        .get_sequences(data("Sequest_test2.fasta"), &first.not_found)
        .unwrap();
    let found: Vec<_> = first.found.iter().chain(&second.found).cloned().collect();
    let found_sequences: Vec<_> = first.sequences.iter().chain(&second.sequences).collect();
    assert_eq!(found.len(), 4);
    assert_eq!(second.not_found.len(), 0);
    assert_eq!(found[2], ("P02666".to_owned(), 0));
    assert_eq!(found_sequences[2], sequences[0]);
    assert_eq!(found[3], ("Q5EEQ7".to_owned(), 2));
    assert_eq!(found_sequences[3], sequences[2]);
}

#[test]
fn sequest_accessions_match_the_class_test() {
    // SequestOutfile_test.cpp, section getACAndACType.
    let file = SequestOutfile::new();
    let cases = [
        (
            ">sp|P02666|CASB_BOVIN Beta-casein precursor - Bos taurus (Bovine).",
            "P02666",
            "SwissProt",
        ),
        (
            ">tr|Q5EEQ7|Q5EEQ7_BOVIN Beta-casein (Fragment) - Bos taurus (Bovine).",
            "Q5EEQ7",
            "SwissProt",
        ),
        ("gi|110174602|gb|DQ660451.1|", "DQ660451.1", "GenBank"),
        ("gi|1655698|emb|Y07752|VCPHEROPH", "Y07752", "EMBL"),
        ("gi|10038695|dbj|BAB12730|", "BAB12730", "DDBJ"),
        ("gi|9628804|ref|NP_043835|", "NP_043835", "NCBI"),
        ("gi|21362794|sp|P58858.0|", "P58858.0", "SwissProt"),
        ("gi|21362794|tr|P58858.0|", "P58858.0", "SwissProt"),
        ("gi|1619818|gnl|PID|d1013471|", "d1013471", "PID"),
        (
            "Q30DX2 Gamma-gliadin/LMW-glutenin chimera Ch7 (Fragment).",
            "Q30DX2",
            "SwissProt",
        ),
        (">P68509|1433F_BOVIN", "P68509", "SwissProt"),
        (">ACBLA (P68509) F_BOVIN", "P68509", "SwissProt"),
    ];
    for (line, accession, accession_type) in cases {
        assert_eq!(
            file.get_ac_and_ac_type(line).unwrap(),
            (accession.to_owned(), accession_type.to_owned()),
            "{line}"
        );
    }
}

#[test]
fn sequest_read_out_header_matches_the_class_test() {
    // SequestOutfile_test.cpp, section readOutHeader.
    let file = SequestOutfile::new();
    assert!(matches!(file.read_out_header(data("a")), Err(Error::Io(_))));
    match file.read_out_header(data("SequestOutfile_headerfile.txt")) {
        Err(Error::Parse { message, .. }) => assert_eq!(message, "No Sequest version found!"),
        other => panic!("expected a parse error, got {other:?}"),
    }
    let header = file
        .read_out_header(data("Sequest.mzXML.13.1.d.out"))
        .unwrap();
    // TOLERANCE_ABSOLUTE(0.0001)
    assert!((header.precursor_mz_value - 866.606).abs() <= 0.0001);
    assert_eq!(header.sequest, "TurboSEQUEST");
    assert_eq!(header.sequest_version, "v.27 (rev. 12)");
    assert_eq!(header.database_type, "amino acids");
    assert_eq!(header.datetime.get(), "2007-01-17 17:29:00");
    assert_eq!(header.charge, 2);
    assert_eq!(header.number_column, Some(0));
    assert_eq!(header.rank_sp_column, Some(1));
    assert_eq!(header.id_column, Some(2));
    assert_eq!(header.mh_column, Some(3));
    assert_eq!(header.delta_cn_column, Some(4));
    assert_eq!(header.xcorr_column, Some(5));
    assert_eq!(header.sp_column, Some(6));
    assert_eq!(header.sf_column, Some(7));
    assert_eq!(header.ions_column, Some(9));
    assert_eq!(header.reference_column, Some(10));
    assert_eq!(header.peptide_column, Some(11));
    assert_eq!(header.score_column, Some(7));
    assert_eq!(header.number_of_columns, 12);
    assert_eq!(header.precursor_mass_type, 0);
    assert_eq!(header.ion_mass_type, 0);
    assert_eq!(header.displayed_peptides, 2);
}

#[test]
fn sequest_protein_hits_accumulate_across_files() {
    // Derived from SequestOutfile.cpp: accessions already among the protein
    // identification's hits are not added again, so loading the same file
    // twice keeps the protein list.
    let file = SequestOutfile::new();
    let mut protein = ProteinIdentification::default();
    file.load(data("SequestOutfile.out"), &mut protein, 1.0, &[], false)
        .unwrap();
    let first: Vec<String> = protein
        .hits
        .iter()
        .map(|hit| hit.accession.clone())
        .collect();
    assert_eq!(first.len(), 14);
    assert_eq!(first[0], "P35574");
    file.load(data("SequestOutfile.out"), &mut protein, 1.0, &[], false)
        .unwrap();
    let second: Vec<String> = protein
        .hits
        .iter()
        .map(|hit| hit.accession.clone())
        .collect();
    assert_eq!(first, second);
    assert_eq!(protein.date_time.as_deref(), Some("2004-03-16 14:32:00"));
}

#[test]
fn sequest_refuses_a_protein_count_past_the_end() {
    // Native refusal: the source reads empty lines past the end of the file.
    let text = std::fs::read_to_string(data("SequestOutfile.out")).unwrap();
    let cut: String = text
        .lines()
        .take(16)
        .map(|line| format!("{line}\n"))
        .collect();
    let path = std::env::temp_dir().join(format!("sequest_cut_{}.out", std::process::id()));
    std::fs::write(&path, cut).unwrap();
    let result = SequestOutfile::new().load(
        &path,
        &mut ProteinIdentification::default(),
        1.0,
        &[],
        false,
    );
    std::fs::remove_file(&path).unwrap();
    match result {
        Err(Error::Parse { line, message }) => {
            assert_eq!(line, 16);
            assert!(message.contains("announces 2 more proteins"), "{message}");
        }
        other => panic!("expected a parse error, got {other:?}"),
    }
}
