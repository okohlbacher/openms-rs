// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/protxml.rs` (`FORMAT/ProtXMLFile.h`).
//!
//! Evidence, strongest first (see `docs/PROTXML_SUPPORT.md` and
//! `tests/data/protxml_provenance.json`):
//!
//! - tier 1, executed differential: the Release C++ `IDFileConverter -test`
//!   of openms4-release-bc9cc12-c19e494-174b576 converted every input below to
//!   idXML on ibminode06, twice, reproduced (`../oracle/F3a-protxml`). Each
//!   output is read back with this crate's idXML reader and compared field by
//!   field with this port's load, and the tool's stderr with this port's log;
//! - tier 1, retained upstream output: `TOPP_IDFileConverter_3`, whose
//!   retained idXML predates `ProteinHit`'s double score (see
//!   `topp_idfileconverter_3_retained_output`);
//! - tier 3: every `START_SECTION` of `ProtXMLFile_test.cpp`, with its literals.
//!
//! The reader needs the `idxml` feature.
#![cfg(feature = "idxml")]

use openms::Error;
use openms::concept::unique_id::UniqueIdGenerator;
use openms::format::idxml;
use openms::format::protxml::{ProtXmlDocument, ProtXmlFile, ReadLimits, SCORE_TYPE};
use openms::identification::{PeptideIdentification, ProteinIdentification};
use std::collections::BTreeSet;

const DATA: &str = "tests/data/protxml";

fn path(name: &str) -> String {
    format!("{DATA}/{name}")
}

fn load(name: &str) -> openms::Result<ProtXmlDocument> {
    ProtXmlFile::new().load(path(name), &mut UniqueIdGenerator::from_seed(1))
}

fn accessions(peptides: &PeptideIdentification, index: usize) -> Vec<&str> {
    peptides.hits[index]
        .evidences
        .iter()
        .map(|evidence| evidence.protein_accession.as_str())
        .collect()
}

// ---------------------------------------------------------------------------
// START_SECTION(void load(const std::string&, ProteinIdentification&, PeptideIdentification&))
// ---------------------------------------------------------------------------

#[test]
fn load_follows_the_class_test() {
    let file = ProtXmlFile::new();
    let mut ids = UniqueIdGenerator::from_seed(7);
    // The class test loads the same file twice to check that members reset.
    for _ in 0..2 {
        let document = file
            .load(path("ProtXMLFile_input_1.protXML"), &mut ids)
            .unwrap();
        let proteins = &document.protein_identification;
        let peptides = &document.peptide_identification;

        // groups
        assert_eq!(proteins.protein_groups.len(), 7);
        assert_eq!(proteins.protein_groups[0].probability, 0.9990);
        assert_eq!(proteins.protein_groups[0].accessions.len(), 1);
        assert_eq!(proteins.protein_groups[3].accessions.len(), 2);
        assert_eq!(
            proteins.protein_groups[3].accessions[0],
            "P01876|IGHA1_HUMAN"
        );
        assert_eq!(
            proteins.protein_groups[3].accessions[1],
            "P01877|IGHA2_HUMAN"
        );
        assert_eq!(proteins.protein_groups[6].probability, 0.2026);
        assert_eq!(proteins.protein_groups[6].accessions.len(), 1);

        assert_eq!(proteins.indistinguishable_groups.len(), 7);
        assert_eq!(proteins.indistinguishable_groups[0].accessions.len(), 1);
        assert_eq!(proteins.indistinguishable_groups[3].accessions.len(), 2);
        assert_eq!(
            proteins.indistinguishable_groups[3].accessions[0],
            "P01876|IGHA1_HUMAN"
        );
        assert_eq!(
            proteins.indistinguishable_groups[3].accessions[1],
            "P01877|IGHA2_HUMAN"
        );
        assert_eq!(proteins.indistinguishable_groups[6].accessions.len(), 1);

        // proteins
        assert_eq!(proteins.hits.len(), 9);
        assert_eq!(proteins.hits[0].accession, "P02787|TRFE_HUMAN");
        assert_eq!(proteins.hits[0].coverage, Some(8.6));
        assert_eq!(proteins.hits[0].score, 0.9990);
        // Indistinguishable, therefore no coverage (the source's -1), but the
        // score is transferred from the leader.
        assert_eq!(proteins.hits[6].accession, "P00739|HPTR_HUMAN");
        assert_eq!(proteins.hits[6].coverage, None);
        assert_eq!(proteins.hits[6].score, 0.2663);
        assert_eq!(proteins.hits[8].accession, "P04217|A1BG_HUMAN");
        assert_eq!(proteins.hits[8].coverage, Some(2.0));
        assert_eq!(proteins.hits[8].score, 0.2026);

        // peptides
        assert_eq!(peptides.hits.len(), 16);
        assert_eq!(peptides.hits[0].sequence.to_string(), "MYLGYEYVTAIR");
        assert_eq!(peptides.hits[0].charge, 2);
        assert_eq!(peptides.hits[0].score, 0.8633);
        let protein_accessions: BTreeSet<&str> = accessions(peptides, 0).into_iter().collect();
        assert_eq!(protein_accessions.len(), 1);
        assert_eq!(
            protein_accessions.into_iter().next(),
            Some("P02787|TRFE_HUMAN")
        );
        // TEST_EQUAL(getMetaValue("is_unique"), true): DataValue(1) == true.
        assert_eq!(peptides.hits[0].metadata["is_unique"].as_i64().unwrap(), 1);
        assert_eq!(
            peptides.hits[0].metadata["is_contributing"]
                .as_i64()
                .unwrap(),
            1
        );
    }
}

// ---------------------------------------------------------------------------
// START_SECTION(void store(...))
// ---------------------------------------------------------------------------

#[test]
fn store_is_not_implemented() {
    let result = ProtXmlFile::new().store(
        "notimplemented.protXML",
        &ProteinIdentification::default(),
        &PeptideIdentification::default(),
        "",
    );
    assert!(matches!(result, Err(Error::Unsupported(_))));
    assert!(!std::path::Path::new("notimplemented.protXML").exists());
}

// ---------------------------------------------------------------------------
// START_SECTION([EXTRA] probability zero protein filtering)
// ---------------------------------------------------------------------------

#[test]
fn probability_zero_proteins_are_dropped() {
    let document = load("ProtXMLFile_input_3.protXML").unwrap();
    let proteins = &document.protein_identification;
    let peptides = &document.peptide_identification;

    // Five groups in the input; the all-subsumable fifth is not inserted.
    assert_eq!(proteins.protein_groups.len(), 4);
    assert_eq!(proteins.protein_groups[0].accessions, ["LEADER_PROT"]);
    assert!((proteins.protein_groups[0].probability - 0.9998).abs() <= 1e-5 * 0.9998);
    assert_eq!(proteins.protein_groups[1].accessions, ["PROT_A", "PROT_B"]);
    assert_eq!(proteins.protein_groups[2].accessions, ["SINGLE_PROT"]);
    assert_eq!(proteins.protein_groups[3].accessions, ["VALID_AFTER"]);
    assert!((proteins.protein_groups[3].probability - 0.95).abs() <= 1e-5 * 0.95);

    assert_eq!(proteins.indistinguishable_groups.len(), 4);
    assert_eq!(
        proteins.indistinguishable_groups[0].accessions,
        ["LEADER_PROT"]
    );
    assert_eq!(
        proteins.indistinguishable_groups[1].accessions,
        ["PROT_A", "PROT_B"]
    );
    assert_eq!(
        proteins.indistinguishable_groups[2].accessions,
        ["SINGLE_PROT"]
    );
    assert_eq!(
        proteins.indistinguishable_groups[3].accessions,
        ["VALID_AFTER"]
    );

    let hits: Vec<(&str, f64)> = proteins
        .hits
        .iter()
        .map(|hit| (hit.accession.as_str(), hit.score))
        .collect();
    let expected = [
        ("LEADER_PROT", 0.9998),
        ("PROT_A", 0.9000),
        ("PROT_B", 0.9000), // inherited from the leader
        ("SINGLE_PROT", 0.8000),
        ("VALID_AFTER", 0.9500),
    ];
    assert_eq!(hits.len(), expected.len());
    for ((accession, score), (want, want_score)) in hits.iter().zip(expected) {
        assert_eq!(*accession, want);
        assert!((score - want_score).abs() <= 1e-5 * want_score);
    }

    let sequences: Vec<String> = peptides
        .hits
        .iter()
        .map(|hit| hit.sequence.to_string())
        .collect();
    assert_eq!(
        sequences,
        [
            "AAAPEPTIDE",
            "BBPEPTIDE",
            "SHAREDPEPTIDE",
            "UNIQUEPEPTIDE",
            "VALIDPEP"
        ]
    );
    assert_eq!(accessions(peptides, 0), ["LEADER_PROT"]);
    assert_eq!(accessions(peptides, 2), ["PROT_A", "PROT_B"]);
    assert_eq!(accessions(peptides, 4), ["VALID_AFTER"]);
}

// ---------------------------------------------------------------------------
// Tier 1: the Release IDFileConverter on every input
// ---------------------------------------------------------------------------

/// The Release build's stderr lines that this port's log reproduces, with the
/// input's name normalised to the oracle's `input.protXML`.
fn our_log(document: &ProtXmlDocument, name: &str) -> Vec<String> {
    document
        .log
        .iter()
        .map(|line| line.replace(&path(name), "input.protXML"))
        .collect()
}

fn expected_log(case: &str) -> Vec<String> {
    std::fs::read_to_string(path(&format!("expected/{case}.stderr.txt")))
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect()
}

/// Modification names the Release build's `ModificationsDB` loads from
/// PSI-MOD (pinned `share/OpenMS/CHEMISTRY/PSI-MOD.obo`, accession given), and
/// which this crate's global registry does not hold: historical PSI-MOD is not
/// bundled (`docs/MODIFICATION_REGISTRY_SUPPORT.md`, "Historical PSI-MOD is not
/// bundled"). They appear only as extra candidates in the source's
/// "not uniquely defined" warnings of these inputs.
const PSI_MOD_ONLY: [(&str, &str); 5] = [
    ("oxidation to L-methionine sulfoxide (M)", "MOD:00719"),
    ("oxidation to L-methionine (R)-sulfoxide (M)", "MOD:00720"),
    ("oxidation to L-methionine (S)-sulfoxide (M)", "MOD:00721"),
    ("N-acetyl-L-alanine (N-term A)", "MOD:00050"),
    ("L-lysine amide (C-term K)", "MOD:00101"),
];

/// The Release log as a registry without PSI-MOD produces it: every
/// [`PSI_MOD_ONLY`] name leaves the candidate list of an ambiguity warning, and
/// a warning left with a single candidate is not emitted (the lookup is then
/// unique). Every other line, and every other candidate, is compared exactly.
fn expected_log_without_psi_mod(case: &str) -> Vec<String> {
    const MARK: &str = "' to represent any of '";
    let mut out = Vec::new();
    for line in expected_log(case) {
        let Some(at) = line.find(MARK) else {
            out.push(line);
            continue;
        };
        let head = &line[..at + MARK.len()];
        let list = line[at + MARK.len()..]
            .strip_suffix("'!")
            .expect("ambiguity warning ends in '!");
        let kept: Vec<&str> = list
            .split(", ")
            .filter(|name| PSI_MOD_ONLY.iter().all(|(psi, _)| psi != name))
            .collect();
        let chosen = head
            .rsplit("Using '")
            .next()
            .and_then(|rest| rest.strip_suffix(MARK))
            .expect("ambiguity warning names its choice");
        assert!(
            PSI_MOD_ONLY.iter().all(|(psi, _)| *psi != chosen),
            "{case}: the source chose a PSI-MOD record, which this comparison cannot model"
        );
        if kept.len() > 1 {
            out.push(format!("{head}{}'!", kept.join(", ")));
        }
    }
    out
}

/// Compare this port's load with the Release C++ output. `IdXMLFile::store`
/// sorts every identification's hits by score (`std::stable_sort`), so the
/// port's hits are sorted the same way first; everything else keeps its
/// order. The run date is compared as the idXML writer spells it: an
/// invalid `DateTime` is `0000-00-00T00:00:00`.
fn assert_matches_release(ours: &ProtXmlDocument, reference: &idxml::IdXmlDocument, case: &str) {
    assert_eq!(reference.protein_identifications.len(), 1, "{case}");
    assert_eq!(reference.peptide_identifications.len(), 1, "{case}");
    let (proteins, theirs) = (
        &ours.protein_identification,
        &reference.protein_identifications[0],
    );
    assert_eq!(proteins.search_engine, theirs.search_engine, "{case}");
    assert_eq!(
        proteins.search_engine_version, theirs.search_engine_version,
        "{case}"
    );
    assert_eq!(
        proteins
            .date_time
            .as_deref()
            .unwrap_or("0000-00-00T00:00:00"),
        theirs.date_time.as_deref().unwrap_or("0000-00-00T00:00:00"),
        "{case}"
    );
    assert_eq!(proteins.score_type, theirs.score_type, "{case}");
    assert_eq!(
        proteins.higher_score_better, theirs.higher_score_better,
        "{case}"
    );
    assert_eq!(
        proteins.search_parameters.database, theirs.search_parameters.database,
        "{case}"
    );
    assert_eq!(
        proteins.search_parameters.digestion_enzyme.to_lowercase(),
        theirs.search_parameters.digestion_enzyme.to_lowercase(),
        "{case}: idXML writes the enzyme name lower-cased"
    );
    assert_eq!(proteins.hits.len(), theirs.hits.len(), "{case}");
    for (index, (ours, theirs)) in proteins.hits.iter().zip(&theirs.hits).enumerate() {
        assert_eq!(ours.accession, theirs.accession, "{case} protein {index}");
        assert_eq!(ours.score, theirs.score, "{case} protein {index}");
        assert_eq!(ours.coverage, theirs.coverage, "{case} protein {index}");
        assert!(ours.metadata.is_empty(), "{case} protein {index}");
    }
    assert_eq!(proteins.protein_groups, theirs.protein_groups, "{case}");
    assert_eq!(
        proteins.indistinguishable_groups, theirs.indistinguishable_groups,
        "{case}"
    );

    let mut peptides = ours.peptide_identification.clone();
    peptides.sort().unwrap();
    let theirs = &reference.peptide_identifications[0];
    assert_eq!(peptides.score_type, theirs.score_type, "{case}");
    assert_eq!(
        peptides.higher_score_better, theirs.higher_score_better,
        "{case}"
    );
    assert_eq!(peptides.hits.len(), theirs.hits.len(), "{case}");
    for (index, (hit, expected)) in peptides.hits.iter().zip(&theirs.hits).enumerate() {
        assert_eq!(
            hit.sequence.to_string(),
            expected.sequence.to_string(),
            "{case} peptide {index}"
        );
        assert_eq!(hit.score, expected.score, "{case} peptide {index}");
        assert_eq!(hit.charge, expected.charge, "{case} peptide {index}");
        let ours: Vec<&str> = hit
            .evidences
            .iter()
            .map(|e| e.protein_accession.as_str())
            .collect();
        let want: Vec<&str> = expected
            .evidences
            .iter()
            .map(|e| e.protein_accession.as_str())
            .collect();
        assert_eq!(ours, want, "{case} peptide {index}");
        assert_eq!(hit.metadata, expected.metadata, "{case} peptide {index}");
    }
}

fn converts_like_the_release_build(input: &str, case: &str) {
    let ours = load(input).unwrap();
    let reference = idxml::load(path(&format!("expected/{case}.idXML"))).unwrap();
    assert_matches_release(&ours, &reference, case);
    assert_eq!(
        our_log(&ours, input),
        expected_log_without_psi_mod(case),
        "{case}"
    );
    assert_eq!(ours.suppressed_log_lines, 0);
}

#[test]
fn release_idfileconverter_3() {
    converts_like_the_release_build("IDFileConverter_3_input.protXML", "IDFileConverter_3_input");
}

#[test]
fn release_class_test_inputs() {
    for n in 1..=3 {
        let case = format!("ProtXMLFile_input_{n}");
        converts_like_the_release_build(&format!("{case}.protXML"), &case);
    }
}

#[test]
fn release_modification_matching() {
    converts_like_the_release_build("f3a_mods.protXML", "f3a_mods");
    // The source's look-up uses the first registry match; for a
    // carbamidomethylated cysteine that is `Carbamidomethyl (N-term)`, so the
    // modification lands on the N-terminus, and an acetylated first residue
    // stays unmodified because `Acetyl (Protein N-term)` has three words.
    let document = load("f3a_mods.protXML").unwrap();
    let sequences: Vec<String> = document
        .peptide_identification
        .hits
        .iter()
        .map(|hit| hit.sequence.to_string())
        .collect();
    assert_eq!(
        sequences,
        [
            "PEPTM(Oxidation)IDEK",
            ".(Carbamidomethyl)CPEPTIDEK",
            "APEPTIDEK",
            "PEPTIDEK.(Amidated)",
            "QPEPTIDEK",
            "PEPTIDEKR",
            ".(Carbamidomethyl)M(Oxidation)CPEPTIDEK",
        ]
    );
}

#[test]
fn release_missing_optional_attributes_and_header() {
    for case in ["f3a_missing", "f3a_empty_group", "f3a_enzyme_case"] {
        converts_like_the_release_build(&format!("{case}.protXML"), case);
    }
    let document = load("f3a_missing.protXML").unwrap();
    assert_eq!(document.protein_identification.hits[0].coverage, None);
    assert_eq!(document.peptide_identification.hits[0].charge, 0);
    let bare = load("f3a_no_header.protXML").unwrap();
    assert_eq!(bare.protein_identification.identifier, "");
    assert_eq!(bare.protein_identification.date_time, None);
    assert_eq!(bare.protein_identification.score_type, "");
    assert_eq!(
        bare.protein_identification
            .search_parameters
            .digestion_enzyme,
        "unknown_enzyme"
    );
    let metadata = &bare.peptide_identification.hits[0].metadata;
    assert_eq!(metadata["is_unique"].as_i64().unwrap(), 0);
    assert_eq!(metadata["is_contributing"].as_i64().unwrap(), 0);
    // The Release build writes an unset run date as `0000-00-00T00:00:00`,
    // which neither its own idXML reader (`DateTime::fromString`) nor this
    // crate's accepts, so this output is checked on its text.
    let release = std::fs::read_to_string(path("expected/f3a_no_header.idXML")).unwrap();
    for line in [
        "enzyme=\"unknown_enzyme\"",
        "<IdentificationRun date=\"0000-00-00T00:00:00\" search_engine=\"\" search_engine_version=\"\"",
        "<ProteinIdentification score_type=\"\" higher_score_better=\"true\"",
        "<ProteinHit id=\"PH_0\" accession=\"BARE\" score=\"0.5\" coverage=\"10.0\"",
        "name=\"protein_group_0\" value=\"0.5,PH_0\"",
        "name=\"indistinguishable_proteins_0\" value=\"0.0,PH_0\"",
        "<PeptideHit score=\"0.9\" sequence=\"BAREPEPTIDE\" charge=\"2\" protein_refs=\"PH_0\"",
        "name=\"is_unique\" value=\"0\"",
        "name=\"is_contributing\" value=\"0\"",
    ] {
        assert!(release.contains(line), "{line}");
    }
    let proteins = &bare.protein_identification;
    assert_eq!(proteins.search_engine, "");
    assert_eq!(proteins.search_engine_version, "");
    assert_eq!(proteins.search_parameters.database, "");
    assert_eq!(proteins.hits.len(), 1);
    assert_eq!(proteins.hits[0].accession, "BARE");
    assert_eq!(proteins.hits[0].score, 0.5);
    assert_eq!(proteins.hits[0].coverage, Some(10.0));
    assert_eq!(proteins.protein_groups.len(), 1);
    assert_eq!(proteins.protein_groups[0].probability, 0.5);
    assert_eq!(proteins.protein_groups[0].accessions, ["BARE"]);
    assert_eq!(proteins.indistinguishable_groups.len(), 1);
    assert_eq!(proteins.indistinguishable_groups[0].probability, 0.0);
    let hit = &bare.peptide_identification.hits[0];
    assert_eq!(hit.sequence.to_string(), "BAREPEPTIDE");
    assert_eq!((hit.score, hit.charge), (0.9, 2));
    assert_eq!(bare.peptide_identification.score_type, "");
    assert!(bare.log.is_empty());
    assert!(expected_log("f3a_no_header").is_empty());
}

/// Which error a refused load must be.
type ErrorKind = fn(&Error) -> bool;

/// The Release build fails these loads (exit 3 or 8); so does the port.
#[test]
fn release_refusals() {
    let cases: [(&str, ErrorKind); 6] = [
        // ElementNotFound: the enzyme 'NoSuchEnzyme' could not be found.
        ("f3a_bad_enzyme", |e| matches!(e, Error::InvalidValue(_))),
        // ParseError: Invalid date time string.
        ("f3a_bad_time", |e| matches!(e, Error::InvalidValue(_))),
        // Xerces NumberFormatException: invalid character encountered.
        ("f3a_charge_text", |e| matches!(e, Error::Parse { .. })),
        // fatalError: Required attribute 'probability' not present!
        ("f3a_missing_required", |e| matches!(e, Error::Parse { .. })),
        // IndexOverflow: the given index was too large: 8 (size = 8).
        ("f3a_pos_high", |e| matches!(e, Error::InvalidRange(_))),
        // IndexOverflow: position 0 wraps to -1 in Size.
        ("f3a_pos_zero", |e| matches!(e, Error::InvalidRange(_))),
    ];
    for (case, kind) in cases {
        let error = load(&format!("{case}.protXML")).unwrap_err();
        assert!(kind(&error), "{case}: {error:?}");
        assert!(
            !std::path::Path::new(&path(&format!("expected/{case}.idXML"))).exists(),
            "{case} has an output"
        );
        assert_eq!(expected_log(case).len().min(1), 1, "{case} left no message");
    }
    let missing = load("f3a_missing_required.protXML").unwrap_err();
    assert!(
        missing
            .to_string()
            .contains("Required attribute 'probability' not present!"),
        "{missing}"
    );
}

/// The Release build accepts these; the port refuses them because they
/// cannot be represented: a charge the source's `(int)` cast wraps from
/// 4294967298 to 2, and a NaN peptide probability.
#[test]
fn unrepresentable_values_are_refused() {
    let release = idxml::load(path("expected/f3a_charge_wrap.idXML")).unwrap();
    assert_eq!(release.peptide_identifications[0].hits[0].charge, 2);
    assert!(matches!(
        load("f3a_charge_wrap.protXML"),
        Err(Error::InvalidValue(_))
    ));
    let nan = std::fs::read_to_string(path("expected/f3a_nsp_nan.idXML")).unwrap();
    assert!(nan.contains("<PeptideHit score=\"NaN\""));
    assert!(matches!(
        load("f3a_nsp_nan.protXML"),
        Err(Error::InvalidValue(_))
    ));
}

// ---------------------------------------------------------------------------
// Tier 1, retained upstream output
// ---------------------------------------------------------------------------

/// `TOPP_IDFileConverter_3` (test-data `topp/CMakeLists.txt:1202-1204`). The
/// retained idXML was written when `ProteinHit` stored its score as a
/// `float`: it spells 0.9990 as `0.999000012874603`, the float nearest 0.999
/// widened to double. At the pin the score is a `double` and the Release
/// build writes `0.999` (see `release_idfileconverter_3`); the upstream
/// FuzzyDiff accepts both. Here the protein scores are compared after
/// rounding both to `f32`, which is exactly the difference between the two
/// builds, and every other field exactly.
#[test]
fn topp_idfileconverter_3_retained_output() {
    let ours = load("IDFileConverter_3_input.protXML").unwrap();
    let mut retained = idxml::load(path("IDFileConverter_3_output.idXML")).unwrap();
    let release = idxml::load(path("expected/IDFileConverter_3_input.idXML")).unwrap();
    let mut widened = 0;
    for (hit, pin) in retained.protein_identifications[0]
        .hits
        .iter_mut()
        .zip(&release.protein_identifications[0].hits)
    {
        assert_eq!(hit.score as f32, pin.score as f32, "{}", hit.accession);
        if hit.score != pin.score {
            widened += 1;
            hit.score = pin.score;
        }
    }
    assert_eq!(widened, 9);
    assert_matches_release(&ours, &retained, "IDFileConverter_3 retained");
}

// ---------------------------------------------------------------------------
// Tier 4: identifiers, streams and ceilings
// ---------------------------------------------------------------------------

#[test]
fn identifiers_come_from_the_callers_generator() {
    let mut ids = UniqueIdGenerator::from_seed(42);
    let document = ProtXmlFile::new()
        .load(path("ProtXMLFile_input_3.protXML"), &mut ids)
        .unwrap();
    let expected = UniqueIdGenerator::from_seed(42).get_unique_id().to_string();
    assert_eq!(document.protein_identification.identifier, expected);
    assert_eq!(document.peptide_identification.identifier, expected);
    assert_eq!(document.protein_identification.score_type, SCORE_TYPE);
    assert_eq!(document.peptide_identification.score_type, SCORE_TYPE);
    assert!(document.peptide_identification.higher_score_better);
}

#[test]
fn read_matches_load_and_names_no_file() {
    let text = std::fs::read(path("f3a_mods.protXML")).unwrap();
    let read = ProtXmlFile::new()
        .read(text.as_slice(), &mut UniqueIdGenerator::from_seed(1))
        .unwrap();
    let loaded = load("f3a_mods.protXML").unwrap();
    assert_eq!(read.protein_identification, loaded.protein_identification);
    assert_eq!(read.peptide_identification, loaded.peptide_identification);
    assert!(read.log[0].starts_with("Non-fatal error while loading '': "));
}

fn wrap(body: &str) -> String {
    format!(
        "<protein_summary><protein_summary_header reference_database=\"d\" sample_enzyme=\"trypsin\">\
         <program_details analysis=\"a\" time=\"2024-01-01T00:00:00\" version=\"v\"/>\
         </protein_summary_header>{body}</protein_summary>"
    )
}

fn read_text(text: &str) -> openms::Result<ProtXmlDocument> {
    ProtXmlFile::new().read(text.as_bytes(), &mut UniqueIdGenerator::from_seed(1))
}

const PEPTIDE: &str = "<peptide peptide_sequence=\"PEPTIDE\" charge=\"2\" nsp_adjusted_probability=\"0.5\" is_nondegenerate_evidence=\"Y\" is_contributing_evidence=\"Y\">";

#[test]
fn shapes_the_source_leaves_undefined_are_refused() {
    // <indistinguishable_protein> before any protein: getHits().back() on an
    // empty vector.
    let orphan = wrap(
        "<protein_group probability=\"1\"><indistinguishable_protein protein_name=\"X\"/></protein_group>",
    );
    assert!(matches!(read_text(&orphan), Err(Error::InvalidValue(_))));
    // <peptide> before any protein: getIndistinguishableProteins().back().
    let peptide = wrap(&format!(
        "<protein_group probability=\"1\">{PEPTIDE}</peptide></protein_group>"
    ));
    assert!(matches!(read_text(&peptide), Err(Error::InvalidValue(_))));
    // <mod_aminoacid_mass> outside a peptide: a null pep_hit_.
    let modification = wrap(
        "<protein_group probability=\"1\"><protein protein_name=\"P\" probability=\"1\" percent_coverage=\"1\">\
         <mod_aminoacid_mass position=\"1\" mass=\"1\"/></protein></protein_group>",
    );
    assert!(matches!(
        read_text(&modification),
        Err(Error::InvalidValue(_))
    ));
    // A non-finite group probability or coverage cannot be stored.
    let group = wrap("<protein_group probability=\"inf\"></protein_group>");
    assert!(matches!(read_text(&group), Err(Error::InvalidValue(_))));
}

#[test]
fn a_prefixed_element_is_not_recognised() {
    // The source matches the qualified name, so `p:protein` is ignored.
    let text = wrap(
        "<protein_group probability=\"1\"><p:protein xmlns:p=\"u\" protein_name=\"P\" probability=\"1\"/></protein_group>",
    );
    let document = read_text(&text).unwrap();
    assert!(document.protein_identification.hits.is_empty());
    assert!(document.protein_identification.protein_groups.is_empty());
}

#[test]
fn ceilings_are_enforced() {
    let text = std::fs::read(path("ProtXMLFile_input_1.protXML")).unwrap();
    let file = ProtXmlFile::new();
    let mut ids = UniqueIdGenerator::from_seed(1);
    let limited = |limits: ReadLimits, ids: &mut UniqueIdGenerator| {
        file.read_with_limits(text.as_slice(), ids, &limits)
    };
    let base = ReadLimits::default();
    // The shared decoder reports its byte ceiling as a parse error.
    assert!(matches!(
        limited(
            ReadLimits {
                max_bytes: 1000,
                ..base
            },
            &mut ids
        ),
        Err(Error::Parse { .. })
    ));
    for limits in [
        ReadLimits {
            max_events: 10,
            ..base
        },
        ReadLimits {
            max_depth: 2,
            ..base
        },
        ReadLimits {
            max_proteins: 8,
            ..base
        },
        ReadLimits {
            max_peptides: 15,
            ..base
        },
    ] {
        assert!(
            matches!(limited(limits, &mut ids), Err(Error::InvalidValue(_))),
            "{limits:?}"
        );
    }
    let exact = ReadLimits {
        max_proteins: 9,
        max_peptides: 16,
        max_log_lines: 1,
        ..base
    };
    let document = limited(exact, &mut ids).unwrap();
    assert_eq!(document.log.len(), 1);
    assert_eq!(document.suppressed_log_lines, 1);
}
