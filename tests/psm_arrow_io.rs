// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `PSMArrowIO` (`FORMAT/PSMArrowIO.h`), the `.idparquet` directory bundle.
//!
//! Two kinds of evidence:
//!
//! * **The upstream class test** `PSMArrowIO_test.cpp` (core `bc9cc12`), every
//!   section, with its fixtures and literals (tier 3).
//! * **C++-written bundles.** `IDMerger_idparquet_input{1,2}.idparquet`,
//!   `IDRipper_4_input.idparquet` and `PercolatorAdapter_idparquet_in.idparquet`
//!   are upstream TOPP fixtures (test data `0cb15f2`, copied unchanged into
//!   `tests/data/arrow_parquet/fixtures/topp/` by lane F1a). The first two were
//!   written by the C++ `IDFileConverter -in IDMerger_1_inputX.idXML -out
//!   IDMerger_idparquet_inputX.idparquet` (`topp/CMakeLists.txt:178` at the
//!   test-data pin), so their import must give back those two idXML files'
//!   identifications, which are transcribed here; the retained
//!   `IDMerger_idparquet_output.idXML` (the C++ import, merged by IDMerger,
//!   which sorts each identification's hits, `IDMerger.cpp:411`) shows the same
//!   hits, scores, charges and protein references. All four are also checked
//!   against the format's own invariants (tier 4).
//! * **Executed differential** (tier 1, below): what the Linux x86_64 Release
//!   build's `PSMArrowIO` imports from those four fixtures and from its own
//!   bundles of the class-test data, field for field.
#![cfg(feature = "arrow")]

use arrow_array::cast::AsArray;
use arrow_array::types::Int32Type;
use arrow_array::{Array, RecordBatch};
use openms::Error;
use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationRecord, ModificationsDB, ResidueModification,
};
use openms::data_structures::datetime::DateTime;
use openms::format::arrow_io_helpers::ParquetWriteConfig;
use openms::format::parquet_file::ParquetFile;
use openms::format::psm_arrow_io::{PROTEIN_GROUPS, PROTEINS, PSMArrowIO, PSMS, SEARCH_PARAMS};
use openms::identification::{
    FlankingResidue, PeptideEvidence, PeptideHit, PeptideIdentification, ProteinIdentification,
};
use openms::metadata::{MetaValue, MetaValueData};
use openms::system::file::TempDir;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

/// The source's `Constants::UserParam::MODIFICATION_DEFINITIONS`.
const MODIFICATION_DEFINITIONS: &str = "modification_definitions";

fn fixture(bundle: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/arrow_parquet/fixtures/topp")
        .join(bundle)
}

fn formula(text: &str) -> EmpiricalFormula {
    text.parse().unwrap()
}

fn config() -> ParquetWriteConfig {
    ParquetWriteConfig::default()
}

/// The class test's `buildMinimalIds`.
fn build_minimal_ids() -> (Vec<ProteinIdentification>, Vec<PeptideIdentification>) {
    let mut prot = ProteinIdentification {
        identifier: "run_1".into(),
        score_type: "score".into(),
        higher_score_better: true,
        primary_ms_run_paths: vec!["/tmp/foo.mzML".into()],
        ..ProteinIdentification::default()
    };
    prot.search_parameters.digestion_enzyme = "Trypsin".into();
    prot.search_parameters
        .metadata
        .insert("extra_features".into(), "COMET:deltaCn,MS:1002049".into());

    let mut pid = PeptideIdentification {
        identifier: "run_1".into(),
        score_type: "score".into(),
        higher_score_better: true,
        rt: Some(123.4),
        mz: Some(567.89),
        ..PeptideIdentification::default()
    };
    pid.set_spectrum_reference("scan=42");
    let mut hit = PeptideHit {
        sequence: AASequence::parse("PEPTIDE").unwrap(),
        charge: 2,
        score: 0.95,
        ..PeptideHit::default()
    };
    hit.metadata.insert("target_decoy".into(), "target".into());
    hit.metadata
        .insert("COMET:deltaCn".into(), MetaValue::try_from(0.5).unwrap());
    hit.evidences.push(PeptideEvidence {
        protein_accession: "sp|P12345|EXAMPLE".into(),
        aa_before: FlankingResidue::Residue('K'),
        aa_after: FlankingResidue::Residue('R'),
        start: Some(42),
        end: Some(48),
    });
    pid.hits.push(hit);
    (vec![prot], vec![pid])
}

/// The class test's `defineMod4b` / `freshRecord4b` record: a tool-defined
/// modification of `origin`, anywhere, with `formula` as its difference.
fn definition(id: &str, origin: char, diff: &str) -> ResidueModification {
    let diff = formula(diff);
    ResidueModification::from_record(ModificationRecord {
        name: id.into(),
        origin: Some(origin),
        diff_mono_mass: diff.mono_mass(),
        diff_formula: diff,
        ..ModificationRecord::default()
    })
    .unwrap()
}

fn read(path: &Path) -> RecordBatch {
    ParquetFile::read_table(path).unwrap()
}

/// The class test's `countSpMetaKey4b`: how many `sp_metavalues` entries of
/// row 0 carry `key`.
fn count_sp_meta_key(search_params: &Path, key: &str) -> usize {
    let table = read(search_params);
    let list = table
        .column_by_name("sp_metavalues")
        .unwrap()
        .as_list::<i32>();
    let entries = list.value(0);
    let names = entries
        .as_struct()
        .column_by_name("name")
        .unwrap()
        .as_string::<i32>();
    (0..names.len()).filter(|&i| names.value(i) == key).count()
}

/// The class test's `columnContains4b`.
fn column_contains(path: &Path, column: &str, needle: &str) -> bool {
    let table = read(path);
    let values = table.column_by_name(column).unwrap().as_string::<i32>();
    (0..values.len()).any(|i| !values.is_null(i) && values.value(i).contains(needle))
}

fn new_dir(tmp: &TempDir, name: &str) -> PathBuf {
    tmp.path().join(format!("{name}.idparquet"))
}

// START_SECTION(([EXTRA] export_then_import_round_trip))
#[test]
fn export_then_import_round_trip() {
    let (prot_ids, pep_ids) = build_minimal_ids();
    let tmp = TempDir::new(false).unwrap();
    let dir = new_dir(&tmp, "round_trip");

    PSMArrowIO::export_to_parquet(&prot_ids, &pep_ids, &dir, true, &config()).unwrap();

    assert!(dir.join("psms.parquet").exists());
    assert!(dir.join("proteins.parquet").exists());
    assert!(dir.join("protein_groups.parquet").exists());
    assert!(dir.join("search_params.parquet").exists());

    let mut prot_ids_in = Vec::new();
    let mut pep_ids_in = Vec::new();
    PSMArrowIO::import_from_parquet(&dir, &mut prot_ids_in, &mut pep_ids_in).unwrap();

    assert_eq!(prot_ids_in.len(), 1);
    // Identifier is synthesized on load per IdXMLFile.cpp:530 parity: the
    // stored "run_1" is replaced, and the identifications are re-stamped in
    // lock-step so both sides agree.
    assert_ne!(prot_ids_in[0].identifier, "");
    assert_ne!(prot_ids_in[0].identifier, "run_1");
    assert_eq!(pep_ids_in[0].identifier, prot_ids_in[0].identifier);
    assert_eq!(prot_ids_in[0].search_parameters.digestion_enzyme, "Trypsin");
    assert_eq!(
        prot_ids_in[0].search_parameters.metadata["extra_features"]
            .as_str()
            .unwrap(),
        "COMET:deltaCn,MS:1002049"
    );

    // Primary MS run path round-trips (spectra_data UserParam).
    assert_eq!(prot_ids_in[0].primary_ms_run_paths, ["/tmp/foo.mzML"]);

    assert_eq!(pep_ids_in.len(), 1);
    assert_eq!(pep_ids_in[0].hits.len(), 1);
    let hit = &pep_ids_in[0].hits[0];
    assert_eq!(hit.sequence.to_string(), "PEPTIDE");
    let delta_cn = hit.metadata["COMET:deltaCn"].as_f64().unwrap();
    assert!((delta_cn - 0.5).abs() <= 1e-5);

    // PeptideEvidence positional fields round-trip (.idparquet schema extension).
    assert_eq!(hit.evidences.len(), 1);
    let ev_in = &hit.evidences[0];
    assert_eq!(ev_in.protein_accession, "sp|P12345|EXAMPLE");
    assert_eq!(ev_in.aa_before, FlankingResidue::Residue('K'));
    assert_eq!(ev_in.aa_after, FlankingResidue::Residue('R'));
    assert_eq!(ev_in.start, Some(42));
    assert_eq!(ev_in.end, Some(48));
}

// START_SECTION(([EXTRA] import_missing_subfile_returns_false))
#[test]
fn import_missing_subfile_is_refused() {
    let (prot_ids, pep_ids) = build_minimal_ids();
    let tmp = TempDir::new(false).unwrap();
    let dir = new_dir(&tmp, "missing");
    PSMArrowIO::export_to_parquet(&prot_ids, &pep_ids, &dir, true, &config()).unwrap();

    // Delete one subfile.
    std::fs::remove_file(dir.join("psms.parquet")).unwrap();

    let mut prot_ids_in = vec![ProteinIdentification::default()];
    let mut pep_ids_in = Vec::new();
    match PSMArrowIO::import_from_parquet(&dir, &mut prot_ids_in, &mut pep_ids_in) {
        // The source's logged message (PSMArrowIO.cpp:115).
        Err(Error::InvalidValue(message)) => assert_eq!(
            message,
            format!(
                "PSMArrowIO: missing required file: {}",
                dir.join("psms.parquet").display()
            )
        ),
        other => panic!("expected the missing-file refusal, got {other:?}"),
    }
    // The outputs are untouched (the source swaps them in only on success).
    assert_eq!(prot_ids_in, vec![ProteinIdentification::default()]);
    assert!(pep_ids_in.is_empty());
}

// START_SECTION(([EXTRA] export_target_is_regular_file_returns_false))
#[test]
fn export_target_is_regular_file_is_refused() {
    let (prot_ids, pep_ids) = build_minimal_ids();
    let tmp = TempDir::new(false).unwrap();
    let path = new_dir(&tmp, "regular");
    // Create a regular file at the target path.
    std::fs::write(&path, "not a directory").unwrap();

    match PSMArrowIO::export_to_parquet(&prot_ids, &pep_ids, &path, true, &config()) {
        // The source's logged message (PSMArrowIO.cpp:37).
        Err(Error::InvalidValue(message)) => assert_eq!(
            message,
            format!(
                "PSMArrowIO: target path exists and is not a directory: {}",
                path.display()
            )
        ),
        other => panic!("expected the not-a-directory refusal, got {other:?}"),
    }
    assert_eq!(std::fs::read(&path).unwrap(), b"not a directory");
}

// START_SECTION(([EXTRA] empty_psms_round_trips))
#[test]
fn empty_psms_round_trips() {
    let (prot_ids, mut pep_ids) = build_minimal_ids();
    pep_ids.clear();
    let tmp = TempDir::new(false).unwrap();
    let dir = new_dir(&tmp, "empty");

    PSMArrowIO::export_to_parquet(&prot_ids, &pep_ids, &dir, true, &config()).unwrap();

    let mut prot_ids_in = Vec::new();
    let mut pep_ids_in = Vec::new();
    PSMArrowIO::import_from_parquet(&dir, &mut prot_ids_in, &mut pep_ids_in).unwrap();
    assert_eq!(pep_ids_in.len(), 0);
    assert_eq!(prot_ids_in.len(), 1);
}

/// The class test's decoy data set: a single decoy PSM (target_decoy=decoy).
fn decoy_ids() -> (Vec<ProteinIdentification>, Vec<PeptideIdentification>) {
    let prot_ids = vec![ProteinIdentification {
        identifier: "run_decoy".into(),
        score_type: "score".into(),
        higher_score_better: true,
        ..ProteinIdentification::default()
    }];
    let mut hit = PeptideHit {
        sequence: AASequence::parse("DECOYSEQ").unwrap(),
        charge: 2,
        score: 0.1,
        ..PeptideHit::default()
    };
    hit.metadata.insert("target_decoy".into(), "decoy".into());
    hit.evidences.push(PeptideEvidence {
        protein_accession: "DECOY_sp|P99999|FAKE".into(),
        ..PeptideEvidence::default()
    });
    let pep_ids = vec![PeptideIdentification {
        identifier: "run_decoy".into(),
        score_type: "score".into(),
        higher_score_better: true,
        hits: vec![hit],
        ..PeptideIdentification::default()
    }];
    (prot_ids, pep_ids)
}

// START_SECTION(([EXTRA] decoy_round_trips))
#[test]
fn decoy_round_trips() {
    // Build a single decoy PSM (target_decoy=decoy) and verify is_decoy round-trips.
    let (prot_ids, pep_ids) = decoy_ids();
    let tmp = TempDir::new(false).unwrap();
    let dir = new_dir(&tmp, "decoy");

    PSMArrowIO::export_to_parquet(&prot_ids, &pep_ids, &dir, true, &config()).unwrap();

    let mut prot_ids_in = Vec::new();
    let mut pep_ids_in = Vec::new();
    PSMArrowIO::import_from_parquet(&dir, &mut prot_ids_in, &mut pep_ids_in).unwrap();

    assert_eq!(pep_ids_in.len(), 1);
    assert_eq!(pep_ids_in[0].hits.len(), 1);
    assert_eq!(
        pep_ids_in[0].hits[0].metadata["target_decoy"]
            .as_str()
            .unwrap(),
        "decoy"
    );
}

/// The class test's multi-rank data set: one identification with three
/// rank-ordered PSMs.
fn multi_rank_ids() -> (Vec<ProteinIdentification>, Vec<PeptideIdentification>) {
    let prot_ids = vec![ProteinIdentification {
        identifier: "run_multirank".into(),
        score_type: "score".into(),
        higher_score_better: true,
        ..ProteinIdentification::default()
    }];
    let mut pid = PeptideIdentification {
        identifier: "run_multirank".into(),
        score_type: "score".into(),
        higher_score_better: true,
        ..PeptideIdentification::default()
    };
    for rank in 1..=3u32 {
        let mut hit = PeptideHit {
            sequence: AASequence::parse(&format!("{}EPTIDE", "P".repeat(rank as usize))).unwrap(),
            charge: 2,
            score: 1.0 / f64::from(rank),
            rank,
            ..PeptideHit::default()
        };
        hit.evidences.push(PeptideEvidence {
            protein_accession: format!("sp|P{rank}|RANK"),
            ..PeptideEvidence::default()
        });
        pid.hits.push(hit);
    }
    (prot_ids, vec![pid])
}

// START_SECTION(([EXTRA] multi_rank_round_trips))
#[test]
fn multi_rank_round_trips() {
    // A PID with three rank-ordered PSMs: all three survive the round trip.
    let (prot_ids, pep_ids) = multi_rank_ids();
    let tmp = TempDir::new(false).unwrap();
    let dir = new_dir(&tmp, "multirank");

    PSMArrowIO::export_to_parquet(&prot_ids, &pep_ids, &dir, true, &config()).unwrap();

    let mut prot_ids_in = Vec::new();
    let mut pep_ids_in = Vec::new();
    PSMArrowIO::import_from_parquet(&dir, &mut prot_ids_in, &mut pep_ids_in).unwrap();

    assert_eq!(pep_ids_in.len(), 1);
    assert_eq!(pep_ids_in[0].hits.len(), 3);
    for (i, hit) in pep_ids_in[0].hits.iter().enumerate() {
        assert_eq!(hit.rank, i as u32 + 1);
    }
    let similar = |a: f64, b: f64| (a - b).abs() <= 1e-5;
    assert!(similar(pep_ids_in[0].hits[0].score, 1.0));
    assert!(similar(pep_ids_in[0].hits[1].score, 0.5));
    assert!(similar(pep_ids_in[0].hits[2].score, 1.0 / 3.0));
}

// START_SECTION(([EXTRA] exportToParquet rejects duplicate ProteinIdentification identifiers))
#[test]
fn export_rejects_duplicate_run_identifiers() {
    let (mut prot_ids, pep_ids) = build_minimal_ids();
    prot_ids.push(ProteinIdentification {
        identifier: "run_1".into(), // duplicate run identifier
        ..ProteinIdentification::default()
    });
    let tmp = TempDir::new(false).unwrap();
    let dir = new_dir(&tmp, "duplicate");

    assert!(matches!(
        PSMArrowIO::export_to_parquet(&prot_ids, &pep_ids, &dir, true, &config()),
        Err(Error::InvalidValue(_))
    ));
    // Rejected before anything is written: the check precedes the directory
    // (PSMArrowIO.cpp:58-62).
    assert!(!dir.exists());
}

// START_SECTION([EXTRA] export/import - the peptidoform column carries the chemistry and the definition restores the name)
#[test]
fn peptidoform_carries_chemistry_and_definition_restores_the_name() {
    let mut registry = ModificationsDB::global().clone();
    registry
        .register_definition(&definition("TestPSM:Adduct", 'K', "C9H11N2O8P"))
        .unwrap();
    let (prots, mut peps) = build_minimal_ids();
    peps[0].hits[0].sequence =
        AASequence::parse_with_registry("AEADNLDDK(TestPSM:Adduct)K", &registry).unwrap();

    let tmp = TempDir::new(false).unwrap();
    let dir = new_dir(&tmp, "adduct");
    PSMArrowIO::export_to_parquet_with_registry(&prots, &peps, &dir, true, &config(), &registry)
        .unwrap();
    assert!(column_contains(
        &dir.join(PSMS),
        "peptidoform",
        "AEADNLDDK[Formula:C9H11N2O8P1|INFO:TestPSM:Adduct]K"
    ));
    assert_eq!(
        count_sp_meta_key(&dir.join(SEARCH_PARAMS), MODIFICATION_DEFINITIONS),
        1
    );

    let mut prots_in = Vec::new();
    let mut peps_in = Vec::new();
    PSMArrowIO::import_from_parquet_with_registry(&dir, &mut prots_in, &mut peps_in, &mut registry)
        .unwrap();
    assert_eq!(peps_in.len(), 1);
    let back = &peps_in[0].hits[0].sequence;
    assert_eq!(back.to_string(), "AEADNLDDK(TestPSM:Adduct)K"); // the name, not a mass bracket
    assert_eq!(back.formula().unwrap(), formula("C54H86N15O28P1"));
    assert_eq!(prots_in.len(), 1);
    assert!(
        prots_in[0]
            .search_parameters
            .metadata
            .contains_key(MODIFICATION_DEFINITIONS)
    );

    // an inbound value (e.g. from an idXML load) is replaced, not duplicated
    let dir2 = new_dir(&tmp, "adduct_again");
    PSMArrowIO::export_to_parquet_with_registry(
        &prots_in,
        &peps_in,
        &dir2,
        true,
        &config(),
        &registry,
    )
    .unwrap();
    assert_eq!(
        count_sp_meta_key(&dir2.join(SEARCH_PARAMS), MODIFICATION_DEFINITIONS),
        1
    );
}

// START_SECTION([EXTRA] import - definitions are registered from search_params.parquet)
#[test]
fn import_registers_definitions_from_search_params() {
    let mut registry = ModificationsDB::global().clone();
    assert!(!registry.has_defined_modification("TestPSM:Fresh"));
    let (mut prots, peps) = build_minimal_ids();
    prots[0].search_parameters.metadata.insert(
        MODIFICATION_DEFINITIONS.into(),
        definition("TestPSM:Fresh", 'K', "C2H2O")
            .to_definition_string()
            .unwrap()
            .into(),
    );

    let tmp = TempDir::new(false).unwrap();
    let dir = new_dir(&tmp, "fresh");
    PSMArrowIO::export_to_parquet_with_registry(&prots, &peps, &dir, true, &config(), &registry)
        .unwrap();
    assert!(!registry.has_defined_modification("TestPSM:Fresh")); // exporting registers nothing
    let mut prots_in = Vec::new();
    let mut peps_in = Vec::new();
    PSMArrowIO::import_from_parquet_with_registry(&dir, &mut prots_in, &mut peps_in, &mut registry)
        .unwrap();
    assert!(registry.has_defined_modification("TestPSM:Fresh"));
}

/// One identification of an `IDMerger_1_inputX.idXML` fixture: `(score,
/// sequence, charge, protein reference)` per hit, in file order.
type Hits = [(f64, &'static str, i32, Option<&'static str>)];

/// Imports a bundle the C++ `IDFileConverter` wrote from an idXML fixture and
/// requires that fixture's identifications back.
fn check_idmerger_input(bundle: &str, date: &str, accessions: [&str; 2], hits: &Hits) {
    let mut proteins = Vec::new();
    let mut peptides = Vec::new();
    PSMArrowIO::import_from_parquet(fixture(bundle), &mut proteins, &mut peptides).unwrap();

    assert_eq!(proteins.len(), 1);
    let run = &proteins[0];
    assert_eq!(run.score_type, "Mascot");
    assert!(run.higher_score_better);
    assert_eq!(run.search_engine, "");
    // synthesizeRunIdentifiers: `unknown` for the empty engine, then the date
    // (ProteinIdentificationArrowIO.cpp:1814-1819).
    assert!(
        run.identifier.starts_with(&format!("unknown_{date}_")),
        "{}",
        run.identifier
    );
    let protein_hits: Vec<(&str, f64)> = run
        .hits
        .iter()
        .map(|h| (h.accession.as_str(), h.score))
        .collect();
    assert_eq!(protein_hits, [(accessions[0], 40.0), (accessions[1], 40.0)]);

    assert_eq!(peptides.len(), 1);
    let peptide = &peptides[0];
    assert_eq!(peptide.identifier, run.identifier);
    assert_eq!(peptide.score_type, "Mascot");
    assert!(peptide.higher_score_better);
    assert_eq!(peptide.mz, Some(623.58));
    assert_eq!(peptide.rt, Some(-1.0));
    let got: Vec<(f64, String, i32, Vec<String>)> = peptide
        .hits
        .iter()
        .map(|h| {
            (
                h.score,
                h.sequence.to_string(),
                h.charge,
                h.evidences
                    .iter()
                    .map(|e| e.protein_accession.clone())
                    .collect(),
            )
        })
        .collect();
    let expected: Vec<(f64, String, i32, Vec<String>)> = hits
        .iter()
        .map(|&(score, sequence, charge, protein)| {
            (
                score,
                sequence.to_owned(),
                charge,
                protein.map(str::to_owned).into_iter().collect(),
            )
        })
        .collect();
    assert_eq!(got, expected);
}

/// `IDMerger_idparquet_input1.idparquet`, written by the C++ build from
/// `IDMerger_1_input1.idXML`.
#[test]
fn cpp_written_idmerger_input1_imports_its_idxml_source() {
    check_idmerger_input(
        "IDMerger_idparquet_input1.idparquet",
        "2006-10-07T06:30:59",
        ["APRT", "APRTT"],
        &[
            (40.5, "APRTPEPTIDER", 2, Some("APRT")),
            (30.5, "APRTPEPTIDERR", 2, Some("APRT")),
            (40.5, "APRTTPEPTIDER", 2, Some("APRTT")),
            (30.5, "APEPTIDER", 2, None),
            (30.5, "APEPTIDERR", 2, None),
        ],
    );
}

/// `IDMerger_idparquet_input2.idparquet`, written by the C++ build from
/// `IDMerger_1_input2.idXML`.
#[test]
fn cpp_written_idmerger_input2_imports_its_idxml_source() {
    check_idmerger_input(
        "IDMerger_idparquet_input2.idparquet",
        "2006-10-07T06:31:08",
        ["BPRT", "BPRTT"],
        &[
            (40.5, "BPRTPEPTIDER", 2, Some("BPRT")),
            (40.5, "BPRTTPEPTIDER", 2, Some("BPRTT")),
            (30.5, "BPRTTPEPTIDERR", 2, Some("BPRTT")),
            (30.5, "BPEPTIDER", 2, None),
        ],
    );
}

/// Every C++-written idparquet fixture imports with the format's invariants:
/// one identification per distinct `peptide_identification_index`, one hit
/// per PSM row, one run per search-parameter row, and every identification
/// stamped with the fresh identifier of the run it was stored under.
#[test]
fn cpp_written_bundles_import_with_the_format_invariants() {
    for bundle in [
        "IDMerger_idparquet_input1.idparquet",
        "IDMerger_idparquet_input2.idparquet",
        "IDRipper_4_input.idparquet",
        "THIRDPARTY/PercolatorAdapter_idparquet_in.idparquet",
    ] {
        let dir = fixture(bundle);
        let psms = read(&dir.join(PSMS));
        let search_params = read(&dir.join(SEARCH_PARAMS));
        let stored_runs: Vec<String> = {
            let column = search_params
                .column_by_name("run_identifier")
                .unwrap()
                .as_string::<i32>();
            (0..column.len())
                .map(|i| column.value(i).to_owned())
                .collect()
        };
        let index = psms
            .column_by_name("peptide_identification_index")
            .unwrap()
            .as_primitive::<Int32Type>();
        let distinct: BTreeSet<i32> = index.iter().flatten().collect();
        let psm_runs = psms
            .column_by_name("run_identifier")
            .unwrap()
            .as_string::<i32>();
        assert!(dir.join(PROTEINS).exists() && dir.join(PROTEIN_GROUPS).exists());

        let mut proteins = Vec::new();
        let mut peptides = Vec::new();
        PSMArrowIO::import_from_parquet(&dir, &mut proteins, &mut peptides).unwrap();

        assert_eq!(proteins.len(), stored_runs.len(), "{bundle}");
        assert_eq!(peptides.len(), distinct.len(), "{bundle}");
        let hits: usize = peptides.iter().map(|p| p.hits.len()).sum();
        assert_eq!(hits, psms.num_rows(), "{bundle}");
        let fresh: BTreeSet<&str> = proteins.iter().map(|p| p.identifier.as_str()).collect();
        for stored in &stored_runs {
            assert!(!fresh.contains(stored.as_str()), "{bundle}: {stored} kept");
        }
        // Rows are grouped by identification in first-seen order, so the
        // first row's run is the first identification's.
        let first_run = stored_runs
            .iter()
            .position(|r| r == psm_runs.value(0))
            .unwrap();
        assert_eq!(
            peptides[0].identifier, proteins[first_run].identifier,
            "{bundle}"
        );
        for peptide in &peptides {
            assert!(fresh.contains(peptide.identifier.as_str()), "{bundle}");
        }
    }
}

/// With the registry form, a failed import leaves the caller's registry
/// unchanged; the source may already have registered the bundle's definitions
/// process-wide when a later table fails.
#[test]
fn failed_import_leaves_the_registry_unchanged() {
    let (mut prots, peps) = build_minimal_ids();
    prots[0].search_parameters.metadata.insert(
        MODIFICATION_DEFINITIONS.into(),
        definition("TestPSM:Staged", 'K', "C2H2O")
            .to_definition_string()
            .unwrap()
            .into(),
    );
    let tmp = TempDir::new(false).unwrap();
    let dir = new_dir(&tmp, "staged");
    let mut registry = ModificationsDB::global().clone();
    PSMArrowIO::export_to_parquet_with_registry(&prots, &peps, &dir, true, &config(), &registry)
        .unwrap();
    // A psms.parquet that is not Parquet fails after the protein side (and
    // its definitions) were read.
    std::fs::write(dir.join(PSMS), b"not parquet").unwrap();
    let mut prots_in = Vec::new();
    let mut peps_in = Vec::new();
    assert!(
        PSMArrowIO::import_from_parquet_with_registry(
            &dir,
            &mut prots_in,
            &mut peps_in,
            &mut registry
        )
        .is_err()
    );
    assert!(!registry.has_defined_modification("TestPSM:Staged"));
    assert!(prots_in.is_empty() && peps_in.is_empty());
}

/// A path that is not a directory is refused with the source's message.
#[test]
fn import_of_a_missing_directory_is_refused() {
    let tmp = TempDir::new(false).unwrap();
    let dir = new_dir(&tmp, "absent");
    let mut proteins = Vec::new();
    let mut peptides = Vec::new();
    match PSMArrowIO::import_from_parquet(&dir, &mut proteins, &mut peptides) {
        // PSMArrowIO.cpp:108.
        Err(Error::InvalidValue(message)) => assert_eq!(
            message,
            format!("PSMArrowIO: not a directory: {}", dir.display())
        ),
        other => panic!("expected the not-a-directory refusal, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// Executed differential against the Linux x86_64 Release build
// ---------------------------------------------------------------------------
//
// `tests/data/psm_arrow_io/cpp_release_records.tsv`, `cpp_release_stderr.txt`
// and the three bundles `tests/data/psm_arrow_io/export_{minimal,decoy,multirank}`
// were produced by the oracle driver `../oracle/integ-F-psmarrowio` (outside
// the repository) linking the Linux x86_64 Release install
// `openms4-release-bc9cc12-c19e494-174b576` on ibminode06, run twice with
// identical records. The driver imports the four C++-written fixtures with
// `PSMArrowIO::importFromParquet` and records every field below; it also
// exports the class test's minimal, decoy and multi-rank sets with
// `PSMArrowIO::exportToParquet` and records their re-import. Doubles are the
// hex of their bits; a run identifier is recorded without its trailing
// `_<unique id>`.

fn release_records() -> Vec<(String, String, String)> {
    let text = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/psm_arrow_io/cpp_release_records.tsv"),
    )
    .unwrap();
    text.lines()
        .map(|line| {
            let mut parts = line.splitn(3, '\t');
            let set = parts.next().unwrap().to_owned();
            let key = parts.next().unwrap().to_owned();
            let value = unescape(parts.next().unwrap_or(""));
            (set, key, value)
        })
        .collect()
}

/// The driver's `esc` undone (`\t`, `\n`, `\\`, `\xNN`).
fn unescape(text: &str) -> String {
    let bytes = text.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'\\' && i + 1 < bytes.len() {
            match bytes[i + 1] {
                b't' => out.push(b'\t'),
                b'n' => out.push(b'\n'),
                b'\\' => out.push(b'\\'),
                b'x' if i + 3 < bytes.len() => {
                    let hex = std::str::from_utf8(&bytes[i + 2..i + 4]).unwrap();
                    out.push(u8::from_str_radix(hex, 16).unwrap());
                    i += 4;
                    continue;
                }
                other => {
                    out.push(b'\\');
                    out.push(other);
                }
            }
            i += 2;
        } else {
            out.push(bytes[i]);
            i += 1;
        }
    }
    String::from_utf8(out).unwrap()
}

fn bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

fn meta_value(value: &MetaValue) -> String {
    match value.data() {
        MetaValueData::Empty => "e:".to_owned(),
        MetaValueData::String(s) => format!("s:{s}"),
        MetaValueData::Integer(i) => format!("i:{i}"),
        MetaValueData::Float(f) => format!("d:{}", bits(*f)),
        MetaValueData::StringList(l) => {
            format!(
                "sl:{}",
                l.iter().map(|s| format!("{s}|")).collect::<String>()
            )
        }
        MetaValueData::IntegerList(l) => {
            format!(
                "il:{}",
                l.iter().map(|i| format!("{i}|")).collect::<String>()
            )
        }
        MetaValueData::FloatList(l) => format!(
            "dl:{}",
            l.iter()
                .map(|f| format!("{}|", bits(*f)))
                .collect::<String>()
        ),
    }
}

fn flank(residue: FlankingResidue) -> char {
    match residue {
        FlankingResidue::Unknown => 'X',
        FlankingResidue::NTerminus => '[',
        FlankingResidue::CTerminus => ']',
        FlankingResidue::Residue(c) => c,
    }
}

fn position(value: Option<usize>) -> String {
    value.map_or_else(|| "-1".to_owned(), |v| v.to_string())
}

type Records = Vec<(String, String, String)>;

/// The driver's `dump`, field for field.
fn records_of(
    set: &str,
    proteins: &[ProteinIdentification],
    peptides: &[PeptideIdentification],
) -> Records {
    let mut out = Vec::new();
    let mut rec = |key: String, value: String| out.push((set.to_owned(), key, value));
    rec("runs".into(), proteins.len().to_string());
    for (i, p) in proteins.iter().enumerate() {
        let k = format!("run{i}");
        let id = match p.identifier.rfind('_') {
            Some(pos) => p.identifier[..pos].to_owned(),
            None => format!("{}|nounderscore", p.identifier),
        };
        rec(format!("{k}.id"), id);
        rec(format!("{k}.engine"), p.search_engine.clone());
        rec(
            format!("{k}.engine_version"),
            p.search_engine_version.clone(),
        );
        let date = p
            .date_time
            .as_deref()
            .and_then(|text| DateTime::parse(text).ok())
            .filter(|date| date.is_valid())
            .and_then(|date| date.format("yyyy-MM-ddThh:mm:ss").ok())
            .unwrap_or_else(|| "invalid".to_owned());
        rec(format!("{k}.date"), date);
        rec(format!("{k}.score_type"), p.score_type.clone());
        let hsb = if p.higher_score_better { "1" } else { "0" };
        rec(format!("{k}.hsb"), hsb.into());
        rec(format!("{k}.threshold"), bits(p.significance_threshold));
        rec(
            format!("{k}.enzyme"),
            p.search_parameters.digestion_enzyme.clone(),
        );
        rec(
            format!("{k}.paths"),
            p.primary_ms_run_paths
                .iter()
                .map(|s| format!("{s}|"))
                .collect(),
        );
        rec(format!("{k}.hits"), p.hits.len().to_string());
        for (j, h) in p.hits.iter().enumerate() {
            rec(
                format!("{k}.hit{j}"),
                format!("{}|{}", h.accession, bits(h.score)),
            );
        }
    }
    rec("peps".into(), peptides.len().to_string());
    for (n, pep) in peptides.iter().enumerate() {
        let k = format!("pep{n}");
        let run = proteins
            .iter()
            .position(|p| p.identifier == pep.identifier)
            .map_or(-1, |i| i as i64);
        rec(format!("{k}.run"), run.to_string());
        rec(
            format!("{k}.rt"),
            pep.rt.map_or_else(|| "none".to_owned(), bits),
        );
        rec(
            format!("{k}.mz"),
            pep.mz.map_or_else(|| "none".to_owned(), bits),
        );
        rec(format!("{k}.score_type"), pep.score_type.clone());
        let hsb = if pep.higher_score_better { "1" } else { "0" };
        rec(format!("{k}.hsb"), hsb.into());
        rec(format!("{k}.ref"), pep.spectrum_reference());
        rec(format!("{k}.hits"), pep.hits.len().to_string());
        for (j, h) in pep.hits.iter().enumerate() {
            let hk = format!("{k}.h{j}");
            rec(
                hk.clone(),
                format!("{}|{}|{}|{}", bits(h.score), h.sequence, h.charge, h.rank),
            );
            let evidences: String = h
                .evidences
                .iter()
                .map(|e| {
                    format!(
                        "{},{},{},{},{};",
                        e.protein_accession,
                        flank(e.aa_before),
                        flank(e.aa_after),
                        position(e.start),
                        position(e.end)
                    )
                })
                .collect();
            rec(format!("{hk}.ev"), evidences);
            // The source keeps a hit's rank as its "rank" meta value
            // (PeptideHit::getRank/setRank, PeptideHit.cpp:152-155 and
            // 341-346: set only when nonzero); this crate keeps it as the
            // `rank` field, so it is recorded among the meta values here.
            let mut metadata = h.metadata.clone();
            if h.rank != 0 && !metadata.contains_key("rank") {
                metadata.insert("rank".into(), i64::from(h.rank).into());
            }
            let metas: String = metadata
                .iter()
                .map(|(key, value)| format!("{key}={};", meta_value(value)))
                .collect();
            rec(format!("{hk}.meta"), metas);
        }
    }
    out
}

fn release_set(set: &str) -> Records {
    release_records()
        .into_iter()
        .filter(|(s, k, _)| s == set && k != "ok" && k != "export_ok")
        .collect()
}

fn assert_records(expected: &Records, actual: &Records) {
    for (e, a) in expected.iter().zip(actual) {
        assert_eq!(e, a);
    }
    assert_eq!(expected.len(), actual.len(), "{}", expected[0].0);
}

/// The four C++-written fixtures import as the Release build imports them,
/// with QPXFile's warnings in the order the Release build logged them.
#[test]
fn fixtures_import_as_the_release_build_imports_them() {
    let records = release_records();
    let mut warnings = Vec::new();
    for (set, bundle) in [
        (
            "IDMerger_idparquet_input1.idparquet",
            "IDMerger_idparquet_input1.idparquet",
        ),
        (
            "IDMerger_idparquet_input2.idparquet",
            "IDMerger_idparquet_input2.idparquet",
        ),
        ("IDRipper_4_input.idparquet", "IDRipper_4_input.idparquet"),
        (
            "PercolatorAdapter_idparquet_in.idparquet",
            "THIRDPARTY/PercolatorAdapter_idparquet_in.idparquet",
        ),
    ] {
        assert!(records.contains(&(set.to_owned(), "ok".to_owned(), "1".to_owned())));
        let mut proteins = Vec::new();
        let mut peptides = Vec::new();
        warnings.extend(
            PSMArrowIO::import_from_parquet(fixture(bundle), &mut proteins, &mut peptides).unwrap(),
        );
        assert_records(&release_set(set), &records_of(set, &proteins, &peptides));
    }
    let stderr = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/psm_arrow_io/cpp_release_stderr.txt"),
    )
    .unwrap();
    // Native difference, documented at the check in qpx_file.rs: the source
    // gives B, Z and X an empty formula (ResidueDB.cpp:236, 239, 245), so it
    // computes an m/z for such a peptidoform (each residue counting as -H2O)
    // and compares it; this crate has no m/z for a residue without chemistry
    // and does not compare. The Release build warned for exactly two such
    // peptidoforms (IDRipper_4_input's BBBB... rows).
    let without_chemistry = |line: &str| {
        let peptidoform = line.split('\'').nth(1).unwrap_or("");
        peptidoform.contains(['B', 'Z', 'X'])
    };
    let all: Vec<&str> = stderr
        .lines()
        .filter(|l| l.starts_with("QPXFile:"))
        .collect();
    assert_eq!(all.iter().filter(|l| without_chemistry(l)).count(), 2);
    let cpp: Vec<&str> = all.into_iter().filter(|l| !without_chemistry(l)).collect();
    let rust: Vec<&str> = warnings
        .iter()
        .map(String::as_str)
        .filter(|l| l.starts_with("QPXFile:"))
        .collect();
    assert_eq!(rust, cpp);
}

/// The Release build's own bundles of the class-test data import here as
/// they import there, and this port's bundles of the same data import to the
/// same records.
#[test]
fn release_bundles_and_this_ports_bundles_import_alike() {
    type Build = fn() -> (Vec<ProteinIdentification>, Vec<PeptideIdentification>);
    let builders: [(&str, Build); 3] = [
        ("minimal", build_minimal_ids),
        ("decoy", decoy_ids),
        ("multirank", multi_rank_ids),
    ];
    let records = release_records();
    let tmp = TempDir::new(false).unwrap();
    for (name, build) in builders {
        let set = format!("export_{name}");
        assert!(records.contains(&(set.clone(), "export_ok".to_owned(), "1".to_owned())));
        let expected = release_set(&set);

        let cpp_bundle = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/psm_arrow_io")
            .join(&set);
        let mut proteins = Vec::new();
        let mut peptides = Vec::new();
        PSMArrowIO::import_from_parquet(&cpp_bundle, &mut proteins, &mut peptides).unwrap();
        assert_records(&expected, &records_of(&set, &proteins, &peptides));

        let (prots, peps) = build();
        let dir = new_dir(&tmp, name);
        PSMArrowIO::export_to_parquet(&prots, &peps, &dir, true, &config()).unwrap();
        let mut proteins = Vec::new();
        let mut peptides = Vec::new();
        PSMArrowIO::import_from_parquet(&dir, &mut proteins, &mut peptides).unwrap();
        assert_records(&expected, &records_of(&set, &proteins, &peptides));
    }
}
