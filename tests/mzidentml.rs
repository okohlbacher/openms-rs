// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Coverage for `FORMAT/MzIdentMLFile.h` and
//! `FORMAT/HANDLERS/MzIdentMLHandler.h`.
//!
//! Every literal taken from `MzIdentMLFile_test.cpp` or from one of the ten
//! unmodified upstream fixtures is transcribed source review (tier 3): the
//! counts, accessions, score types, scores, spectrum references, the
//! `database.fasta`/`MSDB` search databases, `missed_cleavages` 1000, the
//! `Carbamidomethyl (C)` / `Xlink:DTSSP[88] (Protein N-term)` fixed
//! modifications, the `Acetyl (N-term)` variable modification, the 0.5
//! significance threshold, the modification inference cases of issue #5443 and
//! the cross-linking positions, masses, chains and fragment annotations of the
//! two XLMS fixtures. No C++ mzIdentML operation was run, so nothing here is a tier 1
//! differential.
//!
//! The synthetic documents - duplicate ids, dangling references, an
//! unreferenced `SpectrumIdentificationList`, a prefixed namespace, non-ASCII
//! text, entity references and CDATA, the resource ceilings, the C-terminal
//! modification location and the `ProteinDetectionList` - are independently
//! derived from the mzIdentML 1.3.0 schema (tier 4), because no upstream
//! fixture reaches them.
//!
//! All 15 upstream sections are covered, the six cross-linking ones included;
//! `docs/MZIDENTML_SUPPORT.md` carries the section table and the divergences
//! the cross-linking path documents.

#![cfg(feature = "idxml")]

use openms::Error;
use openms::chemistry::ModificationsDB;
use openms::comparison::Tolerance;
use openms::format::mzidentml::{
    self, MzIdentMLDocument, ReadOptions, SCHEMA_VERSION, WriteOptions,
};
use openms::identification::{
    EnzymeTermSpecificity, FlankingResidue, PeptideHit, PeptideIdentification,
};
use openms::system::file::TempDir;
use std::path::PathBuf;

const WHOLE: &str = "mzidentml_whole.mzid";
const MSGF: &str = "mzidentml_msgf_mini.mzid";
const MISSING_LOCATION: &str = "mzidentml_missing_mod_location.mzid";
const THREE_RUNS: &str = "mzidentml_3runs.mzid";
const CROSSLINKING: &str = "mzidentml_crosslinking_v1_3.mzid";
const XLMS_LABELLED: &str = "mzidentml_xlms_labelled.mzid";
const XLMS_UNLABELLED: &str = "mzidentml_xlms_unlabelled.mzid";
const NONCOVALENT: &str = "mzidentml_noncov_assoc_v1_3.mzid";
const EDC: &str = "mzidentml_xlink_edc_v1_3.mzid";
const MULTI_SPECTRA: &str = "mzidentml_multi_spectra_v1_3.mzid";

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn load(name: &str) -> MzIdentMLDocument {
    mzidentml::load(data(name)).expect("upstream mzIdentML fixture loads")
}

fn read(text: &str) -> Result<MzIdentMLDocument, Error> {
    mzidentml::read(text.as_bytes())
}

/// A minimal schema-shaped document with one run, one peptide and one PSM.
fn document(library: &str, results: &str) -> String {
    format!(
        r#"<?xml version="1.0" encoding="UTF-8"?>
<MzIdentML xmlns="http://psidev.info/psi/pi/mzIdentML/1.3" version="1.3.0" id="TEST">
<SequenceCollection>
{library}
</SequenceCollection>
<AnalysisCollection>
 <SpectrumIdentification id="SI" spectrumIdentificationProtocol_ref="SIP" spectrumIdentificationList_ref="SIL">
  <InputSpectra spectraData_ref="SDAT"/>
  <SearchDatabaseRef searchDatabase_ref="SDB"/>
 </SpectrumIdentification>
</AnalysisCollection>
<AnalysisProtocolCollection>
 <SpectrumIdentificationProtocol id="SIP" analysisSoftware_ref="SOF">
  <SearchType><cvParam accession="MS:1001083" cvRef="PSI-MS" name="ms-ms search"/></SearchType>
  <Threshold><cvParam accession="MS:1001494" cvRef="PSI-MS" name="no threshold"/></Threshold>
 </SpectrumIdentificationProtocol>
</AnalysisProtocolCollection>
<DataCollection>
 <Inputs>
  <SearchDatabase location="db.fasta" id="SDB"><DatabaseName><userParam name="db.fasta"/></DatabaseName></SearchDatabase>
  <SpectraData location="run.mzML" id="SDAT"/>
 </Inputs>
 <AnalysisData>
  <SpectrumIdentificationList id="SIL">
{results}
  </SpectrumIdentificationList>
 </AnalysisData>
</DataCollection>
</MzIdentML>
"#
    )
}

/// One `SpectrumIdentificationResult` with a single Mascot-scored PSM.
fn result_with(item_attributes: &str, children: &str) -> String {
    format!(
        r#"   <SpectrumIdentificationResult spectraData_ref="SDAT" spectrumID="scan=1" id="SIR">
    <SpectrumIdentificationItem passThreshold="true" rank="1" peptide_ref="PEP" {item_attributes} id="SII">
{children}
     <cvParam accession="MS:1001171" cvRef="PSI-MS" name="Mascot:score" value="42.5"/>
    </SpectrumIdentificationItem>
   </SpectrumIdentificationResult>"#
    )
}

// ---------------------------------------------------------------------------
// Upstream section: MzIdentMLFile() / ~MzIdentMLFile()
// ---------------------------------------------------------------------------

#[test]
fn adapter_defaults_match_the_source_constructor() {
    // XMLFile("/SCHEMAS/mzIdentML1.3.0.xsd", "1.3.0")
    assert_eq!(SCHEMA_VERSION, "1.3.0");
    let options = ReadOptions::default();
    assert!(options.max_xml_bytes > 0 && options.max_depth > 0);
    assert!(WriteOptions::default().creation_date.is_none());
}

// ---------------------------------------------------------------------------
// Upstream section: void load(...)
// ---------------------------------------------------------------------------

#[test]
fn load_msgf_mini_matches_upstream_literals() {
    let document = load(MSGF);
    let proteins = &document.protein_identifications;
    let peptides = &document.peptide_identifications;
    assert_eq!(proteins.len(), 2);
    assert_eq!(proteins[0].hits.len(), 2);
    assert_eq!(proteins[1].hits.len(), 1);
    assert_eq!(peptides.len(), 5);
    for peptide in peptides {
        assert_eq!(peptide.hits.len(), 1);
    }

    assert_eq!(proteins[0].search_engine, "MS-GF+");
    assert_eq!(proteins[0].search_engine_version, "Beta (v9979)");
    // The upstream test asserts a nonzero date, which only holds because the
    // DOM handler substitutes DateTime::now() for the absent activityDate.
    assert_eq!(proteins[0].date_time, None);
    let parameters = &proteins[0].search_parameters;
    assert_eq!(parameters.database, "database.fasta");
    assert_eq!(parameters.missed_cleavages, 1000);
    assert_eq!(
        parameters.fixed_modifications,
        vec![
            "Carbamidomethyl (C)".to_owned(),
            "Xlink:DTSSP[88] (Protein N-term)".to_owned()
        ]
    );
    assert_eq!(parameters.fragment_tolerance, Tolerance::Absolute(0.0));
    assert_eq!(parameters.precursor_tolerance, Tolerance::Ppm(20.0));
    assert_eq!(parameters.charges, "2-3");
    assert_eq!(parameters.enzyme_specificity, EnzymeTermSpecificity::Full);
    assert_eq!(parameters.digestion_enzyme, "Trypsin");

    assert_eq!(proteins[0].hits[0].accession, "sp|P0A9K9|SLYD_ECOLI");
    assert_eq!(proteins[0].hits[0].sequence, "");
    assert_eq!(proteins[0].hits[1].accession, "sp|P0A786|PYRB_ECOLI");
    assert_eq!(proteins[0].hits[1].sequence, "");

    let expected = [
        (
            "MS-GF:RawScore",
            195.0,
            "LATEFSGNVPVLNAGDGSNQHPTQTLLDLFTIQETQGR",
            "controllerType=0 controllerNumber=1 scan=32805",
        ),
        (
            "MS-GF:RawScore",
            182.0,
            "FLAETDQGPVPVEITAVEDDHVVVDGNHMLAGQNLK",
            "controllerType=0 controllerNumber=1 scan=26090",
        ),
        (
            "MS-GF:RawScore",
            191.0,
            "FLAETDQGPVPVEITAVEDDHVVVDGNHMLAGQNLK",
            "controllerType=0 controllerNumber=1 scan=26157",
        ),
        (
            "MS-GF:RawScore",
            211.0,
            "VGAGPFPTELFDETGEFLC(Carbamidomethyl)K",
            "controllerType=0 controllerNumber=1 scan=15094",
        ),
    ];
    for (index, (score_type, score, sequence, reference)) in expected.iter().enumerate() {
        let peptide = &peptides[index];
        assert_eq!(&peptide.score_type, score_type, "score type {index}");
        assert_eq!(peptide.hits[0].score, *score, "score {index}");
        assert_eq!(
            peptide.hits[0].sequence.to_string(),
            *sequence,
            "sequence {index}"
        );
        assert_eq!(
            peptide.spectrum_reference(),
            *reference,
            "reference {index}"
        );
    }
    // Every PSM keeps its evidence, and the 1-based file positions became
    // 0-based OpenMS positions.
    let evidence = &peptides[0].hits[0].evidences[0];
    assert_eq!(evidence.protein_accession, "sp|P0A786|PYRB_ECOLI");
    assert_eq!(evidence.start, Some(114));
    assert_eq!(evidence.end, Some(151));
    assert_eq!(evidence.aa_before, FlankingResidue::Residue('R'));
    assert_eq!(evidence.aa_after, FlankingResidue::Residue('L'));
}

#[test]
fn load_replaces_rather_than_accumulates() {
    let mut document = load(MSGF);
    mzidentml::load_into(data(MSGF), &mut document).expect("second load replaces");
    assert_eq!(document.protein_identifications.len(), 2);
    assert_eq!(document.peptide_identifications.len(), 5);
}

// ---------------------------------------------------------------------------
// Upstream section: [EXTRA] read mzIdentML Modification without the optional
// location attribute (issue #5443)
// ---------------------------------------------------------------------------

#[test]
fn missing_modification_location_is_inferred_or_skipped() {
    let document = load(MISSING_LOCATION);
    let peptides = &document.peptide_identifications;
    assert_eq!(peptides.len(), 5);
    for peptide in peptides {
        assert!(!peptide.hits.is_empty());
    }
    // 1) N-terminal Acetyl without location: inferred as N-terminal.
    let acetyl = &peptides[0].hits[0].sequence;
    assert_eq!(acetyl.as_str(), "PEPTIDEK");
    assert!(acetyl.n_terminal_modification().is_some());
    assert!(acetyl.c_terminal_modification().is_none());
    assert_eq!(
        acetyl.n_terminal_modification().map(|m| m.name()),
        Some("Acetyl")
    );
    // 2) Oxidation without location: not uniquely terminal, so skipped.
    let oxidation = &peptides[1].hits[0].sequence;
    assert_eq!(oxidation.as_str(), "PEPTIDEM");
    assert!(!oxidation.is_modified());
    assert_eq!(oxidation.to_string(), "PEPTIDEM");
    // 3) Amidated without location: exclusively C-terminal.
    let amidated = &peptides[2].hits[0].sequence;
    assert_eq!(amidated.as_str(), "PEPTIDER");
    assert!(amidated.c_terminal_modification().is_some());
    assert!(amidated.n_terminal_modification().is_none());
    assert_eq!(
        amidated.c_terminal_modification().map(|m| m.name()),
        Some("Amidated")
    );
    // 4) Control: a valid internal location is still parsed.
    assert_eq!(
        peptides[3].hits[0].sequence.to_string(),
        "PEPTIDEC(Carbamidomethyl)K"
    );
    // 5) Acetyl with residues="K" and no location: residue-specific, skipped.
    let acetyl_on_k = &peptides[4].hits[0].sequence;
    assert_eq!(acetyl_on_k.as_str(), "PEPTIKDE");
    assert!(!acetyl_on_k.is_modified());
}

// ---------------------------------------------------------------------------
// Upstream section: void store(...)
// ---------------------------------------------------------------------------

fn store_and_reload(document: &MzIdentMLDocument, name: &str) -> MzIdentMLDocument {
    let directory = TempDir::new(false).expect("temporary directory");
    let path = directory.path().join(name);
    mzidentml::store(&path, document).expect("store writes mzIdentML");
    let reloaded = mzidentml::load(&path).expect("stored mzIdentML loads");
    // A second store of the same document is byte-identical: every id this
    // writer emits is positional, unlike UniqueIdGenerator.
    let again = directory.path().join("again.mzid");
    mzidentml::store(&again, document).expect("second store writes mzIdentML");
    assert_eq!(
        std::fs::read(&path).unwrap(),
        std::fs::read(&again).unwrap(),
        "store is reproducible"
    );
    reloaded
}

#[test]
fn store_round_trip_preserves_whole_document() {
    let original = load(WHOLE);
    let reloaded = store_and_reload(&original, "whole.mzid");
    assert_eq!(
        reloaded.protein_identifications.len(),
        original.protein_identifications.len()
    );
    assert_eq!(
        reloaded.peptide_identifications.len(),
        original.peptide_identifications.len()
    );
    let before = &original.protein_identifications[0];
    let after = &reloaded.protein_identifications[0];
    assert_eq!(after.hits.len(), before.hits.len());
    assert_eq!(after.search_engine, before.search_engine);
    assert_eq!(after.search_engine_version, before.search_engine_version);
    assert_eq!(after.date_time, before.date_time);
    assert_eq!(after.date_time.as_deref(), Some("2006-01-12T12:13:14"));
    assert_eq!(
        after.search_parameters.database,
        before.search_parameters.database
    );
    assert_eq!(
        after.search_parameters.database_version,
        before.search_parameters.database_version
    );
    assert_eq!(
        after.search_parameters.digestion_enzyme,
        before.search_parameters.digestion_enzyme
    );
    assert_eq!(
        after.search_parameters.charges,
        before.search_parameters.charges
    );
    assert_eq!(
        after.search_parameters.mass_type,
        before.search_parameters.mass_type
    );
    assert_eq!(
        after.search_parameters.fragment_tolerance,
        before.search_parameters.fragment_tolerance
    );
    assert_eq!(
        after.search_parameters.precursor_tolerance,
        before.search_parameters.precursor_tolerance
    );
    assert_eq!(
        after.search_parameters.variable_modifications,
        before.search_parameters.variable_modifications
    );
    assert_eq!(
        after
            .search_parameters
            .variable_modifications
            .last()
            .map(String::as_str),
        Some("Acetyl (N-term)")
    );
    assert_eq!(
        after.search_parameters.fixed_modifications,
        before.search_parameters.fixed_modifications
    );
    for (index, (after, before)) in after.hits.iter().zip(&before.hits).enumerate() {
        assert_eq!(after.accession, before.accession, "protein hit {index}");
        assert_eq!(after.sequence, before.sequence, "protein sequence {index}");
    }
    for (index, (after, before)) in reloaded
        .peptide_identifications
        .iter()
        .zip(&original.peptide_identifications)
        .enumerate()
    {
        assert_eq!(after.score_type, before.score_type, "score type {index}");
        assert_eq!(
            after.higher_score_better, before.higher_score_better,
            "orientation {index}"
        );
        assert_eq!(after.mz, before.mz, "m/z {index}");
        assert_eq!(after.rt, before.rt, "RT {index}");
        assert_eq!(
            after.spectrum_reference(),
            before.spectrum_reference(),
            "reference {index}"
        );
        assert_eq!(after.hits.len(), before.hits.len(), "hit count {index}");
        for (hit_index, (after, before)) in after.hits.iter().zip(&before.hits).enumerate() {
            assert_eq!(after.score, before.score, "score {index}/{hit_index}");
            assert_eq!(
                after.sequence, before.sequence,
                "sequence {index}/{hit_index}"
            );
            assert_eq!(after.charge, before.charge, "charge {index}/{hit_index}");
            assert_eq!(
                after.evidences.len(),
                before.evidences.len(),
                "evidence count {index}/{hit_index}"
            );
            for (after, before) in after.evidences.iter().zip(&before.evidences) {
                assert_eq!(after.start, before.start);
                assert_eq!(after.end, before.end);
                assert_eq!(after.aa_before, before.aa_before);
                assert_eq!(after.aa_after, before.aa_after);
            }
        }
    }
}

#[test]
fn store_round_trip_preserves_modified_peptides() {
    let original = load(MSGF);
    let reloaded = store_and_reload(&original, "msgf.mzid");
    assert_eq!(
        reloaded.peptide_identifications[3].hits[0]
            .sequence
            .to_string(),
        "VGAGPFPTELFDETGEFLC(Carbamidomethyl)K"
    );
    // The writer emits the protein-terminal specificity rules the source
    // leaves as a TODO, so a protein N-terminal modification survives.
    assert_eq!(
        reloaded.protein_identifications[0]
            .search_parameters
            .fixed_modifications,
        original.protein_identifications[0]
            .search_parameters
            .fixed_modifications
    );
    assert_eq!(
        reloaded.protein_identifications[0]
            .search_parameters
            .enzyme_specificity,
        EnzymeTermSpecificity::Full
    );
    // Evidence positions must not drift: the source adds one on write without
    // removing one on read, so its own round trip shifts them.
    let evidence = &reloaded.peptide_identifications[0].hits[0].evidences[0];
    assert_eq!(evidence.start, Some(114));
    assert_eq!(evidence.end, Some(151));
}

#[test]
fn store_rejects_a_foreign_extension() {
    let document = load(WHOLE);
    let directory = TempDir::new(false).expect("temporary directory");
    let error = mzidentml::store(directory.path().join("out.idXML"), &document).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    assert!(!directory.path().join("out.idXML").exists());
}

// ---------------------------------------------------------------------------
// Upstream section: [EXTRA] multiple runs
// ---------------------------------------------------------------------------

#[test]
fn three_runs_round_trip() {
    let original = load(THREE_RUNS);
    assert_eq!(original.protein_identifications.len(), 3);
    let reloaded = store_and_reload(&original, "three.mzid");
    assert_eq!(reloaded.protein_identifications.len(), 3);
    for index in 0..3 {
        assert_eq!(
            reloaded.protein_identifications[index].hits.len(),
            original.protein_identifications[index].hits.len(),
            "run {index}"
        );
    }
    assert_eq!(
        reloaded.protein_identifications[0]
            .search_parameters
            .precursor_tolerance,
        Tolerance::Ppm(20.0)
    );
    // Three SpectraData inputs stay distinct across the round trip.
    let locations: Vec<String> = reloaded
        .protein_identifications
        .iter()
        .map(|run| {
            run.metadata
                .get("spectra_data")
                .and_then(|value| value.as_string_list().ok())
                .and_then(|list| list.first().cloned())
                .unwrap_or_default()
        })
        .collect();
    assert_eq!(
        locations,
        vec![
            "/some/path/file1.mzML".to_owned(),
            "/some/path/file2.mzML".to_owned(),
            "/some/path/file3.mzML".to_owned()
        ]
    );
}

// ---------------------------------------------------------------------------
// Upstream section: [EXTRA] thresholds
// ---------------------------------------------------------------------------

#[test]
fn thresholds_round_trip() {
    let mut document = load(WHOLE);
    assert_eq!(document.protein_identifications.len(), 1);
    assert_eq!(
        document.protein_identifications[0].significance_threshold,
        0.5
    );
    assert_eq!(document.peptide_identifications.len(), 3);
    let index = document
        .peptide_identifications
        .iter()
        .position(|p| p.spectrum_reference() == "17")
        .expect("spectrum 17 is present");
    {
        let identification = &document.peptide_identifications[index];
        assert_eq!(identification.hits.len(), 2);
        for hit in &identification.hits {
            assert_eq!(
                hit.metadata
                    .get("pass_threshold")
                    .and_then(|v| v.as_str().ok()),
                Some("false")
            );
        }
    }
    // Add a hit that passes the threshold and has no explicit pass_threshold.
    let mut extra = document.peptide_identifications[index].hits[1].clone();
    extra.metadata.remove("pass_threshold");
    extra.sequence = openms::chemistry::AASequence::parse("TESTER").unwrap();
    extra.score = 0.4;
    document.peptide_identifications[index].hits.push(extra);

    let reloaded = store_and_reload(&document, "thresholds.mzid");
    assert_eq!(reloaded.peptide_identifications.len(), 3);
    let threshold = reloaded.protein_identifications[0].significance_threshold;
    assert_eq!(threshold, 0.5);
    let identification = reloaded
        .peptide_identifications
        .iter()
        .find(|p| p.spectrum_reference() == "17")
        .expect("spectrum 17 survives the round trip");
    assert_eq!(identification.hits.len(), 3);
    for hit in &identification.hits {
        let pass = hit
            .metadata
            .get("pass_threshold")
            .and_then(|v| v.as_str().ok())
            .expect("pass_threshold is written back");
        // q-values: lower is better, so a score below the threshold passes.
        let expected = if hit.score > threshold {
            "false"
        } else {
            "true"
        };
        assert_eq!(pass, expected, "score {}", hit.score);
    }
}

// ---------------------------------------------------------------------------
// Upstream section: [EXTRA] regression test for file loading on example files
// ---------------------------------------------------------------------------

#[test]
fn every_upstream_non_crosslinking_fixture_loads() {
    for name in [WHOLE, MSGF, MISSING_LOCATION, THREE_RUNS] {
        let document = load(name);
        assert!(
            !document.protein_identifications.is_empty(),
            "{name} has runs"
        );
        assert!(
            !document.peptide_identifications.is_empty(),
            "{name} has spectra"
        );
    }
}

// ---------------------------------------------------------------------------
// Upstream section: [EXTRA] compability issues
// ---------------------------------------------------------------------------
// The section is entirely commented out upstream and asserts nothing; the
// conditions its comments enumerate are exercised here instead.

#[test]
fn misplaced_elements_in_a_param_group_are_ignored() {
    // "Misplaced Elements ignored in ParamGroup": an unexpected child of a
    // SpectrumIdentificationResult neither fails the load nor becomes metadata.
    let library = r#"<DBSequence accession="P1" searchDatabase_ref="SDB" id="DBS"><Seq>PEPTIDEK</Seq></DBSequence>
<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>
<PeptideEvidence id="PEV" peptide_ref="PEP" dBSequence_ref="DBS" start="1" end="8" pre="K" post="A" isDecoy="false"/>"#;
    let results = r#"   <SpectrumIdentificationResult spectraData_ref="SDAT" spectrumID="scan=1" id="SIR">
    <SpectrumIdentificationItem passThreshold="true" rank="1" peptide_ref="PEP" chargeState="2" experimentalMassToCharge="500.5" id="SII">
     <PeptideEvidenceRef peptideEvidence_ref="PEV"/>
     <cvParam accession="MS:1001171" cvRef="PSI-MS" name="Mascot:score" value="42.5"/>
    </SpectrumIdentificationItem>
    <Measure id="stray"/>
   </SpectrumIdentificationResult>"#;
    let document = read(&document(library, results)).expect("stray element is ignored");
    let identification = &document.peptide_identifications[0];
    assert_eq!(identification.hits.len(), 1);
    assert!(!identification.metadata.contains_key("Measure"));
}

#[test]
fn a_psm_without_a_recognised_score_yields_no_hit() {
    // "Converting unknown score type to search engine specific score CV":
    // an item with no q-, e- or specific score produces no PeptideHit at all.
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "").replace(
        r#"<cvParam accession="MS:1001171" cvRef="PSI-MS" name="Mascot:score" value="42.5"/>"#,
        r#"<cvParam accession="MS:1001115" cvRef="PSI-MS" name="scan number(s)" value="7"/>"#,
    );
    let document = read(&document(library, &results)).expect("unscored item loads");
    assert!(document.peptide_identifications[0].hits.is_empty());
}

#[test]
fn a_psm_without_peptide_evidence_still_loads() {
    // "PSM without peptide evidences registered in the given search database":
    // the hit exists with no evidence, which is what reading idXML produces.
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let document = read(&document(library, &results)).expect("evidence-free PSM loads");
    let hit = &document.peptide_identifications[0].hits[0];
    assert!(hit.evidences.is_empty());
    assert!(document.protein_identifications[0].hits.is_empty());
}

#[test]
fn an_identification_without_rt_keeps_no_coordinate() {
    // "No RT" / "No MZ": both stay absent rather than becoming zero.
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let document = read(&document(library, &results)).expect("RT-free result loads");
    assert_eq!(document.peptide_identifications[0].rt, None);
    assert_eq!(document.peptide_identifications[0].mz, Some(500.5));
}

#[test]
fn evidence_without_positions_keeps_them_unknown() {
    // "PeptideEvidence without reference to the positional in originating
    // sequence": start and end stay None, and the writer omits both.
    let library = r#"<DBSequence accession="P1" searchDatabase_ref="SDB" id="DBS"/>
<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>
<PeptideEvidence id="PEV" peptide_ref="PEP" dBSequence_ref="DBS"/>"#;
    let results = result_with(
        "chargeState=\"2\" experimentalMassToCharge=\"500.5\"",
        "     <PeptideEvidenceRef peptideEvidence_ref=\"PEV\"/>",
    );
    let document = read(&document(library, &results)).expect("position-free evidence loads");
    let evidence = &document.peptide_identifications[0].hits[0].evidences[0];
    assert_eq!(evidence.start, None);
    assert_eq!(evidence.end, None);
    assert_eq!(evidence.aa_before, FlankingResidue::Unknown);
    let mut text = Vec::new();
    mzidentml::write(&mut text, &document).expect("writes without positions");
    let text = String::from_utf8(text).unwrap();
    assert!(!text.contains(" start="), "{text}");
}

// ---------------------------------------------------------------------------
// Cross-linking (MS:1002494): the six upstream XL sections
// ---------------------------------------------------------------------------

/// A hit's metadata value as text, which is how the upstream sections compare
/// the OpenPepXL user parameters.
fn meta(hit: &PeptideHit, key: &str) -> String {
    hit.metadata
        .get(key)
        .map(ToString::to_string)
        .unwrap_or_default()
}

/// Panics unless the value at `key` parses as a number close to `expected`.
fn meta_number(hit: &PeptideHit, key: &str, expected: f64) {
    let text = meta(hit, key);
    let value: f64 = text
        .parse()
        .unwrap_or_else(|_| panic!("{key} = {text:?} is not a number"));
    assert!(
        (value - expected).abs() < 1e-6,
        "{key} = {value}, expected {expected}"
    );
}

#[test]
fn the_crosslinking_marker_selects_the_crosslinking_path() {
    // The source scans every AdditionalSearchParams for MS:1002494 before it
    // reads anything else; a document that declares it takes the XL path even
    // when no peptide carries a cross-link, and every run is then tagged with
    // the term the writer tests in turn.
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let text = document(library, &results).replace(
        "<SearchType>",
        r#"<AdditionalSearchParams><cvParam accession="MS:1002494" cvRef="PSI-MS" name="crosslinking search"/></AdditionalSearchParams><SearchType>"#,
    );
    let document = read(&text).expect("a cross-linking document is read, not refused");
    let run = &document.protein_identifications[0];
    assert_eq!(
        run.metadata
            .get("SpectrumIdentificationProtocol")
            .map(ToString::to_string)
            .unwrap_or_default(),
        "MS:1002494"
    );
    // OPXLHelper::addPercolatorFeatureList runs on the first run of every
    // cross-linking document, whether or not it carries a cross-link.
    assert_eq!(
        run.search_parameters
            .metadata
            .get("feature_extractor")
            .map(ToString::to_string)
            .unwrap_or_default(),
        "TOPP_PSMFeatureExtractor"
    );
    // No crosslink donor: the source falls back to the linear item path and
    // then still collapses the result per spectrum.
    let identification = &document.peptide_identifications[0];
    assert_eq!(identification.hits.len(), 1);
    assert_eq!(meta(&identification.hits[0], "accessions_beta"), "-");
}

/// A cross-linking document with `library` peptides and one result holding
/// `items`, derived from the mzIdentML 1.3.0 schema rather than any fixture.
fn crosslinking_document(library: &str, items: &str) -> String {
    let results = format!(
        r#"   <SpectrumIdentificationResult spectraData_ref="SDAT" spectrumID="scan=7" id="SIR">
{items}
   </SpectrumIdentificationResult>"#
    );
    document(library, &results).replace(
        "<SearchType>",
        r#"<AdditionalSearchParams><cvParam accession="MS:1002494" cvRef="PSI-MS" name="crosslinking search"/><userParam name="cross_link:mass" type="xsd:double" value="138.0680796"/></AdditionalSearchParams><SearchType>"#,
    )
}

/// One `SpectrumIdentificationItem` of a cross-link group.
fn crosslink_item(id: &str, peptide: &str) -> String {
    format!(
        r#"    <SpectrumIdentificationItem passThreshold="true" rank="1" peptide_ref="{peptide}" chargeState="3" experimentalMassToCharge="500.5" id="{id}">
     <cvParam accession="MS:1002511" cvRef="PSI-MS" name="crosslink spectrum identification item" value="7"/>
     <cvParam accession="MS:1003024" cvRef="PSI-MS" name="OpenPepXL:score" value="0.5"/>
    </SpectrumIdentificationItem>"#
    )
}

#[test]
fn terminal_crosslink_positions_round_trip() {
    // No upstream fixture carries a terminal cross-link: location 0 is the
    // N-terminus and the peptide length plus one the C-terminus, which the
    // reader reports as the adjacent residue plus a terminal specificity.
    let library = r#"<DBSequence accession="P1" searchDatabase_ref="SDB" id="DBS"/>
<Peptide id="PEPA"><PeptideSequence>PEPTIDEK</PeptideSequence>
 <Modification location="0" monoisotopicMassDelta="138.0680796">
  <cvParam accession="XLMOD:02001" cvRef="XLMOD" name="DSS"/>
  <cvParam accession="MS:1002509" cvRef="PSI-MS" name="crosslink donor" value="7"/>
 </Modification>
</Peptide>
<Peptide id="PEPB"><PeptideSequence>KATSIDER</PeptideSequence>
 <Modification location="9" monoisotopicMassDelta="0">
  <cvParam accession="MS:1002510" cvRef="PSI-MS" name="crosslink acceptor" value="7"/>
 </Modification>
</Peptide>
<PeptideEvidence id="PEVA" peptide_ref="PEPA" dBSequence_ref="DBS" start="1" end="8"/>
<PeptideEvidence id="PEVB" peptide_ref="PEPB" dBSequence_ref="DBS" start="11" end="18"/>"#;
    let items = format!(
        "{}\n{}",
        crosslink_item("SIIA", "PEPA"),
        crosslink_item("SIIB", "PEPB")
    );
    let text = crosslinking_document(library, &items);
    let document = read(&text).expect("a terminal cross-link is read");
    let hit = &document.peptide_identifications[0].hits[0];
    assert_eq!(meta(hit, "xl_type"), "cross-link");
    assert_eq!(meta(hit, "xl_pos1"), "0");
    assert_eq!(meta(hit, "xl_term_spec_alpha"), "N_TERM");
    assert_eq!(meta(hit, "xl_pos2"), "7");
    assert_eq!(meta(hit, "xl_term_spec_beta"), "C_TERM");
    assert_eq!(meta(hit, "sequence_beta"), "KATSIDER");
    // Both specificities survive a store and a load.
    let reloaded = store_and_reload(&document, "terminal.mzid");
    let hit = &reloaded.peptide_identifications[0].hits[0];
    assert_eq!(meta(hit, "xl_pos1"), "0");
    assert_eq!(meta(hit, "xl_term_spec_alpha"), "N_TERM");
    assert_eq!(meta(hit, "xl_pos2"), "7");
    assert_eq!(meta(hit, "xl_term_spec_beta"), "C_TERM");
}

#[test]
fn a_loop_link_keeps_both_positions_on_one_chain() {
    // A loop-link is one peptide carrying both halves of the same link, which
    // the source recognises by the donor and the acceptor sharing a value.
    let library = r#"<DBSequence accession="P1" searchDatabase_ref="SDB" id="DBS"/>
<Peptide id="PEPA"><PeptideSequence>PEPTIDEK</PeptideSequence>
 <Modification location="2" residues="E" monoisotopicMassDelta="138.0680796">
  <cvParam accession="XLMOD:02001" cvRef="XLMOD" name="DSS"/>
  <cvParam accession="MS:1002509" cvRef="PSI-MS" name="crosslink donor" value="7"/>
 </Modification>
 <Modification location="8" residues="K" monoisotopicMassDelta="0">
  <cvParam accession="MS:1002510" cvRef="PSI-MS" name="crosslink acceptor" value="7"/>
 </Modification>
</Peptide>
<PeptideEvidence id="PEVA" peptide_ref="PEPA" dBSequence_ref="DBS" start="1" end="8"/>"#;
    let text = crosslinking_document(library, &crosslink_item("SIIA", "PEPA"));
    let document = read(&text).expect("a loop-link is read");
    let hit = &document.peptide_identifications[0].hits[0];
    assert_eq!(meta(hit, "xl_type"), "loop-link");
    assert_eq!(meta(hit, "xl_pos1"), "1");
    assert_eq!(meta(hit, "xl_pos2"), "7");
    // Both link positions become protein coordinates, counting from 1:
    // evidence start 0 plus the link position plus one.
    assert_eq!(meta(hit, "xl_pos1_protein"), "2");
    assert_eq!(meta(hit, "xl_pos2_protein"), "8");
    let reloaded = store_and_reload(&document, "looplink.mzid");
    let hit = &reloaded.peptide_identifications[0].hits[0];
    assert_eq!(meta(hit, "xl_type"), "loop-link");
    assert_eq!(meta(hit, "xl_pos1"), "1");
    assert_eq!(meta(hit, "xl_pos2"), "7");
}

#[test]
fn a_crosslink_group_without_an_experimental_mz_is_skipped() {
    // The source guards the light/heavy split with an emptiness check and then
    // indexes an empty vector when every item's m/z is absent; the group
    // produces no identification here.
    let library = r#"<Peptide id="PEPA"><PeptideSequence>PEPTIDEK</PeptideSequence>
 <Modification location="2" residues="E" monoisotopicMassDelta="138.0680796">
  <cvParam accession="XLMOD:02001" cvRef="XLMOD" name="DSS"/>
  <cvParam accession="MS:1002509" cvRef="PSI-MS" name="crosslink donor" value="7"/>
 </Modification>
</Peptide>"#;
    let item = crosslink_item("SIIA", "PEPA").replace(
        "experimentalMassToCharge=\"500.5\"",
        "experimentalMassToCharge=\"\"",
    );
    let text = crosslinking_document(library, &item);
    let document = read(&text).expect("the group is skipped, not refused");
    assert!(document.peptide_identifications.is_empty());
}

#[test]
fn a_crosslink_user_parameter_is_typed_by_either_attribute() {
    // The source's XL path types a userParam from unitName and its linear path
    // from type; this reader accepts either, so a value from an OpenMS-written
    // file and one from this writer both keep their number.
    let library = r#"<Peptide id="PEPA"><PeptideSequence>PEPTIDEK</PeptideSequence>
 <Modification location="2" residues="E" monoisotopicMassDelta="138.0680796">
  <cvParam accession="XLMOD:02001" cvRef="XLMOD" name="DSS"/>
  <cvParam accession="MS:1002509" cvRef="PSI-MS" name="crosslink donor" value="7"/>
 </Modification>
</Peptide>"#;
    let item = crosslink_item("SIIA", "PEPA").replace(
        "</SpectrumIdentificationItem>",
        r#"     <userParam name="unit_typed" unitName="xsd:double" value="1.5"/>
     <userParam name="type_typed" type="xsd:double" value="2.5"/>
    </SpectrumIdentificationItem>"#,
    );
    let text = crosslinking_document(library, &item);
    let document = read(&text).expect("a typed cross-linking userParam is read");
    let hit = &document.peptide_identifications[0].hits[0];
    for (key, expected) in [("unit_typed", 1.5), ("type_typed", 2.5)] {
        let value = hit
            .metadata
            .get(key)
            .unwrap_or_else(|| panic!("{key} is present"))
            .as_f64()
            .unwrap_or_else(|_| panic!("{key} is typed as a number"));
        assert!((value - expected).abs() < 1e-12, "{key} = {value}");
    }
}

// ---------------------------------------------------------------------------
// Upstream section: [EXTRA] XLMS data labeled cross-linker
// ---------------------------------------------------------------------------

#[test]
fn labelled_crosslinks_load_and_round_trip() {
    let original = load(XLMS_LABELLED);
    let hit = &original.peptide_identifications[1].hits[0];
    assert_eq!(meta(hit, "xl_pos1"), "3");
    assert_eq!(meta(hit, "xl_pos2"), "4");
    assert_eq!(meta(hit, "xl_term_spec_alpha"), "ANYWHERE");
    assert_eq!(meta(hit, "sequence_beta"), "SAVIKTSTR");
    assert_eq!(hit.sequence.to_string(), "FIVKASSGPR");
    assert_eq!(
        original.protein_identifications[0]
            .metadata
            .get("SpectrumIdentificationProtocol")
            .map(ToString::to_string)
            .unwrap_or_default(),
        "MS:1002494"
    );

    let reloaded = store_and_reload(&original, "xlms_labelled.mzid");
    let parameters = &reloaded.protein_identifications[0].search_parameters;
    assert_eq!(parameters.fragment_tolerance, Tolerance::Absolute(0.2));
    assert!(matches!(parameters.precursor_tolerance, Tolerance::Ppm(_)));
    for key in ["cross_link:residue1", "cross_link:residue2"] {
        assert_eq!(
            parameters
                .metadata
                .get(key)
                .map(ToString::to_string)
                .unwrap_or_default(),
            "[K, N-term]"
        );
    }
    let mass: f64 = parameters
        .metadata
        .get("cross_link:mass")
        .map(ToString::to_string)
        .unwrap_or_default()
        .parse()
        .expect("cross_link:mass is a number");
    assert!((mass - 138.0680796).abs() < 1e-6);
    let shift: f64 = parameters
        .metadata
        .get("cross_link:mass_isoshift")
        .map(ToString::to_string)
        .unwrap_or_default()
        .parse()
        .expect("cross_link:mass_isoshift is a number");
    assert!((shift - 12.075321).abs() < 1e-6);
    assert_eq!(
        parameters
            .metadata
            .get("extra_features")
            .map(ToString::to_string)
            .unwrap_or_default(),
        concat!(
            "precursor_mz_error_ppm,OpenPepXL:score,isotope_error,",
            "OpenPepXL:xquest_score,OpenPepXL:xcorr xlink,OpenPepXL:xcorr common,",
            "OpenPepXL:match-odds,OpenPepXL:intsum,OpenPepXL:wTIC,OpenPepXL:TIC,",
            "OpenPepXL:prescore,OpenPepXL:log_occupancy,OpenPepXL:log_occupancy_alpha,",
            "OpenPepXL:log_occupancy_beta,matched_xlink_alpha,matched_xlink_beta,",
            "matched_linear_alpha,matched_linear_beta,ppm_error_abs_sum_linear_alpha,",
            "ppm_error_abs_sum_linear_beta,ppm_error_abs_sum_xlinks_alpha,",
            "ppm_error_abs_sum_xlinks_beta,ppm_error_abs_sum_linear,",
            "ppm_error_abs_sum_xlinks,ppm_error_abs_sum_alpha,ppm_error_abs_sum_beta,",
            "ppm_error_abs_sum,precursor_total_intensity,precursor_target_intensity,",
            "precursor_signal_proportion,precursor_target_peak_count,",
            "precursor_residual_peak_count"
        )
    );

    // One identification per spectrum reference, in the order the merge step's
    // map produces.
    assert_eq!(reloaded.peptide_identifications.len(), 10);
    let identifications = &reloaded.peptide_identifications;
    assert_eq!(identifications[1].rt, identifications[2].rt);
    let rt = identifications[1].rt.expect("a light retention time");
    assert!((rt - 2132.4757).abs() < 1e-3, "{rt}");
    let mz = identifications[1].mz.expect("a light precursor m/z");
    assert!((mz - 721.0845).abs() < 1e-3, "{mz}");
    assert_eq!(
        identifications[1].spectrum_reference(),
        "spectrum=131,spectrum=113"
    );

    assert_eq!(identifications[0].hits.len(), 1);
    assert_eq!(identifications[1].hits.len(), 1);
    assert_eq!(identifications[3].hits.len(), 1);
    let hit = &identifications[1].hits[0];
    assert_eq!(meta(hit, "xl_type"), "cross-link");
    meta_number(hit, "spec_heavy_RT", 2125.5966796875);
    // The upstream literal is 725.109252929687841, the exact decimal
    // expansion of the same double.
    meta_number(hit, "spec_heavy_MZ", 725.109_252_929_687_8);
    assert!((hit.score - -0.190406834856118).abs() < 1e-12);
    assert_eq!(hit.sequence.to_string(), "FIVKASSGPR");
    assert_eq!(meta(hit, "sequence_beta"), "SAVIKTSTR");
    assert_eq!(meta(hit, "xl_pos1"), "3");
    assert_eq!(meta(hit, "xl_pos2"), "4");
    assert_eq!(meta(hit, "xl_term_spec_alpha"), "ANYWHERE");
    assert_eq!(meta(hit, "xl_term_spec_beta"), "ANYWHERE");
    meta_number(hit, "xl_mass", 138.0680796);
    assert_eq!(meta(hit, "xl_mod"), "DSS");
    assert_eq!(hit.peak_annotations[0].annotation, "[alpha|ci$b2]");
    assert_eq!(hit.peak_annotations[0].charge, 1);
    assert_eq!(hit.peak_annotations[1].annotation, "[beta|ci$y2]");
    assert_eq!(hit.peak_annotations[8].annotation, "[alpha|xi$b4]");
    let mono = &identifications[0].hits[0];
    assert_eq!(meta(mono, "xl_type"), "mono-link");
    assert_eq!(meta(mono, "xl_pos1"), "5");
    assert_eq!(meta(mono, "xl_pos2"), "-");
}

// ---------------------------------------------------------------------------
// Upstream section: [EXTRA] XLMS data unlabeled cross-linker
// ---------------------------------------------------------------------------

#[test]
fn unlabelled_crosslinks_round_trip() {
    let original = load(XLMS_UNLABELLED);
    let reloaded = store_and_reload(&original, "xlms_unlabelled.mzid");
    let parameters = &reloaded.protein_identifications[0].search_parameters;
    assert!(matches!(parameters.fragment_tolerance, Tolerance::Ppm(_)));
    assert!(matches!(parameters.precursor_tolerance, Tolerance::Ppm(_)));
    for key in ["cross_link:residue1", "cross_link:residue2"] {
        assert_eq!(
            parameters
                .metadata
                .get(key)
                .map(ToString::to_string)
                .unwrap_or_default(),
            "[K, N-term]"
        );
    }
    let mass: f64 = parameters
        .metadata
        .get("cross_link:mass")
        .map(ToString::to_string)
        .unwrap_or_default()
        .parse()
        .expect("cross_link:mass is a number");
    assert!((mass - 138.0680796).abs() < 1e-6);
    assert_eq!(
        original.protein_identifications[0]
            .metadata
            .get("SpectrumIdentificationProtocol")
            .map(ToString::to_string)
            .unwrap_or_default(),
        "MS:1002494"
    );

    let identifications = &reloaded.peptide_identifications;
    assert_eq!(identifications.len(), 3);
    let rt = identifications[0].rt.expect("a retention time");
    assert!((rt - 2175.3003).abs() < 1e-3, "{rt}");
    let mz = identifications[0].mz.expect("a precursor m/z");
    assert!((mz - 787.740356445313).abs() < 1e-9, "{mz}");
    assert_eq!(
        identifications[0].spectrum_reference(),
        "controllerType=0 controllerNumber=1 scan=2395"
    );

    for identification in identifications {
        assert_eq!(identification.hits.len(), 1);
    }
    assert_eq!(meta(&identifications[0].hits[0], "xl_type"), "mono-link");
    assert_eq!(meta(&identifications[1].hits[0], "xl_type"), "cross-link");
    assert_eq!(meta(&identifications[2].hits[0], "xl_type"), "mono-link");

    let mono = &identifications[0].hits[0];
    assert_eq!(meta(mono, "xl_pos1"), "5");
    assert_eq!(meta(mono, "xl_pos2"), "-");
    assert_eq!(meta(mono, "xl_term_spec_alpha"), "ANYWHERE");
    assert_eq!(meta(mono, "xl_term_spec_beta"), "ANYWHERE");

    let cross = &identifications[1].hits[0];
    assert_eq!(cross.sequence.to_string(), "KNVPIEFPVIDR");
    assert_eq!(meta(cross, "sequence_beta"), "LGCKALHVLFER");
    assert_eq!(meta(cross, "xl_pos1"), "0");
    assert_eq!(meta(cross, "xl_pos2"), "3");
    meta_number(cross, "xl_mass", 138.0680796);
    assert_eq!(meta(cross, "xl_mod"), "DSS");
    assert_eq!(meta(cross, "xl_term_spec_alpha"), "ANYWHERE");
    assert_eq!(meta(cross, "xl_term_spec_beta"), "ANYWHERE");
    assert_eq!(cross.peak_annotations.len(), 5);
    assert_eq!(cross.peak_annotations[0].annotation, "[alpha|ci$y5]");
    assert_eq!(cross.peak_annotations[0].charge, 1);
    assert_eq!(cross.peak_annotations[1].annotation, "[alpha|ci$y7]");
    assert_eq!(cross.peak_annotations[2].annotation, "[beta|ci$y7]");
    assert_eq!(cross.peak_annotations[3].annotation, "[alpha|ci$y8]");
    assert_eq!(cross.peak_annotations[2].charge, 1);
    assert_eq!(cross.peak_annotations[4].charge, 2);

    let loop_free = &identifications[2].hits[0];
    assert_eq!(
        loop_free.sequence.to_string(),
        "VEPSWLGPLFPDK(Xlink:DSS[156])TSNLR"
    );
    assert_eq!(meta(loop_free, "sequence_beta"), "-");
    assert_eq!(meta(loop_free, "xl_pos1"), "12");
    assert_eq!(meta(loop_free, "xl_pos2"), "-");
}

// ---------------------------------------------------------------------------
// Upstream section: [EXTRA] mzIdentML 1.3 crosslinking scores_and_thresholds
// ---------------------------------------------------------------------------

#[test]
fn crosslinking_v1_3_round_trips_at_schema_version_1_3_0() {
    let document = load(CROSSLINKING);
    assert!(!document.protein_identifications.is_empty());
    assert!(!document.peptide_identifications.is_empty());
    let mut total = 0usize;
    for identification in &document.peptide_identifications {
        // Every parsed identification carries a spectrum reference.
        assert!(!identification.spectrum_reference().is_empty());
        total += identification.hits.len();
        if let Some(hit) = identification.hits.first() {
            assert!(!hit.sequence.is_empty());
        }
    }
    assert!(total > 0, "no hits parsed");

    let directory = TempDir::new(false).expect("temporary directory");
    let path = directory.path().join("crosslinking.mzid");
    mzidentml::store(&path, &document).expect("store writes mzIdentML");
    let text = std::fs::read_to_string(&path).expect("the stored document is readable");
    assert!(
        text.lines().any(|line| line.contains("version=\"1.3.0\"")),
        "the default writer path emits version 1.3.0"
    );
    let reloaded = mzidentml::load(&path).expect("the stored document loads");
    assert!(!reloaded.peptide_identifications.is_empty());
    assert!(!reloaded.protein_identifications.is_empty());
}

// ---------------------------------------------------------------------------
// Upstream section: [EXTRA] mzIdentML 1.3 noncovalent association
// ---------------------------------------------------------------------------

#[test]
fn noncovalent_association_resolves_every_sequence() {
    let document = load(NONCOVALENT);
    assert!(!document.protein_identifications.is_empty());
    assert!(!document.peptide_identifications.is_empty());
    let mut total = 0usize;
    for identification in &document.peptide_identifications {
        total += identification.hits.len();
        if let Some(hit) = identification.hits.first() {
            // The modification and cvParam fallbacks resolved.
            assert!(!hit.sequence.is_empty());
        }
    }
    assert!(total > 0, "no hits parsed");
    // Noncovalently associated peptides carry no crosslink donor, so the
    // source reads each item of the result through the linear path. Both
    // candidates of the association survive: the source's merge step keeps
    // only the first hit of every identification, which would drop one.
    let identification = &document.peptide_identifications[0];
    assert_eq!(identification.hits.len(), 2);
    let hit = &identification.hits[0];
    assert_eq!(
        hit.sequence.to_string(),
        "AYALM(Oxidation)TDIHWDDC(Carbamidomethyl)FC(Carbamidomethyl)R"
    );
    assert_eq!(meta(hit, "xl_type"), "");
    assert_eq!(
        identification.hits[1].sequence.to_string(),
        "VHTEC(Carbamidomethyl)C(Carbamidomethyl)HGDLLEC(Carbamidomethyl)ADDR"
    );
}

// ---------------------------------------------------------------------------
// Upstream section: [EXTRA] mzIdentML 1.3 EDC crosslinking
// ---------------------------------------------------------------------------

#[test]
fn edc_crosslinking_reads_linked_and_linear_peptides() {
    // EDC files mix crosslinked and standalone peptides; both must parse.
    let document = load(EDC);
    assert!(!document.protein_identifications.is_empty());
    assert!(!document.peptide_identifications.is_empty());
    let total: usize = document
        .peptide_identifications
        .iter()
        .map(|identification| identification.hits.len())
        .sum();
    assert!(total > 0, "no hits parsed");
    // The fixture embeds its enzyme site pattern in a CDATA section, which is
    // ordinary character data and must not fail the parse.
    let text = std::fs::read_to_string(data(EDC)).expect("the fixture is readable");
    assert!(
        text.contains("<![CDATA["),
        "the fixture still carries CDATA"
    );
}

// ---------------------------------------------------------------------------
// Upstream section: [EXTRA] mzIdentML 1.3 multiple spectra per identification
// ---------------------------------------------------------------------------

#[test]
fn multiple_spectra_keep_distinct_references() {
    let document = load(MULTI_SPECTRA);
    assert!(!document.protein_identifications.is_empty());
    assert!(!document.peptide_identifications.is_empty());
    let references: std::collections::BTreeSet<String> = document
        .peptide_identifications
        .iter()
        .map(PeptideIdentification::spectrum_reference)
        .collect();
    assert!(references.len() > 1, "{references:?}");
}

// ---------------------------------------------------------------------------
// XML character data: entity references and CDATA
// ---------------------------------------------------------------------------

#[test]
fn entity_references_and_cdata_are_character_data() {
    // quick-xml reports every "&...;" as its own event. An arm that ignores
    // them deletes the character they stand for, which is how the sibling
    // Mascot XML reader turned a score of 1.5 into 15.
    let library = r#"<DBSequence accession="P1" searchDatabase_ref="SDB" id="DBS"><Seq>PEP&#84;IDEK<![CDATA[AC]]></Seq></DBSequence>
<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>
<PeptideEvidence id="PEV" peptide_ref="PEP" dBSequence_ref="DBS" start="1" end="8"/>"#;
    let results = r#"   <SpectrumIdentificationResult spectraData_ref="SDAT" spectrumID="scan=1" id="SIR">
    <SpectrumIdentificationItem passThreshold="true" rank="1" peptide_ref="PEP" chargeState="2" experimentalMassToCharge="500.5" id="SII">
     <PeptideEvidenceRef peptideEvidence_ref="PEV"/>
     <cvParam accession="MS:1001171" cvRef="PSI-MS" name="Mascot:score" value="1&#46;5"/>
    </SpectrumIdentificationItem>
   </SpectrumIdentificationResult>"#;
    let document = read(&document(library, results)).expect("entities are character data");
    let hit = &document.peptide_identifications[0].hits[0];
    // 1.5, not 15: the entity reference in the attribute is resolved, not cut.
    assert!((hit.score - 1.5).abs() < 1e-12, "{}", hit.score);
    // Element text keeps both the resolved entity and the CDATA section.
    assert_eq!(
        document.protein_identifications[0].hits[0].sequence,
        "PEPTIDEKAC"
    );
}

#[test]
fn an_external_entity_reference_is_refused() {
    let library = r#"<DBSequence accession="P1" searchDatabase_ref="SDB" id="DBS"><Seq>PEP&external;K</Seq></DBSequence>
<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let error = read(&document(library, &results)).unwrap_err();
    match &error {
        Error::Unsupported(message) => assert!(message.contains("external"), "{message}"),
        other => panic!("expected Unsupported, got {other:?}"),
    }
}

#[test]
fn unsupported_encoding_is_not_silently_interpreted_as_utf8() {
    let text = document("", "").replace("encoding=\"UTF-8\"", "encoding=\"ISO-8859-1\"");
    assert!(text.contains("ISO-8859-1"));
    assert!(matches!(read(&text), Err(Error::Unsupported(_))));
}

#[test]
fn a_doctype_declaration_is_still_refused() {
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let text = document(library, &results).replace(
        "<MzIdentML",
        "<!DOCTYPE MzIdentML [<!ENTITY x \"y\">]>\n<MzIdentML",
    );
    assert!(matches!(read(&text), Err(Error::Unsupported(_))));
}

// ---------------------------------------------------------------------------
// detectVersion
// ---------------------------------------------------------------------------

#[test]
fn detect_version_prefers_the_version_attribute() {
    assert_eq!(mzidentml::detect_version(data(WHOLE)).unwrap(), "1.1.0");
    assert_eq!(mzidentml::detect_version(data(MSGF)).unwrap(), "1.1.0");
    assert_eq!(
        mzidentml::detect_version(data(CROSSLINKING)).unwrap(),
        "1.3.0"
    );
}

#[test]
fn detect_version_falls_back_to_the_namespace_then_the_default() {
    let namespace_only =
        "<?xml version=\"1.0\"?>\n<MzIdentML xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.2\">";
    assert_eq!(
        mzidentml::detect_version_from_reader(namespace_only.as_bytes()).unwrap(),
        "1.2.0"
    );
    assert_eq!(
        mzidentml::detect_version_from_reader("<MzIdentML>".as_bytes()).unwrap(),
        SCHEMA_VERSION
    );
    // Only the first 15 lines are inspected, as TextFile(filename, true, 15).
    let buried = format!("{}<MzIdentML version=\"1.1.0\">", "\n".repeat(40));
    assert_eq!(
        mzidentml::detect_version_from_reader(buried.as_bytes()).unwrap(),
        SCHEMA_VERSION
    );
}

// ---------------------------------------------------------------------------
// Reference resolution: defined behaviour for every broken reference
// ---------------------------------------------------------------------------

#[test]
fn a_duplicate_element_id_is_refused() {
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>
<Peptide id="PEP"><PeptideSequence>PEPTIDER</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let error = read(&document(library, &results)).unwrap_err();
    match &error {
        Error::Parse { message, .. } => assert!(message.contains("duplicate Peptide"), "{message}"),
        other => panic!("expected Parse, got {other:?}"),
    }
}

#[test]
fn a_dangling_peptide_reference_is_refused() {
    let library = r#"<Peptide id="OTHER"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let error = read(&document(library, &results)).unwrap_err();
    match &error {
        Error::Parse { message, .. } => {
            assert!(message.contains("undeclared Peptide"), "{message}");
        }
        other => panic!("expected Parse, got {other:?}"),
    }
}

#[test]
fn a_peptide_evidence_pointing_at_an_undeclared_db_sequence_is_refused() {
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>
<PeptideEvidence id="PEV" peptide_ref="PEP" dBSequence_ref="MISSING"/>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let error = read(&document(library, &results)).unwrap_err();
    match &error {
        Error::Parse { message, .. } => {
            assert!(message.contains("undeclared DBSequence"), "{message}");
        }
        other => panic!("expected Parse, got {other:?}"),
    }
}

#[test]
fn an_unreferenced_spectrum_identification_list_is_refused() {
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let text = document(library, &results).replace(
        r#"<SpectrumIdentificationList id="SIL">"#,
        r#"<SpectrumIdentificationList id="OTHER">"#,
    );
    let error = read(&text).unwrap_err();
    match &error {
        Error::Parse { message, .. } => {
            assert!(message.contains("not referenced"), "{message}");
        }
        other => panic!("expected Parse, got {other:?}"),
    }
}

#[test]
fn two_runs_sharing_one_list_are_refused() {
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let text = document(library, &results).replace(
        "</AnalysisCollection>",
        r#" <SpectrumIdentification id="SI2" spectrumIdentificationProtocol_ref="SIP" spectrumIdentificationList_ref="SIL">
  <InputSpectra spectraData_ref="SDAT"/>
  <SearchDatabaseRef searchDatabase_ref="SDB"/>
 </SpectrumIdentification>
</AnalysisCollection>"#,
    );
    let error = read(&text).unwrap_err();
    match &error {
        Error::Parse { message, .. } => {
            assert!(
                message.contains("reference SpectrumIdentificationList"),
                "{message}"
            );
        }
        other => panic!("expected Parse, got {other:?}"),
    }
}

#[test]
fn dangling_metadata_references_stay_tolerated() {
    // MzIdentMLFile_msgf_mini references SearchProtocol_99 from its second
    // SpectrumIdentification and SearchDB_99 from DBSeq99, and declares
    // neither, so neither may be an error. The unresolved protocol leaves the
    // run without a search engine or parameters; the unresolved
    // searchDatabase_ref on a DBSequence is never consulted at all.
    let document = load(MSGF);
    let second = &document.protein_identifications[1];
    assert_eq!(second.search_engine, "");
    assert_eq!(second.search_engine_version, "");
    assert_eq!(second.search_parameters.missed_cleavages, 0);
    assert!(second.search_parameters.fixed_modifications.is_empty());
    // The database still comes from this run's own SearchDatabaseRef.
    assert_eq!(second.search_parameters.database, "database.fasta");
    assert_eq!(second.hits.len(), 1);
    assert_eq!(second.hits[0].accession, "sp|ABC|xyz");
}

#[test]
fn a_required_element_that_is_absent_is_reported_as_missing() {
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let complete = document(library, &results);
    for (removed, marker) in [
        (
            "<SpectraData location=\"run.mzML\" id=\"SDAT\"/>",
            "SpectraData",
        ),
        (
            "<SpectrumIdentificationProtocol id=\"SIP\" analysisSoftware_ref=\"SOF\">",
            "SpectrumIdentificationProtocol",
        ),
    ] {
        let text = complete.replace(removed, "");
        let error = read(&text).unwrap_err();
        assert!(
            matches!(&error, Error::MissingInformation(message) if message.contains(marker))
                || matches!(&error, Error::Parse { .. }),
            "{marker}: {error:?}"
        );
    }
}

// ---------------------------------------------------------------------------
// Native parser boundaries
// ---------------------------------------------------------------------------

#[test]
fn non_ascii_text_survives_the_round_trip() {
    let library = r#"<DBSequence accession="日本語|P1" searchDatabase_ref="SDB" id="DBS"><Seq>PEPTIDEK</Seq></DBSequence>
<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>
<PeptideEvidence id="PEV" peptide_ref="PEP" dBSequence_ref="DBS" start="1" end="8"/>"#;
    let results = result_with(
        "chargeState=\"2\" experimentalMassToCharge=\"500.5\"",
        "     <PeptideEvidenceRef peptideEvidence_ref=\"PEV\"/>",
    );
    let document = read(&document(library, &results)).expect("non-ASCII document loads");
    assert_eq!(
        document.protein_identifications[0].hits[0].accession,
        "日本語|P1"
    );
    let mut text = Vec::new();
    mzidentml::write(&mut text, &document).expect("non-ASCII document writes");
    let reloaded = mzidentml::read(text.as_slice()).expect("non-ASCII output reloads");
    assert_eq!(
        reloaded.protein_identifications[0].hits[0].accession,
        "日本語|P1"
    );
}

#[test]
fn a_prefixed_namespace_is_resolved() {
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let text = document(library, &results)
        .replace("<MzIdentML xmlns=", "<mzid:MzIdentML xmlns:mzid=")
        .replace("</MzIdentML>", "</mzid:MzIdentML>");
    // Only the root is prefixed here, which leaves the children unbound; a
    // fully prefixed document is what the source cannot read at all.
    assert!(read(&text).is_err());
    let mut all = document(library, &results);
    for tag in [
        "MzIdentML",
        "SequenceCollection",
        "Peptide",
        "PeptideSequence",
        "AnalysisCollection",
        "SpectrumIdentification",
        "InputSpectra",
        "SearchDatabaseRef",
        "AnalysisProtocolCollection",
        "SpectrumIdentificationProtocol",
        "SearchType",
        "Threshold",
        "cvParam",
        "DataCollection",
        "Inputs",
        "SearchDatabase",
        "DatabaseName",
        "userParam",
        "SpectraData",
        "AnalysisData",
        "SpectrumIdentificationList",
        "SpectrumIdentificationResult",
        "SpectrumIdentificationItem",
    ] {
        all = all
            .replace(&format!("<{tag}"), &format!("<mzid:{tag}"))
            .replace(&format!("</{tag}>"), &format!("</mzid:{tag}>"));
    }
    let all = all.replace("<mzid:MzIdentML xmlns=", "<mzid:MzIdentML xmlns:mzid=");
    let document = mzidentml::read(all.as_bytes()).expect("a prefixed document loads");
    assert_eq!(document.peptide_identifications.len(), 1);
}

#[test]
fn a_foreign_root_namespace_is_refused() {
    let text =
        r#"<?xml version="1.0"?><MzIdentML xmlns="http://example.invalid/mzid" version="1.3.0"/>"#;
    assert!(matches!(read(text), Err(Error::Unsupported(_))));
}

#[test]
fn resource_ceilings_refuse_before_allocating() {
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let text = document(library, &results);
    for options in [
        ReadOptions {
            max_xml_bytes: 32,
            ..Default::default()
        },
        ReadOptions {
            max_elements: 3,
            ..Default::default()
        },
        ReadOptions {
            max_depth: 2,
            ..Default::default()
        },
        ReadOptions {
            max_work: 16,
            ..Default::default()
        },
        ReadOptions {
            max_payload_bytes: 512,
            ..Default::default()
        },
    ] {
        assert!(
            mzidentml::read_with_options(text.as_bytes(), &options).is_err(),
            "{options:?}"
        );
    }
    // Zero limits are a caller error, not a parse failure.
    let invalid = ReadOptions {
        max_depth: 0,
        ..Default::default()
    };
    assert!(matches!(
        mzidentml::read_with_options(text.as_bytes(), &invalid),
        Err(Error::InvalidValue(_))
    ));
    // The cross-linking path is metered by the same budget: one item may belong
    // to several groups, so the work of re-scanning its parameters is counted
    // rather than the item alone.
    let options = ReadOptions {
        max_work: 4096,
        ..Default::default()
    };
    assert!(
        mzidentml::load_with_options(data(XLMS_UNLABELLED), &options).is_err(),
        "the cross-linking read draws from the shared work budget"
    );
}

#[test]
fn an_empty_peptide_sequence_is_refused() {
    // The source dereferences the missing text node instead.
    let library = r#"<Peptide id="PEP"><PeptideSequence></PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let error = read(&document(library, &results)).unwrap_err();
    match &error {
        Error::Parse { message, .. } => assert!(message.contains("must not be empty"), "{message}"),
        other => panic!("expected Parse, got {other:?}"),
    }
}

#[test]
fn a_substitution_location_outside_the_peptide_is_refused() {
    // The source writes through as[location - 1] without a bounds check.
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence>
<SubstitutionModification location="99" originalResidue="K" replacementResidue="R"/></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    assert!(matches!(
        read(&document(library, &results)),
        Err(Error::InvalidRange(_))
    ));
}

#[test]
fn a_substitution_without_a_location_replaces_every_occurrence() {
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence>
<SubstitutionModification originalResidue="E" replacementResidue="A"/></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let document = read(&document(library, &results)).expect("substitution applies");
    assert_eq!(
        document.peptide_identifications[0].hits[0]
            .sequence
            .as_str(),
        "PAPTIDAK"
    );
}

#[test]
fn rank_zero_is_read_as_the_first_rank() {
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "")
        .replace("rank=\"1\"", "rank=\"0\"");
    let loaded = read(&document(library, &results)).expect("PMF rank 0 loads");
    assert_eq!(loaded.peptide_identifications[0].hits[0].rank, 0);
    // A negative rank cannot underflow the unsigned rank member here.
    let negative = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "")
        .replace("rank=\"1\"", "rank=\"-3\"");
    assert!(read(&document(library, &negative)).is_err());
}

#[test]
fn the_score_scan_follows_lexicographic_accession_order() {
    // MS-GF+ writes four scores; MS:1002049 (MS-GF:RawScore) wins because it
    // sorts first among the PSM-level statistics, not because of its position.
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with(
        "chargeState=\"2\" experimentalMassToCharge=\"500.5\"",
        r#"     <cvParam accession="MS:1002053" cvRef="PSI-MS" name="MS-GF:EValue" value="1.5"/>
     <cvParam accession="MS:1002049" cvRef="PSI-MS" name="MS-GF:RawScore" value="195"/>"#,
    )
    .replace(
        r#"     <cvParam accession="MS:1001171" cvRef="PSI-MS" name="Mascot:score" value="42.5"/>
"#,
        "",
    );
    let document = read(&document(library, &results)).expect("scored item loads");
    let identification = &document.peptide_identifications[0];
    assert_eq!(identification.score_type, "MS-GF:RawScore");
    assert_eq!(identification.hits[0].score, 195.0);
}

#[test]
fn a_c_terminal_modification_is_written_at_length_plus_one() {
    // The source writes `length`, which its own reader then applies to the
    // last residue; this asserts the schema position and a stable round trip.
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDER</PeptideSequence>
<Modification location="9"><cvParam accession="UNIMOD:2" cvRef="UNIMOD" name="Amidated"/></Modification></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let document = read(&document(library, &results)).expect("C-terminal modification loads");
    let sequence = &document.peptide_identifications[0].hits[0].sequence;
    assert!(sequence.c_terminal_modification().is_some());
    let mut text = Vec::new();
    mzidentml::write(&mut text, &document).expect("writes the C-terminal modification");
    let text = String::from_utf8(text).unwrap();
    assert!(text.contains(r#"<Modification location="9""#), "{text}");
    let reloaded = mzidentml::read(text.as_bytes()).expect("reloads");
    assert_eq!(
        reloaded.peptide_identifications[0].hits[0].sequence,
        *sequence
    );
}

#[test]
fn a_protein_detection_list_appends_to_the_last_run() {
    let library = r#"<DBSequence accession="P1" searchDatabase_ref="SDB" id="DBS"><Seq>PEPTIDEK</Seq></DBSequence>
<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let text = document(library, &results).replace(
        " </AnalysisData>",
        r#"  <ProteinDetectionList id="PDL">
   <ProteinAmbiguityGroup id="PAG">
    <ProteinDetectionHypothesis id="PDH" dBSequence_ref="DBS" passThreshold="true"/>
   </ProteinAmbiguityGroup>
  </ProteinDetectionList>
 </AnalysisData>"#,
    );
    let document = mzidentml::read(text.as_bytes()).expect("protein detection list loads");
    let run = &document.protein_identifications[0];
    assert_eq!(run.hits.len(), 1);
    assert_eq!(run.hits[0].accession, "P1");
    assert_eq!(run.hits[0].sequence, "PEPTIDEK");
}

#[test]
fn peak_annotations_round_trip_through_fragmentation() {
    use openms::identification::PeakAnnotation;
    let mut document = load(WHOLE);
    document.peptide_identifications[0].hits[0].peak_annotations = vec![
        PeakAnnotation {
            mz: 363.908,
            intensity: 208.52,
            charge: 1,
            annotation: "[b2]".into(),
        },
        PeakAnnotation {
            mz: 511.557,
            intensity: 2034.9,
            charge: 1,
            annotation: "[y3-H2O]".into(),
        },
        PeakAnnotation {
            mz: 754.418,
            intensity: 1098.44,
            charge: 2,
            annotation: "not an mzIdentML fragment".into(),
        },
    ];
    let mut text = Vec::new();
    mzidentml::write(&mut text, &document).expect("writes fragment annotations");
    let text = String::from_utf8(text).unwrap();
    assert!(text.contains("<Fragmentation>"), "{text}");
    assert!(text.contains("frag: b ion"), "{text}");
    assert!(text.contains("frag: y ion - H2O"), "{text}");
    // The unmatched annotation is dropped rather than failing the write, as in
    // the source; reading ignores Fragmentation entirely, as it also does.
    let reloaded = mzidentml::read(text.as_bytes()).expect("reloads");
    assert!(
        reloaded.peptide_identifications[0].hits[0]
            .peak_annotations
            .is_empty()
    );
}

#[test]
fn a_peptide_identification_without_hits_cannot_be_written() {
    let mut document = load(WHOLE);
    document.peptide_identifications[0].hits.clear();
    let mut text = Vec::new();
    assert!(matches!(
        mzidentml::write(&mut text, &document),
        Err(Error::MissingInformation(_))
    ));
    assert!(text.is_empty());
}

#[test]
fn an_unlinked_peptide_identification_cannot_be_written() {
    let mut document = load(WHOLE);
    document.peptide_identifications[0].identifier = "nowhere".into();
    let mut text = Vec::new();
    let error = mzidentml::write(&mut text, &document).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
}

#[test]
fn an_evidence_without_a_declared_protein_cannot_be_written() {
    let mut document = load(WHOLE);
    document.peptide_identifications[0].hits[0].evidences[0].protein_accession =
        "UNDECLARED".into();
    let mut text = Vec::new();
    let error = mzidentml::write(&mut text, &document).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
}

#[test]
fn write_ceilings_refuse_before_any_output_is_produced() {
    let document = load(WHOLE);
    let mut text = Vec::new();
    let options = WriteOptions {
        max_output_bytes: 128,
        ..Default::default()
    };
    assert!(mzidentml::write_with_options(&mut text, &document, &options).is_err());
    assert!(text.is_empty());
    let options = WriteOptions {
        max_records: 4,
        ..Default::default()
    };
    assert!(mzidentml::write_with_options(&mut text, &document, &options).is_err());
    assert!(text.is_empty());
    // The record ceiling stops the plan as it grows, not only the emission, so
    // a document far above it is refused without building the whole plan
    // first. A cross-linking document plans two peptides per hit and takes the
    // same ceiling.
    for name in [WHOLE, XLMS_UNLABELLED] {
        let document = load(name);
        let options = WriteOptions {
            max_records: 1,
            ..Default::default()
        };
        assert!(
            mzidentml::write_with_options(&mut text, &document, &options).is_err(),
            "{name}"
        );
        assert!(text.is_empty(), "{name}");
    }
}

#[test]
fn a_caller_registry_resolves_the_modification_names() {
    // The registry is a parameter, not a global, so a caller-owned copy works.
    let registry = ModificationsDB::global();
    let document = mzidentml::load_with_registry(data(MSGF), &ReadOptions::default(), registry)
        .expect("load with an explicit registry");
    assert_eq!(document.peptide_identifications.len(), 5);
}
