// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Coverage for `FORMAT/MzIdentMLFile.h` and
//! `FORMAT/HANDLERS/MzIdentMLHandler.h`.
//!
//! Every literal taken from `MzIdentMLFile_test.cpp` or from one of the four
//! unmodified upstream fixtures is transcribed source review (tier 3): the
//! counts, accessions, score types, scores, spectrum references, the
//! `database.fasta`/`MSDB` search databases, `missed_cleavages` 1000, the
//! `Carbamidomethyl (C)` / `Xlink:DTSSP[88] (Protein N-term)` fixed
//! modifications, the `Acetyl (N-term)` variable modification, the 0.5
//! significance threshold, and the modification inference cases of issue #5443.
//! No C++ was built or run, so nothing here is a tier 1 differential.
//!
//! The synthetic documents - duplicate ids, dangling references, an
//! unreferenced `SpectrumIdentificationList`, a prefixed namespace, non-ASCII
//! text, the resource ceilings, the C-terminal modification location and the
//! `ProteinDetectionList` - are independently derived from the mzIdentML 1.3.0
//! schema (tier 4), because no upstream fixture reaches them.
//!
//! `docs/MZIDENTML_SUPPORT.md` records which upstream sections are ported here
//! and which are not: the six cross-linking sections are not, because the read
//! path they exercise needs `ANALYSIS/XLMS/OPXLHelper.h`, which is unported.
//! What is pinned instead is that such a document is refused explicitly.

#![cfg(feature = "idxml")]

use openms::Error;
use openms::chemistry::ModificationsDB;
use openms::comparison::Tolerance;
use openms::format::mzidentml::{
    self, MzIdentMLDocument, ReadOptions, SCHEMA_VERSION, WriteOptions,
};
use openms::identification::{EnzymeTermSpecificity, FlankingResidue};
use openms::system::file::TempDir;
use std::path::PathBuf;

const WHOLE: &str = "mzidentml_whole.mzid";
const MSGF: &str = "mzidentml_msgf_mini.mzid";
const MISSING_LOCATION: &str = "mzidentml_missing_mod_location.mzid";
const THREE_RUNS: &str = "mzidentml_3runs.mzid";
const CROSSLINKING: &str = "mzidentml_crosslinking_v1_3.mzid";

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
// Cross-linking: the six upstream XL sections are not ported
// ---------------------------------------------------------------------------

#[test]
fn crosslinking_documents_are_refused_explicitly() {
    let error = mzidentml::load(data(CROSSLINKING)).unwrap_err();
    match &error {
        Error::Unsupported(message) => {
            assert!(message.contains("MS:1002494"), "{message}");
            assert!(message.contains("OpenPepXL"), "{message}");
        }
        other => panic!("expected Unsupported, got {other:?}"),
    }
}

#[test]
fn the_crosslinking_marker_is_detected_wherever_it_appears() {
    let library = r#"<Peptide id="PEP"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>"#;
    let results = result_with("chargeState=\"2\" experimentalMassToCharge=\"500.5\"", "");
    let text = document(library, &results).replace(
        "<SearchType>",
        r#"<AdditionalSearchParams><cvParam accession="MS:1002494" cvRef="PSI-MS" name="crosslinking search"/></AdditionalSearchParams><SearchType>"#,
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
}

#[test]
fn a_caller_registry_resolves_the_modification_names() {
    // The registry is a parameter, not a global, so a caller-owned copy works.
    let registry = ModificationsDB::global();
    let document = mzidentml::load_with_registry(data(MSGF), &ReadOptions::default(), registry)
        .expect("load with an explicit registry");
    assert_eq!(document.peptide_identifications.len(), 5);
}
