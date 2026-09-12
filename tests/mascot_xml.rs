// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//
//! Tests for `src/format/mascot_xml.rs`, the port of `FORMAT/MascotXMLFile.h`
//! and the `MascotXMLHandler` it drives.
//!
//! Every `START_SECTION` of the upstream `MascotXMLFile_test.cpp` is
//! reproduced. The upstream section for the four-argument `load` ends by
//! writing the parsed result to idXML and fuzzy-comparing it with the retained
//! `MascotXMLFile_test_out_3.idXML`; here that retained C++ output is read back
//! with this crate's own idXML reader and compared record by record, which is
//! tier-1 differential evidence. See `docs/MASCOT_XML_SUPPORT.md` and
//! `tests/data/mascot_xml_provenance.json`.
#![cfg(feature = "idxml")]

use openms::chemistry::AASequence;
use openms::comparison::Tolerance;
use openms::format::mascot_xml::{
    self as mascot, MascotXmlFile, SpectrumTitleLookup, TitleReferenceFormat,
};
use openms::identification::{PeakMassType, PeptideIdentification};
use openms::{Error, MSExperiment, MSSpectrum, Precursor};
use std::collections::BTreeMap;

const TEST_1: &str = "tests/data/MascotXMLFile_test_1.mascotXML";
const TEST_2: &str = "tests/data/MascotXMLFile_test_2.mascotXML";
const TEST_3: &str = "tests/data/MascotXMLFile_test_3.mascotXML";
const OUT_3: &str = "tests/data/MascotXMLFile_test_out_3.idXML";

/// The absolute tolerance the upstream `FuzzyStringComparator` is configured
/// with for this comparison (`setAcceptableAbsolute(0.0001)`).
const FUZZY_ABSOLUTE: f64 = 0.0001;

fn accessions(identification: &PeptideIdentification, hit: usize) -> Vec<String> {
    let mut values: Vec<String> = identification.hits[hit]
        .evidences
        .iter()
        .map(|e| e.protein_accession.clone())
        .collect();
    values.sort();
    values.dedup();
    values
}

// ---------------------------------------------------------------------------
// START_SECTION((MascotXMLFile()))
// ---------------------------------------------------------------------------

#[test]
fn default_construction() {
    // The upstream section only asserts a non-null pointer; the Rust analogue
    // is that the reader is a stateless value usable straight away.
    let file = MascotXmlFile::new();
    let lookup = SpectrumTitleLookup::new();
    assert!(lookup.is_empty());
    assert_eq!(lookup.len(), 0);
    assert!(lookup.reference_formats().is_empty());
    let result = file.load(TEST_1, &lookup).unwrap();
    assert_eq!(result.peptide_identifications.len(), 3);
    assert_eq!(file, MascotXmlFile);
}

// ---------------------------------------------------------------------------
// START_SECTION((static void initializeLookup(...)))
// ---------------------------------------------------------------------------

#[test]
fn initialize_lookup_reads_the_spectra_and_registers_the_default_formats() {
    let experiment = MSExperiment {
        spectra: vec![MSSpectrum::default()],
        ..Default::default()
    };
    let (lookup, warnings) = MascotXmlFile::initialize_lookup(&experiment, None).unwrap();
    assert!(!lookup.is_empty());
    assert_eq!(lookup.len(), 1);
    // The default spectrum has an empty native ID, so no scan number could be
    // extracted and the source logs exactly that warning.
    assert_eq!(warnings.len(), 1);
    assert!(
        warnings[0].contains("Could not extract scan number"),
        "{warnings:?}"
    );
    assert_eq!(
        lookup.reference_formats(),
        [
            TitleReferenceFormat::ScanNumber,
            TitleReferenceFormat::DtaFileName,
            TitleReferenceFormat::MzThenRt,
        ]
    );
    // Without raw data only the m/z-underscore-RT format is registered, because
    // it is the only one that needs no spectrum to resolve.
    let (lookup, warnings) =
        MascotXmlFile::initialize_lookup(&MSExperiment::default(), None).unwrap();
    assert!(lookup.is_empty());
    assert!(warnings.is_empty());
    assert_eq!(lookup.reference_formats(), [TitleReferenceFormat::MzThenRt]);
    // A caller-supplied regular expression cannot be compiled here.
    assert!(matches!(
        MascotXmlFile::initialize_lookup(&experiment, Some("scan=(?<SCAN>\\d+)")),
        Err(Error::Unsupported(_))
    ));
    // An empty or whitespace-only override is the source default.
    assert!(MascotXmlFile::initialize_lookup(&experiment, Some("")).is_ok());
}

// ---------------------------------------------------------------------------
// START_SECTION((void load(const std::string&, ProteinIdentification&,
//                          PeptideIdentificationList&, SpectrumMetaDataLookup&)))
// ---------------------------------------------------------------------------

#[test]
fn load_mascot_xml_1_0() {
    let lookup = SpectrumTitleLookup::new();
    let result = mascot::load(TEST_1, &lookup).unwrap();
    let search = &result.protein_identification.search_parameters;
    assert_eq!(search.missed_cleavages, 1);
    assert_eq!(search.taxonomy, ". . Eukaryota (eucaryotes)");
    assert_eq!(search.mass_type, PeakMassType::Average);
    assert_eq!(search.database, "MSDB_chordata");
    assert_eq!(search.database_version, "MSDB_chordata_20070910.fasta");
    assert_eq!(search.fragment_tolerance, Tolerance::Absolute(0.2));
    assert_eq!(search.precursor_tolerance, Tolerance::Absolute(1.4));
    assert_eq!(search.charges, "1+, 2+ and 3+");
    assert_eq!(
        search.fixed_modifications,
        [
            "Carboxymethyl (C)",
            "Deamidated (N)",
            "Deamidated (Q)",
            "Guanidinyl (K)"
        ]
    );
    assert_eq!(
        search.variable_modifications,
        ["Acetyl (Protein N-term)", "Biotin (K)", "Carbamyl (K)"]
    );

    let ids = &result.peptide_identifications;
    assert_eq!(ids.len(), 3);
    assert!((ids[0].mz.unwrap() - 789.83).abs() < FUZZY_ABSOLUTE);
    assert!((ids[1].mz.unwrap() - 135.29).abs() < FUZZY_ABSOLUTE);
    assert!((ids[2].mz.unwrap() - 982.58).abs() < FUZZY_ABSOLUTE);

    let protein = &result.protein_identification;
    assert_eq!(protein.hits.len(), 2);
    assert_eq!(protein.hits[0].accession, "AAN17824");
    assert_eq!(protein.hits[1].accession, "GN1736");
    assert_eq!(protein.hits[0].score, 619.0);
    assert_eq!(protein.hits[1].score, 293.0);
    assert_eq!(protein.score_type, "Mascot");
    assert_eq!(protein.search_engine, "Mascot");
    assert_eq!(protein.date_time.as_deref(), Some("2006-03-09 11:31:52"));

    assert!((ids[0].significance_threshold - 31.8621).abs() < 1e-9);
    assert_eq!(ids[0].hits.len(), 2);
    assert_eq!(accessions(&ids[0], 0), ["AAN17824", "GN1736"]);
    assert_eq!(accessions(&ids[0], 1), ["AAN17824"]);
    assert_eq!(accessions(&ids[1], 0), ["GN1736"]);
    assert_eq!(ids[1].hits.len(), 1);
    assert!((ids[0].hits[0].score - 33.85).abs() < 1e-9);
    assert!((ids[0].hits[1].score - 33.12).abs() < 1e-9);
    assert!((ids[1].hits[0].score - 43.9).abs() < 1e-9);
    assert_eq!(ids[0].score_type, "Mascot");
    assert_eq!(ids[1].score_type, "Mascot");

    // Fixed modifications are applied to every matching residue. The upstream
    // literals are PSI-MOD accessions -- LHASGITVTEIPVTATN(MOD:00565)FK(MOD:00445),
    // MRSLGYVAVISAVATDTDK(MOD:00445) and HSK(MOD:00445)LSAK(MOD:00445) -- and
    // this crate's modification table carries UniMod accessions instead, so the
    // same three sequences are asserted residue by residue and then through the
    // UniMod rendering. MOD:00565 is deamidation and MOD:00445 guanidination,
    // which UniMod:7 and UniMod:52 are the same two modifications.
    let sequence = &ids[0].hits[0].sequence;
    assert_eq!(sequence.as_str(), "LHASGITVTEIPVTATNFK");
    assert_eq!(
        sequence.residue_modification(16).unwrap().map(|m| m.name()),
        Some("Deamidated")
    );
    assert_eq!(
        sequence.residue_modification(18).unwrap().map(|m| m.name()),
        Some("Guanidinyl")
    );
    assert!(sequence.n_terminal_modification().is_none());
    assert!(sequence.c_terminal_modification().is_none());
    assert_eq!(
        sequence.to_accession_string(),
        "LHASGITVTEIPVTATN(UniMod:7)FK(UniMod:52)"
    );
    let sequence = &ids[0].hits[1].sequence;
    assert_eq!(sequence.as_str(), "MRSLGYVAVISAVATDTDK");
    assert_eq!(
        sequence.residue_modification(18).unwrap().map(|m| m.name()),
        Some("Guanidinyl")
    );
    assert_eq!(
        sequence.to_accession_string(),
        "MRSLGYVAVISAVATDTDK(UniMod:52)"
    );
    let sequence = &ids[1].hits[0].sequence;
    assert_eq!(sequence.as_str(), "HSKLSAK");
    assert_eq!(
        sequence.residue_modification(2).unwrap().map(|m| m.name()),
        Some("Guanidinyl")
    );
    assert_eq!(
        sequence.residue_modification(6).unwrap().map(|m| m.name()),
        Some("Guanidinyl")
    );
    assert_eq!(
        sequence.to_accession_string(),
        "HSK(UniMod:52)LSAK(UniMod:52)"
    );

    let identifier = &protein.identifier;
    assert_eq!(identifier, "Mascot_2006-03-09 11:31:52");
    for identification in ids {
        assert_eq!(&identification.identifier, identifier);
    }
    // An empty lookup registers no reference format, so no retention time is
    // recovered and the handler reports that once per document.
    assert!(ids.iter().all(|id| id.rt.is_none()));
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.contains("peptide identifications have no retention time value")),
        "{:?}",
        result.warnings
    );
}

#[test]
fn load_mascot_xml_2_1_as_written_by_mascot_server_2_3() {
    let lookup = SpectrumTitleLookup::new();
    let result = mascot::load(TEST_2, &lookup).unwrap();
    let search = &result.protein_identification.search_parameters;
    assert_eq!(search.missed_cleavages, 7);
    assert_eq!(search.taxonomy, "All entries");
    assert_eq!(search.mass_type, PeakMassType::Monoisotopic);
    assert_eq!(search.database, "IPI_human");
    assert_eq!(search.database_version, "ipi.HUMAN.v3.61.fasta");
    assert_eq!(search.fragment_tolerance, Tolerance::Absolute(0.3));
    assert_eq!(search.precursor_tolerance, Tolerance::Absolute(3.0));
    assert_eq!(search.charges, "");
    assert_eq!(search.fixed_modifications, ["Carbamidomethyl (C)"]);
    assert_eq!(
        search.variable_modifications,
        ["Oxidation (M)", "Acetyl (N-term)", "Phospho (Y)"]
    );

    // Not necessarily equal to NumQueries: a query whose first peptide starts
    // at rank 10 contributes no hit, and those empty entries are dropped.
    let ids = &result.peptide_identifications;
    assert_eq!(ids.len(), 1112);
    assert!((ids[0].mz.unwrap() - 304.6967).abs() < FUZZY_ABSOLUTE);
    assert!((ids[1].mz.unwrap() - 314.1815).abs() < FUZZY_ABSOLUTE);
    assert!((ids[1111].mz.unwrap() - 583.7948).abs() < FUZZY_ABSOLUTE);

    let protein = &result.protein_identification;
    assert_eq!(protein.hits.len(), 66);
    assert_eq!(protein.hits[0].accession, "IPI00745872");
    assert_eq!(protein.hits[1].accession, "IPI00908876");
    assert_eq!(protein.hits[0].score, 122.0);
    assert_eq!(protein.hits[1].score, 122.0);
    assert_eq!(protein.score_type, "Mascot");
    assert_eq!(protein.date_time.as_deref(), Some("2011-06-24 19:34:54"));

    assert!((ids[0].significance_threshold - 5.0).abs() < 1e-9);
    assert_eq!(ids[0].hits.len(), 1);
    // The first query's only hit carries no protein evidence, because its
    // accession is attached by the enclosing <protein> and this hit came from
    // the <unassigned> section.
    assert!(ids[0].hits[0].evidences.is_empty());
    // Query 35 is shared by five proteins, which the dedup path collects as
    // five peptide evidences on one hit.
    assert_eq!(
        accessions(&ids[34], 0),
        [
            "IPI00022434",
            "IPI00384697",
            "IPI00745872",
            "IPI00878517",
            "IPI00908876"
        ]
    );

    assert!((ids[0].hits[0].score - 5.34).abs() < 1e-9);
    assert!((ids[49].hits[0].score - 14.83).abs() < 1e-9);
    assert!((ids[49].hits[1].score - 17.5).abs() < 1e-9);
    assert_eq!(ids[0].score_type, "Mascot");
    assert_eq!(ids[1].score_type, "Mascot");

    // VVFIK carries no modification.
    assert_eq!(ids[0].hits[0].sequence, AASequence::parse("VVFIK").unwrap());
    assert_eq!(
        ids[49].hits[0].sequence,
        AASequence::parse("LASYLDK").unwrap()
    );
    // (Acetyl)AAFESDK
    let sequence = &ids[49].hits[1].sequence;
    assert_eq!(sequence.as_str(), "AAFESDK");
    assert_eq!(
        sequence.n_terminal_modification().map(|m| m.name()),
        Some("Acetyl")
    );
    // (Acetyl)GALM(Oxidation)NEIQAAK
    let sequence = &ids[522].hits[0].sequence;
    assert_eq!(sequence.as_str(), "GALMNEIQAAK");
    assert_eq!(
        sequence.n_terminal_modification().map(|m| m.name()),
        Some("Acetyl")
    );
    assert_eq!(
        sequence.residue_modification(3).unwrap().map(|m| m.name()),
        Some("Oxidation")
    );
    // SHY(Phospho)GGSR
    let sequence = &ids[67].hits[0].sequence;
    assert_eq!(sequence.as_str(), "SHYGGSR");
    assert_eq!(
        sequence.residue_modification(2).unwrap().map(|m| m.name()),
        Some("Phospho")
    );

    let identifier = &protein.identifier;
    assert!(!identifier.is_empty());
    for identification in ids {
        assert_eq!(&identification.identifier, identifier);
    }
}

/// One attribute of an idXML element line, or `None` when it is absent.
///
/// `IdXMLFile::store` writes exactly one element per line with plain
/// double-quoted attributes, so the retained output can be read without an XML
/// parser. This keeps the differential independent of this crate's own idXML
/// reader, whose fixed 50-million-unit work budget cannot get through 586
/// modified peptides (see `docs/MASCOT_XML_SUPPORT.md`).
fn attribute<'a>(line: &'a str, name: &str) -> Option<&'a str> {
    let key = format!(" {name}=\"");
    let start = line.find(&key)? + key.len();
    let rest = line.get(start..)?;
    let end = rest.find('"')?;
    rest.get(..end)
}

/// A peptide hit as the retained C++ idXML output records it.
#[derive(Debug)]
struct ReferenceHit {
    score: f64,
    sequence: String,
    charge: i32,
    aa_before: Vec<String>,
    aa_after: Vec<String>,
    accessions: Vec<String>,
    evalue: f64,
}

/// A peptide identification as the retained C++ idXML output records it.
#[derive(Debug)]
struct ReferenceIdentification {
    score_type: String,
    higher_score_better: bool,
    significance_threshold: f64,
    mz: Option<f64>,
    rt: Option<f64>,
    hits: Vec<ReferenceHit>,
}

/// Everything the retained output declares, read line by line.
#[derive(Debug, Default)]
struct Reference {
    protein_accessions: Vec<String>,
    protein_scores: Vec<f64>,
    protein_ids: Vec<String>,
    fixed_modifications: Vec<String>,
    variable_modifications: Vec<String>,
    search: BTreeMap<String, String>,
    run: BTreeMap<String, String>,
    identifications: Vec<ReferenceIdentification>,
}

fn read_reference(path: &str) -> Reference {
    let text = std::fs::read_to_string(path).expect("retained C++ idXML output");
    let mut reference = Reference::default();
    for line in text.lines() {
        let line = line.trim();
        if line.starts_with("<SearchParameters ") {
            for key in [
                "db",
                "db_version",
                "taxonomy",
                "mass_type",
                "charges",
                "enzyme",
                "missed_cleavages",
                "precursor_peak_tolerance",
                "precursor_peak_tolerance_ppm",
                "peak_mass_tolerance",
                "peak_mass_tolerance_ppm",
            ] {
                if let Some(value) = attribute(line, key) {
                    reference.search.insert(key.to_owned(), value.to_owned());
                }
            }
        } else if line.starts_with("<IdentificationRun ") {
            for key in ["date", "search_engine", "search_engine_version"] {
                if let Some(value) = attribute(line, key) {
                    reference.run.insert(key.to_owned(), value.to_owned());
                }
            }
        } else if line.starts_with("<FixedModification ") {
            reference
                .fixed_modifications
                .push(attribute(line, "name").unwrap().to_owned());
        } else if line.starts_with("<VariableModification ") {
            reference
                .variable_modifications
                .push(attribute(line, "name").unwrap().to_owned());
        } else if line.starts_with("<ProteinHit ") {
            reference
                .protein_ids
                .push(attribute(line, "id").unwrap().to_owned());
            reference
                .protein_accessions
                .push(attribute(line, "accession").unwrap().to_owned());
            reference
                .protein_scores
                .push(attribute(line, "score").unwrap().parse().unwrap());
        } else if line.starts_with("<PeptideIdentification ") {
            reference.identifications.push(ReferenceIdentification {
                score_type: attribute(line, "score_type").unwrap().to_owned(),
                higher_score_better: attribute(line, "higher_score_better").unwrap() == "true",
                significance_threshold: attribute(line, "significance_threshold")
                    .unwrap()
                    .parse()
                    .unwrap(),
                mz: attribute(line, "MZ").map(|v| v.parse().unwrap()),
                rt: attribute(line, "RT").map(|v| v.parse().unwrap()),
                hits: Vec::new(),
            });
        } else if line.starts_with("<PeptideHit ") {
            let words = |value: Option<&str>| -> Vec<String> {
                value
                    .unwrap_or("")
                    .split_whitespace()
                    .map(str::to_owned)
                    .collect()
            };
            let refs = words(attribute(line, "protein_refs"));
            let accessions = refs
                .iter()
                .map(|id| {
                    let position = reference
                        .protein_ids
                        .iter()
                        .position(|value| value == id)
                        .expect("protein reference");
                    reference.protein_accessions[position].clone()
                })
                .collect();
            reference
                .identifications
                .last_mut()
                .expect("a PeptideHit inside a PeptideIdentification")
                .hits
                .push(ReferenceHit {
                    score: attribute(line, "score").unwrap().parse().unwrap(),
                    sequence: attribute(line, "sequence").unwrap().to_owned(),
                    charge: attribute(line, "charge").unwrap().parse().unwrap(),
                    aa_before: words(attribute(line, "aa_before")),
                    aa_after: words(attribute(line, "aa_after")),
                    accessions,
                    evalue: f64::NAN,
                });
        } else if line.starts_with("<UserParam ") && attribute(line, "name") == Some("EValue") {
            reference
                .identifications
                .last_mut()
                .unwrap()
                .hits
                .last_mut()
                .unwrap()
                .evalue = attribute(line, "value").unwrap().parse().unwrap();
        }
    }
    reference
}

#[test]
fn load_mascot_xml_2_2_against_the_retained_cpp_idxml_output() {
    // Tier-1 differential: MascotXMLFile_test_out_3.idXML is the output the
    // pinned C++ produced from MascotXMLFile_test_3.mascotXML through
    // MascotXMLFile::load followed by IdXMLFile::store, and the upstream
    // section fuzzy-compares against it with an acceptable absolute deviation
    // of 0.0001. Every record is compared here.
    let lookup = SpectrumTitleLookup::new();
    let result = mascot::load(TEST_3, &lookup).unwrap();
    let reference = read_reference(OUT_3);
    let ours = &result.protein_identification;

    // Search parameters.
    let search = &ours.search_parameters;
    assert_eq!(reference.search["db"], search.database);
    assert_eq!(reference.search["db_version"], search.database_version);
    assert_eq!(reference.search["taxonomy"], search.taxonomy);
    assert_eq!(reference.search["charges"], search.charges);
    assert_eq!(reference.search["mass_type"], "monoisotopic");
    assert_eq!(search.mass_type, PeakMassType::Monoisotopic);
    assert_eq!(
        reference.search["missed_cleavages"],
        search.missed_cleavages.to_string()
    );
    assert_eq!(reference.search["precursor_peak_tolerance_ppm"], "true");
    assert_eq!(reference.search["peak_mass_tolerance_ppm"], "false");
    assert_eq!(
        search.precursor_tolerance,
        Tolerance::Ppm(
            reference.search["precursor_peak_tolerance"]
                .parse()
                .unwrap()
        )
    );
    assert_eq!(
        search.fragment_tolerance,
        Tolerance::Absolute(reference.search["peak_mass_tolerance"].parse().unwrap())
    );
    assert_eq!(search.fixed_modifications, reference.fixed_modifications);
    assert_eq!(
        search.variable_modifications,
        reference.variable_modifications
    );
    // The idXML writer lower-cases the enzyme name; ProteaseDB spells it
    // "Trypsin".
    assert!(
        search
            .digestion_enzyme
            .eq_ignore_ascii_case(&reference.search["enzyme"]),
        "{} vs {}",
        search.digestion_enzyme,
        reference.search["enzyme"]
    );
    // The retained idXML records the run date as an XML dateTime, so only the
    // date/time separator differs from DateTime::get()'s space.
    assert_eq!(
        ours.date_time.as_deref(),
        Some(reference.run["date"].replace('T', " ").as_str())
    );
    assert_eq!(ours.search_engine, reference.run["search_engine"]);
    assert_eq!(
        ours.search_engine_version,
        reference.run["search_engine_version"]
    );

    // Protein hits, in document order.
    assert_eq!(ours.hits.len(), 7);
    assert_eq!(ours.hits.len(), reference.protein_accessions.len());
    for (index, hit) in ours.hits.iter().enumerate() {
        assert_eq!(hit.accession, reference.protein_accessions[index]);
        assert!(
            (hit.score - reference.protein_scores[index]).abs() < FUZZY_ABSOLUTE,
            "protein {index}"
        );
        assert!(hit.sequence.is_empty());
    }

    // Peptide identifications.
    let ours = &result.peptide_identifications;
    assert_eq!(ours.len(), 577);
    assert_eq!(ours.len(), reference.identifications.len());
    let mut total_hits = 0;
    for (index, (ours, reference)) in ours.iter().zip(&reference.identifications).enumerate() {
        assert_eq!(ours.score_type, reference.score_type, "id {index}");
        assert_eq!(
            ours.higher_score_better, reference.higher_score_better,
            "id {index}"
        );
        assert!(
            (ours.significance_threshold - reference.significance_threshold).abs() < FUZZY_ABSOLUTE,
            "id {index}"
        );
        match (ours.mz, reference.mz) {
            (Some(a), Some(b)) => {
                assert!((a - b).abs() < FUZZY_ABSOLUTE, "id {index}: {a} vs {b}")
            }
            (a, b) => panic!("id {index}: m/z {a:?} vs {b:?}"),
        }
        // Neither side has a retention time: the upstream section supplies an
        // empty lookup, so no <pep_scan_title> could be resolved.
        assert_eq!(ours.rt, None, "id {index}");
        assert_eq!(reference.rt, None, "id {index}");
        assert_eq!(ours.hits.len(), reference.hits.len(), "id {index}");
        total_hits += ours.hits.len();
        for (ours, reference) in ours.hits.iter().zip(&reference.hits) {
            assert!(
                (ours.score - reference.score).abs() < FUZZY_ABSOLUTE,
                "id {index}"
            );
            assert_eq!(ours.charge, reference.charge, "id {index}");
            assert_eq!(ours.sequence.to_string(), reference.sequence, "id {index}");
            assert_eq!(
                ours.evidences.len(),
                reference.accessions.len(),
                "id {index}"
            );
            let accessions: Vec<String> = ours
                .evidences
                .iter()
                .map(|e| e.protein_accession.clone())
                .collect();
            assert_eq!(accessions, reference.accessions, "id {index}");
            let before: Vec<String> = ours
                .evidences
                .iter()
                .map(|e| e.aa_before.code().to_string())
                .collect();
            let after: Vec<String> = ours
                .evidences
                .iter()
                .map(|e| e.aa_after.code().to_string())
                .collect();
            assert_eq!(before, reference.aa_before, "id {index}");
            assert_eq!(after, reference.aa_after, "id {index}");
            let evalue = ours.metadata["EValue"].as_f64().unwrap();
            assert!(
                (evalue - reference.evalue).abs() < FUZZY_ABSOLUTE,
                "id {index}: EValue {evalue} vs {}",
                reference.evalue
            );
        }
    }
    assert_eq!(total_hits, 586);
}

// ---------------------------------------------------------------------------
// START_SECTION((void load(const std::string&, ProteinIdentification&,
//                          PeptideIdentificationList&,
//                          std::map<std::string, std::vector<AASequence> >&,
//                          SpectrumMetaDataLookup&)))
// ---------------------------------------------------------------------------

#[test]
fn load_with_a_modified_peptide_map() {
    let mut sequence_1 = AASequence::parse("LHASGITVTEIPVTATNFK").unwrap();
    sequence_1.set_modification(16, "Deamidated").unwrap();
    let mut sequence_2 = AASequence::parse("MRSLGYVAVISAVATDTDK").unwrap();
    sequence_2.set_modification(2, "Phospho").unwrap();
    let mut sequence_3 = AASequence::parse("HSKLSAK").unwrap();
    sequence_3.set_modification(4, "Phospho").unwrap();
    let mut peptides: BTreeMap<String, Vec<AASequence>> = BTreeMap::new();
    peptides.insert(
        "789.83".to_owned(),
        vec![sequence_1.clone(), sequence_2.clone()],
    );
    peptides.insert("135.29".to_owned(), vec![sequence_3.clone()]);

    let lookup = SpectrumTitleLookup::new();
    let result = mascot::load_with_peptides(TEST_1, &peptides, &lookup).unwrap();
    let ids = &result.peptide_identifications;
    assert_eq!(ids.len(), 3);
    assert!((ids[0].mz.unwrap() - 789.83).abs() < FUZZY_ABSOLUTE);
    assert!((ids[1].mz.unwrap() - 135.29).abs() < FUZZY_ABSOLUTE);
    assert!((ids[2].mz.unwrap() - 982.58).abs() < FUZZY_ABSOLUTE);

    let protein = &result.protein_identification;
    assert_eq!(protein.hits.len(), 2);
    assert_eq!(protein.hits[0].accession, "AAN17824");
    assert_eq!(protein.hits[1].accession, "GN1736");
    assert_eq!(protein.hits[0].score, 619.0);
    assert_eq!(protein.hits[1].score, 293.0);
    assert_eq!(protein.score_type, "Mascot");
    assert_eq!(protein.date_time.as_deref(), Some("2006-03-09 11:31:52"));

    assert!((ids[0].significance_threshold - 31.8621).abs() < 1e-9);
    assert_eq!(ids[0].hits.len(), 2);
    assert_eq!(accessions(&ids[0], 0), ["AAN17824", "GN1736"]);
    assert_eq!(accessions(&ids[0], 1), ["AAN17824"]);
    assert_eq!(accessions(&ids[1], 0), ["GN1736"]);
    assert_eq!(ids[1].hits.len(), 1);
    assert!((ids[0].hits[0].score - 33.85).abs() < 1e-9);
    assert!((ids[0].hits[1].score - 33.12).abs() < 1e-9);
    assert!((ids[1].hits[0].score - 43.9).abs() < 1e-9);
    assert_eq!(ids[0].score_type, "Mascot");
    assert_eq!(ids[1].score_type, "Mascot");
    // The supplied sequences replaced the ones the file's fixed modifications
    // would have produced, and the fixed-modification pass was skipped.
    assert_eq!(ids[0].hits[0].sequence, sequence_1);
    assert_eq!(ids[0].hits[1].sequence, sequence_2);
    assert_eq!(ids[1].hits[0].sequence, sequence_3);
}

// ---------------------------------------------------------------------------
// Native tests: title lookup, source quirks, no-panic guarantees and bounds
// ---------------------------------------------------------------------------

fn document(body: &str) -> String {
    format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <mascot_search_results xmlns=\"http://www.matrixscience.com/xmlns/schema/mascot_search_results_2\" majorVersion=\"2\" minorVersion=\"1\">\n\
         {body}\n\
         </mascot_search_results>\n"
    )
}

#[test]
fn the_title_lookup_resolves_the_three_default_formats() {
    let experiment = MSExperiment {
        spectra: vec![
            MSSpectrum {
                native_id: "controllerType=0 controllerNumber=1 scan=818".into(),
                rt: 11.5,
                precursors: vec![Precursor::new(500.5, 2)],
                ..Default::default()
            },
            MSSpectrum {
                native_id: "controllerType=0 controllerNumber=1 scan=673".into(),
                rt: 22.5,
                precursors: vec![Precursor::new(600.5, 3)],
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    let (lookup, warnings) = MascotXmlFile::initialize_lookup(&experiment, None).unwrap();
    assert!(warnings.is_empty());
    for title in [
        "scan=818",
        "Spectrum136 scans:818,",
        "Spectrum3411 scans: 818,",
        "File773 Spectrum198145 scans: 818",
        "6860: Scan 818 (rt=5380.57)",
        "Scan Number: 818",
    ] {
        let meta = lookup.spectrum_meta_data(title, true).unwrap();
        assert_eq!(meta.rt, Some(11.5), "{title}");
        assert_eq!(meta.precursor_mz, Some(500.5), "{title}");
        assert_eq!(meta.scan_number, Some(818), "{title}");
    }
    // A DTA file name: scan 673 and charge 2 come from the name itself.
    let meta = lookup
        .spectrum_meta_data("/path/to/FTAC05_13.673.673.2.dta", true)
        .unwrap();
    assert_eq!(meta.rt, Some(22.5));
    assert_eq!(meta.precursor_mz, Some(600.5));
    // The m/z-underscore-RT form needs no spectrum at all.
    let mut direct = SpectrumTitleLookup::new();
    direct.add_reference_format(TitleReferenceFormat::MzThenRt);
    let meta = direct
        .spectrum_meta_data(
            "575.848571777344_5018.0811_controllerType=0 controllerNumber=1 scan=11515_EcoliMS2small",
            true,
        )
        .unwrap();
    assert!((meta.rt.unwrap() - 5018.0811).abs() < 1e-9);
    assert!((meta.precursor_mz.unwrap() - 575.848571777344).abs() < 1e-9);
    // An unresolvable scan number reports a range error, which the handler
    // downgrades to a warning.
    assert!(matches!(
        lookup.spectrum_meta_data("scan=999999", true),
        Err(Error::InvalidRange(_))
    ));
    // A title matching no format is not an error and yields nothing.
    let meta = lookup
        .spectrum_meta_data("no reference here", true)
        .unwrap();
    assert_eq!(meta, mascot::SpectrumMetaData::default());
}

#[test]
fn a_resolvable_scan_title_sets_the_retention_time() {
    let experiment = MSExperiment {
        spectra: vec![MSSpectrum {
            native_id: "scan=8162".into(),
            rt: 1234.5,
            precursors: vec![Precursor::new(324.7305, 2)],
            ..Default::default()
        }],
        ..Default::default()
    };
    let (lookup, _) = MascotXmlFile::initialize_lookup(&experiment, None).unwrap();
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <hits><hit number=\"1\"><protein accession=\"P1\">\n\
         <prot_score>10</prot_score>\n\
         <peptide query=\"1\"><pep_exp_mz>324.7305</pep_exp_mz><pep_exp_z>2</pep_exp_z>\n\
         <pep_score>14.59</pep_score><pep_seq>IPKFK</pep_seq>\n\
         <pep_scan_title>scan=8162</pep_scan_title></peptide>\n\
         </protein></hit></hits>",
    );
    let result = mascot::read(text.as_bytes(), &lookup).unwrap();
    assert_eq!(result.peptide_identifications.len(), 1);
    assert_eq!(result.peptide_identifications[0].rt, Some(1234.5));
    // A non-empty lookup that resolved nothing is the one condition the source
    // raises as an error rather than a warning.
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <hits><hit number=\"1\"><protein accession=\"P1\">\n\
         <prot_score>10</prot_score>\n\
         <peptide query=\"1\"><pep_exp_mz>324.7305</pep_exp_mz><pep_exp_z>2</pep_exp_z>\n\
         <pep_score>14.59</pep_score><pep_seq>IPKFK</pep_seq>\n\
         <pep_scan_title>nothing recognisable</pep_scan_title></peptide>\n\
         </protein></hit></hits>",
    );
    assert!(matches!(
        mascot::read(text.as_bytes(), &lookup),
        Err(Error::MissingInformation(_))
    ));
}

#[test]
fn a_string_title_supplies_a_retention_time_when_the_query_has_none() {
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <hits><hit number=\"1\"><protein accession=\"P1\">\n\
         <prot_score>10</prot_score>\n\
         <peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz><pep_exp_z>2</pep_exp_z>\n\
         <pep_score>12.0</pep_score><pep_seq>PEPTIDER</pep_seq></peptide>\n\
         </protein></hit></hits>\n\
         <queries><query number=\"1\"><StringTitle>500.0_1234.5</StringTitle></query></queries>",
    );
    let lookup = SpectrumTitleLookup::new();
    let result = mascot::read(text.as_bytes(), &lookup).unwrap();
    assert_eq!(result.peptide_identifications[0].rt, Some(1234.5));
    // An explicit <RTINSECONDS> wins, because it is read into the same slot.
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <hits><hit number=\"1\"><protein accession=\"P1\">\n\
         <prot_score>10</prot_score>\n\
         <peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz><pep_exp_z>2</pep_exp_z>\n\
         <pep_score>12.0</pep_score><pep_seq>PEPTIDER</pep_seq></peptide>\n\
         </protein></hit></hits>\n\
         <queries><query number=\"1\"><RTINSECONDS>99.5</RTINSECONDS>\n\
         <StringTitle>500.0_1234.5</StringTitle></query></queries>",
    );
    let result = mascot::read(text.as_bytes(), &lookup).unwrap();
    assert_eq!(result.peptide_identifications[0].rt, Some(99.5));
}

#[test]
fn a_missing_numqueries_header_is_refused_rather_than_read_out_of_bounds() {
    // The source guard is `peptide_identification_index_ > id_data_.size()`,
    // so index 0 on an empty vector passes it and reads out of bounds.
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date></header>\n\
         <hits><hit number=\"1\"><protein accession=\"P1\">\n\
         <prot_score>10</prot_score>\n\
         <peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz><pep_seq>PEPTIDER</pep_seq></peptide>\n\
         </protein></hit></hits>",
    );
    let lookup = SpectrumTitleLookup::new();
    match mascot::read(text.as_bytes(), &lookup).unwrap_err() {
        Error::Parse { message, .. } => {
            assert!(message.contains("show_header=1"), "{message}");
        }
        other => panic!("unexpected error {other:?}"),
    }
    // A query number one past the end trips the same guard in the source, which
    // compares with `>` instead of `>=`.
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <hits><hit number=\"1\"><protein accession=\"P1\">\n\
         <prot_score>10</prot_score>\n\
         <peptide query=\"2\"><pep_exp_mz>500.0</pep_exp_mz><pep_seq>PEPTIDER</pep_seq></peptide>\n\
         </protein></hit></hits>",
    );
    assert!(matches!(
        mascot::read(text.as_bytes(), &lookup),
        Err(Error::Parse { .. })
    ));
    // Query numbers count from one, so zero is refused instead of underflowing.
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <hits><hit number=\"1\"><protein accession=\"P1\">\n\
         <prot_score>10</prot_score>\n\
         <peptide query=\"0\"><pep_exp_mz>500.0</pep_exp_mz><pep_seq>PEPTIDER</pep_seq></peptide>\n\
         </protein></hit></hits>",
    );
    assert!(matches!(
        mascot::read(text.as_bytes(), &lookup),
        Err(Error::Parse { .. })
    ));
}

#[test]
fn a_hostile_numqueries_is_bounded_rather_than_allocated() {
    let lookup = SpectrumTitleLookup::new();
    // The source calls id_data_.resize() with the converted value, so a huge
    // count commits the memory and a negative one wraps to an enormous size.
    let huge = document("<header><NumQueries>2000000000</NumQueries></header>");
    assert!(matches!(
        mascot::read(huge.as_bytes(), &lookup),
        Err(Error::InvalidValue(_))
    ));
    let negative = document("<header><NumQueries>-1</NumQueries></header>");
    assert!(matches!(
        mascot::read(negative.as_bytes(), &lookup),
        Err(Error::Parse { .. })
    ));
    // An explicit lower ceiling is honoured.
    let text = document("<header><NumQueries>10</NumQueries></header>");
    let limits = mascot::ReadLimits {
        max_queries: 5,
        ..Default::default()
    };
    assert!(matches!(
        mascot::read_with_options(text.as_bytes(), &BTreeMap::new(), &lookup, &limits),
        Err(Error::InvalidValue(_))
    ));
    // Limits above the hard ceilings are refused before anything is read.
    let limits = mascot::ReadLimits {
        max_events: usize::MAX,
        ..Default::default()
    };
    assert!(matches!(
        mascot::read_with_options(text.as_bytes(), &BTreeMap::new(), &lookup, &limits),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn resource_ceilings_bound_bytes_events_depth_and_text() {
    let lookup = SpectrumTitleLookup::new();
    let text = document("<header><NumQueries>1</NumQueries></header>");
    for limits in [
        mascot::ReadLimits {
            max_bytes: 16,
            ..Default::default()
        },
        mascot::ReadLimits {
            max_events: 3,
            ..Default::default()
        },
        mascot::ReadLimits {
            max_depth: 1,
            ..Default::default()
        },
        mascot::ReadLimits {
            max_text_bytes: 1,
            ..Default::default()
        },
    ] {
        assert!(
            mascot::read_with_options(text.as_bytes(), &BTreeMap::new(), &lookup, &limits).is_err(),
            "{limits:?}"
        );
    }
}

#[test]
fn a_warning_element_removes_a_modification_from_the_fixed_list() {
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries>\n\
         <warning number=\"0\">'Phospho (ST)' can only be used as a variable modification; it has been changed</warning>\n\
         </header>\n\
         <fixed_mods><modification identifier=\"1\"><name>Phospho (ST)</name></modification>\n\
         <modification identifier=\"2\"><name>Carbamidomethyl (C)</name></modification></fixed_mods>\n\
         <variable_mods><modification identifier=\"1\"><name>Oxidation (M)</name></modification></variable_mods>",
    );
    let lookup = SpectrumTitleLookup::new();
    let result = mascot::read(text.as_bytes(), &lookup).unwrap();
    let search = &result.protein_identification.search_parameters;
    assert_eq!(search.fixed_modifications, ["Carbamidomethyl (C)"]);
    assert_eq!(search.variable_modifications, ["Oxidation (M)"]);
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.contains("Modification removed as fixed modification")),
        "{:?}",
        result.warnings
    );
}

#[test]
fn mods_and_it_mods_are_read_only_without_the_dedicated_sections() {
    let lookup = SpectrumTitleLookup::new();
    // Without <fixed_mods>/<variable_mods>, the flat lists are used and the
    // specificity groups are expanded.
    let flat = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <search_parameters><MODS>Carboxymethyl (C),Deamidated (NQ)</MODS>\n\
         <IT_MODS>Oxidation (M),Phospho (ST)</IT_MODS></search_parameters>",
    );
    let search = mascot::read(flat.as_bytes(), &lookup)
        .unwrap()
        .protein_identification
        .search_parameters;
    assert_eq!(
        search.fixed_modifications,
        ["Carboxymethyl (C)", "Deamidated (N)", "Deamidated (Q)"]
    );
    assert_eq!(
        search.variable_modifications,
        ["Oxidation (M)", "Phospho (S)", "Phospho (T)"]
    );
    // With the dedicated sections present, the flat lists are ignored.
    let sectioned = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <fixed_mods><modification identifier=\"1\"><name>Carbamidomethyl (C)</name></modification></fixed_mods>\n\
         <variable_mods><modification identifier=\"1\"><name>Oxidation (M)</name></modification></variable_mods>\n\
         <search_parameters><MODS>Carboxymethyl (C)</MODS><IT_MODS>Phospho (Y)</IT_MODS></search_parameters>",
    );
    let search = mascot::read(sectioned.as_bytes(), &lookup)
        .unwrap()
        .protein_identification
        .search_parameters;
    assert_eq!(search.fixed_modifications, ["Carbamidomethyl (C)"]);
    assert_eq!(search.variable_modifications, ["Oxidation (M)"]);
    // The guard is list emptiness, not section presence: an *empty* section
    // contributes no name, so a later <MODS>/<IT_MODS> still populates the
    // list, and the <IT_MODS> fallback expands its specificity groups as it
    // reads rather than at the document end. A second-model review corrected
    // this document's earlier "never both" claim.
    let empty_sections = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <fixed_mods/><variable_mods/>\n\
         <search_parameters><MODS>Carbamidomethyl (C)</MODS>\n\
         <IT_MODS>Phospho (ST)</IT_MODS></search_parameters>",
    );
    let search = mascot::read(empty_sections.as_bytes(), &lookup)
        .unwrap()
        .protein_identification
        .search_parameters;
    assert_eq!(search.fixed_modifications, ["Carbamidomethyl (C)"]);
    assert_eq!(
        search.variable_modifications,
        ["Phospho (S)", "Phospho (T)"]
    );
    // An expansion the modification database does not know is refused, where
    // the source throws ElementNotFound.
    let unknown = document(
        "<header><NumQueries>1</NumQueries></header>\n\
         <search_parameters><MODS>Phospho (Z)</MODS></search_parameters>",
    );
    assert!(matches!(
        mascot::read(unknown.as_bytes(), &lookup),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn a_variable_modification_index_beyond_the_list_is_refused() {
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <variable_mods><modification identifier=\"1\"><name>Oxidation (M)</name></modification></variable_mods>\n\
         <hits><hit number=\"1\"><protein accession=\"P1\"><prot_score>10</prot_score>\n\
         <peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz><pep_score>1.0</pep_score>\n\
         <pep_seq>MPEPTIDE</pep_seq><pep_var_mod_pos>0.90000000.0</pep_var_mod_pos></peptide>\n\
         </protein></hit></hits>",
    );
    let lookup = SpectrumTitleLookup::new();
    match mascot::read(text.as_bytes(), &lookup).unwrap_err() {
        Error::Parse { message, .. } => {
            assert!(
                message.contains("exceeds the declared variable"),
                "{message}"
            );
        }
        other => panic!("unexpected error {other:?}"),
    }
    // A modification marked beyond the sequence is refused too, where the
    // source would index the residue vector out of range.
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <variable_mods><modification identifier=\"1\"><name>Oxidation (M)</name></modification></variable_mods>\n\
         <hits><hit number=\"1\"><protein accession=\"P1\"><prot_score>10</prot_score>\n\
         <peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz><pep_score>1.0</pep_score>\n\
         <pep_seq>M</pep_seq><pep_var_mod_pos>0.001.0</pep_var_mod_pos></peptide>\n\
         </protein></hit></hits>",
    );
    match mascot::read(text.as_bytes(), &lookup).unwrap_err() {
        Error::Parse { message, .. } => {
            assert!(message.contains("beyond the peptide sequence"), "{message}");
        }
        other => panic!("unexpected error {other:?}"),
    }
}

#[test]
fn a_repeated_first_hit_is_collapsed() {
    // Mascot 2.2 repeats the first hit; the source erases the entry at index 1.
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <unassigned>\n\
         <u_peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz><pep_exp_z>2</pep_exp_z>\n\
         <pep_score>12.0</pep_score><pep_seq>PEPTIDER</pep_seq></u_peptide>\n\
         <u_peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz><pep_exp_z>2</pep_exp_z>\n\
         <pep_score>12.0</pep_score><pep_seq>PEPTIDER</pep_seq></u_peptide>\n\
         <u_peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz><pep_exp_z>2</pep_exp_z>\n\
         <pep_score>9.0</pep_score><pep_seq>PEPTIDEK</pep_seq></u_peptide>\n\
         </unassigned>",
    );
    let lookup = SpectrumTitleLookup::new();
    let result = mascot::read(text.as_bytes(), &lookup).unwrap();
    let hits = &result.peptide_identifications[0].hits;
    assert_eq!(hits.len(), 2);
    assert_eq!(hits[0].sequence.as_str(), "PEPTIDER");
    assert_eq!(hits[1].sequence.as_str(), "PEPTIDEK");
}

#[test]
fn identifications_without_a_sequence_are_dropped_with_a_count() {
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>3</NumQueries></header>\n\
         <unassigned>\n\
         <u_peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz><pep_score>1.0</pep_score>\n\
         <pep_seq></pep_seq></u_peptide>\n\
         <u_peptide query=\"2\"><pep_exp_mz>600.0</pep_exp_mz><pep_score>2.0</pep_score>\n\
         <pep_seq>PEPTIDER</pep_seq></u_peptide>\n\
         </unassigned>",
    );
    let lookup = SpectrumTitleLookup::new();
    let result = mascot::read(text.as_bytes(), &lookup).unwrap();
    // Query 1 had a single sequence-less hit and is dropped and counted; query
    // 3 was only padding from NumQueries and is dropped silently.
    assert_eq!(result.peptide_identifications.len(), 1);
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.contains("Removed 1 peptide identifications without sequence")),
        "{:?}",
        result.warnings
    );
}

#[test]
fn malformed_and_non_ascii_documents_are_handled_without_panicking() {
    let lookup = SpectrumTitleLookup::new();
    // Not well formed.
    assert!(mascot::read(b"<mascot_search_results".as_slice(), &lookup).is_err());
    // Unbalanced closing tag.
    let text = "<mascot_search_results majorVersion=\"2\" minorVersion=\"1\"></header></mascot_search_results>";
    assert!(mascot::read(text.as_bytes(), &lookup).is_err());
    // Missing majorVersion attribute.
    assert!(mascot::read(b"<mascot_search_results/>".as_slice(), &lookup).is_err());
    // Non-UTF-8 bytes.
    let mut bytes = document("<header><NumQueries>1</NumQueries></header>").into_bytes();
    bytes.push(0xff);
    assert!(mascot::read(bytes.as_slice(), &lookup).is_err());
    // Non-ASCII text in a protein accession and description survives intact.
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <hits><hit number=\"1\"><protein accession=\"日本語\">\n\
         <prot_desc>Caffeïne — ☕</prot_desc><prot_score>10</prot_score>\n\
         <peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz><pep_score>1.0</pep_score>\n\
         <pep_seq>PEPTIDER</pep_seq></peptide></protein></hit></hits>",
    );
    let result = mascot::read(text.as_bytes(), &lookup).unwrap();
    assert_eq!(result.protein_identification.hits[0].accession, "日本語");
    // A non-ASCII residue letter is not a peptide and is refused.
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <unassigned><u_peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz>\n\
         <pep_score>1.0</pep_score><pep_seq>PEPTÏDER</pep_seq></u_peptide></unassigned>",
    );
    assert!(mascot::read(text.as_bytes(), &lookup).is_err());
    // A file that does not exist.
    assert!(matches!(
        mascot::load("tests/data/does_not_exist.mascotXML", &lookup),
        Err(Error::Io(_))
    ));
}

#[test]
fn a_fractional_protein_score_is_a_conversion_error() {
    // The source reads <prot_score> with toInt32, so a fractional value throws.
    let text = document(
        "<header><NumQueries>1</NumQueries></header>\n\
         <hits><hit number=\"1\"><protein accession=\"P1\">\n\
         <prot_score>619.5</prot_score></protein></hit></hits>",
    );
    let lookup = SpectrumTitleLookup::new();
    assert!(matches!(
        mascot::read(text.as_bytes(), &lookup),
        Err(Error::Parse { .. })
    ));
}

#[test]
fn the_homology_threshold_only_replaces_a_larger_identity_threshold() {
    let build = |homology: &str, identity: &str| {
        document(&format!(
            "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
             <unassigned><u_peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz>\n\
             <pep_score>1.0</pep_score><pep_homol>{homology}</pep_homol>\n\
             <pep_ident>{identity}</pep_ident><pep_seq>PEPTIDER</pep_seq></u_peptide></unassigned>"
        ))
    };
    let lookup = SpectrumTitleLookup::new();
    let threshold = |homology: &str, identity: &str| {
        mascot::read(build(homology, identity).as_bytes(), &lookup)
            .unwrap()
            .peptide_identifications[0]
            .significance_threshold
    };
    // Homology larger than identity: identity wins.
    assert!((threshold("34", "31.8621") - 31.8621).abs() < 1e-9);
    // Homology smaller than identity: homology is kept.
    assert!((threshold("19", "46") - 19.0).abs() < 1e-9);
    // Homology absent (zero): identity wins.
    assert!((threshold("0", "46") - 46.0).abs() < 1e-9);
    // Both thresholds are also recorded on the hit.
    let result = mascot::read(build("34", "31.8621").as_bytes(), &lookup).unwrap();
    let hit = &result.peptide_identifications[0].hits[0];
    assert!((hit.metadata["homology_threshold"].as_f64().unwrap() - 34.0).abs() < 1e-9);
    assert!((hit.metadata["identity_threshold"].as_f64().unwrap() - 31.8621).abs() < 1e-9);
}

#[test]
fn the_file_struct_and_free_functions_agree() {
    let lookup = SpectrumTitleLookup::new();
    let file = MascotXmlFile::new();
    let a = file.load(TEST_1, &lookup).unwrap();
    let b = mascot::load(TEST_1, &lookup).unwrap();
    assert_eq!(a, b);
    let peptides = BTreeMap::new();
    let c = file.load_with_peptides(TEST_1, &peptides, &lookup).unwrap();
    assert_eq!(a, c);
}

// ---------------------------------------------------------------------------
// XML shapes this reader refuses rather than silently mishandling
// ---------------------------------------------------------------------------

/// An entity reference in element text is EXPANDED, as Xerces expands it.
///
/// quick-xml does not expand references inside text: it splits the text at
/// every `&...;` and emits the reference as its own `Event::GeneralRef`, and
/// `BytesText::xml_content()` only decodes and normalises line endings. The
/// catch-all event arm therefore DELETED the reference and concatenated the
/// surrounding fragments, so `<pep_score>1&#46;5</pep_score>` was read as the
/// score **15** — silent corruption of a scientific value. Xerces expands the
/// reference before the C++ handler sees any character data
/// (`XMLFile.cpp:78-94`, `MascotXMLHandler.cpp:600-611`), and
/// so the source's behaviour is expansion, and this reader now matches it.
///
/// A first attempt refused every reference, copying `src/format/mzml.rs` and
/// `src/format/imzml_handler.rs`. That was wrong here: those readers consume
/// only Base64 text, which never contains `&`, whereas a Mascot protein
/// description holding `&amp;` is ordinary valid output, and refusing it failed
/// the whole load. Predefined and numeric references are expanded; a reference
/// that would need a DTD is still an error rather than a silent deletion.
///
/// Found by a second-model review of this port, and the over-correction by a
/// second-model review of the fix. The same original defect was found in the
/// shipped imzML reader by following it there.
#[test]
fn entity_references_are_expanded_and_unresolvable_shapes_refused() {
    let lookup = SpectrumTitleLookup::new();
    let body = |score: &str, sequence: &str| {
        format!(
            "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
             <unassigned><u_peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz>\n\
             <pep_score>{score}</pep_score><pep_seq>{sequence}</pep_seq></u_peptide></unassigned>"
        )
    };
    // The exact input from the finding: it read as 15 before the fix, and now
    // reads the 1.5 the file means, which is what Xerces gives the C++ handler.
    let text = document(&body("1&#46;5", "PEPTIDER"));
    let result = mascot::read(text.as_bytes(), &lookup).unwrap();
    let score = result.peptide_identifications[0].hits[0].score;
    assert!((score - 1.5).abs() < 1e-12, "{score}");
    // Hexadecimal form of the same character reference.
    let text = document(&body("1&#x2E;5", "PEPTIDER"));
    let result = mascot::read(text.as_bytes(), &lookup).unwrap();
    let score = result.peptide_identifications[0].hits[0].score;
    assert!((score - 1.5).abs() < 1e-12, "{score}");
    // All five predefined entities must reach the value, not vanish from it.
    // A peptide sequence is the sharpest probe available: the expanded
    // character is not a residue, so expansion surfaces as an amino-acid error,
    // whereas the old deletion left a clean "PEPTIDER" that parsed
    // successfully. An InvalidValue here therefore proves the character
    // arrived; an Ok would prove it had been dropped.
    for reference in ["&amp;", "&lt;", "&gt;", "&quot;", "&apos;"] {
        let text = document(&body("1.5", &format!("PEP{reference}TIDER")));
        match mascot::read(text.as_bytes(), &lookup) {
            Err(Error::InvalidValue(message)) => {
                assert!(message.contains("amino-acid"), "{reference}: {message}");
            }
            Ok(_) => panic!("{reference} was deleted from the sequence, not expanded"),
            Err(other) => panic!("{reference}: unexpected {other:?}"),
        }
    }
    // And in a numeric field the expanded character is simply the character the
    // file wrote, so the value parses to what the file means.
    let text = document(&body("1&#46;5", "PEPTIDER"));
    assert!(mascot::read(text.as_bytes(), &lookup).is_ok());
    // The undeclared entity that Xerces rejects is refused here too.
    let text = document(&body("1.5", "P&bogus;EPTIDER"));
    assert!(mascot::read(text.as_bytes(), &lookup).is_err());
    // A CDATA section is refused rather than read as plain text.
    let text = document(&body("1.5", "<![CDATA[PEPTIDER]]>"));
    let error = mascot::read(text.as_bytes(), &lookup).unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error:?}");
    // A DTD is refused before any element is read.
    let text = format!(
        "<!DOCTYPE mascot_search_results>\n{}",
        document(&body("1.5", "PEPTIDER"))
    );
    let error = mascot::read(text.as_bytes(), &lookup).unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error:?}");
    // Without the reference the same document reads the score the file means.
    let text = document(&body("1.5", "PEPTIDER"));
    let result = mascot::read(text.as_bytes(), &lookup).unwrap();
    let score = result.peptide_identifications[0].hits[0].score;
    assert!((score - 1.5).abs() < 1e-12, "{score}");
}

/// A document must be one complete element tree, as Xerces requires.
///
/// `check_end_names` only pairs the tags quick-xml actually sees, so a document
/// truncated before its closing root tag used to load successfully — the
/// root-close post-processing simply never ran — and leading junk or a second
/// concatenated document were ignored. Found by a second-model review.
#[test]
fn an_incomplete_or_multi_root_document_is_refused() {
    let lookup = SpectrumTitleLookup::new();
    let parse_error = |text: String| {
        let error = mascot::read(text.as_bytes(), &lookup).unwrap_err();
        assert!(matches!(error, Error::Parse { .. }), "{error:?} for {text}");
    };
    // Truncated before the closing root tag.
    parse_error(
        "<mascot_search_results majorVersion=\"2\" minorVersion=\"1\">\n\
         <NumQueries>1</NumQueries>\n\
         <u_peptide query=\"1\"><pep_seq>PEPTIDER</pep_seq></u_peptide>\n"
            .to_owned(),
    );
    // Character data before the root element.
    parse_error(format!(
        "junk{}",
        document("<header><NumQueries>1</NumQueries></header>")
    ));
    // Character data after it.
    parse_error(format!(
        "{}junk",
        document("<header><NumQueries>1</NumQueries></header>")
    ));
    // Two concatenated documents.
    parse_error(format!(
        "{}{}",
        document("<header><DB>A</DB><NumQueries>1</NumQueries></header>"),
        document("<header><DB>B</DB><NumQueries>1</NumQueries></header>")
    ));
    // White space around the root element is legal, as every fixture has.
    let text = format!(
        "\n  {}\n\n",
        document("<header><NumQueries>1</NumQueries></header>")
    );
    assert!(mascot::read(text.as_bytes(), &lookup).is_ok());
}

/// `<NumQueries>` is an index space, not an allocation.
///
/// The source calls `id_data_.resize(...)` with the declared count, so a
/// five-million-query header commits the whole vector — and repeating the
/// element constructs and drops it again. The count is recorded here and
/// entries are materialised only when a query references them; a *smaller*
/// repeat still truncates, which is the data loss `resize` causes. Measured on
/// the Linux gate node in release mode for this test: 1.59 s before, 0.02 s
/// after.
#[test]
fn a_declared_query_count_is_not_materialised_up_front() {
    let lookup = SpectrumTitleLookup::new();
    let peptide = "<unassigned><u_peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz>\n\
                   <pep_score>1.0</pep_score><pep_seq>PEPTIDER</pep_seq></u_peptide></unassigned>";
    let text = document(&format!(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>5000000</NumQueries></header>\n{peptide}"
    ));
    let result = mascot::read(text.as_bytes(), &lookup).unwrap();
    assert_eq!(result.peptide_identifications.len(), 1);
    // Repeated counts no longer multiply the work.
    let text = document(&format!(
        "<header><Date>2012-03-15T14:20:09Z</Date>\n\
         <NumQueries>5000000</NumQueries><NumQueries>5000000</NumQueries>\n\
         <NumQueries>5000000</NumQueries><NumQueries>5000000</NumQueries></header>\n{peptide}"
    ));
    assert_eq!(
        mascot::read(text.as_bytes(), &lookup)
            .unwrap()
            .peptide_identifications
            .len(),
        1
    );
    // A smaller repeat truncates, so a later reference is out of range.
    let text = document(&format!(
        "<header><NumQueries>5000000</NumQueries><NumQueries>0</NumQueries></header>\n{peptide}"
    ));
    assert!(matches!(
        mascot::read(text.as_bytes(), &lookup),
        Err(Error::Parse { .. })
    ));
    // `std::from_chars` refuses a second '+', so `++1` is not a query count.
    let text = document("<header><NumQueries>++1</NumQueries></header>");
    assert!(matches!(
        mascot::read(text.as_bytes(), &lookup),
        Err(Error::Parse { .. })
    ));
}

/// A reference format that matched but whose value does not convert is an
/// error, not a miss.
///
/// `SpectrumMetaDataLookup::getSpectrumMetaData` returns after the first
/// matching expression, and `toInt32`/`toDouble` throw from inside that block,
/// which `MascotXMLHandler` catches as a warning. Falling through to the next
/// format instead invents a successful association from a different part of the
/// title. Found by a second-model review.
#[test]
fn a_matched_format_whose_value_does_not_convert_does_not_fall_through() {
    let experiment = MSExperiment {
        spectra: vec![MSSpectrum {
            native_id: "scan=1".into(),
            rt: 5.0,
            precursors: vec![Precursor::new(100.0, 1)],
            ..Default::default()
        }],
        ..Default::default()
    };
    let (lookup, _) = MascotXmlFile::initialize_lookup(&experiment, None).unwrap();
    // The scan-number format matches; its digits overflow `toInt32`; the
    // m/z-underscore-RT format is never tried, so RT 12 and m/z 500 are not
    // invented from the leading "500_12".
    let error = lookup
        .spectrum_meta_data("500_12 scan=9223372036854775808", true)
        .unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
    // An overflowing title retention time is a conversion error too, not the
    // infinity Rust's own parser produces.
    let mut direct = SpectrumTitleLookup::new();
    direct.add_reference_format(TitleReferenceFormat::MzThenRt);
    let error = direct
        .spectrum_meta_data(&format!("500_{}", "9".repeat(320)), true)
        .unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
    // Through the handler the conversion error is a warning, the retention time
    // stays unset, and the non-empty lookup then makes that an error overall.
    let text = document(
        "<header><Date>2012-03-15T14:20:09Z</Date><NumQueries>1</NumQueries></header>\n\
         <unassigned><u_peptide query=\"1\"><pep_exp_mz>500.0</pep_exp_mz>\n\
         <pep_score>1.0</pep_score><pep_seq>PEPTIDER</pep_seq>\n\
         <pep_scan_title>500_12 scan=9223372036854775808</pep_scan_title>\n\
         </u_peptide></unassigned>",
    );
    assert!(matches!(
        mascot::read(text.as_bytes(), &lookup),
        Err(Error::MissingInformation(_))
    ));
}

/// Boost's `^` and `$` are line anchors, and scan numbers are 32-bit.
///
/// `initializeLookup` compiles its expressions with the default perl flags, so
/// `^` becomes `syntax_element_start_line` and `$` `syntax_element_end_line`
/// (`basic_regex_parser.hpp`): a wrapped title matches on its later lines and a
/// native ID with a trailing annotation line still yields its scan number.
/// `SpectrumLookup` converts with `toInt32`, so a longer digit run is a
/// conversion failure — reported as the source's `-1` sentinel, i.e. no
/// scan-number entry and a warning. Found by a second-model review.
#[test]
fn the_title_anchors_are_line_anchors_and_scan_numbers_are_32_bit() {
    let mut direct = SpectrumTitleLookup::new();
    direct.add_reference_format(TitleReferenceFormat::MzThenRt);
    let meta = direct
        .spectrum_meta_data("exported spectrum\n500.5_12.5", true)
        .unwrap();
    assert_eq!(meta.precursor_mz, Some(500.5));
    assert_eq!(meta.rt, Some(12.5));
    // Not mid-line, though: `^` is still an anchor.
    let meta = direct.spectrum_meta_data("x 500.5_12.5", true).unwrap();
    assert_eq!(meta, mascot::SpectrumMetaData::default());
    // `$` in the native-ID expression `=(?<SCAN>\d+)$` is a line anchor too.
    let experiment = MSExperiment {
        spectra: vec![MSSpectrum {
            native_id: "scan=818\nannotation".into(),
            rt: 11.5,
            ..Default::default()
        }],
        ..Default::default()
    };
    let (lookup, warnings) = MascotXmlFile::initialize_lookup(&experiment, None).unwrap();
    assert!(warnings.is_empty(), "{warnings:?}");
    assert_eq!(lookup.spectrum(0).unwrap().scan_number, Some(818));
    assert_eq!(
        lookup.spectrum_meta_data("scan=818", true).unwrap().rt,
        Some(11.5)
    );
    // A trailing digit run too long for `toInt32` records no scan number and
    // produces the source's warning.
    let experiment = MSExperiment {
        spectra: vec![MSSpectrum {
            native_id: "scan=2147483648".into(),
            rt: 11.5,
            ..Default::default()
        }],
        ..Default::default()
    };
    let (lookup, warnings) = MascotXmlFile::initialize_lookup(&experiment, None).unwrap();
    assert_eq!(warnings.len(), 1);
    assert!(
        warnings[0].contains("Could not extract scan number"),
        "{warnings:?}"
    );
    assert_eq!(lookup.spectrum(0).unwrap().scan_number, None);
}

/// The scan-number matcher commits to its optional groups, and that cannot
/// change the outcome.
///
/// `[Ss]can( [Nn]umber)?s?[=:]? *(?<SCAN>\d+)` has two greedy optional groups
/// that a real engine would backtrack out of. It never needs to: if ` Number`
/// or the trailing `s` is present and consuming it fails, *not* consuming it
/// requires a digit at the `N` or at the `s` itself, which those characters are
/// not. This asserts the Rust side of that table; the C++ side is Boost's
/// leftmost-match semantics. Note that the leftmost match wins here, unlike the
/// MGF writer's accession table, where the token iterator takes the last.
///
/// Raised as a suspicion by the porting agent and refuted by a second-model
/// review; kept as a test so it is not re-suspected.
#[test]
fn the_scan_number_matcher_needs_no_backtracking() {
    let experiment = MSExperiment {
        spectra: (1..=13)
            .map(|scan| MSSpectrum {
                native_id: format!("scan={scan}"),
                rt: f64::from(scan),
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    };
    let mut lookup = SpectrumTitleLookup::new();
    lookup.read_spectra(&experiment).unwrap();
    lookup.add_reference_format(TitleReferenceFormat::ScanNumber);
    let scan = |title: &str| lookup.spectrum_meta_data(title, false).unwrap().scan_number;
    // Committing to " Number" or to the trailing 's' never loses a match.
    assert_eq!(scan("Scans followed later by 5"), None);
    assert_eq!(scan("Scan Numbers 5"), Some(5));
    assert_eq!(scan("Scan Numberx5"), None);
    assert_eq!(scan("Scan Number scan=7"), Some(7));
    // `[=:]?` precedes ` *`, so a separator after a space is not accepted.
    assert_eq!(scan("Scan Numbers =5"), None);
    // The expression is unanchored, and leftmost wins.
    assert_eq!(scan("xscan=7"), Some(7));
    assert_eq!(scan("scan=1 scan=2"), Some(1));
    assert_eq!(scan("scan=00012"), Some(12));
    assert_eq!(scan("scan=+5"), None);
    assert_eq!(scan("scanx=3"), None);
}
