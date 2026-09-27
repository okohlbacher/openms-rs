// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/percolator_outfile.rs` (`FORMAT/PercolatorOutfile.h`)
//! and the reference path of `src/format/spectrum_meta_data_lookup.rs`.
//!
//! Evidence, strongest first (see `docs/PERCOLATOR_OUTFILE_SUPPORT.md` and
//! `tests/data/percolator_outfile_provenance.json`):
//!
//! - tier 1, retained upstream output: `TOPP_IDFileConverter_12`, `_13` and
//!   `_14` convert `IDFileConverter_12_input.psms` to idXML with
//!   `-score_type qvalue`, `PEP` and `score`; the retained C++ outputs are read
//!   back with this crate's idXML reader and compared record by record;
//! - tier 3: every `START_SECTION` of `PercolatorOutfile_test.cpp`, with its
//!   literals;
//! - tier 3/4: the reference look-up, peptide rewriting and refusal cases,
//!   derived from reading the source, with UniMod masses for the mass shifts.
//!
//! The reader needs no feature; the file is gated on `idxml` because the
//! differential reads the retained outputs with the idXML reader.
#![cfg(feature = "idxml")]

use openms::format::idxml;
use openms::format::percolator_outfile::{
    DEFAULT_REFERENCE_FORMATS, HEADER, PercolatorDocument, PercolatorOutfile, SCORE_TYPE_NAMES,
    ScoreType,
};
use openms::format::spectrum_meta_data_lookup::{
    DEFAULT_SCAN_REGEXP, MetaDataFlags, SpectrumMetaData, SpectrumMetaDataLookup,
};
use openms::identification::PeptideIdentification;
use openms::{Error, MSSpectrum, Precursor};
use std::collections::BTreeSet;

const CLASS_TEST: &str = "tests/data/percolator_outfile/PercolatorOutfile_test.psms";
const TOPP_INPUT: &str = "tests/data/percolator_outfile/IDFileConverter_12_input.psms";

/// `TEST_REAL_SIMILAR`'s default relative tolerance in the class tests.
fn similar(actual: f64, expected: f64) {
    let tolerance = 1e-5 * expected.abs().max(actual.abs());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} is not similar to {expected}"
    );
}

fn read(text: &str, lookup: &mut SpectrumMetaDataLookup) -> openms::Result<PercolatorDocument> {
    PercolatorOutfile::new().read(text.as_bytes(), lookup, ScoreType::QValue)
}

fn row(reference: &str, peptide: &str) -> String {
    format!("{HEADER}\n{reference}\t1.5\t0.01\t0.02\t{peptide}\tP1\n")
}

// ---------------------------------------------------------------------------
// START_SECTION(ScoreType getScoreType(std::string score_type_name))
// ---------------------------------------------------------------------------

#[test]
fn score_type_names_follow_the_class_test() {
    assert_eq!(ScoreType::from_name("qvalue").unwrap(), ScoreType::QValue);
    assert_eq!(ScoreType::from_name("q-value").unwrap(), ScoreType::QValue);
    assert_eq!(
        ScoreType::from_name("PEP").unwrap(),
        ScoreType::PosteriorErrorProbability
    );
    assert_eq!(
        ScoreType::from_name("Posterior Error Probability").unwrap(),
        ScoreType::PosteriorErrorProbability
    );
    assert_eq!(ScoreType::from_name("score").unwrap(), ScoreType::Score);
    // Source extras: "q value", case-insensitivity, the name table.
    assert_eq!("Q VALUE".parse::<ScoreType>().unwrap(), ScoreType::QValue);
    assert_eq!(SCORE_TYPE_NAMES, ["qvalue", "PEP", "score"]);
    assert_eq!(ScoreType::PosteriorErrorProbability.name(), "PEP");
    assert_eq!(ScoreType::default(), ScoreType::QValue);
    assert!(matches!(
        ScoreType::from_name("svm"),
        Err(Error::InvalidValue(_))
    ));
}

// ---------------------------------------------------------------------------
// START_SECTION(void load(...))
// ---------------------------------------------------------------------------

/// The class test's mock-up raw data: three MS2 spectra at RT 2, 3 and 4 with
/// precursor charge `Int(rt)` and m/z `rt * 111.1`, read without native IDs.
fn class_test_lookup() -> SpectrumMetaDataLookup {
    let spectra: Vec<MSSpectrum> = [2.0, 3.0, 4.0]
        .into_iter()
        .map(|rt: f64| {
            let mut spectrum = MSSpectrum::default();
            spectrum.ms_level = 2;
            spectrum.rt = rt;
            spectrum.precursors.push(Precursor {
                charge: rt as i32,
                mz: rt * 111.1,
                ..Precursor::default()
            });
            spectrum
        })
        .collect();
    let mut lookup = SpectrumMetaDataLookup::new();
    // "no native IDs set, so don't parse them"
    assert!(lookup.read_spectra(&spectra, "", false).unwrap().is_empty());
    lookup
}

#[test]
fn load_follows_the_class_test() {
    let mut lookup = class_test_lookup();
    let document = PercolatorOutfile::new()
        .load(CLASS_TEST, &mut lookup, ScoreType::Score)
        .unwrap();
    let proteins = &document.protein_identification;
    let accessions: Vec<&str> = proteins.hits.iter().map(|h| h.accession.as_str()).collect();
    assert_eq!(accessions, ["Protein1", "Protein2", "UniProt_P01834"]);
    assert_eq!(
        proteins.search_parameters.fixed_modifications,
        ["Carbamidomethyl (C)"]
    );

    let peptides = &document.peptide_identifications;
    assert_eq!(peptides.len(), 3);
    for (index, (rt, mz, charge)) in [(2.0, 222.2, 2), (3.0, 333.3, 3), (4.0, 444.4, 4)]
        .into_iter()
        .enumerate()
    {
        similar(peptides[index].rt.unwrap(), rt);
        similar(peptides[index].mz.unwrap(), mz);
        assert_eq!(peptides[index].hits.len(), 1);
        assert_eq!(peptides[index].hits[0].charge, charge);
    }
    similar(peptides[0].hits[0].score, 6.77991);
    similar(peptides[1].hits[0].score, 6.57945);
    similar(peptides[2].hits[0].score, 6.50586);
    assert_eq!(
        peptides[0].hits[0].sequence.to_string(),
        "VDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSK"
    );
    assert_eq!(
        peptides[1].hits[0].sequence.to_string(),
        "VDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSK"
    );
    assert_eq!(
        peptides[2].hits[0].sequence.to_string(),
        "VTLSC(Carbamidomethyl)TGSSSNLGAGYDVHWYQQLPGTAPK"
    );
    similar(
        peptides[0].hits[0].metadata["Percolator_score"]
            .as_f64()
            .unwrap(),
        6.77991,
    );
    assert_eq!(
        peptides[1].hits[0].metadata["Percolator_qvalue"]
            .as_f64()
            .unwrap(),
        0.0
    );
    similar(
        peptides[2].hits[0].metadata["Percolator_PEP"]
            .as_f64()
            .unwrap(),
        1.8014e-14,
    );

    // Source detail beyond the class test: the raw score is "higher is
    // better", the run is Percolator's, and `load` registered the defaults.
    assert_eq!(peptides[0].score_type, "Percolator_score");
    assert!(peptides[0].higher_score_better);
    assert_eq!(peptides[0].identifier, "id");
    assert_eq!(proteins.identifier, "id");
    assert_eq!(proteins.search_engine, "Percolator");
    assert!(proteins.date_time.is_some());
    assert_eq!(lookup.reference_formats().len(), 3);
    for (format, pattern) in lookup
        .reference_formats()
        .iter()
        .zip(DEFAULT_REFERENCE_FORMATS)
    {
        assert_eq!(format.as_str(), pattern);
    }
    assert_eq!(
        document.log,
        ["Created 3 protein hits.\nCreated 3 peptide hits (PSMs)."]
    );
}

// ---------------------------------------------------------------------------
// Tier 1: TOPP_IDFileConverter_12, _13 and _14 (retained C++ outputs)
// ---------------------------------------------------------------------------

fn accession_set(identification: &PeptideIdentification) -> BTreeSet<String> {
    identification.hits[0]
        .evidences
        .iter()
        .map(|e| e.protein_accession.clone())
        .collect()
}

/// `IDFileConverter -in IDFileConverter_12_input.psms -out_type idXML
/// -score_type <name>`: no `mz_file`, so the lookup holds no spectra and
/// every retention time, m/z and charge comes from the Mascot-Percolator
/// references themselves. The retained output is compared record by record;
/// the upstream comparison whitelists only `IdentificationRun date`.
fn converts_like_the_retained_output(score_type: &str, expected: &str) {
    let mut lookup = SpectrumMetaDataLookup::new();
    let ours = PercolatorOutfile::new()
        .load(
            TOPP_INPUT,
            &mut lookup,
            ScoreType::from_name(score_type).unwrap(),
        )
        .unwrap();
    let reference = idxml::load(expected).unwrap();

    assert_eq!(reference.protein_identifications.len(), 1);
    let (proteins, theirs) = (
        &ours.protein_identification,
        &reference.protein_identifications[0],
    );
    assert_eq!(proteins.search_engine, theirs.search_engine);
    assert_eq!(proteins.search_engine_version, theirs.search_engine_version);
    assert_eq!(proteins.score_type, theirs.score_type);
    assert_eq!(proteins.higher_score_better, theirs.higher_score_better);
    assert_eq!(
        proteins.search_parameters.fixed_modifications,
        theirs.search_parameters.fixed_modifications
    );
    assert_eq!(
        proteins.search_parameters.variable_modifications,
        theirs.search_parameters.variable_modifications
    );
    assert_eq!(
        proteins.search_parameters.digestion_enzyme,
        theirs.search_parameters.digestion_enzyme
    );
    assert_eq!(proteins.hits.len(), theirs.hits.len());
    for (ours, theirs) in proteins.hits.iter().zip(&theirs.hits) {
        assert_eq!(ours.accession, theirs.accession);
        assert_eq!(ours.score, theirs.score);
        assert_eq!(ours.sequence, theirs.sequence);
    }

    let (peptides, theirs) = (
        &ours.peptide_identifications,
        &reference.peptide_identifications,
    );
    assert_eq!(peptides.len(), 121);
    assert_eq!(peptides.len(), theirs.len());
    for (index, (ours, theirs)) in peptides.iter().zip(theirs).enumerate() {
        assert_eq!(ours.score_type, theirs.score_type, "id {index}");
        assert_eq!(
            ours.higher_score_better, theirs.higher_score_better,
            "id {index}"
        );
        assert_eq!(ours.rt, theirs.rt, "id {index}");
        assert_eq!(ours.mz, theirs.mz, "id {index}");
        assert_eq!(ours.hits.len(), 1, "id {index}");
        assert_eq!(theirs.hits.len(), 1, "id {index}");
        let (hit, expected) = (&ours.hits[0], &theirs.hits[0]);
        assert_eq!(hit.score, expected.score, "id {index}");
        assert_eq!(hit.charge, expected.charge, "id {index}");
        assert_eq!(
            hit.sequence.to_string(),
            expected.sequence.to_string(),
            "id {index}"
        );
        assert_eq!(hit.evidences.len(), expected.evidences.len(), "id {index}");
        assert_eq!(accession_set(ours), accession_set(theirs), "id {index}");
        for key in ["Percolator_score", "Percolator_qvalue", "Percolator_PEP"] {
            assert_eq!(
                hit.metadata[key].as_f64().unwrap(),
                expected.metadata[key].as_f64().unwrap(),
                "id {index} {key}"
            );
        }
        assert_eq!(hit.metadata.len(), expected.metadata.len(), "id {index}");
    }
    assert!(
        ours.log.is_empty() || ours.log.iter().all(|line| !line.starts_with("Error")),
        "{:?}",
        ours.log
    );
}

#[test]
fn topp_idfileconverter_12_qvalue() {
    converts_like_the_retained_output(
        "qvalue",
        "tests/data/percolator_outfile/IDFileConverter_12_output.idXML",
    );
}

#[test]
fn topp_idfileconverter_13_pep() {
    converts_like_the_retained_output(
        "PEP",
        "tests/data/percolator_outfile/IDFileConverter_13_output.idXML",
    );
}

#[test]
fn topp_idfileconverter_14_score() {
    converts_like_the_retained_output(
        "score",
        "tests/data/percolator_outfile/IDFileConverter_14_output.idXML",
    );
}

// ---------------------------------------------------------------------------
// The spectrum reference path
// ---------------------------------------------------------------------------

/// A reference whose look-up fails keeps what the reference itself supplied
/// before the failure: the MS-GF+ format extracts the charge, then needs the
/// spectrum for RT and m/z, and the empty lookup has none. The source catches
/// the `ElementNotFound`, logs it, and still gives the hit that charge.
#[test]
fn an_unresolvable_reference_keeps_what_it_supplied() {
    let mut lookup = SpectrumMetaDataLookup::new();
    let document = read(&row("sample_SII_7_1_1_3_1", "K.PEPTIDE.A"), &mut lookup).unwrap();
    let peptide = &document.peptide_identifications[0];
    assert_eq!(peptide.hits[0].charge, 3);
    assert_eq!(peptide.rt, None);
    assert_eq!(peptide.mz, None);
    assert_eq!(
        document.log,
        [
            "Error: Could not extract data for spectrum reference 'sample_SII_7_1_1_3_1' from row 1",
            "Created 1 protein hits.\nCreated 1 peptide hits (PSMs).",
            "1 peptide hits without retention time information.",
            "1 peptide hits without mass-to-charge information.",
        ]
    );
}

/// A resolved look-up replaces the whole meta data, including the charge the
/// reference itself carried (source `meta = metadata_[index]`).
#[test]
fn a_resolved_look_up_replaces_the_extracted_values() {
    let mut lookup = class_test_lookup();
    let mut meta = SpectrumMetaData::default();
    lookup
        .add_reference_format(DEFAULT_REFERENCE_FORMATS[0])
        .unwrap();
    lookup
        .meta_data_for_reference(
            "x_SII_2_1_1_9_1",
            &mut meta,
            MetaDataFlags::RT | MetaDataFlags::PRECURSOR_CHARGE,
        )
        .unwrap();
    // INDEX1 = 2 is the second spectrum: RT 3, charge 3, not the 9 of the
    // reference.
    assert_eq!(meta.rt, Some(3.0));
    assert_eq!(meta.precursor_charge, 3);
    assert_eq!(meta.ms_level, 2);

    // Everything requested is in the reference: no look-up at all.
    let mut meta = SpectrumMetaData::default();
    lookup
        .meta_data_for_reference(
            "x_SII_2000_1_1_9_1",
            &mut meta,
            MetaDataFlags::PRECURSOR_CHARGE,
        )
        .unwrap();
    assert_eq!(meta.precursor_charge, 9);
    assert_eq!(meta.rt, None);

    // No format matches: nothing extracted, no error.
    let mut meta = SpectrumMetaData::default();
    lookup
        .meta_data_for_reference("unrelated", &mut meta, MetaDataFlags::ALL)
        .unwrap();
    assert_eq!(meta, SpectrumMetaData::default());
}

/// The X! Tandem format counts from zero and the Mascot format may carry an
/// empty retention time, which then has to be looked up.
#[test]
fn the_x_tandem_and_mascot_formats() {
    let mut lookup = class_test_lookup();
    let tandem = read(&row("file_0_2_1", "K.PEPTIDE.A"), &mut lookup).unwrap();
    let peptide = &tandem.peptide_identifications[0];
    assert_eq!((peptide.rt, peptide.hits[0].charge), (Some(2.0), 2));

    let mut lookup = class_test_lookup();
    let mascot = read(
        &row(
            "query:1;rank:1;spectrum:x_scans=2_file;rt:;mz:500.25;charge:2",
            "K.PEPTIDE.A",
        ),
        &mut lookup,
    )
    .unwrap();
    // RT is empty in the reference, so the look-up by INDEX0 = 2 runs and
    // replaces everything: RT 4, m/z 444.4, charge 4.
    let peptide = &mascot.peptide_identifications[0];
    assert_eq!(peptide.rt, Some(4.0));
    similar(peptide.mz.unwrap(), 444.4);
    assert_eq!(peptide.hits[0].charge, 4);
}

#[test]
fn reference_formats_must_name_a_look_up_group() {
    let mut lookup = SpectrumMetaDataLookup::new();
    assert!(matches!(
        lookup.add_reference_format(r"charge:(?<CHARGE>\d+)"),
        Err(Error::InvalidValue(_))
    ));
    assert!(lookup.add_reference_format(r"(?<SCAN>\d+").is_err());
    assert!(lookup.reference_formats().is_empty());
    lookup.add_reference_format(r"scan=(?<SCAN>\d+)").unwrap();
    assert_eq!(lookup.reference_formats().len(), 1);
}

#[test]
fn spectra_are_indexed_by_scan_number_native_id_and_rt() {
    // Dyadic retention times, so the distance comparisons are exact.
    let spectra: Vec<MSSpectrum> = [(10.0, "scan=5"), (20.0, "scan=7"), (20.5, "other")]
        .into_iter()
        .map(|(rt, id)| {
            let mut spectrum = MSSpectrum::default();
            spectrum.rt = rt;
            spectrum.native_id = id.into();
            spectrum
        })
        .collect();
    let mut lookup = SpectrumMetaDataLookup::new();
    let messages = lookup
        .read_spectra(&spectra, DEFAULT_SCAN_REGEXP, false)
        .unwrap();
    assert_eq!(
        messages,
        [
            "Error: Could not extract scan number from spectrum native ID 'other' using regular expression '=(?<SCAN>\\d+)$'."
        ]
    );
    assert_eq!(lookup.len(), 3);
    assert_eq!(lookup.find_by_scan_number(7).unwrap(), 1);
    assert_eq!(lookup.find_by_native_id("other").unwrap(), 2);
    assert!(lookup.find_by_scan_number(6).is_err());
    assert_eq!(lookup.find_by_index(1, true).unwrap(), 0);
    assert!(lookup.find_by_index(0, true).is_err());
    assert!(lookup.find_by_index(3, false).is_err());
    // Nearest within the tolerance; below wins only when strictly closer.
    assert_eq!(lookup.rt_tolerance, 0.01);
    assert!(lookup.find_by_rt(10.5).is_err());
    lookup.rt_tolerance = 1.0;
    assert_eq!(lookup.find_by_rt(10.5).unwrap(), 0);
    assert_eq!(lookup.find_by_rt(20.25).unwrap(), 2);
    assert_eq!(lookup.find_by_rt(20.125).unwrap(), 1);
    assert_eq!(lookup.find_by_rt(20.5).unwrap(), 2);
    assert!(lookup.find_by_rt(15.0).is_err());
    assert_eq!(lookup.meta_data(1).unwrap().scan_number, Some(7));
    assert!(lookup.meta_data(3).is_err());
    assert!(matches!(
        lookup.read_spectra(&spectra, r"scan=(\d+)", false),
        Err(Error::InvalidValue(_))
    ));
    // A refused read leaves the lookup as it was.
    assert_eq!(lookup.len(), 3);
}

#[test]
fn precursor_retention_times_come_from_the_previous_lower_level() {
    let spectra: Vec<MSSpectrum> = [(1.0, 1), (2.0, 2), (3.0, 2), (4.0, 1), (5.0, 2)]
        .into_iter()
        .map(|(rt, level)| {
            let mut spectrum = MSSpectrum::default();
            spectrum.rt = rt;
            spectrum.ms_level = level;
            if level == 2 {
                spectrum.precursors.push(Precursor::default());
            }
            spectrum
        })
        .collect();
    let mut lookup = SpectrumMetaDataLookup::new();
    assert!(lookup.read_spectra(&spectra, "", true).unwrap().is_empty());
    let precursor_rts: Vec<Option<f64>> = (0..5)
        .map(|i| lookup.meta_data(i).unwrap().precursor_rt)
        .collect();
    assert_eq!(precursor_rts, [None, Some(1.0), Some(1.0), None, Some(4.0)]);
}

// ---------------------------------------------------------------------------
// Peptide rewriting
// ---------------------------------------------------------------------------

fn sequence_of(peptide: &str) -> String {
    let mut lookup = SpectrumMetaDataLookup::new();
    read(&row("psm", peptide), &mut lookup)
        .unwrap()
        .peptide_identifications[0]
        .hits[0]
        .sequence
        .to_string()
}

#[test]
fn flanking_residues_unimod_and_unknown_modifications() {
    assert_eq!(sequence_of("K.PEPTIDE.A"), "PEPTIDE");
    assert_eq!(sequence_of("PEPTIDE"), "PEPTIDE");
    assert_eq!(sequence_of("-.PEPTIDEK.-"), "PEPTIDEK");
    assert_eq!(
        sequence_of("R.PEPC[UNIMOD:4]TIDEM[UNIMOD:35]K.L"),
        "PEPC(Carbamidomethyl)TIDEM(Oxidation)K"
    );

    let mut lookup = SpectrumMetaDataLookup::new();
    let document = read(&row("psm", "K.PEP[unknown]TIDE.A"), &mut lookup).unwrap();
    assert_eq!(
        document.peptide_identifications[0].hits[0]
            .sequence
            .to_string(),
        "PEPTIDE"
    );
    assert_eq!(
        document.log[0],
        "Removing unknown modification(s) from peptide 'PEP[unknown]TIDE'"
    );
}

/// X! Tandem N-terminal mass shifts reported on the first residue. Masses are
/// UniMod's: acetylation 42.010565, carbamidomethylation 57.021464,
/// phosphorylation 79.966331, oxidation 15.994915.
#[test]
fn misassigned_n_terminal_shifts_move_to_the_terminus() {
    // One shift that can be N-terminal becomes terminal.
    assert_eq!(sequence_of("K.S[42.0106]PEPTIDE.A"), ".(Acetyl)SPEPTIDE");
    // ... unless it is carbamidomethylation on C, which stays on the residue.
    assert_eq!(
        sequence_of("K.C[57.0215]PEPTIDE.A"),
        "C(Carbamidomethyl)PEPTIDE"
    );
    // Two shifts, only the first can be N-terminal.
    assert_eq!(
        sequence_of("K.S[42.0106][79.9663]PEPTIDE.A"),
        ".(Acetyl)S(Phospho)PEPTIDE"
    );
    // Two shifts, only the second can be N-terminal.
    assert_eq!(
        sequence_of("K.S[79.9663][42.0106]PEPTIDE.A"),
        ".(Acetyl)S(Phospho)PEPTIDE"
    );
    // A shift no registry entry explains within 0.01 Da is left alone and
    // gains its missing "+".
    assert_eq!(sequence_of("K.S[1.5]PEPTIDE.A"), "S[+1.5]PEPTIDE");
}

// ---------------------------------------------------------------------------
// Refusals
// ---------------------------------------------------------------------------

#[test]
fn malformed_files_are_refused() {
    let mut lookup = SpectrumMetaDataLookup::new();
    // Wrong header (source: Exception::ParseError).
    assert!(matches!(
        read("PSMId\tscore\n", &mut lookup),
        Err(Error::InvalidValue(_))
    ));
    // Empty file.
    assert!(read("", &mut lookup).is_err());
    // Only comments before the header count as nothing.
    let document = read(&format!("# comment\n{HEADER}\n"), &mut lookup).unwrap();
    assert!(document.peptide_identifications.is_empty());
    // A row with fewer than five fields (source reads out of bounds).
    assert!(read(&format!("{HEADER}\npsm\t1\t2\t3\n"), &mut lookup).is_err());
    assert!(read(&format!("{HEADER}\n\npsm\t1\t0\t0\tK.PEP.A\n"), &mut lookup).is_err());
    // A score that does not convert (source: ConversionError) or is not
    // finite (stricter than the source).
    assert!(read(&format!("{HEADER}\npsm\tx\t0\t0\tK.PEP.A\n"), &mut lookup).is_err());
    assert!(read(&format!("{HEADER}\npsm\tnan\t0\t0\tK.PEP.A\n"), &mut lookup).is_err());
    // A peptide shorter than two characters (source reads peptide[len - 2]).
    assert!(read(&format!("{HEADER}\npsm\t1\t0\t0\tK\n"), &mut lookup).is_err());
    // A peptide AASequence cannot parse.
    assert!(
        read(
            &format!("{HEADER}\npsm\t1\t0\t0\tK.PEP(Nonsense)T.A\n"),
            &mut lookup
        )
        .is_err()
    );
}

#[test]
fn accessions_become_sorted_unique_protein_hits() {
    let mut lookup = SpectrumMetaDataLookup::new();
    let text =
        format!("{HEADER}\npsm1\t1\t0\t0\tK.PEPTIDE.A\tb\ta\npsm2\t1\t0\t0\tK.PEPTIDE.A\ta\tC\n");
    let document = read(&text, &mut lookup).unwrap();
    let accessions: Vec<&str> = document
        .protein_identification
        .hits
        .iter()
        .map(|h| h.accession.as_str())
        .collect();
    assert_eq!(accessions, ["C", "a", "b"]);
    let evidences: Vec<&str> = document.peptide_identifications[0].hits[0]
        .evidences
        .iter()
        .map(|e| e.protein_accession.as_str())
        .collect();
    assert_eq!(evidences, ["b", "a"]);
    assert!(
        document
            .protein_identification
            .search_parameters
            .fixed_modifications
            .is_empty()
    );
}
