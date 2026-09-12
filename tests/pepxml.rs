// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "idxml")]

//! Ported from `src/tests/class_tests/openms/source/PepXMLFile_test.cpp` at
//! revision bc9cc12514c768385ce121d6ca4bb710fe1983c4. Every START_SECTION of
//! that file has a counterpart here; see `docs/PEPXML_SUPPORT.md` for the
//! section-by-section mapping and `tests/data/pepxml_provenance.json` for the
//! hashes of the source and of every fixture.

use openms::Error;
use openms::format::pepxml::{
    self, PepXmlDocument, ReadOptions, SpectrumIndex, SpectrumMetaData, WriteOptions,
};
use openms::identification::{FlankingResidue, PeakMassType};

const TEST: &str = include_str!("data/PepXMLFile_test.pepxml");
const EXTENDED: &str = include_str!("data/PepXMLFile_test_extended.pepxml");
const STORE: &str = include_str!("data/PepXMLFile_test_store.pepxml");
const OUT: &str = include_str!("data/PepXMLFile_test_out.pepxml");
const OUT_1: &str = include_str!("data/PepXMLFile_test_out_1.pepxml");
const OUT_MZML: &str = include_str!("data/PepXMLFile_test_out_mzML.pepxml");

fn options(experiment: &str) -> ReadOptions {
    ReadOptions {
        experiment_name: experiment.into(),
        ..Default::default()
    }
}

fn load(text: &str, options: &ReadOptions) -> PepXmlDocument {
    pepxml::read_with_options(text.as_bytes(), options).unwrap()
}

/// The 13 spectra of the upstream `PepXMLFile_test.mzML`, transcribed from its
/// `<spectrum>` identifiers, `MS:1000511` ms level and `MS:1000016` scan start
/// time. The mzML itself is 1.4 MB and this port never reads spectra from the
/// pepXML writer, so only the index it needs is carried here. The file's sha256
/// is recorded in `tests/data/pepxml_provenance.json`.
fn upstream_lookup() -> SpectrumIndex {
    let rows: [(u32, f64); 13] = [
        (1, 0.5927),
        (2, 1.3653),
        (2, 1.719),
        (2, 2.0436),
        (1, 2.3364),
        (2, 2.7843),
        (2, 3.085),
        (2, 3.4018),
        (1, 3.7572),
        (2, 4.4573),
        (2, 4.8288),
        (2, 5.1785),
        (1, 5.4999),
    ];
    let entries = rows
        .iter()
        .enumerate()
        .map(|(index, &(ms_level, rt))| {
            let scan = index as u64 + 1;
            SpectrumMetaData {
                native_id: format!("scan={scan}"),
                rt,
                ms_level,
                scan_number: Some(scan),
            }
        })
        .collect();
    SpectrumIndex::from_entries(entries).unwrap()
}

// ---------------------------------------------------------------------------
// FuzzyStringComparator, as the upstream store sections use it
// ---------------------------------------------------------------------------

/// A number starting at `bytes[start]`, and the offset just past it.
fn scan_number(bytes: &[u8], start: usize) -> Option<(f64, usize)> {
    let mut index = start;
    if matches!(bytes.get(index), Some(b'+' | b'-')) {
        index += 1;
    }
    let digits = index;
    while bytes.get(index).is_some_and(u8::is_ascii_digit) {
        index += 1;
    }
    let integer = index - digits;
    if bytes.get(index) == Some(&b'.') {
        index += 1;
        while bytes.get(index).is_some_and(u8::is_ascii_digit) {
            index += 1;
        }
    }
    if index - digits == 0 || (integer == 0 && index - digits <= 1) {
        return None;
    }
    let mantissa = index;
    if matches!(bytes.get(index), Some(b'e' | b'E')) {
        let mut probe = index + 1;
        if matches!(bytes.get(probe), Some(b'+' | b'-')) {
            probe += 1;
        }
        let exponent = probe;
        while bytes.get(probe).is_some_and(u8::is_ascii_digit) {
            probe += 1;
        }
        if probe > exponent {
            index = probe;
        }
    }
    let text = std::str::from_utf8(&bytes[start..index.max(mantissa)]).ok()?;
    text.parse::<f64>().ok().map(|value| (value, index))
}

/// `FuzzyStringComparator::compareFiles` with the sections' two tolerances.
///
/// Numbers are compared numerically and everything else byte for byte. A number
/// pair is accepted when the absolute difference is within `absolute` OR the
/// ratio of the two is within `ratio`, which is exactly the source's rule: a
/// large relative error is tolerated when the absolute difference is small and
/// a large absolute difference is tolerated when the ratio is small.
fn fuzzy_equal(left: &str, right: &str, absolute: f64, ratio: f64) -> Result<(), String> {
    let (a, b) = (left.as_bytes(), right.as_bytes());
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.len() && j < b.len() {
        match (scan_number(a, i), scan_number(b, j)) {
            (Some((x, next_i)), Some((y, next_j))) => {
                let difference = (x - y).abs();
                let accepted = difference <= absolute
                    || (x != 0.0 && y != 0.0 && x.signum() == y.signum() && {
                        let quotient = x / y;
                        quotient.max(1.0 / quotient) <= ratio
                    });
                if !accepted {
                    return Err(format!("number {x} differs from {y} at byte {i}"));
                }
                i = next_i;
                j = next_j;
            }
            _ => {
                if a[i] != b[j] {
                    let context = |s: &[u8], at: usize| {
                        String::from_utf8_lossy(&s[at.saturating_sub(60)..s.len().min(at + 60)])
                            .into_owned()
                    };
                    return Err(format!(
                        "byte {i} differs\n  left:  {}\n  right: {}",
                        context(a, i),
                        context(b, j)
                    ));
                }
                i += 1;
                j += 1;
            }
        }
    }
    if a[i..].iter().any(|b| !b.is_ascii_whitespace())
        || b[j..].iter().any(|b| !b.is_ascii_whitespace())
    {
        return Err("one input has trailing content".into());
    }
    Ok(())
}

#[test]
fn fuzzy_comparator_matches_the_source_acceptance_rule() {
    // Equal within the absolute tolerance the store sections set.
    assert!(fuzzy_equal("mass=\"1.0\"", "mass=\"1.00000001\"", 1e-7, 1.0 + 1e-7).is_ok());
    // Equal within the relative tolerance, far outside the absolute one.
    assert!(
        fuzzy_equal(
            "mass=\"160.030648985200003\"",
            "mass=\"160.030648199800000\"",
            1e-7,
            1.0 + 1e-7
        )
        .is_ok()
    );
    // Outside both.
    assert!(fuzzy_equal("v=\"1.0\"", "v=\"1.1\"", 1e-7, 1.0 + 1e-7).is_err());
    // Non-numeric text is compared byte for byte.
    assert!(fuzzy_equal("name=\"a\"", "name=\"b\"", 1e-7, 1.0 + 1e-7).is_err());
}

// ---------------------------------------------------------------------------
// START_SECTION(PepXMLFile()) and START_SECTION(~PepXMLFile())
// ---------------------------------------------------------------------------

#[test]
fn default_construction_and_drop() {
    // The source constructor only caches the hydrogen element and clears the
    // parser flags; the destructor is defaulted. The Rust equivalents are the
    // default options and an owned document, which own no external resource.
    let defaults = ReadOptions::default();
    assert!(defaults.experiment_name.is_empty());
    assert!(!defaults.keep_native_spectrum_name);
    assert!(!defaults.parse_unknown_scores);
    assert!(defaults.lookup.is_empty());
    let document = PepXmlDocument::default();
    assert!(document.protein_identifications.is_empty());
    assert!(document.peptide_identifications.is_empty());
    drop(document);
    // Reading twice through the same options leaves no state behind, which is
    // what the source's "load could be called several times" reset covers.
    let first = load(TEST, &options("PepXMLFile_test"));
    let second = load(TEST, &options("PepXMLFile_test"));
    assert_eq!(first, second);
}

// ---------------------------------------------------------------------------
// START_SECTION(void load(..., const SpectrumMetaDataLookup& lookup))
// ---------------------------------------------------------------------------

#[test]
fn load_with_a_spectrum_lookup() {
    let read = ReadOptions {
        lookup: upstream_lookup(),
        ..options("PepXMLFile_test")
    };
    let document = load(TEST, &read);
    assert_eq!(document.peptide_identifications.len(), 18);
    assert_eq!(document.protein_identifications.len(), 2);
    let first = &document.peptide_identifications[0];
    assert!((first.rt.unwrap() - 1.3653).abs() < 1e-9);
    assert!((first.mz.unwrap() - 538.605).abs() < 1e-3);
    // Every query in this fixture carries retention_time_sec, so the lookup is
    // never consulted and the result is identical without it.
    assert_eq!(document, load(TEST, &options("PepXMLFile_test")));
}

// ---------------------------------------------------------------------------
// START_SECTION(void load(..., const std::string& experiment_name = ""))
// ---------------------------------------------------------------------------

#[test]
fn load_two_search_runs_of_one_experiment() {
    let document = load(TEST, &options("PepXMLFile_test"));
    let peptides = &document.peptide_identifications;
    let proteins = &document.protein_identifications;
    assert_eq!(peptides.len(), 18);
    let first = &peptides[0];
    let last = peptides.last().unwrap();

    // Identical for every peptide of the first search run.
    for peptide in &peptides[1..9] {
        assert_eq!(first.identifier, peptide.identifier);
        assert_eq!(first.score_type, peptide.score_type);
        assert_eq!(first.higher_score_better, peptide.higher_score_better);
        assert_eq!(first.significance_threshold, peptide.significance_threshold);
    }

    assert!((first.rt.unwrap() - 1.3653).abs() < 1e-9); // RT of the MS2 spectrum
    assert!((first.mz.unwrap() - 538.605).abs() < 1e-3); // recomputed
    assert_eq!(first.hits.len(), 1);
    let hit = &first.hits[0];
    assert_eq!(hit.sequence.to_string(), ".(Glu->pyro-Glu)ELNKEMAAEKAKAAAG");
    assert_eq!(hit.sequence.as_str(), "ELNKEMAAEKAKAAAG");
    assert_eq!(hit.rank, 0);
    // The score itself is not checked, because the implementation may change.
    assert_eq!(hit.charge, 3);
    assert_eq!(hit.evidences.len(), 3);
    assert_eq!(hit.evidences[0].protein_accession, "ddb000449223");
    assert_eq!(hit.evidences[0].aa_before, FlankingResidue::Residue('R'));
    assert_eq!(hit.evidences[0].aa_after, FlankingResidue::Residue('E'));

    assert!(first.hits[0].sequence.is_modified());
    assert!(first.hits[0].sequence.n_terminal_modification().is_some());
    assert!(first.hits[0].sequence.c_terminal_modification().is_none());

    assert!(peptides[1].hits[0].sequence.is_modified());
    assert!(
        peptides[1].hits[0]
            .sequence
            .n_terminal_modification()
            .is_some()
    );
    assert!(
        peptides[1].hits[0]
            .sequence
            .c_terminal_modification()
            .is_none()
    );

    assert!(peptides[5].hits[0].sequence.is_modified());
    assert!(
        peptides[5].hits[0]
            .sequence
            .n_terminal_modification()
            .is_none()
    );
    assert!(
        peptides[5].hits[0]
            .sequence
            .c_terminal_modification()
            .is_none()
    );

    // A cursory check of a peptide ID from the second search run.
    assert_eq!(last.hits[0].sequence.to_string(), "EISPDTTLLDLQNNDISELR");

    assert_eq!(proteins.len(), 2);
    assert_eq!(proteins[0].identifier, first.identifier);
    assert_eq!(proteins[1].identifier, last.identifier);
    assert!(!proteins[0].identifier.is_empty());
    assert!(!proteins[1].identifier.is_empty());
    assert_ne!(proteins[0].identifier, proteins[1].identifier);
    assert_eq!(proteins[0].search_engine, "X! Tandem (k-score)");
    assert_eq!(proteins[1].search_engine, "SEQUEST");

    assert_eq!(proteins[0].hits.len(), 20);
    let accessions: Vec<&str> = proteins[0]
        .hits
        .iter()
        .map(|hit| hit.accession.as_str())
        .collect();
    // A sample of the accessions that must be present.
    assert!(accessions.contains(&"ddb000449223"));
    assert!(accessions.contains(&"ddb000626346"));
    assert!(accessions.contains(&"rev000409159"));

    let parameters = &proteins[0].search_parameters;
    assert_eq!(parameters.database, "./current.fasta");
    assert_eq!(parameters.mass_type, PeakMassType::Monoisotopic);
    assert_eq!(parameters.digestion_enzyme, "Trypsin");

    assert_eq!(parameters.fixed_modifications.len(), 1);
    assert_eq!(parameters.variable_modifications.len(), 12);
    let variable = &parameters.variable_modifications;
    assert_eq!(variable[0], "Ammonia-loss (N-term C)");
    assert_eq!(variable[1], "Glu->pyro-Glu (N-term E)");
    assert_eq!(variable[2], "Oxidation (M)");
    assert_eq!(variable[3], "Gln->pyro-Glu (N-term Q)");
    assert_eq!(variable[4], "M[+1.0]");
    assert_eq!(variable[5], ".n[+2.0]");
    assert_eq!(variable[6], ".c[+2.0]");
    assert_eq!(variable[7], ".n[+2.5]");
    assert_eq!(variable[8], ".n[+2.5]");
    assert_eq!(variable[9], ".n[-2.5]");
    assert_eq!(variable[10], ".n[+2.5]");
    assert_eq!(variable[11], ".c[+3.4]");

    // A wrong experiment name is a parse error.
    let error = pepxml::read_with_options(TEST.as_bytes(), &options("abcxyz")).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error}");

    // A missing pepXML file is an I/O error.
    let error =
        pepxml::load("this_file_does_not_exist_but_should_be_a_pepXML_file.pepXML").unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error}");
}

#[test]
fn load_without_an_experiment_name_accepts_every_run() {
    // The source's defaulted experiment_name skips the base_name check entirely.
    let document = load(TEST, &ReadOptions::default());
    assert_eq!(document.peptide_identifications.len(), 18);
    assert_eq!(document.protein_identifications.len(), 2);
}

// ---------------------------------------------------------------------------
// START_SECTION([EXTRA] void load(..., const std::string& experiment_name = ""))
// ---------------------------------------------------------------------------

#[test]
fn load_extended_fixture_with_native_names_and_analysis_results() {
    let read = ReadOptions {
        keep_native_spectrum_name: true,
        ..options("PepXMLFile_test")
    };
    let document = load(EXTENDED, &read);
    let peptides = &document.peptide_identifications;
    assert_eq!(peptides.len(), 2);
    let first = &peptides[0];
    let last = peptides.last().unwrap();

    assert_eq!(first.rt, Some(1.3653)); // RT of the MS2 spectrum
    assert!((first.mz.unwrap() - 538.605).abs() < 1e-3); // recomputed
    assert_eq!(first.hits.len(), 1);

    assert_eq!(last.rt, Some(488.652)); // RT of the MS2 spectrum
    assert!((last.mz.unwrap() - 585.316_625_031_9).abs() < 1e-9); // recomputed
    assert_eq!(last.hits.len(), 1);
    assert!(last.metadata.contains_key("swath_assay"));
    assert!(last.metadata.contains_key("status"));
    assert!(last.metadata.contains_key("pepxml_spectrum_name"));
    assert!(!last.experiment_label().is_empty());

    assert_eq!(last.metadata["swath_assay"].to_string(), "EIVLTQSPGTL2:9");
    assert_eq!(last.metadata["status"].to_string(), "target");
    assert_eq!(
        last.metadata["pepxml_spectrum_name"].to_string(),
        "hroest_K120718_SM_OGE10_010_IDA.02552.02552.2"
    );
    assert_eq!(last.experiment_label(), "urine");

    let hit = &last.hits[0];
    assert_eq!(hit.sequence.to_string(), "VVITAPGGNDVK");
    assert_eq!(hit.sequence.as_str(), "VVITAPGGNDVK");
    assert_eq!(hit.rank, 0);
    assert_eq!(hit.charge, 2);

    assert_eq!(hit.analysis_results.len(), 2);
    let peptideprophet = &hit.analysis_results[0];
    assert_eq!(peptideprophet.score_type, "peptideprophet");
    assert!((peptideprophet.main_score - 0.0660).abs() < 1e-9);
    for name in ["fval", "ntt", "empir_irt", "swath_window"] {
        assert!(peptideprophet.sub_scores.contains_key(name), "{name}");
    }
    assert!((peptideprophet.sub_scores["fval"] - 0.7114).abs() < 1e-9);
    assert!((peptideprophet.sub_scores["ntt"] - 2.0).abs() < 1e-9);
    assert!((peptideprophet.sub_scores["empir_irt"] - 79.79).abs() < 1e-9);
    assert!((peptideprophet.sub_scores["swath_window"] - 9.0).abs() < 1e-9);

    let interprophet = &hit.analysis_results[1];
    assert_eq!(interprophet.score_type, "interprophet");
    assert!((interprophet.main_score - 0.93814).abs() < 1e-9);
    assert!(!interprophet.sub_scores.contains_key("fval"));
    assert!(interprophet.sub_scores.contains_key("nss"));
    assert!((interprophet.sub_scores["nrs"] - 10.2137).abs() < 1e-9);

    // InterProphet overwrites PeptideProphet, which overwrote the search score.
    assert_eq!(last.score_type, "InterProphet probability");
    assert!(last.higher_score_better);
    assert!((hit.score - 0.93814).abs() < 1e-9);

    let error = pepxml::read_with_options(EXTENDED.as_bytes(), &options("abcxyz")).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error}");
    let error = pepxml::load_with_options(
        "this_file_does_not_exist_but_should_be_a_pepXML_file.pepXML",
        &read,
    )
    .unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error}");
}

// ---------------------------------------------------------------------------
// START_SECTION(void store(...)) - retained upstream output
// ---------------------------------------------------------------------------

fn store_text(document: &PepXmlDocument, options: &WriteOptions) -> String {
    let mut bytes = Vec::new();
    pepxml::write_with_options(&mut bytes, document, options).unwrap();
    String::from_utf8(bytes).unwrap()
}

#[test]
fn store_reproduces_the_retained_peptideprophet_output() {
    let document = load(STORE, &ReadOptions::default());
    assert_eq!(document.peptide_identifications.len(), 9);
    let written = store_text(
        &document,
        &WriteOptions {
            mz_name: "test".into(),
            peptideprophet_analyzed: true,
            ..Default::default()
        },
    );
    if let Err(message) = fuzzy_equal(&written, OUT, 1e-7, 1.0 + 1e-7) {
        panic!("{message}");
    }
}

#[test]
fn store_reproduces_the_retained_raw_output() {
    let document = load(STORE, &ReadOptions::default());
    let written = store_text(
        &document,
        &WriteOptions {
            mz_name: "test".into(),
            peptideprophet_analyzed: false,
            ..Default::default()
        },
    );
    if let Err(message) = fuzzy_equal(&written, OUT_1, 1e-7, 1.0 + 1e-7) {
        panic!("{message}");
    }
}

// ---------------------------------------------------------------------------
// START_SECTION(void store(..., mz_file = "PepXMLFile_test.mzML", ...))
// ---------------------------------------------------------------------------

#[test]
fn store_with_spectra_metadata_reproduces_the_retained_output() {
    let document = load(STORE, &ReadOptions::default());
    let written = store_text(
        &document,
        &WriteOptions {
            // The source loads this file to build its lookup; this port takes
            // the lookup directly and uses the name only for base_name/raw_data.
            mz_file: "PepXMLFile_test.mzML".into(),
            mz_name: "test".into(),
            peptideprophet_analyzed: true,
            lookup: upstream_lookup(),
            ..Default::default()
        },
    );
    if let Err(message) = fuzzy_equal(&written, OUT_MZML, 1e-7, 1.0 + 1e-7) {
        panic!("{message}");
    }
}

// ---------------------------------------------------------------------------
// START_SECTION([EXTRA] void store(...)) - store/load round trip
// ---------------------------------------------------------------------------

#[test]
fn store_and_reload_the_extended_fixture() {
    let read = ReadOptions {
        keep_native_spectrum_name: true,
        ..options("PepXMLFile_test")
    };
    let document = load(EXTENDED, &read);
    assert_eq!(document.peptide_identifications.len(), 2);
    let first = &document.peptide_identifications[0];
    let last = document.peptide_identifications.last().unwrap();
    assert!((first.mz.unwrap() - 538.605).abs() < 1e-3); // recomputed
    assert!((last.mz.unwrap() - 585.316_625_031_9).abs() < 1e-9); // recomputed

    // peptideprophet_analyzed = false is important here.
    let written = store_text(
        &document,
        &WriteOptions {
            mz_name: "PepXMLFile_test".into(),
            keep_native_spectrum_name: true,
            peptideprophet_analyzed: false,
            ..Default::default()
        },
    );
    let reread = load(&written, &read);
    assert_eq!(
        document.protein_identifications.len(),
        reread.protein_identifications.len()
    );
    assert_eq!(
        document.peptide_identifications.len(),
        reread.peptide_identifications.len()
    );
    assert_eq!(reread.peptide_identifications.len(), 2);
    let first = &reread.peptide_identifications[0];
    let last = reread.peptide_identifications.last().unwrap();

    assert_eq!(first.rt, Some(1.3653)); // RT of the MS2 spectrum
    assert!((first.mz.unwrap() - 538.615_924_863_3).abs() < 1e-6); // recomputed
    assert_eq!(first.hits.len(), 1);

    assert_eq!(last.rt, Some(488.652)); // RT of the MS2 spectrum
    assert!((last.mz.unwrap() - 585.330_421_935_5).abs() < 1e-6); // recomputed
    assert_eq!(last.hits.len(), 1);
    let hit = &last.hits[0];
    assert_eq!(hit.sequence.to_string(), "VVITAPGGNDVK");
    assert_eq!(hit.sequence.as_str(), "VVITAPGGNDVK");
    assert_eq!(hit.rank, 0);
    assert_eq!(hit.charge, 2);

    // The extra attributes survive a store/load round trip.
    assert!(last.metadata.contains_key("swath_assay"));
    assert!(last.metadata.contains_key("status"));
    assert!(last.metadata.contains_key("pepxml_spectrum_name"));
    assert!(!last.experiment_label().is_empty());
    assert_eq!(last.metadata["swath_assay"].to_string(), "EIVLTQSPGTL2:9");
    assert_eq!(last.metadata["status"].to_string(), "target");
    // The writer appends the charge to the retained spectrum name, so the
    // reloaded value carries one more digit than the original.
    assert_eq!(
        last.metadata["pepxml_spectrum_name"].to_string(),
        "hroest_K120718_SM_OGE10_010_IDA.02552.02552.22"
    );
    assert!(
        last.metadata["pepxml_spectrum_name"].to_string()
            == "hroest_K120718_SM_OGE10_010_IDA.02552.02552.22"
    );
    assert_eq!(last.experiment_label(), "urine");

    assert_eq!(hit.analysis_results.len(), 2);
    let peptideprophet = &hit.analysis_results[0];
    assert_eq!(peptideprophet.score_type, "peptideprophet");
    assert!((peptideprophet.main_score - 0.0660).abs() < 1e-9);
    for name in ["fval", "ntt", "empir_irt", "swath_window"] {
        assert!(peptideprophet.sub_scores.contains_key(name), "{name}");
    }
    assert!((peptideprophet.sub_scores["fval"] - 0.7114).abs() < 1e-9);
    assert!((peptideprophet.sub_scores["ntt"] - 2.0).abs() < 1e-9);
    assert!((peptideprophet.sub_scores["empir_irt"] - 79.79).abs() < 1e-9);
    assert!((peptideprophet.sub_scores["swath_window"] - 9.0).abs() < 1e-9);
}

// ---------------------------------------------------------------------------
// START_SECTION(void keepNativeSpectrumName(bool keep))
// ---------------------------------------------------------------------------

#[test]
fn native_spectrum_names_are_only_kept_on_request() {
    // The upstream section is NOT_TESTABLE and points at the store sections,
    // which store and load once with the flag and once without.
    let read = ReadOptions {
        keep_native_spectrum_name: false,
        ..options("PepXMLFile_test")
    };
    let document = load(EXTENDED, &read);
    let last = document.peptide_identifications.last().unwrap();
    assert!(!last.metadata.contains_key("pepxml_spectrum_name"));
    let written = store_text(
        &document,
        &WriteOptions {
            mz_name: "PepXMLFile_test".into(),
            keep_native_spectrum_name: false,
            peptideprophet_analyzed: false,
            ..Default::default()
        },
    );
    let reread = load(&written, &read);
    let last = reread.peptide_identifications.last().unwrap();
    assert!(
        last.metadata
            .get("pepxml_spectrum_name")
            .map(ToString::to_string)
            .as_deref()
            != Some("hroest_K120718_SM_OGE10_010_IDA.02552.02552.22")
    );
    // Without the flag the writer falls back to the reconstructed name.
    assert_eq!(
        last.spectrum_reference(),
        "scan=2",
        "the second identification is the second record"
    );
}

// ---------------------------------------------------------------------------
// START_SECTION(([EXTRA] checking pepxml transformation to reusable
// identifications))
// ---------------------------------------------------------------------------

#[test]
fn pepxml_search_parameters_survive_an_idxml_round_trip() {
    use openms::format::idxml::{self, IdXmlDocument};
    let document = load(STORE, &ReadOptions::default());
    let transfer = IdXmlDocument {
        document_id: String::new(),
        protein_identifications: document.protein_identifications.clone(),
        peptide_identifications: document.peptide_identifications.clone(),
        unreferenced_search_parameters: Vec::new(),
    };
    let mut bytes = Vec::new();
    idxml::write(&mut bytes, &transfer).unwrap();
    let reread = idxml::read(bytes.as_slice()).unwrap();

    let parameters = &document.protein_identifications[0].search_parameters;
    let reread_parameters = &reread.protein_identifications[0].search_parameters;
    assert_eq!(parameters.database, reread_parameters.database);
    assert_eq!(parameters.mass_type, reread_parameters.mass_type);
    assert_eq!(
        parameters.fixed_modifications.len(),
        reread_parameters.fixed_modifications.len()
    );
    assert_eq!(
        parameters.variable_modifications.len(),
        reread_parameters.variable_modifications.len()
    );
    assert!(
        parameters
            .fixed_modifications
            .contains(&reread_parameters.fixed_modifications[0])
    );
    assert!(
        parameters
            .fixed_modifications
            .contains(&"Carbamidomethyl (C)".to_owned())
    );
    for modification in &reread_parameters.variable_modifications {
        assert!(
            parameters.variable_modifications.contains(modification),
            "{modification}"
        );
    }
}

// ---------------------------------------------------------------------------
// Native boundaries: no upstream section covers these
// ---------------------------------------------------------------------------

#[test]
fn non_ascii_paths_and_payloads_do_not_panic() {
    // A recent audit found a reachable panic from byte-slicing a path; every
    // string this module inspects is either constructed here or indexed by char.
    let error = pepxml::load("dir/日本語.pepxml").unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error}");
    let document = PepXmlDocument::default();
    let written = store_text(
        &document,
        &WriteOptions {
            mz_file: "データ/日本語.mzML".into(),
            output_name: "日本語".into(),
            ..Default::default()
        },
    );
    assert!(written.contains("base_name=\"日本語\""));
    // A pepXML whose text, accessions and modification descriptions are
    // non-ASCII still parses, and nothing is sliced at a byte boundary.
    let text = TEST.replace("ddb000449223", "日本語プロテイン");
    let reread = load(&text, &options("PepXMLFile_test"));
    assert_eq!(reread.peptide_identifications.len(), 18);
    assert_eq!(
        reread.peptide_identifications[0].hits[0].evidences[0].protein_accession,
        "日本語プロテイン"
    );
    // A multi-byte experiment name that no run matches is a parse error, not a
    // slice panic.
    let error = pepxml::read_with_options(TEST.as_bytes(), &options("日本語")).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error}");
}

#[test]
fn one_based_indices_of_zero_are_explicit_errors() {
    for (from, to) in [
        (
            "hit_rank=\"1\" peptide=\"ELNKEMAAEKAKAAAG\"",
            "hit_rank=\"0\" peptide=\"ELNKEMAAEKAKAAAG\"",
        ),
        (
            "<mod_aminoacid_mass position=\"1\" mass=\"111.0320\" />",
            "<mod_aminoacid_mass position=\"0\" mass=\"111.0320\" />",
        ),
        (
            "<mod_aminoacid_mass position=\"1\" mass=\"111.0320\" />",
            "<mod_aminoacid_mass position=\"99\" mass=\"111.0320\" />",
        ),
        ("<search_result>", "<search_result search_id=\"0\">"),
        ("assumed_charge=\"3\"", "assumed_charge=\"0\""),
    ] {
        let text = TEST.replacen(from, to, 1);
        assert_ne!(text, TEST, "fixture pattern {from} not found");
        let error =
            pepxml::read_with_options(text.as_bytes(), &options("PepXMLFile_test")).unwrap_err();
        assert!(matches!(error, Error::Parse { .. }), "{to}: {error}");
    }
}

#[test]
fn resource_ceilings_are_checked_before_allocation() {
    for read in [
        ReadOptions {
            max_input_bytes: 1024,
            ..options("PepXMLFile_test")
        },
        ReadOptions {
            max_elements: 10,
            ..options("PepXMLFile_test")
        },
        ReadOptions {
            max_depth: 2,
            ..options("PepXMLFile_test")
        },
        ReadOptions {
            max_identifications: 3,
            ..options("PepXMLFile_test")
        },
        ReadOptions {
            max_protein_hits: 3,
            ..options("PepXMLFile_test")
        },
        ReadOptions {
            max_modifications: 2,
            ..options("PepXMLFile_test")
        },
    ] {
        let error = pepxml::read_with_options(TEST.as_bytes(), &read).unwrap_err();
        assert!(matches!(error, Error::Parse { .. }), "{error}");
    }
    // Zero-valued ceilings are rejected before any input is read.
    let error = pepxml::read_with_options(
        TEST.as_bytes(),
        &ReadOptions {
            max_depth: 0,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error}");
    // The output ceiling is checked while serializing.
    let document = load(STORE, &ReadOptions::default());
    let mut bytes = Vec::new();
    let error = pepxml::write_with_options(
        &mut bytes,
        &document,
        &WriteOptions {
            max_output_bytes: 512,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error}");
}

#[test]
fn malformed_documents_are_rejected_rather_than_partially_loaded() {
    for text in [
        "",
        "<msms_pipeline_analysis/>",
        "<msms_pipeline_analysis date=\"2009-05-25T12:33:22\">",
        "<a/><b/>",
        "<!DOCTYPE x><msms_pipeline_analysis date=\"2009-05-25T12:33:22\"/>",
        "text<msms_pipeline_analysis date=\"2009-05-25T12:33:22\"/>",
    ] {
        let error = pepxml::read(text.as_bytes()).unwrap_err();
        assert!(
            matches!(error, Error::Parse { .. } | Error::Unsupported(_)),
            "{text:?}: {error}"
        );
    }
}

#[test]
fn a_corrupted_date_is_repaired_and_reported() {
    let text = TEST.replace(
        "date=\"2009-05-25T12:33:22\"",
        "date=\"2009:05:25:12:33:22\"",
    );
    let document = load(&text, &options("PepXMLFile_test"));
    assert!(
        document
            .warnings
            .iter()
            .any(|warning| warning.contains("xs:dateTime")),
        "{:?}",
        document.warnings
    );
    assert!(
        document.protein_identifications[0]
            .identifier
            .starts_with("X! Tandem (k-score)_2009-05-25_12:33:22")
    );
    // A date too short for the repair's byte indices is not read out of bounds.
    let text = TEST.replace("date=\"2009-05-25T12:33:22\"", "date=\"x\"");
    let document = load(&text, &options("PepXMLFile_test"));
    assert_eq!(document.peptide_identifications.len(), 18);
    // A multi-byte date whose bytes 4, 7 and 10 all are ':' satisfies the
    // repair's guard, so the three replacements must land on char boundaries.
    // "\u{65e5}x:ab:cd:" occupies bytes 0..2 (the kanji), 3, 4, 5, 6, 7, 8, 9, 10.
    for date in [
        "日x:ab:cd:ef",
        "日本語:05:25:12:33:22",
        "2009-05-25T12:33:22日",
    ] {
        let text = TEST.replace("date=\"2009-05-25T12:33:22\"", &format!("date=\"{date}\""));
        let document = load(&text, &options("PepXMLFile_test"));
        assert_eq!(document.peptide_identifications.len(), 18, "{date}");
    }
}

#[test]
fn unknown_scores_are_recorded_only_on_request() {
    let base = options("PepXMLFile_test");
    let document = load(TEST, &base);
    assert!(
        !document.peptide_identifications[0].hits[0]
            .metadata
            .contains_key("bscore")
    );
    let read = ReadOptions {
        parse_unknown_scores: true,
        ..options("PepXMLFile_test")
    };
    let document = load(TEST, &read);
    let hit = &document.peptide_identifications[0].hits[0];
    assert!((hit.metadata["bscore"].as_f64().unwrap() - 1.0).abs() < 1e-12);
    assert!((hit.metadata["yscore"].as_f64().unwrap() - 1.0).abs() < 1e-12);
}

#[test]
fn the_spectrum_lookup_resolves_by_id_scan_number_and_retention_time() {
    let lookup = upstream_lookup();
    assert_eq!(lookup.len(), 13);
    assert!(!lookup.is_empty());
    assert_eq!(lookup.find_by_native_id("scan=2"), Some(1));
    assert_eq!(lookup.find_by_native_id("scan=99"), None);
    assert_eq!(lookup.find_by_scan_number(12), Some(11));
    assert_eq!(lookup.find_by_scan_number(99), None);
    assert!((lookup.rt_tolerance() - 0.01).abs() < 1e-12);
    assert_eq!(lookup.find_by_rt(1.3653), Some(1));
    assert_eq!(lookup.find_by_rt(1.37), Some(1));
    assert_eq!(lookup.find_by_rt(1.5), None);
    assert_eq!(lookup.find_by_rt(f64::NAN), None);
    assert_eq!(lookup.get(1).unwrap().ms_level, 2);
    assert!(lookup.get(13).is_none());
    let mut lookup = lookup;
    assert!(lookup.set_rt_tolerance(-1.0).is_err());
    assert!(lookup.set_rt_tolerance(1.0).is_ok());
    assert_eq!(lookup.find_by_rt(1.5), Some(1));
    assert!(SpectrumIndex::is_native_id("scan=1"));
    assert!(SpectrumIndex::is_native_id(
        "controllerType=0 controllerNumber=1 scan=1"
    ));
    assert!(!SpectrumIndex::is_native_id(
        "PepXMLFile_test.00002.00002.3"
    ));
    assert!(SpectrumIndex::new().is_empty());
    assert!(
        SpectrumIndex::from_entries(vec![SpectrumMetaData {
            rt: f64::NAN,
            ..Default::default()
        }])
        .is_err()
    );
}

#[test]
fn a_missing_retention_time_is_taken_from_the_lookup() {
    let text = TEST.replace(" retention_time_sec=\"1.3653\"", "");
    // Without spectra the source reports a non-fatal error and leaves the RT.
    let document = load(&text, &options("PepXMLFile_test"));
    assert_eq!(document.peptide_identifications[0].rt, None);
    assert!(
        document
            .warnings
            .iter()
            .any(|warning| warning.contains("no spectra given")),
        "{:?}",
        document.warnings
    );
    let read = ReadOptions {
        lookup: upstream_lookup(),
        ..options("PepXMLFile_test")
    };
    let document = load(&text, &read);
    assert_eq!(document.peptide_identifications[0].rt, Some(1.3653));
}

#[test]
fn decoy_prefixes_annotate_both_peptide_and_protein_hits() {
    use openms::identification::TargetDecoyType;
    let text = TEST.replace(
        "<search_database local_path=\"./current.fasta\" type=\"AA\"/>",
        "<search_database local_path=\"./current.fasta\" type=\"AA\"/>\n\
         <parameter name=\"decoy_search\" value=\"1\"/>\n\
         <parameter name=\"decoy_prefix\" value=\"rev\"/>",
    );
    let document = load(&text, &options("PepXMLFile_test"));
    let hits = &document.protein_identifications[0].hits;
    let decoy = hits
        .iter()
        .find(|hit| hit.accession == "rev000409159")
        .unwrap();
    assert_eq!(decoy.target_decoy_type().unwrap(), TargetDecoyType::Decoy);
    let target = hits
        .iter()
        .find(|hit| hit.accession == "ddb000449223")
        .unwrap();
    assert_eq!(target.target_decoy_type().unwrap(), TargetDecoyType::Target);
    let peptide = document
        .peptide_identifications
        .iter()
        .find(|peptide| peptide.hits[0].sequence.as_str() == "SSLKNYANK")
        .unwrap();
    assert_eq!(
        peptide.hits[0].target_decoy_type().unwrap(),
        TargetDecoyType::Decoy
    );
}

#[test]
fn the_mascot_base_name_workaround_rolls_a_wrong_run_back() {
    // A run summary without base_name defers the experiment check to
    // search_summary, which pops the run again when it does not match.
    let text = TEST.replace(
        "<msms_run_summary base_name=\"PepXMLFile_test\"",
        "<msms_run_summary",
    );
    let document = load(&text, &options("PepXMLFile_test"));
    assert_eq!(document.protein_identifications.len(), 2);
    assert_eq!(document.peptide_identifications.len(), 18);
    assert!(
        document
            .warnings
            .iter()
            .any(|warning| warning.contains("'base_name' attribute")),
        "{:?}",
        document.warnings
    );
    let error = pepxml::read_with_options(text.as_bytes(), &options("something_else")).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error}");
}

#[test]
fn an_unknown_sample_enzyme_does_not_abort_the_load() {
    let text = TEST.replace(
        "<sample_enzyme name=\"trypsin\">",
        "<sample_enzyme name=\"mystery\">",
    );
    let document = load(&text, &options("PepXMLFile_test"));
    let parameters = &document.protein_identifications[0].search_parameters;
    // The enzymatic_search_constraint of this fixture still names trypsin, and
    // the source lets that element overwrite the sample_enzyme.
    assert_eq!(parameters.digestion_enzyme, "Trypsin");
    // Without that element the enzyme stays unknown. The source instead calls
    // ProteaseDB::getEnzyme("mystery") unguarded from 'search_summary' and
    // throws Exception::ElementNotFound there, so it cannot read this file at
    // all; see OpenMS_CPP_ISSUES PEPXML-003.
    let text = text.replace(
        "<enzymatic_search_constraint enzyme=\"trypsin\" max_num_internal_cleavages=\"2\" min_number_termini=\"2\" />",
        "",
    );
    let document = load(&text, &options("PepXMLFile_test"));
    let parameters = &document.protein_identifications[0].search_parameters;
    assert_eq!(parameters.digestion_enzyme, "unknown_enzyme");
    // A 'specificity' element does build a user-defined enzyme, but the source's
    // 'search_summary' handler resets the whole SearchParameters block before
    // the run records it, so that value is unreachable in both ports; see
    // PEPXML-004. The writer's cleavage-site extraction still has to survive
    // such an expression, which arrives from idXML instead.
    let written = store_text(&document, &WriteOptions::default());
    assert!(
        written.contains("<specificity cut=\"\" sense=\"C\"/>"),
        "{written}"
    );
}

#[test]
fn a_custom_cleavage_expression_is_written_without_reading_past_its_end() {
    let mut document = load(STORE, &ReadOptions::default());
    // The shape a user-defined enzyme produces: cut residues and no ')' at all.
    document.protein_identifications[0]
        .search_parameters
        .digestion_enzyme = "user-defined,mystery,KR,P,C".into();
    document.protein_identifications[0]
        .search_parameters
        .digestion_regex = "KR".into();
    let written = store_text(&document, &WriteOptions::default());
    assert!(
        written.contains("<specificity cut=\"KR\" sense=\"C\"/>"),
        "{written}"
    );
    assert!(written.contains("<sample_enzyme name=\"user-defined,mystery,kr,p,c\">"));
}

#[test]
fn modification_declarations_report_their_diagnostics() {
    // mass == massdiff is wrong and the source recomputes the absolute mass.
    let text = TEST.replace(
        "<aminoacid_modification aminoacid=\"M\" massdiff=\"15.9949\" mass=\"147.0354\" variable=\"Y\" />",
        "<aminoacid_modification aminoacid=\"M\" massdiff=\"15.9949\" mass=\"15.9949\" variable=\"Y\" />",
    );
    let document = load(&text, &options("PepXMLFile_test"));
    assert!(
        document
            .warnings
            .iter()
            .any(|warning| warning.contains("mass == massdiff")),
        "{:?}",
        document.warnings
    );
    assert!(
        document.protein_identifications[0]
            .search_parameters
            .variable_modifications
            .contains(&"Oxidation (M)".to_owned())
    );
    // A declaration with neither an amino acid nor a terminus is fatal.
    let text = TEST.replace(
        "<aminoacid_modification aminoacid=\"M\" massdiff=\"15.9949\" mass=\"147.0354\" variable=\"Y\" />",
        "<aminoacid_modification aminoacid=\"\" massdiff=\"15.9949\" mass=\"147.0354\" variable=\"Y\" />",
    );
    let error =
        pepxml::read_with_options(text.as_bytes(), &options("PepXMLFile_test")).unwrap_err();
    assert!(matches!(error, Error::MissingInformation(_)), "{error}");
    // An X!Tandem artefact terminal modification is dropped silently.
    let text = TEST.replace(
        "<terminal_modification terminus=\"c\" protein_terminus=\"c\" massdiff=\"3.4\" mass=\"123.456\" variable=\"Y\" />",
        "<terminal_modification terminus=\"c\" protein_terminus=\"c\" massdiff=\"0.0002\" mass=\"123.456\" variable=\"Y\" />",
    );
    let document = load(&text, &options("PepXMLFile_test"));
    assert_eq!(
        document.protein_identifications[0]
            .search_parameters
            .variable_modifications
            .len(),
        11
    );
}

#[test]
fn preferred_modifications_win_over_a_registry_mass_search() {
    use openms::chemistry::{ModificationsDB, TermSpecificity};
    let registry = ModificationsDB::global();
    let preferred = registry
        .get_modification_handle("Carbamyl", Some('C'), Some(TermSpecificity::Anywhere))
        .unwrap();
    let read = ReadOptions {
        preferred_variable_modifications: vec![preferred],
        ..options("PepXMLFile_test")
    };
    // Carbamyl is +43.005814, far from every declared mass, so a preferred
    // modification list that cannot explain a mass changes nothing.
    let document = load(TEST, &read);
    assert_eq!(
        document.protein_identifications[0]
            .search_parameters
            .variable_modifications[2],
        "Oxidation (M)"
    );
    // A preferred entry named exactly as the declaration's description wins.
    let text = TEST.replace(
        "<aminoacid_modification aminoacid=\"M\" massdiff=\"15.9949\" mass=\"147.0354\" variable=\"Y\" />",
        "<aminoacid_modification aminoacid=\"C\" massdiff=\"43.005814\" mass=\"146.014963\" variable=\"Y\" description=\"Carbamyl (C)\" />",
    );
    let document = load(&text, &read);
    assert!(
        document.protein_identifications[0]
            .search_parameters
            .variable_modifications
            .contains(&"Carbamyl (C)".to_owned())
    );
}

#[test]
fn writing_rejects_records_that_cannot_be_represented() {
    let mut document = load(STORE, &ReadOptions::default());
    document.peptide_identifications[0].score_type = "Posterior Error Probability".into();
    // A posterior error probability adds a peptideprophet analysis result.
    let written = store_text(&document, &WriteOptions::default());
    assert!(written.contains("name=\"Posterior Error Probability\""));
    assert!(written.contains("<analysis_result analysis=\"peptideprophet\">"));
    // A Percolator run without a PEP score is missing information.
    let mut document = load(STORE, &ReadOptions::default());
    document.protein_identifications[0].search_engine = "Percolator".into();
    let error = pepxml::write_with_options(&mut Vec::new(), &document, &WriteOptions::default())
        .unwrap_err();
    assert!(matches!(error, Error::MissingInformation(_)), "{error}");
}

#[test]
fn store_publishes_atomically_and_round_trips_through_a_file() {
    let directory =
        std::env::temp_dir().join(format!("openms-pepxml-{}-{}", std::process::id(), line!()));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("round_trip.pepxml");
    let document = load(STORE, &ReadOptions::default());
    pepxml::store(&path, &document).unwrap();
    let reread = pepxml::load(&path).unwrap();
    assert_eq!(reread.peptide_identifications.len(), 9);
    // The output path's stem becomes the base name when none is supplied.
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.contains("base_name=\"round_trip\""));
    // A failing write leaves the existing file untouched.
    let error = pepxml::store_with_options(
        &path,
        &document,
        &WriteOptions {
            max_output_bytes: 16,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error}");
    assert_eq!(std::fs::read_to_string(&path).unwrap(), text);
    std::fs::remove_dir_all(&directory).ok();
}
