// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `SimpleSearchEngineAlgorithm` against the retained outputs of the upstream
//! TOPP tests `TOPP_SimpleSearchEngine_1` and `_2` (test-data `0cb15f2`,
//! `topp/CMakeLists.txt:3773-3789`), evidence tier 1, compared as decoded
//! idXML under the upstream FuzzyDiff contract (`topp/FuzzyDiff.ini`: a number
//! passes within ratio 1.01 or absolute difference 0.01; text exactly), with
//! the upstream whitelist (`IdentificationRun date`, `SearchParameters db=`).
//! Native checks (tier 4) follow.
#![cfg(all(feature = "mzml", feature = "idxml", feature = "paramxml"))]

use openms::Result;
use openms::analysis::id::simple_search_engine::{
    AnnotatedHit, ExitCodes, MassToleranceUnit, SimpleSearchEngineAlgorithm,
};
use openms::concept::parallel::Threads;
use openms::format::{idxml, paramxml};
use openms::identification::{PeptideIdentification, ProteinIdentification};
use openms::metadata::{MetaInfo, MetaValueData};
use std::path::PathBuf;
use std::sync::Arc;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/simple_search_engine")
        .join(name)
}

fn engine_from_ini(ini: &str) -> Result<SimpleSearchEngineAlgorithm> {
    let tool = paramxml::load(data(ini))?;
    let search = tool.copy("SimpleSearchEngine:1:Search:", true)?;
    SimpleSearchEngineAlgorithm::from_param(&search)
}

/// FuzzyDiff with `ratio 1.01` and `absdiff 0.01`: equal, or within either.
fn fuzzy_eq(a: f64, b: f64) -> bool {
    if a == b || (a - b).abs() <= 0.01 {
        return true;
    }
    let (lo, hi) = if a.abs() < b.abs() { (a, b) } else { (b, a) };
    lo != 0.0 && lo.signum() == hi.signum() && hi / lo <= 1.01
}

fn assert_fuzzy(what: &str, actual: f64, expected: f64) {
    assert!(fuzzy_eq(actual, expected), "{what}: {actual} vs {expected}");
}

fn assert_meta(what: &str, actual: &MetaInfo, expected: &MetaInfo) {
    let keys = |m: &MetaInfo| m.keys().cloned().collect::<Vec<_>>();
    assert_eq!(keys(actual), keys(expected), "{what}: meta keys");
    for (key, value) in expected {
        let got = &actual[key];
        match (value.data(), got.data()) {
            (MetaValueData::Float(e), MetaValueData::Float(a)) => {
                assert_fuzzy(&format!("{what} {key}"), *a, *e)
            }
            (MetaValueData::Integer(e), MetaValueData::Integer(a)) => {
                assert_fuzzy(&format!("{what} {key}"), *a as f64, *e as f64)
            }
            _ => assert_eq!(got.to_string(), value.to_string(), "{what} {key}"),
        }
    }
}

fn run(ini: &str) -> Result<(Vec<ProteinIdentification>, Vec<PeptideIdentification>)> {
    let engine = engine_from_ini(ini)?;
    let (mut proteins, mut peptides) = (Vec::new(), Vec::new());
    let code = engine.search(
        data("SimpleSearchEngine_1.mzML")
            .to_str()
            .expect("UTF-8 path"),
        data("SimpleSearchEngine_1.fasta")
            .to_str()
            .expect("UTF-8 path"),
        &mut proteins,
        &mut peptides,
    )?;
    assert_eq!(code, ExitCodes::ExecutionOk);
    Ok((proteins, peptides))
}

fn compare(ini: &str, expected_file: &str) -> Result<()> {
    let (proteins, peptides) = run(ini)?;
    let expected = idxml::load(data(expected_file))?;

    assert_eq!(proteins.len(), expected.protein_identifications.len());
    let (run, want) = (&proteins[0], &expected.protein_identifications[0]);
    assert_eq!(run.search_engine, want.search_engine);
    assert_eq!(run.score_type, want.score_type);
    assert_eq!(run.higher_score_better, want.higher_score_better);
    let accessions = |p: &ProteinIdentification| {
        p.hits
            .iter()
            .map(|h| h.accession.clone())
            .collect::<Vec<_>>()
    };
    assert_eq!(accessions(run), accessions(want));
    for (hit, want_hit) in run.hits.iter().zip(&want.hits) {
        assert_meta(&hit.accession, &hit.metadata, &want_hit.metadata);
    }
    assert_eq!(run.raw_ms_run_paths, want.raw_ms_run_paths);
    // The tool passes the input as given; the retained run recorded it as a
    // file URI of the bare name.
    let names = |paths: &[String]| {
        paths
            .iter()
            .map(|p| p.rsplit('/').next().unwrap_or(p).to_owned())
            .collect::<Vec<_>>()
    };
    assert_eq!(
        names(&run.primary_ms_run_paths),
        names(&want.primary_ms_run_paths)
    );

    let (sp, want_sp) = (&run.search_parameters, &want.search_parameters);
    assert_eq!(sp.charges, want_sp.charges);
    assert_eq!(sp.mass_type, want_sp.mass_type);
    assert_eq!(sp.fixed_modifications, want_sp.fixed_modifications);
    assert_eq!(sp.variable_modifications, want_sp.variable_modifications);
    assert_eq!(sp.missed_cleavages, want_sp.missed_cleavages);
    assert_eq!(sp.fragment_tolerance, want_sp.fragment_tolerance);
    assert_eq!(sp.precursor_tolerance, want_sp.precursor_tolerance);
    assert_eq!(
        sp.digestion_enzyme.to_lowercase(),
        want_sp.digestion_enzyme.to_lowercase()
    );
    assert_eq!(sp.enzyme_specificity, want_sp.enzyme_specificity);
    assert_meta("search parameters", &sp.metadata, &want_sp.metadata);

    assert_eq!(peptides.len(), expected.peptide_identifications.len());
    for (id, want_id) in peptides.iter().zip(&expected.peptide_identifications) {
        let what = want_id.spectrum_reference();
        assert_eq!(id.spectrum_reference(), what);
        assert_eq!(id.score_type, want_id.score_type);
        assert_eq!(id.higher_score_better, want_id.higher_score_better);
        assert_fuzzy(&what, id.mz.unwrap(), want_id.mz.unwrap());
        assert_fuzzy(&what, id.rt.unwrap(), want_id.rt.unwrap());
        assert_meta(&what, &id.metadata, &want_id.metadata);
        assert_eq!(id.hits.len(), want_id.hits.len(), "{what}");
        for (hit, want_hit) in id.hits.iter().zip(&want_id.hits) {
            assert_eq!(
                hit.sequence.to_string(),
                want_hit.sequence.to_string(),
                "{what}"
            );
            assert_eq!(hit.charge, want_hit.charge, "{what}");
            assert_fuzzy(&format!("{what} score"), hit.score, want_hit.score);
            assert_eq!(hit.evidences, want_hit.evidences, "{what}");
            assert_meta(&what, &hit.metadata, &want_hit.metadata);
            assert_eq!(
                hit.peak_annotations.len(),
                want_hit.peak_annotations.len(),
                "{what} annotations"
            );
            for (a, e) in hit.peak_annotations.iter().zip(&want_hit.peak_annotations) {
                assert_eq!(a.annotation, e.annotation, "{what}");
                assert_eq!(a.charge, e.charge, "{what}");
                assert_fuzzy(&format!("{what} annotation m/z"), a.mz, e.mz);
                assert_fuzzy(
                    &format!("{what} annotation intensity"),
                    a.intensity,
                    e.intensity,
                );
            }
        }
    }
    Ok(())
}

// TOPP_SimpleSearchEngine_1: 0.1 Da fragments, no fixed modification.
#[test]
fn topp_simple_search_engine_1_matches_retained_output() -> Result<()> {
    compare("SimpleSearchEngine_1.ini", "SimpleSearchEngine_1_out.idXML")
}

// TOPP_SimpleSearchEngine_2: 10 ppm fragments, Carbamidomethyl (C) fixed.
#[test]
fn topp_simple_search_engine_2_matches_retained_output() -> Result<()> {
    compare("SimpleSearchEngine_2.ini", "SimpleSearchEngine_2_out.idXML")
}

// Tier 4: the result does not depend on the thread count.
#[test]
fn results_are_identical_across_thread_counts() -> Result<()> {
    let search = |threads: Threads| -> Result<Vec<PeptideIdentification>> {
        let mut engine = engine_from_ini("SimpleSearchEngine_2.ini")?;
        engine.threads = threads;
        let (mut proteins, mut peptides) = (Vec::new(), Vec::new());
        engine.search(
            data("SimpleSearchEngine_1.mzML").to_str().unwrap(),
            data("SimpleSearchEngine_1.fasta").to_str().unwrap(),
            &mut proteins,
            &mut peptides,
        )?;
        for id in &mut peptides {
            id.identifier.clear();
        }
        Ok(peptides)
    };
    let serial = search(Threads::serial())?;
    assert!(!serial.is_empty());
    assert_eq!(search(Threads::from_cli(4))?, serial);
    Ok(())
}

// Tier 4: parameter round trip and source defaults.
#[test]
fn parameters_round_trip_with_source_defaults() -> Result<()> {
    let defaults = SimpleSearchEngineAlgorithm::defaults()?;
    let engine = SimpleSearchEngineAlgorithm::from_param(&defaults)?;
    assert_eq!(engine, SimpleSearchEngineAlgorithm::default());
    assert_eq!(engine.precursor_mass_tolerance_unit, MassToleranceUnit::Ppm);
    assert_eq!(engine.precursor_isotopes, vec![0, 1]);
    assert_eq!(engine.report_top_hits, 1);
    assert_eq!(
        defaults.description("fragment:mass_tolerance_unit")?,
        "Unit of fragment m"
    );
    let mut changed = engine.clone();
    changed.decoys = true;
    changed.fdr_psm = 0.05;
    changed.peptide_motif = "K$".into();
    changed.fragment_mass_tolerance_unit = MassToleranceUnit::Da;
    assert_eq!(
        SimpleSearchEngineAlgorithm::from_param(&changed.to_param()?)?,
        changed
    );
    Ok(())
}

// Tier 4: duplicates are made unique (sorted, as the source's std::set) with a
// warning; restrictions and negative counts are refused.
#[test]
fn parameter_updates_follow_update_members() -> Result<()> {
    let mut param = SimpleSearchEngineAlgorithm::defaults()?;
    param.set_value(
        "modifications:variable",
        openms::param::ParamValue::StringList(vec![
            "Oxidation (M)".into(),
            "Acetyl (N-term)".into(),
            "Oxidation (M)".into(),
        ]),
        "",
        &[],
    )?;
    let (engine, warnings) = SimpleSearchEngineAlgorithm::from_param_with_warnings(&param)?;
    assert_eq!(
        engine.modifications_variable,
        vec!["Acetyl (N-term)".to_owned(), "Oxidation (M)".to_owned()]
    );
    assert!(
        warnings
            .iter()
            .any(|w| w.contains("Duplicate variable modification"))
    );

    let mut bad = SimpleSearchEngineAlgorithm::defaults()?;
    bad.set_value(
        "precursor:mass_tolerance_unit",
        openms::param::ParamValue::String("mDa".into()),
        "",
        &[],
    )?;
    assert!(SimpleSearchEngineAlgorithm::from_param(&bad).is_err());
    let mut negative = SimpleSearchEngineAlgorithm::defaults()?;
    negative.set_value(
        "report:top_hits",
        openms::param::ParamValue::Integer(-1),
        "",
        &[],
    )?;
    assert!(SimpleSearchEngineAlgorithm::from_param(&negative).is_err());
    Ok(())
}

// Tier 4: the source's shortlex tie-break (SimpleSearchEngineAlgorithm.h
// hasBetterScore): score, then modification index, then length, then text.
#[test]
fn has_better_score_is_source_shortlex() {
    let hit = |score: f64, index: usize, sequence: &str| AnnotatedHit {
        sequence: Arc::from(sequence),
        peptide_mod_index: index,
        score,
        prefix_fraction: 0.0,
        suffix_fraction: 0.0,
        mean_error: 0.0,
        isotope_error: 0,
        matched_prefix_ions: 0,
        matched_suffix_ions: 0,
    };
    let better = AnnotatedHit::has_better_score;
    assert!(better(&hit(2.0, 5, "ZZZZ"), &hit(1.0, 0, "A")));
    assert!(better(&hit(1.0, 0, "ZZZZ"), &hit(1.0, 1, "A")));
    assert!(better(&hit(1.0, 0, "ZZ"), &hit(1.0, 0, "AAA")));
    assert!(better(&hit(1.0, 0, "AB"), &hit(1.0, 0, "AC")));
    assert!(!better(&hit(1.0, 0, "AB"), &hit(1.0, 0, "AB")));
}

// Tier 4: decoys double the database; the source's FDR filter then keeps only
// targets at the threshold, and a missing database is reported, not a panic.
#[test]
fn decoy_search_with_fdr_keeps_targets_only() -> Result<()> {
    let mut engine = engine_from_ini("SimpleSearchEngine_2.ini")?;
    engine.decoys = true;
    engine.fdr_psm = 0.5;
    let (mut proteins, mut peptides) = (Vec::new(), Vec::new());
    engine.search(
        data("SimpleSearchEngine_1.mzML").to_str().unwrap(),
        data("SimpleSearchEngine_1.fasta").to_str().unwrap(),
        &mut proteins,
        &mut peptides,
    )?;
    for id in &peptides {
        for hit in &id.hits {
            assert!(
                hit.evidences
                    .iter()
                    .all(|e| !e.protein_accession.starts_with("DECOY_"))
            );
        }
    }
    let missing = engine.search(
        data("SimpleSearchEngine_1.mzML").to_str().unwrap(),
        data("no_such.fasta").to_str().unwrap(),
        &mut proteins,
        &mut peptides,
    );
    assert!(missing.is_err());
    Ok(())
}
