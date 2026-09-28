// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

//! Every section of `IDMergerAlgorithm_test.cpp` at core `bc9cc12`.

use openms::analysis::id::id_merger::{IDMergerAlgorithm, IDMergerParameters};
use openms::chemistry::AASequence;
use openms::identification::{
    FlankingResidue, PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit,
    ProteinIdentification,
};
use openms::metadata::MetaValue;
use openms::param::ParamValue;

fn run(identifier: &str, paths: &[&str], accessions: &[&str]) -> ProteinIdentification {
    ProteinIdentification {
        identifier: identifier.into(),
        primary_ms_run_paths: paths.iter().map(|p| (*p).to_owned()).collect(),
        hits: accessions
            .iter()
            .map(|a| ProteinHit::new(0.0, 1, *a, "").unwrap())
            .collect(),
        ..ProteinIdentification::default()
    }
}

fn peptide(score: f64, sequence: &str, evidence: &[(&str, usize, usize)]) -> PeptideHit {
    let mut hit = PeptideHit::new(score, 1, 1, AASequence::parse(sequence).unwrap()).unwrap();
    for (accession, start, end) in evidence {
        let mut e = PeptideEvidence::new(*accession, *start..=*end).unwrap();
        e.aa_before = FlankingResidue::Residue('A');
        e.aa_after = FlankingResidue::Residue('A');
        hit.evidences.push(e);
    }
    hit
}

fn pep_id(run: &str, hits: &[&PeptideHit], index: i64) -> PeptideIdentification {
    let mut id = PeptideIdentification {
        identifier: run.into(),
        hits: hits.iter().map(|h| (*h).clone()).collect(),
        ..PeptideIdentification::default()
    };
    id.metadata
        .insert("id_merge_index".into(), MetaValue::from(index));
    id
}

#[test]
fn insert_runs_class_test_section_1() {
    let pr1 = run("PR1", &["f1r1.mzML", "f2r1.mzML"], &["A", "B", "C"]);
    let pr2 = run("PR2", &["f1r2.mzML", "f2r2.mzML"], &["A", "D", "E"]);
    let pr3 = run("PR3", &["f1r3.mzML", "f2r3.mzML"], &["D", "E", "F", "G"]);
    let pr4 = run("PR4", &["control.mzML"], &[]); // empty

    let ph0 = peptide(0.0, "AA", &[("A", 1, 3), ("B", 1, 3)]);
    let ph1 = peptide(0.0, "AAA", &[("A", 1, 4)]);
    let ph11 = peptide(0.0, "AAC", &[("C", 1, 4)]);
    let ph2 = peptide(0.0, "AAAA", &[("D", 1, 5), ("E", 1, 5)]);
    let ph3 = peptide(0.0, "AAAAA", &[("D", 1, 6), ("E", 1, 6)]);
    let ph4 = peptide(0.0, "AAAAAA", &[("F", 1, 7)]);
    // Same sequence, different proteins: cannot be recovered while merging.
    let ph5 = peptide(0.1, "AAA", &[("F", 1, 4)]);

    let pes = vec![
        pep_id("PR1", &[&ph0], 0),
        pep_id("PR1", &[&ph1, &ph11], 0),
        pep_id("PR2", &[&ph0], 0), // accessions not in its run: added nonetheless
        pep_id("PR2", &[&ph2, &ph3], 1),
        pep_id("PR3", &[&ph2], 0),
        pep_id("PR3", &[&ph3], 1),
        pep_id("PR3", &[&ph4], 0),
        pep_id("PR3", &[&ph5], 0), // can happen if the second file had no IDs
        pep_id("PR5", &[&ph5], 564), // non-existent run: ignored
    ];
    let prs = vec![pr1, pr2, pr3, pr4];
    let mut ima = IDMergerAlgorithm::new("mymerge", true);
    ima.insert_runs(prs, pes.clone()).unwrap();
    let (prres, peres) = ima.return_results_and_clear();

    assert_eq!(pes.len(), 9);
    assert_eq!(peres.len(), 8);
    assert_eq!(prres.hits.len(), 7);
    assert_eq!(prres.primary_ms_run_paths.len(), 7);
    assert_eq!(peres[2].metadata["id_merge_index"].as_i64().unwrap(), 2);
}

#[cfg(feature = "idxml")]
mod files {
    use super::*;
    use openms::format::idxml;
    use std::path::Path;

    /// Load a class-test file. Both files carry an `id="PI_0"` attribute on
    /// `IdentificationRun` that the IdXML 1.5 schema does not declare; the C++
    /// reader ignores it, the native reader refuses unknown attributes. It is
    /// removed in memory, which is what the source reader does with it.
    fn load(name: &str) -> (Vec<ProteinIdentification>, Vec<PeptideIdentification>) {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name);
        let text = std::fs::read_to_string(path).unwrap().replacen(
            "<IdentificationRun id=\"PI_0\" ",
            "<IdentificationRun ",
            1,
        );
        let d = idxml::read(text.as_bytes()).unwrap();
        (d.protein_identifications, d.peptide_identifications)
    }

    #[test]
    fn insert_runs_class_test_section_2() {
        let (pr1, pe1) = load("newIDMergerTest1.idXML");
        let pe1size = pe1.len();
        let (pr2, pe2) = load("newIDMergerTest2.idXML");
        let pe2size = pe2.len();
        let mut ima = IDMergerAlgorithm::new("mymerge", true);
        ima.insert_runs(pr1, pe1).unwrap();
        ima.insert_runs(pr2, pe2).unwrap();
        assert!(
            ima.insert_runs(
                vec![ProteinIdentification::default()],
                vec![PeptideIdentification::default()]
            )
            .is_err()
        );
        let (prres, peres) = ima.return_results_and_clear();
        assert_eq!(prres.hits.len(), 6);
        assert_eq!(peres.len(), pe1size + pe2size);
        assert!(prres.identifier.starts_with("mymerge"));
        assert_eq!(prres.primary_ms_run_paths.len(), 2);
    }

    #[test]
    fn check_search_setting_consistency_class_test_section_3() {
        let (pr1, pe1) = load("newIDMergerTest1.idXML");
        let (mut pr2, pe2) = load("newIDMergerTest2.idXML");
        // fail with a different db filename
        pr2[0].search_parameters.database = "baz".into();
        let mut ima = IDMergerAlgorithm::new("mymerge", true);
        ima.insert_runs(pr1, pe1).unwrap();
        assert!(ima.insert_runs(pr2.clone(), pe2.clone()).is_err());
        // a Windows path with the right file name is accepted
        pr2[0].search_parameters.database =
            "C:\\foo\\s_pyo_sf370_potato_human_target_decoy_with_contaminants.fasta".into();
        ima.insert_runs(pr2, pe2).unwrap();
        let (prres, _) = ima.return_results_and_clear();
        assert_eq!(prres.hits.len(), 6);
    }
}

#[test]
fn parameters_match_the_source_declarations() {
    let d = IDMergerParameters::defaults().unwrap();
    assert_eq!(
        d.value("annotate_origin").unwrap(),
        &ParamValue::String("true".into())
    );
    assert_eq!(
        d.value("allow_disagreeing_settings").unwrap(),
        &ParamValue::String("false".into())
    );
    let p = IDMergerParameters {
        annotate_origin: false,
        allow_disagreeing_settings: true,
    };
    assert_eq!(
        IDMergerParameters::from_param(&p.to_param().unwrap()).unwrap(),
        p
    );
}

#[test]
fn fixed_identifier_and_reset() {
    // Derived from the source: without a time stamp the identifier is the
    // given one, and returnResultsAndClear leaves an empty, reusable merger.
    let mut ima = IDMergerAlgorithm::new("fixed", false);
    assert_eq!(ima.identifier(), "fixed");
    ima.parameters.annotate_origin = false;
    let mut pep = PeptideIdentification {
        identifier: "R".into(),
        ..PeptideIdentification::default()
    };
    pep.hits.push(peptide(0.0, "PEPTIDE", &[("P1", 0, 6)]));
    let warnings = ima
        .insert_runs(vec![run("R", &[], &["P1", "P1"])], vec![pep])
        .unwrap();
    assert!(warnings.is_empty());
    let (prot, peps) = ima.return_results_and_clear();
    assert_eq!(prot.identifier, "fixed");
    assert_eq!(prot.hits.len(), 1);
    assert_eq!(peps[0].identifier, "fixed");
    assert!(!peps[0].metadata.contains_key("id_merge_index"));
    let (prot, peps) = ima.return_results_and_clear();
    assert!(prot.hits.is_empty() && peps.is_empty());
    assert_eq!(
        ima.insert_runs(Vec::new(), Vec::new()).unwrap(),
        vec!["No ProteinIdentification(Runs) given. Skipping.".to_owned()]
    );
}
