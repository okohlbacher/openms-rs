// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

//! `PeptideProteinResolution_test.cpp` at core `bc9cc12` (section `run`, two
//! input/expected-output idXML pairs compared with `TEST_FILE_SIMILAR`), and
//! invariants derived from the source for `resolve` and the graph search.

use openms::analysis::id::peptide_protein_resolution::{
    ConnectedComponent, PeptideProteinResolution,
};
use openms::chemistry::AASequence;
use openms::identification::{
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinGroup, ProteinHit,
    ProteinIdentification,
};

fn group(probability: f64, accessions: &[&str]) -> ProteinGroup {
    ProteinGroup {
        probability,
        accessions: accessions.iter().map(|a| (*a).to_owned()).collect(),
        ..ProteinGroup::default()
    }
}

fn pep(accessions: &[&str]) -> PeptideIdentification {
    let mut hit = PeptideHit::new(0.5, 1, 2, AASequence::parse("PEPTIDE").unwrap()).unwrap();
    for a in accessions {
        hit.evidences.push(PeptideEvidence::new(*a, 0..=6).unwrap());
    }
    PeptideIdentification {
        hits: vec![hit],
        ..PeptideIdentification::default()
    }
}

fn accessions(p: &PeptideIdentification) -> Vec<String> {
    p.hits[0]
        .evidences
        .iter()
        .map(|e| e.protein_accession.clone())
        .collect()
}

#[test]
fn connected_component_display_matches_the_source_operator() {
    let c = ConnectedComponent {
        prot_grp_indices: [0, 2].into(),
        pep_indices: [1].into(),
    };
    assert_eq!(c.to_string(), "Proteins: 0,2,\nPeptides: 1,");
}

#[test]
fn resolve_keeps_the_source_group_choice() {
    // Derived from the source resolve(): within targets the lower group
    // probability wins; target vs decoy prefers the higher probability.
    let mut protein = ProteinIdentification {
        hits: ["T1", "T2", "D1"]
            .iter()
            .map(|a| ProteinHit::new(0.0, 1, *a, "").unwrap())
            .collect(),
        indistinguishable_groups: vec![
            group(0.9, &["T1"]),
            group(0.4, &["T2"]),
            group(0.4, &["D1"]),
        ],
        ..ProteinIdentification::default()
    };
    protein.hits[2]
        .set_target_decoy_type(openms::identification::TargetDecoyType::Decoy)
        .unwrap();
    let mut peptides = vec![pep(&["T1", "T2"]), pep(&["T1", "D1"]), pep(&["T2", "D1"])];
    PeptideProteinResolution::resolve(&protein, &mut peptides, false, false).unwrap();
    assert_eq!(accessions(&peptides[0]), ["T2"]);
    assert_eq!(accessions(&peptides[1]), ["T1"]);
    // equal probabilities: both best groups are kept...
    assert_eq!(accessions(&peptides[2]), ["T2", "D1"]);
    // ...unless ties are resolved (lowest group index) or targets go first
    let mut tied = vec![pep(&["T2", "D1"])];
    PeptideProteinResolution::resolve(&protein, &mut tied, true, false).unwrap();
    assert_eq!(accessions(&tied[0]), ["T2"]);
    let mut first = vec![pep(&["T2", "D1"])];
    PeptideProteinResolution::resolve(&protein, &mut first, false, true).unwrap();
    assert_eq!(accessions(&first[0]), ["T2"]);

    let mut unknown = vec![pep(&["X"])];
    assert!(PeptideProteinResolution::resolve(&protein, &mut unknown, false, false).is_err());
    assert_eq!(accessions(&unknown[0]), ["X"]);
    protein.indistinguishable_groups.clear();
    assert!(PeptideProteinResolution::resolve(&protein, &mut peptides, false, false).is_err());
}

#[cfg(feature = "idxml")]
mod files {
    use super::*;
    use openms::format::idxml;
    use std::path::Path;

    fn load(name: &str) -> (Vec<ProteinIdentification>, Vec<PeptideIdentification>) {
        let d = idxml::load(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/data")
                .join(name),
        )
        .unwrap();
        (d.protein_identifications, d.peptide_identifications)
    }

    fn groups(groups: &[ProteinGroup]) -> Vec<(f64, Vec<String>)> {
        groups
            .iter()
            .map(|g| (g.probability, g.accessions.clone()))
            .collect()
    }

    fn check(input: &str, expected: &str) {
        let (mut prots, mut peps) = load(input);
        PeptideProteinResolution::run(&mut prots, &mut peps).unwrap();
        let (want_prots, want_peps) = load(expected);
        assert_eq!(prots.len(), want_prots.len());
        for (g, w) in prots.iter().zip(&want_prots) {
            let hits = |p: &ProteinIdentification| {
                p.hits
                    .iter()
                    .map(|h| (h.accession.clone(), h.score))
                    .collect::<Vec<_>>()
            };
            assert_eq!(hits(g), hits(w));
            assert_eq!(
                groups(&g.indistinguishable_groups),
                groups(&w.indistinguishable_groups)
            );
            assert_eq!(groups(&g.protein_groups), groups(&w.protein_groups));
        }
        assert_eq!(peps.len(), want_peps.len());
        for (g, w) in peps.iter().zip(&want_peps) {
            assert_eq!(g.hits.len(), w.hits.len());
            for (gh, wh) in g.hits.iter().zip(&w.hits) {
                assert_eq!(gh.sequence, wh.sequence);
                assert_eq!(gh.evidences, wh.evidences, "{}", gh.sequence);
            }
        }
    }

    #[test]
    fn run_class_test_first_pair() {
        check(
            "PeptideProteinResolution_in.idXML",
            "PeptideProteinResolution_out.idXML",
        );
    }

    #[test]
    fn run_class_test_second_pair() {
        check(
            "PeptideProteinResolution_in2.idXML",
            "PeptideProteinResolution_out2.idXML",
        );
    }
}
