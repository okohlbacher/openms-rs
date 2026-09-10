// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "idxml")]
use openms::analysis::false_discovery_rate::FalseDiscoveryRate;
use openms::analysis::id_filter::{remove_empty_peptide_identifications, update_protein_groups};
use openms::analysis::id_ripper::IDRipper;
use openms::analysis::peptide_indexing::{DecoyRule, PeptideIndexing};
use openms::analysis::protein_inference::BasicProteinInference;
use openms::chemistry::{AASequence, DigestionSpecificity, Protease};
use openms::format::{fasta, idxml};
use openms::identification::{
    PeptideHit, PeptideIdentification, ProteinIdentification, TargetDecoyType,
};
use std::io::Cursor;

#[test]
fn fasta_indexing_inference_fdr_and_origin_partition_preserve_scientific_records() {
    // Independent tiny protein/PSM population. These hand-calculated FDRs test
    // counting conventions; they are not an empirical confidence calibration.
    let database = fasta::read(Cursor::new(b">P1 first target\nAAAKCCCK\n>P2 second target\nAAAKDDDK\n>DECOY_P1 synthetic decoy\nEEEKK\n>UNREFERENCED no evidence\nFFFK\n")).unwrap();
    let mut runs = vec![ProteinIdentification {
        identifier: "run".into(),
        search_engine: "synthetic".into(),
        search_engine_version: "1".into(),
        date_time: Some("2026-09-10T12:00:00".into()),
        primary_ms_run_paths: vec!["first.mzML".into(), "second.mzML".into()],
        ..Default::default()
    }];
    let mut peptides: Vec<_> = [
        ("AAAK", 0.05, 0_i64),
        (
            "C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)K",
            0.01,
            0,
        ),
        ("DDDK", 0.1, 1),
        ("EEEK", 0.5, 1),
    ]
    .into_iter()
    .enumerate()
    .map(|(scan, (sequence, score, origin))| {
        let mut id = PeptideIdentification {
            identifier: "run".into(),
            score_type: "Posterior Error Probability".into(),
            higher_score_better: false,
            hits: vec![PeptideHit::new(score, 7, 2, AASequence::parse(sequence).unwrap()).unwrap()],
            rt: Some(10. * (scan + 1) as f64),
            ..Default::default()
        };
        id.set_spectrum_reference(format!("scan={}", scan + 1));
        id.metadata.insert("id_merge_index".into(), origin.into());
        id
    })
    .collect();
    let indexing = PeptideIndexing {
        decoy_rule: DecoyRule::Prefix("DECOY_".into()),
        enzyme: Some(Protease::Trypsin),
        specificity: Some(DigestionSpecificity::Full),
        write_protein_sequence: true,
        write_protein_description: true,
        ..Default::default()
    }
    .run(&database, &mut runs, &mut peptides)
    .unwrap();
    assert_eq!(
        (
            indexing.target_hits,
            indexing.decoy_hits,
            indexing.non_unique_hits,
            indexing.evidence_count
        ),
        (3, 1, 1, 5)
    );
    assert_eq!(
        peptides[0].hits[0]
            .protein_accessions()
            .into_iter()
            .collect::<Vec<_>>(),
        ["P1", "P2"]
    );
    assert_eq!(peptides[1].hits[0].evidences[0].start, Some(4));
    assert_eq!(peptides[1].hits[0].evidences[0].end, Some(7));
    assert_eq!(
        peptides[3].hits[0].target_decoy_type().unwrap(),
        TargetDecoyType::Decoy
    );
    assert_eq!(runs[0].hits.len(), 3);
    let first_protein = runs[0]
        .hits
        .iter()
        .find(|hit| hit.accession == "P1")
        .unwrap();
    assert_eq!(first_protein.sequence, "AAAKCCCK");
    assert_eq!(first_protein.description(), "first target");
    assert!(peptides[1].hits[0].sequence.is_modified());
    BasicProteinInference::default()
        .run(&mut peptides, &mut runs)
        .unwrap();
    let score = |accession: &str| {
        runs[0]
            .hits
            .iter()
            .find(|h| h.accession == accession)
            .unwrap()
            .score
    };
    assert!((score("P1") - 0.99).abs() < 1e-14);
    assert!((score("P2") - 0.95).abs() < 1e-14);
    assert_eq!(score("DECOY_P1"), 0.5);
    assert_eq!(runs[0].score_type, "Posterior Probability");
    assert_eq!(runs[0].indistinguishable_groups.len(), 3);
    assert_eq!(peptides[1].hits[0].score, 0.01);
    let fdr = FalseDiscoveryRate::default();
    fdr.apply_basic_protein(&mut runs[0], true).unwrap();
    assert_eq!(runs[0].hits.len(), 2);
    // With two higher-scoring targets and one lower decoy, Basic's source
    // pseudocount formula gives target q = (0+1)/(2+1) = 1/3.
    for hit in &runs[0].hits {
        assert!((hit.score - 1. / 3.).abs() < 1e-14);
    }
    let protein_hits = runs[0].hits.clone();
    update_protein_groups(&mut runs[0].indistinguishable_groups, &protein_hits).unwrap();
    assert_eq!(runs[0].indistinguishable_groups.len(), 2);
    fdr.apply_peptides(&mut peptides, false).unwrap();
    remove_empty_peptide_identifications(&mut peptides).unwrap();
    assert_eq!(peptides.len(), 3);
    assert!(peptides.iter().all(|id| id.hits[0].score == 0.));
    runs[0].compute_coverage(&peptides).unwrap();
    assert!(runs[0].hits.iter().all(|hit| hit.coverage == Some(100.)));
    let document = idxml::IdXmlDocument {
        document_id: "protein-workflow".into(),
        protein_identifications: runs,
        peptide_identifications: peptides,
        ..Default::default()
    };
    let mut xml = Vec::new();
    idxml::write(&mut xml, &document).unwrap();
    let restored = idxml::read(Cursor::new(xml)).unwrap();
    assert_eq!(restored, document);
    let split = IDRipper::default()
        .rip(
            &restored.protein_identifications,
            &restored.peptide_identifications,
        )
        .unwrap();
    assert_eq!(split.files.len(), 2);
    assert_eq!(split.files[0].identifier.origin_fullname, "first.mzML");
    assert_eq!(split.files[1].identifier.origin_fullname, "second.mzML");
    assert_eq!(split.files[0].peptide_identifications.len(), 2);
    assert_eq!(split.files[1].peptide_identifications.len(), 1);
    assert_eq!(split.files[0].protein_identifications[0].hits.len(), 2);
    assert_eq!(
        split.files[1].protein_identifications[0].hits[0].accession,
        "P2"
    );
    assert_eq!(split.groups_not_copied, 4);
    assert!(split.skipped_peptide_identifications.is_empty());
}
