// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::analysis::false_discovery_rate::FalseDiscoveryRate;
use openms::analysis::peptide_indexing::{DecoyRule, PeptideIndexing};
use openms::chemistry::{AASequence, DecoyGenerator, Protease, ProteaseDigestion};
use openms::format::fasta::{self, FASTAEntry};
use openms::identification::{
    PeptideHit, PeptideIdentification, ProteinIdentification, TargetDecoyType,
};

fn indexed_decoys() -> (Vec<ProteinIdentification>, Vec<PeptideIdentification>) {
    let mut database = fasta::read(&b">P1 source target\nTESTPEPTIDE\n"[..]).unwrap();
    let target = AASequence::parse(&database[0].sequence).unwrap();
    let decoy = DecoyGenerator::with_seed(17)
        .shuffle(&target, Protease::Trypsin, 1)
        .unwrap()
        .remove(0);
    // Literal upstream multi-variant result; this operation ignores the receiver seed.
    assert_eq!(decoy.as_str().to_owned(), "DIESETEPTTP");
    database.push(FASTAEntry {
        identifier: "DECOY_P1".into(),
        description: "variant 0".into(),
        sequence: decoy.as_str().to_owned(),
    });
    let mut bytes = Vec::new();
    fasta::write(&mut bytes, &database).unwrap();
    let database = fasta::read(bytes.as_slice()).unwrap();
    let mut proteins = vec![ProteinIdentification {
        identifier: "decoy-workflow".into(),
        date_time: Some("2026-09-10T12:00:00".into()),
        search_engine: "synthetic score example".into(),
        ..Default::default()
    }];
    let mut peptides: Vec<_> = [target, decoy]
        .into_iter()
        .enumerate()
        .map(|(i, sequence)| PeptideIdentification {
            identifier: "decoy-workflow".into(),
            score_type: "synthetic score".into(),
            higher_score_better: true,
            hits: vec![PeptideHit::new(20.0 - i as f64 * 10.0, 1, 2, sequence).unwrap()],
            ..Default::default()
        })
        .collect();
    let report = PeptideIndexing {
        decoy_rule: DecoyRule::Prefix("DECOY_".into()),
        enzyme: Some(Protease::Trypsin),
        write_protein_sequence: true,
        ..Default::default()
    }
    .run(&database, &mut proteins, &mut peptides)
    .unwrap();
    assert_eq!((report.target_hits, report.decoy_hits), (1, 1));
    assert_eq!(report.target_and_decoy_hits, 0);
    for (id, expected, accession) in [
        (&peptides[0], TargetDecoyType::Target, "P1"),
        (&peptides[1], TargetDecoyType::Decoy, "DECOY_P1"),
    ] {
        assert_eq!(id.hits[0].target_decoy_type().unwrap(), expected);
        assert_eq!(id.hits[0].evidences.len(), 1);
        assert_eq!(id.hits[0].evidences[0].protein_accession, accession);
    }
    (proteins, peptides)
}

#[test]
fn generated_fasta_decoys_index_and_receive_independent_basic_q_values() {
    let (_, mut peptides) = indexed_decoys();
    FalseDiscoveryRate {
        add_decoy_peptides: true,
        ..Default::default()
    }
    .apply_basic_peptides(&mut peptides)
    .unwrap();
    // Basic conservative formula (D+1)/(T+1): 1/2 then 2/2 at the two thresholds.
    assert_eq!(peptides[0].hits[0].score, 0.5);
    assert_eq!(peptides[1].hits[0].score, 1.0);
    assert_eq!(peptides[0].score_type, "q-value");
    assert!(!peptides[0].higher_score_better);
}

#[test]
fn nonoverlapping_decoys_retain_element_counts_through_native_digestion() {
    let target = AASequence::parse("ACDMKPEPTIDER").unwrap();
    let generator = DecoyGenerator::with_seed(4711);
    let enzyme = Protease::from_name("Trypsin/P").unwrap();
    let mut decoys = generator.shuffle(&target, enzyme, 2).unwrap();
    decoys.push(generator.reverse_protein(&target).unwrap());
    decoys.push(generator.reverse_peptides(&target, enzyme).unwrap());
    let digestion = ProteaseDigestion {
        enzyme,
        missed_cleavages: 0,
        ..Default::default()
    };
    for decoy in decoys {
        assert_eq!(decoy.formula().unwrap(), target.formula().unwrap());
        let products = digestion.digest(&decoy).unwrap();
        assert_eq!(
            products.iter().map(|p| p.sequence.len()).sum::<usize>(),
            target.len()
        );
        assert_eq!(products[0].start, 0);
        assert_eq!(products.last().unwrap().end, target.len());
    }
}

#[cfg(feature = "idxml")]
#[test]
fn generated_database_evidence_and_classification_roundtrip_through_idxml() {
    use openms::format::idxml::{self, IdXmlDocument};
    let (proteins, peptides) = indexed_decoys();
    let document = IdXmlDocument {
        protein_identifications: proteins,
        peptide_identifications: peptides,
        ..Default::default()
    };
    let mut bytes = Vec::new();
    idxml::write(&mut bytes, &document).unwrap();
    assert_eq!(idxml::read(bytes.as_slice()).unwrap(), document);
}
