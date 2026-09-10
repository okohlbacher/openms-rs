// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::analysis::peptide_indexing::{
    DecoyRule, MissingDecoyAction, PeptideIndexing, UnmatchedAction,
};
use openms::chemistry::{AASequence, DigestionSpecificity, Protease};
use openms::format::fasta::FASTAEntry;
use openms::identification::{
    EnzymeTermSpecificity, FlankingResidue, PeptideHit, PeptideIdentification, ProteinHit,
    ProteinIdentification, TargetDecoyType,
};

fn entry(id: &str, sequence: &str) -> FASTAEntry {
    FASTAEntry {
        identifier: id.into(),
        description: format!("description {id}"),
        sequence: sequence.into(),
    }
}
fn run(id: &str) -> ProteinIdentification {
    let mut run = ProteinIdentification {
        identifier: id.into(),
        ..Default::default()
    };
    run.search_parameters.digestion_enzyme = "Trypsin".into();
    run.search_parameters.enzyme_specificity = EnzymeTermSpecificity::Full;
    run
}
fn identification(id: &str, sequence: &str) -> PeptideIdentification {
    PeptideIdentification {
        identifier: id.into(),
        hits: vec![PeptideHit {
            sequence: AASequence::parse(sequence).unwrap(),
            score: 12.5,
            rank: 7,
            ..Default::default()
        }],
        ..Default::default()
    }
}
fn permissive() -> PeptideIndexing {
    PeptideIndexing {
        decoy_rule: DecoyRule::Prefix("DECOY_".into()),
        missing_decoy_action: MissingDecoyAction::Silent,
        unmatched_action: UnmatchedAction::Warn,
        ..Default::default()
    }
}

#[test]
fn independent_pinned_matcher_oracle_121_cases() {
    let fixture = include_str!("data/peptide_indexing_matcher_reference.tsv");
    let mut count = 0;
    for row in fixture.lines().skip(1) {
        let fields: Vec<_> = row.split('\t').collect();
        let indexer = PeptideIndexing {
            max_ambiguities: fields[3].parse().unwrap(),
            max_mismatches: fields[4].parse().unwrap(),
            il_equivalent: fields[5].parse().unwrap(),
            ..Default::default()
        };
        let expected: Vec<usize> = if fields[6] == "-" {
            Vec::new()
        } else {
            fields[6].split(',').map(|s| s.parse().unwrap()).collect()
        };
        assert_eq!(
            indexer.find_matches(fields[1], fields[2]).unwrap(),
            expected,
            "{}",
            fields[0]
        );
        count += 1;
    }
    assert_eq!(count, 121);
}

#[test]
fn source_decoy_inference_regex_and_threshold_cases() {
    // Hand-traced source cases, independently reviewed in fixture provenance.
    let cases: &[(&[&str], &str, bool, bool)] = &[
        (&["Protein1", "DECOY_Protein2"], "DECOY_", true, true),
        (&["Protein1", "DECOYProtein2"], "DECOY", true, true),
        (&["Protein1", "Protein2DECOY_"], "DECOY_", true, false),
        (&["Protein1", "Protein2_DECOY"], "_DECOY", false, true),
        (&["Protein1", "reverse_Protein"], "reverse_", true, true),
        (&["rev_Protein1", "reverse_Protein"], "DECOY_", true, false),
        (&["A", "B", "C", "DECOY_D", "DECOY_E"], "DECOY_", true, true),
        (&["A", "B", "C", "D", "DECOY_E"], "DECOY_", true, false),
        (
            &["DECOY_A", "DECOY_B", "DECOY_C", "DECOY_D", "rev_E"],
            "DECOY_",
            true,
            true,
        ),
        (
            &["DECOY_A", "DECOY_B", "DECOY_C", "rev_D"],
            "DECOY_",
            true,
            false,
        ),
        (&["DECOY_A", "B_DECOY"], "DECOY_", true, false),
        (
            &["DECOY_A", "DECOY_B_DECOY", "C_DECOY", "D_DECOY"],
            "DECOY_",
            true,
            true,
        ),
        (&["DECOY_A", "decoy_B"], "decoy_", true, true),
        (&["A", "reversed_B"], "reverse", true, true),
        (&["A", "B_deco"], "_deco", false, true),
        (&["A", "B_decoyyy"], "_decoyyy", false, true),
        (&["A", "B_decoy_"], "DECOY_", true, false),
    ];
    for (names, affix, prefix, inferred) in cases {
        let database: Vec<_> = names.iter().map(|n| entry(n, "A")).collect();
        let resolved = PeptideIndexing::default()
            .resolve_decoy_rule(&database)
            .unwrap();
        assert_eq!(
            (&resolved.affix[..], resolved.is_prefix, resolved.inferred),
            (*affix, *prefix, *inferred),
            "{names:?}"
        );
    }
    let resolved = PeptideIndexing::default()
        .resolve_decoy_rule(&[entry("DECOY_A", "A"), entry("decoy_B", "A")])
        .unwrap();
    assert!(!resolved.is_decoy("DECOY_A"));
    assert!(resolved.is_decoy("decoy_B"));
    let explicit = PeptideIndexing {
        decoy_rule: DecoyRule::Suffix("_reverse".into()),
        ..Default::default()
    }
    .resolve_decoy_rule(&[])
    .unwrap();
    assert!(explicit.is_decoy("A_reverse"));
    assert!(!explicit.inferred);
}

#[test]
fn source_enzyme_boundary_oracle() {
    // (protein, start, length, enzyme, M/MX clipping, X!Tandem D|P, valid)
    let cases = [
        ("MPEPTIDER", 1, 8, "Trypsin", true, false, true),
        ("MPEPTIDER", 1, 8, "Trypsin", false, false, false),
        ("MAPEPTIDER", 2, 8, "Trypsin", true, false, true),
        ("MAAPEPTIDER", 3, 8, "Trypsin", true, false, false),
        ("AAPEPTIDER", 2, 8, "Trypsin", true, false, false),
        ("MKDPLMMLK", 1, 8, "no cleavage", true, false, true),
        ("MKDPLMMLK", 1, 8, "no cleavage", false, false, false),
        ("KDPLMMLK", 0, 8, "no cleavage", false, false, true),
        ("KDPLMMLK", 0, 2, "no cleavage", false, true, true),
        ("KDPLMMLK", 0, 2, "no cleavage", false, false, false),
        ("ADPEPTIDER", 2, 8, "Trypsin", true, true, true),
        ("ADPEPTIDER", 2, 8, "Trypsin", true, false, false),
        ("AKPEPTIDER", 2, 8, "Trypsin", false, false, false),
        ("AKPEPTIDER", 2, 8, "Trypsin/P", false, false, true),
        ("AKRPEPTIDER", 0, 11, "Trypsin", false, false, true),
    ];
    for (protein, start, length, enzyme, clipping, dp, expected) in cases {
        let config = PeptideIndexing {
            enzyme: Some(Protease::from_name(enzyme).unwrap()),
            allow_nterm_protein_cleavage: clipping,
            max_ambiguities: 0,
            ..permissive()
        };
        let mut proteins = [run("r")];
        if dp {
            proteins[0].search_engine = "XTandem".into();
        }
        let mut peptides = [identification("r", &protein[start..start + length])];
        config
            .run(&[entry("A", protein)], &mut proteins, &mut peptides)
            .unwrap();
        let evidence = &peptides[0].hits[0].evidences;
        assert_eq!(
            !evidence.is_empty(),
            expected,
            "{protein}, {start}, {length}, {enzyme}, {dp}"
        );
        if expected {
            assert_eq!(evidence[0].start, Some(start));
        }
    }
}

#[test]
fn modified_peptides_all_occurrences_order_and_run_reconstruction() {
    let database = [
        entry("B", "ACMKACMK"),
        entry("DECOY_A", "ACMK"),
        entry("C", "PEPTIDER"),
    ];
    let mut proteins = [run("one"), run("two")];
    proteins[0]
        .metadata
        .insert("keep".into(), "run metadata".into());
    proteins[0].hits = vec![ProteinHit::new(99., 4, "B", "ACMKACMK").unwrap()];
    let mut peptides = [
        identification("one", "(Acetyl)AC(Carbamidomethyl)M(Oxidation)K"),
        identification("two", "PEPTIDER"),
    ];
    peptides[0].hits[0]
        .metadata
        .insert("custom".into(), 42_i64.into());
    let original_sequence = peptides[0].hits[0].sequence.clone();
    let config = PeptideIndexing {
        write_protein_sequence: true,
        write_protein_description: true,
        decoy_rule: DecoyRule::Prefix("DECOY_".into()),
        ..Default::default()
    };
    let report = config.run(&database, &mut proteins, &mut peptides).unwrap();
    assert_eq!(report.target_and_decoy_hits, 1);
    assert_eq!(report.target_hits, 1);
    assert_eq!(report.non_unique_hits, 1);
    assert_eq!(report.unique_hits, 1);
    assert_eq!(report.evidence_count, 4);
    let hit = &peptides[0].hits[0];
    assert_eq!(hit.sequence, original_sequence);
    assert_eq!((hit.score, hit.rank), (12.5, 7));
    assert_eq!(hit.metadata["custom"].to_string(), "42");
    assert_eq!(
        hit.target_decoy_type().unwrap(),
        TargetDecoyType::TargetAndDecoy
    );
    assert_eq!(
        hit.evidences
            .iter()
            .map(|e| (e.protein_accession.as_str(), e.start, e.end))
            .collect::<Vec<_>>(),
        vec![
            ("B", Some(0), Some(3)),
            ("B", Some(4), Some(7)),
            ("DECOY_A", Some(0), Some(3))
        ]
    );
    assert_eq!(hit.evidences[0].aa_before, FlankingResidue::NTerminus);
    assert_eq!(hit.evidences[0].aa_after, FlankingResidue::Residue('A'));
    assert_eq!(hit.evidences[1].aa_after, FlankingResidue::CTerminus);
    assert_eq!(
        proteins[0]
            .hits
            .iter()
            .map(|h| &h.accession)
            .collect::<Vec<_>>(),
        vec!["B", "DECOY_A"]
    );
    assert_eq!(proteins[1].hits[0].accession, "C");
    assert_eq!(proteins[0].hits[0].score, 0.0); // source reconstructs fresh hits
    assert_eq!(proteins[0].hits[0].description(), "description B");
    assert_eq!(proteins[0].hits[0].sequence, "ACMKACMK");
    assert_eq!(
        proteins[0].search_parameters.metadata["PeptideIndexer:enzyme"].to_string(),
        "Trypsin"
    );
    assert!(proteins[0].metadata.contains_key("keep"));
}

#[test]
fn repeated_occurrences_in_one_protein_are_unique() {
    let mut proteins = [run("r")];
    let mut peptides = [identification("r", "AA")];
    let report = PeptideIndexing {
        specificity: Some(DigestionSpecificity::None),
        ..permissive()
    }
    .run(&[entry("A", "AAAA")], &mut proteins, &mut peptides)
    .unwrap();
    assert_eq!(report.evidence_count, 3);
    assert_eq!(report.unique_hits, 1);
    assert_eq!(
        peptides[0].hits[0].metadata["protein_references"].to_string(),
        "unique"
    );
}

#[test]
fn source_ambiguities_stop_removal_and_original_sequence_output() {
    let mut proteins = [run("r")];
    let mut peptides = [identification("r", "MLTEAEK")];
    let config = PeptideIndexing {
        max_ambiguities: 1,
        write_protein_sequence: true,
        ..permissive()
    };
    config
        .run(&[entry("A", "*MLT*EAXK")], &mut proteins, &mut peptides)
        .unwrap();
    assert_eq!(peptides[0].hits[0].evidences[0].start, Some(0));
    assert_eq!(peptides[0].hits[0].evidences[0].end, Some(6));
    assert_eq!(proteins[0].hits[0].sequence, "*MLT*EAXK");
    let peptides_b = ["NENE", "NEDE", "DENE", "DEDE"];
    for peptide in peptides_b {
        assert!(
            PeptideIndexing {
                max_ambiguities: 1,
                ..Default::default()
            }
            .find_matches(peptide, "B*EBE*")
            .unwrap()
            .is_empty()
        );
        assert_eq!(
            PeptideIndexing {
                max_ambiguities: 2,
                ..Default::default()
            }
            .find_matches(peptide, "B*EBE*")
            .unwrap(),
            vec![0]
        );
    }
}

#[test]
fn native_case_normalization_and_il_flanks() {
    let config = PeptideIndexing {
        il_equivalent: true,
        specificity: Some(DigestionSpecificity::None),
        ..permissive()
    };
    assert_eq!(config.find_matches("il", "lj").unwrap(), vec![0]);
    assert!(config.find_matches("J", "J").unwrap().is_empty()); // source peptide J stays literal
    let mut proteins = [run("r")];
    let mut peptides = [identification("r", "AC")];
    config
        .run(&[entry("A", "lacj")], &mut proteins, &mut peptides)
        .unwrap();
    let evidence = &peptides[0].hits[0].evidences[0];
    assert_eq!(evidence.aa_before, FlankingResidue::Residue('I'));
    assert_eq!(evidence.aa_after, FlankingResidue::Residue('I'));
}

#[test]
fn missing_decoy_and_unmatched_errors_are_transactional() {
    let mut proteins = [run("r")];
    let mut peptides = [identification("r", "PEPTIDER")];
    let before = (proteins.clone(), peptides.clone());
    assert!(
        PeptideIndexing::default()
            .run(&[entry("A", "PEPTIDER")], &mut proteins, &mut peptides)
            .is_err()
    );
    assert_eq!((proteins.clone(), peptides.clone()), before);
    let config = PeptideIndexing {
        unmatched_action: UnmatchedAction::Error,
        ..permissive()
    };
    assert!(
        config
            .run(&[entry("A", "ACMK")], &mut proteins, &mut peptides)
            .is_err()
    );
    assert_eq!((proteins, peptides), before);
}

#[test]
fn unmatched_warn_and_remove_clear_stale_evidence_and_labels() {
    for action in [UnmatchedAction::Warn, UnmatchedAction::Remove] {
        let config = PeptideIndexing {
            unmatched_action: action,
            ..permissive()
        };
        let mut proteins = [run("r")];
        let mut peptides = [identification("r", "PEPTIDER")];
        peptides[0].hits[0].set_target_decoy_type(TargetDecoyType::Decoy);
        peptides[0].hits[0]
            .evidences
            .push(openms::identification::PeptideEvidence::new("old", 0..=7).unwrap());
        let report = config
            .run(&[entry("A", "ACMK")], &mut proteins, &mut peptides)
            .unwrap();
        assert_eq!(report.unmatched_hits, 1);
        if action == UnmatchedAction::Remove {
            assert!(peptides[0].hits.is_empty());
        } else {
            assert_eq!(
                peptides[0].hits[0].target_decoy_type().unwrap(),
                TargetDecoyType::Unknown
            );
            assert!(peptides[0].hits[0].evidences.is_empty());
            assert_eq!(report.warnings.len(), 1);
        }
    }
}

#[test]
fn unreferenced_proteins_are_retained_only_in_their_own_run() {
    let mut proteins = [run("one"), run("two")];
    for protein in &mut proteins {
        let mut hit = ProteinHit::new(55., 1, "A", "PEPTIDER").unwrap();
        hit.set_target_decoy_type(TargetDecoyType::Target).unwrap();
        protein.hits.push(hit);
    }
    let mut peptides = [identification("one", "PEPTIDER")];
    PeptideIndexing {
        keep_unreferenced_proteins: true,
        ..permissive()
    }
    .run(&[entry("A", "PEPTIDER")], &mut proteins, &mut peptides)
    .unwrap();
    assert_eq!(proteins[0].hits[0].score, 0.);
    assert_eq!(proteins[1].hits[0].score, 55.);
    assert_eq!(
        proteins[1].hits[0].target_decoy_type().unwrap(),
        TargetDecoyType::Unknown
    );
}

#[test]
fn per_run_engine_and_specificity_resolution() {
    let mut proteins = [run("normal"), run("msgf"), run("nonspecific")];
    proteins[1].search_engine = "MS-GF+".into();
    proteins[2].search_parameters.enzyme_specificity = EnzymeTermSpecificity::None;
    let mut peptides = [
        identification("normal", "PEPTIDER"),
        identification("msgf", "PEPTIDER"),
        identification("nonspecific", "EPTIDE"),
    ];
    let report = permissive()
        .run(&[entry("A", "AKPEPTIDER")], &mut proteins, &mut peptides)
        .unwrap();
    assert!(peptides[0].hits[0].evidences.is_empty());
    assert_eq!(peptides[1].hits[0].evidences[0].start, Some(2));
    assert_eq!(peptides[2].hits[0].evidences[0].start, Some(3));
    assert_eq!(report.runs[0].enzyme, Protease::Trypsin);
    assert_eq!(report.runs[1].enzyme, Protease::TrypsinP);
    assert_eq!(report.runs[2].specificity, DigestionSpecificity::None);
    assert_eq!(
        proteins[0].search_parameters.metadata["PeptideIndexer:enzyme"].to_string(),
        "Trypsin"
    );
    assert_eq!(
        proteins[1].search_parameters.metadata["PeptideIndexer:enzyme"].to_string(),
        "Trypsin/P"
    );
}

#[test]
fn input_and_configuration_rejections_preserve_state() {
    for database in [
        vec![],
        vec![entry("A", "ACMK"), entry("A", "ACMK")],
        vec![entry("", "ACMK")],
        vec![entry("A", "A?CMK")],
    ] {
        let mut proteins = [run("r")];
        let mut peptides = [identification("r", "ACMK")];
        let before = (proteins.clone(), peptides.clone());
        assert!(
            permissive()
                .run(&database, &mut proteins, &mut peptides)
                .is_err()
        );
        assert_eq!((proteins, peptides), before);
    }
    for raw in ["A?C", "A-C", "A.C", "A$C", "A\0C", "A(C", "αAC"] {
        assert!(permissive().find_matches("AC", raw).is_err());
    }
    assert!(permissive().find_matches("", "AC").is_err());
    assert_eq!(permissive().find_matches("A*C", "AC").unwrap(), vec![0]);
    assert!(permissive().find_matches("***", "AC").is_err());
    for config in [
        PeptideIndexing {
            max_ambiguities: 11,
            ..permissive()
        },
        PeptideIndexing {
            max_mismatches: 11,
            ..permissive()
        },
        PeptideIndexing {
            max_work: 0,
            ..permissive()
        },
        PeptideIndexing {
            decoy_rule: DecoyRule::Prefix(String::new()),
            ..permissive()
        },
    ] {
        assert!(config.find_matches("A", "A").is_err());
    }
    let mut proteins = [run("r")];
    let mut peptides = [identification("missing", "AC")];
    assert!(
        permissive()
            .run(&[entry("A", "AC")], &mut proteins, &mut peptides)
            .is_err()
    );
    let mut proteins = [run("r"), run("r")];
    let mut peptides = [identification("r", "AC")];
    assert!(
        permissive()
            .run(&[entry("A", "AC")], &mut proteins, &mut peptides)
            .is_err()
    );
}

#[test]
fn resource_bounds_and_output_duplication_are_transactional() {
    for config in [
        PeptideIndexing {
            max_work: 2,
            ..permissive()
        },
        PeptideIndexing {
            max_residues: 2,
            ..permissive()
        },
        PeptideIndexing {
            max_records: 1,
            ..permissive()
        },
        PeptideIndexing {
            max_matches: 1,
            ..permissive()
        },
    ] {
        let mut proteins = [run("r")];
        let mut peptides = [identification("r", "ACMK"), identification("r", "ACMK")];
        let before = (proteins.clone(), peptides.clone());
        assert!(
            config
                .run(&[entry("A", "ACMK")], &mut proteins, &mut peptides)
                .is_err()
        );
        assert_eq!((proteins, peptides), before);
    }
    let config = PeptideIndexing {
        max_matches: 2,
        ..permissive()
    };
    assert!(config.find_matches("AA", "AAAA").is_err());
}

#[test]
fn unknown_settings_fallback_and_unsupported_regex_are_explicit() {
    let mut proteins = [ProteinIdentification {
        identifier: "r".into(),
        ..Default::default()
    }];
    let mut peptides = [identification("r", "ACMK")];
    let report = permissive()
        .run(&[entry("A", "ACMK")], &mut proteins, &mut peptides)
        .unwrap();
    assert_eq!(report.runs[0].enzyme, Protease::Trypsin);
    assert_eq!(report.warnings.len(), 2);
    proteins[0].search_parameters.digestion_regex = "custom".into();
    let before = (proteins.clone(), peptides.clone());
    assert!(
        permissive()
            .run(&[entry("A", "ACMK")], &mut proteins, &mut peptides)
            .is_err()
    );
    assert_eq!((proteins.clone(), peptides.clone()), before);
    PeptideIndexing {
        enzyme: Some(Protease::Trypsin),
        ..permissive()
    }
    .run(&[entry("A", "ACMK")], &mut proteins, &mut peptides)
    .unwrap();
    let config = PeptideIndexing {
        enzyme: Some(Protease::Chymotrypsin),
        il_equivalent: true,
        ..permissive()
    };
    assert!(
        config
            .run(&[entry("A", "ACMK")], &mut proteins, &mut peptides)
            .is_err()
    );
}

#[test]
fn empty_hits_are_a_successful_consistent_cleanup() {
    for empty_record in [false, true] {
        let mut proteins = [run("r")];
        proteins[0]
            .hits
            .push(ProteinHit::new(2., 1, "old", "AC").unwrap());
        let mut peptides = if empty_record {
            vec![PeptideIdentification {
                identifier: "r".into(),
                ..Default::default()
            }]
        } else {
            vec![]
        };
        let report = PeptideIndexing::default()
            .run(&[entry("A", "AC")], &mut proteins, &mut peptides)
            .unwrap();
        assert_eq!(report.peptide_hits, 0);
        assert!(proteins[0].hits.is_empty());
    }
}

#[test]
fn wrapped_search_engines_restore_source_special_cleavages() {
    for wrapper in ["Percolator", "ConsensusID"] {
        let mut proteins = [run("r")];
        proteins[0].search_engine = wrapper.into();
        proteins[0]
            .search_parameters
            .metadata
            .insert("SE:MSGFPLUS".into(), "true".into());
        let mut peptides = [identification("r", "PEPTIDER")];
        let report = permissive()
            .run(&[entry("A", "AKPEPTIDER")], &mut proteins, &mut peptides)
            .unwrap();
        assert_eq!(report.runs[0].enzyme, Protease::TrypsinP);
        assert_eq!(peptides[0].hits[0].evidences[0].start, Some(2));
        assert_eq!(
            proteins[0].search_parameters.metadata["PeptideIndexer:enzyme"].to_string(),
            "Trypsin/P"
        );
        assert!(!proteins[0].metadata.contains_key("PeptideIndexer:enzyme"));
        proteins[0].search_parameters.metadata.remove("SE:MSGFPLUS");
        proteins[0]
            .search_parameters
            .metadata
            .insert("SE:XTANDEM".into(), "true".into());
        let report = permissive()
            .run(&[entry("A", "ADPEPTIDER")], &mut proteins, &mut peptides)
            .unwrap();
        assert!(report.runs[0].allow_random_asp_pro_cleavage);
        assert_eq!(peptides[0].hits[0].evidences[0].start, Some(2));
    }
}

#[test]
fn semi_specificity_and_missing_decoy_warning() {
    let config = PeptideIndexing {
        specificity: Some(DigestionSpecificity::Semi),
        missing_decoy_action: MissingDecoyAction::Warn,
        ..permissive()
    };
    let mut proteins = [run("r")];
    let mut peptides = [
        identification("r", "PEPTIDER"),
        identification("r", "EPTIDE"),
    ];
    let report = config
        .run(&[entry("A", "AAPEPTIDER")], &mut proteins, &mut peptides)
        .unwrap();
    assert_eq!(peptides[0].hits[0].evidences.len(), 1);
    assert!(peptides[1].hits[0].evidences.is_empty());
    assert!(
        report
            .warnings
            .iter()
            .any(|s| s.contains("no peptide hit maps to a decoy"))
    );
    assert_eq!(
        proteins[0]
            .search_parameters
            .metadata
            .iter()
            .filter(|(key, _)| key.starts_with("PeptideIndexer:"))
            .count(),
        10
    );
}
