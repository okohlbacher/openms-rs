// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source: pinned METADATA identification implementations and class tests.

use openms::chemistry::AASequence;
use openms::comparison::Tolerance;
use openms::identification::*;
use openms::kernel::{BaseFeature, ConsensusMap, DataArray, Feature, FeatureMap};
use openms::processing::{Normalizer, SpectrumFilter};
use openms::{Error, MSExperiment, MSSpectrum, Peak1D};
use std::collections::BTreeSet;

fn hit(sequence: &str, score: f64) -> PeptideHit {
    PeptideHit::new(score, 0, 2, AASequence::parse(sequence).unwrap()).unwrap()
}
fn identification(hits: Vec<PeptideHit>) -> PeptideIdentification {
    PeptideIdentification {
        hits,
        ..Default::default()
    }
}
fn evidence(accession: &str, start: usize, end: usize) -> PeptideEvidence {
    PeptideEvidence::new(accession, start..=end).unwrap()
}
fn protein(accession: &str, sequence: &str) -> ProteinHit {
    ProteinHit::new(1.0, 0, accession, sequence).unwrap()
}

#[test]
fn defaults_express_missing_values_without_sentinels() {
    let peptide = PeptideIdentification::default();
    assert_eq!((peptide.rt, peptide.mz), (None, None));
    assert!(peptide.higher_score_better);
    assert!(peptide.is_empty());
    assert_eq!(peptide.significance_threshold, 0.0);
    assert_eq!(PeptideHit::default().rank, 0);
    assert_eq!(ProteinHit::default().coverage, None);
    let search = SearchParameters::default();
    assert_eq!(search.charge_range().unwrap(), None);
    assert_eq!(search.digestion_enzyme, "unknown_enzyme");
    assert_eq!(search.enzyme_specificity, EnzymeTermSpecificity::Unknown);
    assert_eq!(search.fragment_tolerance, Tolerance::Absolute(0.0));
    peptide.validate().unwrap();
    ProteinIdentification::default().validate().unwrap();
}

#[test]
fn evidence_uses_inclusive_zero_based_positions_and_source_marker_order() {
    assert!(evidence("P1", 0, 0).has_valid_limits());
    assert_eq!(evidence("P1", 0, 4).positions().unwrap().count(), 5);
    assert!(!PeptideEvidence::default().has_valid_limits());
    assert!(PeptideEvidence::default().positions().is_err());
    let mut reversed = evidence("P1", 0, 2);
    reversed.start = Some(3);
    assert!(reversed.validate().is_err());
    assert!(FlankingResidue::from_code('a').is_err());
    let mut values: Vec<_> = [']', 'X', '[', 'A']
        .into_iter()
        .map(|c| {
            let mut item = evidence("P1", 0, 4);
            item.aa_before = FlankingResidue::from_code(c).unwrap();
            item
        })
        .collect();
    values.sort();
    assert_eq!(
        values
            .iter()
            .map(|e| e.aa_before.code())
            .collect::<String>(),
        "AX[]"
    );
    let mut missing = evidence("P1", 0, 4);
    missing.start = None;
    assert!(missing < values[0]);
}

#[test]
fn stable_score_sorting_and_best_hit_respect_direction_and_ties() {
    let mut id = identification(vec![hit("PEP", 1.0), hit("AAA", 3.0), hit("CCC", 3.0)]);
    assert_eq!(id.best_hit().unwrap().unwrap().sequence.as_str(), "AAA");
    id.sort().unwrap();
    assert_eq!(
        id.hits
            .iter()
            .map(|h| h.sequence.as_str())
            .collect::<Vec<_>>(),
        ["AAA", "CCC", "PEP"]
    );
    id.higher_score_better = false;
    id.sort().unwrap();
    assert_eq!(
        id.hits
            .iter()
            .map(|h| h.sequence.as_str())
            .collect::<Vec<_>>(),
        ["PEP", "AAA", "CCC"]
    );
    assert_eq!(id.best_hit().unwrap().unwrap().score, 1.0);
    let mut proteins = ProteinIdentification {
        hits: vec![
            protein("P1", "AAA"),
            protein("P2", "AAA"),
            protein("P3", "AAA"),
        ],
        ..Default::default()
    };
    proteins.hits[2].score = 2.0;
    proteins.sort().unwrap();
    assert_eq!(
        proteins
            .hits
            .iter()
            .map(|h| h.accession.as_str())
            .collect::<Vec<_>>(),
        ["P3", "P1", "P2"]
    );
    assert_eq!(proteins.find_hit("P2").unwrap().sequence, "AAA");
}

#[test]
fn invalid_scores_fail_before_sorting_or_filtering() {
    let mut id = identification(vec![hit("PEP", 1.0), hit("AAA", 2.0)]);
    id.hits[1].score = f64::INFINITY;
    let original = id.clone();
    assert!(id.sort().is_err());
    assert!(id.best_hit().is_err());
    assert_eq!(id, original);
    let mut experiment = MSExperiment {
        spectra: vec![MSSpectrum {
            peaks: vec![Peak1D::new(100.0, 100.0)],
            peptide_identifications: vec![id],
            ..Default::default()
        }],
        ..Default::default()
    };
    let original = experiment.clone();
    assert!(
        Normalizer::default()
            .filter_experiment(&mut experiment)
            .is_err()
    );
    assert_eq!(experiment, original);
}

#[test]
fn identity_keys_include_modifications_and_charge() {
    let unmodified = hit("MPEP", 1.0);
    let oxidized = hit("M(Oxidation)PEP", 2.0);
    let mut charged = oxidized.clone();
    charged.charge = 3;
    assert!(!unmodified.same_sequence_and_charge(&oxidized));
    assert!(!oxidized.same_sequence_and_charge(&charged));
    assert_ne!(oxidized.identity_key(), charged.identity_key());
    let mut new_score = oxidized.clone();
    new_score.score = 12.0;
    assert!(oxidized.same_sequence_and_charge(&new_score));
    assert_eq!(oxidized.identity_key(), new_score.identity_key());
}

#[test]
fn accessions_references_and_typed_metadata_remain_attached() {
    let mut candidate = hit("PEP", 1.0);
    candidate.evidences = vec![
        evidence("P2", 0, 2),
        evidence("P1", 4, 6),
        evidence("P2", 3, 5),
        evidence("", 0, 2),
    ];
    assert_eq!(candidate.protein_accessions(), BTreeSet::from(["P1", "P2"]));
    let mut id = identification(vec![candidate, hit("AAA", 2.0)]);
    assert_eq!(id.referencing_hits(&BTreeSet::from(["P2".into()])).len(), 1);
    id.metadata
        .insert("spectrum_reference".into(), 17_i64.into());
    assert_eq!(id.spectrum_reference(), "17");
    id.set_spectrum_reference("scan=25");
    assert_eq!(id.spectrum_reference(), "scan=25");
    id.set_experiment_label("sample A");
    assert_eq!(id.experiment_label(), "sample A");
    id.set_experiment_label("");
    assert!(!id.metadata.contains_key("experiment_label"));
}

#[test]
fn target_decoy_categories_distinguish_mixed_peptide_evidence() {
    let mut peptide = hit("PEP", 1.0);
    assert_eq!(
        peptide.target_decoy_type().unwrap(),
        TargetDecoyType::Unknown
    );
    peptide.set_target_decoy_type(TargetDecoyType::TargetAndDecoy);
    assert!(!peptide.is_decoy().unwrap());
    peptide.set_target_decoy_type(TargetDecoyType::Decoy);
    assert!(peptide.is_decoy().unwrap());
    peptide
        .metadata
        .insert("target_decoy".into(), "TARGET".into());
    assert_eq!(
        peptide.target_decoy_type().unwrap(),
        TargetDecoyType::Target
    );
    let mut protein = protein("P1", "PEP");
    protein
        .set_target_decoy_type(TargetDecoyType::Decoy)
        .unwrap();
    let old = protein.clone();
    assert!(
        protein
            .set_target_decoy_type(TargetDecoyType::TargetAndDecoy)
            .is_err()
    );
    assert_eq!(protein, old);
    protein
        .set_target_decoy_type(TargetDecoyType::Unknown)
        .unwrap();
    assert!(!protein.metadata.contains_key("target_decoy"));
    protein
        .metadata
        .insert("target_decoy".into(), "other".into());
    assert!(protein.validate().is_err());
}

#[test]
fn annotations_and_protein_groups_follow_source_comparison_keys() {
    let mut annotations = vec![
        PeakAnnotation {
            mz: 100.0,
            intensity: 5.0,
            charge: 2,
            annotation: "a1".into(),
        },
        PeakAnnotation {
            mz: 100.0,
            intensity: 8.0,
            charge: 1,
            annotation: "b1".into(),
        },
        PeakAnnotation {
            mz: 100.0,
            intensity: 9.0,
            charge: 1,
            annotation: "a1".into(),
        },
        PeakAnnotation {
            mz: 100.0,
            intensity: 3.0,
            charge: 1,
            annotation: "a1".into(),
        },
    ];
    PeakAnnotation::sort(&mut annotations).unwrap();
    assert_eq!(
        annotations.iter().map(|a| a.intensity).collect::<Vec<_>>(),
        [3.0, 9.0, 8.0, 5.0]
    );
    let group = |p, names: &[&str]| ProteinGroup {
        probability: p,
        accessions: names.iter().map(|s| (*s).into()).collect(),
        ..Default::default()
    };
    let mut groups = vec![
        group(0.5, &["P1"]),
        group(0.9, &["P1", "P2"]),
        group(0.9, &["P2"]),
        group(0.9, &["P1"]),
    ];
    ProteinGroup::sort(&mut groups).unwrap();
    assert_eq!(
        groups
            .iter()
            .map(|g| (g.probability, g.accessions.len(), g.accessions[0].as_str()))
            .collect::<Vec<_>>(),
        [
            (0.9, 1, "P1"),
            (0.9, 1, "P2"),
            (0.9, 2, "P1"),
            (0.5, 1, "P1")
        ]
    );
    let original = groups[0].clone();
    groups[0]
        .float_data_arrays
        .push(DataArray::new("samples", vec![10.0, 20.0, 30.0]));
    groups[0].validate().unwrap(); // Three samples and one accession is valid.
    assert_ne!(groups[0], original);
}

#[test]
fn charge_ranges_include_upstream_goldens_and_signed_boundaries() {
    for (text, expected) in [
        ("1,2,3", (1, 3)),
        ("+2-+5", (2, 5)),
        ("-1,-2,-3", (-3, -1)),
        ("2", (2, 2)),
        (" -3--1 ", (-3, -1)),
        ("2+:5+", (2, 5)),
        ("3-:1-", (-3, -1)),
        ("-2147483648:2147483647", (i32::MIN, i32::MAX)),
        ("2147483648-", (i32::MIN, i32::MIN)),
    ] {
        let params = SearchParameters {
            charges: text.into(),
            ..Default::default()
        };
        assert_eq!(params.charge_range().unwrap(), Some(expected), "{text}");
    }
    for text in [
        "1,",
        "a",
        "3:1",
        "1:2:3",
        "1,2:3",
        "2147483648",
        "--1",
        "2--3",
    ] {
        assert!(
            SearchParameters {
                charges: text.into(),
                ..Default::default()
            }
            .charge_range()
            .is_err(),
            "{text}"
        );
    }
}

#[test]
fn mergeability_preserves_source_database_and_label_rules() {
    let a = SearchParameters {
        database: "/data/proteins.fasta".into(),
        charges: "2,3".into(),
        fixed_modifications: vec!["A".into(), "B".into()],
        ..Default::default()
    };
    let mut b = a.clone();
    b.database = "C:\\data\\proteins.fasta".into();
    b.fixed_modifications = vec!["B".into(), "A".into(), "A".into()];
    b.mass_type = PeakMassType::Average;
    b.missed_cleavages = 99;
    assert!(a.mergeable(&b, "label-free").unwrap());
    b.variable_modifications.push("C".into());
    assert!(!a.mergeable(&b, "label-free").unwrap());
    assert!(a.mergeable(&b, "labeled_MS1").unwrap());
    b.fragment_tolerance = Tolerance::Ppm(0.0);
    assert!(!a.mergeable(&b, "labeled_MS1").unwrap());
}

#[test]
fn coverage_unions_overlaps_duplicates_and_adjacent_inclusive_intervals() {
    let mut proteins = ProteinIdentification {
        identifier: "run-A".into(),
        hits: vec![protein("P1", "ACDEFGHIKL"), protein("P2", "AA")],
        ..Default::default()
    };
    let mut candidate = hit("ACD", 1.0);
    candidate.evidences = vec![
        evidence("P1", 0, 2),
        evidence("P1", 2, 5),
        evidence("P1", 0, 2),
        evidence("P1", 8, 9),
    ];
    let mut id = identification(vec![candidate]);
    id.identifier = "run-B".into(); // Source explicitly does not filter run identifiers.
    proteins.compute_coverage(&[id.clone()]).unwrap();
    assert_eq!(proteins.hits[0].coverage, Some(80.0));
    assert_eq!(proteins.hits[1].coverage, Some(0.0));
    id.hits[0].evidences.push(evidence("P1", 6, 7));
    id.hits[0].evidences.push(evidence("P2", 0, 0));
    proteins.compute_coverage(&[id]).unwrap();
    assert_eq!(proteins.hits[0].coverage, Some(100.0));
    assert_eq!(proteins.hits[1].coverage, Some(50.0));
}

#[test]
fn coverage_matches_independent_bitset_oracle_for_every_small_interval() {
    for length in 1..=12 {
        for seed in 0..9 {
            let mut expected = vec![false; length];
            let mut candidate = hit("A", 1.0);
            for start in 0..length {
                for end in start..length {
                    if (start * 7 + end * 11 + seed) % 9 == 0 {
                        expected[start..=end].fill(true);
                        candidate.evidences.push(evidence("P1", start, end));
                    }
                }
            }
            let mut proteins = ProteinIdentification {
                hits: vec![protein("P1", &"A".repeat(length))],
                ..Default::default()
            };
            proteins
                .compute_coverage(&[identification(vec![candidate])])
                .unwrap();
            assert_eq!(
                proteins.hits[0].coverage,
                Some(100.0 * expected.iter().filter(|&&b| b).count() as f64 / length as f64)
            );
        }
    }
}

#[test]
fn coverage_rejects_unknown_or_out_of_bounds_evidence_atomically() {
    let mut proteins = ProteinIdentification {
        hits: vec![protein("P1", "AAA"), protein("P2", "AAA")],
        ..Default::default()
    };
    proteins.hits[0].coverage = Some(12.0);
    let mut candidate = hit("A", 1.0);
    candidate.evidences = vec![evidence("P1", 0, 0), evidence("P2", 2, 3)];
    let original = proteins.clone();
    assert!(
        proteins
            .compute_coverage(&[identification(vec![candidate.clone()])])
            .is_err()
    );
    assert_eq!(proteins, original); // end == length is the upstream out-of-bounds defect.
    candidate.evidences[1].end = None;
    assert!(
        proteins
            .compute_coverage(&[identification(vec![candidate])])
            .is_err()
    );
    assert_eq!(proteins, original);
    proteins.hits[1].sequence.clear();
    assert!(proteins.compute_coverage(&[]).is_err());
}

#[test]
fn observed_modifications_map_terminals_and_residues_and_apply_skip_sets() {
    let mut candidate = hit(".(Acetyl)AM(Oxidation)C(Carbamidomethyl).(Amidated)", 1.0);
    candidate.evidences = vec![evidence("P1", 2, 4), evidence("P1", 2, 4)];
    let id = identification(vec![candidate]);
    let mut proteins = ProteinIdentification {
        hits: vec![protein("P1", "AAAMCAA"), protein("P2", "AAA")],
        ..Default::default()
    };
    proteins
        .compute_modifications(std::slice::from_ref(&id), &BTreeSet::new())
        .unwrap();
    let observed: Vec<_> = proteins.hits[0]
        .modifications
        .iter()
        .map(|m| (m.position, m.modification.name()))
        .collect();
    assert_eq!(
        observed,
        [
            (2, "Acetyl"),
            (3, "Oxidation"),
            (4, "Amidated"),
            (4, "Carbamidomethyl")
        ]
    );
    proteins
        .compute_modifications(
            &[id],
            &BTreeSet::from(["Carbamidomethyl (C)".into(), "Acetyl".into()]),
        )
        .unwrap();
    assert_eq!(
        proteins.hits[0]
            .modifications
            .iter()
            .map(|m| m.modification.name())
            .collect::<Vec<_>>(),
        ["Oxidation", "Amidated"]
    );
    let old = proteins.hits.clone();
    proteins
        .compute_modifications(&[], &BTreeSet::new())
        .unwrap();
    assert_eq!(proteins.hits, old);
}

#[test]
fn bad_modified_evidence_leaves_all_protein_modifications_unchanged() {
    let mut proteins = ProteinIdentification {
        hits: vec![protein("P1", "AMA")],
        ..Default::default()
    };
    let mut candidate = hit("AM(Oxidation)A", 1.0);
    candidate.evidences = vec![evidence("P1", 0, 0)];
    let before = proteins.clone();
    assert!(
        proteins
            .compute_modifications(&[identification(vec![candidate.clone()])], &BTreeSet::new())
            .is_err()
    );
    assert_eq!(proteins, before);
    candidate.evidences[0].start = None;
    assert!(
        proteins
            .compute_modifications(&[identification(vec![candidate])], &BTreeSet::new())
            .is_err()
    );
    assert_eq!(proteins, before);
}

#[test]
fn spectrum_processing_and_feature_clear_preserve_identification_ownership() {
    let id = identification(vec![hit("PEP", 1.0)]);
    let mut spectrum = MSSpectrum {
        peaks: vec![Peak1D::new(200.0, 4.0), Peak1D::new(100.0, 2.0)],
        peptide_identifications: vec![id.clone()],
        ..Default::default()
    };
    spectrum.sort_by_position().unwrap();
    Normalizer::default()
        .filter_spectrum(&mut spectrum)
        .unwrap();
    spectrum.select(&[0]).unwrap();
    assert_eq!(spectrum.peptide_identifications, std::slice::from_ref(&id));
    spectrum.clear(false);
    assert_eq!(spectrum.peptide_identifications, std::slice::from_ref(&id));
    spectrum.clear(true);
    assert!(spectrum.peptide_identifications.is_empty());
    let mut feature = Feature::default();
    feature.peptide_identifications.push(id.clone());
    let run = ProteinIdentification::default();
    let mut map = FeatureMap {
        features: vec![feature],
        protein_identifications: vec![run.clone()],
        unassigned_peptide_identifications: vec![id.clone()],
        ..Default::default()
    };
    map.validate().unwrap();
    map.clear(false);
    assert_eq!(map.protein_identifications, [run]);
    assert_eq!(map.unassigned_peptide_identifications, [id]);
    map.clear(true);
    assert_eq!(map, FeatureMap::default());
}

#[test]
fn map_and_feature_validation_checks_nested_identifications() {
    let mut invalid = identification(vec![hit("PEP", 1.0)]);
    invalid.mz = Some(f64::INFINITY);
    assert!(
        BaseFeature {
            peptide_identifications: vec![invalid.clone()],
            ..Default::default()
        }
        .validate()
        .is_err()
    );
    assert!(
        FeatureMap {
            unassigned_peptide_identifications: vec![invalid.clone()],
            ..Default::default()
        }
        .validate()
        .is_err()
    );
    assert!(
        ConsensusMap {
            unassigned_peptide_identifications: vec![invalid],
            ..Default::default()
        }
        .validate()
        .is_err()
    );
    let mut run = ProteinIdentification::default();
    run.search_parameters.precursor_tolerance = Tolerance::Ppm(-1.0);
    assert!(
        FeatureMap {
            protein_identifications: vec![run.clone()],
            ..Default::default()
        }
        .validate()
        .is_err()
    );
    assert!(
        ConsensusMap {
            protein_identifications: vec![run],
            ..Default::default()
        }
        .validate()
        .is_err()
    );
}

#[test]
fn peak_file_writers_reject_identification_loss_before_writing_bytes() {
    let spectrum = MSSpectrum {
        peptide_identifications: vec![identification(vec![hit("PEP", 1.0)])],
        ..Default::default()
    };
    let mut bytes = Vec::new();
    assert!(matches!(
        openms::format::dta::write(&mut bytes, &spectrum),
        Err(Error::Unsupported(_))
    ));
    assert!(bytes.is_empty());
    let experiment = MSExperiment {
        spectra: vec![MSSpectrum::default(), spectrum],
        ..Default::default()
    };
    assert!(matches!(
        openms::format::mgf::write(&mut bytes, &experiment),
        Err(Error::Unsupported(_))
    ));
    assert!(bytes.is_empty());
    #[cfg(feature = "mzml")]
    {
        assert!(matches!(
            openms::format::mzml::write(&mut bytes, &experiment),
            Err(Error::Unsupported(_))
        ));
        assert!(bytes.is_empty());
    }
}
