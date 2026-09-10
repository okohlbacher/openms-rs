// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source expectations: pinned PROCESSING/ID/IDFilter and IDFilter_test.cpp/idXML.

use openms::analysis::id_filter::*;
use openms::chemistry::AASequence;
use openms::comparison::Tolerance;
use openms::identification::{
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinGroup, ProteinHit,
    ProteinIdentification,
};
use openms::kernel::{ConsensusFeature, ConsensusMap, DataArray, Feature, FeatureMap};
use openms::metadata::MetaValue;
use std::collections::BTreeSet;

fn hit(sequence: &str, score: f64) -> PeptideHit {
    PeptideHit::new(score, 99, 2, AASequence::parse(sequence).unwrap()).unwrap()
}
fn evidence(accession: &str) -> PeptideEvidence {
    PeptideEvidence {
        protein_accession: accession.into(),
        ..Default::default()
    }
}
fn id(hits: Vec<PeptideHit>) -> PeptideIdentification {
    PeptideIdentification {
        identifier: "run1".into(),
        score_type: "Mascot".into(),
        hits,
        ..Default::default()
    }
}
fn proteins(run: &str, accessions: &[&str]) -> ProteinIdentification {
    ProteinIdentification {
        identifier: run.into(),
        hits: accessions
            .iter()
            .map(|a| ProteinHit::new(0.0, 0, *a, "").unwrap())
            .collect(),
        ..Default::default()
    }
}
fn names(values: &[&str]) -> BTreeSet<String> {
    values.iter().map(|s| (*s).into()).collect()
}
fn scores(id: &PeptideIdentification) -> Vec<f64> {
    id.hits.iter().map(|h| h.score).collect()
}
fn accessions(id: &ProteinIdentification) -> Vec<&str> {
    id.hits.iter().map(|h| h.accession.as_str()).collect()
}
fn annotated(sequence: &str, score: f64, accessions: &[&str]) -> PeptideHit {
    let mut hit = hit(sequence, score);
    hit.evidences = accessions.iter().map(|a| evidence(a)).collect();
    hit
}

// Literal small data from IDFilter_test.idXML; independent of the new idXML reader.
fn source_fixture() -> (Vec<PeptideIdentification>, Vec<ProteinIdentification>) {
    let input = [
        ("MRSLGYVAVISAVATDTDK", 33.85, "Q824A5"),
        ("EGASTDFAALRTFLAEDGK", 12.73, "S53854"),
        ("DLEPGTDYEVTVSTLFGR", 11.79, "S53854"),
        ("LHASGITVTEIPVTATNFK", 34.85, "Q872T5"),
        ("FINFGVNVEVLSRFQTK", 40.0, ""),
        ("MSLLSNMISIVKVGYNAR", 40.0, ""),
        ("THPYGHAIVAGIERYPSK", 39.0, ""),
        ("AITSDFANQAKTVLQNFK", 11.1, ""),
        ("TGCDTWGQGTLVTVSSASTK", 10.93, ""),
        ("TLCHHDATFDNLVWTPK", 10.37, ""),
        ("MSLLSNM(Oxidation)ISIVKVGYNAR", 10.0, ""),
    ];
    let mut peptide = id(input
        .into_iter()
        .map(|(s, score, accession)| {
            let mut hit = hit(s, score);
            if !accession.is_empty() {
                hit.evidences.push(evidence(accession));
            }
            hit
        })
        .collect());
    peptide.sort().unwrap();
    let mut protein = proteins("run1", &["Q824A5", "AAD30739", "S53854", "Q872T5"]);
    for (h, score) in protein.hits.iter_mut().zip([32.3, 28.1, 27.7, 27.4]) {
        h.score = score;
    }
    (vec![peptide], vec![protein])
}

#[test]
fn source_score_top_n_and_dense_rank_expectations() {
    let (original, mut protein) = source_fixture();
    let mut peptides = original.clone();
    filter_peptides_by_score(&mut peptides, 33.0).unwrap();
    assert_eq!(scores(&peptides[0]), [40.0, 40.0, 39.0, 34.85, 33.85]);
    assert_eq!(peptides[0].hits[0].sequence.as_str(), "FINFGVNVEVLSRFQTK");
    filter_peptides_by_score(&mut peptides, 41.0).unwrap();
    assert!(peptides[0].hits.is_empty());
    assert_eq!(peptides[0].score_type, "Mascot");

    let mut peptides = original.clone();
    keep_n_best_peptide_hits(&mut peptides, 3).unwrap();
    assert_eq!(scores(&peptides[0]), [40.0, 40.0, 39.0]);
    let mut peptides = original;
    filter_peptides_by_rank(&mut peptides, 1, Some(5)).unwrap();
    assert_eq!(
        scores(&peptides[0]),
        [40.0, 40.0, 39.0, 34.85, 33.85, 12.73]
    );
    assert!(peptides[0].hits.iter().all(|hit| hit.rank == 99));
    filter_proteins_by_rank(&mut protein, 3, Some(10)).unwrap();
    assert_eq!(accessions(&protein[0]), ["S53854", "Q872T5"]);
}

#[test]
fn best_hits_keep_ties_or_clear_them_and_honor_lower_scores() {
    let (mut peptides, _) = source_fixture();
    keep_best_peptide_hits(&mut peptides, false).unwrap();
    assert_eq!(scores(&peptides[0]), [40.0, 40.0]);
    keep_best_peptide_hits(&mut peptides, true).unwrap();
    assert!(peptides[0].hits.is_empty());

    let mut peptide = id(vec![
        hit("CCC", 2.0),
        hit("AAA", -0.0),
        hit("DDD", 1.0),
        hit("MMM", 0.0),
    ]);
    peptide.higher_score_better = false;
    let mut peptides = vec![peptide];
    filter_peptides_by_score(&mut peptides, 0.0).unwrap();
    keep_best_peptide_hits(&mut peptides, false).unwrap();
    assert_eq!(
        peptides[0]
            .hits
            .iter()
            .map(|h| h.sequence.as_str())
            .collect::<Vec<_>>(),
        ["AAA", "MMM"]
    );
    keep_n_best_peptide_hits(&mut peptides, 1).unwrap();
    keep_best_peptide_hits(&mut peptides, true).unwrap();
    assert_eq!(peptides[0].hits[0].sequence.as_str(), "AAA");
    keep_n_best_peptide_hits(&mut peptides, 0).unwrap();
    assert!(peptides[0].hits.is_empty());
}

#[test]
fn protein_filters_preserve_direction_metadata_and_empty_records_until_cleanup() {
    let mut protein = proteins("run", &["A", "B", "C"]);
    for (hit, score) in protein.hits.iter_mut().zip([2.0, 1.0, 1.0]) {
        hit.score = score;
    }
    protein.higher_score_better = false;
    protein.metadata.insert("sample".into(), "kept".into());
    let mut ids = vec![protein];
    filter_proteins_by_score(&mut ids, 1.0).unwrap();
    assert_eq!(accessions(&ids[0]), ["B", "C"]);
    keep_n_best_protein_hits(&mut ids, 1).unwrap();
    assert_eq!(accessions(&ids[0]), ["B"]);
    assert_eq!(ids[0].metadata["sample"].as_str().unwrap(), "kept");
    filter_proteins_by_accessions(&mut ids, &names(&["B"]), MatchAction::Remove).unwrap();
    assert_eq!(ids.len(), 1);
    remove_empty_protein_identifications(&mut ids).unwrap();
    assert!(ids.is_empty());
}

#[test]
fn length_and_signed_charge_bounds_are_inclusive_without_endpoint_overflow() {
    let mut hits = vec![
        hit("AAA", 1.0),
        hit("C(Carbamidomethyl)AAA", 1.0),
        hit("AAAAA", 1.0),
        hit("AAAA", 1.0),
    ];
    for (hit, charge) in hits.iter_mut().zip([i32::MIN, 0, 3, i32::MAX]) {
        hit.charge = charge;
    }
    let mut ids = vec![id(hits)];
    filter_peptides_by_length(&mut ids, 4, Some(4)).unwrap();
    assert_eq!(ids[0].hits.len(), 2);
    filter_peptides_by_charge(&mut ids, 0, Some(i32::MAX)).unwrap();
    assert_eq!(ids[0].hits.len(), 2);
    filter_peptides_by_charge(&mut ids, i32::MAX, None).unwrap();
    assert_eq!(ids[0].hits[0].charge, i32::MAX);
    filter_peptides_by_length(&mut ids, 0, Some(usize::MAX)).unwrap();
    assert_eq!(ids[0].hits.len(), 1);
    let mut negative = hit("AAA", 0.0);
    negative.charge = i32::MIN;
    let mut ids = vec![id(vec![negative])];
    filter_peptides_by_charge(&mut ids, i32::MIN, Some(i32::MIN)).unwrap();
    assert_eq!(ids[0].hits.len(), 1);
}

#[test]
fn accession_matching_is_any_reference_and_does_not_trim_evidence() {
    let (mut peptides, mut proteins) = source_fixture();
    let accepted = names(&["Q824A5", "Q872T5"]);
    filter_peptides_by_accessions(&mut peptides, &accepted, MatchAction::Keep).unwrap();
    assert_eq!(peptides[0].hits.len(), 2);
    assert_eq!(peptides[0].hits[0].sequence.as_str(), "LHASGITVTEIPVTATNFK");
    filter_proteins_by_accessions(&mut proteins, &accepted, MatchAction::Remove).unwrap();
    assert_eq!(accessions(&proteins[0]), ["AAD30739", "S53854"]);
    let shared = annotated("AAA", 1.0, &["A", "B", ""]);
    let mut ids = vec![id(vec![shared])];
    filter_peptides_by_accessions(&mut ids, &names(&["A"]), MatchAction::Keep).unwrap();
    assert_eq!(ids[0].hits[0].evidences.len(), 3);
    filter_peptides_by_accessions(&mut ids, &names(&[""]), MatchAction::Keep).unwrap();
    assert!(ids[0].hits.is_empty());
}

#[test]
fn modification_and_sequence_filters_include_termini_and_exact_full_ids() {
    let (original, _) = source_fixture();
    assert_eq!(
        extract_peptide_sequences(&original, false).unwrap().len(),
        11
    );
    assert_eq!(
        extract_peptide_sequences(&original, true).unwrap().len(),
        10
    );
    let mut peptides = original.clone();
    filter_peptides_by_modifications(&mut peptides, &BTreeSet::new(), MatchAction::Keep).unwrap();
    assert_eq!(peptides[0].hits.len(), 1);
    let unmodified = extract_peptide_sequences(&peptides, true).unwrap();
    let mut selected = original;
    filter_peptides_by_sequences(&mut selected, &unmodified, true, MatchAction::Keep).unwrap();
    assert_eq!(scores(&selected[0]), [40.0, 10.0]);
    filter_peptides_by_modifications(
        &mut selected,
        &names(&["Oxidation (M)"]),
        MatchAction::Remove,
    )
    .unwrap();
    assert_eq!(scores(&selected[0]), [40.0]);
    let mut terminal = vec![id(vec![
        hit("(Acetyl)PEPTIDER.(Arg-loss)", 1.0),
        hit("M(Oxidation)PEP", 2.0),
    ])];
    filter_peptides_by_modifications(
        &mut terminal,
        &names(&["Acetyl (N-term)", "Oxidation (M)"]),
        MatchAction::Keep,
    )
    .unwrap();
    assert_eq!(terminal[0].hits.len(), 2);
    filter_peptides_by_modifications(
        &mut terminal,
        &names(&["Arg-loss (C-term R)"]),
        MatchAction::Keep,
    )
    .unwrap();
    assert_eq!(terminal[0].hits.len(), 1);
    filter_peptides_by_modifications(&mut terminal, &names(&["Acetyl"]), MatchAction::Keep)
        .unwrap();
    assert!(terminal[0].hits.is_empty()); // source matches full IDs, not aliases
    let modified = extract_peptide_sequences(&peptides, false).unwrap();
    filter_peptides_by_sequences(&mut peptides, &modified, false, MatchAction::Remove).unwrap();
    assert!(peptides[0].hits.is_empty());
}

#[test]
fn source_duplicate_semantics_keep_first_full_hit_or_first_modified_sequence() {
    let mut a = hit("DFPIANGER", 0.3);
    a.charge = 1;
    let mut b = a.clone();
    b.charge = 2;
    let mut c = b.clone();
    c.score = 0.5;
    let mut d = c.clone();
    d.sequence = AASequence::parse("DFPIANGEK").unwrap();
    let mut e = d.clone();
    e.charge = 5;
    let mut ids = vec![id(vec![a, b, c, d.clone(), d.clone(), d, e])];
    remove_duplicate_peptide_hits(&mut ids, DuplicatePolicy::Exact).unwrap();
    assert_eq!(ids[0].hits.len(), 5);
    assert_eq!(ids[0].hits[3].charge, 2);
    assert_eq!(ids[0].hits[4].charge, 5);
    remove_duplicate_peptide_hits(&mut ids, DuplicatePolicy::Sequence).unwrap();
    assert_eq!(ids[0].hits.len(), 2);
    assert_eq!(ids[0].hits[0].score, 0.3); // not the best score

    let mut different_meta = ids[0].hits[0].clone();
    different_meta
        .metadata
        .insert("origin".into(), "second".into());
    ids[0].hits.push(different_meta);
    remove_duplicate_peptide_hits(&mut ids, DuplicatePolicy::Exact).unwrap();
    assert_eq!(ids[0].hits.len(), 3);
    ids[0].hits.push(hit("M(Oxidation)PEP", 1.0));
    ids[0].hits.push(hit("MPEP", 1.0));
    remove_duplicate_peptide_hits(&mut ids, DuplicatePolicy::Sequence).unwrap();
    assert_eq!(ids[0].hits.len(), 4);
}

#[test]
fn exact_duplicates_preserve_distinct_native_analysis_results_and_float_metadata() {
    use openms::identification::AnalysisResult;
    let mut first = hit("AAA", 1.0);
    first
        .metadata
        .insert("numeric".into(), MetaValue::try_from(1.0).unwrap());
    let mut second = first.clone();
    second
        .metadata
        .insert("numeric".into(), MetaValue::try_from(1.0 + 1e-8).unwrap());
    let mut third = first.clone();
    third.analysis_results.push(AnalysisResult {
        score_type: "secondary analysis".into(),
        higher_is_better: false,
        main_score: 0.02,
        sub_scores: Default::default(),
    });
    let mut ids = vec![id(vec![first.clone(), second, third, first])];
    remove_duplicate_peptide_hits(&mut ids, DuplicatePolicy::Exact).unwrap();
    assert_eq!(ids[0].hits.len(), 3);
    assert_eq!(ids[0].hits[2].analysis_results[0].main_score, 0.02);
}

#[test]
fn decoy_and_unique_filters_follow_explicit_source_annotations() {
    let mut hits = vec![hit("AAA", 1.0); 7];
    hits[0]
        .metadata
        .insert("target_decoy".into(), "target".into());
    hits[1]
        .metadata
        .insert("target_decoy".into(), "decoy".into());
    hits[2]
        .metadata
        .insert("target_decoy".into(), "target+decoy".into());
    hits[4].metadata.insert("isDecoy".into(), "true".into());
    hits[5].metadata.insert("isDecoy".into(), "false".into());
    hits[6]
        .metadata
        .insert("target_decoy".into(), "target".into());
    hits[6].metadata.insert("isDecoy".into(), "true".into());
    let mut ids = vec![id(hits)];
    remove_decoy_peptide_hits(&mut ids).unwrap();
    assert_eq!(ids[0].hits.len(), 4);
    assert_eq!(
        ids[0].hits[1].metadata["target_decoy"].as_str().unwrap(),
        "target+decoy"
    );
    for (hit, annotation) in ids[0]
        .hits
        .iter_mut()
        .zip(["non-unique", "unmatched", "", "unique"])
    {
        if !annotation.is_empty() {
            hit.metadata
                .insert("protein_references".into(), annotation.into());
        }
    }
    keep_unique_peptides_per_protein(&mut ids).unwrap();
    assert_eq!(ids[0].hits.len(), 1);
    let mut protein = proteins("run1", &["A", "B"]);
    protein.hits[0]
        .metadata
        .insert("isDecoy".into(), "true".into());
    let mut protein_ids = vec![protein];
    remove_decoy_protein_hits(&mut protein_ids).unwrap();
    assert_eq!(accessions(&protein_ids[0]), ["B"]);
}

#[test]
fn coordinates_and_mass_error_use_observed_mz_with_inclusive_tolerance() {
    let mut ids: Vec<_> = [None, Some(1.0), Some(2.0), Some(2.5), Some(1.5)]
        .into_iter()
        .map(|rt| PeptideIdentification {
            rt,
            ..id(Vec::new())
        })
        .collect();
    filter_peptides_by_rt(&mut ids, 1.0, 1.9).unwrap();
    assert_eq!(
        ids.iter().map(|id| id.rt).collect::<Vec<_>>(),
        [Some(1.0), Some(1.5)]
    );
    let mut ids: Vec<_> = [None, Some(111.1), Some(222.2), Some(225.5), Some(115.5)]
        .into_iter()
        .map(|mz| PeptideIdentification {
            mz,
            ..id(Vec::new())
        })
        .collect();
    filter_peptides_by_mz(&mut ids, 112.0, 223.3).unwrap();
    assert_eq!(
        ids.iter().map(|id| id.mz).collect::<Vec<_>>(),
        [Some(222.2), Some(115.5)]
    );

    let mut a = hit("PEPTIDE", 1.0);
    a.charge = 0; // source interprets unknown charge as +1
    let mut b = a.clone();
    b.charge = 2;
    let mut peptide = id(vec![a.clone(), b]);
    let mz = a.sequence.mz(1).unwrap();
    peptide.mz = Some(mz + 0.125);
    let mut ids = vec![peptide];
    filter_peptides_by_mz_error(&mut ids, Tolerance::Absolute(0.125)).unwrap();
    assert_eq!(ids[0].hits.len(), 1);
    let original = ids.clone();
    let observed = ids[0].mz.unwrap();
    filter_peptides_by_mz_error(&mut ids, Tolerance::Ppm(0.126 / observed * 1e6)).unwrap();
    assert_eq!(ids[0].hits.len(), 1);
    let mut ids = original;
    filter_peptides_by_mz_error(&mut ids, Tolerance::Ppm(0.124 / observed * 1e6)).unwrap();
    assert!(ids[0].hits.is_empty());
}

#[test]
fn cleanup_uses_run_ids_and_union_of_repeated_runs_without_losing_shared_evidence() {
    let mut peptides = vec![
        id(vec![
            annotated("AAA", 1.0, &["A", "B", "A", "GONE"]),
            annotated("CCC", 2.0, &["GONE"]),
        ]),
        PeptideIdentification {
            identifier: "run2".into(),
            ..id(vec![annotated("DDD", 1.0, &["B"])])
        },
        PeptideIdentification {
            identifier: "missing".into(),
            ..id(vec![annotated("MMM", 1.0, &["A"])])
        },
    ];
    let protein_ids = vec![
        proteins("run1", &["A"]),
        proteins("run1", &["B"]),
        proteins("run2", &["A"]),
    ];
    remove_dangling_protein_references(&mut peptides, &protein_ids, false).unwrap();
    assert_eq!(
        peptides[0].hits[0]
            .evidences
            .iter()
            .map(|e| e.protein_accession.as_str())
            .collect::<Vec<_>>(),
        ["A", "B", "A"]
    );
    assert!(peptides[0].hits[1].evidences.is_empty());
    assert!(peptides[1].hits[0].evidences.is_empty());
    assert!(peptides[2].hits[0].evidences.is_empty());
    remove_dangling_protein_references(&mut peptides, &protein_ids, true).unwrap();
    assert_eq!(peptides[0].hits.len(), 1);
    assert_eq!(peptides.len(), 3);
    remove_empty_peptide_identifications(&mut peptides).unwrap();
    assert_eq!(peptides.len(), 1);
    let mut proteins = protein_ids;
    remove_unreferenced_proteins(&mut proteins, &peptides).unwrap();
    assert_eq!(accessions(&proteins[0]), ["A"]);
    assert_eq!(accessions(&proteins[1]), ["B"]);
    assert!(proteins[2].hits.is_empty());
}

#[test]
fn group_cleanup_retains_sample_arrays_and_reports_only_partially_removed_groups() {
    let quantities = DataArray::new("abundances", vec![11.0_f32, 22.0, 33.0]);
    let channels = DataArray::new("channels", vec![1_i32, 2]);
    let filenames = DataArray::new("filenames", vec!["fA".to_string(), "fB".to_string()]);
    let mut groups = vec![
        ProteinGroup {
            accessions: vec!["A".into()],
            probability: 0.1,
            ..Default::default()
        },
        ProteinGroup {
            probability: 0.2,
            accessions: vec!["B".into(), "C".into()],
            float_data_arrays: vec![quantities.clone()],
            integer_data_arrays: vec![channels.clone()],
            string_data_arrays: vec![filenames.clone()],
        },
    ];
    let mut hits = proteins("run1", &["C", "B", "A"]).hits;
    let original = groups.clone();
    assert!(update_protein_groups(&mut groups, &hits).unwrap());
    assert_eq!(groups, original);
    hits.pop();
    assert!(update_protein_groups(&mut groups, &hits).unwrap());
    assert_eq!(groups.len(), 1);
    hits.pop();
    assert!(!update_protein_groups(&mut groups, &hits).unwrap());
    assert_eq!(groups[0].accessions, ["C"]);
    assert_eq!(groups[0].probability, 0.2);
    assert_eq!(groups[0].float_data_arrays, [quantities]);
    assert_eq!(groups[0].integer_data_arrays, [channels]);
    assert_eq!(groups[0].string_data_arrays, [filenames]);
    let mut hits = proteins("run1", &["A", "C", "B"]).hits;
    remove_ungrouped_proteins(&groups, &mut hits).unwrap();
    assert_eq!(hits[0].accession, "C");
    assert_eq!(hits.len(), 1);
}

#[test]
fn consensus_map_cleanup_covers_assigned_and_unassigned_records() {
    let mut feature = ConsensusFeature::default();
    feature.peptide_identifications.push(id(vec![
        annotated("AAA", 2.0, &["A", "DECOY_A"]),
        annotated("CCC", 1.0, &["DECOY_A"]),
    ]));
    let mut map = ConsensusMap {
        features: vec![feature],
        protein_identifications: vec![proteins("run1", &["A", "B"])],
        unassigned_peptide_identifications: vec![id(vec![annotated(
            "DDD",
            1.0,
            &["B", "DECOY_B"],
        )])],
        ..Default::default()
    };
    let original = map.clone();
    remove_dangling_consensus_references(&mut map, false).unwrap();
    assert_eq!(map.features[0].peptide_identifications[0].hits.len(), 2);
    assert!(
        map.features[0].peptide_identifications[0].hits[1]
            .evidences
            .is_empty()
    );
    assert_eq!(
        map.unassigned_peptide_identifications[0].hits[0].evidences[0].protein_accession,
        "B"
    );
    remove_dangling_consensus_references(&mut map, true).unwrap();
    assert_eq!(map.features[0].peptide_identifications[0].hits.len(), 1);
    remove_unreferenced_consensus_proteins(&mut map, true).unwrap();
    assert_eq!(accessions(&map.protein_identifications[0]), ["A", "B"]);
    remove_unreferenced_consensus_proteins(&mut map, false).unwrap();
    assert_eq!(accessions(&map.protein_identifications[0]), ["A"]);
    remove_dangling_consensus_references(&mut map, true).unwrap();
    remove_empty_consensus_identifications(&mut map).unwrap();
    assert!(map.unassigned_peptide_identifications.is_empty());
    assert_eq!(map.features.len(), 1);
    let mut map = original;
    keep_n_best_hits_in_consensus_map(&mut map, 0).unwrap();
    assert!(map.features[0].peptide_identifications[0].hits.is_empty());
    assert!(map.unassigned_peptide_identifications[0].hits.is_empty());
}

#[test]
fn feature_map_adapters_leave_subordinates_and_geometry_attached() {
    let mut subordinate = Feature::default();
    subordinate
        .peptide_identifications
        .push(id(vec![annotated("MMM", 3.0, &["SUB"])]));
    let mut feature = Feature {
        subordinates: vec![subordinate.clone()],
        ..Default::default()
    };
    feature.rt = 4.0;
    feature.peptide_identifications.push(id(vec![
        annotated("AAA", 1.0, &["A"]),
        annotated("DDD", 2.0, &["B"]),
    ]));
    let mut map = FeatureMap {
        features: vec![feature],
        protein_identifications: vec![proteins("run1", &["A", "B", "C", "SUB"])],
        unassigned_peptide_identifications: vec![id(vec![annotated("CCC", 1.0, &["C"])])],
        ..Default::default()
    };
    keep_n_best_hits_in_feature_map(&mut map, 1).unwrap();
    assert_eq!(
        map.features[0].peptide_identifications[0].hits[0]
            .sequence
            .as_str(),
        "DDD"
    );
    remove_unreferenced_feature_proteins(&mut map, true).unwrap();
    assert_eq!(accessions(&map.protein_identifications[0]), ["B", "C"]);
    remove_unreferenced_feature_proteins(&mut map, false).unwrap();
    assert_eq!(accessions(&map.protein_identifications[0]), ["B"]);
    remove_dangling_feature_references(&mut map, true).unwrap();
    remove_empty_feature_identifications(&mut map).unwrap();
    assert!(map.unassigned_peptide_identifications.is_empty());
    assert_eq!(map.features[0].subordinates, [subordinate]);
    assert_eq!(map.features[0].rt, 4.0);
}

#[test]
fn spectrum_ranking_keeps_full_ids_stabilizes_ties_and_checks_compatibility_atomically() {
    let mut ids = vec![
        id(Vec::new()),
        id(vec![hit("AAA", 3.0), hit("CCC", 1.0)]),
        id(vec![hit("DDD", 2.0), hit("MMM", 3.0)]),
        id(Vec::new()),
    ];
    ids[1].set_spectrum_reference("first tie");
    ids[2].set_spectrum_reference("second tie");
    let mut mixed = ids.clone();
    mixed[3].higher_score_better = false;
    let original = mixed.clone();
    assert!(keep_n_best_spectra(&mut mixed, 1).is_err());
    assert_eq!(mixed, original);
    mixed[3].higher_score_better = true;
    mixed[3].score_type.clear();
    let original = mixed.clone();
    assert!(keep_n_best_spectra(&mut mixed, 1).is_err());
    assert_eq!(mixed, original);
    keep_n_best_spectra(&mut ids, 2).unwrap();
    assert_eq!(ids[0].spectrum_reference(), "first tie");
    assert_eq!(ids[1].spectrum_reference(), "second tie");
    assert!(ids.iter().all(|id| id.hits.len() == 2));
    assert_eq!(scores(&ids[1]), [3.0, 2.0]);
    keep_n_best_spectra(&mut ids, 0).unwrap();
    assert!(ids.is_empty());
}

#[test]
fn invalid_parameters_records_and_resource_excess_leave_every_record_unchanged() {
    let mut ids = vec![
        id(vec![hit("AAA", 1.0), hit("CCC", 2.0)]),
        id(vec![hit("DDD", 3.0)]),
    ];
    let original = ids.clone();
    assert!(filter_peptides_by_length(&mut ids, 9, Some(8)).is_err());
    assert!(filter_peptides_by_charge(&mut ids, 2, Some(-1)).is_err());
    assert!(filter_peptides_by_rank(&mut ids, 0, None).is_err());
    assert!(filter_peptides_by_rank(&mut ids, 3, Some(2)).is_err());
    assert!(filter_peptides_by_score(&mut ids, f64::NAN).is_err());
    assert!(filter_peptides_by_rt(&mut ids, 0.0, f64::INFINITY).is_err());
    assert_eq!(ids, original);
    ids[0].mz = Some(100.0); // would remove every first-ID hit
    let original = ids.clone();
    assert!(filter_peptides_by_mz_error(&mut ids, Tolerance::Absolute(0.01)).is_err());
    assert_eq!(ids, original); // second ID lacks m/z
    ids[1].hits[0].score = f64::INFINITY;
    let original = ids.clone();
    assert!(keep_n_best_peptide_hits(&mut ids, 1).is_err());
    assert!(filter_peptides_by_score(&mut ids, 2.0).is_err());
    assert!(remove_duplicate_peptide_hits(&mut ids, DuplicatePolicy::Sequence).is_err());
    assert_eq!(ids, original);
    let mut large = vec![id(vec![hit("AAA", 1.0); 4500])];
    assert!(remove_duplicate_peptide_hits(&mut large, DuplicatePolicy::Exact).is_err());
    assert_eq!(large[0].hits.len(), 4500);
    remove_duplicate_peptide_hits(&mut large, DuplicatePolicy::Sequence).unwrap();
    assert_eq!(large[0].hits.len(), 1);
    let mut a = Feature::default();
    a.peptide_identifications
        .push(id(vec![hit("AAA", 2.0), hit("CCC", 1.0)]));
    let mut b = Feature::default();
    b.peptide_identifications.push(id(vec![hit("DDD", 0.0)]));
    b.peptide_identifications[0].hits[0].score = f64::INFINITY;
    let mut map = FeatureMap::from_features(vec![a, b]);
    let original = map.clone();
    assert!(keep_n_best_hits_in_feature_map(&mut map, 1).is_err());
    assert_eq!(map, original);
}
