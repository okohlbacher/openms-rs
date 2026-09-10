// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::analysis::id_ripper::*;
use openms::chemistry::AASequence;
use openms::identification::{
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinGroup, ProteinHit,
    ProteinIdentification,
};
use openms::metadata::MetaValue;

fn run(name: &str, paths: &[&str]) -> ProteinIdentification {
    ProteinIdentification {
        identifier: name.into(),
        score_type: "protein score".into(),
        search_engine: "engine".into(),
        primary_ms_run_paths: paths.iter().map(|s| (*s).into()).collect(),
        raw_ms_run_paths: vec!["raw/path.raw".into()],
        hits: ["A", "B", "C"]
            .into_iter()
            .enumerate()
            .map(|(i, acc)| ProteinHit::new(i as f64 + 1., 10, acc, "AAAK").unwrap())
            .collect(),
        protein_groups: vec![ProteinGroup {
            probability: 0.9,
            accessions: vec!["A".into(), "B".into()],
            ..Default::default()
        }],
        ..Default::default()
    }
}
fn peptide(run: &str, key: &str, origin: MetaValue, accessions: &[&str]) -> PeptideIdentification {
    let mut hit = PeptideHit::new(0.1, 7, 2, AASequence::parse("AAAK").unwrap()).unwrap();
    hit.evidences = accessions
        .iter()
        .map(|s| PeptideEvidence {
            protein_accession: (*s).into(),
            ..Default::default()
        })
        .collect();
    let mut id = PeptideIdentification {
        identifier: run.into(),
        score_type: "PEP".into(),
        higher_score_better: false,
        hits: vec![hit],
        ..Default::default()
    };
    id.metadata.insert(key.into(), origin);
    id.metadata.insert("sample".into(), "kept".into());
    id
}

#[test]
fn origin_detection_checks_every_record_and_rejects_ambiguous_formats() {
    for (key, expected) in [
        ("file_origin", OriginAnnotationFormat::FileOrigin),
        ("map_index", OriginAnnotationFormat::MapIndex),
        ("id_merge_index", OriginAnnotationFormat::IdMergeIndex),
    ] {
        let a = peptide("run", key, 0_i64.into(), &["A"]);
        assert_eq!(
            detect_origin_annotation_format(std::slice::from_ref(&a)).unwrap(),
            expected
        );
        let mut empty = a.clone();
        empty.hits.clear();
        empty.metadata.remove(key);
        assert!(detect_origin_annotation_format(&[a.clone(), empty]).is_err());
        let mut duplicate = a;
        duplicate.metadata.insert(
            if key == "file_origin" {
                "map_index"
            } else {
                "file_origin"
            }
            .into(),
            "other".into(),
        );
        assert!(detect_origin_annotation_format(&[duplicate]).is_err());
    }
    assert!(detect_origin_annotation_format(&[]).is_err());
    assert!(
        detect_origin_annotation_format(&[peptide("run", "unknown", "x".into(), &["A"])]).is_err()
    );
    assert!(
        detect_origin_annotation_format(&[
            peptide("run", "map_index", 0_i64.into(), &["A"]),
            peptide("run", "id_merge_index", 0_i64.into(), &["A"])
        ])
        .is_err()
    );
}

#[test]
fn file_origin_uses_first_occurrence_indices_and_preserves_input_metadata() {
    let mut p = run("r", &["full/source1", "full/source2"]);
    p.metadata.insert("file_origin".into(), "merged".into());
    p.metadata.insert("operator".into(), "retained".into());
    let proteins = vec![p];
    let peptides = vec![
        peptide("r", "file_origin", "/data/z.raw.mzML".into(), &["B"]),
        peptide("r", "file_origin", "/data/a.mzML".into(), &["A"]),
        peptide(
            "r",
            "file_origin",
            "/data/z.raw.mzML".into(),
            &["C", "A", "B"],
        ),
    ];
    let before = (proteins.clone(), peptides.clone());
    let result = IDRipper::default().rip(&proteins, &peptides).unwrap();
    assert_eq!((proteins, peptides), before);
    assert_eq!(result.files.len(), 2);
    assert_eq!(result.files[0].identifier.file_origin_index, 0);
    assert_eq!(result.files[0].identifier.output_basename, "z.raw");
    assert_eq!(result.files[1].identifier.output_basename, "a");
    assert_eq!(result.files[0].identifier.identification_run_index, None);
    let first = &result.files[0];
    assert_eq!(first.peptide_identifications.len(), 2);
    assert_eq!(
        first.protein_identifications[0]
            .hits
            .iter()
            .map(|h| h.accession.as_str())
            .collect::<Vec<_>>(),
        ["B", "A", "C"]
    );
    assert_eq!(first.protein_identifications[0].hits[0].rank, 10);
    assert_eq!(
        first.protein_identifications[0].metadata["operator"]
            .as_str()
            .unwrap(),
        "retained"
    );
    assert!(
        !first.protein_identifications[0]
            .metadata
            .contains_key("file_origin")
    );
    assert!(
        !first.peptide_identifications[0]
            .metadata
            .contains_key("file_origin")
    );
    assert_eq!(
        first.peptide_identifications[0].metadata["sample"]
            .as_str()
            .unwrap(),
        "kept"
    );
    assert_eq!(
        first.protein_identifications[0].primary_ms_run_paths,
        ["full/source1", "full/source2"]
    );
    assert!(first.protein_identifications[0].protein_groups.is_empty());
    assert_eq!(result.groups_not_copied, 2); // source copyMetaDataOnly omits groups, one copy per output
}

#[test]
fn source_index_modes_reduce_primary_paths_and_sort_numeric_keys() {
    for key in ["map_index", "id_merge_index"] {
        let proteins = [run("r", &["a.mzML", "b.mzML"])];
        let peptides = [
            peptide("r", key, 1_i64.into(), &["B"]),
            peptide("r", key, "0".into(), &["A"]),
            peptide("r", key, MetaValue::try_from(1.0).unwrap(), &["C"]),
        ];
        let result = IDRipper {
            split_identification_runs: true,
            numeric_filenames: false,
        }
        .rip(&proteins, &peptides)
        .unwrap();
        assert_eq!(result.files[0].identifier.file_origin_index, 0);
        assert_eq!(result.files[1].identifier.file_origin_index, 1);
        assert_eq!(result.files[0].identifier.identification_run_index, Some(0));
        assert_eq!(
            result.files[0].protein_identifications[0].primary_ms_run_paths,
            ["a.mzML"]
        );
        assert_eq!(
            result.files[1].protein_identifications[0].primary_ms_run_paths,
            ["b.mzML"]
        );
        assert_eq!(
            result.files[1].protein_identifications[0].raw_ms_run_paths,
            ["raw/path.raw"]
        );
        assert_eq!(result.files[1].peptide_identifications.len(), 2);
    }
}

#[test]
fn same_origin_combines_runs_but_keeps_run_specific_protein_scores() {
    let mut a = run("a", &["same.mzML"]);
    let mut b = run("b", &["same.mzML"]);
    a.hits[0].score = 1.;
    b.hits[0].score = 99.;
    let peptides = [
        peptide("b", "map_index", 0_i64.into(), &["A"]),
        peptide("a", "map_index", 0_i64.into(), &["A"]),
    ];
    let result = IDRipper::default()
        .rip(&[a.clone(), b.clone()], &peptides)
        .unwrap();
    assert_eq!(result.files.len(), 1);
    let runs = &result.files[0].protein_identifications;
    assert_eq!(
        runs.iter()
            .map(|r| r.identifier.as_str())
            .collect::<Vec<_>>(),
        ["b", "a"]
    );
    assert_eq!(runs[0].hits[0].score, 99.);
    assert_eq!(runs[1].hits[0].score, 1.); // source's global accession table would incorrectly use b's hit
    let result = IDRipper {
        split_identification_runs: true,
        numeric_filenames: true,
    }
    .rip(&[a, b], &peptides)
    .unwrap();
    assert_eq!(result.files.len(), 2);
    assert_eq!(result.files[0].protein_identifications[0].identifier, "a");
    assert_eq!(result.files[1].protein_identifications[0].identifier, "b");
}

#[test]
fn basename_collisions_and_same_index_different_paths_fail_explicitly() {
    let proteins = [
        run("a", &["/first/file.mzML"]),
        run("b", &["/second/file.mzML"]),
    ];
    let peptides = [
        peptide("a", "map_index", 0_i64.into(), &["A"]),
        peptide("b", "map_index", 0_i64.into(), &["A"]),
    ];
    for numeric_filenames in [false, true] {
        assert!(
            IDRipper {
                numeric_filenames,
                split_identification_runs: false
            }
            .rip(&proteins, &peptides)
            .is_err()
        );
    }
    assert!(
        IDRipper {
            numeric_filenames: false,
            split_identification_runs: true
        }
        .rip(&proteins, &peptides)
        .is_err()
    );
    let result = IDRipper {
        numeric_filenames: true,
        split_identification_runs: true,
    }
    .rip(&proteins, &peptides)
    .unwrap();
    assert_eq!(result.files.len(), 2);
    assert_ne!(
        result.files[0].identifier.origin_fullname,
        result.files[1].identifier.origin_fullname
    );
}

#[test]
fn empty_and_unreferenced_peptides_remain_available_as_skipped_records() {
    let proteins = [run("r", &["a.mzML"])];
    let mut empty = peptide("r", "map_index", 0_i64.into(), &["A"]);
    empty.hits.clear();
    let no_reference = peptide("r", "map_index", 0_i64.into(), &[]);
    let peptides = [empty, no_reference];
    let result = IDRipper::default().rip(&proteins, &peptides).unwrap();
    assert!(result.files.is_empty());
    assert_eq!(result.skipped_peptide_identifications, peptides);
    assert_eq!(result.groups_not_copied, 0);
}

#[test]
fn unknown_references_duplicate_identities_and_bad_indices_never_return_partial_output() {
    let proteins = [run("r", &["a.mzML"])];
    let valid = peptide("r", "map_index", 0_i64.into(), &["A"]);
    for bad_origin in [
        (-1_i64).into(),
        i64::from(i32::MAX).into(),
        i64::MAX.into(),
        MetaValue::try_from(0.5).unwrap(),
        "oops".into(),
    ] {
        assert!(
            IDRipper::default()
                .rip(
                    &proteins,
                    &[valid.clone(), peptide("r", "map_index", bad_origin, &["A"])]
                )
                .is_err()
        );
    }
    assert!(
        IDRipper::default()
            .rip(
                &proteins,
                &[
                    valid.clone(),
                    peptide("unknown", "map_index", 0_i64.into(), &["A"])
                ]
            )
            .is_err()
    );
    assert!(
        IDRipper::default()
            .rip(
                &proteins,
                &[
                    valid.clone(),
                    peptide("r", "map_index", 0_i64.into(), &["missing", "A"])
                ]
            )
            .is_err()
    );
    assert!(
        IDRipper::default()
            .rip(
                &[proteins[0].clone(), proteins[0].clone()],
                std::slice::from_ref(&valid)
            )
            .is_err()
    );
    let mut duplicate = proteins[0].clone();
    duplicate.hits.push(duplicate.hits[0].clone());
    assert!(IDRipper::default().rip(&[duplicate], &[valid]).is_err());
    assert!(
        IDRipper::default()
            .rip(&proteins, &[peptide("r", "file_origin", "".into(), &["A"])])
            .is_err()
    );
    assert!(
        IDRipper::default()
            .rip(
                &proteins,
                &[peptide("r", "file_origin", 1_i64.into(), &["A"])]
            )
            .is_err()
    );
}
