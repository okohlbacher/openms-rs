// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::Error;
use openms::chemistry::{AASequence, AdductInfo, NASequence};
use openms::identification::graph::*;
use std::collections::{BTreeMap, BTreeSet};

fn parent(g: &mut IdentificationData, text: &str, rna: bool) -> ParentId {
    let mut p = ParentSequence::new(text);
    p.molecule_type = if rna {
        MoleculeType::RNA
    } else {
        MoleculeType::Protein
    };
    p.sequence = "AAAA".into();
    g.register_parent_sequence(p).unwrap()
}
fn peptide(g: &mut IdentificationData, text: &str, parents: &[ParentId]) -> PeptideId {
    let mut p = IdentifiedPeptide::new(AASequence::parse(text).unwrap());
    for &id in parents {
        p.parent_matches.insert(
            id,
            BTreeSet::from([ParentMatch::new(Some(0), Some(text.len() - 1))]),
        );
    }
    g.register_identified_peptide(p).unwrap()
}
fn oligo(g: &mut IdentificationData, text: &str, parents: &[ParentId]) -> OligoId {
    let mut p = IdentifiedOligo::new(NASequence::parse(text).unwrap());
    for &id in parents {
        p.parent_matches.insert(id, BTreeSet::new());
    }
    g.register_identified_oligo(p).unwrap()
}
fn observation(g: &mut IdentificationData, text: &str) -> ObservationId {
    let file = g.register_input_file(InputFile::new("file")).unwrap();
    g.register_observation(Observation::new(text, file))
        .unwrap()
}
fn all_off() -> CleanupOptions {
    CleanupOptions {
        require_observation_match: false,
        require_identified_sequence: false,
        require_parent_match: false,
        require_parent_group: false,
        require_match_group: false,
    }
}
#[test]
fn source_cleanup_literal_counts_four_three_one_and_two_one_one() {
    let mut g = IdentificationData::new().unwrap();
    let protein = parent(&mut g, "protein_1", false);
    let rna = parent(&mut g, "rna_1", true);
    peptide(&mut g, "TEST", &[]);
    let p = peptide(&mut g, "PEPTIDE", &[protein]);
    peptide(&mut g, "EDIT", &[protein]);
    peptide(&mut g, "TESTPEP", &[protein]);
    oligo(&mut g, "ACGU", &[]);
    let r = oligo(&mut g, "UGCA", &[rna]); // empty positions still constitute a parent link
    let c = g
        .register_identified_compound(IdentifiedCompound::new("compound_1"))
        .unwrap();
    let obs = observation(&mut g, "spectrum_1");
    for molecule in [IdentifiedMolecule::from(p), r.into(), c.into()] {
        g.register_observation_match(ObservationMatch::new(molecule, obs))
            .unwrap();
    }
    assert_eq!((g.peptide_count(), g.oligo_count()), (4, 2));
    g.cleanup(CleanupOptions {
        require_observation_match: false,
        ..Default::default()
    })
    .unwrap();
    assert_eq!((g.peptide_count(), g.oligo_count()), (3, 1));
    g.cleanup(CleanupOptions::default()).unwrap();
    assert_eq!((g.peptide_count(), g.oligo_count()), (1, 1));
    assert_eq!(
        (
            g.parent_count(),
            g.compound_count(),
            g.observation_match_count()
        ),
        (2, 1, 3)
    );
}
#[test]
fn all_thirty_two_option_combinations_follow_the_ordered_source_cascade() {
    for bits in 0..32 {
        let o = CleanupOptions {
            require_observation_match: bits & 1 != 0,
            require_identified_sequence: bits & 2 != 0,
            require_parent_match: bits & 4 != 0,
            require_parent_group: bits & 8 != 0,
            require_match_group: bits & 16 != 0,
        };
        let mut g = IdentificationData::new().unwrap();
        let a = parent(&mut g, "a", false);
        let b = parent(&mut g, "b", false);
        let unused = parent(&mut g, "unused", false);
        let rna = parent(&mut g, "rna", true);
        let p = peptide(&mut g, "A", &[a]);
        let q = peptide(&mut g, "C", &[b]);
        let u = peptide(&mut g, "D", &[]);
        let r = oligo(&mut g, "A", &[rna]);
        let c = g
            .register_identified_compound(IdentifiedCompound::new("c"))
            .unwrap();
        let obs = observation(&mut g, "used");
        let empty = observation(&mut g, "unused");
        let adduct = g
            .register_adduct(AdductInfo::parse("M+H;1+").unwrap())
            .unwrap();
        let pm = g
            .register_observation_match(ObservationMatch::new(p, obs))
            .unwrap();
        let um = g
            .register_observation_match(ObservationMatch::new(u, obs))
            .unwrap();
        let mut rm = ObservationMatch::new(r, obs);
        rm.adduct = Some(adduct);
        let rm = g.register_observation_match(rm).unwrap();
        let cm = g
            .register_observation_match(ObservationMatch::new(c, obs))
            .unwrap();
        let pg = g
            .register_parent_group_set(ParentGroupSet {
                groups: vec![ParentGroup::new(BTreeSet::from([a]))],
                ..Default::default()
            })
            .unwrap();
        g.register_observation_match_group(ObservationMatchGroup::new(BTreeSet::from([pm])))
            .unwrap();
        g.cleanup(o).unwrap();
        let q_alive =
            (!o.require_parent_group || !o.require_parent_match) && !o.require_observation_match;
        let u_match = !o.require_parent_match && !o.require_match_group;
        let r_before = !o.require_parent_group || !o.require_parent_match;
        let r_match = r_before && !o.require_match_group;
        let r_alive = r_before && (!o.require_observation_match || r_match);
        assert!(
            g.parent(a).is_ok() && g.peptide(p).is_ok() && g.observation_match(pm).is_ok(),
            "bits {bits}"
        );
        assert_eq!(g.peptide(q).is_ok(), q_alive, "q bits{bits}");
        assert_eq!(
            g.peptide(u).is_ok(),
            !o.require_parent_match && (!o.require_observation_match || u_match)
        );
        assert_eq!(
            g.parent(b).is_ok(),
            !o.require_parent_group && (!o.require_identified_sequence || q_alive)
        );
        assert_eq!(
            g.parent(unused).is_ok(),
            !o.require_parent_group && !o.require_identified_sequence
        );
        assert_eq!(
            g.parent(rna).is_ok(),
            !o.require_parent_group && (!o.require_identified_sequence || r_alive)
        );
        assert_eq!(g.oligo(r).is_ok(), r_alive);
        assert_eq!(g.observation_match(um).is_ok(), u_match);
        assert_eq!(g.observation_match(rm).is_ok(), r_match);
        assert_eq!(g.observation_match(cm).is_ok(), !o.require_match_group);
        assert_eq!(
            g.compound(c).is_ok(),
            !o.require_observation_match || !o.require_match_group
        );
        assert_eq!(
            g.adduct(adduct).is_ok(),
            !o.require_observation_match || r_match
        );
        assert_eq!(g.observation(empty).is_ok(), !o.require_observation_match);
        assert_eq!(
            g.parent_group(pg, 0).unwrap().parent_refs,
            BTreeSet::from([a])
        );
    }
}
#[test]
fn removal_preserves_surviving_ids_repairs_groups_and_exposes_stale_scores() {
    let mut g = IdentificationData::new().unwrap();
    let a = parent(&mut g, "a", false);
    let b = parent(&mut g, "b", false);
    let p = peptide(&mut g, "A", &[a, b]);
    let q = peptide(&mut g, "C", &[a]);
    let obs = observation(&mut g, "s");
    let pm = g
        .register_observation_match(ObservationMatch::new(p, obs))
        .unwrap();
    let qm = g
        .register_observation_match(ObservationMatch::new(q, obs))
        .unwrap();
    let score = g
        .register_score_type(ScoreType::new("score", true))
        .unwrap();
    let pg = g
        .register_parent_group_set(ParentGroupSet {
            label: "keep empty operation".into(),
            groups: vec![ParentGroup {
                parent_refs: BTreeSet::from([a, b]),
                scores: BTreeMap::from([(score, 0.9)]),
            }],
            ..Default::default()
        })
        .unwrap();
    let mut group = ObservationMatchGroup::new(BTreeSet::from([pm, qm]));
    group.result.add_score(score, 0.8, None).unwrap();
    let mg = g.register_observation_match_group(group).unwrap();
    let report = g.remove_parent_sequences_if(|id, _| id == a).unwrap();
    assert_eq!(
        (
            report.removed_parents,
            report.removed_peptides,
            report.removed_observation_matches,
            report.removed_parent_links
        ),
        (1, 1, 1, 2)
    );
    assert!(report.parent_group_scores_may_be_invalid && report.match_group_scores_may_be_invalid);
    assert!(g.parent(a).is_err() && g.peptide(q).is_err() && g.observation_match(qm).is_err());
    assert_eq!(
        g.peptide(p)
            .unwrap()
            .parent_matches
            .keys()
            .copied()
            .collect::<Vec<_>>(),
        [b]
    );
    assert_eq!(g.parent_group(pg, 0).unwrap().scores[&score], 0.9);
    assert_eq!(
        g.observation_match_group(mg).unwrap().result.score(score),
        Some(0.8)
    );
    g.calculate_coverages(true).unwrap();
    assert_eq!(g.parent(b).unwrap().coverage, 0.25);
    let report = g.remove_observation_matches_if(|id, _| id == pm).unwrap();
    assert_eq!(report.removed_match_groups, 1);
    assert!(g.observation_match_group(mg).is_err());
    assert!(g.parent_group_set(pg).unwrap().groups.is_empty());
    assert_eq!(g.parent_group_count(), 0);
    assert!(g.score_type(score).is_ok() && g.input_file_count() == 1);
}
#[test]
fn no_matching_predicate_does_not_trigger_default_cleanup() {
    let mut g = IdentificationData::new().unwrap();
    let p = parent(&mut g, "orphan", false);
    let u = peptide(&mut g, "A", &[]);
    assert_eq!(
        g.remove_parent_sequences_if(|_, _| false).unwrap(),
        CleanupReport::default()
    );
    assert_eq!(
        g.remove_observation_matches_if(|_, _| false).unwrap(),
        CleanupReport::default()
    );
    assert!(g.parent(p).is_ok() && g.peptide(u).is_ok());
}
#[test]
fn fallible_predicate_failure_rolls_back_earlier_selections_in_source_key_order() {
    let mut g = IdentificationData::new().unwrap();
    let z = parent(&mut g, "z", false);
    let a = parent(&mut g, "a", false);
    let mut visited = Vec::new();
    let result = g.try_remove_parent_sequences_if(|id, value| {
        visited.push(value.accession.clone());
        if id == z {
            Err(Error::InvalidValue("stop".into()))
        } else {
            Ok(id == a)
        }
    });
    assert!(result.is_err());
    assert_eq!(visited, ["a", "z"]);
    assert_eq!(g.parent_count(), 2);
    assert!(g.parent(a).is_ok() && g.parent(z).is_ok());
}
#[test]
fn repeated_deletion_reinsertion_frees_slots_and_never_revives_old_ids() {
    let mut g = IdentificationData::with_limits(GraphLimits {
        max_records: 1,
        max_bytes: 8192,
        ..Default::default()
    })
    .unwrap();
    let first = parent(&mut g, "same", false);
    g.remove_parent_sequences_if(|_, _| true).unwrap();
    for _ in 0..1000 {
        let id = parent(&mut g, "same", false);
        assert_ne!(id, first);
        assert!(g.parent(first).is_err());
        g.remove_parent_sequences_if(|_, _| true).unwrap();
        assert!(g.is_empty());
    }
}
#[test]
fn group_key_collisions_keep_first_original_payload_and_invalidate_later_group_id() {
    let mut g = IdentificationData::new().unwrap();
    let a = parent(&mut g, "a", false);
    let b = parent(&mut g, "b", false);
    let p = peptide(&mut g, "A", &[a]);
    let q = peptide(&mut g, "C", &[b]);
    let obs = observation(&mut g, "s");
    let pm = g
        .register_observation_match(ObservationMatch::new(p, obs))
        .unwrap();
    let qm = g
        .register_observation_match(ObservationMatch::new(q, obs))
        .unwrap();
    let score = g
        .register_score_type(ScoreType::new("score", true))
        .unwrap();
    let pg = g
        .register_parent_group_set(ParentGroupSet {
            groups: vec![
                ParentGroup {
                    parent_refs: BTreeSet::from([a]),
                    scores: BTreeMap::from([(score, 1.)]),
                },
                ParentGroup {
                    parent_refs: BTreeSet::from([a, b]),
                    scores: BTreeMap::from([(score, 2.)]),
                },
            ],
            ..Default::default()
        })
        .unwrap();
    let mut first = ObservationMatchGroup::new(BTreeSet::from([pm]));
    first.result.metadata.insert("first".into(), 1.into());
    let first = g.register_observation_match_group(first).unwrap();
    let second = g
        .register_observation_match_group(ObservationMatchGroup::new(BTreeSet::from([pm, qm])))
        .unwrap();
    let report = g.remove_observation_matches_if(|id, _| id == qm).unwrap();
    assert_eq!(
        (report.removed_parent_groups, report.removed_match_groups),
        (1, 1)
    );
    assert_eq!(g.parent_group(pg, 0).unwrap().scores[&score], 1.);
    assert!(g.observation_match_group(first).is_ok() && g.observation_match_group(second).is_err());
    assert_eq!(
        g.register_observation_match_group(ObservationMatchGroup::new(BTreeSet::from([pm])))
            .unwrap(),
        first
    );
}
#[test]
fn empty_groups_removed_but_grouping_operations_history_and_primary_metadata_retained() {
    let mut g = IdentificationData::new().unwrap();
    let score = g.register_score_type(ScoreType::new("s", true)).unwrap();
    let software = g
        .register_processing_software(ProcessingSoftware::new("tool", "1"))
        .unwrap();
    let step = g
        .register_processing_step(ProcessingStep::new(software), None)
        .unwrap();
    g.set_current_processing_step(step).unwrap();
    let pg = g
        .register_parent_group_set(ParentGroupSet {
            groups: vec![ParentGroup::default()],
            ..Default::default()
        })
        .unwrap();
    let mg = g
        .register_observation_match_group(ObservationMatchGroup::default())
        .unwrap();
    let report = g.cleanup(all_off()).unwrap();
    assert_eq!(
        (report.removed_parent_groups, report.removed_match_groups),
        (1, 1)
    );
    assert!(g.observation_match_group(mg).is_err());
    assert!(g.parent_group_set(pg).unwrap().groups.is_empty());
    assert_eq!(
        g.parent_group_set(pg).unwrap().result.steps_and_scores[0].processing_step,
        Some(step)
    );
    assert_eq!(g.current_processing_step(), Some(step));
    assert!(g.score_type(score).is_ok());
}
#[test]
fn copy_and_merge_translate_only_live_sparse_ids() {
    let mut g = IdentificationData::new().unwrap();
    let gone = parent(&mut g, "gone", false);
    let live = parent(&mut g, "live", false);
    let p = peptide(&mut g, "A", &[live]);
    let obs = observation(&mut g, "s");
    let m = g
        .register_observation_match(ObservationMatch::new(p, obs))
        .unwrap();
    g.remove_parent_sequences_if(|id, _| id == gone).unwrap();
    let (copy, t) = g.try_clone_with_translation().unwrap();
    assert!(t.parent(gone).is_err());
    assert_eq!(
        copy.parent(t.parent(live).unwrap()).unwrap(),
        g.parent(live).unwrap()
    );
    assert!(
        copy.observation_match(t.observation_match(m).unwrap())
            .is_ok()
    );
    let mut merged = IdentificationData::new().unwrap();
    let t = merged.merge_from(&g).unwrap();
    assert!(t.parent(gone).is_err());
    assert!(merged.peptide(t.peptide(p).unwrap()).is_ok());
}
#[test]
fn cleanup_work_failure_is_atomic_after_filter_selection() {
    let mut g = IdentificationData::with_limits(GraphLimits {
        max_work: 15_000,
        ..Default::default()
    })
    .unwrap();
    let a = parent(&mut g, "a", false);
    let b = parent(&mut g, "b", false);
    let p = peptide(&mut g, "A", &[a, b]);
    let obs = observation(&mut g, "s");
    g.register_observation_match(ObservationMatch::new(p, obs))
        .unwrap();
    let before = g.record_count();
    let old = g.peptide(p).unwrap().clone();
    let mut called = 0;
    let error = g
        .remove_parent_sequences_if(|id, _| {
            called += 1;
            id == a
        })
        .unwrap_err();
    assert!(error.to_string().contains("work limit"), "{error}");
    assert_eq!(called, 2);
    assert_eq!(g.record_count(), before);
    assert_eq!(g.peptide(p).unwrap(), &old);
    assert!(g.parent(a).is_ok() && g.parent(b).is_ok());
}

#[test]
fn snapshot_allocation_failure_keeps_selected_records_and_ids() {
    let mut g = IdentificationData::with_limits(GraphLimits {
        max_bytes: 3_500,
        ..Default::default()
    })
    .unwrap();
    let a = parent(&mut g, "a", false);
    let b = parent(&mut g, "b", false);
    let error = g.remove_parent_sequences_if(|id, _| id == a).unwrap_err();
    assert!(error.to_string().contains("allocation limit"), "{error}");
    assert_eq!(g.parent_count(), 2);
    assert!(g.parent(a).is_ok() && g.parent(b).is_ok());
}
