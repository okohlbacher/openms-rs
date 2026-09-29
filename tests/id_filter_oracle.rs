// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//
// IDFilter (PROCESSING/ID/IDFilter.h) members ported in wave 4, against the
// executed source.
//
// Expected values are EXECUTED C++ output (tier 1): ../oracle/B8b-idfilter/
// driver.cpp, compiled with libOpenMS's Release flags against the Linux x86_64
// Release install openms4-release-bc9cc12-c19e494-174b576 and run twice on
// ibminode06, both runs byte-identical. tests/data/id_filter_oracle.tsv holds
// its stdout, each line prefixed with the case name; the inputs
// tests/data/IDFilter_test.idXML and IDFilter_test2.idXML are the class-test
// fixtures at bc9cc12, read here by the Rust idXML reader. Every other input is
// built identically in the driver and in this file.
// tests/data/id_filter_provenance.json records the hashes.
#![cfg(feature = "idxml")]

use openms::analysis::id_filter::{
    self, DigestionFilter, HasDecoyAnnotation, HasMaxMetaValue, HasMetaValue, HasMinMetaValue,
    PeptideDigestionFilter,
};
use openms::chemistry::{AASequence, DigestionSpecificity, Protease, ProteaseDigestion};
use openms::format::fasta::FASTAEntry;
use openms::format::idxml;
use openms::identification::graph::{
    IdentificationData, IdentifiedMolecule, IdentifiedPeptide, InputFile, Observation,
    ObservationMatch, ParentId, ParentMatch, ParentSequence, ScoreType, ScoreTypeId,
};
use openms::identification::{
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinGroup, ProteinHit,
    ProteinIdentification,
};
use openms::kernel::{ConsensusFeature, ConsensusMap};
use openms::metadata::{MetaValue, MetaValueData};

const ORACLE: &str = include_str!("data/id_filter_oracle.tsv");
const FIXTURE: &str = include_str!("data/IDFilter_test.idXML");
const FIXTURE2: &str = include_str!("data/IDFilter_test2.idXML");

/// The oracle rows of one case, split on tabs, the case name removed.
fn case(name: &str) -> Vec<Vec<&'static str>> {
    let rows: Vec<Vec<&str>> = ORACLE
        .lines()
        .filter_map(|line| line.strip_prefix(name)?.strip_prefix('\t'))
        .map(|rest| rest.split('\t').collect())
        .collect();
    assert!(!rows.is_empty(), "oracle has no {name} rows");
    rows
}

fn load(text: &str) -> (Vec<ProteinIdentification>, Vec<PeptideIdentification>) {
    let document = idxml::read(text.as_bytes()).unwrap();
    (
        document.protein_identifications,
        document.peptide_identifications,
    )
}

fn number(text: &str) -> f64 {
    text.parse().unwrap()
}

fn sequences(hits: &[PeptideHit]) -> Vec<String> {
    hits.iter().map(|hit| hit.sequence.to_string()).collect()
}

fn hit(sequence: &str) -> PeptideHit {
    PeptideHit {
        sequence: sequence.parse().unwrap(),
        ..PeptideHit::default()
    }
}

#[test]
fn count_hits_and_best_hit() {
    let (proteins, peptides) = load(FIXTURE);
    let mut four = vec![PeptideIdentification::default(); 4];
    four[0].hits.resize(1, PeptideHit::default());
    four[1].hits.resize(3, PeptideHit::default());
    four[3].hits.resize(2, PeptideHit::default());
    for row in case("count_best") {
        match (row[0], row[1]) {
            ("count", "peptides") => {
                assert_eq!(id_filter::count_hits(&peptides).to_string(), row[2])
            }
            ("count", "proteins") => {
                assert_eq!(id_filter::count_hits(&proteins).to_string(), row[2])
            }
            ("count", "four") => assert_eq!(id_filter::count_hits(&four).to_string(), row[2]),
            ("best", "peptide") => {
                let mut work = peptides.clone();
                if row[2] == "1" {
                    work[0].sort().unwrap();
                }
                work[0].higher_score_better = row[4] == "1";
                let best = id_filter::best_hit(&work, row[3] == "1").unwrap();
                assert_eq!(best.is_some(), row[5] == "1", "{row:?}");
                let best = best.unwrap();
                assert_eq!(best.sequence.to_string(), row[6], "{row:?}");
                assert_eq!(best.score, number(row[7]), "{row:?}");
            }
            ("best", "protein") => {
                let best = id_filter::best_hit(&proteins, row[2] == "1")
                    .unwrap()
                    .unwrap();
                assert_eq!(row[3], "1");
                assert_eq!(best.accession, row[4]);
                assert_eq!(best.score, number(row[5]));
            }
            ("best", "several") => {
                let mut several = vec![PeptideIdentification::default(); 3];
                several[1] = peptides[0].clone();
                several[2] = peptides[0].clone();
                several[2].higher_score_better = false;
                several[2].hits[0].score = 99.0;
                let best = id_filter::best_hit(&several, false).unwrap().unwrap();
                assert_eq!(row[2], "1");
                assert_eq!(best.sequence.to_string(), row[3]);
                assert_eq!(best.score, number(row[4]));
                several[2].score_type = "other".into();
                assert_eq!(
                    case("count_best")
                        .iter()
                        .filter(|r| r[1] == "mismatch")
                        .count(),
                    1
                );
                assert!(id_filter::best_hit(&several, false).is_err());
            }
            ("best", "mismatch") => assert_eq!(row[2], "InvalidValue"),
            ("best", "none") => {
                assert_eq!(row[2], "0");
                let none = vec![PeptideIdentification::default(); 2];
                assert!(id_filter::best_hit(&none, false).unwrap().is_none());
            }
            other => panic!("unexpected row {other:?}"),
        }
    }
}

#[test]
fn protein_groups_by_score_keep_the_source_asymmetry() {
    let probabilities = [0.1, 0.5, 0.7, 0.5, 0.9, 0.0, 1.0];
    let rows = case("groups\tgroups");
    assert_eq!(rows.len(), 8);
    for row in rows {
        let mut groups: Vec<ProteinGroup> = probabilities
            .iter()
            .enumerate()
            .map(|(i, &probability)| ProteinGroup {
                probability,
                accessions: vec![format!("P{i}")],
                ..ProteinGroup::default()
            })
            .collect();
        id_filter::filter_groups_by_score(&mut groups, number(row[0]), row[1] == "1").unwrap();
        let kept: Vec<&str> = groups.iter().map(|g| g.accessions[0].as_str()).collect();
        assert_eq!(
            kept,
            row[2..]
                .iter()
                .copied()
                .filter(|s| !s.is_empty())
                .collect::<Vec<_>>(),
            "{row:?}"
        );
    }
}

#[test]
fn rt_predict_p_value_filter() {
    for row in case("rtp\trtp") {
        let (_, mut peptides) = load(FIXTURE2);
        id_filter::filter_peptides_by_rt_predict_p_value(
            &mut peptides,
            "predicted_RT_p_value",
            number(row[0]),
        )
        .unwrap();
        let expected: Vec<&str> = row[1..].iter().copied().filter(|s| !s.is_empty()).collect();
        assert_eq!(sequences(&peptides[0].hits), expected, "{row:?}");
    }
    let (_, peptides) = load(FIXTURE);
    let mut hits = peptides[0].hits.clone();
    for (i, hit) in hits.iter_mut().enumerate() {
        if i % 3 == 0 {
            hit.metadata
                .insert("p".into(), MetaValue::from((i % 2) as i64));
        } else if i % 3 == 1 {
            hit.metadata.insert(
                "p".into(),
                MetaValue::new(MetaValueData::Float(0.1 * i as f64)).unwrap(),
            );
        }
    }
    for row in case("rtp\trtp_meta") {
        let mut ids = vec![PeptideIdentification {
            hits: hits.clone(),
            ..PeptideIdentification::default()
        }];
        id_filter::filter_peptides_by_rt_predict_p_value(&mut ids, "p", number(row[0])).unwrap();
        assert_eq!(sequences(&ids[0].hits), row[1..].to_vec(), "{row:?}");
    }
    // A string value: the source compares its storage reinterpreted as a
    // double; the port refuses and leaves the hits unchanged.
    let mut ids = vec![PeptideIdentification {
        hits: hits.clone(),
        ..PeptideIdentification::default()
    }];
    ids[0].hits[2]
        .metadata
        .insert("p".into(), MetaValue::from("x"));
    let before = ids.clone();
    assert!(id_filter::filter_peptides_by_rt_predict_p_value(&mut ids, "p", 0.5).is_err());
    assert_eq!(ids, before);
}

#[test]
fn regex_removal_on_unmodified_sequences() {
    let (_, peptides) = load(FIXTURE);
    let mut work = peptides.clone();
    for sequence in ["BBBBB", "JJJJJ", "XXXXX"] {
        work[0].hits.push(hit(sequence));
    }
    for row in case("regex\tregex") {
        let mut ids = work.clone();
        id_filter::remove_peptides_with_matching_regex(&mut ids, row[0]).unwrap();
        let expected: Vec<&str> = row[1..].iter().copied().filter(|s| !s.is_empty()).collect();
        assert_eq!(sequences(&ids[0].hits), expected, "regex {:?}", row[0]);
    }
    assert_eq!(case("regex\tregex_invalid")[0][0], "exception");
    let mut ids = work.clone();
    assert!(id_filter::remove_peptides_with_matching_regex(&mut ids, "(").is_err());
    assert_eq!(ids, work);
}

/// The driver's `best_input`: four identifications built from the fixture's
/// eleven hits in file order.
fn best_input(base: &[PeptideHit]) -> Vec<PeptideIdentification> {
    let p0 = PeptideIdentification {
        identifier: "runA".into(),
        score_type: "Mascot".into(),
        higher_score_better: true,
        hits: base.to_vec(),
        ..PeptideIdentification::default()
    };
    let mut p1 = p0.clone();
    for (i, hit) in p1.hits.iter_mut().enumerate() {
        hit.score = 50.0 - base[i].score;
        hit.charge = if i % 2 == 0 { 2 } else { 3 };
    }
    let mut p2 = p0.clone();
    p2.identifier = "runB".into();
    p2.higher_score_better = false;
    let mut p3 = p0.clone();
    p3.identifier = "runB".into();
    p3.hits.truncate(5);
    vec![p0, p1, p2, p3]
}

fn meta_text(hit: &PeptideHit) -> String {
    match hit.metadata.get(id_filter::BEST_PER_PEPTIDE) {
        None => "-".into(),
        Some(value) => format!("i{}", value.as_i64().unwrap()),
    }
}

/// Compares a dump the driver's `dump_peptides` wrote under `tag`.
fn assert_dump(tag: &str, ids: &[PeptideIdentification]) {
    let rows = case(&format!("best_per_peptide\t{tag}"));
    let mut expected_ids = rows.iter().filter(|r| r[0] == "id");
    let mut expected_hits = rows.iter().filter(|r| r[0] == "hit");
    for (p, id) in ids.iter().enumerate() {
        let row = expected_ids
            .next()
            .unwrap_or_else(|| panic!("{tag}: extra id {p}"));
        assert_eq!(row[1], p.to_string(), "{tag}");
        assert_eq!(row[2], id.identifier, "{tag}");
        assert_eq!(row[3], id.hits.len().to_string(), "{tag} id {p}");
        for (i, hit) in id.hits.iter().enumerate() {
            let row = expected_hits.next().unwrap();
            assert_eq!(
                (row[1], row[2]),
                (p.to_string().as_str(), i.to_string().as_str()),
                "{tag}"
            );
            assert_eq!(row[3], hit.sequence.to_string(), "{tag} {p}/{i}");
            assert_eq!(row[4], hit.charge.to_string(), "{tag} {p}/{i}");
            assert_eq!(number(row[5]), hit.score, "{tag} {p}/{i}");
            assert_eq!(row[6], meta_text(hit), "{tag} {p}/{i}");
        }
    }
    assert!(expected_ids.next().is_none(), "{tag}: missing ids");
    assert!(expected_hits.next().is_none(), "{tag}: missing hits");
}

#[test]
fn best_per_peptide_annotation_and_filtering() {
    let (_, peptides) = load(FIXTURE);
    let base = peptides[0].hits.clone();
    let combos = [
        (false, false, 0),
        (true, false, 0),
        (false, true, 0),
        (true, true, 0),
        (false, false, 1),
        (false, false, 3),
        (true, true, 2),
    ];
    for (mods, charges, nr) in combos {
        let tag = format!("{}{}{}", u8::from(mods), u8::from(charges), nr);
        let mut ids = best_input(&base);
        id_filter::annotate_best_per_peptide(&mut ids, mods, charges, nr).unwrap();
        assert_dump(&format!("annotate_{tag}"), &ids);
        let mut ids = best_input(&base);
        id_filter::keep_best_per_peptide(&mut ids, mods, charges, nr).unwrap();
        assert_dump(&format!("keep_{tag}"), &ids);
        let mut ids = best_input(&base);
        id_filter::annotate_best_per_peptide_per_run(&mut ids, mods, charges, nr).unwrap();
        assert_dump(&format!("annotate_run_{tag}"), &ids);
        let mut ids = best_input(&base);
        id_filter::keep_best_per_peptide_per_run(&mut ids, mods, charges, nr).unwrap();
        assert_dump(&format!("keep_run_{tag}"), &ids);

        let input = best_input(&base);
        let mut map = ConsensusMap::default();
        let mut f0 = ConsensusFeature::default();
        f0.peptide_identifications = vec![input[0].clone(), input[1].clone()];
        let mut f1 = ConsensusFeature::default();
        f1.peptide_identifications = vec![input[2].clone()];
        map.features = vec![f0, f1];
        map.unassigned_peptide_identifications = vec![input[3].clone()];
        let mut keep = map.clone();
        id_filter::annotate_best_per_peptide_per_run_in_map(&mut map, mods, charges, nr).unwrap();
        assert_dump(
            &format!("map_annotate_{tag}_f0"),
            &map.features[0].peptide_identifications,
        );
        assert_dump(
            &format!("map_annotate_{tag}_f1"),
            &map.features[1].peptide_identifications,
        );
        assert_dump(
            &format!("map_annotate_{tag}_u"),
            &map.unassigned_peptide_identifications,
        );
        id_filter::keep_best_per_peptide_per_run_in_map(&mut keep, mods, charges, nr).unwrap();
        assert_dump(
            &format!("map_keep_{tag}_f0"),
            &keep.features[0].peptide_identifications,
        );
        assert_dump(
            &format!("map_keep_{tag}_f1"),
            &keep.features[1].peptide_identifications,
        );
        assert_dump(
            &format!("map_keep_{tag}_u"),
            &keep.unassigned_peptide_identifications,
        );
    }
}

#[test]
fn meta_value_predicates_and_move() {
    let (_, peptides) = load(FIXTURE);
    let mut hits = peptides[0].hits.clone();
    for (i, hit) in hits.iter_mut().enumerate() {
        let value = match i % 4 {
            0 => Some(MetaValue::new(MetaValueData::Float(0.5 * i as f64)).unwrap()),
            1 => Some(MetaValue::from(i as i64)),
            2 => Some(MetaValue::from(format!("s{i}"))),
            _ => None,
        };
        if let Some(value) = value {
            hit.metadata.insert("x".into(), value);
        }
        if i == 4 {
            hit.metadata.insert(
                "x".into(),
                MetaValue::new(MetaValueData::Float(2.0000004)).unwrap(),
            );
        }
    }
    let numeric: Vec<PeptideHit> = hits
        .iter()
        .enumerate()
        .filter(|(i, _)| i % 4 != 2)
        .map(|(_, hit)| hit.clone())
        .collect();
    let float = |v: f64| MetaValue::new(MetaValueData::Float(v)).unwrap();
    for row in case("meta_predicates\tmeta") {
        let expected: Vec<&str> = row[1..].iter().copied().filter(|s| !s.is_empty()).collect();
        let kept = match row[0] {
            "has_any" => keep(&hits, |h| {
                HasMetaValue::new("x", MetaValue::default()).matches(h)
            }),
            "has_2.0" => keep(&hits, |h| HasMetaValue::new("x", float(2.0)).matches(h)),
            "has_int5" => keep(&hits, |h| {
                HasMetaValue::new("x", MetaValue::from(5i64)).matches(h)
            }),
            "has_double5" => keep(&hits, |h| HasMetaValue::new("x", float(5.0)).matches(h)),
            "has_s2" => keep(&hits, |h| {
                HasMetaValue::new("x", MetaValue::from("s2")).matches(h)
            }),
            "not_decoy" => {
                let mut k = hits.clone();
                k[0].metadata
                    .insert("target_decoy".into(), MetaValue::from("decoy"));
                k[1].metadata
                    .insert("isDecoy".into(), MetaValue::from("true"));
                k[2].metadata
                    .insert("target_decoy".into(), MetaValue::from("target"));
                k[2].metadata
                    .insert("isDecoy".into(), MetaValue::from("true"));
                k[3].metadata
                    .insert("isDecoy".into(), MetaValue::from("false"));
                id_filter::remove_matching_items(&mut k, |h| HasDecoyAnnotation.matches(h));
                k
            }
            name if name.starts_with("max_") => {
                let predicate = HasMaxMetaValue {
                    key: "x".into(),
                    value: number(&name[4..]),
                };
                let mut k = numeric.clone();
                id_filter::try_keep_matching_items(&mut k, |h| predicate.matches(h)).unwrap();
                k
            }
            name if name.starts_with("min_") => {
                let predicate = HasMinMetaValue {
                    key: "x".into(),
                    value: number(&name[4..]),
                };
                let mut k = numeric.clone();
                id_filter::try_keep_matching_items(&mut k, |h| predicate.matches(h)).unwrap();
                k
            }
            other => panic!("unexpected meta row {other}"),
        };
        assert_eq!(sequences(&kept), expected, "{}", row[0]);
    }
    // The numeric predicates refuse a string value (the source reads garbage).
    let predicate = HasMaxMetaValue {
        key: "x".into(),
        value: 2.0,
    };
    assert!(predicate.matches(&hits[2]).is_err());

    let rows = case("meta_predicates\tmove");
    let mut numbers: Vec<i32> = (0..8).collect();
    let mut target = vec![100];
    id_filter::move_matching_items(&mut numbers, |i| i % 3 == 0, &mut target);
    let text = |v: &[i32]| v.iter().map(i32::to_string).collect::<Vec<_>>();
    assert_eq!(rows[0][0], "kept");
    assert_eq!(text(&numbers), rows[0][1..].to_vec());
    assert_eq!(rows[1][0], "moved");
    assert_eq!(text(&target), rows[1][1..].to_vec());
}

fn keep(hits: &[PeptideHit], predicate: impl FnMut(&PeptideHit) -> bool) -> Vec<PeptideHit> {
    let mut kept = hits.to_vec();
    id_filter::keep_matching_items(&mut kept, predicate);
    kept
}

#[test]
fn unassigned_proteins_are_moved_out_in_source_order() {
    let run = |id: &str, accessions: &[&str]| ProteinIdentification {
        identifier: id.into(),
        hits: accessions
            .iter()
            .map(|a| ProteinHit {
                accession: (*a).into(),
                ..ProteinHit::default()
            })
            .collect(),
        ..ProteinIdentification::default()
    };
    let peptide = |id: &str, accessions: &[&str]| PeptideIdentification {
        identifier: id.into(),
        hits: vec![PeptideHit {
            sequence: "PEPTIDE".parse().unwrap(),
            evidences: accessions
                .iter()
                .map(|a| PeptideEvidence {
                    protein_accession: (*a).into(),
                    ..PeptideEvidence::default()
                })
                .collect(),
            ..PeptideHit::default()
        }],
        ..PeptideIdentification::default()
    };
    let mut map = ConsensusMap::default();
    map.protein_identifications = vec![
        run("run1", &["A", "B", "C", "D", "E"]),
        run("run2", &["F", "G"]),
        run("run1", &["H", "A", "C", "I"]),
        run("run3", &["J"]),
    ];
    let mut f0 = ConsensusFeature::default();
    f0.peptide_identifications = vec![peptide("run1", &["A", "C"])];
    let mut f1 = ConsensusFeature::default();
    f1.peptide_identifications = vec![peptide("run2", &["G"])];
    map.features = vec![f0, f1];
    map.unassigned_peptide_identifications = vec![peptide("run1", &["E"])];
    let result = id_filter::extract_unassigned_proteins(&mut map).unwrap();
    let rows = case("unassigned\tunassigned");
    assert_eq!(rows.len(), result.len());
    for (row, (id, hits)) in rows.iter().zip(&result) {
        assert_eq!(row[0], id);
        let accessions: Vec<&str> = hits.iter().map(|h| h.accession.as_str()).collect();
        assert_eq!(
            accessions,
            row[1..]
                .iter()
                .copied()
                .filter(|s| !s.is_empty())
                .collect::<Vec<_>>()
        );
    }
    let rows = case("unassigned\tremaining");
    for (row, protein) in rows.iter().zip(&map.protein_identifications) {
        assert_eq!(row[0], protein.identifier);
        let accessions: Vec<&str> = protein.hits.iter().map(|h| h.accession.as_str()).collect();
        assert_eq!(
            accessions,
            row[1..]
                .iter()
                .copied()
                .filter(|s| !s.is_empty())
                .collect::<Vec<_>>()
        );
    }
}

#[test]
fn peptide_digestion_filter_counts_as_the_source_tokenizer() {
    let sequences = [
        "DFPIANGER",
        "DFPIAN(Deamidated)GER",
        "DFKPIARN(Deamidated)GER",
        "DFPKIARNGER",
        "K(Dimethyl)FPIAUGR",
        "KRKRKR",
        "PEPTIDEK",
        "A",
        "KPKPRP",
        "RRRR",
    ];
    let bound = |text: &str| -> Option<i32> {
        match text.parse::<i32>().unwrap() {
            -1 => None,
            value => Some(value),
        }
    };
    let rows = case("digestion_peptide\tdigest_pep");
    assert_eq!(rows.len(), 35);
    for row in rows {
        let digestion = ProteaseDigestion {
            enzyme: Protease::from_name(row[0]).unwrap(),
            ..ProteaseDigestion::default()
        };
        let filter = PeptideDigestionFilter::new(&digestion, bound(row[1]), bound(row[2]));
        for (column, sequence) in sequences.iter().enumerate() {
            let hit = hit(sequence);
            assert_eq!(
                filter.matches(&hit).unwrap(),
                row[3 + column] == "1",
                "{} [{}, {}] {sequence}",
                row[0],
                row[1],
                row[2]
            );
        }
    }
}

#[test]
fn digestion_filter_checks_evidence_limits_as_the_source() {
    let fasta = vec![
        FASTAEntry::new("P1", "", "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFK"),
        FASTAEntry::new("P2", "", "MPEPTIDEKAAAKRPPPDPGGR"),
        FASTAEntry::new("P2", "", "APEPTIDEKAAAKRPPPDPGGR"),
    ];
    let evidences: [(&str, i64, i64); 24] = [
        ("P1", 0, 1),
        ("P1", 2, 18),
        ("P1", 2, 19),
        ("P1", 1, 18),
        ("P1", 19, 22),
        ("P1", 24, 27),
        ("P1", 2, 27),
        ("P1", 36, 43),
        ("P1", 37, 43),
        ("P1", 1, 1),
        ("P1", 2, 2),
        ("P1", 0, 43),
        ("P1", 40, 50),
        ("P2", 1, 8),
        ("P2", 0, 8),
        ("P2", 9, 12),
        ("P2", 14, 21),
        ("P2", 13, 21),
        ("P3", 0, 5),
        ("", 0, 5),
        ("P1", -1, 5),
        ("P1", 2, -1),
        ("P1", 0, 0),
        ("P1", 5, 5),
    ];
    let position = |p: i64| usize::try_from(p).ok();
    let rows = case("digestion_evidence\tdigest_ev");
    assert_eq!(rows.len(), 8);
    for row in rows {
        let digestion = ProteaseDigestion {
            enzyme: Protease::Trypsin,
            missed_cleavages: 1,
            specificity: row[2].parse::<DigestionSpecificity>().unwrap(),
            ..ProteaseDigestion::default()
        };
        let filter = DigestionFilter::new(&fasta, &digestion, row[0] == "1", row[1] == "1");
        for (column, &(accession, start, end)) in evidences.iter().enumerate() {
            let evidence = PeptideEvidence {
                protein_accession: accession.into(),
                start: position(start),
                end: position(end),
                ..PeptideEvidence::default()
            };
            assert_eq!(
                filter.matches(&evidence).unwrap(),
                row[3 + column] == "1",
                "{row:?} evidence {column} {accession} {start}-{end}"
            );
        }
    }
    // The source's length `end - start` treats the inclusive end as exclusive:
    // the tryptic peptide WVTFISLLLLFSSAYSR (residues 2..=18) is rejected,
    // 2..=19 accepted.
    let digestion = ProteaseDigestion {
        enzyme: Protease::Trypsin,
        ..ProteaseDigestion::default()
    };
    let filter = DigestionFilter::new(&fasta, &digestion, false, false);
    let peptide = AASequence::parse("WVTFISLLLLFSSAYSR").unwrap();
    assert_eq!(peptide.len(), 17);
    let evidence = PeptideEvidence::new("P1", 2..=18).unwrap();
    assert!(!filter.matches(&evidence).unwrap());
}

/// The driver's `idd_*` graph: three spectra, four peptides (one target, two
/// decoy, one shared parent), seven matches scored under `sc` (higher better)
/// and `low` (lower better), some unscored.
fn identification_graph() -> (IdentificationData, ScoreTypeId, ScoreTypeId) {
    let mut data = IdentificationData::new().unwrap();
    let file = data.register_input_file(InputFile::new("file")).unwrap();
    let mut observation = |id: &str| {
        data.register_observation(Observation::new(id, file))
            .unwrap()
    };
    let s1 = observation("spectrum_1");
    let s2 = observation("spectrum_2");
    let s3 = observation("spectrum_3");
    let sc = data
        .register_score_type(ScoreType::new("sc", true))
        .unwrap();
    let low = data
        .register_score_type(ScoreType::new("low", false))
        .unwrap();
    let pt = data
        .register_parent_sequence(ParentSequence::new("PT"))
        .unwrap();
    let mut decoy = ParentSequence::new("PD");
    decoy.is_decoy = true;
    let pd = data.register_parent_sequence(decoy).unwrap();
    let mut peptide = |sequence: &str, parents: &[ParentId]| {
        let mut value = IdentifiedPeptide::new(sequence.parse().unwrap());
        for &parent in parents {
            value
                .parent_matches
                .entry(parent)
                .or_default()
                .insert(ParentMatch::new(Some(0), Some(sequence.len() - 1)));
        }
        data.register_identified_peptide(value).unwrap()
    };
    let ra = peptide("PEPTIDEK", &[pt]);
    let rb = peptide("AAAK", &[pd]);
    let rc = peptide("PEPTIDE", &[pt, pd]);
    let rd = peptide("GGGR", &[pd]);
    let matches = [
        (ra, s1, Some(10.0), Some(3.0)),
        (rb, s1, Some(12.0), Some(1.0)),
        (rc, s1, None, None),
        (ra, s2, Some(5.0), None),
        (rb, s2, Some(5.0), Some(2.0)),
        (rc, s3, None, None),
        (rd, s3, None, None),
    ];
    for (molecule, observation, sc_value, low_value) in matches {
        let mut value = ObservationMatch::new(molecule, observation);
        if let Some(v) = sc_value {
            value.result.add_score(sc, v, None).unwrap();
        }
        if let Some(v) = low_value {
            value.result.add_score(low, v, None).unwrap();
        }
        data.register_observation_match(value).unwrap();
    }
    (data, sc, low)
}

/// The driver's dump, sorted: the source iterates its graph containers in
/// the heap-address order of the referenced records (the oracle lists
/// spectrum_3 before spectrum_1), an allocator artefact the port does not
/// reproduce; the records must be the same.
fn graph_dump(data: &IdentificationData, sc: ScoreTypeId, low: ScoreTypeId) -> Vec<String> {
    let score = |value: Option<f64>| value.map_or("-".to_owned(), |v| format!("{v}"));
    let mut lines = Vec::new();
    for (_, value) in data.observation_matches() {
        let IdentifiedMolecule::Peptide(peptide) = value.identified_molecule else {
            panic!("peptide matches only");
        };
        lines.push(format!(
            "match\t{}\t{}\t{}\t{}",
            data.observation(value.observation).unwrap().data_id,
            data.peptide(peptide).unwrap().sequence,
            score(value.result.score(sc)),
            score(value.result.score(low)),
        ));
    }
    for (_, parent) in data.parents() {
        lines.push(format!("parent\t{}", parent.accession));
    }
    for (_, peptide) in data.peptides() {
        lines.push(format!("peptide\t{}", peptide.sequence));
    }
    for (_, observation) in data.observations() {
        lines.push(format!("observation\t{}", observation.data_id));
    }
    lines.sort();
    lines
}

#[test]
fn identification_data_filters_match_the_source() {
    for name in [
        "idd_initial",
        "idd_best_sc",
        "idd_best_low",
        "idd_cutoff_sc_10",
        "idd_cutoff_sc_5",
        "idd_cutoff_low_2",
        "idd_decoys",
    ] {
        let mut expected: Vec<String> = case(name).iter().map(|row| row.join("\t")).collect();
        expected.sort();
        let (mut data, sc, low) = identification_graph();
        match name {
            "idd_best_sc" => id_filter::keep_best_match_per_observation(&mut data, sc).unwrap(),
            "idd_best_low" => id_filter::keep_best_match_per_observation(&mut data, low).unwrap(),
            "idd_cutoff_sc_10" => {
                id_filter::filter_observation_matches_by_score(&mut data, sc, 10.0).unwrap()
            }
            "idd_cutoff_sc_5" => {
                id_filter::filter_observation_matches_by_score(&mut data, sc, 5.0).unwrap()
            }
            "idd_cutoff_low_2" => {
                id_filter::filter_observation_matches_by_score(&mut data, low, 2.0).unwrap()
            }
            "idd_decoys" => id_filter::remove_decoys(&mut data).unwrap(),
            _ => {}
        }
        assert_eq!(graph_dump(&data, sc, low), expected, "{name}");
    }
}
