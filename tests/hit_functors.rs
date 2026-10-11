// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of the functor, hash and assignment sections of `PeptideHit_test.cpp`
//! and `ProteinHit_test.cpp` (core SDK bc9cc12), cited by source line. The
//! record sections are ported in `tests/identification.rs`. Literals are tier
//! 3, transcribed from the class tests; values derived from the source instead
//! say so.

use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasher, BuildHasherDefault, DefaultHasher, Hash};

use openms::chemistry::AASequence;
use openms::identification::hits::{SequenceChargeKey, SpectrumMatch};
use openms::identification::{PeakAnnotation, PeptideHit, ProteinHit};
use openms::metadata::MetaInfo;

fn hash_of<T: Hash + ?Sized>(value: &T) -> u64 {
    BuildHasherDefault::<DefaultHasher>::default().hash_one(value)
}

fn peptide(score: f64, rank: u32, charge: i32, sequence: &str) -> PeptideHit {
    PeptideHit::new(score, rank, charge, AASequence::parse(sequence).unwrap()).unwrap()
}

fn protein(score: f64, rank: u32, accession: &str, sequence: &str) -> ProteinHit {
    ProteinHit::new(score, rank, accession, sequence).unwrap()
}

// PeptideHit_test.cpp L226 [PeptideHit::ScoreLess]
#[test]
fn peptide_section_score_less() {
    let a = PeptideHit {
        score: 10.0,
        ..PeptideHit::default()
    };
    let b = PeptideHit {
        score: 20.0,
        ..PeptideHit::default()
    };
    assert!(PeptideHit::score_less(&a, &b));
    assert!(!PeptideHit::score_less(&b, &a));
    assert!(!PeptideHit::score_less(&a, &a));
}

// PeptideHit_test.cpp L238 [PeptideHit::ScoreMore]
#[test]
fn peptide_section_score_more() {
    let a = PeptideHit {
        score: 20.0,
        ..PeptideHit::default()
    };
    let b = PeptideHit {
        score: 10.0,
        ..PeptideHit::default()
    };
    assert!(PeptideHit::score_more(&a, &b));
    assert!(!PeptideHit::score_more(&b, &a));
    assert!(!PeptideHit::score_more(&a, &a));
}

/// `RankLess` and `SequenceLessComparator` have no class-test section; derived
/// from their definitions (PeptideHit.h), `a.getRank() < b.getRank()` and
/// `a.getSequence().toString() < b.getSequence().toString()`.
#[test]
fn peptide_rank_and_sequence_predicates() {
    let first = peptide(1.0, 1, 2, "PEPTIDE");
    let second = peptide(1.0, 2, 2, "AAA");
    assert!(PeptideHit::rank_less(&first, &second));
    assert!(!PeptideHit::rank_less(&second, &first));
    assert!(!PeptideHit::rank_less(&first, &first));
    assert!(PeptideHit::sequence_less(&second, &first));
    assert!(!PeptideHit::sequence_less(&first, &second));
    let oxidised = peptide(1.0, 1, 2, "PEPTM(Oxidation)IDE");
    let plain = peptide(1.0, 1, 2, "PEPTMIDE");
    // "PEPTM(Oxidation)IDE" < "PEPTMIDE" byte-wise, since '(' < 'I'.
    assert!(PeptideHit::sequence_less(&oxidised, &plain));
}

// PeptideHit_test.cpp L669 [EXTRA] PeptideHit::SequenceChargeHash and SequenceChargeEqual
#[test]
fn peptide_section_extra_sequence_charge_hash() {
    let hit1 = peptide(1.0, 1, 2, "PEPTIDE");
    let hit2 = peptide(2.0, 2, 2, "PEPTIDE");
    let hit3 = peptide(1.0, 1, 3, "PEPTIDE");
    let hit4 = peptide(1.0, 1, 2, "PEPTIDER");

    let key = SequenceChargeKey;
    assert_eq!(hash_of(&key(&hit1)), hash_of(&key(&hit2)));
    assert!(key(&hit1) == key(&hit2));
    assert_ne!(hash_of(&key(&hit1)), hash_of(&key(&hit3)));
    assert!(key(&hit1) != key(&hit3));
    assert_ne!(hash_of(&key(&hit1)), hash_of(&key(&hit4)));
    assert!(key(&hit1) != key(&hit4));

    let hit5 = peptide(1.0, 1, 2, "PEPTM(Oxidation)IDE");
    let hit6 = peptide(1.0, 1, 2, "PEPTMIDE");
    assert_ne!(hash_of(&key(&hit5)), hash_of(&key(&hit6)));

    let mut hit_set = HashSet::new();
    hit_set.insert(key(&hit1));
    hit_set.insert(key(&hit2));
    hit_set.insert(key(&hit3));
    hit_set.insert(key(&hit4));
    assert_eq!(hit_set.len(), 3);
    let query = peptide(99.0, 99, 2, "PEPTIDE");
    assert!(hit_set.contains(&key(&query)));
}

/// `std::hash<PeptideHit>` has no class-test section; derived from its
/// "consistent with operator==" contract.
#[test]
fn peptide_hash_follows_equality() {
    let hit1 = peptide(1.0, 1, 2, "PEPTIDE");
    let hit2 = hit1.clone();
    assert_eq!(hash_of(&hit1), hash_of(&hit2));
    let mut hit3 = hit1.clone();
    hit3.score = 2.0;
    assert_ne!(hit1, hit3);
    assert_ne!(hash_of(&hit1), hash_of(&hit3));
    let mut negative_zero = hit1.clone();
    negative_zero.score = 0.0;
    let mut positive_zero = hit1.clone();
    positive_zero.score = -0.0;
    assert_eq!(negative_zero, positive_zero);
    assert_eq!(hash_of(&negative_zero), hash_of(&positive_zero));

    let annotation = PeakAnnotation {
        mz: 100.5,
        intensity: 0.25,
        charge: 1,
        annotation: "y1".into(),
    };
    assert_eq!(hash_of(&annotation), hash_of(&annotation.clone()));
    // `SpectrumMatch` is the same type.
    let alias: SpectrumMatch = hit1.clone();
    assert_eq!(alias, hit1);
}

/// `operator<<` (PeptideHit.cpp:422-427, ProteinHit.cpp:194-198), derived: the
/// score goes through `StringUtils::toStr(double)`, which trims trailing zeros.
#[test]
fn stream_output_of_both_hits() {
    let hit = peptide(0.5, 1, 2, "PEPTIDE");
    assert_eq!(
        hit.to_string(),
        "peptide hit with sequence 'PEPTIDE', charge 2, score 0.5"
    );
    let protein_hit = protein(4.5, 1, "ACC123", "PEPTIDE");
    assert_eq!(
        protein_hit.to_string(),
        "protein hit with accession 'ACC123', score 4.5"
    );
}

/// `PeakAnnotation::writePeakAnnotationsString_` (PeptideHit.cpp:443-457),
/// derived: stable sort by (m/z, charge, annotation, intensity), `toStr` on the
/// numbers, backslash-escaped quoting, `|` between entries.
#[test]
fn peak_annotation_string_is_sorted_and_quoted() {
    let annotations = vec![
        PeakAnnotation {
            mz: 200.25,
            intensity: 0.5,
            charge: 1,
            annotation: "y2".into(),
        },
        PeakAnnotation {
            mz: 100.5,
            intensity: 0.25,
            charge: 2,
            annotation: "[alpha|ci$b1]".into(),
        },
        PeakAnnotation {
            mz: 100.5,
            intensity: 0.75,
            charge: 1,
            annotation: "say \"hi\"".into(),
        },
    ];
    assert_eq!(
        PeakAnnotation::write_peak_annotations_string(&annotations).unwrap(),
        "100.5,0.75,1,\"say \\\"hi\\\"\"|100.5,0.25,2,\"[alpha|ci$b1]\"|200.25,0.5,1,\"y2\""
    );
    assert_eq!(
        PeakAnnotation::write_peak_annotations_string(&[]).unwrap(),
        ""
    );
    let bad = [PeakAnnotation {
        mz: f64::NAN,
        ..annotations[0].clone()
    }];
    assert!(PeakAnnotation::write_peak_annotations_string(&bad).is_err());
}

// ProteinHit_test.cpp L97 ProteinHit& operator=(const MetaInfoInterface& source)
#[test]
fn protein_section_assign_meta_info() {
    // The class test's globals: score 4.4f, rank 3, accession "PROOE34",
    // sequence "ARRAY"; the float score widens to double.
    let score = f64::from(4.4f32);
    let mut hit = protein(score, 3, "PROOE34", "ARRAY");
    // 123.123 is outside the percentage `validate` enforces; the field is
    // public and assigned unchecked, as `setCoverage` is.
    hit.coverage = Some(123.123);
    let mut meta = MetaInfo::new();
    meta.insert("label".into(), 17i32.into());
    hit.assign_meta_info(&meta);
    assert_eq!(hit.score, score);
    assert_eq!(hit.rank, 3);
    assert_eq!(hit.accession, "PROOE34");
    assert_eq!(hit.sequence, "ARRAY");
    assert_eq!(hit.coverage, Some(123.123));
    assert_eq!(hit.metadata.get("label"), Some(&17i32.into()));
}

// ProteinHit_test.cpp L243 [ProteinHit::ScoreLess]
#[test]
fn protein_section_score_less() {
    let a = ProteinHit {
        score: 10.0,
        ..ProteinHit::default()
    };
    let b = ProteinHit {
        score: 20.0,
        ..ProteinHit::default()
    };
    assert!(ProteinHit::score_less(&a, &b));
    assert!(!ProteinHit::score_less(&b, &a));
    assert!(!ProteinHit::score_less(&a, &a));
}

// ProteinHit_test.cpp L255 [ProteinHit::ScoreMore]
#[test]
fn protein_section_score_more() {
    let a = ProteinHit {
        score: 20.0,
        ..ProteinHit::default()
    };
    let b = ProteinHit {
        score: 10.0,
        ..ProteinHit::default()
    };
    assert!(ProteinHit::score_more(&a, &b));
    assert!(!ProteinHit::score_more(&b, &a));
    assert!(!ProteinHit::score_more(&a, &a));
}

/// Derived from ProteinHit.h: equal scores fall back to the accession, and
/// `ProteinHitAccessionHash` looks at the accession only.
#[test]
fn protein_accession_tiebreak_and_accession_hash() {
    let a = protein(1.0, 1, "A", "PEPTIDE");
    let b = protein(1.0, 2, "B", "PROTEIN");
    assert!(ProteinHit::score_less(&a, &b));
    assert!(ProteinHit::score_more(&b, &a));
    let same_accession = protein(9.0, 7, "A", "OTHER");
    assert_eq!(a.accession_hash(), same_accession.accession_hash());
    assert_ne!(a.accession_hash(), b.accession_hash());
    let pointer: &ProteinHit = &a;
    assert_eq!(pointer.accession_hash(), a.accession_hash());
}

// ProteinHit_test.cpp L336 [EXTRA] std::hash<ProteinHit>
#[test]
fn protein_section_extra_hash() {
    let with_coverage = |mut hit: ProteinHit, coverage: f64| {
        hit.coverage = Some(coverage);
        hit
    };
    let hit1 = with_coverage(protein(4.5, 1, "ACC123", "PEPTIDE"), 45.5);
    let hit2 = with_coverage(protein(4.5, 1, "ACC123", "PEPTIDE"), 45.5);
    assert!(hit1 == hit2);
    assert_eq!(hash_of(&hit1), hash_of(&hit2));
    let hit3 = with_coverage(protein(4.5, 1, "ACC456", "PEPTIDE"), 45.5);
    assert_ne!(hash_of(&hit1), hash_of(&hit3));
    let hit4 = with_coverage(protein(5.5, 1, "ACC123", "PEPTIDE"), 45.5);
    assert_ne!(hash_of(&hit1), hash_of(&hit4));
    let hit5 = with_coverage(protein(4.5, 2, "ACC123", "PEPTIDE"), 45.5);
    assert_ne!(hash_of(&hit1), hash_of(&hit5));
    let hit6 = with_coverage(protein(4.5, 1, "ACC123", "PROTEIN"), 45.5);
    assert_ne!(hash_of(&hit1), hash_of(&hit6));
    let hit7 = with_coverage(protein(4.5, 1, "ACC123", "PEPTIDE"), 50.0);
    assert_ne!(hash_of(&hit1), hash_of(&hit7));

    // unordered_set / unordered_map: `ProteinHit` holds floating values and is
    // not `Eq`, so bucket by hash and compare with `==` by hand.
    let mut buckets: HashMap<u64, Vec<(ProteinHit, &str)>> = HashMap::new();
    let mut insert = |hit: &ProteinHit, value: &'static str| {
        let bucket = buckets.entry(hash_of(hit)).or_default();
        match bucket.iter_mut().find(|(key, _)| key == hit) {
            Some(slot) => slot.1 = value,
            None => bucket.push((hit.clone(), value)),
        }
    };
    insert(&hit1, "first");
    insert(&hit2, "first");
    insert(&hit3, "second");
    assert_eq!(buckets.values().map(Vec::len).sum::<usize>(), 2);
    let lookup = |hit: &ProteinHit| {
        buckets
            .get(&hash_of(hit))
            .and_then(|bucket| bucket.iter().find(|(key, _)| key == hit))
            .map(|(_, value)| *value)
    };
    assert_eq!(lookup(&hit1), Some("first"));
    assert_eq!(lookup(&hit2), Some("first"));
    assert_eq!(lookup(&hit3), Some("second"));
}
