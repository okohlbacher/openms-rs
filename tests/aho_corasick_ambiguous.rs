// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

//! `ANALYSIS/ID/AhoCorasickAmbiguous.h`.
//!
//! Every section of the upstream class test
//! `src/tests/class_tests/openms/source/AhoCorasickAmbiguous_test.cpp` at the
//! core pin `bc9cc12`, with its literals. `test_case` is the class test's
//! `testCase`: it collects `needle@query_pos` for every hit of every
//! `nextHits` call and compares the sorted list with the sorted expectation,
//! so multiplicities count and order does not. Tests marked "native" check
//! port-only behaviour.

use openms::analysis::id::aho_corasick_ambiguous::{
    AA, AA_TO_CHAR, ACNode, ACScout, ACTrie, ACTrieState, Bitset, CHAR_TO_AA, Hit, Index,
    MAX_NEEDLE_LENGTH, next_valid_aa,
};

fn test_case(t: &ACTrie, protein: &str, expected: &str, needles: &[&str]) {
    let mut observed = Vec::new();
    let mut state = ACTrieState::default();
    state.set_query(protein);
    while t.next_hits(&mut state).unwrap() {
        for hit in &state.hits {
            observed.push(format!(
                "{}@{}",
                needles[hit.needle_index as usize], hit.query_pos
            ));
        }
    }
    let cleaned: String = expected.chars().filter(|c| !c.is_whitespace()).collect();
    let mut expected: Vec<String> = if cleaned.is_empty() {
        Vec::new()
    } else {
        cleaned.split(',').map(str::to_owned).collect()
    };
    expected.sort();
    observed.sort();
    assert_eq!(observed, expected, "protein {protein}");
}

fn trie(max_aaa: u32, max_mm: u32, needles: &[&str]) -> ACTrie {
    let mut t = ACTrie::new(max_aaa, max_mm);
    t.add_needles_and_compress(needles).unwrap();
    t
}

#[test]
fn construction() {
    let t = ACTrie::new(1, 0);
    assert_eq!(t.size(), 1);
    assert_eq!(t.get_needle_count(), 0);
    assert_eq!(ACTrie::default().get_max_aaa_count(), 0);
    assert_eq!(ACTrie::default().get_max_mm_count(), 0);
}

#[test]
fn add_needle() {
    let mut t = ACTrie::new(1, 2);
    for needle in ["WITHV", "WITHB", "WITHJ", "WITHZ", "WITHX"] {
        t.add_needle(needle).unwrap();
    }
    assert!(t.add_needle("WITH$").is_err());
    assert!(t.add_needle("WITH*").is_err());
    assert_eq!(t.get_needle_count(), 5);
    assert_eq!(t.get_max_aaa_count(), 1);
    assert_eq!(t.get_max_mm_count(), 2);
}

#[test]
fn add_needles() {
    let mut t = ACTrie::new(1, 2);
    for needle in ["WITHV", "WITHB", "WITHJ"] {
        t.add_needle(needle).unwrap();
    }
    assert!(t.add_needles(&["WITHZ", "WITHX", "WITH$"]).is_err());
    assert_eq!(t.get_needle_count(), 5);
}

#[test]
fn next_hits() {
    let needles = ["acd", "adc", "cad", "cda", "dac", "dca"];
    let mut t = trie(0, 0, &needles);
    test_case(
        &t,
        "acdIadcIcadIcdaIdacIdca",
        "acd@0,  adc@4,  cad@8,  cda@12,  dac@16,  dca@20",
        &needles,
    );
    t.set_max_aaa_count(3);
    test_case(
        &t,
        "acdIadcIcadIcdaIdacIdca",
        "acd@0,  adc@4,  cad@8,  cda@12,  dac@16,  dca@20",
        &needles,
    );
    test_case(
        &t,
        "XXX",
        "dac@0,  cad@0,  cda@0,  dca@0,  adc@0,  acd@0",
        &needles,
    );
    test_case(
        &t,
        "aXXX",
        "acd@0, adc@0, dac@1, cad@1, cda@1, dca@1, adc@1, acd@1",
        &needles,
    );
    test_case(&t, "aXBX", "acd@0, cda@1, adc@1", &needles);
    t.set_max_aaa_count(2);
    test_case(&t, "XXX", "", &needles);
    test_case(&t, "aXXX", "acd@0,  adc@0", &needles);
    test_case(&t, "XXXc", "adc@1,  dac@1", &needles);

    let needles = ["eq", "nd", "llll"];
    let t = trie(2, 0, &needles);
    test_case(
        &t,
        "aXXaBBkkZZlllllk",
        "nd@1, nd@4, eq@1, eq@8, llll@10, llll@11   ",
        &needles,
    );

    let needles = ["acd", "adc", "cad", "cda", "dac", "dca"];
    let t = trie(0, 1, &needles);
    test_case(&t, "aaaIIcccIIddd", "", &needles);

    let mut t = trie(0, 3, &needles);
    test_case(
        &t,
        "mmmm",
        "  dac@0,  cad@0,  cda@0,  dca@0,  adc@0,  acd@0, dac@1,  cad@1,  cda@1,  dca@1,  adc@1,  acd@1",
        &needles,
    );
    t.set_max_mm_count(2);
    test_case(&t, "aMMM", "acd@0,  adc@0", &needles);
    t.set_max_aaa_count(1);
    test_case(&t, "aMMB", "  adc@0,  acd@0, cad@1,  acd@1", &needles);

    let needles = ["eq", "nd", "llll"];
    let t = trie(1, 1, &needles);
    test_case(
        &t,
        "aXXaBBkkZZlllllk",
        "nd@0, nd@1, nd@2, nd@3, nd@4, nd@5, eq@0, eq@1, eq@2, eq@7, eq@8, eq@9, llll@9, llll@10, llll@11, llll@12    ",
        &needles,
    );

    let needles = ["acb", "abc", "cda", "bac", "anc", "acn", "dad"];
    let t = trie(0, 0, &needles);
    test_case(&t, "baxyacbIIabcIIbac", "acb@4, abc@9, bac@14", &needles);
    let t = trie(1, 0, &needles);
    test_case(
        &t,
        "baxyacbIIabcIIbac",
        "acb@4, abc@9, bac@0, bac@14, anc@9, acn@4",
        &needles,
    );

    let needles = ["dad", "bax", "bac", "anc", "acn"];
    let t = trie(2, 0, &needles);
    test_case(
        &t,
        "baxyacbIIabcIIbac",
        "dad@0, bax@0, bac@0, bac@14, anc@9, acn@4",
        &needles,
    );

    let needles = ["dady", "baxy", "iibac", "ancii", "yaknif"];
    let t = trie(2, 2, &needles);
    test_case(
        &t,
        "baxyacbIIabcIIbac",
        "dady@0, baxy@0, dady@8, ancii@4, ancii@9, iibac@1, iibac@7, iibac@12, yaknif@3",
        &needles,
    );

    let needles = ["PEPTIDER", "XXXBEBEAR"];
    let t = trie(0, 0, &needles);
    test_case(&t, "PEPTIDERXXXBEBEAR", "PEPTIDER@0, XXXBEBEAR@8", &needles);

    // Offsets in the presence of non-AA characters such as '*'.
    let needles = ["MLTEAEK"];
    let t = trie(0, 0, &needles);
    test_case(&t, "*MLTEAXK*", "", &needles);
    let t = trie(1, 0, &needles);
    test_case(&t, "*MLTEAXK*", "MLTEAEK@1", &needles);

    // A scout does not report hits that do not cover its first ambiguous residue.
    let needles = ["MDDDEADC", "MDD", "DD", "DEADC"];
    let t = trie(4, 0, &needles);
    test_case(
        &t,
        "MBBDEABCRAFG",
        "MDDDEADC@0, MDD@0, DD@1, DD@2, DEADC@3",
        &needles,
    );
}

#[test]
fn set_query() {
    let mut state = ACTrieState::default();
    let q = "PFAINGER";
    state.set_query(q);
    assert_eq!(state.get_query(), q);
    assert!(state.hits.is_empty());
    assert_eq!(state.tree_pos.pos(), 0);
    assert_eq!(state.text_pos_it()[0], q.as_bytes()[0]);
    assert!(state.scouts.is_empty());
    assert_eq!(state.text_pos(), 0);
}

fn check_aa_iterator(expected: &[AA], positions: &[usize], state: &mut ACTrieState) {
    let mut i = 0;
    loop {
        let aa = state.next_valid_aa();
        if !aa.is_valid() {
            break;
        }
        assert_eq!(aa, expected[i]);
        assert_eq!(state.text_pos(), positions[i]);
        i += 1;
    }
    assert_eq!(expected.len(), i);
}

#[test]
fn state_next_valid_aa() {
    let mut state = ACTrieState::default();
    state.set_query("PFBNX");
    let aas: Vec<AA> = b"PFBNX".iter().map(|&c| AA::from_byte(c)).collect();
    check_aa_iterator(&aas, &[1, 2, 3, 4, 5], &mut state);
    state.set_query("?X*B**");
    assert_eq!(state.get_query(), "?X*B**");
    check_aa_iterator(
        &[AA::from_byte(b'X'), AA::from_byte(b'B')],
        &[2, 4],
        &mut state,
    );
    state.set_query("");
    assert!(!state.next_valid_aa().is_valid());
}

#[test]
fn aa_from_index() {
    assert_eq!(AA::from_index(0).unwrap().to_char(), 'A');
    assert_eq!(
        AA::from_index(AA::unambiguous_aa_count() - 1)
            .unwrap()
            .to_char(),
        'V'
    );
    assert!(AA::from_index(AA::unambiguous_aa_count()).is_err());
}

#[test]
fn aa_properties() {
    // The class test's static_asserts and the letter loop.
    let aa = AA::from_byte;
    assert!(!aa(b'?').is_valid());
    assert_eq!(aa(b'?').code(), CHAR_TO_AA[b'?' as usize]);
    assert!(aa(b'G') <= aa(b'B'));
    for c in [b'B', b'J', b'Z', b'X', b'$'] {
        assert!(aa(c).is_ambiguous());
    }
    assert_eq!(aa(b'B').succ(), Some(aa(b'J')));
    assert_eq!((aa(b'B') - aa(b'B')).code(), 0);
    assert_eq!((aa(b'J') - aa(b'B')).code(), 1);
    assert_eq!((aa(b'Z') - aa(b'B')).code(), 2);
    assert_eq!((aa(b'X') - aa(b'B')).code(), 3);
    for c in b'A'..=b'Z' {
        assert!(aa(c).is_valid_for_peptide());
    }
    for c in [b'?', b'$', b' ', b'*', b'3', b'#'] {
        assert!(!aa(c).is_valid_for_peptide());
    }
    assert!(!aa(200).is_valid());
    assert!(!aa(255).is_valid());
    assert_eq!(AA::default(), aa(b'?'));
}

#[test]
fn tables_and_neighbours() {
    // native: the header's table layout and the neighbour ranges the
    // ambiguity expansion relies on (the .cpp's static_asserts).
    let aa = AA::from_byte;
    assert_eq!(aa(b'D').succ(), Some(aa(b'N')));
    assert_eq!(aa(b'E').succ(), Some(aa(b'Q')));
    assert_eq!(aa(b'I').succ(), Some(aa(b'L')));
    assert_eq!(aa(b'Z').succ(), Some(aa(b'X')));
    assert_eq!(aa(b'X').succ(), Some(aa(b'$')));
    assert_eq!(aa(b'$').succ(), Some(aa(b'?')));
    assert_eq!(aa(b'?').succ(), None);
    for (code, &symbol) in AA_TO_CHAR.iter().enumerate().take(26) {
        assert_eq!(CHAR_TO_AA[symbol as usize] as usize, code);
        assert_eq!(
            CHAR_TO_AA[symbol.to_ascii_lowercase() as usize] as usize,
            code
        );
    }
    assert_eq!(
        CHAR_TO_AA.iter().filter(|&&code| code == 27).count(),
        256 - 53
    );
    assert_eq!(aa(b'k').to_char(), 'K');
}

#[test]
fn small_types() {
    // native: Index, Bitset, Hit, ACNode, ACScout.
    assert!(Index::default().is_invalid());
    let mut index = Index::new(3);
    assert!(index.is_valid());
    *index.pos_mut() = 5;
    assert_eq!(index.pos(), 5);
    assert_eq!(Index::from(5u32), index);

    let mut bits = Bitset::default();
    bits.set(0);
    bits.set(3);
    bits.set(40); // ignored
    assert!(bits.test(3) && !bits.test(2) && !bits.test(40));
    assert_eq!(bits.pop_count(), 2);
    bits.clear(0);
    assert_eq!(bits.bits, 8);
    assert_eq!((bits << 1).bits, 16);
    assert_eq!((bits << 32).bits, 0);

    assert_eq!(
        Hit::new(1, 2, 3),
        Hit {
            needle_index: 1,
            needle_length: 2,
            query_pos: 3
        }
    );
    let node = ACNode::default();
    assert_eq!((node.suffix.pos(), node.first_child.pos()), (0, 0));
    assert_eq!(node.edge, AA::INVALID);
    let node = ACNode::new(AA::from_byte(b'K'), 4);
    assert_eq!(node.depth_and_hits.depth, 4);
    assert!(!node.depth_and_hits.has_hit);

    let mut scout = ACScout::new(1, Index::new(2), 1, 0, 3);
    assert_eq!(scout.text_pos(), 1);
    assert_eq!(scout.next_valid_aa("A*C"), AA::from_byte(b'C'));
    assert_eq!(scout.text_pos(), 3);
    let mut pos = 0;
    assert_eq!(next_valid_aa(b"*\0A", &mut pos), AA::INVALID); // NUL ends the query
}

#[test]
fn get_all_hits_matches_next_hits() {
    // native: getAllHits is the concatenation of the nextHits results.
    let needles = ["dady", "baxy", "iibac", "ancii", "yaknif"];
    let t = trie(2, 2, &needles);
    let mut state = ACTrieState::default();
    state.set_query("baxyacbIIabcIIbac");
    t.get_all_hits(&mut state).unwrap();
    let all = state.hits.clone();
    state.set_query("baxyacbIIabcIIbac");
    let mut concatenated = Vec::new();
    while t.next_hits(&mut state).unwrap() {
        concatenated.extend_from_slice(&state.hits);
    }
    // Trailing scout hits after the primary's last hit are kept by
    // getAllHits and dropped by a nextHits loop, as in the source.
    assert!(all.len() >= concatenated.len());
    assert_eq!(&all[..concatenated.len()], &concatenated[..]);
}

#[test]
fn refusals() {
    // native: the source hangs, wraps or reads out of bounds here.
    let mut t = ACTrie::new(0, 0);
    assert!(t.add_needle("").is_err());
    assert!(t.add_needle(&"A".repeat(MAX_NEEDLE_LENGTH + 1)).is_err());
    t.add_needle(&"A".repeat(MAX_NEEDLE_LENGTH)).unwrap();
    let mut state = ACTrieState::default();
    state.set_query("AAAA");
    assert!(t.next_hits(&mut state).is_err()); // not compressed
    t.compress_trie();
    let size = t.size();
    t.compress_trie(); // a second compression is a no-op
    assert_eq!(t.size(), size);
    assert!(t.add_needle("AA").is_err()); // no needles after compression
    assert!(!t.next_hits(&mut state).unwrap());
    let mut unset = ACTrieState::default();
    assert!(t.next_hits(&mut unset).is_err());
    // A refused needle leaves the trie unchanged.
    let mut t = ACTrie::new(0, 0);
    t.add_needle("PEP").unwrap();
    let before = t.size();
    assert!(t.add_needle("QQQ*").is_err());
    assert_eq!(t.size(), before);
    assert_eq!(t.get_needle_count(), 1);
}
