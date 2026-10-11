// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `UniqueIdIndexer_test.cpp` (core `bc9cc12`), one test per
//! `START_SECTION`, literals transcribed (tier 3).
//!
//! The source test shuffles its container with `Math::RandomShuffler{0}`; its
//! assertions hold for any order, so a fixed permutation stands in for it. The
//! `Postcondition` message it checks depends only on the counts, not on the
//! order.

use openms::concept::unique_id_indexer::{UniqueIdIndexer, UniqueIdMap};
use openms::concept::{HasUniqueId, UniqueIdGenerator};

/// The source test's `Dummy`: a `UniqueIdInterface` with a payload.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Dummy {
    dummy: usize,
    unique_id: u64,
}

impl HasUniqueId for Dummy {
    fn unique_id(&self) -> u64 {
        self.unique_id
    }
    fn unique_id_mut(&mut self) -> &mut u64 {
        &mut self.unique_id
    }
}

/// The source test's `DummyVectorIndexed`: elements plus their indexer.
#[derive(Clone, Debug, Default)]
struct DummyVectorIndexed {
    data: Vec<Dummy>,
    indexer: UniqueIdIndexer,
}

impl DummyVectorIndexed {
    fn with_ids(count: usize) -> Self {
        let mut dvi = Self::default();
        for i in 0..count {
            let mut dummy = Dummy {
                dummy: i,
                ..Dummy::default()
            };
            dummy.set_unique_id(10 * i as u64 + 1000);
            dvi.data.push(dummy);
        }
        dvi
    }

    fn unique_id_to_index(&mut self, unique_id: u64) -> Option<usize> {
        self.indexer
            .unique_id_to_index(&self.data, unique_id)
            .ok()
            .flatten()
    }
}

/// A fixed permutation in place of the source's seeded shuffle.
fn permute(data: &mut [Dummy]) {
    data.reverse();
    data.rotate_left(3);
}

#[test]
fn constructor_starts_with_an_empty_map() {
    let indexer = UniqueIdIndexer::new();
    assert!(indexer.unique_id_map().is_empty());
    assert_eq!(indexer, UniqueIdIndexer::default());
}

#[test]
fn unique_id_to_index() {
    let num_uii = 10;
    let mut dvi = DummyVectorIndexed::with_ids(num_uii);
    for i in 0..num_uii {
        assert_eq!(dvi.unique_id_to_index(10 * i as u64 + 1000), Some(i));
    }

    // "shuffling ..."
    permute(&mut dvi.data);
    for i in 0..num_uii {
        let current = dvi.data[i];
        assert_eq!(dvi.unique_id_to_index(current.unique_id()), Some(i));
    }

    dvi.data.pop();
    dvi.data.pop();
    let mut first = Dummy::default();
    first.set_unique_id(12_345_678);
    dvi.data.push(first);
    dvi.data.push(Dummy::default());
    let mut duplicate = Dummy::default();
    duplicate.set_unique_id(12_345_678);
    dvi.data.push(duplicate);
    dvi.data.push(Dummy::default());

    // "shuffling ..."
    permute(&mut dvi.data);
    let error = dvi.indexer.update_unique_id_to_index(&dvi.data);
    assert_eq!(
        error
            .err()
            .map(|error| error.to_string())
            .unwrap_or_default(),
        openms::Error::InvalidValue(
            "Duplicate valid unique ids detected!   RandomAccessContainer has size()==12, \
             num_valid_unique_id==10, uniqueid_to_index_.size()==9"
                .into()
        )
        .to_string()
    );
    // The lookup reports the same error: the source's uniqueIdToIndex only
    // catches std::out_of_range, not the Postcondition.
    assert!(dvi.indexer.unique_id_to_index(&dvi.data, 1).is_err());
}

#[test]
fn resolve_unique_id_conflicts() {
    let mut generator = UniqueIdGenerator::from_seed(0);
    let mut dvi = DummyVectorIndexed::with_ids(10);
    assert_eq!(
        dvi.indexer
            .resolve_unique_id_conflicts(&mut dvi.data, &mut generator),
        0
    );

    // introduce three doubles
    let mut a = Dummy::default();
    a.set_unique_id(1000);
    let mut b = Dummy::default();
    b.set_unique_id(1000 + 30);
    dvi.data.push(a);
    dvi.data.push(b);
    dvi.data.push(b);
    assert!(dvi.indexer.update_unique_id_to_index(&dvi.data).is_err());
    assert_eq!(
        dvi.indexer
            .resolve_unique_id_conflicts(&mut dvi.data, &mut generator),
        3
    );
    assert!(dvi.indexer.update_unique_id_to_index(&dvi.data).is_ok());
    assert_eq!(dvi.indexer.unique_id_map().len(), 13);
    // The first occurrences keep their ids; the later duplicates changed.
    assert_eq!(dvi.data[0].unique_id(), 1000);
    assert_eq!(dvi.data[3].unique_id(), 1030);
    assert_ne!(dvi.data[10].unique_id(), 1000);
    assert_ne!(dvi.data[11].unique_id(), 1030);
    assert_ne!(dvi.data[12].unique_id(), 1030);
}

#[test]
fn resolve_gives_invalid_ids_a_new_id_without_counting_it() {
    let mut generator = UniqueIdGenerator::from_seed(0);
    let mut data = vec![Dummy::default(), Dummy::default()];
    let mut indexer = UniqueIdIndexer::new();
    assert_eq!(
        indexer.resolve_unique_id_conflicts(&mut data, &mut generator),
        0
    );
    assert!(data.iter().all(HasUniqueId::has_valid_unique_id));
    assert_ne!(data[0].unique_id(), data[1].unique_id());
}

#[test]
fn swap() {
    let num_uii = 10;
    let mut dvi = DummyVectorIndexed::with_ids(num_uii);
    assert!(dvi.indexer.update_unique_id_to_index(&dvi.data).is_ok());

    let mut dvi2 = DummyVectorIndexed::default();
    assert_eq!(dvi.indexer.unique_id_map().len(), num_uii);
    assert_eq!(dvi2.indexer.unique_id_map().len(), 0);

    std::mem::swap(&mut dvi, &mut dvi2);
    assert_eq!(dvi.indexer.unique_id_map().len(), 0);
    assert_eq!(dvi2.indexer.unique_id_map().len(), num_uii);

    dvi = dvi2.clone();
    assert_eq!(dvi.indexer.unique_id_map().len(), num_uii);
    dvi.indexer.unique_id_map_mut().clear();
    assert_eq!(dvi.indexer.unique_id_map().len(), 0);

    assert_eq!(dvi.unique_id_to_index(4_321_234_324_124), None);
    assert_eq!(dvi.indexer.unique_id_map().len(), num_uii);

    // The member swap exchanges only the maps.
    let mut empty = UniqueIdIndexer::new();
    dvi.indexer.swap(&mut empty);
    assert_eq!(dvi.indexer.unique_id_map().len(), 0);
    assert_eq!(empty.unique_id_map().len(), num_uii);
}

#[test]
fn stale_entries_are_repaired_on_lookup() {
    let mut dvi = DummyVectorIndexed::with_ids(3);
    let map: &mut UniqueIdMap = dvi.indexer.unique_id_map_mut();
    // A wrong index and an out-of-range index only force a rebuild.
    map.insert(1000, 2);
    map.insert(1010, 99);
    assert_eq!(dvi.unique_id_to_index(1000), Some(0));
    assert_eq!(dvi.unique_id_to_index(1010), Some(1));
    // Invalid ids are never indexed.
    dvi.data.push(Dummy::default());
    assert_eq!(dvi.unique_id_to_index(0), None);
    assert!(!dvi.indexer.unique_id_map().contains_key(&0));
}
