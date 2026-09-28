// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Unique id to index lookup for containers of elements with unique ids:
//! `CONCEPT/UniqueIdIndexer.h`.
//!
//! The source is a CRTP base class that containers such as `FeatureMap` and
//! `ConsensusMap` derive from. It keeps a hash map from unique id to index,
//! refreshed on demand, so that repeated lookups are expected constant time
//! (linear upon failure). Rust has no CRTP mixin with a `mutable` member; the
//! port is [`UniqueIdIndexer`](crate::concept::unique_id_indexer::UniqueIdIndexer), the cached map as a value of its own, and each
//! operation takes the container's elements as a slice of [`HasUniqueId`]
//! values. A container keeps an indexer next to its elements and passes them
//! in, exactly as the source's `getBase_()` reaches the derived container's
//! data.
//!
//! # API mapping
//!
//! | Source (`UniqueIdIndexer.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `template <typename T> class UniqueIdIndexer` | [`UniqueIdIndexer`](crate::concept::unique_id_indexer::UniqueIdIndexer) |
//! | `typedef std::unordered_map<UInt64, Size> UniqueIdMap` | [`UniqueIdMap`](crate::concept::unique_id_indexer::UniqueIdMap) |
//! | `uniqueIdToIndex(UInt64) const` | [`UniqueIdIndexer::unique_id_to_index`](crate::concept::unique_id_indexer::UniqueIdIndexer::unique_id_to_index) |
//! | `updateUniqueIdToIndex() const` | [`UniqueIdIndexer::update_unique_id_to_index`](crate::concept::unique_id_indexer::UniqueIdIndexer::update_unique_id_to_index) |
//! | `resolveUniqueIdConflicts()` | [`UniqueIdIndexer::resolve_unique_id_conflicts`](crate::concept::unique_id_indexer::UniqueIdIndexer::resolve_unique_id_conflicts) |
//! | `swap(UniqueIdIndexer&)` | [`UniqueIdIndexer::swap`](crate::concept::unique_id_indexer::UniqueIdIndexer::swap) (also `std::mem::swap`) |
//! | protected `getBase_()` (both overloads) | the `data` slice argument of each operation |
//! | protected `mutable UniqueIdMap uniqueid_to_index_` | [`UniqueIdIndexer::unique_id_map`](crate::concept::unique_id_indexer::UniqueIdIndexer::unique_id_map), [`UniqueIdIndexer::unique_id_map_mut`](crate::concept::unique_id_indexer::UniqueIdIndexer::unique_id_map_mut) |
//!
//! # Native differences
//!
//! - The source's lookup is `const` and updates its `mutable` map; here it
//!   takes `&mut self`.
//! - `Size(-1)` for "not found" is `None`.
//! - The source's `Exception::Postcondition` for duplicate valid unique ids is
//!   [`crate::Error::InvalidValue`] with the source's message. As in the
//!   source, the map has been updated when the error is reported.
//! - `resolveUniqueIdConflicts` draws new ids from the source's process-wide
//!   `UniqueIdGenerator`; here the caller passes the generator, as everywhere
//!   in `concept::unique_id`.
//! - `FeatureMap::unique_id_to_index` and `ConsensusMap::unique_id_to_index`
//!   in `kernel` rebuild their lookup on every call and do not hold an
//!   indexer; this type is the cached form of the same rule.

use std::collections::HashMap;

use super::unique_id::{HasUniqueId, UniqueId, UniqueIdGenerator};
use crate::{Error, Result};

/// Hash map from unique id to container index: source `UniqueIdMap`.
pub type UniqueIdMap = HashMap<u64, usize>;

/// A cached map from the unique ids of a container's elements to their
/// indices: source `UniqueIdIndexer<T>`.
///
/// # Examples
///
/// ```
/// use openms::concept::unique_id_indexer::UniqueIdIndexer;
/// use openms::concept::UniqueId;
///
/// let elements = [UniqueId(1000), UniqueId(1010), UniqueId(1020)];
/// let mut indexer = UniqueIdIndexer::new();
/// assert_eq!(indexer.unique_id_to_index(&elements, 1010).ok(), Some(Some(1)));
/// assert_eq!(indexer.unique_id_to_index(&elements, 7).ok(), Some(None));
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct UniqueIdIndexer {
    /// Source `uniqueid_to_index_`.
    uniqueid_to_index: UniqueIdMap,
}

impl UniqueIdIndexer {
    /// An indexer with an empty map.
    pub fn new() -> Self {
        Self::default()
    }

    /// The index of the element of `data` with `unique_id`, or `None`: source
    /// `uniqueIdToIndex(unique_id)`.
    ///
    /// The lookup consults the map; if an index is found and that element
    /// indeed has this unique id, it is returned. Otherwise the map is
    /// updated with [`Self::update_unique_id_to_index`] and consulted again.
    /// Subordinate elements are not considered. `data` must be the container
    /// the map was built from; any other slice only costs a rebuild.
    ///
    /// # Errors
    ///
    /// As [`Self::update_unique_id_to_index`], when the rebuild finds
    /// duplicate valid unique ids.
    pub fn unique_id_to_index<E: HasUniqueId>(
        &mut self,
        data: &[E],
        unique_id: u64,
    ) -> Result<Option<usize>> {
        if let Some(&index) = self.uniqueid_to_index.get(&unique_id) {
            if data
                .get(index)
                .is_some_and(|element| element.unique_id() == unique_id)
            {
                return Ok(Some(index));
            }
        }
        self.update_unique_id_to_index(data)?;
        Ok(self.uniqueid_to_index.get(&unique_id).copied())
    }

    /// Updates the map from the unique ids of `data`: source
    /// `updateUniqueIdToIndex()`.
    ///
    /// Every valid unique id is mapped to its element's index; then the
    /// invalid id and every entry whose index is out of range or whose element
    /// no longer carries that id are removed.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`], the source's `Exception::Postcondition`, when
    /// two elements carry the same valid unique id, with the source's message
    /// `Duplicate valid unique ids detected!   RandomAccessContainer has
    /// size()==<n>, num_valid_unique_id==<v>, uniqueid_to_index_.size()==<m>`.
    /// The map keeps its updated state, as in the source.
    pub fn update_unique_id_to_index<E: HasUniqueId>(&mut self, data: &[E]) -> Result<()> {
        let mut num_valid_unique_id = 0usize;
        for (index, element) in data.iter().enumerate() {
            let unique_id = element.unique_id();
            if UniqueId::is_valid(unique_id) {
                self.uniqueid_to_index.insert(unique_id, index);
                num_valid_unique_id += 1;
            }
        }
        self.uniqueid_to_index.remove(&UniqueId::INVALID);
        self.uniqueid_to_index.retain(|&unique_id, &mut index| {
            data.get(index)
                .is_some_and(|element| element.unique_id() == unique_id)
        });
        if self.uniqueid_to_index.len() != num_valid_unique_id {
            return Err(Error::InvalidValue(format!(
                "Duplicate valid unique ids detected!   RandomAccessContainer has size()=={}, \
                 num_valid_unique_id=={num_valid_unique_id}, uniqueid_to_index_.size()=={}",
                data.len(),
                self.uniqueid_to_index.len()
            )));
        }
        Ok(())
    }

    /// Assigns new unique ids to elements whose id is invalid or occurs
    /// earlier in `data`: source `resolveUniqueIdConflicts()`.
    ///
    /// This usually occurs when merging old feature files, which have
    /// sequentially increasing unique ids. The map is rebuilt from scratch: an
    /// element with an invalid id first receives one (`ensureUniqueId`), and
    /// while its id is already taken it draws a new one from `generator`, each
    /// draw counting as one replacement. Subordinate elements are not checked
    /// and may remain non-unique; they are associated with their parent, which
    /// makes identification unique. Returns the number of replacements, not
    /// counting the ids given to elements that had none.
    pub fn resolve_unique_id_conflicts<E: HasUniqueId>(
        &mut self,
        data: &mut [E],
        generator: &mut UniqueIdGenerator,
    ) -> usize {
        let mut invalid_uids = 0;
        self.uniqueid_to_index.clear();
        for (index, element) in data.iter_mut().enumerate() {
            element.ensure_unique_id(generator);
            while self.uniqueid_to_index.contains_key(&element.unique_id()) {
                element.assign_new_unique_id(generator);
                invalid_uids += 1;
            }
            self.uniqueid_to_index.insert(element.unique_id(), index);
        }
        invalid_uids
    }

    /// Exchanges the maps of two indexers: source `swap(UniqueIdIndexer&)`.
    ///
    /// The source notes that only indexers for the same container type can be
    /// swapped; here the maps are plain values and any two can be.
    pub fn swap(&mut self, other: &mut Self) {
        std::mem::swap(&mut self.uniqueid_to_index, &mut other.uniqueid_to_index);
    }

    /// The cached map: the source's protected `uniqueid_to_index_`, which the
    /// upstream class test reaches through a derived accessor.
    pub fn unique_id_map(&self) -> &UniqueIdMap {
        &self.uniqueid_to_index
    }

    /// The cached map, mutably. Any edit is harmless: a wrong or missing entry
    /// only makes the next lookup rebuild the map.
    pub fn unique_id_map_mut(&mut self) -> &mut UniqueIdMap {
        &mut self.uniqueid_to_index
    }
}
