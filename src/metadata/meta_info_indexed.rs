// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The index-keyed half of `MetaInfo` and `MetaInfoInterface`, over a
//! caller-owned [`MetaInfoRegistry`](crate::metadata::meta_info_registry::MetaInfoRegistry).
//!
//! Completes the ports of `src/openms/include/OpenMS/METADATA/MetaInfo.h` and
//! `src/openms/include/OpenMS/METADATA/MetaInfoInterface.h` (core SDK
//! `bc9cc12`). See `docs/META_INFO_REGISTRY_SUPPORT.md`.
//!
//! The source stores meta values in a `std::map<UInt, DataValue>` keyed by the
//! index a process-wide `MetaInfoRegistry` assigns to each name: `setValue(name)`
//! registers the name, the `UInt` overloads address the map directly, and
//! `getKeys` lists entries in **index order**, which is registration order for
//! every name registered after the thirteen reserved ones. This crate stores
//! [`MetaInfo`] by name (`BTreeMap<String, MetaValue>`, name order) and owns no
//! global registry, so the index-keyed members are provided by
//! [`IndexedMetaInfo`](crate::metadata::meta_info_indexed::IndexedMetaInfo), implemented for `MetaInfo`, taking the registry as an
//! argument. `MetaInfoInterface` is composition — every record's `metadata`
//! field — so its `UInt` overloads (`setMetaValue`, `getMetaValue`,
//! `metaValueExists`, `removeMetaValue`, `getKeys`) are the same methods called
//! on that field, and `metaRegistry()` is the caller's registry, like
//! `MetaInfo::registry()`.
//!
//! # Example
//!
//! ```
//! use openms::metadata::meta_info_indexed::IndexedMetaInfo;
//! use openms::metadata::meta_info_registry::MetaInfoRegistry;
//! use openms::metadata::MetaInfo;
//!
//! let mut registry = MetaInfoRegistry::default();
//! let mut meta = MetaInfo::new();
//! meta.set_value_registered(&mut registry, "testname", "x".into())?;
//! meta.set_value_registered(&mut registry, "cluster_id", 4711.into())?;
//! // Name order in the map, registration (index) order through the registry.
//! assert_eq!(meta.keys().collect::<Vec<_>>(), ["cluster_id", "testname"]);
//! assert_eq!(meta.keys_in_index_order(&registry)?, ["cluster_id", "testname"]);
//! assert_eq!(meta.keys_as_indices(&registry)?, [2, 1024]);
//! assert!(meta.exists_at(&registry, 2));
//! # Ok::<(), openms::Error>(())
//! ```

use super::meta_info_registry::MetaInfoRegistry;
use super::{MetaInfo, MetaValue};
use crate::{Error, Result};

/// The `UInt`-index overloads of `MetaInfo` and `MetaInfoInterface`, and the
/// registering name setter, over a caller-owned registry.
///
/// An index means the name the registry holds for it; an index the registry
/// does not hold addresses nothing.
pub trait IndexedMetaInfo {
    /// Set `name` to `value`, registering `name` first when it is new, and
    /// return its index: the source's `setValue(const std::string&, ...)` /
    /// `setMetaValue(name, ...)`, which always registers the name in the
    /// global registry.
    ///
    /// Registering at set time is what gives [`keys_in_index_order`] the
    /// source's insertion-driven order. Setting through the map's own `insert`
    /// stores the value without registering it.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when the registry has no index left
    /// ([`MetaInfoRegistry::register_name`](crate::metadata::meta_info_registry::MetaInfoRegistry::register_name)); nothing is stored then.
    ///
    /// [`keys_in_index_order`]: IndexedMetaInfo::keys_in_index_order
    fn set_value_registered(
        &mut self,
        registry: &mut MetaInfoRegistry,
        name: &str,
        value: MetaValue,
    ) -> Result<u32>;

    /// Set the value of the name registered at `index`: the source's
    /// `setValue(UInt, ...)` / `setMetaValue(UInt, ...)`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `index` is not registered. The source stores
    /// a value under any index, and only a later `getKeys` of names throws on
    /// it; this crate stores by name and so needs the name now. Nothing is
    /// stored then.
    fn set_value_at(
        &mut self,
        registry: &MetaInfoRegistry,
        index: u32,
        value: MetaValue,
    ) -> Result<()>;

    /// The value stored for the name registered at `index`, or `None`: the
    /// source's `getValue(UInt, default_value)` / `getMetaValue(UInt, ...)`,
    /// whose `default_value` is `.unwrap_or(&default)` here.
    fn value_at<'a>(&'a self, registry: &MetaInfoRegistry, index: u32) -> Option<&'a MetaValue>;

    /// Whether a value is stored for the name registered at `index`: the
    /// source's `exists(UInt)` / `metaValueExists(UInt)`.
    fn exists_at(&self, registry: &MetaInfoRegistry, index: u32) -> bool;

    /// Remove the value stored for the name registered at `index`, if any: the
    /// source's `removeValue(UInt)` / `removeMetaValue(UInt)`. Removing an
    /// absent or unregistered index does nothing, as in the source.
    fn remove_value_at(&mut self, registry: &MetaInfoRegistry, index: u32);

    /// The indices of the stored names, ascending: the source's
    /// `getKeys(std::vector<UInt>&)`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a stored name is not registered — possible
    /// here because the map can be filled without the registry, impossible in
    /// the source, whose name setter always registers.
    fn keys_as_indices(&self, registry: &MetaInfoRegistry) -> Result<Vec<u32>>;

    /// The stored names in ascending index order: the source's
    /// `getKeys(std::vector<std::string>&)`, whose order is the index order of
    /// its map. The map's own `keys()` gives name order instead.
    ///
    /// # Errors
    ///
    /// As [`IndexedMetaInfo::keys_as_indices`](crate::metadata::meta_info_indexed::IndexedMetaInfo::keys_as_indices), for a stored name that is not
    /// registered.
    fn keys_in_index_order(&self, registry: &MetaInfoRegistry) -> Result<Vec<String>>;
}

impl IndexedMetaInfo for MetaInfo {
    fn set_value_registered(
        &mut self,
        registry: &mut MetaInfoRegistry,
        name: &str,
        value: MetaValue,
    ) -> Result<u32> {
        let index = registry.register_name(name, "", "")?;
        self.insert(name.to_string(), value);
        Ok(index)
    }

    fn set_value_at(
        &mut self,
        registry: &MetaInfoRegistry,
        index: u32,
        value: MetaValue,
    ) -> Result<()> {
        let name = registry.name_at(index)?;
        self.insert(name.to_string(), value);
        Ok(())
    }

    fn value_at<'a>(&'a self, registry: &MetaInfoRegistry, index: u32) -> Option<&'a MetaValue> {
        registry.name_at(index).ok().and_then(|name| self.get(name))
    }

    fn exists_at(&self, registry: &MetaInfoRegistry, index: u32) -> bool {
        self.value_at(registry, index).is_some()
    }

    fn remove_value_at(&mut self, registry: &MetaInfoRegistry, index: u32) {
        if let Ok(name) = registry.name_at(index) {
            self.remove(name);
        }
    }

    fn keys_as_indices(&self, registry: &MetaInfoRegistry) -> Result<Vec<u32>> {
        let mut indices = Vec::with_capacity(self.len());
        for name in self.keys() {
            let index = registry.index_of(name).ok_or_else(|| {
                Error::InvalidValue(format!(
                    "meta value name '{name}' is not registered in the MetaInfoRegistry"
                ))
            })?;
            indices.push(index);
        }
        indices.sort_unstable();
        Ok(indices)
    }

    fn keys_in_index_order(&self, registry: &MetaInfoRegistry) -> Result<Vec<String>> {
        self.keys_as_indices(registry)?
            .into_iter()
            .map(|index| registry.name_at(index).map(str::to_string))
            .collect()
    }
}
