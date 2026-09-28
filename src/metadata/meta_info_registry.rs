// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Registry which assigns unique integer indices to metadata names.
//!
//! Ports `src/openms/include/OpenMS/METADATA/MetaInfoRegistry.h` and
//! `src/openms/source/METADATA/MetaInfoRegistry.cpp` (core SDK `bc9cc12`). See
//! `docs/META_INFO_REGISTRY_SUPPORT.md`.
//!
//! Registering a new name assigns an index of [`FIRST_REGISTERED_INDEX`](crate::metadata::meta_info_registry::FIRST_REGISTERED_INDEX)
//! (1024) or more, counting up. Indices 1 to 1023 are reserved for fast access
//! and never change; thirteen of them are taken by [`RESERVED_NAMES`](crate::metadata::meta_info_registry::RESERVED_NAMES):
//!
//! | index | name |
//! |---|---|
//! | 1 | `isotopic_range` |
//! | 2 | `cluster_id` |
//! | 3 | `label` |
//! | 4 | `icon` |
//! | 5 | `color` |
//! | 6 | `RT` |
//! | 7 | `MZ` |
//! | 8 | `predicted_RT` |
//! | 9 | `predicted_RT_p_value` |
//! | 10 | `spectrum_reference` |
//! | 11 | `ID` |
//! | 12 | `low_quality` |
//! | 13 | `charge` |
//!
//! The source keeps one process-wide registry behind `MetaInfo::registry()`
//! and serialises every member with `#pragma omp critical`. This crate owns no
//! mutable global state, so [`MetaInfoRegistry`](crate::metadata::meta_info_registry::MetaInfoRegistry) is an ordinary value the
//! caller owns; lookups take `&self` and registration takes `&mut self`, and a
//! registry shared between threads goes behind a `Mutex` or `RwLock`, which is
//! what the source's critical sections amount to.
//!
//! # Example
//!
//! ```
//! use openms::metadata::meta_info_registry::MetaInfoRegistry;
//!
//! let mut registry = MetaInfoRegistry::default();
//! let index = registry.register_name("retention time", "this is just another test", "sec").unwrap();
//! assert_eq!(index, 1024);
//! assert_eq!(registry.unit("retention time").unwrap(), "sec");
//! assert_eq!(registry.index_of("charge"), Some(13));
//! assert_eq!(registry.index_of("unregistered name"), None);
//! ```

use std::collections::BTreeMap;

use crate::{Error, Result};

/// The first index [`MetaInfoRegistry::register_name`](crate::metadata::meta_info_registry::MetaInfoRegistry::register_name) assigns; indices below
/// it are reserved.
pub const FIRST_REGISTERED_INDEX: u32 = 1024;

/// The source's "not registered" answer of `getIndex`, `UInt(-1)`.
/// [`MetaInfoRegistry::index_of`](crate::metadata::meta_info_registry::MetaInfoRegistry::index_of) returns `None` instead, and no name is ever
/// given this index.
pub const UNREGISTERED_INDEX: u32 = u32::MAX;

/// The thirteen names the source's constructor registers, as
/// `(index, name, description)`; every unit is empty.
pub const RESERVED_NAMES: [(u32, &str, &str); 13] = [
    (
        1,
        "isotopic_range",
        "consecutive numbering of the peaks in an isotope pattern. 0 is the monoisotopic peak",
    ),
    (
        2,
        "cluster_id",
        "consecutive numbering of isotope clusters in a spectrum",
    ),
    (3, "label", "label e.g. shown in visualization"),
    (4, "icon", "icon shown in visualization"),
    (
        5,
        "color",
        "color used for visualization e.g. #FF00FF for purple",
    ),
    (6, "RT", "the retention time of an identification"),
    (7, "MZ", "the MZ of an identification"),
    (
        8,
        "predicted_RT",
        "the predicted retention time of a peptide hit",
    ),
    (
        9,
        "predicted_RT_p_value",
        "the predicted RT p-value of a peptide hit",
    ),
    (
        10,
        "spectrum_reference",
        "Reference to a spectrum or feature number",
    ),
    (11, "ID", "Some type of identifier"),
    (
        12,
        "low_quality",
        "Flag which indicates that some entity has a low quality (e.g. a feature pair)",
    ),
    (13, "charge", "Charge of a feature or peak"),
];

/// Name, description and unit of one registered index.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Entry {
    name: String,
    description: String,
    unit: String,
}

/// Registry which assigns unique integer indices to strings
/// (`MetaInfoRegistry`), with a description and a unit per name.
///
/// [`Default`] is the source's constructor: the [`RESERVED_NAMES`](crate::metadata::meta_info_registry::RESERVED_NAMES) registered
/// and the next index [`FIRST_REGISTERED_INDEX`](crate::metadata::meta_info_registry::FIRST_REGISTERED_INDEX). `Clone` is the copy
/// constructor and assignment operator; equality, which the source does not
/// define, compares the whole state.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MetaInfoRegistry {
    next_index: u32,
    name_to_index: BTreeMap<String, u32>,
    entries: BTreeMap<u32, Entry>,
}

impl Default for MetaInfoRegistry {
    fn default() -> Self {
        let mut registry = Self {
            next_index: FIRST_REGISTERED_INDEX,
            name_to_index: BTreeMap::new(),
            entries: BTreeMap::new(),
        };
        for (index, name, description) in RESERVED_NAMES {
            registry.name_to_index.insert(name.to_string(), index);
            registry.entries.insert(
                index,
                Entry {
                    name: name.to_string(),
                    description: description.to_string(),
                    unit: String::new(),
                },
            );
        }
        registry
    }
}

impl MetaInfoRegistry {
    /// Register `name` with `description` and `unit`, and return its index
    /// (`registerName`; the source defaults both texts to empty).
    ///
    /// A name registered before keeps its index, and its stored description
    /// and unit are **not** changed — the two arguments are ignored then, as in
    /// the source. Use [`MetaInfoRegistry::set_description`](crate::metadata::meta_info_registry::MetaInfoRegistry::set_description) and
    /// [`MetaInfoRegistry::set_unit`](crate::metadata::meta_info_registry::MetaInfoRegistry::set_unit) to change them.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when every index up to
    /// [`UNREGISTERED_INDEX`] is taken. The source's unsigned counter would
    /// then wrap around and hand out the sentinel and reserved indices again;
    /// the registry is left unchanged instead.
    pub fn register_name(&mut self, name: &str, description: &str, unit: &str) -> Result<u32> {
        if let Some(index) = self.name_to_index.get(name) {
            return Ok(*index);
        }
        let index = self.next_index;
        if index == UNREGISTERED_INDEX {
            return Err(Error::InvalidRange(format!(
                "MetaInfoRegistry: no index left to register '{name}'"
            )));
        }
        self.name_to_index.insert(name.to_string(), index);
        self.entries.insert(
            index,
            Entry {
                name: name.to_string(),
                description: description.to_string(),
                unit: unit.to_string(),
            },
        );
        self.next_index = index + 1;
        Ok(index)
    }

    /// Set the description of a registered index (`setDescription(UInt, ...)`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an unregistered index ("Unregistered
    /// index!").
    pub fn set_description_at(&mut self, index: u32, description: &str) -> Result<()> {
        self.entry_mut(index)?.description = description.to_string();
        Ok(())
    }

    /// Set the description of a registered name (`setDescription(const
    /// std::string&, ...)`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an unregistered name ("Unregistered
    /// name!").
    pub fn set_description(&mut self, name: &str, description: &str) -> Result<()> {
        let index = self.registered(name, "Unregistered name!")?;
        self.entry_mut(index)?.description = description.to_string();
        Ok(())
    }

    /// Set the unit of a registered index (`setUnit(UInt, ...)`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an unregistered index ("Unregistered
    /// index!").
    pub fn set_unit_at(&mut self, index: u32, unit: &str) -> Result<()> {
        self.entry_mut(index)?.unit = unit.to_string();
        Ok(())
    }

    /// Set the unit of a registered name (`setUnit(const std::string&, ...)`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an unregistered name ("Unregistered
    /// name!").
    pub fn set_unit(&mut self, name: &str, unit: &str) -> Result<()> {
        let index = self.registered(name, "Unregistered name!")?;
        self.entry_mut(index)?.unit = unit.to_string();
        Ok(())
    }

    /// The index of `name`, or `None` when it is not registered (`getIndex`,
    /// which returns [`UNREGISTERED_INDEX`] there).
    pub fn index_of(&self, name: &str) -> Option<u32> {
        self.name_to_index.get(name).copied()
    }

    /// The name registered at `index` (`getName`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an unregistered index ("Unregistered
    /// index!").
    pub fn name_at(&self, index: u32) -> Result<&str> {
        Ok(&self.entry(index)?.name)
    }

    /// The description registered at `index` (`getDescription(UInt)`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an unregistered index ("Unregistered
    /// index!").
    pub fn description_at(&self, index: u32) -> Result<&str> {
        Ok(&self.entry(index)?.description)
    }

    /// The description of `name` (`getDescription(const std::string&)`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an unregistered name ("Unregistered
    /// Name!", capitalised as in the source).
    pub fn description(&self, name: &str) -> Result<&str> {
        let index = self.registered(name, "Unregistered Name!")?;
        Ok(&self.entry(index)?.description)
    }

    /// The unit registered at `index` (`getUnit(UInt)`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an unregistered index ("Unregistered
    /// index!").
    pub fn unit_at(&self, index: u32) -> Result<&str> {
        Ok(&self.entry(index)?.unit)
    }

    /// The unit of `name` (`getUnit(const std::string&)`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an unregistered name ("Unregistered
    /// Name!", capitalised as in the source).
    pub fn unit(&self, name: &str) -> Result<&str> {
        let index = self.registered(name, "Unregistered Name!")?;
        Ok(&self.entry(index)?.unit)
    }

    /// The number of registered names, reserved ones included. Native: the
    /// source exposes no count.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether nothing is registered, which a registry built by [`Default`]
    /// never is. Native, as [`MetaInfoRegistry::len`](crate::metadata::meta_info_registry::MetaInfoRegistry::len).
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    fn registered(&self, name: &str, message: &str) -> Result<u32> {
        self.index_of(name)
            .ok_or_else(|| Error::InvalidValue(format!("{message}: '{name}'")))
    }

    fn entry(&self, index: u32) -> Result<&Entry> {
        self.entries
            .get(&index)
            .ok_or_else(|| unregistered_index(index))
    }

    fn entry_mut(&mut self, index: u32) -> Result<&mut Entry> {
        self.entries
            .get_mut(&index)
            .ok_or_else(|| unregistered_index(index))
    }
}

fn unregistered_index(index: u32) -> Error {
    Error::InvalidValue(format!("Unregistered index!: '{index}'"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exhausting_the_index_space_leaves_the_registry_unchanged() {
        let mut registry = MetaInfoRegistry {
            next_index: UNREGISTERED_INDEX - 1,
            ..MetaInfoRegistry::default()
        };
        assert_eq!(
            registry.register_name("last", "", "").unwrap(),
            UNREGISTERED_INDEX - 1
        );
        let before = registry.clone();
        assert!(matches!(
            registry.register_name("one too many", "", ""),
            Err(Error::InvalidRange(_))
        ));
        assert_eq!(registry, before);
        // A name registered before still answers.
        assert_eq!(
            registry.register_name("last", "", "").unwrap(),
            UNREGISTERED_INDEX - 1
        );
    }
}
