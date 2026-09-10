// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Cross-link modification lookup using the pinned XLMOD vocabulary.
//! The shared OBO loader preserves the source's union of reactive sites; these
//! records do not describe which pairs of sites may chemically react.

use super::{ModificationsDB, OboReadOptions};
use crate::Result;
use std::io::BufRead;
use std::sync::OnceLock;

/// A modification database loaded with the source CrossLinksDB selection rules.
/// The shared global database is immutable. Caller-owned databases support the
/// same checked additions as ModificationsDB through `database_mut`.
#[derive(Debug)]
pub struct CrossLinksDB {
    database: ModificationsDB,
}

impl CrossLinksDB {
    /// Load the bundled XLMOD vocabulary once, without filesystem/network access.
    /// The data's license and provenance are in resources/modifications.
    pub fn global() -> &'static Self {
        static DATABASE: OnceLock<CrossLinksDB> = OnceLock::new();
        DATABASE.get_or_init(|| {
            Self::from_obo(
                include_str!("../../resources/modifications/XLMOD.obo").as_bytes(),
                &OboReadOptions::default(),
            )
            .expect("bundled XLMOD cross-link records must pass the checked parser")
        })
    }

    /// Bounded caller-owned load. Cross-link selection is always enabled;
    /// `options.cross_links_only` is overridden without mutating the caller's
    /// options. Other parser limits retain their supplied values.
    pub fn from_obo(reader: impl BufRead, options: &OboReadOptions) -> Result<Self> {
        let options = OboReadOptions {
            cross_links_only: true,
            ..options.clone()
        };
        Ok(Self {
            database: ModificationsDB::from_obo(reader, &options)?,
        })
    }

    /// Reuse the shared lookup, accession, index, mass-search, and owned-handle API.
    pub fn database(&self) -> &ModificationsDB {
        &self.database
    }

    /// Explicitly extend a caller-owned database using the shared checked API.
    /// As with source addModification, custom additions need not come from XLMOD.
    pub fn database_mut(&mut self) -> &mut ModificationsDB {
        &mut self.database
    }

    /// Full IDs carrying an OBO accession, sorted ascending. This follows source
    /// getAllSearchModifications's PSI-MOD field, which also stores XLMOD IDs;
    /// it is not the general database's UniMod-only search list. Duplicates, if
    /// present in the underlying records, are retained as in the source vector.
    pub fn all_search_modifications(&self) -> Vec<String> {
        let mut names: Vec<_> = self
            .database
            .entries()
            .iter()
            .filter(|record| record.obo_accession().is_some_and(|id| !id.is_empty()))
            .map(|record| record.full_id().to_owned())
            .collect();
        names.sort();
        names
    }
}
