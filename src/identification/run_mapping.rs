// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::{PeptideIdentification, ProteinIdentification, bad};
use crate::concept::constants::user_param::{BASE_NAME, ID_MERGE_INDEX};
use crate::data_structures::list::ListFormat;
use crate::{Error, Result};
use std::borrow::Cow;
use std::collections::BTreeMap;

/// Two-way identification-run/path mapping, including source duplicate failure
/// state and legacy peptide-level path fallback. Only identifiers and primary
/// run paths are consumed; protein hits, groups and other metadata are untouched.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct IdentifierMSRunMapper {
    forward: BTreeMap<String, Vec<String>>,
    reverse: BTreeMap<Vec<String>, String>,
}

impl IdentifierMSRunMapper {
    /// Maximum input run and path descriptors per create call.
    pub const MAX_ITEMS: usize = 1_000_000;
    /// Conservative cumulative owned payload/node allocation per create call.
    pub const MAX_BYTES: usize = 64 * 1024 * 1024;

    pub fn new() -> Self {
        Self::default()
    }
    pub fn from_runs(runs: &[ProteinIdentification]) -> Result<Self> {
        let mut mapping = Self::default();
        mapping.create(runs)?;
        Ok(mapping)
    }

    /// Replace the mapping. A duplicate path-list error deliberately publishes
    /// the complete forward map and only the reverse prefix preceding the first
    /// duplicate, as in C++. This keeps forward source-file resolution usable
    /// after ambiguity. Repeated identifiers overwrite forward paths; repeated
    /// path lists error even when their identifiers are equal.
    ///
    /// Resource errors occur before copying and leave the old mapping unchanged.
    /// Paths are exact, ordered strings: no sorting or filesystem normalization.
    pub fn create(&mut self, runs: &[ProteinIdentification]) -> Result<()> {
        preflight(runs)?;
        let mut next = Self::default();
        for run in runs {
            next.forward
                .insert(run.identifier.clone(), run.primary_ms_run_paths.clone());
        }
        for run in runs {
            if next.reverse.contains_key(&run.primary_ms_run_paths) {
                *self = next;
                return Err(bad(
                    "multiple protein identifications have the same MS run path list",
                ));
            }
            next.reverse
                .insert(run.primary_ms_run_paths.clone(), run.identifier.clone());
        }
        *self = next;
        Ok(())
    }
    pub fn is_empty(&self) -> bool {
        self.forward.is_empty()
    }
    pub fn len(&self) -> usize {
        self.forward.len()
    }
    pub fn has_identifier(&self, identifier: &str) -> bool {
        self.forward.contains_key(identifier)
    }
    /// Empty slice for both absent runs and explicitly empty path lists.
    pub fn ms_run_paths(&self, identifier: &str) -> &[String] {
        self.forward
            .get(identifier)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }
    /// Borrowed lexical traversal replaces the source's copied identifier vector.
    pub fn identifiers(&self) -> impl ExactSizeIterator<Item = &str> {
        self.forward.keys().map(String::as_str)
    }
    pub fn has_run_path(&self, paths: &[String]) -> bool {
        self.reverse.contains_key(paths)
    }
    /// Native Option replaces the source boolean plus mutable string output.
    pub fn try_identifier(&self, paths: &[String]) -> Option<&str> {
        self.reverse.get(paths).map(String::as_str)
    }
    pub fn identifier(&self, paths: &[String]) -> Result<&str> {
        self.try_identifier(paths)
            .ok_or_else(|| bad("MS run paths not found in mapping"))
    }

    /// Select a mapped file using integer id_merge_index, defaulting to zero.
    /// Out-of-range or missing mapping falls back to lenient base_name metadata.
    /// Negative or non-integer indices error when a nonempty mapped run exists,
    /// before that fallback, matching the source unsigned conversion.
    ///
    /// Returned text borrows stored paths/strings where possible. Only a numeric
    /// or list base_name needs bounded source-style formatting and owned text.
    pub fn primary_ms_run_path<'a>(
        &'a self,
        peptide: &'a PeptideIdentification,
    ) -> Result<Cow<'a, str>> {
        let paths = self.ms_run_paths(&peptide.identifier);
        if !paths.is_empty() {
            let index = match peptide.metadata.get(ID_MERGE_INDEX) {
                Some(value) => usize::try_from(value.as_i64()?)
                    .map_err(|_| bad("MS run merge index must be a nonnegative native index"))?,
                None => 0,
            };
            if let Some(path) = paths.get(index) {
                return Ok(Cow::Borrowed(path));
            }
        }
        match peptide.metadata.get(BASE_NAME) {
            Some(value) => value.to_list_text(),
            None => Ok(Cow::Borrowed("")),
        }
    }

    /// Require an unambiguous integer index only for runs with multiple files.
    /// Zero/one-path and unknown runs are exempt, including stale/wrong-type
    /// indices. This is an explicit check, not implicit in path lookup.
    pub fn validate_merge_index(
        &self,
        peptide: &PeptideIdentification,
        psm_index: usize,
    ) -> Result<()> {
        let paths = self.ms_run_paths(&peptide.identifier);
        if paths.len() < 2 {
            return Ok(());
        }
        let diagnostic = |reason: &str| {
            bad(&format!(
                "PSM #{psm_index} in a {}-file run: {reason}",
                paths.len()
            ))
        };
        let value = peptide
            .metadata
            .get(ID_MERGE_INDEX)
            .ok_or_else(|| diagnostic("missing id_merge_index"))?;
        let index = value
            .as_i64()
            .map_err(|_| diagnostic("id_merge_index is not an integer"))?;
        if usize::try_from(index)
            .ok()
            .is_none_or(|index| index >= paths.len())
        {
            return Err(diagnostic("id_merge_index is out of range"));
        }
        Ok(())
    }
}

fn limit() -> Error {
    bad("MS run mapping exceeds its descriptor or byte limit")
}
fn preflight(runs: &[ProteinIdentification]) -> Result<()> {
    let mut count = runs.len();
    if count > IdentifierMSRunMapper::MAX_ITEMS {
        return Err(limit());
    }
    // Two maps, each charged up to two full sparse BTreeMap nodes per input
    // record, including forward overwrites and prefix work before ambiguity.
    let per_run = std::mem::size_of::<(String, Vec<String>)>()
        .checked_mul(48)
        .and_then(|n| n.checked_add(512))
        .ok_or_else(limit)?;
    let mut bytes = runs.len().checked_mul(per_run).ok_or_else(limit)?;
    if bytes > IdentifierMSRunMapper::MAX_BYTES {
        return Err(limit());
    }
    for run in runs {
        count = count
            .checked_add(run.primary_ms_run_paths.len())
            .ok_or_else(limit)?;
        if count > IdentifierMSRunMapper::MAX_ITEMS {
            return Err(limit());
        }
        let descriptors = run
            .primary_ms_run_paths
            .len()
            .checked_mul(2 * std::mem::size_of::<String>())
            .ok_or_else(limit)?;
        bytes = bytes.checked_add(descriptors).ok_or_else(limit)?;
        if bytes > IdentifierMSRunMapper::MAX_BYTES {
            return Err(limit());
        }
        for text in std::iter::once(&run.identifier).chain(&run.primary_ms_run_paths) {
            bytes = text
                .len()
                .checked_mul(2)
                .and_then(|n| bytes.checked_add(n))
                .ok_or_else(limit)?;
            if bytes > IdentifierMSRunMapper::MAX_BYTES {
                return Err(limit());
            }
        }
    }
    Ok(())
}
