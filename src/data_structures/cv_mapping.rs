// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Owned CV mapping records and ordered source container operations.

use std::collections::BTreeSet;

/// Reference declared by a CV mapping file.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CVReference {
    pub name: String,
    pub identifier: String,
}

/// Term selection within a mapping rule. All source value defaults are false.
/// The XML loader separately defaults an omitted `isRepeatable` to true.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CVMappingTerm {
    pub accession: String,
    pub use_term_name: bool,
    pub use_term: bool,
    pub term_name: String,
    pub is_repeatable: bool,
    pub allow_children: bool,
    pub cv_identifier_ref: String,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[repr(u8)]
pub enum RequirementLevel {
    #[default]
    Must = 0,
    Should = 1,
    May = 2,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[repr(u8)]
pub enum CombinationsLogic {
    #[default]
    Or = 0,
    And = 1,
    Xor = 2,
}

/// Scalar fields and the ordered term vector are independently replaceable.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CVMappingRule {
    pub identifier: String,
    pub element_path: String,
    pub requirement_level: RequirementLevel,
    pub scope_path: String,
    pub combinations_logic: CombinationsLogic,
    pub terms: Vec<CVMappingTerm>,
}

/// Rules are replaced directly; reference bulk-set preserves source append behavior.
///
/// The source private map's last value is determined by the ordered references.
/// Its only public query tests identifier membership, so an index of identifiers
/// preserves the complete observable API without another owned copy of every name.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CVMappings {
    pub mapping_rules: Vec<CVMappingRule>,
    references: Vec<CVReference>,
    identifiers: BTreeSet<String>,
}

impl CVMappings {
    pub fn cv_references(&self) -> &[CVReference] {
        &self.references
    }

    pub fn has_cv_reference(&self, identifier: &str) -> bool {
        self.identifiers.contains(identifier)
    }

    /// Append all records, including duplicate identifiers, as source setCVReferences does.
    /// An empty argument leaves existing references unchanged.
    pub fn set_cv_references(&mut self, references: Vec<CVReference>) {
        for reference in &references {
            self.identifiers.insert(reference.identifier.clone());
        }
        self.references.extend(references);
    }

    /// Add a previously unseen identifier. A duplicate is ignored and returns false.
    /// The result replaces the source's process-global warning output.
    pub fn add_cv_reference(&mut self, reference: CVReference) -> bool {
        if !self.identifiers.insert(reference.identifier.clone()) {
            return false;
        }
        self.references.push(reference);
        true
    }

    /// Explicit native replacement of both the ordered references and their index.
    pub fn replace_cv_references(&mut self, references: Vec<CVReference>) {
        let mut next = Self::default();
        next.set_cv_references(references);
        self.references = next.references;
        self.identifiers = next.identifiers;
    }
}
