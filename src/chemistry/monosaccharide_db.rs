// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Built-in ProForma monosaccharide symbols, literal masses and synonyms.
//!
//! The immutable registry embeds the pinned source dataset. It does not discover
//! runtime JSON files or use an installation/environment override.

use crate::{Error, Result};
use std::{collections::BTreeMap, sync::OnceLock};

/// Stored source fields. The formula is text; mass is not recomputed from it.
#[derive(Clone, Debug, PartialEq)]
pub struct Monosaccharide {
    pub symbol: String,
    pub name: String,
    pub mass: f64,
    pub formula: String,
    pub synonyms: Vec<String>,
}

/// Immutable built-in database with source primary/synonym lookup precedence.
#[derive(Debug)]
pub struct MonosaccharideDB {
    records: BTreeMap<String, Monosaccharide>,
    aliases: BTreeMap<String, String>,
}
impl MonosaccharideDB {
    /// One thread-safe shared instance, available without optional features.
    pub fn global() -> &'static Self {
        static DB: OnceLock<MonosaccharideDB> = OnceLock::new();
        DB.get_or_init(|| Self::from_embedded(EMBEDDED))
    }
    /// Exact, case-sensitive lookup. Whitespace and Unicode are not normalized.
    pub fn has_symbol(&self, symbol: &str) -> bool {
        self.aliases.contains_key(symbol)
    }
    /// Return the same borrowed record for a primary symbol and its aliases.
    ///
    /// Like source, all names route through the alias map; an alias inserted by
    /// a later primary can shadow an earlier primary's spelling.
    pub fn get(&self, symbol: &str) -> Option<&Monosaccharide> {
        self.aliases
            .get(symbol)
            .and_then(|primary| self.records.get(primary))
    }
    /// Checked counterpart to `get`; unknown input is not copied into the error.
    pub fn get_or_error(&self, symbol: &str) -> Result<&Monosaccharide> {
        self.get(symbol)
            .ok_or_else(|| Error::InvalidValue("unknown monosaccharide symbol".into()))
    }
    /// Sorted primary symbols only, borrowed from the database.
    pub fn all_symbols(&self) -> Vec<&str> {
        self.records.keys().map(String::as_str).collect()
    }
    pub fn len(&self) -> usize {
        self.records.len()
    }
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
    fn from_embedded(rows: &[EmbeddedMonosaccharide]) -> Self {
        // nlohmann::json uses lexical object-key order. Keep that order even if
        // a future regenerated source resource has a different textual order.
        let ordered: BTreeMap<_, _> = rows.iter().map(|row| (row.symbol, row)).collect();
        let mut result = Self {
            records: BTreeMap::new(),
            aliases: BTreeMap::new(),
        };
        for (symbol, row) in ordered {
            let record = Monosaccharide {
                symbol: symbol.into(),
                name: row.name.into(),
                mass: row.mass,
                formula: row.formula.into(),
                synonyms: row.synonyms.iter().map(|s| (*s).into()).collect(),
            };
            for synonym in &record.synonyms {
                result.aliases.insert(synonym.clone(), symbol.into());
            }
            result.records.insert(symbol.into(), record);
            // Each primary self-mapping follows that entry's synonyms.
            result.aliases.insert(symbol.into(), symbol.into());
        }
        result
    }
}

struct EmbeddedMonosaccharide {
    symbol: &'static str,
    name: &'static str,
    mass: f64,
    formula: &'static str,
    synonyms: &'static [&'static str],
}
const EMBEDDED: &[EmbeddedMonosaccharide] =
    include!("../../resources/monosaccharides/monosaccharides.rs");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lexical_source_alias_precedence_includes_shadowed_primaries() {
        let db = MonosaccharideDB::from_embedded(&[
            EmbeddedMonosaccharide {
                symbol: "B",
                name: "second",
                mass: -0.0,
                formula: "opaque B",
                synonyms: &["A", "shared", "shared"],
            },
            EmbeddedMonosaccharide {
                symbol: "A",
                name: "first",
                mass: -1.0,
                formula: "unparsed A",
                synonyms: &["B", "shared"],
            },
        ]);
        assert_eq!(db.all_symbols(), ["A", "B"]);
        for name in ["A", "B", "shared"] {
            assert_eq!(db.get(name).unwrap().symbol, "B");
        }
        assert_eq!(db.get("B").unwrap().mass.to_bits(), (-0f64).to_bits());
        assert_eq!(db.get("B").unwrap().formula, "opaque B");
        assert_eq!(db.get("B").unwrap().synonyms, ["A", "shared", "shared"]);
        assert_eq!(db.records["A"].mass, -1.0);
        assert_eq!(db.records["A"].formula, "unparsed A");
        assert!(!MonosaccharideDB::from_embedded(&[]).has_symbol("A"));
        assert!(MonosaccharideDB::from_embedded(&[]).is_empty());
    }
}
