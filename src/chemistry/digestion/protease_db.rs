// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Database of the enzymes that digest proteins, ported from
//! `CHEMISTRY/ProteaseDB.h` and the lookup half of its template base
//! `CHEMISTRY/DigestionEnzymeDB.h` at core `bc9cc12`.
//!
//! The source singleton is filled by two providers: the compiled-in
//! `BuiltInProteaseDataProvider` and, when `share/OpenMS/CHEMISTRY/Enzymes.xml`
//! is found, an `EnzymeXMLDataProvider` whose definitions replace the built-in
//! ones by name. At the pin both providers carry the same 33 definitions in the
//! same order (the XML's SHA-256 is recorded in `enzymes.rs`), so this port
//! replays both registration passes over one compiled table to reproduce the
//! source's lookup indices exactly. See `docs/DIGESTION_SUPPORT.md`.

use super::{DigestionEnzymeProtein, enzymes, invalid};
use crate::Result;
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;
use std::sync::OnceLock;

/// Database of the proteases available to protein digestion.
///
/// Source `ProteaseDB` is a process-wide singleton reached through
/// `getInstance()`, immutable once constructed. This port keeps it immutable
/// and compiled: [`ProteaseDB::global`] replaces `getInstance()`, and every
/// enzyme reference it hands out lives for the whole program.
///
/// Names are indexed as the source indexes them: under the enzyme name, under
/// the ASCII lower-cased name, and under every synonym. When two enzymes claim
/// the same key the later registration wins, exactly as the source's
/// `std::map` assignment does.
///
/// Not ported: the public provider constructor
/// `ProteaseDB(std::vector<std::unique_ptr<DigestionEnzymeDataProvider<DigestionEnzymeProtein>>>)`.
/// It needs owned, user-defined [`DigestionEnzymeProtein`] values with an
/// executable cleavage expression, which the compiled enzyme type does not
/// offer yet; `docs/DIGESTION_SUPPORT.md` records the gap.
#[derive(Clone, Copy, Debug, Default)]
pub struct ProteaseDB;

/// The source's three indices after both provider passes, as positions into
/// the compiled table.
#[derive(Debug)]
struct RegistryIndex {
    /// `enzyme_names_`: name, lower-cased name and synonyms.
    names: BTreeMap<String, usize>,
    /// `enzyme_regex_`: cleavage expression, for non-empty expressions only.
    regex: BTreeMap<&'static str, usize>,
    /// `const_enzymes_`, in registration order (see [`ProteaseDB::names`]).
    live: Vec<usize>,
}

impl RegistryIndex {
    /// Replays `DigestionEnzymeDB::addEnzyme_` for one enzyme.
    ///
    /// The source looks the *new* name up in the name index; whatever it finds
    /// (possibly a different enzyme that owns this string as a synonym or as its
    /// lower-cased name) is removed together with all of that enzyme's keys and
    /// its expression entry, even when the expression key has since been
    /// claimed by a third enzyme. Only then is the new enzyme registered.
    fn add(&mut self, table: &[DigestionEnzymeProtein], new: usize) {
        let enzyme = &table[new];
        if let Some(&old) = self.names.get(enzyme.name) {
            self.live.retain(|&index| index != old);
            let old_enzyme = &table[old];
            self.names.remove(old_enzyme.name);
            self.names.remove(&old_enzyme.name.to_ascii_lowercase());
            for synonym in old_enzyme.synonyms {
                self.names.remove(*synonym);
            }
            if !old_enzyme.regex.is_empty() {
                self.regex.remove(old_enzyme.regex);
            }
        }
        self.live.push(new);
        self.names.insert(enzyme.name.to_owned(), new);
        self.names.insert(enzyme.name.to_ascii_lowercase(), new);
        for synonym in enzyme.synonyms {
            self.names.insert((*synonym).to_owned(), new);
        }
        if !enzyme.regex.is_empty() {
            self.regex.insert(enzyme.regex, new);
        }
    }

    fn build(table: &[DigestionEnzymeProtein]) -> Self {
        let mut index = Self {
            names: BTreeMap::new(),
            regex: BTreeMap::new(),
            live: Vec::with_capacity(table.len()),
        };
        // Pass 1: BuiltInProteaseDataProvider. Pass 2: Enzymes.xml, which at
        // the pin defines the same enzymes in the same order and replaces each
        // built-in entry by name.
        for _pass in 0..2 {
            for position in 0..table.len() {
                index.add(table, position);
            }
        }
        index
    }
}

fn registry_index() -> &'static RegistryIndex {
    static INDEX: OnceLock<RegistryIndex> = OnceLock::new();
    INDEX.get_or_init(|| RegistryIndex::build(&enzymes::ENZYMES))
}

impl ProteaseDB {
    /// The process-wide database, replacing the source's `getInstance()`.
    ///
    /// The source constructs its singleton lazily and never frees it; this
    /// port's database is a compiled table and its lookup index is built once,
    /// thread-safely, on first use.
    pub fn global() -> &'static Self {
        &Self
    }

    /// Every enzyme definition of the compiled table, in the pinned
    /// `Enzymes.xml` order.
    ///
    /// This is the storage behind the source's
    /// `beginEnzyme()`/`endEnzyme()` iterator pair. At the pin every definition
    /// is also registered, so this has the same members as [`Self::names`].
    pub fn enzymes(&self) -> &'static [DigestionEnzymeProtein] {
        &enzymes::ENZYMES
    }

    /// Iterates the registered enzymes, replacing `beginEnzyme()`/`endEnzyme()`.
    ///
    /// Iteration follows registration order. The source iterates a
    /// `std::set<const DigestionEnzymeProtein*>`, which orders by heap address;
    /// see [`Self::names`] for why that order is not reproduced.
    pub fn iter(&self) -> impl Iterator<Item = &'static DigestionEnzymeProtein> + '_ {
        registry_index()
            .live
            .iter()
            .map(|&index| &enzymes::ENZYMES[index])
    }

    /// Whether an enzyme is registered under `name`, `getEnzyme`'s key set:
    /// an enzyme name, its ASCII lower-cased form, or a synonym.
    ///
    /// The lookup is otherwise case-sensitive: `"Trypsin"` and `"trypsin"` are
    /// found, `"TRYPSIN"` is not.
    pub fn has_enzyme(&self, name: &str) -> bool {
        registry_index().names.contains_key(name)
    }

    /// Returns the enzyme registered under `name` (supports synonym names).
    ///
    /// Enzymes are registered under their name, the ASCII lower-cased name and
    /// every synonym, as the source's note says: "if unsure use toLowercase".
    ///
    /// # Errors
    ///
    /// [`crate::Error::InvalidValue`] when no enzyme is registered under
    /// `name`, where the source throws `Exception::ElementNotFound`.
    pub fn get_enzyme(&self, name: &str) -> Result<&'static DigestionEnzymeProtein> {
        registry_index()
            .names
            .get(name)
            .map(|&index| &enzymes::ENZYMES[index])
            .ok_or_else(|| invalid(format!("unknown protease {name:?}")))
    }

    /// Whether `enzyme` is one of this database's own entries, the source's
    /// `hasEnzyme(const DigestionEnzymeProtein*)`.
    ///
    /// Source compares pointers, so an equal copy of a registered enzyme is
    /// *not* contained; this compares addresses with [`std::ptr::eq`] for the
    /// same result. Use [`Self::has_enzyme`] with the enzyme's name to ask
    /// whether an equal definition exists.
    pub fn contains(&self, enzyme: &DigestionEnzymeProtein) -> bool {
        self.iter().any(|entry| std::ptr::eq(entry, enzyme))
    }

    /// Returns all enzyme names, without synonyms (`getAllNames`).
    ///
    /// The order is registration order, which at the pin is the order of
    /// `Enzymes.xml`. The source fills this list by walking a set of *pointers*,
    /// so its order is the heap address order of the enzyme objects: an
    /// allocator artefact, not a specified order, that no caller can rely on.
    /// The members are the same.
    pub fn names(&self) -> Vec<&'static str> {
        self.iter().map(|enzyme| enzyme.name).collect()
    }

    /// Whether an enzyme is registered with exactly this cleavage expression
    /// (`hasRegEx`). The text is compared, not compiled.
    pub fn has_regex(&self, regex: &str) -> bool {
        registry_index().regex.contains_key(regex)
    }

    /// Returns the enzyme registered for a cleavage expression
    /// (`getEnzymeByRegEx`).
    ///
    /// Several enzymes can share an expression; the source's index keeps one
    /// per expression, the last one registered, and this returns that same
    /// enzyme. Use [`Self::enzymes_by_regex`] for every enzyme with the
    /// expression.
    ///
    /// # Errors
    ///
    /// [`crate::Error::InvalidValue`] when no enzyme has the expression, where
    /// the source throws `Exception::IllegalArgument` ("was not registered in
    /// Enzyme DB, register first!"). The source's own `@TODO` asks why this
    /// differs from `getEnzyme`'s `ElementNotFound`; both map to the same Rust
    /// variant here.
    pub fn get_enzyme_by_regex(&self, regex: &str) -> Result<&'static DigestionEnzymeProtein> {
        registry_index()
            .regex
            .get(regex)
            .map(|&index| &enzymes::ENZYMES[index])
            .ok_or_else(|| {
                invalid(format!(
                    "Enzyme with regex {regex} was not registered in Enzyme DB, register first!"
                ))
            })
    }

    /// Every registered enzyme whose cleavage expression is exactly `regex`,
    /// in registration order.
    ///
    /// Native addition: the source can return only the one enzyme its index
    /// keeps ([`Self::get_enzyme_by_regex`]).
    ///
    /// # Errors
    ///
    /// [`crate::Error::Unsupported`] when no registered enzyme has the
    /// expression: this is a lookup of pinned expressions, not a regex compiler.
    pub fn enzymes_by_regex(&self, regex: &str) -> Result<Vec<&'static DigestionEnzymeProtein>> {
        let matches: Vec<_> = self.iter().filter(|e| e.regex == regex).collect();
        if matches.is_empty() {
            Err(crate::Error::Unsupported(
                "only cleavage expressions present in the pinned enzyme registry are supported"
                    .into(),
            ))
        } else {
            Ok(matches)
        }
    }

    /// Names of the enzymes available for X! Tandem, those with a non-empty
    /// X! Tandem identifier (`getAllXTandemNames`), in [`Self::names`] order.
    pub fn xtandem_names(&self) -> Vec<&'static str> {
        self.iter()
            .filter(|e| !e.xtandem_id.is_empty())
            .map(|e| e.name)
            .collect()
    }

    /// Names of the enzymes available for Comet (`getAllCometNames`), in
    /// [`Self::names`] order. The source's `-1` "no identifier" sentinel is
    /// `None` here.
    pub fn comet_names(&self) -> Vec<&'static str> {
        self.iter()
            .filter(|e| e.comet_id.is_some())
            .map(|e| e.name)
            .collect()
    }

    /// Names of the enzymes available for MS-GF+ (`getAllMSGFNames`), in
    /// [`Self::names`] order. MS-GF+ numbers enzymes from 0, so identifier 0
    /// counts; the source's `-1` sentinel is `None` here.
    pub fn msgf_names(&self) -> Vec<&'static str> {
        self.iter()
            .filter(|e| e.msgf_id.is_some())
            .map(|e| e.name)
            .collect()
    }

    /// Names of the enzymes available for OMSSA (`getAllOMSSANames`), in
    /// [`Self::names`] order. The source's `-1` sentinel is `None` here.
    pub fn omssa_names(&self) -> Vec<&'static str> {
        self.iter()
            .filter(|e| e.omssa_id.is_some())
            .map(|e| e.name)
            .collect()
    }

    /// Writes the enzyme names as the source's `writeTSV` does: the header
    /// line `OpenMS_AllowedEnzymes`, then one name per line, each line ending
    /// in `\n`.
    ///
    /// Names follow [`Self::names`] order, where the source's follow heap
    /// addresses; the set of lines is the same.
    ///
    /// # Errors
    ///
    /// [`crate::Error::Io`] when writing fails. The source writes through an
    /// unchecked `std::ofstream` and reports nothing.
    pub fn write_tsv(&self, out: &mut impl Write) -> Result<()> {
        out.write_all(b"OpenMS_AllowedEnzymes\n")?;
        for name in self.names() {
            out.write_all(name.as_bytes())?;
            out.write_all(b"\n")?;
        }
        Ok(())
    }

    /// Writes [`Self::write_tsv`]'s text to a file, replacing the source's
    /// `writeTSV(const std::string& filename)`.
    ///
    /// The text is assembled first and written in one call, so a file that
    /// cannot be created leaves nothing behind.
    ///
    /// # Errors
    ///
    /// [`crate::Error::Io`] when the file cannot be created or written. The
    /// source silently writes nothing when its `std::ofstream` fails to open.
    pub fn write_tsv_file(&self, path: impl AsRef<Path>) -> Result<()> {
        let mut text = Vec::new();
        self.write_tsv(&mut text)?;
        std::fs::write(path, text)?;
        Ok(())
    }
}
