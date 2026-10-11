// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The `CHEMISTRY/ModificationsDB.h` members that answer with the source's own
//! residue-matching rule and ordering: `has`, `searchModificationsFast`,
//! `searchModification`, `findModificationIndex`, both
//! `searchModificationsByDiffMonoMassSorted` overloads,
//! `getAllSearchModifications`, `writeTSV`, both `addModification` overloads and
//! `isInstantiated` (core `bc9cc12`).
//!
//! Every expected value in `tests/modifications_db.rs` comes from the executed
//! C++ oracle `../oracle/B8b-dbs`, which builds `ModificationsDB` from the same
//! three providers as [`ModificationsDB::global`]. See
//! `docs/MODIFICATION_REGISTRY_SUPPORT.md`.

use super::{
    GLOBAL_REGISTRY, ModificationsDB, OboReadOptions, ResidueModification, TermSpecificity,
};
use crate::concept::fuzzy_string_comparator::format_g;
use crate::{Error, Result};
use std::cmp::Ordering;
use std::io::Write;
use std::path::Path;
use std::sync::Arc;

/// The source's `residuesMatch_`, with `None` standing for its `'?'` ("no
/// residue given").
///
/// A modification whose origin is not `X` matches its own residue and the
/// wildcards `X`, `.` and `?`. A modification with origin `X` (the terminal
/// records, whose origin this port stores as `None`) matches any residue,
/// except that a *user-defined* one (empty short ID, non-empty full ID) then
/// means a literal `X` and matches only `X` or no residue.
fn residues_match(residue: Option<char>, modification: &ResidueModification) -> bool {
    let residue = residue.unwrap_or('?');
    let origin = modification.origin.unwrap_or('X');
    if origin != 'X' {
        origin == residue || matches!(residue, 'X' | '.' | '?')
    } else {
        let user_defined = modification.name.is_empty() && !modification.full_id.is_empty();
        !(user_defined && residue != '?' && origin != residue)
    }
}

/// The source's residue and terminal-specificity filter; `term == None` is its
/// `NUMBER_OF_TERM_SPECIFICITY` ("any").
fn source_filter(
    modification: &ResidueModification,
    residue: Option<char>,
    term: Option<TermSpecificity>,
) -> bool {
    residues_match(residue, modification) && term.is_none_or(|t| t == modification.term)
}

/// The source's `getTermSpecificityName`, which names `ANYWHERE` `none`
/// (unlike [`TermSpecificity::name`], which serves full-ID construction).
fn term_specificity_name(term: TermSpecificity) -> &'static str {
    match term {
        TermSpecificity::Anywhere => "none",
        TermSpecificity::NTerm => "N-term",
        TermSpecificity::CTerm => "C-term",
        TermSpecificity::ProteinNTerm => "Protein N-term",
        TermSpecificity::ProteinCTerm => "Protein C-term",
    }
}

/// `getAllSearchModifications`' comparator: byte by byte after C `tolower`,
/// the shorter string first when one is a prefix of the other.
///
/// The source compares `tolower(char)` values, and `char` is signed on the
/// x86-64 Linux build, so a byte from 0x80 up sorts *before* ASCII; this
/// compares the same signed values. It is a total order, which the source's
/// `std::sort` requires and Rust's sort relies on.
fn search_name_order(left: &str, right: &str) -> Ordering {
    let key = |byte: u8| byte.to_ascii_lowercase() as i8;
    left.bytes()
        .map(key)
        .cmp(right.bytes().map(key))
        .then_with(|| left.len().cmp(&right.len()))
}

impl ModificationsDB {
    /// Whether the process-wide registry was built (`isInstantiated`).
    ///
    /// Source `initializeModificationsDB` can load other files before the
    /// singleton exists and `isInstantiated` tells a caller whether it is too
    /// late; [`ModificationsDB::global`] takes no files, so this only reports
    /// whether it has been called.
    pub fn is_instantiated() -> bool {
        GLOBAL_REGISTRY.get().is_some()
    }

    /// Whether `name` is a lookup key of any record (`has`): a full ID, short
    /// ID, full name, UniMod or OBO accession, or synonym, matched exactly.
    ///
    /// Unlike [`Self::find`], the UniMod prefix is not normalised.
    pub fn has(&self, name: &str) -> bool {
        self.by_name.contains_key(name)
    }

    /// Resolves `name` in the name index, retrying a lower-case `unimod`
    /// prefix as `UniMod` as the source's searches do.
    fn source_name_indices(&self, name: &str) -> &[usize] {
        if let Some(indices) = self.by_name.get(name) {
            return indices;
        }
        // Source: `mod_name.size() > 6 && toLowered(prefix(mod_name, 6)) ==
        // "unimod"`, then "UniMod" + the rest (Skyline writes `unimod:10`).
        let bytes = name.as_bytes();
        if bytes.len() > 6 && bytes[..6].eq_ignore_ascii_case(b"unimod") {
            let fixed = format!("UniMod{}", &name[6..]);
            if let Some(indices) = self.by_name.get(&fixed) {
                return indices;
            }
        }
        &[]
    }

    /// Returns the modification registered under `name` that matches the
    /// residue and terminal filters, with a flag set when more than one does
    /// (`searchModificationsFast`).
    ///
    /// `residue: None` and `term: None` apply no filter. Residues match with
    /// the source's rule: a non-`X` origin matches its residue and the
    /// wildcards `X`, `.` and `?`; an `X` or terminal origin matches any residue
    /// unless the record is user-defined. A lower-case `unimod` prefix is
    /// retried as `UniMod`.
    ///
    /// Of several matches the source returns the *last* in its name set, which
    /// is ordered by heap address; this returns the last in registry order. The
    /// oracle found the two orders to differ for 6 of 7,806 keys on the Linux
    /// Release build, where the source's pick is an allocator artefact.
    ///
    /// Source logs a warning for an unknown name; this returns `(None, false)`.
    pub fn search_modifications_fast(
        &self,
        name: &str,
        residue: Option<char>,
        term: Option<TermSpecificity>,
    ) -> (Option<&ResidueModification>, bool) {
        let mut found = None;
        let mut count = 0usize;
        for &index in self.source_name_indices(name) {
            let modification = self.entries[index].as_ref();
            if source_filter(modification, residue, term) {
                found = Some(modification);
                count += 1;
            }
        }
        (found, count > 1)
    }

    /// Returns the stored record equal to `modification`, looked up by its
    /// full ID (`searchModification`).
    ///
    /// Equality is [`ResidueModification`]'s complete value equality, the
    /// counterpart of the source's `operator==`. Of several equal records the
    /// first in registry order is returned; the source walks its address-ordered
    /// name set.
    pub fn search_modification(
        &self,
        modification: &ResidueModification,
    ) -> Option<&ResidueModification> {
        self.by_name
            .get(modification.full_id())
            .into_iter()
            .flatten()
            .map(|&index| self.entries[index].as_ref())
            .find(|stored| *stored == modification)
    }

    /// Index of the one record registered under `name`
    /// (`findModificationIndex`); `entries()[index]` is that record.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] with the source's message when no record has
    /// the key ("Modification not found: ...") or more than one does ("More
    /// than one modification with name: ..."), where the source throws
    /// `Exception::ElementNotFound`. The source's third error, a name whose
    /// record is missing from storage, cannot occur here.
    pub fn find_modification_index(&self, name: &str) -> Result<usize> {
        match self.by_name.get(name).map(Vec::as_slice) {
            None | Some([]) => Err(Error::InvalidValue(format!(
                "Modification not found: {name}"
            ))),
            Some([index]) => Ok(*index),
            Some(_) => Err(Error::InvalidValue(format!(
                "More than one modification with name: {name}"
            ))),
        }
    }

    /// Every record whose monoisotopic mass difference lies within
    /// `max_error` of `mass`, closest first
    /// (`searchModificationsByDiffMonoMassSorted`).
    ///
    /// Equal errors keep registry order, as the source's `(error, count)` map
    /// key does. The source's `std::vector<std::string>` overload returns each
    /// record's [`ResidueModification::full_id`]. A record with a zero mass
    /// difference is only considered when `mass` is exactly zero. Residue and
    /// terminal filters use the source's rule, as
    /// [`Self::search_modifications_fast`] describes.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `mass` or `max_error` is not finite or
    /// `max_error` is negative; the source returns an empty list for a
    /// negative or NaN tolerance.
    pub fn search_by_mass_sorted(
        &self,
        mass: f64,
        max_error: f64,
        residue: Option<char>,
        term: Option<TermSpecificity>,
    ) -> Result<Vec<&ResidueModification>> {
        if !mass.is_finite() || !max_error.is_finite() || max_error < 0.0 {
            return Err(Error::InvalidValue(
                "mass and nonnegative tolerance must be finite".into(),
            ));
        }
        let mut hits: Vec<(f64, &ResidueModification)> = self
            .entries
            .iter()
            .map(Arc::as_ref)
            .filter_map(|modification| {
                let error = (modification.diff_mono_mass - mass).abs();
                (error <= max_error
                    && (modification.diff_mono_mass != 0.0 || mass == 0.0)
                    && source_filter(modification, residue, term))
                .then_some((error, modification))
            })
            .collect();
        // Stable: equal errors stay in registry order.
        hits.sort_by(|left, right| left.0.total_cmp(&right.0));
        Ok(hits
            .into_iter()
            .map(|(_, modification)| modification)
            .collect())
    }

    /// Full IDs of every record that can be used in identification searches,
    /// those with a UniMod record ID, sorted case-insensitively
    /// (`getAllSearchModifications`).
    ///
    /// Of two IDs that agree ignoring case, the shorter comes first. IDs equal
    /// ignoring case and length keep registry order; the source's unstable
    /// `std::sort` leaves their order unspecified.
    pub fn all_search_modifications(&self) -> Vec<String> {
        let mut names: Vec<String> = self
            .entries
            .iter()
            .filter(|modification| modification.record_id.is_some_and(|id| id > 0))
            .map(|modification| modification.full_id.clone())
            .collect();
        names.sort_by(|left, right| search_name_order(left, right));
        names
    }

    /// Writes one tab-separated line per record with a header, as the source's
    /// `writeTSV`: `FullId`, `FullName`, `UnimodAccession`, `Origin/AA`,
    /// `TerminusSpecificity` and `DiffMonoMass`.
    ///
    /// The origin of a terminal record is written `X`, the source's stored
    /// origin; anywhere specificity is written `none`, the source's
    /// `getTermSpecificityName`; the mass is written as `std::ostream`'s default
    /// six significant digits (`42.0106`). The executed C++ oracle's file is
    /// reproduced byte for byte for the global registry.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when writing fails; the source's unchecked
    /// `std::ofstream` reports nothing.
    pub fn write_tsv(&self, out: &mut impl Write) -> Result<()> {
        out.write_all(
            b"FullId\tFullName\tUnimodAccession\tOrigin/AA\tTerminusSpecificity\tDiffMonoMass\n",
        )?;
        for modification in &self.entries {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}",
                modification.full_id,
                modification.full_name,
                modification.unimod_accession().unwrap_or_default(),
                modification.origin.unwrap_or('X'),
                term_specificity_name(modification.term),
                format_g(modification.diff_mono_mass),
            )?;
        }
        Ok(())
    }

    /// Writes [`Self::write_tsv`]'s text to a file, replacing
    /// `writeTSV(const std::string& filename)`.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be created or written; the source
    /// silently writes nothing when its stream fails to open.
    pub fn write_tsv_file(&self, path: impl AsRef<Path>) -> Result<()> {
        let mut text = Vec::new();
        self.write_tsv(&mut text)?;
        std::fs::write(path, text)?;
        Ok(())
    }

    /// Adds a modification unless one is already registered under its full ID,
    /// and returns the stored record (both `addModification` overloads).
    ///
    /// Source semantics are kept: the duplicate test looks the *full ID* up in
    /// the whole name index, so any record that has that string as any key
    /// counts, and the first such record in registry order is returned (the
    /// source takes the lowest address). A new record is appended and indexed
    /// under its full ID, short ID, full name and UniMod accession only, the
    /// source's four keys; an absent short ID, full name or accession is indexed
    /// as the empty string, as in the source, so [`Self::has`]`("")` becomes
    /// true. Synonyms and the OBO accession are not indexed by this path.
    ///
    /// The source mutates its `const` singleton; here the registry is a
    /// caller-owned value, and [`ModificationsDB::global`] stays immutable.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the record would exceed the default
    /// [`OboReadOptions`] record or payload ceiling; the registry is then
    /// unchanged.
    pub fn add_modification(
        &mut self,
        modification: ResidueModification,
    ) -> Result<Arc<ResidueModification>> {
        if let Some(&index) = self
            .by_name
            .get(modification.full_id())
            .and_then(|indices| indices.first())
        {
            return Ok(Arc::clone(&self.entries[index]));
        }
        let limits = OboReadOptions::default();
        if self.entries.len() >= limits.max_records {
            return Err(Error::InvalidValue(
                "modification registry record limit exceeded".into(),
            ));
        }
        if modification.payload_bytes()? > limits.max_registry_bytes {
            return Err(Error::InvalidValue(
                "modification definition payload exceeds registry limit".into(),
            ));
        }
        let index = self.entries.len();
        let keys = [
            modification.full_id.clone(),
            modification.name.clone(),
            modification.full_name.clone(),
            modification.unimod_accession().unwrap_or_default(),
        ];
        let stored = Arc::new(modification);
        self.entries.push(Arc::clone(&stored));
        for key in keys {
            let indices = self.by_name.entry(key).or_default();
            // The same key twice (e.g. an empty ID and name) is one set entry
            // in the source; indices are appended in increasing order.
            if indices.last() != Some(&index) {
                indices.push(index);
            }
        }
        Ok(stored)
    }
}
