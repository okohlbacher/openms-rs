// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! PepNovo input file adapter: writes the `PepNovo_PTMs.txt` table of
//! modifications for a PepNovo search.
//!
//! Ports `FORMAT/PepNovoInfile.h`/`.cpp`. See
//! `docs/SEARCH_ENGINE_FILES_SUPPORT.md` for the API mapping.

use crate::chemistry::{ModificationDefinitionsSet, ModificationsDB, TermSpecificity};
use crate::param::value::format_float;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

/// The header line of the PTM table.
pub const HEADER: &str = "#AA\toffset\ttype\tlocations\tsymbol\tPTM\tname";

/// PepNovo input file adapter.
///
/// Two adapters compare equal when their modification definitions do, as in
/// the source's `operator==`; the generated table and key map follow from
/// them.
#[derive(Clone, Debug, Default)]
pub struct PepNovoInfile {
    mods: ModificationDefinitionsSet,
    mods_and_keys: BTreeMap<String, String>,
    ptm_lines: Vec<String>,
}

impl PartialEq for PepNovoInfile {
    fn eq(&self, other: &Self) -> bool {
        self.mods == other.mods
    }
}

impl PepNovoInfile {
    /// An adapter without modifications.
    pub fn new() -> Self {
        Self::default()
    }

    /// Generate the PTM table for the given fixed and variable modifications,
    /// named by their unique identifiers (e.g. `Phospho (S)`).
    ///
    /// The table has the [`HEADER`] line and then one line per fixed and then
    /// per variable modification, each group in name order:
    /// `ORIGIN offset FIXED|OPTIONAL ALL|N_TERM|C_TERM key full-name`, tab
    /// separated. The origin is upper case, or `N_TERM`/`C_TERM` for a
    /// terminal modification of any residue; the offset is the monoisotopic
    /// mass difference as the source's `std::string += double`
    /// (`StringUtils::appendToStr`) writes it, e.g. `79.966330999999997`;
    /// the key is the origin (or
    /// `^`/`$` for any N-/C-terminal residue), a `+` for a non-negative mass,
    /// and the mass rounded half away from zero, e.g. `C+80` or `$+43`.
    ///
    /// The source appends to its table on every call, so a second call
    /// repeats the header; the port replaces the table, as it does the key
    /// map.
    ///
    /// # Errors
    ///
    /// The lookup's error for a name not in [`ModificationsDB::global`], and
    /// [`Error::InvalidValue`] for a protein-terminal modification (the
    /// source's `Exception::InvalidValue` "Invalid term specificity"). On
    /// error the adapter is unchanged.
    pub fn set_modifications(&mut self, fixed: &[&str], variable: &[&str]) -> Result<()> {
        let mods = ModificationDefinitionsSet::from_names(fixed, variable)?;
        let mut mods_and_keys = BTreeMap::new();
        let mut ptm_lines = vec![HEADER.to_owned()];
        for name in mods.fixed_names() {
            ptm_lines.push(ptm_line(&name, false, &mut mods_and_keys)?);
        }
        for name in mods.variable_names() {
            ptm_lines.push(ptm_line(&name, true, &mut mods_and_keys)?);
        }
        self.mods = mods;
        self.mods_and_keys = mods_and_keys;
        self.ptm_lines = ptm_lines;
        Ok(())
    }

    /// The PepNovo keys of the modifications, mapped to their unique
    /// identifiers, e.g. `C+80` to `Phospho (C)` (source `getModifications`).
    pub fn modifications(&self) -> &BTreeMap<String, String> {
        &self.mods_and_keys
    }

    /// The modification definitions last set.
    pub fn modification_definitions(&self) -> &ModificationDefinitionsSet {
        &self.mods
    }

    /// The PTM table, one line per entry.
    pub fn lines(&self) -> &[String] {
        &self.ptm_lines
    }

    /// Write the PTM table, each line followed by a line feed.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when writing fails.
    pub fn write(&self, mut writer: impl Write) -> Result<()> {
        for line in &self.ptm_lines {
            writer.write_all(line.as_bytes())?;
            writer.write_all(b"\n")?;
        }
        Ok(())
    }

    /// Store the PTM table (source `store`).
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be created or written (source
    /// `Exception::UnableToCreateFile`).
    pub fn store(&self, path: impl AsRef<Path>) -> Result<()> {
        let mut buffer = Vec::new();
        self.write(&mut buffer)?;
        std::fs::write(path, buffer)?;
        Ok(())
    }
}

/// Source `handlePTMs_`: one table line, recording its key.
fn ptm_line(
    name: &str,
    variable: bool,
    mods_and_keys: &mut BTreeMap<String, String>,
) -> Result<String> {
    let modification = ModificationsDB::global().get_modification(name, None, None)?;
    let term = modification.term_specificity();
    let mass = modification.diff_mono_mass();
    let kind = if variable { "OPTIONAL" } else { "FIXED" };
    let locations = match term {
        TermSpecificity::CTerm => "C_TERM",
        TermSpecificity::NTerm => "N_TERM",
        TermSpecificity::Anywhere => "ALL",
        other => {
            return Err(Error::InvalidValue(format!(
                "Invalid term specificity {other:?} of '{name}'"
            )));
        }
    };
    let origin_char = modification.origin().unwrap_or('X');
    let mut key = match term {
        TermSpecificity::CTerm => "$".to_owned(),
        TermSpecificity::NTerm => "^".to_owned(),
        _ => String::new(),
    };
    let origin = if term == TermSpecificity::CTerm && origin_char == 'X' {
        "C_TERM".to_owned()
    } else if term == TermSpecificity::NTerm && origin_char == 'X' {
        "N_TERM".to_owned()
    } else {
        key = origin_char.to_string();
        origin_char.to_string()
    };
    if mass >= 0.0 {
        key.push('+');
    }
    key.push_str(&(mass.round() as i32).to_string());
    let line = format!(
        "{}\t{}\t{kind}\t{locations}\t{key}\t{}",
        origin.to_ascii_uppercase(),
        format_float(mass, true),
        modification.full_name()
    );
    mods_and_keys.insert(key, modification.full_id().to_owned());
    Ok(line)
}
