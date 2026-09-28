// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Inspect input file: the comma-separated parameter file Inspect reads, e.g.
//!
//! ```text
//! spectra,dummy4712
//! db,dummy4711
//! protease,Trypsin
//! mod,+57,C,fix,carbamidomethylation
//! PM_tolerance,1.3
//! ```
//!
//! Ports `FORMAT/InspectInfile.h`/`.cpp`, including its modification-line
//! parser [`InspectInfile::handle_ptms`](crate::format::inspect_infile::InspectInfile::handle_ptms), which shares its grammar with
//! `SequestInfile::handlePTMs`. See `docs/SEARCH_ENGINE_FILES_SUPPORT.md`
//! for the API mapping.

use crate::chemistry::EmpiricalFormula;
use crate::data_structures::string_utils::{split, to_float};
use crate::format::file_info::text_format::{ostream_g, to_str};
use crate::format::file_types::{FileType, has_valid_extension};
use crate::format::ptm_xml_file::{PTMXMLFile, PtmInformations};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

/// Inspect input file.
///
/// Unset numeric settings are the source's sentinels: `-1` for the counts
/// and tolerances, `2` for the `blind` and `multicharge` switches (0 false,
/// 1 true, 2 not set). Unset settings are not written.
///
/// Equality is the source's `operator==`, which compares every setting
/// except the database ([`InspectInfile::db`]). A clone copies every
/// setting; the source's copy constructor and assignment operator drop the
/// database.
#[derive(Clone, Debug)]
pub struct InspectInfile {
    spectra: String,
    db: String,
    enzyme: String,
    modifications_per_peptide: i32,
    blind: u32,
    maxptmsize: f32,
    precursor_mass_tolerance: f32,
    peak_mass_tolerance: f32,
    multicharge: u32,
    instrument: String,
    tag_count: i32,
    ptm_name_residues_mass_type: BTreeMap<String, Vec<String>>,
}

impl Default for InspectInfile {
    fn default() -> Self {
        Self {
            spectra: String::new(),
            db: String::new(),
            enzyme: String::new(),
            modifications_per_peptide: -1,
            blind: 2,
            maxptmsize: -1.0,
            precursor_mass_tolerance: -1.0,
            peak_mass_tolerance: -1.0,
            multicharge: 2,
            instrument: String::new(),
            tag_count: -1,
            ptm_name_residues_mass_type: BTreeMap::new(),
        }
    }
}

impl PartialEq for InspectInfile {
    fn eq(&self, other: &Self) -> bool {
        self.spectra == other.spectra
            && self.enzyme == other.enzyme
            && self.modifications_per_peptide == other.modifications_per_peptide
            && self.blind == other.blind
            && self.maxptmsize == other.maxptmsize
            && self.precursor_mass_tolerance == other.precursor_mass_tolerance
            && self.peak_mass_tolerance == other.peak_mass_tolerance
            && self.multicharge == other.multicharge
            && self.instrument == other.instrument
            && self.tag_count == other.tag_count
            && self.ptm_name_residues_mass_type == other.ptm_name_residues_mass_type
    }
}

impl InspectInfile {
    /// An input file with every setting unset.
    pub fn new() -> Self {
        Self::default()
    }

    /// Store the parameters.
    ///
    /// Writes `spectra`, then each set setting in the source order: `db`,
    /// `protease`, `blind`, one `mod,<mass>,<residues>,<type>,<name>` line per
    /// modification in name order, `mods`, `maxptmsize`, `PM_tolerance`,
    /// `IonTolerance`, `multicharge`, `instrument`, `TagCount`. Floats are
    /// written as a C++ stream writes them (six significant digits).
    ///
    /// As in the source, storing lower-cases the stored modification types
    /// and completes `nterm`/`cterm` to `nterminal`/`cterminal`, so
    /// [`InspectInfile::modifications`] reflects that afterwards.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] "invalid file extension, expected 'tsv'" when
    /// the name has a known extension other than `.tsv` (an unknown one,
    /// e.g. `.tmp`, is accepted, as the source's
    /// `FileHandler::hasValidExtension` does), and [`Error::Io`] when the
    /// file cannot be created or written; both are the source's
    /// `Exception::UnableToCreateFile`.
    pub fn store(&mut self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        if !has_valid_extension(&path.to_string_lossy(), FileType::Tsv) {
            return Err(Error::InvalidValue(format!(
                "the file '{}' could not be created: invalid file extension, expected 'tsv'",
                path.display()
            )));
        }
        let mut file = std::fs::File::create(path)?;
        let text = self.render();
        file.write_all(text.as_bytes())?;
        file.flush()?;
        Ok(())
    }

    /// Write the parameters to `writer`; see [`InspectInfile::store`],
    /// without the extension check.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when writing fails.
    pub fn write(&mut self, mut writer: impl Write) -> Result<()> {
        let text = self.render();
        writer.write_all(text.as_bytes())?;
        Ok(())
    }

    fn render(&mut self) -> String {
        let mut out = format!("spectra,{}\n", self.spectra);
        if !self.db.is_empty() {
            out.push_str(&format!("db,{}\n", self.db));
        }
        if !self.enzyme.is_empty() {
            out.push_str(&format!("protease,{}\n", self.enzyme));
        }
        if self.blind != 2 {
            out.push_str(&format!("blind,{}\n", self.blind));
        }
        for (name, entry) in &mut self.ptm_name_residues_mass_type {
            // Source: toLower, then "term" -> "terminal" (fix, opt,
            // cterminal, nterminal).
            entry[2].make_ascii_lowercase();
            if entry[2].ends_with("term") {
                entry[2].push_str("inal");
            }
            out.push_str(&format!(
                "mod,{},{},{},{name}\n",
                entry[1], entry[0], entry[2]
            ));
        }
        if self.modifications_per_peptide > -1 {
            out.push_str(&format!("mods,{}\n", self.modifications_per_peptide));
        }
        if self.maxptmsize >= 0.0 {
            out.push_str(&format!("maxptmsize,{}\n", stream_float(self.maxptmsize)));
        }
        if self.precursor_mass_tolerance >= 0.0 {
            out.push_str(&format!(
                "PM_tolerance,{}\n",
                stream_float(self.precursor_mass_tolerance)
            ));
        }
        if self.peak_mass_tolerance >= 0.0 {
            out.push_str(&format!(
                "IonTolerance,{}\n",
                stream_float(self.peak_mass_tolerance)
            ));
        }
        if self.multicharge != 2 {
            out.push_str(&format!("multicharge,{}\n", self.multicharge));
        }
        if !self.instrument.is_empty() {
            out.push_str(&format!("instrument,{}\n", self.instrument));
        }
        if self.tag_count >= 0 {
            out.push_str(&format!("TagCount,{}\n", self.tag_count));
        }
        out
    }

    /// Parse a modification line into [`InspectInfile::modifications`];
    /// see [`parse_modification_line`] for the grammar. `modifications_file`
    /// is a PTMXML table (e.g. `Inspect_PTMs.xml`), read only when a
    /// modification is not given by a mass.
    ///
    /// # Errors
    ///
    /// As [`parse_modification_line`]. The modifications are cleared
    /// first; a parse error leaves them empty, while a missing or
    /// unreadable modifications file leaves those parsed before it, as in
    /// the source.
    pub fn handle_ptms(
        &mut self,
        modification_line: &str,
        modifications_file: impl AsRef<Path>,
        monoisotopic: bool,
    ) -> Result<()> {
        self.ptm_name_residues_mass_type.clear();
        parse_modification_line(
            modification_line,
            modifications_file.as_ref(),
            monoisotopic,
            &mut self.ptm_name_residues_mass_type,
        )
    }

    /// The modifications: name to `[residues, mass, type]`.
    pub fn modifications(&self) -> &BTreeMap<String, Vec<String>> {
        &self.ptm_name_residues_mass_type
    }

    /// The spectrum file (or directory) to search.
    pub fn spectra(&self) -> &str {
        &self.spectra
    }

    /// Set the spectrum file (or directory) to search.
    pub fn set_spectra(&mut self, spectra: &str) {
        self.spectra = spectra.to_owned();
    }

    /// The `.trie` database to search.
    pub fn db(&self) -> &str {
        &self.db
    }

    /// Set the `.trie` database to search.
    pub fn set_db(&mut self, db: &str) {
        self.db = db.to_owned();
    }

    /// The enzyme (`Trypsin`, `None` or `Chymotrypsin`).
    pub fn enzyme(&self) -> &str {
        &self.enzyme
    }

    /// Set the enzyme.
    pub fn set_enzyme(&mut self, enzyme: &str) {
        self.enzyme = enzyme.to_owned();
    }

    /// The number of modifications permitted in one peptide (-1 unset).
    pub fn modifications_per_peptide(&self) -> i32 {
        self.modifications_per_peptide
    }

    /// Set the number of modifications permitted in one peptide.
    pub fn set_modifications_per_peptide(&mut self, modifications_per_peptide: i32) {
        self.modifications_per_peptide = modifications_per_peptide;
    }

    /// Blind (MS-Alignment) search: 0 no, 1 yes, 2 unset.
    pub fn blind(&self) -> u32 {
        self.blind
    }

    /// Set the blind search switch.
    pub fn set_blind(&mut self, blind: u32) {
        self.blind = blind;
    }

    /// The largest modification (Da) a blind search considers (negative unset).
    pub fn max_ptm_size(&self) -> f32 {
        self.maxptmsize
    }

    /// Set the largest modification a blind search considers (source
    /// `setMaxPTMsize`).
    pub fn set_max_ptm_size(&mut self, maxptmsize: f32) {
        self.maxptmsize = maxptmsize;
    }

    /// The parent mass tolerance in Daltons (negative unset).
    pub fn precursor_mass_tolerance(&self) -> f32 {
        self.precursor_mass_tolerance
    }

    /// Set the parent mass tolerance in Daltons.
    pub fn set_precursor_mass_tolerance(&mut self, precursor_mass_tolerance: f32) {
        self.precursor_mass_tolerance = precursor_mass_tolerance;
    }

    /// How far b and y peaks may shift from their expected masses
    /// (negative unset).
    pub fn peak_mass_tolerance(&self) -> f32 {
        self.peak_mass_tolerance
    }

    /// Set how far b and y peaks may shift.
    pub fn set_peak_mass_tolerance(&mut self, peak_mass_tolerance: f32) {
        self.peak_mass_tolerance = peak_mass_tolerance;
    }

    /// Guess the precursor charge: 0 no, 1 yes, 2 unset.
    pub fn multicharge(&self) -> u32 {
        self.multicharge
    }

    /// Set the multicharge switch.
    pub fn set_multicharge(&mut self, multicharge: u32) {
        self.multicharge = multicharge;
    }

    /// The instrument (e.g. `QTOF`, `ESI-ION-TRAP`).
    pub fn instrument(&self) -> &str {
        &self.instrument
    }

    /// Set the instrument.
    pub fn set_instrument(&mut self, instrument: &str) {
        self.instrument = instrument.to_owned();
    }

    /// The number of tags to generate (-1 unset).
    pub fn tag_count(&self) -> i32 {
        self.tag_count
    }

    /// Set the number of tags to generate.
    pub fn set_tag_count(&mut self, tag_count: i32) {
        self.tag_count = tag_count;
    }
}

/// `std::ostream << float` in the default float field, precision 6.
fn stream_float(value: f32) -> String {
    ostream_g(f64::from(value), 6)
}

fn modification_error(message: &str, modification: &str) -> Error {
    Error::Parse {
        line: 0,
        message: format!("{message} in: {modification}"),
    }
}

/// Source `EmpiricalFormula(const String&)`: an empty text is the empty
/// formula.
fn formula(text: &str) -> Result<EmpiricalFormula> {
    if text.is_empty() {
        Ok(EmpiricalFormula::default())
    } else {
        EmpiricalFormula::parse(text)
    }
}

/// The modification-line parser of `InspectInfile::handlePTMs` (and, with
/// the same grammar, `SequestInfile::handlePTMs`).
///
/// The line lists modifications separated by `:`, each
/// `<mass|composition|name>[,<residues>][,<type>][,<name>]`:
///
/// - a *mass* is a number with an optional leading `+`, trailing `+`, or
///   trailing `-` (a negative mass), kept as written without the signs;
/// - a *name* is looked up in the PTMXML table `modifications_file`, which
///   supplies the composition and the residues (so no residues follow);
/// - anything else is a *composition*, `<formula>` or
///   `<added formula>-<removed formula>`, whose mono- or average mass
///   difference becomes the mass.
///
/// Residues are upper-cased with `*` read as `X`; the type is `FIX` or
/// `OPT` (default `OPT`; any substring of `OPT#FIX#` is accepted, as in
/// the source); an unnamed modification is called `PTM_<n>`, `n` the number
/// of modifications before it. Each entry maps the name to
/// `[residues, mass, type]`, the mass cut after five decimals, or to its
/// first five characters when it has no decimal point (the source's
/// `substr(0, find(".") + 6)` with `npos + 6` wrapping to 5).
///
/// # Errors
///
/// [`Error::Io`] when the modifications file is needed but missing (source
/// `Exception::FileNotFound`, "the file '…' could not be found") or cannot
/// be read, the PTMXML reader's errors, and [`Error::Parse`] with the
/// source's messages, followed by ` in: <modification>`: "There's something
/// wrong with this modification. Aborting!" for a composition that does not
/// parse, "No residues for modification given. Aborting!", "There's
/// something wrong with the type of this modification. Aborting!" for
/// surplus fields, and "There's already a modification with this name.
/// Aborting!". A parse error clears `out`.
pub fn parse_modification_line(
    modification_line: &str,
    modifications_file: &Path,
    monoisotopic: bool,
    out: &mut BTreeMap<String, Vec<String>>,
) -> Result<()> {
    let result = parse_modifications(modification_line, modifications_file, monoisotopic, out);
    if let Err(Error::Parse { .. }) = &result {
        out.clear();
    }
    result
}

fn parse_modifications(
    modification_line: &str,
    modifications_file: &Path,
    monoisotopic: bool,
    out: &mut BTreeMap<String, Vec<String>>,
) -> Result<()> {
    let mut ptm_informations = PtmInformations::new();
    for modification in split(modification_line, ':') {
        if modification.is_empty() {
            continue;
        }
        let mut parts: Vec<&str> = split(modification, ',');
        let first = parts[0];

        // 0 mass, 1 composition, 2 name
        let mut kind: i32 = -1;
        let mut residues = String::new();
        let mut name = String::new();

        let mut mass = first.to_owned();
        if mass.starts_with('+') {
            mass.remove(0);
        }
        if mass.ends_with('+') {
            mass.pop();
        }
        if mass.ends_with('-') {
            mass.pop();
            mass.insert(0, '-');
        }
        if to_float(&mass).is_ok() {
            kind = 0;
        }

        if kind == -1 {
            if ptm_informations.is_empty() {
                if !modifications_file.exists() {
                    return Err(Error::Io(std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        format!(
                            "the file '{}' could not be found",
                            modifications_file.display()
                        ),
                    )));
                }
                ptm_informations = PTMXMLFile::new().load(modifications_file)?;
            }
            if let Some((composition, possible)) = ptm_informations.get(first) {
                mass = composition.clone();
                residues = possible.clone();
                name = first.to_owned();
                kind = 2;
            }
        }

        if kind == -1 {
            mass = first.to_owned();
        }
        if kind == -1 || kind == 2 {
            let formulas = match mass.find('-') {
                Some(pos) => formula(&mass[..pos])
                    .and_then(|add| formula(&mass[pos + 1..]).map(|sub| (add, sub))),
                None => formula(&mass).map(|add| (add, EmpiricalFormula::default())),
            };
            let (add, sub) = formulas.map_err(|_| {
                modification_error(
                    "There's something wrong with this modification. Aborting!",
                    modification,
                )
            })?;
            let difference = if monoisotopic {
                add.mono_mass() - sub.mono_mass()
            } else {
                add.average_mass() - sub.average_mass()
            };
            mass = to_str(difference);
            if kind == -1 {
                kind = 1;
            }
        }

        parts.remove(0);
        if kind < 2 {
            if parts.is_empty() {
                return Err(modification_error(
                    "No residues for modification given. Aborting!",
                    modification,
                ));
            }
            residues = parts.remove(0).replace('*', "X").to_ascii_uppercase();
        }

        let mut kind_of_mod = "OPT".to_owned();
        if let Some(candidate) = parts.first() {
            let candidate = candidate.to_ascii_uppercase();
            if "OPT#FIX#".contains(candidate.as_str()) {
                parts.remove(0);
                kind_of_mod = candidate;
            }
        }
        if parts.len() > 1 {
            return Err(modification_error(
                "There's something wrong with the type of this modification. Aborting!",
                modification,
            ));
        }

        if kind < 2 {
            name = match parts.first() {
                Some(given) => (*given).to_owned(),
                None => format!("PTM_{}", out.len()),
            };
        }

        if out.contains_key(&name) {
            return Err(modification_error(
                "There's already a modification with this name. Aborting!",
                modification,
            ));
        }
        out.insert(name, vec![residues, truncate_mass(&mass), kind_of_mod]);
    }
    Ok(())
}

/// Source `StringUtils::substr(mass, 0, mass.find(".") + 6)`.
fn truncate_mass(mass: &str) -> String {
    let length = match mass.find('.') {
        Some(dot) => dot + 6,
        // npos + 6 wraps around to 5.
        None => 5,
    };
    mass.chars().take(length).collect()
}

#[cfg(test)]
mod tests {
    use super::truncate_mass;

    #[test]
    fn masses_are_cut_as_the_source_cuts_them() {
        assert_eq!(truncate_mass("79.966330999999997"), "79.96633");
        assert_eq!(truncate_mass("10.3"), "10.3");
        assert_eq!(truncate_mass("-16"), "-16");
        assert_eq!(truncate_mass("123456"), "12345");
    }
}
